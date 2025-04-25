module Haplogroupers

# haplogrouper is a julia programme to identify clade-defining variants in phylogenetic trees. 

#module TreeVariants

using Phylo # For Newick tree parsing
using BioSequences
using FASTX
using ArgCheck

# TO-DO: Programme should take as input a msa with a user-specified reference. We make a custom VCF file 
# (everything will be treated as homozygous irrespective of ploidy because we have just one sequence) 
# So the reference is 0/0 and any alternative is 1/1, 2/2, 3/3 etc depending on how many alt alleles we have.

#= """
    Haplogrouper

A struct for analysing phylogenetic trees and identifying 
branch-specific genetic variants.
""" =#
mutable struct Haplogrouper
    tree::Union{RecursiveTree, Nothing}
    leaf_vars::Dict{String, Set{Tuple{Int, DNA, DNA}}}
    leaf_defn_vars::Dict{String, Set{Tuple{Int, DNA, DNA}}}
    node_defn_vars::Dict{String, Set{Tuple{Int, DNA, DNA}}}

    function Haplogrouper()
        new(nothing, 
            Dict{String, Set{Tuple{Int, DNA, DNA}}}(), 
            Dict{String, Set{Tuple{Int, DNA, DNA}}}(),
            Dict{String, Set{Tuple{Int, DNA, DNA}}}())
    end
end
# regexp to match bootstrap:node_height = \)\d+:\d+(?:\.\d+)?
# Wrapper of Phylo.parsenewick()
function parse_newick(h::Haplogrouper, filename::String)
    h.tree = open(Phylo.parsenewick, Phylo.path(filename));
    for i = Phylo.nodeiter(h.tree)
        println(i);
    end
    return nothing
end 

function check_sample_names(h::Haplogrouper, seq_dict::Dict{String, BioSequence{DNAAlphabet{4}}})
    msa_samples = collect(keys(seq_dict))
    tree_leaves = collect(Phylo.getleafnames(h.tree))
    
    msa_only = filter(sample -> sample ∉ tree_leaves, msa_samples)
    tree_only = filter(leaf -> leaf ∉ msa_samples, tree_leaves)
    
    if isempty(msa_only) && isempty(tree_only)
        println("All samples match between MSA and tree ($(length(msa_samples)) samples)")
        return nothing
    else
        if !isempty(msa_only)
            error("Found $(length(msa_only)) samples in MSA that are not in tree: $msa_only")
        end
        if !isempty(tree_only)
            error("Found $(length(tree_only)) leaves in tree that are not in MSA: $tree_only")
        end
    end
end

function check_msa_length(seq_dict::Dict{String, BioSequence{DNAAlphabet{4}}})
    # Check for length consistency
    lengths = [length(seq) for seq in values(seq_dict)]
    if !all(x -> x == lengths[1], lengths)
        for (name, seq) in seq_dict
            println("$name: $(length(seq)) bases")
        end
        error("Not all sequences have the same length in the MSA")
    end
    return nothing
end

function parse_msa(h::Haplogrouper, filename::String, reference_id::String)
    records = collect(FASTA.Reader(open(filename, "r")))
    seq_dict = Dict{String, BioSequence{DNAAlphabet{4}}}()
    for record in records
        seq_dict[identifier(record)] = sequence(LongDNA{4}, record)
    end
    @argcheck haskey(seq_dict, reference_id) "Reference '$reference_id' not found in the MSA file"

    check_msa_length(seq_dict)
    check_sample_names(h, seq_dict) # confirm newick and msa have same samples

    ref_seq = seq_dict[reference_id]
    all_samples = collect(keys(seq_dict))
    println("Samples found: $all_samples")

    for sample in all_samples
        h.leaf_vars[sample] = Set{Tuple{Int, DNA, DNA}}()
    end

    total_variants = 0
    for sample in all_samples
        if sample == reference_id
            continue
        end
        variant_count = 0
        for base in eachindex(ref_seq)
            # Skip positions with N in either reference or sample
            if ref_seq[base] == DNA_N || seq_dict[sample][base] == DNA_N
                continue
            end
            if seq_dict[sample][base] != ref_seq[base]
                variant = (base, ref_seq[base], seq_dict[sample][base])
                push!(h.leaf_vars[sample], variant)
                variant_count += 1
                total_variants += 1
                #println("Variant at position $base in $sample: $(ref_seq[base]) → $(seq_dict[sample][base])")
            end
        end
        #println("Found $variant_count variants in $sample")
    end
    println("Total non-reference variants: $total_variants")
    return nothing
end

# If all samples have non-ref allele then the ref allele is defining for the ref!
function define_reference(h::Haplogrouper, reference_id::String)
    sets_of_positions = [Set(t[1] for t in tuples) for (key, tuples) in h.leaf_vars if key != reference_id]
    return reduce(intersect, sets_of_positions)
end

# This function loops through each sample and finds all unique variants
function define_leaves(h::Haplogrouper)
    for leaf in Phylo.getleafnames(h.tree)
        h.leaf_defn_vars[leaf] = Set{Tuple{Int, DNA, DNA}}()
    end

    # Create a dictionary to track all variants and associated leaves
    variant_to_leaves = Dict{Tuple{Int, DNA, DNA}, Set{String}}()
    for (leaf, variants) in h.leaf_vars
        for variant in variants
            if haskey(variant_to_leaves, variant)
                push!(variant_to_leaves[variant], leaf)
            else
                variant_to_leaves[variant] = Set([leaf])
            end
        end
    end

    for (leaf, variants) in h.leaf_vars
        #println("Processing leaf: $leaf")
        for variant in variants
            # If a variant is associated with only one leaf it is defining
            if length(variant_to_leaves[variant]) == 1
                #println("Defining variant for $leaf: $variant\n")
                push!(h.leaf_defn_vars[leaf], variant)
            end
        end
    end
    
    return nothing
end

# function define_nodes(h::Haplogrouper)
#     internal_nodes = filter(node -> node ∉ Phylo.getleafnames(h.tree), Phylo.getnodenames(h.tree))
#     node_to_descendants = Dict(node => collect_descendant_leaves(h.tree, node) for node in internal_nodes)
    
#     # Create a mapping from node to its parent
#     parent_map = Dict{String, String}()
#     for node in internal_nodes
#         for child in Phylo.getchildren(h.tree, node)
#             parent_map[child] = node
#         end
#     end
    
#     # Initialize node defining variants
#     for node in internal_nodes
#         h.node_defn_vars[node] = Set{Tuple{Int, DNA, DNA}}()
#     end
    
#     # First collect all potential variants for each node
#     node_potential_vars = Dict{String, Set{Tuple{Int, DNA, DNA}}}()
#     for node in internal_nodes
#         println("Processing node: $node")
#         descendants = node_to_descendants[node]
        
#         # Collect all variants from descendants
#         all_descendant_variants = Set{Tuple{Int, DNA, DNA}}()
#         for descendant in descendants
#             union!(all_descendant_variants, h.leaf_vars[descendant])
#         end
        
#         # Find variants present in ALL descendants
#         potential_defining_vars = Set{Tuple{Int, DNA, DNA}}()
#         for variant in all_descendant_variants
#             if all(variant ∈ h.leaf_vars[descendant] for descendant in descendants)
#                 push!(potential_defining_vars, variant)
#             end
#         end
        
#         node_potential_vars[node] = potential_defining_vars
#     end
    
#     # Process nodes from most ancestral to most derived
#     # (We'll do this by considering nodes with fewer descendants first, 
#     # which generally means they're deeper in the tree)
#     sorted_nodes = sort(collect(internal_nodes), by=node -> -length(node_to_descendants[node]))

#     for node in sorted_nodes
#         for variant in node_potential_vars[node]
#             # Check if any ancestor node already has this variant as a defining variant
#             is_already_defined_in_ancestor = false
#             current = node
#             while haskey(parent_map, current)
#                 #parent = parent_map[current]
#                 parent = Phylo.getparent(h.tree, current)
#                 if variant ∈ h.node_defn_vars[parent]
#                     is_already_defined_in_ancestor = true
#                     break
#                 end
#                 current = parent
#             end
            
#             # If not defined in ancestor, add it to this node
#             if !is_already_defined_in_ancestor
#                 push!(h.node_defn_vars[node], variant)
#             end
#         end
#         println("Node $node descendants: $(node_to_descendants[node])")
#         println("$node defining variants: $(h.node_defn_vars[node])\n")
#     end
    
#     return nothing
# end


# Order the internal nodes from root to leaves 
# For each internal node get all of its descendants and the set of all shared variants
# For each variant we consider it defining for this node if all descendants have it AND if no ancestral nodes have already been assigned the variant


function define_nodes(h::Haplogrouper)
    internal_nodes = reverse(filter(node -> node ∉ Phylo.getleafnames(h.tree), Phylo.getnodenames(h.tree)))
    #node_to_descendants = Dict(node => collect_descendant_leaves(h.tree, node) for node in internal_nodes)
    for node in internal_nodes
        h.node_defn_vars[node] = Set{Tuple{Int, DNA, DNA}}()
        descendants = collect_descendant_leaves(h.tree, node)
        
        all_descendant_variants = Set{Tuple{Int, DNA, DNA}}()
        for descendant in descendants
            union!(all_descendant_variants, h.leaf_vars[descendant])
        end

        #all_descendant_variants = reduce(intersect, [h.leaf_vars[k] for k in descendants])
        #println(vals)

        # Find defining variants for this node
        for variant in all_descendant_variants
            # A variant is node-defining if:
            # 1. It appears in ALL descendants of the node
            # 2. It has NOT been placed on an ancestral node
            if all(variant ∈ h.leaf_vars[descendant] for descendant in descendants)
                ancestors = Phylo.getancestors(h.tree, node)
                in_ancestor = false
                for ancestor in ancestors
                    if variant ∈ h.node_defn_vars[ancestor]
                        in_ancestor = true
                        break
                    end
                end
                # If not defined in any ancestor, add it to this node
                if !in_ancestor
                    push!(h.node_defn_vars[node], variant)
                end
            end
        end
        #println("Node $node descendants: $(node_to_descendants[node])")
        println("$node defining variants: $(h.node_defn_vars[node])\n")
    end
    
    return nothing
end

# Helper function to recursively collect all descendant leaves of a node
function collect_descendant_leaves(tree, node)
    children = Phylo.getchildren(tree, node)
    leaves = Set{String}()
    for child in children
        if child in Phylo.getleafnames(tree)
            push!(leaves, child)
        else
            union!(leaves, collect_descendant_leaves(tree, child))
        end
    end
    return leaves
end

msa_path = "/Users/tomroberts/Documents/Oxford_DPhil/second_rotation/chicken_practice/sequence_files/miao_mtdna_seqs_1-61_msa2.fasta"
newick_path = "/Users/tomroberts/Documents/Oxford_DPhil/second_rotation/chicken_practice/full_mitogenome_iqtree/full_mtdna_1-61.treefile"


analyser = Haplogrouper();
parse_newick(analyser, newick_path);
parse_msa(analyser, msa_path, "NC_007235")
define_leaves(analyser);
define_nodes(analyser;)

# iters = collect(Phylo.nodeiter(analyser.tree));
# println("\n", iters[1].name)

# println(Phylo.getnodenames(analyser.tree))




a = [i*i for i in 1:10]
println(a)
println(intersect([9,4,3], [1,4,5]))


my_dict =  Dict("A"=>(1,1), "B"=>(2,2), "C"=>(3,3), "D"=>(1,1))
vals = [my_dict[k] for k in ["A", "B", "C", "D"]]
println(vals)
println(unique(vals))
end