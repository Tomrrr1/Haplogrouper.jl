# Wrapper of Phylo.parsenewick
function parse_newick(filename::String)
    tree = open(Phylo.parsenewick, filename)
    for node in Phylo.nodeiter(tree)
        #println(node)
    end

    return tree
end 

# Recursively collect all descendant leaves of a node
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

# Get a sample's nucleotide at a given position
function get_base_at_pos(tree::RecursiveTree, nucleotide_maps::VariantMap, position::Int64)
    leaf_states = Dict{String, DNA}()
    leaves = Phylo.getleafnames(tree)
    for leaf in leaves
        current_leaf = nucleotide_maps[leaf]
        found = false
        for tup in current_leaf
            if tup[1] == position
                leaf_states[leaf] = tup[3] # [pos, ref, alt]
                found = true
                break
            end
        end
        if !found
            error("No base found for leaf $leaf at position $position")
        end
    end

    return leaf_states
end 

# Confirm that the sample names in the msa match those in the tree
function check_sample_names(tree::RecursiveTree, seq_dict::Dict{String, BioSequence{DNAAlphabet{4}}})
    msa_samples = collect(keys(seq_dict))
    tree_samples = collect(Phylo.getleafnames(tree))
    
    msa_only = filter(sample -> sample ∉ tree_samples, msa_samples)
    tree_only = filter(leaf -> leaf ∉ msa_samples, tree_samples)
    
    if isempty(msa_only) && isempty(tree_only)
        println("\nAll samples match between MSA and tree ($(length(msa_samples)) samples)")
        return true
    else
        if !isempty(msa_only)
            error("Found $(length(msa_only)) samples in MSA that are not in tree: $msa_only")
        end
        if !isempty(tree_only)
            error("Found $(length(tree_only)) leaves in tree that are not in MSA: $tree_only")
        end
    end
end

# Confirm that all sequences in the msa have the same length
function check_msa_length(seq_dict::Dict{String, BioSequence{DNAAlphabet{4}}})
    # Check for length consistency
    lengths = [length(seq) for seq in values(seq_dict)]
    if !all(x -> x == lengths[1], lengths)
        for (name, seq) in seq_dict
            println("$name: $(length(seq)) bases")
        end
        error("Not all sequences have the same length in the MSA")
    end
    
    return true
end