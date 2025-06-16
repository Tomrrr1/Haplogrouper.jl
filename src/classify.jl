function classify_sample(tree::RecursiveTree, all_node_vars::VariantMap, sample_variants::VariantSet, scoring_metric::String="kulczynski")
    internal_nodes = reverse(filter(node -> node ∉ Phylo.getleafnames(tree), Phylo.getnodenames(tree)))
    qualifying_nodes = Tuple{String, VariantSet, Float64, Int}[]
    for n in internal_nodes
        n_name = Phylo.getnodename(tree, n)
        n_vars = if n_name in keys(all_node_vars)
            all_node_vars[n_name]
        else
            VariantSet() # Default empty set
        end

        isempty(n_vars) && continue # Skip empty nodes
        score = if scoring_metric == "jaccard"
            jaccard(sample_variants, n_vars)
        elseif scoring_metric == "kulczynski"
            kulczynski(sample_variants, n_vars)
        end
        push!(qualifying_nodes, (n_name, n_vars, score, length(n_vars)))
    end

    sort!(qualifying_nodes, by = x -> x[3], rev=true)
    highest_score = qualifying_nodes[1][3]
    best_nodes = filter(x -> x[3] == highest_score, qualifying_nodes)
    if length(best_nodes) == 1
        return [[best_nodes[1][1], best_nodes[1][2], best_nodes[1][3]]]
    else
        return [[node[1], node[2], node[3]] for node in best_nodes]
    end
end

function classify_samples_from_msa(tree::RecursiveTree, all_node_vars::VariantMap, msa_file::String, ref_id::String, ancestral_state::Dict{Int, DNA}, scoring_metric::String="kulczynski")
    sample_vars = parse_msa_classify(msa_file, ref_id, ancestral_state)
    assignments = Dict{String, Vector{}}()

    for (sample_id, variants) in sample_vars
        nodes = classify_sample(tree, all_node_vars, variants, scoring_metric)
        assignments[sample_id] = nodes
    end
    
    return assignments
end

"""
    classify(msa_file::String, newick_file::String, scaffold_file::String, ancestral_seq_file::String, ref_id::String, scoring_metric::String = "kulczynski")

Assign new samples to nodes on a previously defined scaffold tree.

# Arguments
- `msa_file::String`: Path to multiple sequence alignment file in FASTA format.
- `newick_file::String`: Path to the tree file in Newick format.
- `scaffold_file::String`: Path to the scaffold tree in TSV format to be used for sample classification. This file must be generated using Haplogrouper.make_scaffold().
- `ancestral_seq_file:: String`: Path to the ancestral sequence file in TSV format. This file must be generated using Haplogrouper.make_scaffold()
- `ref_id::String`: Sequence ID of the reference sample. This must be the same as the reference used when creating the scaffold tree using Haplogrouper.make_scaffold().
- `scoring_metric::String`: A string denoting the scoring method to be used for sample classification. The options are "kulczynski" and "jaccard". The default is "kulczynski".

# Returns
- `nothing`

# Example 
```jldoctest
julia> msa_string = ">ref\\nATCGATCG\\n>s1_new\\nATGGATCT\\n>s2_new\\nATGGATCG\\n>s3_new\\nATCGATCC\\n>s4_new\\nATCGATGC";

julia> tree_string = "((s1,s2),(s3,s4),ref);";

julia> scaffold_string = "Node\\tDescendants\\tDefining_Variants\\nNode 8\\tref,s1,s2,s3,s4\\t\\nNode 3\\ts1,s2\\t3:C>G\\nNode 6\\ts3,s4\\t8:G>C";

julia> ancestral_string = "1\\tA\\n2\\tT\\n3\\tC\\n4\\tG\\n5\\tA\\n6\\tT\\n7\\tC\\n8\\tG\\n";

julia> write("test_sequences.fasta", msa_string);

julia> write("test_tree.nwk", tree_string);

julia> write("test_scaffold.tsv", scaffold_string);

julia> write("test_ancestral.tsv", ancestral_string);

julia> classify("test_sequences.fasta", "test_tree.nwk", "test_scaffold.tsv", "test_ancestral.tsv", "ref", "kulczynski")
Using kulczynski for sample classification.

Assignment results:
s3_new -> Vector{Any}[["Node 6", Set(Tuple{Int64, BioSymbols.DNA, BioSymbols.DNA}[(8, DNA_G, DNA_C)]), 1.0]]
s4_new -> Vector{Any}[["Node 6", Set(Tuple{Int64, BioSymbols.DNA, BioSymbols.DNA}[(8, DNA_G, DNA_C)]), 0.75]]
s1_new -> Vector{Any}[["Node 3", Set(Tuple{Int64, BioSymbols.DNA, BioSymbols.DNA}[(3, DNA_C, DNA_G)]), 0.75]]
s2_new -> Vector{Any}[["Node 3", Set(Tuple{Int64, BioSymbols.DNA, BioSymbols.DNA}[(3, DNA_C, DNA_G)]), 1.0]]
```
See also [`Haplogrouper.make_scaffold`](@ref)
"""
function classify(msa_file::String, newick_file::String, scaffold_file::String, ancestral_seq_file::String, ref_id::String, scoring_metric::String = "kulczynski")
    all_node_vars = read_variants_tsv(scaffold_file)
    tree = parse_newick(newick_file)
    ancestral_seq = read_ancestral_sequence_tsv(ancestral_seq_file)

    println("Using $scoring_metric for sample classification.")
    assignments = classify_samples_from_msa(tree, all_node_vars, msa_file, ref_id, ancestral_seq, scoring_metric)

    println("\nAssignment results:")
    for (sample, node) in assignments
        println("$sample -> $node")
    end

    return nothing
end