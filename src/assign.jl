"""
    assign_sample(tree::RecursiveTree, all_node_vars::Dict{String, Set{Variant}}, sample_variants::Set{Variant})

Assign a new sample to the most appropriate node in the phylogenetic tree based on shared variants.

Parameters:
- `tree`: The phylogenetic tree 
- `all_node_vars`: Dictionary mapping node names to their defining variants
- `sample_variants`: Set of variants found in the sample to be classified

Returns the name of the assigned node.
"""

function assign_sample(tree::RecursiveTree, all_node_vars::Dict{String, Set{Variant}}, sample_variants::Set{Variant})
    internal_nodes = reverse(filter(node -> node ∉ Phylo.getleafnames(tree), Phylo.getnodenames(tree)))
    confidence = 0
    best_node_name = ""
    best_node_vars = Set{Variant}()

    # Loop through all nodes and assign our sample to the one with the highest score
    for n in internal_nodes
        n_name = Phylo.getnodename(tree, n)
        n_vars = if n_name in keys(all_node_vars)
            all_node_vars[n_name]
        else
            Set{Variant}()  # Default empty set
        end

        isempty(n_vars) && continue # Skip empty nodes
        matching_vars = intersect(sample_variants, n_vars) # shared variants [1,2] [1,2,3] => [1,2,3]
        union_vars = union(sample_variants, n_vars)
        score = length(matching_vars) / length(union_vars) # jaccard index
        if score > confidence
            best_node_name = n_name
            best_node_vars = n_vars
            confidence = score
        end
    end
    
    return [best_node_name, best_node_vars, confidence]
end

"""
    assign_samples_from_msa(tree::RecursiveTree, all_node_vars::Dict{String, Set{Variant}}, 
                            msa_file::String, ref_id::String)

Process all samples in an MSA file and assign each to the appropriate node in the phylogenetic tree.

Parameters:
- `tree`: The phylogenetic tree
- `all_node_vars`: Dictionary mapping internal node names to their defining variants
- `msa_file`: Path to the MSA file (FASTA format)
- `ref_id`: ID of the reference sequence in the MSA

Returns a dictionary mapping sample IDs to their assigned nodes.
"""
function assign_samples_from_msa(tree::RecursiveTree, all_node_vars::Dict{String, Set{Variant}}, msa_file::String, ref_id::String)
    # Process the samples to identify their variants. Skip tree check with 'true'
    sample_vars = parse_msa(msa_file, ref_id, tree, true)
    
    # Skip the reference sequence
    delete!(sample_vars, ref_id)
    
    # Assign each sample to a node
    assignments = Dict{String, Vector{}}()
    for (sample_id, variants) in sample_vars
        node = assign_sample(tree, all_node_vars, variants)
        assignments[sample_id] = node
        println("Sample '$sample_id' assigned to node: $node")
    end
    
    return assignments
end