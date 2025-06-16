# Find node defining variants relative to the reconstructed ancestral sequence
function define_nodes(tree::RecursiveTree, nucleotide_maps::VariantMap, reference_id::String)
    ancestral_states = reconstruct_ancestral_seq(tree, nucleotide_maps, reference_id)
    write_ancestral_sequence_tsv(ancestral_states, Phylo.getnodename(tree, Phylo.getroot(tree)), "./ancestral_sequence.tsv")
    node_vars = build_node_vars(tree, ancestral_states)
    
    return node_vars
end

# Get all variants shared by the descendants of each node
function build_node_vars(tree::RecursiveTree, ancestral_states::Dict{String, Dict{Int, DNA}})
    node_vars = VariantMap()
    all_nodes = Phylo.getnodenames(tree)
    root_name = Phylo.getnodename(tree, Phylo.getroot(tree))
    root_states = ancestral_states[root_name]

    # If a 
    for node in all_nodes
        node_vars[node] = VariantSet()
        for (pos, node_base) in ancestral_states[node]
            if !haskey(root_states, pos)
                continue
            end
            root_base = root_states[pos]
            if node_base != root_base
                push!(node_vars[node], (pos, root_base, node_base))
            end
        end
    end

    return node_vars
end