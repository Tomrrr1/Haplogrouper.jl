# Ancestral sequence reconstruction-based variant calling and node definition

function define_nodes_ancestral(tree::RecursiveTree, nucleotide_maps::VariantMap, reference_id::String, prefix::String="")
    ancestral_states = reconstruct_ancestral_seq(tree, nucleotide_maps, reference_id)
    ancestral_filename = isempty(prefix) ? "./ancestral_sequence.tsv" : "./$(prefix)_ancestral_sequence.tsv"
    write_ancestral_sequence_tsv(ancestral_states, Phylo.getnodename(tree, Phylo.getroot(tree)), ancestral_filename)
    node_vars = build_node_vars_ancestral(tree, ancestral_states)
    
    return node_vars
end

function build_node_vars_ancestral(tree::RecursiveTree, ancestral_states::Dict{String, Dict{Int, DNA}})
    node_vars = VariantMap()
    all_nodes = Phylo.getnodenames(tree)
    root_name = Phylo.getnodename(tree, Phylo.getroot(tree))
    root_states = ancestral_states[root_name]

    for node in all_nodes
        node_vars[node] = VariantSet()
        for (pos, node_base) in ancestral_states[node]
            if !haskey(root_states, pos)
                continue
            end
            root_base = root_states[pos]
            # Don't assign - or N as defining variants
            if node_base != root_base && 
               node_base != DNA_Gap && node_base != DNA_N && 
               root_base != DNA_Gap && root_base != DNA_N
                push!(node_vars[node], (pos, root_base, node_base))
            end
        end
    end

    return node_vars
end