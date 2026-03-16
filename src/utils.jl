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