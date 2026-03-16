# Outgroup reference-based variant calling and node definition

function define_nodes_outgroup(tree::RecursiveTree, nucleotide_maps::VariantMap, ref_id::String)
    node_vars = build_node_vars_outgroup(tree, nucleotide_maps, ref_id)
    
    return node_vars
end

function build_node_vars_outgroup(tree::RecursiveTree, nucleotide_maps::VariantMap, ref_id::String)
    node_vars = VariantMap()
    leaf_vars = VariantMap()
    all_nodes = Phylo.getnodenames(tree)

    for (sample_id, sample_data) in nucleotide_maps
        # Skip the outgroup as it has no variants relative to itself
        if sample_id == ref_id
            continue
        end
        leaf_vars[sample_id] = VariantSet()
        
        for (pos, ref_base, sample_base) in sample_data
            # Only include actual variants (differences from outgroup reference)
            if sample_base != ref_base && 
               sample_base != DNA_Gap && sample_base != DNA_N && 
               ref_base != DNA_Gap && ref_base != DNA_N
                push!(leaf_vars[sample_id], (pos, ref_base, sample_base))
            end
        end
    end
    
    # Now build node variants
    for node in all_nodes
        node_vars[node] = VariantSet()
        
        if node in Phylo.getleafnames(tree)
            # Leaf nodes
            if node == ref_id
                continue
            end
            node_vars[node] = copy(leaf_vars[node])
        else
            # Internal nodes (find variants shared by all descendants of the node)
            descendants = collect_descendant_leaves(tree, node)
            valid_descendants = filter(d -> haskey(leaf_vars, d), descendants)
            node_vars[node] = reduce(intersect, [leaf_vars[d] for d in valid_descendants])
        end
    end
    
    return node_vars
end
