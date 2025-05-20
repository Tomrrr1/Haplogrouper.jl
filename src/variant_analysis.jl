# If all samples have non-ref allele then the ref allele is defining for the ref!
function define_reference(leaf_vars::VariantMap, ref_id::String)
    sets_of_positions = [Set(t[1] for t in ts) for (k, ts) in leaf_vars if k != ref_id]
    return reduce(intersect, sets_of_positions)
end

# This function finds unique variants for each sample
function define_leaves(leaf_vars::VariantMap)
    # Create a dictionary to track all variants and associated leaves
    variant_to_leaves = Dict{Variant, Set{String}}()
    for (leaf, variants) in leaf_vars
        for variant in variants
            if haskey(variant_to_leaves, variant)
                push!(variant_to_leaves[variant], leaf)
            else
                variant_to_leaves[variant] = Set([leaf])
            end
        end
    end

    leaf_defn_vars = VariantMap()
    for (leaf, variants) in leaf_vars
        leaf_defn_vars[leaf] = VariantSet()
        for variant in variants
            if length(variant_to_leaves[variant]) == 1
                push!(leaf_defn_vars[leaf], variant)
                println("Defining variant for $leaf: $variant\n")
            end
        end
    end
    
    return leaf_defn_vars
end

# Identify defining variants for internal nodes
function define_nodes(tree::RecursiveTree, leaf_vars::Dict{String, Set{Variant}})
    internal_nodes = reverse(filter(node -> node ∉ Phylo.getleafnames(tree), Phylo.getnodenames(tree)))
    node_defn_vars = VariantMap()
    all_node_vars = VariantMap()
    
    for node in internal_nodes
        node_defn_vars[node] = VariantSet()
        descendants = collect_descendant_leaves(tree, node)
        shared_variants = reduce(intersect, [leaf_vars[d] for d in descendants]) 
        all_node_vars[node] = shared_variants # all variants shared by node descendants

        # Find unique defining variants for this node
        for variant in shared_variants
            ancestors = Phylo.getancestors(tree, node)
            in_ancestor = false
            for ancestor in ancestors
                if variant in node_defn_vars[ancestor]
                    in_ancestor = true
                    break
                end
            end

            if !in_ancestor
                push!(node_defn_vars[node], variant)
            end
        end
    end
    
    return [all_node_vars, node_defn_vars]
end