"""
    output_results_tsv(tree::RecursiveTree, leaf_vars::VariantMap, node_vars::VariantMap, output_file::String)

Write the results of Haplogrouper::define_phylogeny() to a TSV file.

# Arguments
- `tree::RecursiveTree`
- `leaf_vars::VariantMap`: A dictionary of leaves and their corresponding variants relative to a reference.
- `node_vars::VariantMap`: A dictionary of internal nodes and their corresponding variants relative to a reference.
- `output_file::String`: Path to save the output TSV file.

# Returns
- `Nothing`
"""
function output_results_tsv(tree::RecursiveTree, leaf_vars::VariantMap, node_vars::VariantMap, output_file::String)
    open(output_file, "w") do file
        write(file, "Node\tDescendants\tDefining_Variants\n")
        
        # Process internal nodes
        for node in keys(node_vars)
            descendants = collect_descendant_leaves(tree, node)
            descendants_str = join(descendants, ",")
            
            variants = get(node_vars, node, VariantSet())
            variant_strs = ["$(v[1]):$(v[2])>$(v[3])" for v in variants]
            variants_str = join(variant_strs, ",")
            
            write(file, "$node\t$descendants_str\t$variants_str\n")
        end
        
        # Process leaves
        for leaf in keys(leaf_vars)
            variants = get(leaf_vars, leaf, VariantSet())
            variant_strs = ["$(v[1]):$(v[2])>$(v[3])" for v in variants]
            variants_str = join(variant_strs, ",")
            
            write(file, "$leaf\t$leaf\t$variants_str\n")
        end
    end
    
    println("Results written to TSV file: $output_file")
    return nothing
end