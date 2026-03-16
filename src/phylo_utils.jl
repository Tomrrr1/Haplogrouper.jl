# Wrapper of Phylo.parsenewick
function parse_newick(filename::String)
    tree = open(Phylo.parsenewick, filename)

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