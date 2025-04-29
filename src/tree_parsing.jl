# Wrapper of Phylo.parsenewick
function parse_newick(filename::String)
    tree = open(Phylo.parsenewick, Phylo.path(filename))
    for i in Phylo.nodeiter(tree)
        println(i)
    end
    return tree
end 