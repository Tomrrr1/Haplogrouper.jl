# Reconstruct the sequences at each node using Fitch's algorithm
function reconstruct_ancestral_seq(tree::RecursiveTree, nucleotide_maps::Dict{String, Set{Tuple{Int, DNA, DNA}}}, reference_id::String)
    ancestral_states = Dict{String, Dict{Int, DNA}}() # node_name, pos, base
    all_nodes = Phylo.getnodenames(tree)

    for node in all_nodes
        ancestral_states[node] = Dict{Int, DNA}()
    end
    all_positions = Set(t[1] for t in nucleotide_maps[reference_id])
    for position in all_positions
        # Get leaf states for all samples at a given position
        leaf_states = get_base_at_pos(tree, nucleotide_maps, position)
        # Use fitch algorithm to reconstruct the ancestral state at this position
        position_states = fitch(tree, leaf_states)
        for (node, state) in position_states
            ancestral_states[node][position] = state
        end
    end

    return ancestral_states
end 

# Fitch algorithm for ancestral state reconstruction
function fitch(tree::RecursiveTree, leaf_states::Dict{String, DNA})
    node_state_sets = Dict{String, Set{DNA}}()
    leaves = Phylo.getleafnames(tree)
    all_nodes = Phylo.getnodenames(tree)

    # Initialise the leaf states
    for node in all_nodes
        if node in leaves
            node_state_sets[node] = Set([leaf_states[node]])
        end
    end
    node_state_sets = postorder_pass(tree, node_state_sets)
    final_states = topdown_pass(tree, node_state_sets)

    return final_states
end

# Perform the postorder pass of the fitch algorithm
function postorder_pass(tree::RecursiveTree, node_state_sets::Dict{String, Set{DNA}})
    node_postorder = Phylo.traversal(tree, Phylo.postorder)
    internal_node_postorder = filter(node -> node ∉ Phylo.getleaves(tree), node_postorder)

    for node_idx in internal_node_postorder
        node_name = Phylo.getnodename(tree, node_idx)
        children = Phylo.getchildren(tree, node_idx)
        children_names = [Phylo.getnodename(tree, c) for c in children]
        # Assume an internal node a with children b and c.
        if !isempty(children_names)
            # We have the state sets S_b and S_c, which is simply the nucleotide at some position X.
            child_state_sets = [node_state_sets[c] for c in children_names]
            # If the intersection of S_b and S_c is non-empty this becomes S_a (S_b and S_c have the same nucleotide)
            intersect_set = reduce(intersect, child_state_sets) 
            if !isempty(intersect_set)
                    node_state_sets[node_name] = intersect_set
            else # If the intersect is empty, S_a becomes the union of S_b and S_c.
                union_set = reduce(union, child_state_sets)
                node_state_sets[node_name] = union_set
            end
        end
    end

    return node_state_sets
end

# Perform the topdown pass of the fitch algorithms
function topdown_pass(tree::RecursiveTree, node_state_sets::Dict{String, Set{DNA}})
    final_states = Dict{String, DNA}()
    # If the state set of the root contains more than one element, arbitrarily assign any of these characters to the root (we use first()).
    root_name = Phylo.getnodename(tree, Phylo.getroot(tree))
    final_states[root_name] = first(node_state_sets[root_name])
    stack = [root_name]
    while !isempty(stack)
        current_node = pop!(stack)
        current_node_state = final_states[current_node]
        # Let b be a child of node a, and let x_a denote the character assigned to a
        for child_node in Phylo.getchildren(tree, Phylo.getnode(tree, current_node))
            child_name = Phylo.getnodename(tree, child_node)
            child_state_set = node_state_sets[child_name]
            # If x_a is contained in S_b, assign it to b (the parent and child node have the same nucleotide)
            if current_node_state in child_state_set
                final_states[child_name] = current_node_state
            else
                # Otherwise, arbitrarily assign any character from S_b to b
                if !isempty(child_state_set)
                    final_states[child_name] = first(child_state_set)
                end
            end
            # 
            push!(stack, child_name)
        end
    end
    return final_states
end