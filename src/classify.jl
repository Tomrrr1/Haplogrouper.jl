# Helper function to print classification results
function print_classification_results(assignments::Dict{String, Vector{Vector{Any}}}, scoring_metric::String)
    println("Using $scoring_metric for sample classification.")
    println("\nAssignment results:")
    println("=" ^ 80)
    
    for (sample_id, results) in assignments
        println("Sample: $sample_id")
        
        # Show best assignment
        best_result = results[1]
        node_name, shared_variants, unique_to_sample, unique_to_node, score = best_result
        
        println("  Best assignment:")
        println("    Node: $node_name")
        println("    Score: $(round(score, digits=3))")
        
        if node_name != "No classification" && node_name != "No classification due to score < 0.1"
            shared_strs = ["$(pos):$(ref)>$(alt)" for (pos, ref, alt) in shared_variants]
            sample_strs = ["$(pos):$(ref)>$(alt)" for (pos, ref, alt) in unique_to_sample]
            node_strs = ["$(pos):$(ref)>$(alt)" for (pos, ref, alt) in unique_to_node]
            
            println("    Shared variants ($(length(shared_variants))): $(join(shared_strs, ", "))")
            println("    Sample-specific variants ($(length(unique_to_sample))): $(join(sample_strs, ", "))")
            println("    Node-specific variants ($(length(unique_to_node))): $(join(node_strs, ", "))")
        end
        
        # Show second best if available
        if length(results) > 1
            second_result = results[2]
            node_name2, shared_variants2, unique_to_sample2, unique_to_node2, score2 = second_result
            
            println("  Second best assignment:")
            println("    Node: $node_name2")
            println("    Score: $(round(score2, digits=3))")
            
            if node_name2 != "No classification" && node_name2 != "No classification due to score < 0.1"
                shared_strs2 = ["$(pos):$(ref)>$(alt)" for (pos, ref, alt) in shared_variants2]
                sample_strs2 = ["$(pos):$(ref)>$(alt)" for (pos, ref, alt) in unique_to_sample2]
                node_strs2 = ["$(pos):$(ref)>$(alt)" for (pos, ref, alt) in unique_to_node2]
                
                println("    Shared variants ($(length(shared_variants2))): $(join(shared_strs2, ", "))")
                println("    Sample-specific variants ($(length(unique_to_sample2))): $(join(sample_strs2, ", "))")
                println("    Node-specific variants ($(length(unique_to_node2))): $(join(node_strs2, ", "))")
            end
        else
            println("  Second best assignment: None")
        end
        
        println("-" ^ 80)
    end
end

function classify_sample(tree::RecursiveTree, all_node_vars::VariantMap, sample_vars::VariantSet, scoring_metric::String)
    internal_nodes = reverse(filter(node -> node ∉ Phylo.getleafnames(tree), Phylo.getnodenames(tree)))
    qualifying_nodes = Vector{Tuple{String, VariantSet, VariantSet, VariantSet, Float64}}()
    
    # Iterate through each node and calculate the sample assignment score
    for n in internal_nodes
        n_name = Phylo.getnodename(tree, n)
        n_vars = get(all_node_vars, n_name, VariantSet())
        isempty(n_vars) && continue
        
        shared_variants = intersect(sample_vars, n_vars)
        unique_to_sample = setdiff(sample_vars, n_vars)
        unique_to_node = setdiff(n_vars, sample_vars)

        if scoring_metric == "jaccard"
            score = jaccard(sample_vars, n_vars)
        elseif scoring_metric == "kulczynski"
            score = kulczynski(sample_vars, n_vars)
        else 
            error("Unknown scoring metric $scoring_metric. Use 'kulczynski' or 'jaccard'")
        end
        push!(qualifying_nodes, (n_name, shared_variants, unique_to_sample, unique_to_node, score))
    end

    # Sort qualifying nodes by score and return the top 2 nodes for each sample
    sort!(qualifying_nodes, by = x -> x[5], rev=true)
    results = []
    for i in 1:min(5, length(qualifying_nodes))
        node_name, shared, unique_to_sample, unique_to_node, score = qualifying_nodes[i]
        push!(results, [node_name, shared, unique_to_sample, unique_to_node, score])
    end
    
    return results
end

# Method 1: Helper function for ancestral-based classification
function classify_samples_from_msa(
    tree::RecursiveTree, all_node_vars::VariantMap, msa_file::String, ref_id::String, 
    ancestral_state::Dict{Int, DNA}, scoring_metric::String="kulczynski"
    )
    sample_vars = parse_msa_classify(msa_file, ref_id, ancestral_state) # Use ancestral states
    assignments = Dict{String, Vector{Vector{Any}}}()

    for (sample_id, variants) in sample_vars
        nodes = classify_sample(tree, all_node_vars, variants, scoring_metric)
        assignments[sample_id] = nodes
    end
    
    return assignments
end

# Method 2: Helper function for outgroup-based classification  
function classify_samples_from_msa(
    tree::RecursiveTree, all_node_vars::VariantMap, msa_file::String, 
    ref_id::String, scoring_metric::String
    )
    sample_vars = parse_msa_classify(msa_file, ref_id)
    assignments = Dict{String, Vector{Vector{Any}}}()

    for (sample_id, variants) in sample_vars
        nodes = classify_sample(tree, all_node_vars, variants, scoring_metric)
        assignments[sample_id] = nodes
    end
    
    return assignments
end

"""
    classify(msa_file::String, newick_file::String, scaffold_file::String, ancestral_seq_file::String, ref_id::String, scoring_metric::String)

Assign new samples to nodes on a previously annotated phylogenetic tree (using ancestral sequence reconstruction for annotation).

# Arguments
- `msa_file::String`: Path to multiple sequence alignment file in FASTA format.
- `newick_file::String`: Path to the tree file in Newick format.
- `scaffold_file::String`: Path to the scaffold tree in TSV format to be used for sample classification. This file must be generated using Haplogrouper.annotate().
- `ref_id::String`: Sequence ID of the reference sample. This must be the same as the reference used when creating the scaffold tree using Haplogrouper.annotate().
- `ancestral_seq_file:: String`: Path to the ancestral sequence file in TSV format. This file must be generated using Haplogrouper.annotate()
- `scoring_metric::String`: A string denoting the scoring method. Options: "kulczynski" and "jaccard".

# Returns
- `nothing`
"""
function classify(msa_file::String, newick_file::String, scaffold_file::String, ref_id::String, ancestral_seq_file::String, scoring_metric::String)
    all_node_vars = read_variants_tsv(scaffold_file)
    tree = parse_newick(newick_file)
    ancestral_seq = read_ancestral_sequence_tsv(ancestral_seq_file)

    assignments = classify_samples_from_msa(tree, all_node_vars, msa_file, ref_id, ancestral_seq, scoring_metric)
    print_classification_results(assignments, scoring_metric)
    output_assignments_csv(assignments, "haplogrouper_assignments.csv")
    
    return nothing
end

"""
    classify(msa_file::String, newick_file::String, scaffold_file::String, ref_id::String, scoring_metric::String)

Assign new samples to nodes on a previously annotated phylogenetic tree (using outgroup mode for annotation).

# Arguments
- `msa_file::String`: Path to multiple sequence alignment file in FASTA format.
- `newick_file::String`: Path to the tree file in Newick format.
- `scaffold_file::String`: Path to the scaffold tree in TSV format to be used for sample classification. This file must be generated using Haplogrouper.annotate().
- `ref_id::String`: Sequence ID of the reference sample. This must be the same as the reference used when creating the scaffold tree using Haplogrouper.annotate().
- `scoring_metric::String`: A string denoting the scoring method. Options: "kulczynski" and "jaccard".

# Returns
- `nothing`
"""
function classify(msa_file::String, newick_file::String, scaffold_file::String, ref_id::String, scoring_metric::String)
    all_node_vars = read_variants_tsv(scaffold_file)
    tree = parse_newick(newick_file)
    
    assignments = classify_samples_from_msa(tree, all_node_vars, msa_file, ref_id, scoring_metric)
    print_classification_results(assignments, scoring_metric)
    output_assignments_csv(assignments, "haplogrouper_assignments.csv")
    
    return nothing
end