module Haplogrouper

# Import external packages
using ArgCheck
using BioSequences
using DataFrames
using CSV
using FASTX
using Phylo

# Include all module files
include("types.jl")
include("tree_parsing.jl") 
include("msa_parsing.jl")
include("variant_analysis.jl")
include("variant_parsing.jl")
include("output.jl")
include("utils.jl")
include("assign.jl")

export AnalysisResult
export parse_newick
export parse_msa
export define_leaves
export define_nodes
export define_variants
export define_phylogeny
export output_results_tsv
export assign_samples
export assign_samples_from_msa
export read_variants_from_tsv

"""
    define_phylogeny(msa_file::String, newick_file::String, ref_id::String)

Identify defining variants for each branch on a phylogeny

# Arguments
- `msa_file::String`: Path to multiple sequence alignment file in FASTA format
- `newick_file::String`: Path to the Newick tree file
- `ref_id::String`: ID of the reference sequence in the MSA

# Returns
- `AnalysisResult`: Structure containing all analysis results
"""
function define_phylogeny(msa_file::String, newick_file::String, ref_id::String)
    tree = parse_newick(newick_file)
    leaf_vars = parse_msa(msa_file, ref_id, tree, false)
    
    # Find defining variants for leaves and nodes
    leaf_defn_vars = define_leaves(leaf_vars)
    node_vars_result = define_nodes(tree, leaf_vars)
    node_vars = node_vars_result[1]
    node_defn_vars = node_vars_result[2]

    output_results_tsv(tree, leaf_vars, node_vars, "haplogrouper_all_variants.tsv")
    output_results_tsv(tree, leaf_defn_vars, node_defn_vars, "haplogrouper_unique_variants.tsv")
    
    return AnalysisResult(tree, leaf_vars, leaf_defn_vars, node_vars, node_defn_vars)
end


"""
    assign_samples(msa_file::String, newick_file::String, scaffold::String, ref_id::String)

Assign new samples to nodes on a previously defined scaffold tree

# Arguments
- `msa_file::String`: Path to multiple sequence alignment file in FASTA format
- `newick_file::String`: Path to the Newick tree file
- `scaffold::String`: Produced by Haplogrouper::define_phylogeny()
- `ref_id::String`: ID of the reference sequence in the MSA
"""
function assign_samples(msa_file::String, newick_file::String, scaffold::String, ref_id::String)
    all_node_vars = read_variants_from_tsv(scaffold)

    tree = parse_newick(newick_file)
    assignments = assign_samples_from_msa(tree, all_node_vars, msa_file, ref_id)

    println("\nAssignment results:")
    for (sample, node) in assignments
        println("$sample -> $node")
    end

    return nothing
end


end # module