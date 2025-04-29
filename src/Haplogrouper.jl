module Haplogrouper

# Import external packages
using Phylo
using BioSequences
using FASTX
using ArgCheck

# Type definitions
const Variant = Tuple{Int, DNA, DNA}  # position, reference allele, alternate allele
const VariantSet = Set{Variant}
const VariantMap = Dict{String, VariantSet}  # sample/node -> variants

# Include all module files
include("types.jl")
include("tree_parsing.jl") 
include("msa_parsing.jl")
include("variant_analysis.jl")
include("output.jl")
include("utils.jl")

export AnalysisResult
export parse_newick
export parse_msa
export define_leaves
export define_nodes
export define_variants
export define_phylogeny
export output_results_tsv

"""
    define_phylogeny(msa_file::String, newick_file::String, ref_id::String)

Main workflow function to identify defining variants for each branch on a phylogeny

# Arguments
- `msa_file::String`: Path to multiple sequence alignment file in FASTA format
- `newick_file::String`: Path to the Newick tree file
- `ref_id::String`: ID of the reference sequence in the MSA

# Returns
- `AnalysisResult`: Structure containing all analysis results
"""
function define_phylogeny(msa_file::String, newick_file::String, ref_id::String)
    # Parse the newick tree
    tree = parse_newick(newick_file)
    
    # Parse the MSA and identify variants
    leaf_vars = parse_msa(msa_file, ref_id, tree)
    
    # Find defining variants for leaves
    leaf_defn_vars = define_leaves(leaf_vars)
    
    # Find defining variants for internal nodes
    node_defn_vars = define_nodes(tree, leaf_vars)
    
    # Return the analysis result
    return AnalysisResult(tree, leaf_vars, leaf_defn_vars, node_defn_vars)
end

end # module