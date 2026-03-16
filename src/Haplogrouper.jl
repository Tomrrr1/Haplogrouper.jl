module Haplogrouper

using ArgCheck
using BioSequences
using DataFrames
using Documenter
using CSV
using FASTX
using Phylo

include("types.jl")
include("msa_utils.jl")
include("msa.jl")
include("phylo_utils.jl")
include("variant_analysis_ancestral.jl")
include("variant_analysis_outgroup.jl")
include("io.jl")
include("utils.jl")
include("annotate.jl")
include("classify.jl")
include("scoring.jl")
include("fitch.jl")

export annotate
export classify

end

# Mutations are numbered according to a reference sequence, but they are not relative to the reference.
# The variants are relative to the reconstructed ancestral sequence.
# We use a reference to ensure stable coordinates, but the actual variants are relative to the root sequence.
# There are multiple equally probable most parsimonious trees.
# Classification. Given a new sample we identify the variants relative to the ancestral sequence but using the reference
# for coordinates.

# We have two modes: ancestral and outgroup. Ancestral mode involves reconstructing ancestral states 
# using the Fitch algorithm and identifying defining variants based on these states.
# Outgroup mode uses an outgroup to call variants. If all members of a clade have a variant relative
# to the outgroup then the branch ancestral to this clade will have this variant as defining.