module Haplogrouper

using ArgCheck
using BioSequences
using DataFrames
using Documenter
using CSV
using FASTX
using Phylo

include("types.jl")
include("msa_parsing.jl")
include("variant_analysis.jl")
include("io.jl")
include("utils.jl")
include("scaffold.jl")
include("classify.jl")
include("scoring.jl")
include("reconstruct.jl")

export classify
export make_scaffold

end

# Mutations are numbered according to a reference sequence, but they are not relative to the reference.
# The variants are relative to the reconstructed ancestral sequence.
# We use a reference to ensure stable coordinates, but the actual variants are relative to the root sequence.
# There are multiple equally probable most parsimonious trees.
# Classification. Given a new sample we identify the variants relative to the ancestral sequence but using the reference
# for coordinates.