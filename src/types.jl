# Type definitions
const Variant = Tuple{Int, DNA, DNA}  # pos, ref, alt
const VariantSet = Set{Variant} # Set of variant tuples
const VariantMap = Dict{String, VariantSet}  # node -> variants