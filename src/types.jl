# Type definitions
const Variant = Tuple{Int, DNA, DNA}  # pos, ref, alt
const VariantSet = Set{Variant}
const VariantMap = Dict{String, VariantSet}  # node -> variants

struct AnalysisResult
    tree::RecursiveTree
    leaf_vars::VariantMap
    leaf_defn_vars::VariantMap
    node_vars::VariantMap 
    node_defn_vars::VariantMap
end