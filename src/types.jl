# Type definitions
const Variant = Tuple{Int, DNA, DNA}  # position, reference allele, alternate allele
const VariantSet = Set{Variant}
const VariantMap = Dict{String, VariantSet}  # sample/node -> variants

struct AnalysisResult
    tree::RecursiveTree
    leaf_vars::VariantMap
    leaf_defn_vars::VariantMap
    node_defn_vars::VariantMap
end