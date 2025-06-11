# Scoring functions for sample classification

function jaccard(sample_vars::VariantSet, node_vars::VariantSet)
    intersect_vars = intersect(sample_vars, node_vars) # shared variants [1,2] [1,2,3] => [1,2]
    union_vars = union(sample_vars, node_vars) # [1,2] [1,2,3] => [1,2,3]
    
    return length(intersect_vars) / length(union_vars)
end 

function kulczynski(sample_vars::VariantSet, node_vars::VariantSet)
    intersect_vars = intersect(sample_vars, node_vars) # shared variants [1,2] [1,2,3] => [1,2]
    ratio1 = length(intersect_vars) / length(node_vars)
    ratio2 = if isempty(sample_vars) 
        0
    else 
        length(intersect_vars) / length(sample_vars) 
    end
    
    return 0.5 * (ratio1 + ratio2)
end 