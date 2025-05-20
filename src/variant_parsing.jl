"""
    read_variants_from_tsv(filename::String, variants_column::String="Defining_Variants")

Read a TSV file containing all leaf- and node-defining variants into a VariantMap object.

# Arguments
- `filename::String`: Path to the TSV file

# Returns
- `VariantMap`: A dictionary mapping node names to sets of variants
"""
function read_variants_from_tsv(filename::String)
    if !isfile(filename)
        error("File not found: $filename")
    end

    df = CSV.read(filename, DataFrame, delim="\t")
    
    if !("Node" in names(df))
        error("Required column 'Node' not found in $filename")
    end
    if !("Defining_Variants" in names(df))
        error("Required column 'Defining_Variants' not found in $filename")
    end
    
    vars_dict = VariantMap()
    for (i, node_name) in enumerate(df[!, "Node"])
        vars_dict[node_name] = VariantSet()
        vars_str = df[i, "Defining_Variants"] # all variants for a given node
        if ismissing(vars_str) || isempty(vars_str)
            continue
        end
        
        vars_list = split(vars_str, ',')
        for v in vars_list
            m = match(r"^(\d+):([-a-zA-Z])>([-a-zA-Z])$", v) # Match pos:ref>alt
            if isnothing(m)
                @warn "Could not parse variant: $var_str for node: $node_name"
                continue
            end
            
            pos = parse(Int64, m[1])
            ref = LongDNA{4}(m[2])[1] # Extract single DNA character of type DNA
            alt = LongDNA{4}(m[3])[1]
            vars = Variant((pos, ref, alt))
            push!(vars_dict[node_name], vars)
        end
    end
    
    return vars_dict
end