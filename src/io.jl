# Function to output node variants to a tsv file
function output_variants_tsv(tree::RecursiveTree, node_vars::VariantMap, output_file::String)
    open(output_file, "w") do file
        write(file, "Node\tDescendants\tDefining_Variants\n")
        # Process nodes
        for node in keys(node_vars)
            descendants = collect_descendant_leaves(tree, node)
            descendants_str = join(descendants, ",")
            
            variants = get(node_vars, node, VariantSet())
            variant_strs = ["$(v[1]):$(v[2])>$(v[3])" for v in variants]
            variants_str = join(variant_strs, ",")
            
            write(file, "$node\t$descendants_str\t$variants_str\n")
        end
    end

    println("Results written to TSV file: $output_file")
    return nothing
end

function read_variants_tsv(filename::String)
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
            m = match(r"^(\d+):([-a-zA-Z])>([-a-zA-Z])$", String(v)) # Match pos:ref>alt
            if isnothing(m)
                @warn "Could not parse variant: $v for node: $node_name"
                continue
            end
            
            pos = parse(Int64, m[1])
            ref = LongDNA{4}(m[2])[1] # Extract single character of type DNA
            alt = LongDNA{4}(m[3])[1]
            push!(vars_dict[node_name], (pos, ref, alt))
        end
    end
    
    return vars_dict
end

# Function to write the inferred ancestral sequence to a tsv file
function write_ancestral_sequence_tsv(ancestral_states::Dict{String, Dict{Int, DNA}}, root_name::String, outfile::String)
    root_sequence = ancestral_states[root_name]
    open(outfile, "w") do io
        for pos in sort(collect(keys(root_sequence)))
            println(io, "$pos\t$(root_sequence[pos])")
        end
    end

    return nothing
end

# Function to read an ancestral sequence from a tsv file into a Dict{pos, base}
function read_ancestral_sequence_tsv(filename::String)
    ancestral_state = Dict{Int, DNA}()
    open(filename, "r") do io
        for line in eachline(io)
            pos, base = split(line, '\t')
            ancestral_state[parse(Int, pos)] = LongDNA{4}(base)[1]
        end
    end
    
    return ancestral_state
end