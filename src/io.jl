# Function to output node variants to a tsv file
function output_variants_tsv(tree::RecursiveTree, node_vars::VariantMap, output_file::String)
    open(output_file, "w") do file
        write(file, "Node\tDescendants\tDefining_Variants\n")
        # Process nodes
        for node in keys(node_vars)
            descendants = collect_descendant_leaves(tree, node)
            descendants_str = join(descendants, ",")
            
            vars = get(node_vars, node, VariantSet())
            formatted_vars = ["$(v[1]):$(v[2])>$(v[3])" for v in vars]
            vars_str = join(formatted_vars, ",")
            
            write(file, "$node\t$descendants_str\t$vars_str\n")
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
    vars_dict = VariantMap()
    for (i, node_name) in enumerate(df[!, "Node"])
        vars_dict[node_name] = VariantSet()
        vars_str = df[i, "Defining_Variants"] # all variants for a given node
        if ismissing(vars_str) || isempty(vars_str)
            continue
        end
        
        vars_vec = split(vars_str, ',') # Split string at commas
        for var in vars_vec
            m = match(r"^(\d+):([-a-zA-Z])>([-a-zA-Z])$", String(var)) # Match pos:ref>alt
            if isnothing(m)
                @warn "Failed to parse variant: $var for node: $node_name"
                continue
            end
            pos = parse(Int64, m.captures[1])
            ref = LongDNA{4}(m.captures[2])[1] # Extract single character of type DNA
            alt = LongDNA{4}(m.captures[3])[1]
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

# Write sample -> node assignments to CSV
function output_assignments_csv(assignments::Dict{String, Vector{Vector{Any}}}, output_file::String)
    rows = NamedTuple{(:Sample, :Rank, :Node, :Score, :Shared_Count, :SampleSpecific_Count, :NodeSpecific_Count),
                      Tuple{String, Int, String, Float64, Int, Int, Int}}[]
    for (sample_id, results) in sort(collect(assignments); by=first)
        if isempty(results)
            push!(rows, (Sample=sample_id, Rank=1, Node="No classification", Score=0.0,
                         Shared_Count=0, SampleSpecific_Count=0, NodeSpecific_Count=0))
            continue
        end
        for (rank, r) in enumerate(results[1:min(end, 5)])
            node_name, shared_variants, unique_to_sample, unique_to_node, score = r
            push!(rows, (
                Sample=sample_id,
                Rank=rank,
                Node=String(node_name),
                Score=Float64(score),
                Shared_Count=length(shared_variants),
                SampleSpecific_Count=length(unique_to_sample),
                NodeSpecific_Count=length(unique_to_node),
            ))
        end
    end
    df = DataFrame(rows)
    CSV.write(output_file, df)
    println("Assignment results written to CSV file: $output_file")
    return nothing
end