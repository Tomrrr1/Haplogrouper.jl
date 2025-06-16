# Ensure a stable coordinate system anchored on the designated reference
function build_coord_map(seq_dict::Dict{String, BioSequence{DNAAlphabet{4}}}, ref_id::String)
    ref_seq = seq_dict[ref_id]
    ref_coord_map = Vector{Int64}()
    ref_coord = 0
    for base in ref_seq
        if base != DNA_Gap  # Skip gaps in reference
            ref_coord += 1
            push!(ref_coord_map, ref_coord)
        else
            push!(ref_coord_map, ref_coord)  # Gap positions map to previous coordinate
        end
    end

    return ref_coord_map
end

# For each sample, record the base at each position using the reference-derived coordinates
function build_nucleotide_maps(seq_dict::Dict{String, BioSequence{DNAAlphabet{4}}}, ref_id::String)
    nucleotide_maps = VariantMap()
    ref_coord_map = build_coord_map(seq_dict, ref_id)
    ref_seq = seq_dict[ref_id]
    # For each base in each sample record the reference-based position, the reference base and the sample base
    for (sample, _) in seq_dict
        nucleotide_maps[sample] = VariantSet()
        for align_pos in eachindex(ref_seq)
            base = (ref_coord_map[align_pos], ref_seq[align_pos], seq_dict[sample][align_pos])
            push!(nucleotide_maps[sample], base)
        end
    end
    
    return nucleotide_maps
end

# Parse MSA file specific to scaffold mode
function parse_msa_scaffold(filename::String, ref_id::String, tree::RecursiveTree)
    records = collect(FASTA.Reader(open(filename, "r")))
    seq_dict = Dict{String, BioSequence{DNAAlphabet{4}}}()
    for record in records
        seq_dict[identifier(record)] = sequence(LongDNA{4}, record)
    end

    check_msa_length(seq_dict)
    check_sample_names(tree, seq_dict)

    return build_nucleotide_maps(seq_dict, ref_id)
end

# Parse MSA file specific to classification mode
function parse_msa_classify(filename::String, ref_id::String, ancestral_state::Dict{Int, DNA})
    records = collect(FASTA.Reader(open(filename, "r")))
    seq_dict = Dict{String, BioSequence{DNAAlphabet{4}}}()
    for record in records
        seq_dict[identifier(record)] = sequence(LongDNA{4}, record)
    end

    @argcheck haskey(seq_dict, ref_id) "Reference '$ref_id' not found in the MSA file"
    check_msa_length(seq_dict)

    ref_coord_map = build_coord_map(seq_dict, ref_id)
    ref_seq = seq_dict[ref_id]
    sample_vars = VariantMap()
    for (sample, seq) in seq_dict
        sample == ref_id && continue
        sample_vars[sample] = VariantSet()
        for align_pos in eachindex(ref_seq)
            pos = ref_coord_map[align_pos]
            anc_base = get(ancestral_state, pos, nothing)
            sample_base = seq[align_pos]
            if anc_base !== nothing && sample_base != anc_base
                push!(sample_vars[sample], (pos, anc_base, sample_base))
            end
        end
    end

    return sample_vars
end