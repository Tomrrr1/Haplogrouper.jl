# Parse MSA file for tree annotation
function parse_msa_annotate(filename::String, ref_id::String, tree::RecursiveTree)
    records = collect(FASTA.Reader(open(filename, "r")))
    seq_dict = Dict{String, BioSequence{DNAAlphabet{4}}}()
    for record in records
        seq_dict[identifier(record)] = sequence(LongDNA{4}, record)
    end
    check_msa_length(seq_dict)
    check_sample_names(tree, seq_dict)

    return build_nucleotide_maps(seq_dict, ref_id)
end

# Method 1: Parse MSA file for classification (ancestral mode)
function parse_msa_classify(filename::String, ref_id::String, ancestral_state::Dict{Int, DNA})
    records = collect(FASTA.Reader(open(filename, "r")))
    seq_dict = Dict{String, BioSequence{DNAAlphabet{4}}}()
    for record in records
        seq_dict[identifier(record)] = sequence(LongDNA{4}, record)
    end

    if !haskey(seq_dict, ref_id)
        error("Reference '$ref_id' not found in the MSA file")
    end 
    check_msa_length(seq_dict)
    # Find bases that differ from the ancestral sequence
    nucleotide_maps = build_nucleotide_maps(seq_dict, ref_id)
    sample_vars = VariantMap()
    for (sample, tuples) in nucleotide_maps
        sample == ref_id && continue
        sample_vars[sample] = VariantSet()
        for (pos, _, sample_base) in tuples
            anc_base = get(ancestral_state, pos, nothing)
            if anc_base !== nothing && sample_base != anc_base &&
               sample_base != DNA_N && sample_base != DNA_Gap &&
               anc_base != DNA_N && anc_base != DNA_Gap
                push!(sample_vars[sample], (pos, anc_base, sample_base))
            end
        end
    end

    return sample_vars
end

# Method 2: Parse MSA file for classification (outgroup mode)
function parse_msa_classify(filename::String, ref_id::String)
    records = collect(FASTA.Reader(open(filename, "r")))
    seq_dict = Dict{String, BioSequence{DNAAlphabet{4}}}()
    for record in records
        seq_dict[identifier(record)] = sequence(LongDNA{4}, record)
    end

    if !haskey(seq_dict, ref_id)
        error("Reference '$ref_id' not found in the MSA file")
    end 
    check_msa_length(seq_dict)

    nucleotide_maps = build_nucleotide_maps(seq_dict, ref_id)
    sample_vars = VariantMap()
    for (sample, tuples) in nucleotide_maps
        sample == ref_id && continue  # Skip the outgroup itself
        sample_vars[sample] = VariantSet()
        for (pos, ref_base, sample_base) in tuples
            if sample_base != ref_base &&
               sample_base != DNA_N && sample_base != DNA_Gap &&
               ref_base != DNA_N && ref_base != DNA_Gap
                push!(sample_vars[sample], (pos, ref_base, sample_base))
            end
        end
    end

    return sample_vars
end