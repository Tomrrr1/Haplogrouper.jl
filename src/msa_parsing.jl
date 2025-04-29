function check_sample_names(tree::RecursiveTree, seq_dict::Dict{String, BioSequence{DNAAlphabet{4}}})
    msa_samples = collect(keys(seq_dict))
    tree_samples = collect(Phylo.getleafnames(tree))
    
    msa_only = filter(sample -> sample ∉ tree_samples, msa_samples)
    tree_only = filter(leaf -> leaf ∉ msa_samples, tree_samples)
    
    if isempty(msa_only) && isempty(tree_only)
        println("All samples match between MSA and tree ($(length(msa_samples)) samples)")
        return true
    else
        if !isempty(msa_only)
            error("Found $(length(msa_only)) samples in MSA that are not in tree: $msa_only")
        end
        if !isempty(tree_only)
            error("Found $(length(tree_only)) leaves in tree that are not in MSA: $tree_only")
        end
    end
end

function check_msa_length(seq_dict::Dict{String, BioSequence{DNAAlphabet{4}}})
    # Check for length consistency
    lengths = [length(seq) for seq in values(seq_dict)]
    if !all(x -> x == lengths[1], lengths)
        for (name, seq) in seq_dict
            println("$name: $(length(seq)) bases")
        end
        error("Not all sequences have the same length in the MSA")
    end
    return true
end

function parse_msa(filename::String, ref_id::String, tree::RecursiveTree)
    records = collect(FASTA.Reader(open(filename, "r")))
    seq_dict = Dict{String, BioSequence{DNAAlphabet{4}}}()
    for record in records
        seq_dict[identifier(record)] = sequence(LongDNA{4}, record)
    end

    @argcheck haskey(seq_dict, ref_id) "Reference '$ref_id' not found in the MSA file"
    check_msa_length(seq_dict)
    check_sample_names(tree, seq_dict) # confirm newick and msa have same samples

    ref_seq = seq_dict[ref_id]
    all_samples = collect(keys(seq_dict))
    println("Samples found: $all_samples")
    
    # Identify variants for each sample compared to the reference
    leaf_vars = Dict(s => Set{Variant}() for s in all_samples)

    total_variants = 0
    for sample in all_samples
        sample == ref_id && continue
        variant_count = 0
        for base in eachindex(ref_seq)
            if ref_seq[base] == DNA_N || seq_dict[sample][base] == DNA_N
                continue
            end
            if seq_dict[sample][base] != ref_seq[base]
                variant = (base, ref_seq[base], seq_dict[sample][base])
                push!(leaf_vars[sample], variant)
                variant_count += 1
                total_variants += 1
            end
        end
    end
    println("Total non-reference variants: $total_variants")
    
    return leaf_vars
end