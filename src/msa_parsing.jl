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

function parse_msa(filename::String, ref_id::String, tree::RecursiveTree, classify::Bool)
    records = collect(FASTA.Reader(open(filename, "r")))
    seq_dict = Dict{String, BioSequence{DNAAlphabet{4}}}()
    for record in records
        seq_dict[identifier(record)] = sequence(LongDNA{4}, record)
    end

    @argcheck haskey(seq_dict, ref_id) "Reference '$ref_id' not found in the MSA file"
    check_msa_length(seq_dict)
    # Only check tree consistency if not in classification mode
    if !classify
        check_sample_names(tree, seq_dict) # confirm newick and msa have same samples
    end

    ref_seq = seq_dict[ref_id]
    #println(ref_seq)
    all_samples = collect(keys(seq_dict))
    println("Samples found: $all_samples")

    # Create a mapping from alignment positions to reference coordinates
    # This will give us stable coordinates regardless of alignment gaps
    ref_coord_map = Int[]
    ref_coord = 0
    
    for (_, base) in enumerate(ref_seq)
        if base != DNA_Gap  # Skip gaps in reference
            ref_coord += 1
            push!(ref_coord_map, ref_coord)
        else
            push!(ref_coord_map, ref_coord)  # Gap positions map to previous reference coordinate
        end
    end
    
    # Identify variants for each sample compared to the reference
    leaf_vars = Dict(s => VariantSet() for s in all_samples)
    total_variants = 0
    for sample in all_samples
        sample == ref_id && continue
        variant_count = 0
        
        for align_pos in eachindex(ref_seq)
            # Skip positions with N or gaps in either sequence
            if ref_seq[align_pos] == DNA_N || seq_dict[sample][align_pos] == DNA_N
               # || ref_seq[align_pos] == DNA_Gap || seq_dict[sample][align_pos] == DNA_Gap
                continue
            end
            
            if seq_dict[sample][align_pos] != ref_seq[align_pos]
                # Use reference-based coordinate instead of alignment position
                ref_pos = ref_coord_map[align_pos]
                variant = (ref_pos, ref_seq[align_pos], seq_dict[sample][align_pos])
                push!(leaf_vars[sample], variant)
                variant_count += 1
                total_variants += 1
            end
        end
    end
    println("Total non-reference variants: $total_variants")
    
    return leaf_vars
end