# Confirm that the sample names in the msa match those in the tree
function check_sample_names(tree::RecursiveTree, seq_dict::Dict{String, BioSequence{DNAAlphabet{4}}})
    msa_samples = collect(keys(seq_dict))
    tree_samples = collect(Phylo.getleafnames(tree))
    
    msa_only = filter(sample -> sample ∉ tree_samples, msa_samples)
    tree_only = filter(leaf -> leaf ∉ msa_samples, tree_samples)
    
    if isempty(msa_only) && isempty(tree_only)
        println("\nAll samples match between MSA and tree ($(length(msa_samples)) samples)")
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

# Confirm that all sequences in the msa have the same length
function check_msa_length(seq_dict::Dict{String, BioSequence{DNAAlphabet{4}}})
    lengths = [length(seq) for seq in values(seq_dict)]
    if !all(x -> x == lengths[1], lengths)
        for (name, seq) in seq_dict
            println("$name: $(length(seq)) bases")
        end
        error("Not all sequences have the same length in the MSA")
    end
    
    return true
end

# Ensure a stable coordinate system anchored on the designated reference
function build_coord_vec(seq_dict::Dict{String, BioSequence{DNAAlphabet{4}}}, ref_id::String)
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
    ref_coord_map = build_coord_vec(seq_dict, ref_id)
    ref_seq = seq_dict[ref_id]
    # For each base in each sample record the ref-based position, ref base, and sample base
    for (sample, _) in seq_dict
        nucleotide_maps[sample] = VariantSet()
        for align_pos in eachindex(ref_seq)
            base = (ref_coord_map[align_pos], ref_seq[align_pos], seq_dict[sample][align_pos])
            push!(nucleotide_maps[sample], base)
        end
    end
    
    return nucleotide_maps
end