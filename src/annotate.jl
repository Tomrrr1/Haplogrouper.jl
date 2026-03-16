"""
    annotate(msa_file::String, newick_file::String, ref_id::String; method::String="ancestral", prefix::String="")

Identify defining variants for each node on a phylogenetic tree.

# Arguments
- `msa_file::String`: Path to multiple sequence alignment file in FASTA format.
- `newick_file::String`: Path to the tree file in Newick format.
- `ref_id::String`: Sequence ID of the reference sample. In "outgroup" mode this sample must be the outgroup.
- `mode::String`: Mode to use for tree annotation. Options: "ancestral" (default) or "outgroup".
- `prefix::String`: Prefix for output filenames. Default is "".

# Returns
- `nothing`

# Example
```jldoctest
julia> msa_string = ">ref\\nATCGATCG\\n>s1\\nATGGATCT\\n>s2\\nATGGATGG";

julia> tree_string = "((s1, s2), ref);";

julia> write("sequences.fasta", msa_string);

julia> write("tree.nwk", tree_string);

julia> annotate("sequences.fasta", "tree.nwk", "ref")

All samples match between MSA and tree (3 samples)
Using ancestral sequence reconstruction method
Results written to TSV file: haplogrouper_all_variants.tsv
```
"""
function annotate(msa_file::String, newick_file::String, ref_id::String; mode::String="ancestral", prefix::String="")
    tree = parse_newick(newick_file)
    nucleotide_maps = parse_msa_annotate(msa_file, ref_id, tree)
    
    # Generate output filename with prefix
    variants_filename = isempty(prefix) ? "haplogrouper_all_variants.tsv" : "$(prefix)_haplogrouper_all_variants.tsv"
    
    if mode == "ancestral"
        node_vars = define_nodes_ancestral(tree, nucleotide_maps, ref_id, prefix)
        println("Using ancestral sequence reconstruction method")
    elseif mode == "outgroup"
        node_vars = define_nodes_outgroup(tree, nucleotide_maps, ref_id)
        println("Using outgroup reference method with outgroup: $ref_id")
    else
        error("Unknown mode: $mode. Use 'ancestral' or 'outgroup'")
    end
    
    output_variants_tsv(tree, node_vars, variants_filename)
    
    return nothing
end