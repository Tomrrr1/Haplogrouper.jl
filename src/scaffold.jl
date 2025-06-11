"""
    make_scaffold(msa_file::String, newick_file::String, ref_id::String)

Identify defining variants for each node on a phylogenetic tree.

# Arguments
- `msa_file::String`: Path to multiple sequence alignment file in FASTA format.
- `newick_file::String`: Path to the tree file in Newick format.
- `ref_id::String`: Sequence ID of the sample reference sample. 

# Returns
- `nothing`

# Example
```jldoctest
julia> msa_string = ">ref\\nATCGATCG\\n>s1\\nATGGATCT\\n>s2\\nATGGATGG";

julia> tree_string = "((s1, s2), ref);";

julia> write("sequences.fasta", msa_string);

julia> write("tree.nwk", tree_string);

julia> make_scaffold("sequences.fasta", "tree.nwk", "ref")

All samples match between MSA and tree (3 samples)
Results written to TSV file: haplogrouper_all_variants.tsv
```
"""
function make_scaffold(msa_file::String, newick_file::String, ref_id::String)
    tree = parse_newick(newick_file)
    nucleotide_maps = parse_msa_scaffold(msa_file, ref_id, tree)
    node_vars = define_nodes(tree, nucleotide_maps, ref_id)
    output_variants_tsv(tree, node_vars, "haplogrouper_all_variants.tsv")
    
    return nothing
end