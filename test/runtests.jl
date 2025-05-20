using Haplogrouper
using Test

@testset "Haplogrouper.jl" begin
    msa_path = "/Users/tomroberts/Documents/Oxford_DPhil/second_rotation/chicken_practice/sequence_files/miao_mtdna_seqs_1-61_msa2.fasta"
    newick_path = "/Users/tomroberts/Documents/Oxford_DPhil/second_rotation/chicken_practice/full_mitogenome_iqtree/full_mtdna_1-61.treefile"
    define_phylogeny(msa_path, newick_path, "NC_007235")
end
