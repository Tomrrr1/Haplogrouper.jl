# Haplogrouper

[![Build Status](https://github.com/Tomrrr1/Haplogrouper.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Tomrrr1/Haplogrouper.jl/actions/workflows/CI.yml?query=branch%3Amain)

`Haplogrouper.jl` is a Julia package for identifying genetic variants that define clades in phylogenetic trees (i.e. haplogroups). The package also provides functionality for classifying new samples.

## Installation

The package can currently be installed from GitHub using the following commands
```
using Pkg
Pkg.add(url="https://github.com/Tomrrr1/Haplogrouper.jl")
```

## Example usage

The package provides two main functions for *1)* building a scaffold tree, and *2)* classifying new samples based on this tree.

```
using Haplogrouper

# Identify all defining variants of a phylogeny relative to a reference
msa_file = "/path/to/msa"
tree_file = "/path/to/tree"
ref_id = "Ref"

define_phylogeny(msa_file, tree_file, ref_id)
# This function outputs the file "haplogrouper_all_variants.tsv" which will be used by assign_samples()

# Classify new samples 
new_samples_msa = "/path/to/new_msa"
scaffold = "/path/to/scaffold" # File produced by define_phylogeny() 

assign_samples(new_samples_msa, tree_file, scaffold, ref_id)
```
