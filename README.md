# Haplogrouper

[![Build Status](https://github.com/Tomrrr1/Haplogrouper.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Tomrrr1/Haplogrouper.jl/actions/workflows/CI.yml?query=branch%3Amain)

`Haplogrouper.jl` is a Julia package for identifying genetic variants that define clades in phylogenetic trees (i.e. haplogroups). The package also provides functionality for assigning new samples to previously defined haplogroups.


Let different branches share the same variant if they are not correlated. Two different nodes can share the same variant if they are not parents/children.
