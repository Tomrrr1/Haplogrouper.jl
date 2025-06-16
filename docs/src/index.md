# Haplogrouper.jl Documentation

**Haplogrouper** is a [Julia](http://www.julialang.org) package for *(a)* identifying genetic variants that define each node in a phylogenetic tree, and *(b)* classifying new samples based on these identified variants.

Haplogrouper is an alternative to the programme mtPhyl, which, despite its utility, lacks clear documentation, is incompatible with MacOS, and is difficult to obtain (we could only access it through the Wayback Machine).

## Installation

The package is currently available through GitHub and can be installed with the following commands.

```julia
using Pkg
Pkg.add(url="https://github.com/tomrrr1/Haplogrouper.jl")
```

## Main functions
```@docs
Haplogrouper.make_scaffold
Haplogrouper.classify
```