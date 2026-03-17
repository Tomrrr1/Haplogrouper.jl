# Haplogrouper.jl Documentation

**Haplogrouper** is a [Julia](http://www.julialang.org) package for *(a)* annotating phylogenetic trees with genetic variants that are unique to each node, and *(b)* classifying new samples based on these identified variants.

Haplogrouper is an alternative to the programme mtPhyl, which, despite its utility, lacks clear documentation, is incompatible with MacOS, and is difficult to obtain (we can only access it through the Wayback Machine as of 25/07/25).

## Installation

The package is currently available through GitHub and can be installed with the following commands.

```julia
using Pkg
Pkg.add(url="https://github.com/tomrrr1/Haplogrouper.jl")
```

## Using Haplogrouper

Haplogrouper provides two main functions for phylogenetic analysis:

### `annotate()`
Annotates a phylogenetic tree with the genetic variants that are unique to each node. This function supports two modes:
- *Ancestral*: Uses Fitch's algorithm to reconstruct ancestral sequences at internal nodes. The variants defining each node are called relative to the ancestral root sequence.
- *Outgroup*: Variants are called relative to the chosen outgroup/reference sequence.

### `classify()`
Classifies new samples by comparing their variants to the variants at each node of a previously annotated tree. The function uses scoring metrics (Jaccard or Kulczynski) to assign samples to the most appropriate node.

The *Jaccard index* penalises private variants (variants unique to either the sample or node) because it normalises the number of shared variants by the total number of variants across both sets. The Jaccard index is defined as:

```math
\text{Jaccard}(A, B) = \frac{|A \cap B|}{|A \cup B|}
```

Alternatively, the *Kulczynski measure* normalises the number of shared variants by the number of variants in each set individually:

```math
\text{Kulczynski}(A, B) = \frac{1}{2} \left( \frac{|A \cap B|}{|A|} + \frac{|A \cap B|}{|B|} \right)
```

This measure is more robust to the presence of private variants because it accounts for asymmetries in set size.

## Main functions
```@docs
Haplogrouper.annotate
Haplogrouper.classify
```