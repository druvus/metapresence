"""
Metapresence package for evaluating the distribution of aligned reads onto nucleotidic sequences.

While Metapresence was developed for species identification in metagenomic data,
it can be used to assess read distribution on any DNA/RNA sequences from alignment results.
"""

__version__ = "1.0.0"

# We're not importing anything here to avoid circular imports
# The CLI and other modules will import directly from the appropriate submodules