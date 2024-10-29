# SelectiveD

The tool is designed to calculate Patterson's D (ABBA-BABA) test on selected subsets of positions. By comparison of regular D, D synonymous and D nonsynonymous we aim to estimate effects of natural selection on introgressed material.

##Input:

- VCF file (can be multi-species)
- Pylogeny (newick format)
- Name of outgroup for species quartets, used in D test

##Concept:

Script relies on codon structure of provided VCF file. Codons for each species quartet are analyzed independently: the codon only participates in counting ABBA and BABA patters, if it is conservative by two positions in all four species. Synonymous ABBA and BABA patterns are counted only in four fold degenerate codons, nonsynonymous - only in non degenerate codons. In calculation of regular D all positions are used. For time efficiency, after reading a line of VCF file, script handles all possible species quartets. For testing statistical significance of D tests we recommend to run SelectiveD on VCF file divided into n parts, and use the resulting values of D tests for further jacknifing.  

