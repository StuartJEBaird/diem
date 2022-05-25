# diem: Diagnostic index expectation maximisation

The *diem* is a deterministic, likelihood-based genome polarisation algorithm that finds which alleles of genomic markers belong to which side of the barrier. Co-estimates which individuals belong to either side of the barrier and barrier strength. Uses expectation maximisation in likelihood framework.

The *diem* algorithm was implemented in different programming languages, with core functions and accompanying tools stored in this repository. 

1. **diemr** - `R` implementation. The package is available through https://CRAN.R-project.org/package=diemr.
2. **diempy** - `Python` implementation.
3. **diemmca** - `Mathematica` implementation.

When using any of the *diem* algorithm implementations, cite:

> Stuart J. E. Baird, Jan Petružela, Izar Jaroň, Pavel Škrabánek, Natália Martínková. 2022. diemr: Genome polarisation for detecting barriers to geneflow in R. bioRxiv 2022.03.24.485605; doi: https://doi.org/10.1101/2022.03.24.485605.
