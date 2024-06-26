# *diem*: Diagnostic index expectation maximisation

<img
  src="/diemimages/CircleDiagram24MayBitMap.jpg"
  alt="Alt text"
  title=""
  style="display: inline-block; margin: 0 auto; max-width: 300px">
  
*diem* is a genome polarisation algorithm that finds which alleles of genomic markers belong to which side of a barrier, co-estimating which individuals belong on either side. *diem* is a deterministic expectation maximisation algorithm and so uses the likelihood inference framework.

The *diem* algorithm has been implemented in several programming languages, with core functions and accompanying tools stored in this repository. 

1. **diemr** - `R` implementation. The package is available through https://CRAN.R-project.org/package=diemr. ![](https://cranlogs.r-pkg.org/badges/diemr)
2. **diempy** - `Python` implementation.
3. **diemmca** - `Mathematica` implementation. Soon also to be available at https://www.notebookarchive.org

When using any of the *diem* algorithm implementations, cite:

> Stuart J. E. Baird, Jan Petružela, Izar Jaroň, Pavel Škrabánek, Natália Martínková. 2022. Genome polarisation for detecting barriers to geneflow. bioRxiv 2022.03.24.485605; doi: https://doi.org/10.1101/2022.03.24.485605.

> Stuart J. E. Baird, Jan Petružela, Izar Jaroň, Pavel Škrabánek, Natália Martínková. 2023. Genome polarisation for detecting barriers to geneflow. Methods in Ecology and Evolution; doi: https://doi.org/10.1111/2041-210X.14010.
