#' Hybrid index, heterozygosity, error rate
#'
#' Using genotype allele counts, calculates
#' the hybrid index,
#' heterozygosity and error rate in a single individual.
#'
#' @param sCount a numeric vector of length 4 with allele counts for missing data,
#'               homozygots for allele 1, heterozygots, and homozygots for allele 2.
#' @return Returns a named numeric vector with three values: p - hybrid index,
#'         Het - heterozygosity,
#'         Err - error rate.
#' @details Allele counts are genomic state counts multiplied by ploidy. As different
#'   compartments might have different ploidies (e.g. autosomal markers, sex chromosomes,
#'   mitochondrial markers), allele counts should be calculated per compartment and then
#'   summarised to obtain the correct genomic allele counts.   
#'   When all individuals in each compartmenst have the same ploidy, state counts do 
#'   not need to be corrected. 
#' @examples
#' pHetErrOnStateCount(sCount = c(2, 4, 2, 6))
#' #         p       Het       Err
#' # 0.5833333 0.1666667 0.1428571
#' @export
#' @importFrom stats setNames

pHetErrOnStateCount <- function(sCount) {

  # vectorize the data
  if (inherits(sCount, "matrix")) {
    sCount <- c(sCount)
  }
  if (inherits(sCount, "data.frame")) {
    sCount <- unlist(sCount)
  }

  # check data
  if (length(sCount) != 4) {
    stop("sCount length is not 4")
  }
  if (any(!is.numeric(sCount))) {
    stop("sCount contains non numeric data")
  }
  if (any(is.na(sCount))) {
    stop("sCount contains NA")
  }

  # total number of characters with non-missing data
  CallTotal <- sum(sCount * c(0, 1, 1, 1))

  # check if any genotypes were called
  if (CallTotal == 0) {
    Err <- ifelse(sum(sCount) > 0, sCount[1] / sum(sCount), NA)
    res <- c(NA, NA, Err)
  }

  # calculate indices
  if (CallTotal > 0) {
    p <- sum(sCount * c(0, 0, 1, 2)) / (2 * CallTotal)
    Het <- sCount[3] / CallTotal
    Err <- sCount[1] / sum(sCount)
    res <- c(p, Het, Err)
  }

  return(setNames(object = res, nm = c("p", "Het", "Err")))
}
