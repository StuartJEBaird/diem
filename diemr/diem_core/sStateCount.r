#' Count states in a vector
#'
#' Counts genomic states in one sample.
#'
#' @param s character vector with elements "_", "0", "1", "2" representing
#' 			 missing data, homozygots for allele 1, heterozygots, and
#'           homozygots for allele 2. The vector should represent a single individual.
#' @return Numeric vector of length 4 with counts of "_", "0", "1", "2" respectively.
#'
#' @details Summarizes the number of markers that are fixed for an allele in the genome of
#'   one individual. This is used to construct the I4 matrix in \code{\link{diem}}.
#' @seealso \code{\link{emPolarise}} for changing marker polarity.
#' @export
#' @examples
#' genotype <- c("0", "0", "_", "2", "1", "0", "1")
#' sStateCount(genotype)
#' # [1] 1 3 2 1
#'
#' # calculate state counts for a polarised genotype
#' sStateCount(emPolarise(genotype, TRUE))
#' # [1] 1 1 2 3
sStateCount <- function(s) {
  return(c(
    sum(s == "_"),
    sum(s == "0"),
    sum(s == "1"),
    sum(s == "2")
  ))
}
