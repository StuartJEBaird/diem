#' Polarises a marker
#'
#' Changes encodings of genomic markers according to user specification.
#'
#' @param changePolarity logical scalar, indicating whether to leave the marker as is 
#'     (\code{FALSE}) or whether to change its polarity (\code{TRUE}).
#' @param origM character vector of genotypes comprising of _012 encodings.
#' @return Returns a character vector with polarised markers.
#' @note Note that \code{\link{diem}} and \code{\link{importPolarized}} accept also a `U` 
#'   encoding for an unknown or third allele, but \code{emPolarise} requires all `U` to 
#'   be replaced with `_`.
#' @export
#' @seealso \code{\link{diem}} for determining appropriate marker polarity with
#'   respect to a barrier to geneflow.
#' @examples
#' emPolarise(c("0", "0", "1", "2", "2"), TRUE)
#' # [1] "2" "2" "1" "0" "0"
#'
#' emPolarise(c("0", "_", "2", "2", "1"), FALSE)
#' # [1] "0" "_" "2" "2" "1"
emPolarise <- function(origM, changePolarity = TRUE) {
  if (!inherits(origM, "character")) {
    stop("orgiM must be a character vector. It is now ", class(origM))
  }
  if (any(!(origM %in% c("_", "0", "1", "2")))) {
    stop("origM must contain only characters _, 0, 1, 2")
  }
  if (changePolarity) {
    polarisedM <- c("_", "2", "1", "0")[match(origM, c("_", "0", "1", "2"))]
  } else {
    polarisedM <- origM
  }
  return(polarisedM)
}
