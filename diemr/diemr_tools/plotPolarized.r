#' Plots polarized genotypes
#'
#' Plots genotypes that can be optionally polarized.
#' @param genotypes character matrix comprising of _012 encodings.
#' @param HI numeric vector of individual hybrid indices with length equal to number
#'   of rows in \code{genotypes}.
#' @inheritParams diem
#' @param cols vector of four colours, representing missing data, homozygots for 
#'    genotype 0, heterozygots and homozygots for genotype 2.
#' @details To import and polarize genotypes, use the function \code{\link{importPolarized}}.
#'
#'   When using \code{\link{diem}} with argument \code{verbose = TRUE}, hybrid indices,
#'   \code{HI}, can be found in file 'HIwithOptimalPolarities.txt' in folder 'diagnostics'
#'   in the working directory.
#'
#' @return No return value, called for side effects. In the default plot, purple and green 
#'   represent side of the barrier to geneflow encoded as `0` and `2`, respectively, 
#'   yellow shows heterozygots and white missing or undetermined genotypes. Individuals 
#'   are ordered according to the \code{HI}.
#' @importFrom graphics image
#' @importFrom grDevices col2rgb
#' @export
#' @examples
#' gen <- importPolarized(
#'   file = system.file("extdata", "data7x10.txt", package = "diemr"),
#'   changePolarity = c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE),
#'   ChosenInds = 1:7
#' )
#'
#' h <- c(0.625, 0.5, 0.455, 0.455, 0.227, 0.818, 0.292)
#'
#' plotPolarized(genotypes = gen, HI = h)

plotPolarized <- function(genotypes, HI, cols = c("#FFFFFF", "#800080", "#FFE500", "#008080"), ...) {

  nMarkers <- ncol(genotypes)
  nInds <- nrow(genotypes)
  
  # check colors
  areColors <- sapply(cols, function(X) tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE))
  if(sum(areColors) != 4){
  	warning("Argument cols must contain four valid colours. Review the input, now using default colours")
  	print(areColors)
  	cols = c("#FFFFFF", "#800080", "#FFE500", "#008080")
  }
  
  image(
    x = 1:nMarkers,
    y = 1:nInds,
    z = t(matrix(as.numeric(factor(genotypes[order(HI), ], 
                                   levels = c("_", "0", "1", "2"))), 
                 ncol = nMarkers)),
    col = cols,
    xlab = "Markers", ylab = "Individuals", axes = FALSE, useRaster = TRUE
  )
}
