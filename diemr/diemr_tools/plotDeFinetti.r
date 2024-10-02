#' Plot the De Finetti Diagram for Polarized Genotypes
#'
#' This function calculates genotype frequencies from polarized genotypes, ideally
#' imported using the
#' \code{importPolarized} function. It plots individuals onto a ternary De Finetti
#' diagram and includes a curve indicating Hardy-Weinberg equilibrium if specified.
#'
#' @inheritParams plotPolarized
#' @inheritParams importPolarized
#' @param cols A character vector of colors with a length equal to the number of
#'   individuals (rows) in \code{genotypes}.
#' @param HWE Logical indicating whether to plot the curve for Hardy-Weinberg Equilibrium.
#' @param tipLabels A character vector of length 3 with names for the ternary plot vertices.
#' @param ... additional graphical parameters (see \link[graphics]{plot.default}).
#' @details To import and polarize genotypes, use the \link{importPolarized} function. 
#'   Alternatively, the I4 matrix can be used as input for \code{genotypes}.
#'
#' @return No return value; the function is called for its side effects.
#' @importFrom utils modifyList
#' @importFrom graphics lines text
#' @export
#' @examples
#' gen <- importPolarized(
#'   file = system.file("extdata", "data7x10.txt", package = "diemr"),
#'   changePolarity = c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE),
#'   ChosenInds = 1:7
#' )
#'
#' plotDeFinetti(gen, cols = palette.colors(nrow(gen), "Accent"), pch = 19)
plotDeFinetti <- function(genotypes, cols, HWE = TRUE, tipLabels = c("Homozygous 0", "Heterozygous 1", "Homozygous 2"), ...) {
  ######################################
  ## Ternary to cartesian coordinates ##
  ######################################

  ternaryToCartesian <- function(A, B, C) {
    x <- 0.5 * (2 * C + B) / (A + B + C)
    y <- (sqrt(3) / 2) * B / (A + B + C)
    return(data.frame(x = x, y = y))
  }

  userArgs <- list(...)


	if(ncol(genotypes) == 4){
	  nMarkers <- rowSums(genotypes[, 2:4])
	  freq0 <- genotypes[, 2] / nMarkers
	  freq1 <- genotypes[, 3] / nMarkers
	  freq2 <- genotypes[, 4] / nMarkers
	} else {
      nMarkers <- rowSums(matrix(genotypes %in% c("0", "1", "2"), nrow = nrow(genotypes)))
      freq0 <- rowSums(genotypes == "0") / nMarkers
      freq1 <- rowSums(genotypes == "1") / nMarkers
      freq2 <- rowSums(genotypes == "2") / nMarkers
    }

  if (HWE) {
    p <- seq(0, 1, length.out = 100)
    pq <- 2 * p * (1 - p)
    q2 <- (1 - p)^2
    hwe <- ternaryToCartesian(p^2, pq, q2)
  }

  dat <- ternaryToCartesian(freq0, freq1, freq2)

  # plotting arguments
  plottingArgs <- utils::modifyList(list(
    xlim = c(0, 1),
    ylim = c(0, sqrt(3) / 2),
    xlab = "",
    ylab = "",
    pch = 1,
    asp = 1,
    axes = FALSE
  ), userArgs)

  # plotting

  do.call(plot, c(list(x = dat[, 1], y = dat[, 2], type = "n"), plottingArgs))
  if (HWE) {
    lines(hwe, lty = 3, lwd = 0.8, col = "grey")
  }
  plottingArgs$axes <- NULL
  do.call(points, c(list(x = dat[, 1], y = dat[, 2], col = cols), plottingArgs))
  lines(c(0, 1, .5, 0), c(0, 0, sqrt(3) / 2, 0))
  text(0, 0, tipLabels[1], pos = 1, xpd = NA)
  text(1, 0, tipLabels[3], pos = 1, xpd = NA)
  text(0.5, sqrt(3) / 2, tipLabels[2], pos = 3, xpd = NA)
}
