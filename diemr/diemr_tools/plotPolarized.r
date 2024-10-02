#' Plot Polarized Genotypes
#'
#' Plots genotypes that can be optionally polarized.
#' @param genotypes A character matrix comprising of _012 encodings.
#' @param HI A numeric vector of individual hybrid indices with length equal to number
#'   of rows in \code{genotypes}.
#' @param cols A vector of four colors, representing missing data, homozygotes for
#'    genotype 0, heterozygotes and homozygotes for genotype 2.
#' @param ... Additional selected arguments passed to \link[graphics]{image} and
#'    \link[graphics]{axis}.
#' @details To import and polarize genotypes, use the \link{importPolarized} function.
#'
#'   When using \link{diem}, hybrid indices,
#'   \code{HI}, can be found in the file 'HIwithOptimalPolarities.txt'. Alternatively,
#'   calculate \code{HI} from the polarized genotypes as shown in the examples.
#'
#'  By default, the function plots colored tick marks for individuals, changing the 
#'  color at the steepest change in sorted \code{HI}. The second and fourth colors in 
#'  \code{cols} are used for the tick marks. 
#'  * To turn off this feature, use the argument \code{tick = FALSE}.
#'  * To use custom tick mark colors, provide a vector of colors for all individuals
#'  (equal to the number of rows in \code{genotypes}). The **vector of colors must be 
#'  ordered** according to \code{order(HI)}.
#'  * To include individual \code{labels} (e.g., accession numbers), provide a character
#'  vector with the **names in the same order as they are** in the \code{genotypes}. 
#'
#' @return No return value, called for side effects. In the default plot, purple and green
#'   represent sides of the barrier to gene flow encoded as `0` and `2`, respectively,
#'   yellow shows heterozygotes and white represents missing or undetermined genotypes. 
#'   Individuals are ordered according to the sorted \code{HI}.
#' @seealso \link{plotMarkerAxis} to add chromosome information to the x axis.
#' @importFrom graphics image axis
#' @importFrom grDevices col2rgb
#' @importFrom utils modifyList
#' @export
#' @examples
#' gen <- importPolarized(
#'   file = system.file("extdata", "data7x10.txt", package = "diemr"),
#'   changePolarity = c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE),
#'   ChosenInds = 1:7
#' )
#'
#' h <- apply(gen, 1, FUN = function(x) pHetErrOnStateCount(sStateCount(x)))[1, ]
#'
#' plotPolarized(genotypes = gen, HI = h)
#'
#' # Incorrect tick color order
#' plotPolarized(gen, h, col.ticks = c(rep("purple", 5), "green", "purple"), lwd = 3)
#'
#' # Correct tick color order
#' plotPolarized(gen, h, col.ticks = c(rep("purple", 5), "green", "purple")[order(h)], lwd = 3)
#'
#' # Correct individual label order 
#' plotPolarized(gen, h, labels = c(paste("purple", 1:5), "green 1", "purple 6"), ylab = "")
plotPolarized <- function(genotypes, HI, cols = c("#FFFFFF", "#800080", "#FFE500", "#008080"), ...) {
  userArgs <- list(...)

  # check arguments
  nMarkers <- ncol(genotypes)
  nInds <- nrow(genotypes)

  # check colors
  areColors <- sapply(cols, function(X) tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE))
  if (sum(areColors) != 4) {
    warning("Argument cols must contain four valid colors. Review the input, now using default colors")
    cols <- c("#FFFFFF", "#800080", "#FFE500", "#008080")
  }

  # check HI
  if (length(HI) != nInds) stop("Provide HI for all genotypes. Now ", length(HI), " HI for ", nInds, " genotypes.")

  # plotting arguments
  plottingArgs <- utils::modifyList(list(
    # image defaults
    xlab = "Markers",
    ylab = "Individuals",
    axes = FALSE,
    useRaster = TRUE,
    breaks = 0:4,
    # axis defaults
    side = 2,
    at = 1:nInds,
    col.ticks = c(rep(cols[2], which.max(diff(sort(HI)))), rep(cols[4], length(HI) - which.max(diff(sort(HI))))),
    labels = "",
    las = 1
  ), userArgs)
  acceptedImageArgs <- c(
    "zlim", "xlim", "ylim", "add", "xaxs", "yaxs", "xlab", "ylab", "breaks",
    "useRaster", "asp", "cex", "cex.lab", "cex.main", "cex.sub", "axes", "col.axis",
    "cex.axis", "family", "font", "font.axis", "font.lab", "font.main", "font.sub", "lab",
    "xpd"
  )
  acceptedAxisArgs <- c(
    "side", "at", "col.ticks", "labels", "las", "tick", "line",
    "pos", "outer", "font", "lty", "lwd", "lwd.ticks", "hadj", "padj", "gap.axis",
    "xpd"
  )
  imageArgs <- plottingArgs[names(plottingArgs) %in% acceptedImageArgs]
  axisArgs <- plottingArgs[names(plottingArgs) %in% acceptedAxisArgs]
  if(length(axisArgs$labels) > 1) axisArgs$labels <- axisArgs$labels[order(HI)]


  do.call(image, c(list(x = 1:nMarkers, y = 1:nInds, z = t(matrix(
    as.numeric(factor(genotypes[order(HI), ],
      levels = c("_", "0", "1", "2")
    )),
    ncol = nMarkers
  )), col = cols), imageArgs))


  if (length(axisArgs$col.ticks) == 1 && length(axisArgs$labels) == 1) {
    do.call(axis, modifyList(axisArgs, list(labels = rep(axisArgs$labels, length(axisArgs$at)))))
  } else if (length(axisArgs$col.ticks) == 1 && length(axisArgs$labels) > 1) {
    for (i in seq_along(axisArgs$labels)) {
      do.call(axis, modifyList(axisArgs, list(at = axisArgs$at[i], labels = axisArgs$labels[i], col.ticks = axisArgs$col.ticks)))
    }
  } else if (length(axisArgs$col.ticks) > 1 && length(axisArgs$labels) == 1) {
    for (i in seq_along(axisArgs$col.ticks)) {
      do.call(axis, modifyList(axisArgs, list(col.ticks = axisArgs$col.ticks[i], at = axisArgs$at[i], labels = axisArgs$labels)))
    }
  } else if (length(axisArgs$col.ticks) > 1 && length(axisArgs$labels) > 1) {
    for (i in seq_along(axisArgs$col.ticks)) {
      do.call(axis, modifyList(axisArgs, list(col.ticks = axisArgs$col.ticks[i], at = axisArgs$at[i], labels = axisArgs$labels[i])))
    }
  }
}
