#' Add a Marker Axis with Chromosome Names to a Plot of Polarized Genotypes
#'
#' This function adds a marker axis with chromosome names to an existing plot of polarized
#' genotypes. It requires that the plot is already created using \code{plotPolarized}.
#'
#' @param includedSites character. Path to a file with columns \code{CHROM} and \code{POS}.
#' @inheritParams diem
#' @param tickDist numeric. Indicates the spacing of physical tick marks along a chromosome.
#' @details The \code{includedSites} file should ideally be generated by
#' \link[diemr]{\code{vcf2diem}} to ensure congruence between the plotted genotypes and
#' the respective metadata.
#'
#' Tick mark distances within a chromosome are located at \code{tickDist} and formated to
#' multiples of milions.
#' @export
#' @examples
#' \dontrun{
#' # Run this example in a working directory with write permissions
#' myo <- system.file("extdata", "myotis.vcf", package = "diemr")
#' vcf2diem(myo, "myo")
#' inds <- 1:14
#' fit <- diem("myo-001.txt", ChosenInds = inds, ploidy = FALSE)
#' gen <- importPolarized("myo-001.txt", fit$markerPolarity, inds)
#' h <- apply(gen, 1, function(x) pHetErrOnStateCount(sStateCount(x)))[1, ]
#' plotPolarized(gen, h, xlab = "")
#' plotMarkerAxis("myo-includedSites.txt", rep(TRUE, 11), tickDist = 100)
#' }
plotMarkerAxis <- function(includedSites, ChosenSites, tickDist = 1e+06, ...) {
  if (dev.cur() == 1) {
    stop("Plot polarized genotypes with plotPolarized first.")
  }
  axisInfo <- markerAxis(
    includedSites = includedSites,
    ChosenSites = ChosenSites, tickDist = tickDist
  )

  userArgs <- list(...)
  plottingArgs <- utils::modifyList(list(
    # axis defaults
    side = 1,
    las = 1,
    tcl = -.5,
    cex = 1,
    line = 0
  ), userArgs)


  acceptedAxisArgs <- c(
    "side", "col.ticks", "labels", "las", "tick", "line", "tcl", "cex", "cex.axis",
    "pos", "outer", "font", "lty", "lwd", "lwd.ticks", "hadj", "padj", "gap.axis",
    "xpd"
  )
  axisArgs <- plottingArgs[names(plottingArgs) %in% acceptedAxisArgs]

  # ticks of physical distances
  do.call(axis, modifyList(
    axisArgs,
    list(at = axisInfo$ticksPos, labels = axisInfo$ticksPower, cex = axisArgs$cex)
  ))
  # ticks for chromosome separators
  do.call(axis, modifyList(
    axisArgs,
    list(at = axisInfo$CHROMbreaks, tcl = axisArgs$tcl * 2.2, labels = FALSE)
  ))
  # chromosome names
  do.call(mtext, modifyList(
    axisArgs,
    list(
      at = axisInfo$CHROMnamesPos,
      text = axisInfo$CHROMnames,
      line = axisArgs$line + 2.5,
      cex = axisArgs$cex * 1.1
    )
  ))
}
