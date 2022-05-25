#' Model of Diagnostic Marker Based on All Individual State Counts
#'
#' Estimates a diagnostic marker for the state counts of all genomic markers for all
#' individuals. Using the hypothetical, diagnostic marker, calculates individual state
#' counts with respect to their weighted similarity to the diagnostic marker states.
#'
#' @param I4 a matrix or data.frame with 4 numeric columns
#'           representing character state
#' 			 counts for missing data, homozygots for allele 1, heterozygots, and
#'           homozygots for allele 2. Individuals in rows.
#' @param OriginalHI numeric vector of length equal to number of rows in \code{I4},
#'           representing hybrid indices of individuals.
#' @inheritParams diem
#' @param ... parameters to be passed to other functions.
#' @details The \code{OriginalHI} can be calculated with \code{\link{pHetErrOnStateCount}}.
#' @return Matrix with dimensions of I4.
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics points segments par axis box
#' @importFrom utils capture.output
#' @importFrom zoo rollmean
#' @seealso \code{\link{diem}} for utilising the model to determine appropriate marker
#'   polarisation in estimating barriers to geneflow.
#' @examples
#' # state count matrix
#' dat <- matrix(c(0, 0, 1, 3, 1, 2, 2, 0, 2, 1, 4, 1), ncol = 4)
#'
#' # hybrid index calculation, assuming diploid markers
#' HI <- apply(dat * 2, MARGIN = 1, FUN = pHetErrOnStateCount)[1, ]
#'
#' # run model of diagnostics, with the weight of the ideal diagnostic marker being 0.8
#' ModelOfDiagnostic(I4 = dat, OriginalHI = HI, epsilon = 0.8)
#' #      [,1] [,2] [,3] [,4]
#' # [1,]  0.0  5.4  0.4  0.2
#' # [2,]  0.0  0.2  0.0  5.6
#' # [3,]  0.4  4.4  0.8  0.4
#' @export


ModelOfDiagnostic <- function(I4, OriginalHI, epsilon = 0.99, verbose = FALSE, ...) {

  # if the whole matrix of individuals statistics is provided, select only the hybrid index
  # assumes result from pHetErrOnStateCount()
  if (inherits(OriginalHI, "matrix")) {
    OriginalHI <- OriginalHI[1, ]
  }

  # select individuals with available hybrid index
  rows <- !is.na(OriginalHI)

  # if too few individuals have hybrid index, stop
  if (sum(rows) < 2) {
    return(NA)
  }

  if (sum(rows) >= 2) {
    # rescale hybrid index to [0,1]
    RescaledHI <- OriginalHI[rows]
    RescaledHI <- (RescaledHI - min(RescaledHI)) / diff(range(RescaledHI))

    # sorted rescaled hybrid index
    SRHI <- sort(RescaledHI)

    # derivative of the sorted rescaled hybrid index
    WarpSwitchAll <- data.frame(
      rollmean = zoo::rollmean(SRHI, k = 2),
      der = diff(SRHI, lag = 1)
    )

    # Value between the largest change in HI
    WarpSwitch <- WarpSwitchAll$rollmean[which.max(WarpSwitchAll$der)]

    # prepare empty vectors for variables
    RescaledHIsWithNAs <- WarpString <- betaOfRescaledHI <- rep(NA, length(OriginalHI))

    # theoretical diagnostic marker
    RescaledHIsWithNAs[rows] <- RescaledHI
    WarpString[!rows] <- "_"
    # indices of individuals with high/low hybrid index
    high <- which(sapply(RescaledHIsWithNAs > WarpSwitch, isTRUE))
    low <- which(sapply(RescaledHIsWithNAs <= WarpSwitch, isTRUE))
    WarpString[high] <- "2"
    WarpString[low] <- "0"

    # theoretical diagnostic marker in a 4state matrix, weighted by total number of markers
    Nmarkers <- sum(I4[1, ], na.rm = TRUE)
    Warp4 <- Nmarkers * t(apply(matrix(WarpString), MARGIN = 1, FUN = sStateCount))

    # scale hybrid index relative to maximum rate change of the hybrid index
    betaOfRescaledHI[!rows] <- 0.5
    betaOfRescaledHI[high] <- (RescaledHIsWithNAs[high] - WarpSwitch) / (1 - WarpSwitch)
    betaOfRescaledHI[low] <- (WarpSwitch - RescaledHIsWithNAs[low]) / WarpSwitch

    # write diagnostics
    if (inherits(verbose, "character")) {
      folder <- verbose
      verbose <- TRUE
    } else {
      folder <- NA
    }
    if (verbose & is.na(folder)) {
      if (!dir.exists("diagnostics")) dir.create("diagnostics")
      folder <- "diagnostics"
    }
    if (verbose) {
      pdf(paste0(folder, "/SortedRescaledHybridIndex.pdf"), height = 5)
      plot(SRHI,
        xlab = "Inds sorted by rescaled hybrid index", ylab = "Rescaled hybrid index",
        axes = ax <- ifelse(length(SRHI) > 10, TRUE, FALSE), ...
      )
      if (!ax) {
        axis(2, las = 1)
        axis(1, at = 1:length(SRHI), labels = c(1:length(SRHI))[order(RescaledHI)])
        box()
      }
      dev.off()

      pdf(paste0(folder, "/WarpSwitch.pdf"), height = 5)
      plot(WarpSwitchAll, type = "n", xlab = "Hybrid index", ylab = "First derivative of rescaled hybrid index", ...)
      segments(WarpSwitchAll$rollmean, rep(par("usr")[3], nrow(WarpSwitchAll)), WarpSwitchAll$rollmean, WarpSwitchAll$der, ...)
      points(WarpSwitchAll, ...)
      dev.off()

      pdf(paste0(folder, "/RescaledHItobetaOfRescaledHI.pdf"), height = 5)
      plot(RescaledHIsWithNAs, betaOfRescaledHI, type = "n", xlab = "Rescaled hybrid index", ylab = "Weighted rescaled hybrid index", ...)
      segments(RescaledHIsWithNAs, par("usr")[3], RescaledHIsWithNAs, betaOfRescaledHI, ...)
      points(RescaledHIsWithNAs, betaOfRescaledHI, ...)
      dev.off()

      write.table(round(OriginalHI, 4), file = paste0(folder, "/OriginalHI.txt"), sep = "\t", row.names = TRUE)
      cat("######  ModelOfDiagnostics WarpSwitch  #######\n\n", file = paste0(folder, "/WarpSwitch.txt"), append = FALSE)
      cat("WarpSwitch with maximum change in hybrid index: ", round(WarpSwitch, 4), "\n\n", file = paste0(folder, "/WarpSwitch.txt"), append = TRUE)
      cat("WarpString for a diagnostic marker: \n", paste(WarpString, collapse = ""), "\n\n", file = paste0(folder, "/WarpSwitch.txt"), append = TRUE)
      cat("Head of WarpString4:\n", file = paste0(folder, "/WarpSwitch.txt"), append = TRUE)
      capture.output(head(Warp4, 10), file = paste0(folder, "/WarpSwitch.txt"), append = TRUE)
    }


    # update counts relative to the theoretical diagnostic marker
    res <- ((epsilon * betaOfRescaledHI) * Warp4) + ((1 - (epsilon * betaOfRescaledHI)) * I4)
    rownames(res) <- rownames(I4)
    colnames(res) <- colnames(I4)

    return(res)
  }
}
