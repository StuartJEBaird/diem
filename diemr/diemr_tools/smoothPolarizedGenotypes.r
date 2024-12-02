#' Smooth Polarized Genotype States
#'
#' This function smooths polarized genotype states using a Laplace kernel density estimation.
#' It calculates a smoothed version of the genotype states over specified physical extent
#' of genomic content around a site. At each
#' genomic position, the function returns a weighted mode of the genomic state.
#'
#'
#' @inheritParams rank2map
#' @inheritParams plotPolarized
#' @param windows A two-column numeric matrix with indices of start and end positions for 
#'  windows for all markers indicated by \code{ChosenSites}. If \code{windows = NULL}, the function 
#'  calculates the windows using \link{rank2map}.
#' @param ... Additional arguments to be passed to \link{rank2map} if \code{windows = NULL}.
#' @details Ensure that \code{ChosenSites} match those used to import polarized genotypes.
#'
#' The function uses a truncated and scaled Laplace kernel to weight the genotype states 
#' within a window around each marker position,
#' based on physical positions of the markers.
#' 
#' The Laplace kernel weights are calculated for physical positions of the sites centered at
#' the site being smoothed as:
#' \deqn{\frac{10}{19} \exp\left(\frac{-x}{b}\right),}
#' when \eqn{x < 0}, and as:
#' \deqn{\frac{10}{19} \exp\left(\frac{x}{b}\right),}
#' when \eqn{x \geq 0},
#' where \eqn{x} is the position, and \eqn{b} is the scale parameter of the Laplace kernel.
#' The scale parameter is equal to:
#' \deqn{b = \frac{\text{windowSize}}{2 \ln(20)}.}
#' 
#' @importFrom stats aggregate
#' @seealso \link{rank2map}
#' @export
#' @examples
#'  \dontrun{
#'  # Run this example in a working directory with write permissions
#'  myo <- system.file("extdata", "myotis.vcf", package = "diemr")
#'  vcf2diem(myo, "myo")
#'  fit <- diem("myo-001.txt", ChosenInds = 1:14)
#'  gen <- importPolarized("myo-001.txt", changePolarity = fit$markerPolarity, ChosenInds = 1:14)
#'  h <- apply(gen, 1, \(x) pHetErrOnStateCount(sStateCount(x)))[1, ]
#'  gen2 <- smoothPolarizedGenotypes(genotypes = gen, 
#'     includedSites = "myo-includedSites.txt", windowSize = 50)
#'  plotPolarized(gen, h)
#'  plotPolarized(gen2, h)
#'  }
smoothPolarizedGenotypes <- function(genotypes, includedSites, ChosenSites = "all", windows = NULL, windowSize = 250000, ...){
  
  nMarkers <- ncol(genotypes)
  nInds <- nrow(genotypes)
  
  # windows has rank ranges for kernels in rows
  if(is.null(windows)){
    windows <- rank2map(includedSites = includedSites,
      ChosenSites = ChosenSites,
      windowSize = windowSize, ...)
  }
  
  # physical marker positions
  bed <- readIncludedSites(includedSites = includedSites, ChosenSites = ChosenSites)

  # Laplace kernel
  laplaceScale <- windowSize / (2 * log(20))
  centeredPositions <- ceiling(-windowSize / 2):floor(windowSize / 2)
  allLaplaceWeights <- truncatedLaplace(x = centeredPositions,
  	laplaceScale = laplaceScale)
  

  
  gen2 <- matrix(NA, ncol = nMarkers, nrow = nInds)
  
  for(i in 1:nInds){
    for(j in 1:nMarkers){
      whichSites <- windows[j, 1]:windows[j, 2]
      genotypesInWindow <- !is.na(genotypes[i, whichSites])
      if(any(genotypesInWindow)){
      	centeredSitePositions <- bed$POS[whichSites] - bed$POS[j]
      gen2[i, j] <- unbiasedWeightedStateChoice(
      	genomicStates = genotypes[i, whichSites],
      	laplaceWeights = allLaplaceWeights[match(centeredSitePositions, centeredPositions)]
      )
    }
    
    }
  
  }
  return(gen2)
}