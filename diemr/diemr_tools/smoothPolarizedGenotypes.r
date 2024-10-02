#' Smooths Polarized Genotype States
#'
#' This function smooths polarized genotype states using a Laplace kernel density estimation.
#' It calculates a smoothed version of the genotype states over specified windows.
#'
#' @inheritParams rank2map
#' @inheritParams plotPolarized
#' @param windows A two-column numeric matrix with indices of start and end positions for 
#'  windows for all markers indicated by \code{ChosenSites}.
#' @param ... Additional parameters to be passed to \link{rank2map} if \code{windows = NULL}.
#' @details Ensure that \code{ChosenSites} is the same as was used to import polarized genotypes.
#'
#' The function uses a Laplace kernel to weight the genotype states within a window around each marker position,
#' based on physical positions of the markers. The smoothing process accounts for chromosome-level scales.
#' 
#' The Laplace kernel density is calculated as:
#' \deqn{\frac{1}{2b} \exp\left(\frac{-|x - \mu|}{b}\right)}
#' where \eqn{x} is the position, \eqn{\mu} is the center of the kernel, and \eqn{b} is the scale parameter.
#' The scale parameter is \eqn{b = 10^{-4} n}, where \eqn{n} is the position of the last 
#' marker on the chromosome.
#' @importFrom stats weighted.mean
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
smoothPolarizedGenotypes <- function(genotypes, includedSites, ChosenSites = "all", windows = NULL, windowSize = NULL, ...){
  # genotypes are numeric
  genotypes <- suppressWarnings(apply(genotypes, 2, as.numeric))
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
  laplaceDensity <- function(x, mu, b) {
    (1 / (2 * b)) * exp(-abs(x - mu) / b)
  }
  
  # chromosome level scales b
  # https://www.biorxiv.org/content/10.1101/2024.06.03.597101v1.full
  b <- 1e-4 * bed$POS[!duplicated(bed$CHROM, fromLast = TRUE)]
  b <- rep(b, rle(bed$CHROM)$lengths)
  
  gen2 <- matrix(NA, ncol = nMarkers, nrow = nInds)
  
  for(i in 1:nInds){
    for(j in 1:nMarkers){
      whichSites <- windows[j, 1]:windows[j, 2]
      genotypesInWindow <- !is.na(genotypes[i, whichSites])
      if(any(genotypesInWindow)){
      gen2[i,j] <- weighted.mean(genotypes[i, whichSites], 
        weights = laplaceDensity(x = bed$POS[whichSites], mu = bed$POS[j], b = b[j]),
        na.rm = TRUE)
    }
    
    }
  
  }
  return(round(gen2))
}