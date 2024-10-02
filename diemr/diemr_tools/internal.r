# Creates a list of TRUE/FALSE values with elements equal to the number of sites in each
# compartment file

resolveCompartments <- function(files, toBeCompartmentalized, compartmentSizes = NULL) {
  if(missing(files)) stop("Provide paths to diem input files in the `files` argument.")
  if (is.null(compartmentSizes)) {
    compartmentSizes <- unname(sapply(files, FUN = \(x) length(readLines(x))))
  }

  nMarkers <- sum(compartmentSizes)

  if (inherits(toBeCompartmentalized, "character")) {
    toBeCompartmentalized <- lapply(compartmentSizes, FUN = function(x) rep(TRUE, x))
  } else {
    if (inherits(toBeCompartmentalized, "logical")) {
      if (length(toBeCompartmentalized) != nMarkers) {
        compartmentSizes <- unname(sapply(files, FUN = \(x) length(readLines(x))))
        if (length(toBeCompartmentalized) != nMarkers) {
          stop("toBeCompartmentalized does not have the same length (", nMarkers, ") as the number of sites is the files.")
        }
      }

      compartmentLabels <- rep(seq_along(compartmentSizes), compartmentSizes)
      toBeCompartmentalized <- unname(split(toBeCompartmentalized, compartmentLabels))
    } else {
      stop("toBeCompartmentalized must be either `all` or a logical vector with the length equal to the total number of markers in `files`")
    }
  }

  return(toBeCompartmentalized)
}


# Calculates positions of chromosome breaks and physical spacing along the chromosome,
# plots axis with optional graphical parameters

markerAxis <- function(includedSites, ChosenSites, tickDist) {
  bed <- readIncludedSites(includedSites = includedSites, ChosenSites = ChosenSites)

  # chromosome labels
  CHROMbreaks <- c(which(!duplicated(bed$CHROM)) - 0.5, nrow(bed) + 0.5)
  CHROMnamesPos <- zoo::rollmean(CHROMbreaks, k = 2)
  CHROMnames <- unique(bed$CHROM)

  # tick positions

  ticksPos <- integer(0)
  ticksNames <- integer(0)

  for (i in 1:length(CHROMnames)) { # i - chromosome
    pos <- bed$POS[bed$CHROM == CHROMnames[i]] # positions in a chromosome
    if (max(pos) > tickDist) {
      ticksFor <- seq(tickDist, max(pos), by = tickDist)
      for (k in ticksFor) {
        ticksAt <- which(pos > k)[1]
        if (length(ticksAt) > 0) {
          ticksPos <- c(ticksPos, which(bed$CHROM == CHROMnames[i])[ticksAt])
          ticksNames <- c(ticksNames, k)
        }
      }
    }
  }
  ticksPos <- ticksPos - 0.5
  ticksNames <- ticksNames / 1000000

  return(list(
    CHROMbreaks = CHROMbreaks,
    CHROMnamesPos = CHROMnamesPos,
    CHROMnames = CHROMnames,
    ticksPos = ticksPos,
    ticksNames = ticksNames
  ))
}


readIncludedSites <- function(includedSites, ChosenSites = "all"){
  bed <- read.table(includedSites, header = TRUE, sep = "\t")
  if(inherits(ChosenSites, "logical") || inherits(ChosenSites, "numeric")){
    bed <- bed[ChosenSites, ]
  }
  return(bed)
}


# Calculates position of SNPs on a chromosome with respect to the SNP rank in the 
# selected markers
# Original up to the diemr 1.4.1
#' @aliases rank2map
#rank2mapChr <- function(x, windowSize = 3) {
#  n <- length(x)
#  if(is.null(windowSize)){
#    warning("windowSize is NULL. Using windowSize = 1e+07.")
#    windowSize <- 1e+07
#  }
#  halfSize <- windowSize / 2
#  res <- matrix(NA, ncol = 2, nrow = n, dimnames = list(NULL, c("start", "end")))
#  
#  for (i in seq_len(n)) {
#    backwardPos <- i  
#    forwardPos <- i  
#        
#    while (backwardPos > 1 && abs(x[i] - x[backwardPos - 1]) <= halfSize) {
#      backwardPos <- backwardPos - 1
#    }
#    res[i, 1] <- backwardPos
#
#    while (forwardPos < n && abs(x[i] - x[forwardPos + 1]) <= halfSize) {
#      forwardPos <- forwardPos + 1
#    }
#    res[i, 2] <- forwardPos
#  }
#  
#  return(res)
#}
# New from Inchworm by Stuart J.E. Baird
# From diemr 1.4.2
rank2mapChr <- function(x, windowSize) {
  lengthX <- length(x)
  reach <- windowSize / 2
  focus <- 1
  left <- 1
  right <- 1
  lesseq <- x[1] - 1
  ans <- matrix(0, nrow = lengthX, ncol = 2)
  
  for (xfocus in x) {
    if (lesseq > xfocus) stop("unsorted argument")
    
    while (x[right] < xfocus + reach && right < lengthX) {
      right <- right + 1
    }
    
    if (x[right] > xfocus + reach) {
      right <- right - 1
    }
    
    while (x[left] + reach < xfocus && left < focus) {
      left <- left + 1
    }
    
    ans[focus, ] <- c(left, right)
    lesseq <- xfocus
    focus <- focus + 1
  }
  
  return(ans)
}

