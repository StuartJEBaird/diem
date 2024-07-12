# Creates a list of TRUE/FALSE values with elements equal to the number of sites in each
# compartment file

resolveCompartments <- function(files, toBeCompartmentalized, compartmentSizes = NULL) {
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
  bed <- read.table(includedSites, header = TRUE, sep = "\t")[unlist(ChosenSites), ]

  # chromosome labels
  CHROMbreaks <- c(which(!duplicated(bed$CHROM)) - 0.5, nrow(bed) + 0.5)
  CHROMnamesPos <- zoo::rollmean(CHROMbreaks, k = 2)
  CHROMnames <- unique(bed$CHROM)

  # tick positions

  ticksPos <- integer(0)
  ticksPower <- integer(0)

  for (i in 1:length(CHROMnames)) { # i - chromosome
    pos <- bed$POS[bed$CHROM == CHROMnames[i]] # positions in a chromosome
    if (max(pos) > tickDist) {
      ticksFor <- seq(tickDist, max(pos), by = tickDist)
      for (k in ticksFor) {
        ticksAt <- which(pos > k)[1]
        if (length(ticksAt) > 0) {
          ticksPos <- c(ticksPos, which(bed$CHROM == CHROMnames[i])[ticksAt])
          ticksPower <- c(ticksPower, k)
        }
      }
    }
  }
  ticksPos <- ticksPos - 0.5
  ticksPower <- ticksPower / 1000000

  return(list(
    CHROMbreaks = CHROMbreaks,
    CHROMnamesPos = CHROMnamesPos,
    CHROMnames = CHROMnames,
    ticksPos = ticksPos,
    ticksPower = ticksPower
  ))
}
