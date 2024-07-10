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
