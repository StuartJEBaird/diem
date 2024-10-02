#' Imports genomic data polarized according to the specification
#'
#' Reads genotypes from a file and changes marker polarity.
#'
#' @inheritParams diem
#' @param changePolarity A logical vector with length equal to the number of markers.
#' @param verbose Logical whether to show messages on import progress.
#' @param ... Optional numeric vector of \code{compartmentSizes}.
#' @details For details on the input data format, check the \code{file} with
#'   \link{CheckDiemFormat}.
#'
#'   The `changePolarity` argument influences how each marker is imported. Value
#'   `FALSE` means that the marker will be imported as it is saved in the `file`. Value
#'   `TRUE` means that the genotypes encoded as `0` will be imported as `2`, and genotypes
#'   encoded in the `file` as `2` will be imported as `0`.
#' @return Returns a character matrix with rows containing individual genotypes and columns
#'   containing markers.
#' @seealso \link{diem} for determining appropriate marker polarity with
#'   respect to a barrier to gene flow.
#' @export
#' @examples
#' dat <- importPolarized(
#'   files = system.file("extdata", "data7x3.txt", package = "diemr"),
#'   changePolarity = c(FALSE, TRUE, TRUE),
#'   ChosenInds = 1:6,
#'   ChosenSites = "all"
#' )
#' dat
#' #    m1  m2  m3
#' # 1 "0" "1" "2"
#' # 2 "0" "0" "0"
#' # 3 "1" "1" "0"
#' # 4 "1" "2" "0"
#' # 5 "2" "2" "1"
#' # 6 "2" "2" "_"
importPolarized <- function(files, changePolarity, ChosenInds, ChosenSites = "all", nCores = 1, verbose = FALSE, ...) {
  ChosenSites <- resolveCompartments(files = files, toBeCompartmentalized = ChosenSites, ...)
  if(verbose) message("ChosenSites for compartments done ", Sys.time())

  markerLabels <- which(unlist(ChosenSites))
  changePolarity <- resolveCompartments(files = files, toBeCompartmentalized = changePolarity, ...)
  if(verbose) message("changePolarity for compartments done ", Sys.time())

  allCompartments <- parallel::mclapply(1:length(files), mc.cores = nCores,
    FUN = function(i){
    if(sum(ChosenSites[[i]]) == 0){ 
    return(NA)
    }
  
    genotypes <- read.table(file = files[i], as.is = TRUE)
    if ((nchar(genotypes[1, ]) - 1) < max(ChosenInds)) {
      stop("File ", files[[i]], " contains fewer individuals than the maximum index specified in ChosenInds.")
    }
    # select genotypes for chosen individuals and markers
    genotypes <- as.data.frame(strsplit(unlist(genotypes), split = ""))[-1, ][ChosenInds, ChosenSites[[i]], drop = FALSE]
    if (any(grepl("U", genotypes))) {
      genotypes[genotypes == "U"] <- "_"
    }
    # polarise chosen markers
    genotypes <- apply(cbind(as.matrix(changePolarity[[i]][ChosenSites[[i]]]), t(genotypes)),
      MARGIN = 1,
      FUN = function(x) emPolarise(origM = x[2:length(x)], changePolarity = as.logical(x[1]))
    )
    return(genotypes)
    })
  
  if(verbose) message("Importing all compartments done ", Sys.time())

  allCompartments <- Filter(Negate(anyNA), allCompartments)
  allCompartments <- do.call(cbind, allCompartments)
  rownames(allCompartments) <- ChosenInds
  colnames(allCompartments) <- paste0("m", markerLabels)
  
  return(allCompartments)
}
