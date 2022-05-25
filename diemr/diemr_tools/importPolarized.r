#' Imports genomic data polarized according to the specification
#'
#' Reads genotypes from a file and changes marker polarity.
#'
#' @inheritParams diem
#' @param file character vector with a single path to a file with genotypes.
#' @param changePolarity logical vector with length equal to the number of markers.
#' @details For details on the input data format, check the \code{file} with
#'   \code{CheckDiemFormat}.
#'
#'   The `changePolarity` argument influences how each marker is imported. Value
#'   `FALSE` means that the marker will be imported as it is saved in the `file`. Value
#'   `TRUE` means that the genotypes encoded as `0` will be imported as `2`, and genotypes
#'   encoded in the `file` as `2` will be imported as `0`.
#' @return Returns a character matrix with rows containing individual genotypes and columns
#'   containing markers.
#' @seealso \code{\link{diem}} for determining appropriate marker polarity with
#'   respect to a barrier to geneflow.
#' @export
#' @examples
#' dat <- importPolarized(
#'   file = system.file("extdata", "data6x3.txt", package = "diemr"),
#'   changePolarity = c(FALSE, TRUE, TRUE),
#'   ChosenInds = 1:6
#' )
#' dat
#' #    m1  m2  m3
#' # 1 "0" "1" "2"
#' # 2 "0" "0" "0"
#' # 3 "1" "1" "0"
#' # 4 "1" "2" "0"
#' # 5 "2" "2" "1"
#' # 6 "2" "2" "_"

importPolarized <- function(file, changePolarity, ChosenInds) {

  genotypes <- read.table(file = file, as.is = TRUE)
  if ((nchar(genotypes[1, ]) - 1) < max(ChosenInds)) {
    stop("File ", file, " contains fewer individuals than the maximum index specified in ChosenInds.")
  }
  if (nrow(genotypes) != length(changePolarity)) {
    stop("File ", file, " has ", nrow(genotypes), " markers, but changePolarity has length ", length(changePolarity))
  }
  genotypes <- as.data.frame(strsplit(unlist(genotypes), split = ""))[-1, ][ChosenInds, ]
  if (any(grepl("U", genotypes))) {
    genotypes[genotypes == "U"] <- "_"
  }
  # polarise markers
  genotypes <- apply(cbind(as.matrix(changePolarity), t(genotypes)),
    MARGIN = 1,
    FUN = function(x) emPolarise(origM = x[2:length(x)], changePolarity = as.logical(x[1]))
  )
  rownames(genotypes) <- ChosenInds
  colnames(genotypes) <- paste0("m", 1:ncol(genotypes))
  
  return(genotypes)
  
}
