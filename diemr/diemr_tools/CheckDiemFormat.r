#' diem input file checker
#'
#' Checks format of files with genotype data.
#'
#' @inheritParams diem
#' @export
#' @details The input file must have genotypes of one marker for all individuals on one
#'  line. The line must start with a letter "S" and contain only characters
#'  "_" or "U" for unknown genotypes or a third/fourth allele, "0" for homozygots for
#'  allele 1, "1" for heterozygots, and "2" for homozygots for allele 2. Check the 
#'  vignette with \code{browseVignettes(package = "diemr")} for the example of the 
#'  input format.
#'
#'  Ploidies must be given as a list with each element corresponding to a genomic
#'  compartment (aka a file). For each compartment, the numeric vector specifying
#'  ploidies of all individuals chosen for the specific analysis must be given.
#' @return Returns invisible \code{TRUE} if all files are executable by \code{diem}. Exits
#'    with informative error messages otherwise, specifying file names and lines with
#'    potential problems. When too many lines contain problems, the first six are given.
#' @examples
#' # set up input genotypes file names, ploidies and selection of individual samples
#' inputFile <- system.file("extdata", "data6x3.txt", package = "diemr")
#' ploidies <- list(c(2, 1, 2, 2, 2, 1))
#' inds <- 1:6
#'
#' # check input data
#' CheckDiemFormat(files = inputFile, ploidy = ploidies, ChosenInds = inds)
#' #  File check passed: TRUE
#' #  Ploidy check passed: TRUE
CheckDiemFormat <- function(files, ChosenInds, ploidy) {


  ##########################################
  # Checks format in one compartment file
  ##########################################

  compartmentCheck <- function(file, ChosenInds) {
    res <- FALSE
    # filename a character vector
    if (!inherits(file, "character")) {
      stop("The file argument needs to be a character string specifying the path to the input file. Instead, file is ", class(file)[1])
    } else {
      # file exists
      if (!file.exists(file)) {
        stop("File ", file, " cannot be found. A full path to the file might be necessary, or change working directory correspondingly.")
      } else {
        # markers start with a character
        dat <- readLines(file)
        sFormat <- grepl(pattern = "^S", x = dat, ignore.case = FALSE)
        if (any(!sFormat)) {
          stop("Lines ", paste(head(which(!sFormat)), collapse = ", "), " in file ", file, " do not start with a letter 'S'. Prefix 'S' before the genotype string. Check also for invisible characters.")
        } else {
          # number of individuals equal
          nIndividuals <- nchar(dat) - 1
          if (length(unique(nIndividuals)) != 1) {
            return(FALSE)
            Mode <- function(x) {
              ux <- unique(x)
              ux[which.max(tabulate(match(x, ux)))]
            }
            wrongNind <- Mode(nIndividuals) - nIndividuals != 0
            stop("Markers on lines ", paste(head(which(wrongNind)), collapse = ", "), " in file ", file, " were genotyped for a different number of individuals. Make sure the line lengths are the same.")
          } else {
            # maximum index of ChosenInds
            if (max(nIndividuals) < max(ChosenInds)) {
              stop("File ", file, " contains fewer individuals than the maximum index specified in ChosenInds.")
            } else {
              # _012 symbols
              # dat <- sub("^[a-z]", "s", dat)
              fourStateQdata <- grepl("[^S_U012]", dat)
              if (any(fourStateQdata)) {
                stop("File ", file, " contains characters other than _012 on line(s) ", paste(head(which(fourStateQdata)), collapse = ", "))
              } else {
                res <- TRUE
              } # _012 symbols
            } # maximum index of ChosenInds
          } # number of individuals equal
        } # markers start with a character
      } # input file exists
    } # file name as character vector
    return(res)
  }


  for (x in files) {
    res <- compartmentCheck(file = x, ChosenInds = ChosenInds)
  }

  message("File check passed: ", all(unlist(res)))

  res <- FALSE

  # Check ploidies to be a list of length(files) vectors with length of length(ChosenInds)
  if (!inherits(ploidy, "list")) {
    stop("Ploidy must be a list of length ", length(files), " with elements being numeric vectors of length ", length(ChosenInds))
  } else {
    if (length(ploidy) != length(files)) {
      stop("Length of ploidy (", length(ploidy), ") is not equal to the length of files (", length(files), ").")
    } else {
      pLength <- unlist(lapply(ploidy, FUN = function(x) length(x) == length(ChosenInds)))
      if (any(!pLength)) {
        stop("Ploidy for compartment(s) ", which(!pLength), " is not a numeric vector of length ", length(ChosenInds))
      } else {
        if (!all(unlist(ploidy) %in% c(0, 1, 2))) {
          nPloidy <- matrix(unlist(ploidy) %in% c(0, 1, 2), ncol = length(files))
          stop("Ploidy must be 0, 1, or 2. Comparment(s) ", head(which(apply(nPloidy, 2, any))), " contain other characters.")
        } else {
          res <- TRUE
        } # Ploidy 0,1,2
      } # Ploidy in compartments at length of individuals
    } # Ploidy for all compartments
  } # Ploidy is a list

  message("Ploidy check passed: ", res)


  invisible(res)
}
