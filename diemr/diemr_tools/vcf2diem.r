#' Convert vcf files to diem format
#'
#' Reads vcf files and writes genotypes of the most frequent alleles based on
#' chromosome positions to diem format.
#'
#' @param SNP A character vector with a path to the '.vcf' or '.vcf.gz' file, or an \code{vcfR}
#'     object. Diploid data are currently supported.
#' @param filename A character vector with a path where to save the converted genotypes.
#' @param chunk Numeric indicating by how many markers should the result be split into
#'     separate files.
#' @param requireHomozygous A logical or numeric vector indicating whether to require the marker
#' to have at least one or more
#'     homozygous individual(s) for each allele.
#' @inheritParams diem
#'
#' @details Importing vcf files larger than 1GB, and those containing multiallelic
#'    genotypes is not recommended. Instead, use the path to the
#'    vcf file in \code{SNP}. \code{vcf2diem} then reads the file line by line, which is
#'    a preferred solution for data conversion, especially for
#'    very large and complex genomic datasets.
#'
#'    The number of files \code{vcf2diem} creates depends on the \code{chunk} argument
#'    and class of the \code{SNP} object.
#'
#'    * Values of \code{chunk < 100} are interpreted as the number of files into which to
#'    split data in \code{SNP}. For \code{SNP} object of class \code{vcfR}, the number
#'    of markers per file is calculated from the dimensions of \code{SNP}. When class
#'    of \code{SNP} is \code{character}, the number of markers per file is approximated
#'    from a model with a message. If this number of markers per file is inappropriate
#'    for the expected
#'    output, provide the intended number of markers per file in \code{chunk} greater
#'    than 100 (values greater than 10000 are recommended for genomic data).
#'    \code{vcf2diem} will scan the whole input specified in the \code{SNP} file, creating
#'    additional output files until the last line in \code{SNP} is reached.
#'    * Values of \code{chunk >= 100} mean that each output file
#'    in diem format will contain \code{chunk} number of lines with the data in \code{SNP}.
#'
#'    When the vcf file contains markers not informative for genome polarisation,
#'    those are removed and listed in a file ending with *omittedSites.txt* in the
#'    directory specified in the \code{SNP} argument or in the working directory.
#'    The omitted loci are identified by their information in the CHROM and POS columns,
#'    and include the QUAL column data. The last column is an integer specifying
#'    the reason why the respective marker was omitted. The reasons why markers are
#'    not informative for genome polarisation using \code{diem} are:
#'    1. Marker has fewer than 2 alleles representing substitutions.
#'    2. Required homozygous individuals for the 2 most frequent alleles are not present
#'    (optional, controlled
#'    by the \code{requireHomozygous} argument).
#'    3. The second most frequent allele is found only in one heterozygous individual.
#'    4. Dataset is invariant for the most frequent allele.
#'    5. Dataset is invariant for the allele listed as the first ALT in the vcf input.
#'
#'    The CHROM, POS, and QUAL information for loci included in the converted files are
#'    listed in the file ending with *includedSites.txt*. Additional columns show which
#'    allele is
#'    encoded as 0 in its homozygous state and which is encoded as 2.
#' @return No value returned, called for side effects.
#' @importFrom vcfR getFIX extract.gt
#' @importFrom tools file_ext file_path_sans_ext
#' @export
#' @author Natalia Martinkova
#' @author Filip Jagos <521160@mail.muni.cz>
#' @author Jachym Postulka <506194@mail.muni.cz>
#' @examples
#' \dontrun{
#' # vcf2diem will write files to a working directory or a specified folder
#' # make sure the working directory or the folder are at a location with write permission
#' myofile <- system.file("extdata", "myotis.vcf", package = "diemr")
#'
#' vcf2diem(SNP = myofile, filename = "test1")
#' vcf2diem(SNP = myofile, filename = "test2", chunk = 3)
#' }
vcf2diem <- function(SNP, filename, chunk = 1L, requireHomozygous = TRUE, ChosenInds = "all") {
  if (!inherits(SNP, c("character", "vcfR"), which = FALSE)) {
    stop("'SNP' must be either a 'vcfR' object or a 'character' string with path to a vcf file.")
  }
  if (length(SNP) != 1) {
    stop("Provide a single object or path in 'SNP'.")
  }
  if (missing(filename) || !is.character(filename)) {
    stop("Provide a filename in the argument 'filename'.")
  }
  if (inherits(SNP, "character", which = FALSE) && !(tools::file_ext(SNP) %in% c("vcf", "gz"))) {
    stop("Input filename in 'SNP' should have an extension 'vcf' or 'vcf.gz'.")
  }
  if (inherits(SNP, "character", which = FALSE) && !file.exists(SNP)) {
    stop("File ", SNP, " not found.")
  }
  if (!is.numeric(chunk)) {
    stop("'chunk' must be an integer.")
  }
  if (length(chunk) > 1) {
    chunk <- chunk[1]
    warning("Different chunk sizes are not permitted. Using chunk size of ", chunk, " for all files.")
  }
  minHomozygous <- as.integer(requireHomozygous)
  requireHomozygous <- requireHomozygous > 0



  ######################################
  ######################################
  ####  Declare internal functions  ####
  ######################################
  ######################################


  ################################################
  ####  Resolve file structure and chunk size ####
  ################################################

  ResolveOutput <- function(SNP, filename, chunk) {
    origChunk <- chunk
    # use only the first filename
    if (length(filename) > 1) {
      filename <- filename[1]
    }
    # strip file extension and replace with a file number and txt
    fileext <- tools::file_ext(filename)
    filepath <- tools::file_path_sans_ext(filename)
    filename <- filepath
    # file names for loci positions
    lociFiles <- c(
      paste0(filepath, "-omittedSites.txt"),
      paste0(filepath, "-includedSites.txt"),
      paste0(filepath, "-sampleNames.txt")
    )
    # estimate chunk size
    if (chunk < 100) {
      if (any(class(SNP) == "vcfR")) {
        filename <- paste0(
          filename, "-",
          formatC(1:chunk, width = 3, flag = 0),
          ".",
          ifelse(fileext == "", "txt", fileext)
        )
        chunk <- unname(ceiling(nrow(SNP@fix) / chunk))
      } else {
        filesize <- file.size(SNP)[1]
        # model prediction on possible number of markers per chunk
        chunk <- ceiling(10^(-1.2387 + 0.7207 * log10(filesize)) / chunk)
        message("Expecting to include ", chunk, " markers per diem file.\nIf you expect more markers in the file, provide a suitable chunk size.")
      }
    }
    if (chunk >= 100) {
      message("Will include up to ", chunk, " markers per diem file.")
    }
    return(list(filename, chunk, origChunk, lociFiles))
  }


  ###############################################
  ####  Resolve indels and multiple alleles  ####
  ###############################################

  ResolveGenotypes <- function() {
    # attempt to resolve sites with indels in REF and ALT
    INFO[, 4] <- sub("\\*|\\.", "AA", INFO[, 4])
    INFO[, 5] <- sub("\\*|\\.", "AA", INFO[, 5])
    ALLELES <- apply(INFO[, 4:5, drop = FALSE], 1, FUN = \(x) paste(x[1], x[2], sep = ","))
    resolvable <- lapply(strsplit(ALLELES, ",", fixed = TRUE),
      FUN = \(x) which(nchar(x) == 1) - 1 # allele numbers, not indices; removes allele numbers for indels
    )
    resolvable[lengths(resolvable) <= 1] <- NA # sets markers as not resolvable if less than 2 substitutions remain
    reason <- 1

    # remove markers with unresolvable indels
    SNP[indels <- sapply(resolvable, FUN = anyNA, simplify = TRUE), ] <- NA

    # resolve multiallelic markers
    multiallelic <- which(grepl(",", ALLELES) & !indels)
    if (length(multiallelic) > 0) {
      for (i in multiallelic) {
        alleleCounts <- table(unlist(strsplit(SNP, "/|\\|")))
        # check that the site is not invariant
        if (sum(as.character(resolvable[[i]]) %in% names(alleleCounts)) < 2) {
          next
        }
        # select allele numbers with the highest allele counts
        majorAlleles <- names(sort(alleleCounts[names(alleleCounts) %in% as.character(0:9)], decreasing = TRUE)[1:2])
        SNP[i, ] <- gsub(
          pattern = paste0("[", paste(c(0:9)[-(1 + as.numeric(majorAlleles))], collapse = ""), "]"),
          replacement = "\\.",
          x = SNP[i, ]
        )
        SNP[i, ] <- gsub(majorAlleles[1], "A", SNP[i, ])
        SNP[i, ] <- gsub(majorAlleles[2], "B", SNP[i, ])
        SNP[i, ] <- gsub("A", "0", SNP[i, ])
        SNP[i, ] <- gsub("B", "1", SNP[i, ])
      }
    }

    # genotypes from the most frequent substitutions
    SNP[is.na(SNP)] <- "_"
    patterns <- c(".*\\..*", "0.0", "0.1", "1.0", "1.1")
    replacements <- c("_", 0, 1, 1, 2)
    for (i in 1:5) {
      SNP <- sub(patterns[i], replacements[i], SNP, perl = TRUE)
    }

    # identify non-informative markers
    I4 <- t(apply(SNP, MARGIN = 1, FUN = sStateCount))
    nonInformative <- FALSE
    if (requireHomozygous && (I4[, 2] < minHomozygous || I4[, 4] < minHomozygous)) { # homozygous individuals missing
      nonInformative <- TRUE
      reason <- 2
    } else if ((I4[, 3] == 1 && (I4[, 2] == 0 || I4[, 4] == 0))) { # only one heterozygous individual
      nonInformative <- TRUE
      reason <- 3
    } else if (I4[, 2] == 0 && I4[, 3] == 0) { # invariant for the most frequent allele
      nonInformative <- TRUE
      reason <- 4
    } else if (I4[, 4] == 0 && I4[, 3] == 0) { # invariant for the second most frequent allele
      nonInformative <- TRUE
      reason <- 5
    }

    if (any(nonInformative)) {
      cat(paste(c(INFO[nonInformative, c(1:2, 6)], reason), collapse = "\t"),
        file = omittedSites,
        sep = "\n",
        append = TRUE
      )
    } else {
      cat(
        paste(
          c(
            INFO[!nonInformative, c(1:2, 6)],
            strsplit(ALLELES, ",", fixed = TRUE)[[1]][as.numeric(majorAlleles) + 1]
          ),
          collapse = "\t"
        ),
        file = includedSites,
        sep = "\n",
        append = TRUE
      )
    }

    return(SNP[!nonInformative, ])
  }




  #############################
  #############################
  ####   Resolve outputs   ####
  #############################
  #############################

  outputs <- ResolveOutput(SNP = SNP, filename = filename, chunk = chunk)

  filename <- outputs[[1]]
  chunk <- outputs[[2]]
  origChunk <- outputs[[3]]
  omittedSites <- outputs[[4]][1]
  includedSites <- outputs[[4]][2]
  sampleNames <- outputs[[4]][3]

  # initialize loci placement files
  cat("## Reasons for omitting loci:\n## 1 - Marker has fewer than 2 alleles representing substitutions\n## 2 - Required homozygous individuals for the 2 most frequent alleles are not present\n## 3 - The second most frequent allele is found only in one heterozygous individual\n## 4 - Dataset is invariant for the most frequent allele\n## 5 - Dataset is invariant for the allele listed as the first ALT in the vcf input\nCHROM\tPOS\tQUAL\tREASON\n", file = omittedSites, append = FALSE)
  cat("CHROM\tPOS\tQUAL\tallele0\tallele2\n", file = includedSites, append = FALSE)
  cat("sampleNames\n", file = sampleNames, append = FALSE)


  ##############################
  ####  Convert class vcfR  ####
  ##############################

  if (inherits(SNP, "vcfR")) {
    INFO <- vcfR::getFIX(SNP, getINFO = FALSE)
    SNP <- vcfR::extract.gt(SNP)


    SNP <- ResolveGenotypes()


    lims <- c(seq(0, nrow(SNP), by = chunk), nrow(SNP))
    filesizes <- diff(lims)

    for (i in 1:ceiling(nrow(SNP) / chunk)) {
      write.table(
        cbind(
          rep("S", filesizes[i]),
          SNP[(lims[i] + 1):lims[i + 1], ]
        ),
        file = filename[i],
        col.names = FALSE, row.names = FALSE, quote = FALSE, sep = ""
      )
    }
  }


  ###################################
  ####  Convert class character  ####
  ###################################

  if (inherits(SNP, "character")) {
    # initialise file connections
    if (endsWith(SNP, "gz")) {
      infile <- gzcon(file(SNP, open = "rb"))
    } else {
      if (endsWith(SNP, "vcf")) {
        infile <- file(SNP, open = "rt")
      }
    }


    outfile <- file(paste0(filename, "-", formatC(1, width = 3, flag = 0), ".txt"),
      open = "wt"
    )

    on.exit(close(infile))
    on.exit(close(outfile), add = TRUE)

    nLines <- 1
    nFiles <- 2
    Marker <- readLines(infile, n = 1, skipNul = TRUE)

    # skip meta information
    while (substr(Marker, 1, 1) == "#") {
      previousMarker <- Marker
      Marker <- readLines(infile, n = 1)
    }

    # write sample names
    previousMarker <- unlist(strsplit(previousMarker, split = "\t"))
    cat(previousMarker[10:length(previousMarker)], file = sampleNames, sep = "\n", append = TRUE)

    # check ChosenInds
    nInds <- length(previousMarker) - 9
    if (ChosenInds[1] == "all") {
      ChosenInds <- rep(TRUE, nInds)
    }



    if (!(inherits(ChosenInds, "logical") || inherits(ChosenInds, "numeric") || inherits(ChosenInds, "integer")) ||
      any(ChosenInds > nInds) ||
      (inherits(ChosenInds, "logical") && length(ChosenInds) != nInds)) {
      warning("There are ", nInds, " individuals in the vcf file. Converting to diem for all.")
      ChosenInds <- rep(TRUE, nInds)
    }


    # read and resolve genotypes
    while (length(Marker) > 0) { ### does not stop here. why?
      Marker <- matrix(unlist(strsplit(Marker, split = "\t")), nrow = 1)
      INFO <- Marker[, 1:7, drop = FALSE]
      SNP <- Marker[, 10:ncol(Marker), drop = FALSE]
      SNP <- SNP[, ChosenInds, drop = FALSE]
      SNP <- sub(":.+", "", SNP, perl = TRUE)


      SNP <- ResolveGenotypes()

      if (length(SNP) > 0) {
        writeLines(
          paste0(
            "S",
            paste(SNP, collapse = "")
          ),
          con = outfile
        )
        nLines <- nLines + 1
        if (chunk != 1 && nLines > chunk) {
          close(outfile)
          outfile <- file(paste0(filename, "-", formatC(nFiles, width = 3, flag = 0), ".txt"),
            open = "wt"
          )
          message("Done with chunk ", nFiles - 1)
          nFiles <- nFiles + 1
          nLines <- 1
        }
      }
      Marker <- readLines(infile, n = 1)
    }
  }
}
