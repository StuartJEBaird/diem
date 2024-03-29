#' Convert vcf files to diem format
#'
#' Reads vcf files and writes genotypes of the most frequent alleles based on
#' chromosome positions to diem format.
#'
#' @param SNP character vector with a path to the '.vcf' or '.vcf.gz' file, or an \code{vcfR}
#'     object. Diploid data are currently supported.
#' @param filename character vector with a path where to save the converted genotypes.
#' @param chunk numeric indicating by how many markers should the result be split into
#'     separate files. \code{chunk = 1} saves all markers into one file.
#' @param requireHomozygous logical whether to require the marker to have at least one
#'     homozygous individual for each allele.
#' @param ... additional arguments.
#'
#' @details Importing vcf files larger than 1GB is not recommended. The path to the
#'      vcf file in \code{SNP} reads the file line by line, and might be a solution for
#'      very large genomic datasets.
#'
#'      The number of files \code{vcf2diem} creates depends on the \code{chunk} argument
#'      and class of the \code{SNP} object. 
#'      
#'       * When \code{chunk = 1}, one output file will be created.
#'       * Values of \code{chunk < 100} are interpreted as the number of files into which to
#'      split data in \code{SNP}. For \code{SNP} object of class \code{vcfR}, the number
#'      of markers per file is calculated from the dimensions of \code{SNP}. When class
#'      of \code{SNP} is \code{character}, the number of markers per file is approximated
#'      from a model with a message. If this number is inappropriate for the expected
#'      output, provide the intended number of markers per file in \code{chunk} greater
#'      than 100. \code{vcf2diem} will scan the whole input \code{SNP} file, creating
#'      additional output files until the last line in \code{SNP} is reached.
#'       * Values of \code{chunk >= 100} mean that each output file
#'      in diem format will contain \code{chunk} number of lines with the data in \code{SNP}.
#'     
#'      When the vcf file contains markers non-informative for genome polarisation, those 
#'      those are removed and listed in a file *omittedLoci.txt* in the working directory. 
#'      The omitted loci are identified by their information in the CHROM and POS columns.
#'      The CHROM and POS information for loci included in the converted file are in 
#'      *includedLoci.txt*.
#' @return No value returned, called for side effects.
#' @importFrom vcfR getFIX extract.gt
#' @importFrom tools file_ext file_path_sans_ext
#' @export
#' @author Natalia Martinkova
#' @author Filip Jagos <521160@mail.muni.cz>
#' @examples
#' \dontrun{
#' # vcf2diem will write files to a working directory or a specified folder
#' # make sure the working directory or the folder are at a location with write permission
#'   myofile <- system.file("extdata", "myotis.vcf", package = "diemr")
#'   myovcf <- vcfR::read.vcfR(myofile)
#'
#'   vcf2diem(SNP = myofile, filename = "test1")
#'   vcf2diem(SNP = myofile, filename = "test2", chunk = 3)
#'   vcf2diem(SNP = myovcf, filename = "test3")
#'   vcf2diem(SNP = myovcf, filename = "test4", chunk = 3)
#' }
vcf2diem <- function(SNP, filename, chunk = 1L, requireHomozygous = TRUE, ...) {
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



  ######################################
  ######################################
  ####  Declare internal functions  ####
  ######################################
  ######################################


  ################################################
  ####  Resolve file structure and chunk size ####
  ################################################

  ResolveOutput <- function(SNP, filename, chunk) {
    # use only the first filename
    if (length(filename) > 1){
      filename <- filename[1]
    }
    # strip file extension and replace with a file number and txt
    fileext <- tools::file_ext(filename)
    filepath <- tools::file_path_sans_ext(filename)
    filename <- filepath
	# file names for loci positions
    lociFiles <- c(paste0(filepath, "-omittedLoci.txt"),
				   paste0(filepath, "-includedLoci.txt"))
	# estimate chunk size			        
    if (chunk < 100) {
      if (any(class(SNP) == "vcfR")) {
        filename <- paste0(filename, "-",
        formatC(1:chunk, width = nchar(chunk), flag = 0),
        ".",
        ifelse(fileext == "", "txt", fileext)
        )
        chunk <- unname(ceiling(nrow(SNP@fix) / chunk))
      } else {
        filesize <- file.size(SNP)[1]
        # model prediction on possible number of markers per chunk
        chunk <- ceiling(10^(-1.2387 + 0.7207 * log10(filesize)) / chunk)
        message("Expecting to include ", chunk, " markers per diem file.\nIf you expect more markers in the file, provide suitable chunk size.")
      }
    }
    if (chunk >= 100) {
      message("Expecting to include ", chunk, " markers per diem file.")
    }
    return(list(filename, chunk, lociFiles))
  }


  ###############################################
  ####  Resolve indels and multiple alleles  ####
  ###############################################

  ResolveGenotypes <- function() {
    # attempt to resolve sites with indels in ALT
    INFO[, 5] <- sub("\\*", "AA", INFO[, 5])
    resolvable <- lapply(strsplit(INFO[, 5], ",", fixed = TRUE),
      FUN = \(x) which(nchar(x) == 1)
    )
    resolvable[lengths(resolvable) == 0] <- NA

    # remove markers with unresolvable indels
    SNP[nchar(INFO[, 4]) > 1, ] <- NA # REF alleles contain insertion
    SNP[indels <- sapply(resolvable, FUN = anyNA, simplify = TRUE), ] <- NA # ALT alleles do not contain substitutions

    # resolve multiallelic markers
    multiallelic <- which(grepl(",", INFO[, 5]) & !indels)
    for (i in multiallelic) {
      majorAllele <- resolvable[[i]][which.max(sapply(resolvable[[i]],
        FUN = \(allele) sum(unlist(gregexpr(allele, SNP[i, ])) > 0, na.rm = TRUE)
      ))]
      SNP[i, ] <- sub(
        pattern = paste0("[", paste(c(1:9)[-majorAllele], collapse = ""), "]"),
        replacement = "\\.",
        x = SNP[i, ]
      )
      SNP[i, ] <- gsub(majorAllele, "1", SNP[i, ])
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
	nonInformative <- switch(requireHomozygous + 1, (I4[, 2] == 0 & I4[, 4] == 0), (I4[, 2] == 0 | I4[, 4] == 0))  | # only heterozygots, | requiring one homozygous ind for each allele 
					  (I4[, 3] == 1 & I4[, 2] == 0) | (I4[, 3] == 1 & I4[, 4] == 0) | # only one heterozygous individual
					  (I4[, 4] == 0 & I4[, 3] == 0) | # only homozygots for the reference allele
					  (I4[, 2] == 0 & I4[, 3] == 0) # only homozygots for the alternative allele
	
	if(any(nonInformative)){
	  cat(paste(INFO[nonInformative, 1:2], collapse = "\t"),
	    file = omittedLoci, 
	    sep = "\n", 
		append = TRUE
	  )

	} else {
	  cat(paste(INFO[!nonInformative, 1:2], collapse = "\t"),
	    file = includedLoci, 
	    sep = "\n", 
		append = TRUE
	  )
	}
	
    return(SNP[!nonInformative, ])
  }




  #############################
  #############################
  ####  Convert genotypes  ####
  #############################
  #############################

  outputs <- ResolveOutput(SNP = SNP, filename = filename, chunk = chunk)

  filename <- outputs[[1]]
  chunk <- outputs[[2]]
  omittedLoci <- outputs[[3]][1]
  includedLoci <- outputs[[3]][2]

  # initialize loci placement files
  cat("CHROM\tPOS\n", file = omittedLoci, append = FALSE)
  cat("CHROM\tPOS\n", file = includedLoci, append = FALSE)



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


    outfile <- file(paste0(filename, "-01.txt"), open = "wt")

    on.exit(close(infile))
    on.exit(close(outfile), add = TRUE)

    nLines <- 1
    nFiles <- 2
    Marker <- readLines(infile, n = 1)

    # skip meta information
    while (substr(Marker, 1, 1) == "#") {
      Marker <- readLines(infile, n = 1)
    }


    # read and resolve genotypes
    while (length(Marker) > 0) { ### does not stop here. why?
      Marker <- matrix(unlist(strsplit(Marker, split = "\t")), nrow = 1)
      INFO <- Marker[, 1:7, drop = FALSE]
      SNP <- Marker[, 10:ncol(Marker), drop = FALSE]
      SNP <- sub(":.+", "", SNP, perl = TRUE)

      SNP <- ResolveGenotypes()

	  if(length(SNP) > 0){
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
