#' Diagnostic Index Expectation Maximisation
#'
#' Estimates how to assign alleles in a genome to maximise the distinction between two
#' unknown groups of individuals. Using expectation maximisation (EM) in likelihood
#' framework, \code{diem} provides marker
#' polarities for importing data, their likelihood-based diagnostic index and its support
#' for all markers, and hybrid indices for all individuals.
#'
#'
#' @param files character vector with paths to files with genotypes.
#' @param ploidy list of length equal to length of \code{files}. Each element of the list
#'    contains a numeric vector with ploidy numbers for all individuals specified in
#'    the \code{ChosenInds}.
#' @param markerPolarity \code{FALSE} or list of logical vectors.
#' @param ChosenInds numeric vector of indices of individuals to be included in the analysis.
#' @param epsilon numeric, specifying how much the hypothetical diagnostic markers should
#' 			contribute to the likelihood calculations. Must be in \code{[0,1)}, keeping
#'  		tolerance setting of the \code{R} session in mind.
#' @param verbose logical or character with path to directory where run diagnostics will
#'    be saved.
#' @param nCores numeric. Number of cores to be used for parallelisation. Must be
#'     \code{nCores = 1} on Windows.
#' @param maxIterations numeric.
#' @param ... additional arguments.
#' @export
#' @details  Given two alleles of a marker, one allele can belong to one side of a barrier
#'   to geneflow and the other to the other side. Which allele belongs where is a non-trivial
#'   matter. A marker state in an individual can be encoded as 0 if the individual is
#'   homozygous for the first allele, and 2 if the individual is homozygous for the second
#'   allele. Marker polarity determines how the marker will be imported. Marker polarity
#'   equal to \code{FALSE} means that the marker will be imported as-is. A marker with
#'   polarity equal to \code{TRUE} will be imported with states 0 mapped as 2 and states 2
#'   mapped as 0, in effect switching which allele belongs to which side of a barrier to
#'   geneflow.
#'
#'   When \code{markerPolarity = FALSE}, \code{diem} uses random null polarities to
#'   initiate the EM algorithm. To fix the null polarities, \code{markerPolarity} must be
#'   a list of length equal to the length of the \code{files} argument, where each element
#'   in the list is a logical vector of length equal to the number of markers (rows) in
#'   the specific file.
#'
#'   Ploidy needs to be given for each compartment and for each individual. For example,
#'   for a dataset of three diploid mammal males consisting of an autosomal
#'   compartment, an X chromosome
#'   compartment and a Y chromosome compartment, the ploidy list would be
#'   \code{ploidy = list(rep(2, 3), rep(1, 3), rep(1, 3)}. If the dataset consisted of
#'   one male and two females,
#'   ploidy for the sex chromosomes should be vectors reflecting that females have two X
#'   chromosomes, but males only one, and females have no Y chromosomes:
#'   \code{ploidy = list(rep(2, 3), c(1, 2, 2), c(1, 0, 0))}.
#' @note To ensure that the data input format of the genotype files, ploidies and individual
#'   selection are readable for \code{diem}, first use \code{\link{CheckDiemFormat}}.
#'   Fix all errors, and run \code{diem} only once the checks all passed.
#' @seealso \code{\link{CheckDiemFormat}}
#' @return A list including suggested marker polarities, diagnostic indices and support for all
#' 		markers, four genomic state counts matrix for all individuals, and polarity changes 
#'      for the EM iterations.
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
#'
#' # run diem
#' \dontrun{
#' # diem will store temporal files during EM iterations
#' # prior to running diem, set the working directory to a location with write permission
#' fit <- diem(files = inputFile, ChosenInds = inds, ploidy = ploidies, nCores = 1)
#'
#' # run diem with fixed null polarities
#' fit2 <- diem(
#'   files = inputFile, ChosenInds = inds, ploidy = ploidies, nCores = 1,
#'   markerPolarity = list(c(TRUE, FALSE, TRUE))
#' )
#' }
#' @importFrom stats aggregate na.pass
#' @importFrom utils head read.table write.table
#' @importFrom parallel detectCores


diem <- function(files, ploidy = list(2), markerPolarity = FALSE, ChosenInds,
                 epsilon = 0.99999, verbose = FALSE, nCores = parallel::detectCores() - 1,
                 maxIterations = 50, ...) {
  if (is.na(nCores)) nCores <- 1
  if (nCores > parallel::detectCores()) nCores <- parallel::detectCores()


  ###################################
  # prepare diagnostics folder
  ###################################
  if (inherits(verbose, "character")) {
    folder <- verbose
    verbose <- TRUE
    lfolder <- paste0(folder, "/likelihood")
    
    if (!dir.exists(lfolder)) dir.create(lfolder)
    
    logfile <- paste0(folder, "/log.txt")
  } else {
    folder <- getwd()
  }
  if (verbose & folder == getwd()) {
    if (!dir.exists("diagnostics")) dir.create("diagnostics")
    folder <- "diagnostics"
    lfolder <- "diagnostics/likelihood"
    if (!dir.exists(lfolder)) dir.create(lfolder)
    logfile <- paste0(folder, "/log.txt")
  }
  if (!verbose) {
    lfolder <- "likelihood"
    if (!dir.exists(lfolder)) dir.create(lfolder)
  }
  for (i in 1:length(files)) {
    cat("changedMarkers\ttime\titeration\n", file = paste0(lfolder, "/MarkersWithChangedPolarity", i, ".txt"))
  }


  ##############################################
  ##############################################
  ######    declare internal fucntions    ######
  ##############################################
  ##############################################


  ############################
  #  Imports a genotype file
  ############################

  sImport <- function(file) {
    genotypes <- read.table(file = file, as.is = TRUE)
    genotypes <- as.data.frame(strsplit(unlist(genotypes), split = ""))[-1, ][ChosenInds, ]
    if (any(grepl("U", genotypes))) {
      genotypes[genotypes == "U"] <- "_"
    }
    return(genotypes)
  }



  #############################################
  #  Calculates I4, individual state counts
  #############################################

  GetI4ofOneCompartment <- function(file, changePolarity = FALSE, ...) {
    # read file with genotypes, individuals are in rows, markers in columns
    genotypes <- sImport(file = file)
    n <- ncol(genotypes)

    if (length(changePolarity) != n) {
      warning("Marker polarities for compartment ", file, " are not equal to the number of markers (", n, "). Using random polarities.")
      n2 <- floor(n / 2)
      changePolarity <- sample(c(rep(TRUE, n2), rep(FALSE, n - n2)),
        size = n
      )
    }

    # polarise markers
    genotypes <- apply(cbind(changePolarity, t(genotypes)),
      MARGIN = 1,
      FUN = function(x) emPolarise(origM = x[2:length(x)], changePolarity = x[1])
    )

    # calculate state counts
    res <- list(
      I4 = t(apply(genotypes,
        MARGIN = 1,
        FUN = sStateCount
      )),
      changePolarity = changePolarity,
      compartmentSize = n
    )

    return(res)
  }




  ##################################################
  #  Polarises one marker based on DI comparison
  ##################################################

  PolariseAndRankMarker <- function(origM, changePolarity = FALSE, ...) {
    # default polarity
    Polarity <- FALSE


    # indices for original and polarised positions
    indicesM <- as.numeric(factor(origM, levels = c("_", "0", "1", "2")))
    indicesM2 <- as.numeric(factor(origM, levels = c("_", "2", "1", "0")))

    # likelihood for original and polarised marker
    Ans <- c(
      sum(FlatLogI4[indicesM + offsetsM], na.rm = TRUE),
      sum(FlatLogI4[indicesM2 + offsetsM], na.rm = TRUE)
    )
    MaxAns <- max(Ans)

    # reverse order if user wants to change polarity
    if (changePolarity) Ans <- rev(Ans)

    # likelihood-driven polarity
    if (Ans[2] > Ans[1]) Polarity <- TRUE

    # changing marker polarity
    if (!Polarity) polarisedM <- ""
    if (Polarity & !changePolarity) {
      polarisedM <- c("_", "2", "1", "0")[match(origM, c("_", "0", "1", "2"))]
    }
    if (Polarity & changePolarity) polarisedM <- origM

    return(list(
      Polarity = Polarity,
      DI = MaxAns,
      Support = MaxAns - min(Ans),
      polarisedM = polarisedM
    ))
  }





  #####################################
  #   Update all marker polarities
  #####################################

  RelabelCompartment <- function(file, changePolarity = TRUE, compartment, ...) {

    # read file with genotypes
    genotypes <- sImport(file)

    # declare variables
    nMarkers <- ncol(genotypes)
    Polarity <- list()


    # populate changePolarity
    if (length(changePolarity) == 1) {
      warning("changePolarity from within RelabelCompartment is of length 1. Using it for all markers")
      changePolarity <- rep(changePolarity, nMarkers)
    }

    # calculate new polarity based on likelihood
    for (i in 1:nMarkers) {
      Polarity[[i]] <- PolariseAndRankMarker(genotypes[, i], changePolarity = changePolarity[i], ...)
    }

    newPolarity <- unname(unlist(lapply(Polarity, function(x) x[[1]])))

    # calculate change of state counts
    deltaI4 <- matrix(0, ncol = 4, nrow = length(ChosenInds))
    if (sum(newPolarity) == 0) {
      changedMarkers <- data.frame(changedMarkers = NA, time = Sys.time(), iteration = IterationCount)
    }
    if (sum(newPolarity) > 0) {
      changedSC <- lapply(Polarity[newPolarity], function(x) x[[4]])
      changedSC <- matrix(unlist(changedSC), nrow = sum(newPolarity), byrow = TRUE)

      deltaI4 <- t(apply(changedSC, MARGIN = 2, sStateCount))
      changedMarkers <- data.frame(
        changedMarkers = which(newPolarity),
        time = Sys.time(),
        iteration = IterationCount
      )
    }

    write.table(changedMarkers,
      sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE,
      file = paste0(lfolder, "/MarkersWithChangedPolarity", compartment, ".txt")
    )
    write.table(data.frame(
      newPolarity = newPolarity,
      DI = unlist(lapply(Polarity, function(x) x[[2]])),
      Support = unlist(lapply(Polarity, function(x) x[[3]]))
    ),
    file = paste0(lfolder, "/MarkerDiagnostics", IterationCount, "-", compartment, ".txt"), row.names = FALSE, sep = "\t"
    )

    return(list(
      newPolarity = newPolarity,
      deltaI4 = deltaI4
    ))
  }



  ###################################################
  #  Deep check if marker polarity changes repeat
  ###################################################

  CheckDuplicated <- function(changedMarkers) {
    res <- list(
      ActualLimitCycle = logical(),
      iterations = numeric()
    )
    # check if the same set of markers changed polarity before
    markers <- aggregate(changedMarkers[, "changedMarkers"] ~ changedMarkers[, "iteration"],
      FUN = toString, na.action = na.pass
    )
    duplicatedMarkerSets <- duplicated(markers[, 2])
    res[[1]] <- any(duplicatedMarkerSets)

    # check if the overall state is the same as before
    if (res[[1]]) {
      iterations <- which(duplicatedMarkerSets != rev(duplicated(rev(markers[, 2]))))
      if (length(iterations) != 2) stop("Something is very wrong with CheckDuplicated. Infinite loop risk. Run diem with different null polarities")
      res[[1]] <- all(I4changes[[iterations[1]]] == I4changes[[iterations[2]]])
      if (verbose & res[[1]]) {
        cat("\nIterations ", iterations, "called the same marker sets to change and the overall I4 was identical in both cases. Congratulations, the EM converged.\n",
          file = logfile, append = TRUE
        )
      }
      # include previous iterations with identical state in result
      res[[2]] <- iterations
    }

    return(res)
  }





  ####################################################
  #  Summarise marker polarity across compartments
  ####################################################

  summariseMarkersFromCompartments <- function() {
    res <- data.frame(changedMarkers = numeric(), time = character(), iteration = numeric())
    compartments <- rep(0, length(files))
    for (i in 1:length(files)) {
      temp <- read.table(paste0(lfolder, "/MarkersWithChangedPolarity", i, ".txt"),
        sep = "\t", header = TRUE, as.is = TRUE
      )
      temp$changedMarkers <- temp$changedMarkers + sum(compartmentSizes * compartments)
      res <- rbind(res, temp)
      compartments[i] <- 1
    }
    res <- res[order(res$iteration, res$changedMarkers), ]
    return(res)
  }


  summariseDIacrossCompartments <- function(iteration) {
    res <- data.frame(Marker = numeric(), newPolarity = logical(), DI = numeric(), Support = numeric())
    compartments <- rep(0, length(files))
    for (i in 1:length(files)) {
      temp <- read.table(paste0(lfolder, "/MarkerDiagnostics", iteration, "-", i, ".txt"),
        sep = "\t", header = TRUE, as.is = TRUE
      )
      temp$Marker <- 1:nrow(temp)
      temp$Marker <- temp$Marker + sum(compartmentSizes * compartments)
      temp$newPolarity <- unlist(PolarityChanges[[iteration]][i])
      res <- rbind(res, temp)
      compartments[i] <- 1
    }
    return(res)
  }


  ##############################
  #   Small data correction
  ##############################

  SmallDataCorrection <- function() {
    minI4 <- min(I4)
    # return(max(c(0, SmallDataErrorTermGoesTo - minI4)))
    res <- ifelse(test = minI4 > SmallDataErrorTermGoesTo, 
      yes = -1,
      no = ifelse(minI4 == SmallDataErrorTermGoesTo,
        yes = 0,
        no = 1))
    return(res)
  }




  ######################################################
  #   Update state counts with their polarity change   #
  ######################################################

  UpdateM4withDelta <- function(M4, deltaM4) {
    newM4 <- M4
    newM4[, 2] <- newM4[, 2] - deltaM4[, 4] + deltaM4[, 2]
    newM4[, 4] <- newM4[, 4] - deltaM4[, 2] + deltaM4[, 4]
    return(newM4)
  }



  ##########################
  ##########################
  ##   end of functions   ##
  ##########################
  ##########################





  ########################
  #  declare variables   #
  ########################

  IterationCount <- 1
  SmallDataErrorTermGoesTo <- 1
  ActualLimitCycle <- FALSE
  MarkersWithChangedPolarity <- data.frame(
    iteration = numeric(),
    changedMarkers = numeric(),
    time = character()
  )
  nMarkersWithChangedPolarity <- NULL
  PolarityChanges <- I4changes <- HIchanges <- SmallDataI4errorTermChanges <- list()
  if (verbose) {
    message("Starting individual state counts... ", Sys.time())
  }

  # Check ChosenInds
  if (missing(ChosenInds)) {
    warning("ChosenInds is missing, using all individuals")
    ChosenInds <- 1:(nchar(readLines(files[1])[1]) - 1)
  }

  # Check markerPolarity
  if (!inherits(markerPolarity, "list")) {
    warning("markerPolarity is not a list, using random marker polarities")
    markerPolarity <- list(FALSE)[rep(1, length(files))]
  }


  # I4compartments is a list
  I4compartments <- parallel::mclapply(
    X = as.character(1:length(files)), mc.cores = nCores,
    FUN = function(x) {
      GetI4ofOneCompartment(
        file = files[as.numeric(x)],
        changePolarity = markerPolarity[[as.numeric(x)]],
        ChosenInds = ChosenInds, ...
      )
    }
  )
  markerPolarity <- lapply(I4compartments, FUN = function(x) x[[2]])
  changePolarity <- markerPolarity
  compartmentSizes <- unlist(lapply(I4compartments, FUN = function(x) x[[3]]))
  # I4 is a matrix
  I4 <- Reduce("+", lapply(I4compartments, FUN = function(x) x[[1]]))
  nMarkers <- sum(I4[1, ])


  # allele counts accross compartments
  I4compartments <- lapply(I4compartments, FUN = function(x) x[[1]])
  A4compartments <- Map("*", I4compartments, ploidy)



  # small data correction
  I4errorTermDistributor <- compartmentSizes / nMarkers
  A4errorTermDistributor <- Map("*", ploidy, I4errorTermDistributor)
  SmallDataI4errorTerm <- SmallDataCorrection()
  SmallDataA4errorTerm <- Map("*", A4errorTermDistributor, SmallDataI4errorTerm)
  I4 <- I4 + SmallDataI4errorTerm
  A4compartments <- Map("+", A4compartments, Map("*", A4errorTermDistributor, SmallDataI4errorTerm))
  A4 <- Reduce("+", A4compartments)


  ###############################
  # write initial diagnostics
  ###############################
  if (verbose) {
    cat("nMarkers: ", nMarkers, "\nnIndividuals: ", nrow(I4), "\n\n", file = logfile, append = FALSE)
    cat("Top rows in I4 with null polarities:\n", file = logfile, append = TRUE)
    capture.output(I4[1:min(7, length(ChosenInds)), ], file = logfile, append = T)
    cat("\n\nCorrection for small data in I4: ", SmallDataI4errorTerm, "\n\n", file = logfile, append = TRUE)
    cat("Top rows in A4 with null polarities:\n", file = logfile, append = TRUE)
    capture.output(A4[1:min(7, length(ChosenInds)), ], file = logfile, append = T)
    cat("Null marker polarities for each compartment:\n", file = paste0(folder, "/NullMarkerPolarities.txt"), append = FALSE)
    i <- 0
    while (i < length(markerPolarity)) {
      cat(markerPolarity[[i <- i + 1]], "\n", file = paste0(folder, "/NullMarkerPolarities.txt"), append = TRUE)
    }
    message("End individual state counts at ", Sys.time())
  }


  ####################
  # initialise EM
  ####################

  # Hybrid index from allele counts
  OriginalHI <- apply(A4, MARGIN = 1, FUN = pHetErrOnStateCount)[1, ]

  # I4 corrected for epsilon-weighted hybrid index
  V4 <- ModelOfDiagnostic(
    I4 = I4,
    OriginalHI = OriginalHI,
    epsilon = epsilon,
    verbose = verbose, ...
  )

  # DI for state counts
  FlatLogI4 <- apply(V4,
    MARGIN = 1,
    FUN = function(x) x / (sum(x) - 3 * SmallDataI4errorTerm)
  )
  FlatLogI4 <- log(c(FlatLogI4))


  # indices differentiating individuals in DI vector FlatLogI4
  offsetsM <- seq(from = 0, to = 4 * (length(ChosenInds) - 1), by = 4)



  ###############################
  # write diagnostics
  ###############################
  if (verbose) {
    cat("\nTop rows in V4 for epsilon = ", epsilon, " with null polarisation:\n", file = logfile, append = T)
    capture.output(round(V4[1:min(7, length(ChosenInds)), ], 2), file = logfile, append = T)
    cat("\nFlatLogI4 with null polarisation for first 28 values:\n", file = logfile, append = T)
    capture.output(round(FlatLogI4[1:min(28, length(FlatLogI4))], 2), file = logfile, append = TRUE)
    message("Initialised EM at ", Sys.time())
  }



  #########################
  # EM cycle
  #########################

  while (!ActualLimitCycle & (IterationCount <= maxIterations)) {
    message("diadem iteration: ", IterationCount)


    ###################################
    # start UpdateEMstateGivenI4
    ###################################

    I4changes[[IterationCount]] <- I4
    HIchanges[[IterationCount]] <- OriginalHI
    SmallDataI4errorTermChanges[[IterationCount]] <- SmallDataI4errorTerm
    PolarityChanges[[IterationCount]] <- changePolarity

    # start testing marker polarity
    cat("Start testing marker polarity\n")

    res.polarity <- parallel::mclapply(
      X = as.character(1:length(files)), mc.cores = nCores,
      FUN = function(x) {
        RelabelCompartment(
          file = files[as.numeric(x)],
          changePolarity = changePolarity[[as.numeric(x)]], compartment = as.numeric(x), ...
        )
      }
    )

    deltaI4compartments <- lapply(res.polarity, FUN = function(x) x[[2]])
    deltaA4compartments <- Map("*", deltaI4compartments, ploidy)
    deltaI4 <- Reduce("+", deltaI4compartments)

    newPolarity <- lapply(res.polarity, FUN = function(x) x[[1]])


    # check the number of iterations
    MarkersWithChangedPolarity <- summariseMarkersFromCompartments()
    nMarkersWithChangedPolarity <- c(table(MarkersWithChangedPolarity$iteration))
    PotentialLimitCycle <- any(duplicated(nMarkersWithChangedPolarity))
    if (PotentialLimitCycle) {
      ExistingLCcandidate <- CheckDuplicated(MarkersWithChangedPolarity)
      ActualLimitCycle <- ifelse(ExistingLCcandidate[[1]], TRUE, FALSE)
    }


    cat("Potential last iteration:", PotentialLimitCycle, "\nActual last iteration:", ActualLimitCycle, "\n")

    # write diagnostics
    if (verbose) {
      cat("\n\n#########################\n# Iteration ", IterationCount,
        "\n#########################\n\n",
        file = logfile, append = TRUE
      )
      write.table(MarkersWithChangedPolarity,
        file = paste0(folder, "/MarkersWithChangedPolarity.txt"),
        sep = "\t", row.names = FALSE, append = FALSE
      )
    }

    # test polarity change
    ComparePolarity <- lapply(1:length(changePolarity),
      FUN = function(x) {
        data.frame(
          oldPolarity = changePolarity[[x]],
          newPolarity = newPolarity[[x]]
        )
      }
    )
    changePolarity <- lapply(1:length(changePolarity),
      FUN = function(x) rowSums(ComparePolarity[[x]]) == 1
    )


    # update I4, A4
    newI4 <- UpdateM4withDelta(M4 = I4, deltaM4 = deltaI4)
    I4 <- newI4 - SmallDataI4errorTerm
    newA4compartments <- Map(UpdateM4withDelta, M4 = A4compartments, deltaM4 = deltaA4compartments)
    A4compartments <- Map("-", newA4compartments, SmallDataA4errorTerm)

    # check for smallest values in I4 and correct if necessary
    SmallDataI4errorTerm <- SmallDataCorrection()
    SmallDataA4errorTerm <- Map("*", A4errorTermDistributor, SmallDataI4errorTerm)
    I4 <- I4 + SmallDataI4errorTerm
    A4compartments <- Map("+", A4compartments, SmallDataA4errorTerm)
    A4 <- Reduce("+", A4compartments)


    # update FlatLogI4 - DI for state counts
    modelfolder <- paste0("modelDiagnostics", IterationCount)
    dir.create(modelfolder)
    OriginalHI <- apply(A4, MARGIN = 1, FUN = pHetErrOnStateCount)[1, ]
    V4 <- ModelOfDiagnostic(
      I4 = I4,
      OriginalHI = OriginalHI,
      epsilon = epsilon,
      verbose = modelfolder, ...
    )
    FlatLogI4 <- apply(V4,
      MARGIN = 1,
      FUN = function(x) x / sum(x)
    )
    FlatLogI4 <- log(c(FlatLogI4))

    ###################################
    # end UpdateEMstateGivenI4
    ###################################


    if (verbose) {
      cat("\n\nPolarity change at iteration:", IterationCount, "\n", file = logfile, append = TRUE)
      temp <- table(as.data.frame(do.call(rbind, ComparePolarity)))
      capture.output(temp, file = logfile, append = TRUE)
      cat("\n\nPolarities to be changed at the next iteration:\n", file = logfile, append = TRUE)
      temp <- table(do.call(c, changePolarity))
      capture.output(temp, file = logfile, append = TRUE)
      cat("\n\ndeltaI4:\n", file = logfile, append = TRUE)
      capture.output(deltaI4[1:min(7, length(ChosenInds)), ], file = logfile, append = TRUE)
      cat("\n\nUpdated I4:\n", file = logfile, append = TRUE)
      capture.output(I4[1:min(7, length(ChosenInds)), ], file = logfile, append = TRUE)
      cat("\n\nCorrection for small data in I4: ", SmallDataI4errorTerm, "\n\n", file = logfile, append = TRUE)
      cat("\n\nUpdated A4:\n", file = logfile, append = TRUE)
      capture.output(A4[1:min(7, length(ChosenInds)), ], file = logfile, append = TRUE)
      cat("\n\nUpdated V4:\n", file = logfile, append = TRUE)
      capture.output(round(V4[1:min(7, length(ChosenInds)), ], 2), file = logfile, append = TRUE)
      cat("\n\nUpdated FlatLogI4:\n", file = logfile, append = TRUE)
      capture.output(round(FlatLogI4[1:min(28, length(FlatLogI4))], 2), file = logfile, append = TRUE)
    }


    # end one iteration
    IterationCount <- IterationCount + 1
  } # end while for cycle limit



  DI <- summariseDIacrossCompartments(iteration = ExistingLCcandidate[[2]][1] + 1)
  I4 <- I4changes[[ExistingLCcandidate[[2]][1] + 1]] - SmallDataI4errorTermChanges[[ExistingLCcandidate[[2]][1] + 1]]
  rownames(I4) <- ChosenInds
  HI <- HIchanges[[ExistingLCcandidate[[2]][1] + 1]]

  if (verbose) {
    write.table(DI,
      file = paste0(folder, "/MarkerDiagnosticsWithOptimalPolarities.txt"), sep = "\t", row.names = FALSE
    )
    write.table(I4,
      file = paste0(folder, "/I4withOptimalPolarities.txt"), sep = "\t",
      row.names = TRUE, col.names = c("_", "0", "1", "2")
    )
    write.table(HI,
      sep = "\t",
      row.names = ChosenInds, file = paste0(folder, "/HIwithOptimalPolarities.txt")
    )
  }

  # delete folder with help files for likelihood calculations
  if (!verbose) {
    unlink("lfolder", recursive = TRUE)
  }

  return(list(
    markerPolarity = DI$newPolarity,
    epsilon = epsilon,
    DI = DI,
    I4 = I4,
    MarkersWithChangedPolarity = MarkersWithChangedPolarity,
    PolarityChanges = PolarityChanges
  ))
}
