#' @title Run all permutations on the specified multi-generational dataset in
#' RDS format.
#'
#' @description Run a permutation analysis, based on Monte Carlo sampling,
#' for testing the hypothesis that the number of conserved differentially
#' methylated elements (sites, tiles or both), between
#' several generations, is associated to an effect inherited from a treatment
#' and that stochastic effect can be dismissed.
#'
#' @param methylKitRDSFile a string, the name of the RDS file containing the
#' methylKit objet used for the permutation analysis. The RDS file must
#' hold a \code{list} of \code{methylRawList} entries, each
#' \code{methylRawList} contains all the \code{methylRaw} entries related to
#' one generation (first entry = first generation, second entry = second
#' generation, etc..). The number of generations must correspond to the number
#' of entries in the \code{methylKitInfo}.At least 2 generations
#' must be present to do a permutation analysis. More information can be found
#' in the methylKit package
#'
#' @param type One of the "sites","tiles" or "both" strings. Specifies the type
#' of differentially methylated elements should be returned. For
#' retrieving differentially methylated bases type="sites"; for
#' differentially methylated regions type="tiles". Default: "both".
#'
#' @param outputDir a string, the name of the directory that will contain
#' the results of the permutation or \code{NULL}. If the directory does not
#' exist, it will be created. When \code{NULL}, the results of the permutation
#' are not saved. Default: \code{NULL}.
#'
#' @param runObservedAnalysis a \code{logical}, when \code{runObservedAnalysis}
#' = \code{TRUE}, a CpG analysis on the observed dataset is done. Default:
#' \code{TRUE}.
#'
#' @param nbrPermutations, a positive \code{integer}, the total number of
#' permutations that is going to be done. Default: \code{1000}.
#'
#' @param nbrCores a positive \code{integer}, the number of cores to use when
#' processing the analysis. Default: \code{1} and always \code{1} for Windows.
#'
#' @param nbrCoresDiffMeth a positive \code{integer}, the number of cores
#' to use for parallel differential methylation calculations.Parameter
#' used for both sites and tiles analysis. The parameter
#' corresponds to the \code{num.cores} parameter in the package
#' \code{methylKit}.
#' Default: \code{1} and always \code{1} for Windows.
#'
#' @param minReads a positive \code{integer} Bases and regions having lower
#' coverage than this count are discarded. The parameter
#' correspond to the \code{lo.count} parameter in the package \code{methylKit}.
#'
#' @param minMethDiff a positive \code{double} betwwen [0,100], the absolute
#' value of methylation percentage change between cases and controls. The
#' parameter correspond to the \code{difference} parameter in
#' the package \code{methylKit}. Default: \code{10}.
#'
#' @param qvalue a positive \code{double} betwwen [0,1], the cutoff
#' for qvalue of differential methylation statistic. Default: \code{0.01}.
#'
#' @param maxPercReads a \code{double} between [0,100], the percentile of read
#' counts that is going to be used as upper cutoff. Bases ore regions
#' having higher
#' coverage than this percentile are discarded. Parameter used for both CpG
#' sites and tiles analysis. The parameter
#' correspond to the \code{hi.perc} parameter in the package \code{methylKit}.
#' Default: \code{99.9}.
#'
#' @param destrand a \code{logical}, when \code{TRUE} will merge reads on both
#' strands of a CpG dinucleotide to provide better coverage. Only advised
#' when looking at CpG methylation. Parameter used for both CpG
#' sites and tiles analysis.
#' Default: \code{FALSE}.
#'
#' @param minCovBasesForTiles a non-negative \code{integer}, the minimum
#' number of bases to be covered in a given tiling window. The parameter
#' corresponds to the \code{cov.bases} parameter in the package
#' \code{methylKit}.
#' Only used when \code{doingTiles} =
#' \code{TRUE}. Default: \code{0}.
#'
#' @param tileSize a positive \code{integer}, the size of the tiling window.
#' The parameter corresponds to the \code{win.size} parameter in
#' the package \code{methylKit}. Only
#' used when \code{doingTiles} = \code{TRUE}. Default: \code{1000}.
#'
#' @param stepSize a positive \code{integer}, the step size of tiling windows.
#' The parameter corresponds to the \code{stepSize} parameter in
#' the package \code{methylKit}. Only
#' used when \code{doingTiles} = \code{TRUE}. Default: \code{1000}.
#'
#' @param vSeed a \code{integer}, a seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used. Default: \code{-1}.
#'
#' @return TODO
#'
#' @examples
#'
#' ## Path to a methylKit RDS file
#' methylFile <- dir(system.file("extdata", package = "methylInheritance"),
#' pattern = "methylObj_002.RDS", full.names = TRUE)
#'
#' ## Run a permutation analysis
#' \dontrun{runPermutationUsingRDSFile(methylKitRDSFile = methylFile,
#' type = "sites", nbrPermutations = 10, vSeed = 2001)}
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @export
runPermutationUsingRDSFile <- function(methylKitRDSFile,
                                    type = c("both", "sites", "tiles"),
                                    outputDir = NULL,
                                    runObservedAnalysis = TRUE,
                                    nbrPermutations = 1000,
                                    nbrCores = 1,
                                    nbrCoresDiffMeth = 1,
                                    minReads = 10,
                                    minMethDiff = 10,
                                    qvalue = 0.01,
                                    maxPercReads = 99.9,
                                    destrand = FALSE,
                                    minCovBasesForTiles = 0,
                                    tileSize = 1000,
                                    stepSize = 1000,
                                    vSeed = -1) {

    ## Validate that methylKitRDSFile is an existing file
    if (!file.exists(methylKitRDSFile)) {
        stop(paste0("The file \"", methylKitRDSFile, "\" does not exist."))
    }

    ## Extract information from RDS file
    methylInfo <- readRDS(methylKitRDSFile)

    ## Call permutation analysis with the methylInfo object as an parameter
    runPermutationUsingMethylKitInfo(methylInfo, type, outputDir,
                            runObservedAnalysis,
                            nbrPermutations, nbrCores, nbrCoresDiffMeth,
                            minReads, minMethDiff, qvalue, maxPercReads,
                            destrand, minCovBasesForTiles, tileSize, stepSize,
                            vSeed)
}


#' @title Run all permutations on the specified multi-generational dataset in
#' RDS format.
#'
#' @description Run a permutation analysis, based on Monte Carlo sampling,
#' for testing the hypothesis that the number of conserved differentially
#' methylated elements (sites, tiles or both), between
#' several generations, is associated to an effect inherited from a treatment
#' and that stochastic effect can be dismissed.
#'
#' @param methylKitInfo a \code{list} of \code{methylRawList} entries, each
#' \code{methylRawList} contains all the \code{methylRaw} entries related to
#' one generation (first entry = first generation, second entry = second
#' generation, etc..). The number of generations must correspond to the number
#' of entries in the \code{methylKitInfo}.At least 2 generations
#' must be present to do a permutation analysis. More information can be found
#' in the methylKit package.
#'
#' @param type One of the "sites","tiles" or "both" strings. Specifies the type
#' of differentially methylated elements should be returned. For
#' retrieving differentially methylated bases type="sites"; for
#' differentially methylated regions type="tiles". Default: "both".
#'
#' @param outputDir a string, the name of the directory that will contain
#' the results of the permutation or \code{NULL}. If the directory does not
#' exist, it will be created. When \code{NULL}, the results of the permutation
#' are not saved. Default: \code{NULL}.
#'
#' @param runObservedAnalysis a \code{logical}, when \code{runObservedAnalysis}
#' = \code{TRUE}, a CpG analysis on the observed dataset is done. Default:
#' \code{TRUE}.
#'
#' @param nbrPermutations, a positive \code{integer}, the total number of
#' permutations that is going to be done. Default: \code{1000}.
#'
#' @param nbrCores a positive \code{integer}, the number of cores to use when
#' processing the analysis. Default: \code{1} and always \code{1} for Windows.
#'
#' @param nbrCoresDiffMeth a positive \code{integer}, the number of cores
#' to use for parallel differential methylation calculations.Parameter
#' used for both sites and tiles analysis. The parameter
#' corresponds to the \code{num.cores} parameter in the package
#' \code{methylKit}.
#' Default: \code{1} and always \code{1} for Windows.
#'
#' @param minReads a positive \code{integer} Bases and regions having lower
#' coverage than this count are discarded. The parameter
#' correspond to the \code{lo.count} parameter in the package \code{methylKit}.
#'
#' @param minMethDiff a positive \code{double} betwwen [0,100], the absolute
#' value of methylation percentage change between cases and controls. The
#' parameter correspond to the \code{difference} parameter in
#' the methylKit package. Default: \code{10}.
#'
#' @param qvalue a positive \code{double} betwwen [0,1], the cutoff
#' for qvalue of differential methylation statistic. Default: \code{0.01}.
#'
#' @param maxPercReads a \code{double} between [0,100], the percentile of read
#' counts that is going to be used as upper cutoff. Bases ore regions
#' having higher
#' coverage than this percentile are discarded. Parameter used for both CpG
#' sites and tiles analysis. The parameter
#' correspond to the \code{hi.perc} parameter in the package \code{methylKit}.
#' Default: \code{99.9}.
#'
#' @param destrand a \code{logical}, when \code{TRUE} will merge reads on both
#' strands of a CpG dinucleotide to provide better coverage. Only advised
#' when looking at CpG methylation. Parameter used for both CpG
#' sites and tiles analysis.
#' Default: \code{FALSE}.
#'
#' @param minCovBasesForTiles a non-negative \code{integer}, the minimum
#' number of bases to be covered in a given tiling window. The parameter
#' corresponds to the \code{cov.bases} parameter in the package
#' \code{methylKit}.
#' Only used when \code{doingTiles} =
#' \code{TRUE}. Default: \code{0}.
#'
#' @param tileSize a positive \code{integer}, the size of the tiling window.
#' The parameter corresponds to the \code{win.size} parameter in
#' the package \code{methylKit}. Only
#' used when \code{doingTiles} = \code{TRUE}. Default: \code{1000}.
#'
#' @param stepSize a positive \code{integer}, the step size of tiling windows.
#' The parameter corresponds to the \code{stepSize} parameter in
#' the package \code{methylKit}. Only
#' used when \code{doingTiles} = \code{TRUE}. Default: \code{1000}.
#'
#' @param vSeed a \code{integer}, a seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used. Default: \code{-1}.
#'
#' @return TODO
#'
#' @examples
#'
#' ## Load methyl information
#' data(samplesForTransgenerationalAnalysis)
#'
#' ## Run a permutation analysis
#' \dontrun{runPermutationUsingMethylKitInfo(methylKitInfo =
#' samplesForTransgenerationalAnalysis, type = "sites",
#' nbrPermutations = 3, vSeed = 221)}
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom BiocParallel bplapply MulticoreParam SnowParam bptry bpok
#' @importFrom methods new
#' @export
runPermutationUsingMethylKitInfo <- function(methylKitInfo,
                            type = c("both", "sites", "tiles"),
                            outputDir = NULL,
                            runObservedAnalysis = TRUE,
                            nbrPermutations = 1000,
                            nbrCores = 1,
                            nbrCoresDiffMeth = 1,
                            minReads = 10,
                            minMethDiff = 10,
                            qvalue = 0.01,
                            maxPercReads = 99.9,
                            destrand = FALSE,
                            minCovBasesForTiles = 0,
                            tileSize = 1000,
                            stepSize = 1000,
                            vSeed = -1) {

    ## Parameters validation
    validateRunPermutationUsingMethylKitInfo(methylKitInfo = methylKitInfo,
                                type = type, outputDir = outputDir,
                                runObservedAnalysis = runObservedAnalysis,
                                nbrPermutations = nbrPermutations,
                                nbrCores = nbrCores,
                                nbrCoresDiffMeth = nbrCoresDiffMeth,
                                minReads = minReads, minMethDiff = minMethDiff,
                                qvalue = qvalue, maxPercReads = maxPercReads,
                                destrand = destrand,
                                minCovBasesForTiles = minCovBasesForTiles,
                                tileSize = tileSize, stepSize = stepSize,
                                vSeed = vSeed)

    ## Add last slash to path when absent
    if (!is.null(outputDir) &&
            (substr(outputDir, nchar(outputDir), nchar(outputDir)) != "/")) {
        outputDir <- paste0(outputDir, "/")
    }

    ## Set vSeed value when negative seed is given
    if (vSeed <= -1) {
        tSeed <- as.numeric(Sys.time())
        vSeed <- 1e8 * (tSeed - floor(tSeed))
    }
    set.seed(vSeed)

    ## Extract information
    nbGenerations <- length(methylKitInfo)
    nbSamplesByGeneration <- sapply(methylKitInfo, length)
    nbSamples  <- sum(nbSamplesByGeneration)
    allSamples <- unlist(methylKitInfo, recursive = FALSE)

    ## Create all permutations
    permutationSamples <- t(replicate(nbrPermutations, sample(1:nbSamples)))

    ## Create list that will contain all information to run permutation
    finalList <- vector("list", nbrPermutations)

    for (i in 1:nbrPermutations) {
        ## Create list that will contain information for all generations
        ## related to the same permutation analysis
        permutationList <- vector("list", nbGenerations)
        start <- 1
        for (j in 1:nbGenerations) {
            end <- start + nbSamplesByGeneration[j] - 1
            samplePos <- permutationSamples[i, start:end]
            treatment <- methylKitInfo[[j]]@treatment
            newSampleList <- new("methylRawList", allSamples[samplePos],
                                    treatment = treatment)
            permutationList[[j]] <- newSampleList
            start <- end + 1
        }

        finalList[[i]] <- list(sample = permutationList, id = i)
    }

    rm(permutationSamples)

    if (nbrCores == 1) {
        bpParam <- SnowParam()
    } else {
        bpParam <- MulticoreParam(workers = nbrCores)
    }

    redoList <- list()

    if (!is.null(outputDir)) {
        doTiles <- any(type %in% c("tiles", "both"))
        doSites <- any(type %in% c("sites", "both"))
        createOutputDir(outputDir, doingSites = doSites, doingTiles = doTiles)
    }

    ## Call observation analysis
    observedResults <- runObservationUsingMethylKitInfo(methylKitInfo =
                                                            methylKitInfo,
                                    type = type,
                                    outputDir = outputDir,
                                    nbrPermutations = nbrPermutations,
                                    nbrCores = nbrCores,
                                    nbrCoresDiffMeth = nbrCoresDiffMeth,
                                    minReads = minReads,
                                    minMethDiff = minMethDiff,
                                    qvalue = qvalue,
                                    maxPercReads = maxPercReads,
                                    destrand = destrand,
                                    minCovBasesForTiles = minCovBasesForTiles,
                                    tileSize = tileSize,
                                    stepSize = stepSize,
                                    vSeed = vSeed)

    ## Call permutations in parallel mode
    permutationResults <- bplapply(finalList, FUN =
                                        runOnePermutationOnAllGenerations,
                            type = type,
                            outputDir = outputDir,
                            nbrCoresDiffMeth = nbrCoresDiffMeth,
                            minReads = minReads,
                            minMethDiff = minMethDiff,
                            qvalue = qvalue,
                            maxPercReads = maxPercReads,
                            destrand = destrand,
                            minCovBasesForTiles = minCovBasesForTiles,
                            tileSize = tileSize,
                            stepSize = stepSize,
                        BPREDO = redoList,
                        BPPARAM = bpParam)

    ## Create final returned list
    result <- list()
    result[["OBSERVATION"]] <- observedResults
    result[["PERMUTATION"]] <- permutationResults

    return(result)
}


#' @title Run a differentially methylation analysis on each generation
#' present in a dataset.
#'
#' @description Run a differentially methylation analysis on each generation
#' present in a dataset. The number of conserved differentially
#' methylated elements (sites, tile or both) is them calculated. The
#' methylKit package is used to identify the differentially methylated
#' elements.
#'
#' @param methylKitInfo a \code{list} of \code{methylRawList} entries, each
#' \code{methylRawList} contains all the \code{methylRaw} entries related to
#' one generation (first entry = first generation, second entry = second
#' generation, etc..). The number of generations must correspond to the number
#' of entries in the \code{methylKitInfo}.At least 2 generations
#' must be present to calculate the conserved elements. More information can
#' be found in the methylKit package.
#'
#' @param type One of the "sites","tiles" or "both" strings. Specifies the type
#' of differentially methylated elements should be returned. For
#' retrieving differentially methylated bases type="sites"; for
#' differentially methylated regions type="tiles". Default: "both".
#'
#' @param outputDir a string, the name of the directory that will contain
#' the results of the permutation or \code{NULL}. If the directory does not
#' exist, it will be created. When \code{NULL}, the results of the permutation
#' are not saved. Default: \code{NULL}.
#'
#' @param nbrCores a positive \code{integer}, the number of cores to use when
#' processing the analysis. Default: \code{1} and always \code{1} for Windows.
#'
#' @param nbrPermutations, a positive \code{integer}, the total number of
#' permutations that is going to be done. Default: \code{1000}.
#'
#' @param nbrCoresDiffMeth a positive \code{integer}, the number of cores
#' to use for parallel differential methylation calculations.Parameter
#' used for both sites and tiles analysis. The parameter
#' corresponds to the \code{num.cores} parameter in the package
#' \code{methylKit}.
#' Default: \code{1} and always \code{1} for Windows.
#'
#' @param minReads a positive \code{integer} Bases and regions having lower
#' coverage than this count are discarded. The parameter
#' correspond to the \code{lo.count} parameter in the package \code{methylKit}.
#'
#' @param minMethDiff a positive \code{double} betwwen [0,100], the absolute
#' value of methylation percentage change between cases and controls. The
#' parameter correspond to the \code{difference} parameter in
#' the methylKit package. Default: \code{10}.
#'
#' @param qvalue a positive \code{double} betwwen [0,1], the cutoff
#' for qvalue of differential methylation statistic. Default: \code{0.01}.
#'
#' @param maxPercReads a \code{double} between [0,100], the percentile of read
#' counts that is going to be used as upper cutoff. Bases ore regions
#' having higher
#' coverage than this percentile are discarded. Parameter used for both CpG
#' sites and tiles analysis. The parameter
#' correspond to the \code{hi.perc} parameter in the package \code{methylKit}.
#' Default: \code{99.9}.
#'
#' @param destrand a \code{logical}, when \code{TRUE} will merge reads on both
#' strands of a CpG dinucleotide to provide better coverage. Only advised
#' when looking at CpG methylation. Parameter used for both CpG
#' sites and tiles analysis.
#' Default: \code{FALSE}.
#'
#' @param minCovBasesForTiles a non-negative \code{integer}, the minimum
#' number of bases to be covered in a given tiling window. The parameter
#' corresponds to the \code{cov.bases} parameter in the package
#' \code{methylKit}.
#' Only used when \code{doingTiles} =
#' \code{TRUE}. Default: \code{0}.
#'
#' @param tileSize a positive \code{integer}, the size of the tiling window.
#' The parameter corresponds to the \code{win.size} parameter in
#' the package \code{methylKit}. Only
#' used when \code{doingTiles} = \code{TRUE}. Default: \code{1000}.
#'
#' @param stepSize a positive \code{integer}, the step size of tiling windows.
#' The parameter corresponds to the \code{stepSize} parameter in
#' the package \code{methylKit}. Only
#' used when \code{doingTiles} = \code{TRUE}. Default: \code{1000}.
#'
#' @param vSeed a \code{integer}, a seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used. Default: \code{-1}.
#'
#' @return a \code{list} containing the following elements:
#' \itemize{
#' \item \code{SITES} Only present when \code{type} = \code{"sites"} or
#' \code{both}, a \code{list} containing:
#' \itemize{
#' \item\code{i2} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer}, the number of conserved
#' hyper differentially methylated sites between two consecutive generations.
#' The first element represents the intersection of the first and second
#' generations; the second element, the intersection of the second and third
#' generations; etc..
#' \item \code{HYPO} a \code{list} of \code{integer}, the number of conserved
#' hypo differentially methylated sites between two consecutive generations.The
#' first element represents the intersection of the first and second
#' generations; the second element, the intersection of the second and third
#' generations; etc..
#' }
#' \item\code{iAll} a \code{list} containing:
#' \itemize{
#'\item \code{HYPER} a \code{list} of \code{integer}, the number of conserved
#' hyper differentially methylated sites between three or more consecutive
#' generations. The first element represents the intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc..The number of entries depends of the number
#' of generations.
#' \item \code{HYPO} a \code{list} of \code{integer}, the number of conserved
#' hypo differentially methylated sites between three or more consecutive
#' generations. The first element represents the intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc..The number of entries depends of the number of
#' generations.
#' }
#' }
#' \item \code{TILES} Only present when \code{type} = \code{"tiles"} or
#' \code{both}, a \code{list} containing:
#' itemize{
#' \item\code{i2} a \code{list} containing:
#' \itemize{
#' \item \code{HYPER} a \code{list} of \code{integer}, the number of conserved
#' hyper differentially methylated positions between two consecutive
#' generations. The first element represents the intersection of the
#' first and second generations; the second element, the intersection of
#' the second and third generations; etc..
#' \item \code{HYPO} a \code{list} of \code{integer}, the number of conserved
#' hypo differentially methylated positions between two consecutive
#' generations.The first element represents the intersection of the first and
#' second generations; the second element, the intersection of the second
#' and third generations; etc..
#' }
#' \item\code{iAll} a \code{list} containing:
#' \itemize{
#'\item \code{HYPER} a \code{list} of \code{integer}, the number of conserved
#' hyper differentially methylated positions between three or more consecutive
#' generations. The first element represents the intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc..The number of entries depends of the number
#' of generations.
#' \item \code{HYPO} a \code{list} of \code{integer}, the number of conserved
#' hypo differentially methylated positions between three or more consecutive
#' generations. The first element represents the intersection of the first
#' three generations; the second element, the intersection of the first fourth
#' generations; etc..The number of entries depends of the number of
#' generations.
#' }
#' }
#' }
#'
#' @examples
#'
#' ## Load methyl information
#' data(samplesForTransgenerationalAnalysis)
#'
#' ## Run a permutation analysis
#' \dontrun{runObservationUsingMethylKitInfo(methylKitInfo =
#' samplesForTransgenerationalAnalysis, type = "sites",
#' nbrPermutations = 0, vSeed = 221)}
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @export
runObservationUsingMethylKitInfo <- function(methylKitInfo,
                                            type = c("both", "sites", "tiles"),
                                            outputDir = NULL,
                                            nbrPermutations = 0,
                                            nbrCores = 1,
                                            nbrCoresDiffMeth = 1,
                                            minReads = 10,
                                            minMethDiff = 10,
                                            qvalue = 0.01,
                                            maxPercReads = 99.9,
                                            destrand = FALSE,
                                            minCovBasesForTiles = 0,
                                            tileSize = 1000,
                                            stepSize = 1000,
                                            vSeed = -1) {

    ## Parameters validation
    validateRunObservationUsingMethylKitInfo(methylKitInfo = methylKitInfo,
                            type = type, outputDir = outputDir,
                            nbrCores = nbrCores,
                            nbrCoresDiffMeth = nbrCoresDiffMeth,
                            minReads = minReads, minMethDiff = minMethDiff,
                            qvalue = qvalue,
                            maxPercReads = maxPercReads, destrand = destrand,
                            minCovBasesForTiles = minCovBasesForTiles,
                            tileSize = tileSize,
                            stepSize = stepSize, vSeed = vSeed)

    ## Add last slash to path when absent
    if (!is.null(outputDir) &&
        (substr(outputDir, nchar(outputDir), nchar(outputDir)) != "/")) {
        outputDir <- paste0(outputDir, "/")
    }

    ## Set vSeed value when negative seed is given
    if (vSeed <= -1) {
        tSeed <- as.numeric(Sys.time())
        vSeed <- 1e8 * (tSeed - floor(tSeed))
    }
    set.seed(vSeed)

    methylInfo <- list(sample = methylKitInfo, id = 0)

    if (!is.null(outputDir)) {
        doTiles <- any(type %in% c("tiles", "both"))
        doSites <- any(type %in% c("sites", "both"))
        createOutputDir(outputDir, doingSites = doSites, doingTiles = doTiles)
    }

    ## Extract information
    result <- runOnePermutationOnAllGenerations(methylInfoForAllGenerations =
                                                        methylInfo,
                                    type = type, outputDir = outputDir,
                                    nbrCoresDiffMeth = nbrCoresDiffMeth,
                                    minReads = minReads,
                                    minMethDiff = minMethDiff,
                                    qvalue = qvalue,
                                    maxPercReads = maxPercReads,
                                    destrand = destrand,
                                    minCovBasesForTiles = minCovBasesForTiles,
                                    tileSize = tileSize,
                                    stepSize = stepSize)

    return(result)
}

#' @title Run a differentially methylation analysis on each generation
#' present in a dataset.
#'
#' @description Run a differentially methylation analysis on each generation
#' present in a dataset. The number of conserved differentially
#' methylated elements (sites, tile or both) is them calculated. The
#' methylKit package is used to identify the differentially methylated
#' elements.
#'
#' @param methylKitRDSFile a string, the name of the RDS file containing the
#' methylKit objet used for the permutation analysis. The RDS file must
#' hold a \code{list} of \code{methylRawList} entries, each
#' \code{methylRawList} contains all the \code{methylRaw} entries related to
#' one generation (first entry = first generation, second entry = second
#' generation, etc..). The number of generations must correspond to the number
#' of entries in the \code{methylKitInfo}.At least 2 generations
#' must be present to do a permutation analysis. More information can be found
#' in the methylKit package.
#'
#' @param type One of the "sites","tiles" or "both" strings. Specifies the type
#' of differentially methylated elements should be returned. For
#' retrieving differentially methylated bases type="sites"; for
#' differentially methylated regions type="tiles". Default: "both".
#'
#' @param outputDir a string, the name of the directory that will contain
#' the results of the permutation or \code{NULL}. If the directory does not
#' exist, it will be created. When \code{NULL}, the results of the permutation
#' are not saved. Default: \code{NULL}.
#'
#' @param nbrCores a positive \code{integer}, the number of cores to use when
#' processing the analysis. Default: \code{1} and always \code{1} for Windows.
#'
#' @param nbrPermutations, a positive \code{integer}, the total number of
#' permutations that is going to be done. Default: \code{1000}.
#'
#' @param nbrCoresDiffMeth a positive \code{integer}, the number of cores
#' to use for parallel differential methylation calculations.Parameter
#' used for both sites and tiles analysis. The parameter
#' corresponds to the \code{num.cores} parameter in the package
#' \code{methylKit}.
#' Default: \code{1} and always \code{1} for Windows.
#'
#' @param minReads a positive \code{integer} Bases and regions having lower
#' coverage than this count are discarded. The parameter
#' correspond to the \code{lo.count} parameter in the package \code{methylKit}.
#'
#' @param minMethDiff a positive \code{double} betwwen [0,100], the absolute
#' value of methylation percentage change between cases and controls. The
#' parameter correspond to the \code{difference} parameter in
#' the methylKit package. Default: \code{10}.
#'
#' @param qvalue a positive \code{double} betwwen [0,1], the cutoff
#' for qvalue of differential methylation statistic. Default: \code{0.01}.
#'
#' @param maxPercReads a \code{double} between [0,100], the percentile of read
#' counts that is going to be used as upper cutoff. Bases ore regions
#' having higher
#' coverage than this percentile are discarded. Parameter used for both CpG
#' sites and tiles analysis. The parameter
#' correspond to the \code{hi.perc} parameter in the package \code{methylKit}.
#' Default: \code{99.9}.
#'
#' @param destrand a \code{logical}, when \code{TRUE} will merge reads on both
#' strands of a CpG dinucleotide to provide better coverage. Only advised
#' when looking at CpG methylation. Parameter used for both CpG
#' sites and tiles analysis.
#' Default: \code{FALSE}.
#'
#' @param minCovBasesForTiles a non-negative \code{integer}, the minimum
#' number of bases to be covered in a given tiling window. The parameter
#' corresponds to the \code{cov.bases} parameter in the package
#' \code{methylKit}.
#' Only used when \code{doingTiles} =
#' \code{TRUE}. Default: \code{0}.
#'
#' @param tileSize a positive \code{integer}, the size of the tiling window.
#' The parameter corresponds to the \code{win.size} parameter in
#' the package \code{methylKit}. Only
#' used when \code{doingTiles} = \code{TRUE}. Default: \code{1000}.
#'
#' @param stepSize a positive \code{integer}, the step size of tiling windows.
#' The parameter corresponds to the \code{stepSize} parameter in
#' the package \code{methylKit}. Only
#' used when \code{doingTiles} = \code{TRUE}. Default: \code{1000}.
#'
#' @param vSeed a \code{integer}, a seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used. Default: \code{-1}.
#'
#' @return TODO
#'
#' @examples
#'
#' ## Path to a methylKit RDS file
#' methylFile <- dir(system.file("extdata", package = "methylInheritance"),
#' pattern = "methylObj_002.RDS", full.names = TRUE)
#'
#' ## Run a permutation analysis
#' \dontrun{runObservationUsingRDSFile(methylKitRDSFile = methylFile,
#' type = "sites", nbrPermutations = 0, minReads = 8, minMethDiff = 5,
#' vSeed = 2001)}
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @export
runObservationUsingRDSFile <- function(methylKitRDSFile,
                                            type = c("both", "sites", "tiles"),
                                            outputDir = NULL,
                                            nbrPermutations = 0,
                                            nbrCores = 1,
                                            nbrCoresDiffMeth = 1,
                                            minReads = 10,
                                            minMethDiff = 10,
                                            qvalue = 0.01,
                                            maxPercReads = 99.9,
                                            destrand = FALSE,
                                            minCovBasesForTiles = 0,
                                            tileSize = 1000,
                                            stepSize = 1000,
                                            vSeed = -1) {

    ## Validate that methylKitRDSFile is an existing file
    if (!file.exists(methylKitRDSFile)) {
        stop(paste0("The file \"", methylKitRDSFile, "\" does not exist."))
    }

    ## Extract information from RDS file
    methylInfo <- readRDS(methylKitRDSFile)

    ## Call permutation analysis with the methylInfo object as an parameter
    result <- runObservationUsingMethylKitInfo(methylKitInfo = methylInfo,
                                    type = type,
                                    outputDir = outputDir,
                                    nbrPermutations = nbrPermutations,
                                    nbrCores = nbrCores,
                                    nbrCoresDiffMeth = nbrCoresDiffMeth ,
                                    minReads = minReads,
                                    minMethDiff = minMethDiff,
                                    qvalue = qvalue,
                                    maxPercReads = maxPercReads,
                                    destrand = destrand,
                                    minCovBasesForTiles = minCovBasesForTiles,
                                    tileSize = tileSize, stepSize = stepSize,
                                    vSeed = vSeed)

    return(result)
}


#' @title Load all RDS files created by the permutation  and observation
#' analysis
#'
#' @description  Load all RDS files created by the permutation and
#' observation analysis. The function
#' returns an object of \code{class} "methylInheritanceAllResults" that holds
#' all the pertinent information.
#'
#' @param analysisResultsDir a \code{character} string, the path to the
#' directory that contains the analysis results. The path can be the same that
#' for the \code{permutatioNResultsDir} parameter.
#'
#' @param permutationResultsDir a \code{character} string, the path to the
#' directory that contains the permutation results. The path can be the same
#' that for the \code{analysisResultsDir} parameter.
#'
#' @param doingSites a \code{logical}, the data related to differentially
#' methylated sites are loaded when
#' \code{doingSites} = \code{TRUE}. Default: \code{TRUE}.
#'
#' @param doingTiles a \code{logical}, the data related to differentially
#' methylated tiles are loaded when
#' \code{doingTiles} = \code{TRUE}. Default: \code{TRUE}.
#'
#' @return an \code{list} of \code{class} "methylInheritanceAllResults"
#' containing the following elements:
#' \itemize{
#'     \item \code{OBSERVATION}, a \code{list} that contains one or two
#'     \code{list} depending
#'     of the specified parameters. When \code{doingSites} = \code{TRUE}, a
#'     list of called "SITES" is present. When \code{doingTiles} = \code{TRUE},
#'     a list of called "TILES" is present. TODO
#'     \item \code{PERMUTATION}, a \code{list}  contains one or two \code{list}
#'     depending of the specified parameters. When
#'     \code{doingSites} = \code{TRUE}, a
#'     list of called "SITES" is present. When \code{doingTiles} = \code{TRUE},
#'     a list of called "TILES" is present. TODO
#' }
#'
#' @examples
#'
#' ## Get the name of the directory where files are stored
#' filesDir <- dir(system.file("extdata", package = "methylInheritance"),
#' pattern = "TEST", full.names = TRUE)
#'
#' ## Load information from files
#' results <- loadAllRDSResults(analysisResultsDir = filesDir,
#' permutationResultsDir = filesDir, doingSites = TRUE, doingTiles = TRUE)
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @export
loadAllRDSResults <- function(analysisResultsDir,
                                    permutationResultsDir,
                                    doingSites = TRUE, doingTiles = FALSE) {

    ## Add last slash to analysisResultsDIR when absent
    if (!is.null(analysisResultsDir) &&
        (substr(analysisResultsDir, nchar(analysisResultsDir),
                    nchar(analysisResultsDir)) != "/")) {
        analysisResultsDir <- paste0(analysisResultsDir, "/")
    }

    ## Add last slash to permutationResultsDIR when absent
    if (!is.null(permutationResultsDir) &&
        (substr(permutationResultsDir, nchar(permutationResultsDir),
                    nchar(permutationResultsDir)) != "/")) {
        permutationResultsDir <- paste0(permutationResultsDir, "/")
    }

    result<-list()

    ## SITES
    if (doingSites) {
        analysisResults <- readRDS(file = paste0(analysisResultsDir,
                                        "SITES/SITES_observed_results.RDS"))
        analysisStruct <- createDataStructure(interGenerationResult =
                                                    analysisResults)
        result[["OBSERVATION"]][["SITES"]] <- analysisStruct


        filesInDir <- list.files(path = paste0(analysisResultsDir,
                                                                "SITES/"),
                                pattern = "[[:digit:]].RDS", all.files = FALSE,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE,
                                no.. = FALSE)

        sitesPerm <- lapply(filesInDir, FUN = function(x) {readRDS(file = x)})

        t <- lapply(sitesPerm, FUN = function(x) {
                    createDataStructure(interGenerationResult = x)})

        result[["PERMUTATION"]][["SITES"]] <- t
    }

    ## TILES
    if (doingTiles) {
        analysisResults <- readRDS(file = paste0(permutationResultsDir,
                                        "TILES/TILES_observed_results.RDS"))
        analysisStruct <- createDataStructure(interGenerationResult =
                                                    analysisResults)
        result[["OBSERVATION"]][["TILES"]] <- analysisStruct

        filesInDir <- list.files(path = paste0(permutationResultsDir,
                                                            "TILES/"),
                                pattern = "[[:digit:]].RDS", all.files = FALSE,
                                full.names = TRUE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE,
                                no.. = FALSE)

        tilesPerm <- lapply(filesInDir, FUN = function(x) {readRDS(file = x)})

        t <- lapply(tilesPerm, FUN = function(x) {
                    createDataStructure(interGenerationResult = x)})

        result[["PERMUTATION"]][["TILES"]] <- t
    }

    class(result)<-"methylInheritanceAllResults"

    return(result)
}


#' @title TODO
#'
#' @description  TODO
#'
#' @param analysisandPermutationResults TODO
#'
#' @param type One of the "sites" or "tiles" strings. Specifies the type
#' of differentially methylated elements should be returned. For
#' retrieving differentially methylated bases type="sites"; for
#' differentially methylated regions type="tiles". Default: "sites".
#'
#' @param inter TODO
#'
#' @param position TODO
#'
#' @return TODO
#'
#' @examples
#'
#' ## Loading dataset containing all results
#' data(analysisResults)
#'
#' ## Extract information for the intersection between conserved differentially
#' ## methylated sites (type = sites) between the intersection of 2
#' ## generations (inter = i2): F2 and F3 (position = 2)
#' formatForGraph(analysisandPermutationResults = analysisResults,
#' type = "sites", inter="i2", 2)
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @export
formatForGraph <- function(analysisandPermutationResults,
                            type = c("sites", "tiles"),
                            inter=c("i2", "iAll"), position) {

    type <- toupper(type)

    real <- analysisandPermutationResults[["OBSERVED"]][[type]][[inter]]

    dataConserved <- data.frame(TYPE=c("HYPO", "HYPER"),
                                result=c(real[["HYPO"]][[position]],
                                            real[["HYPER"]][[position]]),
                                SOURCE=c("OBSERVED", "OBSERVED"))

    for (i in 1:length(analysisandPermutationResults[["PERMUTATION"]][[type]])) {
        permutation <- analysisandPermutationResults[["PERMUTATION"]][[type]][[i]][[inter]]
        dataConserved <- rbind(dataConserved, data.frame(TYPE=c("HYPO", "HYPER"),
                        result=c(permutation[["HYPO"]][[position]],
                                    permutation[["HYPER"]][[position]]),
                        SOURCE=c("PERMUTATION", "PERMUTATION")))
    }

    return(dataConserved)

}


#' @title Generate a graph for a permutation analysis
#'
#' @description  TODO
#'
#' @param formatForGraphDataFrame TODO
#'
#' @return TODO
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom ggplot2 ggplot geom_text facet_grid theme geom_vline geom_histogram labs aes
#' @export
plotGraph <- function(formatForGraphDataFrame) {

    # Basic graph using dataframe
    # Columns names : TYPE (HYPER or HYPO), result (nbr conseved sites),
    # SOURCE (OBSERVED or PERMUTATION)
    p <- ggplot(data=formatForGraphDataFrame, aes(formatForGraphDataFrame$result)) +
        geom_histogram(col="blue",
                        fill="lightblue",
                        binwidth=2,
                        alpha = .2) +
        labs(title = "") +
        labs(x = "Number of conserved differentially methylated sites",
                y = "Frequency")

    # Split to have one section for HYPER and one for HYPO
    p <- p + facet_grid(.~TYPE)

    # Add vertical line corresponding to the number of conserved elements
    # in the observed results (real results)
    interceptFrame <- subset(formatForGraphDataFrame,
                                formatForGraphDataFrame$SOURCE == "OBSERVED")
    p <- p + geom_vline(aes(xintercept = interceptFrame$result, color="red"),
                        data = interceptFrame, linetype="longdash",
                        show.legend = NA)

    # Calculate the significant level for HYPER AND HYPO
    hypoDataSet <- subset(formatForGraphDataFrame,
                            formatForGraphDataFrame$TYPE == "HYPO")
    hypoTotal <- nrow(hypoDataSet)
    hypoNumber <- interceptFrame[interceptFrame$TYPE == "HYPO",]$result
    signifLevelHypo <- nrow(subset(hypoDataSet,
                                hypoDataSet$result >= hypoNumber))/hypoTotal

    hyperDataSet <- subset(formatForGraphDataFrame,
                            formatForGraphDataFrame$TYPE == "HYPER")
    hyperTotal <- nrow(hyperDataSet)
    hyperNumber <- interceptFrame[interceptFrame$TYPE == "HYPER",]$result
    signifLevelHyper <- nrow(subset(hyperDataSet,
                                hyperDataSet$result >= hyperNumber))/hyperTotal

    # Add significant level value as annotated text
    ann_text <- data.frame(TYPE = c("HYPO", "HYPER"),
                                lab = c(sprintf("%.6f", signifLevelHypo),
                                        sprintf("%.6f", signifLevelHyper)))
    p <- p + geom_text(data = ann_text,
                        aes(x=100, y=100, label=ann_text$lab, size=15,
                                color="red"),
                        inherit.aes=FALSE, parse=FALSE)

    # Add number of observed conserved elements as annotated text
    ann_text_2 <- data.frame(TYPE = c("HYPO", "HYPER"),
                           lab = c(sprintf("%d observed", hypoNumber),
                                    sprintf("%d observed", hyperNumber)))
    p <- p + geom_text(data = ann_text_2,  aes(x = 100, y = 300,
                                                label=ann_text_2$lab, size=15,
                                                color="red"),
                        inherit.aes=FALSE, parse=FALSE)

    p <- p + theme(legend.position="none")

    return(p)

}
