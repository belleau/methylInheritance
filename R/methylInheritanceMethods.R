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
#' pattern = "methylObj_001.RDS", full.names = TRUE)
#'
#' ## Run a permutation analysis
#' \dontrun{runPermutationUsingRDS(methylKitRDSFile = methylFile,
#' type = "sites", nbrPermutations = 10, vSeed = 2001)}
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @export
runPermutationUsingRDSFile <- function(methylKitRDSFile,
                                   type = c("both", "sites", "tiles"),
                                   outputDir = NULL,
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
    validateRunPermutationUsingMethylKitInfo(methylKitInfo, type,
                                outputDir, nbrPermutations,
                                nbrCores, nbrCoresDiffMeth,
                                minReads, minMethDiff, qvalue, maxPercReads,
                                destrand, minCovBasesForTiles, tileSize,
                                stepSize, vSeed)

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

    ## Call permutations in parallel mode
    result <- bplapply(finalList, FUN = runOnePermutationOnAllGenerations,
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

    return(result)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param methylKitRDSFile a string, the name of the file containing the
#' methylKit objet used for the permutation.
#' TODO
#' of files with methylation information for
#' bases or regions in the genome. One \code{list} must contain all files
#' related to the same generation. So, if 3 generations are analyzed, a
#' \code{list} containing 3 \code{list} must be passed. At least 2 generations
#' must be present to do a permutation analysis.
#' The parameter
#' corresponds to the \code{location} parameter in the package \code{methylKit}.
#' A \code{list} of \code{vector} containing
#' \code{0} and \code{1}. The information indicating which files are
#' associated to controls (\code{0}) and which files are cases (\code{1}).
#' One \code{vector} must contain all informations
#' related to the same generation. So, if 3 generations are analyzed, a
#' \code{list} containing 3 \code{vector} must be passed. At least 2
#' generations
#' must be present to do a permutation analysis.
#' The parameter
#' corresponds to the \code{treatment} parameter in the package
#' \code{methylKit}.
#'
#' @param type One of the "sites","tiles" or "both" strings. Specifies the type
#' of differentially methylated elements should be returned. For
#' retrieving differentially methylated bases type="sites"; for
#' differentially methylated regions type="tiles". Default: "both".
#'
#' @param outputDir a string, the name of the directory that will contain
#' the results of the permutation or \code{NULL}. If the directory does not
#' exist, it will be created. When \code{NULL}, the results of the permutation
#' are not saved.
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
#' ##TODO
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom methods new
#' @export
runPermutationUsingRDSTEST <- function(methylKitRDSFile,
                                   type = c("both", "sites", "tiles"),
                                   outputDir = NULL,
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


    ## Extract information from RDS file
    methylInfo <- readRDS(methylKitRDSFile)

    ## Parameters validation
    validateRunPermutationUsingMethylKitInfo(methylInfo, type,
                                   outputDir, nbrPermutations,
                                   nbrCores, nbrCoresDiffMeth,
                                   minReads, minMethDiff, qvalue, maxPercReads,
                                   destrand, minCovBasesForTiles, tileSize,
                                   stepSize, vSeed)

    ## Add last slash to path when absent
    if (substr(outputDir, nchar(outputDir), nchar(outputDir)) != "/") {
        outputDir <- paste0(outputDir, "/")
    }

    ## Set vSeed value when negative seed is given
    if (vSeed <= -1) {
        tSeed <- as.numeric(Sys.time())
        vSeed <- 1e8 * (tSeed - floor(tSeed))
    }
    set.seed(vSeed)

    nbGenerations <- length(methylInfo)

    nbSamplesByGeneration <- sapply(methylInfo, length)
    nbSamples  <- sum(nbSamplesByGeneration)
    allSamples <- unlist(methylInfo, recursive = FALSE)

    permutationSamples <- t(replicate(nbrPermutations, sample(1:nbSamples)))

    finalList <- vector("list", nbrPermutations)

    for (i in 1:nbrPermutations) {

        # Randomnly mixt all samples
        #samples <- sample(allSamples, size=nbSamples, replace = FALSE)

        permutationList <- vector("list", nbGenerations)
        start = 1
        for (j in 1:nbGenerations) {
            end = start + nbSamplesByGeneration[j] - 1
            samplePos <- permutationSamples[i, start:end]
            newSampleList <- new("methylRawList", allSamples[samplePos],
                                 treatment = methylInfo[[1]]@treatment)
            permutationList[[j]] <- newSampleList
            start <- end + 1
        }
        finalList[[i]] <- permutationList
    }

    rm(permutationSamples)

    result <- runOnePermutationOnAllGenerations(finalList[[1]],
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
                                                stepSize = stepSize)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param realAnalysis_output_dir a string, the name of the directory that
#' will contain the results of the permutation. If the directory does not
#' exist, it will be created.
#'
#' @param permutations_output_dir a string, the name of the directory that
#' will contain the results of the permutation. If the directory does not
#' exist, it will be created.
#'
#' @param doingSites a \code{logical}, TODO
#'
#' @param doingTiles a \code{logical}, TODO
#'
#' @return TODO
#'
#' @examples
#'
#' ##TODO
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @export
extractData <- function(realAnalysis_output_dir, permutations_output_dir,
                       doingSites = TRUE, doingTiles=FALSE) {

    ## Add last slash to path when absent
    if (substr(realAnalysis_output_dir, nchar(realAnalysis_output_dir),
                nchar(realAnalysis_output_dir)) != "/") {
        realAnalysis_output_dir <- paste0(realAnalysis_output_dir, "/")
    }

    ## Add last slash to path when absent
    if (substr(permutations_output_dir, nchar(permutations_output_dir),
                nchar(permutations_output_dir)) != "/") {
        permutations_output_dir <- paste0(permutations_output_dir, "/")
    }

    nbr_sites_per_generation <- list()

    ## SITES
    if (doingSites) {
        result <- extractDataFromFile(permutations_output_dir, "SITES")
        nbr_sites_per_generation[["SITES"]] <- result[["SITES"]]
    }

    ## TILES
    if (doingTiles) {
        result <- extractDataFromFile(permutations_output_dir, "TILES")
        nbr_sites_per_generation[["TILES"]] <- result[["TILES"]]
    }

    return(nbr_sites_per_generation)
}
