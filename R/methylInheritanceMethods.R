#' @title TODO
#'
#' @description TODO
#'
#' @param allFilesByGeneration a \code{list} of \code{list}, a \code{list}
#' composed of \code{list}
#' of files with methylation information for
#' bases or regions in the genome. One \code{list} must contain all files
#' related to the same generation. So, if 3 generations are analyzed, a
#' \code{list} containing 3 \code{list} must be passed. At least 2 generations
#' must be present to do a permutation analysis.
#' The parameter
#' corresponds to the \code{location} parameter in the package \code{methylKit}.
#'
#' @param conditionsByGeneration a \code{list} of \code{vector} containing
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
#' @param output_dir a string, the name of the directory that will contain
#' the results of the permutation. If the directory does not exist, it will
#' be created.
#'
#' @param genomeVersion a string, the genome assembly such as hg18, mm9, etc.
#' It can be any string. The parameter
#' correspond to the \code{assembly} parameter in the package \code{methylKit}.
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
#' @param doingSites a \code{logical}, when \code{TRUE} will do the analysis
#' on the CpG dinucleotide sites. Default: \code{TRUE}.
#'
#' @param doingTiles a \code{logical}, when \code{TRUE} will do the analysis
#' on the tiles. Default: \code{FALSE}.
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
#' @importFrom BiocParallel bplapply MulticoreParam SnowParam bptry bpok bpmapply
#' @importFrom utils flush.console write.table
#' @export
runPermutation <- function(allFilesByGeneration, conditionsByGeneration,
                            output_dir,
                            genomeVersion,
                            nbrPermutations = 1000,
                            nbrCores = 1,
                            nbrCoresDiffMeth = 1,
                            doingSites = TRUE,
                            doingTiles = FALSE,
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
    validateRunPermutation(allFilesByGeneration, conditionsByGeneration,
                            output_dir, genomeVersion, nbrPermutations,
                            nbrCores, nbrCoresDiffMeth, doingSites, doingTiles,
                            minReads, minMethDiff, qvalue, maxPercReads,
                            destrand, minCovBasesForTiles, tileSize,
                            stepSize, vSeed)

    ## Add last slash to path when absent
    if (substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") {
        output_dir <- paste0(output_dir, "/")
    }

    ## Set vSeed value when negative seed is given
    if (vSeed <= -1) {
        tSeed <- as.numeric(Sys.time())
        vSeed <- 1e8 * (tSeed - floor(tSeed))
    }
    set.seed(vSeed)

    nbrFiles <- sum(unlist(lapply(allFilesByGeneration, length)))

    nbGenerations <- length(allFilesByGeneration)

    nbSamples  <- sum(length(unlist(allFilesByGeneration)))
    allSamples <- unlist(allFilesByGeneration)

    nbSamplesByGeneration <- sapply(allFilesByGeneration, length)

    finalList <- list()

    for (i in 1:nbrPermutations) {

        # Randomnly mixt all samples
        samples <- sample(allSamples, size=nbSamples, replace = FALSE)

        permutationList <- list()
        start = 1
        for (j in 1:nbGenerations) {
            nbrFilesGeneration <-  nbSamplesByGeneration[j]
            end = start + nbrFilesGeneration - 1
            filesGeneration <- samples[start:end]
            sampleNames <- sapply(filesGeneration, getSampleNameFromFileName)
            infoGeneration <- list(conditions=conditionsByGeneration[[j]],
                                   file.list=as.list(filesGeneration),
                                   designName = paste0("Generation_", j),
                                   sampleNames = sampleNames,
                                   count = i)
            finalList <- append(finalList, list(infoGeneration))
            start <- end + 1
        }
    }

    # Create directories for output files
    createOutputDir(output_dir, nbGenerations, doingSites, doingTiles)

    if (nbrCores == 1) {
        bp_param <- SnowParam()
    } else {
        bp_param <- MulticoreParam(workers = nbrCores)
    }

    result <- bpmapply(FUN = runOnePermutation,
                            finalList,
                            MoreArgs = list(
                            output_dir = output_dir,
                            genomeVersion = genomeVersion,
                            nbrCoresDiffMeth = nbrCoresDiffMeth,
                            doingSites = doingSites,
                            doingTiles = doingTiles,
                            minReads = minReads,
                            minMethDiff = minMethDiff,
                            qvalue = qvalue,
                            maxPercReads = maxPercReads,
                            destrand = destrand,
                            minCovBasesForTiles = minCovBasesForTiles,
                            tileSize = tileSize,
                            stepSize = stepSize),
                            BPPARAM = bp_param)

    result[which(!bpok(result))]
}




#' @title TODO
#'
#' @description TODO
#'
#' @param allFilesByGeneration a \code{list} of \code{list}, a \code{list}
#' composed of \code{list}
#' of files with methylation information for
#' bases or regions in the genome. One \code{list} must contain all files
#' related to the same generation. So, if 3 generations are analyzed, a
#' \code{list} containing 3 \code{list} must be passed. At least 2 generations
#' must be present to do a permutation analysis.
#' The parameter
#' corresponds to the \code{location} parameter in the package \code{methylKit}.
#'
#' @param conditionsByGeneration a \code{list} of \code{vector} containing
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
#' @param output_dir a string, the name of the directory that will contain
#' the results of the permutation. If the directory does not exist, it will
#' be created.
#'
#' @param genomeVersion a string, the genome assembly such as hg18, mm9, etc.
#' It can be any string. The parameter
#' correspond to the \code{assembly} parameter in the package \code{methylKit}.
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
#' @param doingSites a \code{logical}, when \code{TRUE} will do the analysis
#' on the CpG dinucleotide sites. Default: \code{TRUE}.
#'
#' @param doingTiles a \code{logical}, when \code{TRUE} will do the analysis
#' on the tiles. Default: \code{FALSE}.
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
#' @importFrom BiocParallel bplapply MulticoreParam SnowParam bptry bpok bpmapply
#' @importFrom utils flush.console write.table
#' @export
runRealAnalysis <- function(allFilesByGeneration, conditionsByGeneration,
                           output_dir,
                           genomeVersion,
                           nbrCores = 1,
                           nbrCoresDiffMeth = 1,
                           doingSites = TRUE,
                           doingTiles = FALSE,
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
    validateRunRealAnalysis(allFilesByGeneration, conditionsByGeneration,
                           output_dir, genomeVersion,
                           nbrCores, nbrCoresDiffMeth, doingSites, doingTiles,
                           minReads, minMethDiff, qvalue, maxPercReads,
                           destrand, minCovBasesForTiles, tileSize,
                           stepSize, vSeed)

    ## Add last slash to path when absent
    if (substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") {
        output_dir <- paste0(output_dir, "/")
    }

    ## Set vSeed value when negative seed is given
    if (vSeed <= -1) {
        tSeed <- as.numeric(Sys.time())
        vSeed <- 1e8 * (tSeed - floor(tSeed))
    }
    set.seed(vSeed)


    nbrFiles <- sum(unlist(lapply(allFilesByGeneration, length)))

    nbGenerations <- length(allFilesByGeneration)

    nbSamples  <- sum(length(unlist(allFilesByGeneration)))
    allSamples <- unlist(allFilesByGeneration)

    nbSamplesByGeneration <- sapply(allFilesByGeneration, length)

    # Create directories for output files
    createOutputDir(output_dir, nbGenerations, doingSites, doingTiles)

    finalList <- list()

    start = 1
    for (j in 1:nbGenerations) {
        nbrFilesGeneration <-  nbSamplesByGeneration[j]
        end = start + nbrFilesGeneration - 1
        filesGeneration <- allFilesByGeneration[[j]]
        sampleNames <- sapply(filesGeneration, getSampleNameFromFileName)
        infoGeneration <- list(conditions=conditionsByGeneration[[j]],
                               file.list=filesGeneration,
                               designName = paste0("Generation_", j),
                               sampleNames = sampleNames,
                               count = 0)
        finalList <- append(finalList, list(infoGeneration))
        start <- end + 1
    }


    if (nbrCores == 1) {
        bp_param <- SnowParam()
    } else {
        bp_param <- MulticoreParam(workers = nbrCores)
    }

    result <- bpmapply(FUN = runOnePermutation,
                       finalList,
                       MoreArgs = list(
                            output_dir = output_dir,
                            genomeVersion = genomeVersion,
                            nbrCoresDiffMeth = nbrCoresDiffMeth,
                            doingSites = doingSites,
                            doingTiles = doingTiles,
                            minReads = minReads,
                            minMethDiff = minMethDiff,
                            qvalue = qvalue,
                            maxPercReads = maxPercReads,
                            destrand = destrand,
                            minCovBasesForTiles = minCovBasesForTiles,
                            tileSize = tileSize,
                            stepSize = stepSize),
                            BPPARAM = bp_param)

    return(result)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param methylKitRDS a \code{list} of \code{list}, a \code{list}
#' composed of \code{list}
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
#' @param output_dir a string, the name of the directory that will contain
#' the results of the permutation. If the directory does not exist, it will
#' be created.
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
#' @param type One of the "sites","tiles" or "both" strings. Specifies the type
#' of differentially methylated elements should be returned. For
#' retrieving differentially methylated bases type="sites"; for
#' differentially methylated regions type="tiles". Default: "both".
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
#' @importFrom BiocParallel bplapply MulticoreParam SnowParam bptry bpok bpmapply
#' @importFrom utils flush.console write.table
#' @export
runPermutationUsingRDS <- function(methylKitRDS,
                           output_dir,
                           nbrPermutations = 1000,
                           nbrCores = 1,
                           nbrCoresDiffMeth = 1,
                           type = c("both", "sites", "tiles"),
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
    validateRunPermutationethylKitObject(methylKitRDS,
                           output_dir, nbrPermutations,
                           nbrCores, nbrCoresDiffMeth,
                           minReads, minMethDiff, qvalue, maxPercReads,
                           destrand, minCovBasesForTiles, tileSize,
                           stepSize, vSeed)

    ## Add last slash to path when absent
    if (substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") {
        output_dir <- paste0(output_dir, "/")
    }

    ## Set vSeed value when negative seed is given
    if (vSeed <= -1) {
        tSeed <- as.numeric(Sys.time())
        vSeed <- 1e8 * (tSeed - floor(tSeed))
    }
    set.seed(vSeed)

    methylInfo <- readRDS(methylKitRDS)

    nbGenerations <- length(methylInfo)

    nbSamplesByGeneration <- sapply(methylInfo, length)
    nbSamples  <- sum(nbSamplesByGeneration)
    allSamples <- unlist(methylInfo, recursive = FALSE)

    finalList <- vector("list", nbrPermutations)

    for (i in 1:nbrPermutations) {

        # Randomnly mixt all samples
        samples <- sample(allSamples, size=nbSamples, replace = FALSE)

        permutationList <- vector("list", nbGenerations)
        start = 1
        for (j in 1:nbGenerations) {
            end = start + nbSamplesByGeneration[j] - 1
            sampleList <- samples[start:end]
            treatment <- methylInfo[[1]]@treatment
            newSampleList <- new("methylRawList", sampleList,
                                    treatment = treatment)
            permutationList[[j]] <- newSampleList
            start <- end + 1
        }
        finalList[[i]] <- permutationList
    }

    result <- runOnePermutationOnAllGenerations(finalList[[1]],
                      nbrCoresDiffMeth = nbrCoresDiffMeth,
                      type = type,
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


#' @title TODO
#'
#' @description TODO
#'
#' @param allFilesByGeneration a \code{list} of \code{list}, a \code{list}
#' composed of \code{list}
#' of files with methylation information for
#' bases or regions in the genome. One \code{list} must contain all files
#' related to the same generation. So, if 3 generations are analyzed, a
#' \code{list} containing 3 \code{list} must be passed. At least 2 generations
#' must be present to do a permutation analysis.
#' The parameter
#' corresponds to the \code{location} parameter in the package \code{methylKit}.
#'
#' @param conditionsByGeneration a \code{list} of \code{vector} containing
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
#' @param output_dir a string, the name of the directory that will contain
#' the results of the permutation. If the directory does not exist, it will
#' be created.
#'
#' @param genomeVersion a string, the genome assembly such as hg18, mm9, etc.
#' It can be any string. The parameter
#' correspond to the \code{assembly} parameter in the package \code{methylKit}.
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
#' @param doingSites a \code{logical}, when \code{TRUE} will do the analysis
#' on the CpG dinucleotide sites. Default: \code{TRUE}.
#'
#' @param doingTiles a \code{logical}, when \code{TRUE} will do the analysis
#' on the tiles. Default: \code{FALSE}.
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
#' @importFrom BiocParallel bplapply MulticoreParam SnowParam bptry bpok bpmapply
#' @importFrom utils flush.console write.table
#' @export
runPermutationTEST <- function(allFilesByGeneration, conditionsByGeneration,
                            output_dir, genomeVersion,
                            nbrPermutations = 1000,  nbrCores = 1,
                            nbrCoresDiffMeth = 1, doingSites = TRUE,
                            doingTiles = FALSE, minReads = 10,
                            minMethDiff = 10, qvalue = 0.01,
                            maxPercReads = 99.9, destrand = FALSE,
                            minCovBasesForTiles = 0, tileSize = 1000,
                            stepSize = 1000, vSeed = -1) {

    ## Parameters validation
    validateRunPermutation(allFilesByGeneration, conditionsByGeneration,
                           output_dir, genomeVersion, nbrPermutations,
                           nbrCores, nbrCoresDiffMeth, doingSites, doingTiles,
                           minReads, minMethDiff, qvalue, maxPercReads,
                           destrand, minCovBasesForTiles, tileSize,
                           stepSize, vSeed)

    ## Add last slash to path when absent
    if (substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") {
        output_dir <- paste0(output_dir, "/")
    }

    ## Set vSeed value when negative seed is given
    if (vSeed <= -1) {
        tSeed <- as.numeric(Sys.time())
        vSeed <- 1e8 * (tSeed - floor(tSeed))
    }
    set.seed(vSeed)

    nbrFiles <- sum(unlist(lapply(allFilesByGeneration, length)))

    nbGenerations <- length(allFilesByGeneration)

    nbSamples  <- sum(length(unlist(allFilesByGeneration)))
    allSamples <- unlist(allFilesByGeneration)

    nbSamplesByGeneration <- sapply(allFilesByGeneration, length)

    finalList <- vector("list", length(nbrPermutations))

    for (i in 1:nbrPermutations) {

        # Randomnly mixt all samples
        samples <- sample(allSamples, size=nbSamples, replace = FALSE)

        start = 1

        onePermutationList <- vector("list", length(nbGenerations))
        for (j in 1:nbGenerations) {

            nbrFilesGeneration <-  nbSamplesByGeneration[j]
            end = start + nbrFilesGeneration - 1
            filesGeneration <- samples[start:end]
            sampleNames <- sapply(filesGeneration, getSampleNameFromFileName)
            infoGeneration <- list(conditions=conditionsByGeneration[[j]],
                                   file.list=as.list(filesGeneration),
                                   designName = paste0("Generation_", j),
                                   sampleNames = sampleNames,
                                   count = i)
            onePermutationList[[j]] <- infoGeneration
            start <- end + 1
        }

        finalList[[i]] <- onePermutationList
    }

    # Create directories for output files
    createOutputDir(output_dir, nbGenerations, doingSites, doingTiles)


    ## Set worker according to number of threads
    if (nbrCores == 1) {
        bp_param <- SnowParam()
    } else {
        bp_param <- MulticoreParam(workers = nbrCores)
    }

    result <- runTOTO(finalList[[1]], output_dir = output_dir,
            genomeVersion = genomeVersion,
            nbrCoresDiffMeth = nbrCoresDiffMeth,
            doingSites = doingSites,
            doingTiles = doingTiles,
            minReads = minReads,
            minMethDiff = minMethDiff,
            qvalue = qvalue,
            maxPercReads = maxPercReads,
            destrand = destrand,
            minCovBasesForTiles = minCovBasesForTiles,
            tileSize = tileSize,
            stepSize = stepSize)
    #result[which(!bpok(result))]

    result
}

