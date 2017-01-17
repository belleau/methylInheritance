#' @title Parameters validation for the
#' \code{\link{runPermutationUsingMethylKitInfo}} function
#'
#' @description Validation of all parameters needed by the public
#' \code{\link{runPermutationUsingMethylKitInfo}} function.
#'
#' @param methylKitInfo a \code{list} of \code{methylRawList} entries. Each
#' \code{methylRawList} contains all the \code{methylRaw} entries related to
#' one generation. The number of generations must correspond to the number
#' of entries in the \code{methylKitInfo}.At least 2 generations
#' must be present to do a permutation analysis. More information can be found
#' in the Bioconductor methylKit package.
#'
#' @param type One of the "sites","tiles" or "both" strings. Specifies the type
#' of differentially methylated elements should be returned. For
#' retrieving differentially methylated bases type="sites"; for
#' differentially methylated regions type="tiles". Default: "both".
#'
#' @param outputDir a string, the name of the directory that will contain
#' the results of the permutation. If the directory does not exist, it will
#' be created.
#'
#' @param runObservedAnalysis a \code{logical}, when \code{runObservedAnalysis}
#' = \code{TRUE}, a CpG analysis on the observed dataset is done.
#'
#' @param nbrPermutations, a positive \code{integer}, the total number of
#' permutations that is going to be done.
#'
#' @param nbrCores a positive \code{integer}, the number of cores to use when
#' processing the analysis.
#'
#' @param nbrCoresDiffMeth a positive \code{integer}, the number of cores
#' to use for parallel differential methylation calculations.Parameter
#' used for both sites and tiles analysis. The parameter
#' corresponds to the \code{num.cores} parameter in
#' the \code{methylKit} package.
#'
#' @param minReads a positive \code{integer} Bases and regions having lower
#' coverage than this count are discarded. The parameter
#' correspond to the \code{lo.count} parameter in the  \code{methylKit} package.
#'
#' @param minMethDiff a positive \code{double} betwwen [0,100], the absolute
#' value of methylation percentage change between cases and controls. The
#' parameter correspond to the \code{difference} parameter in
#' the  \code{methylKit} package.
#'
#' @param qvalue a positive \code{double} betwwen [0,1], the cutoff
#' for qvalue of differential methylation statistic.
#'
#' @param maxPercReads a \code{double} between [0,100], the percentile of read
#' counts that is going to be used as upper cutoff. Bases ore regions
#' having higher
#' coverage than this percentile are discarded. Parameter used for both CpG
#' sites and tiles analysis. The parameter
#' correspond to the \code{hi.perc} parameter in the  \code{methylKit} package.
#'
#' @param destrand a \code{logical}, when \code{TRUE} will merge reads on both
#' strands of a CpG dinucleotide to provide better coverage. Only advised
#' when looking at CpG methylation. Parameter used for both CpG
#' sites and tiles analysis.
#'
#' @param minCovBasesForTiles a non-negative \code{integer}, the minimum
#' number of bases to be covered in a given tiling window. The parameter
#' corresponds to the \code{cov.bases} parameter in the package
#' \code{methylKit}. Only used when \code{doingTiles} =
#' \code{TRUE}. Default: \code{0}.
#'
#' @param tileSize a positive \code{integer}, the size of the tiling window.
#' The parameter corresponds to the \code{win.size} parameter in
#' the  \code{methylKit} package. Only
#' used when \code{doingTiles} = \code{TRUE}.
#'
#' @param stepSize a positive \code{integer}, the step size of tiling windows.
#' The parameter corresponds to the \code{stepSize} parameter in
#' the  \code{methylKit} package. Only
#' used when \code{doingTiles} = \code{TRUE}.
#'
#' @param vSeed a \code{integer}, a seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used.
#'
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @examples
#'
#' ## Load dataset
#' data(samplesForTransgenerationalAnalysis)
#'
#' ## The function returns 0 when all paramaters are valid
#' methylInheritance:::validateRunPermutationUsingMethylKitInfo(
#' methylKitInfo = samplesForTransgenerationalAnalysis, type = "sites",
#' outputDir = NULL, runObservedAnalysis = TRUE,
#' nbrPermutations = 10000, nbrCores = 1,
#' nbrCoresDiffMeth = 1, minReads = 10, minMethDiff = 25, qvalue = 0.01,
#' maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
#' tileSize = 1000, stepSize = 500, vSeed = 12)
#'
#' ## The function raises an error when at least one paramater is not valid
#' \dontrun{methylInheritance:::validateRunPermutationUsingMethylKitInfo(
#' methylKitInfo = "HI",type = "tiles", outputDir = NULL,
#' runObservedAnalysis = FALSE, nbrPermutations = 10000, nbrCores = 1,
#' nbrCoresDiffMeth = 1, minReads = 10, minMethDiff = 25, qvalue = 0.01,
#' maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
#' tileSize = 1000, stepSize = 500, vSeed = 12)}
#'
#' @author Astrid Deschenes
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
validateRunPermutationUsingMethylKitInfo <- function(methylKitInfo,
                                    type, outputDir, runObservedAnalysis,
                                    nbrPermutations, nbrCores,
                                    nbrCoresDiffMeth,
                                    minReads, minMethDiff, qvalue,
                                    maxPercReads, destrand,
                                    minCovBasesForTiles, tileSize,
                                    stepSize, vSeed) {

    ## Validate that methylKitInfo is a list of methylRawList
    if (class(methylKitInfo) != "list" ||
        !all(sapply(methylKitInfo, class) == "methylRawList")) {
        stop(paste0("methylKitInfo must be a list containing ",
                    "\"methylRawList\" entries; each entry must contain ",
                    "all \"methylRaw\" objects related to one generation"))
    }

    ## Validate that the runObservedAnalysis is a logical
    if (!is.logical(runObservedAnalysis)) {
        stop("runObservedAnalysis must be a logical")
    }

    ## Validate that nbrCores is an positive integer
    if (!(isSingleInteger(nbrCores) || isSingleNumber(nbrCores)) ||
        as.integer(nbrCores) < 1) {
        stop("nbrCores must be a positive integer or numeric")
    }

    ## Validate that nbrCores is set to 1 on Windows system
    if (Sys.info()["sysname"] == "Windows" && as.integer(nbrCores) != 1) {
        stop("nbrCores must be 1 on a Windows system.")
    }

    ## Validate that nbrPermutations is an positive integer
    if (!(isSingleInteger(nbrPermutations) ||
          isSingleNumber(nbrPermutations)) ||
        as.integer(nbrPermutations) < 1) {
        stop("nbrPermutations must be a positive integer or numeric")
    }

    ## Validate all the other parameters
    validateRunAnalysisUsingMethylKitInfo(methylKitInfo = methylKitInfo,
                           type = type, outputDir = outputDir,
                           nbrCores = nbrCores,
                           nbrCoresDiffMeth = nbrCoresDiffMeth,
                           minReads = minReads, minMethDiff = minMethDiff,
                           qvalue = qvalue,
                           maxPercReads = maxPercReads, destrand = destrand,
                           minCovBasesForTiles = minCovBasesForTiles,
                           tileSize = tileSize,
                           stepSize = stepSize, vSeed = vSeed)
}


#' @title Validation of some parameters of the
#' \code{\link{runAnalysisUsingMethylKitInfo}} function
#'
#' @description Validation of some parameters needed by the public
#' \code{\link{runAnalysisUsingMethylKitInfo}} function.
#'
#' @param methylKitInfo a \code{list} of \code{methylRawList} entries. Each
#' \code{methylRawList} contains all the \code{methylRaw} entries related to
#' one generation. The number of generations must correspond to the number
#' of entries in the \code{methylKitInfo}.At least 2 generations
#' must be present to do a permutation analysis. More information can be found
#' in the Bioconductor methylKit package.
#'
#' @param type One of the "sites","tiles" or "both" strings. Specifies the type
#' of differentially methylated elements should be returned. For
#' retrieving differentially methylated bases type="sites"; for
#' differentially methylated regions type="tiles". Default: "both".
#'
#' @param outputDir a string, the name of the directory that will contain
#' the results of the permutation. If the directory does not exist, it will
#' be created.
#'
#' @param nbrCores a positive \code{integer}, the number of cores to use when
#' processing the analysis.
#'
#' @param nbrCoresDiffMeth a positive \code{integer}, the number of cores
#' to use for parallel differential methylation calculations.Parameter
#' used for both sites and tiles analysis. The parameter
#' corresponds to the \code{num.cores} parameter in
#' the \code{methylKit} package.
#'
#' @param minReads a positive \code{integer} Bases and regions having lower
#' coverage than this count are discarded. The parameter
#' correspond to the \code{lo.count} parameter in the  \code{methylKit}
#' package.
#'
#' @param minMethDiff a positive \code{double} betwwen [0,100], the absolute
#' value of methylation percentage change between cases and controls. The
#' parameter correspond to the \code{difference} parameter in
#' the \code{methylKit} package.
#'
#' @param qvalue a positive \code{double} betwwen [0,1], the cutoff
#' for qvalue of differential methylation statistic.
#'
#' @param maxPercReads a \code{double} between [0,100], the percentile of read
#' counts that is going to be used as upper cutoff. Bases ore regions
#' having higher
#' coverage than this percentile are discarded. Parameter used for both CpG
#' sites and tiles analysis. The parameter
#' correspond to the \code{hi.perc} parameter in the  \code{methylKit} package.
#'
#' @param destrand a \code{logical}, when \code{TRUE} will merge reads on both
#' strands of a CpG dinucleotide to provide better coverage. Only advised
#' when looking at CpG methylation. Parameter used for both CpG
#' sites and tiles analysis.
#'
#' @param minCovBasesForTiles a non-negative \code{integer}, the minimum
#' number of bases to be covered in a given tiling window. The parameter
#' corresponds to the \code{cov.bases} parameter in the package
#' \code{methylKit}. Only used when \code{doingTiles} =
#' \code{TRUE}. Default: \code{0}.
#'
#' @param tileSize a positive \code{integer}, the size of the tiling window.
#' The parameter corresponds to the \code{win.size} parameter in
#' the  \code{methylKit} package. Only
#' used when \code{doingTiles} = \code{TRUE}.
#'
#' @param stepSize a positive \code{integer}, the step size of tiling windows.
#' The parameter corresponds to the \code{stepSize} parameter in
#' the  \code{methylKit} package. Only
#' used when \code{doingTiles} = \code{TRUE}.
#'
#' @param vSeed a \code{integer}, a seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used.
#'
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @examples
#'
#' ## Load dataset
#' data(samplesForTransgenerationalAnalysis)
#'
#' ## The function returns 0 when all paramaters are valid
#' methylInheritance:::validateRunAnalysisUsingMethylKitInfo(
#' methylKitInfo = samplesForTransgenerationalAnalysis, type = "sites",
#' outputDir = NULL, nbrCores = 1,
#' nbrCoresDiffMeth = 1, minReads = 10, minMethDiff = 25, qvalue = 0.01,
#' maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
#' tileSize = 1000, stepSize = 500, vSeed = 12)
#'
#' ## The function raises an error when at least one paramater is not valid
#' \dontrun{methylInheritance:::validateRunAnalysisUsingMethylKitInfo(
#' methylKitInfo = samplesForTransgenerationalAnalysis,
#' type = "tiles", outputDir = NULL, nbrCores = 1,
#' nbrCoresDiffMeth = 1, minReads = "HI", minMethDiff = 25, qvalue = 0.01,
#' maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
#' tileSize = 1000, stepSize = 500, vSeed = 12)}
#'
#' @author Astrid Deschenes
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
validateRunAnalysisUsingMethylKitInfo <- function(methylKitInfo,
                                    type, outputDir, nbrCores,
                                    nbrCoresDiffMeth,
                                    minReads, minMethDiff, qvalue,
                                    maxPercReads, destrand,
                                    minCovBasesForTiles, tileSize,
                                    stepSize, vSeed) {

    ## Validate that methylKitInfo is a list of methylRawList
    if (class(methylKitInfo) != "list" ||
        !all(sapply(methylKitInfo, class) == "methylRawList")) {
        stop(paste0("methylKitInfo must be a list containing ",
                    "\"methylRawList\" entries; each entry must contain ",
                    "all \"methylRaw\" objects related to one generation"))
    }

    ## Validate that the output_dir is an not empty string
    if (!is.null(outputDir) && !is.character(outputDir)) {
        stop("output_dir must be a character string or NULL")
    }

    ## Validate that nbrCores is an positive integer
    if (!(isSingleInteger(nbrCores) || isSingleNumber(nbrCores)) ||
        as.integer(nbrCores) < 1) {
        stop("nbrCores must be a positive integer or numeric")
    }

    ## Validate that nbrCores is set to 1 on Windows system
    if (Sys.info()["sysname"] == "Windows" && as.integer(nbrCores) != 1) {
        stop("nbrCores must be 1 on a Windows system.")
    }

    ## Validate that nbrCoresDiffMeth is an positive integer
    if (!(isSingleInteger(nbrCoresDiffMeth) ||
          isSingleNumber(nbrCoresDiffMeth)) ||
        as.integer(nbrCoresDiffMeth) < 1) {
        stop("nbrCoresDiffMeth must be a positive integer or numeric")
    }

    ## Validate that nbrCoresDiffMeth is set to 1 on Windows system
    if (Sys.info()["sysname"] == "Windows" &&
        as.integer(nbrCoresDiffMeth) != 1) {
        stop("nbrCoresDiffMeth must be 1 on a Windows system.")
    }

    ## Validate that minReads is an positive integer
    if (!(isSingleInteger(minReads) || isSingleNumber(minReads)) ||
        as.integer(minReads) < 1) {
        stop("minReads must be a positive integer or numeric")
    }

    ## Validate that minMethDiff is an positive double between [0,100]
    if (!(isSingleNumber(minMethDiff)) ||
        minMethDiff < 0.00 || minMethDiff > 100.00) {
        stop("minMethDiff must be a positive double between [0,100]")
    }

    ## Validate that qvalue is an positive double between [0,1]
    if (!(isSingleNumber(qvalue)) ||
        qvalue < 0.00 || qvalue > 1.00) {
        stop("qvalue must be a positive double between [0,1]")
    }

    ## Validate that maxPercReads is an positive double between [0,100]
    if (!(isSingleNumber(maxPercReads)) ||
        maxPercReads < 0.00 || maxPercReads > 100.00) {
        stop("maxPercReads must be a positive double between [0,100]")
    }

    ## Validate that destrand is a logical
    if (!is.logical(destrand)) {
        stop("destrand must be a logical")
    }

    if (any(type %in% c("both", "tiles"))) {
        ## Validate that minCovBasesForTiles is an positive integer
        if (!(isSingleInteger(minCovBasesForTiles) ||
              isSingleNumber(minCovBasesForTiles)) ||
            as.integer(minCovBasesForTiles) < 0) {
            stop("minCovBasesForTiles must be a positive integer or numeric")
        }

        ## Validate that tileSize is an positive integer
        if (!(isSingleInteger(tileSize) || isSingleNumber(tileSize)) ||
            as.integer(tileSize) < 1) {
            stop("tileSize must be a positive integer or numeric")
        }

        ## Validate that stepSize is an positive integer
        if (!(isSingleInteger(stepSize) || isSingleNumber(stepSize)) ||
            as.integer(stepSize) < 1) {
            stop("stepSize must be a positive integer or numeric")
        }

    }
    ## Validate that vSeed is an integer
    if (!(isSingleInteger(vSeed) || isSingleNumber(vSeed))) {
        stop("vSeed must be an integer or numeric")
    }

    return(0)
}


#' @title TODO
#'
#' @description TODO
#'
#' @param methDiff a S4 \code{methylDiff} class object, a
#' object that holds statistics and locations
#' for differentially methylated regions/bases.
#'
#' @param pDiff a positive \code{double} between \code{0} and \code{100},
#' the cutoff for absolute value of methylation percentage change
#' between test and control.
#'
#' @param qvalue a positive \code{double} inferior to \code{1}, the cutoff
#' for qvalue of differential methylation statistic.
#'
#' @param type One of the \code{"hyper"},\code{"hypo"} or \code{"all"} strings,
#' the string specifies what type of differentially methylated bases/tiles
#' should be treated  For
#' retrieving hyper-methylated tiles/sites \code{type} = \code{"hyper"}; for
#' hypo-methylated \code{type} = \code{"hypo"}. Default: \code{"all"}.
#'
#' @return TODO
#'
#' @examples
#'
#' ## Load permutation results on sites
#' permutationResultsFile <- dir(system.file("extdata",
#' package = "methylInheritance"), pattern = "permutationResultsForSites.RDS",
#' full.names = TRUE)
#' permutationResults <- readRDS(permutationResultsFile)
#'
#' ## Transform result to GRanges
#' resultsGR <- methylInheritance:::getGRangesFromMethylDiff(methDiff =
#' permutationResults, pDiff = 10, qvalue = 0.01, type = "hyper")
#'
#' @author Pascal Belleau
#' @importFrom methylKit getMethylDiff
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @keywords internal
getGRangesFromMethylDiff <- function(methDiff, pDiff, qvalue,
                                        type = c("all", "hyper", "hypo")) {

    methDiffK <- lapply(1:length(methDiff), FUN = function(i, methDiff,
                                                    pDiff, qCut, typeD){
        methK <- getMethylDiff(methDiff[[i]], difference = pDiff,
                                qvalue = qCut, type = typeD)
        GRanges(seqnames = methK$chr, ranges = IRanges(start = methK$start,
                                                        end = methK$end),
                strand = methK$strand, pvalue = methK$pvalue,
                qvalue = methK$qvalue, meth.diff = methK$meth.diff)
    }, methDiff = methDiff, pDiff = pDiff, qCut = qvalue, typeD = type)

    return(methDiffK)
}


#' @title TODO
#'
#' @description TODO
#'
#' @param resultAllGenGR \code{GRanges} from \code{getGRangesFromMethylDiff}, TODO
#'
#' @return \code{list} with 2 elements
#'         i2 list of intersection G1 and G2, G2 and G3, ...
#'         iAll list of intersection G1 and G2 and G3,
#'              G1 and G2 and G3 and G4, ...
#'         TODO
#'
#' @examples
#'
#' ## Load permutation results on sites
#' permutationResultsFile <- dir(system.file("extdata",
#' package = "methylInheritance"), pattern = "permutationResultsForSites.RDS",
#' full.names = TRUE)
#' permutationResults <- readRDS(permutationResultsFile)
#'
#' ## Transform result to GRanges
#' resultsGR <- methylInheritance:::getGRangesFromMethylDiff(methDiff =
#' permutationResults, pDiff = 10, qvalue = 0.01, type = "hyper")
#'
#' ## Extract inter generational conserved sites
#' conservedSitesGR <- interGeneration(resultsGR)
#'
#' @author Pascal Belleau
#' @importFrom GenomicRanges intersect GRanges
#' @importFrom S4Vectors DataFrame values<- values
#' @keywords internal
interGeneration <- function(resultAllGenGR){

    lInter <- list("i2" = list(), "iAll" = list())

    lInter$i2 <- lapply(2:length(resultAllGenGR), FUN = function(i,b){
        upM <- intersect(b[[i-1]][b[[i-1]]$meth.diff > 0],
                            b[[i]][b[[i]]$meth.diff > 0])
        downM <- intersect(b[[i-1]][b[[i-1]]$meth.diff < 0],
                            b[[i]][b[[i]]$meth.diff < 0])
        typeDiff <- DataFrame(typeDiff=rep(1,length(upM)))
        values(upM) <- cbind(values(upM), typeDiff)
        typeDiff <- DataFrame(typeDiff=rep(-1,length(downM)))
        values(downM) <- cbind(values(downM), typeDiff)
        c(upM,downM)
    }, b = resultAllGenGR)

    cur <- lInter$i2[[1]]
    for(i in 3:length(resultAllGenGR)){
        upM <- intersect(cur[cur$typeDiff > 0],
                        resultAllGenGR[[i]][resultAllGenGR[[i]]$meth.diff > 0])
        downM <- intersect(cur[cur$typeDiff < 0],
                            resultAllGenGR[[i]][
                            resultAllGenGR[[i]]$meth.diff < 0])
        typeDiff <- DataFrame(typeDiff=rep(1,length(upM)))
        values(upM) <- cbind(values(upM), typeDiff)
        typeDiff <- DataFrame(typeDiff=rep(-1,length(downM)))
        values(downM) <- cbind(values(downM), typeDiff)

        lInter$iAll[[i-2]] <- c(upM,downM)
        cur <- lInter$iAll[[i-2]]
    }

    return(lInter)
}


#' @title Create directories that will contained the results of the
#' permutations in RDS format.
#'
#' @description Create directories that will contained the results of the
#' permutations in RDS format.
#'
#' @param outputDir a string of \code{character}, the name of the main
#' directory to be created.
#'
#' @param doingSites a \code{logical}, a directory consecrated to contain the
#' results of the permutation analysis for sites is created when
#' \code{doingSites} = \code{TRUE}. Default: \code{TRUE}.
#'
#' @param doingTiles a \code{logical}, a directory consecrated to contain the
#' results of the permutation analysis for tiles is created when
#' \code{doingTiles} = \code{TRUE}. Default: \code{FALSE}.
#'
#' @return \code{0} when all directories are created without problem.
#'
#' @examples
#'
#' ## Create an output directory for SITES only
#' \dontrun{createOutputDir(outputDir = "testSites", doingSites = TRUE,
#' doingTiles = FALSE)}
#'
#' @author Astrid Deschenes
#' @keywords internal
createOutputDir <- function(outputDir, doingSites = TRUE,
                                doingTiles = FALSE) {

    # Create directories for output files
    if (!dir.exists(outputDir)) {
        dir.create(outputDir, showWarnings = TRUE)
    }

    if (doingSites) {
        type <-  "SITES"
        dirName <- paste0(outputDir, type)
        if (!dir.exists(dirName)) {
            dir.create(dirName, showWarnings = TRUE)
        }
    }

    if (doingTiles) {
        type <-  "TILES"
        dirName <- paste0(outputDir, type)
        if (!dir.exists(dirName)) {
            dir.create(dirName, showWarnings = TRUE)
        }
    }

    return(0)
}


#' @title Run the analysis on one permutation dataset using
#' \code{methylKit} package. One permutation dataset
#' includes data for all generations.
#'
#' @description Run CpG site or region analysis using the \code{methylKit}
#' package for each generation present in the dataset. The intersection of
#' conserved elements is obtained for each group of two consecutive
#' generations, as well as, for larger group subset. The output of the
#' analysis is saved in a RDS file when an directory is
#' specified.
#'
#' @param methylInfoForAllGenerations a \code{list} containing the
#' following elements:
#' \itemize{
#' \item \code{sample} a \code{list} of \code{methylRawList} entries, each
#' \code{methylRawList} contains all the \code{methylRaw} entries related to
#' one generation. The number of generations must correspond to the number
#' of entries in the \code{methylKitInfo}. At least 2 generations
#' must be present to do a permutation analysis.
#' \item \code{id} an integer, the permutation id.
#' }
#'
#' @param type One of the "sites","tiles" or "both" strings. Specifies the type
#' of differentially methylated elements should be returned. For
#' retrieving differentially methylated bases type="sites"; for
#' differentially methylated regions type="tiles". Default: "both".
#'
#' @param nbrCoresDiffMeth a positive integer, the number of cores to use for
#' parallel differential methylation calculations.Parameter used for both
#' sites and tiles analysis. The parameter
#' corresponds to the \code{num.cores} parameter in the
#' package \code{methylKit}.
#' Default: \code{1} and always \code{1} for Windows.
#'
#' @param minReads a positive \code{integer} Bases and regions having lower
#' coverage than this count are discarded. The parameter
#' correspond to the \code{lo.count} parameter in the  \code{methylKit} package.
#'
#' @param qvalue a positive \code{double} inferior to \code{1}, the cutoff
#' for qvalue of differential methylation statistic. Default: \code{0.01}.
#'
#' @param maxPercReads a double between [0-100], the percentile of read
#' counts that is going to be used as upper cutoff. Bases ore regions
#' having higher
#' coverage than this percentile are discarded. Parameter used for both CpG
#' sites and tiles analysis. The parameter
#' correspond to the \code{hi.perc} parameter in the  \code{methylKit} package.
#' Default: \code{99.9}.
#'
#' @param minMethDiff a positive integer betwwen [0,100], the absolute value
#' of methylation
#' percentage change between cases and controls. The parameter
#' correspond to the \code{difference} parameter in the
#' package \code{methylKit}.
#' Default: \code{10}.
#'
#' @param destrand a logical, when \code{TRUE} will merge reads on both
#' strands of a CpG dinucleotide to provide better coverage. Only advised
#' when looking at CpG methylation. Parameter used for both
#' sites and tiles analysis.
#' Default: \code{FALSE}.
#'
#' @param minCovBasesForTiles a non-negative integer, the minimum number of
#' bases to be covered in a given tiling window. The parameter
#' corresponds to the \code{cov.bases} parameter in the
#' package \code{methylKit}.
#' Only used when \code{doingTiles} =
#' \code{TRUE}. Default: \code{0}.
#'
#' @param tileSize a positive integer, the size of the tiling window. The
#' parameter corresponds to the \code{win.size} parameter in
#' the  \code{methylKit} package. Only
#' used when \code{doingTiles} = \code{TRUE}. Default: \code{1000}.
#'
#' @param stepSize a positive integer, the step size of tiling windows. The
#' parameter corresponds to the \code{stepSize} parameter in
#' the  \code{methylKit} package. Only
#' used when \code{doingTiles} = \code{TRUE}. Default: \code{1000}.
#'
#' @param doingSites a logical, when \code{TRUE} will do the analysis on the
#' CpG dinucleotide sites. Default: \code{TRUE}.
#'
#' @param doingTiles a logical, when \code{TRUE} will do the analysis on the
#' tiles. Default: \code{FALSE}.
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
#' info <- list(sample = samplesForTransgenerationalAnalysis, id = 100)
#'
#' ## Run a permutation analysis
#' \dontrun{methylInheritance:::runOnePermutationOnAllGenerations(
#' methylInfoForAllGenerations = info, type = "sites", outputDir = NULL,
#' nbrCoresDiffMeth = 1, minReads = 10, minMethDiff = 10, qvalue = 0.01,
#' maxPercReads = 99.9, destrand = FALSE, minCovBasesForTiles = 0,
#' tileSize = 1000, stepSize = 1000)}
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom methylKit filterByCoverage normalizeCoverage unite calculateDiffMeth getMethylDiff getData tileMethylCounts methRead
#' @importFrom GenomicRanges width
#' @keywords internal
runOnePermutationOnAllGenerations <- function(methylInfoForAllGenerations,
                        type = c("both", "sites", "tiles"),
                        outputDir = NULL,
                        nbrCoresDiffMeth = 1,
                        minReads = 10, minMethDiff = 10,
                        qvalue = 0.01, maxPercReads = 99.9,
                        destrand = FALSE, minCovBasesForTiles = 0,
                        tileSize = 1000, stepSize = 1000) {

    doTiles <- any(type %in% c("tiles", "both"))
    doSites <- any(type %in% c("sites", "both"))

    ## Extract info from input list
    methylRawForAllGenerations <- methylInfoForAllGenerations$sample
    id <- methylInfoForAllGenerations$id

    nbrGenerations <- length(methylRawForAllGenerations)

    ## Preparing list that will receive final results
    permutationList <- list()
    if (doTiles) {
        permutationList[["TILES"]] <- list()
    }
    if (doSites) {
        permutationList[["SITES"]] <- list()
    }

    for (i in 1:nbrGenerations) {

        allSamplesForOneGeneration <- methylRawForAllGenerations[[i]]

        ## SITES
        if (doSites) {

            ## Filter sites by coverage
            filtered.sites <- filterByCoverage(allSamplesForOneGeneration,
                                                lo.count = minReads,
                                                lo.perc = NULL,
                                                hi.count = NULL,
                                                hi.perc = maxPercReads)

            ## Normalize coverage
            filtered.sites <- normalizeCoverage(filtered.sites, "median")

            ## Merge all samples to one table
            meth.sites <- unite(filtered.sites, destrand = destrand)

            if (length(meth.sites@.Data[[1]]) == 0) {
                stop("meth.sites IS EMPTY")
            }

            ## Get differentially methylated sites
            permutationList[["SITES"]][[i]] <- calculateDiffMeth(meth.sites,
                                                num.cores = nbrCoresDiffMeth)
        }

        ## TILES
        if (doTiles) {

            ## Summarize methylated base counts over tilling windows
            tiles <- tileMethylCounts(allSamplesForOneGeneration,
                                        win.size = tileSize,
                                        step.size = stepSize,
                                        cov.bases = minCovBasesForTiles)

            ## Filter tiles by coverage
            filtered.tiles <- filterByCoverage(tiles,
                                                lo.count = minReads,
                                                lo.perc = NULL,
                                                hi.count = NULL,
                                                hi.perc = maxPercReads)

            ## Normalize coverage
            filtered.tiles <- normalizeCoverage(filtered.tiles, "median")

            ## Merge all samples to one table
            meth.tiles <- unite(filtered.tiles, destrand = destrand)

            ## Get diff methylated tiles
            permutationList[["TILES"]][[i]] <- calculateDiffMeth(meth.tiles,
                                                num.cores = nbrCoresDiffMeth)
        }
    }

    permutationFinal <- list()

    ## Calculate the number of SITES in the intersection
    if (doSites) {
        permutationFinal[["SITES"]] <- list()
        permutationFinal[["SITES"]][["i2"]] <- list()
        permutationFinal[["SITES"]][["i2"]][["HYPER"]] <- list()
        permutationFinal[["SITES"]][["i2"]][["HYPO"]]  <- list()
        permutationFinal[["SITES"]][["iAll"]][["HYPER"]]  <- list()
        permutationFinal[["SITES"]][["iAll"]][["HYPO"]]   <- list()

        resultGR <- getGRangesFromMethylDiff(permutationList[["SITES"]],
                                        minMethDiff, qvalue, type = "all")

        result <- interGeneration(resultGR)

        if (!is.null(outputDir)) {
            saveInterGenerationResults(outputDir, id, type = "sites", result)
        }

        permutationFinal[["SITES"]][["i2"]][["HYPER"]] <- lapply(result$i2,
                        FUN = function(x) {sum(width(x[x$typeDiff > 0]))})

        permutationFinal[["SITES"]][["i2"]][["HYPO"]]  <- lapply(result$i2,
                        FUN = function(x) {sum(width(x[x$typeDiff < 0]))})


        permutationFinal[["SITES"]][["iAll"]][["HYPER"]] <- lapply(result$iAll,
                        FUN = function(x) {sum(width(x[x$typeDiff > 0]))})

        permutationFinal[["SITES"]][["iAll"]][["HYPO"]]  <- lapply(result$iAll,
                        FUN = function(x) {sum(width(x[x$typeDiff < 0]))})
    }

    ## Calculate the number of TILES in the intersection
    if (doTiles) {
        permutationFinal[["TILES"]] <- list()
        permutationFinal[["TILES"]][["i2"]] <- list()
        permutationFinal[["TILES"]][["i2"]][["HYPER"]] <- list()
        permutationFinal[["TILES"]][["i2"]][["HYPO"]]  <- list()
        permutationFinal[["TILES"]][["iAll"]][["HYPER"]]  <- list()
        permutationFinal[["TILES"]][["iAll"]][["HYPO"]]   <- list()

        resultGR <- getGRangesFromMethylDiff(permutationList[["TILES"]],
                                minMethDiff, qvalue, type = "all")

        result <- interGeneration(resultGR)

        if (!is.null(outputDir)) {
            saveInterGenerationResults(outputDir, id, type = "tiles", result)
        }

        permutationFinal[["TILES"]][["i2"]][["HYPER"]] <- lapply(result$i2,
                            FUN = function(x) {sum(width(x[x$typeDiff > 0]))})

        permutationFinal[["TILES"]][["i2"]][["HYPO"]]  <- lapply(result$i2,
                            FUN = function(x) {sum(width(x[x$typeDiff < 0]))})

        permutationFinal[["TILES"]][["iAll"]][["HYPER"]] <- lapply(result$iAll,
                            FUN = function(x) {sum(width(x[x$typeDiff > 0]))})

        permutationFinal[["TILES"]][["iAll"]][["HYPO"]]  <- lapply(result$iAll,
                            FUN = function(x) {sum(width(x[x$typeDiff < 0]))})
    }

    return(permutationFinal)
}


#' @title Save the result of on CpG site or tile analysis on all generations.
#' The anaysis can come from observed or permutated dataset. Each case is
#' saved with a different extension.
#'
#' @description Save the result of on CpG site or tile analysis on all
#' generations. The results are saved in a RDS file. The anaysis can come
#' from observed or permutated dataset.
#' Each case is saved with a different extension. The files containing the
#' permutation results have the permutation identifiant in their name.
#'
#' @param outputDir a string, the name of the directory that will contain
#' the results of the permutation. The name should end with a slash. The
#' directory should already exists.
#'
#' @param type One of the \code{"sites"} or \code{"tiles"} strings. Specifies
#' the type
#' of differentially methylated elements should be returned. For
#' retrieving differentially methylated bases \code{type} =\code{"sites"}; for
#' differentially methylated regions \code{type} = \code{"tiles"}. Default:
#' \code{"both"}.
#'
#' @param permutationID an \code{integer}, the identifiant of the permutation.
#' When the \code{permutationID} = \code{0}, the results are considered as the
#' observed results and are saved in a file with the "_observed_results.RDS"
#' extension. When the \code{permutationID} != \code{0}, the results are
#' considered as permutation results and are saved in a file with the
#' "_permutation_{permutationID}.RDS" extension.
#'
#' @param type One of the \code{"sites"} or \code{"tiles"} strings. Specifies
#' the type of differentially methylated elements should be saved.
#' Default: \code{"sites"}.
#'
#' @param interGenerationResult a \code{list} that corresponds to the output
#' of the \code{interGeneration} function, the result of on CpG site or tile
#' analysis on all generations.
#'
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @examples
#'
#' ## Load a dataset
#' interGenerationResultFile <- dir(system.file("extdata",
#' package = "methylInheritance"), pattern = "interGenerationResult",
#' full.names = TRUE)
#' interGenerationResult <- readRDS(interGenerationResultFile)
#'
#' ## Save dataset
#' \dontrun{methylInheritance:::saveInterGenerationResults(
#' outputDir = "TEST", permutationID=100, type = "sites",
#' interGenerationResult = interGenerationResult)}
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @keywords internal
saveInterGenerationResults <- function(outputDir, permutationID,
                                        type = c("sites", "tiles"),
                                        interGenerationResult) {

    if (permutationID != 0) {
        ## Save the permutation results
        saveRDS(object = interGenerationResult,
                file = paste0(outputDir,  toupper(type), "/",
                        toupper(type), "_permutation_", permutationID, ".RDS"))
    } else {
        ## Save the observed results
        saveRDS(object = interGenerationResult,
                    file = paste0(outputDir, toupper(type), "/",
                        toupper(type), "_observed_results.RDS"))
    }

    return(0)
}


#' @title TODO
#'
#' @description TODO
#'
#' @param interGenerationResult TODO
#'
#' @param type One of the "sites" or "tiles" strings. Specifies the type
#' of elements that should be returned. For
#' retrieving differentially methylated bases type="sites"; for
#' differentially methylated regions type="tiles". Default: "sites".
#'
#' @param result TODO
#'
#' @return a \code{list} containing:
#'
#' @examples
#'
#' ##TODO
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @keywords internal
createDataStructure <- function(interGenerationResult) {

    result <- list()
    result[["i2"]] <- list()
    result[["i2"]][["HYPER"]] <- list()
    result[["i2"]][["HYPO"]]  <- list()
    result[["iAll"]][["HYPER"]]  <- list()
    result[["iAll"]][["HYPO"]]   <- list()
    result[["i2"]][["HYPER"]] <- lapply(interGenerationResult$i2,
                        FUN = function(x) {sum(width(x[x$typeDiff > 0]))})
    result[["i2"]][["HYPO"]]  <- lapply(interGenerationResult$i2,
                        FUN = function(x) {sum(width(x[x$typeDiff < 0]))})
    result[["iAll"]][["HYPER"]] <- lapply(interGenerationResult$iAll,
                        FUN = function(x) {sum(width(x[x$typeDiff > 0]))})
    result[["iAll"]][["HYPO"]]  <- lapply(interGenerationResult$iAll,
                        FUN = function(x) {sum(width(x[x$typeDiff < 0]))})

    return(result)
}


