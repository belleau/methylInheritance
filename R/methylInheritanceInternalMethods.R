#' @title Run one analysis using \code{methylKit} package
#'
#' @description Run one CpG site or region analysis using the \code{methylKit}
#' package. The output of the analysis is saved in a file in the specified
#' directory.
#'
#' @param info TODO
#'
#' @param output_dir a string or \code{NULL}, the name of the directory that
#' will contain the results of the permutation. If the directory does not
#' exist, it will be created. When \code{NULL}, the results fo the permutation
#' are not saved.
#'
#' @param genomeVersion a string, the genome assembly such as hg18, mm9, etc.
#' It can be any string. The parameter
#' correspond to the \code{assembly} parameter in the \code{methylKit} package.
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
#' @param nbrCoresDiffMeth a positive integer, the number of cores to use for
#' parallel differential methylation calculations.Parameter used for both
#' sites and tiles analysis. The parameter
#' corresponds to the \code{num.cores} parameter in the
#' package \code{methylKit}.
#' Default: \code{1} and always \code{1} for Windows.
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
#' @param debug Default: \code{FALSE}.
#'
#' @return TODO
#'
#' @examples
#'
#' ##TODO
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom methylKit read filterByCoverage normalizeCoverage unite calculateDiffMeth get.methylDiff getData tileMethylCounts methRead
#' @keywords internal
runOnePermutation <- function(info, output_dir,
                                genomeVersion,
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
                                stepSize = 1000) {

    designName  <- info$designName
    count       <- info$count
    file.list   <- info$file.list
    sampleNames <- info$sampleNames
    conditions  <- info$conditions

    ####################################
    ## Prepare data
    ####################################

    myobj <- read(location = file.list,
                                sample.id = as.list(sampleNames),
                                assembly = genomeVersion,
                                context = "CpG",
                                treatment = conditions)

    # myobj <- methRead(location = file.list,
    #                          sample.id = as.list(sampleNames),
    #                          assembly = genomeVersion,
    #                          context = "CpG",
    #                          treatment = conditions)

    ## SITES
    if (doingSites) {
        filtered.sites <- filterByCoverage(myobj,
                                            lo.count = minReads,
                                            lo.perc = NULL,
                                            hi.count = NULL,
                                            hi.perc = maxPercReads)

        filtered.sites <- normalizeCoverage(filtered.sites, "median")

        meth.sites <- unite(filtered.sites, destrand = destrand)

        if (length(meth.sites@.Data[[1]]) == 0) {
            stop("meth.sites IS EMPTY")
        }

        ####################################
        ## Get diff methyl sites
        ####################################

        myDiff.sites <- calculateDiffMeth(meth.sites,
                                            num.cores = nbrCoresDiffMeth)

        myDiff.sites.hyper <- get.methylDiff(myDiff.sites,
                                                difference = minMethDiff,
                                                qvalue = qvalue,
                                                type = "hyper")

        myDiff.sites.hypo  <- get.methylDiff(myDiff.sites,
                                                difference = minMethDiff,
                                                qvalue = qvalue,
                                                type = "hypo")
    }

    ## TILES
    if (doingTiles) {
        tiles <- tileMethylCounts(myobj,
                                    win.size = tileSize,
                                    step.size = stepSize,
                                    cov.bases = minCovBasesForTiles)

        filtered.tiles <- filterByCoverage(tiles,
                                            lo.count = minReads,
                                            lo.perc = NULL,
                                            hi.count = NULL,
                                            hi.perc = maxPercReads)

        filtered.tiles <- normalizeCoverage(filtered.tiles, "median")

        meth.tiles <- unite(filtered.tiles, destrand = destrand)

        ####################################
        ## Get diff methyl tiles
        ####################################
        myDiff.tiles <- calculateDiffMeth(meth.tiles,
                                            num.cores = nbrCoresDiffMeth)

        myDiff.tiles.hyper <- get.methylDiff(myDiff.tiles,
                                                difference = minMethDiff,
                                                qvalue = qvalue,
                                                type = "hyper")

        myDiff.tiles.hypo  <- get.methylDiff(myDiff.tiles,
                                                difference = minMethDiff,
                                                qvalue = qvalue,
                                                type = "hypo")
    }

    #####################
    ## Save data
    #####################

    if (doingSites) {
        printDiffMethylFile(data = myDiff.sites.hypo,
                                output_dir = output_dir, count = count,
                                designName = designName, isHyper = FALSE,
                                isSites = TRUE)

        printDiffMethylFile(data = myDiff.sites.hyper, output_dir = output_dir,
                                count = count,
                                designName = designName, isHyper = TRUE,
                                isSites = TRUE)
    }

    if (doingTiles) {
        printDiffMethylFile(data = myDiff.tiles.hypo, output_dir = output_dir,
                                count = count,
                                designName = designName, isHyper = FALSE,
                                isSites = FALSE)

        printDiffMethylFile(data = myDiff.tiles.hyper, output_dir = output_dir,
                                count = count,
                                designName = designName, isHyper = TRUE,
                                isSites = FALSE)
    }

    warnings()
    return(0)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param data TODO
#'
#' @param output_dir TODO
#'
#' @param count TODO
#'
#' @param designName TODO
#'
#' @param isHyper TODO
#'
#' @param isSites TODO
#'
#' @return TODO
#'
#' @examples
#'
#' ## TODO
#'
#'
#' @author Astrid Deschenes
#' @importFrom methylKit getData
#' @keywords internal
printDiffMethylFile <- function(data, output_dir, count, designName, isHyper,
                                    isSites) {

    dirName <- "/SITES/"
    nameExtension <- ".perbase.txt"
    if (!isSites) {
        dirName <- "/TILES/"
        nameExtension <- ".pertile.txt"
    }

    if (isHyper) {
        nameExtension <- paste0("_hyper", nameExtension)
    } else {

        nameExtension <- paste0("_hypo", nameExtension)
    }

    if(nrow(data) > 0) {
        dataAll <- getData(data)
        dataForFile <- cbind(dataAll[,c(1,2,3)],
                        paste(dataAll[,1], dataAll[,2], dataAll[,3], sep = "."),
                        dataAll[,c(7,4,5,6)])
        colnames(dataForFile)[4]="DMR.ID"
        write.table(dataForFile,
                    paste0(output_dir, dirName, designName, "/", count,
                            nameExtension), quote=F, row.names=F,
                            col.names=F, sep="\t")
    } else {
        write.table(NULL,
                    paste0(output_dir, dirName, designName, "/", count,
                           nameExtension), quote=F, row.names=F,
                    col.names=F, sep="\t")
    }

    return(0)
}

#' @title Extract sample name from file name
#'
#' @description Extract a sample name from the file name. The sample name
#' corresponds to  the file name but without path information.
#'
#' @param fileName a string, the file name used to extract the sample name.
#'
#' @return A sample name extracted from the specified file name.
#'
#' @examples
#'
#' ## Extract sample name from file name
#' methylInheritance:::getSampleNameFromFileName("./test/data001/file_J1.txt")
#'
#' @author Astrid Deschenes
#' @keywords internal
getSampleNameFromFileName <- function(fileName) {
    results <- strsplit(fileName, split="/")[[1]]
    return(results[length(results)])
}


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
#' ## The function returns 0 when all paramaters are valid
#' #methylInheritance:::validateRunPermutationUsingRDS(
#' #allFilesByGeneration = list(list("file01.txt", "file02.txt"),
#' #list("file03.txt", "file04.txt")),
#' #conditionsByGeneration = list(c(0,1), c(0,1)), output_dir = "test",
#' #nbrPermutations = 10000, nbrCores = 1,
#' #nbrCoresDiffMeth = 1, doingSites = TRUE, doingTiles = TRUE,
#' #minReads = 10, minMethDiff = 25, qvalue = 0.01, maxPercReads = 99.9,
#' #destrand = TRUE, minCovBasesForTiles = 10, tileSize = 1000,
#' #stepSize = 500, vSeed = 12)
#'
#' ## The function raises an error when at least one paramater is not valid
#' \dontrun{methylInheritance:::validateRunPermutationUsingRDS(
#' allFilesByGeneration = list(list("file01.txt", "file02.txt"),
#' list("file03.txt", "file04.txt")),
#' conditionsByGeneration = list(c(0,1)), output_dir = "test",
#' nbrPermutations = 10000, nbrCores = 1,
#' nbrCoresDiffMeth = 1, doingSites = TRUE, doingTiles = TRUE,
#' minReads = 10, minMethDiff = 25, qvalue = 0.01, maxPercReads = 99.9,
#' destrand = TRUE, minCovBasesForTiles = 10, tileSize = 1000,
#' stepSize = 500, vSeed = 12)}
#'
#' @author Astrid Deschenes
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
validateRunPermutationUsingMethylKitInfo <- function(methylKitInfo,
                                    type, outputDir,
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

    ## Validate that the output_dir is an not empty string
    if (!is.null(outputDir) && !is.character(outputDir)) {
        stop("output_dir must be a character string or NULL")
    }

    ## Validate that nbrPermutations is an positive integer
    if (!(isSingleInteger(nbrPermutations) ||
          isSingleNumber(nbrPermutations)) ||
        as.integer(nbrPermutations) < 1) {
        stop("nbrPermutations must be a positive integer or numeric")
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


#' @title Parameters validation for the \code{\link{runPermutation}} function
#'
#' @description Validation of all parameters needed by the public
#' #' \code{\link{runPermutation}} function.
#'
#' @param allFilesByGeneration a \code{list} of \code{list}, a \code{list}
#' composed of \code{list}
#' of files with methylation information for
#' bases or regions in the genome. One \code{list} must contain all files
#' related to the same generation. So, if 3 generations are analyzed, a
#' \code{list} containing 3 \code{list} must be passed. At least 2 generations
#' must be present to do a permutation analysis.
#' The parameter
#' corresponds to the \code{location} parameter in the  \code{methylKit} package.
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
#' corresponds to the \code{treatment} parameter in the  \code{methylKit} package.
#'
#' @param output_dir a string, the name of the directory that will contain
#' the results of the permutation. If the directory does not exist, it will
#' be created.
#'
#' @param genomeVersion a string, the genome assembly such as hg18, mm9, etc.
#' It can be any string. The parameter
#' correspond to the \code{assembly} parameter in the  \code{methylKit} package.
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
#' corresponds to the \code{num.cores} parameter in the  \code{methylKit} package.
#'
#' @param doingSites a \code{logical}, when \code{TRUE} will do the analysis
#' on the CpG dinucleotide sites.
#'
#' @param doingTiles a \code{logical}, when \code{TRUE} will do the analysis
#' on the tiles.
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
#' ## The function returns 0 when all paramaters are valid
#' #methylInheritance:::validateRunRealAnalysis(
#' #allFilesByGeneration = list(list("file01.txt", "file02.txt"),
#' #list("file03.txt", "file04.txt")),
#' #conditionsByGeneration = list(c(0,1), c(0,1)), output_dir = "test",
#' #genomeVersion = "hg19", nbrCores = 1,
#' #nbrCoresDiffMeth = 1, doingSites = TRUE, doingTiles = TRUE,
#' #minReads = 10, minMethDiff = 25, qvalue = 0.01, maxPercReads = 99.9,
#' #destrand = TRUE, minCovBasesForTiles = 10, tileSize = 1000,
#' #stepSize = 500, vSeed = 12)
#'
#' ## The function raises an error when at least one paramater is not valid
#' \dontrun{methylInheritance:::validateRunRealAnalysis(
#' allFilesByGeneration = list(list("file01.txt", "file02.txt"),
#' list("file03.txt", "file04.txt")),
#' conditionsByGeneration = list(c(0,1), c(0,1,1)), output_dir = "test",
#' genomeVersion = "hg19", nbrCores = 1,
#' nbrCoresDiffMeth = 1, doingSites = TRUE, doingTiles = TRUE,
#' minReads = 10, minMethDiff = 25, qvalue = 0.01, maxPercReads = 99.9,
#' destrand = TRUE, minCovBasesForTiles = 10, tileSize = 1000,
#' stepSize = 500, vSeed = 12)}
#'
#' @author Astrid Deschenes
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
validateRunRealAnalysis <- function(allFilesByGeneration,
                                    conditionsByGeneration,
                                    output_dir, genomeVersion,
                                    nbrCores, nbrCoresDiffMeth, doingSites,
                                    doingTiles,
                                    minReads, minMethDiff, qvalue,
                                    maxPercReads, destrand,
                                    minCovBasesForTiles, tileSize,
                                    stepSize, vSeed) {

    ## Validate that the allFilesByGeneration is a list
    if (!is.list(allFilesByGeneration)) {
        stop("allFilesByGeneration must be a list")
    }

    ## Validate that the allFilesByGeneration contains at least 2 lists
    if (length(allFilesByGeneration) < 2 ||
        !all(sapply(allFilesByGeneration, is.list))) {
        stop(paste0("allFilesByGeneration must be a list containing at ",
                    "least 2 sub-lists"))
    }

    ## Validate that all entries in allFilesByGeneration are existing files
    for (entry in unlist(allFilesByGeneration)) {
        if (!file.exists(entry)) {
            stop(paste0("The file \"", entry, "\" does not exist."))
        }
    }

    # ## Validate that the conditionsByGeneration is a list
    # if (!is.list(conditionsByGeneration)) {
    #     stop("conditionsByGeneration must be a list")
    # }
    #
    # ## Validate that the conditionsByGeneration contains at least 2 vectors
    # if (length(conditionsByGeneration) < 2 ||
    #     !all(sapply(conditionsByGeneration, is.vector))) {
    #     stop(paste0("conditionsByGeneration must be a list containing at ",
    #                 "least 2 sub-vectors"))
    # }
    #
    # ## Validate that the conditionsByGeneration and AllFilesByGeneration are
    # ## the same length
    # fileLengths <- sapply(allFilesByGeneration, length)
    # conditionLengths <- sapply(conditionsByGeneration, length)
    # if (length(fileLengths) != length(conditionsByGeneration)) {
    #     stop(paste0("The content of the parameters \"allFilesByGeneration\" ",
    #                 "and \"conditionsByGeneration\" must be the same length"))
    # }

    ## Validate that the output_dir is an not empty string
    if (!is.character(output_dir)) {
        stop("output_dir must be a character string")
    }

    ## Validate that the genomeVersion is an not empty string
    if (!is.character(genomeVersion)) {
        stop("genomeVersion must be a character string")
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

    ## Validate that doingSites is a logical
    if (!is.logical(doingSites)) {
        stop("doingSites must be a logical.")
    }

    ## Validate that doingTiles is a logical
    if (!is.logical(doingTiles)) {
        stop("doingTiles must be a logical.")
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
        stop("destrand must be a logical.")
    }

    if (doingTiles) {
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
    ## Validate that nbrPermutations is an positive integer
    if (!(isSingleInteger(vSeed) ||
          isSingleNumber(vSeed)) && as.integer(vSeed) != -1) {
        stop("vSeed must be either -1 or a positive integer or numeric")
    }

    return(0)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param directory a string of \code{character}, TODO
#'
#' @param elementType a sting of \code{character}, TODO
#'
#' @return TODO
#'
#' @examples
#'
#' ## TODO
#'
#' @author Astrid Deschenes
#' @importFrom utils read.table
#' @keywords internal
extractDataFromFile <- function(directory, elementType = c("SITES", "TILES")) {

    ## Set the file extension
    elementPattern <- ".perbase.txt"
    if (elementType == "TILES") {
        elementPattern <- "*.pertile.txt"
    }

    ## Initialize variables
    elements_per_generation <- list()
    elements_per_generation[[elementType]] <- list()

    ## List directories related to methyl diff files
    generationsDir <- list.files(path = paste0(directory, elementType, "/"),
                                    pattern = "Generation_*",
                                    all.files = FALSE,
                                    full.names = FALSE, recursive = FALSE,
                                    ignore.case = FALSE, include.dirs = TRUE,
                                    no.. = FALSE)

    ## Setting the number of generations according to the number of directories
    ## detected
    nbrGenerations <- length(generationsDir)

    ## Two types of files (hyper and hypo) per generation
    nbrExpectedFiles <- nbrGenerations * 2

    ## List all files related to methyl diff for all generations
    sitesFiles <- list.files(path = paste0(directory, elementType),
                                pattern = paste0("*", elementPattern),
                                all.files = FALSE,
                                full.names = FALSE, recursive = TRUE,
                                ignore.case = FALSE, include.dirs = FALSE,
                                no.. = FALSE)

    ## Find all permutations (each permutation has a unique number associated)
    ## that have generated the good number of methyl diff files
    ## for all generations
    id <- sapply(strsplit(sitesFiles, "_hyp"),
                function(x) return(as.integer(strsplit(x[1], "/")[[1]][2])))
    id_table <- table(id)
    id_tabel_subset <- id_table[id_table == nbrExpectedFiles]
    id_final <- as.numeric(names(id_tabel_subset))

    id_final_length <- length(id_final)

    ## Extract the number of conserved sites between each paire of 2 generations
    for (type in c("hypo", "hyper")) {
        genArray <- 1:nbrGenerations
        groupsTwo <- lapply(genArray[-length(genArray)],
                                function(g) return(c(g, g + 1)))
        for (groupTwo in groupsTwo) {
            generation_name <- paste0("Generation_", groupTwo[1], "_and_",
                                            groupTwo[2])
            elements_per_generation[[elementType]][[generation_name]] <- list()
            for (j in id_final){
                fileName01 <- paste0(directory, elementType, "/Generation_",
                                   groupTwo[1], "/", j, "_", type,
                                   elementPattern)
                fileName02 <- paste0(directory, elementType, "/Generation_",
                                     groupTwo[2], "/", j, "_", type,
                                     elementPattern)

                if (file.info(fileName01)$size > 0 &&
                        file.info(fileName02)$size > 0) {
                    sites01 <- read.table(fileName01, stringsAsFactors = F)$V4
                    sites02 <- read.table(fileName02, stringsAsFactors = F)$V4
                    results <- intersect(sites01, sites02)
                    elements_per_generation[[elementType]][[generation_name]][[type]][j] <- length(results)
                } else {
                    elements_per_generation[[elementType]][[generation_name]][[type]][j] <- 0
                }
            }
        }
    }

    return(elements_per_generation)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param methDiff, TODO
#'
#' @param pDiff, TODO
#'
#' @param qCut, TODO
#'
#' @param typeD, TODO
#'
#' @return TODO
#'
#' @examples
#'
#' ## TODO
#'
#' @author Pascal Belleau
#' @importFrom methylKit getMethylDiff
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @keywords internal
getGRangesFromMethylDiff <- function(methDiff, pDiff, qCut, typeD= "all"){

    methDiffK <- lapply(1:length(methDiff), FUN = function(i,
                                                           methDiff,
                                                           pDiff,
                                                           qCut, typeD){
        methK <- getMethylDiff(methDiff[[i]],difference=pDiff,qvalue=qCut, type=typeD)
        GRanges(seqnames = methK$chr, ranges = IRanges(start = methK$start,
                                                       end = methK$end),
                strand = methK$strand, pvalue = methK$pvalue,
                qvalue = methK$qvalue, meth.diff = methK$meth.diff)
    }, methDiff = methDiff, pDiff = pDiff, qCut = qCut, typeD=typeD)

    methDiffK
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
#' ## sum(width(res))
#' ## TODO
#'
#' @author Pascal Belleau
#' @importFrom GenomicRanges intersect GRanges
#' @importFrom S4Vectors DataFrame values
#' @keywords internal
interGeneration <- function(resultAllGenGR){

    lInter <- list("i2" = list(), "iAll" = list())

    lInter$i2 <- lapply(2:length(resultAllGenGR), FUN = function(i,b){
        upM <- intersect(b[[i-1]][b[[i-1]]$meth.diff > 0], b[[i]][b[[i]]$meth.diff > 0])
        downM <- intersect(b[[i-1]][b[[i-1]]$meth.diff < 0], b[[i]][b[[i]]$meth.diff < 0])
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
    lInter
}

#' @title TODO
#'
#' @description TODO
#'
#' @param output_dir a string of \code{character}, TODO
#'
#' @param nbrGenerations a positive \code{integer}, TODO
#'
#' @param doingSites a \code{logical}, TODO
#'
#' @param doingTiles a \code{logical}, TODO
#'
#' @return \code{0} when all directories are created without problem.
#'
#' @examples
#'
#' ## TODO
#'
#' @author Astrid Deschenes
#' @keywords internal
createOutputDir <- function(output_dir, nbrGenerations,
                                doingSites=TRUE, doingTiles=FALSE) {

    # Create directories for output files
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, showWarnings = TRUE)
    }

    if (doingSites) {
        type <-  "SITES"
        dirName = paste0(output_dir, type)
        if (!dir.exists(dirName)) {
            dir.create(dirName, showWarnings = TRUE)
        }
        for (j in 1:nbrGenerations) {
            dirName = paste0(output_dir, type, "/Generation_", j)
            if (!dir.exists(dirName)) {
                dir.create(dirName, showWarnings = TRUE)
            }
        }
    }

    if (doingTiles) {
        type <-  "TILES"
        dirName = paste0(output_dir, type)
        if (!dir.exists(dirName)) {
            dir.create(dirName, showWarnings = TRUE)
        }
        for (j in 1:nbrGenerations) {
            dirName = paste0(output_dir, type, "/Generation_", j)
            if (!dir.exists(dirName)) {
                dir.create(dirName, showWarnings = TRUE)
            }
        }
    }

    return(0)
}

#' @title Run one analysis using \code{methylKit} package
#'
#' @description Run one CpG site or region analysis using the \code{methylKit}
#' package. The output of the analysis is saved in a file in the specified
#' directory.
#'
#' @param info TODO
#'
#' @param output_dir a string, the name of the directory that will contain
#' the results of the permutation. If the directory does not exist, it will
#' be created.
#'
#' @param genomeVersion a string, the genome assembly such as hg18, mm9, etc.
#' It can be any string. The parameter
#' correspond to the \code{assembly} parameter in the \code{methylKit} package.
#'
#' @param minReads a positive \code{integer} Bases and regions having lower
#' coverage than this count are discarded. The parameter
#' correspond to the \code{lo.count} parameter in the  \code{methylKit} package.
#'
#' @param qvalue a positive \code{double} inferior ot \code{1}, the cutoff
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
#' @param nbrCoresDiffMeth a positive integer, the number of cores to use for
#' parallel differential methylation calculations.Parameter used for both
#' sites and tiles analysis. The parameter
#' corresponds to the \code{num.cores} parameter in the
#' package \code{methylKit}.
#' Default: \code{1} and always \code{1} for Windows.
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
#' @param debug Default: \code{FALSE}.
#'
#' @return TODO
#'
#' @examples
#'
#' ##TODO
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom methylKit filterByCoverage normalizeCoverage unite calculateDiffMeth getMethylDiff getData tileMethylCounts methRead
#' @keywords internal
runTOTO <- function(info, output_dir,
                              genomeVersion,
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
                              stepSize = 1000) {

    #print(info)

    nbrGenerations <- length(info)

    permutationList <- list()
    permutationList[["TILES"]] <- list()
    permutationList[["SITES"]] <- list()

    permutationList[["TILES"]][["HYPER"]] <- vector("list", length(nbrGenerations))
    permutationList[["TILES"]][["HYPO"]] <- vector("list", length(nbrGenerations))
    permutationList[["SITES"]][["HYPER"]] <- vector("list", length(nbrGenerations))
    permutationList[["SITES"]][["HYPO"]] <- vector("list", length(nbrGenerations))

    for (i in 1:nbrGenerations) {
        file.list    <- info[[i]][["file.list"]]
        sampleNames  <- info[[i]][["sampleNames"]]
        conditions   <- info[[i]][["conditions"]]
        designName   <- info[[i]][["designName"]]


        myobj <- methRead(location = file.list,
                  sample.id = as.list(sampleNames),
                  assembly = genomeVersion,
                  context = "CpG",
                  treatment = conditions)

        ## SITES
        if (doingSites) {

            filtered.sites <- filterByCoverage(myobj,
                                           lo.count = minReads,
                                           lo.perc = NULL,
                                           hi.count = NULL,
                                           hi.perc = maxPercReads)

            filtered.sites <- normalizeCoverage(filtered.sites, "median")

            meth.sites <- unite(filtered.sites, destrand = destrand)

            if (length(meth.sites@.Data[[1]]) == 0) {
                stop("meth.sites IS EMPTY")
            }

            ####################################
            ## Get diff methyl sites
            ####################################
            myDiff.sites <- calculateDiffMeth(meth.sites,
                                          num.cores = nbrCoresDiffMeth)

            permutationList[["SITES"]][["HYPER"]][[i]] <- getMethylDiff(myDiff.sites,
                                            difference = minMethDiff,
                                            qvalue = qvalue, type = "hyper")

            permutationList[["SITES"]][["HYPO"]][[i]]  <- getMethylDiff(myDiff.sites,
                                            difference = minMethDiff,
                                            qvalue = qvalue, type = "hypo")
        }

        ## TILES
        if (doingTiles) {
            tiles <- tileMethylCounts(myobj, win.size = tileSize,
                                  step.size = stepSize,
                                  cov.bases = minCovBasesForTiles)

            filtered.tiles <- filterByCoverage(tiles,
                                           lo.count = minReads,
                                           lo.perc = NULL,
                                           hi.count = NULL,
                                           hi.perc = maxPercReads)

            filtered.tiles <- normalizeCoverage(filtered.tiles, "median")

            meth.tiles <- unite(filtered.tiles, destrand = destrand)

            ####################################
            ## Get diff methyl tiles
            ####################################
            myDiff.tiles <- calculateDiffMeth(meth.tiles,
                                          num.cores = nbrCoresDiffMeth)

            permutationList[["TILES"]][["HYPER"]][[i]] <- getMethylDiff(myDiff.tiles,
                                            difference = minMethDiff,
                                            qvalue = qvalue, type = "hyper")

            permutationList[["TILES"]][["HYPO"]][[i]]  <- getMethylDiff(myDiff.tiles,
                                            difference = minMethDiff,
                                            qvalue = qvalue, type = "hypo")
        }
    }

    warnings()
    return(permutationList)
}

#' @title Run one permutation using \code{methylKit} package. One permutation
#' includes analysis for all generations associated to the same permutation.
#'
#' @description Run one CpG site or region analysis using the \code{methylKit}
#' package. The output of the analysis is saved in a file in the specified
#' directory.
#'
#' @param methylRawForAllGenerations TODO
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
#' @param qvalue a positive \code{double} inferior ot \code{1}, the cutoff
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
#' @return TODO
#'
#' @examples
#'
#' ##TODO
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom methylKit filterByCoverage normalizeCoverage unite calculateDiffMeth getMethylDiff getData tileMethylCounts methRead
#' @keywords internal
runOnePermutationOnAllGenerations <- function(methylRawForAllGenerations,
                        type = c("both", "sites", "tiles"),
                        outputDir = NULL,
                        nbrCoresDiffMeth = 1,
                        minReads = 10, minMethDiff = 10,
                        qvalue = 0.01, maxPercReads = 99.9,
                        destrand = FALSE, minCovBasesForTiles = 0,
                        tileSize = 1000, stepSize = 1000) {

    doTiles <- any(type %in% c("tiles", "both"))
    doSites <- any(type %in% c("sites", "both"))

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

    ## TODO
    ## add intersection function

    return(permutationList)
}
