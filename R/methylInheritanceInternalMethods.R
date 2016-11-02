#' @title TODO
#'
#' @description TODO
#'
#' @param info TODO
#'
#' @param output_dir a string, the name of the directory that will contain
#' the results of the permutation. If the directory does not exist, it will
#' be created.
#'
#' @param genomeVersion a string, the genome assembly such as hg18, mm9, etc.
#' It can be any string.
#'
#' @param minReads TODO
#'
#' @param minCGs a positive integer, the minimum number of bases to
#' be covered in a given tiling window
#'
#' @param tileSize a positive integer, the size of the tiling window
#'
#' @param stepSize a positive integer, the step size of tiling windows
#'
#' @param minMethDiff a positive integer, the step size of tiling windows
#'
#' @param destrand a logical, when \code{TRUE} will merge reads on both
#' strands of a CpG dinucleotide to provide better coverage. Only advised
#' when looking at CpG methylation.
#' Default: \code{FALSE}.
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
#' @importFrom methylKit read filterByCoverage normalizeCoverage unite calculateDiffMeth get.methylDiff getData tileMethylCounts
#' @keywords internal
runOnePermutation <- function(info, output_dir, genomeVersion, minReads,
                              minCGs, tileSize, stepSize, minMethDiff,
                              destrand = FALSE,
                              doingSites = TRUE,
                              doingTiles = FALSE, debug = FALSE) {

    print(info)

    print("HI")
    designName = info$designName
    count = info$count
    file.list = info$file.list
    sampleNames = info$sampleNames
    conditions = info$conditions

    print(paste0("designName ", designName))
    print(paste0("count ", count))
    print("file.list ")
    print(file.list)
    print(paste0("sampleNames "))
    print(sampleNames)
    print(paste0("conditions "))
    print(conditions)
    print("genomeVersion ")
    print(genomeVersion)
    print(paste0("tileSize ", tileSize))
    print(paste0("stepSize ", stepSize))
    print(paste0("minCGs ", minCGs))
    print(paste0("minReads ", minReads))
    print(paste0("minMethDiff ", minMethDiff))

    ####################################
    ## prepare data
    ####################################

    if (debug) {
        print(paste0(designName, " - ", count))
        flush.console()
    }

    myobj <- methylKit::read(location = file.list,
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
                                       hi.perc = 99.9)
        filtered.sites <- normalizeCoverage(filtered.sites, "median")

        meth.sites <- unite(filtered.sites, destrand = destrand)

        ####################################
        ## Get diff methyl sites
        ####################################
        myDiff.sites <- calculateDiffMeth(meth.sites, num.cores=1)


        myDiff.sites.hyper <- get.methylDiff(myDiff.sites, difference = minMethDiff, qvalue = 0.01, type = "hyper")
        myDiff.sites.hypo  <- get.methylDiff(myDiff.sites, difference = minMethDiff, qvalue = 0.01, type = "hypo")
    }

    ## TILES
    if (doingTiles) {
        tiles <- tileMethylCounts(myobj, win.size = tileSize, stepSize = stepSize, cov.bases=minCGs)
        filtered.tiles <- filterByCoverage(tiles, lo.count = minReads, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
        filtered.tiles <- normalizeCoverage(filtered.tiles, "median")
        meth.tiles <- unite(filtered.tiles, destrand = destrand)


        ####################################
        ## Get diff methyl tiles
        ####################################
        myDiff.tiles <- calculateDiffMeth(meth.tiles, num.cores=1)

        myDiff.tiles.hyper <- get.methylDiff(myDiff.tiles, difference = minMethDiff, qvalue = 0.01, type = "hyper")
        myDiff.tiles.hypo  <- get.methylDiff(myDiff.tiles, difference =  minMethDiff, qvalue = 0.01, type = "hypo")
    }

    if (debug) {
        print(paste0(designName, " - ", count," - Differentially Methylation Analysis"))
        flush.console()
    }

    #####################
    ## Save data
    #####################
    if(nrow(myDiff.sites.hypo) > 0) {
        sitesHypo <- cbind(getData(myDiff.sites.hypo)[,c(1,2,3)],
                paste0(getData(myDiff.sites.hypo)[,1], ".",
                       getData(myDiff.sites.hypo)[,2], ".",
                       getData(myDiff.sites.hypo)[,3]),
                getData(myDiff.sites.hypo)[,c(7,4,5,6)])
        colnames(sitesHypo)[4]="dmr.id"
        write.table(sitesHypo, paste0(output_dir, "/SITES/", designName, "/", count, "_hypo.perbase.txt",sep=""), quote=F, row.names=F, col.names=F,sep="\t")
    }

    if(nrow(myDiff.sites.hyper)>0) {
        sitesHyper <- cbind(getData(myDiff.sites.hyper)[,c(1,2,3)],
                            paste0(getData(myDiff.sites.hyper)[,1], ".",
                                   getData(myDiff.sites.hyper)[,2], ".",
                                   getData(myDiff.sites.hyper)[,3]),
                            getData(myDiff.sites.hyper)[,c(7,4,5,6)])
        colnames(sitesHyper)[4]="dmr.id"
        write.table(sitesHyper, paste0(output_dir, "/SITES/", designName, "/", count, "_hyper.perbase.txt"), quote=F, row.names=F, col.names=F,sep="\t")
    }

    if(nrow(myDiff.tiles.hypo)>0) {
        w=NULL
        w=cbind(getData(myDiff.tiles.hypo)[,c(1,2,3)],
                paste0(getData(myDiff.tiles.hypo)[,1], ".",
                       getData(myDiff.tiles.hypo)[,2], ".",
                       getData(myDiff.tiles.hypo)[,3]),
                getData(myDiff.tiles.hypo)[,c(7,4,5,6)])
        colnames(w)[4]="dmr.id"
        write.table(w, paste0(output_dir, "/TILES/", designName, "/", count, "_hypo.pertile.txt",sep=""), quote=F, row.names=F, col.names=F,sep="\t")
    }

    if(nrow(myDiff.tiles.hyper)>0) {
        w=NULL
        w=cbind(getData(myDiff.tiles.hyper)[,c(1,2,3)],paste0(getData(myDiff.tiles.hyper)[,1], ".",
                                                             getData(myDiff.tiles.hyper)[,2],".",
                                                             getData(myDiff.tiles.hyper)[,3]),getData(myDiff.tiles.hyper)[,c(7,4,5,6)])
        colnames(w)[4]="dmr.id"
        write.table(w, paste0(output_dir, "/TILES/", designName, "/", count, "_hyper.pertile.txt"), quote=F, row.names=F, col.names=F,sep="\t")
    }

    return(0)
}


#' @title Extract sample name from file name
#'
#' @description Extract sample name from file name
#'
#' @param fileName TODO
#'
#' @return TODO
#'
#' @examples
#'
#' ##TODO
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @keywords internal
getSampleNameFromFileName <- function(fileName) {
    results <- strsplit(fileName, split="/")[[1]]
    return(results[length(results)])
}

#' @title Extract sample name from file name
#'
#' @description Extract sample name from file name
#'
#' @param info TODO
#'
#' @param output_dir TODO
#'
#' @param genomeVersion TODO
#'
#' @param minReads TODO
#'
#' @param doingSites TODO
#'
#' @param doingTiles TODO
#'
#' @param debug TODO
#'
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @examples
#'
#' ##TODO
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @keywords internal
validateRunOnePermutation <- function(info, output_dir, genomeVersion,
                                        minReads, minCGs, tileSize,
                                        stepSize, minMethDiff,
                                        destrand,
                                        doingSites, doingTiles, debug) {

    if (!is.list(info)) {
        stop("info must be a list")
    }

    if (is.null(info$count)) {
        stop("info must be a list with an entry called count")
    }

    if (is.null(info$sampleNames)) {
        stop("info must be a list with an entry called sampleNames")
    }

    if (is.null(info$stepSize)) {
        stop("info must be a list with an entry called stepSize")
    }

    if (is.null(info$tileSize)) {
        stop("info must be a list with an entry called tileSize")
    }

    if (is.null(info$file.list)) {
        stop("info must be a list with an entry called file.list")
    }

    if (is.null(info$conditions)) {
        stop("info must be a list with an entry called conditions")
    }

    ## Validate that the genomeVersion is an not empty string
    if (!is.character(genomeVersion)) {
        stop("genomeVersion must be a character string")
    }

    ## Validate that the output_dir directory exist
    charTest <- substr(output_dir, nchar(output_dir), nchar(output_dir))
    if (charTest == "/") {
        output_dir <- substr(output_dir, 1 ,nchar(output_dir)-1)
    }

    if (!file.exists(output_dir)) {
        stop("The directory \'", output_dir, "\' does not exist.")
    }

    if (!(isSingleInteger(minReads) || isSingleNumber(minReads)) ||
        as.integer(minReads) < 1) {
        stop("minReads must be a positive integer or numeric")
    }

    ## Validate that doingSites is a logical
    if (!is.logical(doingSites)) {
        stop("doingSites must be a logical.")
    }

    ## Validate that doingTiles is a logical
    if (!is.logical(doingTiles)) {
        stop("doingTiles must be a logical.")
    }

    ## Validate that debug is a logical
    if (!is.logical(debug)) {
        stop("debug must be a logical.")
    }

    return(0)
}

#' @title Extract sample name from file name
#'
#' @description Extract sample name from file name
#'
#' @param allFilesByGeneration TODO
#'
#' @param conditionsByGeneration TODO
#'
#' @param nbrCores a positive integer, the number of cores to use when
#' processing the analysis. Always \code{1} for Windows.
#'
#' @param minCGs a positive integer, the minimum number of bases to
#' be covered in a given tiling window.
#'
#' @param tileSize a positive integer, the size of the tiling window.
#'
#' @param stepSize a positive integer, the step size of tiling windows.
#'
#' @param minMethDiff a positive integer, the step size of tiling windows.
#'
#' @param destrand a logical, when \code{TRUE} will merge reads on both
#' strands of a CpG dinucleotide to provide better coverage. Only advised
#' when looking at CpG methylation.
#'
#' @param doingSites a logical, when \code{TRUE} will do the analysis on the
#' CpG dinucleotide sites.
#'
#' @param doingTiles a logical, when \code{TRUE} will do the analysis on the
#' tiles.
#'
#' @param vSeed a positive integer, either \code{-1} or a positive integer.
#'
#'
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @examples
#'
#' ##TODO
#'
#' @author Astrid Deschenes
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
validateRunPermutation <- function(allFilesByGeneration,
                                    conditionsByGeneration,
                                    nbrCores, minReads, minCGs, tileSize,
                                    stepSize, minMethDiff,
                                    destrand, genomeVersion,
                                    nbrPermutations, output_dir,
                                    doingTiles,
                                    doingSites, vSeed) {


    ## Validate that nbrCores is an positive integer
    if (!(isSingleInteger(nbrCores) || isSingleNumber(nbrCores)) ||
        as.integer(nbrCores) < 1) {
        stop("nbrCores must be a positive integer or numeric")
    }

    ## Validate that nbrCores is set to 1 on Windows system
    if (Sys.info()["sysname"] == "Windows" && as.integer(nbrCores) != 1) {
        stop("nbrCores must be 1 on a Windows system.")
    }

    ## Validate that minReads is an positive integer
    if (!(isSingleInteger(minReads) || isSingleNumber(minReads)) ||
        as.integer(minReads) < 1) {
        stop("minReads must be a positive integer or numeric")
    }

    ## Validate that destrand is a logical
    if (!is.logical(destrand)) {
        stop("destrand must be a logical.")
    }

    ## Validate that the genomeVersion is an not empty string
    if (!is.character(genomeVersion)) {
        stop("genomeVersion must be a character string")
    }

    ## Validate that nbrPermutations is an positive integer
    if (!(isSingleInteger(nbrPermutations) ||
            isSingleNumber(nbrPermutations)) ||
            as.integer(nbrPermutations) < 1) {
        stop("nbrPermutations must be a positive integer or numeric")
    }

    ## Validate that the output_dir is an not empty string
    if (!is.character(output_dir)) {
        stop("output_dir must be a character string")
    }

    ## Validate that doingSites is a logical
    if (!is.logical(doingSites)) {
        stop("doingSites must be a logical.")
    }

    ## Validate that doingTiles is a logical
    if (!is.logical(doingTiles)) {
        stop("doingTiles must be a logical.")
    }

    ## Validate that nbrPermutations is an positive integer
    if (!(isSingleInteger(vSeed) ||
          isSingleNumber(vSeed)) && as.integer(vSeed) != -1) {
        stop("vSeed must be either -1 or a positive integer or numeric")
    }

    return(0)
}
