#' @title TODO
#'
#' @description TODO
#'
#' @param info TODO
#'
#' @param output_di TODO
#'
#' @param doingSites Default: \code{TRUE}.
#'
#' @param doingTiles Default: \code{FALSE}.
#'
#' @param debug Default: \code{FALSE}.
#'
#' @return
#'
#' @examples
#'
#' ##TODO
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom methylKit read filterByCoverage normalizeCoverage unite calculateDiffMeth get.methylDiff getData tileMethylCounts
#' @keywords internal
runOnePermutation <- function(info, output_dir,
                              doingSites = TRUE,
                              doingTiles = FALSE, debug = FALSE) {

    designName = info$designName
    count = info$count
    file.list = info$file.list
    sampleNames = info$sampleNames
    genomeVersion = info$genomeVersion
    conditions = info$conditions
    mergeStrand = info$mergeStrand
    tileSize = info$tileSize
    stepSize = info$stepSize
    minCGs = info$minCGs
    minReads = info$minReads
    minMethDiff = info$minMethDiff

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
                                treatment = conditions,
                                context = "CpG")

    ## SITES
    if (doingSites) {
        filtered.sites <- filterByCoverage(myobj,
                                       lo.count = minReads,
                                       lo.perc = NULL,
                                       hi.count = NULL,
                                       hi.perc = 99.9)
        filtered.sites <- normalizeCoverage(filtered.sites, "median")

        meth.sites <- unite(filtered.sites, destrand = mergeStrand)

        ####################################
        ## Get diff methyl sites
        ####################################
        myDiff.sites <- calculateDiffMeth(meth.sites, num.cores=1)


        myDiff.sites.hyper <- get.methylDiff(myDiff.sites, difference = minMethDiff, qvalue = 0.01, type = "hyper")
        myDiff.sites.hypo  <- get.methylDiff(myDiff.sites, difference = minMethDiff, qvalue = 0.01, type = "hypo")
    }

    ## TILES
    if (doingTiles) {
        tiles <- tileMethylCounts(myobj, win.size = tileSize, step.size = stepSize, cov.bases=minCGs)
        filtered.tiles <- filterByCoverage(tiles, lo.count = minReads, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
        filtered.tiles <- normalizeCoverage(filtered.tiles, "median")
        meth.tiles <- unite(filtered.tiles, destrand = mergeStrand)


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
}
