#' @title TODO
#'
#' @description TODO
#'
#' @param info TODO
#'
#' @param output_di TODOr
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
runOnePermutation <- function(info,  output_dir,
                              doingSites = TRUE,
                              doingTiles = FALSE, debug = FALSE) {

    designName=info$designName
    count = info$count
    file.list = info$file.list
    sampleNames = info$sampleNames
    genomeVersion = info$info$genomeVersion
    conditions = info$info$conditions
    mergeStrand = info$info$mergeStrand
    tileSize = info$info$tileSize
    stepSize = info$info$stepSize
    minCGs = info$info$minCG
    minReads = info$info$minReads
    minMethDiff = info$info$minMethDiff

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

    filtered.myobj <- filterByCoverage(myobj,
                                       lo.count = minReads,
                                       lo.perc = NULL,
                                       hi.count = NULL,
                                       hi.perc = 99.9)
    filtered.myobj <- normalizeCoverage(filtered.myobj, "median")

    meth <- unite(filtered.myobj, destrand = mergeStrand)

    ## TILES
    tiles <- tileMethylCounts(myobj, win.size = tileSize, step.size = stepSize, cov.bases=minCGs)
    filtered.tiles <- filterByCoverage(tiles, lo.count = minReads, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
    filtered.tiles <- normalizeCoverage(filtered.tiles, "median")
    meth.tiles <- unite(filtered.tiles, destrand = mergeStrand)

    if (debug) {
        print(paste0(designName, " - ", count," - Differentially Methylation Analysis"))
        flush.console()
    }

    ####################################
    ## Get diff methyl sites and tiles
    ####################################
    myDiff <- calculateDiffMeth(meth, num.cores=2)
    myDiff.tiles <- calculateDiffMeth(meth.tiles, num.cores=2)

    myDiff20p.hyper <- get.methylDiff(myDiff, difference = minMethDiff, qvalue = 0.01, type = "hyper")
    myDiff20p.tiles.hyper <- get.methylDiff(myDiff.tiles, difference = minMethDiff, qvalue = 0.01, type = "hyper")
    myDiff20p.hypo <- get.methylDiff(myDiff, difference = minMethDiff, qvalue = 0.01, type = "hypo")
    myDiff20p.tiles.hypo <- get.methylDiff(myDiff.tiles, difference =  minMethDiff, qvalue = 0.01, type = "hypo")

    #####################
    ## Save data
    #####################
    if(nrow(myDiff20p.hypo) > 0) {
        w=cbind(getData(myDiff20p.hypo)[,c(1,2,3)],paste(getData(myDiff20p.hypo)[,1],".",getData(myDiff20p.hypo)[,2],".",getData(myDiff20p.hypo)[,3],sep=""),getData(myDiff20p.hypo)[,c(7,4,5,6)])
        colnames(w)[4]="dmr.id"
        write.table(w, paste0(output_dir, "/", designName, "/", count, "_hypo.perbase.20.percent",sep=""), quote=F, row.names=F, col.names=F,sep="\t")
    }
    if(nrow(myDiff20p.tiles.hypo)>0) {
        w=NULL
        w=cbind(getData(myDiff20p.tiles.hypo)[,c(1,2,3)],paste(getData(myDiff20p.tiles.hypo)[,1],".",getData(myDiff20p.tiles.hypo)[,2],".",getData(myDiff20p.tiles.hypo)[,3],sep=""),getData(myDiff20p.tiles.hypo)[,c(7,4,5,6)])
        colnames(w)[4]="dmr.id"
        write.table(w, paste0(output_dir, "/", designName, "/", count, "_hypo.pertile.20.percent",sep=""), quote=F, row.names=F, col.names=F,sep="\t")
    }
    if(nrow(myDiff20p.hyper)>0) {
        w=NULL
        w=cbind(getData(myDiff20p.hyper)[,c(1,2,3)],paste(getData(myDiff20p.hyper)[,1],".",getData(myDiff20p.hyper)[,2],".",getData(myDiff20p.hyper)[,3],sep=""),getData(myDiff20p.hyper)[,c(7,4,5,6)])
        colnames(w)[4]="dmr.id"
        write.table(w, paste0(output_dir, "/", designName, "/", count, "_hyper.perbase.20.percent",sep=""), quote=F, row.names=F, col.names=F,sep="\t")
    }
    if(nrow(myDiff20p.tiles.hyper)>0) {
        w=NULL
        w=cbind(getData(myDiff20p.tiles.hyper)[,c(1,2,3)],paste(getData(myDiff20p.tiles.hyper)[,1],".",getData(myDiff20p.tiles.hyper)[,2],".",getData(myDiff20p.tiles.hyper)[,3],sep=""),getData(myDiff20p.tiles.hyper)[,c(7,4,5,6)])
        colnames(w)[4]="dmr.id"
        write.table(w, paste0(output_dir, "/", designName, "/", count, "_hyper.pertile.20.percent",sep=""), quote=F, row.names=F, col.names=F,sep="\t")
    }
}
