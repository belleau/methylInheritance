#' @title TODO
#'
#' @description TODO
#'
#' @param allFilesByGeneration TODO
#'
#' @param conditionsByGeneration TODO
#'
#' @param nbrCores a positive integer, the number of cores to use when
#' processing the analysis. Default: \code{1} and always \code{1} for Windows.
#'
#' @param minReads TODO
#'
#' @param minCGs TODO
#'
#' @param tileSize TODO
#'
#' @param stepSize TODO
#'
#' @param minMethDiff TODO
#'
#' @param destrand TODO
#'
#' @param genomeVersion TODO
#'
#' @param nbrPermutations TODO
#'
#' @param output_dir TODO
#'
#'  @param doingSites a logical, when \code{TRUE} will do the analysis on the
#' CpG dinucleotide sites.
#'
#' @param doingTiles a logical, when \code{TRUE} will do the analysis on the
#' tiles.
#'
#' @param vSeed TODO
#'
#'
#' @return TODO
#'
#' @examples
#'
#' ##TODO
#'
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom BiocParallel bplapply MulticoreParam SnowParam bptry bpok
#' @importFrom utils flush.console write.table
#' @export
runPermutation <- function(allFilesByGeneration, conditionsByGeneration,
                            nbrCores = 1,
                            minReads = 10,
                            minCGs, tileSize, stepSize, minMethDiff,
                            destrand,
                            genomeVersion,
                            nbrPermutations = 1000, output_dir, doingTiles,
                            doingSites, vSeed = -1) {

    validateRunPermutation(allFilesByGeneration, conditionsByGeneration,
                           nbrCores, minReads,minCGs, tileSize, stepSize,
                           minMethDiff, destrand, genomeVersion,
                           nbrPermutations, output_dir, doingTiles,
                           doingSites, vSeed)

    nbrFiles <- sum(unlist(lapply(allFilesByGeneration, length)))

    ## Add last slash to path when absent
    if (substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") {
        output_dir <- paste0(output_dir, "/")
    }

    if (vSeed == -1) {
        tSeed <- as.numeric(Sys.time())
        vSeed <- 1e8 * (tSeed - floor(tSeed))
    }
    set.seed(vSeed)

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
            print(is.list(filesGeneration))
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
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, showWarnings = TRUE)
    }
    for (type in c("TILES", "SITES")) {
        dirName = paste0(output_dir, type)
        if (!dir.exists(dirName)) {
            dir.create(dirName, showWarnings = TRUE)
        }
        for (j in 1:nbGenerations) {
          dirName = paste0(output_dir, type, "/Generation_", j)
          if (!dir.exists(dirName)) {
              dir.create(dirName, showWarnings = TRUE)
          }
        }
    }

    if (nbrCores == 1) {
        bp_param <- SnowParam()
    } else {
        bp_param <- MulticoreParam(workers = nbrCores)
    }

    result <- runOnePermutation(finalList[[1]], output_dir = output_dir,
                      genomeVersion = genomeVersion, minReads = minReads,
                      minCGs, tileSize,
                      stepSize, minMethDiff,
                      destrand)
    result <- bptry(bplapply(finalList,
                           FUN = runOnePermutation,
                           output_dir = output_dir,
                           genomeVersion = genomeVersion,
                           minReads = minReads,
                           minCGs =minCGs,
                           tileSize = tileSize,
                           stepSize = stepSize,
                           minMethDiff = minMethDiff,
                           destrand = destrand,
                           BPPARAM = bp_param))


    # if (!all(bpok(result))) {
    #     a <- which(!bpok(result))
    #     print(a)
    #     print(result)
    #     print("e ICI")
    #     print(tail(attr(result[[a[1]]], "traceback")))
    #     print("e 2")
    #     print(tail(attr(result[[a[2]]], "traceback")))
    #     print("e 3")
    #     print(tail(attr(result[[a[3]]], "traceback")))
    #     print("e 4")
    #     print(tail(attr(result[[a[4]]], "traceback")))
    #     print("ALLO")
    #     positions <- grepl("stop", attr(result[[1]], "traceback"))
    #     position <- which(positions)[1] # First error message
    #     msg <- strsplit(attr(result[[1]], "traceback")[position], "[()]")[[1]][2]
    #     msg <- substring(msg, 2, nchar(msg) - 1)
    #
    #     #msg <- lapply(res, function(x) {tail(attr(x, "traceback"))})
    #     print(msg)
    # }
    # print("TOTO")
    result[which(bpok(result))]
}



