#' @title TODO
#'
#' @description TODO
#'
#' @param allFilesByGeneration TODO
#'
#' @param conditionsByGeneration TODO
#'
#' @param nbCores Default = 1.
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
#' @param mergeStrand TODO
#'
#' @param genomeVersion TODO
#'
#' @param nbrPermutations TODO
#'
#' @param output_dir TODO
#'
#' @param doingTiles TODO
#'
#' @param doingSites TODO
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
                           nbCores = 1,
                        minReads, minCGs, tileSize, stepSize, minMethDiff,
                        mergeStrand,
                        genomeVersion,
                        nbrPermutations = 1000, output_dir, doingTiles,
                        doingSites, vSeed = -1) {

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

    nbSamplesByGeneration <- lapply(allFilesByGeneration, length)

    finalList <- list()

    for (i in 1:nbrPermutations) {

        # Randomnly mixt all samples
        samples <- sample(allSamples, size=nbSamples, replace = FALSE)

        permutationList <- list()
        start = 1
        for (j in 1:nbGenerations) {
            nbFiles <-  nbSamplesByGeneration[j]
            end = start + nbrFiles - 1
            filesGeneration <- samples[start:end]
            sampleNames <- sapply(filesGeneration, getSampleNameFromFileName)
            infoGeneration <- list(conditions=conditionsByGeneration[[j]],
                                   file.list=filesGeneration,
                                   designName = paste0("Generation_", j),
                                   sampleNames = sampleNames,
                                   minReads = minReads,
                                   minCGs = minCGs,
                                   tileSize = tileSize,
                                   stepSize = stepSize,
                                   minMethDiff = minMethDiff,
                                   mergeStrand = mergeStrand,
                                   genomeVersion = genomeVersion,
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

    if (nbCores == 1) {
        bp_param <- SnowParam()
    } else {
        bp_param <- MulticoreParam(workers = nbCores)
    }

    res <- bptry(bplapply(finalList,
                          FUN = runOnePermutation,
                          output_dir = output_dir,
                          BPPARAM = bp_param))

    if (!all(bpok(res))) {
        print("ALLO")
        print(res)
        positions <- grepl("stop", attr(res[[1]], "traceback"))
        position <- which(positions)[1] # First error message
        msg <- strsplit(attr(res[[1]], "traceback")[position], "[()]")[[1]][2]
        msg <- substring(msg, 2, nchar(msg) - 1)

        #msg <- lapply(res, function(x) {tail(attr(x, "traceback"))})
        stop(msg)
    }

    res
}



