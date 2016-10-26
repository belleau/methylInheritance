#' @title TODO
#'
#' @description TODO
#'
#' @param allDataByGeneration TODO
#'
#' @param infoByGeneration TODO
#'
#' @param nbCores Default = 1.
#'
#' @param nbrPermutations TODO
#'
#' @param output_dir TODO
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
#' @importFrom parallel mclapply
#' @importFrom utils flush.console write.table
#' @export
runPermutation <- function(allDataByGeneration, infoByGeneration, nbCores = 1,
                        nbrPermutations = 1000, output_dir, vSeed = -1) {

    nbrFiles <- sum(unlist(lapply(allDataByGeneration, length)))


    if (substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") {
        outpout_dir <- paste0(outpout_dir, "/")
    }

    if (vSeed == -1) {
        tSeed <- as.numeric(Sys.time())
        vSeed <- 1e8 * (tSeed - floor(tSeed))
    }
    set.seed(vSeed)

    nbGenerations <- length(allDataByGeneration)

    nbSamples  <- sum(length(unlist(allDataByGeneration)))
    allSamples <- unlist(allDataByGeneration)

    nbSamplesByGeneration <- lapply(allDataByGeneration, length)

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
            infoGeneration <- list(info=infoByGeneration[j],
                                   file.list=filesGeneration,
                                   designName=paste0("Generation_", j),
                                   count = i)
            permutationList <- c(permutationList, infoGeneration)
            start <- end + 1
        }

        finalList <- c(finalList, permutationList)
    }

    if (!dir.exists(output_dir)) {
        dir.create(output_dir, showWarnings = TRUE)
    }
    for (j in 1:nbGenerations) {
        dirName = paste0(output_dir, "/", "Generation_", j)
        if (!dir.exists(dirName)) {
            dir.create(dirName, showWarnings = TRUE)
        }
    }

    res <- mclapply(finalList, runPermutation, mc.preschedule = FALSE,
                    mc.cores = nbCores,
                    mc.set.seed = FALSE)

}
