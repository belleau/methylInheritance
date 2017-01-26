#' @rdname methylInheritanceAllResults
#'
#' @title  Print a methylInheritanceAllResults object
#'
#' @method print methylInheritanceAllResults
#'
#' @description Print a methylInheritanceAllResults object
#'
#' @param x the output object from \code{mergePermutationAndObservation}
#' function, \code{runPermutationUsingRDSFile} function (when
#' \code{runObservationAnalysis} = \code{TRUE} and
#' \code{runPermutationUsingMethylKitInfo} function (when
#' \code{runObservationAnalysis} = \code{TRUE} to be printed
#'
#' @param \ldots arguments passed to or from other methods
#'
#' @export
print.methylInheritanceAllResults <- function(x, ...) {

    nbGenerations <- 0
    if (!is.null(x$OBSERVATION$SITES)) {
        nbGenerations = length(x$OBSERVATION$SITES$i2) + 1
    } else if (!is.null(x$OBSERVATION$SITES)) {
        nbGenerations = length(x$OBSERVATION$TILES$i2) + 1
    }

    nbPermutations <- 0
    if (!is.null(x$PERMUTATION)) {
        nbPermutations <- length(x$PERMUTATION)
    }

    isSites <- FALSE

    if (!is.null(x$OBSERVATION$SITES)) {
        tt <- unlist(x$OBSERVATION$SITES)
        tt.names <- sapply(names(tt), function(x) {strsplit(x, "[.]")})
        tt.analysis <- sapply(tt.names, function(x) {return(x[[1]])})
        tt.types <- sapply(tt.names, function(x) {return(x[[2]])})

        result <- data.frame(SOURCE=rep("OBSERVATION", length(tt)),
                        ELEMENT = rep("SITES", length(tt)),
                        ANALYSIS = tt.analysis, TYPE = tt.types,
                        RESULT=tt, stringsAsFactors = FALSE, row.names = NULL)
        isSites <- TRUE
    }

    if (!is.null(x$OBSERVATION$TILES)) {
        tt <- unlist(x$OBSERVATION$TILES)
        tt.names <- sapply(names(tt), function(x) {strsplit(x, "[.]")})
        tt.analysis <- sapply(tt.names, function(x) {return(x[[1]])})
        tt.types <- sapply(tt.names, function(x) {return(x[[2]])})

        dataTiles <- data.frame(SOURCE=rep("OBSERVATION", length(tt)),
                                ELEMENT = rep("TILES", length(tt)),
                                ANALYSIS = tt.analysis, TYPE = tt.types,
                                RESULT=tt, stringsAsFactors = FALSE,
                                row.names = NULL)

        if (isSites) {
            result <- rbind(result, dataTiles)
        } else {
            result <- dataTiles
        }
    }

    # Print title before printing the content of the object
    cat("Permutation Analysis\n\n")
    cat("Number of Generations: ", nbGenerations, "\n")
    cat("Number of Permutations: ", nbPermutations , "\n\n")
    cat("Observation Results: \n")
    print.data.frame(result)
    invisible(x)
}
