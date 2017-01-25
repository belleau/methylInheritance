###################################################
# Created by Astrid Deschenes
# 2017-01-05
###################################################

###################################################
## Test the methylInheritanceMethods functions
###################################################

DIRECTORY <- system.file("extdata", package = "methylInheritance")

METHYL_OBJ_FILE_01 <- dir(system.file("extdata", package = "methylInheritance"),
                          pattern = "methylObj_001.RDS", full.names = TRUE)

METHYL_OBJ_01 <- readRDS(METHYL_OBJ_FILE_01)

TEST_DIR <- dir(system.file("extdata", package = "methylInheritance"),
                        pattern = "TEST", full.names = TRUE)

data("methylInheritanceResults")

###########################################################
## runPermutationUsingRDSFile() function
###########################################################

## Test when methylKitRDSFile is not a valid RDS file name
test.runPermutationUsingRDSFile_methylKitRDSFile_not_valid <- function() {
    obs <- tryCatch(runPermutationUsingRDSFile(
        methylKitRDSFile = "HI",  outputDir = NULL,
        nbrPermutations = 2, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = FALSE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "The file \"HI\" does not exist."

    message <- paste0(" test.runPermutationUsingRDSFile_methylKitRDSFile_not_valid() ",
                      "- Not valid file for methylKitRDSFile did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when all parameters valid
test.runPermutationUsingRDSFile_good_001 <- function() {
}


###########################################################
# runObservationUsingRDSFile() function
###########################################################

## Test when methylKitRDSFile is not a valid RDS file name
test.runObservationUsingRDSFile_methylKitRDSFile_not_valid <- function() {
    obs <- tryCatch(runObservationUsingRDSFile(
        methylKitRDSFile = "ALLO",  outputDir = NULL,
        nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = FALSE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "The file \"ALLO\" does not exist."

    message <- paste0(" test.runObservationUsingRDSFile_methylKitRDSFile_not_valid() ",
                      "- Not valid file for methylKitRDSFile did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when all parameters valid
test.runObservationUsingRDSFile_good_001 <- function() {
    obs <- tryCatch(runObservationUsingRDSFile(
        methylKitRDSFile = METHYL_OBJ_FILE_01, type = "sites",
        outputDir = NULL, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 5, qvalue = 0.05,
        maxPercReads = 99.9, destrand = FALSE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 200),
        error=conditionMessage)

    exp <- list()
    exp[["OBSERVATION"]] <- list()
    exp[["OBSERVATION"]][["SITES"]] <- list()
    exp[["OBSERVATION"]][["SITES"]][["i2"]] <- list()
    exp[["OBSERVATION"]][["SITES"]][["i2"]][["HYPER"]] <- list(0,0)
    exp[["OBSERVATION"]][["SITES"]][["i2"]][["HYPO"]] <- list(3,3)
    exp[["OBSERVATION"]][["SITES"]][["iAll"]] <- list()
    exp[["OBSERVATION"]][["SITES"]][["iAll"]][["HYPER"]] <- list(0)
    exp[["OBSERVATION"]][["SITES"]][["iAll"]][["HYPO"]] <- list(2)

    message <- paste0(" test.runObservationUsingRDSFile_good_001() ",
                      "- All valid parameters did not generated expected result.")

    checkEquals(obs, exp, msg = message)
}


###########################################################
# extractInfo() function
###########################################################

# Test result when all parameters are good
test.extractInfo_good_01 <- function() {
    obs <- tryCatch(extractInfo(allResults = methylInheritanceResults,
                    type = "sites", inter="i2", 1),
                    error=conditionMessage)

    exp <- data.frame(TYPE = rep(c("HYPO","HYPER"), 21),
                      RESULT = c(2,4,2,4,4,3,1,5,3,3,4,2,0,0,0,1,1,0,6,2,2,5,1,
                            3,2,4,222,67,6,4,183,53,1,6,34,102,2,2,4,3,2,2),
                      SOURCE = c("OBSERVATION", "OBSERVATION",
                                    rep("PERMUTATION", 40)))

    message <- paste0(" test.extractInfo_good_01() ",
                      "- Valid parameters for formatForGraph did not generated expected results.")

    checkEquals(obs, exp, msg = message)
}


###########################################################
## loadAllRDSResults() function
###########################################################

## Test result when all parameters are good
# test.loadAllRDSResults_good_01 <- function() {
#
#     obs <- tryCatch(loadAllRDSResults(analysisResultsDir = TEST_DIR,
#                         permutationResultsDir = TEST_DIR, doingSites = TRUE,
#                         doingTiles = TRUE),
#                     error=conditionMessage)
#
#     exp <- list()
#
#     exp[["OBSERVATION"]] <- list()
#     exp[["OBSERVATION"]][["SITES"]] <- list()
#     exp[["OBSERVATION"]][["SITES"]][["i2"]] <- list(HYPER=list(21, 10), HYPO=list(15, 12))
#     exp[["OBSERVATION"]][["SITES"]][["iAll"]] <- list(HYPER=list(1), HYPO=list(3))
#     i2 <- list(HYPER=list(2000, 3000), HYPO=list(2000, 3000))
#     iAll <- list(HYPER=list(0), HYPO=list(0))
#     exp[["OBSERVATION"]][["TILES"]] <- list()
#     exp[["OBSERVATION"]][["TILES"]][["i2"]] <- list(HYPER=list(2000, 3000), HYPO=list(2000, 3000))
#     exp[["OBSERVATION"]][["TILES"]][["iAll"]] <- list(HYPER=list(0), HYPO=list(0))
#
#     exp[["PERMUTATION"]] <- list()
#     cas_01 <- list()
#     cas_01[["i2"]] <- list(HYPER=list(5, 7), HYPO=list(10, 11))
#     cas_01[["iAll"]] <- list(HYPER=list(2), HYPO=list(3))
#     cas_02 <- list()
#     cas_02[["i2"]] <- list(HYPER=list(8, 9), HYPO=list(4, 7))
#     cas_02[["iAll"]] <- list(HYPER=list(3), HYPO=list(0))
#     cas_03 <- list()
#     cas_03[["i2"]] <- list(HYPER=list(10, 7), HYPO=list(11, 7))
#     cas_03[["iAll"]] <- list(HYPER=list(1), HYPO=list(2))
#     exp[["PERMUTATION"]][["SITES"]] <- list(cas_01, cas_02, cas_03)
#     cas_01 <- list()
#     cas_01[["i2"]] <- list(HYPER=list(0, 0), HYPO=list(1000, 3000))
#     cas_01[["iAll"]] <- list(HYPER=list(0), HYPO=list(1000))
#     cas_02 <- list()
#     cas_02[["i2"]] <- list(HYPER=list(1000, 1000), HYPO=list(1000, 1000))
#     cas_02[["iAll"]] <- list(HYPER=list(0), HYPO=list(0))
#     cas_03 <- list()
#     cas_03[["i2"]] <- list(HYPER=list(0, 3000), HYPO=list(2000, 0))
#     cas_03[["iAll"]] <- list(HYPER=list(0), HYPO=list(0))
#     exp[["PERMUTATION"]][["TILES"]] <- list(cas_01, cas_02, cas_03)
#
#     class(exp) <- "methylInheritanceAllResults"
#
#     message <- paste0(" test.loadAllRDSResults_good_01() ",
#                         "- Valid parameters for loadAllRDSResults() ",
#                         "did not generated expected results.")
#
#     checkEquals(obs, exp, msg = message)
# }


###########################################################
## mergePermutationAndObservation() function
###########################################################

## Test when observationResults is not a list
test.mergePermutationAndObservation_observation_not_list <- function() {
    perm <- list()
    perm[["PERMUTATION"]] <- methylInheritanceResults$PERMUTATION

    obs <- tryCatch(mergePermutationAndObservation(permutationResults = perm,
                                observationResults = "33"),
                    error=conditionMessage)

    exp <- "observationResults must be a list"

    message <- paste0(" test.mergePermutationAndObservation_observation_not_list() ",
                      "- Not a list for observationResults did not generated expected results.")

    checkEquals(obs, exp, msg = message)
}

## Test when permutationResults is not a list
test.mergePermutationAndObservation_permutation_not_list <- function() {
    res <- list()
    res[["OBSERVATION"]] <- methylInheritanceResults$OBSERVATION

    obs <- tryCatch(mergePermutationAndObservation(permutationResults = "allo",
                    observationResults = res),
                    error=conditionMessage)

    exp <- "permutationResults must be a list"

    message <- paste0(" test.mergePermutationAndObservation_permutation_not_list() ",
                      "- Not a list for permutationResults did not generated expected results.")

    checkEquals(obs, exp, msg = message)
}

## Test result when all parameters are good
test.mergePermutationAndObservation_good_01 <- function() {

    observed <- list()
    observed[["OBSERVATION"]] <- list()
    observed[["OBSERVATION"]][["SITES"]] <- list()
    observed[["OBSERVATION"]][["SITES"]][["i2"]] <- list(HYPER=list(21, 10), HYPO=list(15, 12))
    observed[["OBSERVATION"]][["SITES"]][["iAll"]] <- list(HYPER=list(1), HYPO=list(3))
    observed[["OBSERVATION"]][["TILES"]] <- list()
    observed[["OBSERVATION"]][["TILES"]][["i2"]] <- list(HYPER=list(2000, 3000), HYPO=list(2000, 3000))
    observed[["OBSERVATION"]][["TILES"]][["iAll"]] <- list(HYPER=list(0), HYPO=list(0))

    permutated <- list()
    permutated[["PERMUTATION"]] <- list()
    cas_01 <- list()
    cas_01[["i2"]] <- list(HYPER=list(5, 7), HYPO=list(10, 11))
    cas_01[["iAll"]] <- list(HYPER=list(2), HYPO=list(3))
    cas_02 <- list()
    cas_02[["i2"]] <- list(HYPER=list(8, 9), HYPO=list(4, 7))
    cas_02[["iAll"]] <- list(HYPER=list(3), HYPO=list(0))
    cas_03 <- list()
    cas_03[["i2"]] <- list(HYPER=list(10, 7), HYPO=list(11, 7))
    cas_03[["iAll"]] <- list(HYPER=list(1), HYPO=list(2))
    permutated[["PERMUTATION"]][["SITES"]] <- list(cas_01, cas_02, cas_03)
    cas_01 <- list()
    cas_01[["i2"]] <- list(HYPER=list(0, 0), HYPO=list(1000, 3000))
    cas_01[["iAll"]] <- list(HYPER=list(0), HYPO=list(1000))
    cas_02 <- list()
    cas_02[["i2"]] <- list(HYPER=list(1000, 1000), HYPO=list(1000, 1000))
    cas_02[["iAll"]] <- list(HYPER=list(0), HYPO=list(0))
    cas_03 <- list()
    cas_03[["i2"]] <- list(HYPER=list(0, 3000), HYPO=list(2000, 0))
    cas_03[["iAll"]] <- list(HYPER=list(0), HYPO=list(0))
    permutated[["PERMUTATION"]][["TILES"]] <- list(cas_01, cas_02, cas_03)


    obs <- tryCatch(mergePermutationAndObservation(permutationResults = permutated,
                            observationResults = observed),
                    error=conditionMessage)

    exp <- list()
    exp[["PERMUTATION"]] <- permutated[["PERMUTATION"]]
    exp[["OBSERVATION"]] <- observed[["OBSERVATION"]]
    class(exp) <- "methylInheritanceAllResults"

    message <- paste0(" test.mergePermutationAndObservation_good_01() ",
                      "- Valid parameters for mergePermutationAndObservation() ",
                      "did not generated expected results.")

    checkEquals(obs, exp, msg = message)
}

###########################################################
## plotGraph() function
###########################################################

# Test result when all parameters are good
test.plotGraph_good_01 <- function() {

    g <- extractInfo(allResults = methylInheritanceResults,
                                type = "sites", inter="i2", 1)

    obs <- plotGraph(g)

    message <- paste0(" test.plotGraph_good_01() ",
                      "- Valid parameters for plotGraph did not generated expected results.")

    checkTrue("gtable" %in% class(obs), msg = message)
    checkEquals(class(obs[[1]]), "list", msg = message)
    checkEquals(class(obs[[2]]), "data.frame", msg = message)
}
