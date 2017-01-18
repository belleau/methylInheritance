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
    obs <- tryCatch(runPermutationUsingRDSFile(
        methylKitRDSFile = METHYL_OBJ_FILE_01, outputDir = NULL,
        runObservationAnalysis = TRUE,
        nbrPermutations = 2, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 5, qvalue = 0.05,
        maxPercReads = 99.9, destrand = FALSE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    cas_01 <- list()
    cas_01[["SITES"]] <- list()
    cas_01[["SITES"]][["i2"]] <- list()
    cas_01[["SITES"]][["i2"]][["HYPER"]] <- list(9,6)
    cas_01[["SITES"]][["i2"]][["HYPO"]] <- list(2,2)
    cas_01[["SITES"]][["iAll"]] <- list()
    cas_01[["SITES"]][["iAll"]][["HYPER"]] <- list(6)
    cas_01[["SITES"]][["iAll"]][["HYPO"]] <- list(2)
    cas_01[["TILES"]]<-list()
    cas_01[["TILES"]][["i2"]] <- list()
    cas_01[["TILES"]][["i2"]][["HYPER"]] <- list(900,0)
    cas_01[["TILES"]][["i2"]][["HYPO"]] <- list(1000, 2000)
    cas_01[["TILES"]][["iAll"]] <- list()
    cas_01[["TILES"]][["iAll"]][["HYPER"]] <- list(0)
    cas_01[["TILES"]][["iAll"]][["HYPO"]] <- list(100)
    cas_02 <- list()
    cas_02[["SITES"]] <- list()
    cas_02[["SITES"]][["i2"]] <- list()
    cas_02[["SITES"]][["i2"]][["HYPER"]] <- list(0,0)
    cas_02[["SITES"]][["i2"]][["HYPO"]] <- list(1,0)
    cas_02[["SITES"]][["iAll"]] <- list()
    cas_02[["SITES"]][["iAll"]][["HYPER"]] <- list(0)
    cas_02[["SITES"]][["iAll"]][["HYPO"]] <- list(0)
    cas_02[["TILES"]]<-list()
    cas_02[["TILES"]][["i2"]] <- list()
    cas_02[["TILES"]][["i2"]][["HYPER"]] <- list(0,0)
    cas_02[["TILES"]][["i2"]][["HYPO"]] <- list(0,0)
    cas_02[["TILES"]][["iAll"]] <- list()
    cas_02[["TILES"]][["iAll"]][["HYPER"]] <- list(0)
    cas_02[["TILES"]][["iAll"]][["HYPO"]] <- list(0)

    permutated <- list(cas_01, cas_02)

    observed <- list()
    observed[["SITES"]] <- list()
    observed[["SITES"]][["i2"]] <- list()
    observed[["SITES"]][["i2"]][["HYPER"]] <- list(0,0)
    observed[["SITES"]][["i2"]][["HYPO"]] <- list(3,3)
    observed[["SITES"]][["iAll"]] <- list()
    observed[["SITES"]][["iAll"]][["HYPER"]] <- list(0)
    observed[["SITES"]][["iAll"]][["HYPO"]] <- list(2)
    observed[["TILES"]]<-list()
    observed[["TILES"]][["i2"]] <- list()
    observed[["TILES"]][["i2"]][["HYPER"]] <- list(0,0)
    observed[["TILES"]][["i2"]][["HYPO"]] <- list(3800,2000)
    observed[["TILES"]][["iAll"]] <- list()
    observed[["TILES"]][["iAll"]][["HYPER"]] <- list(0)
    observed[["TILES"]][["iAll"]][["HYPO"]] <- list(1000)

    exp <- list()
    exp[["OBSERVATION"]] <- observed
    exp[["PERMUTATION"]] <- permutated

    message <- paste0(" test.runPermutationUsingRDSFile_good_001() ",
                      "- All valid parameters did not generated expected result.")

    checkEquals(obs, exp, msg = message)
}


###########################################################
# runObservationUsingRDSFile() function
###########################################################

## Test when methylKitRDSFile is not a valid RDS file name
test.runObservationUsingRDSFile_methylKitRDSFile_not_valid <- function() {
    obs <- tryCatch(runObservationUsingRDSFile(
        methylKitRDSFile = "ALLO",  outputDir = NULL,
        nbrPermutations = 2, nbrCores = 1, nbrCoresDiffMeth = 1,
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
        outputDir = NULL,
        nbrPermutations = 1, nbrCores = 1, nbrCoresDiffMeth = 1,
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

## Test result when all parameters are good
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
test.loadAllRDSResults_good_01 <- function() {

    obs <- tryCatch(loadAllRDSResults(analysisResultsDir = TEST_DIR,
                        permutationResultsDir = TEST_DIR, doingSites = TRUE,
                        doingTiles = TRUE),
                    error=conditionMessage)

    exp <- list()

    exp[["OBSERVATION"]] <- list()
    exp[["OBSERVATION"]][["SITES"]] <- list()
    exp[["OBSERVATION"]][["SITES"]][["i2"]] <- list(HYPER=list(21, 10), HYPO=list(15, 12))
    exp[["OBSERVATION"]][["SITES"]][["iAll"]] <- list(HYPER=list(1), HYPO=list(3))
    i2 <- list(HYPER=list(2000, 3000), HYPO=list(2000, 3000))
    iAll <- list(HYPER=list(0), HYPO=list(0))
    exp[["OBSERVATION"]][["TILES"]] <- list()
    exp[["OBSERVATION"]][["TILES"]][["i2"]] <- list(HYPER=list(2000, 3000), HYPO=list(2000, 3000))
    exp[["OBSERVATION"]][["TILES"]][["iAll"]] <- list(HYPER=list(0), HYPO=list(0))

    exp[["PERMUTATION"]] <- list()
    cas_01 <- list()
    cas_01[["i2"]] <- list(HYPER=list(5, 7), HYPO=list(10, 11))
    cas_01[["iAll"]] <- list(HYPER=list(2), HYPO=list(3))
    cas_02 <- list()
    cas_02[["i2"]] <- list(HYPER=list(8, 9), HYPO=list(4, 7))
    cas_02[["iAll"]] <- list(HYPER=list(3), HYPO=list(0))
    cas_03 <- list()
    cas_03[["i2"]] <- list(HYPER=list(10, 7), HYPO=list(11, 7))
    cas_03[["iAll"]] <- list(HYPER=list(1), HYPO=list(2))
    exp[["PERMUTATION"]][["SITES"]] <- list(cas_01, cas_02, cas_03)
    cas_01 <- list()
    cas_01[["i2"]] <- list(HYPER=list(0, 0), HYPO=list(1000, 3000))
    cas_01[["iAll"]] <- list(HYPER=list(0), HYPO=list(1000))
    cas_02 <- list()
    cas_02[["i2"]] <- list(HYPER=list(1000, 1000), HYPO=list(1000, 1000))
    cas_02[["iAll"]] <- list(HYPER=list(0), HYPO=list(0))
    cas_03 <- list()
    cas_03[["i2"]] <- list(HYPER=list(0, 3000), HYPO=list(2000, 0))
    cas_03[["iAll"]] <- list(HYPER=list(0), HYPO=list(0))
    exp[["PERMUTATION"]][["TILES"]] <- list(cas_01, cas_02, cas_03)

    class(exp) <- "methylInheritanceAllResults"

    message <- paste0(" test.loadAllRDSResults_good_01() ",
                        "- Valid parameters for loadAllRDSResults() ",
                        "did not generated expected results.")

    checkEquals(obs, exp, msg = message)
}


