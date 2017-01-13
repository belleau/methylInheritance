###################################################
# Created by Astrid Deschenes
# 2017-01-05
###################################################

###################################################
## Test the methylInheritanceMethods functions
###################################################


DIRECTORY <- system.file("extdata", package = "methylInheritance")

METHYL_OBJ_FILE_02 <- dir(system.file("extdata", package = "methylInheritance"),
                       pattern = "methylObj_002.RDS", full.names = TRUE)

METHYL_OBJ_02 <- readRDS(METHYL_OBJ_FILE_02)

data(analysisResults)

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
        methylKitRDSFile = METHYL_OBJ_FILE_02,  outputDir = NULL,
        runObservedAnalysis = TRUE,
        nbrPermutations = 2, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = FALSE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    cas_01 <- list()
    cas_01[["SITES"]] <- list()
    cas_01[["SITES"]][["i2"]] <- list()
    cas_01[["SITES"]][["i2"]][["HYPER"]] <- list(7,8)
    cas_01[["SITES"]][["i2"]][["HYPO"]] <- list(4,4)
    cas_01[["SITES"]][["iAll"]] <- list()
    cas_01[["SITES"]][["iAll"]][["HYPER"]] <- list(3)
    cas_01[["SITES"]][["iAll"]][["HYPO"]] <- list(1)
    cas_01[["TILES"]]<-list()
    cas_01[["TILES"]][["i2"]] <- list()
    cas_01[["TILES"]][["i2"]][["HYPER"]] <- list(0,0)
    cas_01[["TILES"]][["i2"]][["HYPO"]] <- list(1000, 1800)
    cas_01[["TILES"]][["iAll"]] <- list()
    cas_01[["TILES"]][["iAll"]][["HYPER"]] <- list(0)
    cas_01[["TILES"]][["iAll"]][["HYPO"]] <- list(0)
    cas_02 <- list()
    cas_02[["SITES"]] <- list()
    cas_02[["SITES"]][["i2"]] <- list()
    cas_02[["SITES"]][["i2"]][["HYPER"]] <- list(6,9)
    cas_02[["SITES"]][["i2"]][["HYPO"]] <- list(6,4)
    cas_02[["SITES"]][["iAll"]] <- list()
    cas_02[["SITES"]][["iAll"]][["HYPER"]] <- list(2)
    cas_02[["SITES"]][["iAll"]][["HYPO"]] <- list(2)
    cas_02[["TILES"]]<-list()
    cas_02[["TILES"]][["i2"]] <- list()
    cas_02[["TILES"]][["i2"]][["HYPER"]] <- list(0,0)
    cas_02[["TILES"]][["i2"]][["HYPO"]] <- list(0, 0)
    cas_02[["TILES"]][["iAll"]] <- list()
    cas_02[["TILES"]][["iAll"]][["HYPER"]] <- list(0)
    cas_02[["TILES"]][["iAll"]][["HYPO"]] <- list(0)

    permutated <- list(cas_01, cas_02)

    observed <- list()
    observed[["SITES"]] <- list()
    observed[["SITES"]][["i2"]] <- list()
    observed[["SITES"]][["i2"]][["HYPER"]] <- list(9,5)
    observed[["SITES"]][["i2"]][["HYPO"]] <- list(7,7)
    observed[["SITES"]][["iAll"]] <- list()
    observed[["SITES"]][["iAll"]][["HYPER"]] <- list(0)
    observed[["SITES"]][["iAll"]][["HYPO"]] <- list(1)
    observed[["TILES"]]<-list()
    observed[["TILES"]][["i2"]] <- list()
    observed[["TILES"]][["i2"]][["HYPER"]] <- list(0,0)
    observed[["TILES"]][["i2"]][["HYPO"]] <- list(0,2700)
    observed[["TILES"]][["iAll"]] <- list()
    observed[["TILES"]][["iAll"]][["HYPER"]] <- list(0)
    observed[["TILES"]][["iAll"]][["HYPO"]] <- list(0)

    exp <- list()
    exp[["OBSERVATION"]] <- observed
    exp[["PERMUTATION"]] <- permutated

    message <- paste0(" test.runPermutationUsingRDSFile_good_001() ",
                      "- All valid parameters did not generated expected result.")

    checkEquals(obs, exp, msg = message)
}


###########################################################
# runAnalysisUsingRDSFile() function
###########################################################

## Test when methylKitRDSFile is not a valid RDS file name
test.runAnalysisUsingRDSFile_methylKitRDSFile_not_valid <- function() {
    obs <- tryCatch(runAnalysisUsingRDSFile(
        methylKitRDSFile = "ALLO",  outputDir = NULL,
        nbrPermutations = 2, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = FALSE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "The file \"ALLO\" does not exist."

    message <- paste0(" test.runAnalysisUsingRDSFile_methylKitRDSFile_not_valid() ",
                      "- Not valid file for methylKitRDSFile did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}


###########################################################
# formatForGraph() function
###########################################################

## Test result when all parameters are good
test.formatForGraph_good_01 <- function() {
    obs <- tryCatch(formatForGraph(analysisandPermutationResults =
                    analysisResults, type = "sites", inter="i2", 1),
                    error=conditionMessage)

    exp <- data.frame(TYPE = rep(c("HYPO","HYPER"), 21),
                      result = c(2,4,2,4,4,3,1,5,3,3,4,2,0,0,0,1,1,0,6,2,2,5,1,
                                 3,2,4,222,67,6,4,183,53,1,6,34,102,2,2,4,3,2,2),
                      SOURCE = c("OBSERVED", "OBSERVED",
                                    rep("PERMUTATION", 40)))

    message <- paste0(" test.formatForGraph_good_01() ",
                      "- Valid parameters for formatForGraph did not generated expected results.")

    checkEquals(obs, exp, msg = message)
}

