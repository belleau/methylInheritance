###################################################
# Created by Astrid Deschenes
# 2017-01-05
###################################################

###################################################
## Test the methylInheritanceMethods functions
###################################################


DIRECTORY <- system.file("extdata", package = "methylInheritance")

METHYL_OBJ_FILE <- dir(system.file("extdata", package = "methylInheritance"),
                       pattern = "methylObj_001.RDS", full.names = TRUE)

METHYL_OBJ <- readRDS(METHYL_OBJ_FILE)

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
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL,
        nbrPermutations = 2, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = FALSE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    cas_01<-list()
    cas_01[["SITES"]] <- list()
    cas_01[["SITES"]][["i2"]] <- list()
    cas_01[["SITES"]][["i2"]][["HYPER"]] <- list(12,19)
    cas_01[["SITES"]][["i2"]][["HYPO"]] <- list(9,9)
    cas_01[["SITES"]][["iAll"]] <- list()
    cas_01[["SITES"]][["iAll"]][["HYPER"]] <- list(7)
    cas_01[["SITES"]][["iAll"]][["HYPO"]] <- list(2)
    cas_01[["TILES"]]<-list()
    cas_01[["TILES"]][["i2"]] <- list()
    cas_01[["TILES"]][["i2"]][["HYPER"]] <- list(0,0)
    cas_01[["TILES"]][["i2"]][["HYPO"]] <- list(1000, 1800)
    cas_01[["TILES"]][["iAll"]] <- list()
    cas_01[["TILES"]][["iAll"]][["HYPER"]] <- list(0)
    cas_01[["TILES"]][["iAll"]][["HYPO"]] <- list(0)
    cas_02<-list()
    cas_02[["SITES"]] <- list()
    cas_02[["SITES"]][["i2"]] <- list()
    cas_02[["SITES"]][["i2"]][["HYPER"]] <- list(15,15)
    cas_02[["SITES"]][["i2"]][["HYPO"]] <- list(12,10)
    cas_02[["SITES"]][["iAll"]] <- list()
    cas_02[["SITES"]][["iAll"]][["HYPER"]] <- list(4)
    cas_02[["SITES"]][["iAll"]][["HYPO"]] <- list(4)
    cas_02[["TILES"]]<-list()
    cas_02[["TILES"]][["i2"]] <- list()
    cas_02[["TILES"]][["i2"]][["HYPER"]] <- list(1000,0)
    cas_02[["TILES"]][["i2"]][["HYPO"]] <- list(1700, 0)
    cas_02[["TILES"]][["iAll"]] <- list()
    cas_02[["TILES"]][["iAll"]][["HYPER"]] <- list(0)
    cas_02[["TILES"]][["iAll"]][["HYPO"]] <- list(0)

    exp <- list(cas_01, cas_02)

    message <- paste0(" test.runPermutationUsingRDSFile_good_001() ",
                      "- All valid parameters did not generated expected result.")

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

