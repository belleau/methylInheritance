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

###########################################################

## runPermutationUsingRDSFile() function

###########################################################

## Test when methylKitRDSFile is not a valid RDS file name
test.runPermutationUsingRDSFile_methylKitRDSFile_not_valid <- function() {
    obs <- tryCatch(runPermutationUsingRDSFile(
        methylKitRDSFile = "HI",  outputDir = NULL,
        nbrPermutations = 2, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "The file \"HI\" does not exist."

    message <- paste0(" test.runPermutationUsingRDSFile_methylKitRDSFile_not_valid() ",
                      "- Not valid file for methylKitRDSFile did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}
