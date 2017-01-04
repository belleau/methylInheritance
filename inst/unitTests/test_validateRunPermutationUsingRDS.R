###################################################

# Created by Astrid Deschenes

# 2016-12-21

###################################################

###################################################

## Test the methylInheritanceInternalMethods functions

###################################################


DIRECTORY <- system.file("extdata", package = "methylInheritance")

METHYL_OBJ_FILE <- dir(system.file("extdata", package = "methylInheritance"),
                pattern = "methylObj_001.RDS", full.names = TRUE)

###########################################################

## validateRunPermutationUsingRDS() function

###########################################################

## Test when methylKitObject is not a valid RDS file name
test.validateRunPermutationUsingRDS_not_methylKitObject <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
            methylKitRDSFile = "HI",  outputDir = NULL,
            nbrPermutations = 2, nbrCores = 1, nbrCoresDiffMeth = 1,
            minReads = 10, minMethDiff = 10, qvalue = 0.05,
            maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
            tileSize = 1000, stepSize = 100, vSeed = 222),
            error=conditionMessage)

    exp <- "The file \"HI\" does not exist."

    message <- paste0(" test.validateRunPermutationUsingRDS_not_methylKitObject() ",
                      "- Not valid file for methylKitRDSFile did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when outputDir is a number
test.validateRunPermutationUsingRDS_outputDir_as_number <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = 33,
        nbrPermutations = 2, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "output_dir must be a character string or NULL"

    message <- paste0(" test.validateRunPermutationUsingRDS_outputDir_as_number() ",
                      "- Not valid outputDir did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

