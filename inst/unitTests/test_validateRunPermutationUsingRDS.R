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

## Test when nbrPermutations is a string
test.validateRunPermutationUsingRDS_nbrPermutations_as_string <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL,
        nbrPermutations = "TOTO", nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "nbrPermutations must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_nbrPermutations_as_string() ",
                      "- Not valid nbrPermutations did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when nbrCores is zero
test.validateRunPermutationUsingRDS_nbrCores_zero <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL,
        nbrPermutations = 3, nbrCores = 0, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "nbrCores must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_nbrCores_zero() ",
                      "- Not valid nbrCores did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when nbrCores is a negative integer
test.validateRunPermutationUsingRDS_nbrCores_negative <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL,
        nbrPermutations = 3, nbrCores = -1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "nbrCores must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_nbrCores_negative() ",
                      "- Not valid nbrCores did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when nbrCoresDiffMeth is zero
test.validateRunPermutationUsingRDS_nbrCoresDiffMeth_zero <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 0,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "nbrCoresDiffMeth must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_nbrCoresDiffMeth_zero() ",
                      "- Not valid nbrCoresDiffMeth did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when nbrCoresDiffMeth is a negative integer
test.validateRunPermutationUsingRDS_nbrCoresDiffMeth_negative <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = -1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "nbrCoresDiffMeth must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_nbrCoresDiffMeth_negative() ",
                      "- Not valid nbrCoresDiffMeth did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when minReads is zero
test.validateRunPermutationUsingRDS_minReads_zero <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 0, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "minReads must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_minReads_zero() ",
                      "- Not valid minReads did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when minReads is negative
test.validateRunPermutationUsingRDS_minReads_negative <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = -1, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "minReads must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_minReads_negative() ",
                      "- Not valid minReads did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when minMethDiff is negative
test.validateRunPermutationUsingRDS_minMethDiff_negative <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff =-0.1, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "minMethDiff must be a positive double between [0,100]"

    message <- paste0(" test.validateRunPermutationUsingRDS_minMethDiff_negative() ",
                      "- Not valid minMethDiff did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}


## Test when minMethDiff is above 100
test.validateRunPermutationUsingRDS_minMethDiff_above_100 <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 100.1, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "minMethDiff must be a positive double between [0,100]"

    message <- paste0(" test.validateRunPermutationUsingRDS_minMethDiff_above_100() ",
                      "- Not valid minMethDiff did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}
