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

## Test when qvalue is above 1
test.validateRunPermutationUsingRDS_qvalue_above_1 <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 1.01,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "qvalue must be a positive double between [0,1]"

    message <- paste0(" test.validateRunPermutationUsingRDS_qvalue_above_1() ",
                      "- Not valid qvalue did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when qvalue is negative
test.validateRunPermutationUsingRDS_qvalue_negative <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = -0.01,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "qvalue must be a positive double between [0,1]"

    message <- paste0(" test.validateRunPermutationUsingRDS_qvalue_negative() ",
                      "- Not valid qvalue did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when maxPercReads is not a number
test.validateRunPermutationUsingRDS_maxPercReads_not_number <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = "lala", destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "maxPercReads must be a positive double between [0,100]"

    message <- paste0(" test.validateRunPermutationUsingRDS_maxPercReads_not_number() ",
                      "- Not valid maxPercReads did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when maxPercReads is above 100
test.validateRunPermutationUsingRDS_maxPercReads_above_100 <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 100.1, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "maxPercReads must be a positive double between [0,100]"

    message <- paste0(" test.validateRunPermutationUsingRDS_maxPercReads_above_100() ",
                      "- Not valid maxPercReads did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when maxPercReads is negative
test.validateRunPermutationUsingRDS_maxPercReads_negative <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = -0.1, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "maxPercReads must be a positive double between [0,100]"

    message <- paste0(" test.validateRunPermutationUsingRDS_maxPercReads_negative() ",
                      "- Not valid maxPercReads did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}


## Test when destrand is a number
test.validateRunPermutationUsingRDS_destrand_number <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = 20, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "destrand must be a logical"

    message <- paste0(" test.validateRunPermutationUsingRDS_destrand_number() ",
                      "- Not valid destrand did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}


## Test when minCovBasesForTiles is a string and type is both
test.validateRunPermutationUsingRDS_minCovBasesForTiles_string_type_both <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL, type = "both",
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = "ici",
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "minCovBasesForTiles must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_minCovBasesForTiles_string_type_both() ",
                      "- Not valid minCovBasesForTiles did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when minCovBasesForTiles is negative and type is both
test.validateRunPermutationUsingRDS_minCovBasesForTiles_negative_type_both <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL, type = "both",
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = -1,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "minCovBasesForTiles must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_minCovBasesForTiles_negative_type_both() ",
                      "- Not valid minCovBasesForTiles did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when minCovBasesForTiles is a string and type is tiles
test.validateRunPermutationUsingRDS_minCovBasesForTiles_string_type_tiles <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL, type = "tiles",
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = "a",
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "minCovBasesForTiles must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_minCovBasesForTiles_string_type_tiles() ",
                      "- Not valid minCovBasesForTiles did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when minCovBasesForTiles is negative and type is tiles
test.validateRunPermutationUsingRDS_minCovBasesForTiles_negative_type_tiles <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL, type = "tiles",
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = -1,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "minCovBasesForTiles must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_minCovBasesForTiles_negative_type_tiles() ",
                      "- Not valid minCovBasesForTiles did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}


## Test when tileSize is a string and type is both
test.validateRunPermutationUsingRDS_tileSize_string_type_both <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL, type = "both",
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = "yes", stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "tileSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_minCovBasesForTiles_string_type_both() ",
                      "- Not valid tileSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when tileSize is zero and type is both
test.validateRunPermutationUsingRDS_tileSize_zero_type_both <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL, type = "both",
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 0, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "tileSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_tileSize_zero_type_both() ",
                      "- Not valid tileSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when tileSize is negative and type is both
test.validateRunPermutationUsingRDS_tileSize_negative_type_both <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL, type = "both",
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = -1, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "tileSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_tileSize_negative_type_both() ",
                      "- Not valid tileSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when tileSize is a string and type is tiles
test.validateRunPermutationUsingRDS_tileSize_string_type_tiles <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL, type = "tiles",
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = "yes", stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "tileSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_tileSize_string_type_tiles() ",
                      "- Not valid tileSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}


## Test when tileSize is zero and type is tiles
test.validateRunPermutationUsingRDS_tileSize_zero_type_tiles <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL, type = "tiles",
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 0, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "tileSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_tileSize_zero_type_tiles() ",
                      "- Not valid tileSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when tileSize is negative and type is tiles
test.validateRunPermutationUsingRDS_tileSize_negative_type_tiles <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL, type = "tiles",
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = -1, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "tileSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_tileSize_negative_type_tiles() ",
                      "- Not valid tileSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when stepSize is a string and type is tiles
test.validateRunPermutationUsingRDS_stepSize_string_type_tiles <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL, type = "tiles",
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 10000, stepSize = "one", vSeed = 222),
        error=conditionMessage)

    exp <- "stepSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_stepSize_string_type_tiles() ",
                      "- Not valid stepSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when stepSize is zero and type is tiles
test.validateRunPermutationUsingRDS_stepSizee_zero_type_tiles <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL, type = "tiles",
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 10000, stepSize = 0, vSeed = 222),
        error=conditionMessage)

    exp <- "stepSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_stepSizee_zero_type_tiles() ",
                      "- Not valid stepSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when stepSize is negative and type is tiles
test.validateRunPermutationUsingRDS_stepSize_negative_type_tiles <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL, type = "tiles",
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 10000, stepSize = -1, vSeed = 222),
        error=conditionMessage)

    exp <- "stepSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_stepSize_negative_type_tiles() ",
                      "- Not valid stepSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when stepSize is a string and type is both
test.validateRunPermutationUsingRDS_stepSize_string_type_both <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL, type = "both",
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 10000, stepSize = "one", vSeed = 222),
        error=conditionMessage)

    exp <- "stepSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_stepSize_string_type_both() ",
                      "- Not valid stepSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when stepSize is zero and type is both
test.validateRunPermutationUsingRDS_stepSizee_zero_type_both <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL, type = "both",
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 10000, stepSize = 0, vSeed = 222),
        error=conditionMessage)

    exp <- "stepSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_stepSizee_zero_type_both() ",
                      "- Not valid stepSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when stepSize is negative and type is both
test.validateRunPermutationUsingRDS_stepSize_negative_type_both <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL, type = "both",
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 10000, stepSize = -1, vSeed = 222),
        error=conditionMessage)

    exp <- "stepSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_stepSize_negative_type_both() ",
                      "- Not valid stepSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when vSeed is a string
test.validateRunPermutationUsingRDS_vSeed_string <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL, type = "both",
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 10000, stepSize = 100, vSeed = "222"),
        error=conditionMessage)

    exp <- "vSeed must be either -1 or a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_vSeed_string() ",
                      "- Not valid vSeed did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when vSeed is a string
test.validateRunPermutationUsingRDS_vSeed_string <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL, type = "both",
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 10000, stepSize = 100, vSeed = "33"),
        error=conditionMessage)

    exp <- "vSeed must be an integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingRDS_vSeed_string() ",
                      "- Not valid vSeed did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when all parameters valid
test.validateRunPermutationUsingRDS_all_valid_parameters_01 <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingRDS(
        methylKitRDSFile = METHYL_OBJ_FILE,  outputDir = NULL, type = "sites",
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = -3,
        tileSize = -1, stepSize = -2, vSeed = 22),
        error=conditionMessage)

    exp <- 0

    message <- paste0(" test.validateRunPermutationUsingRDS_all_valid_parameters_01() ",
                      "- All valid parameters did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}
