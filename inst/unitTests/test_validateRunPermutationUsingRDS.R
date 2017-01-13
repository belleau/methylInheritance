###################################################

# Created by Astrid Deschenes

# 2016-12-21

###################################################

###################################################
## Test the validateRunPermutationUsingMethylKitInfo function
###################################################

DIRECTORY <- system.file("extdata", package = "methylInheritance")

METHYL_OBJ_FILE <- dir(system.file("extdata", package = "methylInheritance"),
                pattern = "methylObj_001.RDS", full.names = TRUE)

METHYL_OBJ <- readRDS(METHYL_OBJ_FILE)

###########################################################
## validateRunPermutationUsingMethylKitInfo() function
###########################################################

## Test when methylKitInfo is a string
test.validateRunPermutationUsingMethylKitInfo_methylKitInfo_string <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
            methylKitInfo = "HI",  outputDir = NULL, runObservedAnalysis = TRUE,
            nbrPermutations = 2, nbrCores = 1, nbrCoresDiffMeth = 1,
            minReads = 10, minMethDiff = 10, qvalue = 0.05,
            maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
            tileSize = 1000, stepSize = 100, vSeed = 222),
            error=conditionMessage)

    exp <- paste0("methylKitInfo must be a list containing \"methylRawList\" ",
            "entries; each entry must contain all \"methylRaw\" objects ",
            "related to one generation")

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_methylKitInfo_string() ",
                "- Not valid methylKitInfo did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when methylKitInfo is a list of integers
test.validateRunPermutationUsingMethylKitInfo_methylKitInfo_list_of_int <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = list(a=c(1,2), b=c(2,2)), type = "sites", outputDir = NULL,
        runObservedAnalysis = TRUE,
        nbrPermutations = 2, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- paste0("methylKitInfo must be a list containing \"methylRawList\" ",
                  "entries; each entry must contain all \"methylRaw\" objects ",
                  "related to one generation")

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_methylKitInfo_list_of_int() ",
                      "- Not valid methylKitInfo did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when outputDir is a number
test.validateRunPermutationUsingMethylKitInfo_outputDir_as_number <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ, type = "sites", outputDir = 33,
        runObservedAnalysis = TRUE,
        nbrPermutations = 2, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "output_dir must be a character string or NULL"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_outputDir_as_number() ",
                      "- Not valid outputDir did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when runObservedAnalysis is a string
test.validateRunPermutationUsingMethylKitInfo_runObservedAnalysis_string <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, runObservedAnalysis = "allo",
        nbrPermutations = 2, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "runObservedAnalysis must be a logical"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_runObservedAnalysis_string() ",
                      "- Not valid runObservedAnalysis did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when nbrPermutations is a string
test.validateRunPermutationUsingMethylKitInfo_nbrPermutations_as_string <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  type = "sites", outputDir = NULL,
        runObservedAnalysis = TRUE,
        nbrPermutations = "TOTO", nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "nbrPermutations must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_nbrPermutations_as_string() ",
                      "- Not valid nbrPermutations did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when nbrCores is zero
test.validateRunPermutationUsingMethylKitInfo_nbrCores_zero <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ, type = "both", outputDir = NULL,
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 0, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "nbrCores must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_nbrCores_zero() ",
                      "- Not valid nbrCores did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when nbrCores is a negative integer
test.validateRunPermutationUsingMethylKitInfo_nbrCores_negative <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ, type = "both", outputDir = NULL,
        runObservedAnalysis = FALSE,
        nbrPermutations = 3, nbrCores = -1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "nbrCores must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_nbrCores_negative() ",
                      "- Not valid nbrCores did not generated expected message.")


    checkEquals(obs, exp, msg = message)
}

## Test when nbrCoresDiffMeth is zero
test.validateRunPermutationUsingMethylKitInfo_nbrCoresDiffMeth_zero <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ, type = "both", outputDir = NULL,
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 0,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "nbrCoresDiffMeth must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_nbrCoresDiffMeth_zero() ",
                      "- Not valid nbrCoresDiffMeth did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when nbrCoresDiffMeth is a negative integer
test.validateRunPermutationUsingMethylKitInfo_nbrCoresDiffMeth_negative <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ, type = "both", outputDir = NULL,
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = -1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "nbrCoresDiffMeth must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_nbrCoresDiffMeth_negative() ",
                      "- Not valid nbrCoresDiffMeth did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when minReads is zero
test.validateRunPermutationUsingMethylKitInfo_minReads_zero <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ, type = "both", outputDir = NULL,
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 0, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "minReads must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_minReads_zero() ",
                      "- Not valid minReads did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when minReads is negative
test.validateRunPermutationUsingMethylKitInfo_minReads_negative <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ, type = "both", outputDir = NULL,
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = -1, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "minReads must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_minReads_negative() ",
                      "- Not valid minReads did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when minMethDiff is negative
test.validateRunPermutationUsingMethylKitInfo_minMethDiff_negative <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ, type = "both", outputDir = NULL,
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff =-0.1, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "minMethDiff must be a positive double between [0,100]"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_minMethDiff_negative() ",
                      "- Not valid minMethDiff did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}


## Test when minMethDiff is above 100
test.validateRunPermutationUsingMethylKitInfo_minMethDiff_above_100 <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ, type = "both", outputDir = NULL,
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 100.1, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "minMethDiff must be a positive double between [0,100]"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_minMethDiff_above_100() ",
                      "- Not valid minMethDiff did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when qvalue is above 1
test.validateRunPermutationUsingMethylKitInfo_qvalue_above_1 <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  type = "both", outputDir = NULL,
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 1.01,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "qvalue must be a positive double between [0,1]"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_qvalue_above_1() ",
                      "- Not valid qvalue did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when qvalue is negative
test.validateRunPermutationUsingMethylKitInfo_qvalue_negative <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ, type = "both", outputDir = NULL,
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = -0.01,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "qvalue must be a positive double between [0,1]"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_qvalue_negative() ",
                      "- Not valid qvalue did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when maxPercReads is not a number
test.validateRunPermutationUsingMethylKitInfo_maxPercReads_not_number <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ, type = "both", outputDir = NULL,
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = "lala", destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "maxPercReads must be a positive double between [0,100]"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_maxPercReads_not_number() ",
                      "- Not valid maxPercReads did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when maxPercReads is above 100
test.validateRunPermutationUsingMethylKitInfo_maxPercReads_above_100 <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 100.1, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "maxPercReads must be a positive double between [0,100]"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_maxPercReads_above_100() ",
                      "- Not valid maxPercReads did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when maxPercReads is negative
test.validateRunPermutationUsingMethylKitInfo_maxPercReads_negative <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ, type = "both", outputDir = NULL,
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = -0.1, destrand = TRUE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "maxPercReads must be a positive double between [0,100]"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_maxPercReads_negative() ",
                      "- Not valid maxPercReads did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}


## Test when destrand is a number
test.validateRunPermutationUsingMethylKitInfo_destrand_number <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ, type = "both", outputDir = NULL,
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = 20, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "destrand must be a logical"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_destrand_number() ",
                      "- Not valid destrand did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}


## Test when minCovBasesForTiles is a string and type is both
test.validateRunPermutationUsingMethylKitInfo_minCovBasesForTiles_string_type_both <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ, outputDir = NULL, type = "both",
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = "ici",
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "minCovBasesForTiles must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_minCovBasesForTiles_string_type_both() ",
                      "- Not valid minCovBasesForTiles did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when minCovBasesForTiles is negative and type is both
test.validateRunPermutationUsingMethylKitInfo_minCovBasesForTiles_negative_type_both <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, type = "both",
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = -1,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "minCovBasesForTiles must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_minCovBasesForTiles_negative_type_both() ",
                      "- Not valid minCovBasesForTiles did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when minCovBasesForTiles is a string and type is tiles
test.validateRunPermutationUsingMethylKitInfo_minCovBasesForTiles_string_type_tiles <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, type = "tiles",
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = "a",
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "minCovBasesForTiles must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_minCovBasesForTiles_string_type_tiles() ",
                      "- Not valid minCovBasesForTiles did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when minCovBasesForTiles is negative and type is tiles
test.validateRunPermutationUsingMethylKitInfo_minCovBasesForTiles_negative_type_tiles <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, type = "tiles",
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = -1,
        tileSize = 1000, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "minCovBasesForTiles must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_minCovBasesForTiles_negative_type_tiles() ",
                      "- Not valid minCovBasesForTiles did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}


## Test when tileSize is a string and type is both
test.validateRunPermutationUsingMethylKitInfo_tileSize_string_type_both <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, type = "both",
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = "yes", stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "tileSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_minCovBasesForTiles_string_type_both() ",
                      "- Not valid tileSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when tileSize is zero and type is both
test.validateRunPermutationUsingMethylKitInfo_tileSize_zero_type_both <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, type = "both",
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 0, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "tileSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_tileSize_zero_type_both() ",
                      "- Not valid tileSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when tileSize is negative and type is both
test.validateRunPermutationUsingMethylKitInfo_tileSize_negative_type_both <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, type = "both",
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = -1, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "tileSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_tileSize_negative_type_both() ",
                      "- Not valid tileSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when tileSize is a string and type is tiles
test.validateRunPermutationUsingMethylKitInfo_tileSize_string_type_tiles <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, type = "tiles",
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = "yes", stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "tileSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_tileSize_string_type_tiles() ",
                      "- Not valid tileSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}


## Test when tileSize is zero and type is tiles
test.validateRunPermutationUsingMethylKitInfo_tileSize_zero_type_tiles <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, type = "tiles",
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 0, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "tileSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_tileSize_zero_type_tiles() ",
                      "- Not valid tileSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when tileSize is negative and type is tiles
test.validateRunPermutationUsingMethylKitInfo_tileSize_negative_type_tiles <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, type = "tiles",
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = -1, stepSize = 100, vSeed = 222),
        error=conditionMessage)

    exp <- "tileSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_tileSize_negative_type_tiles() ",
                      "- Not valid tileSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when stepSize is a string and type is tiles
test.validateRunPermutationUsingMethylKitInfo_stepSize_string_type_tiles <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, type = "tiles",
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 10000, stepSize = "one", vSeed = 222),
        error=conditionMessage)

    exp <- "stepSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_stepSize_string_type_tiles() ",
                      "- Not valid stepSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when stepSize is zero and type is tiles
test.validateRunPermutationUsingMethylKitInfo_stepSizee_zero_type_tiles <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, type = "tiles",
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 10000, stepSize = 0, vSeed = 222),
        error=conditionMessage)

    exp <- "stepSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_stepSizee_zero_type_tiles() ",
                      "- Not valid stepSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when stepSize is negative and type is tiles
test.validateRunPermutationUsingMethylKitInfo_stepSize_negative_type_tiles <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, type = "tiles",
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 10000, stepSize = -1, vSeed = 222),
        error=conditionMessage)

    exp <- "stepSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_stepSize_negative_type_tiles() ",
                      "- Not valid stepSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when stepSize is a string and type is both
test.validateRunPermutationUsingMethylKitInfo_stepSize_string_type_both <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, type = "both",
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 10000, stepSize = "one", vSeed = 222),
        error=conditionMessage)

    exp <- "stepSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_stepSize_string_type_both() ",
                      "- Not valid stepSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when stepSize is zero and type is both
test.validateRunPermutationUsingMethylKitInfo_stepSizee_zero_type_both <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, type = "both",
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 10000, stepSize = 0, vSeed = 222),
        error=conditionMessage)

    exp <- "stepSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_stepSizee_zero_type_both() ",
                      "- Not valid stepSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when stepSize is negative and type is both
test.validateRunPermutationUsingMethylKitInfo_stepSize_negative_type_both <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, type = "both",
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 10000, stepSize = -1, vSeed = 222),
        error=conditionMessage)

    exp <- "stepSize must be a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_stepSize_negative_type_both() ",
                      "- Not valid stepSize did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when vSeed is a string
test.validateRunPermutationUsingMethylKitInfo_vSeed_string <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, type = "both",
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 10000, stepSize = 100, vSeed = "222"),
        error=conditionMessage)

    exp <- "vSeed must be either -1 or a positive integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_vSeed_string() ",
                      "- Not valid vSeed did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when vSeed is a string
test.validateRunPermutationUsingMethylKitInfo_vSeed_string <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, type = "both",
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = 10,
        tileSize = 10000, stepSize = 100, vSeed = "33"),
        error=conditionMessage)

    exp <- "vSeed must be an integer or numeric"

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_vSeed_string() ",
                      "- Not valid vSeed did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}

## Test when all parameters valid
test.validateRunPermutationUsingMethylKitInfo_all_valid_parameters_01 <- function() {
    obs <- tryCatch(methylInheritance:::validateRunPermutationUsingMethylKitInfo(
        methylKitInfo = METHYL_OBJ,  outputDir = NULL, type = "sites",
        runObservedAnalysis = TRUE,
        nbrPermutations = 3, nbrCores = 1, nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = TRUE, minCovBasesForTiles = -3,
        tileSize = -1, stepSize = -2, vSeed = 22),
        error=conditionMessage)

    exp <- 0

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_all_valid_parameters_01() ",
                      "- All valid parameters did not generated expected message.")

    checkEquals(obs, exp, msg = message)
}
