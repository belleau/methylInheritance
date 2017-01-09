###################################################

# Created by Astrid Deschenes

# 2017-01-09

###################################################

###################################################

## Test the methylInheritanceInternalMethods functions

###################################################

DIRECTORY <- system.file("extdata", package = "methylInheritance")

METHYL_OBJ_FILE <- dir(system.file("extdata", package = "methylInheritance"),
                       pattern = "methylObj_001.RDS", full.names = TRUE)

METHYL_OBJ <- readRDS(METHYL_OBJ_FILE)

###########################################################
## runOnePermutationOnAllGenerations() function
###########################################################

## Test sites when all parameters are valid
test.validateRunPermutationUsingMethylKitInfo_sites_good_01 <- function() {
    ## Extract information
    set.seed(111)
    allSamples <- sample(unlist(METHYL_OBJ, recursive = FALSE), 36, replace = F)
    treatment <- c(0,0,0,0,0,0,1,1,1,1,1,1)
    sampleList01 <- new("methylRawList", allSamples[1:12],
                         treatment = treatment)
    sampleList02 <- new("methylRawList", allSamples[13:24],
                        treatment = treatment)
    sampleList03 <- new("methylRawList", allSamples[25:36],
                        treatment = treatment)
    input <- list(sample = list(sampleList01, sampleList02, sampleList03),
                  id = 1)

    obs <- tryCatch(methylInheritance:::runOnePermutationOnAllGenerations(
        methylInfoForAllGenerations = input, outputDir = NULL, type = "sites",
        nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = FALSE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100),
        error=conditionMessage)

    exp <- list()
    exp[["SITES"]] <- list()
    exp[["SITES"]][["i2"]] <- list()
    exp[["SITES"]][["i2"]][["HYPER"]] <- list(8, 9)
    exp[["SITES"]][["i2"]][["HYPO"]]  <- list(8, 6)
    exp[["SITES"]][["iAll"]][["HYPER"]]  <- list(0)
    exp[["SITES"]][["iAll"]][["HYPO"]]   <- list(0)

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_sites_good_01() ",
                      "- Valid paramters did not generated expected results.")

    checkEquals(obs, exp, msg = message)
}

## Test tiles when all parameters are valid
test.validateRunPermutationUsingMethylKitInfo_tiles_good_01 <- function() {
    ## Extract information
    set.seed(112241)
    allSamples <- sample(unlist(METHYL_OBJ, recursive = FALSE), 36, replace = F)
    treatment <- c(0,0,0,0,0,0,1,1,1,1,1,1)
    sampleList01 <- new("methylRawList", allSamples[1:12],
                        treatment = treatment)
    sampleList02 <- new("methylRawList", allSamples[13:24],
                        treatment = treatment)
    sampleList03 <- new("methylRawList", allSamples[25:36],
                        treatment = treatment)
    input <- list(sample = list(sampleList01, sampleList02, sampleList03),
                  id = 1)

    obs <- tryCatch(methylInheritance:::runOnePermutationOnAllGenerations(
        methylInfoForAllGenerations = input, outputDir = NULL, type = "tiles",
        nbrCoresDiffMeth = 1,
        minReads = 10, minMethDiff = 10, qvalue = 0.05,
        maxPercReads = 99.9, destrand = FALSE, minCovBasesForTiles = 2,
        tileSize = 1000, stepSize = 100),
        error=conditionMessage)

    exp <- list()
    exp[["TILES"]] <- list()
    exp[["TILES"]][["i2"]] <- list()
    exp[["TILES"]][["i2"]][["HYPER"]] <- list(1700, 0)
    exp[["TILES"]][["i2"]][["HYPO"]]  <- list(0, 0)
    exp[["TILES"]][["iAll"]][["HYPER"]]  <- list(0)
    exp[["TILES"]][["iAll"]][["HYPO"]]   <- list(0)

    message <- paste0(" test.validateRunPermutationUsingMethylKitInfo_tiles_good_01() ",
                      "- Valid paramters did not generated expected results.")

    checkEquals(obs, exp, msg = message)
}
