###################################################
# Created by Astrid Deschenes
# 2017-01-25
###################################################

###################################################
## Test the print.rjmcmcNucleosomesMerge.R function
###################################################


### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "methylInheritance" )
}

### }}}

data("methylInheritanceResults")

###########################################################
## print.methylInheritanceAllResults() function
###########################################################

test.print_methylInheritanceAllResults_test_returned_value <- function() {
    result <- print(methylInheritanceResults)

    message <- paste0(" test.print_methylInheritanceAllResults_test_returned_value() ",
                      "- print method did not returned expected value")

    checkEquals(result, methylInheritanceResults, msg = message)
}
