#' methylInheritance: Permutation
#'
#' This package does a permutation analysis, based on Monte Carlo sampling,
#' for testing the hypothesis that the number of conserved differentially
#' methylated elements (sites, tiles or regions), between
#' several generations, is associated to an effect inherited from a treatment
#' and that stochastic effect can be dismissed.
#'
#' @docType package
#'
#' @name methylInheritance-package
#'
#' @aliases methylInheritance-package methylInheritance
#'
#' @author Astrid Deschênes,
#' Pascal Belleau and
#' Arnaud Droit
#'
#' Maintainer:
#' Astrid Deschenes <astrid-louise.deschenes@@crchudequebec.ulaval.ca>
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{runPermutationUsingRDSFile}} { for running a
#'     permutation analysis on the specified multi-generational dataset in
#'     RDS format}
#' }
#'
#' @keywords package
NULL

#' All samples information, formated by \code{methylKit}, in a
#' \code{methylRawList} format (for demo purpose).
#'
#' The object contains a list with 3 entries. Each entry correspond to the
#' information for one generation (first entry = first generation, etc..).
#' There is 12 samples (6 controls and 6 cases) for each generation. Each
#' sample information is stored in a \code{methylRaw} object.
#'
#' This dataset can be
#' used to test the \code{runPermutationUsingRDS} function.
#'
#' @name samplesForTransgenerationalAnalysis
#'
#' @docType data
#'
#' @aliases samplesForTransgenerationalAnalysis
#'
#' @format A \code{methylRawList} containing all samples information. Each
#' entry correspond to the information for one generation (first entry = first
#' generation, etc..). Each sample information is stored in a \code{methylRaw}
#' object.
#'
#' @return A \code{methylRawList} containing all samples information. Each
#' entry correspond to the information for one generation (first entry = first
#' generation, etc..). Each sample information is stored in a \code{methylRaw}
#' object.
#' @seealso
#' \itemize{
#'     \item \code{\link{runPermutationUsingMethylKitInfo}} {for running a
#'     permutation analysis using methylKit info entry}
#' }
#'
#' @usage data(samplesForTransgenerationalAnalysis)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading dataset
#' data(samplesForTransgenerationalAnalysis)
#'
#' ## Run a permutation analysis
#' \dontrun{runPermutationUsingMethylKitInfo(methylKitRDSFile =
#' samplesForTransgenerationalAnalysis, type = "sites", nbrPermutations = 3,
#' vSeed = 2001)}
#'
NULL
