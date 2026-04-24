#' Labor force longitudinal categorical data (Table 1)
#'
#' Five-year employment-status sequences reconstructed from Table 1 in the
#' labor-force example used in Xie and Zimmerman score/Wald antedependence
#' testing work. Category coding is 1 = employed, 2 = unemployed.
#'
#' @format A list with five components:
#'   \describe{
#'     \item{y}{integer matrix of dimension N by 5 containing expanded subject-level sequences}
#'     \item{counts}{data frame with columns Y1, Y2, Y3, Y4, Y5, Count}
#'     \item{n_categories}{number of categories (2)}
#'     \item{time}{integer vector of calendar years (1967:1971)}
#'     \item{status_labels}{character vector c("employed", "unemployed")}
#'   }
#'
#' @source Table 1 (labor-force example) from:
#'   Xie, Y. and Zimmerman, D. L. (2013).
#'   Score and Wald tests for antedependence in categorical longitudinal data.
"labor_force_cat"
