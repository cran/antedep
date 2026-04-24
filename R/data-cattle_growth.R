#' Cattle growth data (Treatments A and B)
#'
#' Longitudinal cattle growth measurements for two treatment groups from
#' Zimmerman and Nunez-Anton antedependence book companion data. This dataset
#' is continuous-response data suitable for Gaussian AD modeling.
#'
#' @format A list with five components:
#'   \describe{
#'     \item{y}{numeric matrix of dimension N by n_time containing all subjects}
#'     \item{y_A}{numeric matrix for Treatment A subjects}
#'     \item{y_B}{numeric matrix for Treatment B subjects}
#'     \item{blocks}{integer vector of length N (1 = Treatment A, 2 = Treatment B)}
#'     \item{time}{integer vector of measurement occasions}
#'   }
#'
#' @source \url{https://homepage.divms.uiowa.edu/~dzimmer/Data-for-AD/cattle_growth_data_Treatment%20A.txt}
#'   and
#'   \url{https://homepage.divms.uiowa.edu/~dzimmer/Data-for-AD/cattle_growth_data_Treatment%20B.txt}
"cattle_growth"
