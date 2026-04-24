#' Cochlear implant speech recognition data
#'
#' Longitudinal speech recognition outcomes for two groups (A/B), including
#' incomplete records, from Zimmerman and Nunez-Anton antedependence book
#' companion data. This dataset is continuous-response data suitable for
#' Gaussian AD modeling.
#'
#' @format A list with six components:
#'   \describe{
#'     \item{y}{numeric matrix of dimension N by n_time containing all subjects}
#'     \item{y_A}{numeric matrix for Group A subjects}
#'     \item{y_B}{numeric matrix for Group B subjects}
#'     \item{blocks}{integer vector of length N (1 = Group A, 2 = Group B)}
#'     \item{group}{character vector of group labels ("A"/"B")}
#'     \item{time}{integer vector of measurement occasions}
#'   }
#'
#' @source \url{https://homepage.divms.uiowa.edu/~dzimmer/Data-for-AD/speech_recognition_data.txt}
"cochlear_implant"
