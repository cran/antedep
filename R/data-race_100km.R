#' 100km race split-time data
#'
#' Split times (in minutes) from a 100km race example with 10 consecutive sections.
#' This is continuous-response longitudinal data suitable for Gaussian AD modeling.
#'
#' @format A list with five components:
#'   \describe{
#'     \item{y}{numeric matrix of dimension N by 10 containing section split times}
#'     \item{age}{numeric vector of subject ages (may include missing values)}
#'     \item{subject}{integer subject identifiers}
#'     \item{time}{integer vector of section indices (1:10)}
#'     \item{section_labels}{character vector of split-time column names}
#'   }
#'
#' @source Combined table extracted from textbook race example data, stored in
#'   \code{data-raw/external/100km_race_combined_extracted.csv}.
"race_100km"
