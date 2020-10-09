#' Simple data set with classifications
#' 
#' Example data for use with CI.classNRI
#' 
#' @format data frame with 3 columns
#' \describe{
#'   \item{ref_class}{The class of the baseline model. Must be a factor}
#'   \item{new_class}{The class  of the new model. Must be a factor}
#'   \item{Outcome}{The outcome of interest (Low or High). Must be a factor}
#' }
"data_class"

#' Simple data set with risk predictions
#' 
#' Example data for use with CI.raplot
#' 
#' @format data frame with 3 columns
#' \describe{
#'   \item{ref}{The prediction from the baseline model}
#'   \item{new}{The prediction from the new model}
#'   \item{outcome}{The outcome of interest (0 or 1)}
#' }
"data_risk"