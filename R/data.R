#' Example dataset for CMP models
#'
#' A simulated dataset demonstrating different outcome types suitable for
#' conditional mixed process modeling.
#'
#' @format A data frame with 500 rows and 8 variables:
#' \describe{
#'   \item{id}{Observation identifier}
#'   \item{x1}{Continuous predictor variable (standardized)}
#'   \item{x2}{Continuous predictor variable (standardized)}
#'   \item{x3}{Binary predictor variable (0/1)}
#'   \item{y_continuous}{Continuous outcome variable}
#'   \item{y_binary}{Binary outcome variable (0/1)}
#'   \item{y_ordered}{Ordered outcome variable (1/2/3)}
#'   \item{y_censored}{Left-censored continuous outcome (0 = censored)}
#' }
#'
#' @details
#' This dataset was simulated to demonstrate various CMP model types:
#' \itemize{
#'   \item \code{y_continuous}: Can be modeled as continuous
#'   \item \code{y_binary}: Can be modeled as probit
#'   \item \code{y_ordered}: Can be modeled as ordered probit
#'   \item \code{y_censored}: Can be modeled as left-censored (Tobit)
#' }
#'
#' The outcomes are correlated through shared unobserved factors, making
#' this suitable for demonstrating the benefits of joint estimation.
#'
#' @examples
#' data(cmp_example)
#' head(cmp_example)
#'
#' # Simple bivariate model
#' \dontrun{
#' model <- cmp(
#'   formula = list(
#'     y_continuous ~ x1 + x2 + x3,
#'     y_binary ~ x1 + x2
#'   ),
#'   data = cmp_example,
#'   indicators = c("continuous", "probit")
#' )
#' summary(model)
#' }
#'
"cmp_example"