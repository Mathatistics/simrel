#' @title Simulation of Linear Model Data
#' @description Simulates univariate, bivariate and multivariate linear model data where users can specify few parameters for simulating data with wide range of properties.
#' @import R6
#' @param  n Number of training samples
#' @param p Number of predictor variables
#' @param q Number of relevant predictor variables
#' @param relpos Position of relevant predictor components
#' @param gamma Decay factor of eigenvalues of predictor variables
#' @param R2 Coefficient of determination
#' @param ntest (Optional) Number of test samples
#' @param muX (Optional) Mean vector of predictor variables
#' @param muY (Optional) Mean vector of response variables
#' @param sim.obj (Optional) Previously fitted simulation object, the parameters will be taken from the object
#' @param lambda.min (Optional) Minimum value that eigenvalues can be
#' @return simrel object (A list)
#' @rdname simrel
#' @export

Simrel <- R6::R6Class(
  "Simrel",
  private = list(
    n = 100,
    p = 10,
    q = 5,
    relpos = c(1, 2, 4),
    R2 = 0.8,
    gamma = 0.7,
    type = "univariate"
  ),
  public = list(
    initialize = function(n, p, q, relpos, R2, gamma, type) {
      if (!missing(n)) private$n <- n
      if (!missing(p)) private$p <- p
      if (!missing(q)) private$q <- q
      if (!missing(relpos)) private$relpos <- relpos
      if (!missing(R2)) private$R2 <- R2
      if (!missing(gamma)) private$gamma <- gamma
      if (!missing(type)) private$type <- type
    }
  )
)

#' @title Multi-Response Simulation of Linear Model Data
#' @description Simulates multivariate linear model data where users can specify few parameters for simulating data with wide range of properties.
#' @import R6
#' @param  n Number of training samples
#' @param p Number of predictor variables
#' @param q Number of relevant predictor variables
#' @param relpos Position of relevant predictor components
#' @param ypos Position of response components while rotation (see details)
#' @param gamma Decay factor of eigenvalues of predictor variables
#' @param R2 Coefficient of determination
#' @param ntest (Optional) Number of test samples
#' @param muX (Optional) Mean vector of predictor variables
#' @param muY (Optional) Mean vector of response variables
#' @param sim.obj (Optional) Previously fitted simulation object, the parameters will be taken from the object
#' @param lambda.min (Optional) Minimum value that eigenvalues can be
#' @return simrel object (A list)
#' @rdname simrel
#' @export

MultiSimrel <- R6::R6Class(
  "MultiSimrel",
  inherit = Simrel,
  private = list(),
  public = list(
    initialize = function(){
      cat("I'm Multiresponse Multivariate Simrel")
    }
  )
)

