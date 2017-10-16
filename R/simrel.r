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
    ..parameters = list(
      n        = 100,
      p        = 10,
      m        = 1,
      q        = 5,
      relpos   = c(1, 2, 4),
      R2       = 0.8,
      gamma    = 0.7,
      eta      = 0,
      rho      = NULL,
      ypos     = NULL,
      type     = "univariate"
    ),
    ..properties = list(
      relpred = NULL,
      irrelpred = NULL,
      lambda = NULL,
      sigma_z = NULL,
      sigma_zinv = NULL,
      sigma_y = NULL,
      sigma_w = NULL,
      sigma_zy = NULL,
      sigma_zw = NULL,
      sigma = NULL,
      rotation_x = NULL,
      rotation_y = NULL,
      beta_z = NULL,
      beta = NULL,
      Rsq = NULL,
      minerror = NULL
    ),
    ..data = list(
      x = NULL,
      y = NULL
    ),
    ..get_cov = function(pos, Rsq, eta = 1, p = p, lambda = lambda){
      out      <- vector("numeric", p)
      alph     <- runif(length(pos), -1, 1)
      out[pos] <- sign(alph) * sqrt(Rsq * abs(alph) / sum(abs(alph)) * lambda[pos] * eta)
      return(out)
    },
    ..get_rotate = function(predPos){
      n    <- length(predPos)
      Qmat <- matrix(rnorm(n ^ 2), n)
      Qmat <- scale(Qmat, scale = FALSE)
      qr.Q(qr(Qmat))
    },
    ..is_pd = function(){
      sigma <- self$get_properties("sigma")
      all(eigen(sigma)$values > 0)
    }
  ),
  active = list(
    parameters = function() private$..parameters,
    lambda = function() {
      gamma <- self$get_parameters("gamma")
      p <- self$get_parameters("p")
      return((exp(-gamma*(1:p)))/(exp(-gamma)))
    },
    eta = function() {
      eta <- self$get_parameters("eta")
      m <- self$get_parameters("m")
      return((exp(-eta * (1:m)))/(exp(-eta)))
    },
    sigma_z = function() diag(self$get_properties("lambda")),
    sigma_zinv = function() {
      diag(1/self$get_properties("lambda"))
    },
    sigma_y = function() return(1),
    predpos = function() {
      type <- private$..parameters$type
      p <- private$..parameters$p
      q <- private$..parameters$q
      relpos <- private$..parameters$relpos
      relpos <- if (!is.list(relpos)) list(relpos) else relpos
      irrelpos <- setdiff(seq_len(p), Reduce(union, relpos))
      out <- lapply(seq_along(relpos), function(i){
        pos      <- relpos[[i]]
        ret      <- c(pos, sample(irrelpos, q[i] - length(pos)))
        irrelpos <<- setdiff(irrelpos, ret)
        return(ret)
      })
      out <- if (type == "univariate") unlist(out) else out
      return(out)
    },
    get_data = function(){
      sigma <- self$get_properties("sigma")
      n <- self$get_parameters("n")
      p <- self$get_parameters("p")
      m <- self$get_parameters("m")
      type <- self$get_parameters("type")
      rotation_x <- self$get_properties("rotation_x")
      rotation_y <- self$get_properties("rotation_y")
      sigma_rot <- chol(sigma)
      train_cal <- matrix(rnorm(n * (p + m), 0, 1), nrow = n) %*% sigma_rot
      Z <- train_cal[, (m + 1):(m + p), drop = F]
      X <- Z %*% t(rotation_x)
      W <- train_cal[, 1:m, drop = F]
      Y <- if (!type == "univariate") W %*% t(rotation_y) else W
      list(y = unname(Y), x = X)
    }
  ),
  public = list(
    initialize = function(...) {
      params <- list(...)
      for (key in names(params)) {
        private$..parameters[[key]] <- params[[key]]
      }
      
      ## Adding Properties to Simrel Object
      self$set_properties("predpos", self$predpos)
      self$set_properties("lambda", self$lambda)
      self$set_properties("sigma_z", self$sigma_z)
      self$set_properties("sigma_zinv", self$sigma_zinv)
      self$set_properties("sigma_y", self$sigma_y)
    },
    get_parameters = function(key) private$..parameters[[key]],
    set_properties = function(key, value) {
      private$..properties[[key]] <- value
    },
    get_properties = function(key) private$..properties[[key]]
  )
)

#' @title Uni-Response Simulation of Linear Model Data
#' @description Simulates univariate linear model data where users can specify few parameters for simulating data with wide range of properties.
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

UniSimrel <- R6::R6Class(
  "UniSimrel",
  inherit = Simrel,
  active = list(
    sigma_zy = function() {
      relpos <- self$get_parameters("relpos")
      R2 <- self$get_parameters("R2")
      p <- self$get_parameters("p")
      lambda <- self$get_properties("lambda")
      private$..get_cov(pos = relpos, Rsq = R2, p = p, lambda = lambda)
    },
    sigma = function() {
      sigma_zy <- self$get_properties("sigma_zy")
      sigma_y <- self$get_properties("sigma_y")
      sigma_z <- self$get_properties("sigma_z")
      rbind(c(sigma_y, t(sigma_zy)), cbind(sigma_zy,  sigma_z))
    },
    irrelpred = function() {
      relpred <- self$get_properties("relpred")
      p <- self$get_parameters("p")
      irrelpred <- setdiff(1:p, relpred)
    },
    rotation_x = function() {
      irrelpred <- self$get_properties("irrelpred")
      relpred <- self$get_properties("relpred")
      p <- self$get_parameters("p")
      out <- diag(p)
      out[irrelpred, irrelpred] <- private$..get_rotate(irrelpred)
      out[relpred, relpred] <- private$..get_rotate(relpred)
      return(out)
    },
    beta_z = function() {
      sigma_zinv <- self$get_properties("sigma_zinv")
      sigma_zy <- self$get_properties("sigma_zy")
      return(sigma_zinv %*% sigma_zy)
    },
    beta = function() {
      rotation_x <- self$get_properties("rotation_x")
      beta_z <- self$get_properties("beta_z")
      return(rotation_x %*% beta_z)
    },
    Rsq = function() {
      beta_z <- self$get_properties("beta_z")
      sigma_zy <- self$get_properties("sigma_zy")
      t(beta_z) %*% sigma_zy
    },
    minerror = function() {
      sigma_y <- self$get_properties("sigma_y")
      Rsq <- self$get_properties("Rsq")
      sigma_y - Rsq
    }
  ),
  public = list(
    initialize = function(...){
      super$initialize(...)
      self$set_properties("sigma_zy", self$sigma_zy)
      self$set_properties("sigma", self$sigma)
      self$set_properties("irrelpred", self$irrelpred)
      self$set_properties("rotation_x", self$rotation_x)
      self$set_properties("beta_z", self$beta_z)
      self$set_properties("beta", self$beta)
      self$set_properties("Rsq", self$Rsq)
      self$set_properties("minerror", self$minerror)
      
      data <- self$get_data
      private$..data[["x"]] <- data$x
      private$..data[["y"]] <- data$y
    },
    data = function() {
      data.frame(x = I(private$..data[["x"]]),
                 y = I(private$..data[["y"]]))
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
  active = list(
    sigma_zw = function() {
      relpos <- self$get_parameters("relpos")
      R2 <- self$get_parameters("R2")
      p <- self$get_parameters("p")
      eta <- self$get_parameters("eta")
      lambda <- self$get_properties("lambda")
      mapply(private$..get_cov, pos = relpos, Rsq = R2, eta = eta, 
             MoreArgs = list(p = p, lambda = lambda))
    },
    sigma = function() {
      sigma_zw <- self$get_properties("sigma_zw")
      sigma_w <- self$get_properties("sigma_w")
      sigma_z <- self$get_properties("sigma_z")
      rbind(c(sigma_w, t(sigma_zw)), cbind(sigma_zw,  sigma_z))
    },
    irrelpred = function() {
      relpred <- self$get_properties("relpred")
      p <- self$get_parameters("p")
      irrelpred <- setdiff(1:p, relpred)
    },
    rotation_x = function() {
      irrelpred <- self$get_properties("irrelpred")
      relpred <- self$get_properties("relpred")
      p <- self$get_parameters("p")
      out <- diag(p)
      out[irrelpred, irrelpred] <- private$..get_rotate(irrelpred)
      out[relpred, relpred] <- private$..get_rotate(relpred)
      return(out)
    },
    beta_z = function() {
      sigma_zinv <- self$get_properties("sigma_zinv")
      sigma_zy <- self$get_properties("sigma_zy")
      return(sigma_zinv %*% sigma_zy)
    },
    beta = function() {
      rotation_x <- self$get_properties("rotation_x")
      beta_z <- self$get_properties("beta_z")
      return(rotation_x %*% beta_z)
    },
    Rsq = function() {
      beta_z <- self$get_properties("beta_z")
      sigma_zy <- self$get_properties("sigma_zy")
      t(beta_z) %*% sigma_zy
    },
    minerror = function() {
      sigma_y <- self$get_properties("sigma_y")
      Rsq <- self$get_properties("Rsq")
      sigma_y - Rsq
    }
  ),
  public = list(
    initialize = function(...){
      super$initialize(...)
      self$set_properties("sigma_zw", self$sigma_zw)
      browser()
      self$set_properties("sigma", self$sigma)
      self$set_properties("irrelpred", self$irrelpred)
      self$set_properties("rotation_x", self$rotation_x)
      self$set_properties("beta_z", self$beta_z)
      self$set_properties("beta", self$beta)
      self$set_properties("Rsq", self$Rsq)
      self$set_properties("minerror", self$minerror)
      
      data <- self$get_data
      private$..data[["x"]] <- data$x
      private$..data[["y"]] <- data$y
    },
    data = function() {
      data.frame(x = I(private$..data[["x"]]),
                 y = I(private$..data[["y"]]))
    }
  )
)

