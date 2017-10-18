## ---- Class: Simrel ----
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
      n          = 100,
      p          = 20,
      q          = 10,
      m          = 1,
      relpos     = c(1, 2, 3),
      ypos       = NULL,
      gamma      = 0.8,
      eta        = 0.3,
      lambda.min = 1e-5,
      rho        = NULL,
      R2         = 0.9,
      ntest      = NULL,
      muX        = NULL,
      muY        = NULL,
      type       = "univariate"
    ),
    ..properties = list(
      relpred    = NULL,
      eigen_x    = NULL,
      eigen_w    = NULL,
      sigma_z    = NULL,
      sigma_zinv = NULL,
      sigma_y    = NULL,
      sigma_w    = NULL,
      sigma_zy   = NULL,
      sigma_zw   = NULL,
      sigma      = NULL,
      rotation_x = NULL,
      rotation_y = NULL,
      beta_z     = NULL,
      beta       = NULL,
      beta0      = NULL,
      Rsq_y      = NULL,
      Rsq_w      = NULL,
      minerror   = NULL
    ),
    ..get_cov = function(pos, Rsq, eta, p, lambda){
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
    ..is_pd = function(sigma_mat){
      all(eigen(sigma_mat)$values > 0)
    },
    ..get_eigen = function(gamma, p) {
      return((exp(-gamma * (1:p))) / (exp(-gamma)))
    },
    ..predpos = function(p, q, relpos) {
      relpos_list <- if (!is.list(relpos)) list(relpos) else relpos
      irrelpos <- setdiff(seq_len(p), Reduce(union, relpos_list))
      out <- lapply(seq_along(relpos_list), function(i){
        pos      <- relpos_list[[i]]
        ret      <- c(pos, sample(irrelpos, q[i] - length(pos)))
        irrelpos <<- setdiff(irrelpos, ret)
        return(ret)
      })
      if(!is.list(relpos)) out <- unlist(out)
      return(out)
    },
    ..get_data = function(n, p, sigma, rotation_x, m = 1,
                          rotation_y = NULL, muX = NULL, muY = NULL){
      sigma_rot <- chol(sigma)
      train_cal <- matrix(rnorm(n * (p + m), 0, 1), nrow = n) %*% sigma_rot
      Z <- train_cal[, (m + 1):(m + p), drop = F]
      X <- Z %*% t(rotation_x)
      W <- train_cal[, 1:m, drop = F]
      Y <- if (all(!is.null(m), m > 1)) W %*% t(rotation_y) else W

      if (!is.null(muX)) X <- sweep(X, 2, muX, "+")
      if (!is.null(muY)) Y <- sweep(Y, 2, muY, "+")

      list(y = unname(Y), x = X)
    }
  ),
  active = list(
    list_properties = function(){
      properties <- private$..properties
      properties <- properties[!sapply(properties, is.null)]
      return(properties)
    },
    list_parameters = function(){
      parameters <- private$..parameters
      parameters <- parameters[!sapply(parameters, is.null)]
      return(parameters)
    }
  ),
  public = list(
    initialize = function(...) {
      self$set_parameters(...)
    },
    base_properties = function() {
      ## Adding Properties to Simrel Object
      self$set_properties("relpred", expression({
        p      <- self$get_parameters("p")
        q      <- self$get_parameters("q")
        relpos <- self$get_parameters("relpos")
        out <- private$..predpos(p, q, relpos)
        if (!is.list(relpos)) out <- unlist(out)
        return(out)
      }))
      self$set_properties("eigen_x", expression({
        p     <- self$get_parameters("p")
        gamma <- self$get_parameters("gamma")
        private$..get_eigen(gamma, p)
      }))
      self$set_properties("eigen_w", expression({
        l <- length(self$get_parameters("ypos"))
        if (l == 0) return(1)
        eta <- self$get_parameters("eta")
        private$..get_eigen(eta, l)
      }))
      self$set_properties("sigma_z", expression({
        diag(self$get_properties("eigen_x"))
      }))
      self$set_properties("sigma_zinv", expression({
        diag(1 / self$get_properties("eigen_x"))
      }))
    },
    get_validated = function() {
      errors <- list()
      warnings <- list()
      params <- self$list_parameters
      list2env(params, envir = environment())
      switch(
        type,
        univariate = {
          req_params <- c("n", "p", "q", "relpos", "R2", "gamma", "type")
          ## --- ensure minimum required parameters
          if (!all(req_params %in% names(params))) {
            stop("Not Enough Parameters to continue")
          }

          ## --- Check muX and muY
          if (exists("muX")) {
            if (length(muX) != p) errors[[1]] <- "Length of muX is not equals to number of predictors (p)"
          }
          if (exists("muY")) {
            if (length(muY) != m) errors[[1]] <- "Length of muY is not equals to number of responses (m)"
          }

          ## --- ensure length params agree
          if (q <= length(relpos)) {
            errors[[2]] <- "Number of position of relevant components should be less than the number of relevant predictors."
          }
          if (q > p) {
            errors[[3]] <- "Relevant number of predictor should be less than or equals to total number of predictors."
          }
          if (!all(c(length(q), length(gamma), length(R2)) == 1)) {
            errors[[6]] <- "Length of q, gamma and R2 should be equal to 1."
          }
          if (!all(relpos %in% 1:p)) {
            errors[[7]] <- paste("Position of relevant component must be between 1 and ", p)
          }
          ## --- Ensure range of parameters
          if (gamma < 0) {
            errors[[4]] <- "Gamma must be greater than zero."
          }
          if (!all(R2 > 0, R2 < 1)) {
            errors[[5]] <- "Coefficient of determination (R2) should lie between 0 and 1."
          }
          ## --- Ensure class of parameters
          if (!is.atomic(relpos)) {
            errors[[6]] <- "Relpos for uniresponse simulation should be a vector."
          }
        },
        multivariate = {
          ## --- Ensure minimum required parameters
          req_params <- c("n", "p", "q", "relpos", "R2",
                          "gamma", "type", "ypos", "m")
          if (any(!req_params %in% names(params))) {
            stop("Not Enough Parameters to continue")
          }

          ## --- Check muX and muY
          if (exists("muX")) {
            if (length(muX) != p) errors[[1]] <- "Length of muX is not equals to number of predictors (p)"
          }
          if (exists("muY")) {
            if (length(muY) != m) errors[[1]] <- "Length of muY is not equals to number of responses (m)"
          }

          ## --- Ensure class of parameters
          if (!is.list(relpos) | !is.list(ypos)) {
            stop("relpos and ypos must be a list")
          }

          if (!all(c(length(q), length(relpos), length(R2)) == length(q))) {
            stop(paste("Length of q, relpos and R2 should be of same length"))
          }

          ## Ensure Lengths of parameters
          if (any(q < sapply(relpos, length))) {
            errors[[2]] <- "Number of position of relevant components should be less than the number of relevant predictors for each response components."
          }
          if (sum(q) > p) {
            errors[[3]] <- "Total number of relevant predictors should be less than or equals to total number of predictors."
          }
          if (!all(unlist(relpos) %in% 1:p)) {
            errors[[4]] <- paste("Position of relevant predictor component of response must be between 1 and ", p)
          }
          if (any(duplicated(unlist(relpos)))) {
            errors[[5]] <- "Relevant position of predictor components must be unique"
          }
          if (!all(unlist(ypos) %in% 1:m)) {
            errors[[6]] <- paste("Position of response component must be between 1 and ", m)
          }
          if (any(duplicated(unlist(ypos)))) {
            errors[[7]] <- "Position of response components must be unique"
          }
          if (length(unlist(ypos)) != m) {
            errors[[8]] <- paste("Position of response components (ypos) must contain all response variabels from 1 to ", m)
          }
          ## --- Ensure range of parameters
          if (gamma < 0) {
            errors[[9]] <- "gamma must be greater than zero."
          }
          if (eta < 0) {
            errors[[9]] <- "eta must be greater than zero."
          }
          if (any(sapply(R2, function(x) !all(x > 0, x < 1)))) {
            errors[[10]] <- "Coefficients of determination (R2) should lie between 0 and 1."
          }

        }
      )
      out <- list(errors = errors, warnings = warnings)
      return(out)
    },
    get_parameters = function(key) {
      if (length(key) == 1) return(private$..parameters[[key]])
      else return(private$..parameters[key])
    },
    set_parameters = function(...) {
      params <- list(...)
      for (key in names(params)) {
        private$..parameters[[key]] <- params[[key]]
      }
    },
    get_properties = function(key) {
      if (length(key) == 1) return(private$..properties[[key]])
      else return(private$..properties[key])
    },
    set_properties = function(key, expression) {
      private$..properties[[key]] <- eval(expression)
    },
    get_data = function(type = c("train", "test", "both")){
      type <- match.arg(type)
      n <- self$get_parameters("n")
      p <- self$get_parameters("p")
      m <- self$get_parameters("m")
      sigma <- self$get_properties("sigma")
      rotation_x <- self$get_properties("rotation_x")
      rotation_y <- self$get_properties("rotation_y")
      muX <- self$get_parameters("muX")
      muY <- self$get_parameters("muY")

      train <- private$..get_data(n, p, sigma, rotation_x, m, rotation_y, muX, muY)
      out <- data.frame(y = I(train$y), x = I(train$x))
      if (type %in% c("test", "both")) {
        ntest <- self$get_parameters("ntest")
        if(is.null(ntest)){
          warning("ntest is not supplied, so n is used for test data")
          ntest <- n
        }
        test <- private$..get_data(ntest, p, sigma, rotation_x, m, rotation_y, muX, muY)
        out <- list(
          train = data.frame(y = I(train$y), x = I(train$x)),
          test = data.frame(y = I(test$y), x = I(test$x))
        )
        return(out)
      }

      return(out)
    }
  )
)


## ---- Class: UniSimrel ----
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
  public = list(
    initialize = function(...){
      super$initialize(...)
      super$set_parameters(type = "univariate")

      ## Validate the parameters
      is_valid <- self$get_validated()
      is_valid <- lapply(is_valid, function(x) x[!sapply(x, is.null)])
      if (length(is_valid$errors)) stop(unlist(is_valid$errors))
      if (length(is_valid$warnings)) warning(unlist(is_valid$warnings))

      ## Adding Properties to Simrel Object
      self$compute_properties()
    },
    add_properties = function(){
      self$set_properties("sigma_y", expression({1}))
      self$set_properties("sigma_zy", expression({
        relpos <- self$get_parameters("relpos")
        R2 <- self$get_parameters("R2")
        p <- self$get_parameters("p")
        lambda <- self$get_properties("eigen_x")
        eta <- self$get_properties("eigen_w")
        private$..get_cov(pos = relpos, Rsq = R2, p = p, lambda = lambda, eta = eta)
      }))
      self$set_properties("sigma", expression({
        sigma_zy <- self$get_properties("sigma_zy")
        sigma_y <- self$get_properties("sigma_y")
        sigma_z <- self$get_properties("sigma_z")
        out <- rbind(c(sigma_y, t(sigma_zy)), cbind(sigma_zy,  sigma_z))
        unname(out)
      }))
      self$set_properties("rotation_x", expression({
        relpred <- self$get_properties("relpred")
        p <- self$get_parameters("p")
        irrelpred <- setdiff(1:p, unlist(relpred))
        out <- diag(p)
        out[irrelpred, irrelpred] <- private$..get_rotate(irrelpred)
        out[relpred, relpred] <- private$..get_rotate(relpred)
        return(out)
      }))
      self$set_properties("beta_z", expression({
        sigma_zinv <- self$get_properties("sigma_zinv")
        sigma_zy <- self$get_properties("sigma_zy")
        return(sigma_zinv %*% sigma_zy)
      }))
      self$set_properties("beta", expression({
        rotation_x <- self$get_properties("rotation_x")
        beta_z <- self$get_properties("beta_z")
        return(rotation_x %*% beta_z)
      }))
      self$set_properties("beta0", expression({
        beta0 <- rep(0, self$get_parameters("m"))
        muX <- self$get_parameters("muX")
        muY <- self$get_parameters("muY")
        if (!is.null(muX)) beta <- self$get_properties("beta")
        if (!is.null(muY)) beta0 <- beta0 + muY
        if (!is.null(muX)) beta0 <- beta0 - t(beta) %*% muX
        return(c(beta0))
      }))
      self$set_properties("Rsq_y", expression({
        beta_z <- self$get_properties("beta_z")
        sigma_zy <- self$get_properties("sigma_zy")
        c(t(beta_z) %*% sigma_zy)
      }))
      self$set_properties("minerror", expression({
        sigma_y <- self$get_properties("sigma_y")
        Rsq <- self$get_properties("Rsq_y")
        c(sigma_y - Rsq)
      }))
    },
    compute_properties = function() {
      super$base_properties()
      self$add_properties()
    }
  )
)


## ---- Class: BiSimrel ----
#' @title Bi-Response Simulation of Linear Model Data
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
BiSimrel <- R6::R6Class(
  "BiSimrel",
  inherit = Simrel,
  public = list(
    initialize = function(...){
      params <- list(...)
      default_params <- list(
        q = c(6, 7, 3),
        m = 2,
        relpos = list(c(1, 2), c(2, 3, 5)),
        ypos = list(1, c(2, 3)),
        R2 = c(0.7, 0.9),
        rho = c(0.7, 0.8),
        eta = 0,
        type = "bivariate"
      )
      default_params[names(params)] <- params
      do.call(super$initialize, default_params)

      ## Validate the parameters
      is_valid <- self$get_validated()
      is_valid <- lapply(is_valid, function(x) x[!sapply(x, is.null)])
      if (length(is_valid$errors)) stop(unlist(is_valid$errors))
      if (length(is_valid$warnings)) warning(unlist(is_valid$warnings))

      ## Adding Properties to Simrel Object
      self$compute_properties()
    },
    add_properties = function() {
      ## Adding Properties to Simrel Object
      self$set_properties("sigma_w", expression({
        eigen_w <- self$get_properties("eigen_w")
        m <- self$get_parameters("m")
        diag(c(eigen_w, rep(1, m - length(eigen_w))))
      }))
      self$set_properties("sigma_zw", expression({
        relpos <- self$get_parameters("relpos")
        m <- self$get_parameters("m")
        R2 <- self$get_parameters("R2")
        p <- self$get_parameters("p")
        lambda <- self$get_properties("eigen_x")
        eta <- self$get_properties("eigen_w")
        out <- mapply(private$..get_cov, pos = relpos, Rsq = R2, eta = eta,
                      MoreArgs = list(p = p, lambda = lambda))
        cbind(out, rep(0, m - length(eta)))
      }))
      self$set_properties("sigma", expression({
        sigma_zw <- self$get_properties("sigma_zw")
        sigma_w <- self$get_properties("sigma_w")
        sigma_z <- self$get_properties("sigma_z")
        out <- cbind(rbind(sigma_w, sigma_zw), rbind(t(sigma_zw), sigma_z))
        unname(out)
      }))
      self$set_properties("rotation_x", expression({
        relpred <- self$get_properties("relpred")
        p <- self$get_parameters("p")
        irrelpred <- setdiff(1:p, unlist(relpred))
        out <- diag(p)
        out[irrelpred, irrelpred] <- private$..get_rotate(irrelpred)
        for (pos in relpred) {
          rotMat         <- private$..get_rotate(pos)
          out[pos, pos] <- rotMat
        }
        return(out)
      }))
      self$set_properties("rotation_y", expression({
        ypos <- self$get_parameters("ypos")
        m <- self$get_parameters("m")
        out <- diag(m)
        for (pos in ypos) {
          rotMat         <- private$..get_rotate(pos)
          out[pos, pos]  <- rotMat
        }
        return(out)
      }))
      self$set_properties("sigma_y", expression({
        rotation_y <- self$get_properties("rotation_y")
        sigma_w <- self$get_properties("sigma_w")
        t(rotation_y) %*% sigma_w %*% rotation_y
      }))
      self$set_properties("sigma_zy", expression({
        rotation_y <- self$get_properties("rotation_y")
        sigma_zw <- self$get_properties("sigma_zw")
        rotation_y %*% t(sigma_zw)
      }))
      self$set_properties("sigma_xy", expression({
        rotation_x <- self$get_properties("rotation_x")
        rotation_y <- self$get_properties("rotation_y")
        sigma_zw <- self$get_properties("sigma_zw")
        rotation_y %*% t(sigma_zw) %*% t(rotation_x)
      }))
      self$set_properties("beta_z", expression({
        sigma_zinv <- self$get_properties("sigma_zinv")
        sigma_zw <- self$get_properties("sigma_zw")
        return(sigma_zinv %*% sigma_zw)
      }))
      self$set_properties("beta", expression({
        rotation_x <- self$get_properties("rotation_x")
        rotation_y <- self$get_properties("rotation_y")
        beta_z <- self$get_properties("beta_z")
        return(rotation_x %*% beta_z %*% t(rotation_y))
      }))
      self$set_properties("beta0", expression({
        beta0 <- rep(0, self$get_parameters("m"))
        muX <- self$get_parameters("muX")
        muY <- self$get_parameters("muY")
        if (!is.null(muX)) beta <- self$get_properties("beta")
        if (!is.null(muY)) beta0 <- beta0 + muY
        if (!is.null(muX)) beta0 <- beta0 - t(beta) %*% muX
        return(c(beta0))
      }))
      self$set_properties("Rsq_w", expression({
        beta_z <- self$get_properties("beta_z")
        sigma_zw <- self$get_properties("sigma_zw")
        sigma_w <- self$get_properties("sigma_w")
        unname(t(beta_z) %*% sigma_zw %*% solve(sigma_w))
      }))
      self$set_properties("Rsq_y", expression({
        rotation_y <- self$get_properties("rotation_y")
        Rsq_w <- self$get_properties("Rsq_w")
        t(rotation_y) %*% Rsq_w %*% rotation_y
      }))
      self$set_properties("minerror", expression({
        sigma_y <- self$get_properties("sigma_y")
        Rsq <- self$get_properties("Rsq_y")
        unname(sigma_y - Rsq)
      }))
    },
    compute_properties = function() {
      super$base_properties()
      self$add_properties()
    }
  )
)


## ---- Class: MultiSimrel ----
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
  public = list(
    initialize = function(...){
      params <- list(...)
      default_params <- list(
        q = c(6, 7),
        m = 3,
        relpos = list(c(1, 2), c(3, 4, 5)),
        ypos = list(1, c(2, 3)),
        R2 = c(0.7, 0.9),
        eta = 0,
        type = "multivariate"
      )
      default_params[names(params)] <- params
      do.call(super$initialize, default_params)

      ## Validate the parameters
      is_valid <- self$get_validated()
      is_valid <- lapply(is_valid, function(x) x[!sapply(x, is.null)])
      if (length(is_valid$errors)) stop(unlist(is_valid$errors))
      if (length(is_valid$warnings)) warning(unlist(is_valid$warnings))

      ## Adding Properties to Simrel Object
      self$compute_properties()
    },
    add_properties = function() {
      ## Adding Properties to Simrel Object
      self$set_properties("sigma_w", expression({
        eigen_w <- self$get_properties("eigen_w")
        m <- self$get_parameters("m")
        diag(c(eigen_w, rep(1, m - length(eigen_w))))
      }))
      self$set_properties("sigma_zw", expression({
        relpos <- self$get_parameters("relpos")
        m <- self$get_parameters("m")
        R2 <- self$get_parameters("R2")
        p <- self$get_parameters("p")
        lambda <- self$get_properties("eigen_x")
        eta <- self$get_properties("eigen_w")
        out <- mapply(private$..get_cov, pos = relpos, Rsq = R2, eta = eta,
                      MoreArgs = list(p = p, lambda = lambda))
        cbind(out, rep(0, m - length(eta)))
      }))
      self$set_properties("sigma", expression({
        sigma_zw <- self$get_properties("sigma_zw")
        sigma_w <- self$get_properties("sigma_w")
        sigma_z <- self$get_properties("sigma_z")
        out <- cbind(rbind(sigma_w, sigma_zw), rbind(t(sigma_zw), sigma_z))
        unname(out)
      }))
      self$set_properties("rotation_x", expression({
        relpred <- self$get_properties("relpred")
        p <- self$get_parameters("p")
        irrelpred <- setdiff(1:p, unlist(relpred))
        out <- diag(p)
        out[irrelpred, irrelpred] <- private$..get_rotate(irrelpred)
        for (pos in relpred) {
          rotMat         <- private$..get_rotate(pos)
          out[pos, pos] <- rotMat
        }
        return(out)
      }))
      self$set_properties("rotation_y", expression({
        ypos <- self$get_parameters("ypos")
        m <- self$get_parameters("m")
        out <- diag(m)
        for (pos in ypos) {
          rotMat         <- private$..get_rotate(pos)
          out[pos, pos]  <- rotMat
        }
        return(out)
      }))
      self$set_properties("sigma_y", expression({
        rotation_y <- self$get_properties("rotation_y")
        sigma_w <- self$get_properties("sigma_w")
        t(rotation_y) %*% sigma_w %*% rotation_y
      }))
      self$set_properties("sigma_zy", expression({
        rotation_y <- self$get_properties("rotation_y")
        sigma_zw <- self$get_properties("sigma_zw")
        rotation_y %*% t(sigma_zw)
      }))
      self$set_properties("sigma_xy", expression({
        rotation_x <- self$get_properties("rotation_x")
        rotation_y <- self$get_properties("rotation_y")
        sigma_zw <- self$get_properties("sigma_zw")
        rotation_y %*% t(sigma_zw) %*% t(rotation_x)
      }))
      self$set_properties("beta_z", expression({
        sigma_zinv <- self$get_properties("sigma_zinv")
        sigma_zw <- self$get_properties("sigma_zw")
        return(sigma_zinv %*% sigma_zw)
      }))
      self$set_properties("beta", expression({
        rotation_x <- self$get_properties("rotation_x")
        rotation_y <- self$get_properties("rotation_y")
        beta_z <- self$get_properties("beta_z")
        return(rotation_x %*% beta_z %*% t(rotation_y))
      }))
      self$set_properties("beta0", expression({
        beta0 <- rep(0, self$get_parameters("m"))
        muX <- self$get_parameters("muX")
        muY <- self$get_parameters("muY")
        if (!is.null(muX)) beta <- self$get_properties("beta")
        if (!is.null(muY)) beta0 <- beta0 + muY
        if (!is.null(muX)) beta0 <- beta0 - t(beta) %*% muX
        return(c(beta0))
      }))
      self$set_properties("Rsq_w", expression({
        beta_z <- self$get_properties("beta_z")
        sigma_zw <- self$get_properties("sigma_zw")
        sigma_w <- self$get_properties("sigma_w")
        unname(t(beta_z) %*% sigma_zw %*% solve(sigma_w))
      }))
      self$set_properties("Rsq_y", expression({
        rotation_y <- self$get_properties("rotation_y")
        Rsq_w <- self$get_properties("Rsq_w")
        t(rotation_y) %*% Rsq_w %*% rotation_y
      }))
      self$set_properties("minerror", expression({
        sigma_y <- self$get_properties("sigma_y")
        Rsq <- self$get_properties("Rsq_y")
        unname(sigma_y - Rsq)
      }))
    },
    compute_properties = function() {
      super$base_properties()
      self$add_properties()
    }
  )
)



## ---- Wrapper Function ----
#' @title Multi-Response Simulation of Linear Model Data
#' @description Simulates multivariate linear model data where users can specify few parameters for simulating data with wide range of properties.
#' @import R6
#' @param  type Number of training samples
#' @param ... All required arguments for different types of simulation
#' @export
simulater <- function(type, ...){
  if (type == "univariate") {
    sobj <- UniSimrel$new(...)
    return(sobj)
  }
  if (type == "multivariate") {
    sobj <- MultiSimrel$new(...)
    return(sobj)
  }
  return(paste(type, "is unknown"))
}
