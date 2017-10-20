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
      eigen_y    = NULL,
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
      rotate_y  <- all(exists("rotation_y"), !is.null(rotation_y))
      sigma_rot <- chol(sigma)
      train_cal <- matrix(rnorm(n * (p + m), 0, 1), nrow = n) %*% sigma_rot
      Z         <- train_cal[, (m + 1):(m + p), drop = F]
      X         <- Z %*% t(rotation_x)
      W         <- train_cal[, 1:m, drop = F]
      Y         <- if (rotate_y) W %*% t(rotation_y) else W

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
      self$set_properties("eigen_y", expression({
        type <- self$get_parameters("type")
        l <- switch(type,
               univariate = 1,
               bivariate = 2,
               multivariate = length(self$get_parameters("ypos")))
        if (l == 1) return(1)
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
        bivariate = {
          n_relpos <- sapply(relpos, length)

          ## --- Ensure minimum required parameters
          req_params <- c("n", "p", "q", "relpos", "R2",
                          "gamma", "type", "rho")
          if (any(!req_params %in% names(params))) {
            stop("Not Enough Parameters to continue")
          }

          ## --- Check muX and muY
          if (exists("muX")) {
            if (length(muX) != p)
              errors[[1]] <- "Length of muX is not equals to number of predictors (p)"
          }
          if (exists("muY")) {
            if (length(muY) != m)
              errors[[2]] <- "Length of muY is not equals to number of responses (m)"
          }

          ## --- Ensure class of parameters
          if (!is.list(relpos)) {
            stop("Position of relevant predictors must be a list")
          }

          compos <- sum(relpos[[1]] %in% relpos[[2]])
          if (compos > 0) {
            cpos1 <- which(relpos[[1]] %in% relpos[[2]])
            cpos2 <- which(relpos[[2]] %in% relpos[[1]])
          }
          if (all(compos == 0, q[3] > 0)) {
            warning("The number of common relevant predictors is set to zero due to orthogonal relevant spaces")
            private$..parameters$q[3] <- 0
          }
          if (all(compos == n_relpos[1], compos == n_relpos[2], q[1] != 0)) {
            private$..parameters$q[3] <- min(q[1:2])
            warnings(paste("The number of common relevant predictors is set to",
                           q[3], "since the relevant spaces are identical.\n"))
          }
          if (all(compos > 0, q[3] == 0)) {
            private$..parameters$q[3] <- compos
            warning(paste("The number of common relevant predictors is set to ",
                  compos, " since the relevant spaces must share common variable(s)\n"))
          }
          if (compos == 0) {
            rho[2] <- rho[1]/sqrt((1 - R2[1]) * (1 - R2[2]))
            if (rho[2] < -0.99999 | rho[2] > 0.9999) {
              stop("Non-PD covariance matrix, choose smaller R2-values\n")
            } else {
              warning(paste("The conditional correlation between y1 and y2 is fixed equal to",
                            round(rho[2], 2), " due to orthogonal relevant spaces \n"))
              private$..parameters$q[3] <- 0
              private$..parameters$rho[2] <- rho[2]
            }
          }

          if (!all(c(length(rho), length(relpos), length(R2)) == 2)) {
            stop(paste("Length of rho, relpos and R2 should be of length 2"))
          }

          ## Ensure Lengths of parameters
          if (length(unlist(sapply(relpos, length))) == 0) {
            errors[[3]] <- "You must have minimum one relevant component for each response."
          }
          if (sum(q) < sum(sapply(relpos, length))) {
            errors[[4]] <- "Number of position of relevant components should be less than the number of relevant predictors for each response components."
          }
          if (sum(q) > p) {
            errors[[4]] <- "Total number of relevant predictors should be less than or equals to total number of predictors."
          }
          if (!all(unlist(relpos) %in% 1:p)) {
            errors[[5]] <- paste("Position of relevant predictor component of response must be between 1 and ", p)
          }

          ## --- Ensure range of parameters
          if (gamma < 0) {
            errors[[7]] <- "gamma must be greater than zero."
          }
          if (any(sapply(R2, function(x) !all(x > 0, x < 1)))) {
            errors[[8]] <- "Coefficients of determination (R2) should lie between 0 and 1."
          }
          if (any(sapply(rho, function(x) !all(x > -1, x < 1)))) {
            errors[[9]] <- "The correlation between responses should lie between -1 and 1."
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
      sigma <- self$get_properties("sigma")
      if (!private$..is_pd(sigma)) stop("\nCorrelation matrix is not positive definite.\n")
      type <- match.arg(type)
      n <- self$get_parameters("n")
      p <- self$get_parameters("p")
      m <- self$get_parameters("m")
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
        eta <- self$get_properties("eigen_y")
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
  private = list(
    ..aux_var = list(
      comn_relpos     = NULL,
      n_compos        = NULL,
      ex_relpos       = NULL,
      n_relpos        = NULL,
      comn_relpos_idx = NULL
    ),
    ..Rfunc = function(alpha1, alpha2, L) {
      t(alpha1) %*% solve(diag(L)) %*% alpha2
    },
    ..a21func = function(a11, a12, a22, R12, l1, l2) {
      l1 / a11 * (R12 - a22 * a12 / l2)
    },
    ..a22func = function(a11, a12, R2, R12, l1, l2) {
      bb    <- R12 * a12 / a11^2 * l1 / l2
      root  <- sqrt(R12^2 * a12^2 / a11^4 * l1^2 / l2^2 -
                    (1 / l2 + l1 / l2 * a12^2 / a11^2) *
                    (l1 / a11^2 * R12^2 - R2))
      denom <- (1 / l2 + l1 / l2 * a12^2 / a11^2)
      w1    <- (bb - root) / denom
      w2    <- (bb + root) / denom
      return(c(w1, w2))
    }
  ),
  public = list(
    initialize = function(...){
      params <- list(...)
      default_params <- list(
        q = c(6, 7, 3),
        m = 2,
        relpos = list(c(1, 2), c(2, 3, 5)),
        R2 = c(0.7, 0.9),
        rho = c(0.7, 0.8),
        eta = 0,
        type = "bivariate"
      )
      if (length(params) != 0) default_params[names(params)] <- params
      do.call(super$initialize, default_params)

      ## Validate the parameters
      is_valid <- self$get_validated()
      is_valid <- lapply(is_valid, function(x) x[!sapply(x, is.null)])
      if (length(is_valid$errors)) stop(unlist(is_valid$errors))
      if (length(is_valid$warnings)) warning(unlist(is_valid$warnings))

      ## Compute auxiliary variables
      self$compute_aux_vars()


      ## Adding Properties to Simrel Object
      self$compute_properties()
    },
    get_aux_var = function(key) private$..aux_var[[key]],
    compute_aux_vars = function() {
      relpos            <- self$get_parameters("relpos")
      private$..aux_var <- within(private$..aux_var, {
        comn_relpos     <- do.call(intersect, relpos)
        n_compos        <- length(comn_relpos)
        ex_relpos       <- append(lapply(relpos, setdiff, comn_relpos),
                                  list(comn_relpos))
        n_relpos        <- sapply(relpos, length)
        comn_relpos_idx <- lapply(relpos, function(x) which(x %in% comn_relpos))
      })
    },
    add_properties = function() {
      ## Adding Properties to Simrel Object
      self$set_properties("relpred", expression({
        p           <- self$get_parameters("p")
        q           <- self$get_parameters("q")
        relpos      <- self$get_parameters("relpos")
        comn_relpos <- self$get_aux_var("comn_relpos")
        ex_relpos   <- self$get_aux_var("ex_relpos")
        out         <- private$..predpos(p, q, ex_relpos)
        if (!is.list(relpos)) out <- unlist(out)
        return(out)
      }))
      self$set_properties("sigma_zy", expression({
        p        <- self$get_parameters("p")
        m        <- self$get_parameters("m")
        relpos   <- self$get_parameters("relpos")
        R2       <- self$get_parameters("R2")
        rho      <- self$get_parameters("rho")
        p        <- self$get_parameters("p")
        lambda   <- self$get_properties("eigen_x")
        eta      <- self$get_properties("eigen_y")
        sigma_zy <- matrix(0, p, 2)
        n_compos <- self$get_aux_var("n_compos")
        n_relpos <- self$get_aux_var("n_relpos")
        cpos     <- self$get_aux_var("comn_relpos_idx")
        if (n_compos == 0) {
          out <- mapply(private$..get_cov, pos = relpos, Rsq = R2, eta = eta,
                        MoreArgs = list(p = p, lambda = lambda))
          R12     <- sum(out[out[, 1] != 0, 1] * out[out[, 2] != 0, 2]/ lambda[relpos[[1]]])
          rho_out <- rho[1] / sqrt(prod(1 - R2 ^ 2))
          if (any(rho_out < -1, rho_out > 1)) {
            stop(paste("Two Response in orthogonal, but highly relevant spaces must be less correlated.",
                       "Choose it less than", rho[1], "or reduce the R2 value."))
          } else {
            rho_out <- c(rho[1], rho_out)
            self$set_parameters(rho = rho_out)
          }
          sigma_zy <- out
        }
        if (any(n_compos == n_relpos)) {
          if (any(n_relpos > n_compos)) {
            id <- which(n_relpos == n_compos)
            id2 <- (1:2)[-id]
            private$..properties$relpred[id] <- vector("list", 1)
            warning(paste("\nThe number of relevant predictors for response", id, "is set to",
                          n_compos, "since its relevant space is a subspace of the other\n"))

            cpos.1 <- cpos[[id]]
            cpos.2 <- cpos[[id2]]
            rho    <- self$get_parameters("rho")
            R2     <- self$get_parameters("R2")
            R12    <- rho[1] - rho[2] * sqrt((1 - R2[1]) * (1 - R2[2]))
            alpha  <- runif(n_compos, -1, 1)
            alpha  <- sign(alpha) * sqrt(R2[id] * abs(alpha)/sum(abs(alpha)) * lambda[relpos[[id]]])

            alpha2sum <- 1
            j <- 1
            while (alpha2sum > R2[id2] & j < 1000) {
              alpha1     <- alpha
              alpha2star <- runif(n_compos - 1, -1, 1)
              alphasum   <- sum(alpha1[-1] * alpha2star/lambda[relpos[[id]][cpos.1[-1]]])
              a21        <- (R12 - alphasum) * lambda[relpos[[id]][cpos.1[1]]]/alpha1[1]
              alpha2star <- c(a21, alpha2star)
              alpha2sum  <- private$..Rfunc(alpha2star, alpha2star, lambda[relpos[[id]]])
              j          <- j + 1
            }
            if (j == 1000) {
              stop("No PD matrix found in 1000 simulations with current parameter setting")
            }
            if ((n_relpos[id2] - n_compos) > 1) {
              alpha2rest <- runif(n_relpos[id2] - n_compos, 0, 1) *
                sample(c(-1, 1), n_relpos[id2] - n_compos, replace = TRUE)
              alpha2rest <- sign(alpha2rest) * sqrt((R2[id2] - c(alpha2sum)) * abs(alpha2rest)/sum(abs(alpha2rest)) * lambda[relpos[[id2]][-cpos.2]])
            } else {
              alpha2rest <- sample(c(-1, 1), 1, replace = FALSE) *
                sqrt((R2[id2] - c(alpha2sum)) * lambda[relpos[[id2]][-cpos.2]])
            }
            alpha_out <- rep(0, p)
            alpha_out[relpos[[id2]][cpos.2]] <- alpha2star
            alpha_out[relpos[[id2]][-cpos.2]] <- alpha2rest

            sigma_zy[relpos[[id]], id] <- alpha1
            sigma_zy[, id2] <- alpha_out
          } else {
            if (n_compos == 1) {
              alpha1 <- sqrt(R2[1] * lambda[relpos[[1]]])
              alpha2 <- sign(rho[1]) * sqrt(R2[2] * lambda[relpos[[2]]])
              rho[2] <- (rho[1] - alpha1 * alpha2/lambda[relpos[[1]]])/
                sqrt((1 - R2[1]) * (1 - R2[2]))
              private$..properties$rho[2] <- rho[2]
              warning(paste(
                "\nThe conditional correlation between y1 and y2 is fixed equal to",
                round(rho[2], 2),
                "due to overlapping and one-dimensional relevant spaces for y1 and y2\n"
              ))
              sigma_zy[relpos[[1]], 1] <- alpha1
              sigma_zy[relpos[[2]], 2] <- alpha2
            }
            if (n_compos == 2) {
              R12 <- rho[1] - rho[2] * sqrt((1 - R2[1]) * (1 - R2[2]))

              R2try <- 2
              tol <- 1
              j <- 1
              while (tol > 0.001 & (j < 5000)) {
                alpha1 <- runif(2, -1, 1)
                alpha1 <- sign(alpha1) *
                  sqrt(R2[1] * abs(alpha1)/
                       sum(abs(alpha1)) * lambda[relpos[[1]]])

                options(warn = -1)
                a22 <- private$..a22func(alpha1[1], alpha1[2], R2[2],
                                R12, lambda[relpos[[2]]][1],
                                lambda[relpos[[2]]][2])
                a22 <- a22[which(abs(a22) == min(abs(a22)))][1]
                if (is.na(a22)) next
                options(warn = 0)

                a21 <- private$..a21func(alpha1[1], alpha1[2], a22, R12,
                                lambda[relpos[[2]]][1],
                                lambda[relpos[[2]]][2])
                alpha2 <- c(a21, a22)
                R2try <- private$..Rfunc(alpha2, alpha2, lambda[relpos[[2]]])
                tol <- abs(R2try - R2[2])
                j <- j + 1
              }
              if (j == 5000) {
                stop("No PD covariance matrix found with current prameter setting\n ")
              }
              sigma_zy[relpos[[1]], 1] <- alpha1
              sigma_zy[relpos[[2]], 2] <- alpha2
            }
            if (n_compos > 2) {
              R12 <- rho[1] - rho[2] * sqrt((1 - R2[1]) * (1 - R2[2]))

              R2try <- 2
              tol   <- 1
              j     <- 1
              while (tol > 0.001 & (j < 5000)) {
                ## Sampling all but two positions
                R2rest.1 <- 1
                R2rest.2 <- 1
                R12rest  <- 1
                k        <- 1
                R1prop   <- runif(1, 0.3, 0.8)
                R2prop   <- runif(1, 0.3, 0.8)
                loop_condition <- any(
                  R2rest.1 > (n_relpos[1] - 2)/n_relpos[1] * R2[1],
                  R2rest.2 > (n_relpos[2] - 2)/n_relpos[2] * R2[2],
                  rho[1] - R12rest < -1,
                  rho[1] - R12rest > 1
                )
                while (loop_condition & (k < 1000)) {
                  alpha1cp <- runif(n_relpos[1] - 2, -1, 1)
                  alpha1cp <- sign(alpha1cp) *
                    sqrt(R1prop * R2[1] *
                         abs(alpha1cp) / sum(abs(alpha1cp)) *
                         lambda[relpos[[1]][-c(1:2)]])
                  R2rest.1 <- sum(alpha1cp * alpha1cp / lambda[relpos[[1]][-c(1:2)]])
                  alpha2cp <- runif(n_relpos[2] - 2, -1, 1)
                  alpha2cp <- sign(alpha2cp) *
                    sqrt(R2prop * R2[2] * abs(alpha2cp) / sum(abs(alpha2cp)) *
                         lambda[relpos[[2]][-c(1:2)]])
                  R2rest.2 <- sum(alpha2cp * alpha2cp / lambda[relpos[[2]][-c(1:2)]])
                  R12rest <- sum(alpha1cp * alpha2cp / lambda[relpos[[2]][-c(1:2)]])
                  k <- k + 1
                }
                alpha1 <- runif(2, -1, 1)
                alpha1 <- sign(alpha1) *
                  sqrt((R2[1] - R2rest.1) *
                       abs(alpha1)/sum(abs(alpha1)) *
                       lambda[relpos[[1]]][1:2])

                options(warn = -1)
                a22 <- private$..a22func(
                  alpha1[1], alpha1[2],
                  R2[2] - R2rest.2,
                  R12 - R12rest,
                  lambda[relpos[[2]]][1],
                  lambda[relpos[[2]]][2]
                )
                a22 <- a22[which(abs(a22) == min(abs(a22)))][1]
                if (is.na(a22)) next
                options(warn = 0)

                a21 <- private$..a21func(
                  alpha1[1], alpha1[2], a22,
                  R12 - R12rest, lambda[relpos[[2]]][1],
                  lambda[relpos[[2]]][2])
                alpha2 <- c(a21, a22, alpha2cp)
                alpha1 <- c(alpha1, alpha1cp)
                R2try  <- private$..Rfunc(alpha2, alpha2, lambda[relpos[[2]]])
                tol    <- abs(R2try - R2[2])
                j      <- j + 1
              }
              if (j == 5000) {
                stop("No PD covariance matrix found with current prameter setting\n ")
              }
              sigma_zy[relpos[[1]], 1] <- alpha1
              sigma_zy[relpos[[2]], 2] <- alpha2
            }
          }
        }
        if (all(n_compos != 0, all(n_compos < n_relpos))) {
          j <- 1
          cp1sum <- cp2sum <- 1
          while (any(cp1sum >= R2[1], cp2sum >= R2[2]) & (j < 1000)) {
            ## cp <- which(relpos[[1]] %in% relpos[[2]])
            cp       <- relpos[[1]][cpos[[1]]]
            R12      <- rho[1] - rho[2] * sqrt((1 - R2[1]) * (1 - R2[2]))
            alpha1cp <- sign(R12) * runif(n_compos, 0, 1)
            alpha2cp <- runif(n_compos, 0, 1)
            cp.sum   <- sum(alpha1cp * alpha2cp/lambda[cp])
            k        <- R12/cp.sum

            sigma_zy[cp, 1] <- alpha1cp <- sqrt(k) * alpha1cp
            sigma_zy[cp, 2] <- alpha2cp <- sqrt(k) * alpha2cp
            cp1sum <- sum(alpha1cp * alpha1cp/lambda[cp])
            cp2sum <- sum(alpha2cp * alpha2cp/lambda[cp])
            j <- j + 1
          }
          if (j == 1000) {
            stop("No PD covariance matrix found with current parameter setting\n")
          }

          id <- which(relpos[[1]] %in% cp)
          ncp <- relpos[[1]][-id]
          alpha1 <- runif(length(ncp), -1, 1)
          sigma_zy[ncp, 1] <- sign(alpha1) *
            sqrt((R2[1] - cp1sum) * abs(alpha1)/
                 sum(abs(alpha1)) * lambda[ncp])

          id <- which(relpos[[2]] %in% cp)
          ncp <- relpos[[2]][-id]
          alpha2 <- runif(length(ncp), -1, 1)
          sigma_zy[ncp, 2] <- alpha2 <- sign(alpha2) *
            sqrt((R2[2] - cp2sum) * abs(alpha2)/
                 sum(abs(alpha2)) * lambda[ncp])
        }
        return(sigma_zy)
      }))
      self$set_properties("sigma_y", expression({
        eigen_y <- self$get_properties("eigen_y")
        rho <- self$get_parameters("rho")
        out <- diag(c(eigen_y, rep(1, 2 - length(eigen_y))))
        out[out != diag(out)] <- rho[1]
        return(out)
      }))
      self$set_properties("sigma_yinv", expression({
        sigma_y <- self$get_properties("sigma_y")
        solve(sigma_y)
      }))
      self$set_properties("sigma", expression({
        sigma_zy <- self$get_properties("sigma_zy")
        sigma_y <- self$get_properties("sigma_y")
        sigma_z <- self$get_properties("sigma_z")
        out <- cbind(rbind(sigma_y, sigma_zy), rbind(t(sigma_zy), sigma_z))
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
      self$set_properties("sigma_xy", expression({
        rotation_x <- self$get_properties("rotation_x")
        sigma_zy <- self$get_properties("sigma_zy")
        t(sigma_zy) %*% t(rotation_x)
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
        t(sigma_zy) %*% beta_z
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
      if (length(params) != 0) default_params[names(params)] <- params
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
        eigen_y <- self$get_properties("eigen_y")
        m <- self$get_parameters("m")
        diag(c(eigen_y, rep(1, m - length(eigen_y))))
      }))
      self$set_properties("sigma_zw", expression({
        relpos <- self$get_parameters("relpos")
        m <- self$get_parameters("m")
        R2 <- self$get_parameters("R2")
        p <- self$get_parameters("p")
        lambda <- self$get_properties("eigen_x")
        eta <- self$get_properties("eigen_y")
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
    class(sobj) <- "simrel"
    return(sobj)
  }
  if (type == "multivariate") {
    sobj <- MultiSimrel$new(...)
    class(sobj) <- "simrel"
    return(sobj)
  }
  return(paste(type, "is unknown"))
}
