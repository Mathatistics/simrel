#' Simulation Plot with ggplot: The true beta, relevant component and eigen structure
#' @keywords simrel-plot, simrelplot, simulation plot, simulation ggplot
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom gridExtra grid.arrange
#' @param obj A simrel object
#' @param ncomp Number of components to plot
#' @param layout A layout matrix of how to layout multiple plots
#' @param print.cov Output estimated covariance structure
#' @param which A character indicating which plot you want as output, it can take \code{TrueBeta}, \code{RelComp} and \code{EstRelComp}
#' @return A list of plots
#' @examples
#' sim.obj <- simrel(n = 100, p = 16, q = c(3, 4, 5),
#'    relpos = list(c(1, 2), c(3, 4), c(5, 7)), m = 5,
#'    ypos = list(c(1, 4), 2, c(3, 5)), type = "multivariate",
#'    R2 = c(0.8, 0.7, 0.9), gamma = 0.8)
#'
#' simrelplot(sim.obj, layout = matrix(c(2, 1, 3, 1), 2))
#'
#' simrelplot(sim.obj, which = c(1, 2))
#'
#' simrelplot(sim.obj, which = c(1, 3), layout = matrix(c(1, 2), 1))
#' @export

simrelplot <- function(simrel_object, ncomp = NULL, which = 1L:3L,
                         layout = NULL, print.cov = FALSE) {
  obj <- simrel_object
  if (!("Simrel" %in% class(obj))) stop("Wrong Simrel Object")
  n <- obj$get_parameters("n")
  nx <- obj$get_parameters("p")
  ny <- obj$get_parameters("m")
  xticks <- 1:nx
  if (is.null(ncomp)) ncomp <- min(n, nx, 20)
  if (is.null(which)) which <- 1L:3L
  if (is.null(layout) & length(which) == 3) layout <- matrix(c(2, 1, 3, 1), 2)
  
  
  plt1 <- expression({
    ## Plot1: True Coefficients Plot
    beta <- `dimnames<-`(
      obj$get_properties("beta"), list(c(1:nx), c(paste0("Y", 1:ny))))
    beta_stk <- `names<-`(
      melt(beta), c("Predictor", "Response", "BetaCoef"))
    plt1 <- ggplot(
      beta_stk,
      aes(Predictor, BetaCoef,
          color = Response,
          group = Response)) +
      geom_hline(yintercept = 0,
                 color = "grey2",
                 linetype = 2) +
      geom_line() +
      geom_point() +
      scale_x_continuous(breaks = xticks) +
      labs(x = "Variable Number",
           y = expression(paste("Regression Coefficients (", beta, ")"))) +
      ggtitle("True Regression Coefficients") +
      theme(legend.position = if (ny == 1) "none" else "bottom")
  })
  
  plt2 <- expression({
    ## Plot 2: Relevant Comoponents Plot
    type <- obj$get_parameters("type")
    idx <- ny
    Sigma <- obj$get_properties("sigma")
    SigmaWZ <- obj$get_properties("sigma_zw")
    covs <- Sigma[-c(1:idx), 1:idx, drop = FALSE]
    covs.sc <- apply(covs, 2, function(x) {
      out <- abs(x)/max(abs(x))
      out[is.nan(out)] <- 0
      return(out)
    })
    covs.dt <- as.data.frame(
      cbind(
        1:ncomp,
        obj$get_properties("eigen_x")[1:ncomp],
        covs.sc[1:ncomp, ]
      ), 1
    )
    names(covs.dt) <- c("comp", "lambda", paste0("W", 1:ny))
    covs.stk <- melt(
      covs.dt, 1:2,
      variable.name = "response",
      value.name = "covariance"
    )
    pjtr <- position_dodge(0.5)
    plt2 <- ggplot(
      covs.stk,
      aes(comp, covariance, color = response, group = response)) +
      geom_bar(data = covs.dt,
               aes(x = comp, y = lambda),
               inherit.aes = FALSE,
               stat = "identity",
               fill = "lightgrey") +
      geom_line(position = pjtr) +
      geom_point(position = pjtr) +
      labs(x = "Components",
           y = "Covariance (absolute value)/ Eigenvalue") +
      theme(legend.position = if (ny == 1) "none" else "bottom") +
      ggtitle("Relevant Components Plot") +
      scale_x_continuous(breaks = xticks) +
      coord_cartesian(xlim = c(1, ncomp)) +
      scale_color_discrete(ifelse(obj$type == "multivariate",
                                  "Response Component", "Response"))
  })
  
  est.covs <- expression({
    dta <- obj$get_data("train")
    X <- scale(dta$x, center = TRUE, scale = FALSE)
    Y <- scale(dta$y, center = TRUE, scale = FALSE)
    
    svdres <- svd(X)
    eigval <- (svdres$d ^ 2)/(n - 1)  #/(svdres$d ^ 2)[1]
    eigval.sc <- eigval/eigval[1]
    
    Z <- X %*% svdres$v
    covs <- t(abs(cov(Y, Z)))
  })
  
  plt3 <- expression({
    ## Plot 3: Estimated Relevant Component Plot
    eval(est.covs)
    covs.sc <- apply(covs, 2, function(x) abs(x)/max(abs(x)))
    covs.dt <- as.data.frame(
      cbind(1:ncomp,
            eigval.sc[1:ncomp],
            covs.sc[1:ncomp, ]),
      1)
    names(covs.dt) <- c("comp", "lambda", paste0("Y", 1:ny))
    covs.stk <- melt(
      covs.dt, 1:2,
      variable.name = "response",
      value.name = "covariance"
    )
    pjtr <- position_dodge(0.5)
    plt3 <- ggplot(
      covs.stk,
      aes(comp, covariance, color = response, group = response)) +
      geom_bar(data = covs.dt,
               aes(x = comp, y = lambda),
               inherit.aes = FALSE,
               stat = "identity",
               fill = "lightgrey") +
      geom_line(position = pjtr) +
      geom_point(position = pjtr) +
      labs(x = "Components",
           y = "Covariance (absolute value)/ Eigenvalue") +
      theme(legend.position = if (ny == 1) "none" else "bottom") +
      ggtitle("Estimated relevant components plot") +
      scale_x_continuous(breaks = xticks) +
      coord_cartesian(xlim = c(1, ncomp)) +
      scale_color_discrete("Response")
  })
  
  plt <- list(
    TrueBeta = plt1,
    RelComp = plt2,
    EstRelComp = plt3
  )
  
  if (length(which) == 1) {
    out <- eval(plt[[which]])
  } else {
    plts <- lapply(which, function(i) eval(plt[[i]]))
    names(plts) <- names(plt)[which]
    
    if (length(which) == 2 & is.null(layout)) layout <- matrix(c(1, 2), 2)
    if (length(which) > length(layout)) layout <- matrix(1:length(which), length(which))
    
    plts$layout_matrix <- layout
    out <- do.call(grid.arrange, plts)
  }
  
  if (print.cov) {
    eval(est.covs)
    covs <- covs[1:min(ncomp, nx), ]
    dimnames(covs) <- list(
      paste0("Component", 1:min(ncomp, nx)),
      paste0("Response", 1:ny)
    )
    attr(out, "covariance") <- covs
    cat("Absolute values of estimated covariances\n")
    print(covs)
  }
  
  if (length(which) == 1) return(out) else return(invisible(out))
}
