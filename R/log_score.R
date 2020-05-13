#' Calculate the log predictive score for a time series with pools of models
#'
#'
#' @title log predictive score with features
#' @param beta p-by-(n-1) matrix
#' @param features T-by-p feature matrix, usually standardized.
#' @param prob T-by-n predictive densities.
#' @param intercept TRUE or FALSE Should intercept be used in feature weights?
#' @param features_select a vector including the numbers of the features to be taken into
#'     consideration. IF is `NULL`, no features. At this time, the intercept must be
#' TRUE.
#' @return
#' @references Geweke & Amisano, (2011) Optimal prediction pools, Journal of Econometrics.
#' @note TODO: log_score_grad(beta, features, prob, intercepts)
#' @author Feng Li
log_score <- function(beta, features, features_select = ncol(features), prob, intercept = T)
{
    ## No features if features_select = NULL
    if(is.null(features_select)){
        features0 <- NULL
        intercept = TRUE
    }else{
        features0<-data.matrix(features[,features_select])
    }

    ## Intercept
    if(intercept) features0 = cbind(rep(1, nrow(prob)), features0)

    num_model <- dim(prob)[2]

    me <- features0 %*% beta
    me[me>709] <- 709 # avoid overflow
    exp_lin = exp(me)

    # Villani, M., Kohn, R., & Giordani, P. (2009). Regression density estimation using
    # smooth adaptive Gaussian mixtures. Journal of Econometrics, 153(2), 155-173.
    deno = matrix (rep((1+rowSums(exp_lin)), num_model-1), ncol = num_model-1)
    w <- exp_lin/ deno # T-by-(n-1)
    w_full = cbind(w, 1 - rowSums(w)) # T-by-n, assuming last is deterministic
    out = sum(log(rowSums(w_full * prob)))
    return(out)
}

#' Gradient of the log score with respect to given models
#'
#' The model gradient
#' @title log_score_grad
#' @param beta p-by-(n-1) matrix
#' @param features T-by-p feature matrix, usually standardized.
#' @param prob T-by-n predictive densities.
#' @param intercept TRUE or FALSE Should intercept be used in feature weights?
#' @param modelcaller which model components should be considered?
#' @return p-by-length(modelcaller) matrix for the gradient.
#' @author Feng Li
log_score_grad <- function(beta, features, features_select = ncol(features), prob,
                           intercept, modelcaller)
{
  if(is.null(features_select)){
      features0<-NULL
      intercept = TRUE
  }else{
    features0<-data.matrix(features[,features_select])
  }
  if(intercept) features0 = cbind(rep(1, nrow(prob)), features0)

  me <- features0 %*% beta # T-by-(n-1)
  me[me>709] <- 709 # avoid overflow

  ex_lin = exp(me) # T-by-(n-1)

  ex_sum = rowSums(ex_lin) + 1 # length-T vector
  ex_p_sum = rowSums(cbind(exp_lin, 1) * prob) # T-by-1

  ## The vectorized version
  grad0 = exp_lin[, modelcaller] * (prob[, modelcaller] * ex_sum -
                                    ex_p_sum) / (ex_sum * ex_p_sum) # T-by-ncol(modelcaller)

  out = apply(grad0, 2, function(x, y) colSums(x*y),y = features0)

  return(out)
}
