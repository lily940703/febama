#' Calculate the log predictive score for a time series with pools of models.
#'
#' The intercept is always included in the model.
#' @title log predictive score with features
#' @param beta list
#' @param features T-by-p feature matrix, usually standardized.
#' @param betaIdx a vector including the numbers of the features to be taken into
#'     consideration. IF is `NULL`, no features. At this time, the intercept must be
#' TRUE.
#' @param prob T-by-n predictive densities.
#' @param sum If TRUE, return the sume of log predictive densities.
#' @return List
#' @references Geweke & Amisano, (2011) Optimal prediction pools, Journal of Econometrics.
#' @author Feng Li
logscore <- function(data, beta, betaIdx, features_used, sum = TRUE)
{
    prob = exp(data$lpd)
    prob[prob == 0] <- 1e-16

    num_models_updated <- length(beta)
    nObs = nrow(prob)

    exp_lin = matrix(NA, nObs, num_models_updated + 1)

    for(model_i in 1:num_models_updated)
    {
        betaCurr = beta[[model_i]]
        betaIdxCurr = betaIdx[[model_i]]
        features_used_curr = features_used[[model_i]]
        features0 = cbind(rep(1, nObs), data$feat[, features_used_curr, drop = FALSE])

        me <- features0[, betaIdxCurr == 1, drop = FALSE] %*% matrix(betaCurr[betaIdxCurr == 1])
        me[me>709] <- 709 # avoid overflow
        exp_lin[, model_i] = exp(me)
    }

    ## Villani, M., Kohn, R., & Giordani, P. (2009). Regression density estimation
    ## using smooth adaptive Gaussian mixtures. Journal of Econometrics, 153(2),
    ## 155-173.
    exp_lin[, num_models_updated + 1] = 1 # assume last model is 1
    weights <- exp_lin/ rowSums(exp_lin) # T-by-n, assuming last is deterministic
    out = log(rowSums(weights * prob))

    if(sum == TRUE)
    {
        out = sum(out)
    }

    return(out)
}

#' Gradient of the log score with respect to given models
#'
#' The model gradient
#' @title logscore_grad
#' @param beta p-by-(n-1) matrix
#' @param features T-by-p feature matrix, usually standardized.
#' @param prob T-by-n predictive densities.
#' @param intercept TRUE or FALSE Should intercept be used in feature weights?
#' @param modelcaller which model components should be considered?
#' @return p-by-length(modelcaller) matrix for the gradient.
#' @author Feng Li
logscore_grad <- function(data, beta, betaIdx, features_used, modelcaller)
{
  if(is.null(betaIdx)){
      features0<-NULL
      intercept = TRUE
  }else{
    features0<-data.matrix(features[,betaIdx])
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
