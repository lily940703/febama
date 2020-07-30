#' Calculate the log predictive score for a time series with pools of models.
#'
#' The intercept is always included in the model.
#' @title log predictive score with features
#' @param data: A list (the same as the output of function \code{lpd_features_multi})
#' \describe{
#'   \item{lpd}{Log probability densities }
#'   \item{feat}{Features}
#' }
#' @param beta: A list of coefficient vectors of the features
#' @param betaIdx: A list of vectors including whether the features to be taken into
#'     consideration. IF is `NULL`, no features. At this time, the intercept must be
#' TRUE.
#' @param features_used: parameter in \item{model_conf}.
#' @param sum: If TRUE, return the sume of log predictive densities.
#' @return the log predictive score
#' @references Geweke & Amisano, (2011) Optimal prediction pools, Journal of Econometrics.
#' @author Feng Li,  Li Li
logscore <- function(data, beta, betaIdx, features_used, sum = TRUE)
{
   if(!is.list(beta))
    {
        beta = betaVec2Lst(beta, betaIdx)
    }
    
    prob = exp(data$lpd)
    prob[prob == 0] <- 1e-16

    num_models_updated <- length(betaIdx)
    nObs = nrow(prob)

    exp_lin = matrix(NA, nObs, num_models_updated + 1)

    for(iComp in 1:num_models_updated)
    {
        betaCurr = beta[[iComp]]
        betaIdxCurr = betaIdx[[iComp]]
        features_used_curr = features_used[[iComp]]
        features0 = cbind(rep(1, nObs), data$feat[, features_used_curr, drop = FALSE])

        me <- features0[, betaIdxCurr == 1, drop = FALSE] %*% matrix(betaCurr[betaIdxCurr == 1])
        me[me>709] <- 709 # avoid overflow
        exp_lin[, iComp] = exp(me)
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
#' @return A list of (number of models -1) vectors for the gradient wrt beta.
#' @author Feng Li, Li Li
logscore_grad <- function(data, beta, betaIdx, features_used, model_update = 1:length(betaIdx))
{

    prob = exp(data$lpd)
    prob[prob == 0] <- 1e-16

    num_models_updated <- ncol(prob) - 1
    nObs = nrow(prob)

    exp_lin = matrix(NA, nObs, num_models_updated + 1)
    exp_lin[, num_models_updated + 1] = 1 # assume last model is 1

    for(iComp in 1:num_models_updated)
    {
        betaCurr = beta[[iComp]]
        betaIdxCurr = betaIdx[[iComp]]
        features_used_curr = features_used[[iComp]]
        features0 = cbind(rep(1, nObs), data$feat[, features_used_curr, drop = FALSE])

        me <- features0[, betaIdxCurr == 1, drop = FALSE] %*% matrix(betaCurr[betaIdxCurr == 1])
        me[me>709] <- 709 # avoid overflow
        exp_lin[, iComp] = exp(me)
    }

    exp_sum = rowSums(exp_lin) # length-T vector
    exp_p_sum = rowSums(exp_lin * prob) # T-by-1

    ## The gradient wrt me=x'beta
    grad0 = exp_lin[, model_update, drop = FALSE] * (prob[, model_update, drop = FALSE] * exp_sum -
                                      exp_p_sum) / (exp_sum * exp_p_sum) # T-by-length(model_caller)

    ## The gradient wrt beta
    out = list()
    iCompdx = 0
    for(iComp in model_update)
    {
        iCompdx = iCompdx + 1
        features_used_curr = features_used[[iComp]]
        betaIdxCurr = betaIdx[[iComp]]
        features0 = cbind(rep(1, nObs), data$feat[, features_used_curr, drop = FALSE])
        out[[iCompdx]] = colSums(grad0[, iCompdx] * features0[, betaIdxCurr == 1, drop = FALSE])
    }
    return(out)
}
