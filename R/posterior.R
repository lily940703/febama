#' Calculate the log posterior for a time series.
#'
#' @param data: A list (the same as the input of function \code{logscore})
#' \describe{
#'   \item{lpd}{Log probability densities }
#'   \item{feat}{Features}
#' }
#' @param beta: A list of coefficient vectors of the features
#' @param betaIdx: A list of vectors including whether the features to be taken into
#'     consideration. IF is `NULL`, no features. At this time, the intercept must be
#' TRUE.
#' @param priArgs: parameter in \item{model_conf}.
#' @param varSelArgs: parameter in \item{model_conf}.
#' @param features_used: parameter in \item{model_conf}.
#' @return the log posterior
#' @export
log_posterior <- function(data, beta, betaIdx, priArgs, varSelArgs, features_used, model_update = 1:length(betaIdx))
{
    beta_list = beta

    ## Special case to allow for vector beta, used in optim MAP
    if(!is.list(beta))
    {
        beta_list = betaVec2Lst(beta, betaIdx)
    }

    ## log prior with conditional
    lpri = log_priors(beta = beta_list[model_update],
                      betaIdx = betaIdx[model_update],
                      varSelArgs = varSelArgs[model_update],
                      priArgs = priArgs[model_update], sum = TRUE)

    ## log score (log likelihood)
    lscore <- logscore(data = data,
                       beta = beta_list,
                       betaIdx = betaIdx,
                       features_used = features_used,
                       sum = TRUE)
    out <- lpri + lscore
    return(out)
}


#' When only optimize the coefficients of a particular model
log_posterior_comp <- function(data, beta_comp, beta, betaIdx, priArgs, varSelArgs, 
                               features_used, model_update)
{
    beta_list = beta
    
    ## Special case to allow for vector beta, used in optim MAP
    if(!is.list(beta))
    {
        beta_list = betaVec2Lst(beta, betaIdx)
    }
    
    ## log prior with conditional
    lpri = log_priors(beta = list(beta_comp),
                      betaIdx = betaIdx[model_update],
                      varSelArgs = varSelArgs[model_update],
                      priArgs = priArgs[model_update], sum = TRUE)
    
    ## log score (log likelihood)
    lscore <- logscore_comp(data = data,
                       beta_comp = beta_comp,
                       beta = beta_list,
                       betaIdx = betaIdx,
                       features_used = features_used,
                       model_update = model_update,
                       sum = TRUE)
    out <- lpri + lscore
    return(out)
}

#' Gradient of the log posterior with respect to given models
#'
#' @return A list of (number of models -1) vectors for the gradient wrt beta.
#' @export
log_posterior_grad <- function(data, beta, betaIdx, priArgs, varSelArgs, features_used,
                               model_update = 1:length(betaIdx), batchRatio = 1)
{
    ## log prior with conditional
    lpri_grad = log_priors_grad(beta = beta[model_update],
                                betaIdx = betaIdx[model_update],
                                varSelArgs = varSelArgs[model_update],
                                priArgs = priArgs[model_update])

    ## log score (log likelihood)
    lscore_grad <- logscore_grad(data = data,
                                 beta = beta,
                                 betaIdx = betaIdx,
                                 features_used = features_used,
                                 model_update = model_update)

    out <- mapply(function(x, y ) x / batchRatio + y,
                  lscore_grad, lpri_grad, SIMPLIFY = FALSE)
    return(out)
}

#' Transform beta (if a vector) into a list consistent with the format of the inputs of 
#' functions \code{logscore, log_priors, log_posterior...} 
#'
#' @return A list of (number of models -1) coefficient vectors of the features
#' @export
betaVec2Lst = function(beta, betaIdx)
{
    beta_list = list()
    a = 1
    for(iComp in 1:length(betaIdx))
    {
        b = a + length(betaIdx[[iComp]]) - 1
        beta_list[[iComp]] = beta[a:b]
        beta_list[[iComp]][betaIdx[[iComp]] == 0] = 0
        a = b + 1
    }

    return(beta_list)
}
