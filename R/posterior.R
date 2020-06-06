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

log_posterior_grad <- function(data, beta, betaIdx, priArgs, varSelArgs, features_used, model_update = 1:length(betaIdx), batchRatio = 1)
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


betaVec2Lst = function(beta, betaIdx)
{
    beta_list = list()
    a = 1
    for(model_i in 1:length(betaIdx))
    {
        b = a + length(betaIdx[[model_i]]) - 1
        beta_list[[model_i]] = beta[a:b]
        a = b + 1
    }

    return(beta_list)
}
