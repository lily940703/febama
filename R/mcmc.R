febama_mcmc <- function(data, model_conf)
{
    ## Extract arguments
    algArgs = model_conf$algArgs
    varSelArgs = model_conf$varSelArgs
    priArgs = model_conf$priArgs

    nIter = algArgs$nIter

    ## Reserve space for beta
    OUT = list()
    OUT[["beta"]] = lapply(model_conf$features,
                           function(x) matrix(nrow=nIter,ncol=length(x)+1))

    ## Initialize variable selection indicators
    OUT[["betaIdx"]] = mapply(function(x, y){
        if(y$init == "all-in")
        {
            x[1, ] = 1
        }
        else if(y$init == "random")
        {
            x[1, ] = c(1, rbinom(length(x) - 1, 1, 0.5)) # intercept is always in.
        }
        return(x)
    }, x = OUT_beta_full, y = varSelArgs, SIMPLIFY = FALSE)

    betaIdx_curr = lapply(OUT_betaIdx_full, function(x) x[1,])
    beta_curr = lapply(betaIdx_curr, function(x) rnorm(length(x)))

    ## Numeric optimization to obtain MAP (Maximum a Posteriori)
    if(algArgs$initOptim == TRUE)
    {
        beta_optim = optim(unlist(beta_curr), fn = log_posterior,
                           data = data,
                           betaIdx = betaIdx_curr,
                           priArgs = priArgs,
                           varSelArgs = varSelArgs,
                           features_used = model_conf$features_used,
                           method = "BFGS",
                           control = list(fnscale = -1, maxit = 10))
        beta_curr = betaVec2Lst(beta_optim$par, betaIdx_curr)
    }

    ## Assign initial values (conditional of variable selections)
    OUT[["beta"]] = mapply(function(x, y){
        x[1, ] = y
        return(x)
    }, x = OUT_beta_full, y = beta_curr, SIMPLIFY = FALSE)



    for (i in 1:nIter)
    {
        beta_betaIdx <-  SGLD_gibbs(data = data,  beta_curr = beta_curr,
                                    betaIdx_curr = betaIdx_curr, model_conf = model_conf)

        beta_curr = beta_betaIdx[["beta"]]
        betaIdx_curr = beta_betaIdx[["betaIdx"]]


        OUT[["beta"]] = mapply(function(x, y){
            x[i, ] = y
            return(x)
        }, x = OUT_beta_full, y = beta_betaIdx[["beta"]], SIMPLIFY = FALSE)

        OUT[["betaIdx"]] = mapply(function(x, y){
            x[i, ] = y
            return(x)
        }, x = OUT_beta_full, y = beta_betaIdx[["beta"]], SIMPLIFY = FALSE)

    }

    return(OUT)
}


SGLD_gibbs <- function(data, beta_curr, betaIdx_curr, model_conf)
{
    nEpoch = model_conf$algArgs$sgld$nEpoch
    burninProp = model_conf$algArgs$sgld$burninProp

    max_batchSize = model_conf$algArgs$sgld$max_batchSize
    batchSize = min(nObs, max_batchSize)
    nBatch = round(nObs/batchSize)

    num_models_updated = ncol(data$lpd) - 1
    nObs = nrow(data$lpd)

    for (model_i in 1:num_models_updated)
    {
        ## 1. propose an update of variable selection indicators
        betaIdx_model_i <- betaIdx_curr[[model_i]]
        beta_model_i = beta_curr[[model_i]]
        nPar = length(betaIdx_model_i)

        betaIdx_prop = betaIdx_curr
        betaIdx_prop[[model_i]] = c(1, rbinom(nPar - 1, 1, prob = 0.5))

        ## 2. conditional on this variable selection indicators, update beta via SGLD
        beta_model_i_sgld = matrix(0, nEpoch * nBatch, nPar)

        for (iIter in 1:(nEpoch * nBatch))
        {
            ## SGLD settings
            if(is.na(stepsize)){ # decay stepsize
                a = model_conf$algArgs$sgld$a
                b = model_conf$algArgs$sgld$b
                gama = model_conf$algArgs$sgld$gama
                stepsize <- a * (b + iIter) ^ (-gama)
            }

            ## Re-split the data into small batches after finish one complete epoch.
            if(iIter %in% seq(1, nEpoch * nBatch, nBatch))
            {
                iBatch = 0
                dataIdxLst = split(sample(1:nObs,nObs),1:nBatch)
            }
            iBatch = iBatch + 1

            data_curr = lapply(data[c("lpd","feat")], function(x) x[dataIdxLst[[iBatch]], ,drop=FALSE])
            batchRatio = length(dataIdxLst[[iBatch]]) / nObs # n/N

            grad_model_i = log_posterior_grad(data = data_curr,
                                              beta = beta_curr,
                                              betaIdx = betaIdx_prop,
                                              priArgs = priArgs,
                                              varSelArgs = varSelArgs,
                                              features_used = features_used,
                                              model_update = model_i,
                                              batchRatio = batchRatio)[[1]]

            ## SGLD
            nPar_new = sum(betaIdx_prop[[model_i]])

            beta_new <- (as.vector(beta_model_i[betaIdx_prop[[model_i]] == 1] + stepsize / 2 * grad_model_i) +
                         as.vector(rmvnorm(1, rep(0, nPar_new), stepsize* diag(nPar_new))))

            beta_model_i_sgld[iIter, betaIdx_prop[[model_i]] == 1] = beta_new

            beta_curr[[model_i]][betaIdx_prop[[model_i]] == 1] = beta_new
            beta_curr[[model_i]][betaIdx_prop[[model_i]] == 0] = 0

        }

        ## Polyak-ruppert averaging
        burnin = 1:(ceiling(burninProp * nEpoch * nBatch))

        beta_prop = beta_curr
        beta_prop[[model_i]] = colMeans(beta_model_i_sgld[-burnin, ])

        ## Metropolis-Hasting accept/reject
        log_post_prop = log_posterior(data = data,
                                      beta = beta_prop,
                                      betaIdx = betaIdx_prop,
                                      priArgs = priArgs,
                                      varSelArgs = varSelArgs,
                                      features_used = features_used,
                                      model_update = model_i)
        log_post_curr = log_posterior(data = data,
                                      beta = beta_curr,
                                      betaIdx = betaIdx_curr,
                                      priArgs = priArgs,
                                      varSelArgs = varSelArgs,
                                      features_used = features_used,
                                      model_update = model_i)


        ## The jump density for the variable selection indicators. TODO: Add adaptive scheme
        logJump.Idx.currATprop <- 1
        logJump.Idx.propATcurr <- 1

        logMHRatio <- (log_post_prop - log_post_curr +
                       logJump.Idx.currATprop - logJump.Idx.propATcurr)


        if(is.na(logMHRatio))
        { ## bad proposal, i.e logJump.currATpropRev = -Inf, or logJump.propATprop = -Inf
            accept_prob_curr <- 0
        }
        else
        {
            accept_prob_curr <- exp(min(0, logMHRatio))
        }

        if(runif(1) < accept_prob_curr) #!is.na(accept.prob.curr)
        {  ## keep the proposal
            beta_curr <- beta_prop
            betaIdx_curr <- betaIdx_prop
        }
        else
        { ## keep the current
        }

    }

    out = list(beta = beta_curr, betaIdx = betaIdx_curr)

    return(out)
}
