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


    ## Loop for SGLD and variable selection (periodically)
    n_VSIter = 10

    SGLD_gibbs(data = data,  beta_curr = beta_curr, betaIdx_curr = betaIdx_curr)

    for (i in 1:nVSIter)
    {
        res_SGLD <- SGLD_gibbs(data = data,  beta = beta_start, betaIdx = betaIdx_start,
                               iter = iter, features_select = features_select, sig = sig)
        I0 <- MH$I
        beta_start <- MH$beta_start
        accept <- MH$accept
        ## 如果I没有被接收，SGLD迭代次数为SGLD_iter_noVS；如果I被接收，SGLD迭代次数为
        ## SGLD_iter不做此步 ，可以注释掉if(){}
        if(accept == 0) {iter = SGLD_iter_noVS}else{iter = SGLD_iter}
    }
    return(list(I = I, B = B, result_all = result_all, acceptance = accept_num/VS_iter))
}


SGLD_gibbs <- function(data, OUT, model_conf)
{
    nEpoch = model_conf$algArgs$sgld$nEpoch
    burninProp = model_conf$algArgs$sgld$burninProp

    max_batchsize = model_conf$algArgs$sgld$max_batchsize
    batchsize = min(nObs, max_batchsize)
    nBatch = round(nObs/batchsize)

    num_models_updated = ncol(data$lpd) - 1
    nObs = nrow(data$lpd)

    for (model_i in 1:num_models_updated)
    {
        ## 1. propose an update of variable selection indicators
        betaIdx_model_i <- OUT[["betaIdx"]][[model_i]][iIter - 1, ]
        nPar = length(betaIdx_model_i)

        betaIdx_prop = c(1, rbinom(nPar - 1, 1, prob = 0.5))

        ## 2. conditional on this variable selection indicators, update beta via SGLD
        beta_sgld = matrix(0, nEpoch * nBatch, nPar)
        beta_sgld[1, ] = beta_model_i

        for (iIter in 2:(nEpoch * nBatch))
        {
            ## SGLD settings
            if(is.null(stepsize)){
                stepsize1 <- a * (b + t) ^ (-gama)
            }else{
                stepsize1 <- stepsize
            }

            ## Re-split the data into small batches after finish one complete epoch.
            if(iIter in seq(2, nEpoch * nBatch, nObs))
            {
                batchIdx = 0
                dataIdxLst = split(sample(1:nObs,nObs),1:nBatch)
            }
            batchIdx = batchIdx + 1

            data_model_i = lapply(data[c("lpd","feat")], function(x) x[dataIdxLst[[batchIdx]], ,drop=FALSE])
            batchRatio = length(dataIdxLst[[batchIdx]]) / nObs

            grad_model_i = log_posterior_grad(data = data_model_i,
                                              beta_curr = beta_model_i,
                                              betaIdx_curr = betaIdx_model_i,
                                              priArgs = priArgs,
                                              varSelArgs = varSelArgs,
                                              features_used = features_used,
                                              model_update = model_conf$features_used,
                                              batchRatio = batchRatio)[[1]]

            ## SGLD
            beta_sgld[iIter, ] <- (beta_model_i + stepsize0 / 2 * grad_model_i +
                                   rmvnorm(1, rep(0,length(beta)), stepsize0* diag(length(beta_model_i))))
        }

        ## Polyak-ruppert averaging
        burnin = 1:(ceiling(burninProp * nIter))
        beta_prop = colMeans(beta_prop[-burnin, ])

        ## Metropolis-Hasting accept/reject
        logpost_prop = log_posterior(data = data,
                                     beta = beta_prop,
                                     betaIdx = betaIdx_prop,
                                     priArgs = priArgs,
                                     varSelArgs = varSelArgs,
                                     features_used = features_used)
        logpost_curr = log_posterior(data = data,
                                     beta = beta_curr,
                                     betaIdx = betaIdx_curr,
                                     priArgs = priArgs,
                                     varSelArgs = varSelArgs,
                                     features_used = features_used)



        logMHRatio <- (logPost.prop - logPost.curr +
                       logJump.currATpropRev - logJump.propATprop +
                       logJump.Idx.currATprop - logJump.Idx.propATcurr)


        if(is.na(logMHRatio))
        { ## bad proposal, i.e logJump.currATpropRev = -Inf, or logJump.propATprop = -Inf
            accept.prob.curr <- 0
        }
        else
        {
            accept.prob.curr <- exp(min(0, logMHRatio))
        }

        if(runif(1) < accept.prob.curr) #!is.na(accept.prob.curr)
        {
            ## keep the proposal
            betaIdx.curr <- betaIdx.prop
            beta.curr <- beta.prop

        }
        else
        { ## keep the current
            betaIdx.prop <- betaIdx.curr
        }

    }
    return(out)
}
