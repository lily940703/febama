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
    nIter = model_conf$algArgs$nIter
    sgldArgs = model_conf$algArgs$sgld

    num_models_updated = length(beta_curr)

    for (model_i in 1:num_models_updated)
    {
        for (iIter in 2:nIter)
        {
            beta_model_i <- OUT[["beta"]][[model_i]][iIter - 1, ]
            betaIdx_model_i <- OUT[["betaIdx"]][[model_i]][iIter - 1, ]

            ## SGLD settings
            if(is.null(stepsize)){
                stepsize1 <- a * (b + t) ^ (-gama)
            }else{
                stepsize1 <- stepsize
            }

            stepsize = 0.1

            ## SGLD
            beta_model_i_new <- beta_model_i + stepsize1 * (gradient_prior(beta_all, I, sig)[,i]/ prior1)
                + stepsize1 * (1/minibatchSize)
                * log_posterior_grad(beta = beta_all, features = features1,
                                     features_select = features_select,
                                     prob= prob1, intercept= intercept)[,i]

                + rmvnorm(1, rep(0,length(beta)), 2*stepsize1* diag(length(beta)))


            if(iIter in VS_interval)
            {
            }
        }

        res[[i]] <- res0
    }

    return(out)
}

MH_step <- function(x, beta0, data, logp = log_posterior,
                    proposal = proposal_I, sig = 10){
    beta_start <- beta0
    rownames(beta_start) <- c("0", which(x == 1))
    xp <- proposal(length(x))
    beta1 <- matrix(NA, nrow = (length(which(xp == 1)) + 1),
                    ncol = dim(beta_start)[2])
    rownames(beta1) <- c("0", which(xp == 1))
    ind <- c("0", which(x==xp & xp==1))
    beta1[ind, ] <- beta_start[ind, ]
    beta1[is.na(beta1)] <-0
    alpha <- min(1, exp(logp(data, beta = beta1, I = xp, prior = prior,
                             logLik = logscore, sig = sig) -
                        logp(data, beta = beta_start, I = x, prior = prior,
                             logLik = logscore, sig = sig)))
    if (runif(1) < alpha){
        accept <- 1
        x <- xp
        beta_start <- beta1
    }else{
        accept <- 0
    }
    return(list(I = x, beta_start = beta_start, accept = accept))
}

proposal_I <- function(n){
    repeat{
        I <- rbinom(n, 1, 0.5)
        if(sum(I != 0) != 0){
            break
        }
    }
    return(I)
}

beta_prepare <- function(res_SGLD_VS){
    beta_pre0 <-list()
    beta_pre <-list()
    for (i in 1 : length(res_SGLD_VS$B)) {
        beta_pre$beta <- res_SGLD_VS$B[[i]]
        beta_pre$features_select <- as.numeric(rownames(res_SGLD_VS$B[[i]]))[-1]
        beta_pre0[[i]] <- beta_pre
    }
    return(beta_pre0)
}
