SGLD_VS <- function(data, logLik, logLik_grad, prior, model_conf)
{
    num_features  <- dim(data$feat)[2]
    num_models <- dim(data$lpd)[2]
    I <- matrix(nrow = num_features , ncol = VS_iter)
    B <- list()
    result_all <- list()
    I[,1] <- rbinom(num_features,1,0.5)
    accept_num <- 1
    iter <- SGLD_iter
    for (i in 1:VS_iter) {
        if(i == 1){
            features_select <- which(I[,i]==1)
            beta_start <- matrix(runif((length(features_select) + 1) * (num_models-1), -10, 10),
                                 ncol = num_models-1)
            rownames(beta_start) <- c("0", features_select)
        }else{
            accept_num <- accept_num +accept
            I[,i] <-I0
            features_select <- which(I[,i]==1)
            beta_start <- beta_start
        }
        res_SGLD <- SGLD_gibbs(data = data, logLik = logscore, logLik_grad = logscore_grad,
                             prior = prior, beta = beta_start, betaIdx = betaIdx,
                             minibatchSize = minibatchSize, stepsize = stepsize,
                             iter = iter, features_select = features_select, sig = sig)
        B[[i]] <- res_SGLD$beta_out
        result_all[[i]] <- res_SGLD
        I0 <- I[,i]
        MH <- MH_step (x = I0, beta0 = B[[i]], data =data,
                       logp = log_posterior, proposal = proposal_I)
        I0 <- MH$I
        beta_start <- MH$beta_start
        accept <- MH$accept
                                        # 如果I没有被接收，SGLD迭代次数为SGLD_iter_noVS；
                                        # 如果I被接收，SGLD迭代次数为SGLD_iter
                                        # 不做此步 ，可以注释掉if(){}
        if(accept == 0) {iter = SGLD_iter_noVS}else{iter = SGLD_iter}
    }
    return(list(I = I, B = B, result_all = result_all, acceptance = accept_num/VS_iter))
}


SGLD_gibbs <- function(data, logLik, logLik_grad, prior, beta, beta_idx, model_conf)
{
    beta_all <- start

    prob <- exp(data$lpd)
    prob[prob == 0] <- 1e-323
    num_models <- dim(prob)[2]
    features <- data$feat

    betaIdx = I

    if (length(prob[,1]) <= 64){
        minibatchSize = 1
    }else {
        minibatchSize = 0.1
    }

    ## Optimize to have good initial values
    ## w_max <- try(optim(
    ##     par = ini,
    ##     fn = log_posterior,
    ##     features = features_y,
    ##     prob = prob,
    ##     ## intercept = intercept,
    ##     intercept = TRUE,
    ##     method = "L-BFGS-B",
    ##     ## method = "BFGS",
    ##     control = list(fnscale = -1)
    ## ))



    res <- list()
    for (i in 1:(num_models-1)) {

        browser()
        prior0 <- log_priors(beta = beta_all, betaIdx = betaIdx, varSelArgs = model_conf$varSelArgs, priArgs = model_conf$priArgs, sum = TRUE)

        logLik0 <- logscore (beta = beta_all, features = features,
                           features_select = features_select,
                           prob = prob, intercept = intercept,
                           sum = TRUE)

        logpost0 <- log_posterior (data, beta_all, I, prior = prior,
                                   logLik = logscore, sig = sig)
        res0 <- list(beta = beta_all[,i], logscore = logLik0, logpost = logpost0,
                     stepsize = NA, prior = prior0)

        if(iter != 0){
            for (t in 1:iter ) {

                beta <- beta_all[,i]
                mini <- sample(1:length(prob[,1]),
                               ceiling(minibatchSize*length(prob[,1])))
                prob1 <- prob[mini,]
                features1 <- features[mini,]

                if(is.null(stepsize)){
                    stepsize1 <- a * (b + t) ^ (-gama)
                }else{
                    stepsize1 <- stepsize
                }

                prior1 <- prior(beta_all, I, sig = sig)

                beta <- (beta + stepsize1 * (gradient_prior(beta_all, I, sig)[,i]/ prior1)
                    + stepsize1 * (1/minibatchSize)
                    * logLik_grad(beta = beta_all, features = features1,
                                  features_select = features_select,
                                  prob= prob1, intercept= intercept)[,i]
                    + rmvnorm(1, rep(0,length(beta)), 2*stepsize1* diag(length(beta)))
                )
                beta_all[,i] <- beta

                prior1 <- prior(beta_all, I, sig = sig)
                logLik1 <- logLik (beta = beta_all, features = features, features_select = features_select,
                                   prob = prob, intercept = intercept)
                logpost1 <- log_posterior (data, beta_all, I, prior = prior,
                                           logLik = logscore, sig = sig)
                res0$beta <- cbind(res0$beta, beta)
                res0$logscore <- c(res0$logscore, logLik1)
                res0$logpost <- c(res0$logpost, logpost1)
                res0$stepsize <- c(res0$stepsize, stepsize1)
                res0$prior <- c(res0$prior, prior1)

                                        # 完成iter次SGLD后，对于每个β，取后10%样本的均值
                                        # 若注释掉if{}， beta_out为每个β的最后一个样本组成的矩阵

                if(t == iter){
                    if(iter < 10) {ind = 1:(iter+1)} else {
                                                       ind =( iter + 1 - floor(0.1 * iter)): (iter+1)}
                    if(dim(res0$beta)[1] == 1){
                        beta_mean <- mean(res0$beta[, ind])

                    }else{
                        beta_mean <- rowMeans(res0$beta[, ind])
                    }
                    beta_all[,i] <- beta_mean
                }
            }
        }
        res[[i]] <- res0
    }
    results <- list(beta=list(), logscore = matrix(ncol = iter+1, nrow = num_models-1),
                    logpost = matrix(ncol = iter+1, nrow = num_models-1),
                    stepsize = matrix(ncol = iter+1, nrow = num_models-1),
                    prior = matrix(ncol = iter+1, nrow = num_models-1),
                    beta_out = beta_all)
    for (j in 1:(num_models-1)) {
        results$beta[[j]] <- res[[j]]$beta
        results$logscore[j,] <-res[[j]]$logscore
        results$logpost[j,] <-res[[j]]$logpost
        results$stepsize[j,] <-res[[j]]$stepsize
        results$prior[j,] <-res[[j]]$prior
    }
    return(results)
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

#' Simple algorithm to optimize model parameters.
#'
#' This can be used as initial values in MCMC or just for simple optimizations
#' @title
#' @param lpd_features
#' @param features_y
#' @return
#' @author Feng Li
optim_beta <- function(lpd_features, features_y = NULL) {
    y_lpd <- lpd_features$lpd

    prob <- exp(y_lpd)
    prob[prob == 0] <- 1e-323

    num_models <- length(lpd_features$lpd[1,])

    ini <-  t(data.matrix(rep(0, num_models-1)))

    w_max <- try(optim(
        par = ini,
        fn = logscore,
        ## gr = gradient,
        features = features_y,
        prob = prob,
        ## intercept = intercept,
        intercept = TRUE,
        method = "L-BFGS-B",
        ## method = "BFGS",
        control = list(fnscale = -1)
    )
    )

    if (w_max$convergence != 0) {
        w_max <- try(optim(
            par = ini,
            fn = logscore,
            ## gr = gradient,
            features = features_y,
            prob = prob,
            ## intercept = intercept,
            intercept = TRUE,
            ## method = "L-BFGS-B",
            method = "BFGS",
            control = list(fnscale = -1)
        )
        )
    }
    beta_optim <- w_max$par
    return(list(beta_optim = beta_optim, logscore = w_max$value))
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
