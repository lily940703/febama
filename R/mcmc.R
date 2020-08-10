#' Optimization to obtain beta and variable selection result by a MCMC method 
#'
#' @param data: A list (the same as the output of function \code{lpd_features_multi})
#' \describe{
#'   \item{lpd}{Log probability densities }
#'   \item{feat}{Features}
#' }
#' @param model_conf
#' @return A list 
#' \describe{
#'   \item{beta}{A list of (number of models -1) matrices of beta in every iteration }
#'   \item{betaIdx}{A list of (number of models -1) matrices of betaIdx in every iteration}
#'   \item{accept_prob}{A list of (number of models -1) matrices of accept probabilities in every iteration}
#' }
#' export
febama_mcmc <- function(data, model_conf)
{
    ## Extract arguments
    algArgs = model_conf$algArgs
    varSelArgs = model_conf$varSelArgs
    priArgs = model_conf$priArgs

    nIter = algArgs$nIter
    num_models_updated = ncol(data$lpd) - 1

    ## Reserve space for OUT beta and initialize variable selection indicators
    OUT = list()
    for(iComp in 1:num_models_updated)
    {
        ## browser()
        ## Determine number of features used.
        nFeat = length(model_conf$features_used[[iComp]]) # default number of features = features_used

        ## nFeat = 0 if only intercept is used.
        if(length(model_conf$features_used[[iComp]]) == 0)
        {
            nFeat = 0
        }

        OUT[["beta"]][[iComp]] = matrix(NA, nIter, nFeat + 1)
        OUT[["betaIdx"]][[iComp]] = matrix(NA, nIter, nFeat + 1)

        ## Initialize variable selection indicators
        if(varSelArgs[[iComp]]$init == "all-in")
        {
            OUT[["betaIdx"]][[iComp]][1, ] = 1
        }else if (varSelArgs[[iComp]]$init == "all-out")
        {
            OUT[["betaIdx"]][[iComp]][1, ] = 0
            OUT[["betaIdx"]][[iComp]][1, 1] = 1
        }else if(varSelArgs[[iComp]]$init == "random")
        {
            OUT[["betaIdx"]][[iComp]][1, ] = c(1, rbinom(nFeat, 1, 0.5)) # intercept is always in.
        }else
        {
            stop("No such init for betaIdx!")
        }

        ## Reserve space for acceptance probabilities with MH corrections
        OUT[["accept_prob"]][[iComp]] = matrix(1, nIter, 1)
    }

    betaIdx_curr = lapply(OUT[["betaIdx"]], function(x) x[1,])
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
                           control = list(fnscale = -1, maxit = 100))
        if( beta_optim$convergence != 0){
            stop("The optimization of initial values is not convergent")
        }
        beta_curr = betaVec2Lst(beta_optim$par, betaIdx_curr)
    }

    ## Assign initial values (conditional of variable selections)
    OUT[["beta"]] = mapply(function(x, y){
        x[1, ] = y
        return(x)
    }, x = OUT[["beta"]], y = beta_curr, SIMPLIFY = FALSE)
                       
    if(nIter > 1){
    for (iIter in 2:nIter)
    { # Loop start with the second iteration. 
      # The first iteration is considered as initial values.
      if(algArgs$algName == "sgld"){
        beta_betaIdx <-  SGLD_gibbs(data = data,  beta_curr = beta_curr,
                                    betaIdx_curr = betaIdx_curr, model_conf = model_conf)
      }else if(algArgs$algName == "MAP"){
        beta_betaIdx <-  MAP_gibbs(data = data,  beta_curr = beta_curr,
                                    betaIdx_curr = betaIdx_curr, model_conf = model_conf)
      }

        ## Extract parameters for loop use
        beta_curr = beta_betaIdx[["beta"]]
        betaIdx_curr = beta_betaIdx[["betaIdx"]]
        accept_prob_curr = beta_betaIdx[["accept_prob"]]

        ## Assign the final output
        OUT[["beta"]] = mapply(function(x, y){
            x[iIter, ] = y
            return(x)
        }, x = OUT[["beta"]], y = beta_betaIdx[["beta"]], SIMPLIFY = FALSE)

        OUT[["betaIdx"]] = mapply(function(x, y){
            x[iIter, ] = y
            return(x)
        }, x = OUT[["betaIdx"]], y = beta_betaIdx[["betaIdx"]], SIMPLIFY = FALSE)

        OUT[["accept_prob"]] = mapply(function(x, y){
            x[iIter, ] = y
            return(x)
        }, x = OUT[["accept_prob"]], y = beta_betaIdx[["accept_prob"]], SIMPLIFY = FALSE)
    
    }
    }
    return(OUT)
}


#' @title Variable selection using MAP with Gibbs
#' @param data: A list (the same as the output of function \code{lpd_features_multi})
#' \describe{
#'   \item{lpd}{Log probability densities }
#'   \item{feat}{Features}
#' }
#' @param beta_curr: A list of coefficient vectors of the features (The initial value of a variable selection)
#' @param betaIdx_curr: A list of vectors including whether the features to be taken into
#'     consideration. (The initial value of a variable selection)
#' @param model_conf
#' @return A list (result of variable selection)
MAP_gibbs <- function(data, beta_curr, betaIdx_curr, model_conf)
{

    priArgs = model_conf$priArgs

    features_used = model_conf$features_used

    num_models_updated = ncol(data$lpd) - 1

    beta_prop = beta_curr
    betaIdx_prop = betaIdx_curr
    accept_prob = list()
    for (iComp in 1:num_models_updated)
    {
        nPar_full = length(betaIdx_curr[[iComp]])

        betaIdx_prop[[iComp]] = c(1, rbinom(nPar_full - 1, 1, prob = 0.5))

        beta_optim = optim(unlist(beta_prop[[iComp]]), fn = log_posterior_comp,
                           data = data,
                           beta = beta_prop,
                           betaIdx = betaIdx_prop,
                           priArgs = priArgs,
                           varSelArgs = varSelArgs,
                           features_used = model_conf$features_used,
                           model_update = iComp,
                           method = "BFGS",
                           control = list(fnscale = -1, maxit = 100))
        if( beta_optim$convergence != 0){
          stop("The optimization of MAP is not convergent")
        }
        
        beta_prop[[iComp]] = beta_optim$par
        beta_prop[[iComp]][betaIdx_prop[[iComp]] == 0] = 0

        ## Metropolis-Hasting accept/reject for variable selection
        if(length(model_conf$varSelArgs[[iComp]]$cand) > 0)
        {
            log_post_prop = log_posterior(data = data,
                                          beta = beta_prop,
                                          betaIdx = betaIdx_prop,
                                          priArgs = priArgs,
                                          varSelArgs = varSelArgs,
                                          features_used = features_used,
                                          model_update = iComp)
            log_post_curr = log_posterior(data = data,
                                          beta = beta_curr,
                                          betaIdx = betaIdx_curr,
                                          priArgs = priArgs,
                                          varSelArgs = varSelArgs,
                                          features_used = features_used,
                                          model_update = iComp)


            ## The jump density for the variable selection indicators. TODO: Add adaptive scheme
            logJump.Idx.currATprop <- 1
            logJump.Idx.propATcurr <- 1

            logMHRatio <- (log_post_prop - log_post_curr +
                           logJump.Idx.currATprop - logJump.Idx.propATcurr)


            if(is.na(logMHRatio))
            { ## bad proposal, i.e logJump.currATpropRev = -Inf, or logJump.propATprop = -Inf
                accept_prob_curr <- 0
            } else{
                accept_prob_curr <- exp(min(0, logMHRatio))
            }

            if(runif(1) < accept_prob_curr) #!is.na(accept.prob.curr)
            {  ## keep the proposal
                beta_curr[[iComp]] <- beta_prop[[iComp]]
                betaIdx_curr[[iComp]] <- betaIdx_prop[[iComp]]
            }
            else
            { ## keep the current
                beta_prop[[iComp]] = beta_curr[[iComp]]
                betaIdx_prop[[iComp]] = betaIdx_curr[[iComp]]
            }
        }
        else
        {
            accept_prob_curr = 1 # SGLD always accepts
            beta_curr[[iComp]] <- beta_prop[[iComp]] ## Accepted
            ## betaIdx_curr unchanged.
        }
        accept_prob[[iComp]] <- accept_prob_curr
        
    }

    out = list(beta = beta_curr, betaIdx = betaIdx_curr, 
               accept_prob = accept_prob)
    return(out)
}



SGLD_gibbs <- function(data, beta_curr, betaIdx_curr, model_conf)
{
  nObs = nrow(data$lpd)
  #nEpoch = model_conf$algArgs$sgld$nEpoch
  stepsize = model_conf$algArgs$sgld$stepsize
  burninProp = model_conf$algArgs$sgld$burninProp
  
  priArgs = model_conf$priArgs
  varSelArgs = model_conf$varSelArgs
  
  features_used = model_conf$features_used
  
  max_batchSize = model_conf$algArgs$sgld$max_batchSize
  batchSize = min(nObs, max_batchSize)
  Batch = round(nObs/batchSize)
  nEpoch = round(200 / nBatch)
  
  num_models_updated = ncol(data$lpd) - 1
  
  beta_prop = beta_curr
  betaIdx_prop = betaIdx_curr
  
  beta_sgld = list()
  for (iComp in 1:num_models_updated)
  {
    nPar_full = length(betaIdx_curr[[iComp]])
    
    ## 1. propose an update of variable selection indicators. Random proposal when
    ## variable selection is enabled: NOT NULL.
    if(length(model_conf$varSelArgs[[iComp]]$cand) > 0)
    {
      betaIdx_prop[[iComp]] = c(1, rbinom(nPar_full - 1, 1, prob = 0.5))
    }
    
    ## 2. conditional on this variable selection indicators, update beta via SGLD
    beta_iComp_sgld = matrix(0, nEpoch * nBatch, nPar_full)
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
        dataIdxLst = suppressWarnings(split(sample(1:nObs,nObs),1:nBatch))
      }
      iBatch = iBatch + 1
      
      data_curr = lapply(data[c("lpd","feat")], function(x) x[dataIdxLst[[iBatch]], ,drop=FALSE])
      batchRatio = length(dataIdxLst[[iBatch]]) / nObs # n/N
      
      grad_iComp = log_posterior_grad(data = data_curr,
                                      beta = beta_prop,
                                      betaIdx = betaIdx_prop,
                                      priArgs = priArgs,
                                      varSelArgs = varSelArgs,
                                      features_used = features_used,
                                      model_update = iComp,
                                      batchRatio = batchRatio)[[1]]
      
      ## SGLD
      nPar1 = sum(betaIdx_prop[[iComp]]) # length of non-zero parameters
      beta_new <- (as.vector(beta_prop[[iComp]][betaIdx_prop[[iComp]] == 1] + stepsize / 2 * grad_iComp) +
                     as.vector(mvtnorm::rmvnorm(1, rep(0, nPar1), stepsize* diag(nPar1))))
      
      beta_iComp_sgld[iIter, betaIdx_prop[[iComp]] == 1] = beta_new
      
      beta_prop[[iComp]][betaIdx_prop[[iComp]] == 1] = beta_new
      betaIdx_prop[[iComp]][betaIdx_prop[[iComp]] == 0] = 0
    }
    
    ## Polyak-Ruppert averaging improve the efficiency of SGLD
    burnin = 1:(ceiling(burninProp * nEpoch * nBatch))
    beta_prop[[iComp]] = colMeans(beta_iComp_sgld[-burnin,, drop = FALSE])
    beta_sgld[[iComp]] = beta_iComp_sgld
    
    ## Metropolis-Hasting accept/reject for variable selection
    if(length(model_conf$varSelArgs[[iComp]]$cand) > 0)
    {
      log_post_prop = log_posterior(data = data,
                                    beta = beta_prop,
                                    betaIdx = betaIdx_prop,
                                    priArgs = priArgs,
                                    varSelArgs = varSelArgs,
                                    features_used = features_used,
                                    model_update = iComp)
      log_post_curr = log_posterior(data = data,
                                    beta = beta_curr,
                                    betaIdx = betaIdx_curr,
                                    priArgs = priArgs,
                                    varSelArgs = varSelArgs,
                                    features_used = features_used,
                                    model_update = iComp)
      
      
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
        beta_curr[[iComp]] <- beta_prop[[iComp]]
        betaIdx_curr[[iComp]] <- betaIdx_prop[[iComp]]
      }
      else
      { ## keep the current
        beta_prop[[iComp]] = beta_curr[[iComp]]
        betaIdx_prop[[iComp]] = betaIdx_curr[[iComp]]
      }
    }
    else
    {
      accept_prob_curr = 1 # SGLD always accepts
      beta_curr[[iComp]] <- beta_prop[[iComp]] ## Accepted
      ## betaIdx_curr unchanged.
    }
    
  }
  
  out = list(beta = beta_curr, betaIdx = betaIdx_curr, 
             accept_prob = accept_prob_curr, beta_sgld = beta_sgld)
  return(out)
}
