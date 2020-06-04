#' This is the prior settings for a variable selection scheme
#'
#' The detailed documentation is in the main setting file for each parameter.
#' @title Model prior settings.
#' @param beta list for models to be updated
#' @param betaIdx same as beta
#' @param varSelArgs list, same as beta length
#' @param priArgs list, same as beta length
#' @param sum logical, if TRUE, sum all conditional priors
#' @return "list" synced
#' @references "list"
#' @export
log_priors <- function(beta, betaIdx, varSelArgs, priArgs, sum = TRUE)
{

    num_models_updated = length(beta)
    out = rapply(priArgs, function (x) 0, how="replace") # retain same out structure as priArgs.

    for(model_i in 1:num_models_updated)
    {
        priArgsCurr <- priArgs[[model_i]]
        betaCurr <- beta[[model_i]]
        betaIdxCurr <- betaIdx[[model_i]] # p-by-1

###----------------------------------------------------------------------------
### Prior for variable selection indicators
###----------------------------------------------------------------------------
        ## intercept as special case. The intercept should always be included.
        varSelCand <- varSelArgs[[model_i]][["cand"]]

        ## Variable section candidates checking. If (1) only intercept is in or (2)
        ## updating candidates is NULL, do not allow for variable selection. We do not
        ## need a prior;
        if(length(betaIdxCurr) == 1 | length(varSelCand) == 0) # intercept is always included.
        {
            candIdx <- NULL
            logDens <- NULL
        }
        else
        {
            if(class(varSelCand) == "character" &&
               tolower(varSelCand) == "2:end")
            {
                candIdx <- 2:length(beta[[model_i]])
            }
            else
            {
                candIdx <- varSelCand
            }

            ## Priors for variable selection indicators
            varSelCandTF <- betaIdxCurr[candIdx]

            if(tolower(priArgsCurr[["betaIdx"]][["type"]]) == "bern") # Bernoulli prior
            {
                ## Note that independent Bernoulli prior can be very informative. See
                ## Scott and Berger 2010 in Annals of Statistics Sec 3.2.
                prob <- priArgs[["prob"]]
                probMat <- matrix(prob, length(candIdx), 1) # TRUE or FALSE of variable selection candidates
                logDens <- sum(dbinom(x = as.numeric(varSelCandTF), size = 1,
                                      prob = probMat, log = TRUE))
            }
            else if(tolower(priArgsCurr[["betaIdx"]][["type"]]) == "beta") # beta prior
            {  ## Scott and Berger 2010 in Annals of Statistics Sec 3.3
                alpha0 = priArgsCurr[["betaIdx"]][["alpha0"]]
                beta0  = priArgsCurr[["betaIdx"]][["beta0"]]

                varSelCand.N <- length(candIdx)
                varSelCandT.N <- sum(varSelCandTF)
                logDens <- (lbeta(alpha0 + varSelCandT.N,
                                  beta0 + varSelCand.N - varSelCandT.N) -
                            lbeta(alpha0, beta0))
            }
        }
        out[[model_i]][["betaIdx"]] <- logDens

###----------------------------------------------------------------------------
### Prior for coefficients
###----------------------------------------------------------------------------
        if(tolower(priArgsCurr[["beta"]][["type"]]) == "cond-mvnorm")
        {
            mean <- priArgsCurr[["beta"]][["mean"]] # mean of density
            covariance <- priArgsCurr[["beta"]][["covariance"]] # Covariates
            shrinkage <- priArgsCurr[["beta"]][["shrinkage"]] # Shrinkage

            ## Normal distribution condition normal The full beta vector is assumed as
            ## normal. Since variable selection is included in the MCMC, The final
            ## proposed beta are those non zeros. We need to using the conditional normal
            ## density See Mardia p. 63

            ## Split the beta vector by selected and non-selected.
            Idx1 <- which(betaIdxCurr == TRUE)
            Idx0 <- which(betaIdxCurr == FALSE)

            Idx0Len <- length(Idx0)
            Idx1Len <- length(Idx1)

            ## The mean vector (recycled if necessary)
            meanVec <- matrix(mean, 1, length(betaCurr))

            if(tolower(covariance) == "identity")
            {
                coVar <- diag(shrinkage, nrow = length(betaCurr))
            }
            else if(tolower(covariance) == "g-prior")
            {
                stop("Not implemented.")
                ## The covariance matrix for the whole beta vector
                ## The covariance matrix for the whole beta vector
                ## X <- X[[model_i]]
                ## coVar0 <- qr.solve(crossprod(X))
                ## coVar0Lst <- lapply(as.vector(rep(NA, nPar),"list"),
                ##                     function(x) coVar0)
                ## coVar <- block.diag(coVar0Lst)
            }

            ## Calculate the log density
            if(Idx0Len == 0L)
            {   ## 1. all are selected or Switch to unconditional prior.
                out[[model_i]][["beta"]] <- dmvnorm(matrix(betaCurr, nrow = 1), meanVec, coVar, log = TRUE)
            }
            else
            {  ## 2. some are selected (at least the intercept is always selected ) Using the conditional prior
                A <- coVar[Idx1, Idx0]%*%solve(coVar[Idx0, Idx0])
                condMean <- meanVec[Idx1] - A%*%meanVec[Idx0]
                condCovar <- coVar[Idx1, Idx1] - A%*%coVar[Idx0, Idx1]

                out[[model_i]][["beta"]] <- dmvnorm(matrix(betaCurr[Idx1], nrow = 1), condMean, condCovar, log = TRUE)
            }
        }
        else
        {
            stop("No such prior.")
        }
    }
###----------------------------------------------------------------------------
### Update the output for prior
###----------------------------------------------------------------------------
    if(sum)
    {
        out = sum(unlist(out))
    }
    return(out)
}


#' Gradient for the model
#'
#' Gradient for priors
#' @param beta list for models to be updated
#' @param betaIdx same as beta
#' @param varSelArgs list, same as beta length
#' @param priArgs list, same as beta length
#' @param chainCaller
#' @param sum logical, if TRUE, sum all conditional priors
#' @return "list". The gradient and Hessian matrix
#' @references NA
#' @export
log_priors_grad <- function(beta, betaIdx, varSelArgs, priArgs)
{
    ## Only update priors for parameters that need to update.  Initial the storage
    ## structure for current log prior
    num_models_updated = length(beta)

    ## Reserve list structure of the gradient and Hessian
    out = beta # retain same out structure as priArgs.

    ## Gradient for the intercept as a special case because its prior is different
    priArgsCurr <- priArgsCurr[[model_i]][["beta"]]
    betaCurr <- beta[[model_i]][1] # the intercepts

    if(tolower(priArgs[["type"]]) == "custom")
    {
        ## Call the any2any() function
        densOutput <- any2any(densArgs = priArgs)

        mean <- priArgsCurr[["beta"]][["mean"]] # mean of density
        covariance <- priArgsCurr[["beta"]][["covariance"]] # Covariates
        shrinkage <- priArgsCurr[["beta"]][["shrinkage"]] # Shrinkage

        ## The output
        out[[model_i]][1] <- -1/(covariance * shrinkage) * (betaCur-mean)
    }

    ## Gradient for beta|I for beta (unconditional)
    priArgs <- priArgs[[model_i]][["beta"]]

    ## Slopes and variable selection indicators(taking away intercept)
    betaCurr <- beta[[model_i]][-1]
    betaIdxCurr <- betaIdx[[model_i]][-1]

    for(model_i in 1:num_models)
    {
        if(length(betaCur) > 0L)
        {
            if(tolower(priArgs[["type"]]) == "cond-mvnorm")
            {
                ## Normal distribution condition normal The full beta vector is assumed as
                ## normal. Since variable selection is included in the MCMC, The final
                ## proposed beta are those non zeros. We need to using the gradient for the
                ## conditional normal density See Mardia p. 63.

                ## Subtract the prior information for the full beta
                mean <- priArgsCurr[["mean"]] # mean of density
                covariance <- priArgsCurr[["covariance"]] # Covariates
                shrinkage <- priArgsCurr[["shrinkage"]] # Shrinkage

                ## Split the beta vector by zero and nonzero.
                Idx1 <- which(betaIdxCurr == TRUE)
                Idx0 <- which(betaIdxCurr == FALSE)
                Idx0Len <- length(Idx0)
                Idx1Len <- length(Idx1)
                betaLen <- length(betaIdxCurr)

                SlopCondGrad <- array(NA, dim(betaIdxCurr))

                ## The mean vector for the whole beta vector (recycled if necessary)
                meanVec <- array(mean, dim(betaIdxCurr))

                ## The covariance matrix for the whole beta vector
                if(tolower(covariance) == "g-prior")
                {
                    coVar0 <- qr.solve(crossprod(X))
                    coVar0Lst <- lapply(
                        as.vector(rep(NA, nPar),"list"),
                        function(x) coVar0)
                    coVar <- block.diag(coVar0Lst)
                }
                else if(tolower(covariance) == "identity")
                {
                    coVar <- diag(betaLen)
                }

                ## The conditional gradient. Consider three situations for
                ## gradient.
                if(Idx0Len == 0)
                {
                    ## 1. all are selected. Switch to unconditional prior.
                    SlopCondGrad[Idx1] <- DensGradHess(
                        B = matrix(beta),
                        mean = matrix(meanVec),
                        covariance = coVar*shrinkage,
                        grad = TRUE, Hess = FALSE)[["grad"]]
                }

                else if(Idx0Len > 0 && Idx0Len < betaLen)
                {
                    ## 2. some are selected (the most common situation)
                    ## The conditional prior
                    A <- coVar[Idx1, Idx0, drop = FALSE]%*%
                        solve(coVar[Idx0, Idx0, drop = FALSE])

                    condMean <- meanVec[Idx1] - A%*%meanVec[Idx0]

                    condCovar <- coVar[Idx1, Idx1, drop = FALSE] -
                        A%*%coVar[Idx0, Idx1, drop = FALSE]

                    ## The conditional gradient
                    SlopCondGrad[Idx1] <- DensGradHess(
                        B = matrix(beta[Idx1]),
                        mean = condMean,
                        covariance = condCovar*shrinkage,
                        grad = TRUE, Hess = FALSE)[["grad"]]
                }
                else
                {
                    ## 3. non are selected, do nothing Switch to unconditional prior.
                }

                out[[model_i]][2:nPar] <- SlopCondGrad
            }
        }
    }

    ## The final gradient output.  The intercept and the conditional gradient
    return(out)
}
