#' Calculate the log joint prior
#'
#' Calculate the log joint prior of  indicator vectors and  coefficient vectors.
#' See \code{model_conf_default} for more information of parameter settings.
#' 
#' @param beta A list of coefficient vectors of the features.
#' @param betaIdx A list of indicator vectors. 
#' @param varSelArgs Variable selection settings.
#' @param priArgs Parameter settings in the priors of \code{beta} and \code{betaIdx}.
#' @param sum Logical, if TRUE, sum all conditional priors.
#' 
#' @return \code{log_priors} returns the value of log joint prior.
#' @export
log_priors <- function(beta, betaIdx, varSelArgs, priArgs, sum = TRUE)
{

    num_models_updated = length(beta)
    out = rapply(priArgs, function (x) 0, how="replace") # retain same out structure as priArgs.

    for(iComp in 1:num_models_updated)
    {
        priArgsCurr <- priArgs[[iComp]]
        betaCurr <- beta[[iComp]]
        betaIdxCurr <- betaIdx[[iComp]] # p-by-1

###----------------------------------------------------------------------------
### Prior for variable selection indicators
###----------------------------------------------------------------------------
        ## intercept as special case. The intercept should always be included.
        varSelCand <- varSelArgs[[iComp]][["cand"]]

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
                candIdx <- 2:length(beta[[iComp]])
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
        out[[iComp]][["betaIdx"]] <- logDens

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
                ## X <- X[[iComp]]
                ## coVar0 <- qr.solve(crossprod(X))
                ## coVar0Lst <- lapply(as.vector(rep(NA, nPar),"list"),
                ##                     function(x) coVar0)
                ## coVar <- block.diag(coVar0Lst)
            }

            ## Calculate the log density
            if(Idx0Len == 0L)
            {   ## 1. all are selected or Switch to unconditional prior.
                out[[iComp]][["beta"]] <- mvtnorm::dmvnorm(matrix(betaCurr, nrow = 1), meanVec, coVar, log = TRUE)
            }
            else
            {  ## 2. some are selected (at least the intercept is always selected ) Using the conditional prior
                A <- coVar[Idx1, Idx0]%*%solve(coVar[Idx0, Idx0])
                condMean <- meanVec[Idx1] - A%*%meanVec[Idx0]
                condCovar <- coVar[Idx1, Idx1] - A%*%coVar[Idx0, Idx1]

                out[[iComp]][["beta"]] <- mvtnorm::dmvnorm(matrix(betaCurr[Idx1], nrow = 1), condMean, condCovar, log = TRUE)
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



# Gradient for priors

log_priors_grad <- function(beta, betaIdx, varSelArgs, priArgs)
{
    num_models_updated = length(betaIdx)
    out = beta

    for(iComp in 1:num_models_updated)
    {
        priArgsCurr <- priArgs[[iComp]]
        betaCurr <- beta[[iComp]]
        betaIdxCurr <- betaIdx[[iComp]] # p-by-1

        ## intercept as special case. The intercept should always be included.
        varSelCand <- varSelArgs[[iComp]][["cand"]]

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
            Idx1 <- which(betaIdxCurr == 1)
            Idx0 <- which(betaIdxCurr == 0)

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
            }

            ## Calculate the gradient w.r.t. log density
            if(Idx0Len == 0L)
            {   ## 1. all are selected or Switch to unconditional prior.
                coVarInv = solve(coVar)
                out[[iComp]] <- -coVarInv %*% (betaCurr - mean)
            }
            else
            {  ## 2. some are selected (at least the intercept is always selected ) Using
                ## the conditional prior
                A <- coVar[Idx1, Idx0]%*%solve(coVar[Idx0, Idx0])
                condMean <- meanVec[Idx1] - A%*%meanVec[Idx0]
                condCovar <- coVar[Idx1, Idx1] - A%*%coVar[Idx0, Idx1]
                coVarInv = solve(condCovar)

                out[[iComp]] <- -coVarInv %*% matrix(betaCurr[Idx1]-meanVec[Idx1])
            }
        }
        else
        {
            stop("No such prior.")
        }
    }

    return(out)
}

