#' This is the prior settings for a variable selection scheme
#'
#' The detailed documentation is in the main setting file for each parameter.
#' @title Model prior settings.
#' @param beta list for models to be updated
#' @param betaIdx same as beta
#' @param varSelArgs
#' @param priArgs
#' @param parUpdate "list"
#' @return "list" synced
#' @references "list"
#' @export
log_priors <- function(beta, betaIdx, varSelArgs, priArgs, sum = TRUE)
{
    num_models_updated = length(beta)
    out = rapply(priArgs, function (x) 0, how="replace") # retain same out structure as priArgs.

    for(model_i in 1:num_models_updated)
    {
        priArgsCur <- priArgs[[model_i]]
###----------------------------------------------------------------------------
### Prior for variable selection indicators
###----------------------------------------------------------------------------
        betaIdx <- betaIdx[[model_i]] # p-by-1
        varSelCand <- varSelArgs[[model_i]][["cand"]]

        ## Variable section candidates checking. If (1) only intercept is in or (2)
        ## updating candidates is NULL, do not allow for variable selection. We do not
        ## need a prior;
        if(length(betaIdx) == 1 | length(varSelCand) == 0) # intercept is always included.
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
            varSelCandTF <- betaIdx[candIdx]

            if(tolower(priArgsCur[["indicators"]][["type"]]) == "bern") # Bernoulli prior
            {
                ## Note that independent Bernoulli prior can be very informative. See
                ## Scott and Berger 2010 in Annals of Statistics Sec 3.2.
                prob <- priArgs[["prob"]]
                probMat <- matrix(prob, length(candIdx), 1) # TRUE or FALSE of variable selection candidates
                logDens <- sum(dbinom(x = as.numeric(varSelCandTF), size = 1,
                                      prob = probMat, log = TRUE))
            }
            else if(tolower(priArgsCur[["indicators"]][["type"]]) == "beta") # beta prior
            {  ## Scott and Berger 2010 in Annals of Statistics Sec 3.3
                alpha0 = priArgsCur[["indicators"]][["alpha"]]
                beta0  = priArgsCur[["indicators"]][["beta"]]

                varSelCand.N <- length(candIdx)
                varSelCandT.N <- sum(varSelCandTF)
                logDens <- (lbeta(alpha + varSelCandT.N,
                                  beta + varSelCand.N - varSelCandT.N) -
                            lbeta(alpha, beta))
            }
        }
        out[[model_i]][["indicators"]] <- logDens

###----------------------------------------------------------------------------
### Prior for coefficients
###----------------------------------------------------------------------------

        ## intercept as special case. The intercept should always be included.
        betaCur <- beta[[model_i]]#intercepts

        if(tolower(priArgsCur[["indicators"]][["type"]]) == "norm")
        {
            mean <- priArgsCur[["indicators"]][["mean"]] # mean of density
            covariance <- priArgsCur[["indicators"]][["covariance"]] # Covariates
            shrinkage <- priArgsCur[["indicators"]][["shrinkage"]] # Shrinkage

            logDens <- dnorm(x = beta[1], mean = mean,
                             sd = sqrt(variance*shrinkage), log = TRUE)
            out[[model_i]][["beta"]][["intercept"]] <- logDens
        }
        else
        {
            stop("No such prior.")
        }

        ## coefficients (conditional on variable selection indicators)
        priArgsCur <- priArgs[[model_i]]

        ## Slopes and variable selection indicators(taking away intercept)
        betaCur <- beta[[model_i]][-1]
        betaIdxNoInt <- betaIdx[[model_i]][-1]

        if(length(betaIdxNoInt) == 0L)
        { ## No covariates at all (only intercept in the model)
            out[[model_i]][["beta"]][["slopes"]] <- NULL
        }
        else
        {
            if(tolower(priArgsCur[["slopes"]][["type"]]) == "cond-mvnorm")
            {
                ## Normal distribution condition normal The full beta vector is assumed as
                ## normal. Since variable selection is included in the MCMC, The final
                ## proposed beta are those non zeros. We need to using the conditional
                ## normal density See Mardia p. 63

                ## Subtract the prior information for the full beta
                betaCur = beta[[model_i]][-1]
                mean <- priArgsCur[["slopes"]][["mean"]] # mean of density
                covariance <- priArgsCur[["slopes"]][["covariance"]] # Covariates
                shrinkage <- priArgsCur[["slopes"]][["shrinkage"]] # Shrinkage

                ## Split the beta vector by selected and non-selected.
                Idx1 <- which(betaIdxNoInt == TRUE)
                Idx0 <- which(betaIdxNoInt == FALSE)

                Idx0Len <- length(Idx0)
                Idx1Len <- length(Idx1)

                betaLen <- length(betaIdxNoInt)
                ## The mean vector (recycled if necessary)
                meanVec <- matrix(mean, 1, betaLen)

                ## The covariance matrix for the whole beta vector
                if(tolower(covariance) == "g-prior")
                {
                    ## The covariance matrix for the whole beta vector
                    X <- X[[model_i]][, -1, drop = FALSE]
                    coVar0 <- qr.solve(crossprod(X))
                    coVar0Lst <- lapply(as.vector(rep(NA, nPar),"list"),
                                        function(x) coVar0)

                    coVar <- block.diag(coVar0Lst)

                }
                else if(tolower(covariance) == "identity")
                {
                    coVar <- diag(betaLen)
                }

                ## Calculate the log density
                if(Idx0Len == 0L)
                {   ## 1. all are selected or Switch to unconditional prior.
                    out[[model_i]][["beta"]][["slopes"]] <- dmvnorm(matrix(betaCur, nrow = 1),
                                                                    meanVec, coVar*shrinkage, log = TRUE)
                }
                else if( Idx0Len == betaLen)
                {   ## 2. non are selected: Switch to unconditional prior.
                    out[[model_i]][["beta"]][["slopes"]] <- dmvnorm(matrix(beta, nrow = 1),
                                                                    meanVec, coVar*shrinkage, log = TRUE)
                }
                else if(Idx0Len > 0 & Idx0Len < betaLen)
                {
                    ## 3. some are selected Using the conditional prior
                    A <- coVar[Idx1, Idx0]%*%solve(coVar[Idx0, Idx0])
                    condMean <- meanVec[Idx1] - A%*%meanVec[Idx0]
                    condCovar <- coVar[Idx1, Idx1] - A%*%coVar[Idx0, Idx1]

                    out[[model_i]][["beta"]][["slopes"]] <- dmvnorm(matrix(beta[Idx1], nrow = 1),
                                                                    condMean, condCovar*shrinkage, log = TRUE)
                }
                else
                {
                    stop("No such conditional priors.")
                }
            }
        }
###----------------------------------------------------------------------------
### Update the output for prior
###----------------------------------------------------------------------------
        logPri[[model_i]] <- out
    }

    if(sum)
    {
        out = sum(unlist(logPri))
    }
    else
    {
        out = logPri
    }

    return(out)
}


#' Gradient for the model
#'
#' Gradient for priors
#' @param X NA
#' @param beta NA
#' @param betaIdx NA
#' @param parLink NA
#' @param varSelArgs NA
#' @param priArgs NA
#' @param chainCaller NA
#' @return "list". The gradient and Hessian matrix
#' @references NA
#' @export
log_priors_grad <- function(X, beta, betaIdx,parLink,
                            varSelArgs, priArgs, chainCaller)
{
    ## Only update priors for parameters that need to update.
    ## Initial the storage structure for current log prior
    CompCaller <- chainCaller[[1]]
    model_i <- chainCaller[[2]]

    ## Reserve list structure of the gradient and Hessian
    gradObsLst <- list()
    HessObsLst <- list()

###----------------------------------------------------------------------------
### Gradient and Hessian for the intercept as a special case
###----------------------------------------------------------------------------
    priArgs <- priArgs[[model_i]][["beta"]][["intercept"]]
    beta <- beta[[model_i]][1, , drop = FALSE] # the intercepts
    link <- parLink[[model_i]]

    if(tolower(priArgs[["type"]]) == "custom")
    {
        ## Call the any2any() function
        densOutput <- any2any(densArgs = priArgs)

        mean <- densOutput$mean
        ## print(densOutput$variance)
        variance <- diag(densOutput$variance, length(beta))
        shrinkage <- priArgs[["output"]][["shrinkage"]]

        ## Gradient and Hessian for the intercept
        GradHessInt <- DensGradHess(B = matrix(beta),
                                    mean = mean,
                                    covariance = variance*shrinkage,
                                    grad = TRUE, Hess = FALSE)

        ## The output
        gradObsLst[["Int"]] <- t(GradHessInt[["grad"]]) # row-vector
                                        # HessObsLst[["Int"]] <- GradHessInt[["Hess"]]
    }
###----------------------------------------------------------------------------
### Gradient for beta|I and Hessian for beta (unconditional)
###----------------------------------------------------------------------------

    priArgs <- priArgs[[model_i]][["beta"]][["slopes"]]
    nPar <- parLink[[model_i]][["nPar"]]

    ## Slopes and variable selection indicators(taking away intercept)
    beta <- beta[[model_i]][-1, , drop = FALSE]
    betaIdxNoInt <- betaIdx[[model_i]][-1, , drop = FALSE]

    X <- X[[model_i]][, -1, drop = FALSE]
    if(length(X) == 0L)
    {
        ## No covariates at all (only intercept in the model)
        gradObsLst[["Slop"]] <- NULL
        HessObsLst[["SlopFull"]] <- NULL
    }
    else
    {
        if(tolower(priArgs[["type"]]) == "cond-mvnorm")
        {
            ## Normal distribution condition normal The full beta vector is assumed
            ## as normal. Since variable selection is included in the MCMC, The
            ## final proposed beta are those non zeros. We need to using the
            ## gradient for the conditional normal density See Mardia p. 63.

            ## Subtract the prior information for the full beta
            mean <- priArgs[["mean"]] # mean of density
            covariance <- priArgs[["covariance"]] # Covariates
            shrinkage <- priArgs[["shrinkage"]] # Shrinkage

            ## Split the beta vector by zero and nonzero.
            Idx1 <- which(betaIdxNoInt == TRUE)
            Idx0 <- which(betaIdxNoInt == FALSE)
            Idx0Len <- length(Idx0)
            Idx1Len <- length(Idx1)
            betaLen <- length(betaIdxNoInt)

            SlopCondGrad <- array(NA, dim(betaIdxNoInt))

            ## The mean vector for the whole beta vector (recycled if necessary)
            meanVec <- array(mean, dim(betaIdxNoInt))

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
                ## 3. non are selected, do nothing
                ##  Switch to unconditional prior.
                ## browser()
                ## SlopCondGrad[Idx1] <- DensGradHess(
                ##     B = beta,
                ##     mean = meanVec,
                ##     covariance = coVar*shrinkage,
                ##     grad = TRUE, Hess = FALSE)[["grad"]]
                ## SlopCondGrad[Idx1] <- NA
            }

            gradObsLst[["Slop"]] <- SlopCondGrad

            ## The unconditional full Hessian matrix
            ## HessObsLst[["SlopFull"]] <- DensGradHess(
            ##     B = matrix(beta),
            ##     mean = matrix(meanVec),
            ##     covariance = coVar*shrinkage,
            ##     grad = FALSE, Hess = TRUE)[["Hess"]]
        }
    }
###----------------------------------------------------------------------------
### The output
###----------------------------------------------------------------------------
    ## The final gradient output.
    ## The intercept and the conditional gradient; The unconditional Hessian

    gradObs <- rbind(gradObsLst[["Int"]], gradObsLst[["Slop"]])

    ## TESTING code to reorganize the hesssian matrix, failed.
    HessObs <- NA
    ## HessObs <- block.diag(HessObsLst)

    ## IdxMat <- matrix(1:15, 5, 3)

    ## IdxSlop <- IdxMat[-1, , drop = FALSE]
    ## IdxInt <- IdxMat[1, ]
    ## HessTest <- matrix(0, 15, 15)

    ## HessTest[IdxSlop, IdxSlop] <- HessObsLst[["SlopFull"]]
    ## HessTest[IdxInt, IdxInt] <- HessObsLst[["Int"]]

    out <- list(gradObs = gradObs, HessObs = HessObs)
    return(out)
}
