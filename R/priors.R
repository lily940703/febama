#' This is the prior settings for a variable selection scheme
#'
#' The detailed documentation is in the main setting file for each parameter.
#' @title Model prior settings.
#' @param Mdl.X "list"
#' @param Mdl.parLink "list"
#' @param Mdl.beta "list"
#' @param Mdl.betaIdx "list"
#' @param Mdl.varSelArgs "list"
#' @param Mdl.priArgs "list"
#' @param parUpdate "list"
#' @param Mdl.logPri  NA
#' @param priCurr "list"
#' @return "list" synced
#' @references "list"
#' @export
log_priors <- function(beta,  Mdl.beta, Mdl.betaIdx, varSelArgs, priArgs,
                       parUpdate)
{
    ## Loop over all updated parameter candidates
###----------------------------------------------------------------------------
### Only update priors for parameters that need to update.
###----------------------------------------------------------------------------
    parUpdateIdx <- names(parUpdate[[CompCaller]])[unlist(parUpdate[[CompCaller]])]
    for(parCaller in parUpdateIdx)
    {
        ## Initial the storage structure for current log prior
        outCurr <-  Mdl.priArgs[[CompCaller]][[parCaller]]

###----------------------------------------------------------------------------
### Prior for variable selection indicators
###----------------------------------------------------------------------------
        betaIdxCurr <- Mdl.betaIdx[[CompCaller]][[parCaller]] # p-yb-lq
        nPar <- Mdl.parLink[[CompCaller]][[parCaller]][["nPar"]]

        ## Variable section candidates checking.
        if(dim(betaIdxCurr)[1] == 1)
        {
            ## Only intercept is in, do not allow for variable selection.
            candIdx <- NULL
        }
        else
        {
            varSelCandConfigCurr <- Mdl.varSelArgs[[CompCaller]][[parCaller]][["cand"]]

            if(class(varSelCandConfigCurr) == "character" &&
               tolower(varSelCandConfigCurr) == "2:end")
            {
                ncolX.ij <- nrow(Mdl.beta[[CompCaller]][[parCaller]])
                candIdx <- 2:ncolX.ij
            }
            else
            {
                candIdx <- varSelCandConfigCurr
            }

        }

        ## Variable selection prior probability.
        if(length(candIdx)>0)
        {
            Mdl.priArgsCurr <- Mdl.priArgs[[CompCaller]][[parCaller]][["indicators"]]
            varSelCandTF <- betaIdxCurr[candIdx, , drop = FALSE]

            if(tolower(Mdl.priArgsCurr[["type"]]) == "bern") # Bernoulli prior
            {
                ## Note that independent Bernoulli prior can be very informative. See
                ## Scott and Berger 2010 in Annals of Statistics Sec 3.2.

                ## The probability is set for each variable involved in the variable
                ## selection procedure via "valSelArgs" The probability is recycled if
                ## necessary.
                prob <- Mdl.priArgsCurr[["prob"]]


                probMat <- matrix(prob, length(candIdx), nPar)
                ## TRUE or FALSE of variable selection candidates

                logDens <- sum(dbinom(x = as.numeric(varSelCandTF), size = 1,
                                      prob = probMat, log = TRUE))
            }
            else if(tolower(Mdl.priArgsCurr[["type"]]) == "beta") # beta prior
            {
                ## Scott and Berger 2010 in Annals of Statistics Sec 3.3
                alpha = Mdl.priArgsCurr[["alpha"]]
                beta  = Mdl.priArgsCurr[["beta"]]

                varSelCand.N <- length(candIdx)
                varSelCandT.N <- sum(varSelCandTF)
                logDens <- (lbeta(alpha + varSelCandT.N,
                                  beta + varSelCand.N - varSelCandT.N) -
                            lbeta(alpha, beta))
            }


        }
        else
        {
            ## No variable selection is used, thus we don't need prior.
            logDens <- NULL
        }
        outCurr[["indicators"]] <- logDens

###----------------------------------------------------------------------------
### Prior for coefficients
###----------------------------------------------------------------------------

### intercept as special case. The intercept should alway be included.
        Mdl.priArgsCurr <- Mdl.priArgs[[CompCaller]][[parCaller]][["beta"]][["intercept"]]

        betaCurr <- Mdl.beta[[CompCaller]][[parCaller]][1,,drop = FALSE]#intercepts

        ## if(is(betaCurr, "try-error")) browser()

        linkCurr <- Mdl.parLink[[CompCaller]][[parCaller]]

        if(tolower(Mdl.priArgsCurr[["type"]]) == "custom")
        {
            ## Call the any2any() function
            densOutput <- any2any(densArgs = Mdl.priArgsCurr, linkArgs = linkCurr)

            mean <- densOutput$mean
            variance <- densOutput$variance

            ## cat(mean, "\t", variance, "\n")

            shrinkage <- Mdl.priArgsCurr[["output"]][["shrinkage"]]

            logDens <- dnorm(x = betaCurr, mean = mean,
                             sd = sqrt(variance*shrinkage), log = TRUE)
            outCurr[["beta"]][["intercept"]] <- logDens
        }

### coefficients (conditional on variable selection indicators)
        Mdl.priArgsCurr <- Mdl.priArgs[[CompCaller]][[parCaller]][["beta"]][["slopes"]]

        ## Slopes and variable selection indicators(taking away intercept)
        betaCurr <- Mdl.beta[[CompCaller]][[parCaller]][-1, , drop = FALSE]
        betaIdxNoIntCurr <- Mdl.betaIdx[[CompCaller]][[parCaller]][-1,,drop = FALSE]


        if(length(betaIdxNoIntCurr) == 0L)
        {
            ## No covariates at all (only intercept in the model)
            outCurr[["beta"]][["slopes"]] <- NULL
        }
        else
        {
            if(tolower(Mdl.priArgsCurr[["type"]]) == "cond-mvnorm")
            {
                ## Normal distribution condition normal The full beta vector is assumed as
                ## normal. Since variable selection is included in the MCMC, The final proposed
                ## beta are those non zeros. We need to using the conditional normal density See
                ## Mardia p. 63

                ## Subtract the prior information for the full beta
                mean <- Mdl.priArgsCurr[["mean"]] # mean of density
                covariance <- Mdl.priArgsCurr[["covariance"]] # Covariates
                shrinkage <- Mdl.priArgsCurr[["shrinkage"]] # Shrinkage

                ## Split the beta vector by selected and non-selected.
                Idx1 <- which(betaIdxNoIntCurr == TRUE)
                Idx0 <- which(betaIdxNoIntCurr == FALSE)

                Idx0Len <- length(Idx0)
                Idx1Len <- length(Idx1)

                betaLen <- length(betaIdxNoIntCurr)

                ## The mean vector (recycled if necessary)
                meanVec <- matrix(mean, 1, betaLen)

                ## The covariance matrix for the whole beta vector
                if(tolower(covariance) == "g-prior")
                {
                    ## The covariance matrix for the whole beta vector
                    X <- Mdl.X[[CompCaller]][[parCaller]][, -1, drop = FALSE]
                    coVar0 <- qr.solve(crossprod(X))
                    coVar0Lst <- lapply(as.vector(rep(NA, nPar),"list"),
                                        function(x) coVar0)

                    coVar <- block.diag(coVar0Lst)

                }
                else if(tolower(covariance) == "identity")
                {
                    coVar <- diag1(betaLen)
                }

                ## Calculate the log density
                if(Idx0Len == 0L)
                {   ## 1. all are selected or
                    ## Switch to unconditional prior.
                    outCurr[["beta"]][["slopes"]] <- dmvnorm(matrix(betaCurr, 1, betaLen,byrow = TRUE),
                                                             meanVec, coVar*shrinkage, log = TRUE)
                }
                else if( Idx0Len == betaLen)
                {   ## 2. non are selected:
                    ## Switch to unconditional prior.
                    outCurr[["beta"]][["slopes"]] <- dmvnorm(matrix(betaCurr, 1, betaLen, byrow = TRUE),
                                                             meanVec, coVar*shrinkage, log = TRUE)
                }
                else if(Idx0Len > 0 & Idx0Len < betaLen)
                {
                    ## 3. some are selected
                    ## Using the conditional prior
                    A <- coVar[Idx1, Idx0]%*%solve(coVar[Idx0, Idx0])
                    condMean <- meanVec[Idx1] - A%*%meanVec[Idx0]
                    condCovar <- coVar[Idx1, Idx1] - A%*%coVar[Idx0, Idx1]

                    outCurr[["beta"]][["slopes"]] <- dmvnorm(matrix(betaCurr[Idx1], 1),
                                                             condMean, condCovar*shrinkage, log = TRUE)
                }
                else
                {
                    stop("Debug me: Unknown situation for conditional priors.")
                }
            }
        }

###----------------------------------------------------------------------------
### Update the output for prior
###----------------------------------------------------------------------------
        Mdl.logPri[[CompCaller]][[parCaller]] <- outCurr
    }

    return(Mdl.logPri)
}


#' Gradient for the model
#'
#' Gradient for priors
#' @param Mdl.X NA
#' @param Mdl.beta NA
#' @param Mdl.betaIdx NA
#' @param Mdl.parLink NA
#' @param Mdl.varSelArgs NA
#' @param Mdl.priArgs NA
#' @param chainCaller NA
#' @return "list". The gradient and Hessian matrix
#' @references NA
#' @export
log_priors_grad <- function(Mdl.X, Mdl.beta, Mdl.betaIdx,Mdl.parLink,
                            Mdl.varSelArgs, Mdl.priArgs, chainCaller)
{
    ## Only update priors for parameters that need to update.
    ## Initial the storage structure for current log prior
    CompCaller <- chainCaller[[1]]
    parCaller <- chainCaller[[2]]

    ## Reserve list structure of the gradient and Hessian
    gradObsLst <- list()
    HessObsLst <- list()

###----------------------------------------------------------------------------
### Gradient and Hessian for the intercept as a special case
###----------------------------------------------------------------------------
    Mdl.priArgsCurr <- Mdl.priArgs[[CompCaller]][[parCaller]][["beta"]][["intercept"]]
    betaCurr <- Mdl.beta[[CompCaller]][[parCaller]][1, , drop = FALSE] # the intercepts
    linkCurr <- Mdl.parLink[[CompCaller]][[parCaller]]

    if(tolower(Mdl.priArgsCurr[["type"]]) == "custom")
    {
        ## Call the any2any() function
        densOutput <- any2any(densArgs = Mdl.priArgsCurr, linkArgs = linkCurr)

        mean <- densOutput$mean
        ## print(densOutput$variance)
        variance <- diag(densOutput$variance, length(betaCurr))
        shrinkage <- Mdl.priArgsCurr[["output"]][["shrinkage"]]

        ## Gradient and Hessian for the intercept
        GradHessInt <- DensGradHess(B = matrix(betaCurr),
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

    Mdl.priArgsCurr <- Mdl.priArgs[[CompCaller]][[parCaller]][["beta"]][["slopes"]]
    nPar <- Mdl.parLink[[CompCaller]][[parCaller]][["nPar"]]

    ## Slopes and variable selection indicators(taking away intercept)
    betaCurr <- Mdl.beta[[CompCaller]][[parCaller]][-1, , drop = FALSE]
    betaIdxNoIntCurr <- Mdl.betaIdx[[CompCaller]][[parCaller]][-1, , drop = FALSE]

    X <- Mdl.X[[CompCaller]][[parCaller]][, -1, drop = FALSE]
    if(length(X) == 0L)
    {
        ## No covariates at all (only intercept in the model)
        gradObsLst[["Slop"]] <- NULL
        HessObsLst[["SlopFull"]] <- NULL
    }
    else
    {
        if(tolower(Mdl.priArgsCurr[["type"]]) == "cond-mvnorm")
        {
            ## Normal distribution condition normal The full beta vector is assumed
            ## as normal. Since variable selection is included in the MCMC, The
            ## final proposed beta are those non zeros. We need to using the
            ## gradient for the conditional normal density See Mardia p. 63.

            ## Subtract the prior information for the full beta
            mean <- Mdl.priArgsCurr[["mean"]] # mean of density
            covariance <- Mdl.priArgsCurr[["covariance"]] # Covariates
            shrinkage <- Mdl.priArgsCurr[["shrinkage"]] # Shrinkage

            ## Split the beta vector by zero and nonzero.
            Idx1 <- which(betaIdxNoIntCurr == TRUE)
            Idx0 <- which(betaIdxNoIntCurr == FALSE)
            Idx0Len <- length(Idx0)
            Idx1Len <- length(Idx1)
            betaLen <- length(betaIdxNoIntCurr)

            SlopCondGrad <- array(NA, dim(betaIdxNoIntCurr))

            ## The mean vector for the whole beta vector (recycled if necessary)
            meanVec <- array(mean, dim(betaIdxNoIntCurr))

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
                    B = matrix(betaCurr),
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
                    B = matrix(betaCurr[Idx1]),
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
                ##     B = betaCurr,
                ##     mean = meanVec,
                ##     covariance = coVar*shrinkage,
                ##     grad = TRUE, Hess = FALSE)[["grad"]]
                ## SlopCondGrad[Idx1] <- NA
            }

            gradObsLst[["Slop"]] <- SlopCondGrad

            ## The unconditional full Hessian matrix
            ## HessObsLst[["SlopFull"]] <- DensGradHess(
            ##     B = matrix(betaCurr),
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
