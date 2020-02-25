## load data
load("BMA_PP.RData")
load("BMA_FF.RData")
load("M4_q1.RData")

#' Calculate the log predictive score for a time series with pools of models
#'
#'
#' @title log predictive score with features
#' @param beta p-by-(n-1) matrix
#' @param features T-by-p feature matrix, usually standardized.
#' @param prob T-by-n predictive densities, currently n=2.
#' @param intercept TRUE or FALSE Should intercept be used in feature weights?
#' @return
#' @references Geweke & Amisano, (2011) Optimal prediction pools, Journal of Econometrics.
#' @note TODO: log_score_grad(beta, features, prob, intercepts)
#' @author Feng Li
log_score<-function(beta, features, prob, intercept){

    if(intercept) features = cbind(rep(1, nrow(prob)), features)

    exp_lin = exp(features%*%beta)

    w <- exp_lin/(1+rowSums(exp_lin)) # T-by-(n-1)

    ## Full (T-by-n) matrix. To keep identification, only first (n-1) are connected with
    ## features. TODO: Common features?
    w_full = cbind(w, 1 - rowSums(w)) # T-by-n

    out = sum(log(rowSums(w_full * prob)))
    return(out)
}


#################################################################################
## performance of five methods
## 1 feature based (42 features)
## 2 optimal pool
## 3 simple average (SA)
## 4 ARIMA
## 5 ETS

for (a in 1:1000) {

    y<-M4_q1[[a]]$x
    log_pred <- PP[[a]]
    features_y <- FF[[a]]

    ## maximizing TODO: change to a better optimization tool.
    set.seed(2019-02-06)
    w_max<-optim(fn=log_score, par=runif(43, min = 0, max = 0),
                 features = features_y,
                 prob = exp(log_pred),
                 intercept = intercept,
                 method="SANN", control = list(fnscale = -1))

    if(w_max$convergence!=0){
        cat("The optimization does not converge in data", a)
    }
    beta_optim<-w_max$par

    ## optimal pool: feature=NULL, intercept =TRUE
    ## w_optim<-optim(fn=log_score, par=runif(43, min = 0, max = 0),
    ##              features = NULL,
    ##              prob = p,
    ##              intercept = TRUE,
    ##              method="SANN", control = list(fnscale = -1))

    ## forecasting
    y_hat_feature<-c()

    feature_y_hat = matrix(nrow = forecast_h, ncol = 42)
    y_new = y
    for (t in 1:forecast_h)
    {
        ## Calculate predictive features. Each individual model provide an one-step-ahead
        ## predictive y values, and recalculate features based on optimized pools. This
        ## will do the forecasting model multiple times (h), consider to simplify it.

        ## ETS model
        ets_fit <- ets(y_new, model = ets_model)
        ets_fore<-forecast(ets_fit, h = 1, level = PI_level)
        ets_fore_mean <- ets_fore$mean
        ets_fore_sd = (ets_fore$lower - ets_fore$mean)/qnorm(1 - PI_level/100)

        ## ARIMA model
        arima_fit <- auto.arima(y_new)
        arima_fore <- forecast(arima_fit, h = forecast_h, level = PI_level)
        arima_fore_mean <- arima_fore$mean
        arima_fore_sd = (ari_fore$lower - ari_fore$mean)/qnorm(1 - PI_level/100)

        ## Update features
        myts <- list(list(x=ts(y_new, frequency = frequency)))
        myfeatures <- THA_features(my)[[1]]$features
        myfeatures <- data.matrix(myfeatures)
        features_y_hat[t, ] <- myfeatures


        ## Update predictive weights
        if(intercept) features = cbind(1, features)
        exp_lin = exp(features%*%beta)
        w <- exp_lin/(1+rowSums(exp_lin)) # T-by-(n-1)
        w_full = cbind(w, 1 - rowSums(w)) # T-by-n


        log_pred_densities[, 1] <- sum(dnorm(y[(t + 1):(t + h)], mean = ets_fore_mean,
                                             sd = ets_fore_sd, log = TRUE))
        log_pred_densities[, 2] <- sum(dnorm(y[(t + 1):(t + h)], mean = arima_fore_mean,
                                             sd = arima_fore_sd, log = TRUE))



        ## Simple version
        y<-c(y,y_hat)

    }

    M4_q1[[a]]$y_hat_feature<-y_hat_feature
    M4_q1[[a]]$y_hat_w<-w_optim*M4_q1[[a]]$ets_fore+(1-w_optim)*M4_q1[[a]]$ari_fore
    M4_q1[[a]]$average<-0.5*(M4_q1[[a]]$ets_fore+M4_q1[[a]]$ari_fore)

    ## mase smape
    ff<-rbind(M4_q1[[a]]$y_hat_feature,M4_q1[[a]]$y_hat_w,M4_q1[[a]]$average,M4_q1[[a]]$ari_fore,M4_q1[[a]]$ets_fore)
    row.names(ff)<-c("y_hat_feature","y_hat_w","average","ari_fore","ets_fore")
    M4_q1[[a]]$ff<-ff

    lentry <- M4_q1[[a]]
    insample <- lentry$x
    ff<- lentry$ff
    frq<-stats::frequency(insample)
    insample <- as.numeric(insample)
    outsample <- as.numeric(lentry$xx)
    masep <- mean(abs(utils::head(insample, -frq) - utils::tail(insample, -frq)))
    repoutsample <- matrix(rep(outsample, each = nrow(ff)), nrow = nrow(ff))
    smape_err <- 200 * abs(ff - repoutsample)/(abs(ff) + abs(repoutsample))
    mase_err <- abs(ff - repoutsample)/masep
    M4_q1[[a]]$mase_err <- rowMeans(mase_err)
    M4_q1[[a]]$smape_err <- rowMeans(smape_err)

    ##log predictive score
    sigma_ari<-M4_q1[[a]]$sigma_ari
    sigma_ets<-M4_q1[[a]]$sigma_ets

    ari_score_sum<-0
    for(j in 1:8){
        ari_score<-dnorm(outsample[j],mean =ff[4,j],sd=sigma_ari )
        ari_score_sum<-ari_score_sum+log(ari_score)
    }
    ets_score_sum<-0
    for(j in 1:8){
        ets_score<-dnorm(outsample[j],mean =ff[5,j],sd=sigma_ets )
        ets_score_sum<-ets_score_sum+log(ets_score)
    }
    ave_score_sum<-0
    for(j in 1:8){
        ave_score<-0.5*dnorm(outsample[j],mean =ff[4,j],sd=sigma_ari)+0.5*dnorm(outsample[j],mean =ff[5,j],sd=sigma_ets)
        ave_score_sum<-ave_score_sum+log(ave_score)
    }
    hat_w_score_sum<-0
    for(j in 1:8){
        hat_w_score<-(1-w_optim)*dnorm(outsample[j],mean =ff[4,j],sd=sigma_ari)+w_optim*dnorm(outsample[j],mean =ff[5,j],sd=sigma_ets)
        hat_w_score_sum<-hat_w_score_sum+log(hat_w_score)
    }
    hat_feature_score_sum<-0
    for(j in 1:8){
        w<-1/(1+exp(-c(1,features_y[(length(features_y[,1])-9+j),])%*%beta_optim))
        hat_feature_score<-(1-w)*dnorm(outsample[j],mean =ff[4,j],sd=sigma_ari)+w*dnorm(outsample[j],mean =ff[5,j],sd=sigma_ets)
        hat_feature_score_sum<-hat_feature_score_sum+log(hat_feature_score)
    }
    score<-c(hat_feature_score_sum,hat_w_score_sum,ave_score_sum,ari_score_sum,ets_score_sum)
    names(score)<-c("y_hat_feature","y_hat_w","average","ari_fore","ets_fore")
    M4_q1[[a]]$score<-score
}


## MASE
mase_err<-c()
for (j in 1:5){
    mase_err0<-c()
    for (i in 1:1000){
        mase_err0<-c(mase_err0,M4_q1[[i]]$mase_err[j])
    }
    mase_err0<-mean(mase_err0)
    mase_err<-c(mase_err,mase_err0)
}


## SMAPE
smape_err<-c()
for (j in 1:5){
    smape_err0<-c()
    for (i in 1:1000){
        smape_err0<-c(smape_err0,M4_q1[[i]]$smape_err[j])
    }
    smape_err0<-mean(smape_err0)
    smape_err<-c(smape_err,smape_err0)
}


score_output<-c()
for (j in 1:5){
    score0<-c()
    for (i in 1:1000){
        if(i == 173){
            score0<-score0
        }else{
            score0<-c(score0,M4_q1[[i]]$score[j])
        }
    }
    score0<-mean(score0)
    score_output<-c(score_output,score0)
}

summery<-rbind(mase_err,smape_err,score_output)
colnames(summery)<-c("y_hat_feature","y_hat_w","average","ari_fore","ets_fore")
rownames(summery)<-c("mase_err","smape_err","log_score")
summery
