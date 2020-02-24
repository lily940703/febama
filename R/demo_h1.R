## feature-based-Bayesian-forecasting-model-averaging
setwd("~/code/febama/")

library(tsfeatures)
## library(M4metaresults)
## library(M4comp2018)
load("data/M4.RData")
## library(M4metalearning)
library(forecast)
library(tseries)

##data(m4 quarterly data)
M4_q<-M4[23001:47000]
set.seed(2019-12-31)
indices <- sample(length(M4_q))
M4_q1 <- M4_q[indices[1:1000]]

## FF: 1000 matrix of features
## PP: 1000 matrix of probability density
## M4_q1: add $ets_fore,$ari_fore, $sigma_ets, $sigma_ari
## save
FF<-list()
PP<-list()

###----------------------------------------------------------------------------
### Model Settings

ets_model = "ANN" # simple exponential smoothing with additive errors
h = 8
PI_level = 90
intercept = FALSE # Do not include intercept in the features
###----------------------------------------------------------------------------
for (a in 1:1000) {
    ## determine ets
    ## ets1<-ets(M4_q1[[a]]$x,model="ZZZ")
    ## ets_model<-ets1$components
    ## # only additive
    ## if(ets_model[1]=="M"){
    ##   ets_model[1]<-"A"
    ## }else{
    ##   ets_model[1]<-ets_model[1]
    ## }
    ## if(ets_model[2]=="M"){
    ##   ets_model[2]<-"A"
    ## }else{
    ##   ets_model[2]<-ets_model[2]
    ## }
    ## if(ets_model[3]=="M"){
    ##   ets_model[3]<-"A"
    ## }else{
    ##   ets_model[3]<-ets_model[3]
    ## }
    ## ets_model<- paste(ets_model[1:3], collapse = "")
    ets2<-ets(M4_q1[[a]]$x,model=ets_model)
    ets_fore<-forecast(ets2, h = h, level = PI_level)
    M4_q1[[a]]$ets_fore<-ets_fore$mean

    ## forecast.ets does not directly provide predictive variance but we could infer from
    ## predict interval.  Remember that PI = pred_mean -/+ qnorm(PI_level)*pred_sd
    ets_fore_sd = (ets_fore$lower - ets_fore$mean)/qnorm(1 - PI_level/100)

    ## determine arima
    arima1<-auto.arima(M4_q1[[a]]$x)
    ari_fore<-forecast(arima1, h = h, level = PI_level)
    M4_q1[[a]]$ari_fore<-ari_fore$mean

    ari_fore_sd = (ari_fore$lower - ari_fore$mean)/qnorm(1 - PI_level/100)

    y<-M4_q1[[a]]$x
    p<-matrix(nrow = length(y),ncol = 2)

    ## p[t,1] of model ets
    sigma_ets<-sqrt(ets2$sigma2)
    M4_q1[[a]]$sigma_ets<-sigma_ets
    for (t in 1:(length(y))) {
        p[t,1]<-dnorm(y[t],mean =ets2$fitted[t],sd=ets_fore_sd)
    }

    ## p[t,2] of model arima
    sigma_ari<-sqrt(arima1$sigma2)
    M4_q1[[a]]$sigma_ari<-sigma_ari
    for (t in 1:(length(y))) {
        p[t,2]<-dnorm(y[t],mean =arima1$fitted[t],sd=ari_fore_sd)
    }
    PP[[a]]<-p

    features_y<-c()
    for (i in 9:length(y)) {
        ## 时间序列长度大于两个周期
        ts<-list(x=ts(y[1:i],frequency = 4))
        ts1<-list(ts)
        myfeatures <- THA_features(ts1)[[1]]$features
        myfeatures<-data.matrix(myfeatures)
        features_y<-rbind(features_y,myfeatures)
    }
    FF[[a]]<-features_y

    ## Time series plot of features
    ## par(mfrow = c(6, 7), mar = c(5, 0, 0, 0))
    ## for(i in 1:42)
    ## {
    ##     plot(features_y[, i], type = "l", col = "red", xlab = colnames(features_y)[i])
    ## }

}


save(PP,file="BMA_PP.RData")
save(FF,file="BMA_FF.RData")
save(M4_q1,file="M4_q1.RData")


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
#' @author Feng Li
log_score<-function(beta, features, prob, intercept){
    if(intercept) features = cbind(1, features)

    exp_lin = exp(features%*%beta)

    w <- exp_lin/(1+rowSums(exp_lin)) # T-by-(n-1)

    ## Full (T-by-n) matrix. To keep identification, only first (n-1) are connected with
    ## features. TODO: Common features?
    w_full = cbind(w, 1 - rowSums(w))

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
    p<-PP[[a]]
    features_y<-FF[[a]]

    ## feature-based method (all 42)
    ##大于两个周期才能算特征，所以前九期权重赋值0.5
    ## log_score1 <- log(0.5*p[1,1]+0.5*p[1,2])+log(0.5*p[2,1]+0.5*p[2,2])+log(0.5*p[3,1]+0.5*p[3,2])+log(0.5*p[4,1]+0.5*p[4,2])+log(0.5*p[5,1]+0.5*p[5,2])+log(0.5*p[6,1]+0.5*p[6,2])+log(0.5*p[7,1]+0.5*p[7,2])+log(0.5*p[8,1]+0.5*p[8,2])+log(0.5*p[9,1]+0.5*p[9,2])
    ## log_score1 = sum(log(rowSums(p[1:9, ])))

    ## maximizing
    set.seed(2019-02-06)
    w_max<-optim(fn=log_score, par=runif(43, min = 0, max = 0),
                 features = features_y,
                 prob = p,
                 intercept = intercept,
                 method="SANN", control = list(fnscale = -1))


    if(w_max$convergence!=0){
        cat("The optimization does not converge in data", a)
    }
    beta_optim<-w_max$par

    ## optimal pool
    y1<-0
    y_score<-function(x){
        for (i in 1:(length(y))) {
            y1<-y1+log(x*p[i,1]+(1-x)*p[i,2])
        }
        return(-y1)
    }
    w_optim<-optimize(y_score,c(0,1),tol = 0.0001)
    w_optim<-w_optim$minimum

    ## forecasting
    y_hat_feature<-c()
    for (i in 1:8) {
        w<-1/(1+exp(-c(1,features_y[length(features_y[,1]),])%*%beta_optim))
        y_hat<-w*(M4_q1[[a]]$ets_fore[i])+(1-w)*(M4_q1[[a]]$ari_fore[i])
        y<-c(y,y_hat)
        ts<-list(x=ts(y[1:length(y)],frequency = 1))
        ts1<-list(ts)
        features <- THA_features(ts1)[[1]]$features
        features<-data.matrix(features)
        features_y<-rbind(features_y,features)
        y_hat_feature<-c(y_hat_feature,y_hat)
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

## #################################################################################
##     ## performance of the feature-based method with different features

##     ## select features to be considered
## feature_select<-list()
## feature_select[[1]]<-10
## feature_select[[2]]<-12
## feature_select[[3]]<-16
## feature_select[[4]]<-17
## feature_select[[5]]<-33
## feature_select[[6]]<-40
## feature_select[[7]]<-c(17,33)
## feature_select[[8]]<-c(12,17,33)
## feature_select[[9]]<-c(10,12,17,33)
## feature_select[[10]]<-c(10,12,16,17,33)
## feature_select[[11]]<-c(10,12,16,17,33,40)

## ## f: feature_select[[f]]
## ## a: M4_q1[[a]]
## for (f in 4:11) {
##     for (a in 1:1000) {
##         y<-M4_q1[[a]]$x
##         p<-PP[[a]]
##         features_y<-FF[[a]]
##         log_score1<-log(0.5*p[1,1]+0.5*p[1,2])+log(0.5*p[2,1]+0.5*p[2,2])+log(0.5*p[3,1]+0.5*p[3,2])+log(0.5*p[4,1]+0.5*p[4,2])+log(0.5*p[5,1]+0.5*p[5,2])+log(0.5*p[6,1]+0.5*p[6,2])+log(0.5*p[7,1]+0.5*p[7,2])+log(0.5*p[8,1]+0.5*p[8,2])+log(0.5*p[9,1]+0.5*p[9,2])
##         log_score<-function(beta){
##             for (i in 10:(length(y))) {
##     ## include intercept term
##                 w<-1/(1+exp(-c(1,features_y[i-9,feature_select[[f]]])%*%beta))
##     ## no intercept term  改四处 optim(...,par,...)，3个w
##     ##w<-1/(1+exp(-features_y[i-9,feature_select[[f]]]%*%beta))
##                 log_score1<-log_score1+log(w*p[i,1]+(1-w)*p[i,2])
##             }
##             return(-log_score1)
##         }

##     ## optimization
##         set.seed(2019-02-14)
##         w_max<-optim(fn=log_score,par=runif((length(feature_select[[f]])+1), 0, 0),method="SANN")
##         if(w_max$convergence!=0){
##             print(a)
##         }
##         beta_optim<-w_max$par

##     ## forcasting
##         y_hat_feature<-c()
##         for (i in 1:8) {
##             w<-1/(1+exp(-c(1,features_y[length(features_y[,1]),feature_select[[f]]])%*%beta_optim))
##             y_hat<-w*(M4_q1[[a]]$ets_fore[i])+(1-w)*(M4_q1[[a]]$ari_fore[i])
##             y<-c(y,y_hat)
##             ts<-list(x=ts(y[1:length(y)],frequency = 4))
##             ts1<-list(ts)
##             features <- THA_features(ts1)[[1]]$features
##             features<-data.matrix(features)
##             features_y<-rbind(features_y,features)
##             y_hat_feature<-c(y_hat_feature,y_hat)
##         }
##         M4_q1[[a]]$y_hat_feature<-y_hat_feature

##     ## mase, smape
##         lentry <- M4_q1[[a]]
##         insample <- lentry$x
##         ff<- lentry$y_hat_feature
##         frq<-stats::frequency(insample)
##         insample <- as.numeric(insample)
##         outsample <- as.numeric(lentry$xx)
##         masep <- mean(abs(utils::head(insample, -frq) - utils::tail(insample, -frq)))
##         smape_err <- 200 * abs(ff - outsample)/(abs(ff) + abs(outsample))
##         mase_err <- abs(ff - outsample)/masep
##         M4_q1[[a]]$mase_err_F <- mean(mase_err)
##         M4_q1[[a]]$smape_err_F <- mean(smape_err)

##     ##log predictive score
##         hat_feature_score_sum<-0
##         for(j in 1:8){
##             w<-1/(1+exp(-c(1,features_y[(length(features_y[,1])-9+j),feature_select[[f]]])%*%beta_optim))
##             hat_feature_score<-(1-w)*dnorm(outsample[j],mean =M4_q1[[a]]$ari_fore[j],sd=M4_q1[[a]]$sigma_ari)+w*dnorm(outsample[j],mean = M4_q1[[a]]$ets_fore[j],sd=M4_q1[[a]]$sigma_ets)
##             hat_feature_score_sum<-hat_feature_score_sum+log(hat_feature_score)
##         }
##         M4_q1[[a]]$score_F<-hat_feature_score_sum
##     }

##     ## The overall predicted performance of the dataset
##     mase_err0<-c()
##     for (i in 1:1000){
##         mase_err0<-c(mase_err0,M4_q1[[i]]$mase_err_F)
##     }
##     mase_err<-mean(mase_err0)

##     smape_err0<-c()
##     for (i in 1:1000){
##         smape_err0<-c(smape_err0,M4_q1[[i]]$smape_err_F)
##     }
##     smape_err<-mean(smape_err0)

##     score0<-c()
##     for (i in 1:1000){
##     ## -Inf in  M4_q1[[173]]$score_F
##     ## delete 173
##         if(i==173){
##             score0<-score0
##         }else{
##             score0<-c(score0,M4_q1[[i]]$score_F)
##         }
##     }
##     score_output<-mean(score0)

##     summery<-cbind(mase_err,smape_err,score_output)
##     colnames(summery)<-c("mase","smape","log score")
##     cat("features:",feature_select[[f]])
##     print(summery)
## }

## ####################################################################################
## ## one step prediction
## for (a in 1:1000) {
##     y<-M4_q1[[a]]$x
##     p<-PP[[a]]
##     features_y<-FF[[a]]

##     ##feature based (feature_select[[8]]:12,17,33)
##     log_score1<-log(0.5*p[1,1]+0.5*p[1,2])+log(0.5*p[2,1]+0.5*p[2,2])+log(0.5*p[3,1]+0.5*p[3,2])+log(0.5*p[4,1]+0.5*p[4,2])+log(0.5*p[5,1]+0.5*p[5,2])+log(0.5*p[6,1]+0.5*p[6,2])+log(0.5*p[7,1]+0.5*p[7,2])+log(0.5*p[8,1]+0.5*p[8,2])+log(0.5*p[9,1]+0.5*p[9,2])
##     log_score<-function(beta){
##         for (i in 10:(length(y))) {
##             w<-1/(1+exp(-features_y[i-9,feature_select[[8]]]%*%beta))
##             log_score1<-log_score1+log(w*p[i,1]+(1-w)*p[i,2])
##         }
##         return(-log_score1)
##     }
##     set.seed(2019-02-14)
##     w_max<-optim(fn=log_score,par=runif(3, 0, 0),method="SANN")
##     beta_optim<-w_max$par
##     w<-1/(1+exp(-features_y[length(features_y[,1]),feature_select[[8]]]%*%beta_optim))
##     y_hat_feature<-w*(M4_q1[[a]]$ets_fore[1])+(1-w)*(M4_q1[[a]]$ari_fore[1])
##     M4_q1[[a]]$y_hat_feature<-y_hat_feature
##     M4_q1[[a]]$w_feature3<-w

##     ##optimal pool
##     y1<-0
##     y_score<-function(x){
##         for (i in 1:(length(y))) {
##             y1<-y1+log(x*p[i,1]+(1-x)*p[i,2])
##         }
##         return(-y1)
##     }
##     w_optim<-optimize(y_score,c(0,1),tol = 0.0001)
##     w_optim<-w_optim$minimum
##     M4_q1[[a]]$y_hat_w<-w_optim*M4_q1[[a]]$ets_fore[1]+(1-w_optim)*M4_q1[[a]]$ari_fore[1]

##     ##SA
##     M4_q1[[a]]$average<-0.5*(M4_q1[[a]]$ets_fore[1]+M4_q1[[a]]$ari_fore[1])

##     ## mase smape
##     ff<-rbind(M4_q1[[a]]$y_hat_feature,M4_q1[[a]]$y_hat_w,M4_q1[[a]]$average,M4_q1[[a]]$ari_fore[1],M4_q1[[a]]$ets_fore[1])
##     row.names(ff)<-c("y_hat_feature","y_hat_w","average","ari_fore","ets_fore")
##     M4_q1[[a]]$ff<-ff

##     lentry <- M4_q1[[a]]
##     insample <- lentry$x
##     ff<- lentry$ff
##     frq<-stats::frequency(insample)
##     insample <- as.numeric(insample)
##     outsample <- as.numeric(lentry$xx)
##     masep <- mean(abs(utils::head(insample, -frq) - utils::tail(insample, -frq)))
##     repoutsample <- matrix(rep(outsample, each = nrow(ff)), nrow = nrow(ff))[,1]
##     smape_err <- 200 * abs(ff - repoutsample)/(abs(ff) + abs(repoutsample))
##     mase_err <- abs(ff - repoutsample)/masep
##     M4_q1[[a]]$mase_err <- mase_err
##     M4_q1[[a]]$smape_err <- smape_err

##     ## log score
##     ari_score<-dnorm(outsample[1],mean =ff[4,1],sd=M4_q1[[a]]$sigma_ari )
##     ets_score<-dnorm(outsample[1],mean =ff[5,1],sd=M4_q1[[a]]$sigma_ets )
##     ave_score<-0.5*dnorm(outsample[1],mean =ff[4,1],sd=M4_q1[[a]]$sigma_ari)+0.5*dnorm(outsample[1],mean =ff[5,1],sd=M4_q1[[a]]$sigma_ets)
##     hat_w_score<-(1-w_optim)*dnorm(outsample[1],mean =ff[4,1],sd=M4_q1[[a]]$sigma_ari)+w_optim*dnorm(outsample[1],mean =ff[5,1],sd=M4_q1[[a]]$sigma_ets)
##     hat_feature_score<-(1-w)*dnorm(outsample[1],mean =ff[4,1],sd=M4_q1[[a]]$sigma_ari)+w*dnorm(outsample[1],mean =ff[5,1],sd=M4_q1[[a]]$sigma_ets)
##     score<-c(hat_feature_score,hat_w_score,ave_score,ari_score,ets_score)
##     names(score)<-c("y_hat_feature","y_hat_w","average","ari_fore","ets_fore")
##     M4_q1[[a]]$score<-score

## }

## mase_err<-c()
## for (j in 1:5){
##     mase_err0<-c()
##     for (i in 1:1000){
##         mase_err0<-c(mase_err0,M4_q1[[i]]$mase_err[j])
##     }
##     mase_err0<-mean(mase_err0)
##     mase_err<-c(mase_err,mase_err0)
## }

## smape_err<-c()
## for (j in 1:5){
##     smape_err0<-c()
##     for (i in 1:1000){
##         smape_err0<-c(smape_err0,M4_q1[[i]]$smape_err[j])
##     }
##     smape_err0<-mean(smape_err0)
##     smape_err<-c(smape_err,smape_err0)
## }

## score_output<-c()
## for (j in 1:5){
##     score0<-c()
##     for (i in 1:1000){
##         score0<-c(score0,log(M4_q1[[i]]$score[j]))
##     }
##     score0<-mean(score0)
##     score_output<-c(score_output,score0)
## }

## summery<-rbind(mase_err,smape_err,score_output)
## colnames(summery)<-c("feature based","optimal pool","SA","ARIMA","ETS")
## rownames(summery)<-c("mase","smape","log score")
## summery
