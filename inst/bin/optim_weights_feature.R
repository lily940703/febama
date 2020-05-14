#################################################################################
## performance of four methods
## 1 optimal pool
## 2 simple average (SA)
## 3 ARIMA
## 4 ETS

require("forecast")
require("M4metalearning")
library(doParallel)
## load historical feature data
#setwd("~/code/febama")
load("E:/time series/R code/feature-based-bayesian-model-averaging/data/historical_lpd_feature_yearly.RData")


#----------------------------------------------------------#
## Function  optim_beta_feature
## param: lpd_feature for a ts with $lpd and $feat from the previous step
## param: features_select (a vector including the numbers of selected features)
## param: intercept = FALSE
## return: a list with beta and the number of the selected feature

optim_beta_feature<-function( lpd_feature, features_select,intercept = FALSE){

  y_lpd <- lpd_feature$lpd
  features_y = lpd_feature$feat

  prob <- exp(y_lpd)
  if(sum(prob==0)!=0){
    cat("The probability predictive densities have 0 value in data",i_ts);
  }
  prob0 <-prob[rowSums(prob==0)!=2,]
  features_y<-features_y[rowSums(prob==0)!=2,]

  ## maximizing TODO: change to a better optimization tool.
  w_max <- try(optim(par = rep(0, 43),
                 fn = logscore,
                 gr = gradient,
                 features = features_y,
                 prob = prob0,
                 features_select = features_select,
                 intercept = intercept,
                 method="BFGS",
                 control = list(fnscale = -1)) #max
               )

  #if(is(w_max, "try-error")) browser()
  beta_optim <- w_max$par
  return(list(beta_optim = beta_optim, value = w_max$value, features_select = features_select))
}

cl <- makeCluster(2)
registerDoParallel(cl)
optimal_beta_feature <- foreach(i_ts = 1:length(data)) %dopar%
  optim_beta_feature(lpd_feature_yearly[[i_ts]],
                     features_select=c(10,12,16,17,40), intercept = TRUE)
stopCluster(cl)


#----------------------------------------------------------#
## Function  lpd_feature
## param: a ts data with $x
## return: lpd_feature (a list including two matrices of probability densities and features)
lpd_feature <- function(data, model_conf) {

  feature_window = model_conf$feature_window
  roll = model_conf$roll
  frequency = model_conf$frequency
  history_burn = model_conf$history_burn
  ets_model = model_conf$ets_model
  forecast_h = model_conf$forecast_h
  train_h = model_conf$train_h
  PI_level = model_conf$PI_level

  ## A single historical data
  y <- data$x
  y1 = scale(y, center = TRUE, scale = TRUE)
  y_mean = attr(y1, "scaled:center")
  y_sd = attr(y1, "scaled:scale")

  y1 = as.numeric(y1)

  ## Calculate historical log predictive density
  log_pred_densities <-
    matrix(nrow = length(y) - history_burn - train_h + 1,
           ncol = 2)

  for (t in (history_burn):(length(y) - train_h))
  {
    ## TODO: We may simplify this to assume the fit and forecast procedure is
    ## invariant with certain time period to save time.
    if(is.null(roll)){
      y01 <- y1[1:t]
    }else if (t < roll){
      y01 <- y1[1:t]
    }else{
      y01 <- y1[(t-roll+1):t]
    }

    ## ETS model
    ets_fit <- ets(y01, model = ets_model)
    ets_fore <- forecast(ets_fit, h = train_h, level = PI_level)
    ets_fore_mean <- ets_fore$mean
    ## forecast.ets does not directly provide predictive variance but we could infer from
    ## predict interval.  Remember that PI = pred_mean -/+ qnorm(PI_level)*pred_sd
    ets_fore_sd = (ets_fore$lower - ets_fore$mean) / qnorm(1 - PI_level /
                                                             100)

    ## ARIMA model
    arima_fit <- auto.arima(y01)
    arima_fore <- forecast(arima_fit, h = train_h, level = PI_level)
    arima_fore_mean <- arima_fore$mean
    arima_fore_sd = (arima_fore$lower - arima_fore$mean) / qnorm(1 - PI_level /
                                                                   100)

    ## To keep numeric stability, we calculate log P(y_pred)
    log_pred_densities[(t - history_burn + 1), 1] <-
      sum(dnorm(
        y1[(t + 1):(t + train_h)],
        mean = ets_fore_mean,
        sd = ets_fore_sd,
        log = TRUE
      ))
    log_pred_densities[(t - history_burn + 1), 2] <-
      sum(dnorm(
        y1[(t + 1):(t + train_h)],
        mean = arima_fore_mean,
        sd = arima_fore_sd,
        log = TRUE
      ))
  }

  ## Calculate historical features
  features_y <-
    matrix(nrow = length(y) - train_h - history_burn + 1,
           ncol = 42)
  myts <- list(list(x = ts(y[1:history_burn], frequency = frequency)))
  myfeatures <- THA_features(myts)[[1]]$features
  names <- colnames(myfeatures)
  colnames(features_y) <- names

  if(is.null(feature_window)){
    for (t in (history_burn):(length(y) - train_h))
    {
      myts <- list(list(x = ts(y[1:t], frequency = frequency)))
      myfeatures <- THA_features(myts)[[1]]$features
      myfeatures <- data.matrix(myfeatures)
      features_y[(t - history_burn + 1),] <- myfeatures
    }
  }else{
    features_y <- features_window( y, window = feature_window,
                                   history_burn, train_h)
  }

  features_y_scaled = scale(features_y, center = TRUE, scale = TRUE)

  lpd_feature <-
    list(lpd = log_pred_densities, feat = features_y_scaled)
  return(lpd_feature)
}

#----------------------------------------------------------#

# delete the features with NaN
feature_clean <- function(lpd_feature){
  for (i in 1:length(lpd_feature)) {
    ind <- which(is.nan(lpd_feature[[i]]$feat), arr.ind = T)[,2]
    ind0 <- ind[!duplicated(ind)]
    lpd_feature[[i]]$feat_mean <- attr(lpd_feature[[i]]$feat, "scaled:center")[-ind0]
    lpd_feature[[i]]$feat_sd <- attr(lpd_feature[[i]]$feat, "scaled:scale")[-ind0]
    lpd_feature[[i]]$feat <- lpd_feature[[i]]$feat[, -ind0]
  }
  return(lpd_feature)
}

#----------------------------------------------------------#
## SGLD + VS

logscore <- function(beta, features, features_select = NULL, prob, intercept){
  if(is.null(features_select)){
    features0<-features
  }else{
    features0<-data.matrix(features[,features_select])
  }
  if(intercept) features0 = cbind(rep(1, nrow(prob)), features0)
  # 避免Inf/Inf
  me <- features0 %*% beta
  me[me>709] <- 709
  exp_lin = exp(me)
  w <- exp_lin/(1+rowSums(exp_lin)) # T-by-(n-1)
  w_full = cbind(w, 1 - rowSums(w)) # T-by-n
  out = sum(log(rowSums(w_full * prob)))
  return(out)
}

gradient <- function(beta, features, features_select = NULL, prob, intercept){
  if(is.null(features_select)){
    features0<-features
  }else{
    features0<-data.matrix(features[,features_select])
  }
  if(intercept) features0 = cbind(rep(1, nrow(prob)), features0)

  ex = exp(features0 %*% beta)
  ex[ex > exp(700)] = exp(700)

  ex_sum = rowSums(ex)
  n = dim(prob)[2]
  gradient0<-function(i, t, p){
    out0 = (
      ex[t,i] * (
        p[t,i] *(1 + ex_sum[t] - ex[t,i])
        - sum(p[t,-c(i,n)] * ex[t,-i]) - p[t,n]
      ) * features0[t,] )/ (
        (sum(ex[t,] * p[t, 1:(n-1)]) + p[t,n] ) * (1 + ex_sum[t])
      )
    return(out0)
  }

  out_t<-c()
  out_i<-c()
  for (i in 1:(n-1)) {
    for (t in 1:length(prob[,1])) {
      out_t<-cbind(out_t, gradient0(i, t, prob))
    }
    out_i<-cbind(out_i, rowSums(out_t))
  }
  return(out_i)
}



SGLD <- function(data, logLik, gradient_logLik, prior, start, minibatchSize = NULL,
                 stepsize = NULL, tol = 1e-5, iter = 5000, samplesize = 0.1,sig = 10,
                 features_select, I, intercept = TRUE, gama = 0.55, a = 0.4, b = 10 ){
  beta <- start

  prob <- exp(data$lpd)
  # when pd=0, delete or assign a minimum?
  prob[prob == 0] <- 1e-323
  features <- data$feat

  prior0 <- prior(beta, I, sig = sig)

  logLik0 <- logLik (beta = beta, features = features,
                     features_select = features_select,
                     prob = prob, intercept = intercept)
  logpost0 <- log_posterior (data, beta, I, prior = prior,
                             logLik = logscore, sig = sig)
  i <- 1
  res <- list(beta = start, logscore = logLik0, logpost = logpost0,
              stepsize = NA, prior = prior0)

  if(! is.null(minibatchSize)){
    minibatchSize = minibatchSize
  } else if (length(prob[,1]) <= 10){
    minibatchSize = 1
  }else {
    minibatchSize = 0.1
  }

  repeat{
    mini <- sample(1:length(prob[,1]),
                   ceiling(minibatchSize*length(prob[,1])))
    prob1 <- prob[mini,]
    features1 <- features[mini,]

    if(is.null(stepsize)){
      stepsize1 <- a * (b + i) ^ (-gama)
    }else{
      stepsize1 <- stepsize
    }

    prior1 <- prior(as.vector(beta), I, sig = sig)

    beta <- (beta + stepsize1 * (as.vector(jacobian(func = dmvnorm , x = as.vector(beta), method="Richardson",
                                                    mean=matrix(0, length(beta), 1),sigma = sig*diag(length(beta))))
                                 / prior1)
             + stepsize1 * (1/minibatchSize) * gradient_logLik(beta = beta, features = features1,
                                                               features_select = features_select,
                                                               prob= prob1, intercept= intercept)
             + mvrnorm(1, rep(0,length(beta)), 2*stepsize1* diag(length(beta)))
    )

    prior1 <- prior(as.vector(beta), I, sig = sig)
    logLik1 <- logLik (beta = beta, features = features, features_select = features_select,
                       prob = prob, intercept = intercept)
    logpost1 <- log_posterior (data, beta, I, prior = prior,
                               logLik = logscore, sig = sig)

    res$beta <- rbind(res$beta, t(beta))
    res$logscore <- rbind(res$logscore, logLik1)
    res$logpost <- rbind(res$logpost, logpost1)
    res$stepsize <- rbind(res$stepsize, stepsize1)
    res$prior <- rbind(res$prior, prior1)

    if (i >= iter)
      break
    i <- i+1
  }

  # n = samplesize * iter
  # beta_SGLD <- tail(res$beta, n)

  # simple average
  #beta_out <- colMeans(beta_SGLD)

  # weighted
  # weights <- tail(res$stepsize, n) / sum(tail(res$stepsize, n))
  # beta_out <- colSums(beta_SGLD * matrix(rep(weights,length(start)), ncol = length(start)))

  return(res)
}


proposal_I <- function(n) rbinom(n, 1, 0.5)

# p(beta | I) * p(I)
prior <- function(beta, I, sig = 10){
  pri <- (
    dmvnorm(as.vector(beta), mean=matrix(0, length(beta), 1), sigma = sig*diag(length(beta)))
    # * 0.5 ^ length(I)
  )
  # avoid prior=0
  if (pri < 1e-323) {pri <- 1e-323}
  return(pri)
}


log_posterior <- function(data, beta, I, prior = prior,
                          logLik = logscore, sig = 10){
  pri <- prior(beta, I, sig = sig)
  features_select <- which(I==1)
  prob <- exp(data$lpd)
  prob[prob == 0] <- 1e-323
  features <- data$feat
  features[is.nan(features)] <- 0
  LS <- logLik(beta = beta, features = features, features_select = features_select,
               prob = prob, intercept = TRUE)
  log_post <- log(pri) + LS
  return(log_post)
}


MH_step <- function(x, beta0, data, logp = log_posterior,
                    proposal = proposal_I, sig = 10){
  beta_start <- as.vector(beta0)
  names(beta_start) <- c("0", which(x == 1))
  xp <- proposal(length(x))
  beta1 <- matrix(NA, nrow = (length(which(xp == 1)) + 1), ncol = 1)
  names(beta1) <- c("0", which(xp == 1))
  ind <- c("0", which(x==xp & xp==1))
  beta1[ind] <- beta_start[ind]
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


SGLD_VS <- function(data, logLik, gradient_logLik, prior, stepsize = NULL,
                    SGLD_iter = 100, VS_iter = 100, minibatchSize = NULL, sig = 10){
  feature_num <- dim(data$feat)[2]
  I <- matrix(nrow = feature_num, ncol = VS_iter)
  B <- list()
  result_all <- list()
  I[,1] <- rbinom(feature_num,1,0.5)
  accept_num <- 0
  for (i in 1:VS_iter) {
    if(i == 1){
      features_select <- which(I[,i]==1)
      beta_start <- runif(length(features_select) + 1, -10, 10)
      names(beta_start) <- c("0", features_select)
    }else{
      accept_num <- accept_num +accept
      I[,i] <-I0
      features_select <- which(I[,i]==1)
      beta_start <- beta_start
    }
    res_SGLD <- SGLD(data = data, logLik = logscore, gradient_logLik = gradient,
                     prior = prior, start = beta_start, I = I[,i],
                     minibatchSize = minibatchSize, stepsize = stepsize,
                     iter = SGLD_iter, features_select = features_select, sig = sig )
    B[[i]] <- res_SGLD$beta
    result_all[[i]] <- res_SGLD
    I0 <- I[,i]
    MH <- MH_step (x = I0, beta0 = t(tail(res_SGLD$beta,1)),
                   data =data, logp = log_posterior, proposal = proposal_I)
    I0 <- MH$I
    beta_start <- MH$beta_start
    accept <- MH$accept
  }
  return(list(I = I, B = B, result_all = result_all, acceptance = accept_num/VS_iter))
}

#Achieve the mean of beta and the selected features in each VS
beta_prepare <- function(samples_SGLD_VS, num){
  beta_samples <- lapply(samples_SGLD_VS$B, tail, n = num)
  beta_mean <- lapply(beta_samples, colMeans)
  beta_pre0 <-list()
  beta_pre <-list()
  for (i in 1 : length(samples_SGLD_VS$B)) {
    beta_pre$beta_mean <- beta_mean[[i]]
    beta_pre$features_select <- as.numeric(names(beta_mean[[i]]))[-1]
    beta_pre0[[i]] <- beta_pre
  }
  return(beta_pre0)
}

#----------------------------------------------------------#
## Function  forecast_results
## param: a ts data with $x, $xx
## param: model_conf
## param: intercept
## param: lpd_feature with $feat, $feat_mean, $feat_sd
## param: beta_pre, a list achieved from beta_prepare function

## return: data[[i_ts]] with $ff_feature, $err_feature

forecast_feature_results <-function(data, model_conf, intercept = T, lpd_feature, beta_pre) {

  #attach(model_conf)
  feature_window = model_conf$feature_window
  roll = model_conf$roll
  frequency = model_conf$frequency
  history_burn = model_conf$history_burn
  ets_model = model_conf$ets_model
  forecast_h = model_conf$forecast_h
  train_h = model_conf$train_h
  PI_level = model_conf$PI_level

  ## forecasting
  y_hat_matrix <- matrix(ncol = forecast_h, nrow = 1)

  y <- data$x
  y01 = scale(y, center = TRUE, scale = TRUE)
  y_mean = attr(y01, "scaled:center")
  y_sd = attr(y01, "scaled:scale")
  y01 = as.numeric(y01)
  y_true = data$xx
  y01_true = as.numeric(scale(y_true, center = y_mean, scale = y_sd))
  y_new = y01
  y_new_nonsd = as.numeric(y)

  features_y = lpd_feature$feat
  features_y_mean = lpd_feature$feat_mean
  features_y_sd = lpd_feature$feat_sd

  lpds = 0
  pred_densities = matrix(NA, forecast_h, 2)
  w_full_mean_h <- c()
  #forecast_h = model_conf$forecast_h
  w_full_all <- list()
  for (t in 1:forecast_h)
  {
    ## NOTE: This part is a recursive process, could not be parallelized.  Calculate
    ## predictive features. Each individual model provide an one-step-ahead predictive y
    ## values, and recalculate features based on optimized pools. This will do the
    ## forecasting model multiple times (h), consider to simplify it.

    ## ETS model
    if(is.null(roll)){
      y_new1 <- y_new
    }else{
      y_new1 <- tail(y_new, roll)
    }
    ets_fit <- ets(y_new1, model = model_conf$ets_model)
    ets_fore <- forecast(ets_fit, h = 1, level = PI_level)
    ets_fore_mean <- ets_fore$mean
    ets_fore_sd = (ets_fore$lower - ets_fore$mean) / qnorm(1 - PI_level /
                                                             100)

    ## ARIMA model
    arima_fit <- auto.arima(y_new1)
    arima_fore <- forecast(arima_fit, h = 1, level = PI_level)
    arima_fore_mean <- arima_fore$mean
    arima_fore_sd = (arima_fore$lower - arima_fore$mean) / qnorm(1 - PI_level /
                                                                   100)

    ## Update features
    y_new_nonsd1 <- tail(y_new_nonsd, feature_window)
    if (!is.null(features_y))
    {
      #用非标准化的数据计算特征
      myts <- list(list(x = ts(y_new_nonsd1, frequency = frequency)))
      myfeatures <- THA_features(myts)[[1]]$features
      myfeatures <- data.matrix(myfeatures)
      myfeatures <- myfeatures[, colnames(myfeatures) %in% colnames(features_y)]
      myfeatures_scaled = scale(t(myfeatures),
                                center = features_y_mean, scale = features_y_sd)
    } else
    {
      myfeatures_scaled = NULL
    }

    ## Update predictive weights
    w_get <- function(beta_pre, myfeatures_scaled){
      myfeatures_scaled <- myfeatures_scaled[, beta_pre$features_select]
      if (intercept){
        myfeatures_scaled = cbind(1, t(myfeatures_scaled))
      }
      exp_lin = exp(myfeatures_scaled %*% beta_pre$beta_mean)
      #  避免Inf
      exp_lin[exp_lin > exp(709)] <- exp(709)
      w <- exp_lin / (1 + rowSums(exp_lin))
      w_full = cbind(w, 1 - rowSums(w))
      return(w_full)
    }
    w_full <- sapply(beta_pre, w_get, myfeatures_scaled = myfeatures_scaled)
    w_full_all[[t]] <- w_full
    w_full_mean <- rowMeans(w_full)
    w_full_mean_h <- cbind(w_full_mean_h, w_full_mean)

    ## The final pooled y
    y_pred_h = sum(cbind(ets_fore_mean, arima_fore_mean) * w_full_mean)
    y_new = c(y_new, y_pred_h)
    y_new_nonsd = c(y_new_nonsd, (y_pred_h * y_sd + y_mean))
    y_hat_matrix[1, t] <- y_pred_h * y_sd + y_mean

    ### The predictive log score
    pred_densities[t, 1] <-
      dnorm(y01_true[t], mean = ets_fore_mean, sd = ets_fore_sd)
    pred_densities[t, 2] <-
      dnorm(y01_true[t], mean = arima_fore_mean, sd = arima_fore_sd)
    lpds0 = log(sum(pred_densities[t, ] * w_full_mean))
    if(lpds0 < log(1e-323)){
      lpds0 = log(1e-323)
    }else{
      lpds0 =lpds0
    }
    lpds = lpds + lpds0
  }
  data$ff_feature <- y_hat_matrix
  colnames(w_full_mean_h) <- seq(1, forecast_h, 1)
  data$w_time_varying <- w_full_mean_h
  data$w_detail <- w_full_all

  ### mase smape
  ff <- data$ff_feature
  insample <- as.numeric(data$x)
  frq <- stats::frequency(insample)
  outsample <- as.numeric(data$xx)
  masep <-
    mean(abs(utils::head(insample,-frq) - utils::tail(insample,-frq)))
  smape_err <- 200 * abs(ff - outsample) / (abs(ff) + abs(outsample))
  mase_err <- abs(ff - outsample) / masep
  mase_err_h <- rowMeans(mase_err)
  smape_err_h <- rowMeans(smape_err)

  data$err_feature <- cbind(lpds, mase_err_h, smape_err_h)
  return(data)
}

#----------------------------------------------------------#
## Function  forecast_performance
## param:data with $err_feature
## return: a matrix performance

load("E:/time series/R code/feature-based-bayesian-model-averaging/data/data_forecast_feature_yearly.RData")
forecast_feature_performance<-function(data){
  # log score
  performance<-c()
  for (i_ts in 1:length(data)) {
    performance<-rbind(performance,data[[i_ts]]$err_feature)
  }
  performance_out1<- sum(performance[,1])

  performance_out2<- colMeans(performance[,2:3])

  performance_out<-cbind(performance_out1,t(performance_out2))
  colnames(performance_out)<-c("Log score","Mase","Smape")
  return(performance_out)
}

forecast_feature_performance(data_forecast_feature_SGLD10)
forecast_feature_performance(data_forecast_feature)

#------------------------------------------------------------------------------------#
## Usage example (long ts)

library(tsfeatures)
library(M4metalearning)
library(forecast)
library(tseries)
library(purrr)
library(ggplot2)
library(M4comp2018)
library(foreach)
library(doParallel)
library(numDeriv)
library("mvtnorm")
library("base")
library(MASS)


# ts

# Daily (length = 8000+)
# $h
# [1] 14
plot(M4[[97406]]$x)
plot(M4[[99458]]$x)

# weekly (length = 2500+)
# $h
# [1] 13
plot(M4[[95004]]$x)
plot(M4[[95224]]$x)


#-------------------------------------------------------------------------#
model_conf1 = list(
  frequency = 1
  , ets_model = "ANN" # simple exponential smoothing with additive errors
  , forecast_h = 14 # Forecasting horizon
  , train_h = 1 # Out of sample training, 1 means rolling forecast
  , history_burn = 1000 # Let the model to start with at least this length of historical data.
  , PI_level = 90 # Predictive Interval level, used to extract out-of-sample variance from forecasting models.
  , roll = NULL # The length of rolling samples, larger than history_burn
  , feature_window =1000 # The length of moving window when computing features
)

data_test <- M4[[97406]]
lpd_feature1 <- lpd_feature(data_test, model_conf1)
lpd_feature1 <- list(lpd_feature1 = lpd_feature1)
lpd_feature_test <- feature_clean(lpd_feature1)[[1]]
num_fea <- length(lpd_feature_test$feat_mean)
set.seed(2020-4-14)
optim_noVS <- SGLD (lpd_feature_test, logLik = logscore, gradient_logLik = gradient,
                    prior = prior, start = runif(num_fea+1, -10, 10), minibatchSize = 0.1,
                    stepsize = 0.1, iter = 500, sig = 10,
                    features_select = seq(1,num_fea,1), I = rep(1,num_fea), intercept = TRUE )
beta_samples <- tail(optim_noVS$beta,100)
beta_mean <- colMeans(beta_samples)
beta_pre <-list()
beta_pre$beta_mean <- beta_mean
beta_pre$features_select <- seq(1,num_fea,1)
beta_pre_noVS <-list()
beta_pre_noVS[[1]] <- beta_pre

set.seed(2020-4-14)
optim <- SGLD_VS(lpd_feature_test, logLik = logscore, gradient_logLik = gradient,
                 prior = prior,  SGLD_iter = 500, VS_iter = 100)
feature_select <- which(rowMeans(optim$I) > 0.5)
I <- as.numeric(rowMeans(optim$I) > 0.5)
beta_pre <- beta_prepare (samples_SGLD_VS = optim, num = 100)

set.seed(2020-4-15)
optim_1 <- SGLD_VS(lpd_feature_test, logLik = logscore, gradient_logLik = gradient,
                   prior = prior,  SGLD_iter = 500, VS_iter = 100)
set.seed(2020-4-16)
optim_2 <- SGLD_VS(lpd_feature_test, logLik = logscore, gradient_logLik = gradient,
                   prior = prior,  SGLD_iter = 500, VS_iter = 100)
set.seed(2020-4-17)
optim_3 <- SGLD_VS(lpd_feature_test, logLik = logscore, gradient_logLik = gradient,
                   prior = prior,  SGLD_iter = 500, VS_iter = 100)
set.seed(2020-4-18)
optim_4 <- SGLD_VS(lpd_feature_test, logLik = logscore, gradient_logLik = gradient,
                   prior = prior,  SGLD_iter = 500, VS_iter = 100)


## Forecasting with time varying weights (rolling samples)
data_test_fore <- forecast_feature_results (data_test, model_conf = model_conf1,
                                            intercept = T,
                                            lpd_feature = lpd_feature_test , beta_pre)
data_test_fore_noVS <- forecast_feature_results (data_test, model_conf1, intercept = T,
                                                 lpd_feature = lpd_feature_test , beta_pre_noVS)


# Analysis------------------------------------------------------------------#
save(model_conf1, data_test,lpd_feature_test,
     optim_noVS, beta_pre_noVS, optim, beta_pre,
     optim_1, optim_2, optim_3, optim_4,
     data_test_fore, data_test_fore_noVS, data_optim,
     file="E:/time series/R code/febama/data/test_daily_SGLD(VS).RData")
load("E:/time series/R code/febama/data/test_daily_SGLD(VS).RData")

# optimization
optim_beta <- optim_beta(lpd_feature_test, features_y = NULL)
# $beta_optim
# [1] 1.062258
#
# $logscore
# [1] 18296.25
w <- exp(optim_beta$beta_optim)/(1+exp(optim_beta$beta_optim))
w_full = cbind(w, 1 - w)
#[1,] 0.7431218 0.2568782

par(mfrow = c(3,1))
plot(optim_noVS$logscore)
abline(h=18296.25, col = "2")
plot(optim_noVS$logpost)
plot(log(optim_noVS$prior))


par(mfrow = c(4,4))
for (i in 1:16) {
  plot(optim$result_all[[i+84]]$logscore, ylab = "log score",
       xlab = paste0('iterations in  no.', i+84, '  VS'),
       ylim = c(18250,18320))
  abline(h=18296.25, col = "2")
}

# error
data_test_fore_noVS$err_feature
data_test_fore$err_feature

# > data_test_fore_noVS$err_feature
#          lpds mase_err_h smape_err_h
# [1,] 6.782679   2.256713     1.02534
# > data_test_fore$err_feature
#          lpds mase_err_h smape_err_h
# [1,] 20.42067   1.740449   0.7920355

data_optim <- forecast_results(data_test, forecast_h=14, optim_beta)
data_optim$mase_err
# Optimal pool 1.602584
# SA           1.661617
# ETS          1.528228
# ARIMA        1.756524
data_optim$smape_err
# Optimal pool 0.7295635
# SA           0.7563267
# ETS          0.6958333
# ARIMA        0.7993238
data_optim$logscore
# Optimal pool  24.44480
# SA            24.14793
# ETS           24.77821
# ARIMA         23.59292

# forecast results
autoplot(data_test$x)+
  xlab("day") + ylab("values")

autoplot(tail(data_test$x,50)) +
  autolayer(data_test$xx, series="true value",size = 1.5) +
  autolayer(ts(as.vector(data_test_fore$ff_feature), start = 13612, end = 13625 ),
            series="SGLD+VS",size = 2) +
  autolayer(ts(as.vector(data_test_fore_noVS$ff_feature), start = 13612, end = 13625 ),
            series="SGLD",size = 1.5) +
  autolayer(ts(as.vector(data_optim$ff[1,]), start = 13612, end = 13625 ),
            series="Optimal pool",size = 1.5) +
  autolayer(ts(as.vector(data_optim$ff[2,]), start = 13612, end = 13625 ),
            series="SA",size = 1.5) +
  autolayer(ts(as.vector(data_optim$ff[3,]), start = 13612, end = 13625 ),
            series="ETS ",size = 1.5) +
  autolayer(ts(as.vector(data_optim$ff[4,]), start = 13612, end = 13625 ),
            series="ARIMA  ",size = 1.5) +
  #ggtitle("") +
  xlab("day") + ylab("values") +
  guides(colour=guide_legend(title="Forecast"))

## whether the features chosen with high probability is stable



optim_1 <- SGLD_VS(lpd_feature_test, logLik = logscore, gradient_logLik = gradient,
                 prior = prior,  SGLD_iter = 500, VS_iter = 100)
optim_2 <- SGLD_VS(lpd_feature_test, logLik = logscore, gradient_logLik = gradient,
                   prior = prior,  SGLD_iter = 500, VS_iter = 100)
optim_3 <- SGLD_VS(lpd_feature_test, logLik = logscore, gradient_logLik = gradient,
                   prior = prior,  SGLD_iter = 500, VS_iter = 100)
optim_4 <- SGLD_VS(lpd_feature_test, logLik = logscore, gradient_logLik = gradient,
                   prior = prior,  SGLD_iter = 500, VS_iter = 100)


par(mfrow = c(1,1))
plot(rowMeans(optim$I), ylab = "The probability of being selected",
     xlab = "Features", type = "l", ylim = c(0,1), lwd = 2)
lines(rowMeans(optim_1$I), col = "5", lwd = 2)
lines(rowMeans(optim_2$I), col = "3", lwd = 2)
lines(rowMeans(optim_3$I), col = "4", lwd = 2)
lines(rowMeans(optim_4$I), col = "6", lwd = 2)

I_mean <- rbind(rowMeans(optim$I), rowMeans(optim_1$I), rowMeans(optim_2$I),
                rowMeans(optim_3$I), rowMeans(optim_4$I))
lines(colMeans(I_mean), col = "2", lwd = 4)
