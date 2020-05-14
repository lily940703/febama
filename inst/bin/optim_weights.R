#' Calculate the log predictive score for a time series with pools of models
#'
#'
#' @title log predictive score with features
#' @param beta p-by-(n-1) matrix
#' @param features T-by-p feature matrix, usually standardized.
#' @param prob T-by-n predictive densities, currently n=2.
#' @param intercept TRUE or FALSE Should intercept be used in feature weights?
#' @param features_select a vector including the numbers of the features to be taken into consideration
#' @return
#' @references Geweke & Amisano, (2011) Optimal prediction pools, Journal of Econometrics.
#' @note TODO: logscore_grad(beta, features, prob, intercepts)
#' @author Feng Li
#'
logscore<-function(beta, features, features_select = NULL, prob, intercept){
  if(is.null(features_select)){
    features0<-features
  }else{
    features0<-data.matrix(features[,features_select])
  }
  if(intercept) features0 = cbind(rep(1, nrow(prob)), features0)
  exp_lin = exp(features0 %*% beta)
  w <- exp_lin/(1+rowSums(exp_lin)) # T-by-(n-1)
  w_full = cbind(w, 1 - rowSums(w)) # T-by-n
  out = sum(log(rowSums(w_full * prob)))
  return(out)
}

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
## Function  optim_beta
## param:lpd_feature with $lpd from the previous step
## param:features_y = NULL(don't take features into consideration)
## return: beta in the optimal weight of ETS

optim_beta <- function(lpd_feature, features_y = NULL) {
  y_lpd <- lpd_feature$lpd
  prob <- exp(y_lpd)
  if (sum(prob == 0) == 2) {
    cat("The probability predictive densities have 0 value")
  }
  prob <- prob[rowSums(prob == 0) != 2, ]

  ## maximizing TODO: change to a better optimization tool.
  ## library("optimx")
  w_max <- try(optim(
    par = 0,
    fn = logscore,
    features = features_y,
    prob = prob,
    ## intercept = intercept,
    intercept = TRUE,
    method = "L-BFGS-B",
    control = list(fnscale = -1)
  )
  )

  if (w_max$convergence != 0) {
    cat("The optimization does not converge in data", i_ts)
  }
  beta_optim <- w_max$par
  return(beta_optim)
}

cl <- makeCluster(2)
registerDoParallel(cl)
optimal_beta <-
  foreach(i_ts = 1:length(data)) %dopar%
  optim_beta(lpd_feature_yearly[[i_ts]], features_y = NULL)
stopCluster(cl)


#----------------------------------------------------------#
## Function  forecast_results
## param: a ts data with $x, $xx
## param:optimal_beta
## param:forecast_h (Forecasting horizon) In yearly data, forecast_h=6.
## return: the ts data with $ff, $mase_err, $smape_err, $logscore

forecast_results<-function(data, forecast_h, optimal_beta){
  library(forecast)
  ## forecasting
  features_y_hat <- matrix(nrow = forecast_h, ncol = 42)
  y_hat_matrix <- matrix(ncol = forecast_h, nrow = 4)
  rownames(y_hat_matrix)<-c("Optimal pool","SA","ETS","ARIMA")
  lpds_all<- matrix(ncol = 1, nrow = 4)

  y <- data$x
  y01 = scale(y, center = TRUE, scale = TRUE)
  y_mean = attr(y01, "scaled:center")
  y_sd = attr(y01, "scaled:scale")
  y01 = as.numeric(y01)
  y_true = data$xx
  y01_true = as.numeric(scale(y_true, center = y_mean, scale = y_sd))
  y_new = y01
  y_new_simple = y01
  y_new_ets = y01
  y_new_arima = y01

  lpds = 0
  lpds_simple = 0
  lpds_ets = 0
  lpds_arima=0
  pred_densities = matrix(NA, forecast_h, 2)
  pred_densities_simple = matrix(NA, forecast_h, 2)
  pred_densities_single = matrix(NA, forecast_h, 2)
  for (t in 1:forecast_h)
  { ## NOTE: This part is a recursive process, could not be parallelized.  Calculate
    ## predictive features. Each individual model provide an one-step-ahead predictive y
    ## values, and recalculate features based on optimized pools. This will do the
    ## forecasting model multiple times (h), consider to simplify it.
    PI_level<-model_conf$PI_level
    ## ETS model
    ets_fit <- ets(y_new, model = model_conf$ets_model)
    ets_fore<-forecast(ets_fit, h = 1, level = PI_level)
    ets_fore_mean <- ets_fore$mean
    ets_fore_sd = (ets_fore$lower - ets_fore$mean)/qnorm(1 -PI_level/100)

    ## ARIMA model
    arima_fit <- auto.arima(y_new)
    arima_fore <- forecast(arima_fit, h = 1, level = PI_level)
    arima_fore_mean <- arima_fore$mean
    arima_fore_sd = (arima_fore$lower - arima_fore$mean)/qnorm(1 - PI_level/100)

    ## ETS model (SA)
    ets_fit_simple <- ets(y_new_simple, model = model_conf$ets_model)
    ets_fore_simple <-forecast(ets_fit_simple, h = 1, level = PI_level)
    ets_fore_mean_simple <- ets_fore_simple$mean
    ets_fore_sd_simple = (ets_fore_simple$lower - ets_fore_simple$mean)/qnorm(1 - PI_level/100)

    ## ARIMA model (SA)
    arima_fit_simple <- auto.arima(y_new_simple)
    arima_fore_simple <- forecast(arima_fit_simple, h = 1, level = PI_level)
    arima_fore_mean_simple <- arima_fore_simple$mean
    arima_fore_sd_simple = (arima_fore_simple$lower - arima_fore_simple$mean)/qnorm(1 - PI_level/100)

    y_pred_h_simple = mean(c(ets_fore_mean_simple, arima_fore_mean_simple))
    y_new_simple = c(y_new_simple, y_pred_h_simple)
    y_hat_matrix[2,t]<-y_pred_h_simple

    ## ETS model (single)
    ets_fit_single <- ets(y_new_ets, model = model_conf$ets_model)
    ets_fore_single <-forecast(ets_fit_single, h = 1, level = PI_level)
    ets_fore_mean_single <- ets_fore_single$mean
    ets_fore_sd_single = (ets_fore_single$lower - ets_fore_single$mean)/qnorm(1 - PI_level/100)
    y_new_ets = c(y_new_ets, ets_fore_mean_single)
    y_hat_matrix[3,t]<-ets_fore_mean_single

    ## ARIMA model (single)
    arima_fit_single <- auto.arima(y_new_arima)
    arima_fore_single <- forecast(arima_fit_single, h = 1, level = PI_level)
    arima_fore_mean_single <- arima_fore_single$mean
    arima_fore_sd_single = (arima_fore_single$lower - arima_fore_single$mean)/qnorm(1 - PI_level/100)
    y_new_arima = c(y_new_arima, arima_fore_mean_single)
    y_hat_matrix[4,t]<-arima_fore_mean_single

    ## Update features

    features_y = NULL
    intercept = TRUE
    if(!is.null(features_y))
    {
      myts <- list(list(x=ts(y_new, frequency = frequency)))
      myfeatures <- THA_features(myts)[[1]]$features
      myfeatures <- data.matrix(myfeatures)
      myfeatures_scaled = scale(myfeatures, center = features_y_mean, scale = features_y_sd)
      # null in myfeatures_scaled
      myfeatures_scaled[is.na( myfeatures_scaled)] <- 0
      ## features_y_hat[t, ] <- myfeatures_scaled
    } else
    {
      myfeatures_scaled = NULL
    }

    ## Update predictive weights
    if(intercept) myfeatures_scaled = cbind(1, myfeatures_scaled)
    exp_lin = exp(myfeatures_scaled %*% optimal_beta)
    w <- exp_lin/(1+rowSums(exp_lin)) # 1-by-(n-1)
    w_full = cbind(w, 1 - rowSums(w)) # 1-by-n
    #feature_w<-c(feature_w,w)

    ## The final pooled y
    y_pred_h = sum(cbind(ets_fore_mean, arima_fore_mean) * w_full)
    y_new = c(y_new, y_pred_h)
    y_hat_matrix[1,t]<-y_pred_h

### The predictive log score

    pred_densities[t, 1] <- dnorm(y01_true[t], mean = ets_fore_mean, sd = ets_fore_sd)
    pred_densities[t, 2] <- dnorm(y01_true[t], mean = arima_fore_mean, sd = arima_fore_sd)
    lpds0 = log(sum(pred_densities[t,] * w_full))
    lpds = lpds + lpds0

    # SA
    pred_densities_simple[t, 1] <- dnorm(y01_true[t], mean = ets_fore_mean_simple, sd = ets_fore_sd_simple)
    pred_densities_simple[t, 2] <- dnorm(y01_true[t], mean = arima_fore_mean_simple, sd = arima_fore_sd_simple)
    lpds_simple0 = log(mean(pred_densities_simple[t, ]))
    lpds_simple = lpds_simple + lpds_simple0

    # Single ets and arima
    pred_densities_single[t, 1] <- dnorm(y01_true[t], mean = ets_fore_mean_single, sd = ets_fore_sd_single)
    pred_densities_single[t, 2] <- dnorm(y01_true[t], mean = arima_fore_mean_single, sd = arima_fore_sd_single)
    lpds_ets0 = log(pred_densities_single[t, 1])
    lpds_ets = lpds_ets + lpds_ets0
    lpds_arima0 = log(pred_densities_single[t, 2])
    lpds_arima = lpds_arima + lpds_arima0
  }
  lpds_all<-rbind(lpds,lpds_simple, lpds_ets, lpds_arima)
  rownames(lpds_all)<-c("Optimal pool","SA","ETS","ARIMA")
  colnames(lpds_all)<-c("log score")
  data$logscore<- lpds_all

### mase smape
  ff<- y_sd * y_hat_matrix + y_mean
  data$ff<-ff
  insample <- as.numeric(data$x)
  frq<-stats::frequency(insample)
  outsample <- as.numeric(data$xx)
  masep <- mean(abs(utils::head(insample, -frq) - utils::tail(insample, -frq)))
  repoutsample <- matrix(rep(outsample, each = nrow(ff)), nrow = nrow(ff))
  smape_err <- 200 * abs(ff - repoutsample)/(abs(ff) + abs(repoutsample))
  mase_err <- abs(ff - repoutsample)/masep
  data$mase_err <- data.matrix(rowMeans(mase_err))
  data$smape_err <- data.matrix(rowMeans(smape_err))
  return(data)
}

cl <- makeCluster(2)
registerDoParallel(cl)
data_forecast <- foreach(i_ts = 1:length(data)) %dopar%
                   forecast_results(data[[i_ts]],forecast_h=6, optimal_beta[[i_ts]])
stopCluster(cl)

save(data_forecast, file="E:/time series/R code/feature-based-bayesian-model-averaging/data/data_forecast_yearly.RData")
#----------------------------------------------------------#
## Function  forecast_performance
## param:data with $mase_err, $smape_err, $logscore
## return: a matrix performance

load("E:/time series/R code/feature-based-bayesian-model-averaging/data/data_forecast_yearly.RData")
forecast_performance<-function(data){
  # log score
  logscore_all<-c()
  for (i_ts in 1:length(data)) {
    # delete the data with -Inf log score
    if(i_ts %in% c(131,311,550,609,686)){
      logscore_all<-logscore_all
    }else{
      logscore_all<-cbind(logscore_all,data[[i_ts]]$logscore)
    }
  }
  logscore_all<-data.matrix(rowSums(logscore_all))
  # mase
  mase_all<-c()
  for (i_ts in 1:length(data)) {
    mase_all<-cbind(mase_all,data[[i_ts]]$mase_err)
    }
  mase_all<-data.matrix(rowMeans(mase_all))
  # smape
  smape_all<-c()
  for (i_ts in 1:length(data)) {
    smape_all<-cbind(smape_all,data[[i_ts]]$smape_err)
  }
  smape_all<-data.matrix(rowMeans(smape_all))
  performance<-cbind(logscore_all,mase_all,smape_all)
  colnames(performance)<-c("Log score","Mase","Smape")
  return(performance)
}

forecast_performance(data_forecast)

##              Log score     Mase    Smape
## Optimal pool -27929.45 3.222006 14.17247
## SA           -28706.42 3.473032 14.56170
## ETS          -37023.01 4.050580 16.05583
## ARIMA        -73666.30 3.425456 15.08570
