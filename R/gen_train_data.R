#测试
#! /usr/bin/env Rscript

## feature-based-Bayesian-forecasting-model-averaging
#setwd("~/code/febama/")

library(tsfeatures)
library(M4metalearning)
library(forecast)
library(tseries)
library(purrr)
library(ggplot2)
library(M4comp2018)
library(foreach)
library(doParallel)

#load("data/M4.rda")

##data(m4 yearly data)
M4_y<-M4[1:23000]

# ##data(m4 quarterly data)
# M4_q<-M4[23001:47000]

set.seed(2020-3-3)
indices <- sample(length(M4_y))
data <- M4_y[indices[1:1000]]
#data<-M4_y


###----------------------------------------------------------------------------
### Model Settings
model_conf = list(
  frequency = 1
  , ets_model = "ANN" # simple exponential smoothing with additive errors
  , forecast_h = 6 # Forecasting horizon
  , train_h = 1 # Out of sample training, 1 means rolling forecast
  , history_burn = 6 # Let the model to start with at least this length of historical data.
  , PI_level = 90 # Predictive Interval level, used to extract out-of-sample variance from forecasting models.
  , intercept = FALSE # Should intercept in the features
)
###----------------------------------------------------------------------------

# ##并行案例
# fun <- function(x){
#  y<-x+1
# }
# # 非并行计算方式，类似于sapply函数的功能
# x <- foreach(x=1:1000,.combine='rbind') %do% fun(x)
# # 启用parallel作为foreach并行计算的后端
# cl <- makeCluster(4)
# registerDoParallel(cl)
# # 并行计算方式
# x <- foreach(x=1:1000,.combine='rbind') %dopar% fun(x)
# stopCluster(cl)


 #attach(model_conf) 

## Function  lpd_feature
## param: a ts data with $x
## return: lpd_feature (a list including two matrices of probability densities and features)
lpd_feature <- function(data) {
  library(tsfeatures)
  library(M4metalearning)
  library(forecast)
  library(tseries)
  library(purrr)
  library(ggplot2)
  library(M4comp2018)
  library(foreach)
  library(doParallel)
  
  frequency = model_conf$frequency
  history_burn = model_conf$history_burn
  ets_model = model_conf$ets_model
  forecast_h = model_conf$forecast_h
  train_h = model_conf$train_h
  PI_level = model_conf$PI_level
  
  ## A single historical data
  y <- data$x
  y01 = scale(y, center = TRUE, scale = TRUE)
  y_mean = attr(y01, "scaled:center")
  y_sd = attr(y01, "scaled:scale")
  
  y01 = as.numeric(y01)
  
  ## Calculate historical log predictive density
  log_pred_densities <-
    matrix(nrow = length(y) - history_burn - train_h + 1,
           ncol = 2)
  
  for (t in (history_burn):(length(y) - train_h))
  {
    ## TODO: We may simplify this to assume the fit and forecast procedure is
    ## invariant with certain time period to save time.
    
    ## ETS model
    ets_fit <- ets(y01[1:t], model = ets_model)
    ets_fore <- forecast(ets_fit, h = train_h, level = PI_level)
    ets_fore_mean <- ets_fore$mean
    ## forecast.ets does not directly provide predictive variance but we could infer from
    ## predict interval.  Remember that PI = pred_mean -/+ qnorm(PI_level)*pred_sd
    ets_fore_sd = (ets_fore$lower - ets_fore$mean) / qnorm(1 - PI_level /
                                                             100)
    
    ## ARIMA model
    arima_fit <- auto.arima(y01[1:t])
    arima_fore <- forecast(arima_fit, h = train_h, level = PI_level)
    arima_fore_mean <- arima_fore$mean
    arima_fore_sd = (arima_fore$lower - arima_fore$mean) / qnorm(1 - PI_level /
                                                                   100)
    
    ## To keep numeric stability, we calculate log P(y_pred)
    log_pred_densities[(t - history_burn + 1), 1] <-
      sum(dnorm(
        y01[(t + 1):(t + train_h)],
        mean = ets_fore_mean,
        sd = ets_fore_sd,
        log = TRUE
      ))
    log_pred_densities[(t - history_burn + 1), 2] <-
      sum(dnorm(
        y01[(t + 1):(t + train_h)],
        mean = arima_fore_mean,
        sd = arima_fore_sd,
        log = TRUE
      ))
  }
  
  ## Calculate historical features
  features_y <-
    matrix(nrow = length(y) - train_h - history_burn + 1,
           ncol = 42)
  
  for (t in (history_burn):(length(y) - train_h))
  {
    myts <- list(list(x = ts(y[1:t], frequency = frequency)))
    myfeatures <- THA_features(myts)[[1]]$features
    myfeatures <- data.matrix(myfeatures)
    features_y[(t - history_burn + 1),] <- myfeatures
  }
  features_y_scaled = scale(features_y, center = TRUE, scale = TRUE)
  
  lpd_feature <-
    list(lpd = log_pred_densities, feat = features_y_scaled)
  return(lpd_feature)
}

 
cl <- makeCluster(2)
registerDoParallel(cl)
lpd_feature_yearly <- foreach(i_ts = 1:length(data)) %dopar% lpd_feature(data[[i_ts]])
stopCluster(cl)

save(lpd_feature_yearly, data, model_conf, file="E:/time series/R code/feature-based-bayesian-model-averaging/data/historical_lpd_feature_yearly.RData")


#--------------------------------------------------------------#
## Calculate the feature with a moving window
y<-M4[[13190]]$x
train_h = 1
history_burn = 6
window = 100
features_window<-function( y, window, history_burn, train_h){
  features_y<-c()
  for (t in (history_burn):(length(y) - train_h))
  {
    if(t <= window){
      myts <-list(list(x=ts(y[1:t], frequency = 1)))
      myfeatures <- THA_features(myts)[[1]]$features
      myfeatures <- data.matrix(myfeatures)
    }else{
      myts <-list(list(x=ts(y[(t-window+1):t], frequency = 1)))
      myfeatures <- THA_features(myts)[[1]]$features
      myfeatures <- data.matrix(myfeatures)
    }
    features_y<-rbind(features_y,myfeatures)
  }
  return(features_y)
}
t1=proc.time()
features_y<-features_window( y, window, history_burn, train_h)
t2=proc.time()
t=t2-t1
print(paste0('执行时间：',t[3][[1]],'秒'))

## Time series plot of features
par(mfrow = c(6, 7), mar = c(5, 0, 0, 0))
for(i in 1:42)
{
  plot(features_y[, i], type = "l", col = "blue", xlab = colnames(features_y)[i])
}
