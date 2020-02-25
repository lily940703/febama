#! /usr/bin/env Rscript

## feature-based-Bayesian-forecasting-model-averaging
setwd("~/code/febama/")

library(tsfeatures)
library(M4metalearning)
library(forecast)
library(tseries)

## library(M4metaresults)

## library(M4comp2018)
load("data/M4.rda")

##data(m4 quarterly data)
M4_q<-M4[23001:47000]

indices <- sample(length(M4_q))

data <- M4_q[indices[1:2]]

## FF: 1000 matrix of features
## PP: 1000 matrix of log predictive probability density
## M4_q1: add $ets_fore,$arima_fore, $sigma_ets, $sigma_ari
## save
feat <- list()
lpred_dens <- list()

###----------------------------------------------------------------------------
### Model Settings
model_conf = list(
  frequency = 4
  , ets_model = "ANN" # simple exponential smoothing with additive errors
  , forecast_h = 8 # Forecasting horizon
  , train_h = 1 # Out of sample training, 1 means rolling forecast
  , history_burn = 16 # Let the model to start with at least this length of historical data.
  , PI_level = 90 # Predictive Interval level, used to extract out-of-sample variance from
              # forecasting models.
  , intercept = FALSE # Should intercept in the features
)
###----------------------------------------------------------------------------

attach(model_conf)
for (i_ts in 1:length(data))
{
  ## TODO: parallel this part

  cat("Processing ", i_ts, "\n")

  ## A single historical data
  y <- data[[i_ts]]$x

  ## Calculate historical log predictive density
  log_pred_densities<-matrix(nrow = length(y) - history_burn - train_h, ncol = 2)
  for (t in (history_burn + 1):(length(y) - train_h))
  {
    ## TODO: We may simplify this to assume the fit and forecast procedure is
    ## invariant with certain time period to save time.

    ## ETS model
    ets_fit <- ets(y[1:(t-train_h)], model = ets_model)
    ets_fore<-forecast(ets_fit, h = train_h, level = PI_level)
    ets_fore_mean <- ets_fore$mean
    ## forecast.ets does not directly provide predictive variance but we could infer from
    ## predict interval.  Remember that PI = pred_mean -/+ qnorm(PI_level)*pred_sd
    ets_fore_sd = (ets_fore$lower - ets_fore$mean)/qnorm(1 - PI_level/100)

    ## ARIMA model
    arima_fit <- auto.arima(y[1:(t - train_h)])
    arima_fore <- forecast(arima_fit, h = train_h, level = PI_level)
    arima_fore_mean <- arima_fore$mean
    arima_fore_sd = (arima_fore$lower - arima_fore$mean)/qnorm(1 - PI_level/100)

    ## To keep numeric stability, we calculate log P(y_pred)
    log_pred_densities[(t - history_burn), 1] <- sum(dnorm(y[(t + 1):(t + train_h)], mean = ets_fore_mean,
                                                           sd = ets_fore_sd, log = TRUE))
    log_pred_densities[(t - history_burn), 2] <- sum(dnorm(y[(t + 1):(t + train_h)], mean = arima_fore_mean,
                                                           sd = arima_fore_sd, log = TRUE))
  }
  lpred_dens[[i_ts]]<-log_pred_densities

  ## Calculate historical features
  features_y <- matrix(nrow = length(y) - train_h - history_burn, ncol = 42)
  for (t in (history_burn + 1):(length(y) - train_h))
  {
    myts <-list(list(x=ts(y[1:(t - train_h)], frequency = frequency)))
    myfeatures <- THA_features(myts)[[1]]$features
    myfeatures <- data.matrix(myfeatures)
    features_y[(t - history_burn), ] <- myfeatures
  }

  features_y_scaled = scale(features_y, center = TRUE, scale = TRUE)

  feat[[i_ts]]<-features_y_scaled

  ## Time series plot of features
  ## par(mfrow = c(6, 7), mar = c(5, 0, 0, 0))
  ## for(i in 1:42)
  ## {
  ##     plot(features_y[, i], type = "l", col = "red", xlab = colnames(features_y)[i])
  ## }

}
save(feat, lpred_dens, data, model_conf, file="data/historical_log_pred_features.RData")
