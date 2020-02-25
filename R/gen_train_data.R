#! /usr/bin/env Rscript

## feature-based-Bayesian-forecasting-model-averaging
setwd("~/code/febama/")

library(tsfeatures)
## library(M4metaresults)
## library(M4comp2018)
load("data/M4.rda")
library(M4metalearning)
library(forecast)
library(tseries)

##data(m4 quarterly data)
M4_q<-M4[23001:47000]
set.seed(2019-12-31)
indices <- sample(length(M4_q))
M4_q1 <- M4_q[indices[1:1000]]

## FF: 1000 matrix of features
## PP: 1000 matrix of log predictive probability density
## M4_q1: add $ets_fore,$ari_fore, $sigma_ets, $sigma_ari
## save
FF<-list()
PP<-list()

###----------------------------------------------------------------------------
### Model Settings
frequency = 4

ets_model = "ANN" # simple exponential smoothing with additive errors
forecast_h = 8
history_burn = 16 # Let the model to start with at least this length of historical data.
PI_level = 90
intercept = FALSE # Do not include intercept in the features
###----------------------------------------------------------------------------
for (a in 1:1000)
{

    ## A single historical data
    y<-M4_q1[[a]]$x

    ## Calculate historical log predictive density
    log_pred_densities<-matrix(nrow = length(y) - h, ncol = 2)
    for (t in (history_burn + 1):(length(y) - h))
    {
        ## TODO: We may simplify this to assume the fit and forecast procedure is
        ## invariant with certain time period to save time.

        ## ETS model
        ets_fit <- ets(y[1:(t-h)], model = ets_model)
        ets_fore<-forecast(ets_fit, h = forecast_h, level = PI_level)
        ets_fore_mean <- ets_fore$mean
        ## forecast.ets does not directly provide predictive variance but we could infer from
        ## predict interval.  Remember that PI = pred_mean -/+ qnorm(PI_level)*pred_sd
        ets_fore_sd = (ets_fore$lower - ets_fore$mean)/qnorm(1 - PI_level/100)

        ## ARIMA model
        arima_fit <- auto.arima(y[1:(t - h)])
        arima_fore <- forecast(arima_fit, h = forecast_h, level = PI_level)
        arima_fore_mean <- arima_fore$mean
        arima_fore_sd = (ari_fore$lower - ari_fore$mean)/qnorm(1 - PI_level/100)

        ## To keep numeric stability, we calculate log P(y_pred)
        log_pred_densities[, 1] <- sum(dnorm(y[(t + 1):(t + h)], mean = ets_fore_mean,
                                             sd = ets_fore_sd, log = TRUE))
        log_pred_densities[, 2] <- sum(dnorm(y[(t + 1):(t + h)], mean = arima_fore_mean,
                                             sd = arima_fore_sd, log = TRUE))



    }
    PP[[a]]<-log_pred_densities

    ## Calculate historical features
    features_y <- matrix(nrow = length(y) - h - history_burn, ncol = 42)
    for (t in (history_burn + 1):(length(y) - h))
    {
        myts <-list(list(x=ts(y[1:(t - h)], frequency = frequency)))
        myfeatures <- THA_features(myts)[[1]]$features
        myfeatures<-data.matrix(myfeatures)
        features_y[t, ] <- myfeatures
    }
    FF[[a]]<-features_y

    ## TODO: Standardize features

    ## Time series plot of features
    ## par(mfrow = c(6, 7), mar = c(5, 0, 0, 0))
    ## for(i in 1:42)
    ## {
    ##     plot(features_y[, i], type = "l", col = "red", xlab = colnames(features_y)[i])
    ## }

}
save(PP, FF, M4_q1,file="data/historical_log_pred_features.RData")
