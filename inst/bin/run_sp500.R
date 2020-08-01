library("tsfeatures")
library("M4metalearning")
library("forecast")
library("tseries")
library("purrr")
library("ggplot2")
library("numDeriv")
library("mvtnorm")
library("base")
library("MASS")
library(rugarch)
library(stochvol)
library(parallel)
library(doParallel)
library(foreach)

## Load R functions,
setwd("E:/time series/R code/febama-master (2)/febama-master")
source("R/features.R")
source("R/models.R")
source("R/mcmc.R")
source("R/priors.R")
source("R/posterior.R")
source("R/logscore.R")
source("R/febama.R")
source("E:/time series/R code/febama_old/R/arrange/op.R")
source("E:/time series/R code/febama_old/R/arrange/logscore.R")

# Data
library(tidyquant)
aapl_prices  <- tq_get("AAPL", get = "stock.prices", from = "1986-01-01", to = "2020-01-08")
aapl_prices_daily <- tq_mutate(data = aapl_prices, select = close, 
                               mutate_fun = periodReturn,
                               period = "daily", type = "log")
n = 250
n_xx = 50
n_x = 200
ts = ts(tail(aapl_prices_daily$daily.returns, n_xx + n_x))
ts_x = ts(ts[1:n_x], frequency = 1, start = 1, end = n_x)
ts_xx = ts(tail(ts, n_xx), frequency = 1, start = n_x+1, end = n_x + n_xx )
daiLogReturn <- list( x = ts_x, xx = ts_xx)
ts.plot(daiLogReturn$x, daiLogReturn$xx, gpars=list(col=c("blue","red")))

###----------------------------------------------------------------------------
### Model config template
###----------------------------------------------------------------------------
num_models = 3

## Default model config template
model_conf_default = list(
  frequency = 1
  #, ets_model = "ANN" 
  , forecast_h = length(daiLogReturn$xx)
  , train_h = 1 
  , history_burn = 100 
  , PI_level = 90 
  , roll = NULL 
  , feature_window = 100 
  , features_used = rep(list(c("ARCH.LM","entropy", "arch_acf", "garch_acf", "arch_r2", "garch_r2", "unitroot_kpss")), num_models - 1)
  , fore_model = c("garch_fore", "tgarch_fore", "auto.arima_fore")
  , lpd_features_parl = list(par = T, ncores = 3)
 
  , varSelArgs = rep(list(list(cand = "2:end", init = "all-in")), num_models - 1)
  
  , priArgs = rep(list(list("beta" = list(type = "cond-mvnorm",
                                          mean = 0, covariance = "identity", shrinkage = 1),
                            "betaIdx" = list(type = "beta", alpha0 = 1, beta0 = 1))), num_models - 1)
  
  , algArgs = list(initOptim = TRUE, 
                   algName = "sgld", 
                   nIter = 1, 
                   "sgld" = list(max_batchSize = 108,
                                 nEpoch = 10,
                                 burninProp = 0.4, 
                                 stepsize = 0.1,
                                 gama = 0.55,
                                 a = 0.4,
                                 b = 10)
  )
  
)

model_conf_curr = model_conf_default
lpd_features0 = lpd_features_multi(data = daiLogReturn, model_conf = model_conf_curr)
lpd_features <- feature_clean(list(lpd_features0))

fe <- lpd_features[[1]]$feat
fm <- lpd_features[[1]]$feat_mean
fs <- lpd_features[[1]]$feat_sd
lpd_features[[1]]$feat<- fe[, unique(unlist(model_conf_curr$features))]
lpd_features[[1]]$feat_mean <- fm[unique(unlist(model_conf_curr$features))]
lpd_features[[1]]$feat_sd <- fs[unique(unlist(model_conf_curr$features))]

out_us = febama_mcmc(data = lpd_features[[1]], model_conf = model_conf_curr)  

model_conf_curr$forecast_h = 20
daiLogReturn$xx = daiLogReturn$xx[1:20]
res_us = forecast_feature_results_multi(ts = daiLogReturn, model_conf = model_conf_curr, 
                                 data = lpd_features[[1]], beta_out = out_us)


out_op = optim_beta(lpd_feature = lpd_features[[1]], features_y = NULL)
res_op = forecast_results_nofea(data = daiLogReturn, model_conf = model_conf_curr, 
                         optimal_beta = out_op)

