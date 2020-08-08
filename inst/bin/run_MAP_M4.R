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

## Load R functions,
setwd("E:/time series/git/febama")
source("R/features.R")
source("R/models.R")
source("R/mcmc.R")
source("R/priors.R")
source("R/posterior.R")
source("R/logscore.R")
source("R/forecast.R")
source("R/Compare/OP.R")

# model setting
num_models = 6
model_conf_default = list(
  frequency = 12
  , ets_model = "ANN" 
  , forecast_h = 18 
  , train_h = 1 
  , history_burn = 36 
  , PI_level = 90 
  , roll = NULL 
  , feature_window = NULL 
  , features_used = rep(list(c("entropy", "arch_acf", "alpha", "beta", "unitroot_kpss")), num_models - 1)
  , fore_model = c("ets_fore",  "naive_fore", "rw_drift_fore","auto.arima_fore",
                   "snaive_fore", "stlm_ar_fore")
  , lpd_features_parl = list(par = F, ncores = 1)
  , varSelArgs = rep(list(list(cand = "2:end", init = "all-in")), num_models - 1)
  
  , priArgs = rep(list(list("beta" = list(type = "cond-mvnorm",
                                          mean = 0, covariance = "identity", shrinkage = 10),
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

# data
library(M4comp2018)
set.seed(2020-0716)
data_test <- M4[sample(c(47001:95000), 10)]

# Do multi-step forecasting 
model_conf_curr = model_conf_default
# Step 1: Compute log predictive density and features
library(foreach)
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
lpd_features0 <- foreach(i_ts = 1:10, .packages = "M4metalearning", .export = model_conf_curr$fore_model) %dopar%
  lpd_features_multi(data = data_test[[i_ts]], model_conf = model_conf_curr)
stopCluster(cl)

lpd_features <- feature_clean(lpd_features0)

for (i in 1:length(lpd_features)) {
  fe <- lpd_features[[i]]$feat
  fm <- lpd_features[[i]]$feat_mean
  fs <- lpd_features[[i]]$feat_sd
  lpd_features[[i]]$feat<- fe[, unique(unlist(model_conf_curr$features))]
  lpd_features[[i]]$feat_mean <- fm[unique(unlist(model_conf_curr$features))]
  lpd_features[[i]]$feat_sd <- fs[unique(unlist(model_conf_curr$features))]
}

# Use MAP to obtain the optimal beta
cl <- makeCluster(3)
registerDoParallel(cl)
out_us <- foreach(i_ts = 1:10, .packages = "M4metalearning", .export = model_conf_curr$fore_model) %dopar%
  febama_mcmc(data = lpd_features[[i_ts]], model_conf = model_conf_curr)
stopCluster(cl)

# Forecast with time-varying weights
cl <- makeCluster(3)
registerDoParallel(cl)
res_us <- foreach(i_ts = 1:10, .packages = "M4metalearning", .export = model_conf_curr$fore_model) %dopar%
  forecast_feature_results_multi(ts = data_test[[i_ts]], model_conf = model_conf_curr, 
                                 data = lpd_features[[i_ts]], beta_out = out_us[[i_ts]])
stopCluster(cl)

# Print average forecasting performance
res_us_comb=c()
for (i in 1:length(res_us)) {
  res_us_comb = rbind(res_us_comb, res_us[[i]]$err_feature)
}
print("The logscore, mase and smape of our method(with features) are")
colMeans(res_us_comb)

# Monthly
lpd_features[[40]]$lpd[240,] = rep(0,6)
lpd_features[[261]]$lpd[141,] = rep(0,6)
lpd_features[[581]]$lpd[1:52,] = rep(0,6)
lpd_features[[676]]$lpd[116,] = rep(0,6)
#Comepare with OP and SA
cl <- makeCluster(3)
registerDoParallel(cl)
out_op <- foreach(i_ts = 1:10, .packages = "M4metalearning", .export = model_conf_curr$fore_model) %dopar%
  optim_beta(lpd_feature = lpd_features[[i_ts]], features_y = NULL)
stopCluster(cl)

cl <- makeCluster(3)
registerDoParallel(cl)
res_op <- foreach(i_ts = 1:10, .packages = "M4metalearning", .export = model_conf_curr$fore_model) %dopar%
  forecast_results_nofea(data = data_test[[i_ts]], model_conf = model_conf_curr, 
                         optimal_beta = out_op[[i_ts]])
stopCluster(cl)

res_ls=c()
res_mase=c()
res_smape=c()
for (i in 1:length(res_op)) {
  res_ls = rbind(res_ls, res_op[[i]]$logscore)
  res_mase = rbind(res_mase, t(res_op[[i]]$mase_err))
  res_smape = rbind(res_smape, t(res_op[[i]]$smape_err))
  
}
print("The logscore, mase and smape of different methods are")
colMeans(res_ls)
colMeans(res_mase)
colMeans(res_smape)
