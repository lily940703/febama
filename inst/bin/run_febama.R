#! /usr/bin/env Rscript

## path = "/gs/home/kangyf/lily/febama"
## setwd(path)
rm(list = ls())
setwd("~/code/febama")

## Load R functions,
## library("febama") # TODO: use a package
source("R/features.R")
source("R/models.R")
source("R/mcmc.R")
source("R/priors.R")
source("R/posterior.R")
source("R/logscore.R")
source("R/febama.R")

# Should recalculate the features and save to path, or load from the saved path.
lpd_features_loc = list("calculate" = FALSE,
                        save_path = "data/lpd_features_yearly.Rdata")

num_models = 3
###----------------------------------------------------------------------------
### Model config template
###----------------------------------------------------------------------------

## Default model config template
model_conf_default = list(
    frequency = 4
  , ets_model = "ANN" # simple exponential smoothing with additive errors
  , forecast_h = 8 # Forecasting horizon
  , train_h = 1 # Out of sample training, 1 means rolling forecast
  , history_burn = 16 # Let the model to start with at least this length of historical data.
  , PI_level = 90 # Predictive Interval level, used to extract out-of-sample variance from forecasting models.
  , roll = NULL # The length of rolling samples, larger than history_burn
  , feature_window = NULL # The length of moving window when computing features
  , features_used = rep(list(c("entropy", "arch_acf", "alpha", "beta", "unitroot_kpss")), num_models - 1)
  , fore_model = c("ets_fore",  "naive_fore", "rw_drift_fore")

    ## Variable selection settings. By default, every model shares the same
    ## settings. Otherwise, write the full list, same applies to priArgs, this would allow
    ## for some models with only intercept.  Variable selection candidates, NULL: no
    ## variable selection use the full covariates provided by $init. ("all-in", "all-out",
    ## "random", or user-input)
  , varSelArgs = rep(list(list(cand = "2:end", init = "all-in")), num_models - 1)

  , priArgs = rep(list(list("beta" = list(type = "cond-mvnorm",
                                          mean = 0, covariance = "identity", shrinkage = 1),
                            "betaIdx" = list(type = "beta", alpha0 = 1, beta0 = 1))), num_models - 1)
  , algArgs = list(initOptim = TRUE, # Use LBFGS to optimize initial values
                   algName = "sgld", # could be NA, results are only based on optimization.
                   nIter = 20, # number of iterations
                   "sgld" = list(max_batchSize = 16,
                                 nEpoch = 3,
                                 burninProp = 0.4, # burnin proportion within each SGLD.
                                 stepsize = NA,
                                 gama = 0.55,
                                 a = 0.4,
                                 b = 10)
                   )

)


## Model without variable selection
model_conf_NoVS = model_conf_default
model_conf_NoVS[["varSelArgs"]] = rep(list(list(cand = NULL, init = "all-in")), num_models - 1)

## Model with only intercept (Bayesian optimal pool)
model_conf_NoFeat = model_conf_default
model_conf_NoFeat[["features_used"]] = rep(list(NULL), num_models - 1)
model_conf_NoFeat[["varSelArgs"]] = rep(list(list(cand = NULL, init = "all-in")), num_models - 1)


###----------------------------------------------------------------------------
### Experiments
###----------------------------------------------------------------------------

## Load current model config
model_conf_curr = model_conf_default
## model_conf_curr = model_conf_NoVS
## model_conf_curr = model_conf_NoFeat
## -------------------------  Experiment  ----------------------------#
## library(foreach)
## library(doParallel)
## cl <- makeCluster(parallel::detectCores())
## registerDoParallel(cl)

## clusterEvalQ(cl,{
    library("tsfeatures")
    library("M4metalearning")
    library("forecast")
    library("tseries")
    library("purrr")
    library("ggplot2")
    # library(M4comp2018)
    library("numDeriv")
    library("mvtnorm")
    library("base")
    library("MASS")
## })

## clusterExport(cl, model_conf_curr$fore_model)

if(lpd_features_loc$calculate == TRUE)
{
    ## Load data
    ## library(M4comp2018)
    message("Loading the time series data.")
    load("data/M4.rda")

    ## set.seed(2020-0503)
    data_test <- M4[sample(c(23001:47000), 10)]

    ## Extract `all 42 features` and given models (model_conf_curr$fore_model)
    message("Extracting LPD and features for given models: ",
            paste(model_conf_curr$fore_model, collapse = ", "))

    lpd_features0 <- lapply(data_test, lpd_features_multi, model_conf=model_conf_curr)
    lpd_features <- feature_clean(lpd_features0)

    save(lpd_features, file = lpd_features_loc$save_path)
    message("LPD and features are saved to: ", lpd_features_loc$save_path)

} else {
    load(lpd_features_loc$save_path)
    message("LPD and features are loaded from: ", lpd_features_loc$save_path)
}

## Extract lpd and candidate features from `model_conf_curr$features`
for (i in 1:length(lpd_features)) {
    fe <- lpd_features[[i]]$feat
    fm <- lpd_features[[i]]$feat_mean
    fs <- lpd_features[[i]]$feat_sd
    lpd_features[[i]]$feat<- fe[, unique(unlist(model_conf_curr$features))]
    lpd_features[[i]]$feat_mean <- fm[unique(unlist(model_conf_curr$features))]
    lpd_features[[i]]$feat_sd <- fs[unique(unlist(model_conf_curr$features))]
}

## Algorithm
OUT = lapply(lpd_features, febama_mcmc, model_conf = model_conf_curr)


## beta_pre <- foreach(i_ts = 1:length(SGLD_VS)) %dopar%
##     beta_prepare(SGLD_VS[[i_ts]])

## fore_feat <- foreach(i_ts = 1:length(data_test)) %dopar%
##     forecast_feature_results_multi(data = data_test[[i_ts]], model_conf_curr = model_conf_curr,
##                                    intercept = T, lpd_feature = lpd_features[[i_ts]],
##                                    beta_pre = beta_pre[[i_ts]])
## perform_feat <- forecast_feature_performance(fore_feat)

## ## Compare
## optim <- foreach(i_ts = 1:length(lpd_feature)) %dopar%
##     optim_beta (lpd_features[[i_ts]], features_y = NULL)

## fore <- foreach(i_ts = 1:length(optim)) %dopar%
##     forecast_results_nofea(data = data_test[[i_ts]],
##                            model_conf_curr = model_conf_curr, optimal_beta = optim[[i_ts]])

## perform <- forecast_performance(fore)

## stopCluster(cl)

## save(data_test, model_conf, lpd_feature, SGLD_VS,
##      beta_pre, fore_feat, perform_feat,
##      optim, fore, perform, file = "test/Q1000.RData")

## ## Visualization

## par(mfrow = c(4,2))
## for (i in 1:4) {
##     plot(SGLD_VS[[1]]$result_all[[i]]$logscore[1,], ylab = "log score",
##          xlab = paste0(i,'th iteration of beta1' ),ylim = c(-18,-15))
##     abline(h=-15.21724, col = "2")

##     plot(SGLD_VS[[1]]$result_all[[i]]$logscore[2,], ylab = "log score",
##          xlab = paste0(i,'th iteration of beta2' ), ylim = c(-18,-15))
##     abline(h=-15.21724, col = "2")
## }

## par(mfrow = c(4,2))
## for (i in 1:4) {
##     plot(SGLD_VS[[5]]$result_all[[i]]$logscore[1,], ylab = "log score",
##          xlab = paste0(i,'th iteration of beta1' ))
##     abline(h=63.32389, col = "2")
##     plot(SGLD_VS[[5]]$result_all[[i]]$logscore[2,], ylab = "log score",
##          xlab = paste0(i,'th iteration of beta2' ))
##     abline(h=63.32389, col = "2")
## }


## par(mfrow = c(4,2))
## for (i in 91:94) {
##     plot(SGLD_VS[[5]]$result_all[[i]]$logscore[1,], ylab = "log score",
##          xlab = paste0(i,'th iteration of beta1' ))
##     abline(h=63.32389, col = "2")
##     plot(SGLD_VS[[5]]$result_all[[i]]$logscore[2,], ylab = "log score",
##          xlab = paste0(i,'th iteration of beta2' ))
##     abline(h=63.32389, col = "2")
## }

## ll <- c()
## for (i in 1:1000) {
##     l <- length(data_test[[i]]$x)
##     ll <- c(ll,l)
## }
## par(mfrow = c(1,1))
## hist(ll, main = "Histogram of length in historical data", ylim = c(0,500))

## rank_ls <- c()
## for (i in 1:1000) {
##     ls <- c(fore_feat[[i]]$err_feature[,"lpds"], fore[[i]]$logscore)
##     names(ls) <- c("Ours", "OP","SA","ets","naive","rw_drift")
##                                         # ascending order
##     ran <- rank(ls)
##     rank_ls <- rbind(rank_ls, ran)
## }

## ls_win <- apply(rank_ls, 2, function(x){sum(x==6)})
## ls_win <- cbind(names(ls_win),data.frame(ls_win), rep("logscore",6))
## colnames(ls_win)=NULL
## rownames(ls_win)=NULL


## rank_mase <- c()
## for (i in 1:1000) {
##     ls <- c(fore_feat[[i]]$err_feature[,"mase_err_h"], fore[[i]]$mase_err)
##     names(ls) <- c("Ours", "OP","SA","ets","naive","rw_drift")
##     ran <- rank(ls, ties.method = "random")
##     rank_mase <- rbind(rank_mase, ran)
## }

## mase_win <- apply(rank_mase, 2, function(x){sum(x==1)})
## mase_win <- cbind(names(mase_win),data.frame(mase_win), rep("mase",6))
## colnames(mase_win)=NULL
## rownames(mase_win)=NULL


## rank_smape <- c()
## for (i in 1:1000) {
##     ls <- c(fore_feat[[i]]$err_feature[,"smape_err_h"], fore[[i]]$smape_err)
##     names(ls) <- c("Ours", "OP","SA","ets","naive","rw_drift")
##     ran <- rank(ls, ties.method = "random")
##     rank_smape <- rbind(rank_smape, ran)
## }

## smape_win <- apply(rank_smape, 2, function(x){sum(x==1)})
## smape_win <- cbind(names(smape_win),data.frame(smape_win), rep("smape",6))
## colnames(smape_win)=NULL
## rownames(smape_win)=NULL


## win <- rbind(data.frame(ls_win), data.frame(mase_win), data.frame(smape_win))

## Q1000_win <- data.frame(method = win[,1], win_num = win[,2], error = win[,3])

## library("ggplot2")
## ggplot(Q1000_win, aes(x=error, y=win_num, fill=method)) + geom_bar(stat='identity')

## ind <- Q1000_win$method %in% c("ets","naive","rw_drift")
## Q1000_win0 <- Q1000_win[!ind,]

## ggplot(Q1000_win0,aes(x=error, y=win_num, fill=method)) + geom_bar(stat='identity')
