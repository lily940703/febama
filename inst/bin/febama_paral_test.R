#! /usr/bin/env Rscript

## path = "/gs/home/kangyf/lily/febama"
## setwd(path)
setwd("~/code/febama")

load("data/M4.rda")

## library(M4comp2018)
set.seed(2020-0503)
data_test <- M4[sample(c(23001:47000), 1000)]
## load("lpd_feature_Y500.RData")

model_conf = list(
    frequency = 4
  , ets_model = "ANN" # simple exponential smoothing with additive errors
  , forecast_h = 8 # Forecasting horizon
  , train_h = 1 # Out of sample training, 1 means rolling forecast
  , history_burn = 16 # Let the model to start with at least this length of historical data.
  , PI_level = 90 # Predictive Interval level, used to extract out-of-sample variance from forecasting models.
  , roll = NULL # The length of rolling samples, larger than history_burn
  , feature_window =NULL # The length of moving window when computing features
    ## , fore_model = list("ets_fore", "auto.arima_fore" )
  , fore_model = list("ets_fore",  "naive_fore", "rw_drift_fore")
)



forecast_feature_results_multi <-function(data, model_conf, intercept = T,
                                          lpd_feature, beta_pre)
{
    ##attach(model_conf)
    feature_window = model_conf$feature_window
    roll = model_conf$roll
    frequency = model_conf$frequency
    history_burn = model_conf$history_burn
    ets_model = model_conf$ets_model
    forecast_h = model_conf$forecast_h
    train_h = model_conf$train_h
    PI_level = model_conf$PI_level
    fore_model = model_conf$fore_model


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
    pred_densities = matrix(NA, forecast_h, length(fore_model))
    colnames(pred_densities) <- unlist(fore_model)
    w_full_mean_h <- c()

    w_full_all <- list()
    for (t in 1:forecast_h)
    {
        ## Update features
        if(!is.null(feature_window)){
            y_new_nonsd1 <- tail(y_new_nonsd, feature_window)
        }else{
            y_new_nonsd1 <- y_new_nonsd
        }

        if (!is.null(features_y))
        {
            ##用非标准化的数据计算特征
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
            exp_lin = exp(matrix(myfeatures_scaled, nrow = 1) %*% beta_pre$beta)
            ##  avoid Inf
            exp_lin[exp_lin > exp(709)] <- exp(709)
            w <- exp_lin / (1 + rowSums(exp_lin))
            w_full = cbind(w, 1 - rowSums(w))
            return(w_full)
        }

        w_full <- sapply(beta_pre, w_get, myfeatures_scaled = myfeatures_scaled)
        w_full_all[[t]] <- w_full
        w_full_mean <- rowMeans(w_full)
        w_full_mean_h <- cbind(w_full_mean_h, w_full_mean)

        ## forecast
        if(is.null(roll)){
            y_new1 <- y_new
        }else{
            y_new1 <- tail(y_new, roll)
        }

        multi_fore <- lapply(fore_model, function(method){
            method_fun <- get(method)
            mean_sd <- method_fun (y_new1, train_h, PI_level)
            return(mean_sd)
        })

        y_pred_multi <- sum (w_full_mean * (sapply(multi_fore, function(mean_sd){
            return(mean_sd[[1]])
        })))

        y_new = c(y_new, y_pred_multi)
        y_new_nonsd = c(y_new_nonsd, (y_pred_multi * y_sd + y_mean))
        y_hat_matrix[1, t] <- y_pred_multi * y_sd + y_mean

        ## The predictive log score
        pd_multi <- sapply(multi_fore, function(mean_sd){
            dnorm(y01_true[t], mean = mean_sd[[1]], sd = mean_sd[[2]], log = F)
        })
        lpds_multi = log(sum(pd_multi * w_full_mean))

        if(lpds_multi < log(1e-323)){
            lpds_multi = log(1e-323)
        }else{
            lpds_multi = lpds_multi
        }
        lpds = lpds + lpds_multi
    }

    data$ff_feature <- y_hat_matrix
    colnames(w_full_mean_h) <- seq(1, forecast_h, 1)
    rownames(w_full_mean_h) <- unlist(fore_model)
    data$w_time_varying <- w_full_mean_h
    data$w_detail <- w_full_all

    ## mase smape
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

forecast_feature_performance<-function(data)
{
    ## log score
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


##--------------------------------------------------------------------------------------#
optim_beta <- function(lpd_feature, features_y = NULL) {
    y_lpd <- lpd_feature$lpd
    prob <- exp(y_lpd)
    prob[prob == 0] <- 1e-323
    num_model <- length(lpd_feature$lpd[1,])
    ini <-  t(data.matrix(rep(0, num_model-1)))

    w_max <- try(optim(
        par = ini,
        fn = log_score,
        ## gr = gradient,
        features = features_y,
        prob = prob,
        ## intercept = intercept,
        intercept = TRUE,
        method = "L-BFGS-B",
        ## method = "BFGS",
        control = list(fnscale = -1)
    )
    )

    if (w_max$convergence != 0) {
        w_max <- try(optim(
            par = ini,
            fn = log_score,
            ## gr = gradient,
            features = features_y,
            prob = prob,
            ## intercept = intercept,
            intercept = TRUE,
            ## method = "L-BFGS-B",
            method = "BFGS",
            control = list(fnscale = -1)
        )
        )
    }
    beta_optim <- w_max$par
    return(list(beta_optim = beta_optim, logscore = w_max$value))
}

## -------------------------  Experiment  ----------------------------#
library(foreach)
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

clusterEvalQ(cl,{
    library(tsfeatures)
    library(M4metalearning)
    library(forecast)
    library(tseries)
    library(purrr)
    library(ggplot2)
    library(M4comp2018)
    library(numDeriv)
    library("mvtnorm")
    library("base")
    library(MASS)})

clusterExport(cl,c("ets_fore", "naive_fore", "rw_drift_fore"))

lpd_feature_Q1000 <- foreach(i_ts = 1:length(data_test)) %dopar%
    lpd_feature_multi(data_test[[i_ts]], model_conf)
lpd_feature_Q1000 <- feature_clean(lpd_feature_Q1000)
                                        # five features
for (i in 1:length(lpd_feature_Q1000)) {
    fe <- lpd_feature_Q1000[[i]]$feat
    fm <- lpd_feature_Q1000[[i]]$feat_mean
    fs <- lpd_feature_Q1000[[i]]$feat_sd
    lpd_feature_Q1000[[i]]$feat<- fe[,colnames(fe) %in% c("entropy", "arch_acf", "alpha", "beta", "unitroot_kpss")]
    lpd_feature_Q1000[[i]]$feat_mean <- fm[names(fm) %in% c("entropy", "arch_acf", "alpha", "beta", "unitroot_kpss")]
    lpd_feature_Q1000[[i]]$feat_sd <- fs[names(fs) %in% c("entropy", "arch_acf", "alpha", "beta", "unitroot_kpss")]
}


SGLD_VS_Q1000 <- foreach(i_ts = 1:length(lpd_feature_Q1000)) %dopar%
    SGLD_VS (data = lpd_feature_Q1000[[i_ts]], logLik = log_score,
             gradient_logLik = gradient_logscore, prior = prior, stepsize = 0.1,
             SGLD_iter = 500, SGLD_iter_noVS = 50, VS_iter = 100,
             minibatchSize = NULL, sig = 10)

beta_pre_Q1000 <- foreach(i_ts = 1:length(SGLD_VS_Q1000)) %dopar%
    beta_prepare(SGLD_VS_Q1000[[i_ts]])

fore_feat_Q1000 <- foreach(i_ts = 1:length(data_test)) %dopar%
    forecast_feature_results_multi(data = data_test[[i_ts]], model_conf = model_conf,
                                   intercept = T, lpd_feature = lpd_feature_Q1000[[i_ts]],
                                   beta_pre = beta_pre_Q1000[[i_ts]])
perform_feat_Q1000 <- forecast_feature_performance(fore_feat_Q1000)

## Compare
optim_Q1000 <- foreach(i_ts = 1:length(lpd_feature_Q1000)) %dopar%
    optim_beta (lpd_feature_Q1000[[i_ts]], features_y = NULL)

fore_Q1000 <- foreach(i_ts = 1:length(optim_Q1000)) %dopar%
    forecast_results_nofea(data = data_test[[i_ts]],
                           model_conf = model_conf, optimal_beta = optim_Q1000[[i_ts]])

perform_Q1000 <- forecast_performance(fore_Q1000)

stopCluster(cl)

save(data_test, model_conf, lpd_feature_Q1000, SGLD_VS_Q1000,
     beta_pre_Q1000, fore_feat_Q1000, perform_feat_Q1000,
     optim_Q1000, fore_Q1000, perform_Q1000, file = "E:/time series/并行代码/test/Q1000.RData")

##-------------------------------------------------------------------------------------#

par(mfrow = c(4,2))
for (i in 1:4) {
    plot(SGLD_VS_Q1000[[1]]$result_all[[i]]$logscore[1,], ylab = "log score",
         xlab = paste0(i,'th iteration of beta1' ),ylim = c(-18,-15))
    abline(h=-15.21724, col = "2")

    plot(SGLD_VS_Q1000[[1]]$result_all[[i]]$logscore[2,], ylab = "log score",
         xlab = paste0(i,'th iteration of beta2' ), ylim = c(-18,-15))
    abline(h=-15.21724, col = "2")
}

par(mfrow = c(4,2))
for (i in 1:4) {
    plot(SGLD_VS_Q1000[[5]]$result_all[[i]]$logscore[1,], ylab = "log score",
         xlab = paste0(i,'th iteration of beta1' ))
    abline(h=63.32389, col = "2")
    plot(SGLD_VS_Q1000[[5]]$result_all[[i]]$logscore[2,], ylab = "log score",
         xlab = paste0(i,'th iteration of beta2' ))
    abline(h=63.32389, col = "2")
}


par(mfrow = c(4,2))
for (i in 91:94) {
    plot(SGLD_VS_Q1000[[5]]$result_all[[i]]$logscore[1,], ylab = "log score",
         xlab = paste0(i,'th iteration of beta1' ))
    abline(h=63.32389, col = "2")
    plot(SGLD_VS_Q1000[[5]]$result_all[[i]]$logscore[2,], ylab = "log score",
         xlab = paste0(i,'th iteration of beta2' ))
    abline(h=63.32389, col = "2")
}

ll <- c()
for (i in 1:1000) {
    l <- length(data_test[[i]]$x)
    ll <- c(ll,l)
}
par(mfrow = c(1,1))
hist(ll, main = "Histogram of length in historical data", ylim = c(0,500))

rank_ls <- c()
for (i in 1:1000) {
    ls <- c(fore_feat_Q1000[[i]]$err_feature[,"lpds"], fore_Q1000[[i]]$logscore)
    names(ls) <- c("Ours", "OP","SA","ets","naive","rw_drift")
                                        # ascending order
    ran <- rank(ls)
    rank_ls <- rbind(rank_ls, ran)
}

ls_win <- apply(rank_ls, 2, function(x){sum(x==6)})
ls_win <- cbind(names(ls_win),data.frame(ls_win), rep("logscore",6))
colnames(ls_win)=NULL
rownames(ls_win)=NULL


rank_mase <- c()
for (i in 1:1000) {
    ls <- c(fore_feat_Q1000[[i]]$err_feature[,"mase_err_h"], fore_Q1000[[i]]$mase_err)
    names(ls) <- c("Ours", "OP","SA","ets","naive","rw_drift")
    ran <- rank(ls, ties.method = "random")
    rank_mase <- rbind(rank_mase, ran)
}

mase_win <- apply(rank_mase, 2, function(x){sum(x==1)})
mase_win <- cbind(names(mase_win),data.frame(mase_win), rep("mase",6))
colnames(mase_win)=NULL
rownames(mase_win)=NULL


rank_smape <- c()
for (i in 1:1000) {
    ls <- c(fore_feat_Q1000[[i]]$err_feature[,"smape_err_h"], fore_Q1000[[i]]$smape_err)
    names(ls) <- c("Ours", "OP","SA","ets","naive","rw_drift")
    ran <- rank(ls, ties.method = "random")
    rank_smape <- rbind(rank_smape, ran)
}

smape_win <- apply(rank_smape, 2, function(x){sum(x==1)})
smape_win <- cbind(names(smape_win),data.frame(smape_win), rep("smape",6))
colnames(smape_win)=NULL
rownames(smape_win)=NULL


win <- rbind(data.frame(ls_win), data.frame(mase_win), data.frame(smape_win))

Q1000_win <- data.frame(method = win[,1], win_num = win[,2], error = win[,3])

library("ggplot2")
ggplot(Q1000_win, aes(x=error, y=win_num, fill=method)) + geom_bar(stat='identity')

ind <- Q1000_win$method %in% c("ets","naive","rw_drift")
Q1000_win0 <- Q1000_win[!ind,]

ggplot(Q1000_win0,aes(x=error, y=win_num, fill=method)) + geom_bar(stat='identity')
