<<<<<<< HEAD
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
#' @note TODO: log_score_grad(beta, features, prob, intercepts)
#' @author Feng Li
#' 
log_score<-function(beta, features, features_select = NULL, prob, intercept){
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
  w_max <- try(optim(par = 0,
                 fn = log_score,
                 features = features_y,
                 prob = prob0,
                 features_select = features_select,
                 intercept = intercept,
                 method="CG",
                 control = list(fnscale = -1)) #max
               )
  
  #if(is(w_max, "try-error")) browser()
  beta_optim <- w_max$par
  return(list(beta_optim=beta_optim,features_select = features_select))
}

cl <- makeCluster(2)
registerDoParallel(cl)
optimal_beta_feature <- foreach(i_ts = 1:length(data)) %dopar% 
  optim_beta_feature(lpd_feature_yearly[[i_ts]], features_select=12, intercept =FALSE)
stopCluster(cl)


#----------------------------------------------------------#
## Function  forecast_results
## param: a ts data with $x, $xx
## param: model_conf
## param: intercept = FALSE
## param: lpd_feature_yearly with $feat
## param: optimal_beta_feature with $beta_optim, $features_select

## return: data[[i_ts]] with $ff_feature, $err_feature

load("E:/time series/R code/feature-based-bayesian-model-averaging/data/data_forecast_yearly.RData")


forecast_feature_results <-
  function(data,
           model_conf,
           intercept = FALSE,
           lpd_feature,
           optimal_beta_feature) {
    library(forecast)
    require("M4metalearning")
    attach(model_conf)
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
    features_y_mean = attr(features_y, "scaled:center")
    features_y_sd = attr(features_y, "scaled:scale")
    
    lpds = 0
    pred_densities = matrix(NA, forecast_h, 2)
    
    for (t in 1:forecast_h)
    {
      ## NOTE: This part is a recursive process, could not be parallelized.  Calculate
      ## predictive features. Each individual model provide an one-step-ahead predictive y
      ## values, and recalculate features based on optimized pools. This will do the
      ## forecasting model multiple times (h), consider to simplify it.
      
      ## ETS model
      ets_fit <- ets(y_new, model = model_conf$ets_model)
      ets_fore <- forecast(ets_fit, h = 1, level = PI_level)
      ets_fore_mean <- ets_fore$mean
      ets_fore_sd = (ets_fore$lower - ets_fore$mean) / qnorm(1 - PI_level /
                                                               100)
      
      ## ARIMA model
      arima_fit <- auto.arima(y_new)
      arima_fore <- forecast(arima_fit, h = 1, level = PI_level)
      arima_fore_mean <- arima_fore$mean
      arima_fore_sd = (arima_fore$lower - arima_fore$mean) / qnorm(1 - PI_level /
                                                                     100)
      
      
      ## Update features
      if (!is.null(features_y))
      {
        #用非标准化的数据计算特征
        myts <- list(list(x = ts(y_new_nonsd, frequency = frequency)))
        myfeatures <- THA_features(myts)[[1]]$features
        myfeatures <- data.matrix(myfeatures)
        myfeatures_scaled = scale(myfeatures, center = features_y_mean, scale = features_y_sd)
        myfeatures_scaled[is.na(myfeatures_scaled)] <- 0
      } else
      {
        myfeatures_scaled = NULL
      }
      
      ## Update predictive weights
      myfeatures_scaled <-
        myfeatures_scaled[, optimal_beta_feature$features_select]
      if (intercept)
        myfeatures_scaled = cbind(1, myfeatures_scaled)
      exp_lin = exp(myfeatures_scaled %*% optimal_beta_feature$beta_optim)
      if (exp_lin == Inf) {
        w <- matrix(1)
      } else{
        w <- exp_lin / (1 + rowSums(exp_lin)) # 1-by-(n-1)
      }
      w_full = cbind(w, 1 - rowSums(w)) # 1-by-n
      
      ## The final pooled y
      y_pred_h = sum(cbind(ets_fore_mean, arima_fore_mean) * w_full)
      y_new = c(y_new, y_pred_h)
      y_new_nonsd = c(y_new_nonsd, (y_pred_h * y_sd + y_mean))
      y_hat_matrix[1, t] <- y_pred_h * y_sd + y_mean
      
      ### The predictive log score
      
      pred_densities[t, 1] <-
        dnorm(y01_true[t], mean = ets_fore_mean, sd = ets_fore_sd)
      pred_densities[t, 2] <-
        dnorm(y01_true[t], mean = arima_fore_mean, sd = arima_fore_sd)
      lpds0 = log(sum(pred_densities[t, ] * w_full))
      lpds = lpds + lpds0
    }
    data$ff_feature <- y_hat_matrix
    
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

cl <- makeCluster(2)
registerDoParallel(cl)
data_forecast_feature <-
  foreach(i_ts = 1:length(data_forecast)) %dopar% 
  forecast_feature_results(
    data_forecast[[i_ts]],
    model_conf,
    intercept = FALSE,
    lpd_feature_yearly[[i_ts]],
    optimal_beta_feature[[i_ts]]
  )
stopCluster(cl)

save(data_forecast_feature, file="E:/time series/R code/feature-based-bayesian-model-averaging/data/data_forecast_feature_yearly.RData")
#----------------------------------------------------------#
## Function  forecast_performance
## param:data with $err_feature
## return: a matrix performance

load("E:/time series/R code/feature-based-bayesian-model-averaging/data/data_forecast_feature_yearly.RData")
forecast_feature_performance<-function(data){ 
  # log score
  performance<-c()
  for (i_ts in 1:length(data_forecast_feature)) {
    performance<-rbind(performance,data_forecast_feature[[i_ts]]$err_feature) 
  }
  performance_out1<- sum(performance[,1])
                                 
  performance_out2<- colMeans(performance[,2:3])
 
  performance_out<-cbind(performance_out1,t(performance_out2))
  colnames(performance_out)<-c("Log score","Mase","Smape")
  return(performance_out)
}

forecast_feature_performance(data_forecast_feature)

#------------------------------------------------------------------------------------# 

