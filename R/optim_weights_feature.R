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

gradient<-function(beta, features, features_select = NULL, prob, intercept){
  if(is.null(features_select)){
    features0<-features
  }else{
    features0<-data.matrix(features[,features_select])
  }
  if(intercept) features0 = cbind(rep(1, nrow(prob)), features0)
  
  ex = exp(features0 %*% beta)
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
  w_max <- try(optim(par = c(0,0),
                 fn = log_score,
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
## Function  SGLD

gra_prior <- function (x, mean, sigma){
  n <- length(x)
  gra <- (
    -((2 * pi)^(-n/2)) * (det(sigma)^(-1/2)) 
    * as.numeric( exp(-0.5 * t(x-mean) %*% solve(sigma) %*% (x-mean)) )
    * solve(sigma) %*% (x-mean)
  )
  return(gra)
}

library(numDeriv)
library(mvtnorm)
library(base)
library(MASS)
SGLD <- function(data, logLik, gradient_logLik, prior, gradient_prior, start, minibatchSize=0.01,
                 stepsize = 0.01, tol = 1e-5, iter = 5000, samplesize = 0.1,
                 features_select, intercept, gama = 0.55, a = 0.15, b = 100 ){
  library(mvtnorm)
  library(MASS)
  beta <- start
  
  prob <- exp(data$lpd)
  # when pd=0, delete or assign a minimum?
  prob[prob == 0] <- 1e-323
  
  features <- data$feat
  logLik0 <- logLik (beta = beta, features = features, features_select = features_select, 
                     prob = prob, intercept = intercept)
  pri0 <- prior(beta, mean=matrix(0, length(beta), 1),sigma = diag(length(beta)))
  post0 <- pri0 * exp(logLik0)
  i <- 1
  res <- list(beta = start, logscore = logLik0, posterior = post0, stepsize = NA)
  
  if(length(prob[,1]) <= 10){
    minibatchSize = 1
  }else if(length(prob[,1]) > 10 & length(prob[,1]) <=100 ){
    minibatchSize = 0.1
  }else {
    minibatchSize = 0.01
  }
  
  repeat{
    mini <- sample(1:length(prob[,1]),
                   ceiling(minibatchSize*length(prob[,1])))
    prob1 <- prob[mini,]
    features1 <- features[mini,]
    
    stepsize1 <- a*(b + i)^(-gama)
    beta <- (beta + stepsize1 * (gradient_prior(beta, mean=matrix(0, length(beta), 1),sigma = diag(length(beta)))
                                 / prior(as.vector(beta), mean=matrix(0, length(beta), 1),sigma = diag(length(beta)))
    )
    + stepsize1 * (1/minibatchSize) * gradient_logLik(beta = beta, features = features1, 
                                                      features_select = features_select, 
                                                      prob= prob1, intercept= intercept)
    + mvrnorm(1, rep(0,length(beta)), 2*stepsize1* diag(length(beta)))
    )
    
    logLik1 <- logLik (beta = beta, features = features, features_select = features_select, 
                       prob = prob, intercept = intercept)
    pri1 <- prior(as.vector(beta), mean=matrix(0, length(beta), 1),sigma = diag(length(beta)))
    post1 <- pri1 * exp(logLik1)
    
    res$beta <- rbind(res$beta, t(beta))
    res$logscore <- rbind(res$logscore, logLik1)
    res$posterior <- rbind(res$posterior, post1)
    res$stepsize <- rbind(res$stepsize, stepsize1)
    
    if (i >= iter)
      break
    i <- i+1
  }
  
  n = samplesize * iter
  beta_SGLD <- tail(res$beta, n)
  
  # simple average
  #beta_out <- colMeans(beta_SGLD)
  
  # weighted
  weights <- tail(res$stepsize, n) / sum(tail(res$stepsize, n))
  beta_out <- colSums(beta_SGLD * matrix(rep(weights,length(start)), ncol = length(start)))
  
  return(list(res = res, beta_optim = beta_out, features_select = features_select ))
}

cl <- makeCluster(2)
registerDoParallel(cl)
set.seed(2020-3-27)
SGLD_beta_feature10 <- foreach(i_ts = 1:length(data)) %dopar% 
  SGLD (data = lpd_feature_yearly[[i_ts]], logLik = log_score , gradient_logLik=gradient,
               prior = dmvnorm, gradient_prior = gra_prior, start = c(0,0), minibatchSize=0.01,
               stepsize = 0.01, tol = 1e-5, iter = 1000, samplesize = 0.1,
               features_select = c(10), intercept = TRUE, gama = 0.55, a = 0.05, b = 10
  )
stopCluster(cl)

## plot samples

library(ggplot2)
library(plotly)

nx <- ny <- 100
xg <- seq(-5, 0, len = nx)
yg <- seq(-5, 0, len = ny)
g <- expand.grid(xg, yg)
z <- c()
for (i in 1:(nx*ny)) {
  z <- rbind(z, log_score (t(g[i,]), lpd_feature0$feat, features_select = c(10), 
                           exp(lpd_feature0$lpd), intercept=TRUE))
}
f_long <- data.frame(x = g[,1], y = g[,2], z = z)

# all samples
beta_sample <- as.data.frame(tail(res1$res$beta,1001))
out <- ggplot(f_long, aes(x, y, z = z)) + 
  geom_raster(aes(fill = z)) + 
  geom_contour(colour = "white", bins = 20) +
  guides(fill = FALSE) +
  geom_path(data = beta_sample, aes(V1, V2, z=0), col = 2, arrow = arrow()) +
  geom_point(data = beta_sample, aes(V1, V2, z=0), size = 1, col = 2) +
  geom_point(data = tail(beta_sample,1), aes(V1, V2, z=0), size = 2, col = 4) +
  geom_point(data = as.data.frame(t(res1$beta)), aes(V1, V2, z=0), size = 2, col = 1) +
  xlab("beta0") +
  ylab("beta1") +
  ggtitle("SGLD")
ggplotly(out)

# tail(,100)
beta_sample <- as.data.frame(tail(res1$res$beta,100))
out <- ggplot(f_long, aes(x, y, z = z)) + 
  geom_raster(aes(fill = z)) + 
  geom_contour(colour = "white", bins = 20) +
  guides(fill = FALSE) +
  geom_path(data = beta_sample, aes(V1, V2, z=0), col = 2, arrow = arrow()) +
  geom_point(data = beta_sample, aes(V1, V2, z=0), size = 1, col = 2) +
  geom_point(data = tail(beta_sample,1), aes(V1, V2, z=0), size = 2, col = 4) +
  geom_point(data = as.data.frame(t(res1$beta)), aes(V1, V2, z=0), size = 2, col = 1) +
  xlab("beta0") +
  ylab("beta1") +
  ggtitle("SGLD for Posterior Sampling (100)")
ggplotly(out)

# path
beta_sample <- as.data.frame(tail(res1$res$beta,2001))
nr <- nrow(beta_sample)
ind <- c(1:10, (1:floor(nr/100))*100, nr)
beta_sample <- beta_sample[ind,]
out <- ggplot(f_long, aes(x, y, z = z)) + 
  geom_raster(aes(fill = z)) + 
  geom_contour(colour = "white", bins = 20) +
  guides(fill = FALSE) +
  geom_path(data = beta_sample, aes(V1, V2, z=0), col = 2, arrow = arrow()) +
  # geom_segment(aes(x = tail(beta_sample$V1,-1), y = tail(beta_sample$V2,-1), z=0,
  #                  xend = head(beta_sample$v1,-1), yend = head(beta_sample$V2,-1), zend=0,
  #                  colour = "segment"), 
  #              data = beta_sample) +
  geom_point(data = beta_sample, aes(V1, V2, z=0), size = 2, col = 2) +
  geom_point(data = tail(beta_sample,1), aes(V1, V2, z=0), size = 2, col = 4) +
  geom_point(data = as.data.frame(t(res1$beta)), aes(V1, V2, z=0), size = 2, col = 1) +
  xlab("beta0") +
  ylab("beta1") +
  ggtitle("The convergence path of SGLD")
ggplotly(out)

# ggplot stepsize
step<-function(gama = 0.55, a = 0.15, b = 100, t){
  steps <- a*(b + t)^(-gama)
  return(steps)
}
x <- seq(1,5000,1)
steps <- step(t = x)
data <- as.data.frame(cbind(x,steps))
ggplot(data, aes(x , steps)) +
  geom_line(color="red")+xlab(" iteration")+ylab("stepsize") 

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
    
    forecast_h = model_conf$forecast_h
    
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
      if (intercept){
        myfeatures_scaled = cbind(1, t(myfeatures_scaled))
      }
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
      if(lpds0 == -Inf){
        lpds0 = log(1e-323)
      }else{
        lpds0 =lpds0 
      }
        
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
data_forecast_feature_SGLD10 <- foreach(i_ts = 1:length(data)) %dopar% 
  forecast_feature_results(
    data[[i_ts]],
    model_conf,
    intercept = TRUE,
    lpd_feature_yearly[[i_ts]],
    SGLD_beta_feature10[[i_ts]]
  )
stopCluster(cl)

save(data_forecast_feature_SGLD, file="E:/time series/R code/febama/data/data_forecast_feature_yearly_SGLD.RData")
load("E:/time series/R code/febama/data/data_forecast_feature_yearly_SGLD.RData")

for (i in 1:1000) {
  if(data_forecast_feature_SGLD[[i]]$err_feature[1] == -Inf)
    print(i)
}

save(data_forecast_feature, file="E:/time series/R code/feature-based-bayesian-model-averaging/data/data_forecast_feature_yearly.RData")
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

