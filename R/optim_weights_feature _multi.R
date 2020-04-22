#' Calculate the log predictive score for a time series with pools of models
#'
#'
#' @title log predictive score with features
#' @param beta p-by-(n-1) matrix
#' @param features T-by-p feature matrix, usually standardized.
#' @param prob T-by-n predictive densities.
#' @param intercept TRUE or FALSE Should intercept be used in feature weights?
#' @param features_select a vector including the numbers of the features to be taken into consideration
#' @return
#' @references Geweke & Amisano, (2011) Optimal prediction pools, Journal of Econometrics.
#' @author Li Li
#' 

#################################################################################

library(tsfeatures)
library(M4metalearning)
library(forecast)
library(tseries)
library(purrr)
library(ggplot2)
library(M4comp2018)
library(foreach)
library(doParallel)
library(numDeriv)
library("mvtnorm")
library("base")
library(MASS)

# ts
data_test <- M4[[96982]]

### Model Settings
model_conf = list(
  frequency = 1
  , ets_model = "ANN" # simple exponential smoothing with additive errors
  , forecast_h = 14 # Forecasting horizon
  , train_h = 1 # Out of sample training, 1 means rolling forecast
  , history_burn = 100 # Let the model to start with at least this length of historical data.
  , PI_level = 90 # Predictive Interval level, used to extract out-of-sample variance from forecasting models.
  , roll = 100 # The length of rolling samples, larger than history_burn
  , feature_window =100 # The length of moving window when computing features
 # , fore_model = list("ets_fore", "auto.arima_fore" )
  , fore_model = list("ets_fore", "auto.arima_fore", "naive_fore", "rw_drift_fore",
                      "snaive_fore", "stlm_ar_fore", "thetaf_fore")
)

# @references FFORMA: Feature-based Forecast Model Averaging
ets_fore <- function(x, train_h, PI_level) {
  ets_fit <- forecast::ets(x, model = "ANN")
  ets_fore <- forecast(ets_fit, h = train_h, level = PI_level)
  ets_fore_mean <- ets_fore$mean
  ets_fore_sd <- (ets_fore$lower - ets_fore$mean) / qnorm(1 - PI_level / 100)
  return(list(ets_fore_mean = as.numeric(ets_fore_mean), 
              ets_fore_sd = as.numeric(ets_fore_sd) ))
}

auto.arima_fore <- function(x, train_h, PI_level) {
  arima_fit <- forecast::auto.arima(x)
  arima_fore <- forecast(arima_fit, h = train_h, level = PI_level)
  arima_fore_mean <- arima_fore$mean
  arima_fore_sd <- (arima_fore$lower - arima_fore$mean) / qnorm(1 - PI_level /100)
  return(list(arima_fore_mean = as.numeric(arima_fore_mean), 
              arima_fore_sd = as.numeric(arima_fore_sd) ))
}


naive_fore <- function(x, train_h, PI_level) {
  naive_fit <- forecast::naive(x, h = train_h, level = PI_level)
  naive_fore_mean <- naive_fit$mean
  naive_fore_sd <- (naive_fit$lower - naive_fit$mean) / qnorm(1 - PI_level /100)
  return(list(naive_fore_mean = as.numeric(naive_fore_mean), 
              naive_fore_sd = as.numeric(naive_fore_sd) ))
}

rw_drift_fore <- function(x, train_h, PI_level) {
  rw_drift_fit <- forecast::rwf(x, drift=TRUE, h = train_h, level = PI_level)
  rw_drift_fore_mean <- rw_drift_fit$mean
  rw_drift_fore_sd <- (rw_drift_fit$lower - rw_drift_fit$mean) / qnorm(1 - PI_level /100)
  return(list(rw_drift_fore_mean = as.numeric(rw_drift_fore_mean), 
              rw_drift_fore_sd = as.numeric(rw_drift_fore_sd) ))
}

snaive_fore <- function(x, train_h, PI_level) {
  snaive_fit <- forecast::snaive(x, h = train_h, level = PI_level)
  snaive_fore_mean <- snaive_fit$mean
  snaive_fore_sd <- (snaive_fit$lower - snaive_fit$mean) / qnorm(1 - PI_level /100)
  return(list(snaive_fore_mean = as.numeric(snaive_fore_mean), 
              snaive_fore_sd = as.numeric(snaive_fore_sd) ))
}

# nnetar_fore <- function(x, train_h, PI_level) {
#   nnetar_fit <- forecast::nnetar(x)
#   nnetar_fore <- forecast(nnetar_fit, h = train_h, level = PI_level)
#   nnetar_fore_mean <- nnetar_fore$mean
#   nnetar_fore_sd <- 
#   return(list(nnetar_fore_mean = as.numeric(nnetar_fore_mean), 
#               nnetar_fore_sd = as.numeric(nnetar_fore_sd) ))
# }

# tbats_fore <- function(x, train_h, PI_level) {
#   tbats_fit <- forecast::tbats(x, use.parallel=FALSE)
#   tbats_fore <- forecast(tbats_fit, h = train_h, level = PI_level)
#   tbats_fore_mean <- tbats_fore$mean

#   #Prediction intervals are not symmetric with respect to point predictions.
#   tbats_fore_sd <- 
#   return(list(tbats_fore_mean = as.numeric(tbats_fore_mean), 
#               tbats_fore_sd = as.numeric(tbats_fore_sd) ))
# }

stlm_ar_fore <- function(x, train_h, PI_level) {
  stlm_ar_fit <- tryCatch({
    forecast::stlm(x, modelfunction = stats::ar)
  }, error = function(e) forecast::auto.arima(x, d=0,D=0))
  stlm_ar_fore <- forecast(stlm_ar_fit, h = train_h, level = PI_level)
  stlm_ar_fore_mean <- stlm_ar_fore$mean
  stlm_ar_fore_sd <- (stlm_ar_fore$lower - stlm_ar_fore$mean) / qnorm(1 - PI_level /100)
  return(list(stlm_ar_fore_mean = as.numeric(stlm_ar_fore_mean), 
              stlm_ar_fore_sd = as.numeric(stlm_ar_fore_sd) ))
}

thetaf_fore <- function(x, train_h, PI_level) {
  thetaf_fit <- forecast::thetaf(x, h = train_h, level = PI_level)
  thetaf_fore_mean <- thetaf_fit$mean
  thetaf_fore_sd <- (thetaf_fit$lower - thetaf_fit$mean) / qnorm(1 - PI_level /100)
  return(list(thetaf_fore_mean = as.numeric(thetaf_fore_mean), 
              thetaf_fore_sd = as.numeric(thetaf_fore_sd) ))
}

#-------------------------------------------------------------------------#
## Compute log predictive density and features


## Function  lpd_feature
## param: a ts data with $x
## param: model_conf
## return: lpd_feature (a list including two matrices of log probability densities and features)

lpd_feature_multi <- function(data, model_conf) {
  
  feature_window = model_conf$feature_window
  roll = model_conf$roll
  frequency = model_conf$frequency
  history_burn = model_conf$history_burn
  ets_model = model_conf$ets_model
  forecast_h = model_conf$forecast_h
  train_h = model_conf$train_h
  PI_level = model_conf$PI_level
  fore_model = model_conf$fore_model
  
  ## A single historical data
  y <- data$x
  y1 = scale(y, center = TRUE, scale = TRUE)
  y_mean = attr(y1, "scaled:center")
  y_sd = attr(y1, "scaled:scale")
  
  y1 = as.numeric(y1)
  
  ## Calculate historical log predictive density
  num_model <- length(fore_model)
  log_pred_densities <-
    matrix(nrow = length(y) - history_burn - train_h + 1,
           ncol = num_model)
  colnames(log_pred_densities) <- unlist(fore_model)
  
  for (t in (history_burn):(length(y) - train_h))
  {
    ## TODO: We may simplify this to assume the fit and forecast procedure is
    ## invariant with certain time period to save time.
    if(is.null(roll)){
      y01 <- y1[1:t]
    }else if (t < roll){
      y01 <- y1[1:t]
    }else{
      y01 <- y1[(t-roll+1):t]
    }
    
    ## To keep numeric stability, we calculate log P(y_pred)
    use_model <- lapply(fore_model, function(method){
      method_fun <- get(method)
      mean_sd <- method_fun (y01, train_h, PI_level)
      return(mean_sd)
    })
    log_pred_den <- lapply(use_model, function(mean_sd){
      lpd <- sum(dnorm(y1[(t + 1):(t + train_h)],
                       mean = mean_sd[[1]],sd = mean_sd[[2]],log = TRUE ))
      return(lpd)
    })
    log_pred_densities[(t - history_burn + 1), ] <- as.numeric(log_pred_den)
  }
  
  ## Calculate historical features
  features_y <-
    matrix(nrow = length(y) - train_h - history_burn + 1,
           ncol = 42)
  myts <- list(list(x = ts(y[1:history_burn], frequency = frequency)))
  myfeatures <- THA_features(myts)[[1]]$features
  names <- colnames(myfeatures) 
  colnames(features_y) <- names
  
  if(is.null(feature_window)){
    for (t in (history_burn):(length(y) - train_h))
    {
      myts <- list(list(x = ts(y[1:t], frequency = frequency)))
      myfeatures <- THA_features(myts)[[1]]$features
      myfeatures <- data.matrix(myfeatures)
      features_y[(t - history_burn + 1),] <- myfeatures
    }
  }else{
    features_y <- features_window( y, window = feature_window, 
                                   history_burn, train_h)
  }
  
  features_y_scaled = scale(features_y, center = TRUE, scale = TRUE)
  
  lpd_feature <-
    list(lpd = log_pred_densities, feat = features_y_scaled)
  return(lpd_feature)
}
features_window<-function( y, window, history_burn, train_h){
  if(window < history_burn){
    stop("The length of window needs to be larger than history_burn")
  }
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


##Delete the features with NaN 
feature_clean <- function(lpd_feature){
  for (i in 1:length(lpd_feature)) {
    ind <- which(is.nan(lpd_feature[[i]]$feat), arr.ind = T)[,2]
    ind0 <- ind[!duplicated(ind)]
    lpd_feature[[i]]$feat_mean <- attr(lpd_feature[[i]]$feat, "scaled:center")[-ind0]
    lpd_feature[[i]]$feat_sd <- attr(lpd_feature[[i]]$feat, "scaled:scale")[-ind0]
    lpd_feature[[i]]$feat <- lpd_feature[[i]]$feat[, -ind0]
  }
  return(lpd_feature)
}

#------------------------------------------------------------------------------------#
# Optimize the parameter beta 

## SGLD + VS

log_score <- function(beta, features, features_select = NULL, prob, intercept = T){
  if(is.null(features_select)){
    features0<-features
  }else{
    features0<-data.matrix(features[,features_select])
  }
  if(intercept) features0 = cbind(rep(1, nrow(prob)), features0)
  num_model <- dim(prob)[2]
  # Avoid Inf/Inf
  me <- features0 %*% beta
  me[me>709] <- 709
  exp_lin = exp(me)
  deno = matrix (rep((1+rowSums(exp_lin)), num_model-1), ncol = num_model-1)
  w <- exp_lin/ deno # T-by-(n-1)
  w_full = cbind(w, 1 - rowSums(w)) # T-by-n
  out = sum(log(rowSums(w_full * prob)))
  return(out)
}

gradient_logscore <- function(beta, features, features_select = NULL, prob, intercept){
  if(is.null(features_select)){
    features0<-features
  }else{
    features0<-data.matrix(features[,features_select])
  }
  if(intercept) features0 = cbind(rep(1, nrow(prob)), features0)
  
  ex = exp(features0 %*% beta)
  ex[ex > exp(700)] = exp(700)
  
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

proposal_I <- function(n) rbinom(n, 1, 0.5)

# p(beta | I) * p(I)
# beta = [beta1, beta2,..., beta(n-1)]
# beta1,...,beta(n-1), iid,  ~ N(0,sigma), sigma = diag{10,...,10}
prior <- function(beta, I, sig = 10){
  pri <- apply(beta, 2, dmvnorm, mean=matrix(0, dim(beta)[1], 1), sigma = sig*diag(dim(beta)[1]))
  pri <- exp(sum(log(pri)))
  # avoid prior=0
  if (pri < 1e-323) {pri <- 1e-323}
  return(pri)
}

# the gradient of the prior can be computed in this way only when
# beta1,...,beta(n-1), iid,  ~ N(0,sigma)
gradient_prior <- function(beta, I, sig){
  coef = -prior(beta, I, sig)
  sigma_inver = solve(sig*diag(dim(beta)[1]))
  mu = matrix(0, dim(beta)[1], 1)
  beta_gra = apply(beta, 2, function(beta){
    sigma_inver %*% (beta-mu)
  })
  return(coef * beta_gra)
}


log_posterior <- function(data, beta, I, prior = prior, 
                          logLik = log_score, sig = 10){
  pri <- prior(beta, I, sig = sig) 
  features_select <- which(I==1)
  prob <- exp(data$lpd)
  prob[prob == 0] <- 1e-323
  features <- data$feat
  LS <- logLik(beta = beta, features = features, features_select = features_select, 
               prob = prob, intercept = TRUE)
  log_post <- log(pri) + LS
  return(log_post)
}

# no gibbs
SGLD <- function(data, logLik, gradient_logLik, prior, start, minibatchSize = NULL,
                 stepsize = NULL, tol = 1e-5, iter = 5000, samplesize = 0.1,sig = 10,
                 features_select, I, intercept = TRUE, gama = 0.55, a = 0.4, b = 10 ){
  beta <- start
 
  prob <- exp(data$lpd)
  prob[prob == 0] <- 1e-323
  num_model <- dim(prob)[2]
  features <- data$feat
  
  prior0 <- prior(beta, I, sig = sig)
  
  logLik0 <- logLik (beta = beta, features = features, 
                     features_select = features_select, 
                     prob = prob, intercept = intercept)
  logpost0 <- log_posterior (data, beta, I, prior = prior, 
                             logLik = log_score, sig = sig)
  i <- 1
  res <- list(beta = list(start), logscore = logLik0, logpost = logpost0, 
              stepsize = NA, prior = prior0)
  
  if(! is.null(minibatchSize)){
    minibatchSize = minibatchSize
  } else if (length(prob[,1]) <= 10){
    minibatchSize = 1
  }else {
    minibatchSize = 0.1
  }
  
  repeat{
    mini <- sample(1:length(prob[,1]),
                   ceiling(minibatchSize*length(prob[,1])))
    prob1 <- prob[mini,]
    features1 <- features[mini,]
    
    if(is.null(stepsize)){
      stepsize1 <- a * (b + i) ^ (-gama)
    }else{
      stepsize1 <- stepsize
    }
    
    prior1 <- prior(beta, I, sig = sig)
    
    beta <- (beta + stepsize1 * (gradient_prior(beta, I, sig)/ prior1)
             + stepsize1 * (1/minibatchSize) * gradient_logLik(beta = beta, features = features1, 
                                                               features_select = features_select, 
                                                               prob= prob1, intercept= intercept)
             + t(mvrnorm(num_model-1, rep(0,dim(beta)[1]), 2*stepsize1* diag(dim(beta)[1])))
    )
    
    prior1 <- prior(beta, I, sig = sig)
    logLik1 <- logLik (beta = beta, features = features, features_select = features_select, 
                       prob = prob, intercept = intercept)
    logpost1 <- log_posterior (data, beta, I, prior = prior, 
                               logLik = log_score, sig = sig)
    
    res$beta[[i+1]] <- beta
    res$logscore <- rbind(res$logscore, logLik1)
    res$logpost <- rbind(res$logpost, logpost1)
    res$stepsize <- rbind(res$stepsize, stepsize1)
    res$prior <- rbind(res$prior, prior1)
    
    if (i >= iter)
      break
    i <- i+1
  }
  
  # n = samplesize * iter
  # beta_SGLD <- tail(res$beta, n)
  
  # simple average
  #beta_out <- colMeans(beta_SGLD)
  
  # weighted
  # weights <- tail(res$stepsize, n) / sum(tail(res$stepsize, n))
  # beta_out <- colSums(beta_SGLD * matrix(rep(weights,length(start)), ncol = length(start)))
  
  return(res)
}

SGLD_gib <- function(data, logLik, gradient_logLik, prior, start, minibatchSize = NULL,
                 stepsize = NULL, tol = 1e-5, iter = 5000, samplesize = 0.1,sig = 10,
                 features_select, I, intercept = TRUE, gama = 0.55, a = 0.4, b = 10 ){
  beta_all <- start
  
  prob <- exp(data$lpd)
  prob[prob == 0] <- 1e-323
  num_model <- dim(prob)[2]
  features <- data$feat
  
  if(! is.null(minibatchSize)){
    minibatchSize = minibatchSize
  } else if (length(prob[,1]) <= 10){
    minibatchSize = 1
  }else {
    minibatchSize = 0.1
  }
  
  res <- list()
  for (i in 1:(num_model-1)) {
    prior0 <- prior(beta_all, I, sig = sig)
    logLik0 <- logLik (beta = beta_all, features = features, 
                       features_select = features_select, 
                       prob = prob, intercept = intercept)
    logpost0 <- log_posterior (data, beta_all, I, prior = prior, 
                               logLik = log_score, sig = sig)
    res0 <- list(beta = beta_all[,i], logscore = logLik0, logpost = logpost0, 
                stepsize = NA, prior = prior0)
    
    if(iter != 0){
    for (t in 1:iter ) {
      
      beta <- beta_all[,i]
      mini <- sample(1:length(prob[,1]),
                     ceiling(minibatchSize*length(prob[,1])))
      prob1 <- prob[mini,]
      features1 <- features[mini,]
      
      if(is.null(stepsize)){
        stepsize1 <- a * (b + t) ^ (-gama)
      }else{
        stepsize1 <- stepsize
      }
      
      prior1 <- prior(beta_all, I, sig = sig)
      
      beta <- (beta + stepsize1 * (gradient_prior(beta_all, I, sig)[,i]/ prior1)
               + stepsize1 * (1/minibatchSize) 
               * gradient_logLik(beta = beta_all, features = features1, 
                                 features_select = features_select, 
                                 prob= prob1, intercept= intercept)[,i]
               + mvrnorm(1, rep(0,length(beta)), 2*stepsize1* diag(length(beta)))
      )
      beta_all[,i] <- beta
      
      prior1 <- prior(beta_all, I, sig = sig)
      logLik1 <- logLik (beta = beta_all, features = features, features_select = features_select, 
                         prob = prob, intercept = intercept)
      logpost1 <- log_posterior (data, beta_all, I, prior = prior, 
                                 logLik = log_score, sig = sig)
      res0$beta <- cbind(res0$beta, beta)
      res0$logscore <- c(res0$logscore, logLik1)
      res0$logpost <- c(res0$logpost, logpost1)
      res0$stepsize <- c(res0$stepsize, stepsize1)
      res0$prior <- c(res0$prior, prior1)
      
      # For every β, after iter iterations of SGLD, take the mean of the last 10% of samples.
      # If comment if{}, beta_out is a matrix of the last sample of each β.
      if(t == iter){
        if(iter < 10) {ind = 1:(iter+1)} else {
          ind =( iter + 1 - floor(0.1 * iter)): (iter+1)}
        beta_mean <- rowMeans(res0$beta[, ind])
        beta_all[,i] <- beta_mean
      }
    }
    }
    res[[i]] <- res0
  }
  results <- list(beta=list(), logscore = matrix(ncol = iter+1, nrow = num_model-1),
                  logpost = matrix(ncol = iter+1, nrow = num_model-1),
                  stepsize = matrix(ncol = iter+1, nrow = num_model-1),
                  prior = matrix(ncol = iter+1, nrow = num_model-1),
                  beta_out = beta_all)
  for (j in 1:(num_model-1)) {
    results$beta[[j]] <- res[[j]]$beta
    results$logscore[j,] <-res[[j]]$logscore
    results$logpost[j,] <-res[[j]]$logpost
    results$stepsize[j,] <-res[[j]]$stepsize
    results$prior[j,] <-res[[j]]$prior
  }
  return(results)
}



MH_step <- function(x, beta0, data, logp = log_posterior, 
                    proposal = proposal_I, sig = 10){
  beta_start <- beta0
  rownames(beta_start) <- c("0", which(x == 1))
  xp <- proposal(length(x))
  beta1 <- matrix(NA, nrow = (length(which(xp == 1)) + 1), 
                  ncol = dim(beta_start)[2])
  rownames(beta1) <- c("0", which(xp == 1))
  ind <- c("0", which(x==xp & xp==1))
  beta1[ind, ] <- beta_start[ind, ]
  beta1[is.na(beta1)] <-0
  alpha <- min(1, exp(logp(data, beta = beta1, I = xp, prior = prior, 
                           logLik = log_score, sig = sig) - 
                        logp(data, beta = beta_start, I = x, prior = prior, 
                             logLik = log_score, sig = sig)))
  if (runif(1) < alpha){
    accept <- 1
    x <- xp
    beta_start <- beta1
  }else{
    accept <- 0
  }
  return(list(I = x, beta_start = beta_start, accept = accept))  
}


SGLD_VS <- function(data, logLik, gradient_logLik, prior, stepsize = NULL, 
                    SGLD_iter = 100, SGLD_iter_noVS = 10, VS_iter = 100, 
                    minibatchSize = NULL, sig = 10){
  feature_num <- dim(data$feat)[2]
  model_num <- dim(data$lpd)[2]
  I <- matrix(nrow = feature_num, ncol = VS_iter)
  B <- list()
  result_all <- list()
  I[,1] <- rbinom(feature_num,1,0.5)
  accept_num <- 1
  iter <- SGLD_iter
  for (i in 1:VS_iter) {
    if(i == 1){
      features_select <- which(I[,i]==1)
      beta_start <- matrix(runif((length(features_select) + 1) * (model_num-1), -10, 10),
                           ncol = model_num-1)
      rownames(beta_start) <- c("0", features_select)
    }else{
      accept_num <- accept_num +accept
      I[,i] <-I0
      features_select <- which(I[,i]==1)
      beta_start <- beta_start
    }
    res_SGLD <- SGLD_gib(data = data, logLik = log_score, gradient_logLik = gradient_logscore, 
                         prior = prior, start = beta_start, I = I[,i],
                         minibatchSize = minibatchSize, stepsize = stepsize,
                         iter = iter, features_select = features_select, sig = sig )
    B[[i]] <- res_SGLD$beta_out
    result_all[[i]] <- res_SGLD
    I0 <- I[,i]
    MH <- MH_step (x = I0, beta0 = B[[i]], data =data, 
                   logp = log_posterior, proposal = proposal_I)
    I0 <- MH$I
    beta_start <- MH$beta_start
    accept <- MH$accept

    # If I is rejected, the iterations of SGLD is SGLD_iter_noVS.
    # If I is accepted, the iterations of SGLD is SGLD_iter.
    # We can ignore this step (the iterations of SGLD is always SGLD_iter)
    # when commenting if(){}
    if(accept == 0) {iter = SGLD_iter_noVS}else{iter = SGLD_iter}
  }
  return(list(I = I, B = B, result_all = result_all, acceptance = accept_num/VS_iter))
}


#Achieve beta and the selected features in each VS
beta_prepare <- function(res_SGLD_VS){
  beta_pre0 <-list()
  beta_pre <-list()
  for (i in 1 : length(res_SGLD_VS$B)) {
    beta_pre$beta <- res_SGLD_VS$B[[i]]
    beta_pre$features_select <- as.numeric(rownames(res_SGLD_VS$B[[i]]))[-1]
    beta_pre0[[i]] <- beta_pre
  }
  return(beta_pre0)
}

#-------------------------------------------------------------------------------------------------#
## Forecasting with time varying weights (rolling samples)
## Function  forecast_results
## param: a ts data with $x, $xx
## param: model_conf
## param: intercept 
## param: lpd_feature with $feat, $feat_mean, $feat_sd
## param: beta_pre, a list achieved from beta_prepare function
## return: data[[i_ts]] with $ff_feature, $err_feature

forecast_feature_results_multi <-function(data, model_conf, intercept = T, 
                                    lpd_feature, beta_pre) {
  
  #attach(model_conf)
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
    y_new_nonsd1 <- tail(y_new_nonsd, feature_window)
    if (!is.null(features_y))
    {
      #Features are calculated using non-standardized data
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
      #  Avoid Inf
      exp_lin[exp_lin > exp(709)] <- exp(709)
      w <- exp_lin / (1 + rowSums(exp_lin))
      w_full = cbind(w, 1 - rowSums(w))
      return(w_full)
    }
    
    w_full <- sapply(beta_pre, w_get, myfeatures_scaled = myfeatures_scaled)
    w_full_all[[t]] <- w_full
    w_full_mean <- rowMeans(w_full)
    w_full_mean_h <- cbind(w_full_mean_h, w_full_mean)
    
    # forecast
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
    
    # The predictive log score
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

#------------------------------------------------------------------------------------------#

## Forecasting with optimal pool

# When taking no feature into consideration, let features = NULL

optim_beta <- function(lpd_feature, features_y = NULL) {
  y_lpd <- lpd_feature$lpd
  prob <- exp(y_lpd)
  prob[prob == 0] <- 1e-323
  num_model <- length(lpd_feature$lpd[1,])
  ini <-  t(data.matrix(rep(0, num_model-1)))
  
  w_max <- try(optim(
    par = ini,
    fn = log_score,
    #gr = gradient,
    features = features_y,
    prob = prob,
    ## intercept = intercept,
    intercept = TRUE,
    method = "L-BFGS-B",
    control = list(fnscale = -1)
  ) 
  )
  
  if (w_max$convergence != 0) {
    stop("The optimization does not converge")
  }
  beta_optim <- w_max$par
  return(list(beta_optim = beta_optim, logscore = w_max$value))
}


## Function  forecast_results
## param: a ts data with $x, $xx
## param: model_conf
## param: optimal_beta
## return: the ts data with $ff, $mase_err, $smape_err, $logscore
forecast_results_nofea<-function(data, model_conf, optimal_beta){
  
  roll = model_conf$roll
  frequency = model_conf$frequency
  forecast_h = model_conf$forecast_h
  PI_level = model_conf$PI_level
  fore_model = model_conf$fore_model
  
  ## forecasting
  y_hat_matrix <- matrix(ncol = forecast_h, nrow = 2 + length(fore_model))
  rownames(y_hat_matrix)<-c("Optimal pool","SA", unlist(fore_model))
  lpds_all<- matrix(ncol = forecast_h, nrow = 2 + length(fore_model))
  rownames(lpds_all)<-c("Optimal pool","SA", unlist(fore_model))
  
  y <- data$x
  y01 = scale(y, center = TRUE, scale = TRUE)
  y_mean = attr(y01, "scaled:center")
  y_sd = attr(y01, "scaled:scale")
  y01 = as.numeric(y01)
  y_true = data$xx
  y01_true = as.numeric(scale(y_true, center = y_mean, scale = y_sd))
  
  # T by 2+n
  y_new <- matrix(rep(y01,2 + length(fore_model)), ncol = 2 + length(fore_model))
  colnames(y_new) <- c("Optimal pool","SA", unlist(fore_model))
  
  pred_densities = matrix(NA, forecast_h, length(fore_model))
  pred_densities_simple = matrix(NA, forecast_h, length(fore_model))
  pred_densities_single = matrix(NA, forecast_h, length(fore_model))
  
  ## get weights (optimal pool)
  exp_lin = exp(optimal_beta$beta_optim)
  w <- exp_lin/(1+rowSums(exp_lin)) # 1-by-(n-1)
  w_full = cbind(w, 1 - rowSums(w)) # 1-by-n
  
  for (t in 1:forecast_h)
  { 
    if(is.null(roll)){
      y_new1 <- y_new
    }else{
      y_new1 <- tail(y_new, roll)
    }
    
    op_fore <- lapply(fore_model, function(method){
      method_fun <- get(method)
      mean_sd <- method_fun (y_new1[,1], train_h, PI_level)
      return(mean_sd)
    })
    y_pred_op <- w_full %*% (sapply(op_fore, function(mean_sd){
      return(mean_sd[[1]])
    }))
    pd_op <- sapply(op_fore, function(mean_sd){
      dnorm(y01_true[t], mean = mean_sd[[1]], sd = mean_sd[[2]], log = F)
    })
    lpds_op = log(sum(pd_op * w_full))
    
    sa_fore <- lapply(fore_model, function(method){
      method_fun <- get(method)
      mean_sd <- method_fun (y_new1[,2], train_h, PI_level)
      return(mean_sd)
    })
    y_pred_sa <- mean(sapply(sa_fore, function(mean_sd){
      return(mean_sd[[1]])
    }))
    pd_sa <- sapply(sa_fore, function(mean_sd){
      dnorm(y01_true[t], mean = mean_sd[[1]], sd = mean_sd[[2]], log = F)
    })
    lpds_sa = log(mean(pd_sa))
    
    y_pred <- c()
    lpd_model <- c()
    for (i in 1 : length(fore_model)) {
      method_fun <- get(fore_model[[i]])
      mean_sd <- method_fun (y_new1[,(i+2)], train_h, PI_level)
      y_pred <- c(y_pred, mean_sd[[1]])
      lpd <- dnorm(y01_true[t], mean = mean_sd[[1]], sd = mean_sd[[2]], log = T)
      lpd_model <- c(lpd_model, lpd)
    }
    y_pred_all <- c(y_pred_op, y_pred_sa, y_pred)
    y_new <- rbind(y_new, y_pred_all)
    lpds <- c(lpds_op, lpds_sa, lpd_model)
    
    y_hat_matrix[,t] <- y_pred_all
    lpds_all[,t] <- lpds
  }    
  data$logscore<- rowSums(lpds_all)
  
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

