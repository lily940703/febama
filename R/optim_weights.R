#' Calculate the log predictive score for a time series with pools of models
#'
#'
#' @title log predictive score with features
#' @param beta p-by-(n-1) matrix
#' @param features T-by-p feature matrix, usually standardized.
#' @param prob T-by-n predictive densities, currently n=2.
#' @param intercept TRUE or FALSE Should intercept be used in feature weights?
#' @return
#' @references Geweke & Amisano, (2011) Optimal prediction pools, Journal of Econometrics.
#' @note TODO: log_score_grad(beta, features, prob, intercepts)
#' @author Feng Li
log_score<-function(beta, features, prob, intercept){

  if(intercept) features = cbind(rep(1, nrow(prob)), features)

  exp_lin = exp(features%*%beta)

  w <- exp_lin/(1+rowSums(exp_lin)) # T-by-(n-1)

  ## Full (T-by-n) matrix. To keep identification, only first (n-1) are connected with
  ## features. TODO: Common features?
  w_full = cbind(w, 1 - rowSums(w)) # T-by-n

  out = sum(log(rowSums(w_full * prob)))
  return(out)
}

#################################################################################
## performance of five methods
## 1 feature based (42 features)
## 2 optimal pool
## 3 simple average (SA)
## 4 ARIMA
## 5 ETS

## load historical feature data
setwd("~/code/febama")
load("data/historical_log_pred_features.RData")
for (i_ts in 1:length(data))
{
  ## Rewrite this loop with palapply

  y <- data[[i_ts]]$x
  y_true = data[[i_ts]]$xx
  y_lpd <- lpred_dens[[i_ts]]
  features_y = feat[[i_ts]]

  ## Historical Features and scale information
  features_y_mean = attr(features_y, "scaled:center")
  features_y_sd = attr(features_y, "scaled:scale")

  ## maximizing TODO: change to a better optimization tool.
  ## library("optimx")
  w_max <- optim(par = runif(1),
                 fn = log_score,
                 ## features = features_y,
                 features = NULL,
                 prob = exp(y_lpd),
                 ## intercept = intercept,
                 intercept = TRUE,
                 method="BFGS",
                 control = list(fnscale = -1))

  if(w_max$convergence!=0){
    cat("The optimization does not converge in data", a)
  }
  beta_optim <- w_max$par

  ## optimal pool: feature=NULL, intercept =TRUE
  ## w_optim<-optim(fn=log_score, par=runif(43, min = 0, max = 0),
  ##              features = NULL,
  ##              prob = p,
  ##              intercept = TRUE,
  ##              method="SANN", control = list(fnscale = -1))

  ## forecasting
  features_y_hat = matrix(nrow = forecast_h, ncol = 42)
  y_new_list = matrix(nrow = forecast_h, 1)
  y_new = y

  for (t in 1:forecast_h)
  { ## NOTE: This part is a recursive process, could not be parallelized.  Calculate
    ## predictive features. Each individual model provide an one-step-ahead predictive y
    ## values, and recalculate features based on optimized pools. This will do the
    ## forecasting model multiple times (h), consider to simplify it.

    ## ETS model
    ets_fit <- ets(y_new, model = ets_model)
    ets_fore<-forecast(ets_fit, h = 1, level = PI_level)
    ets_fore_mean <- ets_fore$mean
    ets_fore_sd = (ets_fore$lower - ets_fore$mean)/qnorm(1 - PI_level/100)

    ## ARIMA model
    arima_fit <- auto.arima(y_new)
    arima_fore <- forecast(arima_fit, h = forecast_h, level = PI_level)
    arima_fore_mean <- arima_fore$mean
    arima_fore_sd = (ari_fore$lower - ari_fore$mean)/qnorm(1 - PI_level/100)

    ## Update features
    if(!is.null(features_y))
    {
      myts <- list(list(x=ts(y_new, frequency = frequency)))
      myfeatures <- THA_features(my)[[1]]$features
      myfeatures <- data.matrix(myfeatures)
      myfeatures_scaled = scale(myfeatures, center = features_y_mean, scale = features_y_sd)
      ## features_y_hat[t, ] <- myfeatures_scaled
    }
    else
    {
      myfeatures_scaled = NULL
    }

    ## Update predictive weights
    if(intercept) my_features_scaled = cbind(1, myfeatures_scaled)

    exp_lin = exp(my_features_scaled%*%beta_optim)
    w <- exp_lin/(1+rowSums(exp_lin)) # T-by-(n-1)
    w_full = cbind(w, 1 - rowSums(w)) # T-by-n

    ## The final pooled y
    y_new = sum(cbind(ets_fore_mean, arima_fore_mean) * w_full)
    y_new_list[t] = y_new

    ## The predictive log score, we need true y here.
    ## log_pred_densities[, 1] <- sum(dnorm(y[(t + 1):(t + h)], mean = ets_fore_mean,
    ##                                      sd = ets_fore_sd, log = TRUE))
    ## log_pred_densities[, 2] <- sum(dnorm(y[(t + 1):(t + h)], mean = arima_fore_mean,
    ##                                      sd = arima_fore_sd, log = TRUE))
    ## out = sum(log(rowSums(w_full * prob)))
  }

}

## MASE and SMAPE
mase_err<-c()
smape_err<-c()
score_output<-c()
for (j in 1:5){
  smape_err0<-c()
  mase_err0<-c()
  score0<-c()
  for (i in 1:1000){
    mase_err0<-c(mase_err0,M4_q1[[i]]$mase_err[j])
    smape_err0<-c(smape_err0,M4_q1[[i]]$smape_err[j])

    if(i == 173){
      score0<-score0
    }else{
      score0<-c(score0,M4_q1[[i]]$score[j])
    }
    score0<-mean(score0)
    score_output<-c(score_output,score0)

  }
  mase_err0<-mean(mase_err0)
  mase_err<-c(mase_err,mase_err0)

  smape_err0<-mean(smape_err0)
  smape_err<-c(smape_err,smape_err0)

}

pred_summary<-rbind(mase_err,smape_err,score_output)
colnames(pred_summary)<-c("y_hat_feature","y_hat_w","average","ari_fore","ets_fore")
rownames(pred_summary)<-c("mase_err","smape_err","log_score")
pred_summary
