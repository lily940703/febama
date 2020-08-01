

log_score <- function(beta, features, features_select = NULL, prob, intercept = T){
  if(is.null(features_select)){
    features0<-features
  }else{
    features0<-data.matrix(features[,features_select])
  }
  if(intercept) features0 = cbind(rep(1, nrow(prob)), features0)
  num_model <- dim(prob)[2]
  me <- features0 %*% beta
  me[me>709] <- 709
  exp_lin = exp(me)
  deno = matrix (rep((1+rowSums(exp_lin)), num_model-1), ncol = num_model-1)
  w <- exp_lin/ deno # T-by-(n-1)
  w_full = cbind(w, 1 - rowSums(w)) # T-by-n
  out = sum(log(rowSums(w_full * prob)))
  return(out)
}

optim_beta <- function(lpd_feature, features_y = NULL) {
  y_lpd <- lpd_feature$lpd
  prob <- exp(y_lpd)
  prob[prob == 0] <- 1e-323
  num_model <- length(lpd_feature$lpd[1,])
  ini <-  t(data.matrix(rep(0, num_model-1)))
  
  w_max <- try(optim(
    par = ini,
    fn = log_score,
    #fn = logscore1,
    features_select = NULL,
    #gr = gradient,
    features = features_y,
    prob = prob,
    ## intercept = intercept,
    intercept = TRUE,
    #method = "L-BFGS-B",
    method = "BFGS",
    control = list(fnscale = -1, maxit = 10000)
  ) 
  )
  if( w_max$convergence != 0){
    stop("The optimization is not convergent")
  }
  beta_optim <- w_max$par
  return(list(beta_optim = beta_optim, logscore = w_max$value))
}





forecast_results_nofea<-function(data, model_conf, optimal_beta){
  
  roll = model_conf$roll
  frequency = model_conf$frequency
  history_burn = model_conf$history_burn
  ets_model = model_conf$ets_model
  forecast_h = model_conf$forecast_h
  train_h = model_conf$train_h
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
  
  # ## get weights (optimal pool) for log_score
  exp_lin = exp(optimal_beta$beta_optim)
  exp_lin[exp_lin > (1e+307)] <- (1e+307)
  w <- exp_lin/(1+rowSums(exp_lin)) # 1-by-(n-1)
  w_full = cbind(w, 1 - rowSums(w)) # 1-by-n
  
  ## get weights (optimal pool) for logscore1
  # beta <- cbind(optimal_beta$beta_optim, rep(1, length(optimal_beta$beta_optim[,1])))
  # exp_lin = exp(beta)
  # w_full <- exp_lin / (rowSums(exp_lin))
  
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
    lpds[lpds < log(1e-323)] = log(1e-323)
    
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

forecast_performance <- function(data){ 
  # log score
  logscore_all<-c()
  for (i_ts in 1:length(data)) {
    logscore_all<-cbind(logscore_all,data[[i_ts]]$logscore) 
  }
  logscore_all<-data.matrix(rowMeans(logscore_all))
  # mase 
  mase_all<-c()
  for (i_ts in 1:length(data)) {
    mase_all<-cbind(mase_all,data[[i_ts]]$mase_err) 
  }
  mase_all<-data.matrix(rowMeans(mase_all))
  # smape 
  smape_all<-c()
  for (i_ts in 1:length(data)) {
    smape_all<-cbind(smape_all,data[[i_ts]]$smape_err) 
  }
  smape_all<-data.matrix(rowMeans(smape_all))
  performance<-cbind(logscore_all,mase_all,smape_all)
  colnames(performance)<-c("Log score","Mase","Smape")
  return(performance)
}
