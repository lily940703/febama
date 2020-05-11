# @references FFORMA: Feature-based Forecast Model Averaging
ets_fore <- function(x, train_h, PI_level) {
  ets_fit <- forecast::ets(x, model = "ANN")
  ets_fore <- forecast(ets_fit, h = train_h, level = PI_level)
  ets_fore_mean <- ets_fore$mean

  ## Obtain the forecasting sd by assuming a Gaussian interval
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
#   #预测区间不关于点预测对称
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
