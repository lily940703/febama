#' @references FFORMA: Feature-based Forecast Model Averaging
#' @models in 'forecast': ets, auto.arima, naive, rw_drift, snaive, stlm_ar, thetaf
#' @details: Obtain the forecasting sd by assuming a Gaussian interval
#' @param x: A \code{ts} object with the input time series
#' @param train_h: The amount of future time steps to forecast
#' @param PI_level: Predictive Interval level, used to extract out-of-sample variance from forecasting models.
#' @return A list including forecasting mean and sd
#' @export

ets_fore <- function(x, train_h, PI_level) {
  ets_fit <- forecast::ets(x, model = "ANN")
  ets_fore <- forecast::forecast(ets_fit, h = train_h, level = PI_level)
  ets_fore_mean <- ets_fore$mean
  ets_fore_sd <- (ets_fore$lower - ets_fore$mean) / qnorm(1 - PI_level / 100)
  return(list(ets_fore_mean = as.numeric(ets_fore_mean),
              ets_fore_sd = as.numeric(ets_fore_sd) ))
}

auto.arima_fore <- function(x, train_h, PI_level) {
  arima_fit <- forecast::auto.arima(x)
  arima_fore <- forecast::forecast(arima_fit, h = train_h, level = PI_level)
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

stlm_ar_fore <- function(x, train_h, PI_level) {
  stlm_ar_fit <- tryCatch({
    forecast::stlm(x, modelfunction = stats::ar)
  }, error = function(e) forecast::auto.arima(x, d=0,D=0))
  stlm_ar_fore <- forecast::forecast(stlm_ar_fit, h = train_h, level = PI_level)
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

#' @models in 'rugarch': garch, egarch, tgarch
#' @details: parameters p = q = 1
#' @param x: A \code{ts} object with the input time series
#' @param train_h: The amount of future time steps to forecast
#' @param PI_level: The following functions don't use it in the calculation.
#  This is used as an input parameter to facilitate calculation with the models above.
#' @return A list including forecasting mean and sd                          
#' @export   
                          
garch_fore <- function(x, train_h,  PI_level) {
  myspec=ugarchspec()
  myfit = ugarchfit(data=x, spec = myspec, solver="hybrid")
  fore = ugarchforecast(myfit, n.ahead=train_h)
  garch_fore_mean <- fitted(fore)
  garch_fore_sd <- sigma(fore)
  return(list(garch_fore_mean = as.numeric(garch_fore_mean),
              garch_fore_sd = as.numeric(garch_fore_sd) ))
}

egarch_fore <- function(x, train_h, PI_level) {
  myspec=ugarchspec(
    variance.model = list(model = "eGARCH", garchOrder = c(1, 1), submodel = NULL, 
                          external.regressors = NULL, variance.targeting = FALSE),
    mean.model = list(armaOrder = c(1, 1), include.mean = TRUE,
                      archm = FALSE, archpow = 1, arfima = FALSE,
                      external.regressors = NULL, archex = FALSE),
    distribution.model = "norm"
  )
  
  myfit = ugarchfit(myspec, data=x, solver="hybrid")
  fore = ugarchforecast(myfit, n.ahead=train_h)
  egarch_fore_mean <- fitted(fore)
  egarch_fore_sd <- sigma(fore)
  return(list(egarch_fore_mean = as.numeric(egarch_fore_mean),
              egarch_fore_sd = as.numeric(egarch_fore_sd) ))
}


tgarch_fore <- function(x, train_h, PI_level) {
  myspec=ugarchspec(
    variance.model = list(model = "fGARCH", garchOrder = c(1, 1), submodel = "TGARCH", 
                          external.regressors = NULL, variance.targeting = FALSE),
    mean.model = list(armaOrder = c(1, 1), include.mean = TRUE,
                      archm = FALSE, archpow = 1, arfima = FALSE,
                      external.regressors = NULL, archex = FALSE),
    distribution.model = "norm"
  )
  
  myfit = ugarchfit(myspec, data=x, solver="hybrid")
  fore = ugarchforecast(myfit, n.ahead=train_h)
  tgarch_fore_mean <- fitted(fore)
  tgarch_fore_sd <- sigma(fore)
  return(list(tgarch_fore_mean = as.numeric(tgarch_fore_mean),
              tgarch_fore_sd = as.numeric(tgarch_fore_sd) ))
}
