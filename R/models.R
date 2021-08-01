#' The forecasting model ETS basd on \code{forecast::ets}
#'
#' @param x The input time series.
#' @param train_h The amount of future time steps to forecast.
#' @param PI_level Predictive Interval level, used to extract out-of-sample variance from forecasting models.
#'
#' @return A list including forecasting mean and sd.
#' @export

ets_fore <- function(x, train_h, PI_level) {
  ets_fit <- forecast::ets(x, model = "AAN")
  ets_fore <- forecast::forecast(ets_fit, h = train_h, level = PI_level)
  ets_fore_mean <- ets_fore$mean
  ets_fore_sd <- (ets_fore$lower - ets_fore$mean) / qnorm(1 - PI_level / 100)
  return(list(ets_fore_mean = as.numeric(ets_fore_mean),
              ets_fore_sd = as.numeric(ets_fore_sd) ))
}

#' The forecasting model ARIMA basd on \code{forecast::auto.arima}
#'
#' @param x The input time series.
#' @param train_h The amount of future time steps to forecast.
#' @param PI_level Predictive Interval level, used to extract out-of-sample variance from forecasting models.
#'
#' @return A list including forecasting mean and sd.
#' @export
auto.arima_fore <- function(x, train_h, PI_level) {
  arima_fit <- forecast::auto.arima(x)
  arima_fore <- forecast::forecast(arima_fit, h = train_h, level = PI_level)
  arima_fore_mean <- arima_fore$mean
  arima_fore_sd <- (arima_fore$lower - arima_fore$mean) / qnorm(1 - PI_level /100)
  return(list(arima_fore_mean = as.numeric(arima_fore_mean),
              arima_fore_sd = as.numeric(arima_fore_sd) ))
}

#' The forecasting model naive basd on \code{forecast::naive}
#'
#' @param x The input time series.
#' @param train_h The amount of future time steps to forecast.
#' @param PI_level Predictive Interval level, used to extract out-of-sample variance from forecasting models.
#'
#' @return A list including forecasting mean and sd.
#' @export
naive_fore <- function(x, train_h, PI_level) {
  naive_fit <- forecast::naive(x, h = train_h, level = PI_level)
  naive_fore_mean <- naive_fit$mean
  naive_fore_sd <- (naive_fit$lower - naive_fit$mean) / qnorm(1 - PI_level /100)
  return(list(naive_fore_mean = as.numeric(naive_fore_mean),
              naive_fore_sd = as.numeric(naive_fore_sd) ))
}

#' The forecasting model random walk basd on \code{forecast::rwf}
#'
#' @param x The input time series.
#' @param train_h The amount of future time steps to forecast.
#' @param PI_level Predictive Interval level, used to extract out-of-sample variance from forecasting models.
#'
#' @return A list including forecasting mean and sd.
#' @export
rw_drift_fore <- function(x, train_h, PI_level) {
  rw_drift_fit <- forecast::rwf(x, drift=TRUE, h = train_h, level = PI_level)
  rw_drift_fore_mean <- rw_drift_fit$mean
  rw_drift_fore_sd <- (rw_drift_fit$lower - rw_drift_fit$mean) / qnorm(1 - PI_level /100)
  return(list(rw_drift_fore_mean = as.numeric(rw_drift_fore_mean),
              rw_drift_fore_sd = as.numeric(rw_drift_fore_sd) ))
}




#' The forecasting model garch basd on \code{rugarch} package
#'
#' @param x The input time series.
#' @param train_h The amount of future time steps to forecast.
#' @param PI_level This function don't use it in the calculation.
#' This is used to make this function consistent with the parameter setting of the top four models.
#'
#' @return A list including forecasting mean and sd.
#' @export
garch_fore <- function(x, train_h,  PI_level) {
  myspec = rugarch::ugarchspec()
  options(warn =-1)
  myfit = rugarch::ugarchfit(data=x, spec = myspec, solver="hybrid")
  fore = rugarch::ugarchforecast(myfit, n.ahead=train_h)
  rm(myfit)
  garch_fore_mean <- fitted(fore)
  garch_fore_sd <- sigma(fore)
  rm(fore)
  return(list(garch_fore_mean = as.numeric(garch_fore_mean),
              garch_fore_sd = as.numeric(garch_fore_sd) ))
}

#' The forecasting model egarch basd on \code{rugarch} package
#'
#' @param x The input time series.
#' @param train_h The amount of future time steps to forecast.
#' @param PI_level This function don't use it in the calculation.
#' This is used to make this function consistent with the parameter setting of the top four models.
#'
#' @return A list including forecasting mean and sd.
#' @export
egarch_fore <- function(x, train_h, PI_level) {
  myspec_e = rugarch::ugarchspec(
    variance.model = list(model = "eGARCH", garchOrder = c(1, 1), submodel = NULL,
                          external.regressors = NULL, variance.targeting = FALSE),
    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,
                      archm = FALSE, archpow = 1, arfima = FALSE,
                      external.regressors = NULL, archex = FALSE),
    distribution.model = "norm"
  )
  options(warn =-1)
  myfit_e = rugarch::ugarchfit(myspec_e, data=x, solver="hybrid")
  fore = rugarch::ugarchforecast(myfit_e, n.ahead=train_h)
  rm(myfit_e)
  egarch_fore_mean <- fitted(fore)
  egarch_fore_sd <- sigma(fore)
  rm(fore)
  return(list(egarch_fore_mean = as.numeric(egarch_fore_mean),
              egarch_fore_sd = as.numeric(egarch_fore_sd) ))
}


#' The forecasting model garch basd on \code{stochvol} package
#'
#' @param x The input time series.
#' @param train_h The amount of future time steps to forecast.
#' @param PI_level This function don't use it in the calculation.
#' This is used to make this function consistent with the parameter setting of the top four models.
#'
#' @return A list including forecasting mean and sd.
#' @export
sv_fore <- function(x, train_h,  PI_level) {
  options(warn =-1)
  myfit = stochvol::svsample(y = x, quiet	= T)
  fore = predict(myfit, train_h)
  sv_fore_mean <- mean(fore$y)
  sv_fore_sd <- sqrt(mean(exp(fore$h)))
  return(list(sv_fore_mean = as.numeric(sv_fore_mean),
              sv_fore_sd = as.numeric(sv_fore_sd) ))
}

