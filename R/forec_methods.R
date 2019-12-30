#' Benchmark methods list
#' 
#' A list of the forecasting methods used for calculating 
#' the point forecasts and prediction intervals with different confidence levels.
#' @import forecast
#' @importFrom stats ar
#' @export
methods_list <- function() {
  methods_list <- list("auto_arima_fun")
  methods_list <- append(methods_list, "ets_fun")
  methods_list <- append(methods_list, "tbats_fun")
  methods_list <- append(methods_list, "stlm_ar_fun")
  methods_list <- append(methods_list, "rw_drift_fun")
  methods_list <- append(methods_list, "thetaf_fun")
  methods_list <- append(methods_list, "naive_fun")
  methods_list <- append(methods_list, "snaive_fun")
  methods_list
}


#' @describeIn methods_list forecast::snaive
#' @param x a \code{ts} object with the input time series
#' @param h the amount of future time steps to forecast
#' @param level the confidence levels for prediction intervals.
#' @export
snaive_fun <- function(x,h,level) {
  model <- forecast::snaive(x, level=level, h)
  forecs <- model$mean
  fitted <- model$fitted
  pil <- model$lower
  piu <- model$upper
  list(forecs=forecs, fitted=fitted, pil=pil, piu=piu)
}

#' @describeIn methods_list forecast::naive
#' @export
naive_fun <- function(x,h,level) {
  model <- forecast::naive(x, level=level, h)
  forecs <- model$mean
  fitted <- model$fitted
  pil <- model$lower
  piu <- model$upper
  list(forecs=forecs, fitted=fitted, pil=pil, piu=piu)
}

#' @describeIn methods_list forecast::auto.arima
#' @export
auto_arima_fun <- function(x, h, level) {
  model <- forecast::auto.arima(x, stepwise=FALSE, approximation=FALSE)
  forecs <- forecast::forecast(model, h=h)$mean
  fitted <- forecast::forecast(model, h=h)$fitted
  pil <- forecast::forecast(model, h=h, level=level)$lower
  piu <- forecast::forecast(model, h=h, level=level)$upper
  list(forecs=forecs, fitted=fitted, pil=pil, piu=piu)
}

#' @describeIn methods_list forecast::ets
#' @export
ets_fun <- function(x, h, level) {
  # for ets method, residuals != x - fitted
  model <- forecast::ets(x)
  forecs <- forecast::forecast(model, h=h)$mean
  fitted <- forecast::forecast(model, h=h)$fitted
  pil <- forecast::forecast(model, h=h, level=level)$lower
  piu <- forecast::forecast(model, h=h, level=level)$upper
  list(forecs=forecs, fitted=fitted, pil=pil, piu=piu)
}

#' @describeIn methods_list forecast::tbats
#' @export
tbats_fun <- function(x, h, level) {
  # for tbats method, residuals != mean - fitted
  model <- forecast::tbats(x, use.parallel=FALSE)
  forecs <- forecast::forecast(model, h=h)$mean
  fitted <- forecast::forecast(model, h=h)$fitted
  pil <- forecast::forecast(model, h=h, level=level)$lower
  piu <- forecast::forecast(model, h=h, level=level)$upper
  list(forecs=forecs, fitted=fitted, pil=pil, piu=piu)
}

#' @describeIn methods_list forecast::stlm with ar modelfunction
#' @export
stlm_ar_fun <- function(x, h, level) {
  model <- tryCatch({
    forecast::stlm(x, modelfunction = stats::ar)
  }, error = function(e) forecast::auto.arima(x, d=0,D=0))
  forecs <- forecast::forecast(model, h=h)$mean
  fitted <- forecast::forecast(model, h=h)$fitted
  pil <- forecast::forecast(model, h=h, level=level)$lower
  piu <- forecast::forecast(model, h=h, level=level)$upper
  list(forecs=forecs, fitted=fitted, pil=pil, piu=piu)
}

#' @describeIn methods_list forecast::rwf
#' @export
rw_drift_fun <- function(x, h, level) {
  model <- forecast::rwf(x, drift=TRUE, h=h, level=level)
  forecs <- model$mean
  fitted <- model$fitted
  pil <- model$lower
  piu <- model$upper
  list(forecs=forecs, fitted=fitted, pil=pil, piu=piu)
}


#' @describeIn methods_list forecast::thetaf
#' @export
thetaf_fun <- function(x, h, level) {
  model <- forecast::thetaf(x, h=h, level=level)
  forecs <- model$mean
  fitted <- model$fitted
  pil <- model$lower
  piu <- model$upper
  colnames(pil) <- colnames(piu) <- paste(level, sep = "", "%")
  list(forecs=forecs, fitted=fitted, pil=pil, piu=piu)
}

