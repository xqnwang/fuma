#' Calculate point forecasts and prediction intervals for a time series dataset
#' 
#' Point forecasts and prediction intervals of all forecasting methods 
#' in \code{methods} are calculated for each time series in \code{dataset}. 
#' @param dataset a list containing the time series. See details for the required format.
#' @param methods a list of strings with the names of the functions that calculate 
#' point forecasts and prediction intervals for time series.
#' @param level the confidence levels for prediction intervals, such as 80, 90.
#' @param parallel logical. If \code{TRUE} then the calculations are conducted in parallel.
#' @param num.cores the specified amount of parallel processes to be used if parallel = TRUE.
#' 
#' @details 
#' \code{dataset} must be a list with each element having the following format:
#' \describe{
#'   \item{x}{a time series object \code{ts} with the historical data.}
#'   \item{h}{the amount of future time steps to forecast.}
#' }
#' \code{methods} a list of strings with the names of the functions 
#' that calculate point forecasts and prediction intervals for time series. 
#' The functions must exist and take as parameters (\code{x}, \code{h}, \code{level}), where
#' \code{x} is the \code{ts} object with the input series, 
#' \code{h} is the amount of future time steps to forecast 
#' and \code{level} denotes confidence levels for prediction intervals. 
#' The output of these functions must be
#' a list with \code{forecs}(point forecasts), \code{fitted}(fitted values), 
#' \code{pil}(lower bounds of prediction intervals) and 
#' \code{piu}(upper bounds of prediction intervals)
#'
#' @return A list with the elements having the following structure
#' \describe{
#'   \item{x}{a time series object \code{ts} with the historical data.}
#'   \item{h}{the amount of future time steps to forecast.}
#'   \item{f}{a matrix with \code{F} rows and \code{n} columns. Each row contains
#'   the fitted values of each method in \code{methods}.}
#'   \item{ff}{a matrix with \code{F} rows and \code{h} columns. Each row contains
#'   the forecasts of each method in \code{methods}.}
#'   \item{lower}{a list with each element being the matrix of lower bounds 
#'   for certain confidence level.}
#'   \item{upper}{a list with each element being the matrix of upper bounds 
#'   for certain confidence level.}
#' }
#'
#' @examples
#' auto_arima_fun <- function(x, h, level) {
#'   model <- forecast::auto.arima(x, stepwise=FALSE, approximation=FALSE)
#'   forecs <- forecast::forecast(model, h=h)$mean
#'   fitted <- forecast::forecast(model, h=h)$fitted
#'   pil <- forecast::forecast(model, h=h, level=level)$lower
#'   piu <- forecast::forecast(model, h=h, level=level)$upper
#'   list(forecs=forecs, fitted=fitted, pil=pil, piu=piu)
#' }
#' ets_fun <- function(x, h, level) {
#'   # for ets method, residuals != x - fitted
#'   model <- forecast::ets(x)
#'   forecs <- forecast::forecast(model, h=h)$mean
#'   fitted <- forecast::forecast(model, h=h)$fitted
#'   pil <- forecast::forecast(model, h=h, level=level)$lower
#'   piu <- forecast::forecast(model, h=h, level=level)$upper
#'   list(forecs=forecs, fitted=fitted, pil=pil, piu=piu)
#' }
#' tbats_fun <- function(x, h, level) {
#'   # for tbats method, residuals != mean - fitted
#'   model <- forecast::tbats(x, use.parallel=FALSE)
#'   forecs <- forecast::forecast(model, h=h)$mean
#'   fitted <- forecast::forecast(model, h=h)$fitted
#'   pil <- forecast::forecast(model, h=h, level=level)$lower
#'   piu <- forecast::forecast(model, h=h, level=level)$upper
#'   list(forecs=forecs, fitted=fitted, pil=pil, piu=piu)
#' }
#' create_example_list <- function() {
#'   methods_list <- list("auto_arima_fun")
#'   methods_list <- append(methods_list, "ets_fun")
#'   methods_list <- append(methods_list, "tbats_fun")
#'   methods_list
#' }
#' methods <- create_example_list()
#' forec_results <- ts_forec(Mcomp::M3[1:5], methods, level = c(80, 90))
#'
#' @references Montero-Manso et al. (2018).
#' @importFrom parallel detectCores makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @export
ts_forec <- function(dataset, methods, level = c(80,90), 
                     parallel = FALSE, num.cores = 2) {
  
  forec_fun <- function(input){
    lapply(input, function (seriesdata) {
      results <- process_methods_list(seriesdata, methods, level)
      method_names <- sapply(results, function (resentry) resentry$method_name)
      f <- t(sapply(results, function (resentry) resentry$f))
      ff <- t(sapply(results, function (resentry) resentry$ff))
      lower <- NULL
      upper <- NULL
      for (i in seq_len(length(level))) {
        level_name <- paste(level[i], sep = "", "%")
        lo <- t(sapply(results, function (resentry) 
          resentry$lower[, colnames(resentry$lower) == level_name]))
        up <- t(sapply(results, function (resentry) 
          resentry$upper[, colnames(resentry$upper) == level_name]))
        row.names(lo) <- row.names(up) <- gsub("_fun", "", method_names)
        lower <- append(lower, list(lo))
        upper <- append(upper, list(up))
      }
      row.names(f) <- row.names(ff) <- gsub("_fun", "", method_names)
      names(lower) <- names(upper) <- paste(level, sep = "", "%")
      
      seriesdata$f <- f
      seriesdata$ff <- ff
      seriesdata$lower <- lower
      seriesdata$upper <- upper
      seriesdata
    })
  }
  
  if (parallel == FALSE){
    ret_list <- forec_fun(dataset)
  } else {
    ncores <- num.cores
    ncores <- ifelse(length(dataset) < ncores, length(dataset), ncores)
    times_sp <- c(rep(floor(length(dataset)/ncores), ncores - 1), 
                  length(dataset) - floor(length(dataset)/ncores) * (ncores - 1))
    data_sp <- split(dataset, rep(seq_len(ncores), times = times_sp))
    cl <- parallel::makeCluster(ncores)
    .env <- c(sapply(methods, function(method) method[1]), "process_methods_list")
    registerDoParallel(cl)
    ret_list <- foreach(data = data_sp, 
                        .combine = base::c, 
                        .export = .env,
                        .packages = 'forecast') %dopar%
      forec_fun(input = data)
    stopCluster(cl)
  }
  ret_list
}

# process each method in methods_list to produce the fitted values, forecasts and prediction intervals
process_methods_list <- function(seriesdata, methods_list, level) {
  lapply(methods_list, function (mentry) {
    method_name <- mentry
    method_fun <- get(mentry)
    forecasts <- tryCatch(method_fun(x=seriesdata$x, h=seriesdata$h, level),
                          error=function(error) {
                            print(error)
                            print(paste("ERROR processing series: ", seriesdata$st))
                            print(paste("The forecast method that produced the error is:",
                                        method_name))
                            #print("Returning snaive forecasts instead")
                            #snaive_fun(seriesdata$x, seriesdata$h, level)
                          })
    list(f=forecasts$fitted, ff=forecasts$forecs, lower=forecasts$pil, upper=forecasts$piu, method_name=method_name)  
  })
}

