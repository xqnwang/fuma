#' Combinate the forecasts and prediction intervals by simple average method
#' 
#' The combinated forecasts and prediction intervals are calculated 
#' by simple average method. 
#' 
#' @param dataset a list of time series. See details for the required format.
#' @param level the confidence level of prediction intervals that are used for combination.
#' @param method.num the number of methods used for averaging. Methods in the first 
#' \code{method.num} rows in \code{ff} will be selected for simple averaging.
#' @param methodname string. the name used to name the combined forecasting results.
#' @param parallel logical. If \code{TRUE} then the calculations are conducted in parallel.
#' @param num.cores the specified amount of parallel processes to be used if parallel = TRUE.
#' 
#' @details 
#' \code{dataset} must be a list with each element having the following format:
#' \describe{
#'   \item{n}{the number of observations in the time series.}
#'   \item{h}{the number of required forecasts.}
#'   \item{period}{interval of the time series. Possible values are 
#'   "YEARLY", "QUARTERLY", "MONTHLY" & "OTHER".}
#'   \item{x}{a time series of length \code{n} (the historical data).}
#'   \item{xx}{a time series of length \code{h} (the future data).}
#'   \item{ff}{a matrix with \code{F} rows and \code{h} columns. Each row contains
#'   the forecasts of each method in \code{methods}.}
#'   \item{lower}{a list with each element being the matrix of lower bounds 
#'   for certain confidence level.}
#'   \item{upper}{a list with each element being the matrix of upper bounds 
#'   for certain confidence level.}
#' }
#' 
#' @return A list with the elements having the following structure
#' \describe{
#'   \item{ff}{a matrix with \code{F+1} rows and \code{h} columns where 
#'   the first \code{F} rows represent the point forecasts of benchmark methods 
#'   and the last row represents the point forecasts by model averaging.}
#'   \item{lower}{a list with one element that contains a matrix. The matrix with 
#'   \code{F+1} rows and \code{h} columns where the first \code{F} rows represent 
#'   the lower bounds of benchmark methods and the last row represents the 
#'   lower bounds by model averaging.}
#'   \item{upper}{a list with one element that contains a matrix. The matrix with 
#'   \code{F+1} rows and \code{h} columns where the first \code{F} rows represent 
#'   the upper bounds of benchmark methods and the last row represents the 
#'   upper bounds by model averaging.}
#' }
#' 
#' @import stats
#' @importFrom parallel detectCores makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @export
simp_forec <- function(dataset, level, method.num, methodname = "comb", 
                       parallel = FALSE, num.cores = 2) {
  
  combforec_fun <- function(input){
    lapply(input, function(lentry){
      index <- seq_len(method.num)
      methodnames <- rownames(lentry$ff)
      
      ff <- data.frame(lentry$ff)[index,]
      lower <- data.frame(lentry$lower[[which(names(lentry$lower) == paste(level, "%", sep = ""))]])[index,]
      upper <- data.frame(lentry$upper[[which(names(lentry$upper) == paste(level, "%", sep = ""))]])[index,]
      lradius <- as.matrix(ff - lower)
      uradius <- as.matrix(upper - ff)
      ma_ff <- apply(ff, 2, mean)
      ma_lradius <- apply(lradius, 2, mean)
      ma_uradius <- apply(uradius, 2, mean)
      ma_lower <- ma_ff - ma_lradius
      ma_upper <- ma_ff + ma_uradius
      
      lentry$ff <- rbind(lentry$ff, ma_ff)
      low <- rbind(lentry$lower[[which(names(lentry$lower) == paste(level, "%", sep = ""))]], ma_lower)
      up <- rbind(lentry$upper[[which(names(lentry$upper) == paste(level, "%", sep = ""))]], ma_upper)
      rownames(lentry$ff) <- rownames(low) <- rownames(up) <- c(methodnames, methodname)
      lentry$lower <- 'names<-' (list(low), paste(level, "%", sep = ""))
      lentry$upper <- 'names<-' (list(up), paste(level, "%", sep = ""))
      
      return(lentry)
    })
  }
  
  if (parallel == FALSE){
    ret_list <- combforec_fun(dataset)
  } else {
    ncores <- num.cores
    ncores <- ifelse(length(dataset) < ncores, length(dataset), ncores)
    times_sp <- c(rep(floor(length(dataset)/ncores), ncores - 1), 
                  length(dataset) - floor(length(dataset)/ncores) * (ncores - 1))
    data_sp <- split(dataset, rep(seq_len(ncores), times = times_sp))
    cl <- parallel::makeCluster(ncores)
    registerDoParallel(cl)
    ret_list <- foreach(data = data_sp, 
                        .combine = base::c, 
                        .packages = 'stats') %dopar%
      combforec_fun(input = data)
    stopCluster(cl)
  }
  
  ret_list
}
