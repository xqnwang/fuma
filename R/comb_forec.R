#' Combinate the forecasts and prediction intervals by assigned weights and threshold ratio
#' 
#' The combinated forecasts and prediction intervals are calculated 
#' by weights given by \code{weightprob}. Besides, the methods used for model averaging 
#' are filtered by threshold ratios for yearly, quarterly and monthly series.
#' 
#' @param dataset a list of time series. See details for the required format.
#' @param weightprob a matrix with \code{n} (number of time series) rows and 
#' \code{F} (number of benchmark methods) matrix. The matrix is applied to assign weights 
#' to each benchmark method for each time series.
#' @param Threshold a vector with 3 threshold ratios that are used for selecting 
#' appropriate methods for yearly, quarterly and monthly series, respectively.
#' @param level the confidence level of prediction intervals that are used for combination.
#' @param combine string representing the interval combination methods. 
#' "mean", "weightedmean" or "median".
#' @param methodname string. the name used to name the combined forecasting results.
#' @param show.methods logical. If \code{TRUE}, then each element of the returned dataset will be 
#' appended to a vector. This vector contains the names of methods that are selected for 
#' model averaging.
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
comb_forec <- function(dataset, weightprob, Threshold, level, 
                       combine = "mean", methodname = "comb", 
                       show.methods = FALSE, parallel = FALSE, num.cores = 2) {

  weightprob <- split(weightprob, row(weightprob))
  comb.fun <- function(L1, L2){
    L1$weightprob <- L2
    return(list(L1))
  }
  dat <- mapply(comb.fun, dataset, weightprob)
  
  combforec_fun <- function(input){
    lapply(input, function(lentry){
      methodnames <- rownames(lentry$ff)
      
      if (toupper(lentry$period) == "YEARLY"){
        if ("snaive" %in% methodnames){
          prob_hat <- lentry$weightprob[-which(methodnames == "snaive")]
        } else {
          prob_hat <- lentry$weightprob
        }
        threshold <- Threshold[1]
      } else if (toupper(lentry$period) == "QUARTERLY"){
        prob_hat <- lentry$weightprob
        threshold <- Threshold[2]
      } else if (toupper(lentry$period) == "MONTHLY"){
        prob_hat <- lentry$weightprob
        threshold <- Threshold[3]
      }
      r <- prob_hat/max(prob_hat)
      index <- which(r >= threshold) # the index of selected models
      wt <- r[index]
      n <- length(index) # the number of selected models
      ma_models <- names(prob_hat)[index]
      
      ff <- data.frame(lentry$ff)[index, ]
      lower <- data.frame(lentry$lower[[which(names(lentry$lower) == paste(level, "%", sep = ""))]])[index, ]
      upper <- data.frame(lentry$upper[[which(names(lentry$upper) == paste(level, "%", sep = ""))]])[index, ]
      lradius <- as.matrix(ff - lower)
      uradius <- as.matrix(upper - ff)
      
      if(combine == "weighted"){
        ma_ff <- apply(ff, 2, function(x) weighted.mean(x, w=wt))
        ma_lradius <- apply(lradius, 2, function(x) weighted.mean(x, w=wt))
        ma_uradius <- apply(uradius, 2, function(x) weighted.mean(x, w=wt))
        ma_lower <- ma_ff - ma_lradius
        ma_upper <- ma_ff + ma_uradius
      }else if(combine == "median"){
        ma_ff <- apply(ff, 2, median)
        ma_lradius <- apply(lradius, 2, median)
        ma_uradius <- apply(uradius, 2, median)
        ma_lower <- ma_ff - ma_lradius
        ma_upper <- ma_ff + ma_uradius
      }else if(combine == "mean"){
        ma_ff <- apply(ff, 2, mean)
        ma_lradius <- apply(lradius, 2, mean)
        ma_uradius <- apply(uradius, 2, mean)
        ma_lower <- ma_ff - ma_lradius
        ma_upper <- ma_ff + ma_uradius
      }
      
      lentry$ff <- rbind(lentry$ff, ma_ff)
      low <- rbind(lentry$lower[[which(names(lentry$lower) == paste(level, "%", sep = ""))]], ma_lower)
      up <- rbind(lentry$upper[[which(names(lentry$upper) == paste(level, "%", sep = ""))]], ma_upper)
      rownames(lentry$ff) <- rownames(low) <- rownames(up) <- c(methodnames, methodname)
      lentry$lower <- 'names<-' (list(low), paste(level, "%", sep = ""))
      lentry$upper <- 'names<-' (list(up), paste(level, "%", sep = ""))
      
      if (show.methods == TRUE){
        lentry$comb.models <- ma_models
      }
      return(lentry)
    })
  }
  
  if (parallel == FALSE){
    ret_list <- combforec_fun(dat)
  } else {
    ncores <- num.cores
    ncores <- ifelse(length(dat) < ncores, length(dat), ncores)
    times_sp <- c(rep(floor(length(dat)/ncores), ncores - 1), 
                  length(dat) - floor(length(dat)/ncores) * (ncores - 1))
    data_sp <- split(dat, rep(seq_len(ncores), times = times_sp))
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
