#' Search the optimal threshold ratio used for model averaging
#' 
#' According to the fitted values of the pre-trained algorithm, 
#' the optimal threshold ratio used for model averaging is searched 
#' from a set of pre-set threshold ratios.
#' 
#' @param dataset a list of time series. See details for the required format.
#' @param fitProb a matrix with \code{n} (number of time series) rows and 
#' \code{F} (number of benchmark methods) matrix. The matrix is applied to assign weights 
#' to each benchmark method for each time series.
#' @param ratio a vector with pre-set threshold ratios. The threshold ratio 
#' quantifies as a number between 0 and 1, where 0 indicates all the benchmark methods 
#' are selected and 1 indicates only the benchmark method with the minimal value of 
#' scores is selected.
#' @param combine string representing the interval combination methods. 
#' "mean", "weightedmean" or "median".
#' @param level the confidence level of prediction intervals that are used for combination.
#' @param parallel logical. If \code{TRUE} then the calculations are conducted in parallel.
#' @param num.cores the specified amount of parallel processes to be used if parallel = TRUE.
#' 
#' @details 
#' \code{dataset} must be a list with each element having the following format:
#' \describe{
#'   \item{n}{the number of observations in the time series.}
#'   \item{h}{the number of required forecasts.}
#'   \item{period}{interval of the time series. Possible values are "YEARLY", "QUARTERLY", "MONTHLY" & "OTHER".}
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
#' @return A list with the following elements
#' \describe{
#'   \item{threshold_y}{numeric, the optimal threshold ratio searched for yearly data.}
#'   \item{threshold_q}{numeric, the optimal threshold ratio searched for quarterly data.}
#'   \item{threshold_m}{numeric, the optimal threshold ratio searched for monthly data.}
#'   \item{msis_y}{a matrix of MSIS values for each pre-set threshold ratio to yearly data.}
#'   \item{msis_q}{a matrix of MSIS values for each pre-set threshold ratio to quarterly data.}
#'   \item{msis_m}{a matrix of MSIS values for each pre-set threshold ratio to monthly data.}
#' }
#' @import magrittr
#' @import stats
#' @importFrom parallel detectCores makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @export
choose.threshold <- function(dataset, fitProb, ratio, 
                             combine = "mean", level, 
                             parallel = FALSE, num.cores = 2) {
  
  ret_list <- lapply(ratio, function(lentry){
    r <- rep(lentry,3)
    combforec <- comb_forec(dataset, weightprob = fitProb, 
                            Threshold = r, level = level, 
                            combine = combine, methodname = "comb", 
                            show.methods = FALSE, parallel = parallel, num.cores = num.cores)
    combscore <- calc_scores(combforec, parallel = parallel, num.cores = num.cores)
    msis_y <- lapply(Filter(function(l) toupper(l$period) == "YEARLY", combscore), function(lentry){
     scores <- lentry$MSIS[[which(names(lentry$MSIS) == paste(level, "%", sep = ""))]]
     return(scores)
    })
    msis_q <- lapply(Filter(function(l) toupper(l$period) == "QUARTERLY", combscore), function(lentry){
      scores <- lentry$MSIS[[which(names(lentry$MSIS) == paste(level, "%", sep = ""))]]
      return(scores)
    })
    msis_m <- lapply(Filter(function(l) toupper(l$period) == "MONTHLY", combscore), function(lentry){
      scores <- lentry$MSIS[[which(names(lentry$MSIS) == paste(level, "%", sep = ""))]]
      return(scores)
    })
    avg_msis_y <- Reduce("+", msis_y)/length(msis_y)
    avg_msis_q <- Reduce("+", msis_q)/length(msis_q)
    avg_msis_m <- Reduce("+", msis_m)/length(msis_m)
    comb_msis_y <- avg_msis_y[which(rownames(avg_msis_y) == "comb"), ] %>% t() %>% data.frame() 
    comb_msis_q <- avg_msis_q[which(rownames(avg_msis_q) == "comb"), ] %>% t() %>% data.frame()
    comb_msis_m <- avg_msis_m[which(rownames(avg_msis_m) == "comb"), ] %>% t() %>% data.frame()
    return(list(msis_y = comb_msis_y, msis_q = comb_msis_q, msis_m = comb_msis_m))
  })
  msis_y <- `rownames<-` (sapply(ret_list, function(lentry){as.numeric(lentry$msis_y)}) 
                          %>% t() %>% data.frame(), paste("r_", ratio, sep = ""))
  msis_q <- `rownames<-` (sapply(ret_list, function(lentry){as.numeric(lentry$msis_q)})
                          %>% t() %>% data.frame(), paste("r_", ratio, sep = ""))
  msis_m <- `rownames<-` (sapply(ret_list, function(lentry){as.numeric(lentry$msis_m)})
                          %>% t() %>% data.frame(), paste("r_", ratio, sep = ""))
  
  threshold_y <- ratio[which.min(apply(msis_y, 1, mean))]
  threshold_q <- ratio[which.min(apply(msis_q, 1, mean))]
  threshold_m <- ratio[which.min(apply(msis_m, 1, mean))]
  
  gen.list <- list()
  gen.list[["threshold_y"]] <- threshold_y
  gen.list[["threshold_q"]] <- threshold_q
  gen.list[["threshold_m"]] <- threshold_m
  gen.list[["msis_y"]] <- msis_y
  gen.list[["msis_q"]] <- msis_q
  gen.list[["msis_m"]] <- msis_m
  return(gen.list)
}
