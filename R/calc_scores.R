#' Calculate various scoring rules for a Time Series Dataset
#'
#' For each series in \code{dataset}, point and interval forecasting performance 
#' are evaluated for all methods in \code{methods} of \code{\link{ts_forec}} 
#' with regard to MASE, sMAPE, MSIS, Spread, Coverage and Upper coverage. 
#' @param dataset the list containing the series. See details for the required format.
#' @param parallel logical. If \code{TRUE} then the calculations are conducted in parallel.
#' @param num.cores the specified amount of parallel processes to be used if parallel = TRUE.
#' 
#' @details 
#' \code{dataset} must be a list with each element having the following format:
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
#' @return A list with the elements having the following structure
#' \describe{
#'   \item{MASE}{a matrix with \code{F} rows and \code{h} columns. Each row contains 
#'   the MASE scores of each method in \code{methods}.
#'   Each column represents the MASE for the j-step point forecasts, where j=1,...,h. }
#'   \item{sMAPE}{a matrix with \code{F} rows and \code{h} columns. Each element 
#'   in the matrix denotes the sMAPE values.}
#'   \item{MSIS}{a list with each element being a matrix for certain confidence level in 
#'   \code{level}. The matrix contains the MSIS values in corresponding forecasting 
#'   horizon for each method.}
#'   \item{Spread}{a list with each element being a matrix for certain confidence level in 
#'   \code{level}. The matrix contains the Spread values in corresponding forecasting 
#'   horizon for each method.}
#'   \item{IfInn}{a list with each element being a matrix for certain confidence level in 
#'   \code{level}. Each element in the matrix measures if the true value lies inside 
#'   the prediction interval.}
#'   \item{IfInu}{a list with each element being a matrix for certain confidence level in 
#'   \code{level}. Each element in the matrix measures if the true value is not larger than 
#'   the upper bound of the prediction interval.}
#' }
#' 
#' @importFrom parallel detectCores makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom stats frequency
#' @export
calc_scores <- function(dataset, parallel = FALSE, num.cores = 2) {
  
  scores_fun <- function(input){
    lapply(input, function(lentry){
      n_methods <- nrow(lentry$ff)
      ff <- lentry$ff
      xx <- 'rownames<-' (matrix(rep(lentry$xx, n_methods), 
                                 nrow = n_methods, byrow = TRUE), 
                          rownames(ff))
      lower <- lentry$lower
      upper <- lentry$upper
      frequency <- stats::frequency(lentry$x)
      scaling <- mean(abs(diff(as.vector(lentry$x), frequency)))
      
      # point forecasting accuracy
      ## MASE
      mase <- abs(xx - ff) / scaling
      ## SMAPE
      smape <- (abs(xx-ff)*200)/(abs(xx)+abs(ff))
      
      # interval forecasting accuracy
      ## MSIS; Spread; Coverage; UpperCoverage
      msis <- spread <- IfInn <- IfInu <- IfInl <- NULL
      for (i in seq_len(length(lower))){
        low <- lower[[i]]
        up <- upper[[i]]
        alpha <- (100 - as.numeric(gsub("%", "", names(lower)[i])))/100
        msis0 <- (up - low +
                    (2 / alpha) * (low - xx) * (low > xx) +
                    (2 / alpha) * (xx - up) * (up < xx)) / scaling
        spread0 <- (up - low) / scaling
        IfInn0 <- ifelse(xx > low & xx < up, 1, 0)
        IfInu0 <- ifelse(xx < up, 1, 0)
        IfInl0 <- ifelse(xx > low, 1, 0)
        msis <- append(msis, list(msis0))
        spread <- append(spread, list(spread0))
        IfInn <- append(IfInn, list(IfInn0))
        IfInu <- append(IfInu, list(IfInu0))
        IfInl <- append(IfInl, list(IfInl0))
      }
      names(msis) <- names(spread) <- names(IfInn) <- names(IfInu) <- names(IfInl) <- names(lower)
      lentry$MASE <- mase
      lentry$sMAPE <- smape
      lentry$MSIS <- msis
      lentry$Spread <- spread
      lentry$IfInn <- IfInn
      lentry$IfInu <- IfInu
      lentry$IfInl <- IfInl
      return(lentry)
    })
  }
  
  if (parallel == FALSE){
    ret_list <- scores_fun(dataset)
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
      scores_fun(input = data)
    stopCluster(cl)
  }
  ret_list
}

