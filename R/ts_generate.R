#' Generate time series with various length by GRATIS.
#'
#' Generate time series with certain seasonal period and length according to 
#' GRATIS proposed by Kang et al. (2019).
#' @param n.ts number of time series to be generated.
#' @param freq seasonal period of time series to be generated. For example, 
#' \code{freq} of the yearly, quarterly, and monthly series are 1, 4, and 12, respectively.
#' @param nComp number of mixing components when simulating time series using MAR models.
#' @param length length vector of time series to be generated. The length of 
#' the generated time series is simply sampled from the length vector.
#' @param h forecasting horizon of the generated time series.
#' @param parallel logical. If TRUE then the generation process is conducted in parallel.
#' @param num.cores the specified amount of parallel processes to be used if parallel = TRUE.
#' @return A list of the generated time series.
#' @examples
#'  x_y <- ts_generate(n.ts = 10, freq = 4, nComp = 2, length = c(20,25,30), h = 8)
#'  x_y$Q1$pars
#'  forecast::autoplot(x_y$Q1$x)
#' @references Kang et al. (2019).
#' @import gratis
#' @export
ts_generate <- function(n.ts = 1, freq = 1, nComp = NULL, length = NULL, 
                        h = 1, parallel = FALSE, num.cores = 2) {
  tsgenerate_fun <- function(input){
    generated.mixture.data <- list()
    for (count in input){
      update_count <- count
      while (count >= update_count) {
        if (is.null(length)){
          n <- 120
        }else{
          n <- sample(length, 1)
        }
        l <- n + h
        if (is.null(nComp)) {
          # set the number of components
          nComp <- sample(1:5, 1)
        } else {
          nComp <- nComp
        }
        pars <- list()
        # seasonal data
        if (freq != 1) {
          if (freq == 4) {
            period <- "QUARTERLY"
          } else if (freq == 12) {
            period <- "MONTHLY"
          } else {
            period <- "OTHER"
          }
          for (i in 1:nComp) {
            pars[[sprintf("pars%d", i)]] <- rnorm(4, 0, 0.5)
          }
          mixture.weights <- rep(NA, nComp)
          for (i in 1:nComp) {
            mixture.weights[i] <- runif(1)
          }
          mixture.weights <- mixture.weights / sum(mixture.weights)
          means.ar.par.list <- lapply(pars, function(x) {
            d <- sample(c(0, 1), 1, prob = c(0.1, 0.9))
            D <- sample(c(0, 1), 1, prob = c(0.6, 0.4))
            c(phi0, pi_coefficients(ar = x[1:2], sar = x[3:4], d = d, D = D, m = freq))
          })
          sigmas.list <- list()
          for (i in 1:nComp) {
            sigmas.list[[i]] <- rep(sigmas[i], l + freq * 10)
          }
          pars$weights <- mixture.weights
          x0 <- rmixnorm_ts(
            n = l + freq * 10,
            means.ar.par.list = means.ar.par.list,
            sigmas.list = sigmas.list,
            weights = mixture.weights
          )
          # allow spikes
          if (runif(1) <= 0.01) {
            x0[sample(1:length(x0), 1)] <- max(x0) * sample(2:5, 1)
          }
          x <- ts(x0[(1 + freq * 10):(n + freq * 10)], frequency = freq)
          xx <- ts(x0[(n + 1 + freq * 10):(l + freq * 10)], frequency = freq)
          if (!(any(is.na(x0)) || (max(abs(x0), na.rm = TRUE) > 1e5))) {
            if (freq == 4) {
              number <- paste0("Q", count)
            } else if (freq == 12) {
              number <- paste0("M", count)
            } else {
              number <- paste0("N", count)
            }
            generated.mixture.data[[number]] <- list()
            generated.mixture.data[[number]]$sn <- number
            generated.mixture.data[[number]]$period<- period
            generated.mixture.data[[number]]$x <- x
            generated.mixture.data[[number]]$xx <- xx
            generated.mixture.data[[number]]$h <- h
            generated.mixture.data[[number]]$n <- n
            generated.mixture.data[[number]]$pars <- pars
            update_count <- update_count + 1
          }
        } else {
          period = "YEARLY"
          # nonseasonal data
          for (i in 1:nComp) {
            pars[[sprintf("pars%d", i)]] <- rnorm(2, 0, 0.5)
          }
          mixture.weights <- rep(NA, nComp)
          for (i in 1:nComp) {
            mixture.weights[i] <- runif(1)
          }
          mixture.weights <- mixture.weights / sum(mixture.weights)
          means.ar.par.list <- lapply(pars, function(x) {
            d <- sample(c(0, 1, 2), 1, prob = c(0.1, 0.6, 0.3))
            c(phi0, pi_coefficients(ar = x[1:2], d = d, m = freq))
          })
          sigmas.list <- list()
          for (i in 1:nComp) {
            sigmas.list[[i]] <- rep(sigmas[i], l + freq * 10)
          }
          pars$weights <- mixture.weights
          x0 <- rmixnorm_ts(
            n = l + freq * 10,
            means.ar.par.list = means.ar.par.list,
            sigmas.list = sigmas.list,
            weights = mixture.weights
          )
          # allow spikes
          if (runif(1) <= 0.01) {
            x0[sample(1:length(x0), 1)] <- max(x0) * sample(2:5, 1)
          }
          x <- ts(x0[(1 + freq * 10):(n + freq * 10)], frequency = freq)
          xx<- ts(x0[(n + 1 + freq * 10):(l + freq * 10)], frequency = freq)
          if (max(abs(x0), na.rm = TRUE) < 1e5) {
            number <- paste0("Y", count)
            generated.mixture.data[[number]] <- list()
            generated.mixture.data[[number]]$sn <- number
            generated.mixture.data[[number]]$period<- period
            generated.mixture.data[[number]]$x <- x
            generated.mixture.data[[number]]$xx <- xx
            generated.mixture.data[[number]]$h <- h
            generated.mixture.data[[number]]$n <- n
            generated.mixture.data[[number]]$pars <- pars
            update_count <- update_count + 1
          }
        }
      }
    }
    return(generated.mixture.data)
  }
  
  index <- seq_len(n.ts)
  sigmas <- sample(c(1:5), 5, replace = TRUE)
  phi0 <- sample(c(0, 3), 1)
  
  if (parallel == FALSE){
    ret_list <- tsgenerate_fun(index)
  } else {
    ncores <- num.cores
    ncores <- ifelse(length(index) < ncores, length(index), ncores)
    times_sp <- c(rep(floor(length(index)/ncores), ncores - 1), 
                  length(index) - floor(length(index)/ncores) * (ncores - 1))
    data_sp <- split(index, rep(seq_len(ncores), times = times_sp))
    
    cl <- parallel::makeCluster(ncores, setup_strategy = "sequential")
    .env <- c("nComp")
    registerDoParallel(cl)
    ret_list <- foreach(data = data_sp, 
                        .combine = base::c, 
                        .export = .env,
                        .packages = 'gratis') %dopar%
      tsgenerate_fun(input = data)
    stopCluster(cl)
  }
  return(ret_list)
}


