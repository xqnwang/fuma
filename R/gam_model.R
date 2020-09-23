#' Train generalized additive models
#'
#' Fit a generalized additive model(GAM) for data to capture 
#' the relationship between the relationship between the X and each column of Y.
#'
#' @param X a dataframe with \code{n} rows and \code{c} columns, 
#' where \code{n} is the number of observations and 
#' \code{c} is the number of explanatory variables.
#' @param Y a dataframe with \code{n} rows and \code{F} columns. 
#' Each column is considered as a response variable to fit the generalized additive model.
#' @param LogY logical. If \code{TRUE} then the logarithmic form of \code{Y} 
#' becomes the response variable.
#' @param k the dimension of the basis used to represent the smooth term 
#' when \code{mgcv} package is used to fit GAM. 
#' If \code{k} is not specified then basis specific defaults are used.
#' @param parallel logical. If \code{TRUE} then the model training is conducted in parallel.
#' @param num.cores the specified amount of parallel processes to be used if parallel = TRUE.
#'
#' @return \code{gam.fun} returns a list with \code{F} generalized additive models 
#' for each response variable.
#'   
#' @examples
#' X <- as.data.frame(state.x77[, "Murder"])
#' Y <- as.data.frame(state.x77[, c("Population", "Illiteracy", "Income", "Frost")])
#' gam.fit <- gam.fun(X, Y, LogY = FALSE, k = 10)
#' @import mgcv
#' @importFrom parallel detectCores makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom Information is.binary
#' @export
gam.fun <- function(X, Y, LogY = TRUE, k = 10, parallel = FALSE, num.cores = 2){ 
  if (LogY == TRUE){
    Y <- log(Y)
  }else{
    Y <- Y
  }
  if (is.null(k)){
    print("Warning: k need to be pre-defined")
  } else {
    k <- rep(k, ncol(X))
    for (i in 1:ncol(X)){
      if(length(unique(X[, i])) < 4){
        k[i] <- 1
      }
    }
    mgcv_fun <- function(input){
      y <- t(input)
      mgcv.model <- NULL
      for (i in 1:ncol(y)){
        f <- CreateGAMFormula(data = X, y = colnames(y)[i], k = k)
        model_fit <- mgcv::gam(f, data = cbind(X, y), method = "REML", select = TRUE) 
        mgcv.model <- append(mgcv.model, list(model_fit))
      }
      names(mgcv.model) <- colnames(y)
      return(mgcv.model)
    }
    if (parallel == FALSE){
      ret_list <- mgcv_fun(t(Y))
    } else {
      ncores <- num.cores
      ncores <- ifelse(ncol(Y) < ncores, ncol(Y), ncores)
      times_sp <- c(rep(floor(ncol(Y)/ncores), ncores - 1), 
                    ncol(Y) - floor(ncol(Y)/ncores) * (ncores - 1))
      sp <- rep(seq_len(ncores), times = times_sp)
      data_sp <- split.data.frame(t(Y), sp)
      cl <- parallel::makeCluster(ncores)
      .env <- "CreateGAMFormula"
      registerDoParallel(cl)
      ret_list <- foreach(data = data_sp, 
                          .combine = base::c, 
                          .export = .env,
                          .packages = 'mgcv') %dopar%
        mgcv_fun(input = data)
      stopCluster(cl)
    }
    return(ret_list)
  }
}


# Generate fomula for the generalized additive model(GAM)
CreateGAMFormula <- function(data, y, k){
  names <- names(data)
  if (length(names)>0){
    for (i in 1:length(names)){
      if (i==1){
        if (is.factor(data[[names[i]]]) | is.character(data[[names[i]]])){
          Formula <- paste0(y," ~", names[i])     
        } else if (Information::is.binary(data[[names[i]]]) | length(unique(data[[names[i]]]))<4){ 
          Formula <- paste0(y," ~", names[i])     
        } else{
          Formula <- paste0(y," ~ s(", names[i],", bs = 'tp', k = ", k[i], ")") 
        }
      } else{
        if (is.factor(data[[names[i]]]) | is.character(data[[names[i]]])){
          Formula <- paste0(Formula, "+ ",names[i])
        } else if (Information::is.binary(data[[names[i]]]) | length(unique(data[[names[i]]]))<4){ 
          Formula <- paste0(Formula, "+ ",names[i])
        } else{
          Formula <- paste0(Formula, "+ s(",names[i],", bs = 'tp', k = ", k[i], ")")  
          
        }
      }
    }
  } 
  return(as.formula(Formula))
}

