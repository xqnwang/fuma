#' Adjusted softmax transformation function
#' 
#' In mathematics, softmax function maps a K-dimensional 
#' arbitrary real vector to another K-dimensional probability vector. 
#' Based on the properties of softmax function, the fitted values 
#' obtained from the pre-trained algorithm for each time series 
#' are transformed into probabilities. According to the adjusted 
#' softmax function, the method with a smaller score value 
#' in the candidate methods for each time series has a 
#' larger probability than the other methods. 
#' Consequently, the appropriate methods for each time series 
#' can be picked based on the obtained probabilities.
#' @param x a dataframe with \code{n} rows and \code{F} columns, where \code{n} 
#' is the number of observations and \code{F} is the number of candidate methods. 
#' @export
softmax.fun <- function(msis){
  mu <- apply(msis, 1, mean)
  sigma <- apply(msis, 1, sd)
  V <- (mu - msis)/sigma
  expV <- exp(V)
  sum_exp <- apply(expV, 1, sum)
  stm <- expV/sum_exp
  colnames(stm) <- colnames(msis)
  return(stm)
}