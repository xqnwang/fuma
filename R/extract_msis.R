#' Extract MSIS values for certain confidence level from a list
#'
#' Extract MSIS values for certain confidence level from a list 
#' that each element contains a list named as \code{MSIS}.
#' @param dataset A list with each element having a \code{MSIS} fields.
#' @param level the confidence level for which MSIS values are extracted.
#'
#' @return A matrix with \code{N} rows and \code{F} columns. Each row contains
#'   the MSIS values of each method in \code{methods}.
#'   
#' @importFrom utils tail 
#' @export
extract_msis <- function(dataset, level) {
  extracted <- t(sapply(dataset, function (lentry) {
    stopifnot("MSIS" %in% names(lentry))
    msis <- lentry$MSIS[[which(names(lentry$MSIS) == paste(level, sep = "", "%"))]]
    msis_h <- as.vector(utils::tail(t(msis),1))
    names(msis_h) <- rownames(msis)
    return(msis_h)
  }))
  return(extracted)
}
