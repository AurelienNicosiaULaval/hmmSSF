#' Predict SSF for new data
#'
#' @param mod Fitted model returned by \code{\link{hmmSSF}}
#' @param new_data Data frame with covariate values used for prediction
#'
#' @return Matrix of predictions in each state
#'
#' @export

predict_ssf <- function(mod, new_data) {
  ssf_MM <- model.matrix(mod$args$ssf_formula, new_data)
  ssf_MM <- ssf_MM[, colnames(ssf_MM) != "(Intercept)", drop = FALSE]
  pred <- exp(ssf_MM %*% mod$par$ssf)
  return(pred)
}
