
#' Print method for hmmSSF model objects
#' @method print hmmSSF
#'
#' @param x Fitted model object, as returned by \code{\link{hmmSSF}}
#' @param ... Unused argument needed for compatibility with generic S3 method
#'
#' @importFrom stats confint
#'
#' @export
print.hmmSSF <- function(x, ...) {
  mod <- x
  cat("Negative log-likelihood:", mod$fit$value, "\n")
  cat("Convergence code:", mod$fit$convergence, "\n\n")

  if (is.null(mod$fit$hessian)) {
    cat("Confidence intervals unavailable without a Hessian.\n\n")
    par <- mod$par
  } else {
    par <- confint(mod, pretty = TRUE)
  }

  cat("SSF model:\n")
  cat(deparse(mod$args$ssf_formula))
  cat("\n")
  print(round(par$ssf, 3))

  cat("\nTPM model:\n")
  cat(deparse(mod$args$tpm_formula))
  cat("\n")
  if(paste(deparse(mod$args$tpm_formula), collapse = " ") == "~1") {
    tpm <- predict_tpm(mod = mod,
                       new_data = mod$args$data[1,],
                       return_CI = FALSE)
    print(round(tpm$mle[,,1], 3))
  } else {
    print(round(par$tpm, 3))
  }
}
