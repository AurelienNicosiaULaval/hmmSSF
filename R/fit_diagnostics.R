#' Diagnose a fitted HMM-SSF
#'
#' Recalculate a numerical gradient of the fitted negative log-likelihood and
#' inspect the Hessian stored in the fitted object. A small gradient and a
#' positive-definite Hessian are local numerical checks, not proof of a global
#' optimum, model adequacy, or identifiability. The gradient depends on the
#' parameterisation and numerical differentiation.
#'
#' @param mod A fitted object returned by \code{\link{hmmSSF}} or
#'   \code{\link{hmmSSF_multistart}}.
#' @param gradient Logical; calculate the numerical gradient? This may be slow.
#' @param eigen_tol Relative tolerance for classifying the Hessian. Eigenvalues
#'   no larger than \code{eigen_tol * max(abs(eigenvalues))} are treated as
#'   numerically non-positive.
#' @return A list containing the optimizer convergence code and message,
#'   objective value, gradient and its maximum absolute element, Hessian
#'   eigenvalues, a Hessian status, and any diagnostic messages. The Hessian
#'   condition number is returned only when the Hessian is positive definite.
#' @export
diagnose_hmmSSF <- function(mod, gradient = TRUE, eigen_tol = sqrt(.Machine$double.eps)) {
  if (!inherits(mod, "hmmSSF")) stop("'mod' must be a fitted hmmSSF object")
  if (!is.logical(gradient) || length(gradient) != 1L || is.na(gradient))
    stop("'gradient' must be TRUE or FALSE")
  if (length(eigen_tol) != 1L || !is.numeric(eigen_tol) ||
      !is.finite(eigen_tol) || eigen_tol < 0)
    stop("'eigen_tol' must be a finite, non-negative number")

  messages <- character()
  grad <- NULL
  if (gradient) {
    grad <- tryCatch({
      inputs <- hmmssf_objective_inputs(mod$args$ssf_formula, mod$args$tpm_formula,
                                        mod$args$data, mod$args$n_states)
      objective <- function(par) hmmssf_checked_nllk(par, inputs, mod$args$data$obs == 1)
      numDeriv::grad(func = objective, x = mod$fit$par)
    }, error = function(e) {
      messages <<- c(messages, paste("Gradient unavailable:", conditionMessage(e)))
      NULL
    })
    if (!is.null(grad) && any(!is.finite(grad)))
      messages <- c(messages, "The numerical gradient contains non-finite values.")
  }
  max_grad <- if (!is.null(grad) && length(grad) > 0L && all(is.finite(grad)))
    max(abs(grad)) else NA_real_

  hessian <- mod$fit$hessian
  eigenvalues <- NULL
  condition <- NA_real_
  status <- "unavailable"
  n_par <- length(mod$fit$par)
  if (is.matrix(hessian) && identical(dim(hessian), c(n_par, n_par)) &&
      all(is.finite(hessian))) {
    eigenvalues <- tryCatch(eigen((hessian + t(hessian)) / 2,
                                  symmetric = TRUE, only.values = TRUE)$values,
                            error = function(e) NULL)
    if (!is.null(eigenvalues) && all(is.finite(eigenvalues))) {
      cutoff <- eigen_tol * max(abs(eigenvalues))
      status <- if (min(eigenvalues) > cutoff) "positive definite" else
        if (min(eigenvalues) < -cutoff) "has negative eigenvalues" else "singular or nearly singular"
      if (status == "positive definite") condition <- max(eigenvalues) / min(eigenvalues)
    }
  }
  if (status == "unavailable") messages <- c(messages, "Hessian diagnostics are unavailable.")
  list(convergence = mod$fit$convergence, optimizer_message = mod$fit$message,
       objective = mod$fit$value, gradient = grad, max_abs_gradient = max_grad,
       hessian_eigenvalues = eigenvalues, hessian_status = status,
       hessian_condition = condition, messages = messages)
}

# Construct the same likelihood inputs as hmmSSF without changing R options.
hmmssf_objective_inputs <- function(ssf_formula, tpm_formula, data, n_states) {
  obs <- subset(data, obs == 1)
  ssf_frame <- stats::model.frame(ssf_formula, data = data, na.action = stats::na.pass)
  tpm_frame <- stats::model.frame(tpm_formula, data = obs, na.action = stats::na.pass)
  ssf_MM <- model.matrix(ssf_formula, ssf_frame)
  ssf_MM <- ssf_MM[, colnames(ssf_MM) != "(Intercept)", drop = FALSE]
  tpm_MM <- model.matrix(tpm_formula, tpm_frame)
  list(ssf_MM = ssf_MM, tpm_MM = tpm_MM, ID = data$ID,
       stratum = data$stratum, n_states = n_states,
       sampling_densities = data$w, n_obs = nrow(obs))
}

# nllk() replaces missing state densities by one. That behavior is useful for
# missing observations but can also hide numerical overflow. Reject overflow
# for finite predictors and unexpected non-finite densities before calling it.
hmmssf_checked_nllk <- function(par, inputs, observed_rows) {
  parameters <- format_par(par, inputs$n_states, colnames(inputs$ssf_MM),
                           colnames(inputs$tpm_MM))
  lp <- inputs$ssf_MM %*% parameters$ssf
  log_limit <- log(.Machine$double.xmax)
  weighted_lp <- lp - log(inputs$sampling_densities)
  if (any(is.infinite(lp)) || any(is.infinite(weighted_lp)) ||
      any(lp[is.finite(lp)] > log_limit) ||
      any(weighted_lp[is.finite(weighted_lp)] > log_limit))
    hmmssf_numerical_error("Numerical overflow in state-dependent densities")
  densities <- state_dens_rcpp(lp, inputs$stratum, inputs$n_states,
                                inputs$sampling_densities, inputs$n_obs)
  observed_lp <- lp[observed_rows, , drop = FALSE]
  if (any(is.finite(observed_lp) & (!is.finite(densities) | densities <= 0)))
    hmmssf_numerical_error("State-dependent density is non-positive or non-finite")
  hmmssf_forward_nllk(densities, parameters$tpm, inputs$tpm_MM, inputs$ID,
                     inputs$stratum, inputs$n_states, check = TRUE)
}

# Mark numerical failures so optimisers may reject a trial point without
# hiding unrelated programming or input errors.
hmmssf_numerical_error <- function(message) {
  stop(structure(list(message = message, call = NULL),
                 class = c("hmmssf_numerical_error", "error", "condition")))
}

# A non-finite trial value is outside the numerically evaluable region. The
# starting point is checked separately, and unrelated errors still propagate.
hmmssf_safe_objective <- function(fn, counter) {
  force(fn)
  force(counter)
  function(par) {
    value <- tryCatch(fn(par), hmmssf_numerical_error = function(e) NA_real_)
    if (length(value) != 1L || !is.finite(value)) {
      counter$rejected <- counter$rejected + 1L
      return(Inf)
    }
    value
  }
}
