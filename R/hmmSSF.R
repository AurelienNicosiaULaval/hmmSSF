#' Fit the HMM-SSF
#'
#' @param ssf_formula Model formula for ssf
#' @param tpm_formula Formula for transition probabilities (default = ~1)
#' @param n_states Number of HMM states
#' @param data Data frame with named columns for ID, stratum, obs, and
#' covariates in ssf_formula
#' @param ssf_par0 Matrix of starting values for SSF parameters, with one row
#' for each covariate and one column for each HMM state.
#' @param tpm_par0 Matrix of starting values for transition probability
#' parameters, with one row for each covariate (including the first row for
#' intercepts) and one column for each off-diagonal transition probability.
#' @param method Optimisation method passed to \code{optim}, defaulting to
#'   \code{"Nelder-Mead"}. See \code{\link{hmmSSF_multistart}} for methods
#'   available with multiple starts.
#' @param optim_opts Additional options to be passed to optim, as a named list
#' (e.g., maxit, trace)
#' @param multistart Logical; use short-run/long-run multi-start optimisation?
#'   The default is \code{FALSE}, which retains the original single-start fit.
#' @param multistart_opts Named list of multi-start settings: \code{n_starts},
#'   \code{n_long}, \code{jitter_sd}, \code{seed}, \code{short_maxit}, and
#'   optionally \code{long_maxit}. See \code{\link{hmmSSF_multistart}}.
#'   Unless overridden, \code{optim_opts$maxit} sets the long-run budget.
#' @param diagnose Logical; calculate and return numerical diagnostics? The
#'   default is \code{FALSE}; calculating a numerical gradient can be slow.
#' @param diagnose_opts Named list of diagnostic settings: \code{gradient}
#'   and \code{eigen_tol}. See \code{\link{diagnose_hmmSSF}}.
#' @return An \code{hmmSSF} object. When requested, \code{multistart} contains
#'   the run history and \code{diagnostics} contains the diagnostic results.
#'   These diagnostics do not establish a global optimum or model validity.
#'
#' @importFrom stats optim
#'
#' @export
#'
#' @useDynLib hmmSSF

hmmSSF <- function(ssf_formula,
                   tpm_formula = ~1,
                   n_states,
                   data,
                   ssf_par0,
                   tpm_par0 = NULL,
                   optim_opts = list(trace = 0, maxit = 5e4),
                   method = "Nelder-Mead",
                   multistart = FALSE,
                   multistart_opts = list(),
                   diagnose = FALSE,
                   diagnose_opts = list()) {

  if (!is.logical(multistart) || length(multistart) != 1L || is.na(multistart))
    stop("'multistart' must be TRUE or FALSE")
  if (!is.logical(diagnose) || length(diagnose) != 1L || is.na(diagnose))
    stop("'diagnose' must be TRUE or FALSE")
  validate_options <- function(x, allowed, name, enabled) {
    if (!is.list(x) || (length(x) &&
        (is.null(names(x)) || anyNA(names(x)) || any(!nzchar(names(x))) ||
         anyDuplicated(names(x)) || any(!names(x) %in% allowed))))
      stop("'", name, "' must be a uniquely named list using supported options")
    if (!enabled && length(x)) stop("'", name, "' requires its feature to be enabled")
  }
  validate_options(multistart_opts,
                   c("n_starts", "n_long", "jitter_sd", "seed", "short_maxit", "long_maxit"),
                   "multistart_opts", multistart)
  validate_options(diagnose_opts, c("gradient", "eigen_tol"), "diagnose_opts", diagnose)

  if (multistart) {
    if (!is.list(optim_opts)) stop("'optim_opts' must be a named list")
    # maxit normally controls the only optimisation. With multiple starts it
    # controls the long runs, unless long_maxit is explicitly supplied.
    if (is.null(multistart_opts$long_maxit) && !is.null(optim_opts$maxit))
      multistart_opts$long_maxit <- optim_opts$maxit
    multi_control <- optim_opts
    multi_control$maxit <- NULL
    mod <- do.call(hmmSSF_multistart,
                   c(list(ssf_formula = ssf_formula, tpm_formula = tpm_formula,
                          n_states = n_states, data = data, ssf_par0 = ssf_par0,
                          tpm_par0 = tpm_par0, method = method,
                          optim_opts = multi_control), multistart_opts))
    if (diagnose)
      mod$diagnostics <- do.call(diagnose_hmmSSF, c(list(mod = mod), diagnose_opts))
    return(mod)
  }

  obs <- subset(data, obs == 1)

  # get ssf model matrix (without intercept)
  options(na.action = 'na.pass')
  ssf_MM <- model.matrix(ssf_formula, data)
  ssf_MM <- ssf_MM[, colnames(ssf_MM) != "(Intercept)", drop = FALSE]

  # get transition probabilities model matrix
  options(na.action = 'na.pass')
  tpm_MM <- model.matrix(tpm_formula, obs)

  # get vector of parameters
  if(is.null(tpm_par0)) {
    # if tpm_par0 not provided, default to -2 for intercept and 0 elsewhere
    tpm_par0 <- matrix(0, nrow = ncol(tpm_MM), ncol = n_states * (n_states - 1))
    tpm_par0[1,] <- - 2
  }
  par <- c(ssf_par0, tpm_par0)

  # get sampling densities for correction
  #sampling_densities <- attr(data, "weights")
  sampling_densities <- data$w

  # optimise negative log likelihood
  fit <- optim(par = par,
               fn = nllk,
               ssf_MM = ssf_MM,
               tpm_MM = tpm_MM,
               sampling_densities = sampling_densities,
               stratum = data$stratum,
               ID = data$ID,
               n_states = n_states,
               n_obs = nrow(obs),
               control = optim_opts,
               hessian = TRUE,
               method = method)

  # unpack fitted parameters
  par <-  format_par(par = fit$par,
                     n_states = n_states,
                     ssf_cov = colnames(ssf_MM),
                     tpm_cov = colnames(tpm_MM))

  # save model formulation in model object
  args <- list(tpm_formula = tpm_formula,
               ssf_formula = ssf_formula,
               data = data,
               n_states = n_states,
               ssf_cov = colnames(ssf_MM),
               tpm_cov = colnames(tpm_MM))

  # returned object
  mod <- list(par = par,
              fit = fit,
              args = args)
  class(mod) <- append("hmmSSF", class(mod))

  if (diagnose)
    mod$diagnostics <- do.call(diagnose_hmmSSF, c(list(mod = mod), diagnose_opts))

  return(mod)
}
