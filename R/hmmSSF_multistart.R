#' Fit an HMM-SSF using short and long optimisation runs
#'
#' Generate perturbations of the supplied starting values on the working
#' parameter scale. Run a short optimisation from each start, then continue
#' the best finite short runs with longer optimisations. All finite short and
#' long results remain eligible for selection. The smallest objective is
#' returned even if its convergence code is nonzero; inspect the run summaries
#' and use \code{\link{diagnose_hmmSSF}} before interpreting the result. A
#' warning is issued if the selected run has a nonzero convergence code.
#' Detectable numerical overflows at trial points are rejected and counted;
#' an invalid initial point causes that run to fail. This does not detect
#' every possible numerical problem.
#' Multiple starts improve the search but do not prove a global optimum.
#'
#' @param ssf_formula,tpm_formula,n_states,data,ssf_par0,tpm_par0 Model inputs
#'   as in \code{\link{hmmSSF}}. Every track must contain at least two
#'   observations and input rows must be ordered by contiguous track and
#'   stratum. Each stratum must begin with the observed row and contain at least
#'   one control. Sampling weights must be finite and positive.
#' @param n_starts Number of short runs, including the supplied start.
#' @param n_long Maximum number of short-run endpoints continued in long runs.
#' @param jitter_sd Standard deviation of independent normal perturbations on
#'   the working scale. A scalar or one value per fitted parameter. The first
#'   start is always the supplied parameter vector, without perturbation.
#' @param seed Integer seed for the starting values. The caller's
#'   \code{.Random.seed} is restored on exit.
#' @param short_maxit,long_maxit Iteration budgets for the short and long runs.
#' @param method Optimisation method: \code{"Nelder-Mead"}, \code{"BFGS"}, or
#'   \code{"CG"}.
#' @param optim_opts Named control list passed to \code{optim}; use the two
#'   iteration-budget arguments instead of supplying \code{maxit} here.
#' @return An \code{hmmSSF} object, with a \code{multistart} element containing
#'   the random seed, all initial values, every optimisation result, a run
#'   summary with rejected trial evaluations, and the selected stage and start.
#'   The Hessian is calculated for the selected result and may be unavailable.
#' @export
hmmSSF_multistart <- function(ssf_formula, tpm_formula = ~1, n_states, data,
                              ssf_par0, tpm_par0 = NULL, n_starts = 10L,
                              n_long = 3L, jitter_sd = 0.5, seed = 1L,
                              short_maxit = 100L, long_maxit = 5000L,
                              method = "Nelder-Mead", optim_opts = list()) {
  positive_integer <- function(x) is.numeric(x) && length(x) == 1L &&
    is.finite(x) && x >= 1 && x == floor(x) && x <= .Machine$integer.max
  if (!positive_integer(n_states) || n_states < 2L) stop("'n_states' must be an integer >= 2")
  for (name in c("n_starts", "n_long", "short_maxit", "long_maxit")) {
    if (!positive_integer(get(name))) stop(paste0("'", name, "' must be a positive integer"))
  }
  if (n_long > n_starts) stop("'n_long' cannot exceed 'n_starts'")
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != floor(seed) || abs(seed) > .Machine$integer.max)
    stop("'seed' must be an integer representable by R")
  if (length(method) != 1L || is.na(method) || !method %in% c("Nelder-Mead", "BFGS", "CG"))
    stop("'method' must be 'Nelder-Mead', 'BFGS', or 'CG'")
  if (!is.list(optim_opts) || (length(optim_opts) &&
      (is.null(names(optim_opts)) || any(!nzchar(names(optim_opts))) || anyDuplicated(names(optim_opts)))))
    stop("'optim_opts' must be a uniquely named list")
  if ("maxit" %in% names(optim_opts)) stop("Use 'short_maxit' and 'long_maxit' instead of optim_opts$maxit")
  if ("fnscale" %in% names(optim_opts) &&
      (length(optim_opts$fnscale) != 1L || !is.finite(optim_opts$fnscale) || optim_opts$fnscale <= 0))
    stop("optim_opts$fnscale must be positive and finite")
  if (!is.data.frame(data) || !all(c("ID", "stratum", "obs", "w") %in% names(data)))
    stop("'data' must contain ID, stratum, obs, and w columns")
  if (anyNA(data$ID) || anyNA(data$stratum) || anyNA(data$obs) ||
      !(is.numeric(data$obs) || is.logical(data$obs)) ||
      !is.numeric(data$w) || any(!is.finite(data$w) | data$w <= 0) ||
      any(!data$obs %in% c(0, 1))) stop("Invalid ID, stratum, obs, or sampling weights")
  if (!(is.factor(data$stratum) ||
        (is.numeric(data$stratum) && all(is.finite(data$stratum)))))
    stop("'stratum' must be a factor or a finite numeric vector")
  if (nrow(data) == 0L) stop("'data' must not be empty")
  if (any(lengths(data[c("ID", "stratum", "obs", "w")]) != nrow(data)))
    stop("ID, stratum, obs, and w must have one value per row")
  # The upstream likelihood assumes consecutive strata and tracks.
  strat <- rle(as.numeric(data$stratum))
  track <- match(data$ID, unique(data$ID))
  ids <- rle(track)
  if (anyDuplicated(strat$values) || anyDuplicated(ids$values))
    stop("Each stratum and each track must form a single contiguous block")
  if (any(strat$lengths < 2L)) stop("Each stratum must contain at least one control")
  last <- cumsum(strat$lengths)
  first <- last - strat$lengths + 1L
  if (any(data$obs[first] != 1) || sum(data$obs) != length(first) ||
      any(track[first] != track[last]))
    stop("Each stratum must have exactly one observed row, listed first, and one track ID")
  if (any(tabulate(track[first]) < 2L))
    stop("Each track must contain at least two observed strata for the upstream likelihood")

  inputs <- hmmssf_objective_inputs(ssf_formula, tpm_formula, data, n_states)
  objective <- function(par) hmmssf_checked_nllk(par, inputs, data$obs == 1)
  n_ssf <- ncol(inputs$ssf_MM)
  n_tpm <- ncol(inputs$tpm_MM)
  if (n_ssf < 1L) stop("The SSF formula must contain at least one non-intercept term")
  if (n_tpm < 1L || any(!is.finite(inputs$tpm_MM)))
    stop("The transition formula must produce at least one column with finite values")
  if (is.null(tpm_par0)) {
    tpm_par0 <- matrix(0, nrow = n_tpm, ncol = n_states * (n_states - 1L))
    tpm_par0[colnames(inputs$tpm_MM) == "(Intercept)", ] <- -2
  }
  if (!is.matrix(ssf_par0) || !is.numeric(ssf_par0) ||
      !identical(dim(ssf_par0), c(n_ssf, as.integer(n_states))) || any(!is.finite(ssf_par0)))
    stop("'ssf_par0' has incorrect dimensions or non-finite values")
  if (!is.matrix(tpm_par0) || !is.numeric(tpm_par0) ||
      !identical(dim(tpm_par0), c(n_tpm, as.integer(n_states * (n_states - 1L)))) || any(!is.finite(tpm_par0)))
    stop("'tpm_par0' has incorrect dimensions or non-finite values")
  initial <- c(ssf_par0, tpm_par0)
  if (!is.numeric(jitter_sd) || !length(jitter_sd) ||
      !length(jitter_sd) %in% c(1L, length(initial)) || any(!is.finite(jitter_sd) | jitter_sd < 0))
    stop("'jitter_sd' must be non-negative and have length 1 or the number of parameters")
  jitter_sd <- rep(jitter_sd, length.out = length(initial))

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  starts <- lapply(seq_len(n_starts), function(i) {
    if (i == 1L) initial else initial + stats::rnorm(length(initial), sd = jitter_sd)
  })
  run_one <- function(par, budget) {
    warnings <- character()
    counter <- new.env(parent = emptyenv())
    counter$rejected <- 0L
    safe_objective <- hmmssf_safe_objective(objective, counter)
    result <- tryCatch(withCallingHandlers({
      initial_value <- objective(par)
      if (length(initial_value) != 1L || !is.finite(initial_value))
        stop("Non-finite objective at the starting values")
      stats::optim(par = par, fn = safe_objective, method = method,
                   control = c(optim_opts, list(maxit = budget)))
    },
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) list(error = conditionMessage(e)))
    result$warnings <- warnings
    result$rejected_evaluations <- counter$rejected
    result
  }
  shorts <- lapply(starts, run_one, budget = short_maxit)
  value_of <- function(x) if (is.null(x$error) && length(x$value) == 1L &&
                              is.finite(x$value) && all(is.finite(x$par))) x$value else Inf
  short_values <- vapply(shorts, value_of, numeric(1))
  eligible <- which(is.finite(short_values))
  if (!length(eligible)) {
    errors <- vapply(shorts, function(x) if (is.null(x$error)) "non-finite objective" else x$error, character(1))
    stop("All short runs failed or returned non-finite objectives. First failure: ", errors[1])
  }
  selected_short <- utils::head(eligible[order(short_values[eligible], eligible)], n_long)
  if (length(selected_short) < n_long)
    warning("Fewer finite short runs than requested long runs; continuing ", length(selected_short), " run(s)")
  longs <- lapply(shorts[selected_short], function(x) run_one(x$par, long_maxit))
  long_values <- vapply(longs, value_of, numeric(1))
  if (!any(is.finite(long_values))) warning("All long runs failed; returning the best finite short run")
  all_values <- c(short_values, long_values)
  # Break exact ties in favour of converged and then long runs.
  all_runs <- c(shorts, longs)
  converged <- vapply(all_runs, function(x)
    length(x$convergence) == 1L && !is.na(x$convergence) && x$convergence == 0L, logical(1))
  is_long <- c(rep(FALSE, length(shorts)), rep(TRUE, length(longs)))
  best <- order(all_values, !converged, !is_long, seq_along(all_values))[1L]
  best_stage <- if (best <= length(shorts)) "short" else "long"
  best_id <- if (best_stage == "short") best else selected_short[best - length(shorts)]
  fit <- if (best_stage == "short") shorts[[best]] else longs[[best - length(shorts)]]
  if (best_stage == "short") warning("The selected result is a short run; inspect the long-run results")
  if (is.null(fit$convergence) || fit$convergence != 0L)
    warning("The selected optimisation did not report convergence; inspect the run summary and diagnostics")
  hessian_error <- NULL
  fit$hessian <- tryCatch(stats::optimHess(par = fit$par, fn = objective, control = optim_opts),
                          error = function(e) { hessian_error <<- conditionMessage(e); NULL })
  if (!is.null(fit$hessian) && any(!is.finite(fit$hessian))) {
    hessian_error <- "Non-finite Hessian"
    fit$hessian <- NULL
  }
  if (is.null(fit$hessian)) warning("Hessian calculation failed; diagnostic and interval estimates may be unavailable")
  summarise <- function(results, stage, ids) {
    data.frame(stage = stage, start = ids,
               objective = vapply(results, value_of, numeric(1)),
               convergence = vapply(results, function(x) if (is.null(x$convergence)) NA_integer_ else as.integer(x$convergence), integer(1)),
               rejected_evaluations = vapply(results, function(x) x$rejected_evaluations, integer(1)),
               message = vapply(results, function(x) if (!is.null(x$error)) x$error else if (!is.null(x$message)) x$message else "", character(1)),
               warnings = vapply(results, function(x) paste(x$warnings, collapse = "; "), character(1)),
               stringsAsFactors = FALSE)
  }
  summary <- rbind(summarise(shorts, "short", seq_along(shorts)),
                   summarise(longs, "long", selected_short))
  args <- list(tpm_formula = tpm_formula, ssf_formula = ssf_formula, data = data,
               n_states = n_states, ssf_cov = colnames(inputs$ssf_MM),
               tpm_cov = colnames(inputs$tpm_MM))
  pars <- format_par(fit$par, n_states, args$ssf_cov, args$tpm_cov)
  mod <- list(par = pars, fit = fit, args = args,
              multistart = list(seed = seed, jitter_sd = jitter_sd, starts = starts,
                                short_runs = shorts, long_runs = longs,
                                summary = summary, selected_stage = best_stage,
                                selected_start = best_id, method = method,
                                short_maxit = short_maxit, long_maxit = long_maxit,
                                optim_opts = optim_opts, hessian_error = hessian_error))
  class(mod) <- c("hmmSSF", "list")
  mod
}
