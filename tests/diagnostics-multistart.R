library(hmmSSF)

# Small deterministic choice sets; these are not simulated animal tracks.
set.seed(12)
n <- 12L
observed <- data.frame(ID = 1, stratum = seq_len(n), obs = 1,
                       x = rnorm(n), y = rnorm(n), w = 1)
controls <- transform(observed, obs = 0, x = 0, y = 0)
dat <- rbind(observed, controls)
dat <- dat[order(dat$stratum, -dat$obs), ]
ssf_start <- matrix(c(1, 0, -1, 0), nrow = 2)

fit_once <- function() {
  hmmSSF_multistart(~x + y, n_states = 2, data = dat, ssf_par0 = ssf_start,
                    n_starts = 4, n_long = 2, seed = 42,
                    short_maxit = 15, long_maxit = 35)
}
set.seed(239)
seed_before <- .Random.seed
warned <- FALSE
fit <- withCallingHandlers(fit_once(), warning = function(w) {
  if (grepl("did not report convergence", conditionMessage(w), fixed = TRUE)) warned <<- TRUE
  invokeRestart("muffleWarning")
})
stopifnot(warned)
stopifnot(identical(seed_before, .Random.seed))
fit2 <- suppressWarnings(fit_once())
stopifnot(identical(fit$multistart$starts, fit2$multistart$starts),
          identical(fit$multistart$summary, fit2$multistart$summary),
          inherits(fit, "hmmSSF"), length(fit$multistart$short_runs) == 4L,
          length(fit$multistart$long_runs) == 2L,
          nrow(fit$multistart$summary) == 6L,
          identical(fit$multistart$starts[[1]], c(ssf_start, -2, -2)))

summary <- fit$multistart$summary
short <- summary[summary$stage == "short", ]
long <- summary[summary$stage == "long", ]
short_order <- order(short$objective, short$start)
stopifnot(identical(long$start, short$start[short_order[1:2]]),
          identical(fit$fit$value, min(summary$objective)),
          fit$multistart$selected_start %in% seq_len(4))

# Verify that the returned parameter vector is evaluated using the package's
# original likelihood and that the gradient matches central differences.
ssf_mm <- model.matrix(~x + y, dat)[, -1, drop = FALSE]
tpm_mm <- model.matrix(~1, dat[dat$obs == 1, , drop = FALSE])
objective <- function(par) nllk(par, ssf_MM = ssf_mm, tpm_MM = tpm_mm,
                                 ID = dat$ID, stratum = dat$stratum,
                                 n_states = 2, sampling_densities = dat$w,
                                 n_obs = n)
stopifnot(isTRUE(all.equal(unname(objective(fit$fit$par)), fit$fit$value)))
# Re-run one continuation independently from the recorded short-run endpoint.
manual_long <- optim(fit$multistart$short_runs[[long$start[1]]]$par,
                     objective, method = "Nelder-Mead", control = list(maxit = 35))
stopifnot(isTRUE(all.equal(manual_long$value, fit$multistart$long_runs[[1]]$value)),
          isTRUE(all.equal(manual_long$par, fit$multistart$long_runs[[1]]$par)))
diagnostic <- diagnose_hmmSSF(fit)
h <- 1e-5
numerical <- vapply(seq_along(fit$fit$par), function(i) {
  a <- b <- fit$fit$par
  a[i] <- a[i] + h
  b[i] <- b[i] - h
  (objective(a) - objective(b)) / (2 * h)
}, numeric(1))
stopifnot(isTRUE(all.equal(unname(diagnostic$gradient), numerical, tolerance = 1e-5)),
          diagnostic$convergence == fit$fit$convergence,
          isTRUE(all.equal(diagnostic$objective, fit$fit$value)))

# Known curvature cases and unavailable Hessian.
test_fit <- fit
np <- length(fit$fit$par)
test_fit$fit$hessian <- diag(seq_len(np))
d <- diagnose_hmmSSF(test_fit, gradient = FALSE)
stopifnot(d$hessian_status == "positive definite", d$hessian_condition == np,
          is.null(d$gradient), is.na(d$max_abs_gradient))
test_fit$fit$hessian <- diag(c(-1, rep(1, np - 1)))
stopifnot(diagnose_hmmSSF(test_fit, gradient = FALSE)$hessian_status == "has negative eigenvalues")
test_fit$fit$hessian <- diag(c(0, rep(1, np - 1)))
stopifnot(diagnose_hmmSSF(test_fit, gradient = FALSE)$hessian_status == "singular or nearly singular")
test_fit$fit$hessian <- NULL
stopifnot(diagnose_hmmSSF(test_fit, gradient = FALSE)$hessian_status == "unavailable")

# Preserve the caller's RNG even when all runs fail.
set.seed(301)
rng <- .Random.seed
err <- try(hmmSSF_multistart(~x + y, n_states = 2, data = dat, ssf_par0 = ssf_start,
                             n_starts = 2, n_long = 1, method = "BFGS", seed = 42,
                             short_maxit = 3, long_maxit = 3,
                             optim_opts = list(ndeps = 0)), silent = TRUE)
stopifnot(inherits(err, "try-error"), identical(rng, .Random.seed))

# Reject malformed inputs before running an optimisation.
expect_error <- function(expr) stopifnot(inherits(try(force(expr), silent = TRUE), "try-error"))
expect_error(hmmSSF_multistart(~x + y, n_states = 2, data = dat, ssf_par0 = ssf_start,
                              n_starts = 1, n_long = 2))
expect_error(hmmSSF_multistart(~x + y, n_states = 2, data = dat, ssf_par0 = ssf_start,
                              optim_opts = list(maxit = 5)))
expect_error(hmmSSF_multistart(~x + y, n_states = 2, data = dat, ssf_par0 = ssf_start,
                              jitter_sd = c(1, 2)))
invalid <- dat
invalid$obs[1] <- 0
expect_error(hmmSSF_multistart(~x + y, n_states = 2, data = invalid, ssf_par0 = ssf_start))

# Diagnostics also accept objects returned by the original single-start fitter.
single <- hmmSSF(~x + y, n_states = 2, data = dat, ssf_par0 = ssf_start,
                 optim_opts = list(maxit = 5))
stopifnot(is.finite(diagnose_hmmSSF(single)$max_abs_gradient))

# Extreme perturbations provoke failed short runs; they are recorded and only
# finite candidates proceed, rather than aborting all remaining runs.
partial <- suppressWarnings(hmmSSF_multistart(~x + y, n_states = 2, data = dat,
  ssf_par0 = ssf_start, n_starts = 5, n_long = 2, seed = 42, jitter_sd = 1000,
  short_maxit = 5, long_maxit = 5))
ps <- partial$multistart$summary
stopifnot(any(ps$stage == "short" & !is.finite(ps$objective)),
          any(nzchar(ps$message[ps$stage == "short" & !is.finite(ps$objective)])),
          all(is.finite(ps$objective[ps$stage == "long"])))

# One SSF covariate, two tracks and a factor transition covariate exercise
# model-matrix dimensions and factor expansion.
dat2 <- dat
dat2$ID <- ifelse(dat2$stratum <= 6, "A", "B")
dat2$group <- factor(ifelse(dat2$stratum %% 2 == 0, "night", "day"))
one_cov <- suppressWarnings(hmmSSF_multistart(~x, tpm_formula = ~group,
  n_states = 2, data = dat2, ssf_par0 = matrix(c(0.5, -0.5), nrow = 1),
  n_starts = 2, n_long = 1, seed = 6, short_maxit = 5, long_maxit = 5))
stopifnot(is.finite(one_cov$fit$value),
          nrow(one_cov$par$ssf) == 1L, nrow(one_cov$par$tpm) == 2L,
          is.finite(diagnose_hmmSSF(one_cov)$max_abs_gradient))

# Regression: the upstream nllk replaces overflow-induced NaN densities with 1
# and can return a finite objective. The multi-start fit must reject that run.
extreme <- dat
extreme$x <- as.numeric(extreme$obs == 1)
extreme$y <- as.numeric(extreme$obs == 1)
overflow <- try(hmmSSF_multistart(~x + y, n_states = 2, data = extreme,
  ssf_par0 = matrix(1000, nrow = 2, ncol = 2), n_starts = 1, n_long = 1,
  jitter_sd = 0, short_maxit = 2, long_maxit = 2), silent = TRUE)
stopifnot(inherits(overflow, "try-error"),
          grepl("overflow", as.character(overflow), fixed = TRUE))

# Main interface: defaults remain single-start and preserve the existing fields.
stopifnot(identical(names(single), c("par", "fit", "args")))
single_diag <- hmmSSF(~x + y, n_states = 2, data = dat, ssf_par0 = ssf_start,
                      optim_opts = list(maxit = 5), diagnose = TRUE,
                      diagnose_opts = list(gradient = FALSE))
stopifnot(identical(single$fit, single_diag$fit),
          identical(single$par, single_diag$par),
          !is.null(single_diag$diagnostics),
          is.null(single_diag$diagnostics$gradient))

# Both options work together, and optim_opts$maxit is the long-run budget.
main_fit <- suppressWarnings(hmmSSF(~x + y, n_states = 2, data = dat,
  ssf_par0 = ssf_start, optim_opts = list(maxit = 35),
  multistart = TRUE, diagnose = TRUE,
  multistart_opts = list(n_starts = 4, n_long = 2, seed = 42,
                         short_maxit = 15)))
stopifnot(identical(main_fit$multistart$starts, fit$multistart$starts),
          isTRUE(all.equal(main_fit$fit$value, fit$fit$value)),
          main_fit$multistart$long_maxit == 35,
          isTRUE(all.equal(main_fit$diagnostics$max_abs_gradient,
                           diagnose_hmmSSF(main_fit)$max_abs_gradient)))

# Explicit long_maxit takes precedence; invalid options cannot be ignored.
override <- suppressWarnings(hmmSSF(~x + y, n_states = 2, data = dat,
  ssf_par0 = ssf_start, optim_opts = list(maxit = 35), multistart = TRUE,
  multistart_opts = list(n_starts = 2, n_long = 1, short_maxit = 5,
                         long_maxit = 7)))
stopifnot(override$multistart$long_maxit == 7)
expect_error(hmmSSF(~x + y, n_states = 2, data = dat, ssf_par0 = ssf_start,
                    multistart_opts = list(seed = 1)))
expect_error(hmmSSF(~x + y, n_states = 2, data = dat, ssf_par0 = ssf_start,
                    diagnose = TRUE, diagnose_opts = list(unknown = TRUE)))
expect_error(hmmSSF(~x + y, n_states = 2, data = dat, ssf_par0 = ssf_start,
                    multistart = TRUE, multistart_opts = list(n_start = 2)))

# A numerically invalid trial is rejected without discarding a valid start.
# Unrelated errors are not converted to objective values.
counter <- new.env(parent = emptyenv())
counter$rejected <- 0L
numerical_error <- getFromNamespace("hmmssf_numerical_error", "hmmSSF")
safe_factory <- getFromNamespace("hmmssf_safe_objective", "hmmSSF")
toy_objective <- function(x) {
  if (abs(x) > 2) numerical_error("outside the numerically evaluable region")
  100 * (x - 1)^2
}
guarded <- safe_factory(toy_objective, counter)
toy_fit <- optim(0, guarded, method = "BFGS")
stopifnot(counter$rejected > 0L, toy_fit$convergence == 0L,
          abs(toy_fit$par - 1) < 1e-3)
other_error <- safe_factory(function(x) stop("unexpected error"), counter)
stopifnot(inherits(try(other_error(0), silent = TRUE), "try-error"))
stopifnot("rejected_evaluations" %in% names(fit$multistart$summary))

# Finite individual terms can still overflow when added in the denominator.
sum_overflow <- data.frame(ID = 1, stratum = rep(1:2, each = 3),
                           obs = rep(c(1, 0, 0), 2), x = 1, w = 1)
make_inputs <- getFromNamespace("hmmssf_objective_inputs", "hmmSSF")
checked_objective <- getFromNamespace("hmmssf_checked_nllk", "hmmSSF")
inputs <- make_inputs(~x, ~1, sum_overflow, 2)
failed <- tryCatch(checked_objective(c(709, 0, -2, -2), inputs, sum_overflow$obs == 1),
                   hmmssf_numerical_error = identity)
stopifnot(inherits(failed, "hmmssf_numerical_error"))

# A one-covariate single-start fit retains the parameter names and shape.
single_cov <- hmmSSF(~0 + x, n_states = 2, data = dat,
                     ssf_par0 = matrix(c(1, -1), nrow = 1),
                     optim_opts = list(maxit = 5))
pred <- predict_ssf(single_cov, dat[1:2, ])
expected <- exp(as.matrix(dat[1:2, "x", drop = FALSE]) %*% single_cov$par$ssf)
stopifnot(identical(dim(pred), c(2L, 2L)),
          isTRUE(all.equal(unname(pred), unname(expected))))

# A missing Hessian still allows printing and prediction without intervals.
no_hessian <- fit
no_hessian$fit$hessian <- NULL
stopifnot(any(grepl("Confidence intervals unavailable", capture.output(print(no_hessian)))),
          all(is.finite(predict_tpm(no_hessian, dat[1, ], return_CI = FALSE)$mle)))
expect_error(confint(no_hessian))

three_states <- suppressWarnings(hmmSSF(~x, n_states = 3, data = dat,
  ssf_par0 = matrix(c(1, 0, -1), nrow = 1), multistart = TRUE,
  multistart_opts = list(n_starts = 1, n_long = 1, short_maxit = 5, long_maxit = 5)))
three_states$fit$hessian <- diag(length(three_states$fit$par))
intervals <- confint(three_states, pretty = TRUE)
stopifnot(nrow(intervals$ssf) == 3L, nrow(intervals$tpm) == 6L,
          length(capture.output(print(three_states))) > 0L)

# Numeric stratum identifiers must be compared at their original precision.
nearby <- dat
nearby$stratum <- 1 + nearby$stratum * .Machine$double.eps
nearby_fit <- suppressWarnings(hmmSSF_multistart(~x + y, n_states = 2,
  data = nearby, ssf_par0 = ssf_start, n_starts = 1, n_long = 1,
  short_maxit = 5, long_maxit = 5))
stopifnot(is.finite(nearby_fit$fit$value))
