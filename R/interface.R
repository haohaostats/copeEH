.scalar <- function(x, name, lower = 0, strict = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      (if (strict) x <= lower else x < lower)) {
    stop(sprintf("`%s` must be one finite number %s %s.", name,
                 if (strict) ">" else ">=", lower), call. = FALSE)
  }
  invisible(x)
}

population_mortality <- function(breaks = 0, hazard) {
  breaks <- .cope_validate_time(breaks)
  if (breaks[1L] != 0 || any(diff(breaks) <= 0)) {
    stop("`breaks` must start at zero and be strictly increasing.", call. = FALSE)
  }
  if (!is.numeric(hazard) || length(hazard) != length(breaks) ||
      anyNA(hazard) || any(!is.finite(hazard)) || any(hazard < 0)) {
    stop("Supply one finite non-negative hazard per break.", call. = FALSE)
  }
  hazard <- as.numeric(hazard)
  cumulative <- c(0, cumsum(diff(breaks) * head(hazard, -1L)))
  H <- function(t) {
    t <- .cope_validate_time(t)
    i <- findInterval(t, breaks)
    cumulative[i] + (t - breaks[i]) * hazard[i]
  }
  list(survival = function(t) exp(-H(t)),
       hazard = function(t) {
         t <- .cope_validate_time(t)
         hazard[findInterval(t, breaks)]
       }, cumhazard = H, breaks = breaks, rates = hazard)
}

fit_cope <- function(data_os, data_pfs, population, degree = 3, eta = 2,
                     tail = c("free", "no_cure"),
                     control = list(maxit = 2000, reltol = 1e-10)) {
  tail <- match.arg(tail)
  for (nm in c("data_os", "data_pfs")) {
    d <- get(nm)
    .validate_endpoint_data(d, nm)
    if (nrow(d) < 1L) stop("Endpoint data must not be empty.", call. = FALSE)
  }
  if (!is.list(population) || !is.function(population$survival) ||
      !is.function(population$hazard)) {
    stop("`population` must supply survival(t) and hazard(t) functions.",
         call. = FALSE)
  }
  if (!is.list(control)) stop("`control` must be a list.", call. = FALSE)
  origin <- .evaluate_prediction_input(population$survival, 0, "survival")
  if (abs(origin - 1) > 1e-8) {
    stop("Population survival must equal one at time zero.", call. = FALSE)
  }
  observed <- sort(unique(c(0, data_os$time, data_pfs$time)))
  sp <- .evaluate_prediction_input(population$survival, observed, "survival")
  if (any(diff(sp) > 1e-10)) {
    stop("Population survival must be non-increasing.", call. = FALSE)
  }
  ans <- fit_cope_eh(data_os, data_pfs, degree, eta,
                   population$survival, population$hazard,
                   penalty_strength = 0, control = control, tail = tail)
  ans$population <- population
  ans$observed_end <- max(observed)
  ans$call <- match.call()
  class(ans) <- "cope_fit"
  if (!ans$converged) warning("Fit did not converge; predictions are unavailable.",
                             call. = FALSE)
  ans
}

coef.cope_fit <- function(object, free = FALSE, ...) {
  coef.cope_eh_fit(object, free = free)
}

predict.cope_fit <- function(object, times, type = c("all", "survival", "state"),
                             ...) {
  if (length(list(...))) stop("Unused prediction arguments; interval estimates are not supported.",
                             call. = FALSE)
  type <- match.arg(type)
  times <- .cope_validate_time(times)
  ord <- order(times)
  s <- .evaluate_prediction_input(object$population$survival, times, "survival")
  if (any(diff(s[ord]) > 1e-10)) stop("Population survival must be non-increasing.",
                                    call. = FALSE)
  out <- predict.cope_eh_fit(object, t = times,
                           population_survival = object$population$survival,
                           population_hazard = object$population$hazard)
  columns <- switch(type, all = names(out), survival = c("time", "os", "pfs"),
                    state = c("time", "probability_pf", "probability_pd",
                              "probability_death"))
  out[columns]
}

summary.cope_fit <- function(object, ...) {
  out <- list(call = object$call, degree = object$degree, eta = object$eta,
              tail = object$tail, converged = object$converged,
              composite_loglik = object$log_likelihood,
              endpoints = data.frame(endpoint = c("OS", "PFS"),
                                     n = c(object$n_os, object$n_pfs),
                                     events = c(object$events_os, object$events_pfs)),
              coefficients = data.frame(k = 0:object$degree,
                OS = unname(object$coefficients[seq_len(object$degree + 1L)]),
                PFS = unname(object$coefficients[object$degree + 1L +
                                                seq_len(object$degree + 1L)])))
  class(out) <- "summary.cope_fit"
  out
}

print.summary.cope_fit <- function(x, ...) {
  cat("COPE-EH: jointly constrained marginal survival\n")
  cat("Degree:", x$degree, "  eta:", x$eta, "  tail:", x$tail, "\n")
  cat("Converged:", x$converged, "  Composite log-likelihood:",
      format(x$composite_loglik, digits = 7), "\n\n")
  print(x$endpoints, row.names = FALSE)
  cat("\nBernstein relative-survival coefficients:\n")
  print(x$coefficients, row.names = FALSE)
  cat("\nPoint estimates only; endpoint dependence is not recovered.\n")
  invisible(x)
}

print.cope_fit <- function(x, ...) {
  print(summary(x))
  invisible(x)
}

plot.cope_fit <- function(x, times = seq(0, max(5, 2 * x$observed_end),
                                        length.out = 301),
                          type = c("survival", "state"),
                          xlab = "Time", ylab = "Probability",
                          col = c("#0072B2", "#D55E00", "#666666"),
                          lty = c(1, 2, 3), lwd = 2, legend = TRUE, ...) {
  type <- match.arg(type)
  times <- sort(unique(.cope_validate_time(times)))
  p <- predict(x, times, type = type)
  graphics::matplot(p$time, as.matrix(p[-1L]), type = "l", ylim = c(0, 1),
                    xlab = xlab, ylab = ylab, col = col, lty = lty,
                    lwd = lwd, bty = "l", ...)
  if (isTRUE(legend)) {
    labels <- if (type == "survival") c("OS", "PFS") else c("PF", "PD", "Death")
    graphics::legend("topright", legend = labels,
                     col = rep_len(col, length(labels)),
                     lty = rep_len(lty, length(labels)), lwd = lwd, bty = "n")
  }
  invisible(p)
}

state_outcomes <- function(object, horizon, utility_pf = 0.8, utility_pd = 0.65,
                            annual_cost_pf = 0, annual_cost_pd = 0,
                            discount = 0.03, step = 0.02) {
  .scalar(horizon, "horizon", strict = TRUE)
  .scalar(step, "step", strict = TRUE)
  for (nm in c("utility_pf", "utility_pd", "annual_cost_pf", "annual_cost_pd",
               "discount")) .scalar(get(nm), nm)
  if (!inherits(object, "cope_fit")) stop("`object` must be a cope_fit.", call. = FALSE)
  times <- sort(unique(c(seq(0, horizon, by = step), horizon)))
  p <- predict(object, times)
  st <- state_time_summary(p)
  ec <- discounted_partitioned_outcomes(p, utility_pf, utility_pd,
                                         annual_cost_pf, annual_cost_pd, discount)
  data.frame(horizon = horizon, discount = discount,
             pf_years = unname(st["pf_time"]), pd_years = unname(st["pd_time"]),
             life_years = unname(st["total_life_years"]),
             qaly = unname(ec["qaly"]), state_cost = unname(ec["state_cost"]))
}

compare_outcomes <- function(treatment, control, wtp, incremental_cost = 0) {
  for (nm in c("treatment", "control")) {
    x <- get(nm)
    needed <- c("horizon", "discount", "life_years", "qaly", "state_cost")
    if (!is.data.frame(x) || nrow(x) != 1L || !all(needed %in% names(x)) ||
        any(!vapply(x[needed], is.numeric, logical(1))) ||
        any(!is.finite(unlist(x[needed])))) {
      stop("Supply single-row results from state_outcomes().", call. = FALSE)
    }
  }
  if (treatment$horizon != control$horizon || treatment$discount != control$discount)
    stop("Treatment and control must use identical horizons and discount rates.", call. = FALSE)
  .scalar(wtp, "wtp", strict = TRUE)
  if (!is.numeric(incremental_cost) || length(incremental_cost) != 1L ||
      !is.finite(incremental_cost)) stop("`incremental_cost` must be finite.", call. = FALSE)
  dq <- treatment$qaly - control$qaly
  dc <- treatment$state_cost - control$state_cost + incremental_cost
  classification <- if (dq > 0 && dc <= 0) "Dominant" else
    if (dq <= 0 && dc > 0) "Dominated" else
    if (dq > 0 && dc > 0) "More effective, more costly" else
    if (dq < 0 && dc < 0) "Less effective, less costly" else "No strict dominance"
  data.frame(delta_ly = treatment$life_years - control$life_years,
             delta_qaly = dq, delta_cost = dc,
             icer = if (dq * dc > 0) dc / dq else NA_real_,
             inmb = wtp * dq - dc, classification = classification)
}

cope_example <- function(arm = c("control", "treatment")) {
  arm <- match.arg(arm)
  env <- new.env(parent = emptyenv())
  utils::data("cope_demo", package = "copeEH", envir = env)
  d <- env$cope_demo
  list(os = d$os[d$os$arm == arm, c("time", "status")],
       pfs = d$pfs[d$pfs$arm == arm, c("time", "status")],
       population = do.call(population_mortality, d$population))
}
