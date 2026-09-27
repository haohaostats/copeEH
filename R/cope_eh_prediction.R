.evaluate_prediction_input <- function(x, t, name) {
  value <- if (is.function(x)) x(t) else x
  if (!is.numeric(value) || length(value) == 0L || anyNA(value) ||
      any(!is.finite(value))) {
    stop(sprintf("`%s` must evaluate to finite numeric values.", name),
         call. = FALSE)
  }
  if (length(value) == 1L) {
    value <- rep(value, length(t))
  }
  if (length(value) != length(t)) {
    stop(sprintf("`%s` must have length 1 or the same length as `t`.", name),
         call. = FALSE)
  }
  as.numeric(value)
}

predict_cope_eh <- function(t, theta_os, theta_pfs, eta,
                            population_survival,
                            population_hazard = NULL) {
  t <- .cope_validate_time(t)
  validate_cope_coefficients(theta_os, theta_pfs)
  eta <- .cope_validate_eta(eta)

  s_pop <- .evaluate_prediction_input(population_survival, t,
                                      "population_survival")
  if (any(s_pop < 0 | s_pop > 1)) {
    stop("`population_survival` must lie in [0, 1].", call. = FALSE)
  }

  r_os <- bernstein_relative_survival(t, theta_os, eta)
  r_pfs <- bernstein_relative_survival(t, theta_pfs, eta)
  s_os <- s_pop * r_os
  s_pfs <- s_pop * r_pfs

  out <- data.frame(
    time = t,
    population_survival = s_pop,
    relative_os = r_os,
    relative_pfs = r_pfs,
    os = s_os,
    pfs = s_pfs,
    probability_pf = s_pfs,
    probability_pd = s_os - s_pfs,
    probability_death = 1 - s_os,
    excess_hazard_os = bernstein_excess_hazard(t, theta_os, eta),
    excess_hazard_pfs = bernstein_excess_hazard(t, theta_pfs, eta),
    check.names = FALSE
  )

  if (!is.null(population_hazard)) {
    h_pop <- .evaluate_prediction_input(population_hazard, t,
                                        "population_hazard")
    if (any(h_pop < 0)) {
      stop("`population_hazard` must be non-negative.", call. = FALSE)
    }
    out$population_hazard <- h_pop
    out$total_hazard_os <- h_pop + out$excess_hazard_os
    out$total_hazard_pfs <- h_pop + out$excess_hazard_pfs
  }

  out
}

diagnose_cope_prediction <- function(prediction, tolerance = 1e-10) {
  required <- c(
    "os", "pfs", "probability_pf", "probability_pd", "probability_death"
  )
  missing_columns <- setdiff(required, names(prediction))
  if (length(missing_columns) > 0L) {
    stop(
      paste("Prediction is missing columns:", paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      is.na(tolerance) || tolerance < 0) {
    stop("`tolerance` must be one non-negative number.", call. = FALSE)
  }

  state_sum <- prediction$probability_pf + prediction$probability_pd +
    prediction$probability_death
  order_violation <- max(prediction$pfs - prediction$os, na.rm = TRUE)
  minimum_state_probability <- min(
    prediction$probability_pf,
    prediction$probability_pd,
    prediction$probability_death,
    na.rm = TRUE
  )
  maximum_state_probability <- max(
    prediction$probability_pf,
    prediction$probability_pd,
    prediction$probability_death,
    na.rm = TRUE
  )
  state_sum_error <- max(abs(state_sum - 1), na.rm = TRUE)

  list(
    valid = isTRUE(order_violation <= tolerance) &&
      isTRUE(minimum_state_probability >= -tolerance) &&
      isTRUE(maximum_state_probability <= 1 + tolerance) &&
      isTRUE(state_sum_error <= tolerance),
    maximum_pfs_minus_os = order_violation,
    minimum_state_probability = minimum_state_probability,
    maximum_state_probability = maximum_state_probability,
    maximum_state_sum_error = state_sum_error
  )
}
