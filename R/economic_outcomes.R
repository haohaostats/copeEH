trapezoid_integral <- function(t, y) {
  if (!is.numeric(t) || !is.numeric(y) || length(t) != length(y) ||
      length(t) < 2L || anyNA(t) || anyNA(y) ||
      any(!is.finite(t)) || any(!is.finite(y))) {
    stop("`t` and `y` must be finite numeric vectors of equal length (>= 2).",
         call. = FALSE)
  }
  if (any(diff(t) <= 0)) {
    stop("`t` must be strictly increasing.", call. = FALSE)
  }
  sum(diff(t) * (y[-length(y)] + y[-1L]) / 2)
}

restricted_mean_survival <- function(t, survival) {
  if (any(survival < 0 | survival > 1)) {
    stop("`survival` must lie in [0, 1].", call. = FALSE)
  }
  trapezoid_integral(t, survival)
}

state_time_summary <- function(prediction) {
  required <- c(
    "time", "os", "probability_pf", "probability_pd", "probability_death"
  )
  missing_columns <- setdiff(required, names(prediction))
  if (length(missing_columns) > 0L) {
    stop(
      paste("Prediction is missing columns:", paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }

  c(
    pf_time = trapezoid_integral(prediction$time, prediction$probability_pf),
    pd_time = trapezoid_integral(prediction$time, prediction$probability_pd),
    death_time = trapezoid_integral(prediction$time,
                                    prediction$probability_death),
    total_life_years = trapezoid_integral(prediction$time, prediction$os)
  )
}

discounted_partitioned_outcomes <- function(
    prediction,
    utility_pf = 0.80,
    utility_pd = 0.65,
    annual_cost_pf = 20000,
    annual_cost_pd = 100000,
    annual_discount_rate = 0.03) {
  required <- c("time", "probability_pf", "probability_pd")
  missing_columns <- setdiff(required, names(prediction))
  if (length(missing_columns) > 0L) {
    stop(
      paste("Prediction is missing columns:", paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }
  inputs <- c(
    utility_pf, utility_pd, annual_cost_pf, annual_cost_pd,
    annual_discount_rate
  )
  if (any(!is.finite(inputs)) || utility_pf < 0 || utility_pd < 0 ||
      annual_cost_pf < 0 || annual_cost_pd < 0 || annual_discount_rate < 0) {
    stop("Utilities, costs, and discount rate must be finite and non-negative.",
         call. = FALSE)
  }

  discount <- (1 + annual_discount_rate)^(-prediction$time)
  qaly_rate <- utility_pf * prediction$probability_pf +
    utility_pd * prediction$probability_pd
  cost_rate <- annual_cost_pf * prediction$probability_pf +
    annual_cost_pd * prediction$probability_pd

  c(
    qaly = trapezoid_integral(prediction$time, discount * qaly_rate),
    state_cost = trapezoid_integral(prediction$time, discount * cost_rate)
  )
}

expand_economic_price_scenarios <- function(
    delta_qaly_truth, state_cost_truth,
    delta_qaly_estimate, state_cost_estimate,
    wtp, price_ratios = c(0.8, 1.0, 1.2)) {
  values <- c(delta_qaly_truth, state_cost_truth, delta_qaly_estimate,
              state_cost_estimate, wtp, price_ratios)
  if (anyNA(values) || any(!is.finite(values)) || wtp <= 0 ||
      any(price_ratios <= 0)) {
    stop("Economic price-scenario inputs must be finite and positive where required.",
         call. = FALSE)
  }
  if (delta_qaly_truth <= 0) {
    stop("Price scenarios require a positive true incremental QALY.",
         call. = FALSE)
  }
  rows <- lapply(price_ratios, function(ratio) {
    acquisition_cost <- ratio * wtp * delta_qaly_truth - state_cost_truth
    if (acquisition_cost < 0) {
      stop(sprintf(
        "The %.1fk price scenario implies a negative acquisition cost.", ratio
      ), call. = FALSE)
    }
    true_delta_cost <- state_cost_truth + acquisition_cost
    estimated_delta_cost <- state_cost_estimate + acquisition_cost
    truth <- c(
      delta_qaly = delta_qaly_truth,
      delta_cost = true_delta_cost,
      inmb = wtp * delta_qaly_truth - true_delta_cost,
      icer = true_delta_cost / delta_qaly_truth
    )
    estimate <- c(
      delta_qaly = delta_qaly_estimate,
      delta_cost = estimated_delta_cost,
      inmb = wtp * delta_qaly_estimate - estimated_delta_cost,
      icer = if (abs(delta_qaly_estimate) > 1e-10) {
        estimated_delta_cost / delta_qaly_estimate
      } else {
        NA_real_
      }
    )
    data.frame(
      price_scenario = paste0(format(ratio, nsmall = 1), "k"),
      price_ratio = ratio,
      incremental_acquisition_cost = acquisition_cost,
      estimand = names(estimate),
      truth = unname(truth), estimate = unname(estimate),
      error = unname(estimate - truth),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
