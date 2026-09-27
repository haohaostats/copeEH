.validate_endpoint_data <- function(data, endpoint_name) {
  if (!is.data.frame(data) || !all(c("time", "status") %in% names(data))) {
    stop(sprintf("`%s` must be a data frame with `time` and `status` columns.",
                 endpoint_name), call. = FALSE)
  }
  if (!is.numeric(data$time) || anyNA(data$time) ||
      any(!is.finite(data$time)) || any(data$time < 0)) {
    stop(sprintf("`%s$time` must contain finite non-negative values.",
                 endpoint_name), call. = FALSE)
  }
  if (!is.numeric(data$status) && !is.integer(data$status) &&
      !is.logical(data$status)) {
    stop(sprintf("`%s$status` must be binary.", endpoint_name), call. = FALSE)
  }
  if (anyNA(data$status) || any(!data$status %in% c(0, 1))) {
    stop(sprintf("`%s$status` must contain only 0 and 1.", endpoint_name),
         call. = FALSE)
  }
  invisible(TRUE)
}

.prepare_endpoint_likelihood <- function(data, degree, eta,
                                         population_survival,
                                         population_hazard,
                                         endpoint_name) {
  .validate_endpoint_data(data, endpoint_name)
  degree <- .cope_validate_degree(degree)
  if (degree < 1L) {
    stop("COPE-EH estimation requires `degree >= 1`.", call. = FALSE)
  }
  eta <- .cope_validate_eta(eta)
  t <- as.numeric(data$time)
  u <- cope_time_map(t, eta)
  basis <- bernstein_basis(u, degree)
  basis_derivative <- bernstein_basis_derivative(u, degree)
  basis_derivative <- basis_derivative * (eta / (t + eta)^2)
  s_pop <- .evaluate_prediction_input(population_survival, t,
                                      "population_survival")
  h_pop <- .evaluate_prediction_input(population_hazard, t,
                                      "population_hazard")

  if (any(s_pop <= 0 | s_pop > 1)) {
    stop("Population survival must lie in (0, 1] at all observed times.",
         call. = FALSE)
  }
  if (any(h_pop < 0)) {
    stop("Population hazard must be non-negative at all observed times.",
         call. = FALSE)
  }

  list(
    time = t,
    status = as.integer(data$status),
    basis = basis,
    basis_derivative = basis_derivative,
    log_population_survival = log(s_pop),
    population_hazard = h_pop,
    endpoint = endpoint_name
  )
}

.endpoint_log_likelihood <- function(theta_free, prepared) {
  theta <- c(1, theta_free)
  relative_survival <- as.vector(prepared$basis %*% theta)
  derivative <- as.vector(prepared$basis_derivative %*% theta)

  if (any(!is.finite(relative_survival)) || any(relative_survival <= 0)) {
    return(-Inf)
  }
  excess_hazard <- -derivative / relative_survival
  total_hazard <- prepared$population_hazard + excess_hazard
  event_rows <- prepared$status == 1L
  if (any(!is.finite(total_hazard[event_rows])) ||
      any(total_hazard[event_rows] <= 0)) {
    return(-Inf)
  }

  log_likelihood <- prepared$log_population_survival + log(relative_survival)
  log_likelihood[event_rows] <- log_likelihood[event_rows] +
    log(total_hazard[event_rows])
  sum(log_likelihood)
}

.endpoint_score <- function(theta_free, prepared) {
  theta <- c(1, theta_free)
  relative_survival <- as.vector(prepared$basis %*% theta)
  derivative <- as.vector(prepared$basis_derivative %*% theta)
  excess_hazard <- -derivative / relative_survival
  total_hazard <- prepared$population_hazard + excess_hazard

  basis_free <- prepared$basis[, -1L, drop = FALSE]
  derivative_free <- prepared$basis_derivative[, -1L, drop = FALSE]
  score_log_survival <- basis_free / relative_survival
  hazard_derivative <-
    -derivative_free / relative_survival +
    basis_free * (derivative / relative_survival^2)
  event_weight <- prepared$status / total_hazard
  colSums(score_log_survival + hazard_derivative * event_weight)
}

build_cope_linear_constraints <- function(
    degree, enforce_order = TRUE, tail = c("free", "no_cure")) {
  degree <- .cope_validate_degree(degree)
  tail <- match.arg(tail)
  if (degree < 1L) {
    stop("COPE-EH estimation requires `degree >= 1`.", call. = FALSE)
  }
  if (tail == "no_cure" && degree < 2L) {
    stop("The no-cure specification requires `degree >= 2`.", call. = FALSE)
  }
  endpoint_parameter_count <- degree - as.integer(tail == "no_cure")
  parameter_count <- 2L * endpoint_parameter_count
  rows <- list()
  rhs <- numeric()
  labels <- character()

  add_constraint <- function(row, boundary, label) {
    rows[[length(rows) + 1L]] <<- row
    rhs[length(rhs) + 1L] <<- boundary
    labels[length(labels) + 1L] <<- label
  }

  for (endpoint in c("os", "pfs")) {
    offset <- if (endpoint == "os") 0L else endpoint_parameter_count

    row <- numeric(parameter_count)
    row[offset + 1L] <- -1
    add_constraint(row, -1, paste0(endpoint, "_theta_1_le_1"))

    if (endpoint_parameter_count >= 2L) {
      for (k in 2L:endpoint_parameter_count) {
        row <- numeric(parameter_count)
        row[offset + k - 1L] <- 1
        row[offset + k] <- -1
        add_constraint(row, 0,
                       paste0(endpoint, "_monotone_", k - 1L, "_", k))
      }
    }

    row <- numeric(parameter_count)
    row[offset + endpoint_parameter_count] <- 1
    add_constraint(row, 0, paste0(endpoint, "_last_free_theta_ge_0"))
  }

  if (isTRUE(enforce_order)) {
    for (k in seq_len(endpoint_parameter_count)) {
      row <- numeric(parameter_count)
      row[k] <- 1
      row[endpoint_parameter_count + k] <- -1
      add_constraint(row, 0, paste0("pfs_le_os_", k))
    }
  }

  ui <- do.call(rbind, rows)
  rownames(ui) <- labels
  colnames(ui) <- c(
    paste0("os_theta_", seq_len(endpoint_parameter_count)),
    paste0("pfs_theta_", seq_len(endpoint_parameter_count))
  )
  names(rhs) <- labels
  list(ui = ui, ci = rhs)
}

default_cope_start <- function(
    degree, version = c("primary", "fallback"),
    tail = c("free", "no_cure")) {
  degree <- .cope_validate_degree(degree)
  version <- match.arg(version)
  tail <- match.arg(tail)
  k <- seq_len(degree)
  fraction <- k / (degree + 1)

  if (version == "primary") {
    theta_os <- exp(-1.0 * fraction)
    theta_pfs <- exp(-2.1 * fraction)
  } else {
    theta_os <- exp(-0.55 * fraction)
    theta_pfs <- exp(-1.45 * fraction)
  }
  if (tail == "no_cure") {
    theta_os <- theta_os[-degree]
    theta_pfs <- theta_pfs[-degree]
  }
  c(theta_os, theta_pfs)
}

.check_strict_feasibility <- function(parameters, constraints,
                                      tolerance = 1e-10) {
  margins <- as.vector(constraints$ui %*% parameters - constraints$ci)
  if (any(!is.finite(margins)) || any(margins <= tolerance)) {
    failing <- rownames(constraints$ui)[margins <= tolerance]
    stop(
      paste(
        "Starting values must be strictly inside every constraint; failing:",
        paste(failing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.cope_penalty <- function(parameters, degree, penalty_strength) {
  if (penalty_strength == 0 || degree < 2L) {
    return(0)
  }
  theta_os <- c(1, parameters[seq_len(degree)])
  theta_pfs <- c(1, parameters[degree + seq_len(degree)])
  penalty_strength * (
    sum(diff(theta_os, differences = 2L)^2) +
      sum(diff(theta_pfs, differences = 2L)^2)
  )
}

.cope_penalty_gradient <- function(parameters, degree, penalty_strength) {
  if (penalty_strength == 0 || degree < 2L) {
    return(rep(0, length(parameters)))
  }
  second_difference_matrix <- diff(diag(degree + 1L), differences = 2L)
  endpoint_gradient <- function(theta_free) {
    theta <- c(1, theta_free)
    as.vector(
      2 * penalty_strength *
        crossprod(second_difference_matrix,
                  second_difference_matrix %*% theta)
    )[-1L]
  }
  c(
    endpoint_gradient(parameters[seq_len(degree)]),
    endpoint_gradient(parameters[degree + seq_len(degree)])
  )
}

.safe_inverse_hessian <- function(hessian, tolerance = 1e-8) {
  if (is.null(hessian) || anyNA(hessian) || any(!is.finite(hessian))) {
    return(NULL)
  }
  hessian <- (hessian + t(hessian)) / 2
  eigenvalues <- tryCatch(
    eigen(hessian, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) numeric()
  )
  if (length(eigenvalues) == 0L || min(eigenvalues) <= tolerance) {
    return(NULL)
  }
  tryCatch(solve(hessian), error = function(e) NULL)
}

fit_cope_eh <- function(data_os, data_pfs, degree, eta,
                        population_survival, population_hazard,
                        penalty_strength = 0,
                        start = NULL,
                        control = list(maxit = 2000, reltol = 1e-10),
                        retry = TRUE,
                        enforce_order = TRUE,
                        tail = c("free", "no_cure")) {
  degree <- .cope_validate_degree(degree)
  eta <- .cope_validate_eta(eta)
  tail <- match.arg(tail)
  if (tail == "no_cure" && degree < 2L) {
    stop("The no-cure specification requires `degree >= 2`.", call. = FALSE)
  }
  endpoint_parameter_count <- degree - as.integer(tail == "no_cure")
  if (!is.numeric(penalty_strength) || length(penalty_strength) != 1L ||
      is.na(penalty_strength) || penalty_strength < 0) {
    stop("`penalty_strength` must be one non-negative number.", call. = FALSE)
  }

  prepared_os <- .prepare_endpoint_likelihood(
    data = data_os,
    degree = degree,
    eta = eta,
    population_survival = population_survival,
    population_hazard = population_hazard,
    endpoint_name = "data_os"
  )
  prepared_pfs <- .prepare_endpoint_likelihood(
    data = data_pfs,
    degree = degree,
    eta = eta,
    population_survival = population_survival,
    population_hazard = population_hazard,
    endpoint_name = "data_pfs"
  )
  constraints <- build_cope_linear_constraints(degree, enforce_order, tail)

  expand_endpoint <- function(parameters) {
    if (tail == "no_cure") c(parameters, 0) else parameters
  }
  expand_all <- function(parameters) {
    c(
      expand_endpoint(parameters[seq_len(endpoint_parameter_count)]),
      expand_endpoint(parameters[
        endpoint_parameter_count + seq_len(endpoint_parameter_count)
      ])
    )
  }

  objective <- function(parameters) {
    ll_os <- .endpoint_log_likelihood(
      expand_endpoint(parameters[seq_len(endpoint_parameter_count)]),
      prepared_os
    )
    ll_pfs <- .endpoint_log_likelihood(
      expand_endpoint(parameters[
        endpoint_parameter_count + seq_len(endpoint_parameter_count)
      ]), prepared_pfs
    )
    if (!is.finite(ll_os) || !is.finite(ll_pfs)) {
      return(.Machine$double.xmax / 100)
    }
    -(ll_os + ll_pfs) +
      .cope_penalty(expand_all(parameters), degree, penalty_strength)
  }
  objective_gradient <- function(parameters) {
    os_parameters <- parameters[seq_len(endpoint_parameter_count)]
    pfs_parameters <- parameters[
      endpoint_parameter_count + seq_len(endpoint_parameter_count)
    ]
    likelihood_gradient <- -c(
      .endpoint_score(expand_endpoint(os_parameters), prepared_os)[
        seq_len(endpoint_parameter_count)
      ],
      .endpoint_score(expand_endpoint(pfs_parameters), prepared_pfs)[
        seq_len(endpoint_parameter_count)
      ]
    )
    full_penalty_gradient <- .cope_penalty_gradient(
      expand_all(parameters), degree, penalty_strength
    )
    kept <- c(
      seq_len(endpoint_parameter_count),
      degree + seq_len(endpoint_parameter_count)
    )
    likelihood_gradient + full_penalty_gradient[kept]
  }

  if (is.null(start)) {
    start <- default_cope_start(degree, "primary", tail)
  }
  if (!is.numeric(start) ||
      length(start) != 2L * endpoint_parameter_count ||
      anyNA(start) || any(!is.finite(start))) {
    stop(sprintf(
      "`start` must contain exactly %d finite values.",
      2L * endpoint_parameter_count
    ), call. = FALSE)
  }
  .check_strict_feasibility(start, constraints)

  run_optimisation <- function(initial) {
    tryCatch(
      stats::constrOptim(
        theta = initial,
        f = objective,
        grad = objective_gradient,
        ui = constraints$ui,
        ci = constraints$ci,
        method = "BFGS",
        control = control
      ),
      error = function(e) {
        list(
          par = rep(NA_real_, 2L * endpoint_parameter_count),
          value = NA_real_,
          convergence = 100L,
          message = conditionMessage(e),
          counts = c("function" = NA_integer_, "gradient" = NA_integer_)
        )
      }
    )
  }

  optimisation <- run_optimisation(start)
  used_fallback <- FALSE
  if (isTRUE(retry) && optimisation$convergence != 0L) {
    fallback <- default_cope_start(degree, "fallback", tail)
    optimisation <- run_optimisation(fallback)
    used_fallback <- TRUE
  }

  converged <- isTRUE(optimisation$convergence == 0L) &&
    all(is.finite(optimisation$par))
  hessian <- NULL
  covariance <- NULL
  constraint_margin <- rep(NA_real_, nrow(constraints$ui))

  if (converged) {
    hessian <- tryCatch(
      stats::optimHess(optimisation$par, objective, objective_gradient),
      error = function(e) NULL
    )
    covariance <- .safe_inverse_hessian(hessian)
    constraint_margin <- as.vector(
      constraints$ui %*% optimisation$par - constraints$ci
    )
    names(constraint_margin) <- rownames(constraints$ui)
  }

  theta_os <- c(
    1,
    expand_endpoint(optimisation$par[seq_len(endpoint_parameter_count)])
  )
  theta_pfs <- c(
    1,
    expand_endpoint(optimisation$par[
      endpoint_parameter_count + seq_len(endpoint_parameter_count)
    ])
  )
  if (converged && isTRUE(enforce_order)) {
    validate_cope_coefficients(theta_os, theta_pfs, tolerance = 1e-7)
  }

  result <- list(
    coefficients = c(
      stats::setNames(theta_os, paste0("os_theta_", 0L:degree)),
      stats::setNames(theta_pfs, paste0("pfs_theta_", 0L:degree))
    ),
    free_coefficients = stats::setNames(
      optimisation$par,
      colnames(constraints$ui)
    ),
    degree = degree,
    eta = eta,
    tail = tail,
    enforce_order = isTRUE(enforce_order),
    penalty_strength = penalty_strength,
    log_likelihood = if (converged) {
      .endpoint_log_likelihood(
        expand_endpoint(optimisation$par[seq_len(endpoint_parameter_count)]),
        prepared_os
      ) +
        .endpoint_log_likelihood(
          expand_endpoint(optimisation$par[
            endpoint_parameter_count + seq_len(endpoint_parameter_count)
          ]), prepared_pfs
        )
    } else {
      NA_real_
    },
    objective = optimisation$value,
    convergence = optimisation$convergence,
    converged = converged,
    message = optimisation$message %||% NULL,
    counts = optimisation$counts,
    used_fallback = used_fallback,
    hessian = hessian,
    covariance = covariance,
    covariance_available = !is.null(covariance),
    minimum_constraint_margin = suppressWarnings(min(constraint_margin,
                                                       na.rm = TRUE)),
    constraint_margin = constraint_margin,
    n_os = nrow(data_os),
    n_pfs = nrow(data_pfs),
    events_os = sum(data_os$status),
    events_pfs = sum(data_pfs$status),
    working_independence = TRUE,
    call = match.call()
  )
  if (!is.finite(result$minimum_constraint_margin)) {
    result$minimum_constraint_margin <- NA_real_
  }
  result$boundary_solution <- isTRUE(
    result$minimum_constraint_margin < 1e-6
  )
  result$regular_asymptotic_inference <- isTRUE(
    result$covariance_available && !result$boundary_solution
  )
  class(result) <- "cope_eh_fit"
  result
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

coef.cope_eh_fit <- function(object, free = FALSE, ...) {
  if (isTRUE(free)) object$free_coefficients else object$coefficients
}

vcov.cope_eh_fit <- function(object, ...) {
  if (is.null(object$covariance)) {
    stop("A positive-definite covariance matrix is not available.",
         call. = FALSE)
  }
  object$covariance
}

print.cope_eh_fit <- function(x, ...) {
  label <- if (isTRUE(x$enforce_order)) {
    "COPE-EH constrained marginal fit"
  } else {
    "Independent Bernstein-EH marginal fit"
  }
  cat(label, "\n")
  cat("  converged:", x$converged, "(code", x$convergence, ")\n")
  cat("  degree:", x$degree, " eta:", x$eta, "\n")
  cat("  tail:", x$tail %||% "free", "\n")
  cat("  OS/PFS events:", x$events_os, "/", x$events_pfs, "\n")
  cat("  log likelihood:", format(x$log_likelihood, digits = 7), "\n")
  cat("  covariance available:", x$covariance_available, "\n")
  invisible(x)
}

predict.cope_eh_fit <- function(object, t, population_survival,
                                population_hazard = NULL,
                                se.fit = FALSE, level = 0.95, ...) {
  if (!isTRUE(object$converged)) {
    stop("Cannot predict from a fit that did not converge.", call. = FALSE)
  }
  degree <- object$degree
  theta_os <- object$coefficients[paste0("os_theta_", 0L:degree)]
  theta_pfs <- object$coefficients[paste0("pfs_theta_", 0L:degree)]
  if (isTRUE(object$enforce_order)) {
    prediction <- predict_cope_eh(
      t = t,
      theta_os = unname(theta_os),
      theta_pfs = unname(theta_pfs),
      eta = object$eta,
      population_survival = population_survival,
      population_hazard = population_hazard
    )
  } else {
    t <- .cope_validate_time(t)
    s_pop <- .evaluate_prediction_input(population_survival, t,
                                        "population_survival")
    r_os <- bernstein_relative_survival(t, unname(theta_os), object$eta)
    r_pfs <- bernstein_relative_survival(t, unname(theta_pfs), object$eta)
    s_os <- s_pop * r_os
    s_pfs <- s_pop * r_pfs
    prediction <- data.frame(
      time = t,
      population_survival = s_pop,
      relative_os = r_os,
      relative_pfs = r_pfs,
      os = s_os,
      pfs = s_pfs,
      probability_pf = s_pfs,
      probability_pd = s_os - s_pfs,
      probability_death = 1 - s_os,
      excess_hazard_os = bernstein_excess_hazard(
        t, unname(theta_os), object$eta
      ),
      excess_hazard_pfs = bernstein_excess_hazard(
        t, unname(theta_pfs), object$eta
      ),
      check.names = FALSE
    )
    if (!is.null(population_hazard)) {
      h_pop <- .evaluate_prediction_input(population_hazard, t,
                                          "population_hazard")
      prediction$population_hazard <- h_pop
      prediction$total_hazard_os <- h_pop + prediction$excess_hazard_os
      prediction$total_hazard_pfs <- h_pop + prediction$excess_hazard_pfs
    }
  }

  if (!isTRUE(se.fit)) {
    return(prediction)
  }
  covariance <- vcov(object)
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must lie strictly between zero and one.", call. = FALSE)
  }
  basis_free <- bernstein_basis(
    cope_time_map(prediction$time, object$eta), degree
  )[, -1L, drop = FALSE]
  if (identical(object$tail, "no_cure")) {
    basis_free <- basis_free[, -ncol(basis_free), drop = FALSE]
  }
  endpoint_parameter_count <- ncol(basis_free)
  weighted_basis <- basis_free * prediction$population_survival
  zeros <- matrix(
    0, nrow = nrow(weighted_basis), ncol = endpoint_parameter_count
  )
  gradient_os <- cbind(weighted_basis, zeros)
  gradient_pfs <- cbind(zeros, weighted_basis)
  gradient_pd <- gradient_os - gradient_pfs
  gradient_death <- -gradient_os

  delta_se <- function(gradient) {
    variance <- rowSums((gradient %*% covariance) * gradient)
    sqrt(pmax(variance, 0))
  }
  z_value <- stats::qnorm(1 - (1 - level) / 2)

  for (quantity in c("os", "pfs", "probability_pd", "probability_death")) {
    gradient <- switch(
      quantity,
      os = gradient_os,
      pfs = gradient_pfs,
      probability_pd = gradient_pd,
      probability_death = gradient_death
    )
    standard_error <- delta_se(gradient)
    prediction[[paste0("se_", quantity)]] <- standard_error
    prediction[[paste0("lower_", quantity)]] <-
      prediction[[quantity]] - z_value * standard_error
    prediction[[paste0("upper_", quantity)]] <-
      prediction[[quantity]] + z_value * standard_error
  }
  prediction
}

fit_independent_bernstein_eh <- function(...) {
  fit_cope_eh(..., enforce_order = FALSE)
}

clip_partitioned_survival <- function(prediction) {
  prediction$pfs <- pmin(prediction$pfs, prediction$os)
  if ("pf" %in% names(prediction)) prediction$pf <- prediction$pfs
  if ("pd" %in% names(prediction)) prediction$pd <- prediction$os - prediction$pfs
  if ("death" %in% names(prediction)) prediction$death <- 1 - prediction$os
  prediction$probability_pf <- prediction$pfs
  prediction$probability_pd <- prediction$os - prediction$pfs
  prediction$probability_death <- 1 - prediction$os
  prediction
}
