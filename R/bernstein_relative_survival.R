.cope_validate_time <- function(t) {
  if (!is.numeric(t) || length(t) == 0L || anyNA(t) ||
      any(!is.finite(t)) || any(t < 0)) {
    stop("`t` must be a non-empty numeric vector of finite, non-negative times.",
         call. = FALSE)
  }
  as.numeric(t)
}

.cope_validate_eta <- function(eta) {
  if (!is.numeric(eta) || length(eta) != 1L || is.na(eta) ||
      !is.finite(eta) || eta <= 0) {
    stop("`eta` must be one finite positive number.", call. = FALSE)
  }
  as.numeric(eta)
}

.cope_validate_degree <- function(degree) {
  if (!is.numeric(degree) || length(degree) != 1L || is.na(degree) ||
      !is.finite(degree) || degree < 0 || degree != floor(degree)) {
    stop("`degree` must be one non-negative integer.", call. = FALSE)
  }
  as.integer(degree)
}

cope_time_map <- function(t, eta) {
  t <- .cope_validate_time(t)
  eta <- .cope_validate_eta(eta)
  t / (t + eta)
}

bernstein_basis <- function(u, degree) {
  degree <- .cope_validate_degree(degree)
  if (!is.numeric(u) || length(u) == 0L || anyNA(u) ||
      any(!is.finite(u)) || any(u < 0) || any(u > 1)) {
    stop("`u` must be a non-empty numeric vector with values in [0, 1].",
         call. = FALSE)
  }

  out <- vapply(
    X = 0L:degree,
    FUN = function(k) stats::dbinom(k, size = degree, prob = u),
    FUN.VALUE = numeric(length(u))
  )
  matrix(out, nrow = length(u), ncol = degree + 1L,
         dimnames = list(NULL, paste0("B", 0L:degree)))
}

bernstein_basis_derivative <- function(u, degree) {
  degree <- .cope_validate_degree(degree)
  if (!is.numeric(u) || length(u) == 0L || anyNA(u) ||
      any(!is.finite(u)) || any(u < 0) || any(u > 1)) {
    stop("`u` must be a non-empty numeric vector with values in [0, 1].",
         call. = FALSE)
  }
  if (degree == 0L) {
    return(matrix(0, nrow = length(u), ncol = 1L,
                  dimnames = list(NULL, "dB0")))
  }

  lower <- bernstein_basis(u, degree - 1L)
  deriv <- matrix(0, nrow = length(u), ncol = degree + 1L)

  for (k in 0L:degree) {
    left <- if (k >= 1L) lower[, k] else 0
    right <- if (k <= degree - 1L) lower[, k + 1L] else 0
    deriv[, k + 1L] <- degree * (left - right)
  }

  colnames(deriv) <- paste0("dB", 0L:degree)
  deriv
}

.cope_validate_theta <- function(theta, name = "theta") {
  if (!is.numeric(theta) || length(theta) < 1L || anyNA(theta) ||
      any(!is.finite(theta))) {
    stop(sprintf("`%s` must be a non-empty finite numeric vector.", name),
         call. = FALSE)
  }
  as.numeric(theta)
}

validate_cope_coefficients <- function(theta_os, theta_pfs,
                                       tolerance = 1e-10) {
  theta_os <- .cope_validate_theta(theta_os, "theta_os")
  theta_pfs <- .cope_validate_theta(theta_pfs, "theta_pfs")

  if (length(theta_os) != length(theta_pfs)) {
    stop("`theta_os` and `theta_pfs` must have the same length.", call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      is.na(tolerance) || tolerance < 0) {
    stop("`tolerance` must be one non-negative number.", call. = FALSE)
  }

  problems <- character()
  if (abs(theta_os[1L] - 1) > tolerance ||
      abs(theta_pfs[1L] - 1) > tolerance) {
    problems <- c(problems, "both coefficient vectors must start at 1")
  }
  if (any(theta_os < -tolerance | theta_os > 1 + tolerance) ||
      any(theta_pfs < -tolerance | theta_pfs > 1 + tolerance)) {
    problems <- c(problems, "all coefficients must lie in [0, 1]")
  }
  if (length(theta_os) > 1L &&
      (any(diff(theta_os) > tolerance) ||
       any(diff(theta_pfs) > tolerance))) {
    problems <- c(problems, "each coefficient vector must be non-increasing")
  }
  if (any(theta_pfs > theta_os + tolerance)) {
    problems <- c(problems, "PFS coefficients must not exceed OS coefficients")
  }

  if (length(problems) > 0L) {
    stop(paste("Invalid COPE-EH coefficients:", paste(problems, collapse = "; ")),
         call. = FALSE)
  }

  invisible(TRUE)
}

bernstein_relative_survival <- function(t, theta, eta) {
  t <- .cope_validate_time(t)
  theta <- .cope_validate_theta(theta)
  eta <- .cope_validate_eta(eta)
  degree <- length(theta) - 1L
  as.vector(bernstein_basis(cope_time_map(t, eta), degree) %*% theta)
}

bernstein_relative_survival_derivative <- function(t, theta, eta) {
  t <- .cope_validate_time(t)
  theta <- .cope_validate_theta(theta)
  eta <- .cope_validate_eta(eta)
  degree <- length(theta) - 1L
  u <- cope_time_map(t, eta)
  d_r_du <- as.vector(bernstein_basis_derivative(u, degree) %*% theta)
  d_u_dt <- eta / (t + eta)^2
  d_r_du * d_u_dt
}

bernstein_excess_hazard <- function(t, theta, eta) {
  rel_survival <- bernstein_relative_survival(t, theta, eta)
  derivative <- bernstein_relative_survival_derivative(t, theta, eta)
  out <- -derivative / rel_survival
  out[rel_survival <= 0] <- Inf
  out
}

bernstein_relative_survival_limit <- function(theta) {
  theta <- .cope_validate_theta(theta)
  theta[length(theta)]
}
