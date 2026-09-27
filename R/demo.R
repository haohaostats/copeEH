run_cope_demo <- function() {
  env <- new.env(parent = emptyenv())
  utils::data("cope_demo", package = "copeEH", envir = env)
  d <- env$cope_demo
  s <- d$settings
  pop <- do.call(population_mortality, d$population)
  fits <- predictions <- outcomes <- truth_outcomes <- list()
  accuracy <- list()
  for (arm in c("control", "treatment")) {
    fits[[arm]] <- fit_cope(d$os[d$os$arm == arm, ], d$pfs[d$pfs$arm == arm, ],
                            pop, degree = s$degree, eta = s$eta)
    tr <- d$truth[d$truth$arm == arm, ]
    predictions[[arm]] <- predict(fits[[arm]], times = tr$time)
    outcomes[[arm]] <- state_outcomes(fits[[arm]], horizon = s$horizon,
                                       utility_pf = s$utility_pf,
                                       utility_pd = s$utility_pd,
                                       annual_cost_pf = s$annual_cost_pf,
                                       annual_cost_pd = s$annual_cost_pd,
                                       discount = s$discount, step = 0.02)
    st <- state_time_summary(tr)
    ec <- discounted_partitioned_outcomes(tr, s$utility_pf, s$utility_pd,
                                           s$annual_cost_pf, s$annual_cost_pd,
                                           s$discount)
    truth_outcomes[[arm]] <- data.frame(horizon = s$horizon, discount = s$discount,
       pf_years = unname(st["pf_time"]), pd_years = unname(st["pd_time"]),
       life_years = unname(st["total_life_years"]), qaly = unname(ec["qaly"]),
       state_cost = unname(ec["state_cost"]))
    for (endpoint in c("os", "pfs")) {
      for (interval in c("observed", "extrapolated")) {
        keep <- if (interval == "observed") tr$time <= s$cutoff else tr$time >= s$cutoff
        times <- tr$time[keep]
        err <- predictions[[arm]][[endpoint]][keep] - tr[[endpoint]][keep]
        accuracy[[length(accuracy) + 1L]] <- data.frame(arm = arm,
          endpoint = toupper(endpoint), interval = interval,
          miae = trapezoid_integral(times, abs(err)) / diff(range(times)),
          max_absolute_error = max(abs(err)))
      }
    }
  }
  economic <- compare_outcomes(outcomes$treatment, outcomes$control,
                                wtp = s$wtp, incremental_cost = s$incremental_cost)
  economic_truth <- compare_outcomes(truth_outcomes$treatment, truth_outcomes$control,
                                      wtp = s$wtp, incremental_cost = s$incremental_cost)
  metrics <- c("pf_years", "pd_years", "life_years", "qaly", "state_cost")
  errors <- do.call(rbind, lapply(names(outcomes), function(arm) {
    estimate <- unlist(outcomes[[arm]][metrics], use.names = FALSE)
    truth <- unlist(truth_outcomes[[arm]][metrics], use.names = FALSE)
    data.frame(arm = arm, outcome = metrics, truth = truth, estimate = estimate,
                 error = estimate - truth, relative_error_percent = 100 * (estimate / truth - 1))
  }))
  structure(list(data = d, fits = fits, predictions = predictions,
                  accuracy = do.call(rbind, accuracy), outcomes = errors,
                  economic = rbind(truth = economic_truth, estimate = economic)),
            class = "cope_demo_result")
}

.demo_km <- function(data) {
  times <- sort(unique(data$time[data$status == 1L]))
  risk <- vapply(times, function(t) sum(data$time >= t), numeric(1))
  events <- vapply(times, function(t) sum(data$time == t & data$status == 1L), numeric(1))
  list(time = c(0, times, max(data$time)),
       survival = c(1, cumprod(1 - events / risk), prod(1 - events / risk)))
}

plot.cope_demo_result <- function(x, type = c("survival", "state"), ...) {
  if (length(list(...))) stop("Additional plot arguments are not supported.", call. = FALSE)
  type <- match.arg(type)
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old))
  graphics::layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE),
                     heights = c(1, 1, 0.15))
  colours <- c(control = "#D55E00", treatment = "#0072B2")
  panel <- 0L
  quantities <- if (type == "survival") c("os", "pfs") else c("probability_pf", "probability_pd")
  for (quantity in quantities) {
    for (arm in c("control", "treatment")) {
      panel <- panel + 1L
      p <- x$predictions[[arm]]
      tr <- x$data$truth[x$data$truth$arm == arm, ]
      endpoint <- switch(quantity, os = "OS", pfs = "PFS",
                           probability_pf = "PF", probability_pd = "PD")
      label <- paste0("(", letters[panel], ") ", endpoint, ": ", arm)
      graphics::par(mar = c(3.7, 3.8, 2.2, 0.8), mgp = c(2.15, 0.55, 0),
                    las = 1, tcl = -0.25, cex.axis = 0.87, family = "sans")
      graphics::plot(p$time, p[[quantity]], type = "n", xlim = c(0, 20),
                     ylim = c(0, 1), xaxs = "i", yaxs = "i", bty = "l",
                     xlab = "Time (years)", ylab = "Probability")
      graphics::rect(x$data$settings$cutoff, 0, 20, 1, col = "#F3F5F7", border = NA)
      if (type == "survival") {
        obs <- x$data[[quantity]]
        km <- .demo_km(obs[obs$arm == arm, ])
        graphics::lines(km$time, km$survival, type = "s", col = "#969696", lwd = 1.2)
      }
      graphics::lines(tr$time, tr[[quantity]], col = "#242424", lty = 2, lwd = 1.8)
      graphics::lines(p$time, p[[quantity]], col = colours[[arm]], lwd = 2.1)
      graphics::abline(v = x$data$settings$cutoff, col = "#A0A0A0", lty = 3)
      graphics::title(main = label, adj = 0, line = 0.65, cex.main = 1.02, font.main = 2)
    }
  }
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  lab <- c("Known truth", "COPE-EH control", "COPE-EH treatment")
  col <- c("#242424", colours)
  lty <- c(2, 1, 1)
  if (type == "survival") {
    lab <- c("Observed KM", lab)
    col <- c("#969696", col)
    lty <- c(1, lty)
  }
  graphics::legend("center", legend = lab, col = col, lty = lty, lwd = 2,
                     horiz = TRUE, bty = "n", text.font = 2, cex = 0.85)
  invisible(x)
}
