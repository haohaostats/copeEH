library(copeEH)

expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}

pop <- population_mortality(c(0, 1, 2), c(0.01, 0.02, 0.04))
stopifnot(isTRUE(all.equal(pop$cumhazard(c(0, 0.5, 1, 2, 3)),
                          c(0, 0.005, 0.01, 0.03, 0.07))))
stopifnot(identical(pop$hazard(c(0, 1, 2)), c(0.01, 0.02, 0.04)))
stopifnot(isTRUE(all.equal(pop$survival(3), exp(-0.07))))
expect_error(population_mortality(c(0, 0), c(0.1, 0.2)))
expect_error(population_mortality(hazard = -1))
expect_error(pop$survival(-1))

set.seed(42)
rng <- .Random.seed
d <- cope_example()
stopifnot(identical(rng, .Random.seed), identical(d$os, cope_example()$os))
rm(".Random.seed", envir = .GlobalEnv)
invisible(cope_example())
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
fit <- fit_cope(d$os, d$pfs, d$population)
stopifnot(inherits(fit, "cope_fit"), fit$converged)
stopifnot(length(coef(fit)) == 8L, length(coef(fit, free = TRUE)) == 6L)
stopifnot(inherits(summary(fit), "summary.cope_fit"))
invisible(capture.output(print(fit)))
t <- c(seq(0, 40, by = 0.1), 100, 1000, 1e6)
p <- predict(fit, t)
stopifnot(all(p$pfs <= p$os + 1e-10), all(p$os <= 1), all(p$pfs >= 0),
          all(diff(p$os) <= 1e-10), all(diff(p$pfs) <= 1e-10),
          all(p$probability_pd >= -1e-10),
          max(abs(p$probability_pf + p$probability_pd + p$probability_death - 1)) < 1e-12)
stopifnot(p$os[1] == 1, p$pfs[1] == 1)
stopifnot(ncol(predict(fit, 0, type = "state")) == 4,
          ncol(predict(fit, c(3, 0, 3), type = "survival")) == 3,
          identical(predict(fit, c(3, 0, 3))$time, c(3, 0, 3)))
expect_error(predict(fit, 1, se.fit = TRUE))
expect_error(predict(fit, -1))
bad <- fit
bad$converged <- FALSE
expect_error(predict(bad, 1))
expect_error(fit_cope(d$os[FALSE, ], d$pfs, d$population))
expect_error(fit_cope(transform(d$os, status = 2), d$pfs, d$population))
expect_error(fit_cope(d$os, d$pfs, list()))
expect_error(fit_cope(d$os, d$pfs, d$population, eta = 0))
expect_error(fit_cope(d$os, d$pfs, d$population, degree = 1, tail = "no_cure"))

nc <- fit_cope(d$os, d$pfs, d$population, tail = "no_cure")
stopifnot(nc$converged, coef(nc)["os_theta_3"] == 0,
          coef(nc)["pfs_theta_3"] == 0)
pn <- predict(nc, seq(0, 40, by = 0.1))
stopifnot(all(pn$pfs <= pn$os + 1e-10))

obj <- fit
obj$degree <- 1L
obj$eta <- 2
obj$tail <- "free"
obj$coefficients <- c(os_theta_0 = 1, os_theta_1 = 1,
                      pfs_theta_0 = 1, pfs_theta_1 = 0)
obj$population <- population_mortality(hazard = 0)
out <- state_outcomes(obj, horizon = 5, discount = 0, step = 0.001,
                       utility_pf = 0.8, utility_pd = 0.6,
                       annual_cost_pf = 10, annual_cost_pd = 20)
pf <- 2 * log(1 + 5 / 2)
stopifnot(abs(out$life_years - 5) < 1e-10,
          abs(out$pf_years - pf) < 1e-6,
          abs(out$pd_years - (5 - pf)) < 1e-6,
          abs(out$qaly - (0.8 * pf + 0.6 * (5 - pf))) < 1e-6,
          abs(out$state_cost - (10 * pf + 20 * (5 - pf))) < 1e-5)
expect_error(state_outcomes(fit, horizon = 0))
expect_error(state_outcomes(fit, horizon = 5, step = 0))
expect_error(state_outcomes(fit, horizon = 5, utility_pf = c(0.7, 0.8)))
partial <- state_outcomes(obj, horizon = 0.015, step = 0.02, discount = 0)
stopifnot(abs(partial$life_years - 0.015) < 1e-12)
same <- compare_outcomes(out, out, wtp = 100000)
stopifnot(same$delta_qaly == 0, same$inmb == 0, is.na(same$icer))
trt <- out
trt$qaly <- trt$qaly - 1
dom <- compare_outcomes(trt, out, wtp = 100000, incremental_cost = 100)
stopifnot(dom$classification == "Dominated", is.na(dom$icer), dom$inmb == -100100)
trt$qaly <- out$qaly + 1
ce <- compare_outcomes(trt, out, wtp = 100000, incremental_cost = 50000)
stopifnot(abs(ce$icer - 50000) < 1e-8, abs(ce$inmb - 50000) < 1e-8)
trt$horizon <- 8
expect_error(compare_outcomes(trt, out, wtp = 100000))

path <- tempfile(fileext = ".rds")
saveRDS(fit, path)
restored <- readRDS(path)
stopifnot(identical(predict(fit, c(0, 5)), predict(restored, c(0, 5))))
unlink(path)
graphics_file <- tempfile(fileext = ".pdf")
grDevices::pdf(graphics_file)
invisible(plot(fit))
invisible(plot(fit, type = "state"))
grDevices::dev.off()
stopifnot(file.info(graphics_file)$size > 0)
unlink(graphics_file)
cat("All copeEH regression checks passed.\n")
