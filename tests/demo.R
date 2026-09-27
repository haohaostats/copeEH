library(copeEH)
data(cope_demo)
stopifnot(nrow(cope_demo$os) == 4000L, nrow(cope_demo$pfs) == 4000L,
          all(cope_demo$pfs$time <= cope_demo$os$time + 1e-10),
          all(cope_demo$truth$pfs <= cope_demo$truth$os))
set.seed(17)
seed <- .Random.seed
d <- cope_example("treatment")
stopifnot(identical(seed, .Random.seed), nrow(d$os) == 2000L)
x <- run_cope_demo()
stopifnot(all(vapply(x$fits, function(f) f$converged, logical(1))),
          nrow(x$accuracy) == 8L, nrow(x$outcomes) == 10L,
          all(is.finite(x$accuracy$miae)),
          all(x$accuracy$miae >= 0),
          max(x$accuracy$miae) < 0.02,
          max(abs(x$outcomes$relative_error_percent[x$outcomes$outcome == "qaly"])) < 1,
          all(x$economic$classification == "More effective, more costly"))
f <- tempfile(fileext = ".pdf")
grDevices::pdf(f, width = 9, height = 6.8)
invisible(plot(x))
invisible(plot(x, type = "state"))
grDevices::dev.off()
stopifnot(file.info(f)$size > 0)
unlink(f)
cat("Fixed teaching dataset and demo checks passed.\n")
