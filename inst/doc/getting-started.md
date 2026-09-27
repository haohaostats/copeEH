# Getting started with copeEH

## 1. Prepare marginal endpoints

For each treatment group prepare two data frames, each with time and status.
Status is 1 for an event and 0 for right censoring. An OS event is death; a PFS
event is progression or death. Original patient identifiers are not required.
Times must be finite and non-negative. Empty datasets are rejected.

The following loads fixed teaching data without changing your RNG state:

```r
library(copeEH)
d <- cope_example(arm = "control")
head(d$os)
head(d$pfs)
```

This example is synthetic, not reconstructed KEYNOTE-024 data. Both endpoints
originate from the same latent patients; the fit uses their marginal records.
There are 2,000 patients per arm with a common ten-year censoring time. No
simulation-generation code is included in the package.

To inspect both arms and their known truth:

```r
data(cope_demo)
result <- run_cope_demo()
result$accuracy
result$outcomes
result$economic
plot(result)
plot(result, type = "state")
```

The plotting methods place panel labels above the upper-left corners, legends
below the panels and use open axes. Shading identifies the unobserved tail.
For portable vector output:

```r
pdf("demo_survival.pdf", width = 9, height = 6.8, useDingbats = FALSE)
plot(result)
dev.off()
```

This matched-model example deliberately uses substantial sample size and
follow-up. It does not establish general superiority or guarantee similarly
small errors in a clinical application. At 10--20 years, MIAE ranges from
0.0014 to 0.0060. QALY errors are below 0.15% in each arm, but control PD-time
error is approximately -6.1%. The true versus estimated incremental QALY is
1.1847 versus 1.1932; the ICER is 105,829 versus 113,422 currency units/QALY.
The INMB is -6,906 versus -16,016 at WTP 100,000. All errors are shown rather
than selecting only favourable outputs. Neither truth curves nor true
coefficients are supplied to the optimiser.

## 2. Specify matched population survival and hazard

```r
pop <- population_mortality(breaks = c(0, 5, 10),
                            hazard = c(0.01, 0.02, 0.04))
pop$survival(c(0, 5, 10, 15))
pop$hazard(c(0, 5, 10, 15))
```

Hazards switch at the specified starts, and the final rate persists indefinitely.
For annual death probabilities use -log1p(-qx), with appropriate attained-age
breaks. For a heterogeneous cohort, construct cohort survival and its matching
survivor-weighted hazard externally. Do not independently average hazards and
survival functions. Custom survival(t) and hazard(t) functions are accepted in
a list but must represent the same process. Their consistency is the caller's
responsibility.

## 3. Fit and inspect

```r
fit <- fit_cope(d$os, d$pfs, population = d$population,
                degree = 3, eta = 2, tail = "free")
summary(fit)
coef(fit)
```

COPE-EH estimates OS and PFS Bernstein coefficients under monotonicity and
coefficient-wise ordering constraints. The composite likelihood treats endpoint
contributions as working independent. It does not identify transition hazards,
recover patient linkage or supply valid joint uncertainty from the inverse
Hessian alone. The public API deliberately reports point estimates only.

The free tail estimates the last relative-survival coefficients. The no_cure
option fixes both to zero and requires degree >= 2. All-cause survival also
contains population mortality; a positive relative-survival tail is not an
immortal all-cause-survival fraction. Eta is a time scale, not a hazard ratio.

## 4. Extrapolate and plot

```r
times <- seq(0, 40, by = 0.1)
p <- predict(fit, times = times)
head(p)
plot(fit, times = times, type = "survival")
plot(fit, times = times, type = "state")
```

PF = PFS, PD = OS - PFS, and death = 1 - OS. These are probabilities obtained
from marginal survival, not fitted transition hazards. Prediction retains the
requested time order; plotting sorts and deduplicates times. The population
functions must remain valid throughout the requested extrapolation horizon.

## 5. Compute restricted state times and economic results

```r
out <- state_outcomes(fit, horizon = 20, utility_pf = 0.8, utility_pd = 0.65,
                      annual_cost_pf = 20000, annual_cost_pd = 100000,
                      discount = 0.03, step = 0.02)
out
```

PF years, PD years and life-years are undiscounted. QALYs and state costs use
the effective annual discount factor (1 + discount)^(-time). The exact horizon
is included in trapezoidal integration. Check numerical sensitivity by reducing
step. The input times and eta must be in years for these annual quantities.

For two actual treatment groups, fit separately and pass both outcomes to
compare_outcomes. The following self-comparison simply demonstrates additional
incremental costs:

```r
compare_outcomes(out, out, wtp = 100000, incremental_cost = 5000)
```

Incremental costs supplied here must already be discounted. The comparison
returns treatment minus control and INMB = WTP * delta QALY - delta cost.
An ICER is NA in dominance quadrants or with a zero increment. A positive ratio
with both increments negative measures savings per QALY lost. Read the quadrant
classification together with the signed increments. These are deterministic
estimates, not cost-effectiveness probabilities.

## 6. Save and restore

```r
saveRDS(fit, "my_cope_fit.rds")
restored <- readRDS("my_cope_fit.rds")
predict(restored, times = c(0, 5, 10), type = "survival")
```

The package mortality helper creates self-contained closures that survive
serialization. Custom functions depending on global variables or external files
require those dependencies to remain available after loading.
