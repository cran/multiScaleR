## ----cran-check, include = FALSE----------------------------------------------
# CRAN limits the whole package check to about 10 minutes, so R CMD check
# stops this vignette here, before the package loads, instead of running its
# code. R CMD build, pkgdown, r-universe, and CI checks with NOT_CRAN=true run
# the full vignette, and the installed vignette always shows the full results.
if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")) &&
    !identical(Sys.getenv("NOT_CRAN"), "true")) {
  knitr::knit_exit()
}

## ----setup, include = FALSE---------------------------------------------------
library(multiScaleR)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6.5,
  fig.height = 4.25,
  warning  = FALSE,
  message  = FALSE
)

# The slice-sampling fit in the first example takes 20 to 40 seconds, so it is
# stored in the package. The text shows the real call. Regenerate the stored
# result with the `bayesian-scales` section of data-raw/generate_vignette_cache.R.
cache_file <- file.path("..", "inst", "extdata", "cache-bayesian-scales.rds")
if (!file.exists(cache_file)) {
  cache_file <- system.file("extdata", "cache-bayesian-scales.rds",
                            package = "multiScaleR")
}
scales_cache <- readRDS(cache_file)

## ----setup-check--------------------------------------------------------------
check_bayes_setup(method = "slice")   # are the packages for this method ready?

## ----prepare-inputs-----------------------------------------------------------
data("landscape_counts")              # counts and one site-level covariate
data("surv_pts")                      # survey locations (an sf object)
pts <- terra::vect(surv_pts)          # convert to a terra SpatVector
land_rast <- terra::rast(system.file("extdata", "landscape.tif",
                                     package = "multiScaleR"))

# bin = TRUE (the default) stores the distance-binned summaries the model uses
inputs <- kernel_prep(pts = pts, raster_stack = land_rast, max_D = 1700,
                      kernel = "gaussian", bin = TRUE, verbose = FALSE)

## ----prior-in-meters----------------------------------------------------------
# Central 95% of the prior for sigma: median 300 m, SD of 1 on the log scale
prior_sigma <- 300 * exp(c(lower = -1.96, median = 0, upper = 1.96) * 1)

# The same values as effective distances (radius holding 90% of kernel weight)
prior_dist <- sapply(prior_sigma, function(s) {
  kernel_dist(sigma = s, kernel = "gaussian")
})
round(rbind(sigma = prior_sigma, effective_distance = prior_dist))

## ----prior-check--------------------------------------------------------------
prepared <- multiScale_bayes(counts ~ site + land1 + land2, inputs,
                             landscape_counts, family = "poisson",
                             priors = list(scale_median = c(land1 = 300,
                                                            land2 = 300)),
                             method = "prepare")

prior <- simulate_bayes_prior(prepared, nsim = 1000, seed = 63)

# Prior for sigma in meters, one column per landscape layer
round(sapply(prior$log_scales, function(x) {
  exp(quantile(x, c(0.025, 0.5, 0.975)))
}))

# Largest count in each simulated data set (the largest observed count is 10)
largest <- vapply(prior$responses, max, numeric(1))
signif(quantile(largest, c(0.025, 0.5, 0.975)), 2)

## ----fit-shared, eval = FALSE-------------------------------------------------
# fit <- multiScale_bayes(counts ~ site + land1 + land2, inputs,
#                         landscape_counts, family = "poisson",
#                         priors = list(scale_median = c(land1 = 300,
#                                                        land2 = 300)),
#                         method = "slice",   # four chains of 2000 iterations
#                         seed = 5211)

## ----load-fit-shared, include = FALSE-----------------------------------------
fit <- scales_cache$fit
precision <- scales_cache$precision

## ----check-shared-call, eval = FALSE------------------------------------------
# precision <- diagnose_bayes_scales(fit)   # needs the posterior package

## ----check-shared-------------------------------------------------------------
summary(fit)$diagnostics$passed      # did the parameter checks pass?
summary(fit)$diagnostics$criteria    # the rule the package applies

# Agreement among chains and effective sample sizes for each log sigma
precision$parameters[c("covariate", "rhat", "ess_bulk", "ess_tail")]

# Each interval limit on the log scale, with its Monte Carlo standard error
precision$quantiles[c("covariate", "probability", "estimate", "mcse")]
precision$passed                     # parameter checks and scale checks together

## ----trace, eval = requireNamespace("bayesplot", quietly = TRUE), fig.height = 3, fig.cap = "***Trace plots for the two scale parameters, one line per chain.***"----
# Chain plots take the stored chains of a slice fit
bayesplot::mcmc_trace(fit$chains, regex_pars = "^scale_beta")

## ----read-shared--------------------------------------------------------------
scales <- posterior_scales(fit)       # sigma in map units (meters here)
scales[c("covariate", "median", "lower", "upper",
         "p_below_min_D", "p_above_max_D")]

# Posterior medians as effective distances
dist90 <- sapply(scales$median, function(s) {
  kernel_dist(sigma = s, kernel = "gaussian")
})
round(setNames(dist90, scales$covariate))

# Regression coefficients on the log scale
coefs <- summary(fit)$coefficients
coefs[c("coefficient", "median", "lower", "upper")]

## ----annual-data--------------------------------------------------------------
annual <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
                              package = "multiScaleR"))
dat <- annual$data            # one row per site and year
inputs_annual <- annual$inputs    # prepared with kernel_prep_by_group()
head(dat, 3)
annual$truth                  # sigma used to simulate each year (m), surveys

## ----group-sd-prior-----------------------------------------------------------
# Prior on tau: half-normal with scale 0.5 (scale_group_sd)
tau <- c(median = qnorm(0.75), upper_95 = qnorm(0.975)) * 0.5
round(rbind(tau = tau, multiplier = exp(tau)), 2)

## ----prepare-annual-----------------------------------------------------------
prepared_annual <- multiScale_bayes(
  count ~ habitat + factor(year) + (1 | site), inputs_annual, dat,
  scales = list(habitat = ~1 + (1 | year)),
  row_id = "site_year",             # column that links rows to landscapes
  family = "nbinom2",               # negative binomial counts
  priors = list(scale_median = 250, scale_group_sd = 0.5),
  method = "prepare"
)
prepared_annual

info <- prepared_annual$inference_metadata
unlist(info[c("target", "method", "population_priors")])
info$uncertainty$intervals

## ----sample-annual, eval = FALSE----------------------------------------------
# bfit <- multiScale_bayes(
#   count ~ habitat + factor(year) + (1 | site), inputs_annual, dat,
#   scales = list(habitat = ~1 + (1 | year)),
#   row_id = "site_year", family = "nbinom2",
#   priors = list(scale_median = 250, scale_group_sd = 0.5),
#   method = "sample",                # NUTS
#   chains = 4, iter = 6000, warmup = 3000, seed = 5211,
#   control = list(adapt_delta = 0.95, max_treedepth = 12)
# )

## ----load-annual--------------------------------------------------------------
bfit <- annual$fits$pooled          # the stored result of the call above

d <- summary(bfit)$diagnostics      # checks computed on the full chains
data.frame(passed = d$passed,
           max_rhat = max(d$parameters$rhat),
           min_bulk_ess = min(d$parameters$ess_bulk),
           min_tail_ess = min(d$parameters$ess_tail),
           divergences = d$divergences,
           depth_hits = d$max_treedepth_hits,
           min_ebfmi = min(d$ebfmi))

## ----read-annual--------------------------------------------------------------
posterior_scales(bfit)[c("group", "median", "lower", "upper", "p_below_min_D")]

summary(bfit)$group_sd$log_scale_group_sd   # tau, spread of log sigma among years
summary(bfit)$group_sd$response_group_sd    # spread of site intercepts, log scale

## ----annual-values, include = FALSE-------------------------------------------
ann <- posterior_scales(bfit)
tau_hat <- summary(bfit)$group_sd$log_scale_group_sd

## ----comparison-plot, fig.cap = "***Posterior sigma by year under three scale models. Crosses mark the values used to simulate the data.***"----
comparison <- annual$comparison     # posterior_scales() of the three stored fits

ggplot2::ggplot(comparison, ggplot2::aes(x = year, y = median, ymin = lower,
                                         ymax = upper, color = model)) +
  ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.6)) +
  ggplot2::geom_point(ggplot2::aes(y = truth), color = "black", shape = 4,
                      size = 3) +
  ggplot2::scale_color_brewer(palette = "Dark2") +
  ggplot2::labs(x = "Year", y = "Sigma (m)", color = "Scale model") +
  ggplot2::theme_classic()

## ----comparison-values, include = FALSE---------------------------------------
cmp <- function(model, year) {
  round(comparison$median[comparison$model == model & comparison$year == year])
}

## ----prior-sensitivity--------------------------------------------------------
annual$sensitivity      # posterior of tau under three priors

## ----predictive-check---------------------------------------------------------
replicated <- simulate(bfit, nsim = 200, seed = 94)   # one column per data set

# Three features of a count data set
features <- function(y) {
  c(mean = mean(y), zero_fraction = mean(y == 0), maximum = max(y))
}
sims <- vapply(replicated, features, numeric(3))

data.frame(observed = features(dat$count),
           lower    = apply(sims, 1, quantile, 0.025),
           median   = apply(sims, 1, median),
           upper    = apply(sims, 1, quantile, 0.975))

## ----quick-reference, eval = FALSE--------------------------------------------
# # 1. Prepare binned inputs and look at the priors
# inputs <- kernel_prep(pts = pts, raster_stack = land_rast, max_D = 1700,
#                       kernel = "gaussian", bin = TRUE, verbose = FALSE)
# prepared <- multiScale_bayes(counts ~ site + land1 + land2, inputs,
#                              landscape_counts, family = "poisson",
#                              priors = list(scale_median = c(land1 = 300,
#                                                             land2 = 300)),
#                              method = "prepare")
# prior <- simulate_bayes_prior(prepared, nsim = 1000, seed = 63)
# kernel_dist(sigma = 300, kernel = "gaussian")   # a sigma as effective distance
# 
# # 2. Fit one scale per layer, check the fit, then read it
# fit <- multiScale_bayes(counts ~ site + land1 + land2, inputs,
#                         landscape_counts, family = "poisson",
#                         priors = list(scale_median = c(land1 = 300,
#                                                        land2 = 300)),
#                         method = "slice", seed = 5211)
# summary(fit)$diagnostics$passed           # parameter checks
# precision <- diagnose_bayes_scales(fit)   # scale checks and MCSE
# bayesplot::mcmc_trace(fit$chains, regex_pars = "^scale_beta")
# posterior_scales(fit)                     # sigma in map units
# summary(fit)$coefficients                 # coefficients on the link scale
# 
# # 3. Let the scale differ among years (method = "prepare" checks it first)
# bfit <- multiScale_bayes(
#   count ~ habitat + factor(year) + (1 | site), inputs_annual, dat,
#   scales = list(habitat = ~1 + (1 | year)),
#   row_id = "site_year", family = "nbinom2",
#   priors = list(scale_median = 250, scale_group_sd = 0.5),
#   method = "sample",
#   chains = 4, iter = 6000, warmup = 3000, seed = 5211,
#   control = list(adapt_delta = 0.95, max_treedepth = 12)
# )
# summary(bfit)$diagnostics                   # R-hat, ESS, divergences, E-BFMI
# posterior_scales(bfit)                      # sigma for each year
# summary(bfit)$group_sd$log_scale_group_sd   # tau, spread of log sigma
# summary(bfit)$group_sd$response_group_sd    # spread of site intercepts
# replicated <- simulate(bfit, nsim = 200, seed = 94)   # predictive check

