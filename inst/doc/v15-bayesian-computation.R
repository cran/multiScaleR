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
library(ggplot2)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6.5,
  fig.height = 4.25,
  warning  = FALSE,
  message  = FALSE
)

# The slice-sampling runs shown below take up to two minutes, so their
# diagnostics, run times, and summaries are stored in the package.
# Regenerate them with the `bayesian-computation` section of
# data-raw/generate_vignette_cache.R.
cache_file <- file.path("..", "inst", "extdata",
                        "cache-bayesian-computation.rds")
if (!file.exists(cache_file)) {
  cache_file <- system.file("extdata", "cache-bayesian-computation.rds",
                            package = "multiScaleR")
}
cache <- readRDS(cache_file)
pilot <- cache$pilot
long  <- cache$long

## ----setup-check--------------------------------------------------------------
setup <- check_bayes_setup("sample")       # also "slice", "laplace", "marginal"
setup$packages[c("package", "minimum")]    # packages NUTS needs, minimum versions

## ----prepare------------------------------------------------------------------
example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
                               package = "multiScaleR"))
dat    <- example$data     # 100 counts from 25 sites over five years
inputs <- example$inputs   # prepared landscape inputs, distances in meters

response <- count ~ habitat + factor(year) + (1 | site)
scales   <- list(habitat = ~1 + (1 | year))   # scale of effect varies by year
priors   <- list(scale_median = 250, scale_group_sd = 0.5)

prepared <- multiScale_bayes(response, inputs, dat, scales = scales,
                             row_id = "site_year", family = "nbinom2",
                             priors = priors, method = "prepare")
prepared

## ----pilot-call, eval = FALSE-------------------------------------------------
# pilot <- multiScale_bayes(response, inputs, dat, scales = scales,
#                           row_id = "site_year", family = "nbinom2",
#                           priors = priors,
#                           method = "slice",   # sampler written in R
#                           chains = 2,         # independent runs
#                           iter = 1000,        # iterations per chain; half are warmup
#                           seed = 5212)
# 
# pilot$diagnostics$passed                # did the sampler checks pass?
# diagnose_bayes_scales(pilot)$parameters # checks on the quantities you report

## ----pilot-show, echo = FALSE-------------------------------------------------
pilot$diagnostics$passed
pilot$precision$parameters[c("quantity", "group", "rhat", "ess_bulk",
                             "ess_tail")]

## ----long-call, eval = FALSE--------------------------------------------------
# slice_long <- multiScale_bayes(response, inputs, dat, scales = scales,
#                                row_id = "site_year", family = "nbinom2",
#                                priors = priors, method = "slice",
#                                chains = 4, iter = 8000, seed = 5213,
#                                cores = 4)       # one chain per worker

## ----nuts-call, eval = FALSE--------------------------------------------------
# nuts_fit <- multiScale_bayes(response, inputs, dat, scales = scales,
#                              row_id = "site_year", family = "nbinom2",
#                              priors = priors, method = "sample",
#                              chains = 4, iter = 6000, warmup = 3000,
#                              seed = 5211,
#                              control = list(adapt_delta = 0.95,
#                                             max_treedepth = 12))

## ----sampler-table, echo = FALSE----------------------------------------------
nuts_fit  <- example$fits$pooled
nuts_info <- example$fit_info["pooled", ]
samplers <- data.frame(
  run = c("Slice pilot", "Slice, long", "NUTS"),
  draws = c(prod(pilot$chain_dim[1:2]), prod(long$chain_dim[1:2]),
            nuts_info$original_draws),
  seconds = c(pilot$seconds, long$seconds, nuts_info$seconds),
  rhat = c(max(pilot$diagnostics$parameters$rhat),
           max(long$diagnostics$parameters$rhat), nuts_info$max_rhat),
  ess = c(min(pilot$diagnostics$parameters$ess_bulk),
          min(long$diagnostics$parameters$ess_bulk), nuts_info$min_bulk_ess),
  passed = c(pilot$diagnostics$passed, long$diagnostics$passed,
             nuts_info$passed)
)
knitr::kable(samplers, digits = c(0, 0, 0, 3, 0, 0), row.names = FALSE,
             col.names = c("Run", "Retained draws", "Seconds", "Largest Rhat",
                           "Smallest bulk ESS", "Passed"))

## ----approximations, eval = FALSE---------------------------------------------
# # One Gaussian curve around the posterior peak
# fast_fit <- multiScale_bayes(response, inputs, dat, scales = scales,
#                              row_id = "site_year", family = "nbinom2",
#                              priors = priors, method = "laplace",
#                              laplace_draws = 2000, seed = 527001)
# 
# # Numerical integration over the spread of scales among years
# marginal_fit <- multiScale_bayes(response, inputs, dat, scales = scales,
#                                  row_id = "site_year", family = "nbinom2",
#                                  priors = priors, method = "marginal",
#                                  laplace_draws = 2000, seed = 527002)
# 
# # The same, reweighted against the full posterior
# corrected_fit <- multiScale_bayes(response, inputs, dat, scales = scales,
#                                   row_id = "site_year", family = "nbinom2",
#                                   priors = priors, method = "marginal",
#                                   laplace_draws = 2000, seed = 527003,
#                                   control = list(correction = "importance",
#                                                  importance_df = 5,
#                                                  importance_min = 65536,
#                                                  importance_max = 131072,
#                                                  importance_refinements = 2))
# 
# fast_fit$diagnostics        # each fit stores its own numerical checks

## ----comparison---------------------------------------------------------------
comparison <- readRDS(system.file("extdata", "bayesian-computation-example.rds",
                                  package = "multiScaleR"))
print(comparison$availability, row.names = FALSE)   # status, run time, reason

## ----sd-table, echo = FALSE---------------------------------------------------
sd_rows <- subset(comparison$quantiles, type == "scale_sd")
sd_rows <- sd_rows[match(c("nuts", "laplace", "marginal", "nested"),
                         sd_rows$method), ]
slice_sd <- long$group_sd$log_scale_group_sd
sd_table <- data.frame(
  computation = c("NUTS", "Slice, long", "Laplace", "Marginal",
                  "Marginal, corrected"),
  checks = c("Passed", "Failed", "Passed", "Failed", "Passed"),
  seconds = c(nuts_info$seconds, long$seconds,
              comparison$availability$seconds),
  median = c(sd_rows$median[1], slice_sd$median, sd_rows$median[2:4]),
  lower  = c(sd_rows$lower[1],  slice_sd$lower,  sd_rows$lower[2:4]),
  upper  = c(sd_rows$upper[1],  slice_sd$upper,  sd_rows$upper[2:4])
)
knitr::kable(sd_table, digits = c(0, 0, 1, 2, 2, 2), row.names = FALSE,
             col.names = c("Computation", "Own checks", "Seconds",
                           "Spread among years", "Lower", "Upper"))

## ----annual-plot, echo = FALSE, fig.height = 3.8, fig.cap = "***Scale of effect in each year from three computations of the same model. Points are posterior medians and bars are 95% intervals.***"----
annual <- subset(comparison$quantiles, type == "annual_log_scale" &
                   method != "marginal")
annual$year <- sub("year=", "", annual$label)
annual$computation <- factor(annual$method,
                             levels = c("nuts", "nested", "laplace"),
                             labels = c("NUTS", "Marginal, corrected",
                                        "Laplace"))
ggplot(annual, aes(x = year, y = exp(median), ymin = exp(lower),
                   ymax = exp(upper), color = computation)) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("black", "#2166AC", "#B2182B")) +
  scale_y_log10() +
  labs(x = "Year", y = "Scale of effect, sigma (m)", color = "Computation") +
  theme_classic()

## ----laplace-check, eval = FALSE----------------------------------------------
# # Probability that each year's scale is 150 m or less, and 350 m or less
# check <- diagnose_bayes_laplace(fast_fit,
#                                 thresholds = list(habitat = c(150, 350)),
#                                 tolerance  = 0.05,     # difference that matters
#                                 candidates = 65536,    # size of the fresh sample
#                                 seed       = 127)
# check$probabilities[, c("group", "threshold", "approximation", "reference",
#                         "status")]

## ----precision-call, eval = FALSE---------------------------------------------
# precision <- diagnose_bayes_scales(slice_long)
# precision$quantiles     # median and interval limits, each with its mcse

## ----precision-show, echo = FALSE---------------------------------------------
q <- subset(long$precision$quantiles, group == "year=2021")
q[c("quantity", "group", "probability", "estimate", "mcse")]

## ----probability-call, eval = FALSE-------------------------------------------
# diagnose_bayes_probabilities(slice_long,
#                              thresholds = list(habitat = c(100, 250)),  # meters
#                              mcse_target = 0.005)$probabilities

## ----probability-show, echo = FALSE-------------------------------------------
p <- subset(cache$probabilities$probabilities, group == "year=2021")
p[c("group", "threshold", "estimate", "mcse", "precision_target_met")]

## ----quick-reference, eval = FALSE--------------------------------------------
# # 1. Check the software and the model
# check_bayes_setup("sample")
# prepared <- multiScale_bayes(response, inputs, dat, scales = scales,
#                              row_id = "site_year", family = "nbinom2",
#                              priors = priors, method = "prepare")
# 
# # 2. Short pilot: time it and read the diagnostics
# pilot <- multiScale_bayes(response, inputs, dat, scales = scales,
#                           row_id = "site_year", family = "nbinom2",
#                           priors = priors, method = "slice",
#                           chains = 2, iter = 1000, seed = 5212)
# pilot$diagnostics$passed                  # FALSE: extend or change sampler
# diagnose_bayes_scales(pilot)$parameters   # Rhat and ESS for reported quantities
# 
# # 3. Full run with NUTS
# nuts_fit <- multiScale_bayes(response, inputs, dat, scales = scales,
#                              row_id = "site_year", family = "nbinom2",
#                              priors = priors, method = "sample",
#                              chains = 4, iter = 6000, warmup = 3000,
#                              seed = 5211,
#                              control = list(adapt_delta = 0.95,
#                                             max_treedepth = 12))
# summary(nuts_fit)$diagnostics             # passed, divergences, Rhat, ESS
# 
# # 4. Precision of what you report (needs the original chains)
# diagnose_bayes_scales(nuts_fit)$quantiles # interval limits with their mcse
# diagnose_bayes_probabilities(nuts_fit,
#                              thresholds = list(habitat = c(100, 250)),
#                              mcse_target = 0.005)$probabilities
# 
# # 5. A fast approximation. Compare it with a sampler before you report it
# fast_fit <- multiScale_bayes(response, inputs, dat, scales = scales,
#                              row_id = "site_year", family = "nbinom2",
#                              priors = priors, method = "laplace",
#                              laplace_draws = 2000, seed = 527001)
# fast_fit$diagnostics                      # numerical checks only
# diagnose_bayes_laplace(fast_fit, thresholds = list(habitat = c(150, 350)),
#                        candidates = 65536, seed = 127)$probabilities
#                                           # "unresolved" means not checked

