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
library(terra)
library(ggplot2)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6.5,
  fig.height = 4.25,
  warning  = FALSE,
  message  = FALSE
)

## ----fit-profile--------------------------------------------------------------
data("landscape_counts")
data("surv_pts")
pts <- vect(surv_pts)
land_rast <- rast(system.file("extdata", "landscape.tif",
                              package = "multiScaleR"))

# bin = TRUE stores the distance-binned summaries the analytic backend needs
inputs <- kernel_prep(pts = pts, raster_stack = land_rast, max_D = 1700,
                      kernel = "gaussian", bin = TRUE, verbose = FALSE)

df    <- data.frame(landscape_counts, inputs$kernel_dat)
start <- glm(counts ~ site + land1 + land2, family = poisson(), data = df)

# Estimate scales and coefficients in one likelihood
fit <- multiScale_optim(start, inputs, backend = "profile_glm",
                        n_cores = 1, verbose = FALSE)

fit$scale_est                                         # sigma and SE, in meters
diagnostics(fit)$profile_glm$uncertainty_available    # TRUE: joint intervals exist

## ----two-intervals------------------------------------------------------------
# Holds the estimated scales fixed
cond  <- predict(fit, interval = "confidence", uncertainty = "conditional")

# Lets scales and coefficients vary together
joint <- predict(fit, interval = "confidence", uncertainty = "joint")

compare <- data.frame(
  expected    = cond$fit,
  cond_lower  = cond$lower,  cond_upper  = cond$upper,
  joint_lower = joint$lower, joint_upper = joint$upper
)
round(head(compare, 5), 2)

## ----two-intervals-plot, fig.cap = "***Expected counts at ten survey points with conditional (black) and joint (blue) 95% confidence intervals.***"----
# Ten survey points spread across the range of expected counts
keep <- order(compare$expected)[seq(5, 95, by = 10)]
long <- data.frame(
  point    = factor(rep(seq_along(keep), 2)),
  expected = rep(compare$expected[keep], 2),
  lower    = c(compare$cond_lower[keep], compare$joint_lower[keep]),
  upper    = c(compare$cond_upper[keep], compare$joint_upper[keep]),
  interval = rep(c("Conditional", "Joint"), each = length(keep))
)

ggplot(long, aes(x = point, y = expected, ymin = lower, ymax = upper,
                 color = interval)) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c(Conditional = "black", Joint = "#2166AC")) +
  labs(x = "Survey point (ordered by expected count)",
       y = "Expected count", color = "Interval") +
  theme_classic()

## ----width-ratio--------------------------------------------------------------
# How much wider is the joint interval at each of the 100 survey points?
ratio <- (joint$upper - joint$lower) / (cond$upper - cond$lower)
round(quantile(ratio, c(0, 0.5, 1)), 2)

## ----scale-inference, eval = requireNamespace("MASS", quietly = TRUE)---------
# One landscape layer, so the model has a single scale of effect
inputs2 <- kernel_prep(pts = pts, raster_stack = land_rast["land2"],
                       max_D = 1700, kernel = "gaussian", bin = TRUE,
                       verbose = FALSE)
df2    <- data.frame(landscape_counts, inputs2$kernel_dat)
start2 <- MASS::glm.nb(counts ~ site + land2, data = df2)

# Estimate sigma, the coefficients, and the negative binomial theta together
fit_nb <- multiScale_optim(start2, inputs2, backend = "joint_nb2",
                           n_cores = 1, verbose = FALSE)

# Refit the model at 21 sigma values between the search limits
scan <- scale_inference(fit_nb, grid_points = 21)
scan$diagnostics$profile_status     # how many regions, and do they touch a limit
round(scan$profile_set[, c("lower", "upper")])   # supported sigma, in meters

## ----expected-versus-observed-------------------------------------------------
annual <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
                              package = "multiScaleR"))
bfit <- annual$fits$pooled          # negative binomial fit, sampled earlier

# Three site-years to predict. Remove the response, which is not an input.
surveys <- annual$data[1:3, ]
surveys$count <- NULL

# Expected count: uncertainty about the mean only
expected <- predict(bfit, surveys, kernel_inputs = annual$inputs,
                    row_id = "site_year", draws = 400)

# Replicated counts: adds observation noise to every posterior draw
replicated <- simulate(bfit, nsim = 400, seed = 11, newdata = surveys,
                       kernel_inputs = annual$inputs, row_id = "site_year")
observed <- t(apply(as.matrix(replicated), 1, quantile,
                    probs = c(0.025, 0.5, 0.975)))

both <- data.frame(
  site_year       = expected$point_id,
  expected_median = expected$median,
  expected_lower  = expected$lower,  expected_upper = expected$upper,
  observed_median = observed[, 2],
  observed_lower  = observed[, 1],   observed_upper = observed[, 3]
)
knitr::kable(both, digits = 1, row.names = FALSE,
             col.names = c("Site and year", "Expected count", "Lower", "Upper",
                           "Recorded count", "Lower", "Upper"))

## ----prediction-record--------------------------------------------------------
record <- attr(expected, "prediction_contract")
record$quantity                         # what the numbers describe
record$uncertainty$observation_noise    # FALSE: observation noise is left out

## ----quick-reference, eval = FALSE--------------------------------------------
# # 1. Fit with a backend that estimates scales and coefficients together
# inputs <- kernel_prep(pts, land_rast, max_D = 1700, kernel = "gaussian",
#                       bin = TRUE)
# start  <- glm(counts ~ site + land1 + land2, family = poisson(),
#               data = data.frame(landscape_counts, inputs$kernel_dat))
# fit    <- multiScale_optim(start, inputs, backend = "profile_glm", n_cores = 1)
# diagnostics(fit)$profile_glm$uncertainty_available   # TRUE: joint intervals exist
# 
# # 2. Scale of effect and its interval
# summary(fit)                   # sigma and effective distance, in map units
# summary(fit, profile = TRUE)   # profile-likelihood limits for sigma
# 
# # One-scale negative binomial fit: every supported region of sigma
# scan <- scale_inference(fit_nb, grid_points = 21)
# scan$diagnostics$profile_status   # "localized", "disconnected", "bound_limited", "unresolved"
# scan$profile_set                  # supported regions, in map units
# 
# # 3. Expected response with two kinds of interval
# predict(fit, interval = "confidence", uncertainty = "conditional")  # scales fixed
# predict(fit, interval = "confidence", uncertainty = "joint")        # scales uncertain
# 
# # 4. Bayesian fit: expected response versus a count you might record
# predict(bfit, newdata, kernel_inputs = new_inputs, row_id = "site_year")
# simulate(bfit, nsim = 400, newdata = newdata, kernel_inputs = new_inputs,
#          row_id = "site_year")                     # adds observation noise
# attr(predict(bfit, newdata, kernel_inputs = new_inputs, row_id = "site_year"),
#      "prediction_contract")$quantity               # what the numbers describe

