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

knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6.5,
  fig.height = 4.25,
  warning  = FALSE,
  message  = FALSE
)

## ----prepare------------------------------------------------------------------
data("landscape_counts")
data("surv_pts")
pts <- vect(surv_pts)
land_rast <- rast(system.file("extdata", "landscape.tif",
                              package = "multiScaleR"))

# bin = TRUE (the default) stores distance-binned summaries for every route
inputs <- kernel_prep(pts = pts, raster_stack = land_rast, max_D = 1700,
                      kernel = "gaussian", bin = TRUE, verbose = FALSE)

# Starting model: its formula and family carry through to the fit
df    <- data.frame(landscape_counts, inputs$kernel_dat)
start <- glm(counts ~ site + land1 + land2, family = poisson(), data = df)

## ----two-routes---------------------------------------------------------------
# Default backend: refit the GLM at every candidate scale
time_refit <- system.time(
  fit_refit <- multiScale_optim(start, inputs, verbose = FALSE)
)

# Analytic backend: scales and coefficients in one likelihood
time_profile <- system.time(
  fit_profile <- multiScale_optim(start, inputs, backend = "profile_glm",
                                  n_cores = 1, verbose = FALSE)
)

# Sigma, its standard error, and the default 95% interval, in meters
est_refit   <- summary(fit_refit)$opt_scale
est_profile <- summary(fit_profile)$opt_scale
scale_table <- data.frame(
  backend   = rep(c("refit", "profile_glm"), each = 2),
  covariate = c(rownames(est_refit), rownames(est_profile)),
  rbind(est_refit, est_profile),
  check.names = FALSE
)
knitr::kable(scale_table, digits = 1, row.names = FALSE,
             col.names = c("Backend", "Covariate", "Sigma (m)", "SE",
                           "Lower", "Upper"))

# Seconds each fit took
round(c(refit = time_refit[["elapsed"]],
        profile_glm = time_profile[["elapsed"]]), 2)

## ----predict-refit------------------------------------------------------------
# refit: standard errors that hold the estimated scales fixed
se_refit <- predict(fit_refit, se.fit = TRUE, uncertainty = "conditional")
round(head(data.frame(expected = se_refit$fit, se = se_refit$se.fit), 3), 2)

## ----predict-refit-joint, error = TRUE----------------------------------------
try({
# refit cannot return a joint interval
predict(fit_refit, interval = "confidence", uncertainty = "joint")
})

## ----predict-profile----------------------------------------------------------
# profile_glm: scales and coefficients vary together
joint <- predict(fit_profile, interval = "confidence", uncertainty = "joint")
round(head(data.frame(expected = joint$fit, se = joint$se.fit,
                      lower = joint$lower, upper = joint$upper), 3), 2)

## ----spline-models------------------------------------------------------------
# Recipe 1: knots and boundaries written out
start_knots <- glm(counts ~ site + land1 +
                     splines::ns(land2, knots = c(-0.5, 0.5),
                                 Boundary.knots = c(-4, 4)),
                   family = poisson(), data = df)

# Recipe 2: R chooses the knots from the data
start_df <- glm(counts ~ site + land1 + splines::ns(land2, df = 3),
                family = poisson(), data = df)

# Written-out knots work with the analytic backend
fit_knots <- multiScale_optim(start_knots, inputs, backend = "profile_glm",
                              n_cores = 1, verbose = FALSE)

## ----spline-df-error, error = TRUE--------------------------------------------
try({
# The df recipe does not
multiScale_optim(start_df, inputs, backend = "profile_glm",
                 n_cores = 1, verbose = FALSE)
})

## ----spline-df-refit----------------------------------------------------------
# The default backend accepts either recipe
fit_df <- multiScale_optim(start_df, inputs, verbose = FALSE)

# Compare the straight-line model with the two curved ones
aic_tab(list(fit_profile, fit_knots, fit_df),
        mod_names = c("linear", "spline, fixed knots", "spline, df = 3"))

## ----mixed-types--------------------------------------------------------------
# Name each covariate and how to summarize its source layer
features <- msr_vars(
  land1  = kernel_var("land1"),                      # mean, scale estimated
  land2  = kernel_var("land2"),                      # mean, scale estimated
  rough3 = surface_var("land3", "sa", radius = 300)  # roughness within 300 m
)

inputs_mix <- kernel_prep(pts = pts, raster_stack = land_rast, max_D = 1700,
                          kernel = "gaussian", scale_vars = features,
                          bin = TRUE, verbose = FALSE)

df_mix    <- data.frame(landscape_counts, inputs_mix$kernel_dat)
start_mix <- glm(counts ~ site + land1 + land2 + rough3,
                 family = poisson(), data = df_mix)

# Either backend accepts the fixed-radius metric
fit_mix_refit   <- multiScale_optim(start_mix, inputs_mix, verbose = FALSE)
fit_mix_profile <- multiScale_optim(start_mix, inputs_mix,
                                    backend = "profile_glm",
                                    n_cores = 1, verbose = FALSE)

round(fit_mix_profile$scale_est, 1)                  # estimated scales (m)
round(coef(summary(fit_mix_profile$opt_mod)), 3)     # regression table

## ----mixed-types-values, include = FALSE--------------------------------------
beta_mix <- coef(summary(fit_mix_profile$opt_mod))

## ----bayes-prepare------------------------------------------------------------
# Validate the formula, rows, priors, and bins. Nothing is sampled.
bayes_check <- multiScale_bayes(counts ~ site + land1 + land2,
                                inputs = inputs, data = landscape_counts,
                                family = "poisson", method = "prepare")
bayes_check

## ----bayes-prepare-error, error = TRUE----------------------------------------
try({
# A term the Bayesian route rejects fails at once
multiScale_bayes(counts ~ site + land1 + splines::ns(land2, df = 3),
                 inputs = inputs, data = landscape_counts,
                 family = "poisson", method = "prepare")
})

## ----quick-reference, eval = FALSE--------------------------------------------
# library(multiScaleR)
# library(terra)
# 
# # 1. Prepare once. Binned inputs serve every route.
# data("landscape_counts"); data("surv_pts")
# pts       <- vect(surv_pts)
# land_rast <- rast(system.file("extdata", "landscape.tif",
#                               package = "multiScaleR"))
# inputs <- kernel_prep(pts = pts, raster_stack = land_rast, max_D = 1700,
#                       kernel = "gaussian", bin = TRUE)
# df     <- data.frame(landscape_counts, inputs$kernel_dat)
# start  <- glm(counts ~ site + land1 + land2, family = poisson(), data = df)
# 
# # 2. Same model, two routes
# fit_refit   <- multiScale_optim(start, inputs)               # any model class
# fit_profile <- multiScale_optim(start, inputs, backend = "profile_glm",
#                                 n_cores = 1)                 # joint intervals
# summary(fit_refit)$opt_scale       # sigma, SE, symmetric 95% interval (m)
# summary(fit_profile)$opt_scale     # sigma, SE, 95% interval on the log scale
# 
# # 3. Prediction uncertainty by route
# predict(fit_refit, se.fit = TRUE, uncertainty = "conditional")       # scales fixed
# predict(fit_profile, interval = "confidence", uncertainty = "joint") # scales uncertain
# 
# # 4. Splines: written-out knots for profile_glm, either recipe for refit
# start_knots <- glm(counts ~ site + land1 +
#                      splines::ns(land2, knots = c(-0.5, 0.5),
#                                  Boundary.knots = c(-4, 4)),
#                    family = poisson(), data = df)
# start_df    <- glm(counts ~ site + land1 + splines::ns(land2, df = 3),
#                    family = poisson(), data = df)
# fit_knots <- multiScale_optim(start_knots, inputs, backend = "profile_glm",
#                               n_cores = 1)
# fit_df    <- multiScale_optim(start_df, inputs)
# aic_tab(list(fit_profile, fit_knots, fit_df),
#         mod_names = c("linear", "spline, fixed knots", "spline, df = 3"))
# 
# # 5. Mixed covariate types: two estimated scales and one fixed-radius metric
# features <- msr_vars(
#   land1  = kernel_var("land1"),
#   land2  = kernel_var("land2"),
#   rough3 = surface_var("land3", "sa", radius = 300)
# )
# inputs_mix <- kernel_prep(pts = pts, raster_stack = land_rast, max_D = 1700,
#                           kernel = "gaussian", scale_vars = features,
#                           bin = TRUE)
# df_mix    <- data.frame(landscape_counts, inputs_mix$kernel_dat)
# start_mix <- glm(counts ~ site + land1 + land2 + rough3,
#                  family = poisson(), data = df_mix)
# fit_mix_refit   <- multiScale_optim(start_mix, inputs_mix)
# fit_mix_profile <- multiScale_optim(start_mix, inputs_mix,
#                                     backend = "profile_glm", n_cores = 1)
# 
# # 6. Check a Bayesian model without fitting it
# bayes_check <- multiScale_bayes(counts ~ site + land1 + land2,
#                                 inputs = inputs, data = landscape_counts,
#                                 family = "poisson", method = "prepare")

