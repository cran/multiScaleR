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
# Draw raster maps with 32 colors rather than terra's default 100. The maps
# look the same, and their stored images are about 40% smaller.
options(terra.pal = terra::map.pal("viridis", 32))
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

## ----load-data----------------------------------------------------------------
data("landscape_counts")   # data frame: counts and one site-level covariate
dat <- landscape_counts

data("surv_pts")           # survey locations (an sf object from the package)
pts <- vect(surv_pts)      # convert to a terra SpatVector

# Three landscape layers on one grid, in meters
land_rast <- rast(system.file("extdata", "landscape.tif",
                              package = "multiScaleR"))

## ----plot-inputs, fig.height = 3.4, fig.cap = "***The continuous land3 layer, alone and with the survey points.***"----
par(mfrow = c(1, 2))
plot(land_rast$land3, main = "land3")
plot(land_rast$land3, main = "Survey points")
points(pts, pch = 19, cex = 0.5, col = "red")
par(mfrow = c(1, 1))

## ----define-vars--------------------------------------------------------------
scale_vars <- msr_vars(
  land1_prop   = kernel_var("land1"),     # proportion of habitat, sigma estimated
  land2_mean   = kernel_var("land2"),     # mean of land2, sigma estimated
  land3_sq_500 = surface_var("land3", metric = "sq",
                             radius = 500)   # roughness of land3 within 500 m
)

scale_vars[, c("covariate", "source", "type", "metric", "radius", "optimize")]

## ----kernel-prep--------------------------------------------------------------
kernel_inputs <- kernel_prep(
  pts          = pts,          # survey locations
  raster_stack = land_rast,    # source layers named in scale_vars
  max_D        = 1700,         # largest distance (m) considered
  kernel       = "gaussian",   # kernel for the two kernel covariates
  scale_vars   = scale_vars,
  verbose      = FALSE
)

# Mean and standard deviation of each covariate, in its own units
round(sapply(kernel_inputs$scl_params, function(x) x), 3)

## ----fit-model----------------------------------------------------------------
df <- data.frame(dat, kernel_inputs$kernel_dat)

# Starting model: formula and family carry through to the optimizer
mod <- glm(counts ~ site + land1_prop + land2_mean + land3_sq_500,
           family = poisson(),
           data   = df)

# Estimate sigma for land1_prop and land2_mean
opt <- multiScale_optim(fitted_mod    = mod,
                        kernel_inputs = kernel_inputs,
                        verbose       = FALSE)

## ----summary------------------------------------------------------------------
summary(opt)

## ----read-scale, include = FALSE----------------------------------------------
est    <- summary(opt)$opt_scale     # sigma, SE, and 95% interval (m)
dist90 <- kernel_dist(opt)           # radius holding 90% of kernel weight (m)

## ----coef-units---------------------------------------------------------------
beta <- coef(opt$opt_mod)              # log scale, per one standard deviation
ci   <- confint.default(opt$opt_mod)   # Wald 95% limits, log scale
sds  <- opt$scl_params$sd              # one standard deviation, in metric units

# A step that is easy to picture, in the units of the metric
step <- c(land3_sq_500 = 0.01)         # 0.01 more roughness, in land3 units
m <- names(step)

effect <- data.frame(
  Covariate  = m,
  SD         = sds[m],
  Step       = step,
  Multiplier = exp(beta[m]  * step / sds[m]),   # change in expected count
  Lower      = exp(ci[m, 1] * step / sds[m]),
  Upper      = exp(ci[m, 2] * step / sds[m])
)
knitr::kable(effect, digits = 3, row.names = FALSE)

## ----sim-truth----------------------------------------------------------------
# Kernel-weighted roughness: cells are weighted by the model's kernel
weighted_vars <- msr_vars(
  land3_sq = surface_var("land3", metric = "sq", weighted = TRUE)
)

# Roughness at a sigma of 300 m, the scale used to simulate the counts
truth_inputs <- kernel_prep(
  pts          = pts,
  raster_stack = land_rast["land3"],
  max_D        = 1500,
  kernel       = "gaussian",
  scale_vars   = weighted_vars,
  sigma        = 300,              # compute the covariate at this sigma
  verbose      = FALSE
)

set.seed(1)
# kernel_dat is standardized, so 0.5 is the effect of one standard deviation
sim_counts <- rpois(100, exp(0.5 + 0.5 * truth_inputs$kernel_dat$land3_sq))

## ----sim-fit------------------------------------------------------------------
free_inputs <- kernel_prep(
  pts          = pts,
  raster_stack = land_rast["land3"],
  max_D        = 1500,             # sigma is searched within this distance
  kernel       = "gaussian",
  scale_vars   = weighted_vars,
  verbose      = FALSE
)

sim_df  <- data.frame(counts = sim_counts, free_inputs$kernel_dat)
sim_mod <- glm(counts ~ land3_sq, family = poisson(), data = sim_df)

sim_opt <- multiScale_optim(fitted_mod    = sim_mod,
                            kernel_inputs = free_inputs,
                            verbose       = FALSE)

# Sigma with its standard error and 95% interval, in meters
sim_scale <- summary(sim_opt)$opt_scale
round(sim_scale, 1)

# Radius that holds 90% of the kernel weight, in meters
sim_dist <- kernel_dist(sim_opt)
round(sim_dist)

## ----sim-profile, fig.height = 3.6, fig.cap = "***AICc across candidate sigma values for kernel-weighted roughness. The dashed line marks the estimate.***"----
# Refit the model at 15 candidate sigma values
plot(profile_sigma(sim_opt, n_pts = 15, verbose = FALSE))

## ----sim-coef-----------------------------------------------------------------
sim_beta <- coef(sim_opt$opt_mod)["land3_sq"]             # per standard deviation
sim_ci   <- confint.default(sim_opt$opt_mod)["land3_sq", ]
sim_sd   <- sim_opt$scl_params$sd["land3_sq"]             # land3 units, fitted sigma

# Change in expected count for 0.01 more roughness
sim_effect <- exp(c(Estimate = unname(sim_beta), Lower = sim_ci[[1]],
                    Upper = sim_ci[[2]]) * 0.01 / sim_sd)
round(sim_effect, 3)

# The value used in the simulation, on the same footing
round(exp(0.5 * 0.01 / truth_inputs$scl_params$sd["land3_sq"]), 3)

## ----sim-marginal, fig.height = 3.6, fig.cap = "***Expected count across the observed range of kernel-weighted roughness at the estimated sigma.***"----
plot_marginal_effects(sim_opt)

## ----project, fig.width = 8, fig.height = 3, fig.cap = "***The covariates of the fitted model, in their own units.***"----
projected <- kernel_scale.raster(raster_stack = land_rast,
                                 multiScaleR  = opt,       # the fitted object
                                 verbose      = FALSE)
names(projected)

plot(projected[[1:3]], nc = 3, mar = c(2, 2, 2, 5))   # room for each legend

## ----standalone-multiscale, fig.width = 8, fig.height = 3, fig.cap = "***Root mean square roughness of land3 at three radii.***"----
sq_multi <- kernel_scale.raster(
  raster_stack = land_rast["land3"],
  scale_vars   = msr_vars(
    sq_200 = surface_var("land3", metric = "sq", radius = 200),
    sq_500 = surface_var("land3", metric = "sq", radius = 500),
    sq_900 = surface_var("land3", metric = "sq", radius = 900)
  ),
  verbose = FALSE
)

# Mean of each surface, in land3 units
round(t(global(sq_multi, "mean", na.rm = TRUE)), 3)

plot(sq_multi, nc = 3, mar = c(2, 2, 2, 5))   # room for each legend

## ----six-metrics, fig.width = 8, fig.height = 5.4, fig.cap = "***The six texture metrics of land3 within 500 m.***"----
six <- kernel_scale.raster(
  raster_stack = land_rast["land3"],
  scale_vars   = msr_vars(
    sa  = surface_var("land3", metric = "sa",  radius = 500),
    sq  = surface_var("land3", metric = "sq",  radius = 500),
    ssk = surface_var("land3", metric = "ssk", radius = 500),
    sku = surface_var("land3", metric = "sku", radius = 500),
    sdq = surface_var("land3", metric = "sdq", radius = 500),
    sdr = surface_var("land3", metric = "sdr", radius = 500)
  ),
  verbose = FALSE
)

plot(six, nc = 3, mar = c(2, 2, 2, 5))   # room for each legend

# Correlation among the metrics at the survey points
round(cor(extract(six, pts, ID = FALSE)), 2)

## ----weighted-vs-hard, fig.height = 3.4, fig.cap = "***Roughness of land3 within a 500 m hard radius, and with Gaussian weights at a sigma of 500 m.***"----
sq_hard <- kernel_scale.raster(
  raster_stack = land_rast["land3"],
  scale_vars   = msr_vars(sq_hard = surface_var("land3", metric = "sq",
                                                radius = 500)),
  verbose      = FALSE
)

sq_weighted <- kernel_scale.raster(
  raster_stack = land_rast["land3"],
  scale_vars   = msr_vars(sq_weighted = surface_var("land3", metric = "sq",
                                                    weighted = TRUE)),
  sigma        = 500,          # kernel width for the weighted metric (m)
  kernel       = "gaussian",
  verbose      = FALSE
)

plot(c(sq_hard, sq_weighted))

## ----quick-reference, eval = FALSE--------------------------------------------
# library(multiScaleR)
# library(terra)
# 
# # 1. Load the data
# data("landscape_counts"); data("surv_pts")
# dat       <- landscape_counts
# pts       <- vect(surv_pts)
# land_rast <- rast(system.file("extdata", "landscape.tif",
#                               package = "multiScaleR"))
# 
# # 2. Define the covariates: two kernel means and one fixed-radius metric
# scale_vars <- msr_vars(
#   land1_prop   = kernel_var("land1"),
#   land2_mean   = kernel_var("land2"),
#   land3_sq_500 = surface_var("land3", metric = "sq", radius = 500)
# )
# 
# # 3. Compute the covariates at the survey points
# kernel_inputs <- kernel_prep(pts = pts, raster_stack = land_rast,
#                              max_D = 1700, kernel = "gaussian",
#                              scale_vars = scale_vars, verbose = FALSE)
# kernel_inputs$scl_params            # mean and SD of each metric, in its units
# 
# # 4. Fit the model and estimate the kernel scales
# df  <- data.frame(dat, kernel_inputs$kernel_dat)
# mod <- glm(counts ~ site + land1_prop + land2_mean + land3_sq_500,
#            family = poisson(), data = df)
# opt <- multiScale_optim(fitted_mod = mod, kernel_inputs = kernel_inputs,
#                         verbose = FALSE)
# summary(opt)                        # sigma, effective distance, regression table
# 
# # 5. Read a metric coefficient in the metric's units
# beta <- coef(opt$opt_mod)
# sds  <- opt$scl_params$sd
# exp(beta["land3_sq_500"] * 0.01 / sds["land3_sq_500"])   # per 0.01 of roughness
# 
# # 6. Estimate the scale of a metric: weight it with the kernel
# weighted_vars <- msr_vars(
#   land3_sq = surface_var("land3", metric = "sq", weighted = TRUE)
# )
# truth_inputs <- kernel_prep(pts = pts, raster_stack = land_rast["land3"],
#                             max_D = 1500, kernel = "gaussian",
#                             scale_vars = weighted_vars, sigma = 300,
#                             verbose = FALSE)
# set.seed(1)
# sim_counts  <- rpois(100, exp(0.5 + 0.5 * truth_inputs$kernel_dat$land3_sq))
# free_inputs <- kernel_prep(pts = pts, raster_stack = land_rast["land3"],
#                            max_D = 1500, kernel = "gaussian",
#                            scale_vars = weighted_vars, verbose = FALSE)
# sim_df  <- data.frame(counts = sim_counts, free_inputs$kernel_dat)
# sim_mod <- glm(counts ~ land3_sq, family = poisson(), data = sim_df)
# sim_opt <- multiScale_optim(fitted_mod = sim_mod, kernel_inputs = free_inputs,
#                             verbose = FALSE)
# summary(sim_opt)$opt_scale                     # sigma, SE, and 95% interval
# kernel_dist(sim_opt)                           # radius holding 90% of the weight
# plot(profile_sigma(sim_opt, n_pts = 15, verbose = FALSE))   # AICc by sigma
# plot_marginal_effects(sim_opt)                 # expected count against roughness
# 
# # 7. Map the covariates of a fitted model, or a metric on its own
# projected <- kernel_scale.raster(raster_stack = land_rast,
#                                  multiScaleR = opt, verbose = FALSE)
# sq_multi  <- kernel_scale.raster(
#   raster_stack = land_rast["land3"],
#   scale_vars   = msr_vars(
#     sq_200 = surface_var("land3", metric = "sq", radius = 200),
#     sq_500 = surface_var("land3", metric = "sq", radius = 500),
#     sq_900 = surface_var("land3", metric = "sq", radius = 900)),
#   verbose = FALSE)
# sq_weighted <- kernel_scale.raster(
#   raster_stack = land_rast["land3"],
#   scale_vars   = msr_vars(sq_weighted = surface_var("land3", metric = "sq",
#                                                     weighted = TRUE)),
#   sigma = 500, kernel = "gaussian", verbose = FALSE)

