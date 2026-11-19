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

## ----fit----------------------------------------------------------------------
data("landscape_counts")   # counts and one site-level covariate
data("surv_pts")           # survey locations
pts <- vect(surv_pts)      # convert to a terra SpatVector

# Three landscape layers on one grid, in meters
land_rast <- rast(system.file("extdata", "landscape.tif",
                              package = "multiScaleR"))

# bin = TRUE stores the distance-binned summaries the fast backend needs
inputs <- kernel_prep(pts = pts, raster_stack = land_rast, max_D = 1700,
                      kernel = "gaussian", bin = TRUE, verbose = FALSE)

df    <- data.frame(landscape_counts, inputs$kernel_dat)
start <- glm(counts ~ site + land1 + land2, family = poisson(), data = df)

# Estimate the scales of effect and the coefficients
fit <- multiScale_optim(start, inputs, backend = "profile_glm",
                        n_cores = 1, verbose = FALSE)

round(fit$scale_est)       # sigma and its standard error, in meters

## ----project-unclamped--------------------------------------------------------
r_unclamped <- kernel_scale.raster(
  raster_stack = land_rast,
  multiScaleR  = fit,       # supplies the fitted kernels and the scaling
  scale_center = TRUE,      # standardize as the model data were standardized
  clamp        = FALSE,     # leave every value as computed
  verbose      = FALSE
)
r_unclamped

## ----predict-unclamped, fig.cap = "***Expected count from the unclamped projection, with survey points.***"----
# Expected count in every raster cell
pred_unclamped <- terra::predict(r_unclamped, fit$opt_mod, type = "response")

plot(pred_unclamped, main = "Expected count, unclamped")
points(pts, pch = 19, cex = 0.5, col = "red")

## ----count-range--------------------------------------------------------------
# Largest expected count on the map and at the survey points
c(map = minmax(pred_unclamped)[2], survey_points = max(fitted(fit$opt_mod)))

## ----compare-ranges-----------------------------------------------------------
covs <- c("land1", "land2")

# Range of each covariate at the survey points (row 1 = min, row 2 = max)
lim <- sapply(fit$opt_mod$model[covs], range)

# Range of each covariate across the whole raster
full <- minmax(r_unclamped[[covs]])

range_tab <- data.frame(
  covariate  = covs,
  survey_min = lim[1, ], survey_max = lim[2, ],
  raster_min = full[1, ], raster_max = full[2, ]
)
knitr::kable(range_tab, digits = 2, row.names = FALSE)

## ----map-outside, fig.cap = "***Cells where at least one covariate lies outside the range at the survey points (yellow), with survey points.***"----
# TRUE where a covariate is below the survey minimum or above the maximum
out_land1 <- r_unclamped$land1 < lim[1, "land1"] |
             r_unclamped$land1 > lim[2, "land1"]
out_land2 <- r_unclamped$land2 < lim[1, "land2"] |
             r_unclamped$land2 > lim[2, "land2"]
outside   <- out_land1 | out_land2

plot(outside, main = "Outside the surveyed range")
points(pts, pch = 19, cex = 0.5, col = "red")

# Percent of cells outside the range: either covariate, then each one
pct_out <- 100 * global(c(outside, out_land1, out_land2), "mean")[, 1]
names(pct_out) <- c("either", "land1", "land2")
round(pct_out, 1)

## ----project-clamped, fig.width = 7.5, fig.height = 3.6, fig.cap = "***Left, expected count from the clamped projection. Right, the clamped prediction divided by the unclamped prediction.***"----
r_clamped <- kernel_scale.raster(
  raster_stack = land_rast,
  multiScaleR  = fit,
  scale_center = TRUE,
  clamp        = TRUE,      # hold values at the surveyed minimum and maximum
  verbose      = FALSE
)
pred_clamped <- terra::predict(r_clamped, fit$opt_mod, type = "response")

# Ratio of the two predictions: 1 means clamping changed nothing
change <- pred_clamped / pred_unclamped

op <- par(mfrow = c(1, 2))
plot(pred_clamped, main = "Expected count, clamped")
plot(change, main = "Clamped / unclamped")
par(op)

# Smallest and largest ratio, and the new largest expected count
round(c(min_ratio = minmax(change)[1], max_ratio = minmax(change)[2],
        max_count = minmax(pred_clamped)[2]), 2)

## ----pct-mx-------------------------------------------------------------------
# Narrower limits: each one moves inward
r_narrow <- kernel_scale.raster(land_rast, multiScaleR = fit,
                                scale_center = TRUE, clamp = TRUE,
                                pct_mx = -0.1, verbose = FALSE)

# Wider limits: each one moves outward
r_wide   <- kernel_scale.raster(land_rast, multiScaleR = fit,
                                scale_center = TRUE, clamp = TRUE,
                                pct_mx = 0.1, verbose = FALSE)

pred_narrow <- terra::predict(r_narrow, fit$opt_mod, type = "response")
pred_wide   <- terra::predict(r_wide,   fit$opt_mod, type = "response")

# For each setting: the limits in force, the percent of cells clamped,
# and the largest expected count
layers <- list(r_narrow, r_clamped, r_wide)
preds  <- list(pred_narrow, pred_clamped, pred_wide)
pct_tab <- data.frame(
  pct_mx      = c(-0.1, 0, 0.1),
  land1_lower = sapply(layers, function(r) minmax(r$land1)[1]),
  land1_upper = sapply(layers, function(r) minmax(r$land1)[2]),
  land2_lower = sapply(layers, function(r) minmax(r$land2)[1]),
  land2_upper = sapply(layers, function(r) minmax(r$land2)[2]),
  pct_clamped = sapply(layers, function(r) {
    changed <- r$land1 != r_unclamped$land1 | r$land2 != r_unclamped$land2
    100 * global(changed, "mean")[1, 1]
  }),
  max_count   = sapply(preds, function(p) minmax(p)[2])
)
knitr::kable(pct_tab, digits = 2, row.names = FALSE)

## ----quick-reference, eval = FALSE--------------------------------------------
# library(multiScaleR)
# library(terra)
# 
# # 1. Fit the model
# data("landscape_counts"); data("surv_pts")
# pts       <- vect(surv_pts)
# land_rast <- rast(system.file("extdata", "landscape.tif",
#                               package = "multiScaleR"))
# inputs <- kernel_prep(pts = pts, raster_stack = land_rast, max_D = 1700,
#                       kernel = "gaussian", bin = TRUE)
# df     <- data.frame(landscape_counts, inputs$kernel_dat)
# start  <- glm(counts ~ site + land1 + land2, family = poisson(), data = df)
# fit    <- multiScale_optim(start, inputs, backend = "profile_glm", n_cores = 1)
# 
# # 2. Project without clamping and predict
# r_unclamped <- kernel_scale.raster(land_rast, multiScaleR = fit,
#                                    scale_center = TRUE, clamp = FALSE)
# pred_unclamped <- terra::predict(r_unclamped, fit$opt_mod, type = "response")
# 
# # 3. Map the cells outside the surveyed range
# covs <- c("land1", "land2")
# lim  <- sapply(fit$opt_mod$model[covs], range)     # survey min and max
# out_land1 <- r_unclamped$land1 < lim[1, "land1"] |
#              r_unclamped$land1 > lim[2, "land1"]
# out_land2 <- r_unclamped$land2 < lim[1, "land2"] |
#              r_unclamped$land2 > lim[2, "land2"]
# outside   <- out_land1 | out_land2
# plot(outside)                                      # where the model extrapolates
# 100 * global(outside, "mean")[1, 1]                # percent of the map
# 
# # 4. Clamp to the surveyed range and predict again
# r_clamped <- kernel_scale.raster(land_rast, multiScaleR = fit,
#                                  scale_center = TRUE, clamp = TRUE)
# pred_clamped <- terra::predict(r_clamped, fit$opt_mod, type = "response")
# plot(pred_clamped / pred_unclamped)                # 1 = unchanged by clamping
# 
# # 5. Narrow or widen the limits
# r_narrow <- kernel_scale.raster(land_rast, multiScaleR = fit,
#                                 scale_center = TRUE, clamp = TRUE,
#                                 pct_mx = -0.1)     # limits move inward
# r_wide   <- kernel_scale.raster(land_rast, multiScaleR = fit,
#                                 scale_center = TRUE, clamp = TRUE,
#                                 pct_mx = 0.1)      # limits move outward

