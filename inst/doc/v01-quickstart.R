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

# Slow steps (optimization, likelihood profiles, raster projection) are stored
# in the package so the vignette builds quickly. The text shows the real calls.
# Regenerate them with the `quickstart` section of
# data-raw/generate_vignette_cache.R.
cache_file <- file.path("..", "inst", "extdata", "cache-quickstart.rds")
if (!file.exists(cache_file)) {
  cache_file <- system.file("extdata", "cache-quickstart.rds",
                            package = "multiScaleR")
}
qs_cache <- readRDS(cache_file)

## ----install, eval = FALSE----------------------------------------------------
# install.packages("multiScaleR")

## ----load-data----------------------------------------------------------------
data("landscape_counts")   # data frame: counts and one site-level covariate
dat <- landscape_counts

data("surv_pts")           # survey locations (an sf object from the package)
pts <- vect(surv_pts)      # convert to a terra SpatVector

# Three landscape layers on one grid, in meters
land_rast <- rast(system.file("extdata", "landscape.tif",
                              package = "multiScaleR"))

## ----plot-data, fig.cap = "***The three landscape layers, and the survey points over `land1`.***"----
plot(land_rast)
plot(land_rast$land1, main = "land1 with survey points")
points(pts, pch = 19, cex = 0.6, col = "red")

## ----kernel-prep--------------------------------------------------------------
kernel_inputs <- kernel_prep(
  pts          = pts,          # survey locations
  raster_stack = land_rast,    # landscape layers to summarize
  max_D        = 1700,         # largest distance (m) considered around a point
  kernel       = "gaussian",   # shape of the distance-decay function
  verbose      = FALSE
)

## ----build-df-----------------------------------------------------------------
df <- data.frame(dat, kernel_inputs$kernel_dat)
head(df, 3)

## ----fit-mod------------------------------------------------------------------
mod <- glm(counts ~ site + land1 + land2,
           family = poisson(),
           data   = df)

## ----optim, eval = FALSE------------------------------------------------------
# opt <- multiScale_optim(fitted_mod    = mod,
#                         kernel_inputs = kernel_inputs)

## ----load-cached-opt, include = FALSE-----------------------------------------
opt <- qs_cache$opt

## ----summary------------------------------------------------------------------
summary(opt)

## ----read-summary-------------------------------------------------------------
est     <- summary(opt)$opt_scale          # sigma, SE, and 95% interval (m)
dist90  <- kernel_dist(opt)                # radius holding 90% of kernel weight (m)
effects <- exp(coef(opt$opt_mod))          # multiplicative change in expected count

round(dist90)
round(effects, 2)

## ----truth-table--------------------------------------------------------------
ci_beta <- confint.default(opt$opt_mod)    # Wald 95% intervals, log scale
check <- data.frame(
  Parameter = c("sigma, land1 (m)", "sigma, land2 (m)",
                "effect of land1", "effect of land2", "effect of site"),
  Truth     = c(250, 500, -0.5, 0.7, 0.3),
  Estimate  = c(est[, "Mean"], coef(opt$opt_mod)[c("land1", "land2", "site")]),
  Lower     = c(est[, "2.5%"], ci_beta[c("land1", "land2", "site"), 1]),
  Upper     = c(est[, "97.5%"], ci_beta[c("land1", "land2", "site"), 2])
)
# Does each 95% interval contain the value used to simulate the data?
check$Covered <- check$Truth >= check$Lower & check$Truth <= check$Upper
knitr::kable(check, digits = 2, row.names = FALSE)

## ----plot-kernel, fig.cap = "***Kernel weight against distance for each landscape covariate.***"----
plot(opt)

## ----plot-kernel-plain, fig.cap = "***The same kernel plots without the text under each panel.***"----
plot(opt, add_label = FALSE)

## ----marginal, fig.cap = "***Estimated effect of each covariate on the expected count.***"----
plot_marginal_effects(opt)

## ----profile-call, eval = FALSE-----------------------------------------------
# prof <- profile_sigma(opt, n_pts = 15)   # 15 sigma values per covariate
# plot(prof)

## ----profile, echo = FALSE, fig.cap = "***AICc across candidate sigma values. The dashed line marks the estimate.***"----
prof <- qs_cache$prof
plot(prof)

## ----project-call, eval = FALSE-----------------------------------------------
# r_scaled <- kernel_scale.raster(raster_stack = land_rast,
#                                 multiScaleR  = opt,     # the fitted object
#                                 scale_center = TRUE,    # same scaling as the model
#                                 clamp        = TRUE)    # stay within the fitted range
# 
# # Expected count in every raster cell
# pred <- terra::predict(r_scaled, opt$opt_mod, type = "response")

## ----project, echo = FALSE, fig.cap = "***Landscape layers smoothed at their estimated scales, in standard deviation units.***"----
r_scaled <- terra::unwrap(qs_cache$r_scaled)
plot(r_scaled)

## ----predict, echo = FALSE, fig.cap = "***Expected count across the landscape, with survey points.***"----
pred <- terra::unwrap(qs_cache$pred)
plot(pred, main = "Expected count")
points(pts, pch = 19, cex = 0.5, col = "red")

## ----full-workflow, eval = FALSE----------------------------------------------
# library(multiScaleR)
# library(terra)
# 
# # 1. Load the data
# data("landscape_counts"); data("surv_pts")
# pts       <- vect(surv_pts)
# land_rast <- rast(system.file("extdata", "landscape.tif",
#                               package = "multiScaleR"))
# 
# # 2. Prepare the kernel inputs
# kernel_inputs <- kernel_prep(pts = pts, raster_stack = land_rast,
#                              max_D = 1700, kernel = "gaussian")
# 
# # 3. Fit a starting model
# df  <- data.frame(landscape_counts, kernel_inputs$kernel_dat)
# mod <- glm(counts ~ site + land1 + land2, family = poisson(), data = df)
# 
# # 4. Estimate the scales of effect
# opt <- multiScale_optim(fitted_mod = mod, kernel_inputs = kernel_inputs)
# 
# # 5. Read the results
# summary(opt)                   # sigma, effective distance, regression table
# summary(opt, profile = TRUE)   # profile-likelihood intervals for sigma (slower)
# kernel_dist(opt)               # radius holding 90% of the kernel weight
# exp(coef(opt$opt_mod))         # multiplicative change in expected count
# plot(opt)                      # kernel weight against distance
# plot(opt, add_label = FALSE)   # the same plot without the text beneath it
# plot_marginal_effects(opt)     # expected count across each covariate
# plot(profile_sigma(opt))       # AICc across sigma: is each scale identified?
# 
# # 6. Map the fitted model
# r_scaled <- kernel_scale.raster(land_rast, multiScaleR = opt,
#                                 scale_center = TRUE, clamp = TRUE)
# pred <- terra::predict(r_scaled, opt$opt_mod, type = "response")

