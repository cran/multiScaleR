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

pkg_extdata <- function(...) {
  src_path <- normalizePath(file.path("..", "inst", "extdata", ...),
                            winslash = "/",
                            mustWork = FALSE)
  if (file.exists(src_path)) {
    return(src_path)
  }

  inst_path <- system.file("extdata", ..., package = "multiScaleR")
  if (nzchar(inst_path) && file.exists(inst_path)) {
    return(inst_path)
  }

  stop("Could not locate required extdata file: ", paste(..., collapse = "/"))
}

## ----install, eval = FALSE----------------------------------------------------
# install.packages("multiScaleR")

## ----load-data----------------------------------------------------------------
data("landscape_counts")
dat <- landscape_counts

data("surv_pts")
pts <- vect(surv_pts)

land_rast <- rast(pkg_extdata("landscape.tif"))

## ----plot-data, fig.cap = "***Landscape rasters with survey locations.***"----
plot(land_rast)
plot(land_rast$land1, main = "land1 with survey points")
points(pts, pch = 19, cex = 0.6, col = 'red')

## ----kernel-prep--------------------------------------------------------------
kernel_inputs <- kernel_prep(
  pts          = pts,
  raster_stack = land_rast,
  max_D        = 1700,
  kernel       = "gaussian",
  verbose      = FALSE
)

## ----build-df-----------------------------------------------------------------
df <- data.frame(dat, kernel_inputs$kernel_dat)

## ----fit-mod------------------------------------------------------------------
mod <- glm(counts ~ site + land1 + land2,
           family = poisson(),
           data   = df)

## ----optim, eval = TRUE-------------------------------------------------------
opt <- multiScaleR::multiScale_optim(fitted_mod    = mod,
                                     kernel_inputs = kernel_inputs)

## ----load-opt, include = FALSE------------------------------------------------
# load(pkg_extdata("opt2.RData"))
# opt <- opt2

## ----summary------------------------------------------------------------------
summary(opt)

## ----plot-kernel, fig.cap = "***Distance-decay of kernel weight for each covariate.***"----
plot(opt)

## ----marginal, fig.cap = "***Marginal effects of each spatial covariate on counts.***"----
plot_marginal_effects(opt)

## ----profile, fig.cap = "***AICc profile across sigma. Red dashed line marks optimized sigma.***"----
prof <- profile_sigma(opt, n_pts = 15, verbose = FALSE)
plot(prof)

## ----project, fig.cap = "***Kernel-smoothed raster layers at optimized scales.***"----
r_scaled <- kernel_scale.raster(raster_stack = land_rast,
                                multiScaleR  = opt,
                                scale_center = TRUE,
                                clamp        = TRUE
                                )
plot(r_scaled)

## ----predict, fig.cap = "***Predicted counts across the landscape.***"--------
pred <- terra::predict(r_scaled, opt$opt_mod, type = "response")
plot(pred, main = "Predicted counts")
plot(surv_pts, add = T, col = 'red')

## ----full-workflow, eval = FALSE----------------------------------------------
# library(multiScaleR)
# library(terra)
# 
# # 1. Load data
# data("landscape_counts"); data("surv_pts")
# pts      <- vect(surv_pts)
# land_rast <- rast(pkg_extdata("landscape.tif"))
# 
# # 2. Prepare kernel inputs
# kernel_inputs <- kernel_prep(pts = pts, raster_stack = land_rast,
#                              max_D = 1700, kernel = "gaussian")
# 
# # 3. Fit initial model
# df  <- data.frame(landscape_counts, kernel_inputs$kernel_dat)
# mod <- glm(counts ~ land1 + land2, family = poisson(), data = df)
# 
# # 4. Optimize scales of effect
# opt <- multiScale_optim(fitted_mod = mod, kernel_inputs = kernel_inputs,
#                         n_cores = 2)
# 
# # 5. Summarize and visualize
# summary(opt)                        # sigma estimates + regression table
# summary(opt, profile = TRUE)        # profile-likelihood CIs (slower)
# plot(opt)                           # kernel decay curves
# plot_marginal_effects(opt)          # covariate effect plots
# prof <- profile_sigma(opt)          # AICc/LL surface across sigma space
# plot(prof)                          # plot profile with optimized sigma marked
# 
# # 6. Project to landscape
# r_scaled <- kernel_scale.raster(land_rast, multiScaleR = opt,
#                                 scale_center = TRUE, clamp = TRUE)
# terra::predict(r_scaled, opt$opt_mod, type = "response")

