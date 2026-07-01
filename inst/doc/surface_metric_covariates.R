## ----setup, include = FALSE---------------------------------------------------
library(multiScaleR)
library(terra)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6.5,
  fig.height = 4.25,
  eval = FALSE,
  warning = FALSE,
  message = FALSE
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

vignette_cache <- readRDS(pkg_extdata("vignette_cache.rds"))

## ----load-data----------------------------------------------------------------
# data("landscape_counts")
# dat <- landscape_counts
# 
# data("surv_pts")
# pts <- vect(surv_pts)
# 
# land_rast <- rast(pkg_extdata("landscape.tif"))

## ----plot-inputs, fig.cap = "***Continuous source surfaces. land3 (right) is smooth; land2 (left) is noisy.***"----
# surfaces <- c(land_rast$land2, land_rast$land3)
# plot(surfaces)
# plot(land_rast$land3, main = "land3 with survey points")
# points(pts, pch = 19, cex = 0.6, col = "red")

## ----define-vars--------------------------------------------------------------
# scale_vars <- msr_vars(
#   land2_mean   = kernel_var("land2"),
#   land3_sq_500 = surface_var("land3", metric = "sq", radius = 500),
#   land3_sa_500 = surface_var("land3", metric = "sa", radius = 500)
# )
# 
# scale_vars

## ----kernel-prep--------------------------------------------------------------
# kernel_inputs <- kernel_prep(
#   pts = pts,
#   raster_stack = land_rast,
#   max_D = 1500,
#   kernel = "gaussian",
#   scale_vars = scale_vars,
#   verbose = FALSE
# )
# 
# head(kernel_inputs$kernel_dat)

## ----fit-model----------------------------------------------------------------
# df <- data.frame(dat, kernel_inputs$kernel_dat)
# 
# mod <- glm(
#   counts ~ site + land2_mean + land3_sq_500 + land3_sa_500,
#   family = poisson(),
#   data = df
# )

## ----optimize-----------------------------------------------------------------
# opt <- multiScale_optim(
#   fitted_mod = mod,
#   kernel_inputs = kernel_inputs,
#   verbose = FALSE
# )

## ----summary, eval = TRUE-----------------------------------------------------
opt <- vignette_cache$surface$opt
summary(opt)

## ----coef-interpret, eval = TRUE----------------------------------------------
cf <- coef(summary(opt$opt_mod))
## Multiplicative change in expected count per 1 SD increase in each covariate
round(exp(cf[, "Estimate"]), 3)

## ----project, fig.cap = "***Projected model covariates at their fitted or fixed scales.***", fig.height = 6, eval = TRUE----
projected <- terra::unwrap(vignette_cache$surface$projected)

plot(projected)
names(projected)

## ----standalone-multiscale, fig.cap = "***RMS roughness (sq) of land3 at three radii. Larger radii fold in broader-scale variation.***", fig.height = 3.5----
# sq_multi <- kernel_scale.raster(
#   raster_stack = land_rast["land3"],
#   scale_vars = msr_vars(
#     sq_200 = surface_var("land3", metric = "sq", radius = 200),
#     sq_500 = surface_var("land3", metric = "sq", radius = 500),
#     sq_900 = surface_var("land3", metric = "sq", radius = 900)
#   ),
#   verbose = FALSE
# )
# 
# plot(sq_multi, main = c("sq - 200 m", "sq - 500 m", "sq - 900 m"))

## ----sa-vs-sq, fig.cap = "***Average roughness (sa) versus RMS roughness (sq) of land3 at 500 m.***", fig.height = 3.5----
# sa_sq <- kernel_scale.raster(
#   raster_stack = land_rast["land3"],
#   scale_vars = msr_vars(
#     sa_500 = surface_var("land3", metric = "sa", radius = 500),
#     sq_500 = surface_var("land3", metric = "sq", radius = 500)
#   ),
#   verbose = FALSE
# )
# 
# plot(sa_sq, main = c("sa - 500 m", "sq - 500 m"))

## ----shape-slope, fig.cap = "***Skewness, kurtosis, RMS slope, and surface area ratio of land3 at 500 m.***", fig.height = 6----
# shape_geom <- kernel_scale.raster(
#   raster_stack = land_rast["land3"],
#   scale_vars = msr_vars(
#     ssk_500 = surface_var("land3", metric = "ssk", radius = 500),
#     sku_500 = surface_var("land3", metric = "sku", radius = 500),
#     sdq_500 = surface_var("land3", metric = "sdq", radius = 500),
#     sdr_500 = surface_var("land3", metric = "sdr", radius = 500)
#   ),
#   verbose = FALSE
# )
# 
# plot(shape_geom,
#      main = c("ssk - 500 m", "sku - 500 m", "sdq - 500 m", "sdr - 500 m"))

## ----optimized-surface-vars, eval = FALSE-------------------------------------
# scale_vars_opt <- msr_vars(
#   land2_mean = kernel_var("land2"),
#   land3_sq   = surface_var("land3", metric = "sq")
# )
# 
# kernel_inputs_opt <- kernel_prep(
#   pts = pts,
#   raster_stack = land_rast,
#   max_D = 1500,
#   scale_vars = scale_vars_opt,
#   verbose = FALSE
# )

## ----weighted-define----------------------------------------------------------
# weighted_vars <- msr_vars(
#   land3_sq_weighted = surface_var("land3", metric = "sq", weighted = TRUE)
# )

## ----weighted-vs-hard, fig.cap = "***Hard-radius versus kernel-weighted RMS roughness of land3 at a 500 m scale.***", fig.height = 3.5----
# sq_hard <- kernel_scale.raster(
#   raster_stack = land_rast["land3"],
#   scale_vars = msr_vars(sq_hard = surface_var("land3", metric = "sq",
#                                               radius = 500)),
#   verbose = FALSE
# )
# 
# sq_weighted <- kernel_scale.raster(
#   raster_stack = land_rast["land3"],
#   scale_vars = weighted_vars,
#   sigma = 500,            # the kernel scale for the weighted metric
#   kernel = "gaussian",
#   verbose = FALSE
# )
# 
# plot(c(sq_hard, sq_weighted), main = c("sq hard radius 500 m", "sq weighted sigma 500"))

## ----quick-reference, eval = FALSE--------------------------------------------
# ## 1. Define covariates: an optimized kernel mean plus fixed-radius roughness
# scale_vars <- msr_vars(
#   land2_mean   = kernel_var("land2"),
#   land3_sq_500 = surface_var("land3", metric = "sq", radius = 500),
#   land3_sa_500 = surface_var("land3", metric = "sa", radius = 500)
# )
# 
# ## 2. Extract and prepare the covariates around each point
# kernel_inputs <- kernel_prep(pts, land_rast, max_D = 1500,
#                              scale_vars = scale_vars, verbose = FALSE)
# 
# ## 3. Fit the starting model
# df  <- data.frame(dat, kernel_inputs$kernel_dat)
# mod <- glm(counts ~ site + land2_mean + land3_sq_500 + land3_sa_500,
#            family = poisson(), data = df)
# 
# ## 4. Optimize the free scale parameter(s)
# opt <- multiScale_optim(fitted_mod = mod, kernel_inputs = kernel_inputs)
# summary(opt)                              # sigma + effective distance, on map units
# 
# ## 5. Project covariates for prediction (centered and scaled like the model)
# projected <- kernel_scale.raster(land_rast, multiScaleR = opt,
#                                  scale_center = TRUE, clamp = TRUE)
# 
# ## Standalone roughness surface, no model required
# sq_500 <- kernel_scale.raster(
#   land_rast["land3"],
#   scale_vars = msr_vars(sq_500 = surface_var("land3", metric = "sq", radius = 500))
# )

