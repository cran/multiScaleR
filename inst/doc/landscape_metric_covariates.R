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

## ----class-level-example------------------------------------------------------
# ## Clumpiness of land-cover class 2, and that class's percentage cover,
# ## both within a fixed 500 m radius
# class_vars <- msr_vars(
#   cover2_clumpy = landscape_var("landcover", metric = "clumpy", radius = 500, class = 2),
#   cover2_pland  = landscape_var("landcover", metric = "pland",  radius = 500, class = 2)
# )
# 
# class_vars

## ----load-data----------------------------------------------------------------
# data("landscape_counts")
# dat <- landscape_counts
# 
# data("surv_pts")
# pts <- vect(surv_pts)
# 
# land_rast <- rast(pkg_extdata("landscape.tif"))

## ----build-landcover----------------------------------------------------------
# landcover <- classify(
#   land_rast$land2,
#   rcl = matrix(
#     c(-Inf, -0.5, 1,
#       -0.5,  0.5, 2,
#        0.5,  Inf, 3),
#     ncol = 3,
#     byrow = TRUE
#   )
# )
# names(landcover) <- "landcover"
# 
# metric_rasters <- c(land_rast$land1, landcover)

## ----plot-inputs, fig.cap = "***Source rasters used to derive model covariates.***"----
# plot(metric_rasters)
# plot(metric_rasters$land1, main = "land1 with survey points")
# points(pts, pch = 19, cex = 0.6, col = "red")

## ----define-vars--------------------------------------------------------------
# scale_vars <- msr_vars(
#   land1_prop = kernel_var("land1"),
#   land1_ed_500 = landscape_var("land1", metric = "ed", radius = 500),
#   landcover_shdi_500 = landscape_var("landcover", metric = "shdi", radius = 500)
# )
# 
# scale_vars

## ----kernel-prep--------------------------------------------------------------
# kernel_inputs <- kernel_prep(
#   pts = pts,
#   raster_stack = metric_rasters,
#   max_D = 1200,
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
#   counts ~ site + land1_prop + land1_ed_500 + landcover_shdi_500,
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
opt <- vignette_cache$landscape$opt
summary(opt)

## ----project, fig.cap = "***Projected derived covariates.***", eval = TRUE----
projected <- terra::unwrap(vignette_cache$landscape$projected)

plot(projected)

## ----projected-names, eval = TRUE---------------------------------------------
names(projected)

## ----standalone-metric-rasters, fig.cap = "***Landscape metric surfaces computed directly from the source rasters at a fixed 500 m radius.***", fig.height = 7, eval = FALSE----
# ## Define several metrics at a fixed 500 m radius
# standalone_vars <- msr_vars(
#   land1_ai    = landscape_var("land1",    metric = "ai",    radius = 500),
#   land1_ed    = landscape_var("land1",    metric = "ed",    radius = 500),
#   land1_lsi   = landscape_var("land1",    metric = "lsi",   radius = 500),
#   land1_pladj = landscape_var("land1",    metric = "pladj", radius = 500),
#   cover_shdi  = landscape_var("landcover", metric = "shdi",  radius = 500),
#   cover_pr    = landscape_var("landcover", metric = "pr",    radius = 500),
#   cover_sidi  = landscape_var("landcover", metric = "sidi",  radius = 500)
# )
# 
# ## Apply directly: no multiScaleR object or sigma values required
# metric_surfaces <- kernel_scale.raster(
#   raster_stack = metric_rasters,
#   scale_vars   = standalone_vars,
#   verbose      = FALSE
# )
# 
# plot(metric_surfaces)

## ----standalone-names, eval = FALSE-------------------------------------------
# names(metric_surfaces)

## ----standalone-single, fig.cap = "***Edge density within a 300 m radius.***", eval = FALSE----
# ed_300 <- kernel_scale.raster(
#   raster_stack = metric_rasters["land1"],
#   scale_vars   = msr_vars(ed_300m = landscape_var("land1", metric = "ed", radius = 300)),
#   verbose      = FALSE
# )
# 
# plot(ed_300, main = "Edge density (m/ha) within 300 m")

## ----standalone-multiscale, fig.cap = "***Edge density at three spatial scales.***", fig.height = 3.5, eval = FALSE----
# ed_multi <- kernel_scale.raster(
#   raster_stack = metric_rasters["land1"],
#   scale_vars   = msr_vars(
#     ed_200m = landscape_var("land1", metric = "ed", radius = 200),
#     ed_500m = landscape_var("land1", metric = "ed", radius = 500),
#     ed_900m = landscape_var("land1", metric = "ed", radius = 900)
#   ),
#   verbose = FALSE
# )
# 
# plot(ed_multi,
#      main = c("ED: 200 m", "ED: 500 m", "ED: 900 m"))

## ----info-theory-surfaces, fig.cap = "***Composition and configuration surfaces for the same landscape.***", fig.height = 3.5, eval = FALSE----
# info_surfaces <- kernel_scale.raster(
#   raster_stack = metric_rasters["landcover"],
#   scale_vars   = msr_vars(
#     H_x    = landscape_var("landcover", metric = "ent",    radius = 500, base = 2),
#     MI     = landscape_var("landcover", metric = "mutinf", radius = 500, base = 2),
#     relMI  = landscape_var("landcover", metric = "relmutinf", radius = 500)
#   ),
#   verbose = FALSE
# )
# 
# plot(info_surfaces, main = c("Marginal entropy (bits)",
#                              "Mutual information (bits)",
#                              "Relative mutual information"))

## ----clumpy-surface, fig.cap = "***Clumpiness and cover of land-cover class 2.***", fig.height = 3.5, eval = FALSE----
# class2_surfaces <- kernel_scale.raster(
#   raster_stack = metric_rasters["landcover"],
#   scale_vars   = msr_vars(
#     clumpy2 = landscape_var("landcover", metric = "clumpy", radius = 500, class = 2),
#     pland2  = landscape_var("landcover", metric = "pland",  radius = 500, class = 2)
#   ),
#   verbose = FALSE
# )
# 
# plot(class2_surfaces, main = c("CLUMPY (class 2)", "PLAND (class 2, %)"))

## ----optimized-landscape-vars, eval = FALSE-----------------------------------
# scale_vars_opt <- msr_vars(
#   land1_prop = kernel_var("land1"),
#   land1_ed = landscape_var("land1", metric = "ed"),
#   landcover_shdi = landscape_var("landcover", metric = "shdi")
# )

