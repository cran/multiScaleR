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

## ----build-landcover----------------------------------------------------------
# Class 1 below 0.45, class 2 from 0.45 to 0.55, class 3 above 0.55
landcover <- classify(
  land_rast$land3,
  rcl = matrix(c(-Inf, 0.45, 1,
                 0.45, 0.55, 2,
                 0.55,  Inf, 3),
               ncol = 3, byrow = TRUE)
)
names(landcover) <- "landcover"

# The source layers for this example, on one grid
metric_rasters <- c(land_rast$land1, land_rast$land2, landcover)

## ----plot-inputs, fig.height = 3.4, fig.cap = "***The binary habitat layer and the three-class land-cover layer.***"----
plot(metric_rasters[[c("land1", "landcover")]])

## ----define-vars--------------------------------------------------------------
scale_vars <- msr_vars(
  land1_prop     = kernel_var("land1"),     # proportion of habitat, sigma estimated
  land2_mean     = kernel_var("land2"),     # mean of land2, sigma estimated
  land1_ed_500   = landscape_var("land1", metric = "ed",
                                 radius = 500),   # edge density within 500 m
  cover_shdi_500 = landscape_var("landcover", metric = "shdi",
                                 radius = 500)    # Shannon diversity within 500 m
)

scale_vars[, c("covariate", "source", "type", "metric", "radius", "optimize")]

## ----kernel-prep--------------------------------------------------------------
kernel_inputs <- kernel_prep(
  pts          = pts,              # survey locations
  raster_stack = metric_rasters,   # source layers named in scale_vars
  max_D        = 1700,             # largest distance (m) considered
  kernel       = "gaussian",       # kernel for the two kernel covariates
  scale_vars   = scale_vars,
  verbose      = FALSE
)

# Mean and standard deviation of each covariate, in its own units
round(sapply(kernel_inputs$scl_params, function(x) x), 3)

## ----fit-model----------------------------------------------------------------
df <- data.frame(dat, kernel_inputs$kernel_dat)

# Starting model: formula and family carry through to the optimizer
mod <- glm(counts ~ site + land1_prop + land2_mean +
             land1_ed_500 + cover_shdi_500,
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

# A step that is easy to picture, in the units of each metric
step <- c(land1_ed_500   = 10,         # 10 m more edge per hectare
          cover_shdi_500 = 0.1)        # 0.1 more Shannon diversity
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
# Edge density within 400 m, the scale used to simulate the counts
truth_inputs <- kernel_prep(
  pts          = pts,
  raster_stack = land_rast["land1"],
  max_D        = 1200,
  scale_vars   = msr_vars(ed_400 = landscape_var("land1", metric = "ed",
                                                 radius = 400)),
  verbose      = FALSE
)

set.seed(1)
# kernel_dat is standardized, so 0.5 is the effect of one standard deviation
sim_counts <- rpois(100, exp(0.5 + 0.5 * truth_inputs$kernel_dat$ed_400))

## ----sim-fit------------------------------------------------------------------
# No radius: the buffer radius becomes a parameter of the model
free_vars <- msr_vars(land1_ed = landscape_var("land1", metric = "ed"))

free_inputs <- kernel_prep(
  pts          = pts,
  raster_stack = land_rast["land1"],
  max_D        = 1200,             # the radius is searched up to this distance
  scale_vars   = free_vars,
  verbose      = FALSE
)

sim_df  <- data.frame(counts = sim_counts, free_inputs$kernel_dat)
sim_mod <- glm(counts ~ land1_ed, family = poisson(), data = sim_df)

sim_opt <- multiScale_optim(fitted_mod    = sim_mod,
                            kernel_inputs = free_inputs,
                            verbose       = FALSE)

# Radius with profile-likelihood limits, in meters
sim_scale <- summary(sim_opt, profile = TRUE)$opt_scale
round(sim_scale, 1)

## ----sim-profile, fig.height = 3.6, fig.cap = "***AICc across candidate radii for edge density. The dashed line marks the estimate.***"----
# Refit the model at 15 candidate radii
plot(profile_sigma(sim_opt, n_pts = 15, verbose = FALSE))

## ----sim-coef-----------------------------------------------------------------
sim_beta <- coef(sim_opt$opt_mod)["land1_ed"]             # per standard deviation
sim_ci   <- confint.default(sim_opt$opt_mod)["land1_ed", ]
sim_sd   <- sim_opt$scl_params$sd["land1_ed"]             # m/ha at the fitted radius

# Change in expected count for 10 m more edge per hectare
sim_effect <- exp(c(Estimate = unname(sim_beta), Lower = sim_ci[[1]],
                    Upper = sim_ci[[2]]) * 10 / sim_sd)
round(sim_effect, 3)

# The value used in the simulation, on the same footing
round(exp(0.5 * 10 / truth_inputs$scl_params$sd["ed_400"]), 3)

## ----sim-marginal, fig.height = 3.6, fig.cap = "***Expected count across the observed range of edge density at the estimated radius.***"----
plot_marginal_effects(sim_opt)

## ----radius-grid--------------------------------------------------------------
# Refit the model at nine radii from 200 to 1,000 m
grid_fit <- multiScale_optim(sim_mod, free_inputs, backend = "refit",
                             radius_grid = seq(200, 1000, by = 100),
                             verbose = FALSE)
radii <- diagnostics(grid_fit)$radius_profile
radii[, c("radius", "log_likelihood", "status", "selected")]

## ----project, fig.height = 5.5, fig.cap = "***The covariates of the fitted model, in their own units.***"----
projected <- kernel_scale.raster(raster_stack = metric_rasters,
                                 multiScaleR  = opt,       # the fitted object
                                 verbose      = FALSE)
names(projected)

plot(projected[[1:4]])

## ----standalone-multiscale, fig.width = 8, fig.height = 3, fig.cap = "***Edge density of land1 (m/ha) at three radii.***"----
ed_multi <- kernel_scale.raster(
  raster_stack = metric_rasters["land1"],
  scale_vars   = msr_vars(
    ed_200 = landscape_var("land1", metric = "ed", radius = 200),
    ed_500 = landscape_var("land1", metric = "ed", radius = 500),
    ed_900 = landscape_var("land1", metric = "ed", radius = 900)
  ),
  verbose = FALSE
)

# Range of each surface, in m of edge per hectare
round(t(global(ed_multi, c("min", "max"), na.rm = TRUE)))

plot(ed_multi, nc = 3, mar = c(2, 2, 2, 5))   # room for each legend

## ----info-theory-surfaces, fig.width = 8, fig.height = 3, fig.cap = "***Marginal entropy, mutual information, and relative mutual information of the land-cover layer within 500 m.***"----
info_surfaces <- kernel_scale.raster(
  raster_stack = metric_rasters["landcover"],
  scale_vars   = msr_vars(
    H_x   = landscape_var("landcover", metric = "ent",
                          radius = 500, base = 2),      # bits
    MI    = landscape_var("landcover", metric = "mutinf",
                          radius = 500, base = 2),      # bits
    relMI = landscape_var("landcover", metric = "relmutinf",
                          radius = 500)                 # 0 to 1
  ),
  verbose = FALSE
)

plot(info_surfaces, nc = 3, mar = c(2, 2, 2, 5))   # room for each legend

## ----clumpy-surface, fig.height = 3.2, fig.cap = "***Percentage cover and clumpiness of land-cover class 2 within 500 m.***"----
class2_surfaces <- kernel_scale.raster(
  raster_stack = metric_rasters["landcover"],
  scale_vars   = msr_vars(
    pland2  = landscape_var("landcover", metric = "pland",
                            radius = 500, class = 2),   # percent of the buffer
    clumpy2 = landscape_var("landcover", metric = "clumpy",
                            radius = 500, class = 2)    # -1 to 1
  ),
  verbose = FALSE
)

plot(class2_surfaces)

## ----quick-reference, eval = FALSE--------------------------------------------
# library(multiScaleR)
# library(terra)
# 
# # 1. Load the data and build a categorical layer
# data("landscape_counts"); data("surv_pts")
# dat       <- landscape_counts
# pts       <- vect(surv_pts)
# land_rast <- rast(system.file("extdata", "landscape.tif",
#                               package = "multiScaleR"))
# landcover <- classify(land_rast$land3,
#                       rcl = matrix(c(-Inf, 0.45, 1,
#                                      0.45, 0.55, 2,
#                                      0.55,  Inf, 3), ncol = 3, byrow = TRUE))
# names(landcover) <- "landcover"
# metric_rasters <- c(land_rast$land1, land_rast$land2, landcover)
# 
# # 2. Define the covariates: two kernel means and two fixed-radius metrics
# scale_vars <- msr_vars(
#   land1_prop     = kernel_var("land1"),
#   land2_mean     = kernel_var("land2"),
#   land1_ed_500   = landscape_var("land1", metric = "ed", radius = 500),
#   cover_shdi_500 = landscape_var("landcover", metric = "shdi", radius = 500)
# )
# 
# # 3. Compute the covariates at the survey points
# kernel_inputs <- kernel_prep(pts = pts, raster_stack = metric_rasters,
#                              max_D = 1700, kernel = "gaussian",
#                              scale_vars = scale_vars, verbose = FALSE)
# kernel_inputs$scl_params            # mean and SD of each metric, in its units
# 
# # 4. Fit the model and estimate the kernel scales
# df  <- data.frame(dat, kernel_inputs$kernel_dat)
# mod <- glm(counts ~ site + land1_prop + land2_mean +
#              land1_ed_500 + cover_shdi_500,
#            family = poisson(), data = df)
# opt <- multiScale_optim(fitted_mod = mod, kernel_inputs = kernel_inputs,
#                         verbose = FALSE)
# summary(opt)                        # sigma, effective distance, regression table
# 
# # 5. Read a metric coefficient in the metric's units
# beta <- coef(opt$opt_mod)
# sds  <- opt$scl_params$sd
# exp(beta["land1_ed_500"] * 10 / sds["land1_ed_500"])   # per 10 m/ha of edge
# 
# # 6. Estimate the radius of a metric: leave radius out
# set.seed(1)
# truth_inputs <- kernel_prep(pts = pts, raster_stack = land_rast["land1"],
#                             max_D = 1200, verbose = FALSE,
#                             scale_vars = msr_vars(
#                               ed_400 = landscape_var("land1", metric = "ed",
#                                                      radius = 400)))
# sim_counts  <- rpois(100, exp(0.5 + 0.5 * truth_inputs$kernel_dat$ed_400))
# free_vars   <- msr_vars(land1_ed = landscape_var("land1", metric = "ed"))
# free_inputs <- kernel_prep(pts = pts, raster_stack = land_rast["land1"],
#                            max_D = 1200, scale_vars = free_vars,
#                            verbose = FALSE)
# sim_df  <- data.frame(counts = sim_counts, free_inputs$kernel_dat)
# sim_mod <- glm(counts ~ land1_ed, family = poisson(), data = sim_df)
# sim_opt <- multiScale_optim(fitted_mod = sim_mod, kernel_inputs = free_inputs,
#                             verbose = FALSE)
# summary(sim_opt, profile = TRUE)$opt_scale     # radius with profile limits
# plot(profile_sigma(sim_opt, n_pts = 15, verbose = FALSE))   # AICc by radius
# grid_fit <- multiScale_optim(sim_mod, free_inputs, backend = "refit",
#                              radius_grid = seq(200, 1000, by = 100))
# diagnostics(grid_fit)$radius_profile        # log-likelihood at each listed radius
# plot_marginal_effects(sim_opt)                 # expected count against m/ha
# 
# # 7. Map the covariates of a fitted model, or a metric on its own
# projected <- kernel_scale.raster(raster_stack = metric_rasters,
#                                  multiScaleR = opt, verbose = FALSE)
# ed_multi  <- kernel_scale.raster(
#   raster_stack = metric_rasters["land1"],
#   scale_vars   = msr_vars(
#     ed_200 = landscape_var("land1", metric = "ed", radius = 200),
#     ed_500 = landscape_var("land1", metric = "ed", radius = 500),
#     ed_900 = landscape_var("land1", metric = "ed", radius = 900)),
#   verbose = FALSE)

