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

# glmmTMB is a suggested package. Without it, the code is shown but not run.
has_glmmTMB <- suppressWarnings(requireNamespace("glmmTMB", quietly = TRUE))

knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6.5,
                      fig.height = 4.25, warning = FALSE, message = FALSE,
                      eval = has_glmmTMB)

# The scale search with a spatial field takes most of a minute, so that fit is
# stored in the package. The text shows the real call. Regenerate the stored
# result with the `spatial-fields` section of data-raw/generate_vignette_cache.R.
cache_file <- file.path("..", "inst", "extdata",
                        "cache-spatial-fields.rds")
if (!file.exists(cache_file)) {
  cache_file <- system.file("extdata", "cache-spatial-fields.rds",
                            package = "multiScaleR")
}
field_cache <- readRDS(cache_file)

# Inline results print a placeholder when glmmTMB is absent and no code ran.
val <- function(x) if (has_glmmTMB) x else "(not run)"

## ----no-glmmTMB-note, echo = FALSE, results = "asis", eval = !has_glmmTMB-----
# cat("> **Note.** `glmmTMB` was not available when this vignette was built, so",
#     "the code is shown without its results.")

## ----simulate-----------------------------------------------------------------
data("surv_pts")
pts <- vect(surv_pts)      # 100 survey points
xy  <- crds(pts)           # their coordinates, in meters
land_rast <- rast(system.file("extdata", "landscape.tif",
                              package = "multiScaleR"))

# Landscape covariates smoothed at the true scales, then standardized
truth  <- kernel_scale.raster(land_rast[[c("land1", "land2")]],
                              sigma = c(250, 500), kernel = "gaussian",
                              verbose = FALSE)
x_true <- scale(as.matrix(extract(truth, pts, ID = FALSE)))

set.seed(1)
# The field: correlation between two points is exp(-distance / 800)
d_mat <- as.matrix(dist(xy))
field <- drop(t(chol(0.6^2 * exp(-d_mat / 800))) %*% rnorm(nrow(xy)))
site  <- rnorm(nrow(xy))   # a covariate measured at the survey point

# Log expected count: intercept, site effect, landscape effects, and the field
log_mu <- 1 + 0.3 * site - 0.5 * x_true[, "land1"] +
  0.7 * x_true[, "land2"] + field
dat <- data.frame(counts = rpois(nrow(xy), exp(log_mu)), site = site)

## ----use-stored-data, include = FALSE-----------------------------------------
# The stored fit used these simulated values. Fall back on the stored copy if
# a platform difference in random number generation changed anything.
if (!isTRUE(all.equal(dat, field_cache$dat))) dat <- field_cache$dat

## ----fit-glm------------------------------------------------------------------
# bin = TRUE stores distance-binned summaries, which speed up every refit
inputs <- kernel_prep(pts = pts, raster_stack = land_rast, max_D = 1700,
                      kernel = "gaussian", bin = TRUE, verbose = FALSE)
df <- data.frame(dat, inputs$kernel_dat)

start_glm <- glm(counts ~ site + land1 + land2, family = poisson(), data = df)
fit_glm   <- multiScale_optim(start_glm, inputs, n_cores = 1, verbose = FALSE)

## ----start-field--------------------------------------------------------------
df$pos <- glmmTMB::numFactor(xy[, 1], xy[, 2])   # coordinates, in meters
df$grp <- factor(rep("all", nrow(df)))           # one field for all points

# exp() here is glmmTMB's exponential covariance, not the exponential function
start_field <- glmmTMB::glmmTMB(
  counts ~ site + land1 + land2 + exp(pos + 0 | grp),
  family = poisson(), data = df,
  start  = list(theta = c(log(0.5), log(500)))   # log field SD, log range (m)
)

## ----fit-field, eval = FALSE--------------------------------------------------
# fit_field <- multiScale_optim(start_field, inputs, n_cores = 1, verbose = FALSE)

## ----load-fit-field, include = FALSE------------------------------------------
fit_field <- field_cache$fit_field
# The stored fit leaves out its copy of the kernel inputs to stay small.
fit_field$kernel_inputs <- inputs
# A glmmTMB fit saved under one glmmTMB version must be rebuilt with up2date()
# before another version can use it. Newer versions add internal model terms.
fit_field$opt_mod <- glmmTMB::up2date(fit_field$opt_mod)

## ----compare------------------------------------------------------------------
s_glm   <- summary(fit_glm)$opt_scale                     # sigma, SE, 95% CI (m)
s_field <- summary(fit_field)$opt_scale
b_glm   <- summary(fit_glm$opt_mod)$coefficients          # log-scale estimates
b_field <- summary(fit_field$opt_mod)$coefficients$cond

compare <- data.frame(
  Parameter = c("sigma, land1 (m)", "sigma, land2 (m)", "intercept",
                "effect of site", "effect of land1", "effect of land2"),
  Truth     = c(250, 500, 1, 0.3, -0.5, 0.7),
  No_field  = c(s_glm[, "Mean"], b_glm[, "Estimate"]),
  SE        = c(s_glm[, "SE"], b_glm[, "Std. Error"]),
  Field     = c(s_field[, "Mean"], b_field[, "Estimate"]),
  SE_field  = c(s_field[, "SE"], b_field[, "Std. Error"]))
knitr::kable(compare, digits = 2, row.names = FALSE,
             col.names = c("Parameter", "Truth", "No field", "SE",
                           "With field", "SE"))

## ----field-par----------------------------------------------------------------
# The field parameters, back-transformed from the log scale
field_par <- exp(glmmTMB::getME(fit_field$opt_mod, "theta"))
names(field_par) <- c("field SD (log scale)", "range (m)")
round(field_par, 2)

# K counts the regression coefficients, the field parameters, and the two sigmas
aic_tab(list(fit_glm, fit_field), mod_names = c("No field", "With field"))

## ----predict-points-----------------------------------------------------------
# Default: includes the estimated field value at each survey point
p_field <- predict(fit_field, se.fit = TRUE, uncertainty = "conditional")

# re.form = NA: sets the field to zero, leaving the covariate effects only
p_zero  <- predict(fit_field, re.form = NA)

round(head(data.frame(observed   = df$counts,
                      with_field = as.numeric(p_field$fit),
                      se         = as.numeric(p_field$se.fit),
                      field_zero = as.numeric(p_zero)), 5), 2)

# The estimated field at the survey points against the simulated one
cor(log(as.numeric(p_field$fit) / as.numeric(p_zero)), field)

## ----predict-grid, fig.width = 7, fig.height = 2.8, fig.cap = "***Predictions on a 225 m grid inside the surveyed area.***"----
# A 225 m grid of new locations inside the surveyed area
grid_xy <- as.matrix(expand.grid(x = seq(800, 3500, by = 225),
                                 y = seq(800, 3500, by = 225)))
grid_inputs <- kernel_prep(pts = vect(grid_xy, type = "points"),
                           raster_stack = land_rast, max_D = 1700,
                           kernel = "gaussian", bin = TRUE, verbose = FALSE)

# One row per new location: site covariate, coordinates, and the field group
grid_dat <- data.frame(site = 0,
                       pos  = glmmTMB::numFactor(grid_xy[, 1], grid_xy[, 2]),
                       grp  = factor("all"),
                       row.names = rownames(grid_inputs$kernel_dat))

with_field <- predict(fit_field, grid_dat, grid_inputs,
                      allow.new.levels = TRUE,          # accept new positions
                      se.fit = TRUE, uncertainty = "conditional")
no_field   <- predict(fit_field, grid_dat, grid_inputs, re.form = NA)

grid_map <- rast(data.frame(
  grid_xy,
  expected   = as.numeric(with_field$fit),
  multiplier = as.numeric(with_field$fit) / as.numeric(no_field),
  se         = as.numeric(with_field$se.fit)), type = "xyz")
plot(grid_map, nc = 3,
     main = c("Expected count", "Field multiplier", "Standard error"))

## ----predict-raster, fig.width = 4.5, fig.height = 3.6, fig.cap = "***Expected count from the covariate effects alone, with the field set to zero.***"----
# The bundled raster has no CRS. Declare a projected one, in meters.
map_rast <- land_rast
crs(map_rast) <- "EPSG:32617"

# Landscape layers at the estimated scales, with site held at zero. The field
# columns are placeholders that re.form = NA ignores.
layers <- kernel_scale.raster(
  map_rast, multiScaleR = fit_field,
  scale_center = TRUE, match_buffer = TRUE, clamp = TRUE, verbose = FALSE,
  site_values = list(site = 0, pos = as.character(df$pos[1]), grp = "all"))
map_zero <- predict(fit_field, layers, re.form = NA)
plot(map_zero, main = "Expected count, field set to zero")
points(pts, pch = 19, cex = 0.4)

## ----quick-reference, eval = FALSE--------------------------------------------
# # 1. Prepare kernel inputs and fit a baseline model without a field
# inputs <- kernel_prep(pts = pts, raster_stack = land_rast, max_D = 1700,
#                       kernel = "gaussian", bin = TRUE)
# df <- data.frame(dat, inputs$kernel_dat)
# start_glm <- glm(counts ~ site + land1 + land2, family = poisson(), data = df)
# fit_glm   <- multiScale_optim(start_glm, inputs, n_cores = 1)
# 
# # 2. Add a spatial field on the survey coordinates and estimate the scales
# df$pos <- glmmTMB::numFactor(xy[, 1], xy[, 2])   # coordinates, in meters
# df$grp <- factor(rep("all", nrow(df)))           # one field for all points
# start_field <- glmmTMB::glmmTMB(
#   counts ~ site + land1 + land2 + exp(pos + 0 | grp),
#   family = poisson(), data = df,
#   start  = list(theta = c(log(0.5), log(500))))  # log field SD, log range (m)
# fit_field <- multiScale_optim(start_field, inputs, n_cores = 1)
# 
# # 3. Read the results
# summary(fit_field)$opt_scale                          # sigma, SE, interval (m)
# summary(fit_field$opt_mod)$coefficients$cond          # log-scale coefficients
# exp(glmmTMB::getME(fit_field$opt_mod, "theta"))       # field SD, range (m)
# aic_tab(list(fit_glm, fit_field), mod_names = c("No field", "With field"))
# 
# # 4. Predict at the survey points, then at new locations
# predict(fit_field, se.fit = TRUE, uncertainty = "conditional")  # with the field
# predict(fit_field, re.form = NA)                                # field at zero
# predict(fit_field, grid_dat, grid_inputs, allow.new.levels = TRUE,
#         se.fit = TRUE, uncertainty = "conditional")   # grid_dat: site, pos, grp
# 
# # 5. Map the covariate effects with the field set to zero
# layers <- kernel_scale.raster(
#   map_rast, multiScaleR = fit_field,             # map_rast has a projected CRS
#   scale_center = TRUE, match_buffer = TRUE, clamp = TRUE,
#   site_values = list(site = 0, pos = as.character(df$pos[1]), grp = "all"))
# map_zero <- predict(fit_field, layers, re.form = NA)

