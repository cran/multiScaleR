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

# The spatial abundance field section shows its calls without running them,
# because those functions need specific versions of unmarked and TMB. Its
# printed table comes from the `distance-removal` section of
# data-raw/generate_vignette_cache.R.
cache_file <- file.path("..", "inst", "extdata",
                        "cache-distance-removal.rds")
if (!file.exists(cache_file)) {
  cache_file <- system.file("extdata", "cache-distance-removal.rds",
                            package = "multiScaleR")
}
field <- readRDS(cache_file)$field

## ----load-data----------------------------------------------------------------
data("surv_pts")                 # 100 survey locations (an sf object)
pts <- vect(surv_pts)            # convert to a terra SpatVector

# Binary habitat layer on a 20 m grid, in meters
land1 <- rast(system.file("extdata", "landscape.tif",
                          package = "multiScaleR"))[["land1"]]

# Step 1 of every analysis: extract the raster around each survey point
kernel_inputs <- kernel_prep(pts = pts, raster_stack = land1,
                             max_D = 1700, kernel = "gaussian",
                             verbose = FALSE)

## ----simulate-----------------------------------------------------------------
# Habitat smoothed at the true scale, then standardized (mean 0, SD 1)
land1_250 <- kernel_scale.raster(raster_stack = land1, sigma = 250,
                                 kernel = "gaussian", verbose = FALSE)
hab_true  <- as.numeric(scale(terra::extract(land1_250, pts, ID = FALSE)[, 1]))

sim <- sim_dat_gdr(
  density          = exp(log(400) + 0.5 * hab_true),  # birds per square km
  distance_scale   = 40,                # half-normal detection scale (m)
  removal          = 0.2,               # chance per minute of a first detection
  dist_breaks      = c(0, 25, 50, 75),  # three bands within a 75 m radius
  interval_lengths = c(2, 3, 5),        # a 10-minute count in three intervals
  interval_units   = "minutes",
  site_covs        = kernel_inputs$kernel_dat,   # covariates for the model
  event_keys       = data.frame(survey = paste0("pt", 1:100)),
  seed             = 2
)
umf <- sim$observed              # the unmarkedFrameGDR to analyze

head(umf@yDistance, 3)           # birds by distance band
head(umf@yRemoval, 3)            # the same birds by time interval
sim$truth[1:3, c("expected_density", "expected_abundance",
                 "realized_abundance", "realized_observed")]

## ----histories----------------------------------------------------------------
# Where each column of the record matrix sits: 2 time intervals by 3 distance
# bands (0-25 m, 25-75 m, and beyond the 75 m radius)
cell_map <- expand.grid(interval = 1:2, distance_bin = 1:3)
cell_map$column <- paste0("cell", 1:6)

# One row per bird, with a 1 where the bird was recorded
records <- matrix(0, nrow = 3, ncol = 6, dimnames = list(NULL, cell_map$column))
records[1, "cell1"] <- 1                 # interval 1, 0-25 m
records[2, "cell4"] <- 1                 # interval 2, 25-75 m
records[3, c("cell5", "cell2")] <- 1     # beyond 75 m first, then 0-25 m

linked <- prepare_gdr_histories(
  histories    = records,
  history_keys = data.frame(survey = c("a", "a", "b"),
                            bird   = c("one", "two", "one")),
  event_keys   = data.frame(survey = c("a", "b", "c")),   # every survey done
  event_id = "survey", individual_id = "bird", cell_map = cell_map,
  dist_breaks = c(0, 25, 75, Inf), radius = 75, units_in = "m",
  interval_lengths = c(5, 5), interval_units = "minutes",
  surveyed = c(TRUE, TRUE, TRUE)         # survey c found no birds
)

linked$first[, c(".history_row", "first_interval", "first_distance_bin",
                 "status")]
linked$yDistance
linked$yRemoval

# The two tables go straight into an unmarked frame
toy_umf <- unmarked::unmarkedFrameGDR(
  yDistance = linked$yDistance, yRemoval = linked$yRemoval,
  dist.breaks = c(0, 25, 75), unitsIn = "m", period.lengths = c(5, 5)
)

## ----validate-----------------------------------------------------------------
checked <- validate_unmarked_data(
  object = umf, event_keys = sim$event_keys, event_id = "survey",
  interval_units = "minutes",
  effort = rep(10, 100)               # minutes surveyed at each point
)
checked$valid                         # TRUE: no errors were found
nrow(checked$findings)                # number of errors and warnings

## ----intervals----------------------------------------------------------------
umf_obs <- umf
# Who ran each count: a covariate recorded once per survey
unmarked::siteCovs(umf_obs)$observer <- factor(rep(c("Ann", "Bob"), times = 50))

expanded <- expand_gdr_intervals(
  object = umf_obs, event_keys = sim$event_keys, event_id = "survey",
  interval_units  = "minutes",
  site_covariates = c(rem_observer = "observer")   # new name = source column
)
head(unmarked::obsCovs(expanded$data), 4)

## ----fit-mod------------------------------------------------------------------
mod <- unmarked::gdistremoval(
  lambdaformula   = ~ land1,      # density depends on the landscape covariate
  removalformula  = ~ 1,          # one removal probability
  distanceformula = ~ 1,          # one distance scale
  data     = umf,
  keyfun   = "halfnorm",          # half-normal detection function
  mixture  = "P",                 # Poisson abundance
  output   = "density",           # report density instead of abundance
  unitsOut = "kmsq",              # birds per square km
  K        = 60,                  # largest abundance considered at a point
  starts   = c(log(400), 0, log(40), qlogis(0.2))  # link-scale starting values
)

## ----optim--------------------------------------------------------------------
opt <- multiScale_optim(fitted_mod    = mod,
                        kernel_inputs = kernel_inputs,
                        n_cores       = 1,
                        verbose       = FALSE)

## ----summary------------------------------------------------------------------
summary(opt)

## ----components---------------------------------------------------------------
fit <- opt$opt_mod                 # the model refitted at the estimated scale

dens <- unmarked::predict(fit, type = "lambda")   # one row per survey point
dsig <- unmarked::predict(fit, type = "dist")     # one row per survey point
prem <- unmarked::predict(fit, type = "rem")      # one row per point and interval

# First survey point
rbind(density = dens[1, ], distance = dsig[1, ], removal = prem[1, ])

# Turn the two detection components into probabilities for a 10-minute count
p_near  <- 2 * dsig$Predicted[1]^2 / 75^2 *
  (1 - exp(-75^2 / (2 * dsig$Predicted[1]^2)))    # average over the 75 m circle
p_heard <- 1 - (1 - prem$Predicted[1])^10         # detected at least once
round(c(distance = p_near, removal = p_heard, both = p_near * p_heard), 2)

## ----truth-table--------------------------------------------------------------
est     <- summary(opt)$opt_scale                    # sigma and 95% limits (m)
b       <- unmarked::coef(fit)                       # link-scale coefficients
ci_lam  <- unmarked::confint(fit, type = "lambda")   # Wald limits, log scale
ci_dist <- unmarked::confint(fit, type = "dist")     # log scale
ci_rem  <- unmarked::confint(fit, type = "rem")      # logit scale

check <- data.frame(
  Parameter = c("Scale of effect, sigma (m)", "Density at average habitat",
                "Effect of land1 (log scale)", "Distance scale (m)",
                "Removal probability per minute"),
  Truth     = c(250, 400, 0.5, 40, 0.2),
  Estimate  = c(est[1, "Mean"], exp(b[1]), b[2], exp(b[3]), plogis(b[4])),
  Lower     = c(est[1, "2.5%"], exp(ci_lam[1, 1]), ci_lam[2, 1],
                exp(ci_dist[1, 1]), plogis(ci_rem[1, 1])),
  Upper     = c(est[1, "97.5%"], exp(ci_lam[1, 2]), ci_lam[2, 2],
                exp(ci_dist[1, 2]), plogis(ci_rem[1, 2]))
)
# Does each 95% interval contain the value used to simulate the data?
check$Covered <- check$Truth >= check$Lower & check$Truth <= check$Upper
knitr::kable(check, digits = 2, row.names = FALSE)

## ----plot-kernel, fig.cap = "***Kernel weight against distance for `land1`.***"----
plot(opt)

## ----marginal, fig.cap = "***Estimated density across the observed range of kernel-weighted `land1`.***"----
plot_marginal_effects(opt, type = "lambda",
                      ylab = "Density (birds per square km)")

## ----map, fig.cap = "***Predicted density, in birds per square kilometer, with survey points.***"----
# Smooth the raster at the estimated scale, standardized as in the model
r_scaled <- kernel_scale.raster(raster_stack = land1, multiScaleR = opt,
                                scale_center = TRUE, clamp = TRUE,
                                verbose = FALSE)

# Density for every raster cell, with its standard error and 95% limits
dens_map <- terra::predict(r_scaled, fit, type = "lambda")
names(dens_map)

plot(dens_map$Predicted, main = "Density (birds per square km)")
points(pts, pch = 19, cex = 0.5)

## ----refitter, eval = FALSE---------------------------------------------------
# # A second set of starting values, named like the model coefficients
# cold <- unmarked::coef(mod)
# cold[] <- c(log(300), 0, log(30), qlogis(0.4))
# 
# refit <- make_gdr_refitter(model = mod, cold_starts = list(cold = cold),
#                            audit_dir = tempfile("gdr-audit-"))
# 
# # The same search as before, with every candidate scale fitted from both starts
# opt_checked <- multiScale_optim(fitted_mod = mod, kernel_inputs = kernel_inputs,
#                                 refit_fn = refit, n_cores = 1)

## ----field-call, eval = FALSE-------------------------------------------------
# # Template: an independent random intercept for each block (TMB engine)
# template <- unmarked::gdistremoval(
#   lambdaformula = ~ land1 + (1 | block), removalformula = ~ 1,
#   distanceformula = ~ 1, data = umf_blocks, keyfun = "halfnorm",
#   mixture = "P", output = "density", unitsOut = "kmsq", K = 80, engine = "TMB",
#   starts = c(log(400), 0, log(40), qlogis(0.2), log(0.5))
# )
# 
# # Extra starting values, named by parameter block and position
# cold <- c("beta_lambda[1]"   = log(300),      # log density at average habitat
#           "beta_lambda[2]"   = 0,             # land1 effect
#           "beta_dist[1]"     = log(30),       # log distance scale
#           "beta_rem[1]"      = qlogis(0.4),   # logit removal probability
#           "lsigma_lambda[1]" = log(0.3))      # log SD of the block effects
# settings <- list(method = "BFGS", se = TRUE, threads = 1L,
#                  control = list(maxit = 2500L, reltol = 1e-10))
# 
# field_fit <- fit_gdr_response_field(
#   model  = template,
#   coords = block_xy,        # one row per survey: the center of its block
#   crs    = 32617,           # EPSG code of the projected coordinates
#   range  = 1500,            # correlation range, in map units (m)
#   cold_starts = list(cold = cold), tmb_settings = settings,
#   audit_dir   = tempfile("gdr-field-audit-")
# )

## ----field-results, echo = FALSE----------------------------------------------
knitr::kable(field$compare, digits = c(0, 0, 2, 2, 1, 1), row.names = FALSE,
             col.names = c("Model", "Density at average habitat",
                           "Effect of land1", "SD of block effects",
                           "Log likelihood", "AIC"))

## ----quick-reference, eval = FALSE--------------------------------------------
# library(multiScaleR)
# library(terra)
# 
# # 1. Prepare the kernel inputs
# data("surv_pts")
# pts   <- vect(surv_pts)
# land1 <- rast(system.file("extdata", "landscape.tif",
#                           package = "multiScaleR"))[["land1"]]
# kernel_inputs <- kernel_prep(pts = pts, raster_stack = land1,
#                              max_D = 1700, kernel = "gaussian")
# 
# # 2. Simulate a survey (or build the frame from your own records)
# land1_250 <- kernel_scale.raster(raster_stack = land1, sigma = 250,
#                                  kernel = "gaussian")
# hab_true  <- as.numeric(scale(terra::extract(land1_250, pts, ID = FALSE)[, 1]))
# sim <- sim_dat_gdr(density = exp(log(400) + 0.5 * hab_true),
#                    distance_scale = 40, removal = 0.2,
#                    dist_breaks = c(0, 25, 50, 75),
#                    interval_lengths = c(2, 3, 5), interval_units = "minutes",
#                    site_covs = kernel_inputs$kernel_dat,
#                    event_keys = data.frame(survey = paste0("pt", 1:100)),
#                    seed = 2)
# umf <- sim$observed
# 
# # 3. Check the data
# validate_unmarked_data(object = umf, event_keys = sim$event_keys,
#                        event_id = "survey", interval_units = "minutes",
#                        effort = rep(10, 100))$valid
# 
# # 4. Fit a starting model, then estimate the scale of effect
# mod <- unmarked::gdistremoval(lambdaformula = ~ land1, removalformula = ~ 1,
#                               distanceformula = ~ 1, data = umf,
#                               keyfun = "halfnorm", mixture = "P",
#                               output = "density", unitsOut = "kmsq", K = 60,
#                               starts = c(log(400), 0, log(40), qlogis(0.2)))
# opt <- multiScale_optim(fitted_mod = mod, kernel_inputs = kernel_inputs,
#                         n_cores = 1)
# 
# # 5. Read the results
# summary(opt)                                  # sigma, effective distance, model
# fit <- opt$opt_mod
# unmarked::predict(fit, type = "lambda")       # density, birds per square km
# unmarked::predict(fit, type = "dist")         # distance scale, m
# unmarked::predict(fit, type = "rem")          # removal probability per minute
# plot(opt)                                     # kernel weight against distance
# plot_marginal_effects(opt, type = "lambda")   # density across the covariate
# 
# # 6. Map density
# r_scaled <- kernel_scale.raster(raster_stack = land1, multiScaleR = opt,
#                                 scale_center = TRUE, clamp = TRUE)
# dens_map <- terra::predict(r_scaled, fit, type = "lambda")

