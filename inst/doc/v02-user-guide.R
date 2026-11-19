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

# Three fits are stored in the package so the vignette builds quickly. The text
# shows the real calls.
# - N-mixture model (about five minutes): inst/extdata/cache-user-guide-nmix.rds
# - Zero-inflated and occupancy models (about 20 s together):
#   inst/extdata/cache-user-guide.rds
# Regenerate them with the `nmix` and `user-guide` sections of
# data-raw/generate_vignette_cache.R.
# Every other model in this vignette runs when it is built.
stored_file <- function(name) {
  path <- file.path("..", "inst", "extdata", name)
  if (!file.exists(path)) {
    path <- system.file("extdata", name, package = "multiScaleR")
  }
  path
}
opt_nmix_stored <- readRDS(stored_file("cache-user-guide-nmix.rds"))
guide_cache     <- readRDS(stored_file("cache-user-guide.rds"))

## ----kernel-shapes, fig.height = 5.5, fig.cap = "***Weight against distance for the four kernels, each drawn with sigma = 300 m. The dashed line marks the radius that holds 90% of the kernel weight.***"----
library(ggplot2)   # adds a title to each panel
library(cowplot)   # arranges the four panels on one page

# plot_kernel() draws a kernel without a fitted model.
# add_label = FALSE leaves out the text under each panel.
k_gaus  <- plot_kernel(sigma = 300, kernel = "gaussian", add_label = FALSE)
k_exp   <- plot_kernel(sigma = 300, kernel = "exp",      add_label = FALSE)
k_expow <- plot_kernel(sigma = 300, kernel = "expow",    add_label = FALSE,
                       shape = 5)         # shape of the exponential power curve
k_fixed <- plot_kernel(sigma = 300, kernel = "fixed",    add_label = FALSE)

plot_grid(k_gaus  + ggtitle("Gaussian"),
          k_exp   + ggtitle("Negative exponential"),
          k_expow + ggtitle("Exponential power, shape 5"),
          k_fixed + ggtitle("Fixed buffer"))

## ----kernel-dist--------------------------------------------------------------
# Radius (m) of the circle holding 90% of the weight, sigma = 300 m
kernel_dist(sigma = 300, kernel = "gaussian")
kernel_dist(sigma = 300, kernel = "exp")

# The same Gaussian kernel with 95% of the weight
kernel_dist(sigma = 300, kernel = "gaussian", prob = 0.95)

# Line-based definition: distance holding 90% of the weight along a transect
kernel_dist(sigma = 300, kernel = "gaussian", dimension = "1d")

## ----load-data----------------------------------------------------------------
data("landscape_counts")   # counts and one site-level covariate, 100 rows
dat <- landscape_counts

data("surv_pts")           # survey locations (an sf object from the package)
pts <- vect(surv_pts)      # convert to a terra SpatVector

# Three landscape layers on one 20 m grid
land_rast <- rast(system.file("extdata", "landscape.tif",
                              package = "multiScaleR"))

## ----kernel-prep--------------------------------------------------------------
kernel_inputs <- kernel_prep(
  pts          = pts,          # survey locations
  raster_stack = land_rast,    # landscape layers to summarize
  max_D        = 1700,         # largest distance (m) searched around a point
  kernel       = "gaussian",   # shape of the distance-decay function
  verbose      = FALSE
)
kernel_inputs

## ----lean-prep----------------------------------------------------------------
# Same preparation, keeping only the distance bins
lean_inputs <- kernel_prep(pts = pts, raster_stack = land_rast,
                           max_D = 1700, kernel = "gaussian",
                           store_cell_data = FALSE, verbose = FALSE)

# Size of each object in megabytes
c(full = as.numeric(object.size(kernel_inputs)),
  lean = as.numeric(object.size(lean_inputs))) / 1e6

## ----start-model--------------------------------------------------------------
df <- data.frame(dat, kernel_inputs$kernel_dat)

# Starting model with all three landscape layers
mod_full <- glm(counts ~ site + land1 + land2 + land3,
                family = poisson(),
                data   = df)

## ----optim-full---------------------------------------------------------------
opt_full <- multiScale_optim(fitted_mod    = mod_full,
                             kernel_inputs = kernel_inputs,
                             verbose       = FALSE)   # no progress messages

## ----ram, eval = FALSE--------------------------------------------------------
# # Estimated memory for three workers, and the largest safe worker count
# estimate_multiscale_ram(kernel_inputs, fitted_mod = mod_full, n_cores = 3)

## ----summary-full-------------------------------------------------------------
summary(opt_full)

## ----optim-reduced------------------------------------------------------------
# Starting model without land3 (the Quick Start model)
mod <- glm(counts ~ site + land1 + land2,
           family = poisson(),
           data   = df)

opt <- multiScale_optim(fitted_mod    = mod,
                        kernel_inputs = kernel_inputs,
                        verbose       = FALSE)

opt$scale_est                      # sigma and its standard error, in meters
round(kernel_dist(opt))            # radius holding 90% of the kernel weight (m)
round(exp(coef(opt$opt_mod)), 2)   # multiplicative change in expected count

## ----profile-interval---------------------------------------------------------
# Delta-method limits (the default), in meters
round(summary(opt)$opt_scale)

# Profile-likelihood limits. The standard error column does not change
round(summary(opt, profile = TRUE)$opt_scale)

## ----plot-kernels, fig.cap = "***Fitted kernels for `land1` and `land2`.***"----
plot(opt)

## ----plot-kernels-99, fig.cap = "***The same kernels with the radius that holds 99% of the weight, and no text under the panels.***"----
# prob changes the share of weight. dimension = "1d" gives the line-based distance
plot(opt, prob = 0.99, add_label = FALSE)

## ----marginal, fig.cap = "***Expected count across each covariate in the reduced model.***"----
plot_marginal_effects(opt)

## ----profile, fig.cap = "***AICc across candidate sigma values for each covariate in the full model. The dashed line marks the estimate.***"----
prof_full <- profile_sigma(opt_full, n_pts = 15,   # 15 sigma values per covariate
                           verbose = FALSE)
plot(prof_full)

## ----start-strategy-----------------------------------------------------------
# One kernel_prep() call per kernel
fixed_inputs <- kernel_prep(pts = pts, raster_stack = land_rast,
                            max_D = 1700, kernel = "fixed", verbose = FALSE)

# One start, the default for the smooth kernels
opt_fixed_single <- multiScale_optim(mod, fixed_inputs, verbose = FALSE,
                                     start_strategy = "single")

# Default for a fixed kernel: scan sigma values, run short trial fits, then
# optimize once from the best of them
set.seed(1)                    # the trial starts are jittered at random
opt_fixed <- multiScale_optim(mod, fixed_inputs, verbose = FALSE)

# Buffer radius (m) and log-likelihood from each strategy
data.frame(start  = c("single", "screen"),
           land1  = c(opt_fixed_single$scale_est["land1", "Mean"],
                      opt_fixed$scale_est["land1", "Mean"]),
           land2  = c(opt_fixed_single$scale_est["land2", "Mean"],
                      opt_fixed$scale_est["land2", "Mean"]),
           logLik = c(logLik(opt_fixed_single$opt_mod),
                      logLik(opt_fixed$opt_mod)))

## ----compare-kernels----------------------------------------------------------
# Prepare the same data with the two remaining kernels
exp_inputs   <- kernel_prep(pts = pts, raster_stack = land_rast,
                            max_D = 1700, kernel = "exp", verbose = FALSE)
expow_inputs <- kernel_prep(pts = pts, raster_stack = land_rast,
                            max_D = 1700, kernel = "expow", verbose = FALSE)

# Refit the reduced model with each kernel
opt_exp   <- multiScale_optim(mod, exp_inputs,   verbose = FALSE)
opt_expow <- multiScale_optim(mod, expow_inputs, verbose = FALSE)

# Rank the four kernels by AICc. opt_fixed is the screened fit from above
aic_tab(list(opt, opt_exp, opt_expow, opt_fixed),
        mod_names = c("Gaussian", "Negative exponential",
                      "Exponential power", "Fixed buffer"))

## ----compare-structures-------------------------------------------------------
# Candidate starting models
mod_land  <- glm(counts ~ land1 + land2, family = poisson(), data = df)
mod_land1 <- glm(counts ~ land1,         family = poisson(), data = df)
mod_land2 <- glm(counts ~ land2,         family = poisson(), data = df)
mod_site  <- glm(counts ~ site,          family = poisson(), data = df)

# Estimate the scales of effect for each landscape model
opt_land  <- multiScale_optim(mod_land,  kernel_inputs, verbose = FALSE)
opt_land1 <- multiScale_optim(mod_land1, kernel_inputs, verbose = FALSE)
opt_land2 <- multiScale_optim(mod_land2, kernel_inputs, verbose = FALSE)

model_set <- list(opt_full, opt, opt_land, opt_land1, opt_land2, mod_site)
model_names <- c("site + land1 + land2 + land3", "site + land1 + land2",
                 "land1 + land2", "land1", "land2", "site")

aic_tab(model_set, mod_names = model_names)   # ranked by AICc
bic_tab(model_set, mod_names = model_names)   # ranked by BIC

## ----negative-binomial--------------------------------------------------------
# Negative binomial starting model for the running example
mod_nb <- MASS::glm.nb(counts ~ site + land1 + land2, data = df)
opt_nb <- multiScale_optim(mod_nb, kernel_inputs, verbose = FALSE)

opt_nb$scale_est          # sigma (m), nearly identical to the Poisson fit
opt_nb$opt_mod$theta      # dispersion: large values mean little overdispersion

# Is the extra parameter worth it?
aic_tab(list(opt, opt_nb), mod_names = c("Poisson", "Negative binomial"))

## ----zinb-fit, eval = FALSE---------------------------------------------------
# set.seed(12345)
# zi_rast <- sim_rast(user_seed = 999, dim = 500, resolution = 30)
# zi_rast <- subset(zi_rast, c("bin1", "cont1"))   # one binary, one continuous layer
# 
# # 400 random points at least 1,000 m inside the raster edge
# interior <- buffer(vect(ext(zi_rast)), -1000)
# zi_pts   <- spatSample(zi_rast[[1]], size = 400, method = "random",
#                        ext = ext(interior), as.points = TRUE)
# 
# # Covariates weighted at the true scales: 400 m for bin1, 200 m for cont1
# truth <- kernel_prep(pts = zi_pts, raster_stack = zi_rast, max_D = 1500,
#                      kernel = "gaussian", sigma = c(400, 200),
#                      verbose = FALSE)$kernel_dat
# 
# # Excess zeros depend on cont1. Counts depend on bin1
# excess_zero <- rbinom(400, size = 1, prob = plogis(-0.25 - 1 * truth$cont1))
# abundance   <- rnbinom(400, mu = exp(0.5 + 0.75 * truth$bin1), size = 2)
# zi_dat      <- data.frame(cnt = ifelse(excess_zero == 1, 0, abundance))
# 
# # Prepare the inputs without telling the package the true scales
# zi_inputs <- kernel_prep(pts = zi_pts, raster_stack = zi_rast, max_D = 1500,
#                          kernel = "gaussian", verbose = FALSE)
# zi_df <- data.frame(zi_dat, zi_inputs$kernel_dat)
# 
# # Count model | zero-inflation model
# zinb_mod <- pscl::zeroinfl(cnt ~ bin1 | cont1, dist = "negbin", data = zi_df)
# zinb_opt <- multiScale_optim(zinb_mod, zi_inputs, verbose = FALSE)

## ----zinb-load, include = FALSE-----------------------------------------------
zinb_opt <- guide_cache$zinb_opt

## ----zinb-read----------------------------------------------------------------
round(summary(zinb_opt)$opt_scale)      # sigma (m): truth is 400 and 200
round(coef(zinb_opt$opt_mod), 2)        # truth: 0.5, 0.75, -0.25, -1

## ----binomial-----------------------------------------------------------------
df$present <- as.integer(df$counts > 0)    # 1 when any animal was counted
table(df$present)

mod_bin <- glm(present ~ site + land1 + land2,
               family = binomial(link = "logit"),
               data   = df)
opt_bin <- multiScale_optim(mod_bin, kernel_inputs, verbose = FALSE)

round(summary(opt_bin)$opt_scale)          # sigma, SE, and 95% limits (m)
round(exp(coef(opt_bin$opt_mod)), 2)       # odds ratios

## ----nmix-data, fig.height = 3.5, fig.cap = "***The two simulated layers used for the `unmarked` examples.***"----
library(unmarked)

# Four simulated layers on a 10 m grid. Keep one binary and one continuous
um_rast <- subset(sim_rast(user_seed = 999, dim = 250), c("bin1", "cont2"))
plot(um_rast)

# 75 sites, 5 surveys each, individual detection probability 0.5
nmix_sim <- sim_dat_unmarked(alpha = 1,                # abundance intercept (log)
                             beta  = c(0.75, -0.75),   # effects of bin1, cont2
                             sigma = c(75, 150),       # true scales (m)
                             kernel = "gaussian", type = "count",
                             n_points = 75, n_surv = 5, det = 0.5,
                             raster_stack = um_rast, max_D = 550,
                             user_seed = 555)

nmix_inputs <- kernel_prep(pts = nmix_sim$pts, raster_stack = um_rast,
                           max_D = 550, kernel = "gaussian", verbose = FALSE)

# Counts (sites by surveys) with the kernel-weighted site covariates
nmix_umf <- unmarkedFramePCount(y = nmix_sim$y,
                                siteCovs = nmix_inputs$kernel_dat)

# Starting model: ~ detection ~ abundance
nmix_mod <- pcount(~ 1 ~ bin1 + cont2, data = nmix_umf, K = 100)

## ----nmix-optim, eval = FALSE-------------------------------------------------
# opt_nmix <- multiScale_optim(fitted_mod    = nmix_mod,
#                              kernel_inputs = nmix_inputs)

## ----nmix-load, include = FALSE-----------------------------------------------
opt_nmix <- opt_nmix_stored
# The stored fit came from this call under other object names. Show the call
# as it is written above.
opt_nmix$call <- quote(multiScale_optim(fitted_mod = nmix_mod,
                                        kernel_inputs = nmix_inputs))

## ----nmix-summary-------------------------------------------------------------
summary(opt_nmix)

# Detection is reported on the logit scale. Convert it to a probability
plogis(coef(opt_nmix$opt_mod, type = "det"))

## ----nmix-compare-------------------------------------------------------------
# K in this table counts parameters: 3 abundance, 1 detection, 2 sigma
aic_tab(list(opt_nmix), mod_names = "bin1 + cont2")

# K in pcount() is the upper limit of the abundance sum. Check that the
# log-likelihood does not change when you raise it
c(K_100 = logLik(nmix_mod), K_200 = logLik(update(nmix_mod, K = 200)))

## ----occupancy, eval = FALSE--------------------------------------------------
# # 250 sites, 5 surveys each. Occupancy effects on the logit scale
# occ_sim <- sim_dat_unmarked(alpha = 0.25, beta = c(-0.75, 1),
#                             sigma = c(225, 75),      # true scales (m)
#                             kernel = "gaussian", type = "occ",
#                             n_points = 250, n_surv = 5, det = 0.5,
#                             raster_stack = um_rast, max_D = 800,
#                             user_seed = 555)
# 
# occ_inputs <- kernel_prep(pts = occ_sim$pts, raster_stack = um_rast,
#                           max_D = 800, kernel = "gaussian", verbose = FALSE)
# occ_umf <- unmarkedFrameOccu(y = occ_sim$y, siteCovs = occ_inputs$kernel_dat)
# 
# # Starting model: ~ detection ~ occupancy
# occ_mod <- occu(~ 1 ~ bin1 + cont2, data = occ_umf)
# opt_occ <- multiScale_optim(occ_mod, occ_inputs, verbose = FALSE)

## ----occupancy-load, include = FALSE------------------------------------------
opt_occ <- guide_cache$opt_occ

## ----occupancy-read-----------------------------------------------------------
round(summary(opt_occ)$opt_scale)                 # sigma (m): truth is 225 and 75
round(coef(opt_occ$opt_mod), 2)                   # logit scale
plogis(coef(opt_occ$opt_mod, type = "det"))       # detection probability

## ----project, fig.cap = "***Layers for prediction from the reduced count model, in standard deviation units.***"----
r_scaled <- kernel_scale.raster(raster_stack = land_rast,
                                multiScaleR  = opt,     # the fitted object
                                scale_center = TRUE,    # same scaling as the model
                                clamp        = TRUE,    # stay within the fitted range
                                verbose      = FALSE)
plot(r_scaled)

## ----predict, fig.cap = "***Expected count across the landscape, with survey points.***"----
# Expected count in every raster cell
pred <- terra::predict(r_scaled, opt$opt_mod, type = "response")
plot(pred, main = "Expected count")
points(pts, pch = 19, cex = 0.5, col = "red")

## ----predict-unmarked, fig.height = 5.5, fig.cap = "***Expected abundance from the N-mixture model, with its standard error and 95% confidence limits.***"----
um_scaled <- kernel_scale.raster(raster_stack = um_rast, multiScaleR = opt_nmix,
                                 scale_center = TRUE, clamp = TRUE,
                                 verbose = FALSE)

# type = "state" is abundance for pcount() and occupancy probability for occu()
predict_state <- function(model, data, ...) {
  unmarked::predict(model, type = "state", newdata = data, ...)
}
abundance_map <- terra::predict(um_scaled, opt_nmix$opt_mod,
                                fun = predict_state, na.rm = TRUE)
plot(abundance_map)

## ----quick-reference, eval = FALSE--------------------------------------------
# library(multiScaleR)
# library(terra)
# 
# # 1. Understand the kernel
# plot_kernel(sigma = 300, kernel = "gaussian", add_label = FALSE)
# kernel_dist(sigma = 300, kernel = "gaussian")      # radius with 90% of the weight
# 
# # 2. Prepare the data
# data("landscape_counts"); data("surv_pts")
# dat       <- landscape_counts
# pts       <- vect(surv_pts)
# land_rast <- rast(system.file("extdata", "landscape.tif",
#                               package = "multiScaleR"))
# kernel_inputs <- kernel_prep(pts = pts, raster_stack = land_rast,
#                              max_D = 1700, kernel = "gaussian")
# 
# # 3. Fit a starting model and estimate the scales
# df  <- data.frame(dat, kernel_inputs$kernel_dat)
# mod <- glm(counts ~ site + land1 + land2, family = poisson(), data = df)
# opt <- multiScale_optim(fitted_mod = mod, kernel_inputs = kernel_inputs)
# 
# # 4. Read and check the fit
# summary(opt)                                # sigma, distance, coefficients, warnings
# summary(opt, profile = TRUE)                # profile-likelihood limits for sigma
# kernel_dist(opt)                            # effective distance (map units)
# plot(opt)                                   # fitted kernels
# plot(opt, prob = 0.99, add_label = FALSE)   # 99% radius, no text under panels
# plot_marginal_effects(opt)                  # expected response across covariates
# plot(profile_sigma(opt, n_pts = 15))        # is each scale identified?
# 
# # 5. Guard against a stalled search (a fixed kernel screens starts by default)
# fixed_inputs <- kernel_prep(pts = pts, raster_stack = land_rast,
#                             max_D = 1700, kernel = "fixed")
# set.seed(1)
# opt_fixed <- multiScale_optim(mod, fixed_inputs)
# 
# # 6. Compare models (each sigma counts as a parameter)
# mod_land <- glm(counts ~ land1 + land2, family = poisson(), data = df)
# opt_land <- multiScale_optim(mod_land, kernel_inputs)
# aic_tab(list(opt, opt_land), mod_names = c("site + land1 + land2",
#                                            "land1 + land2"))
# 
# # 7. Other families: same workflow, different starting model
# mod_nb     <- MASS::glm.nb(counts ~ site + land1 + land2, data = df)
# df$present <- as.integer(df$counts > 0)
# mod_bin    <- glm(present ~ site + land1 + land2,
#                   family = binomial(link = "logit"), data = df)
# opt_nb     <- multiScale_optim(mod_nb,  kernel_inputs)
# opt_bin    <- multiScale_optim(mod_bin, kernel_inputs)
# # pscl::zeroinfl(cnt ~ bin1 | cont1, dist = "negbin", data = zi_df)
# # unmarked::pcount(~ 1 ~ bin1 + cont2, data = nmix_umf, K = 100)
# # unmarked::occu(~ 1 ~ bin1 + cont2, data = occ_umf)
# 
# # 8. Map the fitted model
# r_scaled <- kernel_scale.raster(raster_stack = land_rast, multiScaleR = opt,
#                                 scale_center = TRUE, clamp = TRUE)
# pred <- terra::predict(r_scaled, opt$opt_mod, type = "response")

