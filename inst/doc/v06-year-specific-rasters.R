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

## ----annual-data--------------------------------------------------------------
data("surv_pts")
pts  <- vect(surv_pts)     # 100 survey locations as a terra SpatVector
land <- rast(system.file("extdata", "landscape.tif",
                         package = "multiScaleR"))

# The bundled objects carry no coordinate reference system (CRS).
# Assign the same projected CRS, in meters, to both.
crs(pts)  <- "EPSG:26915"
crs(land) <- "EPSG:26915"

# Year 1: the binary habitat layer (1 = habitat, 0 = not habitat)
habitat_1 <- land$land1
names(habitat_1) <- "habitat"

# Year 2: the same map after habitat is cleared from one square
habitat_2 <- habitat_1
cleared   <- ext(1000, 2600, 1400, 3000)
habitat_2[cells(habitat_2, cleared)] <- 0

# One row per observation, with the year of its survey: 50 points per year
set.seed(93)
obs <- data.frame(year = factor(sample(rep(c("year1", "year2"), each = 50))))

## ----plot-maps, fig.width = 7.5, fig.height = 3.6, fig.cap = "***Habitat (yellow) in each year, with the points surveyed that year. The red square in year 2 marks the cleared area.***"----
op <- par(mfrow = c(1, 2))
plot(habitat_1, main = "Year 1", legend = FALSE)
points(pts[obs$year == "year1", ], pch = 19, cex = 0.5, col = "red")
plot(habitat_2, main = "Year 2", legend = FALSE)
points(pts[obs$year == "year2", ], pch = 19, cex = 0.5, col = "red")
lines(vect(cleared), col = "red", lwd = 2)   # outline of the cleared square
par(op)

## ----simulate-counts----------------------------------------------------------
# Kernel-weighted habitat at the true scale (sigma = 250 m) on each map
true_1 <- kernel_scale.raster(habitat_1, sigma = 250, kernel = "gaussian",
                              verbose = FALSE)
true_2 <- kernel_scale.raster(habitat_2, sigma = 250, kernel = "gaussian",
                              verbose = FALSE)

# Read each point's value from the map of its own year
hab <- ifelse(obs$year == "year1",
              extract(true_1, pts)[, 2],
              extract(true_2, pts)[, 2])
hab_z <- (hab - mean(hab)) / sd(hab)    # standardize across all 100 rows

# Poisson counts: habitat lowers the expected count, year 2 raises it
obs$count <- rpois(nrow(obs),
                   exp(0.25 - 0.5 * hab_z + 0.3 * (obs$year == "year2")))

## ----grouped-preparation------------------------------------------------------
prepared <- kernel_prep_by_group(
  pts             = pts,                  # one point per observation
  raster_stacks   = list(year1 = habitat_1, year2 = habitat_2),
  group           = obs$year,             # which map each observation uses
  max_D           = 1700,                 # largest distance (m) considered
  kernel          = "gaussian",
  store_cell_data = FALSE,                # keep only the distance-binned summaries
  verbose         = FALSE
)

head(prepared$kernel_dat, 3)              # starting covariate values
table(prepared$raster_group)              # observations assigned to each map

## ----grouped-fit--------------------------------------------------------------
model_data <- data.frame(obs, prepared$kernel_dat)

# Starting model: one habitat slope, plus a shift in the baseline for year 2
start <- glm(count ~ habitat + year, family = poisson(), data = model_data)

# Estimate the shared scale of effect and the coefficients
fit <- multiScale_optim(start, prepared, n_cores = 1, verbose = FALSE)

est     <- summary(fit)$opt_scale         # sigma, SE, and 95% interval (m)
dist90  <- kernel_dist(fit)               # radius holding 90% of kernel weight (m)
ci_beta <- confint.default(fit$opt_mod)   # Wald 95% intervals, log scale

check <- data.frame(
  Parameter = c("sigma, habitat (m)", "effect of habitat", "effect of year 2"),
  Truth     = c(250, -0.5, 0.3),
  Estimate  = c(est[, "Mean"], coef(fit$opt_mod)[c("habitat", "yearyear2")]),
  Lower     = c(est[, "2.5%"], ci_beta[c("habitat", "yearyear2"), 1]),
  Upper     = c(est[, "97.5%"], ci_beta[c("habitat", "yearyear2"), 2])
)
knitr::kable(check, digits = 2, row.names = FALSE)

# Did every prepared observation enter the fitted model?
unlist(diagnostics(fit)$sample_size[c("fitted_n", "prepared_n")])

## ----annual-projection, fig.width = 7.5, fig.height = 3.6, fig.cap = "***Expected count across the landscape in each year, on a common color scale.***"----
surface_1 <- kernel_scale.raster(habitat_1, multiScaleR = fit,
                                 scale_center = TRUE,
                                 site_values  = list(year = "year1"),
                                 verbose = FALSE)
surface_2 <- kernel_scale.raster(habitat_2, multiScaleR = fit,
                                 scale_center = TRUE,
                                 site_values  = list(year = "year2"),
                                 verbose = FALSE)

# Expected count in every raster cell, for each year
pred_1 <- terra::predict(surface_1, fit$opt_mod, type = "response")
pred_2 <- terra::predict(surface_2, fit$opt_mod, type = "response")

# Draw both years on one color scale
zlim <- range(minmax(c(pred_1, pred_2)))
op <- par(mfrow = c(1, 2))
plot(pred_1, main = "Year 1", range = zlim)
plot(pred_2, main = "Year 2", range = zlim)
par(op)

## ----annual-ratio, fig.cap = "***Expected count in year 2 divided by expected count in year 1. The red square marks the cleared area.***"----
# How many times higher is the expected count in year 2?
ratio <- pred_2 / pred_1
plot(ratio, main = "Year 2 / year 1")
lines(vect(cleared), col = "red", lwd = 2)   # outline of the cleared square

round(minmax(ratio)[, 1], 2)

## ----quick-reference, eval = FALSE--------------------------------------------
# library(multiScaleR)
# library(terra)
# 
# # 1. Align the annual maps: one grid, one projected CRS, matching layer names
# data("surv_pts")
# pts  <- vect(surv_pts)
# land <- rast(system.file("extdata", "landscape.tif",
#                          package = "multiScaleR"))
# crs(pts)  <- "EPSG:26915"
# crs(land) <- "EPSG:26915"
# habitat_1 <- land$land1
# names(habitat_1) <- "habitat"
# habitat_2 <- habitat_1
# cleared   <- ext(1000, 2600, 1400, 3000)
# habitat_2[cells(habitat_2, cleared)] <- 0
# 
# # 2. Label each observation with its year (and add your response, here `count`)
# set.seed(93)
# obs <- data.frame(year = factor(sample(rep(c("year1", "year2"), each = 50))))
# 
# # 3. Assign each observation to its map and prepare pooled covariates
# prepared <- kernel_prep_by_group(
#   pts = pts, raster_stacks = list(year1 = habitat_1, year2 = habitat_2),
#   group = obs$year, max_D = 1700, kernel = "gaussian",
#   store_cell_data = FALSE
# )
# 
# # 4. Fit one scale of effect and one habitat coefficient across years
# model_data <- data.frame(obs, prepared$kernel_dat)
# start <- glm(count ~ habitat + year, family = poisson(), data = model_data)
# fit   <- multiScale_optim(start, prepared, n_cores = 1)
# summary(fit)                           # shared sigma, distance, coefficients
# kernel_dist(fit)                       # radius holding 90% of the kernel weight
# diagnostics(fit)$sample_size           # rows prepared and rows fitted
# 
# # 5. Map the fitted model onto each year
# surface_1 <- kernel_scale.raster(habitat_1, multiScaleR = fit,
#                                  scale_center = TRUE,
#                                  site_values  = list(year = "year1"))
# surface_2 <- kernel_scale.raster(habitat_2, multiScaleR = fit,
#                                  scale_center = TRUE,
#                                  site_values  = list(year = "year2"))
# pred_1 <- terra::predict(surface_1, fit$opt_mod, type = "response")
# pred_2 <- terra::predict(surface_2, fit$opt_mod, type = "response")
# plot(pred_2 / pred_1)                  # change in expected count between years

