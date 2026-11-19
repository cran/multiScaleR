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

## ----model-call, eval = FALSE-------------------------------------------------
# forecast_fit <- multiScale_bayes(
#   count ~ habitat + (1 | site) + (1 | year),   # site and year baselines
#   inputs = annual$inputs, data = annual$data,
#   scales = list(habitat = ~1 + (1 | year)),    # scale of effect varies by year
#   family = "nbinom2", row_id = "site_year",
#   priors = list(scale_median = 250),           # prior median scale, in meters
#   method = "sample", chains = 4, iter = 8000, warmup = 4000, seed = 630127,
#   control = list(adapt_delta = 0.99, max_treedepth = 12)
# )

## ----load-example-------------------------------------------------------------
# Survey data and prepared landscape inputs for 2021 to 2025
annual <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
                              package = "multiScaleR"))

# The stored fit, with future surveys and their prepared landscape inputs
example <- readRDS(system.file("extdata", "bayesian-prediction-example.rds",
                               package = "multiScaleR"))
forecast_fit <- example$fit

summary(forecast_fit)$diagnostics$passed   # TRUE: sampling checks passed

## ----known-surveys------------------------------------------------------------
# Three surveyed site-years. Leave out the response, which is not an input.
surveys <- annual$data[1:3, c("site", "year", "site_year")]

# Expected count: uncertainty about the mean only
expected <- predict(forecast_fit, surveys, kernel_inputs = annual$inputs,
                    row_id = "site_year")

# Replicated counts: 1,000 simulated surveys, with observation noise
replicated <- simulate(forecast_fit, nsim = 1000, seed = 11, newdata = surveys,
                       kernel_inputs = annual$inputs, row_id = "site_year")
recorded <- t(apply(as.matrix(replicated), 1, quantile,
                    probs = c(0.025, 0.5, 0.975)))

known <- data.frame(expected[c("point_id", "median", "lower", "upper")],
                    recorded[, c(2, 1, 3)])
knitr::kable(known, digits = 1, row.names = FALSE,
             col.names = c("Site and year", "Expected count", "Lower", "Upper",
                           "Recorded count", "Lower", "Upper"))

## ----future-scales------------------------------------------------------------
years <- data.frame(year = c("2025", "2026"))    # 2025 surveyed, 2026 not

# Draw a deviation for the new year from the estimated spread among years
scale_new <- predict_scales(forecast_fit, years, new_levels = "sample",
                            seed = 801)

# Reference: set the new year's deviation to zero
scale_ref <- predict_scales(forecast_fit, years[2, , drop = FALSE],
                            new_levels = "zero")

knitr::kable(rbind(scale_new, scale_ref)[c("group", "distribution", "median",
                                           "lower", "upper", "p_below_min_D",
                                           "p_above_max_D")],
             digits = c(0, 0, 0, 0, 0, 2, 3), row.names = FALSE,
             col.names = c("Year", "Distribution", "Median (m)", "Lower",
                           "Upper", "Below 100 m", "Above 1,200 m"))

## ----new-surveys--------------------------------------------------------------
new_surveys   <- example$newdata         # two known sites, two new, two years
future_inputs <- example$kernel_inputs   # prepared landscape for these 8 rows
head(new_surveys, 4)

## ----new-predictions----------------------------------------------------------
# Expected counts, drawing an effect for each new site and year
future_mean <- predict(forecast_fit, new_surveys,
                       kernel_inputs = future_inputs, row_id = "site_year",
                       new_scale_levels = "sample",
                       new_response_levels = "sample", seed = 802)

# The same rows with new effects set to zero (a typical site and year)
future_typical <- predict(forecast_fit, new_surveys,
                          kernel_inputs = future_inputs, row_id = "site_year",
                          new_scale_levels = "zero",
                          new_response_levels = "zero")

future <- data.frame(future_mean[c("point_id", "median", "lower", "upper")],
                     future_typical[c("lower", "upper")])
knitr::kable(future, digits = 1, row.names = FALSE,
             col.names = c("Site and year", "Expected count", "Lower", "Upper",
                           "Typical lower", "Typical upper"))

attr(future_mean, "new_groups")   # which new groups were drawn

## ----totals-------------------------------------------------------------------
# 1,000 replicate data sets for the eight future surveys
future_counts <- simulate(forecast_fit, nsim = 1000, seed = 803,
                          newdata = new_surveys, kernel_inputs = future_inputs,
                          row_id = "site_year", new_scale_levels = "sample",
                          new_response_levels = "sample")

# Range of counts each survey could record (first four shown)
by_survey <- t(apply(as.matrix(future_counts), 1, quantile,
                     probs = c(0.025, 0.5, 0.975)))
rownames(by_survey) <- new_surveys$site_year
head(by_survey, 4)

# Total over the eight surveys: sum within each replicate, then summarize
total <- quantile(colSums(future_counts), c(0.025, 0.5, 0.975))
total

# Not an interval for the total: the survey limits added together
colSums(by_survey)

## ----map-raster---------------------------------------------------------------
# Rebuild the example's 2025 habitat map (60 x 60 cells of 100 m)
set.seed(521)
blank <- rast(nrows = 60, ncols = 60, xmin = 0, xmax = 6000,
              ymin = 0, ymax = 6000, crs = "EPSG:26915")
xy <- xyFromCell(blank, seq_len(ncell(blank)))
for (j in 1:5) {            # the five yearly maps share one random stream
  habitat_2025 <- setValues(blank, sin(xy[, 1] / 400 + j / 4) +
                              cos(xy[, 2] / 550) + rnorm(nrow(xy), sd = 0.2))
}
names(habitat_2025) <- "habitat"

## ----map-predict--------------------------------------------------------------
# A grid of locations every 200 m across the surveyed area
grid     <- expand.grid(x = seq(1400, 4600, by = 200),
                        y = seq(1400, 4600, by = 200))
grid_pts <- vect(grid, geom = c("x", "y"), crs = "EPSG:26915")

# Prepare the landscape around each location, with the settings used in fitting
grid_inputs <- kernel_prep(pts = grid_pts, raster_stack = habitat_2025,
                           max_D = 1200, kernel = "gaussian", sigma = 250,
                           nbins = 48, store_cell_data = FALSE, verbose = FALSE)

# Each location is a new site in a surveyed year
grid_data <- data.frame(cell = grid_inputs$binned$point_ids,
                        site = paste0("cell_", seq_len(nrow(grid))),
                        year = "2025")

# Expected count at a typical site (site effect set to zero)
grid_pred <- predict(forecast_fit, grid_data, kernel_inputs = grid_inputs,
                     row_id = "cell", new_response_levels = "zero")

# Put the posterior median and the ratio of the interval limits on a raster
map <- rast(cbind(grid, median = grid_pred$median,
                  ratio = grid_pred$upper / grid_pred$lower), type = "xyz")

## ----map-plot, fig.width = 7, fig.height = 3.4, fig.cap = "***Expected count in 2025 at a typical site (left) and the ratio of its upper to its lower 95% limit (right). Points are the survey sites.***"----
# The 25 survey sites sit on a regular 800 m grid
sites <- vect(expand.grid(x = seq(1400, 4600, by = 800),
                          y = seq(1400, 4600, by = 800)),
              geom = c("x", "y"), crs = "EPSG:26915")
plot(map, main = c("Expected count", "Upper limit / lower limit"),
     mar = c(2, 2.5, 2, 5), pax = list(cex.axis = 0.7),
     fun = function() points(sites, pch = 19, cex = 0.5))

## ----quick-reference, eval = FALSE--------------------------------------------
# # 0. Stored fit, surveys to predict, and their prepared landscape inputs
# forecast_fit  <- example$fit
# surveys       <- annual$data[1:3, c("site", "year", "site_year")]
# years         <- data.frame(year = c("2025", "2026"))
# new_surveys   <- example$newdata
# future_inputs <- example$kernel_inputs
# 
# # 1. Known sites and years: expected count, then counts you might record
# expected   <- predict(forecast_fit, surveys, kernel_inputs = annual$inputs,
#                       row_id = "site_year")
# replicated <- simulate(forecast_fit, nsim = 1000, seed = 11, newdata = surveys,
#                        kernel_inputs = annual$inputs, row_id = "site_year")
# 
# # 2. Scale of effect in a surveyed year and a new year
# scale_new <- predict_scales(forecast_fit, years, new_levels = "sample",
#                             seed = 801)                  # new-year forecast
# scale_ref <- predict_scales(forecast_fit, years[2, , drop = FALSE],
#                             new_levels = "zero")         # typical-year reference
# 
# # 3. Expected counts at new sites and in new years
# future_mean <- predict(forecast_fit, new_surveys,
#                        kernel_inputs = future_inputs, row_id = "site_year",
#                        new_scale_levels = "sample",
#                        new_response_levels = "sample", seed = 802)
# attr(future_mean, "new_groups")                          # what was drawn
# 
# # 4. Counts you might record, and a total across surveys
# future_counts <- simulate(forecast_fit, nsim = 1000, seed = 803,
#                           newdata = new_surveys, kernel_inputs = future_inputs,
#                           row_id = "site_year", new_scale_levels = "sample",
#                           new_response_levels = "sample")
# quantile(colSums(future_counts), c(0.025, 0.5, 0.975))   # interval for the total
# 
# # 5. Map: prepare a grid of locations and predict at a typical site
# grid_inputs <- kernel_prep(pts = grid_pts, raster_stack = habitat_2025,
#                            max_D = 1200, kernel = "gaussian", sigma = 250,
#                            nbins = 48, store_cell_data = FALSE, verbose = FALSE)
# grid_pred <- predict(forecast_fit, grid_data, kernel_inputs = grid_inputs,
#                      row_id = "cell", new_response_levels = "zero")

