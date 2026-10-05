## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 6, fig.height = 4.25)
library(multiScaleR)
library(terra)

## ----annual-data--------------------------------------------------------------
set.seed(93)
habitat_1 <- rast(nrows = 35, ncols = 40, xmin = 0, xmax = 800,
                  ymin = 0, ymax = 700, crs = "EPSG:26915")
xy <- xyFromCell(habitat_1, seq_len(ncell(habitat_1)))
values(habitat_1) <- as.integer(
  sin(xy[, 1] / 80) + cos(xy[, 2] / 105) > 0
)
names(habitat_1) <- "habitat"
habitat_2 <- habitat_1
habitat_2[xy[, 1] > 350 & xy[, 1] < 500 &
          xy[, 2] > 250 & xy[, 2] < 440] <- 0

observations <- data.frame(
  year = factor(rep(c("year1", "year2"), each = 30)),
  x = runif(60, 120, 680), y = runif(60, 120, 580)
)
rownames(observations) <- paste0("site_", seq_len(nrow(observations)))
points <- sf::st_as_sf(observations, coords = c("x", "y"), crs = 26915)

## ----grouped-preparation------------------------------------------------------
prepared <- kernel_prep_by_group(
  pts = points,
  raster_stacks = list(year1 = habitat_1, year2 = habitat_2),
  group = observations$year,
  max_D = 100,
  kernel = "gaussian",
  bin = TRUE,
  store_cell_data = FALSE,
  verbose = FALSE
)
head(prepared$kernel_dat)
table(prepared$raster_group)

## ----grouped-fit--------------------------------------------------------------
observations$survived <- rbinom(
  nrow(observations), 1,
  plogis(-0.4 + 0.9 * prepared$kernel_dat$habitat +
           0.3 * (observations$year == "year2"))
)
model_data <- cbind(observations, prepared$kernel_dat)
initial_model <- glm(survived ~ habitat + year,
                     family = binomial(), data = model_data,
                     na.action = na.fail)
fit <- multiScale_optim(initial_model, prepared,
                        n_cores = 1, verbose = FALSE)
fit$scale_est
coef(fit$opt_mod)
diagnostics(fit)$sample_size

## ----annual-projection, eval=FALSE--------------------------------------------
# surface_year1 <- kernel_scale.raster(habitat_1, multiScaleR = fit,
#                                      scale_center = TRUE, verbose = FALSE)
# surface_year2 <- kernel_scale.raster(habitat_2, multiScaleR = fit,
#                                      scale_center = TRUE, verbose = FALSE)
# predict_year <- function(surface, year_value) {
#   terra::predict(surface, fit$opt_mod, type = "response",
#                  fun = function(model, data, ...) {
#                    data$year <- factor(year_value,
#                                        levels = levels(observations$year))
#                    predict(model, newdata = data, ...)
#                  })
# }
# pred_year1 <- predict_year(surface_year1, "year1")
# pred_year2 <- predict_year(surface_year2, "year2")

## ----quick-reference, eval=FALSE----------------------------------------------
# # 1. Align the annual maps on one projected grid with matching layer names.
# maps <- list(year1 = habitat_1, year2 = habitat_2)
# 
# # 2. Assign each observation to its map; bins and scaling are pooled.
# prepared <- kernel_prep_by_group(points, maps, observations$year,
#                                   max_D = 100, bin = TRUE,
#                                   store_cell_data = FALSE)
# 
# # 3. Fit the response model appropriate for the sampling design.
# model_data <- cbind(observations, prepared$kernel_dat)
# initial_model <- glm(survived ~ habitat + year,
#                      family = binomial(), data = model_data)
# fit <- multiScale_optim(initial_model, prepared)
# summary(fit)               # shared habitat scale and model coefficients
# diagnostics(fit)$sample_size

