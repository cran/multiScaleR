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
library(multiScaleR)
library(terra)

# The step-selection backend fits its models with the survival package, which
# multiScaleR lists under Suggests. Code chunks run only when it is installed.
has_survival <- requireNamespace("survival", quietly = TRUE)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6.5,
  fig.height = 4.25,
  warning  = FALSE,
  message  = FALSE,
  eval     = has_survival
)

## ----simulate-landscape, fig.cap = "***The simulated forest layer. Forest is 1 and open ground is 0.***"----
set.seed(321)

# A 200 x 200 grid of 20 m cells in a projected CRS, so map units are meters
noise <- rast(nrows = 200, ncols = 200, xmin = 0, xmax = 4000,
              ymin = 0, ymax = 4000, crs = "EPSG:32617")
values(noise) <- rnorm(ncell(noise))

# Smooth the noise and cut it at zero to get patches of forest
forest <- focal(noise, w = focalMat(noise, 100, "Gauss"), fun = "sum",
                na.rm = TRUE) > 0
forest <- forest * 1                   # logical to numeric 0/1
names(forest) <- "forest"
plot(forest)

## ----simulate-steps-----------------------------------------------------------
# The landscape as the animal perceives it: forest weighted by the true kernel
true_sigma <- 150
perceived  <- focal(forest, w = focalMat(forest, true_sigma, "Gauss"),
                    fun = "sum", na.rm = TRUE)

n_sets <- 400                          # choice sets (one per observed step)
n_alt  <- 5                            # candidate end points per choice set

# Each choice set starts at a random location away from the map edge.
# Candidates leave it in a random direction with a random step length (m).
start_x <- rep(runif(n_sets, 1000, 3000), each = n_alt)
start_y <- rep(runif(n_sets, 1000, 3000), each = n_alt)
step_m  <- rgamma(n_sets * n_alt, shape = 2, scale = 100)
angle   <- runif(n_sets * n_alt, 0, 2 * pi)

steps <- data.frame(
  step_id = rep(seq_len(n_sets), each = n_alt),   # choice-set ID
  x       = start_x + step_m * cos(angle),
  y       = start_y + step_m * sin(angle),
  sl_km   = step_m / 1000                         # step length in km
)

# Perceived forest at each end point, standardized across all end points
f <- extract(perceived, as.matrix(steps[, c("x", "y")]))[, 1]
z <- (f - mean(f)) / sd(f)

# True selection: +1 per standard deviation of forest, -3 per km of step.
# The animal takes exactly one candidate per choice set.
score <- exp(1 * z - 3 * steps$sl_km)
steps$case_ <- 0
for (s in seq_len(n_sets)) {
  rows <- which(steps$step_id == s)
  steps$case_[sample(rows, 1, prob = score[rows])] <- 1
}
head(steps, 5)

## ----fit----------------------------------------------------------------------
library(survival)

# Step end points as a terra SpatVector, in the row order of `steps`
pts <- vect(steps, geom = c("x", "y"), crs = "EPSG:32617")

# max_D is the largest distance (m) considered around a point
inputs <- kernel_prep(pts = pts, raster_stack = forest, max_D = 600,
                      kernel = "gaussian", bin = TRUE, verbose = FALSE)

# Add the standardized, kernel-weighted forest covariate to the step table
df <- data.frame(steps, inputs$kernel_dat)

# Conditional logistic regression: compare steps only within a choice set
mod <- clogit(case_ ~ forest + sl_km + strata(step_id),
              data = df, model = TRUE)

# Estimate sigma and the selection coefficients in one likelihood
opt <- multiScale_optim(fitted_mod = mod, kernel_inputs = inputs,
                        backend = "conditional_logit",
                        n_cores = 1)            # this backend runs serially

## ----read-scale---------------------------------------------------------------
opt$scale_est               # sigma in meters. This backend reports no SE
kernel_dist(opt)            # radius holding 90% of the kernel weight (m)

b <- coef(opt$opt_mod)      # selection coefficients on the log scale
round(cbind(estimate = b, confint(opt$opt_mod)), 2)
round(exp(b), 2)            # relative selection strength
opt$scl_params$sd           # one standard deviation of kernel-weighted forest

## ----read-values, include = FALSE---------------------------------------------
sig  <- opt$scale_est["forest", "Mean"]
d90  <- kernel_dist(opt)["forest", "Mean"]
rss  <- exp(b)
sd_f <- opt$scl_params$sd[["forest"]]

## ----plot-kernel, fig.cap = "***Kernel weight against distance from a step's end point.***", fig.width = 4.5, fig.height = 4----
plot(opt)

## ----plot-effects, fig.cap = "***Relative selection strength across each covariate.***", fig.show = "hold", out.width = "48%", fig.width = 4, fig.height = 3.6----
# link = TRUE exponentiates the curves, giving ratios in place of log ratios
plot_marginal_effects(opt, link = TRUE)

## ----predict------------------------------------------------------------------
p_set  <- predict(opt)                   # probability within each choice set
p_risk <- predict(opt, type = "risk")    # relative selection score, exp(link)

# The first choice set, with standardized forest at the estimated sigma
first <- data.frame(steps[1:5, c("step_id", "case_", "sl_km")],
                    forest      = round(model.frame(opt$opt_mod)$forest[1:5], 2),
                    probability = round(p_set[1:5], 3),
                    score       = round(p_risk[1:5], 2))
first

## ----cross-validate-----------------------------------------------------------
# Four folds. All rows of a choice set stay together
folds <- grouped_cv_folds(df$step_id, k = 4)

# Reference model: step length only, no landscape covariate
cv_move <- cross_validate_clogit(case_ ~ sl_km + strata(step_id),
                                 data = df, folds = folds)

# Forest model: sigma is estimated again within every training fold
cv_forest <- cross_validate_clogit(case_ ~ forest + sl_km + strata(step_id),
                                   data = df, folds = folds,
                                   kernel_inputs = inputs)

cv_tab(list(movement = cv_move, forest = cv_forest),
       metric = "neg_log_score")[, c("Modnames", "Score", "Delta", "N_scored")]

# Sigma (m) estimated in each training fold
round(sapply(cv_forest$training, function(fold) fold$scale_est$Mean))

## ----quick-reference, eval = FALSE--------------------------------------------
# library(multiScaleR); library(terra); library(survival)
# 
# # 1. Prepare kernel inputs at every step end point (used and available)
# pts    <- vect(steps, geom = c("x", "y"), crs = "EPSG:32617")
# inputs <- kernel_prep(pts = pts, raster_stack = forest, max_D = 600,
#                       kernel = "gaussian", bin = TRUE)
# df     <- data.frame(steps, inputs$kernel_dat)
# 
# # 2. Fit the starting model, then estimate the scale of effect
# mod <- clogit(case_ ~ forest + sl_km + strata(step_id), data = df, model = TRUE)
# opt <- multiScale_optim(fitted_mod = mod, kernel_inputs = inputs,
#                         backend = "conditional_logit", n_cores = 1)
# 
# # 3. Read, plot, and predict
# opt$scale_est                            # sigma (m), no SE from this backend
# kernel_dist(opt)                         # radius holding 90% of kernel weight
# exp(coef(opt$opt_mod))                   # relative selection strength
# plot(opt)                                # kernel weight against distance
# plot_marginal_effects(opt, link = TRUE)  # relative selection curves
# predict(opt)                             # probability within each choice set
# predict(opt, type = "risk")              # relative selection score
# 
# # 4. Score on held-out choice sets
# folds     <- grouped_cv_folds(df$step_id, k = 4)
# cv_move   <- cross_validate_clogit(case_ ~ sl_km + strata(step_id),
#                                    data = df, folds = folds)
# cv_forest <- cross_validate_clogit(case_ ~ forest + sl_km + strata(step_id),
#                                    data = df, folds = folds,
#                                    kernel_inputs = inputs)
# cv_tab(list(movement = cv_move, forest = cv_forest), metric = "neg_log_score")

