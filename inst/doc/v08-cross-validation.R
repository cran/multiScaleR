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

knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6.5,
                      fig.height = 4.25, warning = FALSE, message = FALSE)
options(scipen = 10)   # print large values without scientific notation

## ----prepare------------------------------------------------------------------
data("landscape_counts")   # counts and one site-level covariate
data("surv_pts")           # survey locations
pts  <- vect(surv_pts)     # convert to a terra SpatVector
land <- rast(system.file("extdata", "landscape.tif", package = "multiScaleR"))

# bin = TRUE stores the distance-binned summaries used by the default
# fitting backend of cross_validate()
inputs <- kernel_prep(pts = pts, raster_stack = land, max_D = 1700,
                      kernel = "gaussian", bin = TRUE, verbose = FALSE)

# Response data and kernel-weighted covariates, one row per survey point
dat <- data.frame(landscape_counts, inputs$kernel_dat)

## ----spatial-folds------------------------------------------------------------
# 1,000 m blocks assigned to three folds. Coordinates come from terra.
folds <- spatial_cv_folds(crds(pts), block_size = 1000, k = 3, crs = 32617)
folds

## ----fold-map, fig.height = 5, fig.cap = "***Survey points over `land1`, colored by the fold that holds them out (black, fold 1; blue, fold 2; red, fold 3). Dashed lines mark the 1,000 m blocks.***"----
# One row per survey point: coordinates, block, and assessment fold
where <- attr(folds, "fold_info")$assignment
fold_col <- c(fold1 = "black", fold2 = "#2166AC", fold3 = "#B2182B")

plot(land$land1, col = gray(c(0.95, 0.75)), legend = FALSE)
abline(v = 1000 * 1:3, h = 1000 * 1:3, lty = 2)           # block edges
points(where$x, where$y, pch = 19, col = fold_col[where$fold])

## ----buffer-folds-------------------------------------------------------------
# Drop training blocks that have a point within 320 m of an assessment point
folds_buffer <- spatial_cv_folds(crds(pts), block_size = 1000, k = 3,
                                 buffer = 320, crs = 32617)
folds_buffer

## ----grouped-folds------------------------------------------------------------
# Twelve records from six sites, two visits each
site_id <- rep(c("A", "B", "C", "D", "E", "F"), each = 2)
grouped <- grouped_cv_folds(site_id, k = 3, ids = paste0("visit", 1:12))
grouped$fold1$assessment    # both visits to a site are held out together

## ----random-folds-------------------------------------------------------------
set.seed(1)
ids  <- rownames(dat)
part <- sample(rep(1:3, length.out = nrow(dat)))   # random fold label per row
folds_random <- lapply(1:3, function(k) {
  list(analysis = ids[part != k], assessment = ids[part == k])
})
names(folds_random) <- paste0("fold", 1:3)         # fold names are required

## ----cross-validate-----------------------------------------------------------
cv <- cross_validate(counts ~ site + land1 + land2,
                     data          = dat,
                     kernel_inputs = inputs,
                     folds         = folds,
                     family        = "poisson",
                     keep_selection = TRUE)   # needed later for select_cv()
cv

## ----fold-scales--------------------------------------------------------------
# Sigma (m) estimated from the analysis rows of each fold
round(t(sapply(cv$training, function(fold) {
  setNames(fold$scale_est$Mean, rownames(fold$scale_est))
})))

## ----reference-models---------------------------------------------------------
# Intercept only: predicts the training mean count everywhere
null_cv <- cross_validate_glm(counts ~ 1, dat, folds, family = "poisson",
                              keep_selection = TRUE)
# Site covariate only: no landscape information
site_cv <- cross_validate_glm(counts ~ site, dat, folds, family = "poisson",
                              keep_selection = TRUE)
candidates <- list(landscape = cv, site_only = site_cv, null = null_cv)

## ----cv-tab-------------------------------------------------------------------
tab_rmse <- cv_tab(candidates, metric = "rmse")
tab_nls  <- cv_tab(candidates, metric = "neg_log_score")

tab_rmse[, c("Modnames", "Rank", "Score", "Delta", "N_scored", "Warnings")]
tab_nls[, c("Modnames", "Rank", "Score", "Delta", "N_scored", "Warnings")]

## ----cv-tab-fold--------------------------------------------------------------
by_fold <- cv_tab(candidates, metric = "rmse", by_fold = TRUE)
by_fold[, c("Fold", "Modnames", "Rank", "Score", "Delta")]

## ----complete-----------------------------------------------------------------
cv$complete        # TRUE: every fold was fitted and scored
cv$fold_results    # rows used for fitting and assessment, and the fold status

## ----design-comparison--------------------------------------------------------
# The same landscape model, scored on random rows and on the buffered blocks
cv_random <- cross_validate(counts ~ site + land1 + land2, dat, inputs,
                            folds_random, family = "poisson")
cv_buffer <- cross_validate(counts ~ site + land1 + land2, dat, inputs,
                            folds_buffer, family = "poisson")

rmse_of <- function(x) x$summary$estimate[x$summary$metric == "rmse"]
data.frame(design   = c("Random rows", "Spatial blocks", "Blocks, 320 m buffer"),
           rmse     = round(c(rmse_of(cv_random), rmse_of(cv), rmse_of(cv_buffer)), 2),
           complete = c(cv_random$complete, cv$complete, cv_buffer$complete))

## ----spline-------------------------------------------------------------------
library(splines)

# A natural spline in land1 that chooses its knots from the data
df_cv <- cross_validate(counts ~ site + ns(land1, df = 3) + land2,
                        dat, inputs, folds, family = "poisson")
df_cv$complete                # FALSE: no fold was scored
df_cv$errors$fold1$message    # the reason, recorded for each fold

# The message names the fix
curve_cv <- cross_validate(counts ~ site + ns(land1, df = 3) + land2,
                           dat, inputs, folds, family = "poisson",
                           backend = "refit")
curve_cv$complete
cv_tab(list(linear = cv, curved = curve_cv),
       metric = "rmse")[, c("Modnames", "Rank", "Score", "Delta")]

## ----nb2, eval = requireNamespace("MASS", quietly = TRUE)---------------------
# Negative binomial versions of the site-only and landscape models
site_nb <- cross_validate_glm(counts ~ site, dat, folds, family = "nbinom2")
land_nb <- cross_validate(counts ~ site + land1 + land2, dat, inputs, folds,
                          family = "nbinom2")

count_models <- list(site_poisson = site_cv, site_nb = site_nb,
                     landscape_poisson = cv, landscape_nb = land_nb)
cv_tab(count_models,
       metric = "neg_log_score")[, c("Modnames", "Rank", "Score", "Delta", "Warnings")]

# Theta estimated in each fold
round(rbind(site_nb      = sapply(site_nb$distribution, function(fold) fold$theta),
            landscape_nb = sapply(land_nb$distribution, function(fold) fold$theta)), 1)

## ----select-cv----------------------------------------------------------------
selected <- select_cv(candidates, criterion = "AICc")
selected$selection$chosen          # model chosen in each fold
selected$selection$tables$fold1    # training AICc table for the first fold

## ----mixed-cv, eval = requireNamespace("glmmTMB", quietly = TRUE)-------------
set.seed(10)
site        <- factor(rep(sprintf("site%02d", 1:18), each = 8))   # 18 sites
habitat     <- rnorm(length(site))                       # standardized habitat
site_effect <- rnorm(nlevels(site), sd = 0.7)[site]      # one value per site
surveys <- data.frame(site, habitat, visit = rep(1:8, nlevels(site)))
surveys$count <- rpois(nrow(surveys), exp(0.9 + 0.4 * habitat + site_effect))

# New sites: all eight visits to a site are held out together
new_sites <- grouped_cv_folds(surveys$site, k = 3, ids = rownames(surveys))

habitat_new   <- cross_validate_mixed(count ~ habitat + (1 | site), surveys,
                                      new_sites, family = "poisson",
                                      prediction = "new_group")
intercept_new <- cross_validate_mixed(count ~ 1 + (1 | site), surveys,
                                      new_sites, family = "poisson",
                                      prediction = "new_group")
mixed_tab <- cv_tab(list(habitat = habitat_new, intercept = intercept_new),
                    metric = "neg_log_score")
mixed_tab[, c("Modnames", "Rank", "Score", "Delta", "Unit")]

## ----mixed-cv-later, eval = requireNamespace("glmmTMB", quietly = TRUE)-------
# Later visits: train on visits 1 to 6, predict visits 7 and 8 at every site
ids <- rownames(surveys)
later_visits <- list(later = list(analysis   = ids[surveys$visit <= 6],
                                  assessment = ids[surveys$visit > 6]))
habitat_later <- cross_validate_mixed(count ~ habitat + (1 | site), surveys,
                                      later_visits, family = "poisson",
                                      prediction = "existing_group")
habitat_later$summary[, c("metric", "estimate", "n_assessment")]

## ----clogit, eval = requireNamespace("survival", quietly = TRUE)--------------
library(survival)   # provides strata()

i <- 1:72
choices <- data.frame(
  step    = (i - 1) %/% 3 + 1,          # choice set label, 1 to 24
  used    = as.numeric(i %% 3 == 1),    # the first step of each set was used
  habitat = sin(i * 1.7)                # covariate unrelated to use
)
choice_folds <- grouped_cv_folds(choices$step, k = 3)   # whole sets held out

choice_models <- list(
  null    = cross_validate_clogit(used ~ strata(step), choices, choice_folds),
  habitat = cross_validate_clogit(used ~ habitat + strata(step), choices,
                                  choice_folds))
cv_tab(choice_models,
       metric = "neg_log_score")[, c("Modnames", "Rank", "Score", "Delta", "Unit")]

## ----quick-reference, eval = FALSE--------------------------------------------
# library(multiScaleR)
# library(terra)
# 
# # 1. Prepare the data
# data("landscape_counts"); data("surv_pts")
# pts  <- vect(surv_pts)
# land <- rast(system.file("extdata", "landscape.tif", package = "multiScaleR"))
# inputs <- kernel_prep(pts = pts, raster_stack = land, max_D = 1700,
#                       kernel = "gaussian", bin = TRUE, verbose = FALSE)
# dat <- data.frame(landscape_counts, inputs$kernel_dat)
# 
# # 2. Build folds that match the prediction question
# folds <- spatial_cv_folds(crds(pts), block_size = 1000, k = 3, crs = 32617)
# folds_buffer <- spatial_cv_folds(crds(pts), block_size = 1000, k = 3,
#                                  buffer = 320, crs = 32617)
# site_id <- rep(c("A", "B", "C", "D", "E", "F"), each = 2)   # group per row
# grouped <- grouped_cv_folds(site_id, k = 3, ids = paste0("visit", 1:12))
# 
# # 3. Cross-validate: scales and coefficients are re-estimated in every fold
# cv <- cross_validate(counts ~ site + land1 + land2, data = dat,
#                      kernel_inputs = inputs, folds = folds,
#                      family = "poisson", keep_selection = TRUE)
# cv                      # pooled scores, lower is better
# cv$complete             # TRUE: every fold was fitted and scored
# cv$training             # scales and coefficients estimated in each fold
# 
# # 4. Reference models on the same folds
# null_cv <- cross_validate_glm(counts ~ 1, dat, folds, family = "poisson",
#                               keep_selection = TRUE)
# site_cv <- cross_validate_glm(counts ~ site, dat, folds, family = "poisson",
#                               keep_selection = TRUE)
# candidates <- list(landscape = cv, site_only = site_cv, null = null_cv)
# 
# # 5. Compare
# cv_tab(candidates, metric = "rmse")                   # pooled, in counts
# cv_tab(candidates, metric = "neg_log_score")          # whole distribution
# cv_tab(candidates, metric = "rmse", by_fold = TRUE)   # is the ranking consistent?
# 
# # 6. Score an AICc selection rule
# selected <- select_cv(candidates, criterion = "AICc")
# selected$selection$chosen                             # model chosen per fold
# 
# # 7. Repeated surveys: choose the prediction target, then compare candidates
# new_sites <- grouped_cv_folds(surveys$site, k = 3, ids = rownames(surveys))
# habitat_new <- cross_validate_mixed(count ~ habitat + (1 | site), surveys,
#                                     new_sites, family = "poisson",
#                                     prediction = "new_group")   # unseen sites
# intercept_new <- cross_validate_mixed(count ~ 1 + (1 | site), surveys,
#                                       new_sites, family = "poisson",
#                                       prediction = "new_group")
# cv_tab(list(habitat = habitat_new, intercept = intercept_new),
#        metric = "neg_log_score")                    # joint score per survey

