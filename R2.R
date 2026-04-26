# =========================================================
# Building a Random Forest Model
# =========================================================

setwd("D:/")

suppressPackageStartupMessages({
  library(data.table)
  library(rio)
  library(ranger)
})

## =========================
## Read optimal hyperparameters (CSV)
## =========================
best_csv <- "D:/outputs_CAP_mlr3_WEIGHT_DENSITY_FULL_V2/best_params.csv"
if (!file.exists(best_csv)) stop("Optimal parameter file not found: ", best_csv)

best_dt <- fread(best_csv)
if (nrow(best_dt) < 1) stop("best_params.csv is empty: ", best_csv)
best <- as.list(best_dt[1])

# num.trees.code -> {300, 500, 800}
trees_map <- c(300L, 500L, 800L)
best_code <- suppressWarnings(as.integer(best[["regr.ranger_custom.num.trees.code"]]))
if (!is.finite(best_code)) best_code <- 2L
best_code <- max(1L, min(3L, best_code))
best_num_trees <- trees_map[best_code]

best_mtry    <- suppressWarnings(as.integer(best[["regr.ranger_custom.mtry"]]))
best_mns     <- suppressWarnings(as.integer(best[["regr.ranger_custom.min.node.size"]]))
best_depth   <- suppressWarnings(as.integer(best[["regr.ranger_custom.max.depth"]]))
best_sf      <- suppressWarnings(as.numeric(best[["regr.ranger_custom.sample.fraction"]]))
best_split   <- as.character(best[["regr.ranger_custom.splitrule"]])
best_replace <- as.logical(best[["regr.ranger_custom.replace"]])
best_respect <- as.character(best[["regr.ranger_custom.respect.unordered.factors"]])
best_seed    <- suppressWarnings(as.integer(best[["regr.ranger_custom.seed"]]))

best_regfac  <- suppressWarnings(as.numeric(best[["regr.ranger_custom.regularization.factor"]]))
if (!is.finite(best_regfac)) best_regfac <- NA_real_

# Fallback defaults
if (!is.finite(best_mtry))  best_mtry  <- 4L
if (!is.finite(best_mns))   best_mns   <- 7L
if (!is.finite(best_depth)) best_depth <- 10L
if (!is.finite(best_sf))    best_sf    <- 0.6
if (!is.finite(best_seed))  best_seed  <- 124754L
if (!best_respect %in% c("ignore","order","partition")) best_respect <- "order"
if (!best_split %in% c("variance","extratrees","maxstat","beta","poisson")) best_split <- "variance"
if (!isTRUE(best_replace) && !identical(best_replace, FALSE)) best_replace <- TRUE

cat("[Applying Optimal Hyperparameters]\n",
    "num.trees.code =", best_code, "-> num.trees =", best_num_trees, "\n",
    "mtry =", best_mtry, "\n",
    "min.node.size =", best_mns, "\n",
    "max.depth =", best_depth, "\n",
    "sample.fraction =", best_sf, "\n",
    "splitrule =", best_split, "\n",
    "replace =", best_replace, "\n",
    "respect.unordered.factors =", best_respect, "\n",
    "seed =", best_seed, "\n",
    "regularization.factor =", best_regfac, "\n")

## =========================
## Load modeling data and train model
## =========================
data <- import("APrate_modeling_data.csv")
setDT(data)
set.seed(124754)

# 7:3 split
train_idx <- sample(nrow(data), round(0.7 * nrow(data)))
train <- copy(data[train_idx, ])
test  <- copy(data[-train_idx, ])

## ---------- Standardization (numeric predictors only; based on training set) ----------
predictors <- setdiff(colnames(train), "CAP")
predictors <- setdiff(predictors, c("lon","lat"))

is_num <- sapply(train[, ..predictors], function(x) is.numeric(x) || is.integer(x))
num_predictors     <- predictors[is_num]
non_num_predictors <- predictors[!is_num]

centers <- sapply(train[, ..num_predictors], function(x) mean(x, na.rm = TRUE))
scales  <- sapply(train[, ..num_predictors], function(x) stats::sd(x, na.rm = TRUE))
scales[is.na(scales) | scales == 0] <- 1

std_apply <- function(DT, cols, centers, scales) {
  for (nm in cols) {
    cn <- centers[[nm]]
    sc <- scales[[nm]]
    suppressWarnings(DT[, (nm) := as.numeric(get(nm))])
    DT[, (nm) := (get(nm) - cn) / sc]
  }
  invisible(NULL)
}

# Standardize training set
std_apply(train, num_predictors, centers, scales)

# Standardize test set using training parameters
for (nm in num_predictors) suppressWarnings(test[, (nm) := as.numeric(get(nm))])
std_apply(test, num_predictors, centers, scales)

## ---------- Train ranger random forest (quantreg=TRUE) ----------
fit_args <- list(
  formula     = CAP ~ .,
  data        = train[, c("CAP", predictors), with = FALSE],
  num.trees   = best_num_trees,
  mtry        = best_mtry,
  min.node.size = best_mns,
  max.depth   = best_depth,
  sample.fraction = best_sf,
  splitrule   = best_split,
  replace     = best_replace,
  seed        = best_seed,
  respect.unordered.factors = best_respect,
  importance  = "permutation",
  na.action   = "na.omit",
  quantreg    = TRUE
)

if (is.finite(best_regfac)) fit_args$regularization.factor <- best_regfac
fit_args <- fit_args[names(fit_args) %in% names(formals(ranger::ranger))]

fit.ranger <- do.call(ranger::ranger, fit_args)

cat("=== Variable Importance (Permutation) ===\n")
imp <- sort(fit.ranger$variable.importance, decreasing = TRUE)

for (i in 1:length(imp)) {
  cat(names(imp)[i], ":", imp[i], "\n")
}

op <- par(no.readonly = TRUE)
barplot(
  imp[order(imp)[1:min(30, length(imp))]],
  horiz = TRUE,
  las = 1,
  cex.names = 0.7,
  main = "Variable importance (ranger)",
  xlim = c(0, max(imp[order(imp)[1:min(30, length(imp))]]) * 1.2),
  xaxt = "n"
)
axis(1, at = seq(0, max(imp[order(imp)[1:min(30, length(imp))]]) * 1.2, by = 0.05))
par(op)

## ---------- Train/Test prediction and evaluation ----------
train_pred <- predict(fit.ranger, data = train[, ..predictors])$predictions
test_pred  <- predict(fit.ranger, data = test[,  ..predictors])$predictions

fwrite(as.data.table(train),      "train_std.csv")
fwrite(as.data.table(train_pred), "train_pred.csv")
fwrite(as.data.table(test),       "test_std.csv")
fwrite(as.data.table(test_pred),  "test_pred.csv")

## ---------- Evaluation metrics ----------
data_rf_R2   <- cor(test$CAP, test_pred, use = "complete.obs")^2
data_rf_MAE  <- mean(abs(test$CAP - test_pred), na.rm = TRUE)
data_rf_RMSE <- sqrt(mean((test$CAP - test_pred)^2, na.rm = TRUE))

# Model Efficiency Coefficient (MEC)
mean_observed <- mean(test$CAP, na.rm = TRUE)
mec_numerator <- sum((test$CAP - test_pred)^2, na.rm = TRUE)
mec_denominator <- sum((test$CAP - mean_observed)^2, na.rm = TRUE)
data_rf_MEC <- 1 - (mec_numerator / mec_denominator)

# Mean Error (ME)
data_rf_ME <- mean(test_pred - test$CAP, na.rm = TRUE)

cat(sprintf("R2=%.6f  MAE=%.6f  RMSE=%.6f\n", data_rf_R2, data_rf_MAE, data_rf_RMSE))
cat(sprintf("MEC=%.6f  ME=%.6f\n", data_rf_MEC, data_rf_ME))

plot(test$CAP, test_pred,
     xlab = "Observed",
     ylab = "Predicted",
     main = "Observed vs Predicted")
abline(lm(test_pred ~ test$CAP), col = "red")
