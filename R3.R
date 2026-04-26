# =========================================================
# SHAP
# =========================================================

## =========================
## Environment & Paths
## =========================
setwd("D:/")

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(rio)
  library(ranger)
  library(foreign)
  library(openxlsx)
  library(ggplot2)
  library(scales)
  library(grid)
  library(fastshap)
})

# Define output path
output_dir <- "D:/train_SHAP/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

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

## ---------- Standardization preparation ----------
predictors <- setdiff(colnames(train), "CAP")
predictors <- setdiff(predictors, c("lon","lat"))

is_num <- sapply(train[, ..predictors], function(x) is.numeric(x) || is.integer(x))
num_predictors     <- predictors[is_num]
non_num_predictors <- predictors[!is_num]

centers <- sapply(train[, ..num_predictors], function(x) mean(x, na.rm = TRUE))
scales_ <- sapply(train[, ..num_predictors], function(x) stats::sd(x, na.rm = TRUE))
scales_[is.na(scales_) | scales_ == 0] <- 1

std_apply <- function(DT, cols, centers, scales_) {
  for (nm in cols) {
    cn <- centers[[nm]]
    sc <- scales_[[nm]]
    suppressWarnings(DT[, (nm) := as.numeric(get(nm))])
    DT[, (nm) := (get(nm) - cn) / sc]
  }
  invisible(NULL)
}

# Keep raw data copy before scaling
train_raw <- copy(train)
test_raw  <- copy(test)

# Standardize training set
std_apply(train, num_predictors, centers, scales_)

# Standardize test set using training parameters
for (nm in num_predictors) suppressWarnings(test[, (nm) := as.numeric(get(nm))])
std_apply(test, num_predictors, centers, scales_)

## ---------- Train ranger random forest ----------
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
  na.action   = "na.omit"
)

# regularization.factor only if supported
if (is.finite(best_regfac)) fit_args$regularization.factor <- best_regfac

# Remove unsupported parameters
fit_args <- fit_args[names(fit_args) %in% names(formals(ranger::ranger))]

fit.ranger <- do.call(ranger::ranger, fit_args)

## ---------- Variable importance ----------
imp <- sort(fit.ranger$variable.importance, decreasing = TRUE)
print(imp)

op <- par(no.readonly = TRUE)
barplot(
  imp[1:min(30, length(imp))],
  las = 2, cex.names = 0.7,
  main = "Variable importance (ranger)"
)
par(op)

## ---------- Train/Test prediction ----------
train_pred <- predict(fit.ranger, data = train[, ..predictors])$predictions
test_pred  <- predict(fit.ranger, data = test[,  ..predictors])$predictions

# Save outputs
fwrite(as.data.table(train), file.path(output_dir, "train.csv"))
fwrite(data.table(train_pred = train_pred), file.path(output_dir, "train_pred.csv"))
fwrite(as.data.table(test), file.path(output_dir, "test.csv"))
fwrite(data.table(test_pred = test_pred), file.path(output_dir, "test_pred.csv"))
fwrite(
  data.table(variable = names(imp), importance = imp),
  file.path(output_dir, "ranger_importance.csv")
)

# Export scaling parameters
scaler_params <- data.table(
  feature = names(centers),
  center  = as.numeric(centers),
  scale   = as.numeric(scales_)
)
fwrite(scaler_params, file.path(output_dir, "scaler_params.csv"))

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

plot(
  test$CAP, test_pred,
  xlab = "Observed",
  ylab = "Predicted",
  main = "Observed vs Predicted"
)
abline(lm(test_pred ~ test$CAP), col = "red")

## =========================
## SHAP Calculation
## =========================
cat("\nStarting SHAP calculation...\n")

# -------- 1. Select samples for SHAP --------
set.seed(2026)
shap_n <- min(2000, nrow(train))
shap_idx <- sample(seq_len(nrow(train)), shap_n)

# Export X_shap in original units
X_shap <- copy(train_raw[shap_idx, ..predictors])
X_shap[, CAP_true := train_raw$CAP[shap_idx]]

# Keep standardized version for fastshap
X_shap_std <- copy(train[shap_idx, ..predictors])

# fastshap input should be data.frame
X_shap_df <- as.data.frame(X_shap_std)

# -------- 2. Prediction wrapper --------
pred_wrapper <- function(object, newdata) {
  predict(object, data = as.data.frame(newdata))$predictions
}

# -------- 3. Compute SHAP --------
set.seed(2026)
nsim_shap <- 50

shap_mat <- fastshap::explain(
  object       = fit.ranger,
  X            = as.data.frame(train[, ..predictors]),
  newdata      = X_shap_df,
  pred_wrapper = pred_wrapper,
  nsim         = nsim_shap,
  adjust       = TRUE
)

shap_dt <- as.data.table(shap_mat)
setnames(shap_dt, colnames(X_shap_df))

# -------- 4. Global importance mean(|SHAP|) --------
shap_imp <- data.table(
  variable = names(shap_dt),
  mean_abs_shap = sapply(shap_dt, function(x) mean(abs(x), na.rm = TRUE))
)
setorder(shap_imp, -mean_abs_shap)

## =========================
## Save SHAP results to Excel
## =========================
metrics_out <- data.table(
  target = "CAP",
  shap_n = shap_n,
  nsim_shap = nsim_shap,
  R2   = data_rf_R2,
  MAE  = data_rf_MAE,
  RMSE = data_rf_RMSE,
  MEC  = data_rf_MEC,
  ME   = data_rf_ME
)

xlsx_out <- file.path(output_dir, "TrainOut_CAP_SHAP.xlsx")

wb <- createWorkbook()

addWorksheet(wb, "metrics")
writeData(wb, "metrics", metrics_out)

addWorksheet(wb, "X_shap")
writeData(wb, "X_shap", X_shap)

addWorksheet(wb, "SHAP_values")
writeData(wb, "SHAP_values", shap_dt)

addWorksheet(wb, "SHAP_importance")
writeData(wb, "SHAP_importance", shap_imp)

saveWorkbook(wb, xlsx_out, overwrite = TRUE)
cat("SHAP Excel file saved:", xlsx_out, "\n")