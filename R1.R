# =========================================================
# cross-validation + Parameter tuning (Bayesian)
# + Sampling density correction
# =========================================================

# =========================
# 0) Environment & Paths
# =========================
setwd("D:/")
model_csv <- "APrate_modeling_data.csv"
out_dir <- "outputs_CAP_mlr3_WEIGHT_DENSITY_FULL_V2"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

set.seed(124754)

suppressPackageStartupMessages({
  library(rio)
  
  library(mlr3)
  library(mlr3pipelines)
  library(mlr3tuning)
  library(paradox)
  library(bbotk)
  
  library(ranger)
  library(R6)
  
  library(spatstat.geom)
  library(spatstat.explore)
})

# =========================
# 0.1) Tuner: prefer mbo, otherwise random_search
# =========================
get_tuner_compat <- function() {
  keys <- tryCatch(mlr3tuning::mlr_tuners$keys(), error = function(e) character(0))
  if ("mbo" %in% keys) return(list(tuner = tnr("mbo"), mode = "mbo"))
  
  if (requireNamespace("mlr3mbo", quietly = TRUE)) {
    ok <- tryCatch({
      suppressPackageStartupMessages(library(mlr3mbo))
      TRUE
    }, error = function(e) FALSE)
    if (ok) {
      keys2 <- tryCatch(mlr3tuning::mlr_tuners$keys(), error = function(e) character(0))
      if ("mbo" %in% keys2) return(list(tuner = tnr("mbo"), mode = "mbo"))
    }
  }
  list(tuner = tnr("random_search"), mode = "random_search")
}

# =========================
# 0.2) Terminator: compatible with n_evals / n_evaluations
# =========================
make_terminator_evals <- function(n_eval) {
  tryCatch(
    trm("evals", n_evals = n_eval),
    error = function(e1) {
      tryCatch(
        trm("evals", n_evaluations = n_eval),
        error = function(e2) stop(
          "Failed to create trm('evals') terminator.\n",
          "e1: ", conditionMessage(e1), "\n",
          "e2: ", conditionMessage(e2)
        )
      )
    }
  )
}

# =========================
# A) Sampling density -> weights
#    (bw.CvL + density.ppp(at='points'))
# =========================
calc_density_weights <- function(df, lon_col = "lon", lat_col = "lat", expand = 0.02) {
  if (!(lon_col %in% names(df) && lat_col %in% names(df))) {
    stop("calc_density_weights: Missing lon/lat columns in data.")
  }
  
  lon <- df[[lon_col]]
  lat <- df[[lat_col]]
  ok <- is.finite(lon) & is.finite(lat)
  
  w <- rep(NA_real_, nrow(df))
  sigma_used <- NA_real_
  
  if (sum(ok) < 10) {
    warning("Too few valid lon/lat points; density cannot be estimated. Weight set to 1.")
    w[ok] <- 1
    return(list(weight = w, sigma = sigma_used))
  }
  
  lon2 <- lon[ok]
  lat2 <- lat[ok]
  
  xr <- range(lon2); yr <- range(lat2)
  dx <- diff(xr); dy <- diff(yr)
  win <- owin(
    xrange = xr + c(-1, 1) * dx * expand,
    yrange = yr + c(-1, 1) * dy * expand
  )
  
  pp <- ppp(x = lon2, y = lat2, window = win)
  sig <- bw.CvL(pp)
  sigma_used <- as.numeric(sig)
  
  lam <- density.ppp(pp, sigma = sig, at = "points", edge = TRUE)
  lam <- as.numeric(lam)
  lam[!is.finite(lam) | lam <= 0] <- NA_real_
  
  ww <- 1 / lam
  ww[!is.finite(ww)] <- NA_real_
  
  ww <- ww / mean(ww, na.rm = TRUE)  # normalize mean(w)=1
  
  w[ok] <- ww
  list(weight = w, sigma = sigma_used)
}

# Manual weighted RMSE (for verification)
rmse_weighted_manual <- function(obs, pred, w) {
  ok <- is.finite(obs) & is.finite(pred) & is.finite(w) & w > 0
  if (sum(ok) < 3) return(NA_real_)
  sqrt(sum(w[ok] * (obs[ok] - pred[ok])^2) / sum(w[ok]))
}

# =========================
# 1) Load data & split train/test
# =========================
data <- import(model_csv)
stopifnot("CAP" %in% names(data))

if (!all(c("lon", "lat") %in% names(data))) {
  stop("lon and lat columns are required for sampling density correction, but were not found.")
}

idx <- sample.int(nrow(data), round(0.7 * nrow(data)))
train_dt <- data[idx, , drop = FALSE]
test_dt  <- data[-idx, , drop = FALSE]

# Calculate weights (density estimated within each dataset)
w_train_obj <- calc_density_weights(train_dt, "lon", "lat", expand = 0.02)
w_test_obj  <- calc_density_weights(test_dt,  "lon", "lat", expand = 0.02)

train_dt$weight <- w_train_obj$weight
test_dt$weight  <- w_test_obj$weight

cat("bw.CvL sigma(train) =", w_train_obj$sigma, "\n")
cat("bw.CvL sigma(test)  =", w_test_obj$sigma, "\n")
cat("train weight NA:", sum(!is.finite(train_dt$weight)), "/", nrow(train_dt), "\n")
cat("test  weight NA:", sum(!is.finite(test_dt$weight)),  "/", nrow(test_dt),  "\n")

# lon/lat excluded from modeling
train_dt <- train_dt[, setdiff(names(train_dt), c("lon", "lat")), drop = FALSE]
test_dt  <- test_dt[,  setdiff(names(test_dt),  c("lon", "lat")), drop = FALSE]

task_train <- TaskRegr$new(id = "CAP_train", backend = train_dt, target = "CAP")
task_test  <- TaskRegr$new(id = "CAP_test",  backend = test_dt,  target = "CAP")

# Weight column roles
task_train$set_col_roles("weight", roles = c("weights_learner", "weights_measure"))
task_test$set_col_roles("weight",  roles = c("weights_learner", "weights_measure"))

p <- length(task_train$feature_names)
cat("Train:", task_train$nrow, "Test:", task_test$nrow, "p:", p, "\n")

# =========================
# 2) Custom Ranger learner
# =========================
LearnerRegrRangerCustom <- R6Class(
  "LearnerRegrRangerCustom",
  inherit = LearnerRegr,
  public = list(
    initialize = function() {
      
      ps_learner <- ps(
        num.trees.code = p_int(1L, 3L, tags = "train"),
        mtry = p_int(1L, p, tags = "train"),
        
        min.node.size = p_int(5L, 10L, tags = "train"),
        max.depth = p_int(5L, 20L, tags = "train"),
        sample.fraction = p_dbl(0.4, 0.8, tags = "train"),
        
        regularization.factor = p_dbl(0.01, 1.0, tags = "train"),
        splitrule = p_fct(levels = c("variance", "maxstat"), tags = "train"),
        replace = p_lgl(tags = "train"),
        respect.unordered.factors = p_fct(levels = c("ignore", "order", "partition"), tags = "train"),
        seed = p_int(1L, .Machine$integer.max, tags = "train")
      )
      
      super$initialize(
        id = "regr.ranger_custom",
        feature_types = c("logical","integer","numeric","factor","ordered","character"),
        predict_types = "response",
        param_set = ps_learner,
        packages = "ranger",
        properties = c("missings", "weights")
      )
      
      self$param_set$values <- list(
        num.trees.code = 2L,
        replace = TRUE,
        respect.unordered.factors = "order",
        seed = 99999L,
        splitrule = "variance",
        regularization.factor = 0.5,
        
        min.node.size = 7L,
        max.depth = 10L,
        sample.fraction = 0.6
      )
    }
  ),
  private = list(
    .train = function(task) {
      dt <- task$data()
      yname <- task$target_names
      feat_names <- task$feature_names
      pp <- length(feat_names)
      
      pv <- self$param_set$values
      
      trees_map <- c(300L, 500L, 800L)
      code <- as.integer(pv$num.trees.code)
      code <- max(1L, min(3L, code))
      num.trees <- trees_map[code]
      
      mtry <- as.integer(pv$mtry)
      mtry <- max(1L, min(pp, mtry))
      
      mns <- suppressWarnings(as.integer(pv$min.node.size))
      if (length(mns) != 1L || !is.finite(mns)) mns <- 7L
      mns <- max(5L, min(10L, mns))
      
      sf <- suppressWarnings(as.numeric(pv$sample.fraction))
      if (length(sf) != 1L || !is.finite(sf)) sf <- 0.6
      sf <- max(0.4, min(0.8, sf))
      
      md <- suppressWarnings(as.integer(pv$max.depth))
      if (length(md) != 1L || !is.finite(md)) md <- 10L
      md <- max(5L, min(20L, md))
      
      w <- task$data(cols = "weight")$weight
      ok_w <- is.finite(w) & (w > 0)
      
      dt2 <- dt[ok_w, , drop = FALSE]
      w2  <- w[ok_w]
      
      if (nrow(dt2) < 20) stop("Too few valid weighted samples: ", nrow(dt2))
      
      args <- list(
        dependent.variable.name = yname,
        data = dt2,
        num.trees = num.trees,
        mtry = mtry,
        
        min.node.size = mns,
        max.depth = md,
        sample.fraction = sf,
        
        splitrule = as.character(pv$splitrule),
        replace = isTRUE(pv$replace),
        respect.unordered.factors = as.character(pv$respect.unordered.factors),
        seed = as.integer(pv$seed),
        na.action = "na.omit",
        
        case.weights = w2
      )
      
      if (!is.null(pv$regularization.factor) && pv$splitrule != "maxstat") {
        rf <- suppressWarnings(as.numeric(pv$regularization.factor))
        if (length(rf) == 1L && is.finite(rf)) {
          rf <- max(0.0, min(1.0, rf))
          args$regularization.factor <- rf
        }
      }
      
      self$model <- do.call(ranger::ranger, args)
    },
    
    .predict = function(task) {
      newdata <- task$data(cols = task$feature_names)
      pred <- predict(self$model, data = newdata)$predictions
      PredictionRegr$new(task = task, response = pred)
    }
  )
)

lrn_ranger <- LearnerRegrRangerCustom$new()

# =========================
# 3) Pipeline: scale + learner
# =========================
graph <- po("scale") %>>% po("learner", lrn_ranger)
glrn  <- GraphLearner$new(graph)

# =========================
# 4) Auto-detect prefix + search space
# =========================
ids <- glrn$param_set$ids()
nt_id <- ids[grepl("num\\.trees\\.code$", ids)]
if (length(nt_id) != 1L) {
  stop("Cannot uniquely identify num.trees.code. Current ids:\n", paste(ids, collapse = ", "))
}
prefix <- sub("num\\.trees\\.code$", "", nt_id)
pn <- function(x) paste0(prefix, x)

mtry_center <- max(1L, round(p / 3))
mtry_low    <- max(1L, floor(mtry_center * 0.5))
mtry_up     <- min(p,  ceiling(mtry_center * 2.0))

par_list <- list(
  p_int(lower = 1L, upper = 3L),
  p_int(lower = mtry_low, upper = mtry_up),
  p_int(lower = 5L, upper = 10L),
  p_int(lower = 5L, upper = 20L),
  p_dbl(lower = 0.4, upper = 0.8),
  p_dbl(lower = 0.01, upper = 1.0),
  p_fct(levels = c("variance"))
)

id_vec <- c(
  pn("num.trees.code"),
  pn("mtry"),
  pn("min.node.size"),
  pn("max.depth"),
  pn("sample.fraction"),
  pn("regularization.factor"),
  pn("splitrule")
)

for (i in seq_along(par_list)) par_list[[i]]$id <- id_vec[i]
names(par_list) <- id_vec
search_space <- ParamSet$new(par_list)

# =========================
# 5) 5-fold CV + tuning
# =========================
resampling <- rsmp("cv", folds = 5)

measure <- msr("regr.rmse")

tuner_obj <- get_tuner_compat()
tuner <- tuner_obj$tuner
n_eval <- if (tuner_obj$mode == "mbo") 70 else 200
terminator <- make_terminator_evals(n_eval)

at <- AutoTuner$new(
  learner = glrn,
  resampling = resampling,
  measure = measure,
  search_space = search_space,
  terminator = terminator,
  tuner = tuner,
  store_models = FALSE
)

cat("Start tuning... tuner =", tuner_obj$mode, "evals =", n_eval, "\n")
at$train(task_train)

best_params <- at$tuning_result$learner_param_vals
cat("\nBest hyperparameters (CV-RMSE, weighted by weights_measure):\n")
print(best_params)

# Save best hyperparameters
saveRDS(best_params, file.path(out_dir, "best_params.rds"))
write.csv(as.data.frame(best_params), file.path(out_dir, "best_params.csv"), row.names = FALSE)
cat("Saved best_params to:", file.path(out_dir, "best_params.rds"), "\n")

# Save tuning archive
stringify_list_cols <- function(df) {
  if (!is.data.frame(df)) df <- as.data.frame(df)
  for (nm in names(df)) {
    if (is.list(df[[nm]])) {
      df[[nm]] <- vapply(df[[nm]], function(x) {
        if (is.null(x)) return(NA_character_)
        if (length(x) == 0) return("")
        if (is.atomic(x)) return(paste(as.character(x), collapse = ";"))
        paste(capture.output(dput(x)), collapse = "")
      }, character(1))
    }
  }
  df
}

arch <- at$tuning_instance$archive
archive_df <- NULL
if (!is.null(arch$data)) archive_df <- arch$data
if (is.null(archive_df)) archive_df <- tryCatch(arch$as_data_table(), error = function(e) NULL)
if (is.null(archive_df)) archive_df <- tryCatch(arch$as_data_frame(), error = function(e) NULL)

if (is.null(archive_df)) {
  warning("Cannot extract archive data, skipping tuning_archive.csv")
} else {
  archive_df <- as.data.frame(archive_df)
  archive_df <- stringify_list_cols(archive_df)
  utils::write.csv(archive_df, file.path(out_dir, "tuning_archive.csv"), row.names = FALSE)
}

# =========================
# 6) Independent test evaluation + save
# =========================
pred_test <- at$predict(task_test)

test_rmse_mlr3 <- pred_test$score(measure)

truth <- pred_test$truth
resp  <- pred_test$response
w_test <- task_test$data(cols = "weight")$weight
w_test_used <- w_test[pred_test$row_ids]

test_wrmse_manual <- rmse_weighted_manual(truth, resp, w_test_used)

test_r2 <- if (sum(is.finite(truth) & is.finite(resp)) >= 3) {
  cor(truth, resp, use = "complete.obs")^2
} else NA_real_

cat("\nIndependent test: RMSE(mlr3) =", test_rmse_mlr3,
    " | wRMSE(manual) =", test_wrmse_manual,
    " | R2 =", test_r2, "\n")

cv_best <- tryCatch(at$tuning_result$perf, error = function(e) NA_real_)
cv_best_num <- suppressWarnings(as.numeric(cv_best))[1]
if (!is.finite(cv_best_num)) cv_best_num <- NA_real_

test_rmse_num <- suppressWarnings(as.numeric(test_rmse_mlr3))[1]
if (!is.finite(test_rmse_num)) test_rmse_num <- NA_real_

test_wrmse_num <- suppressWarnings(as.numeric(test_wrmse_manual))[1]
if (!is.finite(test_wrmse_num)) test_wrmse_num <- NA_real_

test_r2_num <- suppressWarnings(as.numeric(test_r2))[1]
if (!is.finite(test_r2_num)) test_r2_num <- NA_real_

stats_df <- data.frame(
  metric = c("CV_best_RMSE_mlr3", "Test_RMSE_mlr3", "Test_wRMSE_manual", "Test_R2"),
  value  = c(cv_best_num, test_rmse_num, test_wrmse_num, test_r2_num)
)

write.csv(stats_df, file.path(out_dir, "model_performance_stats.csv"), row.names = FALSE)
write.csv(
  data.frame(CAP_obs = truth, CAP_pred = resp, weight = w_test_used),
  file.path(out_dir, "test_pred.csv"),
  row.names = FALSE
)

saveRDS(at, file.path(out_dir, "autotuner_final.rds"))
cat("\nAll results saved to:", out_dir, "\n")