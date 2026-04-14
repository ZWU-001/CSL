rm(list=ls())

## Required packages
library(gamlss)
library(nnls)
library(mgcv)
library(nnet)
library(randomForest)
library(glmnet)
library(nhanesA)
library(dplyr)
library(tidyr)
library(labelled)
library(ggplot2)

## Split conformal quantile threshold
conformal_threshold <- function(scores, alpha) {
  m <- length(scores)
  k <- ceiling((m + 1) * (1 - alpha))
  sort(scores, decreasing = FALSE)[k]
}

## Base learner 1: linear model
fit_lm <- function(X, y) {
  lm(y ~ ., data = data.frame(X, y = y))
}
pred_lm <- function(fit, newX) {
  as.numeric(predict(fit, newdata = newX))
}

## Base learner 2: GAM
fit_gam <- function(X, y) {
  mgcv::gam(
    y ~ s(x1) + s(x2) + s(x3) + s(x4) +
      s(x5) + s(x6) + s(x7) + s(x8) +
      s(x9) + s(x10) + s(x11) + s(x12) +
      factor(x13) + factor(x14) + factor(x15) +
      factor(x16) + factor(x17) + factor(x18),
    data = data.frame(X, y = y),
    family = gaussian()
  )
}
pred_gam <- function(fit, newX) {
  as.numeric(predict(fit, newdata = newX, type = "response"))
}

## Base learner 3: neural network
fit_nnet <- function(X, y) {
  nnet::nnet(y ~ ., data = data.frame(X, y = y), size = 5,
             linout = TRUE, trace = FALSE)
}
pred_nnet <- function(fit, newX) {
  as.numeric(predict(fit, newdata = newX, type = "raw"))
}

## Base learner 4: random forest
fit_rf <- function(X, y) {
  randomForest::randomForest(x = X, y = y)
}
pred_rf <- function(fit, newX) {
  as.numeric(predict(fit, newdata = newX))
}

## Base learner 5: LASSO
fit_lasso <- function(X, y) {
  Xmm <- model.matrix(~ ., data = X)[, -1]
  glmnet::cv.glmnet(Xmm, y)
}

pred_lasso <- function(fit, newX) {
  Xmm_new <- model.matrix(~ ., data = newX)[, -1]
  as.numeric(predict(fit, newx = Xmm_new, s = "lambda.min"))
}

## Wrapper for fitting the first five base learners
fit_base_manual <- function(k, X, y) {
  if (k == 1) return(fit_lm(X, y))
  if (k == 2) return(fit_gam(X, y))
  if (k == 3) return(fit_nnet(X, y))
  if (k == 4) return(fit_rf(X, y))
  if (k == 5) return(fit_lasso(X, y))
  stop("k must be 1..5 for base learners.")
}

## Wrapper for predicting from the first five base learners
pred_base_manual <- function(k, fit, newX) {
  if (k == 1) return(pred_lm(fit, newX))
  if (k == 2) return(pred_gam(fit, newX))
  if (k == 3) return(pred_nnet(fit, newX))
  if (k == 4) return(pred_rf(fit, newX))
  if (k == 5) return(pred_lasso(fit, newX))
  stop("k must be 1..5 for base learners.")
}

## Super Learner library
SL.library <- c("SL.lm", "SL.gam", "SL.nnet", "SL.randomForest", "SL.lasso", "SL.gamlss")

# Parameters
alpha = 0.1
V = 10

## ------------------------------------------------------------
## Data extraction: NHANES 2021-2023
## ------------------------------------------------------------
demo_raw <- nhanes('DEMO_L')     
bmx_raw  <- nhanes('BMX_L')      
lab_raw  <- nhanes('BIOPRO_L')   
tri_raw  <- nhanes('TRIGLY_L')
kiq_raw  <- nhanes("KIQ_U_L")    
mcq_raw  <- nhanes("MCQ_L")
dia_raw  <- nhanes("DIQ_L")
ghb_raw  <- nhanes("GHB_L")
bpq_raw  <- nhanes("BPQ_L")
rx_raw   <- nhanes("RXQ_RX_L")
bpxo_raw <- nhanes("BPXO_L")

## Identify participants with kidney conditions
kiq_flag <- kiq_raw %>%
  transmute(
    SEQN,
    KIQ025_chr = to_character(KIQ025)
  ) %>%
  mutate(
    KIQ025_chr = case_when(
      KIQ025_chr == "Yes" ~ "Yes",
      KIQ025_chr == "No" ~ "No",
      TRUE ~ NA_character_
    )
  )

seqn_kiq_yes <- kiq_flag %>%
  filter(KIQ025_chr == "Yes") %>%
  pull(SEQN)

## Merge the NHANES modules
df_full <- demo_raw %>%
  left_join(bmx_raw, by = "SEQN") %>%
  left_join(lab_raw, by = "SEQN") %>%
  left_join(tri_raw, by = "SEQN") %>%
  left_join(mcq_raw, by = "SEQN") %>%
  left_join(dia_raw, by = "SEQN") %>%
  left_join(ghb_raw, by = "SEQN") %>%
  left_join(bpq_raw, by = "SEQN") %>%
  left_join(rx_raw, by = "SEQN") %>%
  left_join(bpxo_raw, by = "SEQN") %>%
  select(
    SEQN, 
    RIAGENDR, RIDAGEYR, RIDRETH3, 
    BMXWT, BMXHT, BMXBMI, RIDEXPRG,
    MCQ160L, DIQ010, LBXGH, BPQ020, RXQ033, 
    BMXWAIST, BPXOSY1, 
    contains("SCR"),  
    contains("SUA"),  
    contains("SAL"),  
    contains("STP"),  
    matches("CH"),
    matches("STR"),
    matches("SCA"),
    matches("SPH"),
    matches("ALP"),
    matches("ALT"),
    matches("AST")
  ) %>%
  filter(!is.na(BMXWT))

## Remove participants with kidney conditions
df_full <- df_full %>%
  filter(!(SEQN %in% seqn_kiq_yes))

## Recode selected categorical variables
df_full <- df_full %>%
  mutate(
    RIAGENDR_chr = to_character(RIAGENDR),
    RIDEXPRG_chr = to_character(RIDEXPRG),
    MCQ160L_chr  = to_character(MCQ160L),
    DIQ010_chr   = to_character(DIQ010),
    BPQ020_chr   = to_character(BPQ020),
    RXQ033_chr   = to_character(RXQ033),
    RIDRETH3_chr = to_character(RIDRETH3)
  ) %>%
  mutate(
    ## Pregnancy status:
    ## men and non-pregnant/uncertain women are coded as "No"
    RIDEXPRG_chr = case_when(
      RIAGENDR_chr == "Male" ~ "No",
      RIAGENDR_chr == "Female" & is.na(RIDEXPRG_chr) ~ "No",
      RIDEXPRG_chr == "The participant was not pregnant at exam" ~ "No",
      RIDEXPRG_chr == "Yes, positive lab pregnancy test or self-reported pregnant at exam" ~ "Yes",
      RIDEXPRG_chr == "Cannot ascertain if the participant is pregnant at exam" ~ "No",
      TRUE ~ RIDEXPRG_chr
    ),
    
    ## Liver condition
    MCQ160L_chr = case_when(
      MCQ160L_chr == "Yes" ~ "Yes",
      MCQ160L_chr == "No" ~ "No",
      MCQ160L_chr %in% c("Refused", "Don't know") ~ NA_character_,
      TRUE ~ NA_character_
    ),
    
    ## Diabetes status
    DIQ010_chr = case_when(
      DIQ010_chr == "Yes" ~ "Yes",
      DIQ010_chr == "No" ~ "No",
      DIQ010_chr == "Borderline" ~ "Borderline",
      DIQ010_chr %in% c("Refused", "Don't know") ~ NA_character_,
      TRUE ~ NA_character_
    ),
    
    ## High blood pressure
    BPQ020_chr = case_when(
      BPQ020_chr == "Yes" ~ "Yes",
      BPQ020_chr == "No" ~ "No",
      BPQ020_chr %in% c("Refused", "Don't know") ~ NA_character_,
      TRUE ~ NA_character_
    ),
    
    ## Prescription medication use
    RXQ033_chr = case_when(
      RXQ033_chr == "Yes" ~ "Yes",
      RXQ033_chr == "No" ~ "No",
      RXQ033_chr %in% c("Refused", "Don't know") ~ NA_character_,
      TRUE ~ NA_character_
    )
  )

## Exclude pregnant participants
df_processed <- df_full %>%
  filter(RIDEXPRG_chr != "Yes")

## Select final variables
sl_data <- df_processed %>%
  select(
    LBXSAL, LBXSTP, LBXSCR, LBXSUA, LBXSTR, LBXSCH, LBXSCA, LBXSPH,
    RIAGENDR, RIDAGEYR, BMXBMI,
    MCQ160L_chr, DIQ010_chr, LBXGH, BPQ020_chr, RXQ033_chr,
    BMXWAIST, BPXOSY1, RIDRETH3_chr
  ) %>%
  drop_na()

## Construct final analysis data set
## outcome is modelled on the log scale
dataset <- data.frame(
  y_raw = sl_data$LBXSCR,
  y   = log(sl_data$LBXSCR),
  x1  = sl_data$LBXSAL,
  x2  = sl_data$LBXSTP,
  x3  = sl_data$LBXSUA,
  x4  = sl_data$LBXSTR,
  x5  = sl_data$LBXSCH,
  x6  = sl_data$LBXSCA,
  x7  = sl_data$LBXSPH,
  x8  = sl_data$RIDAGEYR,
  x9  = sl_data$BMXBMI,
  x10 = sl_data$LBXGH,
  x11 = sl_data$BMXWAIST,
  x12 = sl_data$BPXOSY1,
  x13 = factor(as.integer(sl_data$RIAGENDR), levels = c(1, 2), labels = c("Male", "Female")),
  x14 = factor(sl_data$MCQ160L_chr, levels = c("Yes", "No")),
  x15 = factor(sl_data$DIQ010_chr, levels = c("Yes", "No", "Borderline")),
  x16 = factor(sl_data$BPQ020_chr, levels = c("Yes", "No")),
  x17 = factor(sl_data$RXQ033_chr, levels = c("Yes", "No")),
  x18 = factor(sl_data$RIDRETH3_chr)
)

## Optional sensitivity analysis:
## remove observations with serum creatinine > 1.5 mg/dL
# dataset_clean <- dataset[dataset$y_raw <= 1.5, c("y_raw", "y", paste0("x", 1:18))]
# cat("Original sample size:", nrow(dataset), "\n")
# cat("Cleaned sample size :", nrow(dataset_clean), "\n")
# cat("Number removed      :", sum(dataset$y_raw > 1.5, na.rm = TRUE), "\n")
# dataset <- dataset_clean

## Predictor names
all_x_vars <- paste0("x", 1:18)

## ------------------------------------------------------------
## Train / calibration / test split
## ------------------------------------------------------------
set.seed(66666)

n <- nrow(dataset)
split_indices <- sample(1:n)

n3 <- floor(0.1 * n)
n1 <- floor((n - n3) * 0.8)
n2 <- n - n3 - n1

idx1 <- split_indices[1:n1]               # Training indices
idx2 <- split_indices[(n1 + 1):(n1 + n2)] # Calibration indices
idx3 <- split_indices[(n1 + n2 + 1):n]    # Testing indices

training_set <- dataset[idx1, ]
calibration_set <- dataset[idx2, ]
testing_set <- dataset[idx3, ]

X_train <- training_set[, all_x_vars, drop = FALSE]
y_train <- training_set$y
n_train <- length(y_train)

X_calib <- calibration_set[, all_x_vars, drop = FALSE]
y_calib <- calibration_set$y

X_test <- testing_set[, all_x_vars, drop = FALSE]
y_test <- testing_set$y
y_test_raw <- testing_set$y_raw
n_test <- nrow(X_test)

## Objects for storing final results
L_final <- rep(NA, n_test)
U_final <- rep(NA, n_test)
y_pred_plot <- rep(NA_real_, n_test)

## Record indices where WMV yields multiple disjoint regions
multi_interval_idx <- integer(0)
multi_interval_count <- 0

## ------------------------------------------------------------
## Super Learner weight estimation via V-fold CV
## ------------------------------------------------------------
fold_id <- sample(rep(1:V, length.out = n_train))
folds_list <- lapply(1:V, function(v) which(fold_id == v))

## Cross-validated predictions for the first five learners
Z_base <- matrix(NA, nrow = n_train, ncol = 5)
colnames(Z_base) <- c("lm", "gam", "nnet", "rf", "lasso")

for (v in 1:V) {
  idx_te <- folds_list[[v]]
  idx_tr <- setdiff(1:n_train, idx_te)
  
  X_tr <- X_train[idx_tr, , drop = FALSE]
  y_tr <- y_train[idx_tr]
  X_te <- X_train[idx_te, , drop = FALSE]
  
  for (k in 1:5) {
    fit_k_cv <- fit_base_manual(k, X_tr, y_tr)
    Z_base[idx_te, k] <- pred_base_manual(k, fit_k_cv, X_te)
  }
}

## Cross-validated predictions for GAMLSS
Z_gamlss <- rep(NA, n_train)

for (v in 1:V) {
  idx_te <- folds_list[[v]]
  idx_tr <- setdiff(1:n_train, idx_te)
  
  data_fold_train <- data.frame(
    X_train[idx_tr, , drop = FALSE],
    y = y_train[idx_tr]
  )
  
  data_fold_test <- X_train[idx_te, , drop = FALSE]
  
  fit_g <- gamlss(
    formula        = y ~ .,
    sigma.formula  = ~ .,
    data           = data_fold_train,
    family         = NO(),
    trace          = TRUE
  )
  
  converged <- fit_g$converged
  if (!converged) warning("GAMLSS did not converge in fold ")
  
  Z_gamlss[idx_te] <- predict(fit_g, newdata = data_fold_test, what = "mu", type = "response")
}

## Combine all cross-validated predictions
Z_all <- cbind(Z_base, Z_gamlss)

## Estimate Super Learner weights by non-negative least squares
w_hat <- nnls(Z_all, y_train)$x
w_hat <- w_hat / sum(w_hat)

print(w_hat)


## Active learners: positive weight
## Dominant learner: weight greater than 0.5
active_idx <- which(w_hat > 0)
dominant_idx <- which(w_hat > 0.5)

## ------------------------------------------------------------
## Split-CSL performance on the testing set
## ------------------------------------------------------------
if (length(dominant_idx) == 1){
  
  ## Case 1: a dominant learner exists
  k <- dominant_idx
  cat("\nDominant learner detected, k =", k, " weight =", w_hat[k], "\n")
  
  if (k %in% 1:5){
    
    ## Base learners 1--5:
    ## use absolute residuals as non-conformity scores
    fit_k <- fit_base_manual(k, X_train, y_train)
    yhat_cal <- pred_base_manual(k, fit_k, X_calib)
    q_k <- conformal_threshold(abs(y_calib - yhat_cal), alpha)
    
    for (i in 1:n_test) {
      newX_i <- X_test[i, , drop = FALSE]
      yhat_i <- pred_base_manual(k, fit_k, newX_i)
      
      y_pred_plot[i] <- yhat_i
      L_final[i] <- yhat_i - q_k
      U_final[i] <- yhat_i + q_k
    }
  } else {
    
    ## GAMLSS learner:
    ## use absolute quantile residuals as non-conformity scores
    data_train_gamlss <- data.frame(X_train, y = y_train)
    
    fit_g <- gamlss(
      formula       = y ~ .,
      sigma.formula = ~ .,
      data          = data_train_gamlss,
      family        = NO(),
      trace         = TRUE
    )
    
    print(paste("GAMLSS converged:", fit_g$converged))
    
    mu_cal <- predict(fit_g, newdata = X_calib, what = "mu", type = "response")
    sg_cal <- predict(fit_g, newdata = X_calib, what = "sigma", type = "response")
    
    F_cal <- pNO(y_calib, mu = mu_cal, sigma = sg_cal)
    qres_cal <- qnorm(F_cal)
    
    q_gamlss <- conformal_threshold(abs(qres_cal), alpha)
    
    p_lower <- pnorm(-q_gamlss)
    p_upper <- pnorm(q_gamlss)
    
    for (i in 1:n_test) {
      newX_i <- X_test[i, , drop = FALSE]
      
      mu_i <- predict(fit_g, newdata = newX_i, what = "mu", type = "response")
      sg_i <- predict(fit_g, newdata = newX_i, what = "sigma", type = "response")
      
      y_pred_plot[i] <- mu_i
      L_final[i] <- qNO(p_lower, mu = mu_i, sigma = sg_i)
      U_final[i] <- qNO(p_upper, mu = mu_i, sigma = sg_i)
    }
  }
} else{
  
  ## Case 2: no dominant learner
  ## construct learner-specific split conformal intervals and combine by WMV
  cat("\nNo dominant learner. Construct split CP intervals for active learners...\n")
  cat("Active learners:", active_idx, "\n")
  
  intervals_k <- matrix(NA, nrow = 2, ncol = length(SL.library))
  base_intervals <- list()
  
  for (k in active_idx){
    
    if (k %in% 1:5) {
      
      fit_k <- fit_base_manual(k, X_train, y_train)
      yhat_cal <- pred_base_manual(k, fit_k, X_calib)
      q_k <- conformal_threshold(abs(y_calib - yhat_cal), alpha)
      
      base_intervals[[k]] <- list(fit = fit_k, q = q_k)
      
    } else {
      
      data_train_gamlss <- data.frame(X_train, y = y_train)
      
      fit_g <- gamlss(
        formula       = y ~ .,
        sigma.formula = ~ .,
        data          = data_train_gamlss,
        family        = NO(),
        trace         = TRUE
      )
      
      print(paste("GAMLSS converged:", fit_g$converged))
      
      mu_cal <- predict(fit_g, newdata = X_calib, what = "mu", type = "response")
      sg_cal <- predict(fit_g, newdata = X_calib, what = "sigma", type = "response")
      
      F_cal <- pNO(y_calib, mu = mu_cal, sigma = sg_cal)
      qres_cal <- qnorm(F_cal)
      
      q_gamlss <- conformal_threshold(abs(qres_cal), alpha)
      
      base_intervals[[k]] <- list(fit = fit_g, q = q_gamlss)
    }
  }
  
  for (i in 1:n_test){
    
    print(i)
    
    newX_i <- X_test[i, , drop = FALSE]
    pred_i_weighted <- 0
    
    for (k in active_idx){
      
      if (k %in% 1:5) {
        
        fit_k <- base_intervals[[k]]$fit
        q_k <- base_intervals[[k]]$q
        
        pred_i <- pred_base_manual(k, fit_k, newX_i)
        pred_i_weighted <- pred_i_weighted + w_hat[k] * pred_i
        
        intervals_k[1, k] <- pred_i - q_k
        intervals_k[2, k] <- pred_i + q_k
        
      } else{
        
        fit_g <- base_intervals[[k]]$fit
        q_gamlss <- base_intervals[[k]]$q
        
        mu_i <- predict(fit_g, newdata = newX_i, what = "mu", type = "response")
        sg_i <- predict(fit_g, newdata = newX_i, what = "sigma", type = "response")
        
        pred_i_weighted <- pred_i_weighted + w_hat[k] * mu_i
        
        p_lower <- pnorm(-q_gamlss)
        p_upper <- pnorm(q_gamlss)
        
        intervals_k[1, k] <- qNO(p_lower, mu = mu_i, sigma = sg_i)
        intervals_k[2, k] <- qNO(p_upper, mu = mu_i, sigma = sg_i)
      }
    }
    
    y_pred_plot[i] <- pred_i_weighted
    
    ## Combine learner-specific intervals by weighted majority vote
    endpoints <- sort(unique(c(intervals_k[1, active_idx], intervals_k[2, active_idx])))
    
    range_extension <- 0.1 * (max(endpoints) - min(endpoints))
    endpoints <- c(min(endpoints) - range_extension,
                   endpoints,
                   max(endpoints) + range_extension)
    
    test_points <- seq(min(endpoints), max(endpoints), by = 0.0001)
    vote_scores <- numeric(length(test_points))
    
    for (j in seq_along(test_points)) {
      y_grid <- test_points[j]
      in_interval <- numeric(length(SL.library))
      
      for (k in active_idx) {
        if (y_grid >= intervals_k[1, k] && y_grid <= intervals_k[2, k]) {
          in_interval[k] <- w_hat[k]
        }
      }
      vote_scores[j] <- sum(in_interval)
    }
    
    above_threshold <- vote_scores > 0.5
    
    if (any(above_threshold)) {
      
      regions <- rle(above_threshold)
      
      ## Record cases where WMV produces multiple disjoint accepted regions
      n_regions <- sum(regions$values)
      if (n_regions > 1) {
        multi_interval_idx <- c(multi_interval_idx, i)
        multi_interval_count <- multi_interval_count + 1
      }
      
      ## Retain the accepted region
      ## needs further implementation if multiple accepted regions occur
      start_idx <- 1
      best_region <- NULL
      max_width <- -Inf
      
      for (u in seq_along(regions$values)) {
        if (regions$values[u]) {
          region_start <- start_idx
          region_end <- start_idx + regions$lengths[u] - 1
          region_width <- test_points[region_end] - test_points[region_start]
          
          if (region_width > max_width) {
            max_width <- region_width
            best_region <- list(start = region_start, end = region_end)
          }
        }
        start_idx <- start_idx + regions$lengths[u]
      }
      
      L_final[i] <- test_points[best_region$start]
      U_final[i] <- test_points[best_region$end]
      
    } else {
      L_final[i] <- NA
      U_final[i] <- NA
    }
  }
}

## ------------------------------------------------------------
## Back-transform intervals and predictions to the original scale
## ------------------------------------------------------------
y_pred_plot_raw <- exp(y_pred_plot)
L_final_raw <- exp(L_final)
U_final_raw <- exp(U_final)

## Empirical coverage and mean interval width on the test set
cover_test <- (y_test_raw >= L_final_raw) & (y_test_raw <= U_final_raw)
width_test <- U_final_raw - L_final_raw

cat("\n--- APPLICATION (Split CSL / WMV) ---\n")
cat("Test size:", n_test, "\n")
cat("Coverage:", round(mean(cover_test, na.rm = TRUE), 3), "\n")
cat("Mean width:", round(mean(width_test, na.rm = TRUE), 4), "\n")

cat("\nTest points with MORE THAN ONE prediction interval after WMV combination:\n")
print(multi_interval_idx)
cat("Total count:", multi_interval_count, "\n")

## ------------------------------------------------------------
## Ordered test-set interval plot
## ------------------------------------------------------------
plot_test_df <- data.frame(
  id     = seq_len(n_test),
  y_true = y_test_raw,
  y_pred = y_pred_plot_raw,
  L      = L_final_raw,
  U      = U_final_raw
) %>%
  filter(!is.na(y_true), !is.na(y_pred), !is.na(L), !is.na(U)) %>%
  arrange(y_true) %>%
  mutate(order_id = row_number()) %>%
  mutate(
    width = U - L,
    mid   = (L + U) / 2
  )

print(plot_test_df)

p_all_test <- ggplot(plot_test_df, aes(x = order_id)) +
  geom_linerange(
    aes(ymin = L, ymax = U),
    color = "blue",
    alpha = 0.35,
    linewidth = 0.3
  ) +
  geom_point(
    aes(y = y_pred, shape = "Point prediction", color = "Point prediction"),
    size = 1.4,
    stroke = 0.7
  ) +
  geom_point(
    aes(y = y_true, shape = "True value", color = "True value"),
    size = 1.4,
    alpha = 0.85
  ) +
  scale_shape_manual(
    values = c("Point prediction" = 4, "True value" = 16)
  ) +
  scale_color_manual(
    values = c("Point prediction" = "red", "True value" = "darkturquoise")
  ) +
  labs(
    x = "Ordered Testing Point",
    y = "Serum creatinine (LBXSCR)",
    shape = NULL,
    color = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

print(p_all_test)










