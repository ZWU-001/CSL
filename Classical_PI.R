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

## linear model
fit_lm <- function(X, y) {
  lm(y ~ ., data = data.frame(X, y = y))
}

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

idx1 <- split_indices[1:n1]               
idx2 <- split_indices[(n1 + 1):(n1 + n2)] 
idx3 <- split_indices[(n1 + n2 + 1):n]    

training_set <- dataset[idx1, ]
calibration_set <- dataset[idx2, ]
testing_set <- dataset[idx3, ]

## training and calibration are combined
data_full <- rbind(training_set, calibration_set)

X_train <- data_full[, all_x_vars, drop = FALSE]
y_train <- data_full$y

X_test <- testing_set[, all_x_vars, drop = FALSE]
y_test_raw <- testing_set$y_raw

fit_lm_model <- fit_lm(X_train, y_train)

lm_pred <- predict(
  fit_lm_model,
  newdata = X_test,
  interval = "prediction",
  level = 0.90
)

L_log <- lm_pred[, "lwr"]
U_log <- lm_pred[, "upr"]

L <- exp(L_log)
U <- exp(U_log)

coverage_lm <- mean(y_test_raw >= L & y_test_raw <= U)
mean_width_lm <- mean(U - L)

round(coverage_lm, 3)
round(mean_width_lm, 4)
















