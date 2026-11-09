rm(list=ls())

library(readxl)
library(dplyr)
library(ggplot2)

## Set working directory
setwd("/Users/jameswu/Desktop/CSL_application/Github")

## Load the xlsx file
# files <- "results_full_csl_profiles_full_data.xlsx"

## With outliers removed
files <- "results_full_csl_profiles_outlier_removed.xlsx"

short_var_labels <- c(
  x1  = "LBXSAL",
  x2  = "LBXSTP",
  x3  = "LBXSUA",
  x4  = "LBXSTR",
  x5  = "LBXSCH",
  x6  = "LBXSCA",
  x7  = "LBXSPH",
  x8  = "RIDAGEYR",
  x9  = "BMXBMI",
  x10 = "LBXGH",
  x11 = "BMXWAIST",
  x12 = "BPXOSY1"
)

var_order <- c("x8", "x9", "x11", "x12", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x10")
label_order <- short_var_labels[var_order]  # ordered short labels

obs_profile_plot_df_all <- lapply(files, function(f) {
  read_xlsx(f) %>%
    mutate(
      value   = as.numeric(value),
      pred    = as.numeric(pred),
      L       = as.numeric(L),
      U       = as.numeric(U),
      width   = as.numeric(width),
      var     = as.character(var),
      label   = as.character(label),
      profile = as.character(profile)
    )
}) %>%
  bind_rows() %>%
  mutate(
    # Replace label with short abbreviation, ordered
    short_label = factor(short_var_labels[var], levels = label_order),
    profile = factor(profile, levels = c(
      "High observed profile",
      "Middle observed profile",
      "Low observed profile"
    ))
  )

p_obs_profiles_full_facet_all <- ggplot(
  obs_profile_plot_df_all,
  aes(x = value)
) +
  geom_ribbon(
    aes(ymin = L, ymax = U, fill = "Prediction interval"),
    alpha = 0.32
  ) +
  geom_line(
    aes(y = L, colour = "Lower/upper bound"),
    linewidth = 0.9
  ) +
  geom_line(
    aes(y = U, colour = "Lower/upper bound"),
    linewidth = 0.9
  ) +
  geom_line(
    aes(y = pred, colour = "Point prediction"),
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  facet_grid(profile ~ short_label, scales = "free_x") +
  scale_fill_manual(
    values = c("Prediction interval" = "#E8A6A1")
  ) +
  scale_colour_manual(
    values = c(
      "Lower/upper bound" = "blue",
      "Point prediction"  = "red"
    )
  ) +
  labs(
    x      = NULL,
    y      = "Prediction interval for serum creatinine (LBXSCR)",
    fill   = NULL,
    colour = NULL
  ) +
  theme_gray(base_size = 12) +
  coord_cartesian(ylim = c(0.3, 1.9)) +
  scale_y_continuous(breaks = c(0.5, 1.0, 1.5)) +
  theme(
    strip.text       = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90", colour = "grey60", linewidth = 0.6),
    panel.grid.major = element_line(linewidth = 0.45, colour = "grey80"),
    panel.grid.minor = element_line(linewidth = 0.25, colour = "grey88"),
    axis.title       = element_text(face = "bold"),
    axis.text        = element_text(colour = "black"),
    panel.spacing    = unit(0.6, "lines"),
    legend.position  = "bottom"
  )

print(p_obs_profiles_full_facet_all)

