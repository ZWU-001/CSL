rm(list=ls())

library(readxl)
library(dplyr)
library(ggplot2)

## Set working directory
setwd("/Users/jameswu/Desktop/CSL_application/Github")

## Load the xlsx file
# files <- "Results_Full_Task1_original.xlsx"

## With outliers removed
files <- "Results_Full_Task1_removal.xlsx"

plot_test_df <- lapply(files, read_xlsx) %>%
  bind_rows() %>%
  mutate(
    y_true = as.numeric(y_true),
    y_pred = as.numeric(y_pred),
    L      = as.numeric(L),
    U      = as.numeric(U)
  ) %>%
  filter(!is.na(y_true), !is.na(y_pred), !is.na(L), !is.na(U)) %>%
  mutate(
    y_pred = round(y_pred, 6),
    L      = round(L, 5),
    U      = round(U, 5)
  ) %>%
  arrange(y_true) %>%
  mutate(
    order_id = row_number(),
    width    = U - L,
    mid      = (L + U) / 2
  )

print(plot_test_df)

## Coverage and mean width
coverage   <- mean(plot_test_df$y_true >= plot_test_df$L &
                     plot_test_df$y_true <= plot_test_df$U)
mean_width <- mean(plot_test_df$width)

cat(sprintf("Coverage:   %.3f\n", coverage))
cat(sprintf("Mean width: %.4f\n", mean_width))

p_all_test <- ggplot(plot_test_df, aes(x = order_id)) +
  geom_linerange(
    aes(ymin = L, ymax = U),
    color     = "blue",
    alpha     = 0.35,
    linewidth = 0.3
  ) +
  geom_point(
    aes(y = y_pred, shape = "Point prediction", color = "Point prediction"),
    size   = 1.4,
    stroke = 0.7
  ) +
  geom_point(
    aes(y = y_true, shape = "True value", color = "True value"),
    size  = 1.4,
    alpha = 0.85
  ) +
  scale_shape_manual(
    values = c("Point prediction" = 4, "True value" = 16)
  ) +
  scale_color_manual(
    values = c("Point prediction" = "red", "True value" = "darkturquoise")
  ) +
  labs(
    x     = "Ordered Testing Point",
    y     = "Serum creatinine (LBXSCR)",
    shape = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(p_all_test)