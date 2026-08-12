
library(tidyfst)
pkg_load(tidyverse,fs,arrow,DBI,duckdb)

con<- dbConnect(
  duckdb(),
  dbdir = "data/sel_data.duckdb"
)

dbListTables(con)
tbl(con,"review_citation") %>% 
  print(width = Inf)

tbl(con,"review_citation") %>% 
  transmute(cnp = citation_normalized_percentile,pr_duration) %>% 
  collect() -> overall_test

# =============================================================================
# Model selection: review duration (pr_duration) and normalized citation
# impact (cnp)
#
# Goal: compare a small set of candidate functional forms and identify which
# one best describes the conditional mean of cnp given pr_duration.
#
# Selection strategy (no cross-validation):
#   1. BIC on the full analysis sample  -> primary criterion
#   2. A single random 80/20 hold-out   -> one honest out-of-sample RMSE
#   3. Adjusted R-squared and effect sizes -> practical relevance
#
# IMPORTANT INTERPRETATION LIMIT
# The data contain only two variables. No journal, field, year or author
# information is available, so every result below is an UNADJUSTED
# statistical association. It must not be read as a causal effect of review
# duration on citation impact.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Packages
# -----------------------------------------------------------------------------

library(tidyfst)

pkg_load(
  dplyr,
  tidyr,
  purrr,
  tibble,
  ggplot2,
  splines,
  segmented,
  broom,
  gt,
  conflicted
)

# Force the dplyr versions of the generic verbs. Several loaded packages
# define select() and filter() as well.
conflict_prefer("select", "dplyr", quiet = TRUE)
conflict_prefer("filter", "dplyr", quiet = TRUE)
conflict_prefer("lag", "dplyr", quiet = TRUE)

set.seed(20260801)


# -----------------------------------------------------------------------------
# 2. Output location
# -----------------------------------------------------------------------------

# The project root is expected to already contain a data folder. Only the
# which_model subfolder is created here.
if (!dir.exists("data")) {
  stop(
    "No 'data' folder found in the working directory. ",
    "Set the working directory to the project root before running this script."
  )
}

output_dir <- file.path("data", "which_model")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

out_path <- function(...) file.path(output_dir, ...)


# -----------------------------------------------------------------------------
# 3. Analysis sample
# -----------------------------------------------------------------------------

# overall_test is assumed to be already in the environment:
#   cnp          numeric, normalized citation impact
#   pr_duration  numeric, review duration in days

stopifnot(exists("overall_test"))

raw_n <- nrow(overall_test)

# Extremely long review durations are rare and would dominate the curvature of
# any flexible model. The upper 0.5% is trimmed rather than winsorized so that
# no artificial pile-up is created at the boundary.
duration_cap <- quantile(overall_test$pr_duration, 0.995, na.rm = TRUE)

model_data <- overall_test %>%
  filter(
    !is.na(cnp),
    !is.na(pr_duration),
    pr_duration > 0,
    pr_duration <= duration_cap
  ) %>%
  mutate(log_duration = log(pr_duration))

analysis_n <- nrow(model_data)

message(
  "Rows supplied: ", format(raw_n, big.mark = ","),
  " | rows analysed: ", format(analysis_n, big.mark = ","),
  " | duration cap: ", round(duration_cap, 1), " days"
)


# -----------------------------------------------------------------------------
# 4. Single hold-out split
# -----------------------------------------------------------------------------

# One random split only. This is deliberately not repeated cross-validation:
# with millions of rows the sampling error of a 20% hold-out RMSE is already
# negligible relative to the differences between candidate models.

holdout_flag <- runif(analysis_n) < 0.2

train_data <- model_data[!holdout_flag, ]
test_data  <- model_data[holdout_flag, ]

message(
  "Training rows: ", format(nrow(train_data), big.mark = ","),
  " | hold-out rows: ", format(nrow(test_data), big.mark = ",")
)


# -----------------------------------------------------------------------------
# 5. Break points for the segmented candidates
# -----------------------------------------------------------------------------

# Break point locations are estimated on a modest subsample. Locating a change
# point does not require millions of observations, and this is by far the most
# expensive step if run on the full data. Once located, the break points are
# treated as fixed and the piecewise model is refitted on the full sample with
# ordinary lm(), which is cheap.

break_sample <- model_data %>% slice_sample(n = min(150000L, analysis_n))

estimate_breaks <- function(n_breaks) {
  base_fit <- lm(cnp ~ pr_duration, data = break_sample)
  
  start_values <- as.numeric(quantile(
    break_sample$pr_duration,
    probs = seq_len(n_breaks) / (n_breaks + 1)
  ))
  
  attempt <- try(
    segmented(
      base_fit,
      seg.Z = ~pr_duration,
      psi = start_values,
      control = seg.control(n.boot = 0, it.max = 30)
    ),
    silent = TRUE
  )
  
  if (inherits(attempt, "try-error")) {
    return(NULL)
  }
  
  list(
    psi = as.numeric(attempt$psi[, "Est."]),
    se  = as.numeric(attempt$psi[, "St.Err"])
  )
}

break_results <- map(1:3, estimate_breaks)
names(break_results) <- paste0("segmented_", 1:3, "bp")


# -----------------------------------------------------------------------------
# 6. Candidate model formulas
# -----------------------------------------------------------------------------

# Piecewise linear terms are built explicitly as hinge functions so that the
# model is an ordinary lm() and can be fitted on the full sample.
hinge_terms <- function(psi) {
  paste0("pmax(pr_duration - ", round(psi, 4), ", 0)", collapse = " + ")
}

candidate_formulas <- list(
  linear        = cnp ~ pr_duration,
  log_linear    = cnp ~ log_duration,
  quadratic     = cnp ~ pr_duration + I(pr_duration^2),
  spline_df3    = cnp ~ ns(pr_duration, df = 3),
  spline_df4    = cnp ~ ns(pr_duration, df = 4),
  spline_df5    = cnp ~ ns(pr_duration, df = 5)
)

for (nm in names(break_results)) {
  br <- break_results[[nm]]
  if (is.null(br)) next
  candidate_formulas[[nm]] <- as.formula(
    paste("cnp ~ pr_duration +", hinge_terms(br$psi))
  )
}

candidate_labels <- c(
  linear          = "Linear",
  log_linear      = "Log-linear in duration",
  quadratic       = "Quadratic",
  spline_df3      = "Natural spline, 3 df",
  spline_df4      = "Natural spline, 4 df",
  spline_df5      = "Natural spline, 5 df",
  segmented_1bp   = "Piecewise linear, 1 break",
  segmented_2bp   = "Piecewise linear, 2 breaks",
  segmented_3bp   = "Piecewise linear, 3 breaks"
)


# -----------------------------------------------------------------------------
# 7. Fit every candidate twice: full sample and training sample
# -----------------------------------------------------------------------------

rmse <- function(observed, predicted) {
  sqrt(mean((observed - predicted)^2))
}

mae <- function(observed, predicted) {
  mean(abs(observed - predicted))
}

evaluate_candidate <- function(model_name) {
  form <- candidate_formulas[[model_name]]
  
  full_fit  <- lm(form, data = model_data)
  train_fit <- lm(form, data = train_data)
  
  holdout_pred <- predict(train_fit, newdata = test_data)
  full_summary <- summary(full_fit)
  
  tibble(
    model          = model_name,
    label          = candidate_labels[[model_name]],
    n_parameters   = length(coef(full_fit)),
    bic            = BIC(full_fit),
    aic            = AIC(full_fit),
    adj_r_squared  = full_summary$adj.r.squared,
    holdout_rmse   = rmse(test_data$cnp, holdout_pred),
    holdout_mae    = mae(test_data$cnp, holdout_pred)
  )
}

model_metrics <- map_dfr(names(candidate_formulas), evaluate_candidate)


# -----------------------------------------------------------------------------
# 8. Comparison table
# -----------------------------------------------------------------------------

linear_bic  <- model_metrics$bic[model_metrics$model == "linear"]
linear_rmse <- model_metrics$holdout_rmse[model_metrics$model == "linear"]

comparison_table <- model_metrics %>%
  mutate(
    delta_bic          = bic - min(bic),
    bic_gain_vs_linear = linear_bic - bic,
    rmse_gain_pct      = 100 * (linear_rmse - holdout_rmse) / linear_rmse,
    is_best_bic        = bic == min(bic),
    is_best_holdout    = holdout_rmse == min(holdout_rmse)
  ) %>%
  arrange(bic) %>%
  mutate(bic_rank = row_number()) %>%
  select(
    bic_rank,
    model,
    label,
    n_parameters,
    bic,
    delta_bic,
    bic_gain_vs_linear,
    aic,
    adj_r_squared,
    holdout_rmse,
    holdout_mae,
    rmse_gain_pct,
    is_best_bic,
    is_best_holdout
  )

best_model_name <- comparison_table$model[1]
best_model_label <- comparison_table$label[1]
best_fit <- lm(candidate_formulas[[best_model_name]], data = model_data)

message("Lowest BIC: ", best_model_label)
message(
  "Lowest hold-out RMSE: ",
  comparison_table$label[which.min(comparison_table$holdout_rmse)]
)


# -----------------------------------------------------------------------------
# 9. Practical-relevance check
# -----------------------------------------------------------------------------

# With this sample size almost anything is statistically significant, so the
# question that matters is how much predictive accuracy actually improves and
# how large the fitted change in cnp is across the observed duration range.

duration_grid <- tibble(
  pr_duration = seq(
    min(model_data$pr_duration),
    max(model_data$pr_duration),
    length.out = 400
  )
) %>%
  mutate(log_duration = log(pr_duration))

grid_prediction <- duration_grid %>%
  mutate(fitted_cnp = predict(best_fit, newdata = duration_grid))

fitted_range <- diff(range(grid_prediction$fitted_cnp))
cnp_sd <- sd(model_data$cnp)

effect_summary <- tibble(
  quantity = c(
    "Rows analysed",
    "Duration cap (days)",
    "Best model by BIC",
    "Adjusted R-squared of best model",
    "Hold-out RMSE improvement over linear (%)",
    "Fitted cnp range across durations",
    "Standard deviation of cnp",
    "Fitted range as share of one SD"
  ),
  value = c(
    format(analysis_n, big.mark = ","),
    format(round(duration_cap, 1), nsmall = 1),
    best_model_label,
    format(round(comparison_table$adj_r_squared[1], 5), nsmall = 5),
    format(round(comparison_table$rmse_gain_pct[1], 3), nsmall = 3),
    format(round(fitted_range, 4), nsmall = 4),
    format(round(cnp_sd, 4), nsmall = 4),
    format(round(fitted_range / cnp_sd, 3), nsmall = 3)
  )
)

print(effect_summary, n = Inf)


# -----------------------------------------------------------------------------
# 10. Break points and stage slopes, if a piecewise model is competitive
# -----------------------------------------------------------------------------

# Piecewise models are the only candidates whose parameters translate directly
# into a sentence such as "before day X impact declines by Y per 10 days".

stage_slope_table <- NULL
break_point_table <- NULL

best_piecewise <- comparison_table %>%
  filter(grepl("^segmented", model)) %>%
  slice(1)

if (nrow(best_piecewise) == 1) {
  pw_name <- best_piecewise$model
  pw_psi <- break_results[[pw_name]]$psi
  pw_se  <- break_results[[pw_name]]$se
  pw_fit <- lm(candidate_formulas[[pw_name]], data = model_data)
  
  break_point_table <- tibble(
    break_index = seq_along(pw_psi),
    break_day = round(pw_psi, 1),
    ci_lower_day = round(pw_psi - 1.96 * pw_se, 1),
    ci_upper_day = round(pw_psi + 1.96 * pw_se, 1)
  )
  
  # Cumulative sums of the hinge coefficients give the slope in each stage.
  pw_coef <- coef(pw_fit)
  slope_terms <- pw_coef[-1]
  stage_slopes <- cumsum(slope_terms)
  
  stage_bounds <- c(min(model_data$pr_duration), pw_psi,
                    max(model_data$pr_duration))
  
  stage_slope_table <- tibble(
    stage = seq_along(stage_slopes),
    from_day = round(head(stage_bounds, -1), 1),
    to_day = round(stage_bounds[-1], 1),
    slope_per_day = as.numeric(stage_slopes),
    change_per_10_days = 10 * as.numeric(stage_slopes)
  ) %>%
    mutate(
      direction = if_else(slope_per_day < 0, "decreasing", "increasing")
    )
  
  print(break_point_table, n = Inf)
  print(stage_slope_table, n = Inf)
}


# -----------------------------------------------------------------------------
# 11. Presentation table
# -----------------------------------------------------------------------------

display_table <- comparison_table %>%
  mutate(
    selected = case_when(
      is_best_bic & is_best_holdout ~ "Best on both criteria",
      is_best_bic ~ "Lowest BIC",
      is_best_holdout ~ "Lowest hold-out RMSE",
      TRUE ~ ""
    )
  ) %>%
  select(
    Rank = bic_rank,
    Model = label,
    Parameters = n_parameters,
    BIC = bic,
    `Delta BIC` = delta_bic,
    `Adj. R2` = adj_r_squared,
    `Hold-out RMSE` = holdout_rmse,
    `Hold-out MAE` = holdout_mae,
    `RMSE gain vs linear (%)` = rmse_gain_pct,
    Note = selected
  )

model_gt <- display_table %>%
  gt() %>%
  tab_header(
    title = "Which model best describes review duration and citation impact?",
    subtitle = paste0(
      "Ranked by BIC on ", format(analysis_n, big.mark = ","),
      " observations; hold-out RMSE from a single random 20% split"
    )
  ) %>%
  fmt_number(columns = c(BIC, `Delta BIC`), decimals = 0) %>%
  fmt_number(
    columns = c(`Hold-out RMSE`, `Hold-out MAE`),
    decimals = 5
  ) %>%
  fmt_number(columns = `Adj. R2`, decimals = 5) %>%
  fmt_number(columns = `RMSE gain vs linear (%)`, decimals = 3) %>%
  tab_source_note(
    source_note = paste(
      "Only two variables are available, so all results are unadjusted",
      "associations and carry no causal interpretation."
    )
  ) %>%
  tab_options(table.font.size = px(13))

print(model_gt)


# -----------------------------------------------------------------------------
# 12. Figures
# -----------------------------------------------------------------------------

# Binned means give a model-free picture of the relationship that the fitted
# curve can be judged against.
binned_means <- model_data %>%
  mutate(duration_bin = cut(pr_duration, breaks = 60)) %>%
  group_by(duration_bin) %>%
  summarise(
    mid_duration = mean(pr_duration),
    mean_cnp = mean(cnp),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  filter(n_obs >= 500)

fit_plot <- ggplot() +
  geom_point(
    data = binned_means,
    aes(x = mid_duration, y = mean_cnp),
    colour = "grey45", size = 1.6
  ) +
  geom_line(
    data = grid_prediction,
    aes(x = pr_duration, y = fitted_cnp),
    colour = "#1f78b4", linewidth = 1
  ) +
  labs(
    title = paste("Selected model:", best_model_label),
    subtitle = "Grey points are binned observed means; blue line is the fitted curve",
    x = "Review duration (days)",
    y = "Normalized citation impact (cnp)"
  ) +
  theme_minimal(base_size = 12)

# Break points are marked with a dashed rule and labelled with the estimated
# day, so the figure can be read without consulting the tables. Labels are
# staggered vertically to stay legible when two break points sit close
# together.
if (!is.null(break_point_table)) {
  label_top <- max(binned_means$mean_cnp)
  label_bottom <- min(binned_means$mean_cnp)
  label_span <- label_top - label_bottom
  
  break_labels <- break_point_table %>%
    mutate(
      label_y = label_top - label_span *
        0.06 * ((break_index - 1) %% 3),
      label_text = paste0("Break: day ", round(break_day, 0)),
      ci_text = paste0(
        "95% CI ", round(ci_lower_day, 0), "-", round(ci_upper_day, 0)
      )
    )
  
  fit_plot <- fit_plot +
    geom_vline(
      data = break_labels,
      aes(xintercept = break_day),
      linetype = "dashed", colour = "#d95f02", linewidth = 0.6
    ) +
    geom_label(
      data = break_labels,
      aes(
        x = break_day,
        y = label_y,
        label = paste0(label_text, "\n", ci_text)
      ),
      colour = "#a34700",
      fill = "white",
      label.size = 0.25,
      size = 3.1,
      lineheight = 1.05,
      hjust = -0.05
    )
  
  # The caption records which model the break points come from, because the
  # plotted curve is the BIC-selected model and that is not necessarily the
  # piecewise model shown by the dashed rules.
  fit_plot <- fit_plot +
    labs(
      caption = paste0(
        "Dashed rules: break points from ", best_piecewise$label,
        " (estimated on a 150,000-row subsample)."
      )
    )
}

bic_plot <- comparison_table %>%
  mutate(label = reorder(label, -delta_bic)) %>%
  ggplot(aes(x = delta_bic, y = label)) +
  geom_col(fill = "#4c9f70") +
  labs(
    title = "Model ranking by BIC",
    subtitle = "Difference from the lowest BIC; zero marks the selected model",
    x = "Delta BIC",
    y = NULL
  ) +
  theme_minimal(base_size = 12)


# -----------------------------------------------------------------------------
# 13. Export
# -----------------------------------------------------------------------------

gtsave(model_gt, out_path("model_comparison.html"))
gtsave(model_gt, out_path("model_comparison.rtf"))

write.csv(comparison_table, out_path("model_comparison.csv"), row.names = FALSE)
write.csv(effect_summary, out_path("effect_summary.csv"), row.names = FALSE)

if (!is.null(break_point_table)) {
  write.csv(
    break_point_table, out_path("break_points.csv"), row.names = FALSE
  )
  write.csv(
    stage_slope_table, out_path("stage_slopes.csv"), row.names = FALSE
  )
}

ggsave(out_path("selected_model_fit.png"), fit_plot,
       width = 8, height = 5.2, dpi = 300)
ggsave(out_path("model_ranking_bic.png"), bic_plot,
       width = 8, height = 5.2, dpi = 300)

saveRDS(best_fit, out_path("best_model.rds"))

message("All outputs written to ", output_dir)
