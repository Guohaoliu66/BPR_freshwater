# ============================================================
# PBR group-interaction analysis for richness, FRic, and SES-FRic
# Outputs:
# - Figure4.png saved to ../000_output_file
# ============================================================

# 1. Packages & Settings --------------------------------------
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(patchwork)
library(performance)
library(purrr)
library(readr)
library(readxl)
library(tibble)

input_file  <- file.path("..", "000_input_file", "All_data.xlsx")
output_dir  <- file.path("..", "000_output_file")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

trophic_reference   <- "Consumer"
ecosystem_reference <- "River"

use_original_complete_case_filter <- TRUE
use_random_quadratic_slope        <- FALSE


# 2. Helper Functions -----------------------------------------
z_within_dataset <- function(x) {
  x_sd <- sd(x, na.rm = TRUE)
  if (is.na(x_sd) || x_sd == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

set_reference <- function(x, reference) {
  x <- factor(x)
  if (!reference %in% levels(x)) {
    stop(paste0("Reference level '", reference, "' was not found. Available levels: ", paste(levels(x), collapse = ", ")))
  }
  relevel(x, ref = reference)
}

format_p <- function(p_value) {
  if (is.null(p_value) || is.na(p_value) || length(p_value) == 0) return("p = NA")
  if (p_value < 0.001) return("p < 0.001")
  paste0("p = ", formatC(p_value, format = "f", digits = 3))
}

extract_lrt <- function(model_test) {
  if (is.null(model_test) || nrow(model_test) < 2) {
    return(tibble(chi_square = NA_real_, degrees_freedom = NA_real_, p_value = NA_real_))
  }
  
  p_col <- grep("Pr\\(", colnames(model_test), ignore.case = TRUE, value = TRUE)
  chisq_col <- grep("Chisq|Chi", colnames(model_test), ignore.case = TRUE, value = TRUE)
  df_col <- grep("Chi Df|Df", colnames(model_test), ignore.case = TRUE, value = TRUE)
  
  p_val <- if (length(p_col) > 0) model_test[2, p_col[1]] else NA_real_
  chisq_val <- if (length(chisq_col) > 0) model_test[2, chisq_col[1]] else NA_real_
  df_val <- if (length(df_col) > 0) model_test[2, df_col[1]] else NA_real_
  
  tibble(
    chi_square      = as.numeric(chisq_val),
    degrees_freedom = as.numeric(df_val),
    p_value         = as.numeric(p_val)
  )
}

fixed_design_matrix <- function(model, newdata) {
  fixed_formula <- nobars(formula(model))
  fixed_terms   <- delete.response(terms(fixed_formula))
  design_matrix <- model.matrix(fixed_terms, newdata)
  beta_names    <- names(fixef(model))
  design_matrix[, beta_names, drop = FALSE]
}

contrast_statistics <- function(contrast, beta, variance_matrix) {
  estimate       <- drop(contrast %*% beta)
  standard_error <- sqrt(drop(contrast %*% variance_matrix %*% contrast))
  p_value        <- 2 * pnorm(abs(estimate / standard_error), lower.tail = FALSE)
  
  c(
    estimate        = estimate,
    standard_error  = standard_error,
    confidence_low  = estimate - 1.96 * standard_error,
    confidence_high = estimate + 1.96 * standard_error,
    p_value         = p_value
  )
}


# 3. Data Preparation -----------------------------------------
dat_raw <- read_xlsx(input_file) %>%
  mutate(
    Richness      = as.numeric(Richness),
    FRic          = as.numeric(FRic),
    SES_FRic      = as.numeric(SES_FRic),
    TP            = as.numeric(TP),
    bio1_1500     = as.numeric(bio1_1500),
    bio12_1500    = as.numeric(bio12_1500),
    AL            = readr::parse_number(as.character(AL)),
    Trophic_level = trimws(as.character(Trophic_level)),
    Ecosystem     = trimws(as.character(Ecosystem)),
    Dataset       = trimws(as.character(Dataset))
  )

if (any(dat_raw$TP < 0, na.rm = TRUE)) {
  stop("Negative TP values were detected.")
}

dat_core <- dat_raw %>%
  mutate(logTP = log(TP + 1)) %>%
  filter(
    !is.na(Richness), !is.na(FRic), !is.na(SES_FRic),
    !is.na(TP), !is.na(logTP),
    !is.na(Trophic_level), Trophic_level != "",
    !is.na(Ecosystem), Ecosystem != "",
    !is.na(Dataset), Dataset != ""
  )

if (use_original_complete_case_filter) {
  dat_core <- dat_core %>%
    filter(
      !is.na(bio1_1500), !is.na(bio12_1500), !is.na(bio4_1500),
      !is.na(bio15_1500), !is.na(Elevation), !is.na(AL)
    )
}

logTP_center <- mean(dat_core$logTP, na.rm = TRUE)

dat <- dat_core %>%
  group_by(Dataset) %>%
  mutate(
    Richness_z = z_within_dataset(Richness),
    FRic_z     = z_within_dataset(FRic)
  ) %>%
  ungroup() %>%
  mutate(
    Dataset       = factor(Dataset),
    Trophic_level = set_reference(Trophic_level, trophic_reference),
    Ecosystem     = set_reference(Ecosystem, ecosystem_reference),
    logTP_c       = logTP - logTP_center
  )


# 4. Model Fitting --------------------------------------------
random_effect_string <- if (use_random_quadratic_slope) {
  "(1 + logTP_c + I(logTP_c^2) || Dataset)"
} else {
  "(1 + logTP_c || Dataset)"
}

fit_group_models <- function(data, response, group_var, random_effect, analysis_name) {
  model_data <- data %>%
    filter(!is.na(.data[[response]]), !is.na(logTP_c), !is.na(.data[[group_var]]), !is.na(Dataset)) %>%
    droplevels()
  
  base_formula      <- as.formula(paste0(response, " ~ logTP_c + I(logTP_c^2) + ", group_var, " + ", random_effect))
  linear_formula    <- as.formula(paste0(response, " ~ logTP_c + I(logTP_c^2) + ", group_var, " + logTP_c:", group_var, " + ", random_effect))
  quadratic_formula <- as.formula(paste0(response, " ~ logTP_c + I(logTP_c^2) + ", group_var, " + I(logTP_c^2):", group_var, " + ", random_effect))
  full_formula      <- as.formula(paste0(response, " ~ (logTP_c + I(logTP_c^2)) * ", group_var, " + ", random_effect))
  
  model_control <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  
  base_model      <- lmer(base_formula,      data = model_data, REML = FALSE, control = model_control)
  linear_model    <- lmer(linear_formula,    data = model_data, REML = FALSE, control = model_control)
  quadratic_model <- lmer(quadratic_formula, data = model_data, REML = FALSE, control = model_control)
  full_model      <- lmer(full_formula,      data = model_data, REML = FALSE, control = model_control)
  
  list(
    analysis_name   = analysis_name,
    response        = response,
    group_var       = group_var,
    data            = model_data,
    base_model      = base_model,
    linear_model    = linear_model,
    quadratic_model = quadratic_model,
    full_model      = full_model,
    overall_test    = anova(base_model, full_model),
    linear_test     = anova(quadratic_model, full_model),
    quadratic_test  = anova(linear_model, full_model)
  )
}

fits <- list(
  richness_trophic   = fit_group_models(dat, "Richness_z", "Trophic_level", random_effect_string, "richness_trophic"),
  fric_trophic       = fit_group_models(dat, "FRic_z",     "Trophic_level", random_effect_string, "fric_trophic"),
  ses_trophic        = fit_group_models(dat, "SES_FRic",   "Trophic_level", random_effect_string, "ses_trophic"),
  richness_ecosystem = fit_group_models(dat, "Richness_z", "Ecosystem",     random_effect_string, "richness_ecosystem"),
  fric_ecosystem     = fit_group_models(dat, "FRic_z",     "Ecosystem",     random_effect_string, "fric_ecosystem"),
  ses_ecosystem      = fit_group_models(dat, "SES_FRic",   "Ecosystem",     random_effect_string, "ses_ecosystem")
)


# 5. In-Memory Calculation of Annotation Metrics --------------
model_comparison_row <- function(fit) {
  linear_result    <- extract_lrt(fit$linear_test)
  quadratic_result <- extract_lrt(fit$quadratic_test)
  
  r2_val <- tryCatch({
    r2_res <- performance::r2_nakagawa(fit$full_model)
    list(m = as.numeric(r2_res$R2_marginal), c = as.numeric(r2_res$R2_conditional))
  }, error = function(e) list(m = NA_real_, c = NA_real_))
  
  tibble(
    analysis               = fit$analysis_name,
    conditional_R2         = r2_val$c,
    linear_interaction_p   = linear_result$p_value,
    quadratic_interaction_p = quadratic_result$p_value
  )
}

group_curve_table <- function(fit) {
  data         <- fit$data
  group_var    <- fit$group_var
  group_levels <- levels(data[[group_var]])
  
  map_dfr(group_levels, function(current_group) {
    group_rows <- as.character(data[[group_var]]) == current_group
    group_data <- data[group_rows, ]
    
    tibble(
      analysis       = fit$analysis_name,
      group          = current_group,
      n_observations = nrow(group_data),
      n_datasets     = n_distinct(group_data$Dataset)
    )
  })
}

table_model_comparison <- map_dfr(fits, model_comparison_row)
table_group_curves     <- map_dfr(fits, group_curve_table)


# 6. Prediction and Plotting Functions ------------------------
fixed_prediction_ci <- function(model, newdata, level = 0.95) {
  design_matrix   <- fixed_design_matrix(model, newdata)
  beta            <- fixef(model)
  variance_matrix <- as.matrix(vcov(model))[names(beta), names(beta), drop = FALSE]
  
  predicted      <- drop(design_matrix %*% beta)
  standard_error <- sqrt(rowSums((design_matrix %*% variance_matrix) * design_matrix))
  critical_value <- qnorm(1 - (1 - level) / 2)
  
  bind_cols(
    newdata,
    tibble(
      predicted = predicted,
      conf.low  = predicted - critical_value * standard_error,
      conf.high = predicted + critical_value * standard_error
    )
  )
}

make_interaction_plot <- function(
    fit,
    y_label,
    panel_title,
    logTP_center,
    model_table,
    curves_table,
    add_zero_line = FALSE) {
  model        <- fit$full_model
  data         <- fit$data
  response     <- fit$response
  group_var    <- fit$group_var
  group_levels <- levels(data[[group_var]])
  
  prediction_grid <- map_dfr(group_levels, function(current_group) {
    group_rows <- as.character(data[[group_var]]) == current_group
    group_data <- data[group_rows, ]
    
    newdata <- data.frame(
      logTP_c = seq(min(group_data$logTP_c, na.rm = TRUE), max(group_data$logTP_c, na.rm = TRUE), length.out = 200)
    )
    newdata[[group_var]] <- factor(rep(current_group, nrow(newdata)), levels = group_levels)
    newdata
  })
  
  predictions <- fixed_prediction_ci(model, prediction_grid) %>%
    mutate(logTP = logTP_c + logTP_center)
  
  m_info <- model_table %>% filter(analysis == fit$analysis_name)
  c_info <- curves_table %>% filter(analysis == fit$analysis_name)
  
  r2c_str <- if (nrow(m_info) > 0 && !is.na(m_info$conditional_R2)) {
    formatC(m_info$conditional_R2, digits = 3, format = "f")
  } else "N/A"
  
  p_lin <- if (nrow(m_info) > 0) format_p(m_info$linear_interaction_p) else "p = NA"
  p_qua <- if (nrow(m_info) > 0) format_p(m_info$quadratic_interaction_p) else "p = NA"
  
  group_lines <- map_chr(group_levels, function(g) {
    g_row <- c_info %>% filter(group == g)
    n_val <- if (nrow(g_row) > 0) format(g_row$n_observations, big.mark = ",") else "N/A"
    k_val <- if (nrow(g_row) > 0) as.character(g_row$n_datasets) else "N/A"
    sprintf("%s: n = %s, K = %s", g, n_val, k_val)
  })
  
  annotation_label <- paste0(
    paste(group_lines, collapse = "\n"),
    "\nTP × group ", p_lin,
    "\nTP² × group ", p_qua,
    "\nR² = ", r2c_str
  )
  
  x_range      <- range(data$logTP, na.rm = TRUE)
  annotation_x <- x_range[2] - 0.02 * diff(x_range)
  
  group_colors <- setNames(
    colorRampPalette(c("#2166AC", "#B2182B"))(length(group_levels)),
    group_levels
  )
  
  zero_line_layer <- if (add_zero_line) {
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.6
    )
  } else {
    NULL
  }
  
  ggplot() +
    zero_line_layer +
    geom_point(
      data = data,
      aes(x = logTP, y = .data[[response]], color = .data[[group_var]]),
      alpha = 0.035, size = 0.45
    ) +
    geom_ribbon(
      data = predictions,
      aes(x = logTP, ymin = conf.low, ymax = conf.high, fill = .data[[group_var]], group = .data[[group_var]]),
      alpha = 0.15, color = NA
    ) +
    geom_line(
      data = predictions,
      aes(x = logTP, y = predicted, color = .data[[group_var]], group = .data[[group_var]]),
      linewidth = 1.25
    ) +
    annotate(
      "text", x = annotation_x, y = Inf, label = annotation_label,
      hjust = 1, vjust = 1.15, size = 3.3, lineheight = 1.1
    ) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    coord_cartesian(ylim = c(-2.5, 3.5)) +
    scale_y_continuous(breaks = seq(-2, 3, by = 1)) +
    labs(
      title = panel_title,
      x     = "log(TP + 1)",
      y     = y_label,
      color = NULL,
      fill  = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title      = element_text(size = 14)
    )
}


# 7. Generate Panels ------------------------------------------
p_trophic_rich <- make_interaction_plot(
  fits$richness_trophic,
  "Z-Richness",
  "Richness: TP × trophic level",
  logTP_center,
  table_model_comparison,
  table_group_curves
)

p_trophic_fric <- make_interaction_plot(
  fits$fric_trophic,
  "Z-FRic",
  "FRic: TP × trophic level",
  logTP_center,
  table_model_comparison,
  table_group_curves
)

p_trophic_ses <- make_interaction_plot(
  fits$ses_trophic,
  "SES-FRic",
  "SES-FRic: TP × trophic level",
  logTP_center,
  table_model_comparison,
  table_group_curves,
  add_zero_line = TRUE
)

p_ecosystem_rich <- make_interaction_plot(
  fits$richness_ecosystem,
  "Z-Richness",
  "Richness: TP × ecosystem",
  logTP_center,
  table_model_comparison,
  table_group_curves
)

p_ecosystem_fric <- make_interaction_plot(
  fits$fric_ecosystem,
  "Z-FRic",
  "FRic: TP × ecosystem",
  logTP_center,
  table_model_comparison,
  table_group_curves
)

p_ecosystem_ses <- make_interaction_plot(
  fits$ses_ecosystem,
  "SES-FRic",
  "SES-FRic: TP × ecosystem",
  logTP_center,
  table_model_comparison,
  table_group_curves,
  add_zero_line = TRUE
)

final_plot <-
  (p_trophic_rich | p_trophic_fric | p_trophic_ses) /
  (p_ecosystem_rich | p_ecosystem_fric | p_ecosystem_ses) +
  plot_annotation(tag_levels = "a")


# 8. Save Plot ------------------------------------------------
ggsave(
  file.path(output_dir, "Figure4.png"),
  plot = final_plot,
  width = 12,
  height = 10,
  dpi = 600
)