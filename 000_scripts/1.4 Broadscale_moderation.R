

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(lme4)
  library(lmerTest)
  library(performance)
  library(ggplot2)
  library(patchwork)
  library(openxlsx)
  library(tibble)
})

get_script_directory <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(p)) return(dirname(p))
  }
  getwd()
}

setwd(get_script_directory())

input_file <- file.path("..", "000_input_file", "All_data.xlsx")
output_dir <- file.path("..", "000_output_file")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_file)) {
  stop("Input file not found: ", normalizePath(input_file, mustWork = FALSE))
}

response_metadata <- tribble(
  ~response,     ~response_label,
  "Richness_z",  "Taxonomic richness",
  "FRic_z",      "Functional richness"
)

moderator_metadata <- tribble(
  ~moderator,                      ~moderator_label,              ~moderator_type,
  "Temperature_z",                 "Annual mean temperature",     "continuous",
  "Temperature_seasonality_z",     "Temperature seasonality",     "continuous",
  "Precipitation_z",                "Annual precipitation",        "continuous",
  "Precipitation_seasonality_z",   "Precipitation seasonality",   "continuous",
  "Elevation_z",                    "Elevation",                   "continuous",
  "Landuse_class",                 "Land use",                    "factor"
)

continuous_curve_values <- c(-1, 0, 1)
continuous_curve_labels <- c("Low (-1 SD)", "Mean (0)", "High (+1 SD)")
continuous_palette <- c("Low (-1 SD)" = "#2166AC", "Mean (0)" = "#4D4D4D", "High (+1 SD)" = "#B2182B")
response_colors <- c("Functional richness" = "#C51B32", "Taxonomic richness" = "#2878B5")
random_effect_text <- "(1 + logTP_c || Dataset)"
fit_control <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

z_within_group <- function(x) {
  x <- as.numeric(x)
  keep <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(keep)) return(out)
  s <- stats::sd(x[keep])
  out[keep] <- if (!is.finite(s) || s == 0) 0 else as.numeric(scale(x[keep]))
  out
}

z_global <- function(x) {
  x <- as.numeric(x)
  keep <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(keep)) return(out)
  s <- stats::sd(x[keep])
  if (!is.finite(s) || s == 0) stop("Continuous moderator has 0 variance.")
  out[keep] <- as.numeric(scale(x[keep]))
  out
}

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

safe_r2 <- function(model) {
  res <- tryCatch(
    suppressWarnings(performance::r2_nakagawa(model)),
    error = function(e) NULL
  )
  
  m_r2 <- NA_real_
  c_r2 <- NA_real_
  
  if (!is.null(res)) {
    m_r2 <- unname(res$R2_marginal)
    c_r2 <- unname(res$R2_conditional)
  }
  
  if (is.na(c_r2) && requireNamespace("MuMIn", quietly = TRUE)) {
    res_mumin <- tryCatch(
      suppressWarnings(MuMIn::r.squaredGLMM(model)),
      error = function(e) NULL
    )
    if (!is.null(res_mumin)) {
      m_r2 <- unname(res_mumin[1, "R2m"])
      c_r2 <- unname(res_mumin[1, "R2c"])
    }
  }
  
  if (is.na(c_r2) && is.finite(m_r2) && lme4::isSingular(model, tol = 1e-4)) {
    c_r2 <- m_r2
  }
  
  c(marginal = m_r2, conditional = c_r2)
}

extract_lrt <- function(m1, m2) {
  comp <- suppressMessages(anova(m1, m2, test = "Chisq"))
  df_col <- intersect(c("Df", "Chi Df"), names(comp))[1]
  tibble(chisq = comp$Chisq[2], df = comp[[df_col]][2], p_value = comp$`Pr(>Chisq)`[2])
}

fixed_effect_table <- function(model) {
  tab <- as.data.frame(summary(model)$coefficients)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  crit <- if ("df" %in% names(tab)) stats::qt(0.975, df = tab$df) else rep(stats::qnorm(0.975), nrow(tab))
  p_vals <- if ("Pr(>|t|)" %in% names(tab)) tab[["Pr(>|t|)"]] else rep(NA_real_, nrow(tab))
  
  tibble(
    term = tab$term,
    estimate = tab$Estimate,
    standard_error = tab$`Std. Error`,
    degrees_freedom = tab$df,
    statistic = tab$`t value`,
    p_value = p_vals,
    confidence_low = tab$Estimate - crit * tab$`Std. Error`,
    confidence_high = tab$Estimate + crit * tab$`Std. Error`
  )
}

fixed_prediction_ci <- function(model, newdata, confidence_level = 0.95) {
  fixed_terms <- stats::delete.response(stats::terms(lme4::nobars(stats::formula(model))))
  dm <- stats::model.matrix(fixed_terms, data = newdata)
  beta <- lme4::fixef(model)
  beta_names <- names(beta)
  
  dm <- dm[, beta_names, drop = FALSE]
  vcov_mat <- as.matrix(stats::vcov(model))[beta_names, beta_names, drop = FALSE]
  
  pred <- as.numeric(dm %*% beta)
  pred_se <- sqrt(pmax(rowSums((dm %*% vcov_mat) * dm), 0))
  crit <- stats::qnorm(1 - (1 - confidence_level) / 2)
  
  newdata %>%
    mutate(
      predicted = pred,
      standard_error = pred_se,
      conf_low = pred - crit * pred_se,
      conf_high = pred + crit * pred_se
    )
}

draw_key_ci_line <- function(data, params, size) {
  col <- if (length(data$colour) == 0L || is.na(data$colour)) "black" else data$colour
  grid::grobTree(
    grid::rectGrob(width = grid::unit(0.82, "npc"), height = grid::unit(0.82, "npc"),
                   gp = grid::gpar(col = NA, fill = scales::alpha(col, 0.16))),
    grid::segmentsGrob(x0 = grid::unit(0.14, "npc"), x1 = grid::unit(0.86, "npc"),
                       y0 = grid::unit(0.50, "npc"), y1 = grid::unit(0.50, "npc"),
                       gp = grid::gpar(col = col, lwd = 2.4, lineend = "butt"))
  )
}

dat_raw <- read_xlsx(input_file)

dat <- dat_raw %>%
  mutate(
    across(c(Richness, FRic, TP, bio1_1500, bio4_1500, bio12_1500, bio15_1500, Elevation), as.numeric),
    logTP = log(TP + 1),
    Dataset = factor(Dataset),
    Landuse_class = factor(Landuse_class)
  ) %>%
  filter(!is.na(Dataset), is.finite(TP), is.finite(logTP)) %>%
  group_by(Dataset) %>%
  mutate(
    Richness_z = z_within_group(Richness),
    FRic_z = z_within_group(FRic)
  ) %>%
  ungroup() %>%
  mutate(
    logTP_c = logTP - mean(logTP, na.rm = TRUE),
    Temperature_z = z_global(bio1_1500),
    Temperature_seasonality_z = z_global(bio4_1500),
    Precipitation_z = z_global(bio12_1500),
    Precipitation_seasonality_z = z_global(bio15_1500),
    Elevation_z = z_global(Elevation)
  )

if ("Forest" %in% levels(dat$Landuse_class)) {
  dat$Landuse_class <- stats::relevel(dat$Landuse_class, ref = "Forest")
}

global_logTP_mean <- mean(dat$logTP, na.rm = TRUE)

Richness_best <- lmer(
  Richness_z ~ logTP + I(logTP^2) + Landuse_class + bio1_1500 + bio12_1500 + bio15_1500 + Elevation + (1 + logTP || Dataset),
  REML = FALSE, data = dat, na.action = stats::na.omit, control = fit_control
)

FRic_best <- lmer(
  FRic_z ~ logTP + I(logTP^2) + Landuse_class + bio12_1500 + bio4_1500 + bio15_1500 + Elevation + (1 + logTP || Dataset),
  REML = FALSE, data = dat, na.action = stats::na.omit, control = fit_control
)

extract_additive_model <- function(model, model_name, response_label) {
  r2 <- safe_r2(model)
  mf <- stats::model.frame(model)
  aic_val <- stats::AIC(model)
  is_sing <- lme4::isSingular(model, tol = 1e-4)
  
  fixed_effect_table(model) %>%
    mutate(
      model = model_name,
      response = response_label,
      n_observations = nrow(mf),
      n_datasets = dplyr::n_distinct(mf$Dataset),
      marginal_R2 = r2["marginal"],
      conditional_R2 = r2["conditional"],
      AIC = aic_val,
      singular_model = is_sing,
      .before = 1
    )
}

additive_model_table <- bind_rows(
  extract_additive_model(Richness_best, "Richness_best", "Taxonomic richness"),
  extract_additive_model(FRic_best, "FRic_best", "Functional richness")
)

fit_moderator_models <- function(data, response, response_label, moderator, moderator_label, moderator_type) {
  analysis_data <- data %>%
    select(all_of(c(response, moderator, "logTP", "logTP_c", "Dataset"))) %>%
    drop_na() %>%
    droplevels()
  
  if (n_distinct(analysis_data$Dataset) < 2L) stop("Fewer than 2 datasets for ", response_label, " × ", moderator_label)
  if (moderator_type == "factor" && nlevels(analysis_data[[moderator]]) < 2L) stop("Fewer than 2 factor levels for ", moderator_label)
  
  f_base <- as.formula(paste0(response, " ~ logTP_c + I(logTP_c^2) + ", moderator, " + ", random_effect_text))
  f_lin  <- as.formula(paste0(response, " ~ logTP_c + I(logTP_c^2) + ", moderator, " + logTP_c:", moderator, " + ", random_effect_text))
  f_quad <- as.formula(paste0(response, " ~ logTP_c + I(logTP_c^2) + ", moderator, " + I(logTP_c^2):", moderator, " + ", random_effect_text))
  f_full <- as.formula(paste0(response, " ~ (logTP_c + I(logTP_c^2)) * ", moderator, " + ", random_effect_text))
  
  m_base <- lmer(f_base, data = analysis_data, REML = FALSE, control = fit_control)
  m_lin  <- lmer(f_lin,  data = analysis_data, REML = FALSE, control = fit_control)
  m_quad <- lmer(f_quad, data = analysis_data, REML = FALSE, control = fit_control)
  m_full <- lmer(f_full, data = analysis_data, REML = FALSE, control = fit_control)
  
  t_overall <- extract_lrt(m_base, m_full)
  t_lin     <- extract_lrt(m_quad, m_full)
  t_quad    <- extract_lrt(m_lin,  m_full)
  r2        <- safe_r2(m_full)
  analysis_id <- paste(response, moderator, sep = "__")
  
  model_test <- tibble(
    analysis = analysis_id,
    response = response,
    response_label = response_label,
    moderator = moderator,
    moderator_label = moderator_label,
    moderator_type = moderator_type,
    n_observations = nrow(stats::model.frame(m_full)),
    n_datasets = n_distinct(analysis_data$Dataset),
    marginal_R2 = r2["marginal"],
    conditional_R2 = r2["conditional"],
    base_AIC = stats::AIC(m_base),
    full_AIC = stats::AIC(m_full),
    delta_AIC_full_minus_base = stats::AIC(m_full) - stats::AIC(m_base),
    overall_interaction_chisq = t_overall$chisq,
    overall_interaction_df = t_overall$df,
    overall_interaction_p = t_overall$p_value,
    linear_interaction_chisq = t_lin$chisq,
    linear_interaction_df = t_lin$df,
    linear_interaction_p = t_lin$p_value,
    quadratic_interaction_chisq = t_quad$chisq,
    quadratic_interaction_df = t_quad$df,
    quadratic_interaction_p = t_quad$p_value,
    singular_full_model = lme4::isSingular(m_full, tol = 1e-4)
  )
  
  interaction_coefficients <- fixed_effect_table(m_full) %>%
    filter(grepl(":", term, fixed = TRUE)) %>%
    mutate(
      component = if_else(grepl("I\\(logTP_c\\^2\\)", term), "Quadratic interaction", "Linear interaction"),
      analysis = analysis_id,
      response = response,
      response_label = response_label,
      moderator = moderator,
      moderator_label = moderator_label,
      moderator_type = moderator_type,
      .before = 1
    )
  
  list(
    analysis = analysis_id,
    data = analysis_data,
    base_model = m_base,
    linear_interaction_model = m_lin,
    quadratic_interaction_model = m_quad,
    full_model = m_full,
    model_test = model_test,
    interaction_coefficients = interaction_coefficients
  )
}

moderator_combos <- tidyr::crossing(response_metadata, moderator_metadata)
moderator_fits <- pmap(moderator_combos, function(response, response_label, moderator, moderator_label, moderator_type) {
  message("Fitting: ", response_label, " × ", moderator_label)
  fit_moderator_models(dat, response, response_label, moderator, moderator_label, moderator_type)
})
names(moderator_fits) <- map_chr(moderator_fits, "analysis")

model_test_table <- map_dfr(moderator_fits, "model_test") %>%
  group_by(response) %>%
  mutate(
    overall_interaction_p_FDR = p.adjust(overall_interaction_p, method = "BH"),
    FDR_supported = overall_interaction_p_FDR < 0.05
  ) %>%
  ungroup() %>%
  arrange(match(response, response_metadata$response), overall_interaction_p_FDR)

interaction_coefficient_table <- map_dfr(moderator_fits, "interaction_coefficients") %>%
  left_join(
    select(model_test_table, analysis, overall_interaction_p, overall_interaction_p_FDR, FDR_supported),
    by = "analysis"
  ) %>%
  arrange(match(response, response_metadata$response), moderator_label, component, term)

fdr_supported_table <- model_test_table %>%
  filter(FDR_supported) %>%
  arrange(match(response, response_metadata$response), overall_interaction_p_FDR)

Table_S2 <- model_test_table %>%
  transmute(
    Response = response_label,
    Moderator = moderator_label,
    N = n_observations,
    K = n_datasets,
    `Marginal R2` = round(marginal_R2, 3),
    `Conditional R2` = round(conditional_R2, 3),
    `AIC base` = round(base_AIC, 1),
    `AIC full` = round(full_AIC, 1),
    `Delta AIC (full - base)` = round(delta_AIC_full_minus_base, 1),
    `LRT chi-square` = round(overall_interaction_chisq, 2),
    df = overall_interaction_df,
    P = vapply(overall_interaction_p, format_p, character(1)),
    `FDR-adjusted P` = vapply(overall_interaction_p_FDR, format_p, character(1))
  )

make_prediction_data <- function(fit_object, moderator_type) {
  model <- fit_object$full_model
  analysis_data <- fit_object$data
  test_row <- fit_object$model_test
  moderator <- test_row$moderator
  
  if (moderator_type == "continuous") {
    prediction_grid <- tidyr::expand_grid(
      logTP_c = seq(min(analysis_data$logTP_c, na.rm = TRUE), max(analysis_data$logTP_c, na.rm = TRUE), length.out = 200),
      moderator_value = continuous_curve_values
    ) %>%
      mutate(curve_level = factor(continuous_curve_labels[match(moderator_value, continuous_curve_values)], levels = continuous_curve_labels))
    prediction_grid[[moderator]] <- prediction_grid$moderator_value
  } else {
    lvls <- levels(analysis_data[[moderator]])
    prediction_grid <- map_dfr(lvls, function(lvl) {
      sub_data <- filter(analysis_data, .data[[moderator]] == lvl)
      tibble(
        logTP_c = seq(min(sub_data$logTP_c, na.rm = TRUE), max(sub_data$logTP_c, na.rm = TRUE), length.out = 200),
        curve_level = factor(lvl, levels = lvls)
      )
    })
    prediction_grid[[moderator]] <- prediction_grid$curve_level
  }
  
  fixed_prediction_ci(model, prediction_grid) %>%
    mutate(
      logTP = logTP_c + global_logTP_mean,
      analysis = fit_object$analysis,
      response = test_row$response,
      response_label = test_row$response_label,
      moderator = test_row$moderator,
      moderator_label = test_row$moderator_label,
      moderator_type = test_row$moderator_type,
      .before = 1
    )
}

prediction_data_all <- map_dfr(fdr_supported_table$analysis, function(id) {
  fit <- moderator_fits[[id]]
  make_prediction_data(fit, fit$model_test$moderator_type)
})

make_moderation_plot <- function(fit_object, corrected_test_row) {
  p_data <- filter(prediction_data_all, analysis == corrected_test_row$analysis)
  r2_txt <- if (is.finite(corrected_test_row$marginal_R2)) sprintf("%.3f", corrected_test_row$marginal_R2) else "NA"
  annot_txt <- paste0("Marginal R² = ", r2_txt, "\nFDR-adjusted P = ", format_p(corrected_test_row$overall_interaction_p_FDR))
  
  if (corrected_test_row$moderator_type == "continuous") {
    line_cols <- continuous_palette
  } else {
    lvls <- levels(p_data$curve_level)
    line_cols <- stats::setNames(grDevices::colorRampPalette(c("#2166AC", "#4D4D4D", "#B2182B"))(length(lvls)), lvls)
  }
  
  ggplot() +
    geom_ribbon(data = p_data, aes(x = logTP, ymin = conf_low, ymax = conf_high, fill = curve_level, group = curve_level),
                alpha = 0.16, colour = NA, show.legend = FALSE) +
    geom_line(data = p_data, aes(x = logTP, y = predicted, colour = curve_level, group = curve_level),
              linewidth = 1.25, lineend = "round", key_glyph = draw_key_ci_line) +
    annotate("text", x = Inf, y = Inf, label = annot_txt, hjust = 1.04, vjust = 1.12, size = 3.4, lineheight = 1.08) +
    scale_colour_manual(values = line_cols, name = NULL, drop = FALSE) +
    scale_fill_manual(values = line_cols, name = NULL, drop = FALSE) +
    coord_cartesian(ylim = c(-2, 2)) +
    scale_y_continuous(breaks = seq(-2, 2, by = 1), expand = expansion(mult = c(0.02, 0.03))) +
    labs(title = paste0("TP × ", corrected_test_row$moderator_label), x = "log(TP + 1)",
         y = if (corrected_test_row$response == "Richness_z") "Z-Richness" else "Z-FRic") +
    guides(fill = "none", colour = guide_legend(title = NULL, nrow = 1, byrow = TRUE,
                                                keywidth = grid::unit(0.68, "cm"), keyheight = grid::unit(0.55, "cm"))) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 12.5, hjust = 0.5, margin = margin(b = 5)),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10.5, colour = "black"),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.text = element_text(size = 9.5, face = "bold"),
      legend.key = element_blank(),
      plot.margin = margin(8, 8, 8, 8)
    )
}

prediction_plot_list <- map(seq_len(nrow(fdr_supported_table)), function(i) {
  row <- fdr_supported_table[i, ]
  make_moderation_plot(moderator_fits[[row$analysis]], row)
})
names(prediction_plot_list) <- fdr_supported_table$analysis

if (length(prediction_plot_list) > 0L) {
  prediction_figure <- patchwork::wrap_plots(prediction_plot_list, ncol = 2) +
    patchwork::plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 16, face = "bold"))
  
  print(prediction_figure)
  
  pred_rows <- ceiling(length(prediction_plot_list) / 2)
  ggsave(
    filename = file.path(output_dir, "FigureS3.png"),
    plot = prediction_figure,
    width = 9, height = max(5, 3.6 * pred_rows), dpi = 600, bg = "white", limitsize = FALSE
  )
}

effect_order <- c(
  "Annual mean temperature", "Temperature seasonality", "Annual precipitation",
  "Precipitation seasonality", "Elevation", "Land use: Mixed vs Forest", "Land use: Human vs Forest"
)
effect_levels <- rev(effect_order)

forest_data <- interaction_coefficient_table %>%
  mutate(
    component_clean = factor(
      if_else(grepl("Quadratic", component, ignore.case = TRUE), "Quadratic interaction", "Linear interaction"),
      levels = c("Linear interaction", "Quadratic interaction")
    ),
    landuse_level = if_else(moderator == "Landuse_class", sub(":.*$", "", sub(".*Landuse_class", "", term)), NA_character_),
    effect_label = case_when(
      moderator != "Landuse_class" ~ moderator_label,
      grepl("Mixed", landuse_level, ignore.case = TRUE) ~ "Land use: Mixed vs Forest",
      grepl("Human", landuse_level, ignore.case = TRUE) ~ "Land use: Human vs Forest",
      TRUE ~ paste0("Land use: ", landuse_level, " vs Forest")
    ),
    response_label = factor(response_label, levels = c("Functional richness", "Taxonomic richness")),
    FDR_status = factor(if_else(FDR_supported, "P_adj < 0.05", "P_adj >= 0.05"), levels = c("P_adj >= 0.05", "P_adj < 0.05")),
    effect_label = factor(effect_label, levels = effect_levels)
  ) %>%
  filter(!is.na(effect_label), !is.na(component_clean), !is.na(response_label), !is.na(FDR_status),
         is.finite(estimate), is.finite(confidence_low), is.finite(confidence_high)) %>%
  mutate(
    y_base = as.numeric(effect_label),
    y_offset = case_when(
      response_label == "Functional richness" ~ 0.13,
      response_label == "Taxonomic richness" ~ -0.13,
      TRUE ~ 0
    ),
    y_position = y_base + y_offset
  )

interaction_forest <- ggplot(forest_data, aes(colour = response_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, colour = "grey55") +
  geom_segment(aes(x = confidence_low, xend = confidence_high, y = y_position, yend = y_position), linewidth = 0.75) +
  geom_point(aes(x = estimate, y = y_position, shape = FDR_status), size = 3.1, stroke = 1.1) +
  facet_wrap(~ component_clean, nrow = 1, scales = "free_x") +
  scale_y_continuous(breaks = seq_along(effect_levels), labels = effect_levels, expand = expansion(add = c(0.55, 0.55))) +
  scale_colour_manual(values = response_colors, drop = FALSE) +
  scale_shape_manual(values = c("P_adj >= 0.05" = 1, "P_adj < 0.05" = 16),
                     labels = c(expression(italic(P)[adj] >= 0.05), expression(italic(P)[adj] < 0.05)), drop = FALSE) +
  labs(x = "Interaction coefficient (95% CI)", y = NULL, colour = NULL, shape = "Overall interaction test") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 13, colour = "black"),
    axis.text.y = element_text(size = 11, colour = "black"),
    axis.text.x = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 12, margin = margin(t = 8)),
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.justification = "center",
    legend.margin = margin(l = 8),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    panel.spacing = grid::unit(1.4, "cm"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 40)
  ) +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2, title.position = "top"))

print(interaction_forest)

ggsave(
  filename = file.path(output_dir, "Figure6.png"),
  plot = interaction_forest,
  width = 10, height = 5, units = "in", dpi = 600, bg = "white", limitsize = FALSE
)

openxlsx::write.xlsx(
  list(
    Final_additive_models = additive_model_table,
    Interaction_tests = model_test_table,
    Interaction_coefficients = interaction_coefficient_table,
    FDR_supported_models = fdr_supported_table,
    Predictions_FDR_supported = prediction_data_all,
    Table_S2 = Table_S2
  ),
  file = file.path(output_dir, "Supplementary_broadscale_result.xlsx"),
  overwrite = TRUE
)

openxlsx::write.xlsx(
  Table_S2,
  file = file.path(output_dir, "Table S2.xlsx"),
  overwrite = TRUE
)
