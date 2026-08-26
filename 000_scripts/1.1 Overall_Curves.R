# ==============================================================================
# Overall TP relationships for taxonomic richness, FRic, and SES-FRic
# ==============================================================================

library(readxl)
library(dplyr)
library(readr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(tidyr)
library(purrr)
library(ggeffects)
library(patchwork)

# ==============================================================================
# 1. File paths and helper functions
# ==============================================================================
get_script_directory <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    current_path <- rstudioapi::getActiveDocumentContext()$path
    
    if (nzchar(current_path)) {
      return(dirname(current_path))
    }
  }
  
  getwd()
}

script_dir <- get_script_directory()
setwd(script_dir)

input_file <- file.path("..", "000_input_file", "All_data.xlsx")
output_dir <- file.path("..", "000_output_file")

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

z_within_group <- function(x) {
  standard_deviation <- sd(x, na.rm = TRUE)

  if (is.na(standard_deviation) || standard_deviation == 0) {
    return(rep(0, length(x)))
  }

  as.numeric(scale(x))
}

# ==============================================================================
# 2. Read and prepare the data
# ==============================================================================

dat <- read_xlsx(input_file) %>%
  mutate(
    Richness = as.numeric(Richness),
    FRic = as.numeric(FRic),
    SES_FRic = as.numeric(SES_FRic),
    TP = as.numeric(TP),
    logTP = log(TP + 1),
    bio1_1500 = as.numeric(bio1_1500),
    bio12_1500 = as.numeric(bio12_1500),
    AL = readr::parse_number(as.character(AL))
  ) %>%
  filter(
    !is.na(Richness),
    !is.na(FRic),
    !is.na(SES_FRic),
    !is.na(TP),
    !is.na(logTP),
    !is.na(Trophic_level),
    !is.na(Ecosystem),
    !is.na(Dataset),
    !is.na(bio1_1500),
    !is.na(bio12_1500),
    !is.na(bio4_1500),
    !is.na(bio15_1500),
    !is.na(Elevation),
    !is.na(AL)
  ) %>%
  group_by(Dataset) %>%
  mutate(
    Richness_z = z_within_group(Richness),
    FRic_z = z_within_group(FRic)
  ) %>%
  ungroup() %>%
  mutate(
    Dataset = as.factor(Dataset),
    Trophic_level = as.factor(Trophic_level),
    Ecosystem = as.factor(Ecosystem)
  ) %>%
  droplevels()

# Richness and observed FRic are standardized within each dataset.
# SES-FRic is not standardized again because it is already a standardized
# effect size derived from the richness-constrained null model.

# ==============================================================================
# 3. Fit linear and quadratic mixed-effects models
# ==============================================================================

# Taxonomic richness
mod_rich_lin <- lmer(
  Richness_z ~ logTP +
    (1 + logTP || Dataset),
  data = dat,
  REML = FALSE
)

mod_rich_quad <- lmer(
  Richness_z ~ logTP + I(logTP^2) +
    (1 + logTP || Dataset),
  data = dat,
  REML = FALSE
)

# Observed functional richness
mod_fric_lin <- lmer(
  FRic_z ~ logTP +
    (1 + logTP || Dataset),
  data = dat,
  REML = FALSE
)

mod_fric_quad <- lmer(
  FRic_z ~ logTP + I(logTP^2) +
    (1 + logTP || Dataset),
  data = dat,
  REML = FALSE
)

# Standardized effect size of functional richness
mod_ses_lin <- lmer(
  SES_FRic ~ logTP +
    (1 + logTP || Dataset),
  data = dat,
  REML = FALSE
)

mod_ses_quad <- lmer(
  SES_FRic ~ logTP + I(logTP^2) +
    (1 + logTP || Dataset),
  data = dat,
  REML = FALSE
)


extract_quad_coefs <- function(mod) {
  cf <- summary(mod)$coefficients
  terms <- c("logTP", "I(logTP^2)")
  cf[terms, c("Estimate", "Std. Error")]
}

results <- rbind(
  Richness = extract_quad_coefs(mod_rich_quad),
  FRic     = extract_quad_coefs(mod_fric_quad),
  SES_FRic = extract_quad_coefs(mod_ses_quad)
)

print(results)


# Compare the linear and quadratic models
print(anova(mod_rich_lin, mod_rich_quad))
print(anova(mod_fric_lin, mod_fric_quad))
print(anova(mod_ses_lin, mod_ses_quad))

# ==============================================================================
# 4. Format coefficient labels
# ==============================================================================

get_math_labels_list <- function(model) {
  coefficients <- summary(model)$coefficients

  beta_linear <- coefficients["logTP", "Estimate"]
  se_linear <- coefficients["logTP", "Std. Error"]
  beta_quadratic <- coefficients["I(logTP^2)", "Estimate"]
  se_quadratic <- coefficients["I(logTP^2)", "Std. Error"]

  format_coefficient <- function(value, standard_error) {
    value_scientific <- formatC(value, format = "e", digits = 1)
    value_parts <- strsplit(value_scientific, "e")[[1]]

    se_scientific <- formatC(standard_error, format = "e", digits = 1)
    se_parts <- strsplit(se_scientific, "e")[[1]]

    sprintf(
      "%.1f %%*%% 10^%d ~ (SE == %.1f %%*%% 10^%d)",
      as.numeric(value_parts[1]),
      as.integer(value_parts[2]),
      as.numeric(se_parts[1]),
      as.integer(se_parts[2])
    )
  }

  c(
    line1 = paste0(
      "italic(beta)[1]==",
      format_coefficient(beta_linear, se_linear)
    ),
    line2 = paste0(
      "italic(beta)[2]==",
      format_coefficient(beta_quadratic, se_quadratic)
    )
  )
}

# ==============================================================================
# 5. Define color palettes
# ==============================================================================

blue_palette <- function(n) {
  colorRampPalette(
    c("#9ECAE1", "#6BAED6", "#3182BD", "#08519C")
  )(n)
}

red_palette <- function(n) {
  colorRampPalette(
    c("#FEE0D2", "#FC9272", "#DE2D26", "#A50F15")
  )(n)
}

purple_palette <- function(n) {
  colorRampPalette(
    c("#DADAEB", "#9E9AC8", "#756BB1", "#54278F")
  )(n)
}

# ==============================================================================
# 6. Create a panel with dataset-specific and overall predictions
# ==============================================================================

make_gradient_plot <- function(
    data,
    response,
    y_label,
    model,
    line_palette_function,
    main_line_color,
    ribbon_fill_color = main_line_color,
    add_zero_line = FALSE,
    text_color = "black") {

  plot_data <- data %>%
    mutate(Dataset = as.factor(Dataset))

  tp_range <- plot_data %>%
    group_by(Dataset) %>%
    summarise(
      minimum_logTP = min(logTP, na.rm = TRUE),
      maximum_logTP = max(logTP, na.rm = TRUE),
      .groups = "drop"
    )

  individual_predictions <- map_dfr(
    seq_len(nrow(tp_range)),
    function(index) {
      prediction_data <- data.frame(
        Dataset = tp_range$Dataset[index],
        logTP = seq(
          tp_range$minimum_logTP[index],
          tp_range$maximum_logTP[index],
          length.out = 50
        )
      )

      prediction_data$predicted <- predict(
        model,
        newdata = prediction_data,
        re.form = NULL
      )

      prediction_data
    }
  )

  overall_predictions <- as.data.frame(
    ggeffects::ggpredict(
      model,
      terms = "logTP [all]"
    )
  )

  dataset_levels <- levels(plot_data$Dataset)
  color_values <- line_palette_function(length(dataset_levels))
  names(color_values) <- dataset_levels

  x_range <- range(plot_data$logTP, na.rm = TRUE)
  x_position <- x_range[2] - 0.05 * diff(x_range)
  label_list <- get_math_labels_list(model)

  plot_object <- ggplot() +
    geom_point(
      data = plot_data,
      aes(
        x = logTP,
        y = .data[[response]],
        color = Dataset
      ),
      alpha = 0.02,
      size = 1
    ) +

    # Draw the color-matched confidence interval first so that the
    # dataset-specific prediction lines remain visible above it.
    geom_ribbon(
      data = overall_predictions,
      aes(
        x = x,
        ymin = conf.low,
        ymax = conf.high
      ),
      inherit.aes = FALSE,
      fill = ribbon_fill_color,
      alpha = 0.12
    ) +

    # Conditional prediction for each dataset
    geom_line(
      data = individual_predictions,
      aes(
        x = logTP,
        y = predicted,
        group = Dataset,
        color = Dataset
      ),
      inherit.aes = FALSE,
      linewidth = 0.3,
      alpha = 0.55
    ) +

    # Overall fixed-effect prediction
    geom_line(
      data = overall_predictions,
      aes(
        x = x,
        y = predicted
      ),
      inherit.aes = FALSE,
      color = main_line_color,
      linewidth = 1.4
    ) +
    annotate(
      "text",
      x = x_position,
      y = 4.2,
      label = label_list["line1"],
      parse = TRUE,
      size = 4.5,
      hjust = 1,
      color = text_color
    ) +
    annotate(
      "text",
      x = x_position,
      y = 3.3,
      label = label_list["line2"],
      parse = TRUE,
      size = 4.5,
      hjust = 1,
      color = text_color
    ) +
    scale_color_manual(values = color_values) +
    coord_cartesian(ylim = c(-4, 4)) +
    scale_y_continuous(breaks = seq(-4, 4, by = 2)) +
    labs(
      x = "log(TP + 1)",
      y = y_label
    ) +
    theme_classic(base_size = 18) +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      plot.margin = grid::unit(c(1, 1, 1, 1), "lines")
    )

  if (add_zero_line) {
    plot_object <- plot_object +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "grey40",
        linewidth = 0.6
      )
  }

  plot_object
}

# ==============================================================================
# 7. Create the three panels
# ==============================================================================

p_rich <- make_gradient_plot(
  data = dat,
  response = "Richness_z",
  y_label = "Z-Richness",
  model = mod_rich_quad,
  line_palette_function = blue_palette,
  main_line_color = "#2166AC"
)

p_fric <- make_gradient_plot(
  data = dat,
  response = "FRic_z",
  y_label = "Z-FRic",
  model = mod_fric_quad,
  line_palette_function = red_palette,
  main_line_color = "#B2182B"
)

p_ses <- make_gradient_plot(
  data = dat,
  response = "SES_FRic",
  y_label = "SES-FRic",
  model = mod_ses_quad,
  line_palette_function = purple_palette,
  main_line_color = "#54278F",
  add_zero_line = TRUE
)

# Combine panels a-c
final_plot <- p_rich | p_fric | p_ses

final_plot <- final_plot +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 24))

# ==============================================================================
# 8. Display the figure
# ==============================================================================

print(final_plot)

# ==============================================================================
# 9. Save the figure
# This section can be run independently after changing the dimensions.
# ==============================================================================


ggsave(file.path(output_dir, "Figure3.png"), plot = final_plot, width = 15, height = 7, dpi = 600)





# ============================================================
# Organism-specific TP–diversity relationships
# Richness, observed FRic, and SES-FRic
# ============================================================


organism_order <- c(
  "Aquatic Plants",
  "Benthic Diatoms",
  "Phytoplankton",
  "Fish",
  "Benthic Macroinvertebrates",
  "Zooplankton"
)

organism_labels <- c(
  "Aquatic Plants" = "Aquatic Plants",
  "Benthic Diatoms" = "Diatom",
  "Phytoplankton" = "Phytoplankton",
  "Fish" = "Fish",
  "Benthic Macroinvertebrates" = "Macroinvertebrates",
  "Zooplankton" = "Zooplankton"
)

dat <- dat %>%
  mutate(
    Dataset = as.character(Dataset),
    Organism = factor(Organism, levels = organism_order)
  ) %>%
  filter(!is.na(Organism))


# ------------------------------------------------------------
# 1. Quadratic mixed-effects model
# ------------------------------------------------------------

fit_quad_mixed <- function(data, response) {
  
  # Backticks allow response names containing special characters
  form1 <- as.formula(
    paste0(
      "`", response, "`",
      " ~ logTP + I(logTP^2) + (1 + logTP || Dataset)"
    )
  )
  
  lmer(
    form1,
    data = data,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5),
      check.conv.singular = .makeCC(
        action = "ignore",
        tol = 1e-4
      ),
      check.conv.grad = .makeCC(
        action = "ignore",
        tol = 1e-3
      )
    )
  )
}


# ------------------------------------------------------------
# 2. Format model coefficients
# ------------------------------------------------------------

get_math_labels_two_lines <- function(model) {
  
  coefs <- summary(model)$coefficients
  
  b1 <- coefs["logTP", "Estimate"]
  s1 <- coefs["logTP", "Std. Error"]
  
  b2 <- coefs["I(logTP^2)", "Estimate"]
  s2 <- coefs["I(logTP^2)", "Std. Error"]
  
  f_fmt <- function(val, se, idx) {
    
    b_sci <- formatC(val, format = "e", digits = 1)
    b_p   <- strsplit(b_sci, "e")[[1]]
    
    s_sci <- formatC(se, format = "e", digits = 1)
    s_p   <- strsplit(s_sci, "e")[[1]]
    
    sprintf(
      paste0(
        "italic(beta)[%s]==%.1f %%*%% 10^%d",
        "~(SE==%.1f %%*%% 10^%d)"
      ),
      idx,
      as.numeric(b_p[1]),
      as.integer(b_p[2]),
      as.numeric(s_p[1]),
      as.integer(s_p[2])
    )
  }
  
  paste0(
    "atop(",
    f_fmt(b1, s1, "1"),
    ",",
    f_fmt(b2, s2, "2"),
    ")"
  )
}


# ------------------------------------------------------------
# 3. Colour palettes
# ------------------------------------------------------------

make_base_palette <- function(response) {
  
  if (response == "Richness_z") {
    
    # Blue: taxonomic richness
    c("#08306B", "#4292C6", "#C6DBEF")
    
  } else if (response == "FRic_z") {
    
    # Red: observed FRic
    c("#99000D", "#EF3B2C", "#FCBBA1")
    
  } else if (response == "SES_FRic") {
    
    # Purple: SES-FRic
    c("#54278F", "#756BB1", "#CBC9E2")
    
  } else {
    
    c("#444444", "#888888", "#DDDDDD")
  }
}


# ------------------------------------------------------------
# 4. Plotting function
# ------------------------------------------------------------

make_organism_plot <- function(
    data,
    response,
    y_lab,
    ncol = 3,
    add_zero_line = FALSE,
    y_limits = c(-4, 4.5)
) {
  
  # Retain only observations available for the focal response
  plot_data <- data %>%
    filter(
      !is.na(.data[[response]]),
      !is.na(logTP),
      !is.na(Dataset),
      !is.na(Organism)
    ) %>%
    droplevels()
  
  group_sym <- sym("Organism")
  
  base_pal <- make_base_palette(response)
  main_col <- base_pal[1]
  
  groups <- levels(droplevels(plot_data$Organism))
  
  dataset_color_vec <- c()
  
  for (g in groups) {
    
    g_dat <- plot_data %>%
      filter(Organism == g)
    
    g_datasets <- sort(unique(as.character(g_dat$Dataset)))
    n_ds <- length(g_datasets)
    
    pal <- colorRampPalette(base_pal)(max(n_ds, 3))[seq_len(n_ds)]
    names(pal) <- g_datasets
    
    dataset_color_vec <- c(dataset_color_vec, pal)
  }
  
  nested_res <- plot_data %>%
    group_by(!!group_sym) %>%
    nest() %>%
    mutate(
      model = map(
        data,
        ~ fit_quad_mixed(.x, response)
      ),
      
      pred_ind = map2(
        data,
        model,
        ~ {
          tp_range <- .x %>%
            group_by(Dataset) %>%
            summarise(
              min_t = min(logTP, na.rm = TRUE),
              max_t = max(logTP, na.rm = TRUE),
              .groups = "drop"
            )
          
          ind_list <- lapply(
            seq_len(nrow(tp_range)),
            function(j) {
              
              df_sub <- data.frame(
                Dataset = tp_range$Dataset[j],
                logTP = seq(
                  tp_range$min_t[j],
                  tp_range$max_t[j],
                  length.out = 80
                )
              )
              
              tryCatch(
                {
                  df_sub$pred <- predict(
                    .y,
                    newdata = df_sub,
                    re.form = NULL,
                    allow.new.levels = TRUE
                  )
                  
                  df_sub
                },
                error = function(e) NULL
              )
            }
          )
          
          bind_rows(Filter(Negate(is.null), ind_list))
        }
      ),
      
      pred_overall = map(
        model,
        ~ as.data.frame(
          ggpredict(.x, terms = "logTP [all]")
        )
      ),
      
      label_str = map_chr(
        model,
        get_math_labels_two_lines
      )
    )
  
  df_points <- nested_res %>%
    select(!!group_sym, data) %>%
    unnest(data) %>%
    mutate(
      col = dataset_color_vec[as.character(Dataset)]
    )
  
  df_ind <- nested_res %>%
    select(!!group_sym, pred_ind) %>%
    unnest(pred_ind) %>%
    mutate(
      col = dataset_color_vec[as.character(Dataset)]
    )
  
  df_overall <- nested_res %>%
    select(!!group_sym, pred_overall) %>%
    unnest(pred_overall) %>%
    mutate(col = main_col)
  
  df_labels <- nested_res %>%
    ungroup() %>%
    mutate(
      x = map_dbl(data, ~ max(.x$logTP, na.rm = TRUE)),
      y = y_limits[2] - 0.7,
      col = main_col
    )
  
  p <- ggplot()
  
  # Add the null expectation only for SES-FRic
  if (add_zero_line) {
    p <- p +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.5,
        colour = "grey55"
      )
  }
  
  p +
    geom_point(
      data = df_points,
      aes(
        x = logTP,
        y = .data[[response]],
        color = col
      ),
      alpha = 0.20,
      size = 0.5
    ) +
    geom_line(
      data = df_ind,
      aes(
        x = logTP,
        y = pred,
        color = col,
        group = Dataset
      ),
      linewidth = 0.6,
      alpha = 0.6
    ) +
    geom_ribbon(
      data = df_overall,
      aes(
        x = x,
        ymin = conf.low,
        ymax = conf.high,
        fill = col
      ),
      alpha = 0.15
    ) +
    geom_line(
      data = df_overall,
      aes(
        x = x,
        y = predicted,
        color = col
      ),
      linewidth = 1.3,
      alpha = 0.8
    ) +
    geom_text(
      data = df_labels,
      aes(
        x = x,
        y = y,
        label = label_str,
        color = col
      ),
      parse = TRUE,
      hjust = 1.05,
      size = 2.8
    ) +
    scale_color_identity() +
    scale_fill_identity() +
    facet_wrap(
      ~ Organism,
      ncol = ncol,
      labeller = as_labeller(organism_labels)
    ) +
    coord_cartesian(ylim = y_limits) +
    scale_y_continuous(
      breaks = seq(
        ceiling(y_limits[1]),
        floor(y_limits[2]),
        by = 2
      )
    ) +
    labs(
      x = expression(log(TP + 1)),
      y = y_lab
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}


# ------------------------------------------------------------
# 5. Generate the three panels
# ------------------------------------------------------------

p_organism_rich <- make_organism_plot(
  data = dat,
  response = "Richness_z",
  y_lab = "Z-Richness",
  y_limits = c(-4, 4.5)
)

p_organism_fric <- make_organism_plot(
  data = dat,
  response = "FRic_z",
  y_lab = "Z-FRic",
  y_limits = c(-4, 4.5)
)

# SES_FRic is used directly and is NOT Z-transformed
p_organism_ses <- make_organism_plot(
  data = dat,
  response = "SES_FRic",
  y_lab = "SES-FRic",
  add_zero_line = TRUE,
  y_limits = c(-4, 4.5)
)


# ------------------------------------------------------------
# 6. Combine the plots
# ------------------------------------------------------------

p_organism_combined <-
  p_organism_rich /
  p_organism_fric /
  p_organism_ses +
  plot_annotation(tag_levels = "a")

print(p_organism_combined)


ggsave(file.path(output_dir, "FigureS2.png"), plot = p_organism_combined, width = 10, height = 20, dpi = 600)
