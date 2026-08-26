# ==============================================================================
# Site-level SES-FRic with richness-controlled null assemblages
# ==============================================================================
#
# Purpose
#   1. Calculate observed functional richness (FRic) for every site.
#   2. Remove the intrinsic effect of taxonomic richness using null models.
#   3. Use r1 as the primary null model:
#        - preserve the number of species with complete traits at each site;
#        - sample species without replacement from the dataset-specific pool;
#        - weight sampling by species incidence across sites.
#   4. Use r0 as a sensitivity analysis:
#        - preserve site richness;
#        - sample all species equiprobably.
#   5. Export SES, empirical P values, null-distribution diagnostics, trait
#      coverage diagnostics, processing logs, and merge diagnostics.
#
# Input structure for every raw .xlsx file
#   sheet "species": first column = Site; remaining columns = species
#   sheet "traits":  first column = Species; remaining columns = traits
#
# Main formula
#   SES_FRic = (FRic_observed - mean(FRic_null)) / sd(FRic_null)
#
# Important
#   - The response used in downstream models should be SES_FRic directly.
#   - Do not Z-standardize SES_FRic again within datasets, because SES = 0 is
#     the biologically meaningful null expectation.
# ==============================================================================

# ==============================================================================
# 0. Packages
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(cluster)
  library(ape)
  library(geometry)
})

# ==============================================================================
# 1. Paths and analysis settings
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

raw_data_path <- file.path("..", "000_input_file", "raw_data")
all_data_file <- file.path("..", "000_input_file", "All_data.xlsx")
output_dir <- file.path("..", "000_output_file")
output_file <- file.path(
  output_dir,
  "All_data_with_SES_FRic_null_models.xlsx"
)

# Optional RDS containing every randomized FRic value. This can be large.
save_null_fric_draws_rds <- FALSE
null_draws_rds_file <- file.path(
  output_dir,
  "SES_FRic_null_draws.rds"
)

# Column names used to merge calculated results into All_data.xlsx.
dataset_column_all_data <- "Dataset"
site_column_all_data <- "Site"

# Primary analysis and sensitivity analysis.
null_models_to_run <- c("r1", "r0")
primary_null_model <- "r1"

# Formal analysis: use at least 999 randomizations.
n_iterations <- 999L

# Retain at most four PCoA axes. The actual number is also limited by the
# richness of the smallest eligible assemblage within each dataset.
n_axes_max <- 4L

# Sites with fewer usable species are reported but receive no SES-FRic.
minimum_species_with_traits <- 3L

# Require almost all randomized FRic values to be valid.
minimum_valid_null_fraction <- 0.95

# Numerical threshold below which null SD is treated as zero.
null_sd_tolerance <- 1e-12

# Flag null distributions with substantial skewness.
null_skewness_warning_threshold <- 1.0

# Trait coverage is always exported. By default, low-coverage sites are flagged
# but retained. Set exclude_sites_below_trait_coverage to TRUE for a strict
# coverage-filtered sensitivity analysis.
trait_coverage_flag_threshold <- 0.80
exclude_sites_below_trait_coverage <- FALSE

# Reproducible base seed. A deterministic dataset-model-specific seed is
# generated from this value so results do not depend on file processing order.
random_seed <- 321L

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(raw_data_path)) {
  stop(
    "Raw-data directory not found: ",
    normalizePath(raw_data_path, mustWork = FALSE)
  )
}

if (!file.exists(all_data_file)) {
  stop(
    "All_data.xlsx not found: ",
    normalizePath(all_data_file, mustWork = FALSE)
  )
}

if (length(null_models_to_run) == 0L ||
    !all(null_models_to_run %in% c("r0", "r1"))) {
  stop("null_models_to_run must contain only 'r0' and/or 'r1'.")
}

if (!primary_null_model %in% null_models_to_run) {
  stop("primary_null_model must be included in null_models_to_run.")
}

if (n_iterations < 2L) {
  stop("n_iterations must be at least 2.")
}

if (n_axes_max < 1L) {
  stop("n_axes_max must be at least 1.")
}

if (minimum_species_with_traits < 2L) {
  stop("minimum_species_with_traits must be at least 2.")
}

if (minimum_valid_null_fraction <= 0 ||
    minimum_valid_null_fraction > 1) {
  stop("minimum_valid_null_fraction must be in (0, 1].")
}

if (trait_coverage_flag_threshold < 0 ||
    trait_coverage_flag_threshold > 1) {
  stop("trait_coverage_flag_threshold must be in [0, 1].")
}

# ==============================================================================
# 2. General helper functions
# ==============================================================================

normalize_identifier <- function(x) {
  x <- stringr::str_trim(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}

dataset_name_from_file <- function(file) {
  tools::file_path_sans_ext(basename(file))
}

seed_from_text <- function(base_seed, text) {
  integer_values <- utf8ToInt(enc2utf8(as.character(text)))

  if (length(integer_values) == 0L) {
    return(as.integer(base_seed))
  }

  weighted_sum <- sum(
    (seq_along(integer_values) * as.double(integer_values)) %% 2147483646
  )

  as.integer((as.double(base_seed) + weighted_sum) %% 2147483646 + 1)
}

clean_trait_column <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    return(as.numeric(x))
  }

  x_character <- stringr::str_trim(as.character(x))
  x_character[x_character %in% c("", "NA", "NaN", "NULL")] <- NA_character_

  x_numeric <- suppressWarnings(as.numeric(x_character))
  observed <- !is.na(x_character)

  if (all(!observed | !is.na(x_numeric))) {
    return(x_numeric)
  }

  factor(x_character)
}

is_informative_trait <- function(x) {
  values <- x[!is.na(x)]
  length(values) > 0L && dplyr::n_distinct(values) > 1L
}

adjusted_sample_skewness <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)

  if (n < 3L) {
    return(NA_real_)
  }

  x_sd <- stats::sd(x)

  if (!is.finite(x_sd) || x_sd <= 0) {
    return(NA_real_)
  }

  x_mean <- mean(x)

  n / ((n - 1) * (n - 2)) * sum(((x - x_mean) / x_sd)^3)
}

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  as.numeric(
    stats::quantile(
      x,
      probs = probability,
      names = FALSE,
      type = 8,
      na.rm = TRUE
    )
  )
}

# ==============================================================================
# 3. Construct one dataset-specific functional trait space
# ==============================================================================

make_trait_space <- function(
    trait_table,
    community_pa,
    n_axes_max,
    minimum_species_with_traits) {

  trait_columns_input <- setdiff(names(trait_table), "Species")

  if (length(trait_columns_input) == 0L) {
    stop("The traits sheet contains no trait columns.")
  }

  traits_clean <- trait_table %>%
    mutate(
      Species = normalize_identifier(Species),
      across(all_of(trait_columns_input), clean_trait_column)
    ) %>%
    filter(!is.na(Species)) %>%
    distinct(Species, .keep_all = TRUE)

  informative_before_missing_filter <- vapply(
    traits_clean[trait_columns_input],
    is_informative_trait,
    logical(1)
  )

  informative_columns <- trait_columns_input[
    informative_before_missing_filter
  ]

  if (length(informative_columns) == 0L) {
    stop(
      "No informative trait columns remain after removing empty or constant traits."
    )
  }

  species_before_complete_case_filter <- traits_clean$Species

  traits_clean <- traits_clean %>%
    select(Species, all_of(informative_columns)) %>%
    filter(if_all(all_of(informative_columns), ~ !is.na(.x)))

  # A trait may become constant after incomplete species are excluded.
  informative_after_missing_filter <- vapply(
    traits_clean[informative_columns],
    is_informative_trait,
    logical(1)
  )

  informative_columns <- informative_columns[
    informative_after_missing_filter
  ]

  if (length(informative_columns) == 0L) {
    stop(
      "No informative trait columns remain after excluding incomplete species."
    )
  }

  traits_clean <- traits_clean %>%
    select(Species, all_of(informative_columns))

  species_pool <- intersect(colnames(community_pa), traits_clean$Species)
  species_pool <- species_pool[
    colSums(community_pa[, species_pool, drop = FALSE]) > 0
  ]

  if (length(species_pool) < minimum_species_with_traits) {
    stop(
      "Fewer than ",
      minimum_species_with_traits,
      " occurring species have complete trait information."
    )
  }

  traits_for_distance <- traits_clean %>%
    filter(Species %in% species_pool) %>%
    slice(match(species_pool, Species)) %>%
    column_to_rownames("Species")

  gower_distance <- cluster::daisy(
    traits_for_distance,
    metric = "gower"
  )

  pcoa_result <- ape::pcoa(
    gower_distance,
    correction = "cailliez"
  )

  coordinate_matrix <- pcoa_result$vectors.cor

  if (is.null(coordinate_matrix) || ncol(coordinate_matrix) == 0L) {
    coordinate_matrix <- pcoa_result$vectors
  }

  if (is.null(coordinate_matrix) || ncol(coordinate_matrix) == 0L) {
    stop("PCoA produced no coordinate axes.")
  }

  coordinate_matrix <- as.matrix(coordinate_matrix)
  rownames(coordinate_matrix) <- rownames(traits_for_distance)

  finite_axes <- vapply(
    seq_len(ncol(coordinate_matrix)),
    function(axis_index) {
      axis_values <- coordinate_matrix[, axis_index]
      all(is.finite(axis_values)) && stats::var(axis_values) > 1e-14
    },
    logical(1)
  )

  coordinate_matrix <- coordinate_matrix[, finite_axes, drop = FALSE]

  if (ncol(coordinate_matrix) == 0L) {
    stop("PCoA produced no positive, variable coordinate axes.")
  }

  community_with_traits <- community_pa[
    ,
    rownames(coordinate_matrix),
    drop = FALSE
  ]

  richness_with_traits <- rowSums(community_with_traits)
  richness_for_axis_selection <- richness_with_traits[
    richness_with_traits >= minimum_species_with_traits
  ]

  if (length(richness_for_axis_selection) == 0L) {
    stop(
      "No site contains at least ",
      minimum_species_with_traits,
      " species with complete trait information."
    )
  }

  richness_axis_limit <- min(richness_for_axis_selection) - 1L

  n_axes <- min(
    ncol(coordinate_matrix),
    n_axes_max,
    richness_axis_limit
  )

  if (!is.finite(n_axes) || n_axes < 1L) {
    stop("No PCoA axis can be retained under the site-richness constraint.")
  }

  trait_space <- coordinate_matrix[, seq_len(n_axes), drop = FALSE]

  centered_coordinates <- scale(
    coordinate_matrix,
    center = TRUE,
    scale = FALSE
  )
  axis_inertia <- colSums(centered_coordinates^2)
  total_inertia <- sum(axis_inertia)

  trait_space_quality <- if (is.finite(total_inertia) && total_inertia > 0) {
    sum(axis_inertia[seq_len(n_axes)]) / total_inertia
  } else {
    NA_real_
  }

  list(
    trait_space = trait_space,
    community_with_traits = community_with_traits,
    informative_trait_names = informative_columns,
    n_traits_input = length(trait_columns_input),
    n_traits_used = length(informative_columns),
    n_species_traits_before_complete_filter = length(
      species_before_complete_case_filter
    ),
    n_species_pool = length(species_pool),
    n_axes_available = ncol(coordinate_matrix),
    n_axes = n_axes,
    trait_space_quality = trait_space_quality
  )
}

# ==============================================================================
# 4. FRic calculation
# ==============================================================================

calculate_fric <- function(species, trait_space) {
  species <- unique(intersect(species, rownames(trait_space)))

  if (length(species) == 0L) {
    return(NA_real_)
  }

  points_all_species <- trait_space[species, , drop = FALSE]
  dimensions <- ncol(points_all_species)

  if (dimensions == 1L) {
    # At least two species are needed to define a one-dimensional range.
    if (length(species) < 2L) {
      return(NA_real_)
    }

    unique_values <- unique(as.numeric(points_all_species[, 1]))

    # Multiple species with the same trait coordinate have zero FRic.
    if (length(unique_values) < 2L) {
      return(0)
    }

    return(diff(range(unique_values, na.rm = TRUE)))
  }

  # The observed species count must be greater than the retained dimension.
  if (length(species) <= dimensions) {
    return(NA_real_)
  }

  # Remove duplicate functional coordinates only after checking the original
  # species count. Duplicate or affinely dependent points imply zero volume,
  # not a missing value.
  points <- unique(as.data.frame(points_all_species))
  points <- as.matrix(points)

  if (nrow(points) <= dimensions) {
    return(0)
  }

  centered_points <- sweep(
    points,
    MARGIN = 2,
    STATS = points[1, ],
    FUN = "-"
  )

  affine_rank <- qr(
    centered_points,
    tol = 1e-10
  )$rank

  if (affine_rank < dimensions) {
    return(0)
  }

  hull <- tryCatch(
    geometry::convhulln(
      points,
      options = "FA Qt Qbb"
    ),
    error = function(e) NULL
  )

  if (is.null(hull) ||
      is.null(hull$vol) ||
      !is.finite(hull$vol)) {
    return(NA_real_)
  }

  as.numeric(hull$vol)
}

calculate_observed_fric <- function(
    community_with_traits,
    trait_space,
    minimum_species_with_traits) {

  apply(community_with_traits, 1, function(site_row) {
    present_species <- colnames(community_with_traits)[site_row > 0]

    if (length(present_species) < minimum_species_with_traits) {
      return(NA_real_)
    }

    calculate_fric(present_species, trait_space)
  })
}

# ==============================================================================
# 5. Richness-controlled null references
# ==============================================================================

calculate_null_reference <- function(
    community_with_traits,
    trait_space,
    richness_values,
    n_iterations,
    null_model,
    minimum_species_with_traits,
    minimum_valid_null_fraction,
    null_sd_tolerance,
    null_skewness_warning_threshold,
    simulation_seed) {

  species_pool <- colnames(community_with_traits)
  incidence <- colSums(community_with_traits > 0)

  sampling_probability <- if (null_model == "r1") {
    incidence
  } else {
    NULL
  }

  minimum_valid <- ceiling(
    n_iterations * minimum_valid_null_fraction
  )

  unique_richness <- sort(unique(as.integer(richness_values)))
  summary_rows <- vector("list", length(unique_richness))
  null_draws <- setNames(
    vector("list", length(unique_richness)),
    as.character(unique_richness)
  )

  set.seed(simulation_seed)

  for (richness_index in seq_along(unique_richness)) {
    site_richness <- unique_richness[richness_index]
    richness_key <- as.character(site_richness)

    if (is.na(site_richness) ||
        site_richness < minimum_species_with_traits ||
        site_richness > length(species_pool)) {

      null_draws[[richness_key]] <- numeric(0)

      summary_rows[[richness_index]] <- tibble(
        Richness_with_traits = site_richness,
        FRic_null_mean = NA_real_,
        FRic_null_median = NA_real_,
        FRic_null_sd = NA_real_,
        FRic_null_q025 = NA_real_,
        FRic_null_q975 = NA_real_,
        FRic_null_skewness = NA_real_,
        FRic_null_zero_fraction = NA_real_,
        Null_iterations_valid = 0L,
        Null_iterations_invalid = n_iterations,
        Null_valid_fraction = 0,
        Null_distribution_skewed = NA,
        Null_reference_status = "invalid_or_ineligible_richness"
      )

      next
    }

    null_fric_raw <- replicate(n_iterations, {
      randomized_species <- sample(
        species_pool,
        size = site_richness,
        replace = FALSE,
        prob = sampling_probability
      )

      calculate_fric(randomized_species, trait_space)
    })

    null_fric_raw <- as.numeric(null_fric_raw)
    null_fric <- null_fric_raw[is.finite(null_fric_raw)]
    null_draws[[richness_key]] <- null_fric

    valid_iterations <- length(null_fric)
    invalid_iterations <- n_iterations - valid_iterations
    valid_fraction <- valid_iterations / n_iterations

    null_mean <- if (valid_iterations > 0L) {
      mean(null_fric)
    } else {
      NA_real_
    }

    null_median <- if (valid_iterations > 0L) {
      stats::median(null_fric)
    } else {
      NA_real_
    }

    null_sd <- if (valid_iterations > 1L) {
      stats::sd(null_fric)
    } else {
      NA_real_
    }

    null_skewness <- adjusted_sample_skewness(null_fric)

    zero_fraction <- if (valid_iterations > 0L) {
      mean(null_fric == 0)
    } else {
      NA_real_
    }

    null_status <- dplyr::case_when(
      valid_iterations < minimum_valid ~ "too_few_valid_null_values",
      is.na(null_sd) ~ "null_sd_missing",
      null_sd <= null_sd_tolerance ~ "null_sd_zero_or_near_zero",
      TRUE ~ "ok"
    )

    summary_rows[[richness_index]] <- tibble(
      Richness_with_traits = site_richness,
      FRic_null_mean = null_mean,
      FRic_null_median = null_median,
      FRic_null_sd = null_sd,
      FRic_null_q025 = safe_quantile(null_fric, 0.025),
      FRic_null_q975 = safe_quantile(null_fric, 0.975),
      FRic_null_skewness = null_skewness,
      FRic_null_zero_fraction = zero_fraction,
      Null_iterations_valid = valid_iterations,
      Null_iterations_invalid = invalid_iterations,
      Null_valid_fraction = valid_fraction,
      Null_distribution_skewed = if_else(
        is.finite(null_skewness),
        abs(null_skewness) > null_skewness_warning_threshold,
        NA
      ),
      Null_reference_status = null_status
    )
  }

  list(
    summary = bind_rows(summary_rows),
    draws = null_draws
  )
}

calculate_empirical_null_statistics <- function(
    observed_value,
    null_values) {

  null_values <- null_values[is.finite(null_values)]
  n_null <- length(null_values)

  if (!is.finite(observed_value) || n_null == 0L) {
    return(tibble(
      FRic_empirical_p_lower = NA_real_,
      FRic_empirical_p_upper = NA_real_,
      FRic_empirical_p_two_sided = NA_real_,
      FRic_null_percentile = NA_real_,
      FRic_probit_effect = NA_real_
    ))
  }

  p_lower <- (1 + sum(null_values <= observed_value)) / (n_null + 1)
  p_upper <- (1 + sum(null_values >= observed_value)) / (n_null + 1)
  p_two_sided <- min(1, 2 * min(p_lower, p_upper))

  percentile <- (
    sum(null_values < observed_value) +
      0.5 * sum(null_values == observed_value) +
      0.5
  ) / (n_null + 1)

  # Avoid infinite qnorm values at numerical boundaries.
  percentile_for_probit <- min(
    max(percentile, 0.5 / (n_null + 1)),
    1 - 0.5 / (n_null + 1)
  )

  tibble(
    FRic_empirical_p_lower = p_lower,
    FRic_empirical_p_upper = p_upper,
    FRic_empirical_p_two_sided = p_two_sided,
    FRic_null_percentile = percentile,
    FRic_probit_effect = stats::qnorm(percentile_for_probit)
  )
}

# ==============================================================================
# 6. Process one raw dataset
# ==============================================================================

process_one_dataset <- function(file) {
  dataset_name <- dataset_name_from_file(file)
  message("Processing SES-FRic: ", dataset_name)

  community_raw <- readxl::read_excel(file, sheet = "species")
  trait_raw <- readxl::read_excel(file, sheet = "traits")

  if (ncol(community_raw) < 2L) {
    stop("The species sheet has fewer than two columns.")
  }

  if (ncol(trait_raw) < 2L) {
    stop("The traits sheet has fewer than two columns.")
  }

  names(community_raw)[1] <- "Site"
  names(trait_raw)[1] <- "Species"

  community_raw$Site <- normalize_identifier(community_raw$Site)
  trait_raw$Species <- normalize_identifier(trait_raw$Species)

  species_names <- normalize_identifier(names(community_raw)[-1])

  if (anyNA(species_names) || anyDuplicated(species_names)) {
    stop("Species names in the species sheet are missing or duplicated.")
  }

  names(community_raw)[-1] <- species_names

  if (anyNA(community_raw$Site) || anyDuplicated(community_raw$Site)) {
    stop("Site identifiers in the species sheet are missing or duplicated.")
  }

  community_numeric <- community_raw %>%
    select(-Site) %>%
    mutate(
      across(
        everything(),
        ~ suppressWarnings(as.numeric(as.character(.x)))
      )
    )

  community_matrix <- as.matrix(community_numeric)
  community_matrix[is.na(community_matrix)] <- 0

  # Presence-absence conversion: positive value = present; otherwise absent.
  community_pa <- ifelse(community_matrix > 0, 1L, 0L)
  community_pa <- matrix(
    community_pa,
    nrow = nrow(community_matrix),
    ncol = ncol(community_matrix),
    dimnames = list(
      community_raw$Site,
      colnames(community_matrix)
    )
  )

  occurring_species <- colSums(community_pa) > 0
  community_pa <- community_pa[, occurring_species, drop = FALSE]

  if (ncol(community_pa) < minimum_species_with_traits) {
    stop(
      "Fewer than ",
      minimum_species_with_traits,
      " species occur in this dataset."
    )
  }

  trait_result <- make_trait_space(
    trait_table = trait_raw,
    community_pa = community_pa,
    n_axes_max = n_axes_max,
    minimum_species_with_traits = minimum_species_with_traits
  )

  community_with_traits <- trait_result$community_with_traits
  trait_space <- trait_result$trait_space

  richness_raw <- rowSums(community_pa)
  richness_with_traits <- rowSums(community_with_traits)

  trait_species_coverage <- ifelse(
    richness_raw > 0,
    richness_with_traits / richness_raw,
    NA_real_
  )

  coverage_below_threshold <- is.finite(trait_species_coverage) &
    trait_species_coverage < trait_coverage_flag_threshold

  coverage_is_eligible <- if (exclude_sites_below_trait_coverage) {
    !coverage_below_threshold
  } else {
    rep(TRUE, length(coverage_below_threshold))
  }

  observed_fric <- calculate_observed_fric(
    community_with_traits = community_with_traits,
    trait_space = trait_space,
    minimum_species_with_traits = minimum_species_with_traits
  )

  site_base <- tibble(
    Dataset = dataset_name,
    Site = rownames(community_pa),
    Richness_raw_PA = as.integer(richness_raw),
    Richness_with_traits = as.integer(richness_with_traits),
    Trait_species_coverage = as.numeric(trait_species_coverage),
    Trait_coverage_below_threshold = coverage_below_threshold,
    Trait_coverage_eligible = coverage_is_eligible,
    FRic_observed_SES_space = as.numeric(observed_fric),
    Species_pool_size = trait_result$n_species_pool,
    Site_pool_fraction = Richness_with_traits / trait_result$n_species_pool,
    Traits_input = trait_result$n_traits_input,
    Traits_used = trait_result$n_traits_used,
    Trait_names_used = paste(
      trait_result$informative_trait_names,
      collapse = " | "
    ),
    PCoA_axes_available = trait_result$n_axes_available,
    PCoA_axes = trait_result$n_axes,
    Trait_space_quality = trait_result$trait_space_quality
  )

  richness_values_for_null <- site_base %>%
    filter(
      Richness_with_traits >= minimum_species_with_traits,
      Trait_coverage_eligible
    ) %>%
    pull(Richness_with_traits)

  if (length(richness_values_for_null) == 0L) {
    stop("No eligible site remains for null-model analysis.")
  }

  site_model_results <- vector("list", length(null_models_to_run))
  null_reference_results <- vector("list", length(null_models_to_run))
  processing_log_results <- vector("list", length(null_models_to_run))
  null_draw_results <- if (save_null_fric_draws_rds) {
    setNames(vector("list", length(null_models_to_run)), null_models_to_run)
  } else {
    NULL
  }

  for (model_index in seq_along(null_models_to_run)) {
    current_null_model <- null_models_to_run[model_index]
    current_seed <- seed_from_text(
      random_seed,
      paste(dataset_name, current_null_model, sep = "::")
    )

    message(
      "  Null model ",
      current_null_model,
      ": ",
      n_iterations,
      " iterations per richness value"
    )

    null_result <- calculate_null_reference(
      community_with_traits = community_with_traits,
      trait_space = trait_space,
      richness_values = richness_values_for_null,
      n_iterations = n_iterations,
      null_model = current_null_model,
      minimum_species_with_traits = minimum_species_with_traits,
      minimum_valid_null_fraction = minimum_valid_null_fraction,
      null_sd_tolerance = null_sd_tolerance,
      null_skewness_warning_threshold = null_skewness_warning_threshold,
      simulation_seed = current_seed
    )

    null_reference <- null_result$summary %>%
      mutate(
        Dataset = dataset_name,
        Null_model = current_null_model,
        Null_seed = current_seed,
        Null_iterations_requested = n_iterations,
        Species_pool_size = trait_result$n_species_pool,
        PCoA_axes = trait_result$n_axes,
        Trait_space_quality = trait_result$trait_space_quality,
        .before = 1
      )

    empirical_statistics <- purrr::map2_dfr(
      site_base$FRic_observed_SES_space,
      site_base$Richness_with_traits,
      function(observed_value, site_richness) {
        richness_key <- as.character(site_richness)
        null_values <- null_result$draws[[richness_key]]

        if (is.null(null_values)) {
          null_values <- numeric(0)
        }

        calculate_empirical_null_statistics(
          observed_value,
          null_values
        )
      }
    )

    site_results_current <- site_base %>%
      left_join(
        null_reference %>%
          select(
            Richness_with_traits,
            FRic_null_mean,
            FRic_null_median,
            FRic_null_sd,
            FRic_null_q025,
            FRic_null_q975,
            FRic_null_skewness,
            FRic_null_zero_fraction,
            Null_iterations_valid,
            Null_iterations_invalid,
            Null_valid_fraction,
            Null_distribution_skewed,
            Null_reference_status
          ),
        by = "Richness_with_traits"
      ) %>%
      bind_cols(empirical_statistics) %>%
      mutate(
        SES_FRic = if_else(
          Trait_coverage_eligible &
            Null_reference_status == "ok" &
            is.finite(FRic_observed_SES_space),
          (FRic_observed_SES_space - FRic_null_mean) / FRic_null_sd,
          NA_real_
        ),
        SES_status = case_when(
          Richness_raw_PA == 0L ~ "no_species_present",
          Richness_with_traits == 0L ~ "no_species_with_complete_traits",
          Richness_with_traits < minimum_species_with_traits ~
            "too_few_species_with_complete_traits",
          !Trait_coverage_eligible ~ "trait_coverage_below_threshold",
          !is.finite(FRic_observed_SES_space) ~
            "observed_fric_not_estimable",
          is.na(Null_reference_status) ~ "null_reference_missing",
          Null_reference_status != "ok" ~ Null_reference_status,
          !is.finite(SES_FRic) ~ "ses_not_estimable",
          TRUE ~ "ok"
        ),
        Analysis_ready = SES_status == "ok",
        Coverage_sensitivity_ready = SES_status == "ok" &
          !Trait_coverage_below_threshold,
        Null_model = current_null_model,
        Null_seed = current_seed,
        Null_iterations_requested = n_iterations
      )

    site_model_results[[model_index]] <- site_results_current
    null_reference_results[[model_index]] <- null_reference

    processing_log_results[[model_index]] <- tibble(
      Dataset = dataset_name,
      File = basename(file),
      Null_model = current_null_model,
      Processing_status = "ok",
      Error_message = NA_character_,
      Sites = nrow(community_pa),
      Species_occurring = ncol(community_pa),
      Species_with_complete_traits = trait_result$n_species_pool,
      Traits_input = trait_result$n_traits_input,
      Traits_used = trait_result$n_traits_used,
      PCoA_axes_available = trait_result$n_axes_available,
      PCoA_axes = trait_result$n_axes,
      Trait_space_quality = trait_result$trait_space_quality,
      Trait_coverage_median = stats::median(
        trait_species_coverage,
        na.rm = TRUE
      ),
      Trait_coverage_minimum = suppressWarnings(
        min(trait_species_coverage, na.rm = TRUE)
      ),
      Sites_below_trait_coverage_threshold = sum(
        coverage_below_threshold,
        na.rm = TRUE
      ),
      Sites_with_SES = sum(site_results_current$SES_status == "ok"),
      Sites_with_skewed_null = sum(
        site_results_current$Null_distribution_skewed %in% TRUE,
        na.rm = TRUE
      ),
      Null_iterations = n_iterations,
      Null_seed = current_seed
    )

    if (save_null_fric_draws_rds) {
      null_draw_results[[current_null_model]] <- null_result$draws
    }
  }

  list(
    site_results_all_models = bind_rows(site_model_results),
    null_reference = bind_rows(null_reference_results),
    processing_log = bind_rows(processing_log_results),
    null_draws = null_draw_results
  )
}

# ==============================================================================
# 7. Process all dataset-specific raw files
# ==============================================================================

raw_files <- list.files(
  raw_data_path,
  pattern = "\\.xlsx$",
  full.names = TRUE,
  ignore.case = TRUE
)

raw_files <- raw_files[!grepl("^~\\$", basename(raw_files))]
raw_files <- sort(raw_files)

if (length(raw_files) == 0L) {
  stop("No .xlsx files were found in the raw-data directory.")
}

site_result_list <- list()
null_reference_list <- list()
processing_log_list <- list()
all_null_draws <- if (save_null_fric_draws_rds) list() else NULL

for (file_index in seq_along(raw_files)) {
  current_file <- raw_files[file_index]
  current_dataset <- dataset_name_from_file(current_file)

  current_result <- tryCatch(
    process_one_dataset(current_file),
    error = function(error_condition) error_condition
  )

  if (inherits(current_result, "error")) {
    message(
      "Failed: ",
      current_dataset,
      " -> ",
      conditionMessage(current_result)
    )

    processing_log_list[[length(processing_log_list) + 1L]] <- tibble(
      Dataset = current_dataset,
      File = basename(current_file),
      Null_model = NA_character_,
      Processing_status = "error",
      Error_message = conditionMessage(current_result),
      Sites = NA_integer_,
      Species_occurring = NA_integer_,
      Species_with_complete_traits = NA_integer_,
      Traits_input = NA_integer_,
      Traits_used = NA_integer_,
      PCoA_axes_available = NA_integer_,
      PCoA_axes = NA_integer_,
      Trait_space_quality = NA_real_,
      Trait_coverage_median = NA_real_,
      Trait_coverage_minimum = NA_real_,
      Sites_below_trait_coverage_threshold = NA_integer_,
      Sites_with_SES = NA_integer_,
      Sites_with_skewed_null = NA_integer_,
      Null_iterations = n_iterations,
      Null_seed = NA_integer_
    )
  } else {
    site_result_list[[length(site_result_list) + 1L]] <-
      current_result$site_results_all_models

    null_reference_list[[length(null_reference_list) + 1L]] <-
      current_result$null_reference

    processing_log_list[[length(processing_log_list) + 1L]] <-
      current_result$processing_log

    if (save_null_fric_draws_rds) {
      all_null_draws[[current_dataset]] <- current_result$null_draws
    }
  }
}

ses_all_models <- bind_rows(site_result_list)
null_reference_all <- bind_rows(null_reference_list)
processing_log <- bind_rows(processing_log_list)

if (nrow(ses_all_models) == 0L) {
  stop(
    "SES-FRic could not be calculated for any dataset. Check Processing_log."
  )
}

if (anyDuplicated(ses_all_models[c("Dataset", "Site", "Null_model")])) {
  stop("Dataset-Site-Null_model keys are duplicated in the SES results.")
}

# ==============================================================================
# 8. Create primary site table and append sensitivity-model columns
# ==============================================================================

ses_by_site <- ses_all_models %>%
  filter(Null_model == primary_null_model)

if (nrow(ses_by_site) == 0L) {
  stop("No results were produced for primary_null_model.")
}

if (anyDuplicated(ses_by_site[c("Dataset", "Site")])) {
  stop("Dataset-Site keys are duplicated in the primary SES results.")
}

sensitivity_models <- setdiff(
  null_models_to_run,
  primary_null_model
)

sensitivity_columns <- c(
  "FRic_null_mean",
  "FRic_null_sd",
  "FRic_null_skewness",
  "Null_distribution_skewed",
  "FRic_empirical_p_two_sided",
  "FRic_null_percentile",
  "FRic_probit_effect",
  "SES_FRic",
  "SES_status"
)

for (sensitivity_model in sensitivity_models) {
  sensitivity_table <- ses_all_models %>%
    filter(Null_model == sensitivity_model) %>%
    select(
      Dataset,
      Site,
      all_of(sensitivity_columns)
    )

  sensitivity_names <- paste0(
    sensitivity_columns,
    "_",
    sensitivity_model
  )

  names(sensitivity_table)[
    match(sensitivity_columns, names(sensitivity_table))
  ] <- sensitivity_names

  ses_by_site <- ses_by_site %>%
    left_join(
      sensitivity_table,
      by = c("Dataset", "Site")
    )
}

# ==============================================================================
# 9. Merge SES-FRic with All_data.xlsx
# ==============================================================================

all_data <- readxl::read_excel(all_data_file)

required_merge_columns <- c(
  dataset_column_all_data,
  site_column_all_data
)

missing_merge_columns <- setdiff(
  required_merge_columns,
  names(all_data)
)

if (length(missing_merge_columns) > 0L) {
  stop(
    "The following merge columns are absent from All_data.xlsx: ",
    paste(missing_merge_columns, collapse = ", "),
    ". Change dataset_column_all_data or site_column_all_data at the top ",
    "of the script."
  )
}

all_data_for_join <- all_data %>%
  mutate(
    .Dataset_join = normalize_identifier(
      .data[[dataset_column_all_data]]
    ),
    .Site_join = normalize_identifier(
      .data[[site_column_all_data]]
    )
  )

ses_for_join <- ses_by_site %>%
  mutate(
    .Dataset_join = normalize_identifier(Dataset),
    .Site_join = normalize_identifier(Site)
  ) %>%
  select(-Dataset, -Site)

# If this script is rerun on a previously merged workbook, replace old SES
# output columns instead of creating .x and .y duplicates.
calculated_columns <- setdiff(
  names(ses_for_join),
  c(".Dataset_join", ".Site_join")
)

existing_calculated_columns <- intersect(
  names(all_data_for_join),
  calculated_columns
)

if (length(existing_calculated_columns) > 0L) {
  message(
    "Replacing existing calculated SES columns in All_data: ",
    paste(existing_calculated_columns, collapse = ", ")
  )

  all_data_for_join <- all_data_for_join %>%
    select(-all_of(existing_calculated_columns))
}

merged_data <- all_data_for_join %>%
  left_join(
    ses_for_join,
    by = c(".Dataset_join", ".Site_join")
  ) %>%
  select(-.Dataset_join, -.Site_join)

unmatched_all_data_keys <- all_data_for_join %>%
  distinct(.Dataset_join, .Site_join) %>%
  anti_join(
    ses_for_join %>% distinct(.Dataset_join, .Site_join),
    by = c(".Dataset_join", ".Site_join")
  ) %>%
  transmute(
    Dataset = .Dataset_join,
    Site = .Site_join
  )

unmatched_raw_keys <- ses_for_join %>%
  distinct(.Dataset_join, .Site_join) %>%
  anti_join(
    all_data_for_join %>% distinct(.Dataset_join, .Site_join),
    by = c(".Dataset_join", ".Site_join")
  ) %>%
  transmute(
    Dataset = .Dataset_join,
    Site = .Site_join
  )

merge_summary <- tibble(
  Measure = c(
    "Rows in original All_data",
    "Rows in merged data",
    "Rows matched to raw SES results",
    "Rows with estimable primary SES_FRic",
    "Rows passing trait-coverage sensitivity threshold",
    "Unique unmatched Dataset-Site keys in All_data",
    "Unique raw Dataset-Site keys absent from All_data"
  ),
  Value = c(
    nrow(all_data),
    nrow(merged_data),
    sum(!is.na(merged_data$SES_status)),
    sum(is.finite(merged_data$SES_FRic)),
    sum(merged_data$Coverage_sensitivity_ready %in% TRUE, na.rm = TRUE),
    nrow(unmatched_all_data_keys),
    nrow(unmatched_raw_keys)
  )
)

# ==============================================================================
# 10. Richness-dependence diagnostics
# ==============================================================================

calculate_richness_diagnostic <- function(data) {
  data_valid <- data %>%
    filter(
      SES_status == "ok",
      is.finite(SES_FRic),
      is.finite(Richness_with_traits)
    )

  if (nrow(data_valid) < 3L ||
      dplyr::n_distinct(data_valid$Richness_with_traits) < 2L) {
    return(tibble(
      Sites_used = nrow(data_valid),
      Pearson_r_SES_richness = NA_real_,
      Spearman_rho_SES_richness = NA_real_
    ))
  }

  tibble(
    Sites_used = nrow(data_valid),
    Pearson_r_SES_richness = suppressWarnings(
      stats::cor(
        data_valid$SES_FRic,
        data_valid$Richness_with_traits,
        method = "pearson"
      )
    ),
    Spearman_rho_SES_richness = suppressWarnings(
      stats::cor(
        data_valid$SES_FRic,
        data_valid$Richness_with_traits,
        method = "spearman"
      )
    )
  )
}

richness_diagnostics_by_dataset <- ses_all_models %>%
  group_by(Dataset, Null_model) %>%
  group_modify(~ calculate_richness_diagnostic(.x)) %>%
  ungroup()

richness_diagnostics_overall <- ses_all_models %>%
  group_by(Null_model) %>%
  group_modify(~ calculate_richness_diagnostic(.x)) %>%
  ungroup() %>%
  mutate(Dataset = "__OVERALL__", .before = 1)

richness_diagnostics <- bind_rows(
  richness_diagnostics_overall,
  richness_diagnostics_by_dataset
)

# ==============================================================================
# 11. Export analysis settings and outputs
# ==============================================================================

analysis_settings <- tibble(
  Setting = c(
    "Primary null model",
    "Null models run",
    "Iterations per dataset-model-richness combination",
    "Maximum PCoA axes",
    "Minimum species with complete traits",
    "Minimum valid null fraction",
    "Null SD tolerance",
    "Null skewness warning threshold",
    "Trait coverage flag threshold",
    "Exclude sites below trait coverage threshold",
    "Base random seed",
    "Save null FRic draws as RDS",
    "Raw data directory",
    "All_data input file",
    "Main output file"
  ),
  Value = as.character(c(
    primary_null_model,
    paste(null_models_to_run, collapse = ", "),
    n_iterations,
    n_axes_max,
    minimum_species_with_traits,
    minimum_valid_null_fraction,
    null_sd_tolerance,
    null_skewness_warning_threshold,
    trait_coverage_flag_threshold,
    exclude_sites_below_trait_coverage,
    random_seed,
    save_null_fric_draws_rds,
    normalizePath(raw_data_path, mustWork = FALSE),
    normalizePath(all_data_file, mustWork = FALSE),
    normalizePath(output_file, mustWork = FALSE)
  ))
)

output_sheets <- list(
  Merged_data = merged_data,
  SES_by_site = ses_by_site,
  SES_all_models = ses_all_models,
  Null_reference = null_reference_all,
  Richness_diagnostics = richness_diagnostics,
  Processing_log = processing_log,
  Merge_summary = merge_summary,
  Unmatched_All_data = unmatched_all_data_keys,
  Unmatched_raw = unmatched_raw_keys,
  Analysis_settings = analysis_settings
)

excel_max_rows <- 1048576L
too_large_sheets <- names(output_sheets)[
  vapply(output_sheets, nrow, integer(1)) + 1L > excel_max_rows
]

if (length(too_large_sheets) > 0L) {
  stop(
    "The following sheets exceed the Excel row limit: ",
    paste(too_large_sheets, collapse = ", "),
    ". Export these tables as CSV or RDS instead."
  )
}

writexl::write_xlsx(
  output_sheets,
  path = output_file
)

if (save_null_fric_draws_rds) {
  saveRDS(
    all_null_draws,
    file = null_draws_rds_file,
    compress = "xz"
  )
}

message(
  "Completed. Main output written to: ",
  normalizePath(output_file, mustWork = FALSE)
)

if (save_null_fric_draws_rds) {
  message(
    "Null FRic draws written to: ",
    normalizePath(null_draws_rds_file, mustWork = FALSE)
  )
}

print(merge_summary)
