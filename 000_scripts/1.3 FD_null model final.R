# ==============================================================================
# Integrated Analysis of Functional Space Shifts across TP Classes:
# 1. Functional Trait-Space Net Change (SES from Incidence-Weighted Null Model)
# 2. Functional-Centroid (Incidence-based CWM) Shifts (SES from Trait-Label Permutation)
# ==============================================================================

# ==============================================================================
# 0. Packages
# ==============================================================================

required_packages <- c(
  "readxl",
  "writexl",
  "dplyr",
  "tibble",
  "tidyr",
  "purrr",
  "stringr",
  "cluster",
  "ape",
  "geometry",
  "ggplot2",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install these packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(cluster)
  library(ape)
  library(geometry)
  library(ggplot2)
  library(patchwork)
})

# ==============================================================================
# 1. Paths and Settings
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
mapping_file <- file.path("..", "000_input_file", "trait_mapping.xlsx")
output_dir <- file.path("..", "000_output_file")

# ------------------------------------------------------------------------------
# THE ONLY SWITCH NEEDED BETWEEN TEST AND FORMAL ANALYSIS
# ------------------------------------------------------------------------------
run_calculation <- FALSE
test_mode <- FALSE  # TRUE = quick workflow test; FALSE = formal 999 iterations

mode_suffix <- if (test_mode) "_TEST" else ""

output_file <- file.path(
  output_dir,
  paste0("FD_null_model", mode_suffix, ".xlsx")
)

draws_rds_file <- file.path(
  output_dir,
  paste0("Integrated_functional_draws", mode_suffix, ".rds")
)


figure_combined_png <- file.path(
  output_dir,
  paste0("Figure_integrated_net_and_centroid_shifts", mode_suffix, ".png")
)


save_draws_rds <- TRUE
resume_from_checkpoints <- TRUE

# TP grouping parameter (top and bottom 20%)
tp_group_fraction <- 0.20
null_model <- "r1" # Incidence-weighted null model

# Permutations & Integration
n_iterations <- if (test_mode) 49L else 999L
n_mc <- if (test_mode) 500L else 2000L
n_axes_max <- 3L
minimum_axes <- 1L

# Strict unified quality control thresholds
minimum_site_richness <- 3L
minimum_species_pool <- 4L
minimum_sites_per_group <- 3L
minimum_valid_fraction <- 0.80
null_sd_tolerance <- 1e-12
random_seed <- 123L
analysis_version <- "integrated_v2_centroid_traitshuffle"

minimum_valid_draws <- max(1L, ceiling(n_iterations * minimum_valid_fraction))

dimension_levels <- c("Habitat", "Resource", "Size")
comparison_levels <- c("Low-Mid", "Mid-High", "Low-High")
tp_group_levels <- c("Low", "Mid", "High")
metric_levels <- c("Loss", "Gain", "Turnover", "Net")

comparison_definitions <- list(
  "Low-Mid"  = c("Low", "Mid"),
  "Mid-High" = c("Mid", "High"),
  "Low-High" = c("Low", "High")
)

broad_dims <- list(
  Size = c("Size"),
  Resource = c(
    "Resource_acquisition",
    "Feeding",
    "Trophic",
    "Trophic_strategy"
  ),
  Habitat = c(
    "Habitat",
    "Life_form",
    "Attachment",
    "Guild",
    "Habit",
    "Buoyancy"
  )
)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

settings_tag <- sprintf(
  "%s_%s_f%02d_n%d_mc%d_a%d_r%d_g%d",
  analysis_version,
  null_model,
  round(100 * tp_group_fraction),
  n_iterations,
  n_mc,
  n_axes_max,
  minimum_site_richness,
  minimum_sites_per_group
)

# ==============================================================================
# 2. Trait Mapping & Helper Functions
# ==============================================================================

required_mapping_columns <- c("Organism", "Trait_category", "Original_trait")

read_trait_mapping <- function(file) {
  if (!file.exists(file)) {
    stop("Trait-mapping file not found: ", normalizePath(file, FALSE))
  }
  available_sheets <- readxl::excel_sheets(file)
  mapping_sheet <- if ("Trait_mapping" %in% available_sheets) "Trait_mapping" else available_sheets[1]
  mapping <- readxl::read_excel(file, sheet = mapping_sheet)
  
  missing_columns <- setdiff(required_mapping_columns, names(mapping))
  if (length(missing_columns) > 0L) {
    stop("Trait mapping is missing columns: ", paste(missing_columns, collapse = ", "))
  }
  
  mapping %>%
    transmute(
      Organism = stringr::str_trim(as.character(Organism)),
      Trait_category = stringr::str_trim(as.character(Trait_category)),
      Original_trait = stringr::str_trim(as.character(Original_trait))
    ) %>%
    filter(
      !is.na(Organism), Organism != "",
      !is.na(Trait_category), Trait_category != "",
      !is.na(Original_trait), Original_trait != ""
    ) %>%
    distinct()
}

normalize_identifier <- function(x) {
  x <- stringr::str_trim(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}

clean_trait_column <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
  x_character <- stringr::str_trim(as.character(x))
  x_character[x_character %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x_numeric <- suppressWarnings(as.numeric(x_character))
  observed <- !is.na(x_character)
  if (all(!observed | !is.na(x_numeric))) return(x_numeric)
  factor(x_character)
}

is_informative_trait <- function(x) {
  values <- x[!is.na(x)]
  length(values) > 1L && dplyr::n_distinct(values) > 1L
}

dataset_name_from_file <- function(file) {
  tools::file_path_sans_ext(basename(file))
}

get_organism <- function(file) {
  dataset_name <- dataset_name_from_file(file)
  code <- substr(dataset_name, 4, 5)
  dplyr::case_when(
    code == "BD" ~ "Diatoms",
    code == "BM" ~ "Macroinvertebrates",
    code == "FI" ~ "Fish",
    code == "PP" ~ "Phytoplankton",
    code == "ZP" ~ "Zooplankton",
    code == "AP" ~ "Aquatic plants",
    TRUE ~ "Unknown"
  )
}

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(stats::quantile(x, probability, names = FALSE, na.rm = TRUE))
}

clamp_value <- function(x, lower, upper) {
  pmin(pmax(x, lower), upper)
}

safe_one_sample_test <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L || stats::sd(x) <= null_sd_tolerance) return(NA_real_)
  stats::t.test(x, mu = 0)$p.value
}

seed_from_text <- function(base_seed, text) {
  text_values <- utf8ToInt(enc2utf8(text))
  if (length(text_values) == 0L) return(as.integer(base_seed))
  text_component <- sum((seq_along(text_values) + 31L) * text_values)
  as.integer(
    (as.double(base_seed) + as.double(text_component)) %%
      (.Machine$integer.max - 1L) + 1L
  )
}

empty_log_row <- function(dataset, organism, dimension = NA_character_, status = "skipped", message = NA_character_) {
  tibble(
    Dataset = dataset,
    Organism = organism,
    Dimension = dimension,
    Processing_status = status,
    Message = message,
    Raw_sites = NA_integer_,
    Eligible_sites = NA_integer_,
    Sites_per_TP_group = NA_integer_,
    Traits_requested = NA_integer_,
    Traits_used = NA_integer_,
    Trait_names_used = NA_character_,
    Species_in_PCoA = NA_integer_,
    PCoA_axes = NA_integer_,
    PCoA_variance_retained = NA_real_
  )
}

# ==============================================================================
# 3. Geometric Convex-Hull, Centroid & Null Statistics Functions
# ==============================================================================

build_hull <- function(species, trait_space) {
  species <- unique(intersect(species, rownames(trait_space)))
  if (length(species) == 0L) return(NULL)
  
  points <- as.matrix(trait_space[species, , drop = FALSE])
  points <- points[apply(points, 1, function(x) all(is.finite(x))), , drop = FALSE]
  points <- unique(points)
  dimensions <- ncol(points)
  
  if (nrow(points) == 0L) return(NULL)
  
  if (dimensions == 1L || nrow(points) <= dimensions) {
    lower <- min(points[, 1])
    upper <- max(points[, 1])
    vol <- max(1e-6, upper - lower)
    return(list(
      type = "interval",
      dimensions = 1L,
      points = points[, 1, drop = FALSE],
      lower = lower,
      upper = upper,
      volume = vol
    ))
  }
  
  hull <- tryCatch(
    geometry::convhulln(points, options = "FA"),
    error = function(e) {
      jittered <- points + matrix(stats::rnorm(length(points), sd = 1e-8), nrow = nrow(points))
      tryCatch(geometry::convhulln(jittered, options = "FA"), error = function(e2) NULL)
    }
  )
  
  if (is.null(hull) || is.null(hull$vol) || !is.finite(hull$vol) || hull$vol <= 0) {
    lower <- min(points[, 1])
    upper <- max(points[, 1])
    return(list(
      type = "interval",
      dimensions = 1L,
      points = points[, 1, drop = FALSE],
      lower = lower,
      upper = upper,
      volume = max(1e-6, upper - lower)
    ))
  }
  
  list(
    type = "convex_hull",
    dimensions = dimensions,
    points = points,
    hull = hull,
    volume = as.numeric(hull$vol)
  )
}

calculate_overlap_volume <- function(hull_a, hull_b, n_points) {
  if (is.null(hull_a) || is.null(hull_b)) return(NA_real_)
  if (hull_a$dimensions != hull_b$dimensions) return(NA_real_)
  
  if (hull_a$dimensions == 1L) {
    overlap <- max(0, min(hull_a$upper, hull_b$upper) - max(hull_a$lower, hull_b$lower))
    return(clamp_value(overlap, 0, min(hull_a$volume, hull_b$volume)))
  }
  
  bounds_a <- rbind(apply(hull_a$points, 2, min), apply(hull_a$points, 2, max))
  bounds_b <- rbind(apply(hull_b$points, 2, min), apply(hull_b$points, 2, max))
  volume_a_box <- prod(bounds_a[2, ] - bounds_a[1, ])
  volume_b_box <- prod(bounds_b[2, ] - bounds_b[1, ])
  
  selected_bounds <- if (volume_a_box <= volume_b_box) bounds_a else bounds_b
  lower_bounds <- selected_bounds[1, ]
  upper_bounds <- selected_bounds[2, ]
  ranges <- upper_bounds - lower_bounds
  bounding_volume <- prod(ranges)
  
  if (!is.finite(bounding_volume) || bounding_volume <= 0) return(NA_real_)
  
  random_points <- matrix(stats::runif(n_points * hull_a$dimensions), nrow = n_points, ncol = hull_a$dimensions)
  random_points <- sweep(random_points, 2, ranges, "*")
  random_points <- sweep(random_points, 2, lower_bounds, "+")
  
  inside_a <- tryCatch(geometry::inhulln(hull_a$hull, random_points), error = function(e) rep(FALSE, n_points))
  inside_b <- tryCatch(geometry::inhulln(hull_b$hull, random_points), error = function(e) rep(FALSE, n_points))
  
  overlap <- bounding_volume * mean(inside_a & inside_b)
  clamp_value(overlap, 0, min(hull_a$volume, hull_b$volume))
}

calculate_volume_metrics <- function(group_hulls, n_points) {
  purrr::map_dfr(
    names(comparison_definitions),
    function(comparison_name) {
      group_names <- comparison_definitions[[comparison_name]]
      hull_ref <- group_hulls[[group_names[1]]]
      hull_tgt <- group_hulls[[group_names[2]]]
      
      if (is.null(hull_ref) || is.null(hull_tgt)) return(tibble())
      
      overlap <- calculate_overlap_volume(hull_ref, hull_tgt, n_points)
      if (!is.finite(overlap)) return(tibble())
      
      ref_vol <- hull_ref$volume
      tgt_vol <- hull_tgt$volume
      union_vol <- ref_vol + tgt_vol - overlap
      
      if (!is.finite(union_vol) || union_vol <= 1e-14) return(tibble())
      
      loss <- clamp_value((ref_vol - overlap) / union_vol, 0, 1)
      gain <- clamp_value((tgt_vol - overlap) / union_vol, 0, 1)
      
      tibble(
        Comparison = comparison_name,
        Reference_group = group_names[1],
        Target_group = group_names[2],
        Reference_volume = ref_vol,
        Target_volume = tgt_vol,
        Overlap_volume = overlap,
        Union_volume = union_vol,
        Stable = clamp_value(overlap / union_vol, 0, 1),
        Loss = loss,
        Gain = gain,
        Turnover = loss + gain,
        Net = gain - loss
      )
    }
  )
}

calculate_centroid_statistics <- function(site_centroids, group_table) {
  selected_groups <- group_table %>% filter(Selected_for_TP_comparison)
  group_centroids <- lapply(tp_group_levels, function(g) {
    g_sites <- selected_groups$Site[selected_groups$TP_group == g]
    colMeans(site_centroids[g_sites, , drop = FALSE])
  })
  names(group_centroids) <- tp_group_levels
  
  group_sizes <- vapply(tp_group_levels, function(g) sum(selected_groups$TP_group == g), integer(1))
  
  pairwise <- purrr::imap_dfr(
    comparison_definitions,
    function(group_names, comp_name) {
      diff_vec <- group_centroids[[group_names[2]]] - group_centroids[[group_names[1]]]
      tibble(
        Comparison = comp_name,
        Reference_group = group_names[1],
        Target_group = group_names[2],
        Centroid_distance = sqrt(sum(diff_vec^2))
      )
    }
  )
  
  overall_centroid <- colMeans(site_centroids[selected_groups$Site, , drop = FALSE])
  omnibus <- sum(vapply(tp_group_levels, function(g) {
    diff_vec <- group_centroids[[g]] - overall_centroid
    group_sizes[[g]] * sum(diff_vec^2)
  }, numeric(1)))
  
  centroid_table <- bind_rows(lapply(tp_group_levels, function(g) {
    cur <- as.list(group_centroids[[g]])
    names(cur) <- colnames(site_centroids)
    as_tibble(cur) %>% mutate(TP_group = g, Sites_in_group = group_sizes[[g]], .before = 1)
  }))
  
  list(pairwise = pairwise, omnibus = omnibus, group_centroids = centroid_table)
}


# Randomly reassign complete multivariate trait-space coordinates among species.
# Community composition and TP groups remain unchanged.
# Entire coordinate vectors are permuted jointly, preserving correlations among axes
# and the geometry of the original functional trait space.
permute_trait_labels <- function(trait_space) {
  trait_space <- as.matrix(trait_space)
  if (nrow(trait_space) <= 1L) return(trait_space)

  perm_index <- sample.int(nrow(trait_space), size = nrow(trait_space), replace = FALSE)
  permuted_space <- trait_space[perm_index, , drop = FALSE]

  # Re-attach original species labels so that the randomly selected trait vector
  # is assigned to each fixed species identity in the community matrix.
  rownames(permuted_space) <- rownames(trait_space)
  colnames(permuted_space) <- colnames(trait_space)

  permuted_space
}

summarize_metric_against_null <- function(observed_val, null_vec, n_requested) {
  null_vec <- null_vec[is.finite(null_vec)]
  n_valid <- length(null_vec)
  
  if (!is.finite(observed_val) || n_valid < minimum_valid_draws) {
    return(tibble(
      Observed = observed_val, Null_mean = NA_real_, Null_sd = NA_real_,
      Null_q025 = NA_real_, Null_q975 = NA_real_, SES = NA_real_,
      Observed_minus_null = NA_real_, P_upper = NA_real_, P_lower = NA_real_,
      P_two_sided = NA_real_, Valid_draws = n_valid, Null_status = "insufficient_valid_draws"
    ))
  }
  
  null_m <- mean(null_vec)
  null_s <- stats::sd(null_vec)
  p_up <- (1 + sum(null_vec >= observed_val)) / (n_valid + 1)
  p_low <- (1 + sum(null_vec <= observed_val)) / (n_valid + 1)
  p_two <- min(1, 2 * min(p_up, p_low))
  
  if (!is.finite(null_s) || null_s <= null_sd_tolerance) {
    return(tibble(
      Observed = observed_val, Null_mean = null_m, Null_sd = null_s,
      Null_q025 = safe_quantile(null_vec, 0.025), Null_q975 = safe_quantile(null_vec, 0.975),
      SES = if (abs(observed_val - null_m) <= 1e-12) 0 else NA_real_,
      Observed_minus_null = observed_val - null_m, P_upper = p_up, P_lower = p_low,
      P_two_sided = p_two, Valid_draws = n_valid, Null_status = "degenerate_null_distribution"
    ))
  }
  
  tibble(
    Observed = observed_val,
    Null_mean = null_m,
    Null_sd = null_s,
    Null_q025 = safe_quantile(null_vec, 0.025),
    Null_q975 = safe_quantile(null_vec, 0.975),
    SES = (observed_val - null_m) / null_s,
    Observed_minus_null = observed_val - null_m,
    P_upper = p_up,
    P_lower = p_low,
    P_two_sided = p_two,
    Valid_draws = n_valid,
    Null_status = "ok"
  )
}

assign_tp_groups <- function(site_tp, randomize_ties = FALSE) {
  site_tp <- site_tp %>%
    transmute(Site = normalize_identifier(Site), TP = as.numeric(TP)) %>%
    filter(!is.na(Site), is.finite(TP))
  
  n_sites <- nrow(site_tp)
  sites_per_group <- floor(n_sites * tp_group_fraction)
  
  if (sites_per_group < minimum_sites_per_group || 3L * sites_per_group > n_sites) {
    return(NULL)
  }
  
  tie_breaker <- if (randomize_ties) stats::runif(n_sites) else rank(site_tp$Site, ties.method = "first")
  sorted <- site_tp[order(site_tp$TP, tie_breaker), , drop = FALSE]
  
  middle_start <- floor((n_sites - sites_per_group) / 2) + 1L
  low_idx <- seq_len(sites_per_group)
  mid_idx <- middle_start:(middle_start + sites_per_group - 1L)
  high_idx <- (n_sites - sites_per_group + 1L):n_sites
  
  sorted$TP_group <- NA_character_
  sorted$TP_group[low_idx] <- "Low"
  sorted$TP_group[mid_idx] <- "Mid"
  sorted$TP_group[high_idx] <- "High"
  
  sorted %>%
    mutate(
      TP_group = factor(TP_group, levels = tp_group_levels),
      Selected_for_TP_comparison = !is.na(TP_group)
    ) %>%
    arrange(match(Site, site_tp$Site))
}

# ==============================================================================
# 4. Integrated Dataset Processing Function
# ==============================================================================

process_dataset_integrated <- function(file, trait_map) {
  dataset <- dataset_name_from_file(file)
  organism <- get_organism(file)
  message("Processing: ", dataset, " (", organism, ")")
  
  sheet_names <- readxl::excel_sheets(file)
  required_sheets <- c("species", "traits", "environment")
  missing_sheets <- setdiff(required_sheets, sheet_names)
  
  if (length(missing_sheets) > 0L) {
    return(list(
      vol_results = tibble(), vol_draws = tibble(),
      centroid_results = tibble(), centroid_draws = tibble(),
      omnibus_results = tibble(), site_centroids = tibble(),
      tp_groups = tibble(), group_centroids = tibble(),
      log = empty_log_row(dataset, organism, status = "failed", message = "Missing sheets")
    ))
  }
  
  community_raw <- readxl::read_excel(file, sheet = "species")
  traits_raw <- readxl::read_excel(file, sheet = "traits")
  environment_raw <- readxl::read_excel(file, sheet = "environment")
  
  names(community_raw)[1] <- "Site"
  names(traits_raw)[1] <- "Species"
  names(environment_raw)[1] <- "Site"
  
  community_raw$Site <- normalize_identifier(community_raw$Site)
  traits_raw$Species <- normalize_identifier(traits_raw$Species)
  environment_raw$Site <- normalize_identifier(environment_raw$Site)
  
  community_numeric <- community_raw %>%
    select(-Site) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(as.character(.x)))))
  community_matrix <- as.matrix(community_numeric)
  community_matrix[!is.finite(community_matrix)] <- 0
  community_matrix <- 1L * (community_matrix > 0)
  rownames(community_matrix) <- community_raw$Site
  
  environment_clean <- environment_raw %>%
    transmute(Site = normalize_identifier(Site), TP = suppressWarnings(as.numeric(TP))) %>%
    filter(!is.na(Site), is.finite(TP), Site %in% rownames(community_matrix)) %>%
    distinct(Site, .keep_all = TRUE)
  
  organism_mapping <- trait_map %>% filter(Organism == organism)
  
  vol_res_list <- list()
  vol_draws_list <- list()
  cwm_res_list <- list()
  cwm_draws_list <- list()
  omni_res_list <- list()
  site_cnt_list <- list()
  tp_grp_list <- list()
  grp_cnt_list <- list()
  log_list <- list()
  
  for (dimension_name in names(broad_dims)) {
    message("  Dimension: ", dimension_name)
    
    requested_traits <- organism_mapping %>%
      filter(Trait_category %in% broad_dims[[dimension_name]]) %>%
      pull(Original_trait) %>%
      unique()
    
    actual_traits <- intersect(requested_traits, names(traits_raw))
    if (length(actual_traits) == 0L) {
      log_list[[dimension_name]] <- empty_log_row(dataset, organism, dimension_name, "skipped", "No mapped traits")
      next
    }
    
    traits_dim <- traits_raw %>%
      select(Species, all_of(actual_traits)) %>%
      filter(!is.na(Species)) %>%
      mutate(across(all_of(actual_traits), clean_trait_column)) %>%
      distinct(Species, .keep_all = TRUE)
    
    inf_traits <- actual_traits[vapply(traits_dim[actual_traits], is_informative_trait, logical(1))]
    if (length(inf_traits) == 0L) {
      log_list[[dimension_name]] <- empty_log_row(dataset, organism, dimension_name, "skipped", "No informative traits")
      next
    }
    
    traits_dim <- traits_dim %>%
      select(Species, all_of(inf_traits)) %>%
      filter(if_all(all_of(inf_traits), ~ !is.na(.x)))
    
    inf_traits <- inf_traits[vapply(traits_dim[inf_traits], is_informative_trait, logical(1))]
    if (length(inf_traits) == 0L) {
      log_list[[dimension_name]] <- empty_log_row(dataset, organism, dimension_name, "skipped", "Constant traits after filtering")
      next
    }
    
    initial_species <- intersect(colnames(community_matrix), traits_dim$Species)
    if (length(initial_species) < minimum_species_pool) {
      log_list[[dimension_name]] <- empty_log_row(dataset, organism, dimension_name, "skipped", "Too few trait-complete species")
      next
    }
    
    comm_initial <- community_matrix[environment_clean$Site, initial_species, drop = FALSE]
    usable_richness <- rowSums(comm_initial)
    eligible_sites <- names(usable_richness)[usable_richness >= minimum_site_richness]
    
    site_tp <- environment_clean %>%
      filter(Site %in% eligible_sites) %>%
      arrange(match(Site, eligible_sites))
    
    observed_groups <- assign_tp_groups(site_tp, randomize_ties = FALSE)
    if (is.null(observed_groups)) {
      log_list[[dimension_name]] <- empty_log_row(dataset, organism, dimension_name, "skipped", "TP grouping failed")
      next
    }
    
    sites_per_group <- sum(observed_groups$TP_group == "Low", na.rm = TRUE)
    comm_eligible <- community_matrix[site_tp$Site, initial_species, drop = FALSE]
    species_pool <- colnames(comm_eligible)[colSums(comm_eligible) > 0]
    
    if (length(species_pool) < minimum_species_pool) {
      log_list[[dimension_name]] <- empty_log_row(dataset, organism, dimension_name, "skipped", "Too few occurring species")
      next
    }
    
    trait_data <- traits_dim %>%
      filter(Species %in% species_pool) %>%
      slice(match(species_pool, Species)) %>%
      select(all_of(inf_traits)) %>%
      as.data.frame(check.names = FALSE)
    rownames(trait_data) <- species_pool
    
    comm_dim <- comm_eligible[, species_pool, drop = FALSE]
    
    obs_species_groups <- purrr::map(tp_group_levels, function(g) {
      g_sites <- observed_groups$Site[observed_groups$TP_group == g & !is.na(observed_groups$TP_group)]
      pres <- colSums(comm_dim[g_sites, , drop = FALSE]) > 0
      names(pres)[pres]
    })
    names(obs_species_groups) <- tp_group_levels
    
    gower_dist <- cluster::daisy(trait_data, metric = "gower")
    pcoa_res <- ape::pcoa(gower_dist, correction = "cailliez")
    coords <- pcoa_res$vectors.cor
    if (is.null(coords) || ncol(coords) == 0L) coords <- pcoa_res$vectors
    coords <- as.matrix(coords)
    rownames(coords) <- rownames(trait_data)
    
    valid_axes <- which(apply(coords, 2, function(x) all(is.finite(x)) && stats::sd(x) > null_sd_tolerance))
    if (length(valid_axes) < minimum_axes) {
      log_list[[dimension_name]] <- empty_log_row(dataset, organism, dimension_name, "failed", "No valid PCoA axes")
      next
    }
    
    axes_used <- valid_axes[seq_len(min(length(valid_axes), n_axes_max))]
    trait_space <- coords[, axes_used, drop = FALSE]
    colnames(trait_space) <- paste0("Axis", seq_len(ncol(trait_space)))
    
    all_axis_var <- sum(apply(coords[, valid_axes, drop = FALSE], 2, var))
    used_axis_var <- sum(apply(trait_space, 2, var))
    pcoa_var_retained <- if (is.finite(all_axis_var) && all_axis_var > 0) used_axis_var / all_axis_var else NA_real_
    
    site_richness <- rowSums(comm_dim)
    site_centroids <- (comm_dim %*% trait_space) / site_richness
    site_centroids[!is.finite(site_centroids)] <- NA_real_
    
    obs_group_hulls <- purrr::map(obs_species_groups, build_hull, trait_space = trait_space)
    obs_vol_metrics <- calculate_volume_metrics(obs_group_hulls, n_points = n_mc)
    obs_centroid_stats <- calculate_centroid_statistics(site_centroids, observed_groups)
    
    dim_seed <- seed_from_text(random_seed, paste(dataset, dimension_name))
    set.seed(dim_seed)
    
    species_inc <- colSums(comm_dim)[rownames(trait_space)]
    sampling_prob <- if (null_model == "r1") species_inc / sum(species_inc) else rep(1 / nrow(trait_space), nrow(trait_space))
    
    null_vol_list <- vector("list", n_iterations)
    null_pairwise_cnt_list <- vector("list", n_iterations)
    null_omni_cnt_list <- vector("list", n_iterations)
    
    for (iter in seq_len(n_iterations)) {
      null_groups <- lapply(obs_species_groups, function(o_set) {
        sample(rownames(trait_space), size = length(o_set), replace = FALSE, prob = sampling_prob)
      })
      null_hulls <- purrr::map(null_groups, build_hull, trait_space = trait_space)
      null_vol <- calculate_volume_metrics(null_hulls, n_points = n_mc)
      if (nrow(null_vol) > 0L) {
        null_vol_list[[iter]] <- null_vol %>% mutate(Iteration = iter, .before = 1)
      }
      
      # Panel b null model: keep observed community composition and TP groups fixed,
      # but randomly reassign complete multivariate trait-space coordinates among species.
      # This preserves site richness, species composition, taxonomic turnover and TP
      # grouping, while breaking the species-trait association.
      permuted_trait_space <- permute_trait_labels(trait_space)

      permuted_site_centroids <- (comm_dim %*% permuted_trait_space) / site_richness
      permuted_site_centroids[!is.finite(permuted_site_centroids)] <- NA_real_

      perm_stats <- calculate_centroid_statistics(
        permuted_site_centroids,
        observed_groups
      )
      
      null_pairwise_cnt_list[[iter]] <- perm_stats$pairwise %>%
        transmute(Iteration = iter, Comparison, Null_centroid_distance = Centroid_distance)
      null_omni_cnt_list[[iter]] <- tibble(Iteration = iter, Null_omnibus_separation = perm_stats$omnibus)
    }
    
    null_vol_all <- bind_rows(null_vol_list)
    null_pairwise_cnt_all <- bind_rows(null_pairwise_cnt_list)
    null_omni_cnt_all <- bind_rows(null_omni_cnt_list)
    
    vol_summary <- purrr::map_dfr(comparison_levels, function(comp) {
      obs_row <- obs_vol_metrics %>% filter(Comparison == comp)
      purrr::map_dfr(metric_levels, function(met) {
        obs_val <- if (nrow(obs_row) == 1L) obs_row[[met]][1] else NA_real_
        null_vec <- null_vol_all %>% filter(Comparison == comp) %>% pull(all_of(met))
        summarize_metric_against_null(obs_val, null_vec, n_iterations) %>%
          mutate(Comparison = comp, Metric = met, .before = 1)
      })
    }) %>%
      mutate(Dataset = dataset, Organism = organism, Dimension = dimension_name, PCoA_axes = ncol(trait_space), .before = 1)
    
    centroid_summary <- purrr::map_dfr(comparison_levels, function(comp) {
      obs_val <- obs_centroid_stats$pairwise %>% filter(Comparison == comp) %>% pull(Centroid_distance)
      null_vec <- null_pairwise_cnt_all %>% filter(Comparison == comp) %>% pull(Null_centroid_distance)
      summarize_metric_against_null(obs_val, null_vec, n_iterations) %>%
        mutate(Comparison = comp, Metric = "Centroid_distance", .before = 1)
    }) %>%
      mutate(Dataset = dataset, Organism = organism, Dimension = dimension_name, .before = 1)
    
    omni_summary <- summarize_metric_against_null(obs_centroid_stats$omnibus, null_omni_cnt_all$Null_omnibus_separation, n_iterations) %>%
      mutate(Dataset = dataset, Organism = organism, Dimension = dimension_name, Statistic = "Omnibus separation", .before = 1)
    
    vol_res_list[[dimension_name]] <- vol_summary
    vol_draws_list[[dimension_name]] <- null_vol_all %>% mutate(Dataset = dataset, Organism = organism, Dimension = dimension_name, .before = 1)
    cwm_res_list[[dimension_name]] <- centroid_summary
    cwm_draws_list[[dimension_name]] <- null_pairwise_cnt_all %>% mutate(Dataset = dataset, Organism = organism, Dimension = dimension_name, .before = 1)
    omni_res_list[[dimension_name]] <- omni_summary
    
    site_cnt_out <- as_tibble(site_centroids, rownames = "Site") %>%
      mutate(Dimension_richness = site_richness[Site]) %>%
      left_join(observed_groups %>% select(Site, TP, TP_group, Selected_for_TP_comparison), by = "Site") %>%
      mutate(Dataset = dataset, Organism = organism, Dimension = dimension_name, .before = 1)
    for (ax in paste0("Axis", seq_len(n_axes_max))) if (!ax %in% names(site_cnt_out)) site_cnt_out[[ax]] <- NA_real_
    site_cnt_list[[dimension_name]] <- site_cnt_out
    
    grp_cnt_out <- obs_centroid_stats$group_centroids %>%
      mutate(Dataset = dataset, Organism = organism, Dimension = dimension_name, .before = 1)
    for (ax in paste0("Axis", seq_len(n_axes_max))) if (!ax %in% names(grp_cnt_out)) grp_cnt_out[[ax]] <- NA_real_
    grp_cnt_list[[dimension_name]] <- grp_cnt_out
    
    tp_grp_list[[dimension_name]] <- observed_groups %>%
      mutate(Dataset = dataset, Organism = organism, Dimension = dimension_name, .before = 1)
    
    log_list[[dimension_name]] <- empty_log_row(dataset, organism, dimension_name, "ok", "Completed") %>%
      mutate(
        Raw_sites = nrow(community_matrix),
        Eligible_sites = nrow(site_tp),
        Sites_per_TP_group = sites_per_group,
        Traits_requested = length(requested_traits),
        Traits_used = length(inf_traits),
        Trait_names_used = paste(inf_traits, collapse = " | "),
        Species_in_PCoA = nrow(trait_space),
        PCoA_axes = ncol(trait_space),
        PCoA_variance_retained = pcoa_var_retained
      )
  }
  
  list(
    vol_results = bind_rows(vol_res_list),
    vol_draws = bind_rows(vol_draws_list),
    centroid_results = bind_rows(cwm_res_list),
    centroid_draws = bind_rows(cwm_draws_list),
    omnibus_results = bind_rows(omni_res_list),
    site_centroids = bind_rows(site_cnt_list),
    tp_groups = bind_rows(tp_grp_list),
    group_centroids = bind_rows(grp_cnt_list),
    log = bind_rows(log_list)
  )
}

# ==============================================================================
# 5. Execution & Strict Unified Cohort Alignment
# ==============================================================================

if (run_calculation) {
  if (!dir.exists(raw_data_path)) stop("Raw-data directory not found.")
  trait_map <- read_trait_mapping(mapping_file)
  dataset_files <- list.files(raw_data_path, pattern = "\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
  dataset_files <- dataset_files[!grepl("^~\\$", basename(dataset_files))]
  
  all_results <- list()
  for (file in dataset_files) {
    dataset <- dataset_name_from_file(file)
    checkpoint_file <- file.path(checkpoint_dir, paste0(make.names(dataset), "_", settings_tag, ".rds"))
    
    if (resume_from_checkpoints && file.exists(checkpoint_file)) {
      message("Loading checkpoint: ", dataset)
      current_res <- readRDS(checkpoint_file)
    } else {
      current_res <- tryCatch(
        process_dataset_integrated(file, trait_map),
        error = function(e) {
          organism <- get_organism(file)
          list(
            vol_results = tibble(), vol_draws = tibble(), centroid_results = tibble(),
            centroid_draws = tibble(), omnibus_results = tibble(), site_centroids = tibble(),
            tp_groups = tibble(), group_centroids = tibble(),
            log = empty_log_row(dataset, organism, status = "failed", message = conditionMessage(e))
          )
        }
      )
      saveRDS(current_res, checkpoint_file)
    }
    all_results[[dataset]] <- current_res
  }
  
  vol_results_all <- bind_rows(lapply(all_results, `[[`, "vol_results"))
  centroid_results_all <- bind_rows(lapply(all_results, `[[`, "centroid_results"))
  omnibus_results_all <- bind_rows(lapply(all_results, `[[`, "omnibus_results"))
  site_centroids_all <- bind_rows(lapply(all_results, `[[`, "site_centroids"))
  tp_groups_all <- bind_rows(lapply(all_results, `[[`, "tp_groups"))
  group_centroids_all <- bind_rows(lapply(all_results, `[[`, "group_centroids"))
  processing_log_all <- bind_rows(lapply(all_results, `[[`, "log"))
  
  # Common cohort filtering: ensures identical K across panels
  inclusion <- processing_log_all %>%
    filter(Processing_status == "ok") %>%
    distinct(Dataset, Organism, Dimension) %>%
    left_join(
      vol_results_all %>%
        group_by(Dataset, Organism, Dimension) %>%
        summarise(
          Vol_valid_SES = sum(Null_status == "ok" & is.finite(SES) & Metric == "Net"),
          .groups = "drop"
        ),
      by = c("Dataset", "Organism", "Dimension")
    ) %>%
    left_join(
      centroid_results_all %>%
        group_by(Dataset, Organism, Dimension) %>%
        summarise(
          Centroid_valid_SES = sum(Null_status == "ok" & is.finite(SES)),
          .groups = "drop"
        ),
      by = c("Dataset", "Organism", "Dimension")
    ) %>%
    mutate(
      Vol_valid_SES = replace_na(Vol_valid_SES, 0L),
      Centroid_valid_SES = replace_na(Centroid_valid_SES, 0L),
      Included_unified_analysis = (Vol_valid_SES == length(comparison_levels) &
                                     Centroid_valid_SES == length(comparison_levels)),
      Exclusion_reason = case_when(
        Included_unified_analysis ~ NA_character_,
        Vol_valid_SES < length(comparison_levels) ~ "Volume net change null distribution incomplete",
        Centroid_valid_SES < length(comparison_levels) ~ "Centroid displacement trait-permutation incomplete",
        TRUE ~ "Not included"
      )
    )
  
  valid_dataset_dimensions <- inclusion %>%
    filter(Included_unified_analysis) %>%
    select(Dataset, Organism, Dimension)
  
  vol_results_complete <- vol_results_all %>% inner_join(valid_dataset_dimensions, by = c("Dataset", "Organism", "Dimension"))
  centroid_results_complete <- centroid_results_all %>% inner_join(valid_dataset_dimensions, by = c("Dataset", "Organism", "Dimension"))
  
  vol_results_global <- vol_results_complete %>%
    group_by(Dimension, Comparison, Metric) %>%
    summarise(
      K_datasets = n_distinct(Dataset),
      Mean_observed = mean(Observed, na.rm = TRUE),
      Mean_null = mean(Null_mean, na.rm = TRUE),
      Mean_SES = mean(SES, na.rm = TRUE),
      SD_SES = stats::sd(SES, na.rm = TRUE),
      SE_SES = SD_SES / sqrt(K_datasets),
      Median_SES = stats::median(SES, na.rm = TRUE),
      Q25_SES = safe_quantile(SES, 0.25),
      Q75_SES = safe_quantile(SES, 0.75),
      Proportion_p_lt_0_05 = mean(P_two_sided < 0.05, na.rm = TRUE),
      P_mean_SES = safe_one_sample_test(SES),
      .groups = "drop"
    ) %>%
    mutate(
      T_critical = stats::qt(0.975, df = K_datasets - 1L),
      CI95_SES_low = Mean_SES - T_critical * SE_SES,
      CI95_SES_high = Mean_SES + T_critical * SE_SES,
      P_mean_SES_FDR = p.adjust(P_mean_SES, method = "BH")
    )
  
  centroid_results_global <- centroid_results_complete %>%
    group_by(Dimension, Comparison) %>%
    summarise(
      K_datasets = n_distinct(Dataset),
      Mean_observed_distance = mean(Observed, na.rm = TRUE),
      Mean_null_distance = mean(Null_mean, na.rm = TRUE),
      Mean_SES = mean(SES, na.rm = TRUE),
      SD_SES = stats::sd(SES, na.rm = TRUE),
      SE_SES = SD_SES / sqrt(K_datasets),
      Median_SES = stats::median(SES, na.rm = TRUE),
      Q25_SES = safe_quantile(SES, 0.25),
      Q75_SES = safe_quantile(SES, 0.75),
      Proportion_upper_p_lt_0_05 = mean(P_upper < 0.05, na.rm = TRUE),
      P_mean_SES = safe_one_sample_test(SES),
      .groups = "drop"
    ) %>%
    mutate(
      T_critical = stats::qt(0.975, df = K_datasets - 1L),
      CI95_SES_low = Mean_SES - T_critical * SE_SES,
      CI95_SES_high = Mean_SES + T_critical * SE_SES,
      P_mean_SES_FDR = p.adjust(P_mean_SES, method = "BH")
    )
  
  export_list <- list(
    Volume_complete = vol_results_complete %>% mutate(across(where(is.factor), as.character)),
    Volume_global = vol_results_global %>% mutate(across(where(is.factor), as.character)),
    Centroid_complete = centroid_results_complete %>% mutate(across(where(is.factor), as.character)),
    Centroid_global = centroid_results_global %>% mutate(across(where(is.factor), as.character)),
    Omnibus_results = omnibus_results_all %>% mutate(across(where(is.factor), as.character)),
    Site_centroids = site_centroids_all %>% mutate(across(where(is.factor), as.character)),
    Group_centroids = group_centroids_all %>% mutate(across(where(is.factor), as.character)),
    TP_groups = tp_groups_all %>% mutate(across(where(is.factor), as.character)),
    Inclusion = inclusion %>% mutate(across(where(is.factor), as.character)),
    Processing_log = processing_log_all %>% mutate(across(where(is.factor), as.character))
  )
  writexl::write_xlsx(export_list, output_file)
  
  if (save_draws_rds) {
    vol_draws_all <- bind_rows(lapply(all_results, `[[`, "vol_draws"))
    centroid_draws_all <- bind_rows(lapply(all_results, `[[`, "centroid_draws"))
    saveRDS(list(vol_draws = vol_draws_all, centroid_draws = centroid_draws_all), draws_rds_file)
  }
} else {
  if (!file.exists(output_file)) stop("Saved result workbook not found.")
  vol_results_complete <- readxl::read_excel(output_file, sheet = "Volume_complete")
  vol_results_global <- readxl::read_excel(output_file, sheet = "Volume_global")
  centroid_results_complete <- readxl::read_excel(output_file, sheet = "Centroid_complete")
  centroid_results_global <- readxl::read_excel(output_file, sheet = "Centroid_global")
}

# ==============================================================================
# 6. Combined Visualization (Pure Forest Plots: Tagged via plot_annotation)
# ==============================================================================

vol_results_complete <- vol_results_complete %>%
  mutate(
    Dimension = factor(Dimension, levels = dimension_levels),
    Comparison = factor(Comparison, levels = comparison_levels)
  )

vol_results_global <- vol_results_global %>%
  mutate(
    Dimension = factor(Dimension, levels = dimension_levels),
    Comparison = factor(Comparison, levels = comparison_levels)
  )

centroid_results_complete <- centroid_results_complete %>%
  mutate(
    Dimension = factor(Dimension, levels = dimension_levels),
    Comparison = factor(Comparison, levels = comparison_levels)
  )

centroid_results_global <- centroid_results_global %>%
  mutate(
    Dimension = factor(Dimension, levels = dimension_levels),
    Comparison = factor(Comparison, levels = comparison_levels)
  )

# Calculate unified K labels
dimension_k <- vol_results_complete %>%
  transmute(Dataset, Dimension = as.character(Dimension)) %>%
  distinct() %>%
  count(Dimension, name = "K_datasets")

dimension_legend_labels <- vapply(
  dimension_levels,
  function(d) {
    cur_k <- dimension_k$K_datasets[dimension_k$Dimension == d]
    if (length(cur_k) == 0L) cur_k <- 0L
    paste0(d, " (K = ", cur_k[1], ")")
  },
  character(1)
)

dimension_colors <- c(Habitat = "#1B9E77", Resource = "#E67E22", Size = "#4C78A8")

# Unified theme: left-aligned titles and clean forest plot structure
base_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 15, hjust = 0.5, margin = margin(b = 10)),
    axis.title = element_text(size = 14),
    axis.title.x = element_blank(), 
    axis.text = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    plot.margin = margin(10, 14, 10, 10)
  )

# ------------------------------------------------------------------------------
# Panel a: Net Trait-Space Change (Forest Plot Only: Y in [-2, 2])
# ------------------------------------------------------------------------------
net_plot_data_global <- vol_results_global %>% filter(Metric == "Net")

plot_a_net <- ggplot(
  net_plot_data_global,
  aes(x = Comparison, y = Mean_SES, colour = Dimension)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.55) +
  
  geom_point(
    data = vol_results_complete %>% filter(Metric == "Net"),
    aes(x = Comparison, y = SES, colour = Dimension),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.55),
    alpha = 0.3,      
    size = 1.5,      
    shape = 16,
    show.legend = FALSE 
  ) +
  # ============================

geom_errorbar(
  aes(ymin = CI95_SES_low, ymax = CI95_SES_high),
  width = 0,
  linewidth = 0.95,
  position = position_dodge(width = 0.55)
) +
  geom_point(
    shape = 18,
    size = 4.2,
    position = position_dodge(width = 0.55)
  ) +
  coord_cartesian(ylim = c(-2, 2)) +
  scale_y_continuous(breaks = seq(-2, 2, by = 1)) +
  scale_colour_manual(
    values = dimension_colors,
    breaks = dimension_levels,
    labels = dimension_legend_labels,
    drop = FALSE
  ) +
  labs(
    title = "Functional trait space",
    y = "Net trait space change (SES)",
    colour = "Trait dimension"
  ) +
  base_theme

# ------------------------------------------------------------------------------
# Panel b: Centroid Displacement (Forest Plot Only: Y in [-3, 3])
# ------------------------------------------------------------------------------
plot_b_centroid <- ggplot(
  centroid_results_global,
  aes(x = Comparison, y = Mean_SES, colour = Dimension)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.55) +
  
  geom_point(
    data = centroid_results_complete,
    aes(x = Comparison, y = SES, colour = Dimension),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.55),
    alpha = 0.3,
    size = 1.5,
    shape = 16,
    show.legend = FALSE
  ) +
  # ============================

geom_errorbar(
  aes(ymin = CI95_SES_low, ymax = CI95_SES_high),
  width = 0,
  linewidth = 0.95,
  position = position_dodge(width = 0.55)
) +
  geom_point(
    shape = 18,
    size = 4.2,
    position = position_dodge(width = 0.55)
  ) +
  coord_cartesian(ylim = c(-5, 5)) +
  scale_y_continuous(breaks = seq(-5, 5, by = 2.5)) +
  scale_colour_manual(
    values = dimension_colors,
    breaks = dimension_levels,
    labels = dimension_legend_labels,
    drop = FALSE
  ) +
  labs(
    title = "Community functional centroid",
    y = "Centroid displacement (SES)",
    colour = "Trait dimension"
  ) +
  base_theme
# ------------------------------------------------------------------------------
# Combine a and b using patchwork with plot_annotation(tag_levels = "a")
# ------------------------------------------------------------------------------
combined_figure <- (plot_a_net + plot_b_centroid) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(size = 18)
  )

print(combined_figure)

# Save figures
ggsave(file.path(output_dir, "Figure5.png"), combined_figure, width = 11, height = 5, dpi = 300)
