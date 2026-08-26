library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(readxl)
library(patchwork)

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

out_dir <- file.path("..", "000_output_file", "Figure2")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# --- Figure 2a ---
dat_var <- read.csv(file.path("..", "000_input_file", "dataset_var.csv"))

dat_rob <- dat_var %>%
  filter(!is.na(X), !is.na(Y)) %>%
  st_as_sf(coords = c("X", "Y"), crs = 4326) %>%
  st_transform(crs = "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs")

world_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(continent != "Antarctica") %>%
  st_transform(crs = "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs")

color_palette <- c("Producer" = "#1f77b4", "Consumer" = "#ff7f0e")
shape_palette <- c("Lake" = 17, "River" = 16)

Fig2a <- ggplot() +
  geom_sf(data = world_map, fill = "grey95", color = "grey80", linewidth = 0.2) +
  geom_sf(data = dat_rob, aes(color = Group, shape = Ecosystem), size = 2.3, alpha = 1) +
  scale_color_manual(values = color_palette, name = "Organism") +
  scale_shape_manual(values = shape_palette, name = "Ecosystem") +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = c(0.15, 0.25),
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 7),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.3, "cm"),       
    legend.spacing.y = unit(0.05, "cm"),     
    legend.margin = margin(0, 0, 0, 0),     
    legend.box.margin = margin(-5, -5, -5, -5)
  ) +
  guides(
    color = guide_legend(byrow = TRUE),
    shape = guide_legend(byrow = TRUE)
  )

ggsave(file.path(out_dir, "Fig2a.png"), plot = Fig2a, width = 7, height = 4.5, units = "in", dpi = 600, bg = "white")

# --- Figure 2b ---
dat <- read_xlsx(file.path("..", "000_input_file", "All_data.xlsx")) %>%
  mutate(logTP = log(TP + 1))

dat_dataset <- dat %>%
  distinct(Dataset, Ecosystem, Trophic_level)

eco_summary <- dat_dataset %>%
  count(Category = Ecosystem) %>%
  mutate(Type = "Ecosystem")

trophic_summary <- dat_dataset %>%
  count(Category = Trophic_level) %>%
  mutate(Type = "Trophic level")

combined_summary <- bind_rows(eco_summary, trophic_summary) %>%
  mutate(
    Category = factor(Category, levels = rev(c("River", "Lake", "Producer", "Consumer"))),
    Type = factor(Type, levels = c("Ecosystem", "Trophic level"))
  )

fig2b <- ggplot(combined_summary, aes(x = n, y = Category, fill = Type)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = n), hjust = -0.15, size = 3.5, color = "black") +
  facet_grid(Type ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = c("Ecosystem" = "#547ac0", "Trophic level" = "steelblue3")) +
  labs(x = "Number of datasets", y = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  expand_limits(x = max(combined_summary$n) * 1.15)

ggsave(file.path(out_dir, "Fig2b.png"), plot = fig2b, width = 4.5, height = 4, units = "in", dpi = 600, bg = "white")

# --- Figure 2c ---
dat_plot <- dat %>%
  mutate(Landuse_class = factor(Landuse_class, levels = c("Human", "Mixed", "Forest")))

custom_theme <- theme_minimal() +
  theme(
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 8),
    axis.title = element_text(size = 9),
    plot.background = element_rect(fill = "white", color = NA)
  )

bar_color <- "#1f77b4"

p1 <- ggplot(dat_plot, aes(x = logTP)) +
  geom_histogram(bins = 12, fill = bar_color, color = "white", linewidth = 0.2) +
  labs(x = "Total phosphorus (log)", y = "Number of sites") +
  custom_theme

p2 <- ggplot(dat_plot, aes(x = bio1_1500)) +
  geom_histogram(bins = 12, fill = bar_color, color = "white", linewidth = 0.2) +
  labs(x = "Annual mean temperature (°C)", y = "Number of sites") +
  custom_theme

p3 <- ggplot(dat_plot, aes(x = bio12_1500)) +
  geom_histogram(bins = 12, fill = bar_color, color = "white", linewidth = 0.2) +
  labs(x = "Annual mean precipitation (mm)", y = "Number of sites") +
  custom_theme

p4 <- ggplot(dat_plot, aes(x = Landuse_class)) +
  geom_bar(fill = bar_color, color = "white", linewidth = 0.2, width = 0.55) +
  labs(x = "Land-use class", y = "Number of sites") +
  custom_theme

Fig2c <- (p1 + p2 + p3 + p4) + 
  plot_layout(nrow = 1, widths = c(1, 1, 1, 0.95)) & 
  theme(plot.margin = margin(10, 5, 10, 5))

ggsave(file.path(out_dir, "Fig2c.png"), plot = Fig2c, width = 11, height = 3, units = "in", dpi = 600, bg = "white")