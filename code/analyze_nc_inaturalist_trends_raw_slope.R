# North Carolina iNaturalist bird trends vs. BBS, 2015-2025
#
# Put this script and the three CSV files in the same folder, set that folder
# as the R working directory, then run:
#   source("analyze_nc_inaturalist_trends_raw_slope.R")
# Requires tidyverse: install.packages("tidyverse")
# Outputs are saved under ./nc_inat_trend_outputs_raw_slope/.

library(tidyverse)

setwd("C:/Z_School/school/HurlbertLab/lab meeting")


start_year <- 2015L
end_year <- 2025L
years_expected <- start_year:end_year
output_dir <- "nc_inat_trend_outputs_NEW"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

birds <- read_csv("north_carolina_bird_observations_by_year.csv", show_col_types = FALSE)
stats <- read_csv("nc_inat_stats_2015_2025.csv", show_col_types = FALSE)
bbs <- read_csv("north_carolina_bbs_2015_2025_trend.csv", show_col_types = FALSE) %>%
  transmute(Species, BBS_APC_percent = as.numeric(Trend))

bird_species <- setdiff(names(birds), "year")

# Each dplyr::selected denominator is taken from the SAME year as the species count.
# The proxy is observations made by users with 10+ bird observations; the file
# does not contain the number of those users. bird_observers is the available
# annual bird-contributor count; a separate lifetime-ever count is not present.
metric_denominators <- c(
  total_iNat_observations = "total_observations",
  total_iNat_users = "total_inat_users",
  total_bird_observations = "bird_observations",
  annual_bird_observers = "bird_observers",
  superuser_observation_proxy =
    "observations_by_users_with_10_or_more_bird_observations"
)

missing_stats <- setdiff(unique(unname(metric_denominators)), names(stats))
if (length(missing_stats) > 0) {
  stop("Missing stats columns: ", paste(missing_stats, collapse = ", "))
}
if (!all(years_expected %in% birds$year) || !all(years_expected %in% stats$year)) {
  stop("Both input files must include every year from 2015 through 2025.")
}
missing_bbs <- setdiff(bird_species, bbs$Species)
if (length(missing_bbs) > 0) {
  stop("No BBS trend found for: ", paste(missing_bbs, collapse = ", "))
}

# A row per species-year with the matching annual iNaturalist denominators.
annual_counts <- birds %>%
  filter(year %in% years_expected) %>%
  pivot_longer(
    cols = all_of(bird_species), names_to = "Species", values_to = "inat_count"
  ) %>%
  inner_join(stats %>% filter(year %in% years_expected), by = "year")

# Calculate each fraction separately for each metric and species-year:
# standardized_rate = species' bird count in year / dplyr::selected metric in that year.
standardized <- imap_dfr(metric_denominators, function(denom_col, metric_name) {
  annual_counts %>%
    transmute(
      year, Species, inat_count,
      standardization = metric_name,
      denominator_name = denom_col,
      denominator_value = .data[[denom_col]],
      standardized_rate = if_else(
        is.na(.data[[denom_col]]) | .data[[denom_col]] <= 0,
        NA_real_, inat_count / .data[[denom_col]]
      )
    )
})

# Save annual standardized fractions in wide form: one column per metric.
standardized_wide <- standardized %>%
  dplyr::select(year, Species, standardized_rate, standardization) %>%
  pivot_wider(names_from = standardization, values_from = standardized_rate) %>%
  arrange(Species, year)
write_csv(standardized_wide,
          file.path(output_dir, "standardized_species_counts_2015_2025.csv"))

# Fit a raw-scale straight-line regression for each species x standardization:
#   standardized_rate = intercept + slope * (year - 2015)
# No log transformation is used. The slope is an absolute change in standardized
# fraction per year, not a percent change. Require all 11 annual fractions.
trend_fits <- standardized %>%
  group_by(Species, standardization) %>%
  group_modify(~ {
    d <- .x %>% filter(year %in% years_expected) %>% arrange(year)
    if (nrow(d) != length(years_expected) ||
        !identical(as.integer(d$year), as.integer(years_expected)) ||
        any(is.na(d$standardized_rate))) {
      stop(paste0(
        "Cannot fit the full-period raw-scale trend for ", .y$Species,
        " / ", .y$standardization,
        ": expected 11 non-missing annual rates for 2015-2025."
      ))
    }
    fit <- lm(standardized_rate ~ I(year - start_year), data = d)
    tibble(
      iNat_raw_slope_per_year = unname(coef(fit)[["I(year - start_year)"]]),
      iNat_raw_intercept_2015 = unname(coef(fit)[["(Intercept)"]]),
      n_years = nrow(d),
      period = "2015-2025"
    )
  }) %>%
  ungroup() %>%
  left_join(bbs, by = "Species")

# For a meaningful slope correlation, express the BBS trend on each standardized
# series' raw-fraction scale. The BBS file gives an annual percent trend (APC),
# not a raw fraction/year slope. First create the BBS-implied path anchored at
# that species-metric's actual 2015 iNaturalist fraction, then fit a raw linear
# slope to that path. This makes BBS and iNaturalist slopes have matching units.
first_rates <- standardized %>%
  filter(year == start_year) %>%
  dplyr::select(Species, standardization, first_iNat_rate_2015 = standardized_rate)

bbs_raw_slopes <- trend_fits %>%
  dplyr::select(Species, standardization, BBS_APC_percent) %>%
  left_join(first_rates, by = c("Species", "standardization")) %>%
  group_by(Species, standardization) %>%
  group_modify(~ {
    if (is.na(.x$first_iNat_rate_2015[1]) || .x$first_iNat_rate_2015[1] <= 0) {
      stop(paste0("Cannot anchor BBS trajectory for ", .y$Species,
                  " / ", .y$standardization, ": 2015 iNaturalist fraction is not positive."))
    }
    bbs_path <- tibble(
      year = years_expected,
      bbs_implied_rate = .x$first_iNat_rate_2015[1] *
        (1 + .x$BBS_APC_percent[1] / 100)^(years_expected - start_year)
    )
    bbs_fit <- lm(bbs_implied_rate ~ I(year - start_year), data = bbs_path)
    tibble(
      BBS_raw_slope_per_year = unname(coef(bbs_fit)[["I(year - start_year)"]]),
      BBS_implied_rate_2025 = tail(bbs_path$bbs_implied_rate, 1)
    )
  }) %>%
  ungroup()

trend_fits <- trend_fits %>%
  left_join(bbs_raw_slopes, by = c("Species", "standardization")) %>%
  mutate(metric_note = if_else(
    standardization == "superuser_observation_proxy",
    "Proxy: denominator is observation count from 10+ bird-observation users, not user count.",
    ""
  ))
write_csv(trend_fits, file.path(output_dir, "inat_and_bbs_raw_slope_estimates.csv"))

# Correlate raw-scale annual slopes across species, separately by standardization.
# Both slope columns are in standardized-fraction units per year.
correlation_summary <- trend_fits %>%
  group_by(standardization) %>%
  group_modify(~ tibble(
    n_species = sum(complete.cases(.x$iNat_raw_slope_per_year,
                                   .x$BBS_raw_slope_per_year)),
    pearson_r = cor(.x$iNat_raw_slope_per_year, .x$BBS_raw_slope_per_year,
                    use = "complete.obs", method = "pearson"),
    spearman_rho = cor(.x$iNat_raw_slope_per_year, .x$BBS_raw_slope_per_year,
                       use = "complete.obs", method = "spearman")
  )) %>%
  ungroup() %>%
  arrange(desc(pearson_r))
write_csv(correlation_summary, file.path(output_dir, "raw_slope_correlation_summary.csv"))

# Prepare plot values. Plot raw annual fractions (not indices), the raw-scale
# iNaturalist linear fit, and the BBS APC trajectory anchored to the first
# iNaturalist fraction in 2015. The BBS curve is implied by its period APC;
# annual observed BBS counts/indices were not included in the supplied CSV.
plot_data <- standardized %>%
  left_join(trend_fits %>%
              dplyr::select(Species, standardization, iNat_raw_slope_per_year,
                     iNat_raw_intercept_2015, BBS_APC_percent),
            by = c("Species", "standardization")) %>%
  left_join(first_rates, by = c("Species", "standardization")) %>%
  mutate(
    iNat_fitted_rate = iNat_raw_intercept_2015 +
      iNat_raw_slope_per_year * (year - start_year),
    BBS_implied_rate = first_iNat_rate_2015 *
      (1 + BBS_APC_percent / 100)^(year - start_year)
  )

make_trend_plot <- function(d, metric_name, facet_species = TRUE) {
  p <- ggplot(d, aes(x = year)) +
    geom_line(aes(y = standardized_rate, color = "iNaturalist annual standardized fractions"),
              linewidth = 0.55, alpha = 0.8) +
    geom_point(aes(y = standardized_rate, color = "iNaturalist annual standardized fractions"),
               size = 1.2, alpha = 0.85) +
    geom_line(aes(y = iNat_fitted_rate, color = "iNaturalist raw-scale linear fit"),
              linewidth = 0.9) +
    geom_line(aes(y = BBS_implied_rate, color = "BBS APC trajectory, anchored at 2015 iNaturalist value"),
              linewidth = 0.9, linetype = "dashed") +
    scale_color_manual(values = c(
      "iNaturalist annual standardized fractions" = "#8BAE9D",
      "iNaturalist raw-scale linear fit" = "#147D64",
      "BBS APC trajectory, anchored at 2015 iNaturalist value" = "#B64926"
    )) +
    scale_x_continuous(breaks = years_expected, minor_breaks = NULL) +
    labs(
      title = paste("Raw-scale iNaturalist vs. BBS trends:", gsub("_", " ", metric_name)),
      subtitle = "iNaturalist slope is fitted to annual fractions; BBS APC path starts at the 2015 iNaturalist fraction",
      x = "Year", y = "Standardized fraction (species count / same-year metric)", color = NULL,
      caption = "Green: linear regression on raw fractions. Dashed orange: BBS annual percent trend applied to the iNaturalist 2015 value. BBS annual survey values were not supplied."
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  if (facet_species) {
    p <- p + facet_wrap(~ Species, scales = "free_y")
  } else {
    p <- p + labs(subtitle = paste(unique(d$Species),
      "| raw fractions; BBS trajectory anchored to 2015 iNaturalist value"))
  }
  p
}

for (metric_name in unique(plot_data$standardization)) {
  metric_data <- plot_data %>% filter(standardization == metric_name)
  safe_metric <- gsub("[^A-Za-z0-9_-]", "_", metric_name)

  overview <- make_trend_plot(metric_data, metric_name, facet_species = TRUE)
  ggsave(file.path(output_dir, paste0("raw_slope_overlay_", safe_metric, ".png")),
         overview, width = 13, height = 8, dpi = 300)

  for (species_name in unique(metric_data$Species)) {
    species_data <- metric_data %>% filter(Species == species_name)
    species_plot <- make_trend_plot(species_data, metric_name, facet_species = FALSE)
    safe_species <- gsub("[^A-Za-z0-9_-]", "_", species_name)
    ggsave(file.path(output_dir,
                     paste0("raw_slope_overlay_", safe_metric, "_", safe_species, ".png")),
           species_plot, width = 8, height = 5, dpi = 300)
  }
}

cat("\nCompleted. Outputs saved in: ", normalizePath(output_dir), "\n", sep = "")
cat("iNaturalist slopes are fitted to raw standardized fractions, not log-transformed values.\n")
cat("BBS raw-scale slopes are derived from the BBS APC trajectory anchored at each 2015 iNaturalist fraction.\n")
cat("\nRaw-slope correlation summary across species:\n")
print(correlation_summary)


#making a big panel with all species and standardization:

# Each panel is one species x standardization combination.
# plot_data must already contain:
# standardized_rate, iNat_fitted_rate, BBS_implied_rate,
# Species, standardization, and year.
library(ggh4x)

all_metrics <- names(metric_denominators)
all_species <- bird_species

plot_data_grid <- plot_data %>%
  mutate(
    Species = factor(Species, levels = all_species),
    standardization = factor(standardization, levels = all_metrics)
  ) %>%
  arrange(Species, standardization, year)

combined_grid <- ggplot(plot_data_grid, aes(x = year)) +
  geom_line(
    aes(y = standardized_rate,
        color = "iNaturalist annual standardized fractions"),
    linewidth = 0.4, alpha = 0.75
  ) +
  geom_point(
    aes(y = standardized_rate,
        color = "iNaturalist annual standardized fractions"),
    size = 0.7, alpha = 0.8
  ) +
  geom_line(
    aes(y = iNat_fitted_rate,
        color = "iNaturalist raw-scale linear fit"),
    linewidth = 0.75
  ) +
  geom_line(
    aes(y = BBS_implied_rate,
        color = "BBS trend anchored at 2015 iNaturalist fraction"),
    linewidth = 0.75, linetype = "dashed"
  ) +
  ggh4x::facet_grid2(
    rows = vars(Species),
    cols = vars(standardization),
    scales = "free_y",
    independent = "y"
  ) +
  scale_color_manual(values = c(
    "iNaturalist annual standardized fractions" = "#8BAE9D",
    "iNaturalist raw-scale linear fit" = "#147D64",
    "BBS trend anchored at 2015 iNaturalist fraction" = "#B64926"
  )) +
  scale_x_continuous(
    breaks = c(2015, 2020, 2025),
    minor_breaks = NULL
  ) +
  labs(
    title = "iNaturalist standardized fractions and BBS trend trajectories",
    subtitle = paste(
      "Rows are species; columns are standardizations.",
      "Each count is divided by that metric's value in the same year."
    ),
    x = "Year",
    y = "Standardized fraction (species count / same-year metric)",
    color = NULL,
    caption = paste(
      "Green: raw-scale linear fit to iNaturalist fractions.",
      "Dashed orange: BBS 2015-2025 trend anchored at the 2015 iNaturalist fraction.",
      "The BBS file provides a period trend, not annual observed BBS values."
    )
  ) +
  theme_minimal(base_size = 9) +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(size = 8, face = "bold"),
    strip.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.spacing = grid::unit(0.35, "lines")
  )

ggsave(
  file.path(output_dir, "raw_slope_all_species_by_standardization.png"),
  plot = combined_grid,
  width = 18, height = 16, dpi = 300,
  limitsize = FALSE
)

# PDF is convenient for zooming into the full panel grid.
ggsave(
  file.path(output_dir, "raw_slope_all_species_by_standardization.pdf"),
  plot = combined_grid,
  width = 18, height = 16,
  limitsize = FALSE
)


#now adding a column for raw Inat trends

# Create an additional, unstandardized raw-count series for each species.
# Fit its raw-scale linear trend across 2015-2025, and create a BBS APC-implied
# path anchored at that species' first iNaturalist raw count in 2015.
raw_count_plot_data <- annual_counts %>%
  select(year, Species, inat_count) %>%
  group_by(Species) %>%
  group_modify(~ {
    d <- .x %>% arrange(year)
    if (nrow(d) != length(years_expected) ||
        !identical(as.integer(d$year), as.integer(years_expected)) ||
        any(is.na(d$inat_count))) {
      stop(paste0("Raw-count panel requires 11 non-missing years for ", .y$Species))
    }
    raw_fit <- lm(inat_count ~ I(year - start_year), data = d)
    d %>% mutate(
      standardization = "raw_iNat_counts",
      standardized_rate = inat_count,
      iNat_fitted_rate = unname(coef(raw_fit)[["(Intercept)"]]) +
        unname(coef(raw_fit)[["I(year - start_year)"]]) * (year - start_year)
    )
  }) %>%
  ungroup() %>%
  left_join(bbs, by = "Species") %>%
  group_by(Species) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    first_iNat_rate_2015 = standardized_rate[year == start_year][1],
    BBS_implied_rate = first_iNat_rate_2015 *
      (1 + BBS_APC_percent / 100)^(year - start_year)
  ) %>%
  ungroup()

# One combined grid: species are rows; the first column is raw iNaturalist
# counts, followed by the standardized metrics. Each cell is one species x
# measure combination. Every panel has an independent y-axis because raw counts
# and standardized fractions have different magnitudes. Standardized panels
# use the year-specific denominators; the raw-count panel does not divide.
all_metrics <- c("raw_iNat_counts", names(metric_denominators))
all_species <- bird_species
plot_data_grid <- bind_rows(plot_data, raw_count_plot_data) %>%
  mutate(
    Species = factor(Species, levels = all_species),
    standardization = factor(standardization, levels = all_metrics)
  ) %>%
  arrange(Species, standardization, year)

combined_grid <- ggplot(plot_data_grid, aes(x = year)) +
  geom_line(aes(y = standardized_rate,
                color = "iNaturalist annual values (raw counts or standardized fractions)"),
            linewidth = 0.4, alpha = 0.75) +
  geom_point(aes(y = standardized_rate,
                 color = "iNaturalist annual values (raw counts or standardized fractions)"),
             size = 0.7, alpha = 0.8) +
  geom_line(aes(y = iNat_fitted_rate,
                color = "iNaturalist raw-scale linear fit"),
            linewidth = 0.75) +
  geom_line(aes(y = BBS_implied_rate,
                color = "BBS APC trajectory anchored at 2015 iNaturalist value"),
            linewidth = 0.75, linetype = "dashed") +
  ggh4x::facet_grid2(
    rows = vars(Species),
    cols = vars(standardization),
    scales = "free_y",
    independent = "y",
    labeller = labeller(standardization = as_labeller(c(
      raw_iNat_counts = "Raw iNaturalist counts",
      total_iNat_observations = "Per total iNat observations",
      total_iNat_users = "Per total iNat users",
      total_bird_observations = "Per bird observations",
      annual_bird_observers = "Per annual bird observers",
      superuser_observation_proxy = "Per superuser-observation proxy"
    )))
  ) +
  scale_color_manual(values = c(
    "iNaturalist annual values (raw counts or standardized fractions)" = "#8BAE9D",
    "iNaturalist raw-scale linear fit" = "#147D64",
    "BBS APC trajectory anchored at 2015 iNaturalist value" = "#B64926"
  )) +
  scale_x_continuous(breaks = c(2015, 2020, 2025), minor_breaks = NULL) +
  labs(
    title = "iNaturalist annual values and BBS trend trajectories",
    subtitle = paste(
      "Rows are species; first column is raw counts, remaining columns are standardizations.",
      "Each standardized annual count is divided by that metric's value in the same year."
    ),
    x = "Year",
    y = "Raw count (first column) or standardized fraction (other columns)",
    color = NULL,
    caption = paste(
      "Solid green: raw-scale linear fit to annual iNaturalist values.",
      "Dashed orange: supplied BBS 2015-2025 APC trajectory anchored to the first 2015 iNaturalist value in each panel.",
      "BBS annual observations were not supplied."
    )
  ) +
  theme_minimal(base_size = 9) +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(size = 8, face = "bold"),
    strip.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.spacing = grid::unit(0.35, "lines")
  )

ggsave(
  filename = file.path(output_dir, "raw_slope_species_by_measure_including_raw_counts.png"),
  
  plot = combined_grid,
  width = 18,
  height = 16,
  dpi = 300,
  limitsize = FALSE
)

# Also save a PDF version, which is easier to zoom when viewing the full grid.
ggsave(
  filename = file.path(output_dir, "raw_slope_species_by_measure_including_raw_counts.pdf"),
  plot = combined_grid,
  width = 18,
  height = 16,
  limitsize = FALSE
)

cat("\nCompleted. Outputs saved in: ", normalizePath(output_dir), "\n", sep = "")
cat("iNaturalist slopes are fitted to raw standardized fractions, not log-transformed values.\n")
cat("BBS raw-scale slopes are derived from the BBS APC trajectory anchored at each 2015 iNaturalist fraction.\n")
cat("\nRaw-slope correlation summary across species:\n")
print(correlation_summary)


#ok now sorting by bbs decling to not declining? and diff colors for each standardization


# Create an additional, unstandardized raw-count series for each species.
# Fit its raw-scale linear trend across 2015-2025, and create a BBS APC-implied
# path anchored at that species' first iNaturalist raw count in 2015.
raw_count_plot_data <- annual_counts %>%
  select(year, Species, inat_count) %>%
  group_by(Species) %>%
  group_modify(~ {
    d <- .x %>% arrange(year)
    if (nrow(d) != length(years_expected) ||
        !identical(as.integer(d$year), as.integer(years_expected)) ||
        any(is.na(d$inat_count))) {
      stop(paste0("Raw-count panel requires 11 non-missing years for ", .y$Species))
    }
    raw_fit <- lm(inat_count ~ I(year - start_year), data = d)
    d %>% mutate(
      standardization = "raw_iNat_counts",
      standardized_rate = inat_count,
      iNat_fitted_rate = unname(coef(raw_fit)[["(Intercept)"]]) +
        unname(coef(raw_fit)[["I(year - start_year)"]]) * (year - start_year)
    )
  }) %>%
  ungroup() %>%
  left_join(bbs, by = "Species") %>%
  group_by(Species) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    first_iNat_rate_2015 = standardized_rate[year == start_year][1],
    BBS_implied_rate = first_iNat_rate_2015 *
      (1 + BBS_APC_percent / 100)^(year - start_year)
  ) %>%
  ungroup()

# One combined grid: species are rows; the first column is raw iNaturalist
# counts, followed by the standardized metrics. Species are ordered by the BBS
# trend from most declining (most negative APC) to least declining / increasing.
# Every panel has an independent y-axis. Each metric has its own color for both
# annual iNaturalist values and the fitted iNaturalist line; BBS is black dashed.
all_metrics <- c("raw_iNat_counts", names(metric_denominators))
species_order <- bbs %>%
  arrange(BBS_APC_percent) %>%
  pull(Species)

plot_data_grid <- bind_rows(plot_data, raw_count_plot_data) %>%
  mutate(
    Species = factor(Species, levels = species_order),
    standardization = factor(standardization, levels = all_metrics)
  ) %>%
  arrange(Species, standardization, year)

metric_colors <- c(
  raw_iNat_counts = "#555555",
  total_iNat_observations = "#1B9E77",
  total_iNat_users = "#D95F02",
  total_bird_observations = "#7570B3",
  annual_bird_observers = "#E7298A",
  superuser_observation_proxy = "#66A61E",
  `BBS trend (anchored at 2015 iNaturalist value)` = "#222222"
)

combined_grid <- ggplot(plot_data_grid, aes(x = year)) +
  geom_line(aes(y = standardized_rate, color = standardization,
                group = interaction(Species, standardization)),
            linewidth = 0.4, alpha = 0.75) +
  geom_point(aes(y = standardized_rate, color = standardization),
             size = 0.7, alpha = 0.8) +
  geom_line(aes(y = iNat_fitted_rate, color = standardization,
                group = interaction(Species, standardization)),
            linewidth = 0.75) +
  geom_line(aes(y = BBS_implied_rate,
                color = "BBS trend (anchored at 2015 iNaturalist value)"),
            linewidth = 0.75, linetype = "dashed", show.legend = TRUE) +
  ggh4x::facet_grid2(
    rows = vars(Species),
    cols = vars(standardization),
    scales = "free_y",
    independent = "y",
    labeller = labeller(standardization = as_labeller(c(
      raw_iNat_counts = "Raw iNaturalist counts",
      total_iNat_observations = "Per total iNat observations",
      total_iNat_users = "Per total iNat users",
      total_bird_observations = "Per bird observations",
      annual_bird_observers = "Per annual bird observers",
      superuser_observation_proxy = "Per superuser-observation proxy"
    )))
  ) +
  scale_color_manual(
    values = metric_colors,
    breaks = c(all_metrics, "BBS trend (anchored at 2015 iNaturalist value)"),
    labels = c(
      "Raw iNaturalist counts",
      "Per total iNat observations",
      "Per total iNat users",
      "Per bird observations",
      "Per annual bird observers",
      "Per superuser-observation proxy",
      "BBS trend (anchored at 2015 iNaturalist value)"
    )
  ) +
  scale_x_continuous(breaks = c(2015, 2020, 2025), minor_breaks = NULL) +
  labs(
    title = "iNaturalist annual values and BBS trend trajectories",
    subtitle = paste(
      "Species rows ordered from most declining to least declining by BBS APC.",
      "First column is raw counts; remaining columns are standardizations using the same-year denominator."
    ),
    x = "Year",
    y = "Raw count (first column) or standardized fraction (other columns)",
    color = NULL,
    caption = paste(
      "Each iNaturalist metric has its own color; raw counts are gray.",
      "Black dashed: supplied BBS 2015-2025 APC trajectory anchored to the 2015 iNaturalist value in each panel.",
      "BBS annual observations were not supplied."
    )
  ) +
  theme_minimal(base_size = 9) +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(size = 8, face = "bold"),
    strip.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.spacing = grid::unit(0.35, "lines")
  )

ggsave(
  filename = file.path(output_dir, "raw_slope_species_by_measure_including_raw_counts.png"),
  
  plot = combined_grid,
  width = 18,
  height = 16,
  dpi = 300,
  limitsize = FALSE
)

# Also save a PDF version, which is easier to zoom when viewing the full grid.
ggsave(
  filename = file.path(output_dir, "raw_slope_species_by_measure_including_raw_counts.pdf"),
  plot = combined_grid,
  width = 18,
  height = 16,
  limitsize = FALSE
)

cat("\nCompleted. Outputs saved in: ", normalizePath(output_dir), "\n", sep = "")
cat("iNaturalist slopes are fitted to raw standardized fractions, not log-transformed values.\n")
cat("BBS raw-scale slopes are derived from the BBS APC trajectory anchored at each 2015 iNaturalist fraction.\n")
cat("\nRaw-slope correlation summary across species:\n")
print(correlation_summary)






