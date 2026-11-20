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
output_dir <- "nc_inat_trend_outputs_raw_slope"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

birds <- read_csv("north_carolina_bird_observations_by_year.csv", show_col_types = FALSE)
stats <- read_csv("north_carolina_inaturalist_stats_by_year.csv", show_col_types = FALSE)
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
