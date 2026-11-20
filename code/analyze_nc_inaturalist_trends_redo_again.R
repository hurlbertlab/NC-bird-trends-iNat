# North Carolina iNaturalist bird trends vs. BBS, 2015-2025
#
# Run from the folder containing the three CSV files:
#   source("analyze_nc_inaturalist_trends.R")
# Requires tidyverse: install.packages("tidyverse")
# Writes tables and plots to ./nc_inat_trend_outputs/.

library(tidyverse)

setwd("C:/Z_School/school/HurlbertLab/lab meeting")


start_year <- 2015L
end_year <- 2025L
years_expected <- start_year:end_year
output_dir <- "nc_inat_trend_outputs_nov13"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


birds <- read_csv("north_carolina_bird_observations_by_year.csv", show_col_types = FALSE)
stats <- read_csv("north_carolina_inaturalist_stats_by_year.csv", show_col_types = FALSE)
bbs <- read_csv("north_carolina_bbs_2015_2025_trend.csv", show_col_types = FALSE) %>%
  transmute(Species, BBS_APC_percent = as.numeric(Trend))

bird_species <- setdiff(names(birds), "year")

# Each denominator is a yearly total shared across all species.
# The final item is only a proxy: the source file has a count of observations
# made by 10+ bird-observation users, not the number of those users.
metric_denominators <- c(
  total_iNat_observations = "total_observations",
  total_iNat_users = "total_inat_users",
  total_bird_observations = "bird_observations",
  annual_bird_observers = "bird_observers",
  superuser_observation_proxy =
    "observations_by_users_with_10_or_more_bird_observations"
)

required_stats <- unique(unname(metric_denominators))
missing_stats <- setdiff(required_stats, names(stats))
if (length(missing_stats) > 0) {
  stop("Missing columns in stats CSV: ", paste(missing_stats, collapse = ", "))
}

# Validate the requested comparison period and species match.
if (!all(years_expected %in% birds$year) || !all(years_expected %in% stats$year)) {
  stop("Both input files must include every year from 2015 through 2025.")
}
missing_bbs <- setdiff(bird_species, bbs$Species)
if (length(missing_bbs) > 0) {
  stop("No BBS trend found for: ", paste(missing_bbs, collapse = ", "))
}

# Long-format annual counts joined to annual denominators; restrict explicitly
# to the 11 years 2015-2025 (inclusive).
annual_counts <- birds %>%
  filter(year %in% years_expected) %>%
  pivot_longer(
    cols = all_of(bird_species), names_to = "Species", values_to = "inat_count"
  ) %>%
  inner_join(stats %>% filter(year %in% years_expected), by = "year")

# One rate column for each available standardization.
standardized <- imap_dfr(metric_denominators, function(denom_col, metric_name) {
  annual_counts %>%
    transmute(
      year, Species, inat_count,
      standardization = metric_name,
      denominator = denom_col,
      denominator_value = .data[[denom_col]],
      standardized_rate = if_else(
        is.na(.data[[denom_col]]) | .data[[denom_col]] <= 0,
        NA_real_, inat_count / .data[[denom_col]]
      )
    )
})

# Save standardized annual values as one row per species-year, one column per metric.
standardized_wide <- standardized %>%
  dplyr::select(year, Species, standardization, standardized_rate) %>%
  pivot_wider(names_from = standardization, values_from = standardized_rate) %>%
  arrange(Species, year)
write_csv(standardized_wide,
          file.path(output_dir, "standardized_species_counts_2015_2025.csv"))

# Fit EXACTLY ONE log-linear trend per species x standardization across the
# full 2015-2025 period: log(rate) = intercept + slope * year.
# This is not a set of year-to-year slopes. APC = 100*(exp(slope)-1).
# Require all 11 positive annual rates so every estimate uses the same period.
trend_fits <- standardized %>%
  group_by(Species, standardization) %>%
  group_modify(~ {
    d <- .x %>% filter(year %in% years_expected) %>% arrange(year)
    if (nrow(d) != length(years_expected) ||
        !identical(as.integer(d$year), as.integer(years_expected)) ||
        any(is.na(d$standardized_rate)) || any(d$standardized_rate <= 0)) {
      stop(paste0(
        "Cannot estimate a full-period log-linear trend for ", .y$Species,
        " / ", .y$standardization,
        ": expected 11 positive, non-missing annual rates for 2015-2025."
      ))
    }
    fit <- lm((standardized_rate) ~ year, data = d)
    slope <- unname(coef(fit)[["year"]])
    tibble(
      iNat_log_slope = slope,
      iNat_APC_percent = 100*(slope),
      n_years = nrow(d),
      period = "2015-2025"
    )
  }) %>%
  ungroup() %>%
  left_join(bbs, by = "Species") %>%
  mutate(
    metric_note = if_else(
      standardization == "superuser_observation_proxy",
      "Proxy only: denominator is observations by 10+ bird-observation users, not a user count.",
      ""
    )
  )
write_csv(trend_fits, file.path(output_dir, "inat_and_bbs_trend_estimates.csv"))

# Correlation across species between the single iNaturalist APC estimate for
# 2015-2025 and the BBS APC for the same period. This is not a correlation of
# annual time-series values. Observer count distinctions unavailable in these
# files: bird_observers supports annual bird observers, but not lifetime-ever
# bird contributors as a separate measure.
correlation_summary <- trend_fits %>%
  group_by(standardization) %>%
  group_modify(~ {
    tibble(
      n_species = sum(complete.cases(.x$iNat_APC_percent, .x$BBS_APC_percent)),
      pearson_r = cor(.x$iNat_APC_percent, .x$BBS_APC_percent,
                      use = "complete.obs", method = "pearson"),
      spearman_rho = cor(.x$iNat_APC_percent, .x$BBS_APC_percent,
                         use = "complete.obs", method = "spearman")
    )
  }) %>%
  ungroup() %>%
  arrange(desc(pearson_r))
write_csv(correlation_summary, file.path(output_dir, "correlation_summary.csv"))

correlation_summary2 <- trend_fits %>%
  group_by(standardization) %>%
  summarize(cor = cor(iNat_APC_percent, BBS_APC_percent, method ="pearson"))


par(mfrow= c(2,3))
for (s in unique(trend_fits$standardization)) {
  tmp = filter(trend_fits, standardization == s)
  plot(tmp$iNat_APC_percent, tmp$BBS_APC_percent, main = s)
}






# Plot annual standardized values and the fitted trend for each metric.
# The annual rate is calculated explicitly for every species-year as:
#   iNat species observations in that year / that year's dplyr::selected denominator.
# Each iNaturalist series and the BBS-implied trend is indexed to 100 in 2015.
# One overview graph is saved per standardization (species are separate facets
# with free y-scales). Individual species PNGs are also saved for easier review.
# BBS supplies a single overall 2015-2025 APC, not annual survey values, so its
# plotted line is the trajectory implied by that APC, not observed BBS data.
annual_plot_data <- standardized %>%
  left_join(trend_fits %>%
              dplyr::select(Species, standardization, iNat_log_slope,
                     BBS_APC_percent),
            by = c("Species", "standardization")) %>%
  group_by(Species, standardization) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    baseline_rate_2015 = standardized_rate[year == start_year][1],
    iNat_observed_index = 100 * standardized_rate / baseline_rate_2015,
    iNat_fitted_index = 100 * exp(iNat_log_slope * (year - start_year)),
    BBS_trend_index = 100 * (1 + BBS_APC_percent / 100)^(year - start_year)
  ) %>%
  ungroup()

if (any(is.na(annual_plot_data$baseline_rate_2015) |
        annual_plot_data$baseline_rate_2015 <= 0)) {
  stop("Cannot index a standardized species series to 2015: its 2015 rate is missing or non-positive.")
}

make_trend_plot <- function(d, metric_name, facet_species = TRUE) {
  p <- ggplot(d, aes(x = year)) +
    geom_line(aes(y = iNat_observed_index,
                  color = "iNaturalist annual standardized values"),
              linewidth = 0.45, alpha = 0.7) +
    geom_point(aes(y = iNat_observed_index,
                   color = "iNaturalist annual standardized values"),
               size = 1.1, alpha = 0.8) +
    geom_line(aes(y = iNat_fitted_index,
                  color = "iNaturalist fitted trend (2015-2025)"),
              linewidth = 0.9) +
    geom_line(aes(y = BBS_trend_index,
                  color = "BBS trend-implied trajectory (2015-2025)"),
              linewidth = 0.9, linetype = "dashed") +
    scale_color_manual(values = c(
      "iNaturalist annual standardized values" = "#8BAE9D",
      "iNaturalist fitted trend (2015-2025)" = "#147D64",
      "BBS trend-implied trajectory (2015-2025)" = "#B64926"
    )) +
    scale_x_continuous(breaks = years_expected, minor_breaks = NULL) +
    labs(
      title = paste("iNaturalist vs. BBS:", gsub("_", " ", metric_name)),
      subtitle = "Annual iNaturalist counts are divided by that same year's dplyr::selected metric; all indices start at 100 in 2015",
      x = "Year", y = "Index (2015 = 100)", color = NULL,
      caption = paste(
        "Green line is one log-linear iNaturalist trend fitted across 2015-2025.",
        "Dashed orange line is implied by the supplied BBS period trend."
      )
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  if (facet_species) {
    p <- p + facet_wrap(~ Species, scales = "free_y")
  } else {
    p <- p + labs(subtitle = paste(unique(d$Species),
      "| annual iNaturalist values use year-specific denominator; index = 100 in 2015"))
  }
  p
}

for (metric_name in unique(annual_plot_data$standardization)) {
  metric_data <- annual_plot_data %>% filter(standardization == metric_name)
  safe_metric <- gsub("[^A-Za-z0-9_-]", "_", metric_name)

  # One overview PNG per standardization; each species gets its own free-scale facet.
  overview <- make_trend_plot(metric_data, metric_name, facet_species = TRUE)
  ggsave(file.path(output_dir, paste0("trend_overlay_", safe_metric, ".png")),
         overview, width = 13, height = 8, dpi = 300)

  # Also save one separate graph per species for visual comparison.
  for (species_name in unique(metric_data$Species)) {
    species_data <- metric_data %>% filter(Species == species_name)
    species_plot <- make_trend_plot(species_data, metric_name, facet_species = FALSE)
    safe_species <- gsub("[^A-Za-z0-9_-]", "_", species_name)
    ggsave(file.path(output_dir,
                     paste0("trend_overlay_", safe_metric, "_", safe_species, ".png")),
           species_plot, width = 8, height = 5, dpi = 300)
  }
}

cat("\nCompleted. Outputs saved in: ", normalizePath(output_dir), "\n", sep = "")
cat("\nEach iNaturalist APC is one log-linear slope fit over all 11 years, 2015-2025.\n")
cat("BBS trend overlay is the trajectory implied by its single period trend; annual BBS observations were not supplied.\n")
cat("\nCorrelation summary across species:\n")
print(correlation_summary)
