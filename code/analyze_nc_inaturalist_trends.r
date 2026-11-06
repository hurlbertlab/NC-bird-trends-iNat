# Standardize North Carolina iNaturalist bird counts, estimate trends, and
# compare them with the supplied BBS annual trends (2015-2025).
#
# Put this script in the same folder as the three CSV files, set that folder as
# the R working directory, and run: source("analyze_nc_inaturalist_trends.R")
# Requires: tidyverse (install.packages("tidyverse") if needed).
# Outputs are written to ./nc_inat_trend_outputs/.

library(tidyverse)

start_year <- 2015L
end_year <- 2025L
output_dir <- "nc_inat_trend_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

birds_file <- "north_carolina_bird_observations_by_year.csv"
stats_file <- "north_carolina_inaturalist_stats_by_year.csv"
bbs_file <- "north_carolina_bbs_2015_2025_trend.csv"

birds <- read_csv(birds_file, show_col_types = FALSE)
stats <- read_csv(stats_file, show_col_types = FALSE)
bbs <- read_csv(bbs_file, show_col_types = FALSE) %>%
  mutate(Trend = as.numeric(Trend))

# The stats file has bird_observers (annual bird-observation contributors),
# but does not distinguish lifetime-ever bird observers from users who
# contributed a bird observation in the current year. These two requested
# descriptions therefore cannot be independently estimated from this file.
# The alias below is retained to make the duplicate definition explicit.
#
# The stats file also has observations_by_users_with_10_or_more_bird_observations,
# which is an OBSERVATION count, not a count of superusers. Its use below is
# only a clearly labeled proxy; the exact superuser-user denominator is absent.

bird_species <- setdiff(names(birds), "year")
needed_stats <- c(
  "total_observations", "total_inat_users", "bird_observations",
  "bird_observers",
  "observations_by_users_with_10_or_more_bird_observations"
)
missing_stats <- setdiff(needed_stats, names(stats))
if (length(missing_stats) > 0) {
  stop("Missing required columns in stats CSV: ", paste(missing_stats, collapse = ", "))
}

# One denominator per available standardization. The two observer metrics
# requested in the prompt map to the same annual bird_observers field.
metric_denominators <- c(
  total_iNat_observations = "total_observations",
  total_iNat_users = "total_inat_users",
  total_bird_observations = "bird_observations",
  annual_bird_observers = "bird_observers",
  bird_contributors_same_year = "bird_observers", # duplicate alias, not independent
  superuser_observation_proxy =
    "observations_by_users_with_10_or_more_bird_observations"
)

# This is intentionally NA: no count of superuser USERS is present in the input.
missing_superuser_metric <- "superuser_users_unavailable"

birds_long <- birds %>%
  pivot_longer(
    cols = all_of(bird_species),
    names_to = "Species",
    values_to = "inat_bird_count"
  )

analysis <- birds_long %>%
  inner_join(stats, by = "year") %>%
  filter(year >= start_year, year <= end_year) %>%
  mutate(!!missing_superuser_metric := NA_real_)

# Make a long table of rates for every available metric. Do not replace a zero
# denominator with an arbitrary value: those rates are undefined (NA).
standardized <- imap_dfr(metric_denominators, function(denominator, metric) {
  analysis %>%
    transmute(
      year, Species, inat_bird_count,
      standardization = metric,
      denominator_name = denominator,
      denominator_value = .data[[denominator]],
      standardized_rate = if_else(
        is.na(.data[[denominator]]) | .data[[denominator]] <= 0,
        NA_real_,
        inat_bird_count / .data[[denominator]]
      )
    )
})

# Add an explicit unavailable row for the exact number of active superusers.
# This preserves the requested metric in the output without mislabeling the
# available observation-count proxy as a user count.
unavailable_superusers <- analysis %>%
  transmute(
    year, Species, inat_bird_count,
    standardization = missing_superuser_metric,
    denominator_name = "not_in_input_data",
    denominator_value = NA_real_,
    standardized_rate = NA_real_
  )

standardized_all <- bind_rows(standardized, unavailable_superusers)

# Wide annual output: one standardized column per metric/species-year row.
standardized_wide <- standardized_all %>%
  select(year, Species, standardized_rate, standardization) %>%
  pivot_wider(names_from = standardization, values_from = standardized_rate) %>%
  arrange(Species, year)
write_csv(standardized_wide, file.path(output_dir, "standardized_species_counts_2015_2025.csv"))

# Estimate an annual percent change (APC) for each species x standardization
# using a log-linear model: log(rate) = intercept + slope * year.
# APC = 100 * (exp(slope) - 1), directly comparable in form to the BBS Trend.
# All analyzed iNaturalist rates must be positive for this model; no pseudocount
# is added. If a series contains zero rates, the script records the number of
# usable years and fits only positive rates (inspect n_years in the output).
inat_trends <- standardized %>%
  filter(!is.na(standardized_rate), standardized_rate > 0) %>%
  group_by(Species, standardization) %>%
  group_modify(~ {
    fit <- lm(log(standardized_rate) ~ year, data = .x)
    slope <- unname(coef(fit)[["year"]])
    tibble(
      inat_APC_percent = 100 * (exp(slope) - 1),
      n_years = nrow(.x),
      first_year = min(.x$year),
      last_year = max(.x$year)
    )
  }) %>%
  ungroup()

# Include the unavailable exact superuser-user metric as an explicit row.
unavailable_trends <- tidyr::crossing(
  Species = bird_species,
  standardization = missing_superuser_metric
) %>%
  mutate(
    inat_APC_percent = NA_real_, n_years = 0L,
    first_year = NA_integer_, last_year = NA_integer_
  )

trend_comparison <- inat_trends %>%
  full_join(unavailable_trends, by = c("Species", "standardization")) %>%
  left_join(bbs %>% select(Species, BBS_APC_percent = Trend), by = "Species") %>%
  mutate(
    metric_note = case_when(
      standardization == "annual_bird_observers" ~
        "Uses annual bird_observers; lifetime-ever and current-year contributor counts are not separately available.",
      standardization == "bird_contributors_same_year" ~
        "Duplicate alias of annual_bird_observers; same source field.",
      standardization == "superuser_observation_proxy" ~
        "Uses observations by 10+ bird-observation users, not a count of those users.",
      standardization == missing_superuser_metric ~
        "Unavailable: input has no count of active 10+ bird-observation users.",
      TRUE ~ ""
    )
  )
write_csv(trend_comparison, file.path(output_dir, "inat_and_bbs_trend_estimates.csv"))

# Correlate species-level iNaturalist APC estimates with BBS APCs. The
# observer alias is excluded from the ranking because it duplicates a metric;
# the proxy is included but separately labeled. Pearson and Spearman use the
# same complete species set for each available standardization.
correlation_summary <- trend_comparison %>%
  filter(
    !standardization %in% c("bird_contributors_same_year", missing_superuser_metric),
    !is.na(inat_APC_percent), !is.na(BBS_APC_percent)
  ) %>%
  group_by(standardization) %>%
  group_modify(~ {
    pearson <- cor.test(.x$inat_APC_percent, .x$BBS_APC_percent,
                        method = "pearson", exact = FALSE)
    spearman <- cor.test(.x$inat_APC_percent, .x$BBS_APC_percent,
                         method = "spearman", exact = FALSE)
    tibble(
      n_species = nrow(.x),
      pearson_r = unname(pearson$estimate),
      pearson_p = pearson$p.value,
      spearman_rho = unname(spearman$estimate),
      spearman_p = spearman$p.value
    )
  }) %>%
  ungroup() %>%
  arrange(desc(pearson_r))
write_csv(correlation_summary, file.path(output_dir, "correlation_summary.csv"))

# Build annual iNaturalist index (2015 = 100) and the BBS-implied trajectory
# (2015 = 100, using the published annual trend). Facet each figure by species.
# Each available standardization gets an overlay plot; the unavailable exact
# superuser-user metric is skipped because it has no denominator data.
plot_data <- standardized %>%
  filter(!is.na(standardized_rate), standardized_rate > 0) %>%
  group_by(Species, standardization) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    inat_index_2015_100 = 100 * standardized_rate / standardized_rate[year == start_year][1]
  ) %>%
  ungroup() %>%
  left_join(bbs %>% select(Species, BBS_APC_percent = Trend), by = "Species") %>%
  mutate(bbs_index_2015_100 = 100 * (1 + BBS_APC_percent / 100)^(year - start_year))

for (metric in unique(plot_data$standardization)) {
  metric_data <- plot_data %>% filter(standardization == metric)
  p <- ggplot(metric_data, aes(x = year)) +
    geom_line(aes(y = inat_index_2015_100, color = "iNaturalist standardized index"),
              linewidth = 0.8) +
    geom_point(aes(y = inat_index_2015_100, color = "iNaturalist standardized index"),
               size = 1.3) +
    geom_line(aes(y = bbs_index_2015_100, color = "BBS trend-implied index"),
              linewidth = 0.8, linetype = "dashed") +
    facet_wrap(~ Species, scales = "free_y") +
    scale_color_manual(values = c(
      "iNaturalist standardized index" = "#147D64",
      "BBS trend-implied index" = "#B64926"
    )) +
    scale_x_continuous(breaks = seq(start_year, end_year, by = 2)) +
    labs(
      title = paste("North Carolina bird trends:", gsub("_", " ", metric)),
      subtitle = "Both series indexed to 100 in 2015; BBS line uses its 2015-2025 annual trend",
      x = "Year", y = "Index (2015 = 100)", color = NULL,
      caption = "iNaturalist series = species count divided by the selected annual denominator."
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    )
  safe_metric <- gsub("[^A-Za-z0-9_-]", "_", metric)
  ggsave(
    filename = file.path(output_dir, paste0("overlay_", safe_metric, ".png")),
    plot = p, width = 13, height = 8, dpi = 300
  )
}

# Report the key comparability check in the R console.
cat("\nOutputs saved to:", normalizePath(output_dir), "\n")
cat("\nCorrelation summary (available distinct denominators plus proxy):\n")
print(correlation_summary)
cat("\nReminder: total_iNat_users and total_iNat_observations are common annual\n",
    "denominators shared by all species. Such common scaling shifts every species'\n",
    "log-linear APC by the same amount, so cross-species correlations with BBS\n",
    "can be identical across standardizations.\n", sep = "")
