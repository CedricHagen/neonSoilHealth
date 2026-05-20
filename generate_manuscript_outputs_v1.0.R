cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(viridis)
  library(scales)
  library(gridExtra)
  library(maps)
  library(broom)
})

setwd("~/Desktop/neonSoilHealth/manuscript_submission/")
dir.create("manuscript_outputs", showWarnings = FALSE)
dir.create("manuscript_outputs/tables", showWarnings = FALSE)
dir.create("manuscript_outputs/figures", showWarnings = FALSE)
dir.create("manuscript_outputs/data", showWarnings = FALSE)

cat("Output directories created.\n")

cat("Loading dataset...\n")
neon_plfa <- read.csv("neon_soil_health_2026-05-18.csv", stringsAsFactors = FALSE)
neon_plfa$collectDate <- as.Date(neon_plfa$collectDate)

cat("Loaded", nrow(neon_plfa), "samples from", length(unique(neon_plfa$siteID)), "sites\n")

cat("\n=== Calculating Gram-positive:Gram-negative ratio ===\n")

neon_plfa$gram_positive <- rowSums(neon_plfa[, c(
  "i14To0ScaledConcentration",
  "i15To0ScaledConcentration",
  "aC15To0ScaledConcentration",
  "i16To0ScaledConcentration",
  "i17To0ScaledConcentration",
  "c17To0AnteisoScaledConcentration"
)], na.rm = TRUE)

neon_plfa$gram_negative <- rowSums(neon_plfa[, c(
  "c16To1n7ScaledConcentration",
  "c18To1n11ScaledConcentration"
)], na.rm = TRUE)

neon_plfa$GP_GN_ratio <- ifelse(neon_plfa$gram_negative > 0,
                                 neon_plfa$gram_positive / neon_plfa$gram_negative,
                                 NA)

cat("GP:GN ratio calculated for", sum(!is.na(neon_plfa$GP_GN_ratio)), "samples\n")

cat("\n=== Generating clean dataset (without chemistry) ===\n")

chemistry_cols <- c(
  "organicCPercent", "nitrogenPercent", "CN_ratio",
  "organicd13C", "d15N", "biomass_per_gC_nmol_gC",
  "plotID.chem", "collectDateTime.chem", "collectDate.chem",
  "year.chem", "month.chem"
)

neon_plfa_clean <- neon_plfa %>%
  select(-any_of(chemistry_cols))

write.csv(neon_plfa_clean,
          "manuscript_outputs/data/neon_plfa_synthesis_v1.0.csv",
          row.names = FALSE)

cat("✓ Saved: neon_plfa_synthesis_v1.0.csv\n")

cat("\n=== Generating supplementary tables ===\n")

site_metadata <- neon_plfa_clean %>%
  group_by(siteID) %>%
  summarize(
    domain = first(domainID),
    latitude = round(first(decimalLatitude), 4),
    longitude = round(first(decimalLongitude), 4),
    elevation_m = round(first(elevation), 1),
    soil_order = first(soilOrder),
    soil_suborder = first(soilSuborder),
    n_samples = n(),
    n_years = n_distinct(year),
    years_sampled = paste(sort(unique(year)), collapse=", "),
    first_sample = min(collectDate, na.rm = TRUE),
    last_sample = max(collectDate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(domain, siteID)

write.csv(site_metadata, "manuscript_outputs/tables/Table1_site_metadata.csv", row.names = FALSE)

calc_summary <- function(x) {
  x <- x[!is.na(x)]
  data.frame(
    n = length(x),
    mean = round(mean(x), 2),
    sd = round(sd(x), 2),
    min = round(min(x), 2),
    q25 = round(quantile(x, 0.25), 2),
    median = round(median(x), 2),
    q75 = round(quantile(x, 0.75), 2),
    max = round(max(x), 2)
  )
}

summary_stats <- bind_rows(
  calc_summary(neon_plfa_clean$microbial_biomass_nmol_g) %>%
    mutate(metric = "Microbial biomass (nmol/g)", .before = 1),
  calc_summary(neon_plfa_clean$fungal_plfa_nmol_g) %>%
    mutate(metric = "Fungal PLFA (nmol/g)", .before = 1),
  calc_summary(neon_plfa_clean$bacterial_plfa_nmol_g) %>%
    mutate(metric = "Bacterial PLFA (nmol/g)", .before = 1),
  calc_summary(neon_plfa_clean$FB_ratio) %>%
    mutate(metric = "F:B ratio", .before = 1),
  calc_summary(neon_plfa_clean$stress_index) %>%
    mutate(metric = "Stress index", .before = 1),
  calc_summary(neon_plfa_clean$GP_GN_ratio) %>%
    mutate(metric = "GP:GN ratio", .before = 1)
)

write.csv(summary_stats, "manuscript_outputs/tables/Table2_summary_statistics.csv", row.names = FALSE)

temporal_matrix <- neon_plfa_clean %>%
  group_by(siteID, year) %>%
  summarize(n_samples = n(), .groups = "drop") %>%
  pivot_wider(names_from = year, values_from = n_samples, values_fill = 0)
write.csv(temporal_matrix, "manuscript_outputs/tables/temporal_coverage_matrix.csv", row.names = FALSE)

completeness_by_site <- neon_plfa_clean %>%
  group_by(siteID, domainID) %>%
  summarize(
    n_total = n(),
    pct_biomass = round(100 * sum(!is.na(microbial_biomass_nmol_g)) / n(), 1),
    pct_fungal = round(100 * sum(!is.na(fungal_plfa_nmol_g)) / n(), 1),
    pct_bacterial = round(100 * sum(!is.na(bacterial_plfa_nmol_g)) / n(), 1),
    pct_FB = round(100 * sum(!is.na(FB_ratio)) / n(), 1),
    pct_stress = round(100 * sum(!is.na(stress_index)) / n(), 1),
    pct_GPGN = round(100 * sum(!is.na(GP_GN_ratio)) / n(), 1),
    .groups = "drop"
  )
write.csv(completeness_by_site, "manuscript_outputs/tables/completeness_by_site.csv", row.names = FALSE)

horizon_summary <- neon_plfa_clean %>%
  group_by(horizon) %>%
  summarize(
    n_samples = n(),
    pct_total = round(100 * n() / nrow(neon_plfa_clean), 1),
    n_sites = n_distinct(siteID),
    mean_biomass = round(mean(microbial_biomass_nmol_g, na.rm = TRUE), 2),
    sd_biomass = round(sd(microbial_biomass_nmol_g, na.rm = TRUE), 2),
    mean_FB = round(mean(FB_ratio, na.rm = TRUE), 3),
    sd_FB = round(sd(FB_ratio, na.rm = TRUE), 3),
    mean_stress = round(mean(stress_index, na.rm = TRUE), 2),
    sd_stress = round(sd(stress_index, na.rm = TRUE), 2),
    mean_GPGN = round(mean(GP_GN_ratio, na.rm = TRUE), 2),
    sd_GPGN = round(sd(GP_GN_ratio, na.rm = TRUE), 2),
    .groups = "drop"
  )
write.csv(horizon_summary, "manuscript_outputs/tables/horizon_summary.csv", row.names = FALSE)

domain_summary <- neon_plfa_clean %>%
  group_by(domainID) %>%
  summarize(
    n_samples = n(),
    n_sites = n_distinct(siteID),
    mean_biomass = round(mean(microbial_biomass_nmol_g, na.rm = TRUE), 2),
    sd_biomass = round(sd(microbial_biomass_nmol_g, na.rm = TRUE), 2),
    mean_FB = round(mean(FB_ratio, na.rm = TRUE), 3),
    sd_FB = round(sd(FB_ratio, na.rm = TRUE), 3),
    mean_stress = round(mean(stress_index, na.rm = TRUE), 2),
    sd_stress = round(sd(stress_index, na.rm = TRUE), 2),
    mean_GPGN = round(mean(GP_GN_ratio, na.rm = TRUE), 2),
    sd_GPGN = round(sd(GP_GN_ratio, na.rm = TRUE), 2),
    .groups = "drop"
  )
write.csv(domain_summary, "manuscript_outputs/tables/domain_summary.csv", row.names = FALSE)

cat("✓ All tables generated\n")


cat("\n=== Generating Figures ===\n")
theme_set(theme_bw(base_size = 11))

cat("Generating Figure 1a: Site Map...\n")

north_america <- map_data("world", region = c("USA", "Canada", "Mexico", "Puerto Rico"))

fig1a <- ggplot() +
  geom_polygon(data = north_america,
               aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(data = site_metadata,
             aes(x = longitude, y = latitude, size = n_samples, fill = domain),
             shape = 21, alpha = 0.8, stroke = 0.5) +
  scale_size_continuous(range = c(2, 10), name = "N samples",
                        breaks = c(100, 200, 300, 400, 500)) +
  scale_fill_viridis_d(option = "turbo", name = "Domain") +
  coord_fixed(1.3, xlim = c(-180, -50), ylim = c(14, 72)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.box = "vertical",
        panel.background = element_rect(fill = "aliceblue"),
        plot.margin = margin(10, 40, 10, 10))  # Extra right margin for legend

ggsave("manuscript_outputs/figures/Figure1a_site_map.png", fig1a,
       width = 12, height = 7, dpi = 300)
ggsave("manuscript_outputs/figures/Figure1a_site_map.pdf", fig1a,
       width = 12, height = 7)

cat("✓ Figure 1a saved\n")


cat("Generating Figure 1b: Temporal Coverage...\n")

temporal_long <- temporal_matrix %>%
  pivot_longer(-siteID, names_to = "year", values_to = "n_samples") %>%
  mutate(year = as.numeric(year))

fig1b <- ggplot(temporal_long, aes(x = year, y = siteID, fill = n_samples)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "plasma", name = "N samples",
                       breaks = c(0, 50, 100, 200, 300, 400, 500)) +
  scale_x_continuous(breaks = 2017:2024) +
  labs(x = "Year", y = "Site") +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_text(size = 6),
        panel.grid = element_blank())

ggsave("manuscript_outputs/figures/Figure1b_temporal_heatmap.png", fig1b,
       width = 8, height = 10, dpi = 300)
ggsave("manuscript_outputs/figures/Figure1b_temporal_heatmap.pdf", fig1b,
       width = 8, height = 10)

cat("✓ Figure 1b saved\n")

cat("Generating Figure 1c: Horizon Distribution...\n")

fig1c <- ggplot(horizon_summary, aes(x = "", y = n_samples, fill = horizon)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  scale_fill_manual(values = c("M" = "#8B4513", "O" = "#2E8B57"),
                    labels = c("M" = paste0("Mineral (", horizon_summary$pct_total[horizon_summary$horizon == "M"], "%)"),
                               "O" = paste0("Organic (", horizon_summary$pct_total[horizon_summary$horizon == "O"], "%)"))) +
  labs(fill = "Horizon") +
  theme_void()

ggsave("manuscript_outputs/figures/Figure1c_horizon_distribution.png", fig1c,
       width = 6, height = 5, dpi = 300)

cat("✓ Figure 1c saved\n")

cat("Generating Figure 2: Metric Distributions...\n")

fig2a <- ggplot(neon_plfa_clean, aes(x = log10(microbial_biomass_nmol_g))) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  labs(x = "log10(Biomass) [nmol/g]", y = "Count") +
  theme_minimal()

fig2b <- ggplot(neon_plfa_clean, aes(x = FB_ratio)) +
  geom_histogram(bins = 50, fill = "darkgreen", alpha = 0.7) +
  labs(x = "Fungal:Bacterial Ratio", y = "Count") +
  theme_minimal()

fig2c <- ggplot(neon_plfa_clean %>% filter(stress_index < 5),
                aes(x = stress_index)) +
  geom_histogram(bins = 50, fill = "coral", alpha = 0.7) +
  labs(x = "Stress Index", y = "Count") +
  theme_minimal()

fig2d <- ggplot(neon_plfa_clean %>% filter(!is.na(GP_GN_ratio) & GP_GN_ratio < 10),
                aes(x = GP_GN_ratio)) +
  geom_histogram(bins = 50, fill = "purple", alpha = 0.7) +
  labs(x = "Gram+:Gram- Ratio", y = "Count") +
  theme_minimal()

fig2_combined <- grid.arrange(fig2a, fig2b, fig2c, fig2d, ncol = 2, nrow = 2)

ggsave("manuscript_outputs/figures/Figure2_distributions.png", fig2_combined,
       width = 8, height = 10, dpi = 300)
ggsave("manuscript_outputs/figures/Figure2_distributions.pdf", fig2_combined,
       width = 8, height = 10)

cat("✓ Figure 2 saved\n")

cat("Generating Figure 3: Latitudinal Gradients...\n")

site_averages <- neon_plfa_clean %>%
  group_by(siteID) %>%
  summarize(
    latitude = first(decimalLatitude),
    mean_biomass = mean(microbial_biomass_nmol_g, na.rm = TRUE),
    mean_FB = mean(FB_ratio, na.rm = TRUE),
    mean_stress = mean(stress_index, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

add_lm_stats <- function(data, x_col, y_col) {
  df <- data.frame(x = data[[x_col]], y = data[[y_col]]) %>%
    filter(!is.na(x) & !is.na(y))

  if (nrow(df) < 3) return(NULL)

  fit <- lm(y ~ x, data = df)
  fit_summary <- summary(fit)

  r2 <- round(fit_summary$r.squared, 3)
  pval <- fit_summary$coefficients["x", "Pr(>|t|)"]
  pval_text <- if (pval < 0.001) "p < 0.001" else paste0("p = ", round(pval, 3))
  n <- nrow(df)

  list(
    label = paste0("N = ", n, "\nR² = ", r2, "\n", pval_text),
    r2 = r2,
    pval = pval,
    n = n
  )
}

stats_biomass <- add_lm_stats(site_averages, "latitude", "mean_biomass")

fig3a <- ggplot(site_averages, aes(x = latitude, y = log10(mean_biomass))) +
  geom_point(size = 3, alpha = 0.7, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "darkblue", fill = "lightblue") +
  annotate("text", x = 20, y = 2.9,
           label = stats_biomass$label, hjust = 0, vjust = 1, size = 3.5,
           fontface = "bold") +
  labs(x = "Latitude (°N)", y = "log10(Mean Biomass) [nmol/g]") +
  theme_minimal() +
  theme(plot.margin = margin(10, 10, 10, 10))

stats_FB <- add_lm_stats(site_averages, "latitude", "mean_FB")

fig3b <- ggplot(site_averages, aes(x = latitude, y = mean_FB)) +
  geom_point(size = 3, alpha = 0.7, color = "darkgreen") +
  geom_smooth(method = "lm", se = TRUE, color = "darkgreen", fill = "lightgreen") +
  annotate("text", x = 20, y = 0.29,
           label = stats_FB$label, hjust = 0, vjust = 1, size = 3.5,
           fontface = "bold") +
  labs(x = "Latitude (°N)", y = "Mean F:B Ratio") +
  theme_minimal() +
  theme(plot.margin = margin(10, 10, 10, 10))

stats_stress <- add_lm_stats(site_averages, "latitude", "mean_stress")

fig3c <- ggplot(site_averages, aes(x = latitude, y = mean_stress)) +
  geom_point(size = 3, alpha = 0.7, color = "coral") +
  geom_smooth(method = "lm", se = TRUE, color = "red4", fill = "pink") +
  annotate("text", x = 20, y = 2.6,
           label = stats_stress$label, hjust = 0, vjust = 1, size = 3.5,
           fontface = "bold") +
  labs(x = "Latitude (°N)", y = "Mean Stress Index") +
  theme_minimal() +
  theme(plot.margin = margin(10, 10, 10, 10))

fig3_combined <- grid.arrange(fig3a, fig3b, fig3c, ncol = 1)

ggsave("manuscript_outputs/figures/Figure3_latitudinal_gradients.png", fig3_combined,
       width = 8, height = 11, dpi = 300)
ggsave("manuscript_outputs/figures/Figure3_latitudinal_gradients.pdf", fig3_combined,
       width = 8, height = 11)

lat_stats <- data.frame(
  metric = c("log10(Biomass)", "F:B Ratio", "Stress Index"),
  n_sites = c(stats_biomass$n, stats_FB$n, stats_stress$n),
  r_squared = c(stats_biomass$r2, stats_FB$r2, stats_stress$r2),
  p_value = c(stats_biomass$pval, stats_FB$pval, stats_stress$pval)
)
write.csv(lat_stats, "manuscript_outputs/tables/latitudinal_gradient_stats.csv", row.names = FALSE)

cat("✓ Figure 3 saved\n")


cat("Generating Figure 4: Domain patterns...\n")

fig4a <- ggplot(neon_plfa_clean, aes(x = domainID, y = log10(microbial_biomass_nmol_g))) +
  geom_boxplot(fill = "steelblue", alpha = 0.6, outlier.size = 0.5) +
  labs(x = "NEON Domain", y = "log10(Biomass) [nmol/g]") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig4b <- ggplot(neon_plfa_clean, aes(x = domainID, y = FB_ratio)) +
  geom_boxplot(fill = "darkgreen", alpha = 0.6, outlier.size = 0.5) +
  labs(x = "NEON Domain", y = "F:B Ratio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig4c <- ggplot(neon_plfa_clean %>% filter(stress_index < 5),
                aes(x = domainID, y = stress_index)) +
  geom_boxplot(fill = "coral", alpha = 0.6, outlier.size = 0.5) +
  labs(x = "NEON Domain", y = "Stress Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig4_combined <- grid.arrange(fig4a, fig4b, fig4c, ncol = 1)

ggsave("manuscript_outputs/figures/Figure4_domain_patterns.png", fig4_combined,
       width = 10, height = 12, dpi = 300)
ggsave("manuscript_outputs/figures/Figure4_domain_patterns.pdf", fig4_combined,
       width = 10, height = 12)

cat("✓ Figure 4 saved\n")


cat("Generating Figure 5: Data Quality...\n")

plfa_completeness <- data.frame(
  metric = c("Biomass", "Fungal PLFA", "Bacterial PLFA", "F:B Ratio", "Stress Index", "GP:GN Ratio"),
  pct_complete = c(
    100 * sum(!is.na(neon_plfa_clean$microbial_biomass_nmol_g)) / nrow(neon_plfa_clean),
    100 * sum(!is.na(neon_plfa_clean$fungal_plfa_nmol_g)) / nrow(neon_plfa_clean),
    100 * sum(!is.na(neon_plfa_clean$bacterial_plfa_nmol_g)) / nrow(neon_plfa_clean),
    100 * sum(!is.na(neon_plfa_clean$FB_ratio)) / nrow(neon_plfa_clean),
    100 * sum(!is.na(neon_plfa_clean$stress_index)) / nrow(neon_plfa_clean),
    100 * sum(!is.na(neon_plfa_clean$GP_GN_ratio)) / nrow(neon_plfa_clean)
  )
)

fig5a <- ggplot(plfa_completeness, aes(x = reorder(metric, pct_complete), y = pct_complete)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = paste0(round(pct_complete, 1), "%")),
            hjust = -0.1, size = 4) +
  coord_flip() +
  ylim(0, 105) +
  labs(x = "", y = "% Complete") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank())

ggsave("manuscript_outputs/figures/Figure5a_completeness.png", fig5a,
       width = 8, height = 5, dpi = 300)
ggsave("manuscript_outputs/figures/Figure5a_completeness.pdf", fig5a,
       width = 8, height = 5)

site_year_means <- neon_plfa_clean %>%
  group_by(siteID, year) %>%
  summarize(
    mean_biomass = mean(microbial_biomass_nmol_g, na.rm = TRUE),
    mean_FB = mean(FB_ratio, na.rm = TRUE),
    mean_stress = mean(stress_index, na.rm = TRUE),
    mean_GPGN = mean(GP_GN_ratio, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

temporal_stability <- site_year_means %>%
  group_by(siteID) %>%
  summarize(
    n_years = n(),
    cv_biomass = 100 * sd(mean_biomass, na.rm = TRUE) / mean(mean_biomass, na.rm = TRUE),
    cv_FB = 100 * sd(mean_FB, na.rm = TRUE) / mean(mean_FB, na.rm = TRUE),
    cv_stress = 100 * sd(mean_stress, na.rm = TRUE) / mean(mean_stress, na.rm = TRUE),
    cv_GPGN = 100 * sd(mean_GPGN, na.rm = TRUE) / mean(mean_GPGN, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_years >= 3)

write.csv(temporal_stability,
          "manuscript_outputs/tables/temporal_stability_cv.csv",
          row.names = FALSE)

cv_summary <- data.frame(
  metric = c("Biomass", "F:B Ratio", "Stress Index", "GP:GN Ratio"),
  n_sites = c(
    sum(!is.na(temporal_stability$cv_biomass)),
    sum(!is.na(temporal_stability$cv_FB)),
    sum(!is.na(temporal_stability$cv_stress)),
    sum(!is.na(temporal_stability$cv_GPGN))
  ),
  median_cv = c(
    round(median(temporal_stability$cv_biomass, na.rm = TRUE), 1),
    round(median(temporal_stability$cv_FB, na.rm = TRUE), 1),
    round(median(temporal_stability$cv_stress, na.rm = TRUE), 1),
    round(median(temporal_stability$cv_GPGN, na.rm = TRUE), 1)
  ),
  mean_cv = c(
    round(mean(temporal_stability$cv_biomass, na.rm = TRUE), 1),
    round(mean(temporal_stability$cv_FB, na.rm = TRUE), 1),
    round(mean(temporal_stability$cv_stress, na.rm = TRUE), 1),
    round(mean(temporal_stability$cv_GPGN, na.rm = TRUE), 1)
  )
)

write.csv(cv_summary,
          "manuscript_outputs/tables/cv_summary_stats.csv",
          row.names = FALSE)

cat("Temporal stability: n =", nrow(temporal_stability), "sites with ≥3 years\n")
cat("Median CVs: Biomass =", cv_summary$median_cv[1], "%, F:B =", cv_summary$median_cv[2],
    "%, Stress =", cv_summary$median_cv[3], "%, GP:GN =", cv_summary$median_cv[4], "%\n")

calc_icc <- function(data, value_col) {
  model_data <- data %>%
    select(siteID, year, value = all_of(value_col)) %>%
    filter(!is.na(value))

  if(nrow(model_data) < 10) return(NA)

  site_means <- model_data %>%
    group_by(siteID) %>%
    summarize(site_mean = mean(value, na.rm = TRUE), .groups = "drop")

  grand_mean <- mean(model_data$value, na.rm = TRUE)
  between_var <- var(site_means$site_mean)

  within_var <- model_data %>%
    left_join(site_means, by = "siteID") %>%
    mutate(dev = (value - site_mean)^2) %>%
    pull(dev) %>%
    mean(na.rm = TRUE)

  icc <- between_var / (between_var + within_var)
  return(icc)
}

icc_biomass <- calc_icc(site_year_means, "mean_biomass")
icc_FB <- calc_icc(site_year_means, "mean_FB")
icc_stress <- calc_icc(site_year_means, "mean_stress")
icc_GPGN <- calc_icc(site_year_means, "mean_GPGN")

icc_results <- data.frame(
  metric = c("Biomass", "F:B Ratio", "Stress Index", "GP:GN Ratio"),
  ICC = round(c(icc_biomass, icc_FB, icc_stress, icc_GPGN), 3)
)

write.csv(icc_results,
          "manuscript_outputs/tables/icc_results.csv",
          row.names = FALSE)

cat("ICCs: Biomass =", round(icc_biomass, 3), ", F:B =", round(icc_FB, 3),
    ", Stress =", round(icc_stress, 3), ", GP:GN =", round(icc_GPGN, 3), "\n")

cv_long <- temporal_stability %>%
  select(siteID, n_years, cv_biomass, cv_FB, cv_stress, cv_GPGN) %>%
  pivot_longer(cols = starts_with("cv_"),
               names_to = "metric",
               values_to = "cv") %>%
  mutate(metric = recode(metric,
                         "cv_biomass" = "Biomass",
                         "cv_FB" = "F:B Ratio",
                         "cv_stress" = "Stress Index",
                         "cv_GPGN" = "GP:GN Ratio"))

fig5b <- ggplot(cv_long, aes(x = metric, y = cv)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7, outlier.alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  geom_hline(yintercept = c(20, 30, 40), linetype = "dashed",
             color = "gray60", linewidth = 0.3) +
  labs(x = "", y = "Coefficient of Variation (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("manuscript_outputs/figures/Figure5b_temporal_stability.png", fig5b,
       width = 8, height = 5, dpi = 300)
ggsave("manuscript_outputs/figures/Figure5b_temporal_stability.pdf", fig5b,
       width = 8, height = 5)

cat("✓ Figure 5 saved\n")


cat("Generating Figure 6: PLFA Metric Relationships with statistics...\n")

get_stats_label <- function(data, x_var, y_var) {
  df <- data.frame(x = data[[x_var]], y = data[[y_var]]) %>%
    filter(!is.na(x) & !is.na(y) & is.finite(x) & is.finite(y))

  if (nrow(df) < 3) return(list(label = "Insufficient data", n = 0))

  fit <- lm(y ~ x, data = df)
  fit_summary <- summary(fit)

  r2 <- round(fit_summary$r.squared, 3)
  pval <- fit_summary$coefficients["x", "Pr(>|t|)"]
  pval_text <- if (pval < 0.001) "p < 0.001" else paste0("p = ", round(pval, 3))

  list(
    label = paste0("N = ", nrow(df), "\nR² = ", r2, "\n", pval_text),
    n = nrow(df),
    r2 = r2,
    pval = pval
  )
}

data_6a <- neon_plfa_clean %>%
  mutate(log_biomass = log10(microbial_biomass_nmol_g)) %>%
  filter(!is.na(log_biomass) & !is.na(FB_ratio) & is.finite(log_biomass))

stats_6a <- get_stats_label(data_6a, "log_biomass", "FB_ratio")

fig6a <- ggplot(data_6a, aes(x = log_biomass, y = FB_ratio)) +
  geom_point(alpha = 0.2, size = 0.8, color = "steelblue") +
  geom_smooth(method = "lm", color = "darkblue", se = TRUE, fill = "lightblue") +
  annotate("text", x = 0.8, y = 2.1,
           label = stats_6a$label, hjust = 0, vjust = 1, size = 3.5,
           fontface = "bold") +
  labs(x = "log10(Biomass) [nmol/g]", y = "F:B Ratio") +
  theme_minimal()

data_6b <- neon_plfa_clean %>%
  mutate(log_bacterial = log10(bacterial_plfa_nmol_g),
         log_fungal = log10(fungal_plfa_nmol_g + 0.1)) %>%
  filter(!is.na(log_bacterial) & !is.na(log_fungal) &
         is.finite(log_bacterial) & is.finite(log_fungal))

stats_6b <- get_stats_label(data_6b, "log_bacterial", "log_fungal")

fig6b <- ggplot(data_6b, aes(x = log_bacterial, y = log_fungal)) +
  geom_point(alpha = 0.2, size = 0.8, color = "darkgreen") +
  geom_smooth(method = "lm", color = "darkgreen", se = TRUE, fill = "lightgreen") +
  annotate("text", x = 0.2, y = 2.7,
           label = stats_6b$label, hjust = 0, vjust = 1, size = 3.5,
           fontface = "bold") +
  labs(x = "log10(Bacterial PLFA) [nmol/g]",
       y = "log10(Fungal PLFA + 0.1) [nmol/g]") +
  theme_minimal()

data_6c <- neon_plfa_clean %>%
  filter(!is.na(stress_index) & !is.na(FB_ratio) & stress_index < 5)

stats_6c <- get_stats_label(data_6c, "FB_ratio", "stress_index")

fig6c <- ggplot(data_6c, aes(x = FB_ratio, y = stress_index)) +
  geom_point(alpha = 0.2, size = 0.8, color = "coral") +
  geom_smooth(method = "lm", color = "red4", se = TRUE, fill = "pink") +
  annotate("text", x = 0.05, y = 4.7,
           label = stats_6c$label, hjust = 0, vjust = 1, size = 3.5,
           fontface = "bold") +
  labs(x = "F:B Ratio", y = "Stress Index") +
  theme_minimal()

fig6_combined <- grid.arrange(fig6a, fig6b, fig6c, ncol = 1)

ggsave("manuscript_outputs/figures/Figure6_plfa_relationships.png", fig6_combined,
       width = 8, height = 11, dpi = 300)
ggsave("manuscript_outputs/figures/Figure6_plfa_relationships.pdf", fig6_combined,
       width = 8, height = 11)

fig6_stats <- data.frame(
  panel = c("(a) F:B vs Biomass", "(b) Fungal vs Bacterial", "(c) Stress vs F:B"),
  n_samples = c(stats_6a$n, stats_6b$n, stats_6c$n),
  r_squared = c(stats_6a$r2, stats_6b$r2, stats_6c$r2),
  p_value = c(stats_6a$pval, stats_6b$pval, stats_6c$pval)
)
write.csv(fig6_stats, "manuscript_outputs/tables/figure6_relationship_stats.csv", row.names = FALSE)

cat("✓ Figure 6 saved\n")


cat("\n=== Saving R workspace ===\n")
save.image("manuscript_outputs/neon_plfa_analysis_workspace.RData")
cat("✓ Saved workspace\n")