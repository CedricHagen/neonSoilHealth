
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)


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
            hjust = -0.1, size = 3.8) + 
  coord_flip() +
  ylim(0, 105) +
  labs(x = "", y = "% Complete") +
  theme_minimal(base_size = 14) +  
  theme(
    axis.title = element_text(size = 13, face = "bold"),  
    axis.text.y = element_text(size = 12),  
    axis.text.x = element_text(size = 11),  
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 20, 10, 10)  
  )

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
  geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +  
  geom_hline(yintercept = c(20, 30, 40), linetype = "dashed",
             color = "gray60", linewidth = 0.3) +
  labs(x = "", y = "Coefficient of Variation (%)") +
  theme_minimal(base_size = 14) +  
  theme(
    axis.title = element_text(size = 13, face = "bold"),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),  
    axis.text.y = element_text(size = 11),  
    plot.margin = margin(10, 10, 10, 20)  
  )


fig5a_labeled <- arrangeGrob(
  fig5a,
  top = textGrob("a", x = 0.02, hjust = 0,
                 gp = gpar(fontface = "bold", fontsize = 16))
)

fig5b_labeled <- arrangeGrob(
  fig5b,
  top = textGrob("b", x = 0.02, hjust = 0,
                 gp = gpar(fontface = "bold", fontsize = 16))
)

fig5_combined <- arrangeGrob(
  fig5a_labeled,
  fig5b_labeled,
  ncol = 2,  
  widths = c(1, 1)  
)

ggsave("manuscript_outputs/figures/Figure5_combined.png",
       fig5_combined,
       width = 14, height = 6, dpi = 300)

ggsave("manuscript_outputs/figures/Figure5_combined.pdf",
       fig5_combined,
       width = 14, height = 6)
