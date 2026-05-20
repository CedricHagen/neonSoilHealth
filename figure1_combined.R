
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)
library(viridis)
library(maps)

north_america <- map_data("world", region = c("USA", "Canada", "Mexico", "Puerto Rico"))

fig1a <- ggplot() +
  geom_polygon(data = north_america,
               aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(data = site_metadata,
             aes(x = longitude, y = latitude, size = n_samples, fill = domain),
             shape = 21, alpha = 0.8, stroke = 0.5) +
  scale_size_continuous(range = c(2, 10), name = "Number of samples",
                        breaks = c(100, 200, 300, 400, 500)) +
  scale_fill_viridis_d(option = "turbo", name = "Domain") +
  coord_fixed(1.3, xlim = c(-180, -50), ylim = c(14, 72)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 14) +  
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 11),  
    axis.title = element_text(size = 13),  
    axis.text = element_text(size = 11),  
    panel.background = element_rect(fill = "aliceblue"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(ncol = 2, override.aes = list(size = 4)))


temporal_long <- temporal_matrix %>%
  pivot_longer(-siteID, names_to = "year", values_to = "n_samples") %>%
  mutate(year = as.numeric(year))

fig1b <- ggplot(temporal_long, aes(x = year, y = siteID, fill = n_samples)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "plasma", name = "Number of samples",
                       breaks = c(0, 50, 100, 200, 300, 400, 500)) +
  scale_x_continuous(breaks = 2017:2024) +
  labs(x = "Year", y = "Site") +
  theme_minimal(base_size = 14) +  
  theme(
    axis.text.y = element_text(size = 7),  
    axis.text.x = element_text(size = 11),  
    axis.title = element_text(size = 13, face = "bold"),  
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 11),  
    panel.grid = element_blank()
  )


fig1c <- ggplot(horizon_summary, aes(x = "", y = n_samples, fill = horizon)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 1) +
  coord_polar("y") +
  scale_fill_manual(
    values = c("M" = "#8B4513", "O" = "#2E8B57"),
    labels = c(
      "M" = paste0("Mineral\n(", horizon_summary$pct_total[horizon_summary$horizon == "M"], "%)"),
      "O" = paste0("Organic\n(", horizon_summary$pct_total[horizon_summary$horizon == "O"], "%)")
    )
  ) +
  labs(fill = "Horizon") +
  theme_void(base_size = 14) +  
  theme(
    legend.title = element_text(size = 13, face = "bold"),  
    legend.text = element_text(size = 12)  
  )


layout_matrix <- rbind(
  c(1, 1, 1, 1), 
  c(2, 2, 3, 3) 
)

fig1a_labeled <- arrangeGrob(
  fig1a,
  top = textGrob("a", x = 0.02, hjust = 0,
                 gp = gpar(fontface = "bold", fontsize = 16))
)

fig1b_labeled <- arrangeGrob(
  fig1b,
  top = textGrob("b", x = 0.02, hjust = 0,
                 gp = gpar(fontface = "bold", fontsize = 16))
)

fig1c_labeled <- arrangeGrob(
  fig1c,
  top = textGrob("c", x = 0.02, hjust = 0,
                 gp = gpar(fontface = "bold", fontsize = 16))
)

fig1_combined <- arrangeGrob(
  fig1a_labeled,
  fig1b_labeled,
  fig1c_labeled,
  layout_matrix = layout_matrix,
  heights = c(1.2, 1)  
)

ggsave("manuscript_outputs/figures/Figure1_combined.png",
       fig1_combined,
       width = 14, height = 11, dpi = 300)

ggsave("manuscript_outputs/figures/Figure1_combined.pdf",
       fig1_combined,
       width = 14, height = 11)
