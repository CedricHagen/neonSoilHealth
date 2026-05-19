source("plfa_soil_health_functions.R")

# Install if needed (recommended to do once, not repeatedly)
# install.packages(c("neonUtilities","withr","dplyr","jsonlite"))

sites <- c("STER", "RMNP")
start_date <- as.Date("2018-01-01")
end_date   <- Sys.Date()

res <- run_plfa_soil_health(
  sites = sites,
  start_date = start_date,
  end_date = end_date,
  release = "current",
  include_provisional = FALSE,
  package = "basic",
  data_dir = "neon_downloads",
  quiet = FALSE
)

df <- res$data
str(df)
summary(df$microbial_biomass_nmol_g)
