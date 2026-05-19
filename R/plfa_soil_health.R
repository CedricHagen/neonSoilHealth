# R/plfa_soil_health.R
# NEON Soil Health Explorer - backend functions
# Updated: 2026-04-01 (robust PLFA site metadata parsing + cache validation)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(lubridate)
  library(jsonlite)
})

NEON_API_BASE <- "https://data.neonscience.org/api/v0"
DP_PLFA <- "DP1.10104.001"
DP_SOIL_PERIODIC <- "DP1.10086.001"
DP_MEGAPIT <- "DP1.00096.001"

`%||%` <- function(x, y) if (!is.null(x)) x else y

# -------- Package checks (production-friendly) --------
require_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required R packages: ", paste(missing, collapse = ", "),
      "\nInstall with: install.packages(c(", paste0('"', missing, '"', collapse = ", "), "))",
      call. = FALSE
    )
  }
}

# -------- Robust HTTP fetch (fixes macOS/SSL readLines issues) --------
fetch_url_text <- function(url) {
  # 1) curl (best)
  if (requireNamespace("curl", quietly = TRUE)) {
    h <- curl::new_handle()
    curl::handle_setopt(h, useragent = "R NEON Soil Health Explorer")
    res <- curl::curl_fetch_memory(url, handle = h)
    return(rawToChar(res$content))
  }
  
  # 2) httr (fallback)
  if (requireNamespace("httr", quietly = TRUE)) {
    r <- httr::GET(url, httr::user_agent("R NEON Soil Health Explorer"))
    httr::stop_for_status(r)
    return(httr::content(r, as = "text", encoding = "UTF-8"))
  }
  
  # 3) base readLines (last resort)
  paste(readLines(url, warn = FALSE), collapse = "\n")
}

neon_api_get <- function(path) {
  url <- if (startsWith(path, "http")) path else paste0(NEON_API_BASE, path)
  txt <- tryCatch(fetch_url_text(url), error = function(e) NULL)
  if (is.null(txt) || !nzchar(txt)) stop("NEON API request failed: ", url, call. = FALSE)
  jsonlite::fromJSON(txt, simplifyDataFrame = FALSE)
}

get_neon_release_choices <- function() {
  res <- neon_api_get("/releases")
  rel <- NULL
  
  if (!is.null(res$data$release)) rel <- res$data$release
  if (is.null(rel) && !is.null(res$data)) {
    if (is.list(res$data) && length(res$data) > 0 && !is.null(res$data[[1]]$release)) {
      rel <- vapply(res$data, function(x) x$release, character(1))
    }
  }
  
  rel <- unique(na.omit(rel))
  rel <- rel[grepl("^RELEASE-\\d{4}$", rel)]
  rel <- sort(rel)
  
  unique(c(rel, "current"))
}

default_release <- function(release_choices) {
  yrs <- suppressWarnings(as.integer(sub("^RELEASE-", "", release_choices)))
  if (all(is.na(yrs))) return("current")
  release_choices[which.max(yrs)]
}

# -------- NEW: robust extraction of site codes from possibly nested API structures --------
extract_neon_site_codes <- function(x) {
  if (is.null(x)) return(character(0))
  
  # already clean
  if (is.character(x)) return(x)
  if (is.factor(x)) return(as.character(x))
  
  # list cases (common when API structure changes)
  if (is.list(x)) {
    
    # If it's a list of site objects: each element may have $siteCode
    codes <- purrr::map_chr(x, function(el) {
      if (is.null(el)) return(NA_character_)
      
      if (is.character(el) || is.factor(el)) {
        elc <- as.character(el)
        # return first string that looks like a site code
        hit <- elc[str_detect(elc, "^[A-Z]{4}$")]
        if (length(hit) > 0) return(hit[[1]])
        # sometimes element itself is a single code
        if (length(elc) == 1) return(elc[[1]])
        return(NA_character_)
      }
      
      if (is.list(el)) {
        sc <- el$siteCode %||% el$siteID %||% el$code %||% el$site
        if (is.character(sc) || is.factor(sc)) {
          scc <- as.character(sc)
          hit <- scc[str_detect(scc, "^[A-Z]{4}$")]
          if (length(hit) > 0) return(hit[[1]])
        }
      }
      
      NA_character_
    })
    
    codes <- codes[!is.na(codes)]
    if (length(codes) > 0) return(codes)
    
    # Fallback: search recursively for 4-letter site codes
    flat <- unlist(x, recursive = TRUE, use.names = FALSE)
    flat <- as.character(flat)
    flat[str_detect(flat, "^[A-Z]{4}$")]
  } else {
    # unknown type
    as.character(x)
  }
}

extract_site_codes_from_product <- function(prod) {
  raw <- prod$data$siteCodes %||%
    prod$data$product$siteCodes %||%
    prod$data$sites %||%
    prod$data$product$sites
  
  codes <- extract_neon_site_codes(raw)
  codes <- unique(na.omit(codes))
  codes <- codes[str_detect(codes, "^[A-Z]{4}$")]
  codes <- sort(unique(codes))
  
  # sanity check: if absurdly large, something is wrong
  if (length(codes) > 200) {
    # try a stricter recursive extract from whole product payload
    flat <- unlist(prod, recursive = TRUE, use.names = FALSE)
    flat <- as.character(flat)
    codes2 <- unique(flat[str_detect(flat, "^[A-Z]{4}$")])
    codes2 <- sort(unique(codes2))
    if (length(codes2) <= 200 && length(codes2) > 0) codes <- codes2
  }
  
  codes
}

validate_plfa_site_table <- function(df) {
  if (is.null(df) || !is.data.frame(df)) return(FALSE)
  if (!all(c("siteID", "label") %in% names(df))) return(FALSE)
  
  # coerce for validation
  sid <- df$siteID
  if (is.list(sid)) return(FALSE)
  sid <- as.character(sid)
  sid <- sid[!is.na(sid)]
  if (length(sid) == 0) return(FALSE)
  if (!all(str_detect(sid, "^[A-Z]{4}$"))) return(FALSE)
  
  # NEON has < 200 sites; PLFA subset much smaller
  if (nrow(df) > 200) return(FALSE)
  
  TRUE
}

# -------- PLFA site table (robust, cached, never returns empty) --------
get_plfa_site_table <- function(cache_dir = "cache", force_refresh = FALSE) {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  cache_file <- file.path(cache_dir, "plfa_site_table.rds")
  
  if (file.exists(cache_file) && !force_refresh) {
    out <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (validate_plfa_site_table(out)) {
      # light cleanup
      out$siteID <- as.character(out$siteID)
      out$label <- as.character(out$label)
      out <- out %>% filter(str_detect(siteID, "^[A-Z]{4}$")) %>% distinct(siteID, .keep_all = TRUE)
      return(out)
    }
  }
  
  prod <- neon_api_get(paste0("/products/", DP_PLFA))
  site_codes <- extract_site_codes_from_product(prod)
  
  if (length(site_codes) == 0) {
    stop("Failed to extract PLFA site codes from NEON product metadata (DP1.10104.001).", call. = FALSE)
  }
  
  sites <- purrr::map(site_codes, function(sc) {
    s <- tryCatch(neon_api_get(paste0("/sites/", sc)), error = function(e) NULL)
    
    if (is.null(s) || is.null(s$data)) {
      return(tibble::tibble(
        siteID = as.character(sc),
        siteName = NA_character_,
        domainID = NA_character_,
        domainName = NA_character_,
        stateCode = NA_character_,
        stateName = NA_character_,
        latitude = NA_real_,
        longitude = NA_real_,
        siteType = NA_character_
      ))
    }
    
    d <- s$data
    tibble::tibble(
      siteID = as.character(d$siteCode %||% d$siteID %||% sc),
      siteName = d$siteName %||% NA_character_,
      domainID = d$domainCode %||% d$domainID %||% NA_character_,
      domainName = d$domainName %||% NA_character_,
      stateCode = d$stateCode %||% NA_character_,
      stateName = d$stateName %||% NA_character_,
      latitude = suppressWarnings(as.numeric(d$siteLatitude %||% d$siteLatitudeDecimal %||% d$latitude)),
      longitude = suppressWarnings(as.numeric(d$siteLongitude %||% d$siteLongitudeDecimal %||% d$longitude)),
      siteType = d$siteType %||% NA_character_
    )
  }) %>% bind_rows()
  
  sites <- sites %>%
    filter(str_detect(siteID, "^[A-Z]{4}$")) %>%
    distinct(siteID, .keep_all = TRUE) %>%
    mutate(
      domain_num = suppressWarnings(as.integer(str_extract(domainID, "\\d+"))),
      domain_num = ifelse(is.na(domain_num), 999L, domain_num),
      domainID_show = ifelse(is.na(domainID) | !nzchar(domainID), "D??", domainID),
      siteName_show = ifelse(is.na(siteName) | !nzchar(siteName), "Unknown site name", siteName),
      label = sprintf("%s — %s (%s)", domainID_show, siteID, siteName_show)
    ) %>%
    arrange(domain_num, siteID)
  
  saveRDS(sites, cache_file)
  sites
}

# -------- NEON download helpers --------
download_neon_product <- function(dpID, site, release,
                                  startdate = NULL, enddate = NULL,
                                  package = "basic",
                                  check.size = FALSE,
                                  cache_dir = "cache/neon_downloads") {
  
  require_pkgs(c("neonUtilities", "withr"))
  
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  
  args <- list(
    dpID = dpID,
    site = site,
    release = release,
    package = package,
    check.size = check.size
  )
  
  if (!is.null(startdate) && nzchar(startdate)) args$startdate <- startdate
  if (!is.null(enddate) && nzchar(enddate)) args$enddate <- enddate
  
  withr::with_dir(cache_dir, {
    do.call(neonUtilities::loadByProduct, args)
  })
}

standardize_collect_datetime <- function(df) {
  if (!is.data.frame(df)) return(df)
  
  dt <- NULL
  if ("collectDate" %in% names(df)) dt <- df$collectDate
  if (is.null(dt) && "collectDateTime" %in% names(df)) dt <- df$collectDateTime
  if (is.null(dt) && "startDateTime" %in% names(df)) dt <- df$startDateTime
  
  if (is.null(dt)) {
    df$collectDateTime <- as.POSIXct(NA)
    df$collectDate <- as.Date(NA)
    return(df)
  }
  
  if (inherits(dt, "POSIXct")) {
    df$collectDateTime <- dt
  } else {
    df$collectDateTime <- suppressWarnings(lubridate::ymd_hms(dt, tz = "UTC"))
    if (all(is.na(df$collectDateTime))) df$collectDateTime <- suppressWarnings(lubridate::ymd_hm(dt, tz = "UTC"))
    if (all(is.na(df$collectDateTime))) df$collectDateTime <- suppressWarnings(lubridate::ymd(dt, tz = "UTC"))
  }
  
  df$collectDate <- as.Date(df$collectDateTime)
  df$year <- lubridate::year(df$collectDateTime)
  df$month <- lubridate::month(df$collectDateTime)
  df
}

infer_horizon_from_sampleID <- function(sampleID) {
  ifelse(grepl("-O-", sampleID), "O",
         ifelse(grepl("-M-", sampleID), "M", NA_character_))
}

best_col <- function(df, candidates) {
  candidates <- candidates[candidates %in% names(df)]
  if (length(candidates) == 0) return(NA_character_)
  candidates[[1]]
}

lipid_col <- function(df, base_conc_name) {
  scaled <- sub("Concentration$", "ScaledConcentration", base_conc_name)
  cand <- c(scaled, base_conc_name)
  best_col(df, cand)
}

# -------- Pull + merge data products --------
pull_plfa_table <- function(site, release, startdate = NULL, enddate = NULL, cache_dir = "cache/neon_downloads") {
  microb <- download_neon_product(DP_PLFA, site = site, release = release,
                                  startdate = startdate, enddate = enddate,
                                  package = "basic", check.size = FALSE,
                                  cache_dir = cache_dir)
  
  if (!is.null(microb$sme_scaledMicrobialBiomass)) {
    df <- microb$sme_scaledMicrobialBiomass
    df$.plfa_table_used <- "sme_scaledMicrobialBiomass"
  } else if (!is.null(microb$sme_microbialBiomass)) {
    df <- microb$sme_microbialBiomass
    df$.plfa_table_used <- "sme_microbialBiomass"
  } else {
    stop("Could not find sme_scaledMicrobialBiomass or sme_microbialBiomass in DP1.10104.001 download.", call. = FALSE)
  }
  
  df <- standardize_collect_datetime(df)
  
  if (!("siteID" %in% names(df))) df$siteID <- site
  if (!("sampleID" %in% names(df))) stop("PLFA table missing sampleID.", call. = FALSE)
  
  if (!("horizon" %in% names(df))) df$horizon <- infer_horizon_from_sampleID(df$sampleID)
  df
}

pull_soil_periodic_tables <- function(site, release, startdate = NULL, enddate = NULL, cache_dir = "cache/neon_downloads") {
  soil <- download_neon_product(DP_SOIL_PERIODIC, site = site, release = release,
                                startdate = startdate, enddate = enddate,
                                package = "basic", check.size = FALSE,
                                cache_dir = cache_dir)
  
  core <- soil$sls_soilCoreCollection %||% tibble::tibble()
  moist <- soil$sls_soilMoisture %||% tibble::tibble()
  chem <- soil$sls_soilChemistry %||% tibble::tibble()
  
  core <- standardize_collect_datetime(core)
  moist <- standardize_collect_datetime(moist)
  chem <- standardize_collect_datetime(chem)
  
  list(core = core, moisture = moist, chemistry = chem)
}

pull_megapit_soil_type <- function(site, release, cache_dir = "cache/neon_downloads") {
  out <- list(soilType = NA_character_, soilOrder = NA_character_, soilSuborder = NA_character_,
              soilGreatGroup = NA_character_, soilSubgroup = NA_character_)
  
  mega <- tryCatch(
    download_neon_product(DP_MEGAPIT, site = site, release = release,
                          package = "basic", check.size = FALSE, cache_dir = cache_dir),
    error = function(e) NULL
  )
  if (is.null(mega)) return(out)
  
  tbls <- mega[names(mega) != "citation"]
  
  pick <- NULL
  for (nm in names(tbls)) {
    df <- tbls[[nm]]
    if (!is.data.frame(df)) next
    nms <- names(df)
    if (any(grepl("soilTaxon|taxon|soilOrder|soilSuborder|soilGreatGroup|soilSubgroup", nms, ignore.case = TRUE))) {
      pick <- df
      break
    }
  }
  if (is.null(pick) || !is.data.frame(pick) || nrow(pick) == 0) return(out)
  
  nms <- names(pick)
  get_first <- function(regex) {
    hit <- nms[grepl(regex, nms, ignore.case = TRUE)]
    if (length(hit) == 0) return(NA_character_)
    val <- pick[[hit[[1]]]][[1]]
    if (is.null(val)) return(NA_character_)
    as.character(val)
  }
  
  soilOrder <- get_first("^soil.*order$|soilOrder|taxonOrder")
  soilSuborder <- get_first("suborder")
  soilGreatGroup <- get_first("greatgroup|greatGroup")
  soilSubgroup <- get_first("subgroup|subGroup")
  
  parts <- c(
    if (!is.na(soilOrder) && nzchar(soilOrder)) paste0("Order: ", soilOrder) else NULL,
    if (!is.na(soilSuborder) && nzchar(soilSuborder)) paste0("Suborder: ", soilSuborder) else NULL,
    if (!is.na(soilGreatGroup) && nzchar(soilGreatGroup)) paste0("Great group: ", soilGreatGroup) else NULL,
    if (!is.na(soilSubgroup) && nzchar(soilSubgroup)) paste0("Subgroup: ", soilSubgroup) else NULL
  )
  soilType <- if (length(parts) == 0) NA_character_ else paste(parts, collapse = " | ")
  
  list(
    soilType = soilType,
    soilOrder = soilOrder,
    soilSuborder = soilSuborder,
    soilGreatGroup = soilGreatGroup,
    soilSubgroup = soilSubgroup
  )
}

merge_plfa_soil <- function(plfa, soil_tables) {
  core <- soil_tables$core
  moist <- soil_tables$moisture
  chem <- soil_tables$chemistry
  
  soilTemp_col <- best_col(core, c("soilTemp", "soilTemperature"))
  soilMoist_col <- best_col(moist, c("soilMoisture", "soilWaterContent", "volumetricWaterContent", "vwc", "waterContent"))
  
  core_keep <- core %>%
    select(any_of(c("sampleID", "horizon", "plotID", "collectDateTime", "collectDate", "year", "month", soilTemp_col))) %>%
    rename(soilTemp = !!soilTemp_col)
  
  moist_keep <- moist %>%
    select(any_of(c("sampleID", "horizon", "plotID", "collectDateTime", "collectDate", "year", "month", soilMoist_col))) %>%
    rename(soilMoisture = !!soilMoist_col)
  
  chem_keep <- chem %>%
    select(any_of(c("sampleID", "horizon", "plotID", "collectDateTime", "collectDate", "year", "month",
                    "organicCPercent", "nitrogenPercent", "organicd13C", "d15N"))) %>%
    distinct(sampleID, .keep_all = TRUE)
  
  out <- plfa %>%
    left_join(core_keep, by = "sampleID", suffix = c("", ".core")) %>%
    left_join(moist_keep, by = "sampleID", suffix = c("", ".moist")) %>%
    left_join(chem_keep, by = "sampleID", suffix = c("", ".chem"))
  
  if ("horizon.core" %in% names(out) && "horizon" %in% names(out)) {
    out <- out %>% mutate(horizon = coalesce(horizon, horizon.core))
  }
  
  out
}

# -------- Metric calculations --------
compute_soil_health_metrics <- function(df) {
  total_col <- best_col(df, c("totalLipidScaledConcentration", "totalLipidConcentration"))
  c18_col <- best_col(df, c("c18To0ScaledConcentration", "c18To0Concentration"))
  
  if (is.na(total_col)) {
    df$microbial_biomass_nmol_g <- NA_real_
  } else if (is.na(c18_col)) {
    df$microbial_biomass_nmol_g <- suppressWarnings(as.numeric(df[[total_col]]))
  } else {
    total <- suppressWarnings(as.numeric(df[[total_col]]))
    c18 <- suppressWarnings(as.numeric(df[[c18_col]]))
    df$microbial_biomass_nmol_g <- ifelse(is.na(c18), total, total - c18)
    df$microbial_biomass_nmol_g <- pmax(df$microbial_biomass_nmol_g, 0)
  }
  
  df$total_microbial_biomass <- df$microbial_biomass_nmol_g
  
  fungal_cand <- c(
    lipid_col(df, "cis18To2n912Concentration"),
    lipid_col(df, "trans18To2n912Concentration")
  )
  fungal_cand <- fungal_cand[!is.na(fungal_cand)]
  fungal_col <- if (length(fungal_cand) > 0) fungal_cand[[1]] else NA_character_
  
  bacterial_base <- c(
    "i15To0Concentration",
    "aC15To0Concentration",
    "i16To0Concentration",
    "i17To0Concentration",
    "c17To0AnteisoConcentration",
    "cis16To1n9Concentration",
    "c18To1n11Concentration",
    "cyclo17To0Concentration",
    "cyclo19To0Concentration",
    "lipid10Methyl16To0Concentration",
    "lipid10Methyl17To0Concentration",
    "lipid10Methyl18To0Concentration"
  )
  bacterial_cols <- purrr::map_chr(bacterial_base, ~ lipid_col(df, .x))
  bacterial_cols <- bacterial_cols[!is.na(bacterial_cols)]
  
  cyclo_cols <- c(lipid_col(df, "cyclo17To0Concentration"),
                  lipid_col(df, "cyclo19To0Concentration"))
  cyclo_cols <- cyclo_cols[!is.na(cyclo_cols)]
  precursor_cols <- c(lipid_col(df, "cis16To1n9Concentration"),
                      lipid_col(df, "c18To1n11Concentration"))
  precursor_cols <- precursor_cols[!is.na(precursor_cols)]
  
  df$fungal_plfa_nmol_g <- if (!is.na(fungal_col)) suppressWarnings(as.numeric(df[[fungal_col]])) else NA_real_
  
  if (length(bacterial_cols) > 0) {
    df$bacterial_plfa_nmol_g <- df %>%
      select(all_of(bacterial_cols)) %>%
      mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
      rowSums(na.rm = TRUE)
    df$bacterial_plfa_nmol_g <- ifelse(df$bacterial_plfa_nmol_g == 0, NA_real_, df$bacterial_plfa_nmol_g)
  } else {
    df$bacterial_plfa_nmol_g <- NA_real_
  }
  
  df$FB_ratio <- ifelse(
    is.na(df$fungal_plfa_nmol_g) | is.na(df$bacterial_plfa_nmol_g) | df$bacterial_plfa_nmol_g <= 0,
    NA_real_,
    df$fungal_plfa_nmol_g / df$bacterial_plfa_nmol_g
  )
  
  if (length(cyclo_cols) > 0 && length(precursor_cols) > 0) {
    cyclo_sum <- df %>%
      select(all_of(cyclo_cols)) %>%
      mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
      rowSums(na.rm = TRUE)
    
    prec_sum <- df %>%
      select(all_of(precursor_cols)) %>%
      mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
      rowSums(na.rm = TRUE)
    
    df$stress_index <- ifelse(prec_sum > 0, cyclo_sum / prec_sum, NA_real_)
  } else {
    df$stress_index <- NA_real_
  }
  
  df$organicCPercent <- suppressWarnings(as.numeric(df$organicCPercent))
  df$nitrogenPercent <- suppressWarnings(as.numeric(df$nitrogenPercent))
  df$CN_ratio <- ifelse(!is.na(df$organicCPercent) & !is.na(df$nitrogenPercent) & df$nitrogenPercent > 0,
                        df$organicCPercent / df$nitrogenPercent,
                        NA_real_)
  
  df$biomass_per_gC_nmol_gC <- ifelse(!is.na(df$organicCPercent) & df$organicCPercent > 0,
                                      df$microbial_biomass_nmol_g * 100 / df$organicCPercent,
                                      NA_real_)
  
  attr(df, "marker_cols_used") <- list(
    fungal_col = fungal_col,
    bacterial_cols = bacterial_cols,
    cyclo_cols = cyclo_cols,
    precursor_cols = precursor_cols,
    total_col = total_col,
    c18_col = c18_col
  )
  
  df
}

compute_trends <- function(df, metrics, group_vars = c("siteID", "horizon")) {
  require_pkgs(c("broom"))
  
  df <- df %>% filter(!is.na(year))
  out <- list()
  
  for (m in metrics) {
    if (!(m %in% names(df))) next
    dd <- df %>%
      select(all_of(c(group_vars, "year", m))) %>%
      mutate(value = suppressWarnings(as.numeric(.data[[m]]))) %>%
      filter(!is.na(value)) %>%
      group_by(across(all_of(group_vars)), year) %>%
      summarize(value = mean(value, na.rm = TRUE), .groups = "drop")
    
    if (nrow(dd) == 0) next
    
    res <- dd %>%
      group_by(across(all_of(group_vars))) %>%
      group_modify(function(dat, keys) {
        if (dplyr::n_distinct(dat$year) < 3) {
          return(tibble::tibble(
            metric = m,
            slope_per_year = NA_real_,
            p_value = NA_real_,
            r2 = NA_real_,
            n_years = dplyr::n_distinct(dat$year)
          ))
        }
        fit <- lm(value ~ year, data = dat)
        sm <- summary(fit)
        tibble::tibble(
          metric = m,
          slope_per_year = unname(coef(fit)[["year"]]),
          p_value = unname(coef(sm)[2, 4]),
          r2 = sm$r.squared,
          n_years = dplyr::n_distinct(dat$year)
        )
      }) %>%
      ungroup()
    
    out[[m]] <- res
  }
  
  bind_rows(out)
}

run_neon_soil_health <- function(sites, release,
                                 start_date = NULL, end_date = NULL,
                                 cache_dir = "cache/neon_downloads") {
  
  diagnostics <- list(
    release_requested = release,
    sites_requested = sites
  )
  
  all <- list()
  
  for (s in sites) {
    plfa <- pull_plfa_table(s, release, startdate = start_date, enddate = end_date, cache_dir = cache_dir)
    soil <- pull_soil_periodic_tables(s, release, startdate = start_date, enddate = end_date, cache_dir = cache_dir)
    
    merged <- merge_plfa_soil(plfa, soil)
    
    st <- pull_megapit_soil_type(s, release, cache_dir = cache_dir)
    merged$soilType <- st$soilType
    merged$soilOrder <- st$soilOrder
    merged$soilSuborder <- st$soilSuborder
    merged$soilGreatGroup <- st$soilGreatGroup
    merged$soilSubgroup <- st$soilSubgroup
    
    merged <- compute_soil_health_metrics(merged)
    
    merged <- merged %>%
      mutate(
        siteID = coalesce(siteID, s),
        horizon = coalesce(horizon, infer_horizon_from_sampleID(sampleID))
      )
    
    all[[s]] <- merged
  }
  
  out <- bind_rows(all)
  
  diagnostics$n_rows <- nrow(out)
  diagnostics$n_sites_in_output <- dplyr::n_distinct(out$siteID)
  
  marker <- NULL
  if (length(all) > 0) marker <- attr(all[[1]], "marker_cols_used")
  diagnostics$marker_cols_used <- marker %||% list()
  
  diagnostics$FB_all_na <- all(is.na(out$FB_ratio))
  diagnostics$stress_all_na <- all(is.na(out$stress_index))
  diagnostics$soc_all_na <- all(is.na(out$organicCPercent))
  
  list(data = out, diagnostics = diagnostics)
}