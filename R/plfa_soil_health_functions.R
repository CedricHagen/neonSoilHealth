# plfa_soil_health_functions.R
# ============================================================
# NEON PLFA-based soil health indicators across sites + time
#
# Data products used:
#   DP1.10104.001  Soil Microbe Biomass (PLFA)
#   DP1.10086.001  Soil physical and chemical properties, periodic
#
# Outputs (per sample):
#   - microbial_biomass_nmol_g (NEON correctedTotLipidConc)
#   - fungal_plfa_nmol_g
#   - bacterial_plfa_nmol_g
#   - fb_ratio
#   - stress_index_cy_pre
#   - soc_g_g (SOC as g C / g dry soil)
#   - biomass_per_soc_nmol_gC
# ============================================================

check_required_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing packages: ", paste(missing, collapse = ", "),
      "\nInstall with: install.packages(c(", paste(sprintf('"%s"', missing), collapse = ", "), "))",
      call. = FALSE
    )
  }
}

pick_first_existing <- function(x, candidates) {
  # x: character vector (typically names(df))
  for (cnd in candidates) {
    if (cnd %in% x) return(cnd)
  }
  return(NA_character_)
}

row_sum_any_of <- function(df, cols) {
  cols2 <- intersect(cols, names(df))
  if (length(cols2) == 0) return(rep(NA_real_, nrow(df)))
  m <- as.matrix(df[, cols2, drop = FALSE])
  suppressWarnings(storage.mode(m) <- "numeric")
  out <- rowSums(m, na.rm = TRUE)
  # If all were NA in a row, rowSums gives 0. Convert those to NA to avoid false zeros.
  all_na <- apply(is.na(df[, cols2, drop = FALSE]), 1, all)
  out[all_na] <- NA_real_
  out
}

as_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) return(as.Date(x))
  if (is.character(x)) {
    # handle "YYYY-MM-DD" or "YYYY-MM-DDTHH:MM:SSZ"
    return(as.Date(substr(x, 1, 10)))
  }
  as.Date(x)
}

format_ym <- function(d) {
  d <- as_date_safe(d)
  format(d, "%Y-%m")
}

# --- NEON site metadata for map + UI choices (API call) -------
get_neon_sites_api <- function() {
  check_required_pkgs(c("jsonlite"))
  url <- "https://data.neonscience.org/api/v0/sites"
  x <- jsonlite::fromJSON(url, flatten = TRUE)
  
  dat <- x$data
  if (is.null(dat)) stop("NEON sites API returned no 'data' field.", call. = FALSE)
  
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  
  code_col <- pick_first_existing(names(dat), c("siteCode", "siteID", "site"))
  name_col <- pick_first_existing(names(dat), c("siteName", "name"))
  lat_col  <- pick_first_existing(names(dat), c("latitude", "decimalLatitude", "location.latitude", "location.decimalLatitude"))
  lon_col  <- pick_first_existing(names(dat), c("longitude", "decimalLongitude", "location.longitude", "location.decimalLongitude"))
  dom_col  <- pick_first_existing(names(dat), c("domainCode", "domainID"))
  state_col <- pick_first_existing(names(dat), c("stateCode", "state"))
  
  if (any(is.na(c(code_col, name_col, lat_col, lon_col)))) {
    stop("Could not standardize NEON site metadata from API response. Columns seen: ",
         paste(names(dat), collapse = ", "), call. = FALSE)
  }
  
  out <- data.frame(
    siteCode = dat[[code_col]],
    siteName = dat[[name_col]],
    latitude = dat[[lat_col]],
    longitude = dat[[lon_col]],
    domainCode = if (!is.na(dom_col)) dat[[dom_col]] else NA_character_,
    stateCode  = if (!is.na(state_col)) dat[[state_col]] else NA_character_,
    stringsAsFactors = FALSE
  )
  
  # Drop rows with missing coordinates
  out <- out[!is.na(out$latitude) & !is.na(out$longitude), ]
  out
}

# --- Download NEON data (uses neonUtilities::loadByProduct) ---
download_neon_products <- function(
    sites,
    start_date,
    end_date,
    release = "current",
    include_provisional = FALSE,
    package = "basic",
    data_dir = "neon_downloads",
    token = Sys.getenv("NEON_TOKEN", unset = NA_character_),
    quiet = TRUE
) {
  check_required_pkgs(c("neonUtilities", "withr"))
  
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  
  start_ym <- format_ym(start_date)
  end_ym   <- format_ym(end_date)
  
  # loadByProduct downloads into working directory; isolate into data_dir
  res <- withr::with_dir(data_dir, {
    if (!quiet) message("Downloading DP1.10104.001 (Soil Microbe Biomass / PLFA)...")
    microb <- neonUtilities::loadByProduct(
      dpID = "DP1.10104.001",
      site = sites,
      startdate = start_ym,
      enddate = end_ym,
      package = package,
      release = release,
      include.provisional = include_provisional,
      check.size = FALSE,
      token = token
    )
    
    if (!quiet) message("Downloading DP1.10086.001 (Soil physical/chemical properties)...")
    soil <- neonUtilities::loadByProduct(
      dpID = "DP1.10086.001",
      site = sites,
      startdate = start_ym,
      enddate = end_ym,
      package = package,
      release = release,
      include.provisional = include_provisional,
      check.size = FALSE,
      token = token
    )
    
    list(microb = microb, soil = soil, start_ym = start_ym, end_ym = end_ym)
  })
  
  res
}

# --- Identify SOC field + units from variables_10086 ----------
pick_soc_field_from_variables <- function(soil_vars, cn_df) {
  check_required_pkgs(c("dplyr"))
  
  if (is.null(soil_vars) || !is.data.frame(soil_vars)) return(list(field = NA_character_, units = NA_character_))
  
  # Prefer "organic carbon" fields; avoid isotope fields.
  vars <- soil_vars
  vars$desc_low <- tolower(vars$description)
  
  candidates <- vars |>
    dplyr::filter(
      grepl("organic carbon|total organic carbon|%c|percent c|carbon concentration", desc_low) &
        !grepl("delta|isotope|d13|13c|ratio", desc_low)
    ) |>
    dplyr::pull(fieldName)
  
  candidates <- intersect(candidates, names(cn_df))
  
  # Priority order (if present)
  priority <- c(
    "soilTotalOrganicCarbon",
    "soilOrganicCarbon",
    "totalOrganicCarbon",
    "percentOrganicCarbon",
    "soilPercentOrganicCarbon",
    "percentC",
    "soilTotalCarbon" # fallback; not always organic-only
  )
  
  field <- NA_character_
  for (p in priority) {
    if (p %in% candidates) { field <- p; break }
  }
  if (is.na(field) && length(candidates) > 0) field <- candidates[[1]]
  
  units <- NA_character_
  if (!is.na(field)) {
    units <- vars$units[match(field, vars$fieldName)]
    if (length(units) == 0) units <- NA_character_
    units <- units[[1]]
  }
  
  list(field = field, units = units)
}

convert_soc_to_g_per_g <- function(value, units) {
  # Convert SOC to g C / g dry soil
  if (is.na(units) || is.null(units)) return(rep(NA_real_, length(value)))
  
  if (units %in% c("gramsPerKilogram")) return(value / 1000)
  if (units %in% c("percent", "percentByMass")) return(value / 100)
  if (units %in% c("gramsPerGram")) return(value)
  
  # Unknown units: return NA
  rep(NA_real_, length(value))
}

# --- Extract SOC table from soil list (DP1.10086.001) ---------
extract_soc <- function(soil_list) {
  check_required_pkgs(c("dplyr"))
  
  # Best guess: table name includes carbonNitrogen (per NEON docs)
  cn_name <- names(soil_list)[grepl("carbonNitrogen", names(soil_list), ignore.case = TRUE)][1]
  if (is.na(cn_name)) {
    # fallback: search any table with likely carbon columns
    cn_name <- names(soil_list)[grepl("carbon|nitrogen|cn", names(soil_list), ignore.case = TRUE)][1]
  }
  if (is.na(cn_name)) return(NULL)
  
  cn <- soil_list[[cn_name]]
  if (!is.data.frame(cn)) return(NULL)
  
  soil_vars <- soil_list$variables_10086
  soc_info <- pick_soc_field_from_variables(soil_vars, cn)
  soc_field <- soc_info$field
  soc_units <- soc_info$units
  
  if (is.na(soc_field)) return(NULL)
  
  # Join key preference
  join_key <- pick_first_existing(names(cn), c("sampleID", "cnSampleID", "sampleCode"))
  if (is.na(join_key)) return(NULL)
  
  cn2 <- cn |>
    dplyr::filter(!is.na(.data[[soc_field]]))
  
  if (nrow(cn2) == 0) return(NULL)
  
  # Prefer acidTreatment == Y for carbonate soils if column exists (NEON: organic C is reported; acid runs appear separately)
  if ("acidTreatment" %in% names(cn2)) {
    cn2 <- cn2 |>
      dplyr::mutate(.acid_priority = dplyr::if_else(.data$acidTreatment %in% c("Y", "Yes", "TRUE", "T"), 1L, 0L)) |>
      dplyr::group_by(.data[[join_key]]) |>
      dplyr::filter(.acid_priority == max(.acid_priority, na.rm = TRUE)) |>
      dplyr::ungroup()
  }
  
  # Average replicates if present
  out <- cn2 |>
    dplyr::group_by(.data[[join_key]]) |>
    dplyr::summarise(
      soc_raw = mean(.data[[soc_field]], na.rm = TRUE),
      .groups = "drop"
    )
  
  out$soc_units <- soc_units
  out$soc_g_g <- convert_soc_to_g_per_g(out$soc_raw, out$soc_units)
  
  # Keep join key standardized
  names(out)[names(out) == join_key] <- "join_id"
  
  out
}

# --- Compute PLFA indicators ---------------------------------
compute_plfa_indicators <- function(core_microb_df) {
  # Biomarkers are implemented using common PLFA conventions.
  # The code is robust to NEON column-name differences by checking
  # for multiple possible field names and using those that exist.
  
  df <- core_microb_df
  
  # Fungal markers (commonly: 18:2 ω6,9; often also 18:1 ω9)
  fungal_markers <- c(
    "c18To2n912ScaledConcentration", "cis18To2n912ScaledConcentration",
    "c18To2w6w9ScaledConcentration", "cis18To2w6w9ScaledConcentration",
    "c18To1n9ScaledConcentration", "cis18To1n9ScaledConcentration"
  )
  
  # AMF marker often used: 16:1 ω5c (sometimes represented as 16:1 cis11)
  amf_markers <- c("c16To1Cis11ScaledConcentration", "cis16To1Cis11ScaledConcentration")
  
  fungal_plfa <- row_sum_any_of(df, unique(c(fungal_markers, amf_markers)))
  
  # Bacterial markers (branched + cyclo + monoenoic ω7)
  gram_pos <- c(
    "i14To0ScaledConcentration",
    "i15To0ScaledConcentration", "a15To0ScaledConcentration",
    "i16To0ScaledConcentration",
    "i17To0ScaledConcentration", "a17To0ScaledConcentration",
    "tenMe16To0ScaledConcentration", "tenMe17To0ScaledConcentration", "tenMe18To0ScaledConcentration"
  )
  gram_neg <- c(
    "c16To1n7ScaledConcentration", "cis16To1n7ScaledConcentration",
    "c18To1n7ScaledConcentration", "cis18To1n7ScaledConcentration",
    "cyclo17To0ScaledConcentration", "cyclo19To0ScaledConcentration"
  )
  bacterial_plfa <- row_sum_any_of(df, unique(c(gram_pos, gram_neg)))
  
  # F:B ratio
  fb <- fungal_plfa / bacterial_plfa
  fb[is.infinite(fb)] <- NA_real_
  
  # Stress index: (cy17 + cy19) / (16:1ω7 + 18:1ω7)
  cy <- row_sum_any_of(df, c("cyclo17To0ScaledConcentration", "cyclo19To0ScaledConcentration"))
  pre <- row_sum_any_of(df, c(
    "c16To1n7ScaledConcentration", "cis16To1n7ScaledConcentration",
    "c18To1n7ScaledConcentration", "cis18To1n7ScaledConcentration"
  ))
  stress <- cy / pre
  stress[is.infinite(stress)] <- NA_real_
  
  df$fungal_plfa_nmol_g <- fungal_plfa
  df$bacterial_plfa_nmol_g <- bacterial_plfa
  df$fb_ratio <- fb
  df$stress_index_cy_pre <- stress
  
  df
}

# --- Main wrapper: download + merge + compute ----------------
run_plfa_soil_health <- function(
    sites,
    start_date,
    end_date,
    release = "current",
    include_provisional = FALSE,
    package = "basic",
    data_dir = "neon_downloads",
    token = Sys.getenv("NEON_TOKEN", unset = NA_character_),
    quiet = TRUE
) {
  check_required_pkgs(c("dplyr"))
  
  dl <- download_neon_products(
    sites = sites,
    start_date = start_date,
    end_date = end_date,
    release = release,
    include_provisional = include_provisional,
    package = package,
    data_dir = data_dir,
    token = token,
    quiet = quiet
  )
  
  microb <- dl$microb
  soil   <- dl$soil
  
  if (is.null(microb$sme_scaledMicrobialBiomass)) {
    stop("DP1.10104.001: Could not find table 'sme_scaledMicrobialBiomass' in downloaded object.", call. = FALSE)
  }
  if (is.null(soil$sls_soilCoreCollection)) {
    stop("DP1.10086.001: Could not find table 'sls_soilCoreCollection' in downloaded object.", call. = FALSE)
  }
  
  mdf <- microb$sme_scaledMicrobialBiomass
  
  # NEON recommended correction for post-Nov 2021 (c18:0 contamination correction)
  # (Matches NEON tutorial logic)
  pre <- is.na(mdf$c18To0ScaledConcentration)
  mdf$microbial_biomass_nmol_g <- NA_real_
  mdf$microbial_biomass_nmol_g[pre] <- mdf$totalLipidScaledConcentration[pre]
  mdf$microbial_biomass_nmol_g[!pre] <- mdf$totalLipidScaledConcentration[!pre] - mdf$c18To0ScaledConcentration[!pre]
  
  # Merge with soil core metadata on sampleID
  core <- soil$sls_soilCoreCollection
  if (!("sampleID" %in% names(core)) || !("sampleID" %in% names(mdf))) {
    stop("Expected 'sampleID' in both soil core collection and microbial biomass tables.", call. = FALSE)
  }
  
  core_microb <- merge(core, mdf, by = "sampleID")
  
  # Identify a date column from the merged table
  date_col <- pick_first_existing(names(core_microb), c("collectDate", "collectDateTime", "startDateTime"))
  if (is.na(date_col)) {
    stop("Could not find a collection date column (collectDate/collectDateTime/startDateTime) after merging.", call. = FALSE)
  }
  core_microb$collectDate_std <- as_date_safe(core_microb[[date_col]])
  
  # Filter by site + exact date range
  core_microb <- core_microb |>
    dplyr::filter(.data$siteID %in% sites) |>
    dplyr::filter(.data$collectDate_std >= as_date_safe(start_date),
                  .data$collectDate_std <= as_date_safe(end_date))
  
  # QA filter if present
  if ("analysisResultsQF" %in% names(core_microb)) {
    core_microb <- core_microb |>
      dplyr::filter(grepl("OK", .data$analysisResultsQF))
  }
  
  # Compute indicator metrics
  core_microb <- compute_plfa_indicators(core_microb)
  
  # SOC extraction + join
  soc_tbl <- extract_soc(soil)
  if (!is.null(soc_tbl)) {
    # Find a common join field
    join_candidates <- c("sampleID", "cnSampleID", "sampleCode")
    join_in_core <- intersect(join_candidates, names(core_microb))
    join_id_col <- NA_character_
    
    # soc_tbl has "join_id" (whatever it was)
    # attempt to map: if soc was keyed by sampleID, it likely matches core$sampleID
    # If keyed by cnSampleID, core may also contain cnSampleID.
    if ("join_id" %in% names(soc_tbl)) {
      # pick a join candidate that exists in core and plausibly matches join_id
      # try sampleID first
      if ("sampleID" %in% names(core_microb) && grepl("^A000000", soc_tbl$join_id[1]) == FALSE) {
        join_id_col <- "sampleID"
      }
      # More robust: if core has cnSampleID and join_id looks like A000000...
      if (is.na(join_id_col) && "cnSampleID" %in% names(core_microb)) join_id_col <- "cnSampleID"
      if (is.na(join_id_col) && "sampleCode" %in% names(core_microb)) join_id_col <- "sampleCode"
      if (is.na(join_id_col)) join_id_col <- "sampleID"
    }
    
    # Build a temp key in core to join against soc_tbl$join_id
    core_microb$join_id <- core_microb[[join_id_col]]
    
    core_microb <- dplyr::left_join(core_microb, soc_tbl, by = "join_id")
    
    # Biomass / SOC (nmol PLFA per g C)
    core_microb$biomass_per_soc_nmol_gC <- core_microb$microbial_biomass_nmol_g / core_microb$soc_g_g
    core_microb$biomass_per_soc_nmol_gC[is.infinite(core_microb$biomass_per_soc_nmol_gC)] <- NA_real_
  } else {
    core_microb$soc_g_g <- NA_real_
    core_microb$biomass_per_soc_nmol_gC <- NA_real_
  }
  
  # Final tidy-ish output with standard names
  out <- core_microb |>
    dplyr::rename(
      collectDate = collectDate_std
    ) |>
    dplyr::select(
      siteID, plotID, sampleID, collectDate,
      # biomass
      microbial_biomass_nmol_g,
      # indicators
      fungal_plfa_nmol_g, bacterial_plfa_nmol_g, fb_ratio,
      stress_index_cy_pre,
      soc_g_g,
      biomass_per_soc_nmol_gC,
      # useful context if present
      dplyr::any_of(c("horizon", "horizonType", "nlcdClass", "plotType", "decimalLatitude", "decimalLongitude"))
    )
  
  list(
    data = out,
    download_months = list(start_ym = dl$start_ym, end_ym = dl$end_ym),
    release = release,
    include_provisional = include_provisional,
    package = package
  )
}
