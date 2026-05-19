# R/app.R
# NEON Soil Health Explorer
# Updated: 2026-04-01

suppressPackageStartupMessages({
  library(shiny)
  library(dplyr)
  library(ggplot2)
  library(DT)
  library(leaflet)
  library(stringr)
})

source("plfa_soil_health.R", local = TRUE)

require_pkgs(c(
  "neonUtilities","withr","jsonlite",
  "dplyr","tidyr","purrr","stringr","lubridate",
  "ggplot2","DT","leaflet","broom","tibble"
))

`%||%` <- function(a, b) if (!is.null(a)) a else b

DRIVER_VARS <- c(
  "soil temp" = "soilTemp",
  "soil moisture" = "soilMoisture",
  "microbial biomass" = "microbial_biomass_nmol_g",
  "year" = "year",
  "organic C %" = "organicCPercent",
  "N %" = "nitrogenPercent",
  "organic d13C" = "organicd13C",
  "d15N" = "d15N",
  "CN ratio" = "CN_ratio",
  "fungal PLFA" = "fungal_plfa_nmol_g",
  "bacterial PLFA" = "bacterial_plfa_nmol_g",
  "F:B ratio" = "FB_ratio",
  "stress index" = "stress_index"
)

METRICS_CORE <- c(
  "microbial_biomass_nmol_g",
  "fungal_plfa_nmol_g",
  "bacterial_plfa_nmol_g",
  "FB_ratio",
  "stress_index",
  "biomass_per_gC_nmol_gC",
  "organicCPercent",
  "nitrogenPercent",
  "CN_ratio"
)

coerce_to_char_vec <- function(x) {
  if (is.null(x)) return(character(0))
  if (is.list(x)) {
    return(vapply(x, function(z) {
      if (is.null(z) || length(z) == 0) return(NA_character_)
      as.character(z[[1]])
    }, character(1)))
  }
  as.character(x)
}

message_plot <- function(msg) {
  msg <- as.character(msg)
  if (length(msg) == 0 || is.na(msg[[1]])) msg <- "No plot to display."
  ggplot() +
    annotate("text", x = 0, y = 0, label = msg[[1]], hjust = 0, vjust = 1, size = 4) +
    xlim(-0.1, 1) + ylim(-0.1, 1) +
    theme_void()
}

horizon_class <- function(h) {
  h0 <- toupper(trimws(as.character(h)))
  dplyr::case_when(
    is.na(h0) ~ NA_character_,
    stringr::str_detect(h0, "^M") ~ "Mineral",
    stringr::str_detect(h0, "MIN") ~ "Mineral",
    stringr::str_detect(h0, "^O") ~ "Organic",
    stringr::str_detect(h0, "ORG") ~ "Organic",
    TRUE ~ h0
  )
}

lm_stats_tbl <- function(df) {
  df <- df %>% filter(is.finite(x) & is.finite(y))
  if (!is.data.frame(df) || nrow(df) < 3 || length(unique(df$x)) < 2) {
    return(tibble::tibble(
      n = nrow(df),
      slope = NA_real_,
      intercept = NA_real_,
      r2 = NA_real_,
      p_value = NA_real_
    ))
  }
  m <- stats::lm(y ~ x, data = df)
  s <- summary(m)
  tibble::tibble(
    n = nrow(df),
    slope = unname(coef(m)[["x"]]),
    intercept = unname(coef(m)[["(Intercept)"]]),
    r2 = unname(s$r.squared),
    p_value = unname(s$coefficients["x","Pr(>|t|)"])
  )
}

ui <- fluidPage(
  titlePanel("NEON Soil Health Explorer"),
  sidebarLayout(
    sidebarPanel(
      tags$p("Download NEON PLFA (DP1.10104.001) + soil properties (DP1.10086.001), compute soil health indicators, and explore patterns across sites, time, and drivers."),
      
      uiOutput("release_ui"),
      
      selectizeInput(
        "sites",
        "Sites to download/analyze (PLFA sites only):",
        choices = NULL,
        selected = NULL,
        multiple = TRUE,
        options = list(placeholder = "Select one or more sites…")
      ),
      
      dateRangeInput(
        "dates",
        "Date range (collection date filter):",
        start = as.Date("2015-01-01"),
        end = Sys.Date()
      ),
      
      actionButton("run", "Run analysis", class = "btn-primary"),
      
      tags$hr(),
      
      selectizeInput(
        "display_site",
        "Display site (plots/tabs):",
        choices = NULL,
        selected = NULL,
        multiple = FALSE,
        options = list(placeholder = "Run analysis first…")
      ),
      
      tags$hr(),
      
      downloadButton("download_data", "Download processed dataset (CSV)"),
      br(), br(),
      downloadButton("download_plots", "Download plots (ZIP)"),
      
      tags$hr(),
      verbatimTextOutput("status_text")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Map & Summary",
                 leafletOutput("site_map", height = 520),
                 tags$hr(),
                 uiOutput("soil_type_ui"),
                 tags$hr(),
                 DTOutput("summary_table")
        ),
        
        tabPanel("Key metrics (time series)",
                 plotOutput("ts_biomass", height = 260),
                 plotOutput("ts_fb", height = 260),
                 plotOutput("ts_stress", height = 260),
                 plotOutput("ts_biomassC", height = 260)
        ),
        
        tabPanel("Timepoint summaries",
                 tags$p("Boxplots summarize within-timepoint variability across sampled plots/subplots. Points removed for readability."),
                 plotOutput("tp_biomass", height = 320),
                 plotOutput("tp_stress", height = 320)
        ),
        
        tabPanel("Seasonality",
                 plotOutput("season_biomass", height = 300),
                 plotOutput("season_fb", height = 300),
                 plotOutput("season_stress", height = 300)
        ),
        
        tabPanel("Trends",
                 tags$p("Linear trend estimates require ≥ 3 unique years. Empty results often mean insufficient temporal coverage."),
                 DTOutput("trends_table")
        ),
        
        tabPanel("Drivers",
                 tags$p("Select a response + explanatory variable and horizon. A single large plot is shown with a linear fit and model statistics."),
                 fluidRow(
                   column(
                     4,
                     selectInput(
                       "drivers_response",
                       "Response metric:",
                       choices = as.list(DRIVER_VARS),
                       selected = "microbial_biomass_nmol_g"
                     )
                   ),
                   column(
                     4,
                     selectInput(
                       "drivers_explan",
                       "Explanatory variable:",
                       choices = as.list(DRIVER_VARS),
                       selected = "soilTemp"
                     )
                   ),
                   column(
                     4,
                     radioButtons(
                       "drivers_horizon",
                       "Horizon:",
                       choices = list("Mineral" = "Mineral", "Organic" = "Organic", "Both" = "Both"),
                       selected = "Mineral",
                       inline = TRUE
                     )
                   )
                 ),
                 tags$hr(),
                 uiOutput("drivers_soiltype_banner"),
                 plotOutput("drivers_grid", height = 520),
                 DTOutput("drivers_fit_table")
        ),
        
        tabPanel("Diagnostics",
                 DTOutput("diagnostics_table")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  rv <- reactiveValues(
    site_table = NULL,
    result = NULL,
    diagnostics = NULL,
    status = "Ready."
  )
  
  output$release_ui <- renderUI({
    rel <- tryCatch(get_neon_release_choices(), error = function(e) "current")
    if (length(rel) == 0) rel <- "current"
    selectInput("release", "NEON release:", choices = rel, selected = default_release(rel))
  })
  
  # --- Load site table on startup ---
  observe({
    rv$status <- "Loading NEON site metadata for PLFA sites…"
    output$status_text <- renderText(rv$status)
    
    tbl <- tryCatch(get_plfa_site_table(cache_dir = "cache", force_refresh = FALSE), error = function(e) NULL)
    if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0) {
      tbl <- tryCatch(get_plfa_site_table(cache_dir = "cache", force_refresh = TRUE), error = function(e) NULL)
    }
    
    rv$site_table <- tbl
    
    if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0) {
      rv$status <- "No NEON PLFA sites loaded (metadata issue). Check internet access, or install the 'curl' R package."
      updateSelectizeInput(session, "sites", choices = character(0), selected = character(0), server = TRUE)
      output$status_text <- renderText(rv$status)
      output$site_map <- renderLeaflet({
        leaflet() %>% addProviderTiles(providers$CartoDB.Positron)
      })
      return()
    }
    
    # defensive coercion
    tbl$siteID <- coerce_to_char_vec(tbl$siteID)
    tbl$label <- coerce_to_char_vec(tbl$label)
    
    tbl <- tbl %>%
      filter(!is.na(siteID) & str_detect(siteID, "^[A-Z]{4}$")) %>%
      distinct(siteID, .keep_all = TRUE)
    
    # build choices safely
    vals <- as.character(tbl$siteID)
    labs <- as.character(tbl$label)
    n <- min(length(vals), length(labs))
    vals <- vals[seq_len(n)]
    labs <- labs[seq_len(n)]
    choices <- stats::setNames(vals, labs)
    
    updateSelectizeInput(session, "sites", choices = choices, selected = character(0), server = TRUE)
    
    rv$status <- paste0("Loaded ", nrow(tbl), " PLFA sites. Select sites and click Run.")
    output$status_text <- renderText(rv$status)
    
    output$site_map <- renderLeaflet({
      pts <- tbl %>% filter(!is.na(latitude) & !is.na(longitude))
      
      m <- leaflet() %>% addProviderTiles(providers$CartoDB.Positron)
      
      if (nrow(pts) > 0) {
        m <- m %>%
          addCircleMarkers(
            data = pts,
            lng = ~longitude, lat = ~latitude,
            radius = 5,
            popup = ~paste0("<b>", siteID, "</b><br/>", siteName, "<br/>", domainID),
            label = ~siteID,
            stroke = FALSE, fillOpacity = 0.8
          ) %>%
          fitBounds(
            lng1 = min(pts$longitude, na.rm = TRUE),
            lat1 = min(pts$latitude, na.rm = TRUE),
            lng2 = max(pts$longitude, na.rm = TRUE),
            lat2 = max(pts$latitude, na.rm = TRUE)
          )
      }
      m
    })
  })
  
  # --- Run analysis ---
  observeEvent(input$run, {
    if (is.null(input$release) || !nzchar(input$release)) {
      rv$status <- "ERROR: Please select a NEON release."
      output$status_text <- renderText(rv$status)
      showNotification(rv$status, type = "error", duration = 6)
      return()
    }
    if (is.null(input$sites) || length(input$sites) == 0) {
      rv$status <- "ERROR: Please select at least one site."
      output$status_text <- renderText(rv$status)
      showNotification(rv$status, type = "error", duration = 6)
      return()
    }
    
    start_date <- as.character(input$dates[1])
    end_date <- as.character(input$dates[2])
    
    rv$status <- "Starting…"
    output$status_text <- renderText(rv$status)
    
    withProgress(message = "Running NEON Soil Health Explorer…", value = 0, {
      
      incProgress(0.05, detail = "Preparing request…")
      
      sites <- input$sites
      rel <- input$release
      
      incProgress(0.10, detail = paste0("Downloading + processing ", length(sites), " site(s)…"))
      
      all_results <- list()
      all_diags <- list()
      
      for (i in seq_along(sites)) {
        s <- sites[[i]]
        incProgress(0.80 / length(sites), detail = paste0("Site ", s, " (", i, "/", length(sites), "): downloading + merging…"))
        
        out <- tryCatch(
          run_neon_soil_health(
            sites = c(s),
            release = rel,
            start_date = start_date,
            end_date = end_date,
            cache_dir = "cache/neon_downloads"
          ),
          error = function(e) e
        )
        
        if (inherits(out, "error")) {
          rv$status <- paste0("ERROR at site ", s, ": ", conditionMessage(out))
          output$status_text <- renderText(rv$status)
          showNotification(rv$status, type = "error", duration = 10)
          rv$result <- NULL
          rv$diagnostics <- NULL
          return()
        }
        
        all_results[[s]] <- out$data
        all_diags[[s]] <- out$diagnostics
      }
      
      incProgress(0.05, detail = "Finalizing outputs…")
      combined <- bind_rows(all_results)
      
      if (is.null(combined) || !is.data.frame(combined) || nrow(combined) == 0) {
        rv$status <- "ERROR: Analysis completed but returned 0 rows (no data after filtering). Try expanding the date range."
        output$status_text <- renderText(rv$status)
        showNotification(rv$status, type = "error", duration = 10)
        rv$result <- NULL
        rv$diagnostics <- all_diags
        return()
      }
      
      rv$result <- combined
      rv$diagnostics <- all_diags
      
      disp_choices <- unique(na.omit(combined$siteID))
      disp_choices <- disp_choices[order(disp_choices)]
      
      updateSelectizeInput(
        session, "display_site",
        choices = disp_choices,
        selected = disp_choices[[1]],
        server = TRUE
      )
      
      rv$status <- paste0("Done. Processed ", nrow(combined), " rows across ", length(unique(combined$siteID)), " sites.")
      output$status_text <- renderText(rv$status)
      
      incProgress(1, detail = "Complete.")
    })
  })
  
  display_df <- reactive({
    req(rv$result)
    req(input$display_site)
    rv$result %>% filter(siteID == input$display_site)
  })
  
  # --- Map zoom to display site ---
  observeEvent(input$display_site, {
    req(rv$site_table)
    req(input$display_site)
    
    st <- rv$site_table %>% filter(siteID == input$display_site) %>% slice(1)
    if (nrow(st) == 0) return()
    if (is.na(st$latitude[[1]]) || is.na(st$longitude[[1]])) return()
    
    leafletProxy("site_map") %>%
      clearGroup("selected") %>%
      addCircleMarkers(
        data = st,
        lng = ~longitude, lat = ~latitude,
        radius = 10, color = "black", weight = 2,
        fillOpacity = 0.9,
        group = "selected",
        popup = ~paste0("<b>", siteID, "</b><br/>", siteName, "<br/>", domainID)
      ) %>%
      setView(lng = st$longitude[[1]], lat = st$latitude[[1]], zoom = 9)
  })
  
  output$soil_type_ui <- renderUI({
    req(display_df())
    d <- display_df()
    soilType <- unique(na.omit(d$soilType))
    if (length(soilType) == 0) {
      tags$p(tags$strong("Soil type (Megapit): "), "Not available for this site (or megapit pull failed).")
    } else {
      tags$p(tags$strong("Soil type (Megapit): "), soilType[[1]])
    }
  })
  
  output$drivers_soiltype_banner <- renderUI({
    req(display_df())
    d <- display_df()
    soilType <- unique(na.omit(d$soilType))
    if (length(soilType) == 0) {
      tags$div(style="padding:8px; background:#f7f7f7; border:1px solid #ddd;",
               tags$strong("Soil type (Megapit): "), "Not available.")
    } else {
      tags$div(style="padding:8px; background:#f7f7f7; border:1px solid #ddd;",
               tags$strong("Soil type (Megapit): "), soilType[[1]])
    }
  })
  
  output$summary_table <- renderDT({
    req(display_df())
    d <- display_df()
    
    summ <- d %>%
      summarize(
        siteID = first(siteID),
        n_samples = n(),
        n_years = n_distinct(year, na.rm = TRUE),
        horizons = paste(sort(unique(na.omit(horizon))), collapse = ", "),
        biomass_median = median(microbial_biomass_nmol_g, na.rm = TRUE),
        FB_median = median(FB_ratio, na.rm = TRUE),
        stress_median = median(stress_index, na.rm = TRUE),
        organicC_median = median(organicCPercent, na.rm = TRUE),
        N_median = median(nitrogenPercent, na.rm = TRUE)
      )
    
    datatable(summ, options = list(dom = "t"), rownames = FALSE)
  })
  
  ts_plot <- function(df, y, ylab) {
    ggplot(df, aes(x = collectDate, y = .data[[y]])) +
      geom_line(na.rm = TRUE) +
      geom_point(size = 1.4, alpha = 0.6, na.rm = TRUE) +
      facet_wrap(~horizon, scales = "free_y") +
      labs(x = "Collection date", y = ylab) +
      theme_bw()
  }
  
  season_plot <- function(df, y, ylab) {
    ggplot(df, aes(x = factor(month), y = .data[[y]])) +
      geom_boxplot(na.rm = TRUE, outlier.alpha = 0.25) +
      facet_wrap(~horizon, scales = "free_y") +
      labs(x = "Month", y = ylab) +
      theme_bw()
  }
  
  timepoint_box <- function(df, y, ylab) {
    df2 <- df %>% mutate(tp = as.factor(collectDate))
    ggplot(df2, aes(x = tp, y = .data[[y]])) +
      geom_boxplot(na.rm = TRUE, outlier.alpha = 0.25) +
      facet_wrap(~horizon, scales = "free_y") +
      labs(x = "Timepoint (date)", y = ylab) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  output$ts_biomass <- renderPlot({
    req(display_df())
    ts_plot(display_df(), "microbial_biomass_nmol_g", "Microbial biomass (nmol PLFA / g soil; c18:0-corrected)")
  })
  output$ts_fb <- renderPlot({
    req(display_df())
    ts_plot(display_df(), "FB_ratio", "F:B ratio (fungal / bacterial)")
  })
  output$ts_stress <- renderPlot({
    req(display_df())
    ts_plot(display_df(), "stress_index", "Stress index (cyclopropyl / precursor)")
  })
  output$ts_biomassC <- renderPlot({
    req(display_df())
    ts_plot(display_df(), "biomass_per_gC_nmol_gC", "Biomass per g C (nmol PLFA / g C)")
  })
  
  output$tp_biomass <- renderPlot({
    req(display_df())
    timepoint_box(display_df(), "microbial_biomass_nmol_g", "Microbial biomass (nmol/g; c18:0-corrected)")
  })
  output$tp_stress <- renderPlot({
    req(display_df())
    timepoint_box(display_df(), "stress_index", "Stress index")
  })
  
  output$season_biomass <- renderPlot({
    req(display_df())
    season_plot(display_df(), "microbial_biomass_nmol_g", "Microbial biomass (nmol/g; c18:0-corrected)")
  })
  output$season_fb <- renderPlot({
    req(display_df())
    season_plot(display_df(), "FB_ratio", "F:B ratio")
  })
  output$season_stress <- renderPlot({
    req(display_df())
    season_plot(display_df(), "stress_index", "Stress index")
  })
  
  output$trends_table <- renderDT({
    req(display_df())
    d <- display_df()
    tr <- compute_trends(d, metrics = METRICS_CORE, group_vars = c("siteID", "horizon"))
    datatable(tr, options = list(pageLength = 25), rownames = FALSE)
  })
  

  build_drivers_outputs <- reactive({
    d <- tryCatch(display_df(), error = function(e) NULL)
    if (is.null(d) || !is.data.frame(d) || nrow(d) == 0) {
      return(list(
        plot = message_plot("Run analysis and select a display site to view Drivers."),
        stats = tibble::tibble()
      ))
    }
    
    resp <- input$drivers_response
    expl <- input$drivers_explan
    hz_sel <- input$drivers_horizon
    
    if (is.null(resp) || !nzchar(resp) || !(resp %in% names(d))) {
      return(list(plot = message_plot("Selected response is not available for this site."), stats = tibble::tibble()))
    }
    if (is.null(expl) || !nzchar(expl) || !(expl %in% names(d))) {
      return(list(plot = message_plot("Selected explanatory variable is not available for this site."), stats = tibble::tibble()))
    }
    
    d <- d %>% mutate(horizonClass = horizon_class(horizon))
    
    if (hz_sel %in% c("Mineral","Organic")) {
      d <- d %>% filter(horizonClass == hz_sel)
    } # Both = no filter
    
    # Numeric pairing
    df <- tibble::tibble(
      horizonClass = d$horizonClass,
      x = suppressWarnings(as.numeric(d[[expl]])),
      y = suppressWarnings(as.numeric(d[[resp]]))
    ) %>%
      filter(is.finite(x) & is.finite(y))
    
    if (!is.data.frame(df) || nrow(df) < 3) {
      return(list(
        plot = message_plot("Not enough non-missing paired data for this site / horizon (try widening date range)."),
        stats = tibble::tibble()
      ))
    }
    
    # Labels
    nice_name <- function(var) {
      nm <- names(DRIVER_VARS)[match(var, DRIVER_VARS)]
      if (is.na(nm) || !nzchar(nm)) var else nm
    }
    xlab <- nice_name(expl)
    ylab <- nice_name(resp)
    
    # Identity case
    if (identical(resp, expl)) {
      p <- ggplot(df, aes(x = x, y = y)) +
        geom_point(alpha = 0.6, na.rm = TRUE) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        labs(x = xlab, y = ylab,
             title = "Identity plot (response = explanatory)") +
        theme_bw()
      
      st <- tibble::tibble(
        horizon = hz_sel,
        n = nrow(df),
        slope = NA_real_,
        intercept = NA_real_,
        r2 = NA_real_,
        p_value = NA_real_,
        note = "Response and explanatory are identical; showing 1:1 line (no regression)."
      )
      
      return(list(plot = p, stats = st))
    }
    
    # Fit models
    if (hz_sel == "Both") {
      # per-horizon stats + overall
      st_h <- df %>%
        group_by(horizonClass) %>%
        group_modify(~lm_stats_tbl(.x)) %>%
        ungroup() %>%
        rename(horizon = horizonClass)
      
      st_all <- lm_stats_tbl(df) %>%
        mutate(horizon = "Overall")
      
      st <- bind_rows(st_all, st_h) %>%
        mutate(
          equation = ifelse(is.na(slope), NA_character_,
                            paste0("y = ", signif(intercept, 4), " + ", signif(slope, 4), "x")),
          p_value = as.numeric(p_value),
          r2 = as.numeric(r2)
        )
      
      # Overall annotation
      ann <- st %>% filter(horizon == "Overall") %>% slice(1)
      ann_txt <- paste0(
        "Overall fit:\n",
        "n=", ann$n, "\n",
        "slope=", signif(ann$slope, 4), "\n",
        "R²=", signif(ann$r2, 4), "\n",
        "p=", signif(ann$p_value, 4)
      )
      
      xpos <- unname(stats::quantile(df$x, 0.05, na.rm = TRUE))
      ypos <- unname(stats::quantile(df$y, 0.95, na.rm = TRUE))
      
      p <- ggplot(df, aes(x = x, y = y, color = horizonClass)) +
        geom_point(alpha = 0.6, na.rm = TRUE) +
        geom_smooth(method = "lm", se = FALSE, na.rm = TRUE) +
        annotate("text", x = xpos, y = ypos, label = ann_txt, hjust = 0, vjust = 1, size = 4) +
        labs(x = xlab, y = ylab, color = "Horizon") +
        theme_bw()
      
      return(list(plot = p, stats = st))
    } else {
      st <- lm_stats_tbl(df) %>%
        mutate(
          horizon = hz_sel,
          equation = ifelse(is.na(slope), NA_character_,
                            paste0("y = ", signif(intercept, 4), " + ", signif(slope, 4), "x"))
        )
      
      ann_txt <- paste0(
        "Fit (", hz_sel, "):\n",
        "n=", st$n, "\n",
        "slope=", signif(st$slope, 4), "\n",
        "R²=", signif(st$r2, 4), "\n",
        "p=", signif(st$p_value, 4)
      )
      
      xpos <- unname(stats::quantile(df$x, 0.05, na.rm = TRUE))
      ypos <- unname(stats::quantile(df$y, 0.95, na.rm = TRUE))
      
      p <- ggplot(df, aes(x = x, y = y)) +
        geom_point(alpha = 0.6, na.rm = TRUE) +
        geom_smooth(method = "lm", se = FALSE, na.rm = TRUE) +
        annotate("text", x = xpos, y = ypos, label = ann_txt, hjust = 0, vjust = 1, size = 4) +
        labs(x = xlab, y = ylab) +
        theme_bw()
      
      return(list(plot = p, stats = st))
    }
  })
  
  output$drivers_grid <- renderPlot({
    out <- build_drivers_outputs()
    out$plot
  })
  
  output$drivers_fit_table <- renderDT({
    out <- build_drivers_outputs()
    st <- out$stats
    if (is.null(st) || !is.data.frame(st) || nrow(st) == 0) {
      return(datatable(tibble::tibble(message = "No model statistics to display."), options = list(dom = "t"), rownames = FALSE))
    }
    datatable(
      st,
      options = list(dom = "t"),
      rownames = FALSE
    )
  })
  
  output$diagnostics_table <- renderDT({
    req(rv$diagnostics)
    dd <- rv$diagnostics
    
    rows <- lapply(names(dd), function(site) {
      x <- dd[[site]]
      tibble::tibble(
        site = site,
        release_requested = as.character(x$release_requested %||% NA),
        n_rows = as.numeric(x$n_rows %||% NA),
        FB_all_na = as.logical(x$FB_all_na %||% NA),
        stress_all_na = as.logical(x$stress_all_na %||% NA),
        soc_all_na = as.logical(x$soc_all_na %||% NA),
        fungal_col = as.character(x$marker_cols_used$fungal_col %||% NA),
        n_bacterial_cols = length(x$marker_cols_used$bacterial_cols %||% character(0)),
        cyclo_cols = paste(x$marker_cols_used$cyclo_cols %||% character(0), collapse = ", "),
        precursor_cols = paste(x$marker_cols_used$precursor_cols %||% character(0), collapse = ", ")
      )
    })
    
    datatable(bind_rows(rows), options = list(pageLength = 25), rownames = FALSE)
  })
  
  output$download_data <- downloadHandler(
    filename = function() paste0("neon_soil_health_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$result)
      write.csv(rv$result, file, row.names = FALSE)
    }
  )
  
  output$download_plots <- downloadHandler(
    filename = function() paste0("neon_soil_health_plots_", input$display_site, "_", Sys.Date(), ".zip"),
    content = function(file) {
      req(display_df())
      d <- display_df()
      
      tmpdir <- tempfile("plots_")
      dir.create(tmpdir)
      
      plots <- list(
        ts_biomass = ts_plot(d, "microbial_biomass_nmol_g", "Microbial biomass (nmol/g; c18:0-corrected)"),
        ts_fb = ts_plot(d, "FB_ratio", "F:B ratio"),
        ts_stress = ts_plot(d, "stress_index", "Stress index"),
        ts_biomassC = ts_plot(d, "biomass_per_gC_nmol_gC", "Biomass per g C (nmol/g C)"),
        season_biomass = season_plot(d, "microbial_biomass_nmol_g", "Microbial biomass"),
        season_fb = season_plot(d, "FB_ratio", "F:B ratio"),
        season_stress = season_plot(d, "stress_index", "Stress index")
      )
      
      for (nm in names(plots)) {
        ggsave(
          filename = file.path(tmpdir, paste0(nm, ".png")),
          plot = plots[[nm]],
          width = 10, height = 6, dpi = 300
        )
      }
      
      oldwd <- getwd()
      setwd(tmpdir)
      on.exit(setwd(oldwd), add = TRUE)
      files <- list.files(tmpdir, pattern = "\\.png$", full.names = FALSE)
      utils::zip(zipfile = file, files = files)
    }
  )
}

shinyApp(ui, server)