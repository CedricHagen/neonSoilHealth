# NEON Soil Health Explorer - User Guide

A Shiny application and R workflow to download, process, and visualize NEON soil health indicators across sites and time, with a focus on microbial community lipid biomarkers (PLFA) and associated soil chemistry, temperature, moisture, and soil type.

## Purpose

This tool is designed to support:

- **Cross-site comparisons** of microbial biomass and microbial community indicators
- **Temporal patterns** (seasonality, sampling-event summaries, and trends)
- **Driver screening** via simple linear regressions (response vs explanatory variables)
- **Exportable, publication-friendly outputs** (processed dataset + plot ZIP)

---

## What This Tool Does

### Inputs

From the UI you choose:
- NEON Data Release (e.g., RELEASE-2026)
- NEON sites (only sites with PLFA data are offered)
- Date range
- Click **Run** (backend downloads + processing start)

### Outputs (Core Metrics)

Computed per sample (and summarized in plots/tables):

- **Total microbial biomass (corrected)** (nmol lipids / g soil)
- **Fungal PLFA** (nmol/g)
- **Bacterial PLFA** (nmol/g)
- **F:B ratio** (fungal:bacterial ratio)
- **Stress index** (cyclopropyl : precursor ratio)
- **Microbial biomass / SOC** (biomass normalized by organic C fraction)
- **Soil chemistry**: organic C %, N %, CN ratio, organic δ¹³C, δ¹⁵N
- **Soil temperature and soil moisture** at sampling (when available)
- **Soil type / classification** from NEON Megapit pedon descriptions (static per site)

### Visualizations

- **Time series** of metrics by date (site-level and/or horizon subsets)
- **Seasonality** (month / bout-style summaries where possible)
- **Timepoint summary** (boxplots per sampling event)
- **Trends** (per-year slope, p-value, and R² where enough years exist)
- **Drivers tab**: one large regression plot at a time with:
  - Response selector (limited list)
  - Explanatory selector (limited list)
  - Horizon selector (Mineral vs Organic)
  - Regression line + R² + p-value
  - 1:1 line if response == explanatory
- **Map** of NEON sites (with automatic zoom/highlight to selected sites)

### Exports

- Download the processed dataset (CSV)
- Export all plots as a ZIP (PNGs)

---

## Repository Structure

Typical layout:

```
neonSoilHealth/
  README.md
  R/
    app.R
    plfa_soil_health.R
  cache/
    neon_downloads/        # raw NEON downloads (cached)
    processed/             # processed datasets (cached)
    plots/                 # temporary plot files for ZIP export
```

---

## Requirements

### R Packages

The app uses (at minimum):

- `shiny`, `shinyWidgets`
- `dplyr`, `tidyr`, `stringr`, `purrr`, `readr`
- `ggplot2`
- `lubridate`
- `leaflet`
- `broom` (for regression summaries)
- `neonUtilities` (for NEON downloads)

If you are deploying (shinyapps.io / Posit Connect), also ensure:

- `rmarkdown` (optional, only if you later add report generation)
- `zip` (or `utils::zip()` availability on server)

---

## Running the App Locally

From the project root in R:

```r
setwd("path/to/neonSoilHealth")

# Install packages if needed
# install.packages(c("shiny","dplyr","tidyr","ggplot2","leaflet","broom","neonUtilities", ...))

shiny::runApp("R")
```

---

## NEON Data Products Used

This tool pulls from NEON released data products via `neonUtilities::loadByProduct()`:

### Soil microbe biomass (PLFA)
- **DP1.10104.001**

### Soil physical & chemical properties, periodic
- **DP1.10086.001**
  - Used for:
    - Soil core collection metadata (incl. collection timing and horizon info)
    - Soil moisture at sampling (where provided)
    - `sls_soilChemistry` table (organicCPercent, nitrogenPercent, organicd13C, d15N, etc.)

### Soil type / pedon classification (Megapit)
- The app uses NEON Megapit pedon descriptions to attach a site-level soil classification/soil type string
- This is treated as static per site and repeated across rows for that site in the output dataset

**Note:** Some sites may not have Megapit descriptors available. The tool will keep `soil_type` as NA in those cases.

---

## Site List Behavior

- The site selector is populated with terrestrial NEON sites that have PLFA data available in the chosen release
- The list is sorted by domain (D01, D02, …) then by site code
- The app does not auto-select any site at startup (you choose)

---

## Processing and QA/QC

### 1. Join Keys

Most tables are merged using NEON sample identifiers:
- The PLFA table is merged to soil tables using `sampleID` (or the appropriate NEON linking key present in the downloaded tables)

### 2. Collection Date Handling

The pipeline standardizes a single column named:
- `collectDate` (Date)

If the upstream table uses a different field name (e.g., `collectDateTime` or `startDateTime`), the pipeline converts it and stores the standardized `collectDate`.

### 3. Filtering

The pipeline uses NEON QA/QC flags where available. Typical behavior:
- Keep rows that pass lab QA/QC (e.g., `analysisResultsQF` indicates acceptable status)
- Drop rows missing essential identifiers (`siteID`, `sampleID`, `collectDate`)

Diagnostics in the app report how many rows were retained and which marker columns were used.

---

## Metric Calculations

### A) Total Microbial Biomass (Corrected; Excludes 18:0)

**Goal:** Prevent known contamination issues from inflating biomass.

NEON's workflow notes that C18:0 (stearic acid; "c18To0") contamination can occur due to the PLFA method materials. When `c18To0ScaledConcentration` is reported, it should be subtracted from the total lipid signal to obtain a reliable biomass estimate.

This app implements the same logic:

- If `c18To0ScaledConcentration` is missing/NA (older data periods):
  - `microbial_biomass_total = totalLipidScaledConcentration`
- If `c18To0ScaledConcentration` is present:
  - `microbial_biomass_total = totalLipidScaledConcentration - c18To0ScaledConcentration`

**Units:** nmol lipids / g dry soil (as reported by NEON for lipid scaled concentrations)

**Output column:** `microbial_biomass_total_nmol_g` (this is the corrected total)

### B) Fungal PLFA Biomass

Fungal biomass is estimated from a fungal biomarker PLFA (most commonly 18:2ω6,9; naming varies by dataset release).

**Output:** `fungal_plfa_nmol_g`

If the fungal marker column is not present for a given release/site subset, the tool reports this in Diagnostics and `fungal_plfa_nmol_g` may be NA.

### C) Bacterial PLFA Biomass

Bacterial biomass is computed as the sum of bacterial biomarker PLFAs (a conservative set that is consistently reported in NEON PLFA outputs).

Because NEON column availability can differ by lab era/release, the tool:
- Dynamically detects which bacterial marker columns exist
- Sums those that are present

**Output:** `bacterial_plfa_nmol_g`

Diagnostics lists the exact bacterial marker columns used for your run.

### D) Fungal:Bacterial Ratio (F:B)

Computed as:

```
FB_ratio = fungal_plfa_nmol_g / bacterial_plfa_nmol_g
```

**Rules:**
- If either fungal or bacterial biomass is missing, the ratio is NA
- If bacterial biomass is zero, ratio is NA (avoid divide-by-zero artifacts)

**Output:** `FB_ratio`

### E) Stress Index

The app implements a cyclopropyl:precursor stress index, conceptually tracking membrane shifts often associated with nutrient/physiological stress.

Computed as:

```
stress_index = (cy17:0 + cy19:0) / (16:1ω7c + 18:1ω7c)
```

Where available. As with other PLFA-based metrics:
- Column availability varies across releases/sites
- If required markers are missing, stress is NA
- Diagnostics reports which stress-related columns were found and used

**Output:** `stress_index`

### F) Soil Organic C, N, Isotopes, and CN Ratio

Pulled from DP1.10086.001 → `sls_soilChemistry`:

- `organicCPercent`
- `nitrogenPercent`
- `organicd13C`
- `d15N`

**Additional derived column:**

```
CN_ratio = organicCPercent / nitrogenPercent
```

**Outputs:**
- `organicCPercent`
- `nitrogenPercent`
- `organicd13C`
- `d15N`
- `CN_ratio`

**Note:** No interpolation is performed between sampling events/years. If a chemistry value is missing for a sample date, it remains NA.

### G) Biomass / SOC

This is computed only when `organicCPercent` is present.

Convert C percent to a fraction:

```
C_frac = organicCPercent / 100
```

Then:

```
biomass_per_SOC = microbial_biomass_total_nmol_g / C_frac
```

**Output:** `biomass_per_SOC_nmol_per_gC`

**Notes:**
- If `organicCPercent` is NA or 0, this value is NA

---

## Trends Tab (Slopes, p-values, R²)

Trends are computed using simple linear regression:

1. Aggregate to annual values by site × horizon (e.g., mean or median per year)
2. Fit: `metric ~ year`

The app only reports slopes/p-values/R² when there are enough data to support it. If a site has only 1–2 distinct years, trend stats will remain blank/NA.

**Outputs:**
- `slope_per_year`
- `p_value`
- `r2`
- Plus metadata (siteID, horizon, metric name)

---

## Drivers Tab (Interactive Regression Plot)

The Drivers tab is designed for exploratory screening.

### Controls
- **Horizon toggle:** Mineral vs Organic
- **Response variable selector**
- **Explanatory variable selector**

### Allowed Variables (Intentionally Limited)

The app only exposes these options:
- Soil temperature
- Soil moisture
- Microbial biomass
- Year
- Organic C %
- N %
- Organic δ¹³C
- δ¹⁵N
- CN ratio
- Fungal PLFA
- Bacterial PLFA
- F:B ratio
- Stress index

### Outputs

- One large plot at a time (readable)
- Linear fit line (or 1:1 line if response == explanatory)
- R² and p-value shown on the plot panel

---

## Soil Type (Megapit) Display

The app pulls a site-level soil classification/soil type string from NEON Megapit pedon descriptions and displays it at the top of the Drivers tab.

In the processed dataset, the column `soil_type` is repeated across all rows for a site (because soil type does not change at the timescales of periodic sampling).

---

## Caching and Performance

To make repeated use fast and API-friendly:

- Downloads are cached under `cache/neon_downloads/`
- Processed merged results can be cached under `cache/processed/`

If you suspect a corrupted download or want to force re-download:
- Delete the relevant cache subfolders (or use the app's "Clear cache" control if enabled)

---

## Diagnostics

The Diagnostics panel is intended to make the pipeline transparent and publishable:

- Release requested
- Sites requested
- Number of rows downloaded / processed
- Whether key outputs are all NA (e.g., FB_ratio)
- Exact PLFA marker columns used for fungal, bacterial, and stress metrics
- SOC/chemistry availability stats

This is especially useful because NEON PLFA column presence can vary across analytical eras.

---

## Data Licensing / Access Notes

Historically, NEON data products have been distributed under CC0 (public domain dedication). NEON has also announced a transition to required login/API tokens and an updated data license (CC BY 4.0) expected in Summer 2026.

If you encounter authentication errors, you may need to configure a NEON API token for neonUtilities downloads.

---

## Recommended Citation

When publishing outputs derived from this tool:

- Cite NEON for each data product used (DP1.10104.001, DP1.10086.001, and Megapit resources if used)
- Include the release name (e.g., RELEASE-2026) and date accessed
- Cite relevant literature supporting PLFA interpretation and stress indices where appropriate

---

## Contact / Contributions

If you extend this tool (new metrics, mixed-effects models, advanced drivers, etc.), keep Diagnostics updated so that the pipeline remains transparent and reproducible.

For questions or contributions, please see the main [README.md](../README.md) file.

---

## Reference Notes

These are the main references that informed this guide:

- NEON tutorial describing C18:0 contamination and the correction: subtract `c18To0ScaledConcentration` from `totalLipidScaledConcentration` where applicable
- NEON tutorial documenting use of DP1.10104.001 and DP1.10086.001 together, and showing merge on sampleID
- NEON User Guide for Soil Microbe Biomass (DP1.10104.001): PLFA as a proxy for viable biomass; mentions applications including fungal:bacterial ratios and stress ratios
- NEON metagenomics tutorial (DP1.10086.001) confirming `sls_soilChemistry` variables: organicCPercent, nitrogenPercent, organicd13C, d15N, CNratio
- NEON Megapit pedon "Megapit Details" pages contain soil classification/type strings (Taxonomy/classification, map unit, etc.) that can be attached as site-level soil type
- NEON Data Usage Policy (CC0) and NEON announcement of planned licensing/authentication changes (Summer 2026)
