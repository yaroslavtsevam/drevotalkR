# drevotalkR

**`drevotalkR`** is a fast, lightweight R package for processing raw telemetry and SD-card backups from **TreeTalker** devices (`TT+` v2/v3, `TT-Soil`, `TT-Reference` / `TTR`, `TTG`, `TT-Cloud`), transforming biophysical signals into **ecosystem service indicators**, and generating publication-ready `ggplot2` visualizations with switchable **English** (`lang = "en"`, default) and **Russian** (`lang = "ru"`) labels.

## Key Design Features

- **Minimal dependency footprint**: requires only `data.table` and `ggplot2` (plus base R packages).
- **High speed + standard outputs**: all tabular, grouped, and rolling-window computations are powered internally by `data.table`, while **every user-facing function returns a pure `base::data.frame`**.
- **Self-contained scientific algorithms** (with explicit source attribution in code & help pages):
  - Savitzky–Golay polynomial smoothing (`sg_filter`, adapted from `prospectr::savitzkyGolay` / Savitzky & Golay, 1964)
  - Solar altitude & azimuth calculation (`get_sun_position`, adapted from `suncalc` / Meeus, 1998)
  - Thermal dissipation zero-flow $\Delta T_{\max}$ estimation (`calc_dTmax` methods `"default"`, `"pd"`, `"mw"`, `"dr"`, adapted from `TREX::tdm_dt.max`, Peters et al., 2021)
  - Hampel + MAD radial growth spike filter (`clean_growth_series`, adapted from `pracma::hampel` & `forecast::tsclean`)
  - Bounded small-gap interpolation (`fill_small_gaps`, adapted from `baytrends::fillMissing` & `zoo::na.approx`)

## Installation from GitHub

```r
# install.packages("remotes")
remotes::install_github("yaroslavtsevam/drevotalkR")
```

Or install locally from the repository folder:

```r
remotes::install_local("drevotalkR")
```

## Quick Start Workflow

```r
library(drevotalkR)

# 1. Scrape from Cloud servers or import local SD-card backups
raw_list <- tt_scrape(
  id = c("21030001"),
  import_folder = "path/to/backup_folder"
)

# 2. Run full biophysical processing (returns a pure data.frame)
bio_df <- tt_process(
  raw_list,
  step_time = 3600,
  grow_calc = TRUE,
  extrap_dates = TRUE,
  dtmax_method = "mw" # "default", "pd", "mw", or "dr"
)

# 3. Transform into Ecosystem Services (cooling, transpiration, carbon, stability)
es_df <- tt_calc_ecosystem_services(
  bio_df,
  ref_species = "TTR",
  wood_density = 0.50,
  carbon_fraction = 0.50
)

# Daily or seasonal summary table (pure data.frame)
es_daily  <- tt_summarise_services(es_df, period = "day")
es_season <- tt_summarise_services(es_df, period = "season")

# 4. Generate standard plots (English by default, switchable to Russian via lang = "ru")
tt_plot_timeseries(es_df, var = "Fdh", group_by = "Species", lang = "en")
tt_plot_diurnal(es_df, var = "ES_dT", group_by = "Species", lang = "ru")
tt_plot_ecoservices(es_df, service = "cooling", hot_only = TRUE, lang = "en")
tt_plot_waterbalance(es_df, group_by = "Species", unit = "L", lang = "en")
tt_plot_scatter(es_df, x = "VPD", y = "Fdh", group_by = "Species", lang = "en")
```

## License

GPL (>= 3)

