#' Transform biophysical TreeTalker series into Ecosystem Service indicators
#'
#' Calculates hourly ecosystem service metrics from processed TreeTalker data:
#' \enumerate{
#'   \item \bold{Microclimate regulation (urban cooling & humidifying)}:
#'     canopy temperature anomaly (`ES_dT`, \eqn{^\circ}C), relative humidity
#'     anomaly (`ES_dRH`, \%), VPD buffering (`ES_dVPD`, kPa), and 12-band
#'     canopy spectral transmittance ratios (`Alb_450`..`Alb_860`) relative to
#'     open-area reference sensors (`ref_species = "TTR"` or `ref_data`).
#'   \item \bold{Hydrological regulation}: whole-tree water uptake (`Flux`, L/h),
#'     crown-projected transpiration (`Flux_mm`, mm/h), and canopy conductance
#'     (`Gc`).
#'   \item \bold{Carbon sequestration & biomass increment}: basal area increment
#'     (`dBA_cm2`, \eqn{\text{cm}^2}), dry woody biomass increment (`dBiomass_kg`,
#'     kg), sequestered carbon (`dCarbon_kg`, kg C), and \eqn{\text{CO}_2}
#'     equivalent (`CO2_eq_kg`, kg \eqn{\text{CO}_2}).
#'   \item \bold{Biomechanical stability}: stem tilt anomaly (`phi_d`, degrees)
#'     and wind sway oscillation factor (`Osc_fct`, g).
#' }
#'
#' @param data A processed `data.frame` from [tt_process()] (or a raw list from
#'   [tt_types_to_list()], which will be passed to [tt_process()] automatically).
#' @param ref_species Character vector of species codes identifying open-area
#'   reference stations (`TT-Reference`, default `"TTR"`).
#' @param ref_id Optional character/numeric vector of `TT_ID` values acting as
#'   open-area reference stations.
#' @param ref_data Optional external meteorological `data.frame` with columns
#'   `Time` and `Tair_Ref` (and optionally `RH_Ref`, `VPD_Ref`).
#' @param wood_density Basic wood density in \eqn{\text{g cm}^{-3}} (or
#'   \eqn{\text{t m}^{-3}}, default `0.50`). Can also be a named numeric vector
#'   by `Species`.
#' @param carbon_fraction Carbon mass fraction of dry woody biomass (default `0.50`).
#' @param form_factor Stem form factor for volume calculation (default `0.45`).
#' @param tree_height_m Default tree height in meters when `Height` column is
#'   absent (default `15`).
#' @return A pure `base::data.frame` enriched with ecosystem service columns.
#' @export
tt_calc_ecosystem_services <- function(data, ref_species = "TTR", ref_id = NULL,
                                       ref_data = NULL, wood_density = 0.50,
                                       carbon_fraction = 0.50, form_factor = 0.45,
                                       tree_height_m = 15) {
  if (is.list(data) && !is.data.frame(data)) {
    data <- tt_process(data)
  }
  if (is.null(data) || nrow(data) == 0L) return(to_df(data))
  dt <- data.table::as.data.table(data)

  if (!"Time" %in% names(dt) && "Timestamp" %in% names(dt)) {
    dt[, Time := floor_datetime(Timestamp, 3600)]
  }
  if (!"VPD" %in% names(dt) && all(c("Tair", "RH") %in% names(dt))) {
    dt[, VPD := (0.6108 * exp((17.27 * Tair) / (Tair + 237.3))) * (1 - RH / 100)]
  }

  # 1. Microclimate & Radiation Interception Services relative to Reference (TTR)
  is_ref <- rep(FALSE, nrow(dt))
  if ("Species" %in% names(dt) && !is.null(ref_species)) {
    is_ref <- is_ref | (as.character(dt$Species) %in% as.character(ref_species))
  }
  if ("TT_ID" %in% names(dt) && !is.null(ref_id)) {
    is_ref <- is_ref | (as.character(dt$TT_ID) %in% as.character(ref_id))
  }

  spec_a_cols <- c(
    "AS7262_450_A", "AS7262_500_A", "AS7262_550_A", "AS7262_570_A", "AS7262_600_A", "AS7262_650_A",
    "AS7263_610_A", "AS7263_680_A", "AS7263_730_A", "AS7263_760_A", "AS7263_810_A", "AS7263_860_A"
  )
  present_spec <- intersect(spec_a_cols, names(dt))

  if (!is.null(ref_data) && is.data.frame(ref_data) && nrow(ref_data) > 0L) {
    ref_dt <- data.table::as.data.table(ref_data)
    if (!"Time" %in% names(ref_dt) && "Timestamp" %in% names(ref_dt)) {
      ref_dt[, Time := floor_datetime(Timestamp, 3600)]
    }
    if ("Tair" %in% names(ref_dt) && !"Tair_Ref" %in% names(ref_dt)) {
      data.table::setnames(ref_dt, "Tair", "Tair_Ref")
    }
    if ("RH" %in% names(ref_dt) && !"RH_Ref" %in% names(ref_dt)) {
      data.table::setnames(ref_dt, "RH", "RH_Ref")
    }
    if (!"VPD_Ref" %in% names(ref_dt) && all(c("Tair_Ref", "RH_Ref") %in% names(ref_dt))) {
      ref_dt[, VPD_Ref := (0.6108 * exp((17.27 * Tair_Ref) / (Tair_Ref + 237.3))) * (1 - RH_Ref / 100)]
    }
    keep_ref <- intersect(c("Time", "Tair_Ref", "RH_Ref", "VPD_Ref"), names(ref_dt))
    ref_summary <- unique(ref_dt[, ..keep_ref], by = "Time")
    dt <- merge(dt, ref_summary, by = "Time", all.x = TRUE, sort = FALSE)
  } else if (any(is_ref)) {
    ref_sub <- dt[is_ref]
    ref_summary <- ref_sub[, .(
      Tair_Ref = mean(as.numeric(Tair), na.rm = TRUE),
      RH_Ref   = mean(as.numeric(RH), na.rm = TRUE),
      VPD_Ref  = mean(as.numeric(VPD), na.rm = TRUE)
    ), by = Time]

    if (length(present_spec) > 0L) {
      spec_ref <- ref_sub[, lapply(.SD, function(x) mean(as.numeric(x), na.rm = TRUE)),
                          by = Time, .SDcols = present_spec]
      data.table::setnames(spec_ref, present_spec, sub("_A$", "_Ref", present_spec))
      ref_summary <- merge(ref_summary, spec_ref, by = "Time", all = TRUE, sort = FALSE)
    }
    # Remove any pre-existing *_Ref columns before merging
    old_ref_cols <- intersect(names(ref_summary)[-1L], names(dt))
    if (length(old_ref_cols) > 0L) dt[, (old_ref_cols) := NULL]
    dt <- merge(dt, ref_summary, by = "Time", all.x = TRUE, sort = FALSE)
  } else {
    # Fallback: use site hourly 90th percentile of Tair (upper canopy / unshaded proxy)
    # when no dedicated TTR station is present
    if (!"Tair_Ref" %in% names(dt) && "Tair" %in% names(dt)) {
      dt[, `:=`(
        Tair_Ref = stats::quantile(as.numeric(Tair), probs = 0.90, na.rm = TRUE),
        RH_Ref   = stats::quantile(as.numeric(RH), probs = 0.10, na.rm = TRUE),
        VPD_Ref  = stats::quantile(as.numeric(VPD), probs = 0.90, na.rm = TRUE)
      ), by = Time]
    }
  }

  if (all(c("Tair", "Tair_Ref") %in% names(dt))) {
    dt[, ES_dT := as.numeric(Tair) - as.numeric(Tair_Ref)]
  }
  if (all(c("RH", "RH_Ref") %in% names(dt))) {
    dt[, ES_dRH := as.numeric(RH) - as.numeric(RH_Ref)]
  }
  if (all(c("VPD", "VPD_Ref") %in% names(dt))) {
    dt[, ES_dVPD := as.numeric(VPD) - as.numeric(VPD_Ref)]
  }

  # Spectral canopy transmittance ratios (fixes Alb_760 bug from original ttfull.R)
  wavelengths <- c("450", "500", "550", "570", "600", "650",
                   "610", "680", "730", "760", "810", "860")
  sensor_pref <- c(rep("AS7262_", 6L), rep("AS7263_", 6L))
  for (i in seq_along(wavelengths)) {
    wl <- wavelengths[i]
    col_a <- paste0(sensor_pref[i], wl, "_A")
    col_ref <- paste0(sensor_pref[i], wl, "_Ref")
    col_alb <- paste0("Alb_", wl)
    if (all(c(col_a, col_ref) %in% names(dt))) {
      dt[, (col_alb) := {
        r_val <- as.numeric(get(col_a)) / as.numeric(get(col_ref))
        r_val[is.infinite(r_val) | r_val < 0 | r_val > 5] <- NA_real_
        r_val
      }]
    }
  }

  # 2. Carbon Sequestration & Woody Biomass Increment Service
  if (!"DBH" %in% names(dt)) dt[, DBH := 60]
  if (!"growth" %in% names(dt) && "growt_sens" %in% names(dt)) {
    dt <- data.table::as.data.table(calc_growth(dt))
  }
  if (!"growth" %in% names(dt)) dt[, growth := 0]

  r0_cm_vec <- resolve_diameter_cm(dt[["DBH"]]) / 2
  inc_cm_vec <- ifelse(is.na(dt[["growth"]]), 0, pmax(0, as.numeric(dt[["growth"]])) / 10)
  dt[, `:=`(
    BA_cm2  = pi * ((r0_cm_vec + inc_cm_vec)^2),
    dBA_cm2 = pi * ((r0_cm_vec + inc_cm_vec)^2 - (r0_cm_vec^2))
  )]

  # Resolve species-specific or scalar wood density (t m^-3 = g cm^-3)
  if (!is.null(names(wood_density)) && "Species" %in% names(dt)) {
    wd_vec <- unname(wood_density[as.character(dt$Species)])
    wd_vec[is.na(wd_vec)] <- mean(wood_density, na.rm = TRUE)
  } else {
    wd_vec <- as.numeric(wood_density)[1L]
  }
  h_vec <- if ("Height" %in% names(dt)) {
    ifelse(is.na(dt$Height) | dt$Height <= 0, tree_height_m, as.numeric(dt$Height))
  } else {
    tree_height_m
  }

  # Biomass expansion factor BEF = 1.25 (stem + branches)
  bef <- 1.25
  dt[, `:=`(
    Biomass_kg  = (BA_cm2 / 10000) * h_vec * form_factor * bef * wd_vec * 1000,
    dBiomass_kg = (dBA_cm2 / 10000) * h_vec * form_factor * bef * wd_vec * 1000
  )]
  dt[, `:=`(
    Carbon_kg  = Biomass_kg * carbon_fraction,
    dCarbon_kg = dBiomass_kg * carbon_fraction,
    CO2_eq_kg  = dBiomass_kg * carbon_fraction * (44 / 12)
  )]

  if ("Time" %in% names(dt) && "TT_ID" %in% names(dt)) {
    data.table::setorder(dt, TT_ID, Time)
  }
  to_df(dt)
}

#' Summarise TreeTalker Ecosystem Services by day, month, or season
#'
#' Aggregates hourly ecosystem service metrics into daily, monthly, or seasonal
#' tables per tree (`TT_ID`), `Species`, or `Site` using `data.table` and
#' returns a pure `base::data.frame`.
#'
#' @param data A `data.frame` returned by [tt_calc_ecosystem_services()].
#' @param period Aggregation time step: `"day"` (default), `"month"`, or
#'   `"season"` (entire monitoring period).
#' @param by Character vector of grouping columns (default
#'   `c("TT_ID", "Species", "Site")`, intersected with available columns).
#' @param hot_threshold Reference air temperature threshold (\eqn{^\circ}C) for
#'   evaluating peak heat-wave canopy cooling (default `25`).
#' @return A pure `base::data.frame` of aggregated ecosystem services.
#' @export
tt_summarise_services <- function(data, period = c("day", "month", "season"),
                                  by = c("TT_ID", "Species", "Site"),
                                  hot_threshold = 25) {
  period <- match.arg(period)
  if (is.null(data) || nrow(data) == 0L) return(data.frame())

  if (!"ES_dT" %in% names(data) && !"dCarbon_kg" %in% names(data)) {
    data <- tt_calc_ecosystem_services(data)
  }
  dt <- data.table::as.data.table(data)
  if (all(c("TT_ID", "Time") %in% names(dt))) {
    data.table::setorder(dt, TT_ID, Time)
    dt <- unique(dt, by = c("TT_ID", "Time"))
    dt[, dt_step_h := {
      d_h <- diff(as.numeric(Time)) / 3600
      d_pos <- d_h[d_h > 0 & d_h <= 6]
      if (length(d_pos) > 0L) stats::median(d_pos, na.rm = TRUE) else 1
    }, by = TT_ID]
  } else {
    dt[, dt_step_h := 1]
  }

  if (!"Date" %in% names(dt) && "Time" %in% names(dt)) {
    dt[, Date := as.Date(Time, tz = "UTC")]
  }
  if (!"Month" %in% names(dt) && "Time" %in% names(dt)) {
    dt[, Month := as.POSIXlt(Time, tz = "UTC")$mon + 1L]
  }

  grp_cols <- intersect(by, names(dt))
  if (period == "day" && "Date" %in% names(dt)) {
    grp_cols <- c(grp_cols, "Date")
  } else if (period == "month" && "Month" %in% names(dt)) {
    grp_cols <- c(grp_cols, "Month")
  }

  safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(as.numeric(x), na.rm = TRUE)
  safe_sum  <- function(x) if (all(is.na(x))) NA_real_ else sum(as.numeric(x), na.rm = TRUE)
  safe_max  <- function(x) if (all(is.na(x))) NA_real_ else max(as.numeric(x), na.rm = TRUE)
  safe_min  <- function(x) if (all(is.na(x))) NA_real_ else min(as.numeric(x), na.rm = TRUE)
  safe_sd   <- function(x) if (sum(!is.na(x)) < 2L) NA_real_ else stats::sd(as.numeric(x), na.rm = TRUE)

  # Ensure target columns exist so aggregation expressions stay uniform
  req_num <- c("Flux", "Flux_mm", "Fdh", "Gc", "ES_dT", "ES_dRH", "ES_dVPD",
               "Tair_Ref", "Tair", "NDVI", "StWC", "growth", "dBA_cm2",
               "dBiomass_kg", "dCarbon_kg", "CO2_eq_kg", "Osc_fct", "phi_d")
  for (rc in setdiff(req_num, names(dt))) dt[, (rc) := NA_real_]

  res <- dt[, .(
    n_obs          = .N,
    Flux_L_sum     = safe_sum(Flux * dt_step_h),
    Flux_mm_sum    = safe_sum(Flux_mm * dt_step_h),
    Fdh_mean       = safe_mean(Fdh),
    Fdh_max        = safe_max(Fdh),
    Gc_mean        = safe_mean(Gc),
    ES_dT_mean     = safe_mean(ES_dT),
    ES_dT_hot_mean = {
      t_ref <- ifelse(!is.na(Tair_Ref), Tair_Ref, Tair)
      hot_idx <- !is.na(t_ref) & t_ref >= hot_threshold
      if (any(hot_idx)) safe_mean(ES_dT[hot_idx]) else NA_real_
    },
    ES_dT_min      = safe_min(ES_dT),
    ES_dRH_mean    = safe_mean(ES_dRH),
    ES_dVPD_mean   = safe_mean(ES_dVPD),
    NDVI_mean      = safe_mean(NDVI),
    StWC_mean      = safe_mean(StWC),
    growth_mm_max  = safe_max(growth),
    dBA_cm2_max    = safe_max(dBA_cm2),
    dBiomass_kg    = safe_max(dBiomass_kg),
    dCarbon_kg     = safe_max(dCarbon_kg),
    CO2_eq_kg      = safe_max(CO2_eq_kg),
    Osc_fct_mean   = safe_mean(Osc_fct),
    Osc_fct_max    = safe_max(Osc_fct),
    phi_d_sd       = safe_sd(phi_d)
  ), by = grp_cols]

  to_df(res)
}

