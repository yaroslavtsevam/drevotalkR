#' Calculate daily zero-flow temperature difference (dTmax)
#'
#' Computes the zero-flow baseline temperature difference \eqn{\Delta T_{\max}}
#' for thermal dissipation probes using standard daily maxima or the three
#' methods from the `TREX` methodology (`"pd"`, `"mw"`, `"dr"`).
#'
#' @details
#' **Source attribution**: Methods `"pd"` (Pre-Dawn), `"mw"` (Moving Window),
#' and `"dr"` (Double Regression) are adapted from `TREX::tdm_dt.max` (Richard
#' Peters, Christoforos Pappas et al., 2021, *Methods in Ecology and Evolution*,
#' 12(5), 942-950, GPL-3) and implemented on `data.table`.
#'
#' @param data A `data.frame` or `data.table` containing at least `TT_ID`,
#'   `Time` (or `Timestamp`), and `dT` (or `Theat_1`, `Tref_1`, `Theat_0`, `Tref_0`).
#' @param method Character string specifying the \eqn{\Delta T_{\max}} method:
#'   `"default"` (daily maximum), `"pd"` (pre-dawn 00:00--08:00 maximum with low
#'   VPD filter), `"mw"` (rolling maximum over `window_days`), or `"dr"` (double
#'   regression over `window_days`).
#' @param window_days Integer window length in days for `"mw"` and `"dr"`
#'   (default `7`).
#' @param vpd_threshold Maximum pre-dawn `VPD` (kPa) for `"pd"` filtering when
#'   `VPD` is present (default `0.5`).
#' @return A pure `base::data.frame` with columns `dT`, `dTmax`, and
#'   method-specific baselines (`dTmax_pd`, `dTmax_mw`, `dTmax_dr`).
#' @references
#' Peters, R. L., Pappas, C., Hurley, A. G., Poyatos, R., Flo, V., Zweifel, R.,
#' Steppe, K., & Fonti, P. (2021). Assimilate, process and analyse thermal
#' dissipation sap flow data using the TREX R package. *Methods in Ecology and
#' Evolution*, 12(5), 942-950.
#' @export
calc_dTmax <- function(data, method = c("default", "pd", "mw", "dr"),
                       window_days = 7L, vpd_threshold = 0.5) {
  method <- match.arg(method)
  if (is.null(data) || nrow(data) == 0L) return(to_df(data))
  dt <- data.table::as.data.table(data)

  ts_col <- intersect(c("Time", "Datetime", "Timestamp"), names(dt))[1L]
  if (is.na(ts_col)) stop("Data must contain a 'Time' or 'Timestamp' column.")
  if (!"TT_ID" %in% names(dt)) dt[, TT_ID := "1"]

  if (!"dT" %in% names(dt)) {
    if (all(c("Theat_1", "Tref_1", "Theat_0", "Tref_0") %in% names(dt))) {
      dt[, dTon := as.numeric(Theat_1) - as.numeric(Tref_1)]
      dt[, dToff := as.numeric(Theat_0) - as.numeric(Tref_0)]
      dt[, dT := ifelse(!is.na(dToff) & abs(dToff) < 15, dTon - dToff, dTon)]
    } else if (all(c("Theat_1", "Tref_1") %in% names(dt))) {
      dt[, dT := as.numeric(Theat_1) - as.numeric(Tref_1)]
    } else {
      stop("Data must contain 'dT' or thermistor columns ('Theat_1', 'Tref_1').")
    }
  }
  # Filter physically implausible dT values
  dt[dT <= 0 | dT > 35, dT := NA_real_]

  dt[, Date := as.Date(get(ts_col), tz = "UTC")]
  dt[, Hour := as.POSIXlt(get(ts_col), tz = "UTC")$hour]
  data.table::setorderv(dt, c("TT_ID", ts_col))

  # 1. Default daily max per TT_ID & Date
  daily_dt <- dt[, .(
    dTmax_def = if (any(!is.na(dT))) max(dT, na.rm = TRUE) else NA_real_,
    dTmax_pd = {
      pd_sub <- Hour <= 8L & !is.na(dT)
      if ("VPD" %in% names(.SD) && any(pd_sub & !is.na(VPD) & VPD <= vpd_threshold)) {
        max(dT[pd_sub & !is.na(VPD) & VPD <= vpd_threshold], na.rm = TRUE)
      } else if (any(pd_sub)) {
        max(dT[pd_sub], na.rm = TRUE)
      } else if (any(!is.na(dT))) {
        max(dT, na.rm = TRUE)
      } else {
        NA_real_
      }
    }
  ), by = .(TT_ID, Date)]

  w <- max(1L, as.integer(window_days))
  daily_dt[, `:=`(
    dTmax_pd = fill_small_gaps(dTmax_pd, max_gap = w * 2L, method = "linear")
  ), by = TT_ID]

  # 2. Moving window (mw) and Double regression (dr) on daily series
  daily_dt[, dTmax_mw := {
    v <- dTmax_pd
    if (all(is.na(v))) v <- dTmax_def
    n_v <- length(v)
    k_w <- min(w, n_v)
    mw_vec <- data.table::frollapply(v, k_w, FUN = function(z) {
      if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE)
    }, align = "center", fill = NA_real_)
    # Edge fill
    fill_small_gaps(
      ifelse(is.na(mw_vec), v, mw_vec),
      max_gap = w * 2L, method = "linear"
    )
  }, by = TT_ID]

  daily_dt[, dTmax_dr := {
    v <- dTmax_mw
    n_v <- length(v)
    k_w <- min(w, n_v)
    mean_vec <- data.table::frollmean(v, n = k_w, align = "center", na.rm = TRUE)
    mean_vec <- ifelse(is.na(mean_vec), v, mean_vec)
    # Filter raw daily maxima below the first regression baseline and recompute
    v_filt <- ifelse(!is.na(dTmax_pd) & dTmax_pd >= mean_vec, dTmax_pd, mean_vec)
    dr_vec <- data.table::frollmean(v_filt, n = k_w, align = "center", na.rm = TRUE)
     ifelse(is.na(dr_vec), v_filt, dr_vec)
  }, by = TT_ID]

  dt <- merge(dt, daily_dt, by = c("TT_ID", "Date"), all.x = TRUE, sort = FALSE)
  data.table::setorderv(dt, c("TT_ID", ts_col))

  dt[, dTmax := switch(
    method,
    default = dTmax_def,
    pd = ifelse(is.na(dTmax_pd), dTmax_def, dTmax_pd),
    mw = ifelse(is.na(dTmax_mw), dTmax_def, dTmax_mw),
    dr = ifelse(is.na(dTmax_dr), dTmax_def, dTmax_dr)
  )]
  # Ensure dTmax >= dT so ((dTmax - dT)/dT) is non-negative
  dt[!is.na(dT) & !is.na(dTmax) & dTmax < dT, dTmax := dT]

  to_df(dt)
}

#' Calculate sap flux density and whole-tree transpiration
#'
#' Applies Granier (1985) and Doerfler (2019) calibrations to TreeTalker thermal
#' dissipation probe signals and scales sap flux to whole-tree (`Flux`, L/h) and
#' crown-projected (`Flux_mm`, mm/h) water use.
#'
#' @param data A `data.frame` containing thermistor temperatures (`Tref_0`,
#'   `Theat_0`, `Tref_1`, `Theat_1`) or `dT`, plus `TT_ID` and `Time`/`Timestamp`.
#' @param dtmax_method Method for \eqn{\Delta T_{\max}}: `"default"`, `"pd"`,
#'   `"mw"`, or `"dr"` (see [calc_dTmax()]).
#' @param window_days Rolling window in days for `"mw"` and `"dr"` (default `7`).
#' @return A pure `base::data.frame` with added sap flow columns (`dTon`,
#'   `dToff`, `dT`, `dTmax`, `Fd`, `Fdh`, `SFD_Gra`, `SFD_Doe`, `Flux`,
#'   `Flux_mm`, `Gc`).
#' @references
#' Granier, A. (1985). Une nouvelle methode pour la mesure du flux de seve brute
#' dans le tronc des arbres. *Annales des Sciences Forestieres*, 42(2), 193-200.
#' @export
calc_sap_flux_density <- function(data, dtmax_method = c("default", "pd", "mw", "dr"),
                                  window_days = 7L) {
  dtmax_method <- match.arg(dtmax_method)
  if (is.null(data) || nrow(data) == 0L) return(to_df(data))

  df_dt <- calc_dTmax(data, method = dtmax_method, window_days = window_days)
  dt <- data.table::as.data.table(df_dt)

  # Granier index K = (dTmax - dT) / dT
  dt[, `:=`(
    Fd = {
      k_val <- pmax(0, (dTmax - dT) / dT)
      fd_val <- 0.00011899 * (k_val^1.231)
      fd_val[is.infinite(fd_val) | fd_val < 0 | fd_val > 0.001] <- NA_real_
      fd_val
    }
  )]
  # Fdh in L m^-2 h^-1; SFD_Gra in cm^3 cm^-2 h^-1; SFD_Doe in L m^-2 h^-1
  dt[, `:=`(
    Fdh = Fd * 3600 * 1000,
    SFD_Gra = Fd * 360,
    SFD_Doe = {
      k_val <- pmax(0, (dTmax - dT) / dT)
      doe <- 15.68 * (k_val^1.19) * 10
      doe[is.infinite(doe) | doe < 0] <- NA_real_
      doe
    }
  )]

  # Whole-tree sap flux (L h^-1) and crown-projected transpiration (mm h^-1)
  if (!"DBH" %in% names(dt)) dt[, DBH := 60]
  diam_cm_vec <- resolve_diameter_cm(dt[["DBH"]])
  dt[, `:=`(
    Flux = Fdh * ((diam_cm_vec / 100)^1.8777) * 0.755,
    Flux_mm = {
      flux_l <- Fdh * ((diam_cm_vec / 100)^1.8777) * 0.755
      crown_area_m2 <- pmax(1, pi * (0.12 * diam_cm_vec)^2)
      flux_l / crown_area_m2
    }
  )]

  # Canopy stomatal conductance proxy Gc = Fdh / VPD
  if ("VPD" %in% names(dt)) {
    dt[, Gc := ifelse(!is.na(VPD) & VPD > 0.05, Fdh / VPD, NA_real_)]
  } else {
    dt[, Gc := NA_real_]
  }

  to_df(dt)
}

#' Calculate stem volumetric water content (StWC)
#'
#' Computes temperature-compensated volumetric wood moisture content (`StWC`, %)
#' from TreeTalker capacitive frequency sensor (`MoistSensHz`) and reference
#' wood temperature (`Tref_0`) using species-specific calibrations for v3 and v2
#' sensors.
#'
#' @param data A `data.frame` with `MoistSensHz`, `Tref_0`, and optionally
#'   `Species`.
#' @param version Sensor hardware version: `"v3"` (packet `4D`) or `"v2"`
#'   (packet `45`).
#' @return A pure `base::data.frame` with column `StWC` (%).
#' @export
calc_wood_moist <- function(data, version = c("v3", "v2")) {
  version <- match.arg(version)
  if (is.null(data) || nrow(data) == 0L) return(to_df(data))
  dt <- data.table::as.data.table(data)

  if (!"Species" %in% names(dt)) dt[, Species := "Unknown"]
  hz <- as.numeric(dt[["MoistSensHz"]])
  tref <- as.numeric(dt[["Tref_0"]])
  sp <- as.character(dt[["Species"]])

  if (version == "v3") {
    # Temperature-corrected dimensionless moisture index
    ec_t <- -0.09 * tref - 6.7e-4 * (hz - 20 * tref) + 18.4
    stwc <- ifelse(
      grepl("Picea|Larix|\u0415\u043b\u044c|\u041b\u0438\u0441\u0442\u0432\u0435\u043d\u043d\u0438\u0446\u0430", sp, ignore.case = TRUE),
      (ec_t - 7.8017) / 0.0556,
      ifelse(
        grepl("Quercus|Acer|Fagus|\u0414\u0443\u0431|\u041a\u043b\u0435\u043d|\u041a\u043b\u0451\u043d|\u0411\u0443\u043a", sp, ignore.case = TRUE),
        (ec_t - 8.1112) / 0.0366,
        ifelse(
          grepl("Pinus|\u0421\u043e\u0441\u043d\u0430", sp, ignore.case = TRUE),
          (ec_t - 7.0398) / 0.0725,
          (ec_t - 8.0592) / 0.0454
        )
      )
    )
  } else {
    # v2 calibration
    stwc <- ifelse(
      grepl("Pinus|Picea|Larix|Abies|\u0421\u043e\u0441\u043d\u0430|\u0415\u043b\u044c|\u041b\u0438\u0441\u0442\u0432\u0435\u043d\u043d\u0438\u0446\u0430", sp, ignore.case = TRUE),
      152.7 - 0.0085 * hz + 0.14 * tref,
      138.4 - 0.0078 * hz + 0.12 * tref
    )
  }
  stwc[stwc <= 0 | stwc > 100 | is.infinite(stwc)] <- NA_real_
  dt[, StWC := stwc]
  to_df(dt)
}

#' Calculate soil moisture and water potential (pF) for TT-Soil sensors
#'
#' Computes soil volumetric water content (`Wsoil`, `Wchern`, `Wdern`, `Wash`)
#' and matric potential (`pFchern`, `pFdern`) from packet `55` frequency and
#' temperature channels.
#'
#' @param data A `data.frame` of TT-Soil (`55`) records with `MoistSensHz` and
#'   `Tref_0`.
#' @return A pure `base::data.frame` with soil moisture and pF columns.
#' @export
calc_soil_moist <- function(data) {
  if (is.null(data) || nrow(data) == 0L) return(to_df(data))
  dt <- data.table::as.data.table(data)

  hz <- as.numeric(dt[["MoistSensHz"]])
  tsoil <- as.numeric(dt[["Tref_0"]])

  # Calibrations from ttfull.R for chernozem, sod-podzolic (dern), and ash/urban soils
  dt[, `:=`(
    Wchern = pmax(0, pmin(100, 146.5 - 0.0076 * hz + 0.15 * tsoil)),
    Wdern  = pmax(0, pmin(100, 132.0 - 0.0069 * hz + 0.12 * tsoil)),
    Wash   = pmax(0, pmin(100, 139.2 - 0.0072 * hz + 0.13 * tsoil))
  )]
  dt[, `:=`(
    Wsoil   = Wdern,
    pFchern = ifelse(Wchern > 1, 6.5 - 1.35 * log10(Wchern), NA_real_),
    pFdern  = ifelse(Wdern > 1, 6.2 - 1.40 * log10(Wdern), NA_real_)
  )]
  to_df(dt)
}

#' Calculate stem radial growth and stem radius from IR distance sensor
#'
#' Converts raw infrared distance counts (`growt_sens`) to distance in mm,
#' cleans outliers using [clean_growth_series()], and computes cumulative radial
#' stem increment (`growth`, mm) and current stem radius (`R`, mm).
#'
#' @param data A `data.frame` containing `TT_ID`, `Time` (or `Timestamp`), and
#'   `growt_sens`.
#' @return A pure `base::data.frame` with columns `growth_cmV`, `growth_cmVb`,
#'   `growth` (mm), and `R` (mm).
#' @export
calc_growth <- function(data) {
  if (is.null(data) || nrow(data) == 0L) return(to_df(data))
  dt <- data.table::as.data.table(data)

  if (!"growt_sens" %in% names(dt)) {
    dt[, `:=`(growth_cmV = NA_real_, growth_cmVb = NA_real_, growth = NA_real_, R = NA_real_)]
    return(to_df(dt))
  }
  if (!"TT_ID" %in% names(dt)) dt[, TT_ID := "1"]
  if (!"DBH" %in% names(dt)) dt[, DBH := 60]

  # Polynomial conversion from ADC counts to distance (mm)
  gs <- as.numeric(dt[["growt_sens"]])
  gs[gs <= 100 | gs >= 120000] <- NA_real_
  dist_mm <- (1.132e-9 * (gs^2) - 2.299e-4 * gs + 14.37) * 10
  dt[, growth_cmV := dist_mm]

  dt[, `:=`(
    growth_cmVb = clean_growth_series(growth_cmV, k = 7L, t0 = 3)
  ), by = TT_ID]

  # Radial growth is decrease in bark-to-sensor distance relative to initial baseline
  dt[, growth := {
    g_clean <- growth_cmVb
    valid <- !is.na(g_clean)
    if (sum(valid) >= 3L) {
      init_dist <- stats::median(head(g_clean[valid], 10L), na.rm = TRUE)
      inc <- init_dist - g_clean
      inc_smooth <- sg_filter(inc, m = 0L, p = 2L, w = min(15L, sum(valid)))
      pmax(0, inc_smooth)
    } else {
      rep(NA_real_, .N)
    }
  }, by = TT_ID]

  r0_mm_vec <- resolve_diameter_cm(dt[["DBH"]]) * 5
  dt[, R := ifelse(!is.na(growth), r0_mm_vec + growth, r0_mm_vec)]

  to_df(dt)
}

#' Full biophysical processing pipeline for TreeTalker data
#'
#' Processes raw typed tables from [tt_types_to_list()] (`mydata_4D`, `mydata_45`,
#' `mydata_49`, and `mydata_55`), performs sensor calibrations (temperature,
#' battery voltage, VPD, sap flow, wood/soil moisture, radial growth, 12-band
#' spectra, 13 vegetation indices, and 3D stem tilt angles) using `data.table`
#' internally, and returns a unified pure `base::data.frame`.
#'
#' @param mydata_list Output list from [tt_types_to_list()] (or [tt_scrape()]).
#' @param step_time Time grid resolution in seconds (default `3600`).
#' @param grow_calc Logical; compute filtered radial growth series (default `TRUE`).
#' @param extrap_dates Logical; run [tt_extrapolate_dates()] to repair broken
#'   RTC timestamps (default `FALSE`).
#' @param tz_offset Optional numeric timezone offset in hours from UTC.
#' @param dtmax_method Method for \eqn{\Delta T_{\max}}: `"default"`, `"pd"`,
#'   `"mw"`, or `"dr"`.
#' @param lat_default Default latitude (decimal degrees) for solar altitude
#'   filtering of spectral sensors when `lat` is absent from metadata (default `55.75`).
#' @param lon_default Default longitude (decimal degrees, default `37.62`).
#' @return A unified pure `base::data.frame` containing all calibrated sensor
#'   channels and biophysical variables.
#' @export
tt_process <- function(mydata_list, step_time = 3600, grow_calc = TRUE,
                       extrap_dates = FALSE, tz_offset = NULL,
                       dtmax_method = c("default", "pd", "mw", "dr"),
                       lat_default = 55.75, lon_default = 37.62) {
  dtmax_method <- match.arg(dtmax_method)
  if (is.data.frame(mydata_list)) {
    mydata_list <- tt_types_to_list(mydata_list)
  }

  # Helper to process main stem/soil sensor tables (4D, 55, 45)
  proc_main <- function(df_in, ver = "v3", is_soil = FALSE) {
    if (is.null(df_in) || nrow(df_in) == 0L) return(data.table::data.table())
    if (isTRUE(extrap_dates)) {
      df_in <- tt_extrapolate_dates(df_in, step_sec = step_time, tz_offset = tz_offset)
    }
    dt <- data.table::as.data.table(df_in)
    if (!"Time" %in% names(dt)) {
      dt[, Time := floor_datetime(Timestamp, step_sec = step_time)]
    }

    if (ver == "v3") {
      dt[, `:=`(
        Bat_mV  = 2 * 1100 * as.numeric(adc_Vbat) / as.numeric(adc_bandgap),
        Tref_0  = conv_analog_t_v3(Tref_0),
        Theat_0 = conv_analog_t_v3(Theat_0),
        Tref_1  = conv_analog_t_v3(Tref_1),
        Theat_1 = conv_analog_t_v3(Theat_1),
        Tair    = as.numeric(Tair) / 10,
        RH      = pmin(100, pmax(0, as.numeric(RH) / 10))
      )]
    } else {
      # v2 (packet 45)
      dt[, `:=`(
        Bat_mV  = 2 * (1100 * 131072 / as.numeric(adc_bandgap)),
        Tref_0  = tt_norm_analog_t(Tref_0, adc_bandgap) / 10,
        Theat_0 = tt_norm_analog_t(Theat_0, adc_bandgap) / 10,
        Tref_1  = tt_norm_analog_t(Tref_1, adc_bandgap) / 10,
        Theat_1 = tt_norm_analog_t(Theat_1, adc_bandgap) / 10,
        Tair    = as.numeric(Tair) / 10,
        RH      = pmin(100, pmax(0, as.numeric(RH)))
      )]
    }
    dt[, Volt := Bat_mV / 1000]
    dt[Tair < -50 | Tair > 65, Tair := NA_real_]
    dt[Tref_0 < -45 | Tref_0 > 60, Tref_0 := NA_real_]
    dt[Tref_1 < -45 | Tref_1 > 60, Tref_1 := NA_real_]
    dt[Theat_0 < -45 | Theat_0 > 75, Theat_0 := NA_real_]
    dt[Theat_1 < -45 | Theat_1 > 75, Theat_1 := NA_real_]

    # Vapor Pressure Deficit VPD (kPa)
    dt[, VPD := (0.6108 * exp((17.27 * Tair) / (Tair + 237.3))) * (1 - RH / 100)]
    dt[VPD < 0, VPD := 0]

    # 3D Accelerometer / Clinometer inclinations (degrees) and wind sway
    rad2deg <- 180 / pi
    gx <- as.numeric(dt[["gx_mean"]])
    gy <- as.numeric(dt[["gy_mean"]])
    gz <- as.numeric(dt[["gz_mean"]])
    dt[, `:=`(
      Theta   = atan2(gx, sqrt(gy^2 + gz^2)) * rad2deg,
      Psi     = atan2(gy, sqrt(gx^2 + gz^2)) * rad2deg,
      Phi     = atan2(sqrt(gx^2 + gy^2), gz) * rad2deg,
      Accel   = sqrt(gx^2 + gy^2 + gz^2) / 16384,
      Osc_fct = (abs(as.numeric(gx_sd)) + abs(as.numeric(gy_sd)) + abs(as.numeric(gz_sd))) / 16384
    )]
    dt[, phi_d := Phi - mean(Phi, na.rm = TRUE), by = TT_ID]

    if (isTRUE(is_soil)) {
      dt <- data.table::as.data.table(calc_soil_moist(dt))
    } else {
      dt <- data.table::as.data.table(calc_sap_flux_density(dt, dtmax_method = dtmax_method))
      dt <- data.table::as.data.table(calc_wood_moist(dt, version = ver))
      if (isTRUE(grow_calc)) {
        dt <- data.table::as.data.table(calc_growth(dt))
      }
    }
    dt
  }

  dt_4d <- proc_main(mydata_list$mydata_4D, ver = "v3", is_soil = FALSE)
  dt_55 <- proc_main(mydata_list$mydata_55, ver = "v3", is_soil = TRUE)
  dt_45 <- proc_main(mydata_list$mydata_45, ver = "v2", is_soil = FALSE)

  dt_main <- data.table::rbindlist(list(dt_4d, dt_55, dt_45), fill = TRUE, use.names = TRUE)
  if (nrow(dt_main) > 0L) {
    dt_main <- unique(dt_main, by = c("TT_ID", "Time"))
  }

  # Process spectral packet 49
  dt_49 <- data.table::as.data.table(mydata_list$mydata_49)
  if (!is.null(dt_49) && nrow(dt_49) > 0L) {
    if (isTRUE(extrap_dates)) {
      dt_49 <- data.table::as.data.table(
        tt_extrapolate_dates(dt_49, step_sec = step_time, tz_offset = tz_offset, align_by_tair = FALSE)
      )
    }
    if (!"Time" %in% names(dt_49)) {
      dt_49[, Time := floor_datetime(Timestamp, step_sec = step_time)]
    }

    # Filter negative or saturated spectral counts
    spec_raw_cols <- c(
      "AS7262_450", "AS7262_500", "AS7262_550", "AS7262_570", "AS7262_600", "AS7262_650",
      "AS7263_610", "AS7263_680", "AS7263_730", "AS7263_760", "AS7263_810", "AS7263_860"
    )
    for (sc in intersect(spec_raw_cols, names(dt_49))) {
      dt_49[get(sc) <= 0 | get(sc) >= 65000, (sc) := NA_real_]
    }

    # Apply laboratory linear calibration equations (from ttfull.R / ttlight.R)
    # Fixes channel shift bug in ttlight.R where 760/810/860 were misaligned
    dt_49[, `:=`(
      AS7262_450_A = -32.79   + 0.3518  * AS7262_450,
      AS7262_500_A = -29.88   + 0.2383  * AS7262_500,
      AS7262_550_A = -54.022  + 0.2549  * AS7262_550,
      AS7262_570_A = -57.768  + 0.2564  * AS7262_570,
      AS7262_600_A = -48.88   + 0.2414  * AS7262_600,
      AS7262_650_A = -24.167  + 0.2909  * AS7262_650,
      AS7263_610_A = -58.4081 + 0.1651  * AS7263_610,
      AS7263_680_A = -59.5665 + 0.2323  * AS7263_680,
      AS7263_730_A = -25.19   + 0.1977  * AS7263_730,
      AS7263_760_A = -19.22   + 0.1626  * AS7263_760,
      AS7263_810_A = -6.71    + 0.10557 * AS7263_810,
      AS7263_860_A =  8.00    + 0.08825 * AS7263_860
    )]
    spec_a_cols <- paste0(spec_raw_cols, "_A")
    for (ac in spec_a_cols) dt_49[get(ac) <= 0, (ac) := NA_real_]

    # Solar altitude filter (retain daytime measurements where sun is above horizon)
    lat_vec <- if ("lat" %in% names(dt_49)) ifelse(is.na(dt_49$lat), lat_default, as.numeric(dt_49$lat)) else lat_default
    lon_vec <- if ("lon" %in% names(dt_49)) ifelse(is.na(dt_49$lon), lon_default, as.numeric(dt_49$lon)) else lon_default
    sun_pos <- get_sun_position(dt_49$Time, lat_vec, lon_vec)
    dt_49[, `:=`(altitude = sun_pos$altitude, azimuth = sun_pos$azimuth)]

    # Vegetation Indices (using calibrated channels)
    b450 <- dt_49$AS7262_450_A
    b500 <- dt_49$AS7262_500_A
    b550 <- dt_49$AS7262_550_A
    b650 <- dt_49$AS7262_650_A
    b680 <- dt_49$AS7263_680_A
    b730 <- dt_49$AS7263_730_A
    b760 <- dt_49$AS7263_760_A
    b810 <- dt_49$AS7263_810_A
    b860 <- dt_49$AS7263_860_A

    dt_49[, `:=`(
      NDVI  = (b810 - b680) / (b810 + b680),
      NDVIm = (b760 - b650) / (b760 + b650),
      NDVIr = (b810 - b650) / (b810 + b650),
      NDRE  = (b760 - b730) / (b760 + b730),
      MSAVI = (2 * b810 + 1 - sqrt(pmax(0, (2 * b810 + 1)^2 - 8 * (b810 - b680)))) / 2,
      ReCI  = (b810 / b730) - 1,
      PSRI  = (b680 - b500) / b760,
      RVI   = b810 / b680,
      TVI   = 0.5 * (120 * (b760 - b550) - 200 * (b680 - b550)),
      GNDVI = (b760 - b550) / (b760 + b550),
      EVI   = 2.5 * (b810 - b680) / (b810 + 6 * b680 - 7.5 * b450 + 1),
      VARI  = (b550 - b650) / (b550 + b650 - b450),
      OSAVI = (b860 - b680) / (b860 + b680 + 0.16)
    )]
    for (idx_col in c("NDVI", "NDVIm", "NDVIr", "NDRE", "GNDVI", "VARI", "OSAVI")) {
      dt_49[get(idx_col) < -1 | get(idx_col) > 1, (idx_col) := NA_real_]
    }
    # Mask nighttime measurements where solar altitude <= 0
    idx_all_cols <- c("NDVI", "NDVIm", "NDVIr", "NDRE", "MSAVI", "ReCI", "PSRI",
                      "RVI", "TVI", "GNDVI", "EVI", "VARI", "OSAVI")
    dt_49[!is.na(altitude) & altitude <= 0, (c(spec_a_cols, idx_all_cols)) := NA_real_]

    spec_keep <- c(
      "TT_ID", "Time", spec_raw_cols, spec_a_cols,
      "integration_T", "gain", "altitude", "azimuth",
      "NDVI", "NDVIm", "NDVIr", "NDRE", "MSAVI", "ReCI", "PSRI", "RVI",
      "TVI", "GNDVI", "EVI", "VARI", "OSAVI"
    )
    spec_keep <- intersect(spec_keep, names(dt_49))
    dt_49_sub <- unique(dt_49[, ..spec_keep], by = c("TT_ID", "Time"))

    if (nrow(dt_main) > 0L) {
      dt_main <- merge(dt_main, dt_49_sub, by = c("TT_ID", "Time"), all = TRUE, sort = FALSE)
    } else {
      dt_main <- dt_49_sub
    }
  }

  if (nrow(dt_main) > 0L) {
    dt_main[, `:=`(
      Date  = as.Date(Time, tz = "UTC"),
      DOY   = as.POSIXlt(Time, tz = "UTC")$yday + 1L,
      Hour  = as.POSIXlt(Time, tz = "UTC")$hour,
      Month = as.POSIXlt(Time, tz = "UTC")$mon + 1L
    )]
    data.table::setorder(dt_main, TT_ID, Time)
  }

  to_df(dt_main)
}

#' @rdname tt_process
#' @export
ttfull_calc <- tt_process

#' Legacy single-module TreeTalker wrappers returning pure data.frames
#'
#' Backward-compatible wrappers (`ttbattery`, `ttgranier`, `ttStWC`, `ttgrowth`,
#' `ttlight`, `ttmeteo`, `ttstability`) implemented on `data.table` and
#' returning pure `base::data.frame` objects.
#'
#' @param mydata A list from [tt_types_to_list()] or a processed `data.frame`.
#' @return A pure `base::data.frame`.
#' @name tt_legacy_modules
NULL

#' @rdname tt_legacy_modules
#' @export
ttbattery <- function(mydata) {
  df <- if (is.data.frame(mydata) && "Bat_mV" %in% names(mydata)) to_df(mydata) else tt_process(mydata, grow_calc = FALSE)
  keep <- intersect(c("TT_ID", "Species", "Site", "Time", "Date", "Bat_mV", "Volt", "adc_Vbat", "adc_bandgap"), names(df))
  df[, keep, drop = FALSE]
}

#' @rdname tt_legacy_modules
#' @export
ttgranier <- function(mydata) {
  df <- if (is.data.frame(mydata) && "Fd" %in% names(mydata)) to_df(mydata) else tt_process(mydata, grow_calc = FALSE)
  keep <- intersect(c("TT_ID", "Species", "Site", "DBH", "Time", "Date", "DOY", "Hour",
                      "Tref_0", "Theat_0", "Tref_1", "Theat_1", "dTon", "dToff", "dT", "dTmax",
                      "Fd", "Fdh", "SFD_Gra", "SFD_Doe", "Flux", "Flux_mm", "Gc"), names(df))
  df[, keep, drop = FALSE]
}

#' @rdname tt_legacy_modules
#' @export
ttStWC <- function(mydata) {
  df <- if (is.data.frame(mydata) && "StWC" %in% names(mydata)) to_df(mydata) else tt_process(mydata, grow_calc = FALSE)
  keep <- intersect(c("TT_ID", "Species", "Site", "Time", "Date", "MoistSensHz", "Tref_0", "StWC"), names(df))
  df[, keep, drop = FALSE]
}

#' @rdname tt_legacy_modules
#' @export
ttgrowth <- function(mydata) {
  df <- if (is.data.frame(mydata) && "growth" %in% names(mydata)) to_df(mydata) else tt_process(mydata, grow_calc = TRUE)
  keep <- intersect(c("TT_ID", "Species", "Site", "DBH", "Time", "Date", "DOY",
                      "growt_sens", "growth_cmV", "growth_cmVb", "growth", "R"), names(df))
  df[, keep, drop = FALSE]
}

#' @rdname tt_legacy_modules
#' @export
ttlight <- function(mydata) {
  df <- if (is.data.frame(mydata) && "NDVI" %in% names(mydata)) to_df(mydata) else tt_process(mydata, grow_calc = FALSE)
  keep <- intersect(c("TT_ID", "Species", "Site", "Time", "Date", "DOY", "Hour",
                      grep("^AS726", names(df), value = TRUE),
                      "NDVI", "NDVIm", "NDVIr", "NDRE", "MSAVI", "ReCI", "PSRI",
                      "RVI", "TVI", "GNDVI", "EVI", "VARI", "OSAVI"), names(df))
  df[, keep, drop = FALSE]
}

#' @rdname tt_legacy_modules
#' @export
ttmeteo <- function(mydata) {
  df <- if (is.data.frame(mydata) && "VPD" %in% names(mydata)) to_df(mydata) else tt_process(mydata, grow_calc = FALSE)
  keep <- intersect(c("TT_ID", "Species", "Site", "Time", "Date", "DOY", "Hour", "Tair", "RH", "VPD"), names(df))
  df[, keep, drop = FALSE]
}

#' @rdname tt_legacy_modules
#' @export
ttstability <- function(mydata) {
  df <- if (is.data.frame(mydata) && "Phi" %in% names(mydata)) to_df(mydata) else tt_process(mydata, grow_calc = FALSE)
  keep <- intersect(c("TT_ID", "Species", "Site", "Time", "Date", "DOY",
                      "gx_mean", "gy_mean", "gz_mean", "gx_sd", "gy_sd", "gz_sd",
                      "Theta", "Psi", "Phi", "phi_d", "Accel", "Osc_fct"), names(df))
  df[, keep, drop = FALSE]
}
