#' @import data.table
#' @import ggplot2
#' @importFrom stats setNames
#' @importFrom utils head tail
NULL

# Register non-standard evaluation variables for R CMD check
utils::globalVariables(c(
  ".", ".N", ".SD", ".I", ".GRP", ":=",
  "TT_ID", "Server_time", "Time", "Date", "Datetime", "DOY", "Hour", "Month",
  "Species", "Site", "DBH", "Name", "V1", "V2", "V3", "V4", "V5", "V6", "V7",
  "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18",
  "V19", "V20", "V21", "V22", "V23", "Dev_Ty", "Rec_Nr", "Timestamp", "Time_s",
  "Tref_0", "Theat_0", "Tref_1", "Theat_1", "growt_sens", "adc_bandgap",
  "Bits", "RH", "Tair", "gz_mean", "gz_sd", "gy_mean", "gy_sd", "gx_mean",
  "gx_sd", "MoistSensHz", "adc_Vbat", "Bat_mV", "Volt", "VPD", "dTon", "dToff",
  "dT", "dTmax", "dTmax_pd", "dTmax_mw", "dTmax_dr", "Fd", "Fdh", "Flux",
  "Flux_mm", "SFD_Gra", "SFD_Doe", "Gc", "StWC", "Wsoil", "Wchern", "Wdern",
  "Wash", "pFchern", "pFdern", "growth", "growth_cmV", "growth_cmVb", "R",
  "Theta", "Psi", "Phi", "phi_d", "Accel", "Osc_fct", "NDVI", "NDRE", "MSAVI",
  "ReCI", "PSRI", "GNDVI", "EVI", "VARI", "OSAVI", "AS7262_450", "AS7262_500",
  "AS7262_550", "AS7262_570", "AS7262_600", "AS7262_650", "AS7263_610",
  "AS7263_680", "AS7263_730", "AS7263_760", "AS7263_810", "AS7263_860",
  "AS7262_450_A", "AS7262_500_A", "AS7262_550_A", "AS7262_570_A",
  "AS7262_600_A", "AS7262_650_A", "AS7263_610_A", "AS7263_680_A",
  "AS7263_730_A", "AS7263_760_A", "AS7263_810_A", "AS7263_860_A",
  "Tair_Ref", "RH_Ref", "VPD_Ref", "ES_dT", "ES_dRH", "ES_dVPD",
  "Alb_450", "Alb_500", "Alb_550", "Alb_570", "Alb_600", "Alb_650",
  "Alb_610", "Alb_680", "Alb_730", "Alb_760", "Alb_810", "Alb_860",
  "BA_cm2", "dBA_cm2", "Biomass_kg", "dBiomass_kg", "Carbon_kg", "dCarbon_kg",
  "CO2_eq_kg", "ischarged", "charge_cycle", "altitude", "azimuth",
  "val", "mean_val", "sd_val", "se_val", "cum_val", "group_var", "facet_var",
  "x_val", "y_val", "..keep_cols", "..keep_ref", "..spec_keep",
  "Flux_L_sum", "Group", "Time_hr", "dTmax_def", "mean_flux", "n_obs",
  "dt_step_h"
))

#' Convert a table to a pure base::data.frame
#'
#' Internal helper ensuring that all user-facing functions return a pure
#' `data.frame` without `data.table` or `tbl_df` classes.
#'
#' @param x A `data.table`, `data.frame`, or `list`.
#' @return A pure `base::data.frame`.
#' @keywords internal
to_df <- function(x) {
  if (is.null(x)) return(data.frame())
  if (data.table::is.data.table(x)) {
    data.table::setDF(x)
    return(x)
  }
  as.data.frame(x, stringsAsFactors = FALSE)
}

#' Convert DBH vector to stem diameter in cm without row-wise threshold jumps
#'
#' Resolves whether a `DBH` vector represents stem circumference (cm) or stem
#' diameter (cm) across the dataset (using the stand median when
#' `dbh_unit = "auto"`), avoiding a 10x basal-area jump at 40 cm within a stand.
#'
#' @param dbh Numeric vector of `DBH` values from device metadata.
#' @param dbh_unit Unit of `dbh`: `"auto"` (default; circumference if stand
#'   median > 40, else diameter), `"circumference_cm"`, or `"diameter_cm"`.
#' @return Numeric vector of stem diameters in centimeters.
#' @keywords internal
resolve_diameter_cm <- function(dbh, dbh_unit = c("auto", "circumference_cm", "diameter_cm")) {
  dbh_unit <- match.arg(dbh_unit)
  dbh_num <- as.numeric(dbh)
  dbh_num[is.na(dbh_num) | dbh_num <= 0] <- 60
  if (dbh_unit == "auto") {
    med_val <- stats::median(dbh_num, na.rm = TRUE)
    dbh_unit <- if (!is.na(med_val) && med_val <= 40) "diameter_cm" else "circumference_cm"
  }
  if (dbh_unit == "circumference_cm") dbh_num / pi else dbh_num
}

#' Fast vectorized POSIXct floor rounding
#'
#' Replaces `lubridate::floor_date()` using integer division of epoch seconds.
#'
#' @param x A `POSIXct` vector.
#' @param step_sec Step size in seconds (default `3600` for 1 hour).
#' @return A `POSIXct` vector rounded down to the nearest multiple of `step_sec`.
#' @export
floor_datetime <- function(x, step_sec = 3600) {
  if (!inherits(x, "POSIXct")) {
    x <- as.POSIXct(x, tz = "UTC")
  }
  tz <- attr(x, "tzone")
  if (is.null(tz) || identical(tz, "")) tz <- "UTC"
  as.POSIXct(floor(as.numeric(x) / step_sec) * step_sec,
             origin = "1970-01-01", tz = tz)
}

#' Savitzky-Golay smoothing filter
#'
#' Computes Savitzky-Golay polynomial smoothing or derivatives for a numeric
#' vector without requiring the external `prospectr` or `signal` packages.
#'
#' @details
#' **Source attribution**: Adapted from `prospectr::savitzkyGolay` (Antoine
#' Stevens and Leonardo Ramirez-Lopez, GPL-3) and the original convolution
#' least-squares formulation of Savitzky & Golay (1964), *Analytical Chemistry*,
#' 36(8), 1627-1639. Uses Vandermonde pseudoinverse weights (`qr.solve`) applied
#' via `stats::filter`, with linear edge extrapolation.
#'
#' @param x Numeric vector to smooth.
#' @param m Differentiation order (`0` = smoothing only).
#' @param p Polynomial order (default `3`).
#' @param w Window size (must be odd and `>= p + 1`, default `11`).
#' @return Numeric vector of the same length as `x`.
#' @references
#' Savitzky, A., & Golay, M. J. E. (1964). Smoothing and Differentiation of Data
#' by Simplified Least Squares Procedures. *Analytical Chemistry*, 36(8),
#' 1627-1639.
#' @export
sg_filter <- function(x, m = 0L, p = 3L, w = 11L) {
  x <- as.numeric(x)
  n <- length(x)
  w <- as.integer(w)
  p <- as.integer(p)
  m <- as.integer(m)
  if (w %% 2L == 0L) w <- w + 1L
  if (w <= p) w <- p + (1L + (p %% 2L == 1L))
  if (n < w || sum(!is.na(x)) < w) return(x)

  # Interpolate internal NAs temporarily for convolution
  na_idx <- is.na(x)
  x_work <- x
  if (any(na_idx) && sum(!na_idx) >= 2L) {
    x_work <- stats::approx(seq_len(n)[!na_idx], x[!na_idx],
                            xout = seq_len(n), rule = 2)$y
  }

  half_w <- (w - 1L) %/% 2L
  k <- seq.int(-half_w, half_w)
  # Vandermonde design matrix S_{i,j} = k_i^j for j = 0..p
  S <- outer(k, 0L:p, `^`)
  # Pseudoinverse (S^T S)^{-1} S^T
  G <- qr.solve(t(S) %*% S, t(S))
  # Filter coefficients for m-th derivative
  coeffs <- factorial(m) * G[m + 1L, ]

  # Pad edges by reflection to preserve length n
  left_pad <- 2 * x_work[1L] - x_work[seq.int(half_w + 1L, 2L, by = -1L)]
  right_pad <- 2 * x_work[n] - x_work[seq.int(n - 1L, n - half_w, by = -1L)]
  x_pad <- c(left_pad, x_work, right_pad)

  # stats::filter applies rev(coeffs) in convolution mode
  filtered <- as.numeric(stats::filter(x_pad, rev(coeffs), sides = 1L))
  out <- filtered[seq.int(w, w + n - 1L)]
  out[na_idx] <- NA_real_
  out
}

#' Calculate solar position (altitude and azimuth)
#'
#' Computes solar altitude and azimuth for given timestamps and geographic
#' coordinates without requiring the `suncalc` package.
#'
#' @details
#' **Source attribution**: Adapted from `suncalc::getSunlightPosition` (Benoit
#' Thieurmel & Achraf Elmarhraoui, GPL-2; originally `suncalc.js` by Vladimir
#' Agafonkin) based on astronomical formulas by Jean Meeus (1998), *Astronomical
#' Algorithms*, 2nd ed., Willmann-Bell (Chapters 25 and 47).
#'
#' @param date A `POSIXct` vector of timestamps.
#' @param lat Numeric latitude in decimal degrees (North positive).
#' @param lon Numeric longitude in decimal degrees (East positive).
#' @return A pure `data.frame` with columns `date`, `lat`, `lon`, `altitude`
#'   (radians), `azimuth` (radians), and `altitude_deg` (degrees above horizon).
#' @references
#' Meeus, J. (1998). *Astronomical Algorithms* (2nd ed.). Willmann-Bell.
#' @export
get_sun_position <- function(date, lat, lon) {
  if (!inherits(date, "POSIXct")) {
    date <- as.POSIXct(date, tz = "UTC")
  }
  rad <- pi / 180
  # Days since J2000.0 (2000-01-01 12:00 UTC)
  d <- as.numeric(date) / 86400 - 0.5 + 2440588 - 2451545

  # Solar mean anomaly
  M <- rad * (357.5291 + 0.98560028 * d)
  # Equation of center
  C <- rad * (1.9148 * sin(M) + 0.0200 * sin(2 * M) + 0.0003 * sin(3 * M))
  # Ecliptic longitude
  L <- M + C + rad * 102.9372 + pi
  # Obliquity of the Earth
  e <- rad * 23.4397

  # Declination and right ascension
  dec <- asin(sin(L) * sin(e))
  ra <- atan2(sin(L) * cos(e), cos(L))

  # Sidereal time and hour angle
  lw <- rad * (-lon)
  phi <- rad * lat
  theta <- rad * (280.16 + 360.9856235 * d) - lw
  H <- theta - ra

  altitude <- asin(sin(phi) * sin(dec) + cos(phi) * cos(dec) * cos(H))
  azimuth <- atan2(sin(H), cos(H) * sin(phi) - tan(dec) * cos(phi))

  data.frame(
    date = date,
    lat = lat,
    lon = lon,
    altitude = as.numeric(altitude),
    azimuth = as.numeric(azimuth),
    altitude_deg = as.numeric(altitude / rad),
    stringsAsFactors = FALSE
  )
}

#' Interpolate short gaps in a time series
#'
#' Fills missing (`NA`) runs of length at most `max_gap` using linear or spline
#' interpolation while leaving longer gaps as `NA`.
#'
#' @details
#' **Source attribution**: Adapted from the bounded gap-interpolation logic in
#' `baytrends::fillMissing` (Rebecca Murphy et al., GPL-3) and `zoo::na.approx`
#' (Achim Zeileis & Gabor Grothendieck, GPL-2), implemented in base R using
#' `rle(is.na(x))` and `stats::approx` / `stats::spline`.
#'
#' @param x Numeric vector with potential `NA` values.
#' @param max_gap Maximum contiguous `NA` run length to interpolate (default `7`).
#' @param method Interpolation method: `"linear"` (default) or `"spline"`.
#' @return Numeric vector with short gaps filled.
#' @export
fill_small_gaps <- function(x, max_gap = 7L, method = c("linear", "spline")) {
  method <- match.arg(method)
  x <- as.numeric(x)
  n <- length(x)
  na_mask <- is.na(x)
  if (!any(na_mask) || sum(!na_mask) < 2L) return(x)

  runs <- rle(na_mask)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L

  # Identify which NA positions belong to runs > max_gap or leading/trailing NAs
  long_na <- rep(FALSE, n)
  for (i in seq_along(runs$lengths)) {
    if (runs$values[i]) {
      is_edge <- (starts[i] == 1L) || (ends[i] == n)
      if (runs$lengths[i] > max_gap || is_edge) {
        long_na[starts[i]:ends[i]] <- TRUE
      }
    }
  }

  valid_idx <- which(!na_mask)
  if (method == "linear") {
    interp <- stats::approx(valid_idx, x[valid_idx], xout = seq_len(n), rule = 1)$y
  } else {
    interp <- stats::spline(valid_idx, x[valid_idx], xout = seq_len(n))$y
  }
  out <- x
  fill_pos <- na_mask & !long_na
  out[fill_pos] <- interp[fill_pos]
  out
}

#' Clean radial growth sensor series (Hampel + MAD filter)
#'
#' Removes spikes and step artifacts from TreeTalker IR distance sensor series
#' without requiring `pracma`, `forecast`, or `changepoint`.
#'
#' @details
#' **Source attribution**: Combines the sliding Hampel median filter from
#' `pracma::hampel` (Hans W. Borchers, GPL-3; Davies & Gather, 1993) and
#' robust MAD residual outlier replacement from `forecast::tsclean` (Rob J.
#' Hyndman, GPL-3).
#'
#' @param x Numeric vector of raw distance/growth measurements.
#' @param k Half-window size for running median (window `2*k + 1`, default `7`).
#' @param t0 Threshold multiplier of MAD (default `3`).
#' @return Cleaned numeric vector.
#' @references
#' Davies, L., & Gather, U. (1993). The identification of multiple outliers.
#' *Journal of the American Statistical Association*, 88(423), 782-792.
#' @export
clean_growth_series <- function(x, k = 7L, t0 = 3) {
  x <- as.numeric(x)
  n <- length(x)
  valid <- !is.na(x)
  if (sum(valid) < 5L) return(x)

  w <- 2L * as.integer(k) + 1L
  if (w > sum(valid)) w <- max(3L, 2L * (sum(valid) %/% 4L) + 1L)

  x_filled <- fill_small_gaps(x, max_gap = 24L, method = "linear")
  med_idx <- !is.na(x_filled)
  out <- x
  if (sum(med_idx) >= w) {
    med_vals <- stats::runmed(x_filled[med_idx], k = w, endrule = "keep")
    abs_dev <- abs(x_filled[med_idx] - med_vals)
    # Local MAD scale estimate (1.4826 * running median of |x - med|)
    local_mad <- 1.4826 * stats::runmed(abs_dev, k = w, endrule = "keep")
    global_mad <- stats::median(abs_dev, na.rm = TRUE) * 1.4826
    scale_est <- ifelse(local_mad > 1e-6, local_mad, max(global_mad, 1e-6))

    outliers <- abs_dev > (t0 * scale_est)
    cleaned_sub <- x_filled[med_idx]
    cleaned_sub[outliers] <- med_vals[outliers]
    out[med_idx] <- cleaned_sub
  }
  out[!valid] <- NA_real_
  out
}

#' Bilingual label dictionary (EN / RU)
#'
#' Internal lookup table for plot titles, axis labels, and legends with English
#' (`"en"`) as default and Russian (`"ru"`) as optional switch.
#'
#' @param key String identifier of the label.
#' @param lang Language code: `"en"` (default) or `"ru"`.
#' @return Localized character string.
#' @keywords internal
tt_i18n <- function(key, lang = c("en", "ru")) {
  lang <- match.arg(lang)
  dict <- list(
    time = c(en = "Date / Time", ru = "\u0414\u0430\u0442\u0430 / \u0412\u0440\u0435\u043c\u044f"),
    hour = c(en = "Hour of Day", ru = "\u0427\u0430\u0441 \u0441\u0443\u0442\u043e\u043a"),
    doy = c(en = "Day of Year (DOY)", ru = "\u0414\u0435\u043d\u044c \u0433\u043e\u0434\u0430 (DOY)"),
    month = c(en = "Month", ru = "\u041c\u0435\u0441\u044f\u0446"),
    species = c(en = "Species", ru = "\u041f\u043e\u0440\u043e\u0434\u0430"),
    site = c(en = "Site", ru = "\u041e\u0431\u044a\u0435\u043a\u0442"),
    tt_id = c(en = "Device ID (TT_ID)", ru = "\u041d\u043e\u043c\u0435\u0440 \u0443\u0441\u0442\u0440\u043e\u0439\u0441\u0442\u0432\u0430 (TT_ID)"),
    Bat_mV = c(en = "Battery Voltage (mV)", ru = "\u041d\u0430\u043f\u0440\u044f\u0436\u0435\u043d\u0438\u0435 \u0430\u043a\u043a\u0443\u043c\u0443\u043b\u044f\u0442\u043e\u0440\u0430 (\u043c\u0412)"),
    Volt = c(en = "Battery Voltage (V)", ru = "\u041d\u0430\u043f\u0440\u044f\u0436\u0435\u043d\u0438\u0435 \u0430\u043a\u043a\u0443\u043c\u0443\u043b\u044f\u0442\u043e\u0440\u0430 (\u0412)"),
    Tair = c(en = "Air Temperature (\u00B0C)", ru = "\u0422\u0435\u043c\u043f\u0435\u0440\u0430\u0442\u0443\u0440\u0430 \u0432\u043e\u0437\u0434\u0443\u0445\u0430 (\u00B0C)"),
    RH = c(en = "Relative Humidity (%)", ru = "\u041e\u0442\u043d\u043e\u0441\u0438\u0442\u0435\u043b\u044c\u043d\u0430\u044f \u0432\u043b\u0430\u0436\u043d\u043e\u0441\u0442\u044c (%)"),
    VPD = c(en = "Vapor Pressure Deficit, VPD (kPa)", ru = "\u0414\u0435\u0444\u0438\u0446\u0438\u0442 \u0443\u043f\u0440\u0443\u0433\u043e\u0441\u0442\u0438 \u043f\u0430\u0440\u0430, VPD (\u043a\u041f\u0430)"),
    dT = c(en = "Probe Temperature Difference, dT (\u00B0C)", ru = "\u0420\u0430\u0437\u043d\u043e\u0441\u0442\u044c \u0442\u0435\u043c\u043f\u0435\u0440\u0430\u0442\u0443\u0440 \u0437\u043e\u043d\u0434\u043e\u0432, dT (\u00B0C)"),
    Fd = c(en = "Sap Flux Density, Fd (m3 m-2 s-1)", ru = "\u041f\u043b\u043e\u0442\u043d\u043e\u0441\u0442\u044c \u0441\u043e\u043a\u043e\u0434\u0432\u0438\u0436\u0435\u043d\u0438\u044f, Fd (\u043c3 \u043c-2 \u0441-1)"),
    Fdh = c(en = "Sap Flux Density, Fdh (L m-2 h-1)", ru = "\u041f\u043b\u043e\u0442\u043d\u043e\u0441\u0442\u044c \u0441\u043e\u043a\u043e\u0434\u0432\u0438\u0436\u0435\u043d\u0438\u044f, Fdh (\u043b \u043c-2 \u0447-1)"),
    Flux = c(en = "Whole-Tree Sap Flux (L h-1)", ru = "\u0421\u043e\u043a\u043e\u0434\u0432\u0438\u0436\u0435\u043d\u0438\u0435 \u0434\u0435\u0440\u0435\u0432\u0430 (\u043b \u0447-1)"),
    Flux_mm = c(en = "Canopy Transpiration (mm h-1)", ru = "\u0422\u0440\u0430\u043d\u0441\u043f\u0438\u0440\u0430\u0446\u0438\u044f \u043a\u0440\u043e\u043d\u044b (\u043c\u043c \u0447-1)"),
    Gc = c(en = "Canopy Conductance, Gc (L m-2 h-1 kPa-1)", ru = "\u0423\u0441\u0442\u044c\u0438\u0447\u043d\u0430\u044f \u043f\u0440\u043e\u0432\u043e\u0434\u0438\u043c\u043e\u0441\u0442\u044c, Gc (\u043b \u043c-2 \u0447-1 \u043a\u041f\u0430-1)"),
    StWC = c(en = "Stem Volumetric Water Content, StWC (%)", ru = "\u041e\u0431\u044a\u0451\u043c\u043d\u0430\u044f \u0432\u043b\u0430\u0436\u043d\u043e\u0441\u0442\u044c \u0434\u0440\u0435\u0432\u0435\u0441\u0438\u043d\u044b \u0441\u0442\u0432\u043e\u043b\u0430, StWC (%)"),
    Wsoil = c(en = "Soil Volumetric Water Content (%)", ru = "\u041e\u0431\u044a\u0451\u043c\u043d\u0430\u044f \u0432\u043b\u0430\u0436\u043d\u043e\u0441\u0442\u044c \u043f\u043e\u0447\u0432\u044b (%)"),
    growth = c(en = "Cumulative Radial Growth (mm)", ru = "\u041d\u0430\u043a\u043e\u043f\u043b\u0435\u043d\u043d\u044b\u0439 \u0440\u0430\u0434\u0438\u0430\u043b\u044c\u043d\u044b\u0439 \u043f\u0440\u0438\u0440\u043e\u0441\u0442 (\u043c\u043c)"),
    Phi = c(en = "Stem Inclination Angle, Phi (\u00B0)", ru = "\u0423\u0433\u043e\u043b \u043d\u0430\u043a\u043b\u043e\u043d\u0430 \u0441\u0442\u0432\u043e\u043b\u0430, Phi (\u00B0)"),
    phi_d = c(en = "Inclination Deviation, dPhi (\u00B0)", ru = "\u041e\u0442\u043a\u043b\u043e\u043d\u0435\u043d\u0438\u0435 \u0443\u0433\u043b\u0430 \u043d\u0430\u043a\u043b\u043e\u043d\u0430, dPhi (\u00B0)"),
    Osc_fct = c(en = "Wind Sway Oscillation Factor (g)", ru = "\u0424\u0430\u043a\u0442\u043e\u0440 \u0432\u0435\u0442\u0440\u043e\u0432\u043e\u0433\u043e \u0440\u0430\u0441\u043a\u0430\u0447\u0438\u0432\u0430\u043d\u0438\u044f (g)"),
    NDVI = c(en = "Canopy Transmitted NDVI", ru = "\u041f\u043e\u0434\u043a\u0440\u043e\u043d\u043e\u0432\u044b\u0439 \u0438\u043d\u0434\u0435\u043a\u0441 NDVI"),
    NDRE = c(en = "Canopy Red-Edge Index (NDRE)", ru = "\u0418\u043d\u0434\u0435\u043a\u0441 \u043a\u0440\u0430\u0441\u043d\u043e\u0433\u043e \u043a\u0440\u0430\u044f (NDRE)"),
    ES_dT = c(en = "Canopy Cooling Effect, dT (\u00B0C)", ru = "\u041e\u0445\u043b\u0430\u0436\u0434\u0430\u044e\u0449\u0438\u0439 \u044d\u0444\u0444\u0435\u043a\u0442 \u043a\u0440\u043e\u043d\u044b, dT (\u00B0C)"),
    ES_dRH = c(en = "Canopy Humidifying Effect, dRH (%)", ru = "\u0423\u0432\u043b\u0430\u0436\u043d\u044f\u044e\u0449\u0438\u0439 \u044d\u0444\u0444\u0435\u043a\u0442 \u043a\u0440\u043e\u043d\u044b, dRH (%)"),
    ES_dVPD = c(en = "Canopy VPD Buffering, dVPD (kPa)", ru = "\u0411\u0443\u0444\u0435\u0440\u0438\u0437\u0430\u0446\u0438\u044f VPD \u043a\u0440\u043e\u043d\u043e\u0439, dVPD (\u043a\u041f\u0430)"),
    dCarbon_kg = c(en = "Sequestered Carbon Increment (kg C)", ru = "\u0414\u0435\u043f\u043e\u043d\u0438\u0440\u043e\u0432\u0430\u043d\u043d\u044b\u0439 \u0443\u0433\u043b\u0435\u0440\u043e\u0434 \u043f\u0440\u0438\u0440\u043e\u0441\u0442\u0430 (\u043a\u0433 C)"),
    CO2_eq_kg = c(en = "Sequestered CO2 Equivalent (kg CO2)", ru = "\u042d\u043a\u0432\u0438\u0432\u0430\u043b\u0435\u043d\u0442 \u043f\u043e\u0433\u043b\u043e\u0449\u0451\u043d\u043d\u043e\u0433\u043e CO2 (\u043a\u0433 CO2)"),
    cum_flux = c(en = "Cumulative Transpiration (L)", ru = "\u041d\u0430\u043a\u043e\u043f\u043b\u0435\u043d\u043d\u0430\u044f \u0442\u0440\u0430\u043d\u0441\u043f\u0438\u0440\u0430\u0446\u0438\u044f (\u043b)"),
    cum_flux_mm = c(en = "Cumulative Transpiration (mm)", ru = "\u041d\u0430\u043a\u043e\u043f\u043b\u0435\u043d\u043d\u0430\u044f \u0442\u0440\u0430\u043d\u0441\u043f\u0438\u0440\u0430\u0446\u0438\u044f (\u043c\u043c)"),
    title_timeseries = c(en = "TreeTalker Time Series", ru = "\u0412\u0440\u0435\u043c\u0435\u043d\u043d\u043e\u0439 \u0440\u044f\u0434 TreeTalker"),
    title_diurnal = c(en = "Diurnal Cycle (Mean \u00B1 SE)", ru = "\u0421\u0443\u0442\u043e\u0447\u043d\u044b\u0439 \u0445\u043e\u0434 (\u0421\u0440\u0435\u0434\u043d\u0435\u0435 \u00B1 SE)"),
    title_cooling = c(en = "Ecosystem Service: Microclimate Cooling", ru = "\u042d\u043a\u043e\u0441\u0438\u0441\u0442\u0435\u043c\u043d\u044b\u0439 \u0441\u0435\u0440\u0432\u0438\u0441: \u041e\u0445\u043b\u0430\u0436\u0434\u0435\u043d\u0438\u0435 \u043c\u0438\u043a\u0440\u043e\u043a\u043b\u0438\u043c\u0430\u0442\u0430"),
    title_transpiration = c(en = "Ecosystem Service: Water Regulation (Transpiration)", ru = "\u042d\u043a\u043e\u0441\u0438\u0441\u0442\u0435\u043c\u043d\u044b\u0439 \u0441\u0435\u0440\u0432\u0438\u0441: \u0420\u0435\u0433\u0443\u043b\u044f\u0446\u0438\u044f \u0432\u043e\u0434\u043d\u043e\u0433\u043e \u0431\u0430\u043b\u0430\u043d\u0441\u0430"),
    title_carbon = c(en = "Ecosystem Service: Carbon Sequestration", ru = "\u042d\u043a\u043e\u0441\u0438\u0441\u0442\u0435\u043c\u043d\u044b\u0439 \u0441\u0435\u0440\u0432\u0438\u0441: \u0414\u0435\u043f\u043e\u043d\u0438\u0440\u043e\u0432\u0430\u043d\u0438\u0435 \u0443\u0433\u043b\u0435\u0440\u043e\u0434\u0430"),
    title_stability = c(en = "Ecosystem Service: Biomechanical Wind Stability", ru = "\u042d\u043a\u043e\u0441\u0438\u0441\u0442\u0435\u043c\u043d\u044b\u0439 \u0441\u0435\u0440\u0432\u0438\u0441: \u0412\u0435\u0442\u0440\u043e\u0443\u0441\u0442\u043e\u0439\u0447\u0438\u0432\u043e\u0441\u0442\u044c"),
    title_waterbalance = c(en = "Cumulative Tree Water Balance", ru = "\u041a\u0443\u043c\u0443\u043b\u044f\u0442\u0438\u0432\u043d\u044b\u0439 \u0432\u043e\u0434\u043d\u044b\u0439 \u0431\u0430\u043b\u0430\u043d\u0441 \u0434\u0435\u0440\u0435\u0432\u044c\u0435\u0432")
  )
  if (!key %in% names(dict)) return(key)
  unname(dict[[key]][lang])
}
