#' Calibrate TT+ v2 thermistor reading using internal ADC bandgap reference
#'
#' Vectorized pure-R implementation of the TreeTalker v2 thermistor lookup-table
#' bandgap normalization (`ADC_256x_TBL`), replacing the C++ (`tt.cpp`) and
#' row-by-row R loops without requiring `Rcpp`.
#'
#' @param temp_tenths Numeric vector of raw temperature values in tenths of a
#'   degree Celsius (as transmitted in packet `45`).
#' @param bandgap Numeric vector of internal reference `adc_bandgap` counts.
#' @return Numeric vector of bandgap-corrected temperatures in tenths of a
#'   degree Celsius (divide by 10 for degrees Celsius).
#' @export
tt_norm_analog_t <- function(temp_tenths, bandgap) {
  adc_tbl <- c(
    14289, 14199, 14105, 14009, 13911, 13809, 13705, 13599, 13489, 13378,
    13263, 13146, 13025, 12903, 12778, 12651, 12521, 12389, 12255, 12119,
    11981, 11842, 11701, 11558, 11413, 11267, 11119, 10970, 10820, 10669,
    10516, 10363, 10209, 10054,  9899,  9743,  9587,  9431,  9275,  9119,
     8963,  8808,  8653,  8499,  8345,  8192,  8039,  7887,  7736,  7586,
     7438,  7290,  7145,  7000,  6857,  6716,  6577,  6439,  6303,  6168,
     6036,  5905,  5776,  5649,  5523,  5400,  5279,  5159,  5042,  4927,
     4814,  4703,  4593,  4486,  4381,  4279,  4178,  4079,  3982,  3887,
     3794,  3705,  3617,  3531,  3448,  3366,  3286,  3208,  3132,  3057,
     2984,  2913,  2843,  2775,  2709,  2644,  2580,  2519,  2458,  2400,
     2342,  2287,  2232,  2179,  2127,  2077,  2028,  1980,  1933,  1888,
     1843,  1800,  1757,  1716,  1675,  1636,  1597,  1560,  1524,  1488,
     1453
  )
  t_grid <- seq(-200, 1000, by = 10)

  temp_tenths <- as.numeric(temp_tenths)
  bandgap <- as.numeric(bandgap)
  n <- max(length(temp_tenths), length(bandgap))
  if (n == 0L) return(numeric(0))
  temp_tenths <- rep_len(temp_tenths, n)
  bandgap <- rep_len(bandgap, n)

  valid <- !is.na(temp_tenths) & !is.na(bandgap) & bandgap > 0
  out <- rep(NA_real_, n)
  if (!any(valid)) return(out)

  # Step 1: Degrees2ADC piecewise linear interpolation
  adc_raw <- stats::approx(t_grid, adc_tbl, xout = temp_tenths[valid], rule = 2)$y
  # Step 2: Bandgap scaling ((1100 << 17) / 3300 = 43690.67)
  adc_scaled <- adc_raw * 43691 / bandgap[valid]
  # Step 3: ADC2Degrees inverse piecewise linear lookup
  out[valid] <- stats::approx(rev(adc_tbl), rev(t_grid), xout = adc_scaled, rule = 2)$y
  out
}

#' @rdname tt_norm_analog_t
#' @export
Normalize <- tt_norm_analog_t

#' Convert TT+ v3 thermistor ADC counts to temperature in Celsius
#'
#' Applies the cubic calibration polynomial for TreeTalker v3 thermistors.
#'
#' @param dn Numeric vector of raw digital number (ADC) readings.
#' @return Numeric vector of temperatures in degrees Celsius.
#' @export
conv_analog_t_v3 <- function(dn) {
  dn <- as.numeric(dn)
  127.6 - 0.006045 * dn + 1.26e-7 * (dn^2) - 1.15e-12 * (dn^3)
}

#' Reconstruct and align corrupted TreeTalker RTC timestamps
#'
#' Detects battery recharge events, identifies records with reset/corrupted
#' internal clocks (e.g., year 1970 or before 2015), estimates the device
#' logging interval (`3600` or `5400` s), extrapolates timestamps within each
#' charge cycle, and optionally aligns unanchored series by cross-correlating
#' air temperature (`Tair`) with reference devices at the same site.
#'
#' @param data A `data.frame` of TreeTalker records (e.g., `mydata_4D` or
#'   `mydata_45`).
#' @param step_sec Default logging interval in seconds when it cannot be
#'   inferred from valid timestamps (default `3600`).
#' @param tz_offset Optional numeric UTC offset in hours (e.g. `3` for MSK) or
#'   `NULL` to keep UTC / infer from `lon` (`round(lon / 15)`).
#' @param align_by_tair Logical; if `TRUE`, attempts cross-correlation (`ccf`)
#'   alignment of `Tair` for charge cycles lacking any valid timestamp anchor.
#' @return A pure `base::data.frame` with reconstructed `Timestamp` and `Time`
#'   columns.
#' @export
tt_extrapolate_dates <- function(data, step_sec = 3600, tz_offset = NULL,
                                 align_by_tair = TRUE) {
  if (is.null(data) || nrow(data) == 0L) return(to_df(data))
  dt <- data.table::as.data.table(data)

  ts_col <- if ("Timestamp" %in% names(dt)) "Timestamp" else if ("Time" %in% names(dt)) "Time" else NULL
  if (is.null(ts_col)) return(to_df(dt))

  if (!inherits(dt[[ts_col]], "POSIXct")) {
    ts_num <- suppressWarnings(as.numeric(dt[[ts_col]]))
    dt[, (ts_col) := as.POSIXct(ts_num, origin = "1970-01-01", tz = "UTC")]
  }

  # Detect battery charge cycles per TT_ID
  vbat_col <- intersect(c("adc_Vbat", "Bat_mV", "Volt"), names(dt))[1L]
  if (!is.na(vbat_col) && "TT_ID" %in% names(dt)) {
    thresh <- if (vbat_col == "Volt") 0.25 else 120
    dt[, ischarged := {
      dv <- c(0, diff(as.numeric(get(vbat_col))))
      dv[is.na(dv)] <- 0
      as.integer(dv > thresh)
    }, by = TT_ID]
    dt[, charge_cycle := cumsum(ischarged), by = TT_ID]
  } else {
    dt[, `:=`(ischarged = 0L, charge_cycle = 0L)]
  }

  min_valid <- as.POSIXct("2016-01-01 00:00:00", tz = "UTC")
  max_valid <- Sys.time() + 86400 * 30

  # Reconstruct timestamps within each (TT_ID, charge_cycle)
  dt[, (ts_col) := {
    ts_vec <- get(ts_col)
    srv_vec <- if ("Server_time" %in% names(.SD)) Server_time else rep(as.POSIXct(NA, tz = "UTC"), .N)
    # Use Server_time where RTC Timestamp is invalid but Server_time is valid
    valid_rtc <- !is.na(ts_vec) & ts_vec >= min_valid & ts_vec <= max_valid
    valid_srv <- !is.na(srv_vec) & srv_vec >= min_valid & srv_vec <= max_valid
    ts_work <- ts_vec
    ts_work[!valid_rtc & valid_srv] <- srv_vec[!valid_rtc & valid_srv]
    valid_now <- !is.na(ts_work) & ts_work >= min_valid & ts_work <= max_valid

    if (any(valid_now) && !all(valid_now)) {
      # Estimate step size (3600 vs 5400) from valid consecutive timestamps
      diffs <- diff(as.numeric(ts_work[valid_now]))
      diffs_pos <- diffs[diffs > 600 & diffs < 10800]
      est_step <- if (length(diffs_pos) > 0L) {
        med_d <- stats::median(diffs_pos, na.rm = TRUE)
        if (abs(med_d - 5400) < abs(med_d - 3600)) 5400 else 3600
      } else {
        step_sec
      }
      anchor_idx <- which(valid_now)[1L]
      anchor_time <- ts_work[anchor_idx]
      idx_seq <- seq_len(.N)
      extrap <- anchor_time + (idx_seq - anchor_idx) * est_step
      ts_work[!valid_now] <- extrap[!valid_now]
    }
    ts_work
  }, by = .(TT_ID, charge_cycle)]

  # Optional cross-correlation alignment by Tair for completely unanchored cycles
  if (isTRUE(align_by_tair) && "Tair" %in% names(dt)) {
    valid_rows <- !is.na(dt[[ts_col]]) & dt[[ts_col]] >= min_valid & !is.na(dt[["Tair"]])
    if (any(valid_rows) && any(!valid_rows)) {
      ref_dt <- dt[valid_rows, .(Tair_ref = mean(as.numeric(Tair), na.rm = TRUE)),
                   by = .(Time_hr = floor_datetime(get(ts_col), step_sec))]
      data.table::setorder(ref_dt, Time_hr)
      if (nrow(ref_dt) >= 24L) {
        dt[, (ts_col) := {
          ts_vec <- get(ts_col)
          if (all(is.na(ts_vec) | ts_vec < min_valid) && .N >= 12L && sum(!is.na(Tair)) >= 12L) {
            t_sub <- as.numeric(Tair)
            t_sub <- fill_small_gaps(t_sub, max_gap = 12L)
            cc <- suppressWarnings(stats::ccf(ref_dt$Tair_ref, t_sub, lag.max = nrow(ref_dt) - 1L,
                                              plot = FALSE, na.action = stats::na.pass))
            if (!is.null(cc$acf) && any(!is.na(cc$acf))) {
              best_lag <- cc$lag[which.max(cc$acf)]
              ts_vec <- ref_dt$Time_hr[1L] + (as.integer(best_lag) + seq_len(.N) - 1L) * step_sec
            }
          }
          ts_vec
        }, by = .(TT_ID, charge_cycle)]
      }
    }
  }

  # Apply timezone offset if specified or inferable from longitude
  if (is.null(tz_offset) && "lon" %in% names(dt) && any(!is.na(dt$lon))) {
    tz_offset <- round(stats::median(as.numeric(dt$lon), na.rm = TRUE) / 15)
  }
  if (!is.null(tz_offset) && is.finite(tz_offset) && tz_offset != 0) {
    dt[, Time := get(ts_col) + as.numeric(tz_offset) * 3600]
  } else {
    dt[, Time := get(ts_col)]
  }
  dt[, Time := floor_datetime(Time, step_sec = step_sec)]

  to_df(dt)
}

#' @rdname tt_extrapolate_dates
#' @export
extrapolate_tt_date <- tt_extrapolate_dates

