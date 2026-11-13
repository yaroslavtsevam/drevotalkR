#' Plot TreeTalker biophysical or ecosystem service time series
#'
#' Generates a `ggplot2` time-series plot for any TreeTalker sensor channel or
#' ecosystem service metric, supporting combined (`"all_in_one"`) or faceted
#' (`"split"`) views and switchable English (`"en"`, default) or Russian (`"ru"`)
#' labels.
#'
#' @param data A processed `data.frame` from [tt_process()] or
#'   [tt_calc_ecosystem_services()].
#' @param var Name of the numeric column to plot (default `"Fdh"`).
#' @param group_by Grouping column for color/facets (default `"TT_ID"`).
#' @param mode Plot layout: `"all_in_one"` (single panel) or `"split"`
#'   (faceted by `group_by`).
#' @param smooth Logical; overlay a smoothed trend line (default `FALSE`).
#' @param hline Optional numeric value for a horizontal reference line (e.g.,
#'   `3500` mV for battery threshold or `0` for `ES_dT`).
#' @param lang Language for axis and title labels: `"en"` (default) or `"ru"`.
#' @return A `ggplot` object.
#' @export
tt_plot_timeseries <- function(data, var = "Fdh", group_by = "TT_ID",
                               mode = c("all_in_one", "split"),
                               smooth = FALSE, hline = NULL,
                               lang = c("en", "ru")) {
  mode <- match.arg(mode)
  lang <- match.arg(lang)
  if (is.list(data) && !is.data.frame(data)) data <- tt_process(data)
  if (is.null(data) || nrow(data) == 0L || !var %in% names(data)) {
    stop(sprintf("Column '%s' not found in data.", var))
  }

  dt <- data.table::as.data.table(data)
  ts_col <- intersect(c("Time", "Datetime", "Timestamp", "Date"), names(dt))[1L]
  if (is.na(ts_col)) stop("No timestamp column ('Time' or 'Timestamp') found in data.")
  if (!group_by %in% names(dt)) dt[, (group_by) := "TreeTalker"]

  dt_plot <- dt[!is.na(get(ts_col)) & !is.na(get(var)),
                .(Time = get(ts_col),
                  val = as.numeric(get(var)),
                  group_var = as.character(get(group_by)))]
  df_plot <- to_df(dt_plot)

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Time, y = val, color = group_var)) +
    ggplot2::geom_line(alpha = if (smooth) 0.45 else 0.85, linewidth = 0.6) +
    ggplot2::labs(
      title = paste0(tt_i18n("title_timeseries", lang), ": ", tt_i18n(var, lang)),
      x = tt_i18n("time", lang),
      y = tt_i18n(var, lang),
      color = tt_i18n(tolower(group_by), lang)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")

  if (isTRUE(smooth) && nrow(df_plot) >= 10L) {
    p <- p + ggplot2::geom_smooth(se = FALSE, linewidth = 1.0, method = "loess", formula = y ~ x)
  }
  if (!is.null(hline) && is.finite(hline)) {
    p <- p + ggplot2::geom_hline(yintercept = hline, linetype = "dashed", color = "firebrick", linewidth = 0.6)
  }
  if (mode == "split") {
    p <- p + ggplot2::facet_wrap(~ group_var, scales = "free_y")
  }
  p
}

#' Plot diurnal cycle (Mean +/- SE) of biophysical or ecosystem service variables
#'
#' Computes hourly diurnal profiles (`0..23` h) with mean and standard error
#' (\eqn{\pm \text{SE}}) ribbons using `data.table` and renders a publication-ready
#' `ggplot2` chart without requiring the `openair` package.
#'
#' @param data A processed `data.frame` from [tt_process()] or
#'   [tt_calc_ecosystem_services()].
#' @param var Numeric column name to plot (default `"ES_dT"`).
#' @param group_by Grouping variable for color/fill (default `"Species"`).
#' @param facet_by Optional column name for faceting (e.g., `"Month"` or `"Site"`).
#' @param lang Language for labels: `"en"` (default) or `"ru"`.
#' @return A `ggplot` object.
#' @export
tt_plot_diurnal <- function(data, var = "ES_dT", group_by = "Species",
                            facet_by = NULL, lang = c("en", "ru")) {
  lang <- match.arg(lang)
  if (is.list(data) && !is.data.frame(data)) data <- tt_calc_ecosystem_services(data)
  if (is.null(data) || nrow(data) == 0L) stop("Input data is empty.")
  if (!var %in% names(data)) {
    data <- tt_calc_ecosystem_services(data)
    if (!var %in% names(data)) stop(sprintf("Column '%s' not found in data.", var))
  }

  dt <- data.table::as.data.table(data)
  if (!"Hour" %in% names(dt)) {
    ts_col <- intersect(c("Time", "Datetime", "Timestamp"), names(dt))[1L]
    dt[, Hour := as.POSIXlt(get(ts_col), tz = "UTC")$hour]
  }
  if (!group_by %in% names(dt)) {
    group_by <- if ("TT_ID" %in% names(dt)) "TT_ID" else "Group"
    if (!group_by %in% names(dt)) dt[, Group := "All"]
  }

  if (!is.null(facet_by) && facet_by %in% names(dt)) {
    dt_sum <- dt[!is.na(get(var)) & !is.na(Hour), .(
      mean_val = mean(as.numeric(get(var)), na.rm = TRUE),
      sd_val   = stats::sd(as.numeric(get(var)), na.rm = TRUE),
      n_obs    = sum(!is.na(get(var)))
    ), by = .(Hour, group_var = as.character(get(group_by)), facet_var = as.character(get(facet_by)))]
  } else {
    dt_sum <- dt[!is.na(get(var)) & !is.na(Hour), .(
      mean_val = mean(as.numeric(get(var)), na.rm = TRUE),
      sd_val   = stats::sd(as.numeric(get(var)), na.rm = TRUE),
      n_obs    = sum(!is.na(get(var)))
    ), by = .(Hour, group_var = as.character(get(group_by)))]
  }
  dt_sum[, se_val := ifelse(is.na(sd_val) | n_obs <= 1L, 0, sd_val / sqrt(n_obs))]
  df_sum <- to_df(dt_sum)

  p <- ggplot2::ggplot(df_sum, ggplot2::aes(x = Hour, y = mean_val, color = group_var, fill = group_var)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                         alpha = 0.2, color = NA) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_x_continuous(breaks = seq(0, 23, by = 3)) +
    ggplot2::labs(
      title = paste0(tt_i18n("title_diurnal", lang), ": ", tt_i18n(var, lang)),
      x = tt_i18n("hour", lang),
      y = tt_i18n(var, lang),
      color = tt_i18n(tolower(group_by), lang),
      fill = tt_i18n(tolower(group_by), lang)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")

  if (var %in% c("ES_dT", "ES_dRH", "ES_dVPD", "phi_d")) {
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40")
  }
  if (!is.null(facet_by) && "facet_var" %in% names(df_sum)) {
    p <- p + ggplot2::facet_wrap(~ facet_var)
  }
  p
}

#' Plot comparative Ecosystem Service indicators across species or sites
#'
#' Creates comparative boxplot/jitter charts for the four primary TreeTalker
#' ecosystem services: canopy cooling (`"cooling"`), water regulation
#' (`"transpiration"`), carbon sequestration (`"carbon"`), and biomechanical
#' wind stability (`"stability"`).
#'
#' @param data A `data.frame` from [tt_calc_ecosystem_services()] or
#'   [tt_process()].
#' @param service Ecosystem service type: `"cooling"` (default),
#'   `"transpiration"`, `"carbon"`, or `"stability"`.
#' @param group_by Grouping column on the x-axis (default `"Species"`).
#' @param hot_only Logical; for `"cooling"`, restrict to hot hours (`Tair_Ref >= 25` C,
#'   default `FALSE`).
#' @param lang Language for labels: `"en"` (default) or `"ru"`.
#' @return A `ggplot` object.
#' @export
tt_plot_ecoservices <- function(data,
                                service = c("cooling", "transpiration", "carbon", "stability"),
                                group_by = "Species", hot_only = FALSE,
                                lang = c("en", "ru")) {
  service <- match.arg(service)
  lang <- match.arg(lang)
  if (!"ES_dT" %in% names(data) && !"dCarbon_kg" %in% names(data)) {
    data <- tt_calc_ecosystem_services(data)
  }
  dt <- data.table::as.data.table(data)
  if (!group_by %in% names(dt)) {
    group_by <- if ("TT_ID" %in% names(dt)) "TT_ID" else "Group"
    if (!group_by %in% names(dt)) dt[, Group := "TreeTalker"]
  }

  if (service == "cooling") {
    if (isTRUE(hot_only) && "Tair_Ref" %in% names(dt)) {
      dt_sub <- dt[!is.na(Tair_Ref) & Tair_Ref >= 25]
      if (nrow(dt_sub) > 0L) dt <- dt_sub
    }
    var_col <- "ES_dT"
    title_key <- "title_cooling"
    df_plot <- to_df(dt[!is.na(get(var_col)), .(group_var = as.character(get(group_by)),
                                                val = as.numeric(get(var_col)))])
  } else if (service == "transpiration") {
    # Daily transpiration per tree (L/day)
    sum_df <- tt_summarise_services(to_df(dt), period = "day", by = unique(c("TT_ID", group_by)))
    dt_sum <- data.table::as.data.table(sum_df)
    var_col <- "Flux"
    title_key <- "title_transpiration"
    df_plot <- to_df(dt_sum[!is.na(Flux_L_sum), .(group_var = as.character(get(group_by)),
                                                  val = as.numeric(Flux_L_sum))])
  } else if (service == "carbon") {
    sum_df <- tt_summarise_services(to_df(dt), period = "season", by = unique(c("TT_ID", group_by)))
    dt_sum <- data.table::as.data.table(sum_df)
    var_col <- "dCarbon_kg"
    title_key <- "title_carbon"
    df_plot <- to_df(dt_sum[!is.na(dCarbon_kg), .(group_var = as.character(get(group_by)),
                                                  val = as.numeric(dCarbon_kg))])
  } else {
    var_col <- "Osc_fct"
    title_key <- "title_stability"
    df_plot <- to_df(dt[!is.na(get(var_col)), .(group_var = as.character(get(group_by)),
                                                val = as.numeric(get(var_col)))])
  }

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = group_var, y = val, fill = group_var)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1.2) +
    ggplot2::labs(
      title = tt_i18n(title_key, lang),
      x = tt_i18n(tolower(group_by), lang),
      y = tt_i18n(var_col, lang),
      fill = tt_i18n(tolower(group_by), lang)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "none")

  if (service == "cooling") {
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "firebrick")
  }
  p
}

#' Plot cumulative tree water balance (transpiration) over time
#'
#' Computes and plots cumulative whole-tree (`unit = "L"`) or crown-projected
#' (`unit = "mm"`) transpiration trajectories per tree or species.
#'
#' @param data A processed `data.frame` containing `Time`, `TT_ID`, and `Flux`
#'   (or `Flux_mm`).
#' @param group_by Grouping column for color (default `"Species"`).
#' @param unit Transpiration unit: `"L"` (liters per tree, default) or `"mm"`
#'   (millimeters per crown projection area).
#' @param lang Language for labels: `"en"` (default) or `"ru"`.
#' @return A `ggplot` object.
#' @export
tt_plot_waterbalance <- function(data, group_by = "Species",
                                 unit = c("L", "mm"), lang = c("en", "ru")) {
  unit <- match.arg(unit)
  lang <- match.arg(lang)
  if (is.list(data) && !is.data.frame(data)) data <- tt_process(data)
  dt <- data.table::as.data.table(data)

  flux_col <- if (unit == "L") "Flux" else "Flux_mm"
  y_key <- if (unit == "L") "cum_flux" else "cum_flux_mm"
  if (!flux_col %in% names(dt)) stop(sprintf("Column '%s' not found in data.", flux_col))
  if (!group_by %in% names(dt)) {
    group_by <- if ("TT_ID" %in% names(dt)) "TT_ID" else "Group"
    if (!group_by %in% names(dt)) dt[, Group := "TreeTalker"]
  }

  ts_col <- intersect(c("Time", "Datetime", "Timestamp"), names(dt))[1L]
  data.table::setorderv(dt, c("TT_ID", ts_col))
  dt <- unique(dt, by = c("TT_ID", ts_col))

  # Compute mean hourly flux per group and cumulative sum scaled by time step (hours)
  dt_grp <- dt[!is.na(get(ts_col)) & !is.na(get(flux_col)), .(
    mean_flux = mean(as.numeric(get(flux_col)), na.rm = TRUE)
  ), by = .(Time = get(ts_col), group_var = as.character(get(group_by)))]
  data.table::setorder(dt_grp, group_var, Time)
  dt_grp[, dt_step_h := {
    d_h <- diff(as.numeric(Time)) / 3600
    d_pos <- d_h[d_h > 0 & d_h <= 6]
    if (length(d_pos) > 0L) stats::median(d_pos, na.rm = TRUE) else 1
  }, by = group_var]
  dt_grp[, cum_val := cumsum(ifelse(is.na(mean_flux), 0, mean_flux * dt_step_h)), by = group_var]

  df_plot <- to_df(dt_grp)
  ggplot2::ggplot(df_plot, ggplot2::aes(x = Time, y = cum_val, color = group_var)) +
    ggplot2::geom_line(linewidth = 1.0) +
    ggplot2::labs(
      title = tt_i18n("title_waterbalance", lang),
      x = tt_i18n("time", lang),
      y = tt_i18n(y_key, lang),
      color = tt_i18n(tolower(group_by), lang)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

#' Plot eco-physiological scatter relationships with regression fit
#'
#' Plots bivariate relationships between TreeTalker variables (e.g., `VPD` vs
#' `Fdh`, `NDVI` vs `Fdh`, `Tair_Ref` vs `ES_dT`) with linear regression lines
#' and overall \eqn{R^2} in the subtitle.
#'
#' @param data A processed `data.frame`.
#' @param x Name of x-axis variable (default `"VPD"`).
#' @param y Name of y-axis variable (default `"Fdh"`).
#' @param group_by Grouping column for color (default `"Species"`).
#' @param lang Language for labels: `"en"` (default) or `"ru"`.
#' @return A `ggplot` object.
#' @export
tt_plot_scatter <- function(data, x = "VPD", y = "Fdh", group_by = "Species",
                            lang = c("en", "ru")) {
  lang <- match.arg(lang)
  if (is.list(data) && !is.data.frame(data)) data <- tt_process(data)
  dt <- data.table::as.data.table(data)
  if (!all(c(x, y) %in% names(dt))) {
    stop(sprintf("Columns '%s' and/or '%s' not found in data.", x, y))
  }
  if (!group_by %in% names(dt)) {
    group_by <- if ("TT_ID" %in% names(dt)) "TT_ID" else "Group"
    if (!group_by %in% names(dt)) dt[, Group := "TreeTalker"]
  }

  dt_sub <- dt[!is.na(get(x)) & !is.na(get(y)),
               .(x_val = as.numeric(get(x)),
                 y_val = as.numeric(get(y)),
                 group_var = as.character(get(group_by)))]
  df_plot <- to_df(dt_sub)

  r2_val <- if (nrow(df_plot) >= 3L && stats::sd(df_plot$x_val) > 0 && stats::sd(df_plot$y_val) > 0) {
    round(summary(stats::lm(y_val ~ x_val, data = df_plot))$r.squared, 3)
  } else {
    NA_real_
  }

  ggplot2::ggplot(df_plot, ggplot2::aes(x = x_val, y = y_val, color = group_var)) +
    ggplot2::geom_point(alpha = 0.45, size = 1.5) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.9) +
    ggplot2::labs(
      title = paste0(tt_i18n(y, lang), " vs. ", tt_i18n(x, lang)),
      subtitle = if (!is.na(r2_val)) paste0("R\u00B2 = ", r2_val) else NULL,
      x = tt_i18n(x, lang),
      y = tt_i18n(y, lang),
      color = tt_i18n(tolower(group_by), lang)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

#' Classic TreeTalker diagnostic plot wrappers
#'
#' Backward-compatible plotting functions (`ttbatteryplot`, `ttgranierplot`,
#' `ttStWCplot`, `ttgrowthplot`, `ttlightplot`, `ttstabilityplot`) with
#' switchable English (`lang = "en"`, default) and Russian (`lang = "ru"`) labels.
#'
#' @param mydata A processed `data.frame` or list from [tt_types_to_list()].
#' @param mode Layout mode: `"all_in_one"` or `"split"`.
#' @param lang Language code: `"en"` (default) or `"ru"`.
#' @return A `ggplot` object.
#' @name tt_legacy_plots
NULL

#' @rdname tt_legacy_plots
#' @export
ttbatteryplot <- function(mydata, mode = c("all_in_one", "split"), lang = c("en", "ru")) {
  df <- ttbattery(mydata)
  tt_plot_timeseries(df, var = "Bat_mV", group_by = "TT_ID", mode = mode, hline = 3500, lang = lang)
}

#' @rdname tt_legacy_plots
#' @export
ttgranierplot <- function(mydata, mode = c("all_in_one", "split"), lang = c("en", "ru")) {
  df <- ttgranier(mydata)
  tt_plot_timeseries(df, var = "Fdh", group_by = "TT_ID", mode = mode, lang = lang)
}

#' @rdname tt_legacy_plots
#' @export
ttStWCplot <- function(mydata, mode = c("all_in_one", "split"), lang = c("en", "ru")) {
  df <- ttStWC(mydata)
  tt_plot_timeseries(df, var = "StWC", group_by = "TT_ID", mode = mode, lang = lang)
}

#' @rdname tt_legacy_plots
#' @export
ttgrowthplot <- function(mydata, mode = c("all_in_one", "split"), lang = c("en", "ru")) {
  df <- ttgrowth(mydata)
  tt_plot_timeseries(df, var = "growth", group_by = "TT_ID", mode = mode, smooth = TRUE, lang = lang)
}

#' @rdname tt_legacy_plots
#' @export
ttlightplot <- function(mydata, mode = c("all_in_one", "split"), lang = c("en", "ru")) {
  df <- ttlight(mydata)
  tt_plot_timeseries(df, var = "NDVI", group_by = "TT_ID", mode = mode, lang = lang)
}

#' @rdname tt_legacy_plots
#' @export
ttstabilityplot <- function(mydata, mode = c("all_in_one", "split"), lang = c("en", "ru")) {
  df <- ttstability(mydata)
  tt_plot_timeseries(df, var = "Phi", group_by = "TT_ID", mode = mode, lang = lang)
}

