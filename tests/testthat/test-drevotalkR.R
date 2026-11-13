test_that("Ported algorithms (Savitzky-Golay, SunPosition, GapFill, Hampel) work and return pure data.frame/numeric", {
  x <- sin(seq(0, 2 * pi, length.out = 50))
  x_smooth <- sg_filter(x, m = 0, p = 3, w = 11)
  expect_type(x_smooth, "double")
  expect_length(x_smooth, 50)
  expect_lt(max(abs(x - x_smooth)), 0.05)

  dates <- as.POSIXct("2024-06-21 12:00:00", tz = "UTC") + (0:5) * 3600
  sun <- get_sun_position(dates, lat = 55.75, lon = 37.62)
  expect_identical(class(sun), "data.frame")
  expect_true(all(c("date", "altitude", "azimuth", "altitude_deg") %in% names(sun)))
  expect_gt(sun$altitude_deg[1], 30)

  gapped <- c(1, 2, NA, NA, 5, 6, rep(NA, 10), 17, 18)
  filled <- fill_small_gaps(gapped, max_gap = 3, method = "linear")
  expect_equal(filled[3:4], c(3, 4))
  expect_true(all(is.na(filled[7:16])))

  growth_raw <- c(seq(50, 48, length.out = 20), 120, seq(47.9, 47, length.out = 20))
  growth_clean <- clean_growth_series(growth_raw, k = 5, t0 = 3)
  expect_lt(growth_clean[21], 60)
})

test_that("Packet parsing, deduplication, and processing return pure base::data.frame objects", {
  # Build synthetic 48-hour dataset for 1 tree (Tilia) + 1 reference station (TTR) + 1 soil (55) + 1 TTG (61)
  base_ts <- as.numeric(as.POSIXct("2024-07-01 00:00:00", tz = "UTC")) + (0:47) * 3600
  srv_str <- format(as.POSIXct(base_ts, origin = "1970-01-01", tz = "UTC"), "%d.%m.%y %H:%M:%S")
  hr <- (0:47) %% 24
  diurnal <- sin((hr - 6) * pi / 12)

  # Tree 1 (52000001, Tilia, 4D + 49)
  tair_tree <- round((22 + 5 * diurnal) * 10)
  rh_tree   <- round((65 - 15 * diurnal) * 10)
  # Ref station (52000099, TTR, 4D + 49) - warmer and drier in daytime
  tair_ref  <- round((24 + 7 * diurnal) * 10)
  rh_ref    <- round((55 - 18 * diurnal) * 10)

  # Invert conv_analog_t_v3 (~ 16000 DN is ~31 C, 17500 DN is ~22 C)
  tref0_dn  <- rep(17500, 48)
  theat0_dn <- rep(17500, 48)
  tref1_dn  <- rep(17500, 48)
  # Lower dT during midday (higher sap flow), higher dT at night (~16200 vs 16800)
  theat1_dn <- round(16400 + 350 * pmax(0, diurnal))
  grow_dn   <- round(seq(20000, 19500, length.out = 48))

  lines_tree_4d <- sprintf(
    "%s,SN_%02d;52000001;4D;%d;%d;%d;%d;1050;0;%d;%d;16000;15;500;10;300;12;%d;%d;12500;2000",
    srv_str, 1:48, base_ts, tref0_dn, theat0_dn, grow_dn, rh_tree, tair_tree, tref1_dn, theat1_dn
  )
  lines_ref_4d <- sprintf(
    "%s,SN_%02d;52000099;4D;%d;%d;%d;%d;1050;0;%d;%d;16000;5;100;5;100;5;%d;%d;12500;2000",
    srv_str, 1:48, base_ts, tref0_dn, theat0_dn, grow_dn, rh_ref, tair_ref, tref1_dn, theat1_dn
  )
  lines_tree_49 <- sprintf(
    "%s,SN_%02d;52000001;49;%d;800;600;1500;2400;3200;3500;500;600;900;850;700;550;50;2",
    srv_str, 1:48, base_ts
  )
  lines_ref_49 <- sprintf(
    "%s,SN_%02d;52000099;49;%d;1600;1400;2800;3800;4800;5000;1100;1200;1600;1500;1400;1200;50;2",
    srv_str, 1:48, base_ts
  )
  lines_soil_55 <- sprintf(
    "%s,SN_%02d;52000055;55;%d;17500;17500;20000;1050;0;600;200;16000;5;100;5;100;5;17500;17500;14000;2000",
    srv_str[1:10], 1:10, base_ts[1:10]
  )
  lines_ttg_61 <- sprintf(
    "%s,SN_%02d;52000061;61;%d;1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17",
    srv_str[1:5], 1:5, base_ts[1:5]
  )

  raw_lines <- c(lines_tree_4d, lines_ref_4d, lines_tree_49, lines_ref_49, lines_soil_55, lines_ttg_61)
  raw_dt <- data.table::fread(text = raw_lines, sep = ";", header = FALSE, fill = TRUE)
  raw_dt[, V2 := as.character(V2)]

  meta <- data.frame(
    TT_ID = c("52000001", "52000099", "52000055"),
    Species = c("Tilia cordata", "TTR", "Soil"),
    DBH = c(75, 60, 60),
    Site = c("UrbanPark", "UrbanPark", "UrbanPark"),
    lat = c(55.75, 55.75, 55.75),
    lon = c(37.62, 37.62, 37.62),
    stringsAsFactors = FALSE
  )
  raw_merged <- merge(raw_dt, meta, by.x = "V2", by.y = "TT_ID", all.x = TRUE, sort = FALSE)
  v_cols <- paste0("V", 1:21)
  other_cols <- setdiff(names(raw_merged), v_cols)
  data.table::setcolorder(raw_merged, c(v_cols, other_cols))

  parsed <- tt_types_to_list(raw_merged, split_by_55 = TRUE)
  expect_type(parsed, "list")
  for (tbl in parsed) {
    expect_identical(class(tbl), "data.frame")
  }
  expect_equal(nrow(parsed$mydata_4D), 96L)
  expect_equal(nrow(parsed$mydata_49), 96L)
  expect_equal(nrow(parsed$mydata_55), 10L)
  expect_equal(nrow(parsed$mydata_61), 5L)

  # Full biophysical processing
  proc_df <- tt_process(parsed, grow_calc = TRUE, dtmax_method = "mw")
  expect_identical(class(proc_df), "data.frame")
  expect_true(all(c("Fdh", "Flux", "Flux_mm", "StWC", "growth", "NDVI", "Phi", "Osc_fct") %in% names(proc_df)))

  # Ecosystem services calculation
  es_df <- tt_calc_ecosystem_services(proc_df, ref_species = "TTR")
  expect_identical(class(es_df), "data.frame")
  expect_true(all(c("ES_dT", "ES_dRH", "ES_dVPD", "Alb_730", "Alb_760", "dBA_cm2", "dBiomass_kg", "dCarbon_kg", "CO2_eq_kg") %in% names(es_df)))

  # Verify Alb_730 and Alb_760 are both present and distinct (bug fix verification)
  tilia_rows <- es_df[es_df$Species == "Tilia cordata" & !is.na(es_df$Alb_730), ]
  expect_gt(nrow(tilia_rows), 0L)
  expect_false(identical(tilia_rows$Alb_730, tilia_rows$Alb_760))
  # Verify canopy cooling (Tilia is cooler than TTR)
  expect_lt(mean(tilia_rows$ES_dT, na.rm = TRUE), 0)

  # Summarise services
  sum_day <- tt_summarise_services(es_df, period = "day")
  sum_season <- tt_summarise_services(es_df, period = "season")
  expect_identical(class(sum_day), "data.frame")
  expect_identical(class(sum_season), "data.frame")
  expect_gt(nrow(sum_day), 0L)

  # Test plots in EN (default) and RU
  p_ts_en <- tt_plot_timeseries(es_df, var = "Fdh", lang = "en")
  p_ts_ru <- tt_plot_timeseries(es_df, var = "Fdh", lang = "ru")
  p_di_en <- tt_plot_diurnal(es_df, var = "ES_dT", lang = "en")
  p_es_ru <- tt_plot_ecoservices(es_df, service = "cooling", lang = "ru")
  p_wb_en <- tt_plot_waterbalance(es_df, unit = "L", lang = "en")
  p_sc_en <- tt_plot_scatter(es_df, x = "VPD", y = "Fdh", lang = "en")

  expect_s3_class(p_ts_en, "ggplot")
  expect_s3_class(p_ts_ru, "ggplot")
  expect_s3_class(p_di_en, "ggplot")
  expect_s3_class(p_es_ru, "ggplot")
  expect_s3_class(p_wb_en, "ggplot")
  expect_s3_class(p_sc_en, "ggplot")
  expect_match(p_ts_en$labels$y, "Sap Flux Density")
  expect_match(p_ts_ru$labels$y, "Плотность сокодвижения")
})

test_that("Logger terminal message cleaning, TT+ v2 (45) battery, and RTC date extrapolation work", {
  # Simulate a raw serial/UART download from an individual TreeTalker (like 219B0781.txt & 8123B048_backup.txt)
  raw_logger_dump <- c(
    ":ok",
    "TT+ rel. 4.0c",
    "wdt interval=1075 ms",
    "SPI",
    "I2c",
    "SN=219B0781",
    "last record @3BC000",
    "measure period=3600",
    "init done",
    ">",
    "Please insert:",
    "1 -> GO",
    "6 -> Read Memory",
    ">219B0781;1;45;1628247600;230;242;57754;41758;17;53;249;-4124;0;172;0;-257;0;230;290;10396",
    "219B0781;2;49;1628247600;3116;3199;4561;7919;12189;13061;10562;9352;4540;6669;5382;7518;50;3",
    "18.08.21 12:00:00,219B0781;3;45;1628251200;235;247;57433;41758;17;53;255;-4124;0;172;0;-256;0;235;296;10642",
    "219B0781;4;49;1628251200;2469;1229;4701;7108;10097;9929;7023;4662;4317;4370;4393;3016;50;3",
    "OK",
    "AT_COMMAND_NOT_FOUND",
    "RAKwireless RAK3172",
    "+EVT:JOINED",
    ">"
  )

  cleaned <- tt_clean_logger_lines(raw_logger_dump)
  expect_length(cleaned, 4L)
  expect_true(all(grepl("^[0-9]{2}\\.[0-9]{2}\\.[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2},", cleaned)))

  # Write to a temp file and test tt_read_backup + tt_types_to_list + tt_process
  tmp_dir <- file.path(tempdir(), "tt_bak_test")
  dir.create(tmp_dir, showWarnings = FALSE)
  tmp_file <- file.path(tmp_dir, "219B0781.txt")
  writeLines(raw_logger_dump, tmp_file)

  bak_df <- tt_read_backup(tmp_file, metadata = data.frame(TT_ID = "219B0781", Species = "Quercus robur", DBH = 80))
  expect_identical(class(bak_df), "data.frame")
  expect_equal(nrow(bak_df), 4L)
  expect_true(all(bak_df$V2 == "219B0781"))

  parsed_v2 <- tt_types_to_list(bak_df)
  expect_equal(nrow(parsed_v2$mydata_45), 2L)
  expect_equal(nrow(parsed_v2$mydata_49), 2L)
  expect_true(all(parsed_v2$mydata_45$TT_ID == "219B0781"))

  proc_v2 <- tt_process(parsed_v2, grow_calc = FALSE)
  expect_identical(class(proc_v2), "data.frame")
  # Verify Bug #1 fix: Bat_mV is in millivolts (> 3000 mV) and Volt is in Volts (> 3 V)
  expect_true(all(proc_v2$Bat_mV > 3000 & proc_v2$Bat_mV < 9000))
  expect_true(all(proc_v2$Volt > 3 & proc_v2$Volt < 9))

  # Test tt_extrapolate_dates with broken RTC timestamps (year 1970)
  ts_good <- as.POSIXct("2024-07-01 00:00:00", tz = "UTC") + (0:23) * 3600
  df_clock <- data.frame(
    TT_ID = rep(c("T_REF", "T_BROKEN"), each = 24),
    Timestamp = c(ts_good, ts_good[1:12], rep(as.POSIXct("1970-01-01 00:05:00", tz = "UTC"), 12)),
    Server_time = as.POSIXct(NA, tz = "UTC"),
    Tair = rep(20 + 5 * sin((0:23) * pi / 12), 2),
    adc_Vbat = 2000,
    stringsAsFactors = FALSE
  )
  repaired <- tt_extrapolate_dates(df_clock, step_sec = 3600)
  expect_identical(class(repaired), "data.frame")
  expect_true(all(repaired$Time >= as.POSIXct("2024-07-01 00:00:00", tz = "UTC")))
  unlink(tmp_dir, recursive = TRUE)
})
