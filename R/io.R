#' Clean raw TreeTalker download lines from logger and terminal messages
#'
#' Filters raw text lines extracted from individual TreeTalker loggers (via
#' serial/UART cable or LoRa terminal) or TT-Cloud backups, removing firmware
#' boot headers (`TT+ rel. 4.0c`, `SPI`, `I2c`, `SN=...`, `init done`, `>`),
#' interactive menu prompts (`1 -> GO`, `6 -> Read Memory`), and AT/LoRaWAN
#' status messages (`OK`, `AT_COMMAND_NOT_FOUND`, `RAKwireless`, `+EVT:JOINED`).
#' Also normalizes missing server timestamp prefixes per line.
#'
#' @param lines Character vector of raw lines read from a TreeTalker file.
#' @param add_server_time Logical; if `TRUE` (default), prepends a placeholder
#'   server timestamp `"01.01.70 00:00:00,"` to any valid data line that lacks
#'   a `"dd.mm.yy hh:mm:ss,"` prefix.
#' @return Character vector of cleaned, semicolon-delimited telemetry lines.
#' @export
tt_clean_logger_lines <- function(lines, add_server_time = TRUE) {
  if (is.null(lines) || length(lines) == 0L) return(character(0))
  lines <- iconv(as.character(lines), from = "latin1", to = "UTF-8", sub = "")
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0L) return(character(0))

  # Strip leading terminal prompts (e.g. ":ok", ">") glued to the start of a data line
  lines <- sub("^(:ok\\s*|[:>\\s]+)", "", lines, ignore.case = TRUE)
  lines <- trimws(lines)

  # Keep only lines where the 3rd semicolon-delimited field is a valid TT packet code
  valid_pkt <- grepl("^([^;]+;){2}\\s*(45|49|4B|4C|4D|55|60|61)\\s*;", lines, ignore.case = TRUE)
  lines <- lines[valid_pkt]
  if (length(lines) == 0L) return(character(0))

  if (isTRUE(add_server_time)) {
    has_srv <- grepl("^[0-9]{2}\\.[0-9]{2}\\.[0-9]{2}\\s+[0-9]{2}:[0-9]{2}:[0-9]{2}\\s*,", lines)
    if (any(!has_srv)) {
      lines[!has_srv] <- paste0("01.01.70 00:00:00,", lines[!has_srv])
    }
  }
  lines
}

#' Normalize V1/V2 order (Rec_Nr vs TT_ID) and merge device metadata
#'
#' Internal helper that detects whether `TT_ID` is in `V1` (standard `ttcloud.txt`
#' and direct logger dumps like `219B0781;1;45;...`) or `V2` (synthetic tables
#' `SN_01;52000001;4D;...`), normalizes `V2` to `TT_ID`, and merges `metadata`.
#'
#' @param dt A `data.table` of raw records.
#' @param metadata Optional `data.frame`, `data.table`, or file path to metadata.
#' @return A `data.table` with normalized `V1` (`[Server_time,]Rec_Nr`), `V2` (`TT_ID`),
#'   and merged metadata columns placed after `V1..Vn`.
#' @keywords internal
merge_tt_metadata <- function(dt, metadata = NULL) {
  if (is.null(dt) || nrow(dt) == 0L) return(dt)
  dt <- data.table::as.data.table(dt)
  if (!all(c("V1", "V2") %in% names(dt))) return(dt)

  v1_chr <- as.character(dt[["V1"]])
  has_comma <- grepl(",", v1_chr, fixed = TRUE)
  srv_part <- ifelse(has_comma, sub(",.*$", "", v1_chr), NA_character_)
  tok1 <- trimws(ifelse(has_comma, sub("^.*,", "", v1_chr), v1_chr))
  tok2 <- trimws(as.character(dt[["V2"]]))

  # In direct TT dumps (e.g. "219B0781;1;45;...") and ttcloud ("dd.mm.yy hh:mm:ss,219B0781;1;45;..."),
  # tok1 is the 8-char device serial (TT_ID) and tok2 is the hex record counter (Rec_Nr).
  # In pre-swapped tables ("dd.mm.yy hh:mm:ss,SN_01;52000001;4D;..."), tok1 starts with SN_ or tok2 is the 8-char ID.
  tok1_is_id <- all(!grepl("^SN_", head(tok1, 20L), ignore.case = TRUE)) &&
    all(nchar(head(tok1[nzchar(tok1)], 20L)) >= 6L) &&
    (any(nchar(head(tok2[nzchar(tok2)], 20L)) <= 5L) ||
       length(unique(head(tok1, 100L))) < length(unique(head(tok2, 100L))))

  if (isTRUE(tok1_is_id)) {
    dt[, `:=`(
      V1 = ifelse(has_comma, paste0(srv_part, ",", tok2), tok2),
      V2 = tok1
    )]
  } else {
    dt[, V2 := tok2]
  }

  if (!is.null(metadata)) {
    meta_dt <- if (is.character(metadata) && length(metadata) == 1L && file.exists(metadata)) {
      data.table::fread(metadata, showProgress = FALSE)
    } else {
      data.table::as.data.table(metadata)
    }
    if ("id" %in% names(meta_dt) && !"TT_ID" %in% names(meta_dt)) {
      data.table::setnames(meta_dt, "id", "TT_ID")
    }
    if ("TT_ID" %in% names(meta_dt)) {
      meta_dt[, TT_ID := trimws(as.character(TT_ID))]
      meta_dt <- unique(meta_dt, by = "TT_ID")
      # Remove any already-present metadata columns in dt before merging
      dup_meta <- intersect(setdiff(names(meta_dt), "TT_ID"), names(dt))
      if (length(dup_meta) > 0L) dt[, (dup_meta) := NULL]
      dt <- merge(dt, meta_dt, by.x = "V2", by.y = "TT_ID", all.x = TRUE, sort = FALSE)
      v_cols <- grep("^V[0-9]+$", names(dt), value = TRUE)
      v_cols <- v_cols[order(as.integer(sub("^V", "", v_cols)))]
      other_cols <- setdiff(names(dt), v_cols)
      data.table::setcolorder(dt, c(v_cols, other_cols))
    }
  }
  dt
}

#' Parse raw TreeTalker records into a list of typed data.frames
#'
#' Splits raw TreeTalker telemetry records by hexadecimal packet identifier
#' (`45`, `49`, `4B`, `4C`, `4D`, `55`, `60`, `61`), assigns standard column
#' names, converts timestamps and numeric sensor channels using `data.table`,
#' and returns a named list of pure `base::data.frame` tables.
#'
#' @param ttdata A `data.frame` or `data.table` of raw TreeTalker records (at
#'   least 4 columns: server/record header, `TT_ID`, packet type `Dev_Ty`,
#'   device `Timestamp`, followed by sensor fields).
#' @param split_by_55 Logical; if `TRUE`, retains `mydata_55` (TT-Soil v3) as a
#'   separate table (default `TRUE`). If `FALSE`, merges `55` packets with `4D`.
#' @return A named list of pure `data.frame` objects: `mydata_45`, `mydata_49`,
#'   `mydata_4B`, `mydata_4C`, `mydata_4D`, `mydata_55`, `mydata_60`,
#'   `mydata_61`.
#' @export
tt_types_to_list <- function(ttdata, split_by_55 = TRUE) {
  if (is.null(ttdata) || nrow(ttdata) == 0L) {
    empty_df <- data.frame()
    return(list(
      mydata_45 = empty_df, mydata_49 = empty_df, mydata_4B = empty_df,
      mydata_4C = empty_df, mydata_4D = empty_df, mydata_55 = empty_df,
      mydata_60 = empty_df, mydata_61 = empty_df
    ))
  }

  dt <- data.table::as.data.table(ttdata)
  known_meta <- c("Species", "DBH", "Name", "Site", "SiteIndex", "Corniferous",
                  "DTT", "installation_height", "lat", "lon", "Dist", "Height",
                  "CloudIndex", "TTver", "CloudName")
  v_raw <- grep("^V[0-9]+$", names(dt), value = TRUE)
  if (length(v_raw) >= 4L) {
    meta_cols <- setdiff(names(dt), v_raw)
  } else {
    meta_cols <- intersect(names(dt), known_meta)
    raw_cols <- setdiff(names(dt), meta_cols)
    std_v <- paste0("V", seq_along(raw_cols))
    data.table::setnames(dt, raw_cols, std_v)
  }

  # Normalize V1 (Rec_Nr) and V2 (TT_ID) order
  dt <- merge_tt_metadata(dt, metadata = NULL)

  # Parse server timestamp from V1 (format: "dd.mm.yy hh:mm:ss,Rec_Nr")
  v1_chr <- as.character(dt[["V1"]])
  has_comma <- grepl(",", v1_chr, fixed = TRUE)
  srv_str <- ifelse(has_comma, sub(",.*$", "", v1_chr), NA_character_)
  rec_str <- ifelse(has_comma, sub("^.*,", "", v1_chr), v1_chr)
  srv_time <- suppressWarnings(as.POSIXct(srv_str, format = "%d.%m.%y %H:%M:%S", tz = "UTC"))

  dt[, `:=`(
    Server_time = srv_time,
    V1 = rec_str,
    V2 = as.character(V2),
    V3 = toupper(trimws(as.character(V3)))
  )]

  # Helper to extract and format a specific packet type
  extract_packet <- function(sub_dt, col_names, num_indices, ts_col = "Timestamp") {
    if (is.null(sub_dt) || nrow(sub_dt) == 0L) {
      res <- as.data.frame(setNames(
        replicate(length(col_names), logical(0), simplify = FALSE),
        col_names
      ))
      return(res)
    }
    sub_dt <- data.table::copy(sub_dt)
    req_v <- paste0("V", seq_along(col_names))
    missing_v <- setdiff(req_v, names(sub_dt))
    for (mv in missing_v) sub_dt[, (mv) := NA]

    keep_cols <- c("Server_time", req_v, intersect(meta_cols, names(sub_dt)))
    sub_dt <- sub_dt[, ..keep_cols]
    data.table::setnames(sub_dt, req_v, col_names)

    for (idx in num_indices) {
      cname <- col_names[idx]
      if (cname %in% names(sub_dt)) {
        sub_dt[, (cname) := suppressWarnings(as.numeric(as.character(get(cname))))]
      }
    }
    if (ts_col %in% names(sub_dt)) {
      ts_num <- suppressWarnings(as.numeric(sub_dt[[ts_col]]))
      sub_dt[, (ts_col) := as.POSIXct(ts_num, origin = "1970-01-01", tz = "UTC")]
      # Fallback to Server_time when RTC Timestamp is missing or before 2015
      bad_ts <- is.na(sub_dt[[ts_col]]) | sub_dt[[ts_col]] < as.POSIXct("2015-01-01", tz = "UTC")
      if (any(bad_ts) && "Server_time" %in% names(sub_dt)) {
        sub_dt[bad_ts & !is.na(Server_time) & Server_time >= as.POSIXct("2015-01-01", tz = "UTC"),
               (ts_col) := Server_time]
      }
    }
    to_df(sub_dt)
  }

  # 1. Packet 49: Spectral sensor (AS7262 VIS + AS7263 NIR)
  cols_49 <- c(
    "Rec_Nr", "TT_ID", "Dev_Ty", "Timestamp",
    "AS7263_610", "AS7263_680", "AS7263_730", "AS7263_760", "AS7263_810", "AS7263_860",
    "AS7262_450", "AS7262_500", "AS7262_550", "AS7262_570", "AS7262_600", "AS7262_650",
    "integration_T", "gain"
  )
  mydata_49 <- extract_packet(dt[V3 == "49"], cols_49, 4:18)

  # 2. Packet 4D (TT+ v3) & Packet 55 (TT-Soil v3)
  cols_4d <- c(
    "Rec_Nr", "TT_ID", "Dev_Ty", "Timestamp",
    "Tref_0", "Theat_0", "growt_sens", "adc_bandgap", "Bits",
    "RH", "Tair", "gz_mean", "gz_sd", "gy_mean", "gy_sd", "gx_mean", "gx_sd",
    "Tref_1", "Theat_1", "MoistSensHz", "adc_Vbat"
  )
  if (isTRUE(split_by_55)) {
    mydata_4D <- extract_packet(dt[V3 == "4D"], cols_4d, 4:21)
    mydata_55 <- extract_packet(dt[V3 == "55"], cols_4d, 4:21)
  } else {
    mydata_4D <- extract_packet(dt[V3 %in% c("4D", "55")], cols_4d, 4:21)
    mydata_55 <- extract_packet(dt[0L], cols_4d, 4:21)
  }

  # 3. Packet 45: TT+ v2
  cols_45 <- c(
    "Rec_Nr", "TT_ID", "Dev_Ty", "Timestamp",
    "Tref_0", "Theat_0", "growt_sens", "adc_bandgap", "Bits",
    "RH", "Tair", "gz_mean", "gz_sd", "gy_mean", "gy_sd", "gx_mean", "gx_sd",
    "Tref_1", "Theat_1", "MoistSensHz"
  )
  mydata_45 <- extract_packet(dt[V3 == "45"], cols_45, 4:20)

  # 4. Packet 4B: TT-Cloud status
  cols_4b <- c(
    "Rec_Nr", "TT_ID", "Dev_Ty", "Timestamp",
    "Acc_Rec", "Sent_Rec", "Pag_Op", "Bandgap",
    "Bat_mV", "GSM_field", "Nev_ce", "N_err", "Not_Usd"
  )
  mydata_4B <- extract_packet(dt[V3 == "4B"], cols_4b, 4:13)

  # 5. Packet 4C: TT-Cloud RSSI table
  dt_4c <- dt[V3 == "4C"]
  mydata_4C <- if (nrow(dt_4c) > 0L) to_df(dt_4c) else data.frame()

  # 6. Packets 60 & 61: TTG clinometer/accelerometer
  cols_60 <- paste0("C", 1:21)
  cols_60[1:4] <- c("Rec_Nr", "TT_ID", "Dev_Ty", "Timestamp")
  mydata_60 <- extract_packet(dt[V3 == "60"], cols_60, 4:21)
  mydata_61 <- extract_packet(dt[V3 == "61"], cols_60, 4:21)

  list(
    mydata_45 = mydata_45,
    mydata_49 = mydata_49,
    mydata_4B = mydata_4B,
    mydata_4C = mydata_4C,
    mydata_4D = mydata_4D,
    mydata_55 = mydata_55,
    mydata_60 = mydata_60,
    mydata_61 = mydata_61
  )
}

#' @rdname tt_types_to_list
#' @export
tttypes_to_list <- tt_types_to_list

#' Remove duplicate TreeTalker telemetry records
#'
#' Fast deduplication of raw or parsed TreeTalker tables using `data.table`,
#' returning pure `data.frame` objects.
#'
#' @param data A `data.frame`, `data.table`, or named `list` of tables.
#' @return Deduplicated `data.frame` (or list of `data.frame`s).
#' @export
tt_remove_duplicates <- function(data) {
  dedup_one <- function(df) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(to_df(df))
    dt <- data.table::as.data.table(df)
    key_cols <- intersect(c("TT_ID", "Dev_Ty", "Timestamp", "V2", "V3", "V4"), names(dt))
    if (length(key_cols) >= 2L) {
      dt <- unique(dt, by = key_cols)
    } else {
      dt <- unique(dt)
    }
    to_df(dt)
  }
  if (is.list(data) && !is.data.frame(data)) {
    return(lapply(data, dedup_one))
  }
  dedup_one(data)
}

#' @rdname tt_remove_duplicates
#' @export
remove_duplicates <- tt_remove_duplicates

#' Read raw TreeTalker SD-card or serial logger backup files
#'
#' Reads semicolon-delimited raw telemetry files from a backup folder, cleans
#' out any firmware boot banners, serial menu prompts, or AT/LoRaWAN messages
#' via [tt_clean_logger_lines()], normalizes `TT_ID` and `Rec_Nr` columns,
#' optionally merges device metadata (`TT_ID`, `Species`, `DBH`, `Site`), and
#' returns a pure `data.frame`.
#'
#' @param backup_folder Path to directory containing raw `.txt`, `.csv`, or
#'   `.dat` files, or a character vector of file paths.
#' @param metadata Optional `data.frame` or path to a tab/semicolon-separated
#'   metadata file containing at least column `TT_ID` (or `id`).
#' @return A pure `data.frame` of raw TreeTalker records ready for
#'   [tt_types_to_list()].
#' @export
tt_read_backup <- function(backup_folder, metadata = NULL) {
  if (length(backup_folder) == 1L && dir.exists(backup_folder)) {
    files <- list.files(backup_folder, pattern = "\\.(txt|csv|dat)$",
                        full.names = TRUE, ignore.case = TRUE)
    files <- files[!grepl("tt_Descriptive|metadata|serial_numbers", basename(files), ignore.case = TRUE)]
  } else {
    files <- backup_folder[file.exists(backup_folder)]
  }
  if (length(files) == 0L) return(data.frame())

  dt_list <- lapply(files, function(f) {
    raw_lines <- readLines(f, warn = FALSE)
    clean_lines <- tt_clean_logger_lines(raw_lines, add_server_time = TRUE)
    if (length(clean_lines) == 0L) return(NULL)
    tryCatch(
      data.table::fread(text = clean_lines, sep = ";", fill = TRUE, header = FALSE, showProgress = FALSE),
      error = function(e) NULL
    )
  })
  dt <- data.table::rbindlist(dt_list, fill = TRUE, use.names = TRUE)
  if (nrow(dt) == 0L) return(data.frame())

  dt <- merge_tt_metadata(dt, metadata = metadata)
  to_df(dt)
}

#' @rdname tt_read_backup
#' @export
tt_backup_import <- tt_read_backup

#' Restore missing server timestamps in raw cloud/SD backup logs
#'
#' Strips any logger/terminal messages using [tt_clean_logger_lines()],
#' reconstructs the `"dd.mm.yy HH:MM:SS,..."` server header in raw backup
#' files from the internal device epoch timestamp (`V4`), and optionally writes
#' the restored file to disk.
#'
#' @param input_data Path to a raw backup file, a character vector of raw lines,
#'   or a raw `data.frame`.
#' @param output_file Optional file path to write the restored semicolon-separated
#'   backup text.
#' @return A pure `data.frame` of restored raw records.
#' @export
tt_restore_cloud_backup <- function(input_data, output_file = NULL) {
  if (is.character(input_data) && length(input_data) == 1L && file.exists(input_data)) {
    clean_lines <- tt_clean_logger_lines(readLines(input_data, warn = FALSE), add_server_time = FALSE)
    if (length(clean_lines) == 0L) return(data.frame())
    dt <- data.table::fread(text = clean_lines, sep = ";", fill = TRUE, header = FALSE, showProgress = FALSE)
  } else if (is.character(input_data)) {
    clean_lines <- tt_clean_logger_lines(input_data, add_server_time = FALSE)
    if (length(clean_lines) == 0L) return(data.frame())
    dt <- data.table::fread(text = clean_lines, sep = ";", fill = TRUE, header = FALSE, showProgress = FALSE)
  } else {
    dt <- data.table::as.data.table(input_data)
  }
  if (nrow(dt) == 0L) return(data.frame())

  v1_chr <- as.character(dt[["V1"]])
  already_has_ts <- grepl("^[0-9]{2}\\.[0-9]{2}\\.[0-9]{2}\\s+[0-9]{2}:[0-9]{2}:[0-9]{2}\\s*,", v1_chr)
  ts_num <- suppressWarnings(as.numeric(dt[["V4"]]))
  hex_mask <- is.na(ts_num) & grepl("^[0-9A-Fa-f]{8}$", trimws(as.character(dt[["V4"]])))
  if (any(hex_mask)) {
    ts_num[hex_mask] <- suppressWarnings(as.numeric(strtoi(trimws(as.character(dt[["V4"]][hex_mask])), base = 16L)))
  }
  ts_posix <- as.POSIXct(ifelse(is.na(ts_num), 0, ts_num), origin = "1970-01-01", tz = "UTC")
  ts_fmt <- format(ts_posix, "%d.%m.%y %H:%M:%S")
  dt[, V1 := ifelse(already_has_ts, v1_chr, paste0(ts_fmt, ",", v1_chr))]

  if (!is.null(output_file)) {
    data.table::fwrite(dt, file = output_file, sep = ";", col.names = FALSE, quote = FALSE, na = "")
  }
  to_df(dt)
}

#' @rdname tt_restore_cloud_backup
#' @export
restore_cloud_backup <- tt_restore_cloud_backup

#' Scrape and assemble TreeTalker telemetry from cloud servers and local backups
#'
#' Downloads `ttcloud.txt` streams from standard Altervista endpoints or custom
#' URLs, cleans out any logger/terminal messages, optionally combines them with
#' local backup files and site metadata, removes duplicates, and parses packets
#' into typed `data.frame` tables.
#'
#' @param id Character or numeric vector of TreeTalker Cloud IDs (e.g., `"21030001"`).
#' @param custom_url Optional character vector of full URLs or local file paths.
#' @param import_folder Optional path to a local folder containing backup files
#'   and/or `tt_Descriptive.txt`.
#' @param metadata Optional `data.frame` with device metadata (`TT_ID`, `Species`,
#'   `DBH`, `Site`, `lat`, `lon`).
#' @param split_by_55 Logical; separate TT-Soil (`55`) packets from `4D` (default `TRUE`).
#' @return A named list of pure `data.frame` objects (output of [tt_types_to_list()]).
#' @export
tt_scrape <- function(id = NULL, custom_url = NULL, import_folder = NULL,
                      metadata = NULL, split_by_55 = TRUE) {
  urls <- character(0)
  if (!is.null(id) && length(id) > 0L) {
    urls <- c(
      paste0("http://naturetalkers.altervista.org/", id, "/ttcloud.txt"),
      paste0("http://ittn.altervista.org/", id, "/ttcloud.txt")
    )
  }
  if (!is.null(custom_url)) {
    urls <- c(urls, custom_url)
  }

  fetched <- lapply(urls, function(u) {
    tryCatch({
      raw_lines <- readLines(u, warn = FALSE)
      clean_lines <- tt_clean_logger_lines(raw_lines, add_server_time = TRUE)
      if (length(clean_lines) == 0L) return(NULL)
      data.table::fread(text = clean_lines, sep = ";", fill = TRUE, header = FALSE, showProgress = FALSE)
    }, error = function(e) NULL, warning = function(w) NULL)
  })
  dt_cloud <- data.table::rbindlist(fetched, fill = TRUE, use.names = TRUE)

  dt_local <- data.table::data.table()
  if (!is.null(import_folder) && dir.exists(import_folder)) {
    desc_file <- file.path(import_folder, "tt_Descriptive.txt")
    if (is.null(metadata) && file.exists(desc_file)) {
      metadata <- tryCatch(data.table::fread(desc_file, showProgress = FALSE), error = function(e) NULL)
    }
    bak_df <- tt_read_backup(import_folder, metadata = NULL)
    if (nrow(bak_df) > 0L) dt_local <- data.table::as.data.table(bak_df)
  }

  dt_all <- data.table::rbindlist(list(dt_cloud, dt_local), fill = TRUE, use.names = TRUE)
  if (nrow(dt_all) == 0L) return(tt_types_to_list(data.frame(), split_by_55 = split_by_55))

  dt_all <- merge_tt_metadata(dt_all, metadata = metadata)
  dedup_df <- tt_remove_duplicates(dt_all)
  tt_types_to_list(dedup_df, split_by_55 = split_by_55)
}

#' @rdname tt_scrape
#' @export
ttscrape <- tt_scrape

#' Subset TreeTalker data by device IDs and date range
#'
#' Filters a single `data.frame` or a named list of TreeTalker tables by
#' `TT_ID` (or `Species`/`Site`) and timestamp range without modifying the
#' global environment.
#'
#' @param mydata A `data.frame` or a named list of `data.frame`s from
#'   [tt_types_to_list()].
#' @param id_list Optional character/numeric vector of `TT_ID` values to keep.
#' @param start_date Optional start date/time (`POSIXct` or `"YYYY-MM-DD"`).
#' @param end_date Optional end date/time (`POSIXct` or `"YYYY-MM-DD"`).
#' @return Subsetted `data.frame` or list of `data.frame`s.
#' @export
tt_subset <- function(mydata, id_list = NULL, start_date = NULL, end_date = NULL) {
  sub_one <- function(df) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(to_df(df))
    dt <- data.table::as.data.table(df)
    if (!is.null(id_list) && "TT_ID" %in% names(dt)) {
      dt <- dt[as.character(TT_ID) %in% as.character(id_list)]
    }
    ts_col <- intersect(c("Time", "Datetime", "Timestamp", "Server_time"), names(dt))[1L]
    if (!is.na(ts_col)) {
      if (!is.null(start_date)) {
        s_dt <- as.POSIXct(start_date, tz = "UTC")
        dt <- dt[get(ts_col) >= s_dt]
      }
      if (!is.null(end_date)) {
        e_dt <- as.POSIXct(end_date, tz = "UTC")
        dt <- dt[get(ts_col) <= e_dt]
      }
    }
    to_df(dt)
  }

  if (is.list(mydata) && !is.data.frame(mydata)) {
    return(lapply(mydata, sub_one))
  }
  sub_one(mydata)
}

#' @rdname tt_subset
#' @export
ttsubset <- tt_subset

#' Export processed TreeTalker tables to CSV/TSV files
#'
#' Writes a processed `data.frame` or list of `data.frame`s to disk using fast
#' `data.table::fwrite()`.
#'
#' @param mydata A `data.frame` or list of `data.frame`s.
#' @param out_dir Output directory path (default `tempdir()`).
#' @param sep Field delimiter (default `";"`).
#' @return Invisibly returns a character vector of written file paths.
#' @export
tt_export <- function(mydata, out_dir = tempdir(), sep = ";") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (is.data.frame(mydata)) {
    out_path <- file.path(out_dir, "drevotalkR_export.csv")
    data.table::fwrite(mydata, file = out_path, sep = sep)
    return(invisible(out_path))
  }
  written <- character(0)
  for (nm in names(mydata)) {
    tbl <- mydata[[nm]]
    if (is.data.frame(tbl) && nrow(tbl) > 0L) {
      out_path <- file.path(out_dir, paste0(nm, ".csv"))
      data.table::fwrite(tbl, file = out_path, sep = sep)
      written <- c(written, out_path)
    }
  }
  invisible(written)
}

#' @rdname tt_export
#' @export
ttexport <- tt_export
