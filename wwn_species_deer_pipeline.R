#!/usr/bin/env Rscript

# End-to-end helper for:
# 1) Species-level tree-ring standardization (detrend + ARSTAN outputs)
# 2) Species establishment/die-off plots
# 3) Standardized deer series + lag correlations vs RWI, VPDMIN/VPDMAX, PDSI
# 4) Direct comparison: PDSI lag 2 vs best deer-VPD relationship
#
# Usage:
#   Rscript wwn_species_deer_pipeline.R [input_dir] [output_dir]

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else "."
output_dir <- if (length(args) >= 2) args[[2]] else file.path(".", "wwn_species_deer_outputs")
max_lag <- 2L
sampling_year <- 2024L
species_zoom_window_years <- 120L

# Royal palette aligned with existing project scripts.
COL_RWI <- "#7F0000"
COL_DEER <- "#A67C52"
COL_VPD <- "#4169E1"
COL_PPT <- "#0B3C8A"
COL_PDSI <- "#CC7722"
COL_FIT <- "#2F3E56"

# Keep local libraries first.
dir.create(".Rlibs", showWarnings = FALSE)
.libPaths(c(normalizePath(".Rlibs"), .libPaths()))

suppressPackageStartupMessages({
  library(dplR)
  library(readxl)
  library(ggplot2)
})

canon <- function(x) gsub("[^a-z0-9]+", "", tolower(x))

ww_title <- function(subject, method) {
  paste0("Wesselman Woods - ", subject, " (", method, ")")
}

pick_col <- function(nms, patterns) {
  for (p in patterns) {
    idx <- which(grepl(p, nms, fixed = TRUE))
    if (length(idx) > 0) return(idx[1])
  }
  NA_integer_
}

to_year_df <- function(x) {
  x <- as.data.frame(x, check.names = FALSE)
  yrs <- suppressWarnings(as.integer(rownames(x)))
  out <- data.frame(year = yrs, x, check.names = FALSE)
  out <- out[!is.na(out$year), , drop = FALSE]
  out <- out[order(out$year), , drop = FALSE]
  rownames(out) <- NULL
  out
}

estimate_ar1 <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  out <- tryCatch(as.numeric(acf(x, lag.max = 1, plot = FALSE)$acf[2]), error = function(e) NA_real_)
  if (!is.finite(out)) NA_real_ else out
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 4) return(list(r = NA_real_, p = NA_real_, n = n, r2 = NA_real_))
  ct <- tryCatch(cor.test(x[ok], y[ok], method = method), error = function(e) NULL)
  if (is.null(ct)) return(list(r = NA_real_, p = NA_real_, n = n, r2 = NA_real_))
  r <- as.numeric(ct$estimate)
  list(r = r, p = as.numeric(ct$p.value), n = n, r2 = r^2)
}

lag_vec <- function(x, k) {
  if (k <= 0) return(x)
  c(rep(NA_real_, k), head(x, -k))
}

zfun <- function(x) {
  x <- as.numeric(x)
  if (sum(is.finite(x)) < 2) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

clean_series_id <- function(x) {
  x <- gsub("\\s+", "", x)
  x <- gsub("[^A-Za-z0-9]", "", x)
  x[x == ""] <- "SERIES"
  toupper(x)
}

parse_loose_tucson <- function(path) {
  lines <- readLines(path, warn = FALSE)
  series_vec <- character(0)
  year_vec <- integer(0)
  value_vec <- numeric(0)
  skipped <- 0L

  for (line in lines) {
    line <- trimws(line)
    if (line == "" || startsWith(line, "#")) next

    tokens <- strsplit(line, "\\s+")[[1]]
    if (length(tokens) < 2) {
      skipped <- skipped + 1L
      next
    }

    sid <- NA_character_
    decade <- NA_integer_
    vals <- character(0)

    # Pattern A: <ID> <YYYY> <values...>
    if (length(tokens) >= 3 && grepl("^-?\\d{4}$", tokens[2])) {
      y <- suppressWarnings(as.integer(tokens[2]))
      if (!is.na(y) && y >= 1500 && y <= 2200) {
        sid <- tokens[1]
        decade <- y
        vals <- tokens[-c(1, 2)]
      }
    }

    # Pattern B: <ID><YYYY> <values...> or <ID>_<YYYY> <values...>
    if (is.na(decade)) {
      m <- regexec("^(.+?)[_]?([12][0-9]{3})$", tokens[1])
      mm <- regmatches(tokens[1], m)[[1]]
      if (length(mm) == 3) {
        y <- suppressWarnings(as.integer(mm[3]))
        if (!is.na(y) && y >= 1500 && y <= 2200) {
          sid <- mm[2]
          decade <- y
          vals <- tokens[-1]
        }
      }
    }

    if (is.na(decade) || length(vals) == 0) {
      skipped <- skipped + 1L
      next
    }

    sid <- clean_series_id(sid)

    for (i in seq_along(vals)) {
      v <- suppressWarnings(as.numeric(vals[i]))
      if (!is.finite(v)) next
      if (v %in% c(999, -9999)) break
      if (v <= 0) next
      series_vec <- c(series_vec, sid)
      year_vec <- c(year_vec, decade + i - 1L)
      value_vec <- c(value_vec, v)
    }
  }

  if (length(value_vec) == 0) stop("Loose parser could not recover values from ", basename(path))

  d <- data.frame(series = series_vec, year = year_vec, value = value_vec, stringsAsFactors = FALSE)
  k <- paste(d$series, d$year, sep = "::")
  dup <- sum(duplicated(k))
  d <- d[!duplicated(k), , drop = FALSE]

  years <- sort(unique(d$year))
  ids <- sort(unique(d$series))
  mat <- matrix(NA_real_, nrow = length(years), ncol = length(ids), dimnames = list(as.character(years), ids))
  mat[cbind(match(d$year, years), match(d$series, ids))] <- d$value

  list(
    rwl = as.data.frame(mat, check.names = FALSE),
    note = paste0("loose_parser used; skipped_lines=", skipped, "; duplicate_values_dropped=", dup)
  )
}

safe_read_rwl <- function(path) {
  out <- tryCatch(read.rwl(path, verbose = FALSE), error = function(e) e)
  if (!inherits(out, "error")) {
    return(list(rwl = as.data.frame(out, check.names = FALSE), method = "read.rwl", note = ""))
  }
  loose <- parse_loose_tucson(path)
  list(rwl = loose$rwl, method = "loose_parser", note = paste("read.rwl failed:", conditionMessage(out), "|", loose$note))
}

sanitize_rwl <- function(rwl) {
  rwl <- as.data.frame(rwl, check.names = FALSE)
  years <- suppressWarnings(as.integer(rownames(rwl)))
  keep <- !is.na(years)
  rwl <- rwl[keep, , drop = FALSE]
  years <- years[keep]
  if (length(years) == 0 || ncol(rwl) == 0) stop("No usable years/columns after read.")

  ord <- order(years)
  rwl <- rwl[ord, , drop = FALSE]
  rownames(rwl) <- as.character(years[ord])

  if (any(duplicated(rownames(rwl)))) rwl <- rwl[!duplicated(rownames(rwl)), , drop = FALSE]

  colnames(rwl) <- make.unique(clean_series_id(colnames(rwl)), sep = "_dup")
  rwl[] <- lapply(rwl, function(v) {
    v <- suppressWarnings(as.numeric(v))
    v[!is.finite(v)] <- NA_real_
    v[v <= 0] <- NA_real_
    v
  })

  keep_col <- colSums(!is.na(rwl)) > 0
  rwl <- rwl[, keep_col, drop = FALSE]
  keep_row <- rowSums(!is.na(rwl)) > 0
  rwl <- rwl[keep_row, , drop = FALSE]

  if (nrow(rwl) == 0 || ncol(rwl) == 0) stop("No usable values after sanitize.")
  rwl
}

make_wide <- function(chron_list, col_name) {
  pieces <- lapply(names(chron_list), function(sp) {
    x <- chron_list[[sp]]
    if (!(col_name %in% names(x))) return(NULL)
    data.frame(year = as.integer(rownames(x)), value = as.numeric(x[[col_name]]), stringsAsFactors = FALSE)
  })
  names(pieces) <- names(chron_list)
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (length(pieces) == 0) return(data.frame(year = integer(0)))

  out <- NULL
  for (sp in names(pieces)) {
    d <- pieces[[sp]]
    names(d)[2] <- sp
    out <- if (is.null(out)) d else merge(out, d, by = "year", all = TRUE, sort = TRUE)
  }
  out
}

prewhiten_series <- function(x, max_order = 5L) {
  x <- as.numeric(x)
  idx <- which(is.finite(x))
  out <- rep(NA_real_, length(x))
  if (length(idx) < 8) return(list(residual = out, order = NA_integer_))

  x0 <- x[idx]
  order_max <- min(max_order, max(1L, floor(length(x0) / 3)))
  fit <- tryCatch(ar(x0, aic = TRUE, order.max = order_max, method = "yw"), error = function(e) NULL)
  if (is.null(fit) || is.null(fit$resid)) return(list(residual = out, order = NA_integer_))

  resid_vals <- as.numeric(fit$resid)
  if (length(resid_vals) < length(x0)) resid_vals <- c(rep(NA_real_, length(x0) - length(resid_vals)), resid_vals)
  out[idx[seq_along(resid_vals)]] <- resid_vals
  list(residual = out, order = fit$order)
}

process_species_file <- function(file, species_name, species_out_dir) {
  prefix <- substr(paste0(toupper(species_name), "XXX"), 1, 3)

  tryCatch({
    read_out <- safe_read_rwl(file)
    rwl <- sanitize_rwl(read_out$rwl)

    # Per-series first/last year (establishment/die-off proxy).
    rwl_years <- suppressWarnings(as.integer(rownames(rwl)))
    series_first <- vapply(rwl, function(v) {
      yy <- rwl_years[is.finite(v)]
      if (length(yy) == 0) NA_integer_ else min(yy)
    }, integer(1))
    series_last <- vapply(rwl, function(v) {
      yy <- rwl_years[is.finite(v)]
      if (length(yy) == 0) NA_integer_ else max(yy)
    }, integer(1))
    life_df <- data.frame(
      species = species_name,
      series_id = names(series_first),
      first_year = as.integer(series_first),
      last_year = as.integer(series_last),
      stringsAsFactors = FALSE
    )
    life_df <- life_df[is.finite(life_df$first_year) & is.finite(life_df$last_year), , drop = FALSE]

    write.csv(to_year_df(rwl), file.path(species_out_dir, paste0(species_name, "_rwl_cleaned_matrix.csv")), row.names = FALSE, na = "")
    tryCatch(write.rwl(rwl, file.path(species_out_dir, paste0(species_name, "_cleaned.rwl"))), error = function(e) NULL)

    rwl <- fill.internal.NA(rwl)
    rwi <- detrend(rwl = rwl, method = "Spline", nyrs = NULL, f = 0.5, pos.slope = FALSE)
    crn_std <- chron(rwi, prefix = prefix, biweight = TRUE, prewhiten = FALSE)
    crn_pw <- chron(rwi, prefix = prefix, biweight = TRUE, prewhiten = TRUE)
    crn_ars <- tryCatch(
      chron.ars(rwi, biweight = TRUE, maxLag = 10, firstAICmin = TRUE, verbose = FALSE, prewhitenMethod = c("ar.yw", "arima.CSS-ML")),
      error = function(e) e
    )

    write.csv(to_year_df(rwi), file.path(species_out_dir, paste0(species_name, "_rwi.csv")), row.names = FALSE, na = "")
    write.csv(to_year_df(crn_std), file.path(species_out_dir, paste0(species_name, "_chron_standard.csv")), row.names = FALSE, na = "")
    write.csv(to_year_df(crn_pw), file.path(species_out_dir, paste0(species_name, "_chron_prewhiten.csv")), row.names = FALSE, na = "")

    arstan_error <- ""
    arstan_obj <- NULL
    if (inherits(crn_ars, "error")) {
      arstan_error <- conditionMessage(crn_ars)
    } else {
      arstan_obj <- crn_ars
      write.csv(to_year_df(crn_ars), file.path(species_out_dir, paste0(species_name, "_chron_arstan.csv")), row.names = FALSE, na = "")
    }

    yrs <- suppressWarnings(as.integer(rownames(rwl)))
    row_out <- data.frame(
      species = species_name,
      status = "ok",
      input_file = normalizePath(file, winslash = "/", mustWork = FALSE),
      read_method = read_out$method,
      n_series = ncol(rwl),
      first_year = min(yrs, na.rm = TRUE),
      last_year = max(yrs, na.rm = TRUE),
      ar1_standard = estimate_ar1(crn_pw$std),
      ar1_residual = estimate_ar1(crn_pw$res),
      arstan_error = arstan_error,
      read_note = read_out$note,
      error = "",
      stringsAsFactors = FALSE
    )

    list(row = row_out, life_df = life_df, prewhite = crn_pw, arstan = arstan_obj)
  }, error = function(e) {
    row_out <- data.frame(
      species = species_name,
      status = "error",
      input_file = normalizePath(file, winslash = "/", mustWork = FALSE),
      read_method = NA_character_,
      n_series = NA_integer_,
      first_year = NA_integer_,
      last_year = NA_integer_,
      ar1_standard = NA_real_,
      ar1_residual = NA_real_,
      arstan_error = NA_character_,
      read_note = "",
      error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
    list(row = row_out, life_df = NULL, prewhite = NULL, arstan = NULL)
  })
}

# --------------------
# A) Species ARSTAN run + establishment plots
# --------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
species_dir <- file.path(output_dir, "species")
dir.create(species_dir, recursive = TRUE, showWarnings = FALSE)

all_rwl <- sort(list.files(input_dir, pattern = "\\.rwl$", ignore.case = TRUE, full.names = TRUE))
all_stems <- tools::file_path_sans_ext(basename(all_rwl))

# Species files are uppercase alpha-numeric stems (e.g., ACER1, QUER, LIST).
is_species <- grepl("^[A-Z][A-Z0-9]*$", all_stems)
species_files <- all_rwl[is_species]
if (length(species_files) == 0) {
  warning("No species-like .rwl files matched; falling back to all .rwl files.")
  species_files <- all_rwl
}

message("Species files to process: ", length(species_files))

species_summary <- list()
prewhite_list <- list()
arstan_list <- list()
lifespan_rows <- list()

for (i in seq_along(species_files)) {
  f <- species_files[i]
  sp <- tools::file_path_sans_ext(basename(f))
  sp_out <- file.path(species_dir, sp)
  dir.create(sp_out, recursive = TRUE, showWarnings = FALSE)

  message(sprintf("[%d/%d] %s", i, length(species_files), basename(f)))

  res <- process_species_file(f, sp, sp_out)
  species_summary[[length(species_summary) + 1L]] <- res$row
  if (!is.null(res$life_df) && nrow(res$life_df) > 0) lifespan_rows[[length(lifespan_rows) + 1L]] <- res$life_df
  if (!is.null(res$prewhite)) prewhite_list[[sp]] <- res$prewhite
  if (!is.null(res$arstan)) arstan_list[[sp]] <- res$arstan
}

if (length(species_summary) > 0) {
  species_summary_df <- do.call(rbind, species_summary)
  write.csv(species_summary_df, file.path(species_dir, "species_batch_summary.csv"), row.names = FALSE, na = "")
}

if (length(prewhite_list) > 0) {
  write.csv(make_wide(prewhite_list, "std"), file.path(species_dir, "species_standard_chronologies_wide.csv"), row.names = FALSE, na = "")
  write.csv(make_wide(prewhite_list, "res"), file.path(species_dir, "species_residual_chronologies_wide.csv"), row.names = FALSE, na = "")
  write.csv(make_wide(prewhite_list, "samp.depth"), file.path(species_dir, "species_sample_depth_wide.csv"), row.names = FALSE, na = "")
}

if (length(arstan_list) > 0) {
  write.csv(make_wide(arstan_list, "ars"), file.path(species_dir, "species_arstan_chronologies_wide.csv"), row.names = FALSE, na = "")
}

# Establishment/die-off outputs requested by Jim.
if (length(lifespan_rows) > 0) {
  lifespan_df <- do.call(rbind, lifespan_rows)
  lifespan_df$terminal_year <- pmin(as.integer(lifespan_df$last_year), sampling_year)
  lifespan_df$censored_at_sampling <- as.integer(lifespan_df$last_year) >= sampling_year
  lifespan_df$death_year <- ifelse(lifespan_df$censored_at_sampling, NA_integer_, lifespan_df$terminal_year)
  lifespan_df$death_year <- as.integer(lifespan_df$death_year)
  write.csv(lifespan_df, file.path(species_dir, "species_series_lifespan.csv"), row.names = FALSE, na = "")

  counts_list <- lapply(split(lifespan_df, lifespan_df$species), function(df_sp) {
    if (nrow(df_sp) == 0) return(NULL)
    yr_start <- min(df_sp$first_year, na.rm = TRUE)
    yr_end <- sampling_year
    if (!is.finite(yr_start) || !is.finite(yr_end) || yr_start > yr_end) return(NULL)
    yr_seq <- seq(yr_start, yr_end)
    est_counts <- tabulate(match(df_sp$first_year, yr_seq), nbins = length(yr_seq))
    die_match <- match(df_sp$death_year, yr_seq)
    die_match <- die_match[is.finite(die_match)]
    die_counts <- tabulate(die_match, nbins = length(yr_seq))
    active_counts <- vapply(yr_seq, function(y) sum(df_sp$first_year <= y & df_sp$terminal_year >= y, na.rm = TRUE), integer(1))
    alive_sampling <- sum(df_sp$censored_at_sampling, na.rm = TRUE)
    data.frame(
      species = unique(df_sp$species),
      year = yr_seq,
      established = est_counts,
      died_off = die_counts,
      active_series = active_counts,
      alive_at_sampling = ifelse(yr_seq == sampling_year, alive_sampling, 0L),
      stringsAsFactors = FALSE
    )
  })
  counts_list <- Filter(Negate(is.null), counts_list)

  if (length(counts_list) > 0) {
    species_counts <- do.call(rbind, counts_list)
    write.csv(species_counts, file.path(species_dir, "species_establishment_dieoff_by_year.csv"), row.names = FALSE, na = "")

    events_df <- rbind(
      data.frame(species = species_counts$species, year = species_counts$year, event = "Established", count = species_counts$established, stringsAsFactors = FALSE),
      data.frame(species = species_counts$species, year = species_counts$year, event = "Died off", count = -species_counts$died_off, stringsAsFactors = FALSE),
      data.frame(species = species_counts$species, year = species_counts$year, event = "Alive at sampling", count = species_counts$alive_at_sampling, stringsAsFactors = FALSE)
    )

    events_xmin <- min(species_counts$year, na.rm = TRUE)
    events_xmax <- max(species_counts$year, na.rm = TRUE)
    zoom_start <- max(events_xmin, sampling_year - species_zoom_window_years)
    zoom_break <- if ((sampling_year - zoom_start) <= 120) 10 else 20

    p_events_base <- ggplot(events_df, aes(x = year, y = count, fill = event)) +
      geom_col(width = 0.9) +
      facet_wrap(~species, scales = "free_y", ncol = 3) +
      geom_vline(xintercept = sampling_year, color = COL_VPD, linetype = "dashed", linewidth = 0.45) +
      scale_fill_manual(values = c("Established" = COL_PPT, "Died off" = COL_PDSI, "Alive at sampling" = COL_VPD)) +
      labs(
        title = ww_title("Species Establishment and Die-Off", "Series Timing by Year"),
        subtitle = paste0(
          "Positive bars = start year; negative bars = inferred die-off year. ",
          "Series ending in ", sampling_year, " are treated as alive at sampling."
        ),
        x = "Year",
        y = "Series count"
      ) +
      theme_minimal(base_size = 11)

    p_events_full <- p_events_base +
      scale_x_continuous(breaks = seq(events_xmin, events_xmax, by = 20))
    ggsave(file.path(species_dir, "species_establishment_dieoff_faceted_fullrange.png"), p_events_full, width = 14.5, height = 9.6, dpi = 320)

    p_events_zoom <- p_events_base +
      coord_cartesian(xlim = c(zoom_start, sampling_year)) +
      scale_x_continuous(breaks = seq(zoom_start, sampling_year, by = zoom_break))
    ggsave(file.path(species_dir, "species_establishment_dieoff_faceted.png"), p_events_zoom, width = 14.5, height = 9.6, dpi = 320)

    p_active_base <- ggplot(species_counts, aes(x = year, y = active_series)) +
      geom_line(linewidth = 0.9, color = COL_VPD) +
      facet_wrap(~species, scales = "free_y", ncol = 3) +
      geom_vline(xintercept = sampling_year, color = COL_PDSI, linetype = "dashed", linewidth = 0.45) +
      labs(
        title = ww_title("Species Establishment Status", "Active Series Through Time"),
        subtitle = paste0("Counts are right-censored at sampling year ", sampling_year, "."),
        x = "Year",
        y = "Active series"
      ) +
      theme_minimal(base_size = 11)

    p_active_full <- p_active_base +
      scale_x_continuous(breaks = seq(events_xmin, events_xmax, by = 20))
    ggsave(file.path(species_dir, "species_active_series_through_time_fullrange.png"), p_active_full, width = 14.5, height = 9.6, dpi = 320)

    p_active_zoom <- p_active_base +
      coord_cartesian(xlim = c(zoom_start, sampling_year)) +
      scale_x_continuous(breaks = seq(zoom_start, sampling_year, by = zoom_break))
    ggsave(file.path(species_dir, "species_active_series_through_time.png"), p_active_zoom, width = 14.5, height = 9.6, dpi = 320)
  }
}
# --------------------
# B) Deer standardization + lag analysis
# --------------------
deer_dir <- file.path(output_dir, "deer")
dir.create(deer_dir, recursive = TRUE, showWarnings = FALSE)

deer_file <- file.path(input_dir, "wwn_deer.xlsx")
if (!file.exists(deer_file)) stop("Missing deer workbook: ", deer_file)

deer_raw <- read_excel(deer_file, sheet = 1)
deer_names <- canon(names(deer_raw))

year_idx <- pick_col(deer_names, c("year", "yr"))
deer_idx <- pick_col(deer_names, c("deerpopulation", "deerpop", "deer"))
rwi_idx <- pick_col(deer_names, c("rwi"))

if (is.na(year_idx)) stop("Could not find YEAR column in wwn_deer.xlsx")
if (is.na(deer_idx)) stop("Could not find deer population column in wwn_deer.xlsx")

deer_df <- data.frame(
  year = as.integer(as.numeric(deer_raw[[year_idx]])),
  deer_population = as.numeric(deer_raw[[deer_idx]]),
  stringsAsFactors = FALSE
)
if (!is.na(rwi_idx)) deer_df$RWI <- as.numeric(deer_raw[[rwi_idx]])
deer_df <- deer_df[is.finite(deer_df$year), , drop = FALSE]
deer_df <- deer_df[order(deer_df$year), , drop = FALSE]

# Load PRISM monthly file (latest by name), aggregate to annual means.
prism_candidates <- sort(list.files(input_dir, pattern = "^PRISM_.*vpdmin.*vpdmax.*\\.csv$", ignore.case = TRUE, full.names = TRUE))
if (length(prism_candidates) == 0) stop("No PRISM file with vpdmin/vpdmax found in input_dir.")
prism_file <- prism_candidates[length(prism_candidates)]

prism_lines <- readLines(prism_file, warn = FALSE)
header_idx <- which(grepl("^\\s*\\ufeff?[Dd]ate\\s*,", prism_lines))[1]
if (is.na(header_idx)) {
  prism_line_norm <- tolower(gsub("\"", "", prism_lines))
  header_idx <- which(grepl("date\\s*,\\s*ppt\\s*,\\s*tmin\\s*,\\s*tmean\\s*,\\s*tmax", prism_line_norm))[1]
}
if (is.na(header_idx)) stop("Could not find PRISM header row containing date,ppt,tmin,tmean,tmax.")

prism_raw <- read.csv(prism_file, skip = header_idx - 1, stringsAsFactors = FALSE, check.names = FALSE)
prism_names <- canon(names(prism_raw))

date_idx <- pick_col(prism_names, c("date"))
vpdmin_idx <- pick_col(prism_names, c("vpdmin"))
vpdmax_idx <- pick_col(prism_names, c("vpdmax"))

if (is.na(date_idx) || is.na(vpdmin_idx) || is.na(vpdmax_idx)) {
  stop("Could not find Date/VPDMIN/VPDMAX columns in PRISM file: ", basename(prism_file))
}

prism_date_raw <- as.character(prism_raw[[date_idx]])
prism_date_digits <- gsub("[^0-9]", "", prism_date_raw)

prism_ann <- data.frame(
  year = suppressWarnings(as.integer(substr(prism_date_digits, 1, 4))),
  month = suppressWarnings(as.integer(substr(prism_date_digits, 5, 6))),
  VPDMIN = as.numeric(prism_raw[[vpdmin_idx]]),
  VPDMAX = as.numeric(prism_raw[[vpdmax_idx]]),
  stringsAsFactors = FALSE
)

prism_month_counts <- aggregate(month ~ year, data = prism_ann, FUN = function(z) length(unique(z[is.finite(z)])))
names(prism_month_counts)[2] <- "PRISM_month_n"
prism_incomplete_years <- prism_month_counts$year[is.finite(prism_month_counts$PRISM_month_n) & prism_month_counts$PRISM_month_n < 12]

prism_ann <- aggregate(cbind(VPDMIN, VPDMAX) ~ year, data = prism_ann, FUN = function(z) mean(z, na.rm = TRUE))
prism_ann <- merge(prism_ann, prism_month_counts, by = "year", all.x = TRUE, sort = TRUE)
prism_ann <- prism_ann[is.na(prism_ann$PRISM_month_n) | prism_ann$PRISM_month_n == 12, , drop = FALSE]

vpdmin_vals <- as.numeric(prism_raw[[vpdmin_idx]])
vpdmax_vals <- as.numeric(prism_raw[[vpdmax_idx]])
neg_vpdmin_mask <- is.finite(vpdmin_vals) & vpdmin_vals < 0
neg_vpdmax_mask <- is.finite(vpdmax_vals) & vpdmax_vals < 0
vpd_order_mask <- is.finite(vpdmin_vals) & is.finite(vpdmax_vals) & vpdmin_vals > vpdmax_vals

climate_qc <- data.frame(
  check = c(
    "Negative monthly VPDMIN",
    "Negative monthly VPDMAX",
    "Monthly rows with VPDMIN > VPDMAX",
    "Incomplete PRISM years excluded from annual deer regressions"
  ),
  n_flagged = c(
    sum(neg_vpdmin_mask),
    sum(neg_vpdmax_mask),
    sum(vpd_order_mask),
    length(prism_incomplete_years)
  ),
  years = c(
    if (sum(neg_vpdmin_mask) == 0) "" else paste(sort(unique(as.integer(substr(prism_date_digits[neg_vpdmin_mask], 1, 4)))), collapse = ", "),
    if (sum(neg_vpdmax_mask) == 0) "" else paste(sort(unique(as.integer(substr(prism_date_digits[neg_vpdmax_mask], 1, 4)))), collapse = ", "),
    if (sum(vpd_order_mask) == 0) "" else paste(sort(unique(as.integer(substr(prism_date_digits[vpd_order_mask], 1, 4)))), collapse = ", "),
    if (length(prism_incomplete_years) == 0) "" else paste(sort(unique(prism_incomplete_years)), collapse = ", ")
  ),
  note = c(
    "Vapor pressure deficit should not be negative",
    "Vapor pressure deficit should not be negative",
    "Minimum VPD should not exceed maximum VPD",
    if (length(prism_incomplete_years) == 0) "All PRISM deer years had 12 months available" else paste0("Dropped incomplete PRISM years: ", paste(sort(unique(prism_incomplete_years)), collapse = ", "))
  ),
  stringsAsFactors = FALSE
)

# PDSI source priority:
# 1) pdsi.csv with Date + "PDSI avg" (or similar),
# 2) chronology_climate.xlsx fallback.
pdsi_ann <- NULL
pdsi_source_file <- "not found"

pdsi_csv_path <- file.path(input_dir, "pdsi.csv")
if (file.exists(pdsi_csv_path)) {
  pd <- read.csv(pdsi_csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  pd_names <- canon(names(pd))
  pd_date_idx <- pick_col(pd_names, c("date", "year", "yr"))
  pd_pdsi_idx <- pick_col(pd_names, c("pdsiavg", "pdsi"))

  if (!is.na(pd_date_idx) && !is.na(pd_pdsi_idx)) {
    date_raw <- gsub("[^0-9]", "", as.character(pd[[pd_date_idx]]))
    pdsi_ann <- data.frame(
      year = suppressWarnings(as.integer(substr(date_raw, 1, 4))),
      PDSI = as.numeric(pd[[pd_pdsi_idx]]),
      stringsAsFactors = FALSE
    )
    pdsi_ann <- pdsi_ann[is.finite(pdsi_ann$year), , drop = FALSE]
    if (nrow(pdsi_ann) > 0) {
      pdsi_ann <- aggregate(PDSI ~ year, data = pdsi_ann, FUN = function(z) mean(z, na.rm = TRUE))
      pdsi_source_file <- basename(pdsi_csv_path)
    } else {
      pdsi_ann <- NULL
    }
  }
}

if (is.null(pdsi_ann)) {
  pdsi_xlsx_path <- file.path(input_dir, "chronology_climate.xlsx")
  if (file.exists(pdsi_xlsx_path)) {
    cc <- read_excel(pdsi_xlsx_path, sheet = 1)
    cc_names <- canon(names(cc))
    cc_year_idx <- pick_col(cc_names, c("year", "yr"))
    cc_pdsi_idx <- pick_col(cc_names, c("pdsi"))
    if (!is.na(cc_year_idx) && !is.na(cc_pdsi_idx)) {
      pdsi_ann <- data.frame(
        year = as.integer(as.numeric(cc[[cc_year_idx]])),
        PDSI = as.numeric(cc[[cc_pdsi_idx]]),
        stringsAsFactors = FALSE
      )
      pdsi_ann <- pdsi_ann[is.finite(pdsi_ann$year), , drop = FALSE]
      if (nrow(pdsi_ann) > 0) {
        pdsi_ann <- aggregate(PDSI ~ year, data = pdsi_ann, FUN = function(z) mean(z, na.rm = TRUE))
        pdsi_source_file <- basename(pdsi_xlsx_path)
      } else {
        pdsi_ann <- NULL
      }
    }
  }
}

# Merge deer + RWI + annual VPD + PDSI
deer_merge <- merge(deer_df, prism_ann, by = "year", all.x = TRUE, sort = TRUE)
if (!is.null(pdsi_ann)) deer_merge <- merge(deer_merge, pdsi_ann, by = "year", all.x = TRUE, sort = TRUE)

write.csv(climate_qc, file.path(deer_dir, "climate_qc_flags.csv"), row.names = FALSE, na = "")

# Deer standardization and prewhitening
deer_merge$deer_z <- zfun(deer_merge$deer_population)
prew <- prewhiten_series(deer_merge$deer_z, max_order = 5L)
deer_merge$deer_resid <- prew$residual

deer_merge$RWI_z <- if ("RWI" %in% names(deer_merge)) zfun(deer_merge$RWI) else NA_real_
deer_merge$VPDMIN_z <- zfun(deer_merge$VPDMIN)
deer_merge$VPDMAX_z <- zfun(deer_merge$VPDMAX)
if ("PDSI" %in% names(deer_merge)) deer_merge$PDSI_z <- zfun(deer_merge$PDSI)

write.csv(deer_merge, file.path(deer_dir, "deer_standardized_timeseries.csv"), row.names = FALSE, na = "")

# Lag correlations: standardized + prewhitened deer vs RWI, VPDMIN, VPDMAX, PDSI
predictors <- c("RWI", "VPDMIN", "VPDMAX", "PDSI")
predictors <- predictors[predictors %in% names(deer_merge)]
responses <- c("deer_z", "deer_resid")
methods <- c("pearson", "spearman")

lag_rows <- list()
k <- 1L
for (resp in responses) {
  y <- deer_merge[[resp]]
  for (pred in predictors) {
    x0 <- as.numeric(deer_merge[[pred]])
    for (lg in 0:max_lag) {
      x <- lag_vec(x0, lg)
      for (meth in methods) {
        ct <- safe_cor(y, x, method = meth)
        lag_rows[[k]] <- data.frame(
          response = resp,
          predictor = pred,
          lag = lg,
          method = meth,
          r = ct$r,
          r2 = ct$r2,
          p = ct$p,
          n = ct$n,
          stringsAsFactors = FALSE
        )
        k <- k + 1L
      }
    }
  }
}

lag_df <- if (length(lag_rows) > 0) do.call(rbind, lag_rows) else data.frame()
write.csv(lag_df, file.path(deer_dir, "deer_standardized_lag_correlations.csv"), row.names = FALSE, na = "")

if (nrow(lag_df) > 0) {
  best <- lag_df[is.finite(lag_df$r), , drop = FALSE]
  if (nrow(best) > 0) {
    best <- best[order(-abs(best$r), best$p), , drop = FALSE]
    write.csv(head(best, 20), file.path(deer_dir, "deer_standardized_lag_top20.csv"), row.names = FALSE, na = "")
  }
}

# Scatter panels with explicit R2 labels (addresses Jim's request).
plot_predictor_colors <- c(RWI = COL_RWI, VPDMIN = COL_VPD, VPDMAX = COL_PPT, PDSI = COL_PDSI)
for (resp in responses) {
  best_rows <- lag_df[
    lag_df$response == resp & lag_df$method == "pearson" & is.finite(lag_df$r),
    , drop = FALSE
  ]
  if (nrow(best_rows) == 0) next

  best_rows <- do.call(rbind, lapply(split(best_rows, best_rows$predictor), function(dfp) dfp[which.max(abs(dfp$r)), , drop = FALSE]))
  if (nrow(best_rows) == 0) next

  panel_rows <- list()
  ann_rows <- list()
  j <- 1L

  for (ii in seq_len(nrow(best_rows))) {
    pred <- best_rows$predictor[ii]
    lg <- best_rows$lag[ii]
    x <- lag_vec(as.numeric(deer_merge[[pred]]), lg)
    y <- as.numeric(deer_merge[[resp]])
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 4) next

    panel_name <- paste0(pred, " (lag ", lg, ")")
    panel_rows[[j]] <- data.frame(
      predictor_value = x[ok],
      response_value = y[ok],
      panel = panel_name,
      predictor = pred,
      stringsAsFactors = FALSE
    )

    ct <- safe_cor(y, x, method = "pearson")
    ann_rows[[j]] <- data.frame(
      panel = panel_name,
      label = paste0("r = ", sprintf("%.3f", ct$r), "\nR2 = ", sprintf("%.3f", ct$r2), "\np = ", sprintf("%.4f", ct$p), "\nn = ", ct$n),
      stringsAsFactors = FALSE
    )
    j <- j + 1L
  }

  if (length(panel_rows) == 0) next

  plot_df <- do.call(rbind, panel_rows)
  ann_df <- do.call(rbind, ann_rows)

  p_scatter <- ggplot(plot_df, aes(x = predictor_value, y = response_value, color = predictor)) +
    geom_point(alpha = 0.9, size = 2.4) +
    geom_smooth(method = "lm", se = FALSE, color = COL_FIT, linewidth = 0.9) +
    facet_wrap(~panel, scales = "free_x") +
    scale_color_manual(values = plot_predictor_colors[intersect(names(plot_predictor_colors), unique(plot_df$predictor))]) +
    geom_text(
      data = ann_df,
      aes(label = label),
      x = -Inf,
      y = Inf,
      hjust = -0.05,
      vjust = 1.1,
      color = COL_FIT,
      size = 3.15,
      inherit.aes = FALSE
    ) +
    labs(
      title = ww_title("Deer Population Predictors", paste0("Best-Lag Regression Panels (", resp, ")")),
      subtitle = "Each panel shows best lag by absolute Pearson r with explicit R2",
      x = "Lagged predictor value",
      y = resp,
      color = "Predictor"
    ) +
    theme_minimal(base_size = 11)

  ggsave(
    file.path(deer_dir, paste0("deer_best_lag_scatter_r2_", resp, ".png")),
    p_scatter,
    width = 11.5,
    height = 6.8,
    dpi = 320
  )
}

# Direct comparison for Evan: PDSI lag 2 vs best VPD lag (same response + method)
comparison_rows <- list()
if ("PDSI" %in% predictors && any(c("VPDMIN", "VPDMAX") %in% predictors) && nrow(lag_df) > 0) {
  cidx <- 1L
  for (resp in responses) {
    for (meth in methods) {
      pdsi_l2 <- lag_df[
        lag_df$response == resp & lag_df$method == meth & lag_df$predictor == "PDSI" & lag_df$lag == 2,
        , drop = FALSE
      ]
      vpd_rows <- lag_df[
        lag_df$response == resp & lag_df$method == meth & lag_df$predictor %in% c("VPDMIN", "VPDMAX") & is.finite(lag_df$r),
        , drop = FALSE
      ]

      if (nrow(pdsi_l2) == 1 && nrow(vpd_rows) > 0) {
        best_vpd <- vpd_rows[which.max(abs(vpd_rows$r)), , drop = FALSE]
        winner <- if (abs(best_vpd$r) > abs(pdsi_l2$r)) "VPD stronger" else if (abs(best_vpd$r) < abs(pdsi_l2$r)) "PDSI lag2 stronger" else "Tie"
        comparison_rows[[cidx]] <- data.frame(
          response = resp,
          method = meth,
          pdsi_lag = 2,
          pdsi_r = pdsi_l2$r,
          pdsi_r2 = pdsi_l2$r2,
          pdsi_p = pdsi_l2$p,
          best_vpd_predictor = best_vpd$predictor,
          best_vpd_lag = best_vpd$lag,
          best_vpd_r = best_vpd$r,
          best_vpd_r2 = best_vpd$r2,
          best_vpd_p = best_vpd$p,
          delta_abs_r = abs(best_vpd$r) - abs(pdsi_l2$r),
          delta_r2 = best_vpd$r2 - pdsi_l2$r2,
          winner = winner,
          stringsAsFactors = FALSE
        )
        cidx <- cidx + 1L
      }
    }
  }
}

if (length(comparison_rows) > 0) {
  comp_df <- do.call(rbind, comparison_rows)
  write.csv(comp_df, file.path(deer_dir, "pdsi_lag2_vs_vpd_comparison.csv"), row.names = FALSE, na = "")

  plot_df <- rbind(
    data.frame(
      response = comp_df$response,
      method = comp_df$method,
      model = "PDSI lag 2",
      abs_r = abs(comp_df$pdsi_r),
      model_type = "PDSI",
      stringsAsFactors = FALSE
    ),
    data.frame(
      response = comp_df$response,
      method = comp_df$method,
      model = paste0(comp_df$best_vpd_predictor, " lag ", comp_df$best_vpd_lag),
      abs_r = abs(comp_df$best_vpd_r),
      model_type = "VPD",
      stringsAsFactors = FALSE
    )
  )

  p_comp <- ggplot(plot_df, aes(x = model, y = abs_r, fill = model_type)) +
    geom_col(width = 0.72) +
    facet_grid(response ~ method) +
    scale_fill_manual(values = c("PDSI" = COL_PDSI, "VPD" = COL_VPD)) +
    labs(
      title = ww_title("Deer Climate Comparison", "PDSI Lag 2 vs Best VPD Relationship"),
      subtitle = "Bars show absolute correlation strength (same response series and method)",
      x = "Model",
      y = "|correlation|",
      fill = "Driver"
    ) +
    theme_minimal(base_size = 11)
  ggsave(file.path(deer_dir, "pdsi_lag2_vs_vpd_comparison.png"), p_comp, width = 10.5, height = 6.5, dpi = 320)
}

# Regression-first panels requested: readable scatter + directional path + full stats.
build_driver_panel <- function(response_name, lag_table, data_in) {
  target_predictors <- c("RWI", "PDSI", "VPDMIN", "VPDMAX")
  target_predictors <- target_predictors[target_predictors %in% unique(lag_table$predictor)]
  if (length(target_predictors) == 0) return(NULL)

  rows <- list()
  stats_rows <- list()
  j <- 1L

  for (pred in target_predictors) {
    if (pred == "PDSI") {
      cand <- lag_table[
        lag_table$response == response_name & lag_table$method == "pearson" &
          lag_table$predictor == pred & lag_table$lag == 2,
        , drop = FALSE
      ]
      if (nrow(cand) == 0) {
        tmp <- lag_table[
          lag_table$response == response_name & lag_table$method == "pearson" & lag_table$predictor == pred,
          , drop = FALSE
        ]
        if (nrow(tmp) > 0) cand <- tmp[which.max(abs(tmp$r)), , drop = FALSE]
      }
    } else {
      tmp <- lag_table[
        lag_table$response == response_name & lag_table$method == "pearson" & lag_table$predictor == pred,
        , drop = FALSE
      ]
      cand <- if (nrow(tmp) > 0) tmp[which.max(abs(tmp$r)), , drop = FALSE] else tmp
    }

    if (nrow(cand) == 0) next

    lg <- cand$lag[1]
    x <- lag_vec(as.numeric(data_in[[pred]]), lg)
    y <- as.numeric(data_in[[response_name]])
    yr <- as.integer(data_in$year)
    ok <- is.finite(x) & is.finite(y) & is.finite(yr)
    if (sum(ok) < 4) next

    panel_name <- paste0(pred, " (lag ", lg, ")")

    rows[[j]] <- data.frame(
      year = yr[ok],
      predictor_value = x[ok],
      response_value = y[ok],
      panel = panel_name,
      predictor = pred,
      stringsAsFactors = FALSE
    )

    stats_rows[[j]] <- data.frame(
      response = response_name,
      predictor = pred,
      lag = lg,
      r = cand$r[1],
      r2 = cand$r2[1],
      p = cand$p[1],
      n = cand$n[1],
      panel = panel_name,
      label = paste0(
        "r = ", sprintf("%.3f", cand$r[1]),
        "\nR2 = ", sprintf("%.3f", cand$r2[1]),
        "\np = ", sprintf("%.4f", cand$p[1]),
        "\nn = ", cand$n[1]
      ),
      stringsAsFactors = FALSE
    )
    j <- j + 1L
  }

  if (length(rows) == 0) return(NULL)
  list(
    plot_df = do.call(rbind, rows),
    stats_df = do.call(rbind, stats_rows)
  )
}

for (resp in responses) {
  panel_obj <- build_driver_panel(resp, lag_df, deer_merge)
  if (is.null(panel_obj)) next

  plot_df <- panel_obj$plot_df
  stats_df <- panel_obj$stats_df

  write.csv(
    stats_df,
    file.path(deer_dir, paste0("deer_vpd_pdsi_regression_stats_", resp, ".csv")),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    stats_df,
    file.path(deer_dir, paste0("deer_rwi_vpd_pdsi_regression_stats_", resp, ".csv")),
    row.names = FALSE,
    na = ""
  )

  p_reg <- ggplot(plot_df, aes(x = predictor_value, y = response_value)) +
    geom_path(
      aes(group = panel),
      color = "gray70",
      alpha = 0.55,
      linewidth = 0.45,
      arrow = grid::arrow(type = "open", length = grid::unit(0.08, "inches"))
    ) +
    geom_point(aes(color = year), size = 2.5, alpha = 0.95) +
    geom_smooth(method = "lm", se = FALSE, color = COL_FIT, linewidth = 0.95) +
    facet_wrap(~panel, scales = "free_x") +
    scale_color_gradient(low = COL_PPT, high = COL_PDSI) +
    geom_text(
      data = stats_df,
      aes(label = label),
      x = -Inf,
      y = Inf,
      hjust = -0.05,
      vjust = 1.08,
      color = COL_FIT,
      size = 3.15,
      inherit.aes = FALSE
    ) +
    labs(
      title = ww_title("Deer vs RWI, VPD, and PDSI", paste0("Regression Panels (", resp, ")")),
      subtitle = "Scatter + regression line + directional year path; stats shown per panel (r, R2, p, n)",
      x = "Lagged predictor value",
      y = resp,
      color = "Year"
    ) +
    theme_minimal(base_size = 11)

  ggsave(
    file.path(deer_dir, paste0("deer_vpd_pdsi_regression_panels_", resp, ".png")),
    p_reg,
    width = 11.8,
    height = 6.8,
    dpi = 320
  )
  ggsave(
    file.path(deer_dir, paste0("deer_rwi_vpd_pdsi_regression_panels_", resp, ".png")),
    p_reg,
    width = 11.8,
    height = 6.8,
    dpi = 320
  )
}
# Run summary text for fast handoff.
summary_lines <- c(
  paste0("Run timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Input directory: ", normalizePath(input_dir, winslash = "/", mustWork = FALSE)),
  paste0("Output directory: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE)),
  "",
  paste0("Species .rwl files processed: ", length(species_files)),
  paste0("Deer AR prewhitening order (AIC-selected): ", ifelse(is.na(prew$order), "NA", prew$order)),
  paste0("PRISM source used: ", basename(prism_file)),
  paste0("PDSI source used: ", pdsi_source_file),
  paste0("Incomplete PRISM years excluded: ", if (length(prism_incomplete_years) == 0) "none" else paste(sort(unique(prism_incomplete_years)), collapse = ", ")),
  "",
  "Key outputs:",
  "- species/species_residual_chronologies_wide.csv",
  "- species/species_establishment_dieoff_faceted.png",
  "- species/species_establishment_dieoff_faceted_fullrange.png",
  "- species/species_active_series_through_time.png",
  "- species/species_active_series_through_time_fullrange.png",
  "- deer/deer_standardized_lag_correlations.csv (includes r2)",
  "- deer/deer_best_lag_scatter_r2_deer_resid.png",
  "- deer/climate_qc_flags.csv",
  "- deer/pdsi_lag2_vs_vpd_comparison.csv"
)
writeLines(summary_lines, con = file.path(output_dir, "run_summary.txt"))

message("Done. Outputs written to: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))



