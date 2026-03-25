# Wesselman Woods integrated handoff analysis
# Climate (PRISM + MODIS), deer browsing, and stand-age dynamics
# Built for reuse by future students.

if (dir.exists(".Rlibs")) {
  .libPaths(c(normalizePath(".Rlibs"), .libPaths()))
}

req_pkgs <- c("ggplot2", "dplyr", "dplR", "readxl", "zoo", "SPEI", "patchwork")
missing_pkgs <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(paste(
    "Missing required packages:", paste(missing_pkgs, collapse = ", "),
    "\nInstall these packages and re-run."
  ))
}

library(ggplot2)
library(dplyr)
library(dplR)
library(readxl)
library(zoo)
library(SPEI)
library(patchwork)

save_dpi <- 420
COL_RWI <- '#7F0000'
COL_DEER <- '#A67C52'
COL_VPD <- '#4169E1'
COL_SPEI <- '#DAA520'
COL_PPT <- '#0B3C8A'
COL_ETA <- '#2E8B57'
COL_NDVI <- '#39D353'
COL_PDSI <- '#CC7722'
COL_TEMP <- '#7FA9D9'
COL_FIT <- '#2F3E56'
COL_EVENT <- '#B22222'
COL_HEAT_NEG <- COL_PPT
COL_HEAT_POS <- COL_PDSI

theme_pub <- function(base_size = 15) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = 'bold'),
      plot.title.position = "plot",
      plot.margin = margin(14, 22, 18, 14),
      panel.spacing = grid::unit(0.9, "lines"),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10))
    )
}

save_standard_wavelet_plot <- function(years, values, file_out, main_txt, subtitle_txt, crn_label = "RWI") {
  keep <- is.finite(years) & is.finite(values)
  years <- years[keep]
  values <- values[keep]

  if (length(values) < 32) {
    return(FALSE)
  }

  wave_obj <- dplR::morlet(y1 = values, x1 = years, p2 = NULL, dj = 0.25, siglvl = 0.95)

  png(file_out, width = 1600, height = 1200, res = 220)
  par(oma = c(0, 0, 4.8, 0))
  dplR::wavelet.plot(
    wave_obj,
    useRaster = NA,
    x.lab = "Year",
    period.lab = "Period (years)",
    crn.lab = crn_label
  )
  mtext(main_txt, side = 3, outer = TRUE, line = 2.0, cex = 1.3, font = 2)
  mtext(subtitle_txt, side = 3, outer = TRUE, line = 0.8, cex = 0.84)
  dev.off()

  TRUE
}

# ==========================
# USER SETTINGS
# ==========================
short_climate_file <- "pca_climate_merged.csv"          # MODIS-era table (2000+)
short_climate_fallback_file <- "chronology_climate.xlsx"
pdsi_file <- "pdsi.csv"
prism_file <- "PRISM_ppt_tmin_tmean_tmax_tdmean_vpdmin_vpdmax_stable_4km_189501_202401_37.9810_-87.5064.csv"
use_latest_prism_file <- FALSE
prism_pattern <- "^PRISM_.*stable_4km_.*\\.csv$"
prism_latitude <- 37.9810

rwl_long_file <- "dated8.rwl"
rwl_bulk_file <- "BulkClean.rwl"
deer_file <- "wwn_deer.xlsx"
deer_sheet <- 1
ndvi_file <- "ee-chart.csv"
enso_file <- "enso_nino34_noaa.txt"
enso_annual_file <- "enso_annual_nino34.csv"

response_var <- "RWI"
year_var <- "year"
max_lag <- 2
sampling_year <- 2024

hydro_marker_quantile <- 0.85
documented_hydro_events_file <- 'documented_flood_drought_years.csv'

out_dir <- "wwn_handoff_outputs"
plot_dir <- file.path(out_dir, "plots")
dir.create(out_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

# ==========================
# HELPER FUNCTIONS
# ==========================

ww_title <- function(subject, analysis) {
  paste0("Wesselman Woods - ", subject, " (", analysis, ")")
}

canon_name <- function(x) {
  y <- tolower(x)
  y <- trimws(gsub("\ufeff", "", y, fixed = TRUE))
  y <- gsub("\\s*\\(.*\\)$", "", y)
  y <- gsub("[^a-z0-9]", "", y)
  y
}

pick_col <- function(clean_names, candidates) {
  for (nm in candidates) {
    idx <- which(clean_names == nm)
    if (length(idx) > 0) return(idx[1])
  }
  NA_integer_
}

read_table_by_ext <- function(path, excel_sheet = 1L) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    return(as.data.frame(read_excel(path, sheet = excel_sheet), stringsAsFactors = FALSE, check.names = FALSE))
  }
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

load_pdsi_annual <- function(pdsi_path, fallback_path = NULL, fallback_sheet = 1L) {
  parse_pdsi <- function(df, year_candidates, pdsi_candidates) {
    nm <- canon_name(names(df))
    year_idx <- pick_col(nm, year_candidates)
    pdsi_idx <- pick_col(nm, pdsi_candidates)
    if (is.na(year_idx) || is.na(pdsi_idx)) return(NULL)

    yr_raw <- as.character(df[[year_idx]])
yr <- suppressWarnings(as.integer(yr_raw))
yr_digits <- gsub("[^0-9]", "", yr_raw)
yr4 <- suppressWarnings(as.integer(substr(yr_digits, 1, 4)))
use_yr4 <- is.finite(yr4) & (!is.finite(yr) | yr < 1000 | yr > 3000)
yr[use_yr4] <- yr4[use_yr4]

    pdsi <- suppressWarnings(as.numeric(df[[pdsi_idx]]))
    out <- data.frame(Year = yr, PDSI = pdsi, stringsAsFactors = FALSE)
    out <- out[is.finite(out$Year), , drop = FALSE]
    if (nrow(out) == 0) return(NULL)
    out %>%
      group_by(Year) %>%
      summarise(PDSI = if (all(is.na(PDSI))) NA_real_ else mean(PDSI, na.rm = TRUE), .groups = "drop")
  }

  result <- list(data = NULL, source = "not found")
  if (file.exists(pdsi_path)) {
    pd <- read.csv(pdsi_path, check.names = FALSE, stringsAsFactors = FALSE)
    parsed <- parse_pdsi(pd, c("date", "year", "yr"), c("pdsiavg", "pdsi"))
    if (!is.null(parsed)) {
      result$data <- parsed
      result$source <- basename(pdsi_path)
      return(result)
    }
  }

  if (!is.null(fallback_path) && file.exists(fallback_path)) {
    fb <- read_table_by_ext(fallback_path, excel_sheet = fallback_sheet)
    parsed <- parse_pdsi(fb, c("year", "yr", "date"), c("pdsi", "pdsiavg"))
    if (!is.null(parsed)) {
      result$data <- parsed
      result$source <- basename(fallback_path)
      return(result)
    }
  }

  result
}

zfun <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

build_label_positions <- function(df, group_col, x_col, y_col, x_frac = 0.03, y_frac = 0.05) {
  groups <- split(df, df[[group_col]])
  bind_rows(lapply(groups, function(dfp) {
    x_rng <- range(dfp[[x_col]], na.rm = TRUE)
    y_rng <- range(dfp[[y_col]], na.rm = TRUE)

    x_span <- diff(x_rng)
    y_span <- diff(y_rng)
    if (!is.finite(x_span) || x_span == 0) x_span <- max(0.1, abs(x_rng[1]) * 0.08, abs(x_rng[2]) * 0.08, na.rm = TRUE)
    if (!is.finite(y_span) || y_span == 0) y_span <- max(0.1, abs(y_rng[1]) * 0.08, abs(y_rng[2]) * 0.08, na.rm = TRUE)

    data.frame(
      group_id = unique(dfp[[group_col]])[1],
      label_x = x_rng[2] - x_frac * x_span,
      label_y = y_rng[2] - y_frac * y_span,
      stringsAsFactors = FALSE
    )
  }))
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

cor_ci_fisher <- function(r, n, alpha = 0.05) {
  if (!is.finite(r) || !is.finite(n) || n <= 3) return(c(NA_real_, NA_real_))
  r <- max(min(r, 0.999999), -0.999999)
  z <- 0.5 * log((1 + r) / (1 - r))
  se <- 1 / sqrt(n - 3)
  zc <- qnorm(1 - alpha / 2)
  lo <- z - zc * se
  hi <- z + zc * se
  c((exp(2 * lo) - 1) / (exp(2 * lo) + 1), (exp(2 * hi) - 1) / (exp(2 * hi) + 1))
}
sig_code <- function(p) {
  out <- rep("", length(p))
  out[!is.na(p) & p < 0.10] <- "."
  out[!is.na(p) & p < 0.05] <- "*"
  out[!is.na(p) & p < 0.01] <- "**"
  out[!is.na(p) & p < 0.001] <- "***"
  out
}

heatmap_sig_label <- function(p) {
  out <- sig_code(p)
  out[is.na(p)] <- ""
  out
}

fmt_p <- function(p) {
  p <- as.numeric(p)
  out <- rep(NA_character_, length(p))
  out[is.na(p)] <- "p = NA"
  out[!is.na(p) & p < 0.001] <- "p < 0.001"
  keep <- !is.na(p) & p >= 0.001
  out[keep] <- paste0("p = ", format(round(p[keep], 3), nsmall = 3))
  out
}

safe_min <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

safe_cor <- function(x, y, method = "pearson") {
  keep <- complete.cases(x, y)
  n <- sum(keep)
  if (n < 4) {
    return(data.frame(r = NA_real_, r2 = NA_real_, p = NA_real_, n = n))
  }
  xk <- x[keep]
  yk <- y[keep]
  if (sd(xk) == 0 || sd(yk) == 0) {
    return(data.frame(r = NA_real_, r2 = NA_real_, p = NA_real_, n = n))
  }
  ct <- suppressWarnings(cor.test(xk, yk, method = method))
  r_val <- unname(ct$estimate)
  data.frame(r = r_val, r2 = r_val^2, p = ct$p.value, n = n)
}

make_qc_row <- function(check, n_flagged = 0L, years = integer(0), note = "", action = "") {
  data.frame(
    check = check,
    n_flagged = as.integer(n_flagged),
    years = if (length(years) == 0) "" else paste(sort(unique(years)), collapse = ", "),
    note = note,
    action = action,
    stringsAsFactors = FALSE
  )
}

safe_read_rwl <- function(path) {
  x <- try(read.tucson(path), silent = TRUE)
  if (!inherits(x, "try-error")) return(x)

  message("read.tucson failed for ", path, "; attempting fallback parser.")
  lines <- readLines(path, warn = FALSE)
  if (length(lines) == 0) stop("RWL file is empty: ", path)

  tok <- strsplit(trimws(lines), "\\s+")
  rows <- vector("list", length(tok))
  k <- 1
  for (tt in tok) {
    if (length(tt) < 3) next
    id <- tt[1]
    decade <- suppressWarnings(as.integer(tt[2]))
    vals <- suppressWarnings(as.numeric(tt[-c(1, 2)]))
    if (is.na(decade) || length(vals) == 0) next
    yrs <- decade + seq_along(vals) - 1
    vals[vals == 999] <- NA_real_
    rows[[k]] <- data.frame(series = id, year = yrs, rw = vals, stringsAsFactors = FALSE)
    k <- k + 1
  }

  long <- bind_rows(rows)
  if (nrow(long) == 0) stop("Fallback parser produced no rows for: ", path)

  long <- long %>%
    filter(!is.na(year)) %>%
    group_by(series, year) %>%
    summarise(rw = if (all(is.na(rw))) NA_real_ else mean(rw, na.rm = TRUE), .groups = "drop")

  years <- sort(unique(long$year))
  series <- sort(unique(long$series))
  mat <- matrix(NA_real_, nrow = length(years), ncol = length(series),
                dimnames = list(as.character(years), series))
  rr <- match(long$year, years)
  cc <- match(long$series, series)
  mat[cbind(rr, cc)] <- long$rw
  as.data.frame(mat, stringsAsFactors = FALSE)
}

clean_rwl_for_detrend <- function(rwl, min_obs = 8) {
  rwl_num <- as.data.frame(lapply(rwl, function(x) suppressWarnings(as.numeric(x))))
  rownames(rwl_num) <- rownames(rwl)
  rwl_num[rwl_num <= 0] <- NA
  keep <- colSums(!is.na(rwl_num)) >= min_obs
  rwl_num <- rwl_num[, keep, drop = FALSE]
  if (ncol(rwl_num) == 0) stop("No usable ring-width series after filtering.")
  rwl_num
}
make_lag_corr <- function(df, year_col, response, predictors, max_lag = 2, method = "pearson") {
  d <- df[order(df[[year_col]]), ]
  rows <- list()
  k <- 1
  for (v in predictors) {
    for (lg in 0:max_lag) {
      xv <- if (lg == 0) d[[v]] else dplyr::lag(d[[v]], n = lg)
      ct <- safe_cor(d[[response]], xv, method = method)
      rows[[k]] <- data.frame(
        variable = v,
        lag = lg,
        method = method,
        r = ct$r,
        r2 = ct$r2,
        p = ct$p,
        n = ct$n,
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }
  bind_rows(rows)
}

make_matrix_corr <- function(df, vars, method = "pearson") {
  rows <- list()
  k <- 1
  for (v1 in vars) {
    for (v2 in vars) {
      ct <- safe_cor(df[[v1]], df[[v2]], method = method)
      rows[[k]] <- data.frame(
        var_x = v1,
        var_y = v2,
        method = method,
        r = ct$r,
        r2 = ct$r2,
        p = ct$p,
        n = ct$n,
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }
  bind_rows(rows)
}

plot_lag_heatmap <- function(tab, title_txt, subtitle_txt, file_out) {
  d <- tab %>%
    mutate(
      lag_f = paste0("Lag ", lag),
      variable = factor(variable, levels = rev(unique(variable))),
      label = heatmap_sig_label(p)
    )

  p <- ggplot(d, aes(lag_f, variable, fill = r)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = label), size = 3.2, fontface = "bold") +
    scale_fill_gradient2(low = COL_HEAT_NEG, mid = 'white', high = COL_HEAT_POS, midpoint = 0, na.value = 'gray90') +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "Predictor lag",
      y = "Climate variable",
      fill = "Correlation"
    ) +
    theme_pub(base_size = 15) +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1),
      axis.text.y = element_text(margin = margin(r = 6))
    )

  plot_width <- max(10.8, 0.78 * length(unique(d$variable)) + 7.2)
  plot_height <- max(6.4, 0.62 * length(unique(d$variable)) + 2.8)
  ggsave(file_out, p, width = plot_width, height = plot_height, dpi = save_dpi)
}

plot_matrix_heatmap <- function(tab, title_txt, subtitle_txt, file_out) {
  d <- tab %>%
    mutate(
      var_x = factor(var_x, levels = unique(var_x)),
      var_y = factor(var_y, levels = rev(unique(var_y))),
      label = heatmap_sig_label(p)
    )

  p <- ggplot(d, aes(var_x, var_y, fill = r)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = label), size = 3.0, fontface = "bold") +
    scale_fill_gradient2(low = COL_HEAT_NEG, mid = 'white', high = COL_HEAT_POS, midpoint = 0, na.value = 'gray90') +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = NULL,
      y = NULL,
      fill = "Correlation"
    ) +
    theme_pub(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      axis.text.y = element_text(margin = margin(r = 6))
    )

  plot_width <- max(9.6, 0.68 * length(unique(d$var_x)) + 6.2)
  plot_height <- max(7.2, 0.58 * length(unique(d$var_y)) + 3.2)
  ggsave(file_out, p, width = plot_width, height = plot_height, dpi = save_dpi)
}

climate_driver_colors <- function(vars) {
  base_cols <- c(
    "RWI" = COL_RWI,
    "RWI_std" = COL_RWI,
    "RWI_res" = COL_RWI,
    "PPT" = COL_PPT,
    "PPT_mm" = COL_PPT,
    "VPDMIN" = COL_VPD,
    "VPDMAX" = COL_VPD,
    "SPEI" = COL_SPEI,
    "SPEI6_gs" = COL_SPEI,
    "SPEI6_min" = COL_SPEI,
    "ENSO_NINO34" = COL_SPEI,
    "PDSI" = COL_PDSI,
    "TMEAN" = COL_TEMP,
    "TMEAN_C" = COL_TEMP,
    "TMIN_C" = COL_TEMP,
    "TMAX_C" = COL_TEMP,
    "EVT" = COL_ETA,
    "ETA" = COL_ETA,
    "NDVI" = COL_NDVI,
    "ClimateComposite" = COL_DEER
  )
  out <- setNames(rep(COL_VPD, length(vars)), vars)
  hit <- intersect(names(base_cols), vars)
  out[hit] <- base_cols[hit]
  out
}

fmt_num <- function(x, digits = 3) {
  if (is.na(x)) return("NA")
  format(round(x, digits), nsmall = digits)
}

corr_symbol <- function(method) {
  ifelse(tolower(method) == "spearman", "rho", "r")
}

corr_power_label <- function(method) {
  ifelse(tolower(method) == "spearman", "rho^2", "R2")
}

build_corr_label <- function(est, est2, p, n, method) {
  method <- tolower(method)
  paste0(
    corr_symbol(method), " = ", vapply(est, fmt_num, character(1)),
    "\n", corr_power_label(method), " = ", vapply(est2, fmt_num, character(1)),
    "\n", vapply(p, fmt_p, character(1)),
    "\nn = ", n
  )
}

build_best_lag_regression_panels <- function(
  df,
  year_col,
  response,
  predictors,
  lag_table,
  title_txt,
  subtitle_txt,
  file_out,
  stats_out,
  ncol = 3,
  lag_selection_method = "pearson",
  stat_method = "pearson"
) {
  d <- df[order(df[[year_col]]), , drop = FALSE]
  lag_src <- lag_table[lag_table$method == lag_selection_method & lag_table$variable %in% predictors, , drop = FALSE]

  if (nrow(lag_src) == 0) return(invisible(NULL))

  panel_rows <- list()
  stats_rows <- list()
  k <- 1

  for (v in predictors) {
    vv <- lag_src[lag_src$variable == v & is.finite(lag_src$r), , drop = FALSE]
    if (nrow(vv) == 0) next

    best <- vv[which.max(abs(vv$r)), , drop = FALSE]
    lg <- as.integer(best$lag[1])

    xv <- if (lg == 0) d[[v]] else dplyr::lag(d[[v]], n = lg)
    yv <- d[[response]]
    yr <- d[[year_col]]

    keep <- complete.cases(xv, yv, yr)
    if (sum(keep) < 4) next

    panel_name <- paste0(v, " (lag ", lg, ")")

    cp <- safe_cor(yv, xv, method = "pearson")
    cs <- safe_cor(yv, xv, method = "spearman")
    primary_stats <- if (tolower(stat_method) == "spearman") cs else cp

    panel_rows[[k]] <- data.frame(
      year = yr[keep],
      variable = v,
      panel = panel_name,
      predictor = xv[keep],
      response_value = yv[keep],
      stringsAsFactors = FALSE
    )

    stats_rows[[k]] <- data.frame(
      variable = v,
      lag = lg,
      method = stat_method,
      lag_selection_method = lag_selection_method,
      r = primary_stats$r,
      r2 = primary_stats$r2,
      p = primary_stats$p,
      n = primary_stats$n,
      pearson_r = cp$r,
      pearson_r2 = cp$r2,
      pearson_p = cp$p,
      spearman_rho = cs$r,
      spearman_rho2 = cs$r2,
      spearman_p = cs$p,
      panel = panel_name,
      stringsAsFactors = FALSE
    )
    k <- k + 1
  }

  if (length(panel_rows) == 0) return(invisible(NULL))

  panel_df <- bind_rows(panel_rows)
  stats_df <- bind_rows(stats_rows)

  stats_df$label <- build_corr_label(stats_df$r, stats_df$r2, stats_df$p, stats_df$n, stats_df$method)
  label_pos <- build_label_positions(panel_df, "panel", "predictor", "response_value", x_frac = 0.03, y_frac = 0.07)
  names(label_pos)[names(label_pos) == "group_id"] <- "panel"
  stats_df <- left_join(stats_df, label_pos, by = "panel")

  write.csv(stats_df, stats_out, row.names = FALSE)

  col_map <- climate_driver_colors(unique(panel_df$variable))

  p <- ggplot(panel_df, aes(predictor, response_value)) +
    geom_path(
      aes(group = panel),
      color = "gray70",
      alpha = 0.55,
      linewidth = 0.45,
      arrow = grid::arrow(type = "open", length = grid::unit(0.08, "inches"))
    ) +
    geom_point(aes(color = variable), alpha = 0.95, size = 2.5) +
    geom_smooth(method = "lm", se = FALSE, color = COL_FIT, linewidth = 0.95) +
    facet_wrap(~panel, scales = "free_x", ncol = ncol) +
    scale_color_manual(values = col_map) +
    geom_text(
      data = stats_df,
      aes(x = label_x, y = label_y, label = label),
      hjust = 1,
      vjust = 1,
      inherit.aes = FALSE,
      color = COL_FIT,
      size = 2.9,
      lineheight = 0.94
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "Lagged predictor value",
      y = response,
      color = "Predictor"
    ) +
    theme_pub(base_size = 15) +
    theme(strip.text = element_text(face = "bold"))

  plot_width <- max(14.6, 3.55 * min(ncol, nrow(stats_df)) + 2.2)
  plot_height <- max(8.8, 2.7 * ceiling(nrow(stats_df) / ncol))
  ggsave(file_out, p, width = plot_width, height = plot_height, dpi = save_dpi)

  invisible(stats_df)
}

make_marker_bands <- function(years, half_width = 0.5) {
  yrs <- sort(unique(as.integer(years[is.finite(years)])))
  if (length(yrs) == 0) {
    return(data.frame(xmin = numeric(0), xmax = numeric(0), marker_year = integer(0)))
  }
  data.frame(
    xmin = yrs - half_width,
    xmax = yrs + half_width,
    marker_year = yrs,
    stringsAsFactors = FALSE
  )
}

plot_faceted_timeseries_with_markers <- function(
  df,
  year_col,
  variables,
  title_txt,
  subtitle_txt,
  file_out,
  marker_years = integer(0),
  ncol = 2
) {
  vars <- variables[variables %in% names(df)]
  if (length(vars) == 0) return(invisible(NULL))

  ts_rows <- bind_rows(lapply(vars, function(v) {
    data.frame(
      Year = df[[year_col]],
      variable = v,
      value = as.numeric(df[[v]]),
      stringsAsFactors = FALSE
    )
  }))

  ts_rows <- ts_rows[is.finite(ts_rows$Year), , drop = FALSE]
  if (nrow(ts_rows) == 0) return(invisible(NULL))

  band_df <- make_marker_bands(intersect(marker_years, ts_rows$Year))
  col_map <- climate_driver_colors(vars)

  p <- ggplot(ts_rows, aes(Year, value, color = variable)) +
    geom_rect(
      data = band_df,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "grey80",
      alpha = 0.35
    ) +
    geom_line(linewidth = 0.95, na.rm = TRUE) +
    geom_point(size = 1.1, alpha = 0.75, na.rm = TRUE) +
    facet_wrap(~variable, scales = "free_y", ncol = ncol) +
    scale_color_manual(values = col_map, guide = "none") +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "Year",
      y = "Observed value"
    ) +
    theme_pub(base_size = 15) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  plot_height <- max(8.0, 2.5 * ceiling(length(vars) / ncol))
  ggsave(file_out, p, width = 13.8, height = plot_height, dpi = save_dpi)

  invisible(ts_rows)
}

build_all_lag_regression_panels <- function(
  df,
  year_col,
  response,
  predictors,
  max_lag,
  title_txt,
  subtitle_txt,
  file_out,
  stats_out,
  ncol = 4,
  stat_method = "pearson"
) {
  d <- df[order(df[[year_col]]), , drop = FALSE]
  panel_rows <- list()
  stats_rows <- list()
  k <- 1L

  for (v in predictors) {
    if (!(v %in% names(d))) next
    for (lg in 0:max_lag) {
      xv <- if (lg == 0) d[[v]] else dplyr::lag(d[[v]], n = lg)
      yv <- d[[response]]
      yr <- d[[year_col]]

      keep <- complete.cases(xv, yv, yr)
      if (sum(keep) < 4) next

      cp <- safe_cor(yv, xv, method = "pearson")
      cs <- safe_cor(yv, xv, method = "spearman")
      primary_stats <- if (tolower(stat_method) == "spearman") cs else cp
      panel_name <- paste0(v, " (lag ", lg, ")")

      panel_rows[[k]] <- data.frame(
        year = yr[keep],
        variable = v,
        lag = lg,
        panel = panel_name,
        predictor = xv[keep],
        response_value = yv[keep],
        stringsAsFactors = FALSE
      )

      stats_rows[[k]] <- data.frame(
        response = response,
        variable = v,
        lag = lg,
        method = stat_method,
        r = primary_stats$r,
        r2 = primary_stats$r2,
        p = primary_stats$p,
        n = primary_stats$n,
        pearson_r = cp$r,
        pearson_r2 = cp$r2,
        pearson_p = cp$p,
        spearman_rho = cs$r,
        spearman_rho2 = cs$r2,
        spearman_p = cs$p,
        panel = panel_name,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }

  if (length(panel_rows) == 0) return(invisible(NULL))

  panel_df <- bind_rows(panel_rows)
  stats_df <- bind_rows(stats_rows)
  stats_df$label <- build_corr_label(stats_df$r, stats_df$r2, stats_df$p, stats_df$n, stats_df$method)
  label_pos <- build_label_positions(panel_df, "panel", "predictor", "response_value", x_frac = 0.03, y_frac = 0.07)
  names(label_pos)[names(label_pos) == "group_id"] <- "panel"
  stats_df <- left_join(stats_df, label_pos, by = "panel")

  write.csv(stats_df, stats_out, row.names = FALSE)

  col_map <- climate_driver_colors(unique(panel_df$variable))
  p <- ggplot(panel_df, aes(predictor, response_value)) +
    geom_path(
      aes(group = panel),
      color = "gray72",
      alpha = 0.45,
      linewidth = 0.42,
      arrow = grid::arrow(type = "open", length = grid::unit(0.07, "inches"))
    ) +
    geom_point(aes(color = variable), alpha = 0.92, size = 2.1) +
    geom_smooth(method = "lm", se = FALSE, color = COL_FIT, linewidth = 0.9) +
    facet_wrap(~panel, scales = "free_x", ncol = ncol) +
    scale_color_manual(values = col_map) +
    geom_text(
      data = stats_df,
      aes(x = label_x, y = label_y, label = label),
      hjust = 1,
      vjust = 1,
      inherit.aes = FALSE,
      color = COL_FIT,
      size = 2.55,
      lineheight = 0.94
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "Lagged predictor value",
      y = response,
      color = "Predictor"
    ) +
    theme_pub(base_size = 14) +
    theme(strip.text = element_text(face = "bold"))

  plot_height <- max(9.4, 2.45 * ceiling(nrow(stats_df) / ncol))
  plot_width <- max(16.2, 3.4 * min(ncol, nrow(stats_df)) + 2.4)
  ggsave(file_out, p, width = plot_width, height = plot_height, dpi = 320)

  invisible(stats_df)
}
# ==========================
# INPUT DISCOVERY
# ==========================

if (use_latest_prism_file || prism_file == "") {
  prism_candidates <- list.files(".", pattern = prism_pattern, full.names = TRUE)
  if (length(prism_candidates) == 0) {
    stop("No PRISM files found with pattern: ", prism_pattern)
  }
  prism_mtime <- file.info(prism_candidates)$mtime
  prism_file <- prism_candidates[which.max(prism_mtime)]
}

short_climate_active_file <- short_climate_file
if (!file.exists(short_climate_active_file) && file.exists(short_climate_fallback_file)) {
  short_climate_active_file <- short_climate_fallback_file
}
if (!file.exists(short_climate_active_file)) {
  stop("Missing short climate file. Tried: ", short_climate_file, " and ", short_climate_fallback_file)
}
if (!file.exists(prism_file)) stop("Missing file: ", prism_file)
if (!file.exists(rwl_long_file)) stop("Missing file: ", rwl_long_file)
if (!file.exists(rwl_bulk_file)) stop("Missing file: ", rwl_bulk_file)
if (!file.exists(deer_file)) stop("Missing file: ", deer_file)

pdsi_loaded <- load_pdsi_annual(pdsi_file, fallback_path = short_climate_active_file, fallback_sheet = 1L)
pdsi_ann <- pdsi_loaded$data
pdsi_source_file <- pdsi_loaded$source
if (is.null(pdsi_ann)) {
  warning("No usable PDSI column found in pdsi.csv or fallback short climate table.")
}
ndvi_source_file <- "not found"

# ==========================
# A) MODIS-ERA SHORT CLIMATE TABLE (2000+)
# ==========================

short_df <- read_table_by_ext(short_climate_active_file, excel_sheet = 1L)
short_names <- trimws(names(short_df))
names(short_df) <- short_names
short_file_used <- basename(short_climate_active_file)

if (!(year_var %in% names(short_df))) {
  y_idx <- pick_col(canon_name(names(short_df)), c("year", "yr"))
  if (is.na(y_idx)) stop("Could not find year column in short climate file.")
  names(short_df)[y_idx] <- year_var
}
if (!(response_var %in% names(short_df))) stop("Could not find RWI column in short climate file.")

short_df[[year_var]] <- as.integer(short_df[[year_var]])
short_df <- short_df[order(short_df[[year_var]]), ]

short_names_canon <- canon_name(names(short_df))
evt_idx <- pick_col(short_names_canon, c("evt", "eta", "eta_mm", "et", "actualevapotranspiration"))
if (!is.na(evt_idx) && !"EVT" %in% names(short_df)) {
  names(short_df)[evt_idx] <- "EVT"
}

if (!is.null(pdsi_ann)) {
  pdsi_short <- pdsi_ann %>%
    transmute(!!year_var := as.integer(Year), PDSI_csv = PDSI)
  short_df <- left_join(short_df, pdsi_short, by = year_var)
  if ("PDSI" %in% names(short_df)) {
    short_df$PDSI <- dplyr::coalesce(short_df$PDSI_csv, short_df$PDSI)
  } else {
    short_df$PDSI <- short_df$PDSI_csv
  }
  short_df$PDSI_csv <- NULL
}

if (file.exists(ndvi_file)) {
  ndvi_raw <- read.csv(ndvi_file, check.names = FALSE, stringsAsFactors = FALSE)
  ndvi_names <- canon_name(names(ndvi_raw))
  ndvi_year_idx <- pick_col(ndvi_names, c("year", "yr"))
  ndvi_val_idx <- pick_col(ndvi_names, c("ndvi", "meanndvi", "ndviavg"))
  if (!is.na(ndvi_year_idx) && !is.na(ndvi_val_idx)) {
    ndvi_ann <- data.frame(
      year = as.integer(suppressWarnings(as.numeric(ndvi_raw[[ndvi_year_idx]]))),
      NDVI_file = as.numeric(ndvi_raw[[ndvi_val_idx]]),
      stringsAsFactors = FALSE
    )
    ndvi_ann <- ndvi_ann[is.finite(ndvi_ann$year), , drop = FALSE]
    ndvi_ann <- ndvi_ann %>%
      group_by(year) %>%
      summarise(NDVI_file = if (all(is.na(NDVI_file))) NA_real_ else mean(NDVI_file, na.rm = TRUE), .groups = "drop")
    names(ndvi_ann)[1] <- year_var
    short_df <- left_join(short_df, ndvi_ann, by = year_var)
    if ("NDVI" %in% names(short_df)) {
      short_df$NDVI <- dplyr::coalesce(short_df$NDVI_file, short_df$NDVI)
    } else {
      short_df$NDVI <- short_df$NDVI_file
    }
    short_df$NDVI_file <- NULL
    ndvi_source_file <- basename(ndvi_file)
  }
}

short_predictors <- intersect(c("EVT", "PPT", "TMEAN", "PDSI", "SPEI", "NDVI"), names(short_df))
if (length(short_predictors) == 0) stop("No expected short-period climate predictors found.")

# ==========================
# B) PRISM MONTHLY -> ANNUAL TABLE (1895+)
# ==========================

prism_lines <- readLines(prism_file, warn = FALSE)
header_idx <- which(grepl("^\\s*\\ufeff?[Dd]ate\\s*,", prism_lines))[1]
if (is.na(header_idx)) {
  prism_line_norm <- tolower(gsub("\"", "", prism_lines))
  header_idx <- which(grepl("date\\s*,\\s*ppt\\s*,\\s*tmin\\s*,\\s*tmean\\s*,\\s*tmax", prism_line_norm))[1]
}
if (is.na(header_idx)) stop("Could not find PRISM header row containing date,ppt,tmin,tmean,tmax.")

prism_raw <- read.csv(prism_file, skip = header_idx - 1, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(prism_raw) < 5) stop("PRISM file has too few columns.")

prism_clean_names <- canon_name(names(prism_raw))

date_idx <- pick_col(prism_clean_names, c("date"))
if (is.na(date_idx)) date_idx <- 1
ppt_idx <- pick_col(prism_clean_names, c("ppt"))
tmin_idx <- pick_col(prism_clean_names, c("tmin"))
tmean_idx <- pick_col(prism_clean_names, c("tmean"))
tmax_idx <- pick_col(prism_clean_names, c("tmax"))
tdmean_idx <- pick_col(prism_clean_names, c("tdmean", "tedmean"))
vpdmin_idx <- pick_col(prism_clean_names, c("vpdmin"))
vpdmax_idx <- pick_col(prism_clean_names, c("vpdmax"))

needed <- c(ppt_idx, tmin_idx, tmean_idx, tmax_idx)
if (any(is.na(needed))) {
  stop("PRISM parser could not find one or more required columns: Date, ppt, tmin, tmean, tmax. Found columns: ", paste(names(prism_raw), collapse = ", "))
}

to_num <- function(x) suppressWarnings(as.numeric(x))
prism <- data.frame(
  Date = as.character(prism_raw[[date_idx]]),
  PPT_in = to_num(prism_raw[[ppt_idx]]),
  TMIN_F = to_num(prism_raw[[tmin_idx]]),
  TMEAN_F = to_num(prism_raw[[tmean_idx]]),
  TMAX_F = to_num(prism_raw[[tmax_idx]]),
  stringsAsFactors = FALSE
)

prism$TDMEAN_F <- if (!is.na(tdmean_idx)) to_num(prism_raw[[tdmean_idx]]) else NA_real_
prism$VPDMIN <- if (!is.na(vpdmin_idx)) to_num(prism_raw[[vpdmin_idx]]) else NA_real_
prism$VPDMAX <- if (!is.na(vpdmax_idx)) to_num(prism_raw[[vpdmax_idx]]) else NA_real_

prism$Date <- trimws(prism$Date)
prism$Date <- sub("^([0-9]{4})-([0-9]{2})$", "\\1-\\2-01", prism$Date)
prism$Date <- sub("^([0-9]{4})([0-9]{2})$", "\\1-\\2-01", prism$Date)
prism$Date <- as.Date(prism$Date)
prism <- prism[!is.na(prism$Date), ]
prism <- prism[order(prism$Date), ]

prism$Year <- as.integer(format(prism$Date, "%Y"))
prism$Month <- as.integer(format(prism$Date, "%m"))

# Unit conversions for physically meaningful SPEI input.
prism$PPT_mm <- prism$PPT_in * 25.4
prism$TMIN_C <- (prism$TMIN_F - 32) * 5 / 9
prism$TMEAN_C <- (prism$TMEAN_F - 32) * 5 / 9
prism$TMAX_C <- (prism$TMAX_F - 32) * 5 / 9

prism$SPEI6 <- NA_real_
if (all(!is.na(prism$PPT_mm)) && all(!is.na(prism$TMEAN_C))) {
  start_vec <- c(prism$Year[1], prism$Month[1])
  ppt_ts <- ts(prism$PPT_mm, start = start_vec, frequency = 12)
  tmean_ts <- ts(prism$TMEAN_C, start = start_vec, frequency = 12)

  pet_ts <- tryCatch(
    SPEI::thornthwaite(Tave = tmean_ts, lat = prism_latitude),
    error = function(e) rep(NA_real_, length(ppt_ts))
  )

  if (!all(is.na(pet_ts))) {
    wb_ts <- ppt_ts - pet_ts
    spei_obj <- tryCatch(SPEI::spei(wb_ts, scale = 6), error = function(e) NULL)
    if (!is.null(spei_obj)) {
      prism$SPEI6 <- as.numeric(spei_obj$fitted)
    }
  }
}

prism_qc_rows <- list(
  make_qc_row(
    "Negative monthly PPT_mm",
    sum(is.finite(prism$PPT_mm) & prism$PPT_mm < 0),
    prism$Year[is.finite(prism$PPT_mm) & prism$PPT_mm < 0],
    "Physical precipitation values should not be negative",
    if (sum(is.finite(prism$PPT_mm) & prism$PPT_mm < 0) > 0) "Inspect source values" else "No issue found"
  ),
  make_qc_row(
    "Negative monthly VPDMIN",
    sum(is.finite(prism$VPDMIN) & prism$VPDMIN < 0),
    prism$Year[is.finite(prism$VPDMIN) & prism$VPDMIN < 0],
    "Vapor pressure deficit should not be negative",
    if (sum(is.finite(prism$VPDMIN) & prism$VPDMIN < 0) > 0) "Inspect source values" else "No issue found"
  ),
  make_qc_row(
    "Negative monthly VPDMAX",
    sum(is.finite(prism$VPDMAX) & prism$VPDMAX < 0),
    prism$Year[is.finite(prism$VPDMAX) & prism$VPDMAX < 0],
    "Vapor pressure deficit should not be negative",
    if (sum(is.finite(prism$VPDMAX) & prism$VPDMAX < 0) > 0) "Inspect source values" else "No issue found"
  ),
  make_qc_row(
    "Monthly rows with VPDMIN > VPDMAX",
    sum(is.finite(prism$VPDMIN) & is.finite(prism$VPDMAX) & prism$VPDMIN > prism$VPDMAX),
    prism$Year[is.finite(prism$VPDMIN) & is.finite(prism$VPDMAX) & prism$VPDMIN > prism$VPDMAX],
    "Minimum VPD should not exceed maximum VPD",
    if (sum(is.finite(prism$VPDMIN) & is.finite(prism$VPDMAX) & prism$VPDMIN > prism$VPDMAX) > 0) "Inspect source values" else "No issue found"
  )
)

annual_prism_full <- prism %>%
  group_by(Year) %>%
  summarise(
    PRISM_month_n = n_distinct(Month),
    PPT_mm = sum(PPT_mm, na.rm = TRUE),
    TMIN_C = mean(TMIN_C, na.rm = TRUE),
    TMEAN_C = mean(TMEAN_C, na.rm = TRUE),
    TMAX_C = mean(TMAX_C, na.rm = TRUE),
    VPDMIN = mean(VPDMIN, na.rm = TRUE),
    VPDMAX = mean(VPDMAX, na.rm = TRUE),
    SPEI6_gs = mean(SPEI6[Month %in% 4:10], na.rm = TRUE),
    SPEI6_min = safe_min(SPEI6),
    .groups = "drop"
  )

if (all(is.na(prism$SPEI6))) {
  annual_prism_full$SPEI6_gs <- NA_real_
  annual_prism_full$SPEI6_min <- NA_real_
}

incomplete_prism_rows <- annual_prism_full %>% filter(PRISM_month_n < 12)
prism_qc_rows[[length(prism_qc_rows) + 1L]] <- make_qc_row(
  "Incomplete PRISM years excluded from annual regressions",
  nrow(incomplete_prism_rows),
  incomplete_prism_rows$Year,
  if (nrow(incomplete_prism_rows) == 0) {
    "All annual PRISM years had 12 months available"
  } else {
    paste(paste0(incomplete_prism_rows$Year, " (", incomplete_prism_rows$PRISM_month_n, " months)"), collapse = "; ")
  },
  if (nrow(incomplete_prism_rows) == 0) "No issue found" else "Dropped from annual PRISM summaries"
)

annual_prism <- annual_prism_full %>% filter(PRISM_month_n == 12)
prism_incomplete_years <- incomplete_prism_rows$Year

if (!is.null(pdsi_ann)) {
  annual_prism <- left_join(annual_prism, pdsi_ann, by = "Year")
} else if (!("PDSI" %in% names(annual_prism))) {
  annual_prism$PDSI <- NA_real_
}

prism_qc <- bind_rows(prism_qc_rows)
write.csv(prism_qc, file.path(out_dir, "climate_qc_flags.csv"), row.names = FALSE)

enso_ann <- NULL
enso_source_file <- "not found"
if (file.exists(enso_file)) {
  enso_raw <- read.table(enso_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  enso_names <- canon_name(names(enso_raw))
  enso_year_idx <- pick_col(enso_names, c("yr", "year"))
  enso_mon_idx <- pick_col(enso_names, c("mon", "month"))
  enso_nino34_idx <- pick_col(enso_names, c("nino34"))
  enso_anom_idx <- NA_integer_
  if (!is.na(enso_nino34_idx) && enso_nino34_idx < ncol(enso_raw)) {
    next_name <- enso_names[enso_nino34_idx + 1]
    if (grepl("^anom", next_name)) enso_anom_idx <- enso_nino34_idx + 1L
  }
  if (is.na(enso_anom_idx)) {
    anom_candidates <- which(grepl("^anom", enso_names))
    if (length(anom_candidates) > 0) enso_anom_idx <- tail(anom_candidates, 1)
  }

  if (!is.na(enso_year_idx) && !is.na(enso_mon_idx) && !is.na(enso_anom_idx)) {
    enso_ann <- data.frame(
      Year = as.integer(as.numeric(enso_raw[[enso_year_idx]])),
      Month = as.integer(as.numeric(enso_raw[[enso_mon_idx]])),
      ENSO_NINO34 = as.numeric(enso_raw[[enso_anom_idx]]),
      stringsAsFactors = FALSE
    )
    enso_ann <- enso_ann[is.finite(enso_ann$Year) & is.finite(enso_ann$Month), , drop = FALSE]
    if (nrow(enso_ann) > 0) {
      enso_ann <- enso_ann %>%
        group_by(Year) %>%
        summarise(ENSO_NINO34 = mean(ENSO_NINO34, na.rm = TRUE), .groups = "drop")
      enso_source_file <- basename(enso_file)
      write.csv(enso_ann, file.path(out_dir, "enso_annual_nino34.csv"), row.names = FALSE)
    } else {
      enso_ann <- NULL
    }
  }
} else if (file.exists(enso_annual_file)) {
  enso_raw <- read.csv(enso_annual_file, check.names = FALSE, stringsAsFactors = FALSE)
  enso_names <- canon_name(names(enso_raw))
  enso_year_idx <- pick_col(enso_names, c("year", "yr"))
  enso_value_idx <- pick_col(enso_names, c("ensonino34", "nino34", "anom"))

  if (!is.na(enso_year_idx) && !is.na(enso_value_idx)) {
    enso_ann <- data.frame(
      Year = as.integer(as.numeric(enso_raw[[enso_year_idx]])),
      ENSO_NINO34 = as.numeric(enso_raw[[enso_value_idx]]),
      stringsAsFactors = FALSE
    )
    enso_ann <- enso_ann[is.finite(enso_ann$Year), , drop = FALSE]
    if (nrow(enso_ann) > 0) {
      enso_ann <- enso_ann %>%
        group_by(Year) %>%
        summarise(ENSO_NINO34 = if (all(is.na(ENSO_NINO34))) NA_real_ else mean(ENSO_NINO34, na.rm = TRUE), .groups = "drop")
      enso_source_file <- basename(enso_annual_file)
      write.csv(enso_ann, file.path(out_dir, "enso_annual_nino34.csv"), row.names = FALSE)
    } else {
      enso_ann <- NULL
    }
  }
}

# ==========================
# C) TREE-RING CHRONOLOGIES
# ==========================
grow_long <- clean_rwl_for_detrend(safe_read_rwl(rwl_long_file))
long_det <- detrend(grow_long, method = "Spline", f = 0.5, pos.slope = FALSE)
long_crn <- chron(long_det, prefix = "WWN", biweight = TRUE, prewhiten = TRUE)

long_df <- data.frame(
  Year = as.integer(rownames(long_crn)),
  RWI_std = as.numeric(long_crn[, "std"]),
  RWI_res = as.numeric(long_crn[, "res"])
)

long_merge <- inner_join(long_df, annual_prism, by = "Year")
if (!is.null(enso_ann)) {
  long_merge <- left_join(long_merge, enso_ann, by = "Year")
  short_df <- left_join(short_df, enso_ann %>% rename(year = Year), by = "year")
}

long_comp_candidates <- intersect(c("PPT_mm", "PDSI", "SPEI6_gs", "TMEAN_C"), names(long_merge))
long_comp_parts <- lapply(long_comp_candidates, function(v) {
  zv <- zfun(long_merge[[v]])
  if (v == "TMEAN_C") zv <- -zv
  zv
})
if (length(long_comp_parts) > 0) {
  long_comp_mat <- do.call(cbind, long_comp_parts)
  long_merge$ClimateComposite <- rowMeans(long_comp_mat, na.rm = TRUE)
  long_merge$ClimateComposite[is.nan(long_merge$ClimateComposite)] <- NA_real_
} else {
  long_merge$ClimateComposite <- NA_real_
}

write.csv(long_merge, file.path(out_dir, "long_prism_rwi_annual_merge.csv"), row.names = FALSE)

# ==========================
# D) PRISM LONG-PERIOD CORRELATIONS (PEARSON + SPEARMAN)
# ==========================

long_predictors <- c("PPT_mm", "TMEAN_C", "TMIN_C", "TMAX_C", "VPDMIN", "VPDMAX", "SPEI6_gs", "SPEI6_min", "PDSI")
long_predictors <- long_predictors[sapply(long_predictors, function(v) v %in% names(long_merge) && !all(is.na(long_merge[[v]])))]

long_lag_pear <- make_lag_corr(long_merge, "Year", "RWI_res", long_predictors, max_lag = max_lag, method = "pearson")
long_lag_spear <- make_lag_corr(long_merge, "Year", "RWI_res", long_predictors, max_lag = max_lag, method = "spearman")

write.csv(long_lag_pear, file.path(out_dir, "long_prism_lag_corr_pearson.csv"), row.names = FALSE)
write.csv(long_lag_spear, file.path(out_dir, "long_prism_lag_corr_spearman.csv"), row.names = FALSE)

plot_lag_heatmap(
  long_lag_pear,
  ww_title("PRISM Climate vs RWI", "Long-Period Lag Correlations (Pearson)"),
  "Tile fill = correlation; stars show significance (. <0.10, * <0.05, ** <0.01, *** <0.001)",
  file.path(plot_dir, "01_long_prism_lag_heatmap_pearson.png")
)

plot_lag_heatmap(
  long_lag_spear,
  ww_title("PRISM Climate vs RWI", "Long-Period Lag Correlations (Spearman)"),
  "Tile fill = correlation; stars show significance (. <0.10, * <0.05, ** <0.01, *** <0.001)",
  file.path(plot_dir, "02_long_prism_lag_heatmap_spearman.png")
)

build_best_lag_regression_panels(
  df = long_merge,
  year_col = "Year",
  response = "RWI_res",
  predictors = long_predictors,
  lag_table = long_lag_pear,
  title_txt = ww_title("RWI vs PRISM Climate", "Best-Lag Regression Panels"),
  subtitle_txt = "Panels use strongest absolute Pearson lag per predictor; arrows show time direction",
  file_out = file.path(plot_dir, "02b_long_prism_best_lag_regression_panels.png"),
  stats_out = file.path(out_dir, "long_prism_best_lag_regression_stats.csv"),
  ncol = 3
)
if (!is.null(enso_ann) && "ENSO_NINO34" %in% names(long_merge) && !all(is.na(long_merge$ENSO_NINO34))) {
  enso_lag_pear <- make_lag_corr(long_merge, "Year", "RWI_res", c("ENSO_NINO34"), max_lag = max_lag, method = "pearson")
  enso_lag_spear <- make_lag_corr(long_merge, "Year", "RWI_res", c("ENSO_NINO34"), max_lag = max_lag, method = "spearman")
  write.csv(enso_lag_pear, file.path(out_dir, "rwi_vs_enso_lag_corr_pearson.csv"), row.names = FALSE)
  write.csv(enso_lag_spear, file.path(out_dir, "rwi_vs_enso_lag_corr_spearman.csv"), row.names = FALSE)

  build_best_lag_regression_panels(
    df = long_merge,
    year_col = "Year",
    response = "RWI_res",
    predictors = c("ENSO_NINO34"),
    lag_table = enso_lag_pear,
    title_txt = ww_title("RWI vs ENSO (Nino3.4)", "Best-Lag Regression Panel"),
    subtitle_txt = "Annual mean Nino3.4 anomaly from NOAA CPC; panel uses strongest absolute Pearson lag",
    file_out = file.path(plot_dir, "03c_rwi_enso_best_lag_regression_panel.png"),
    stats_out = file.path(out_dir, "rwi_enso_best_lag_regression_stats.csv"),
    ncol = 1
  )

  enso_wave_df <- long_merge %>%
    dplyr::select(Year, RWI_res, ENSO_NINO34) %>%
    filter(complete.cases(.))

  old_enso_wavelet_file <- file.path(plot_dir, "03d_rwi_enso_wavelet_coherence.png")
  if (file.exists(old_enso_wavelet_file)) {
    file.remove(old_enso_wavelet_file)
  }

  save_standard_wavelet_plot(
    years = enso_wave_df$Year,
    values = enso_wave_df$RWI_res,
    file_out = file.path(plot_dir, "03d_rwi_enso_wavelet_transform.png"),
    main_txt = ww_title("RWI Residual Chronology", "Wavelet Transform"),
    subtitle_txt = "Morlet wavelet on ENSO-overlap years; 2-7 year bands indicate ENSO-scale periodicity.",
    crn_label = "RWI residual"
  )
}

# Long-period composite overlay
long_series <- c("RWI_res", "PPT_mm", "TMEAN_C", "SPEI6_gs", "VPDMAX")
long_series <- long_series[sapply(long_series, function(v) v %in% names(long_merge) && !all(is.na(long_merge[[v]])))]

long_z_rows <- lapply(long_series, function(v) {
  data.frame(Year = long_merge$Year, z = zfun(long_merge[[v]]), series = v, stringsAsFactors = FALSE)
})
long_z <- bind_rows(long_z_rows)

p_long_comp <- ggplot(long_z, aes(Year, z, color = series)) +
  geom_line(linewidth = 1.0, alpha = 0.9) +
  scale_color_manual(values = c('RWI_res' = COL_RWI, 'PPT_mm' = COL_PPT, 'TMEAN_C' = COL_TEMP, 'SPEI6_gs' = COL_SPEI, 'VPDMAX' = COL_VPD)) +
  labs(
    title = ww_title("RWI + PRISM Climate", "Long-Period Composite z-score Overlay"),
    subtitle = "Legend explains each line. SPEI uses monthly water balance (mm, deg C) with 6-month scale.",
    x = "Year",
    y = "Standardized value (z)",
    color = "Line legend"
  ) +
  theme_pub(base_size = 15)

ggsave(file.path(plot_dir, "03_long_prism_composite_zscore_overlay.png"), p_long_comp, width = 12.8, height = 6.4, dpi = save_dpi)

if ("ClimateComposite" %in% names(long_merge) && !all(is.na(long_merge$ClimateComposite))) {
  long_comp_ts <- bind_rows(
    data.frame(Year = long_merge$Year, z = zfun(long_merge$RWI_res), series = "RWI_res"),
    data.frame(Year = long_merge$Year, z = zfun(long_merge$ClimateComposite), series = "ClimateComposite")
  )

  p_long_comp_overlay <- ggplot(long_comp_ts, aes(Year, z, color = series)) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = c("RWI_res" = COL_RWI, "ClimateComposite" = COL_DEER)) +
    labs(
      title = ww_title("RWI + Climate Composite", "Long-Period Time Series"),
      subtitle = "ClimateComposite = mean z-score across available long-term climate variables (TMEAN_C sign reversed)",
      x = "Year",
      y = "Standardized value (z)",
      color = "Line legend"
    ) +
    theme_pub(base_size = 15)

  ggsave(file.path(plot_dir, "03f_long_prism_climate_composite_timeseries.png"), p_long_comp_overlay, width = 11.8, height = 5.8, dpi = save_dpi)
}


# Flood and drought marker years from climate signals
marker_base <- long_merge %>%
  transmute(
    Year = Year,
    RWI_std = RWI_res,
    PPT_mm = PPT_mm,
    SPEI6_gs = if ("SPEI6_gs" %in% names(long_merge)) SPEI6_gs else NA_real_,
    VPDMIN = if ("VPDMIN" %in% names(long_merge)) VPDMIN else NA_real_,
    VPDMAX = if ("VPDMAX" %in% names(long_merge)) VPDMAX else NA_real_
  )

spei_z <- if (all(is.na(marker_base$SPEI6_gs))) rep(NA_real_, nrow(marker_base)) else zfun(marker_base$SPEI6_gs)
flood_mat <- cbind(zfun(marker_base$PPT_mm), spei_z, zfun(marker_base$VPDMIN))
drought_mat <- cbind(-zfun(marker_base$PPT_mm), -spei_z, zfun(marker_base$VPDMAX))

marker_base$flood_score <- rowMeans(flood_mat, na.rm = TRUE)
marker_base$drought_score <- rowMeans(drought_mat, na.rm = TRUE)
marker_base$flood_score[is.nan(marker_base$flood_score)] <- NA_real_
marker_base$drought_score[is.nan(marker_base$drought_score)] <- NA_real_

flood_thr <- as.numeric(quantile(marker_base$flood_score, probs = hydro_marker_quantile, na.rm = TRUE))
drought_thr <- as.numeric(quantile(marker_base$drought_score, probs = hydro_marker_quantile, na.rm = TRUE))

marker_base$is_flood <- !is.na(marker_base$flood_score) & marker_base$flood_score >= flood_thr
marker_base$is_drought <- !is.na(marker_base$drought_score) & marker_base$drought_score >= drought_thr

marker_base$event_type <- dplyr::case_when(
  marker_base$is_flood & marker_base$is_drought ~ "Flood + Drought marker",
  marker_base$is_flood ~ "Flood marker",
  marker_base$is_drought ~ "Drought marker",
  TRUE ~ "No marker"
)

marker_events <- marker_base %>%
  filter(event_type != "No marker") %>%
  mutate(event_source = "Climate-derived")

marker_band_years <- sort(unique(marker_events$Year))

write.csv(marker_base, file.path(out_dir, "flood_drought_marker_scores.csv"), row.names = FALSE)

# Optional documented event overlay file
# Expected columns: Year, EventType, Source (optional)
documented_events <- data.frame(stringsAsFactors = FALSE)
if (file.exists(documented_hydro_events_file)) {
  doc <- read.csv(documented_hydro_events_file, stringsAsFactors = FALSE)
  if (all(c("Year", "EventType") %in% names(doc))) {
    doc$Year <- as.integer(doc$Year)
    evt <- tolower(trimws(as.character(doc$EventType)))
    doc$EventTypeLabel <- ifelse(grepl("flood|wet", evt), "Documented flood",
                                 ifelse(grepl("drought|dry", evt), "Documented drought", "Documented event"))
    doc$SourceLabel <- if ("Source" %in% names(doc)) as.character(doc$Source) else "Documented"
    documented_events <- doc[!is.na(doc$Year), c("Year", "EventTypeLabel", "SourceLabel")]
  }
} else {
  write.csv(
    data.frame(Year = integer(0), EventType = character(0), Source = character(0), stringsAsFactors = FALSE),
    file.path(out_dir, "documented_flood_drought_years_template.csv"),
    row.names = FALSE
  )
}

documented_marker_rows <- data.frame(stringsAsFactors = FALSE)
if (nrow(documented_events) > 0) {
  documented_marker_rows <- documented_events %>%
    transmute(
      Year = Year,
      RWI_std = NA_real_,
      PPT_mm = NA_real_,
      SPEI6_gs = NA_real_,
      VPDMIN = NA_real_,
      VPDMAX = NA_real_,
      flood_score = NA_real_,
      drought_score = NA_real_,
      is_flood = grepl("flood", tolower(EventTypeLabel)),
      is_drought = grepl("drought", tolower(EventTypeLabel)),
      event_type = EventTypeLabel,
      event_source = SourceLabel
    )
  marker_band_years <- sort(unique(c(marker_band_years, documented_events$Year)))
}

marker_events_export <- bind_rows(marker_events, documented_marker_rows) %>%
  arrange(Year, event_source)

write.csv(marker_events_export, file.path(out_dir, "flood_drought_marker_years.csv"), row.names = FALSE)

max_rwi <- suppressWarnings(max(marker_base$RWI_std, na.rm = TRUE))
if (!is.finite(max_rwi)) max_rwi <- 1

marker_points <- marker_events %>%
  transmute(
    Year = Year,
    y = max_rwi + 0.15,
    event_type = event_type,
    source = event_source
  )

if (nrow(documented_events) > 0) {
  marker_points <- bind_rows(
    marker_points,
    documented_events %>% transmute(Year = Year, y = max_rwi + 0.35, event_type = EventTypeLabel, source = SourceLabel)
  )
}

line_events <- marker_events_export %>%
  distinct(Year, event_type, event_source, .keep_all = TRUE)
flood_lines <- line_events %>% filter(grepl("Flood", event_type))
drought_lines <- line_events %>% filter(grepl("Drought", event_type))

p_marker_years <- ggplot(marker_base, aes(Year, RWI_std)) +
  geom_line(color = COL_RWI, linewidth = 1.05) +
  geom_vline(data = flood_lines, aes(xintercept = Year, linetype = "Flood marker year"), color = COL_PPT, alpha = 0.35, linewidth = 0.55, show.legend = TRUE) +
  geom_vline(data = drought_lines, aes(xintercept = Year, linetype = "Drought marker year"), color = COL_PDSI, alpha = 0.35, linewidth = 0.55, show.legend = TRUE) +
  geom_point(data = marker_points, aes(Year, y, shape = event_type, color = source), size = 2.6, inherit.aes = FALSE) +
  scale_linetype_manual(values = c("Flood marker year" = "dotted", "Drought marker year" = "dashed"), name = "Line legend") +
  labs(
    title = ww_title("Flood and Drought Marker Years", "RWI Timeline with Marker Years"),
    subtitle = paste0("Marker years are climate-derived (top ", round((1 - hydro_marker_quantile) * 100), "% signal years). Add documented events via documented_flood_drought_years.csv"),
    x = "Year",
    y = "RWI (std chronology)",
    color = "Marker source",
    shape = "Marker type"
  ) +
  expand_limits(y = max_rwi + 0.55) +
  theme_pub(base_size = 15)

ggsave(file.path(plot_dir, "03a_flood_drought_marker_years_rwi.png"), p_marker_years, width = 13.0, height = 6.4, dpi = save_dpi)

score_long <- bind_rows(
  marker_base %>% transmute(Year = Year, score = flood_score, score_type = "Flood score"),
  marker_base %>% transmute(Year = Year, score = drought_score, score_type = "Drought score")
)
threshold_df <- data.frame(score_type = c("Flood score", "Drought score"), threshold = c(flood_thr, drought_thr), stringsAsFactors = FALSE)

marker_score_points <- bind_rows(
  marker_base %>% filter(is_flood) %>% transmute(Year = Year, score = flood_score, score_type = "Flood score"),
  marker_base %>% filter(is_drought) %>% transmute(Year = Year, score = drought_score, score_type = "Drought score")
)

p_marker_scores <- ggplot(score_long, aes(Year, score, color = score_type)) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  geom_hline(data = threshold_df, aes(yintercept = threshold, color = score_type), linetype = "dashed", linewidth = 0.7, show.legend = FALSE) +
  geom_point(data = marker_score_points, aes(Year, score, color = score_type), size = 2.0, show.legend = FALSE) +
  labs(
    title = ww_title("Flood and Drought Marker Years", "Climate Marker Score Time Series"),
    subtitle = "Dashed lines are marker thresholds; points are selected marker years",
    x = "Year",
    y = "Marker score (z composite)",
    color = "Score series"
  ) +
  theme_pub(base_size = 15)

ggsave(file.path(plot_dir, "03b_flood_drought_marker_scores.png"), p_marker_scores, width = 12.2, height = 5.8, dpi = save_dpi)

long_ts_vars <- c("PPT_mm", "TMIN_C", "TMEAN_C", "TMAX_C", "VPDMIN", "VPDMAX", "SPEI6_gs", "SPEI6_min", "PDSI", "ENSO_NINO34")
long_ts_vars <- long_ts_vars[long_ts_vars %in% names(long_merge)]

plot_faceted_timeseries_with_markers(
  df = long_merge,
  year_col = "Year",
  variables = long_ts_vars,
  title_txt = ww_title("Climatology Climate Variables", "Faceted Time Series with Marker Years"),
  subtitle_txt = "Each panel shows one climate variable. Light gray bands mark climate-derived flood and drought marker years.",
  file_out = file.path(plot_dir, "03e_climatology_climate_timeseries_marker_bands.png"),
  marker_years = marker_band_years,
  ncol = 2
)

build_all_lag_regression_panels(
  df = long_merge,
  year_col = "Year",
  response = "RWI_res",
  predictors = long_predictors,
  max_lag = max_lag,
  title_txt = ww_title("RWI vs Climatology Climate", "All Lag Regression Panels"),
  subtitle_txt = "Every panel shows one climate-variable lag (0 to 2 years) against residual RWI with Pearson fit and sample size.",
  file_out = file.path(plot_dir, "02c_climatology_all_lag_regression_panels.png"),
  stats_out = file.path(out_dir, "long_climatology_all_lag_regression_stats.csv"),
  ncol = 4
)

build_all_lag_regression_panels(
  df = long_merge,
  year_col = "Year",
  response = "RWI_res",
  predictors = long_predictors,
  max_lag = max_lag,
  title_txt = ww_title("RWI vs Climatology Climate", "All Lag Scatter Panels (Spearman)"),
  subtitle_txt = "Every panel shows one climate-variable lag (0 to 2 years) against residual RWI with the same scatter style and Spearman rank statistics.",
  file_out = file.path(plot_dir, "02d_climatology_all_lag_spearman_panels.png"),
  stats_out = file.path(out_dir, "long_climatology_all_lag_regression_stats_spearman.csv"),
  ncol = 4,
  stat_method = "spearman"
)
# ==========================
# E) MODIS-ERA SHORT CORRELATIONS + COMPOSITES
# ==========================

short_base_vars <- unique(c(year_var, response_var, short_predictors, "ENSO_NINO34"))
short_base_vars <- short_base_vars[short_base_vars %in% names(short_df)]
short_core <- short_df %>% dplyr::select(all_of(short_base_vars))

comp_candidates <- intersect(c("EVT", "ETA", "PPT", "PDSI", "SPEI", "TMEAN"), names(short_core))
comp_parts <- lapply(comp_candidates, function(v) {
  zv <- zfun(short_core[[v]])
  if (v == "TMEAN") zv <- -zv
  zv
})
if (length(comp_parts) > 0) {
  comp_mat <- do.call(cbind, comp_parts)
  short_core$ClimateComposite <- rowMeans(comp_mat, na.rm = TRUE)
  short_core$ClimateComposite[is.nan(short_core$ClimateComposite)] <- NA_real_
} else {
  short_core$ClimateComposite <- NA_real_
}

short_profile_predictors <- intersect(
  c("ClimateComposite", short_predictors, "ENSO_NINO34"),
  names(short_core)
)
short_profile_predictors <- short_profile_predictors[
  sapply(short_profile_predictors, function(v) !all(is.na(short_core[[v]])))
]
short_core_vars <- c(response_var, short_profile_predictors)
short_core <- short_core %>% dplyr::select(all_of(c(year_var, short_core_vars)))

short_mat_pear <- make_matrix_corr(short_core, short_core_vars, method = "pearson")
short_mat_spear <- make_matrix_corr(short_core, short_core_vars, method = "spearman")

write.csv(short_mat_pear, file.path(out_dir, "short_modis_corr_matrix_pearson.csv"), row.names = FALSE)
write.csv(short_mat_spear, file.path(out_dir, "short_modis_corr_matrix_spearman.csv"), row.names = FALSE)

plot_matrix_heatmap(
  short_mat_pear,
  ww_title("MODIS-Era Climate + Growth", "Correlation Matrix (Pearson)"),
  "Includes all tested short-window inputs, with ClimateComposite and ENSO_NINO34 shown when available; stars show significance.",
  file.path(plot_dir, "04_short_modis_corr_heatmap_pearson.png")
)

plot_matrix_heatmap(
  short_mat_spear,
  ww_title("MODIS-Era Climate + Growth", "Correlation Matrix (Spearman)"),
  "Includes all tested short-window inputs, with ClimateComposite and ENSO_NINO34 shown when available; stars show significance.",
  file.path(plot_dir, "05_short_modis_corr_heatmap_spearman.png")
)

short_lag_pear <- make_lag_corr(short_core, year_var, response_var, short_profile_predictors, max_lag = max_lag, method = "pearson")
short_lag_spear <- make_lag_corr(short_core, year_var, response_var, short_profile_predictors, max_lag = max_lag, method = "spearman")

write.csv(short_lag_pear, file.path(out_dir, "short_modis_lag_corr_pearson.csv"), row.names = FALSE)
write.csv(short_lag_spear, file.path(out_dir, "short_modis_lag_corr_spearman.csv"), row.names = FALSE)

plot_lag_heatmap(
  short_lag_pear,
  ww_title("MODIS-Era Climate vs RWI", "Lag Correlations (Pearson)"),
  "Tile fill = correlation for every tested short-window input, including ClimateComposite and ENSO_NINO34 when available.",
  file.path(plot_dir, "06_short_modis_lag_heatmap_pearson.png")
)

plot_lag_heatmap(
  short_lag_spear,
  ww_title("MODIS-Era Climate vs RWI", "Lag Correlations (Spearman)"),
  "Tile fill = rank correlation for every tested short-window input, including ClimateComposite and ENSO_NINO34 when available.",
  file.path(plot_dir, "07_short_modis_lag_heatmap_spearman.png")
)

# RWI vs each short predictor (best-lag regression panels with full stats)
build_best_lag_regression_panels(
  df = short_core,
  year_col = year_var,
  response = response_var,
  predictors = short_profile_predictors,
  lag_table = short_lag_pear,
  title_txt = ww_title("RWI vs MODIS-Era Predictors", "Best-Lag Regression Panels"),
  subtitle_txt = "Panels use the strongest absolute Pearson lag for every tested short-window input; arrows show time direction.",
  file_out = file.path(plot_dir, "08_short_modis_rwi_predictor_panels.png"),
  stats_out = file.path(out_dir, "short_modis_best_lag_regression_stats.csv"),
  ncol = 3
)

plot_faceted_timeseries_with_markers(
  df = short_core,
  year_col = year_var,
  variables = short_profile_predictors,
  title_txt = ww_title("MODIS-Era Climate and Productivity Variables", "Faceted Time Series with Marker Years"),
  subtitle_txt = "Each panel shows one tested short-window input, including ClimateComposite and ENSO_NINO34 when available. Gray bands mark marker years.",
  file_out = file.path(plot_dir, "09b_modis_climate_timeseries_marker_bands.png"),
  marker_years = marker_band_years,
  ncol = 2
)

build_all_lag_regression_panels(
  df = short_core,
  year_col = year_var,
  response = response_var,
  predictors = short_profile_predictors,
  max_lag = max_lag,
  title_txt = ww_title("RWI vs MODIS-Era Predictors", "All Lag Regression Panels"),
  subtitle_txt = "Every panel shows one tested short-window predictor-lag combination (0 to 2 years) against RWI with Pearson fit and sample size.",
  file_out = file.path(plot_dir, "08b_short_modis_all_lag_regression_panels.png"),
  stats_out = file.path(out_dir, "short_modis_all_lag_regression_stats.csv"),
  ncol = 4
)

build_all_lag_regression_panels(
  df = short_core,
  year_col = year_var,
  response = response_var,
  predictors = short_profile_predictors,
  max_lag = max_lag,
  title_txt = ww_title("RWI vs MODIS-Era Predictors", "All Lag Scatter Panels (Spearman)"),
  subtitle_txt = "Every panel shows one tested short-window predictor-lag combination (0 to 2 years) against RWI with the same scatter style and Spearman rank statistics.",
  file_out = file.path(plot_dir, "08c_short_modis_all_lag_spearman_panels.png"),
  stats_out = file.path(out_dir, "short_modis_all_lag_regression_stats_spearman.csv"),
  ncol = 4,
  stat_method = "spearman"
)

short_comp_ts <- bind_rows(
  data.frame(year = short_core[[year_var]], z = zfun(short_core[[response_var]]), series = "RWI"),
  data.frame(year = short_core[[year_var]], z = zfun(short_core$ClimateComposite), series = "ClimateComposite")
)

p_short_comp_ts <- ggplot(short_comp_ts, aes(year, z, color = series)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c('RWI' = COL_RWI, 'ClimateComposite' = COL_DEER)) +
  labs(
    title = ww_title("RWI + Climate Composite", "MODIS-Era Time Series"),
    subtitle = "ClimateComposite = mean z-score across available climate variables (TMEAN sign reversed)",
    x = "Year",
    y = "Standardized value (z)",
    color = "Line legend"
  ) +
  theme_pub(base_size = 15)

ggsave(file.path(plot_dir, "09_short_modis_composite_timeseries.png"), p_short_comp_ts, width = 11.5, height = 5.8, dpi = save_dpi)

comp_cor_p <- safe_cor(short_core[[response_var]], short_core$ClimateComposite, method = "pearson")
comp_cor_s <- safe_cor(short_core[[response_var]], short_core$ClimateComposite, method = "spearman")

p_short_comp_scatter <- ggplot(short_core, aes(ClimateComposite, .data[[response_var]])) +
  geom_point(color = COL_RWI, alpha = 0.9, size = 2.6) +
  geom_smooth(method = "lm", se = TRUE, color = COL_FIT, linewidth = 1.0) +
  annotate(
    "text", x = Inf, y = Inf,
    label = paste0(
      "Pearson r = ", format(round(comp_cor_p$r, 3), nsmall = 3), " (", fmt_p(comp_cor_p$p), ")",
      "\nSpearman rho = ", format(round(comp_cor_s$r, 3), nsmall = 3), " (", fmt_p(comp_cor_s$p), ")"
    ),
    hjust = 1.05, vjust = 1.1, size = 3.5
  ) +
  labs(
    title = ww_title("RWI vs Climate Composite", "MODIS-Era Composite Regression"),
    subtitle = "Composite combines all available short-period climate inputs",
    x = "Climate composite (z-average)",
    y = "RWI"
  ) +
  theme_pub(base_size = 15)

ggsave(file.path(plot_dir, "10_short_modis_composite_vs_rwi.png"), p_short_comp_scatter, width = 9.8, height = 6.2, dpi = save_dpi)

write.csv(short_core, file.path(out_dir, "short_modis_with_composite.csv"), row.names = FALSE)

# ==========================
# F) DEER BROWSING + CLIMATE + GROWTH
# ==========================

deer_raw <- read_excel(deer_file, sheet = deer_sheet)
clean_deer_names <- canon_name(names(deer_raw))

deer_year_idx <- pick_col(clean_deer_names, c("year", "yr"))
deer_rwi_idx <- pick_col(clean_deer_names, c("rwi"))
deer_pop_idx <- pick_col(clean_deer_names, c("deerpopulation", "deerpop", "deer", "deercount", "deerabundance"))

if (is.na(deer_year_idx)) stop("Could not find year column in deer workbook.")
if (is.na(deer_pop_idx)) stop("Could not find deer population column in deer workbook.")

deer_df <- data.frame(
  year = as.integer(as.numeric(deer_raw[[deer_year_idx]])),
  deer_population = as.numeric(deer_raw[[deer_pop_idx]]),
  stringsAsFactors = FALSE
)

if (!is.na(deer_rwi_idx)) {
  deer_df$deer_sheet_rwi <- as.numeric(deer_raw[[deer_rwi_idx]])
}

extra_numeric_idx <- which(sapply(deer_raw, is.numeric))
extra_numeric_idx <- setdiff(extra_numeric_idx, c(deer_year_idx, deer_pop_idx, deer_rwi_idx))

if (length(extra_numeric_idx) > 0) {
  for (i in extra_numeric_idx) {
    deer_df[[clean_deer_names[i]]] <- as.numeric(deer_raw[[i]])
  }
}

prism_by_year <- annual_prism %>% rename(year = Year)

deer_merge <- deer_df %>%
  left_join(short_core, by = c("year" = year_var)) %>%
  left_join(prism_by_year, by = "year")

if (all(c("PDSI.x", "PDSI.y") %in% names(deer_merge))) {
  deer_merge$PDSI <- dplyr::coalesce(deer_merge$PDSI.y, deer_merge$PDSI.x)
  deer_merge$PDSI.x <- NULL
  deer_merge$PDSI.y <- NULL
} else if ("PDSI.x" %in% names(deer_merge)) {
  names(deer_merge)[names(deer_merge) == "PDSI.x"] <- "PDSI"
} else if ("PDSI.y" %in% names(deer_merge)) {
  names(deer_merge)[names(deer_merge) == "PDSI.y"] <- "PDSI"
}

if (!is.null(enso_ann)) {
  deer_merge <- left_join(deer_merge, enso_ann %>% rename(year = Year), by = "year")
}

write.csv(deer_merge, file.path(out_dir, "deer_climate_growth_merged.csv"), row.names = FALSE)

# Deer + RWI time series
if ("RWI" %in% names(deer_merge)) {
  deer_ts <- bind_rows(
    data.frame(year = deer_merge$year, z = zfun(deer_merge$deer_population), series = "Deer population"),
    data.frame(year = deer_merge$year, z = zfun(deer_merge$RWI), series = "RWI")
  )

  p_deer_ts <- ggplot(deer_ts, aes(year, z, color = series)) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = c('Deer population' = COL_DEER, 'RWI' = COL_RWI)) +
    labs(
      title = ww_title("Deer Browsing + Growth", "Time-Series Overlay"),
      subtitle = "Both series are standardized (z) for direct visual comparison",
      x = "Year",
      y = "Standardized value (z)",
      color = "Line legend"
    ) +
    theme_pub(base_size = 15)

  ggsave(file.path(plot_dir, "11_deer_rwi_timeseries_overlay.png"), p_deer_ts, width = 11.2, height = 5.6, dpi = save_dpi)

  deer_rwi_p <- safe_cor(deer_merge$RWI, deer_merge$deer_population, method = "pearson")
  deer_rwi_s <- safe_cor(deer_merge$RWI, deer_merge$deer_population, method = "spearman")

  p_deer_scatter <- ggplot(deer_merge, aes(deer_population, RWI)) +
    geom_point(color = COL_RWI, size = 2.6, alpha = 0.9) +
    geom_smooth(method = "lm", se = TRUE, color = COL_FIT, linewidth = 1.0) +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste0(
        "Pearson r = ", format(round(deer_rwi_p$r, 3), nsmall = 3), " (", fmt_p(deer_rwi_p$p), ")",
        "\nSpearman rho = ", format(round(deer_rwi_s$r, 3), nsmall = 3), " (", fmt_p(deer_rwi_s$p), ")"
      ),
      hjust = 1.05, vjust = 1.1, size = 3.5
    ) +
    labs(
      title = ww_title("RWI vs Deer Population", "Correlation Regression"),
      subtitle = "25-year overlap period",
      x = "Deer population",
      y = "RWI"
    ) +
    theme_pub(base_size = 15)

  ggsave(file.path(plot_dir, "12_deer_population_vs_rwi_scatter.png"), p_deer_scatter, width = 9.5, height = 6.2, dpi = save_dpi)
  if ("VPDMIN" %in% names(deer_merge)) {
    vpd_2005 <- deer_merge$VPDMIN[deer_merge$year == 2005]
    vpd_2005 <- vpd_2005[!is.na(vpd_2005)]
    vpd_2005 <- if (length(vpd_2005) > 0) vpd_2005[1] else NA_real_

    fac_long <- bind_rows(
      data.frame(
        year = deer_merge$year,
        VPDMIN = deer_merge$VPDMIN,
        response = deer_merge$deer_population,
        panel = "VPDMIN vs Deer Population",
        stringsAsFactors = FALSE
      ),
      data.frame(
        year = deer_merge$year,
        VPDMIN = deer_merge$VPDMIN,
        response = deer_merge$RWI,
        panel = "VPDMIN vs RWI",
        stringsAsFactors = FALSE
      )
    )
    fac_long <- fac_long[complete.cases(fac_long), ]
    fac_long$obs_group <- ifelse(fac_long$year == 2005, "2005 drought year", "Other years")

    panel_label <- function(dfp) {
      cp <- safe_cor(dfp$response, dfp$VPDMIN, method = "pearson")
      cs <- safe_cor(dfp$response, dfp$VPDMIN, method = "spearman")
      slope_p <- NA_real_
      keep <- complete.cases(dfp$response, dfp$VPDMIN)
      if (sum(keep) >= 4 && sd(dfp$VPDMIN[keep]) > 0) {
        slope_p <- tryCatch(summary(lm(response ~ VPDMIN, data = dfp[keep, ]))$coefficients[2, 4], error = function(e) NA_real_)
      }
      paste0(
        "Pearson r = ", ifelse(is.na(cp$r), "NA", format(round(cp$r, 3), nsmall = 3)),
        " (", fmt_p(cp$p), ")",
        "\nSpearman rho = ", ifelse(is.na(cs$r), "NA", format(round(cs$r, 3), nsmall = 3)),
        " (", fmt_p(cs$p), ")",
        "\nSlope ", fmt_p(slope_p),
        "\nn = ", cp$n
      )
    }

    panel_ann <- bind_rows(lapply(split(fac_long, fac_long$panel), function(dfp) {
      data.frame(panel = unique(dfp$panel), label = panel_label(dfp), stringsAsFactors = FALSE)
    }))

    p_vpd_faceted <- ggplot(fac_long, aes(VPDMIN, response)) +
      geom_point(aes(color = obs_group), alpha = 0.92, size = 2.3) +
      geom_smooth(method = "lm", se = TRUE, color = COL_FIT, linewidth = 0.8) +
      facet_wrap(~panel, scales = "free_y", ncol = 2) +
      geom_text(
        data = panel_ann,
        aes(x = Inf, y = Inf, label = label),
        inherit.aes = FALSE,
        hjust = 1.05,
        vjust = 1.1,
        size = 3.1
      ) +
      scale_color_manual(values = c("2005 drought year" = COL_EVENT, "Other years" = COL_DEER)) +
      labs(
        title = ww_title("VPDMIN Sensitivity", "Faceted Scatter (Deer vs Canopy Growth)"),
        subtitle = "2005 drought observation is highlighted in red; dashed line marks the 2005 VPDMIN value when available",
        x = "VPDMIN",
        y = "Response value",
        color = "Observation group"
      ) +
      theme_pub(base_size = 15)

    if (!is.na(vpd_2005)) {
      p_vpd_faceted <- p_vpd_faceted +
        geom_vline(xintercept = vpd_2005, linetype = "dashed", color = COL_EVENT, linewidth = 0.8)
    }

    ggsave(file.path(plot_dir, "12b_vpdmin_deer_rwi_faceted_scatter_2005.png"), p_vpd_faceted, width = 12.4, height = 6.4, dpi = save_dpi)
  }
}

# Deer vs climate drivers
possible_drivers <- unique(c(short_profile_predictors, "PDSI", "PPT_mm", "TMEAN_C", "SPEI6_gs", "VPDMAX", "VPDMIN"))
possible_drivers <- possible_drivers[possible_drivers %in% names(deer_merge)]

if (length(possible_drivers) > 0) {
  deer_corr_rows <- list()
  k <- 1
  for (drv in possible_drivers) {
    for (meth in c("pearson", "spearman")) {
      ct <- safe_cor(deer_merge$deer_population, deer_merge[[drv]], method = meth)
      deer_corr_rows[[k]] <- data.frame(driver = drv, method = meth, r = ct$r, r2 = ct$r2, p = ct$p, n = ct$n, stringsAsFactors = FALSE)
      k <- k + 1
    }
  }
  deer_driver_corr <- bind_rows(deer_corr_rows)
  write.csv(deer_driver_corr, file.path(out_dir, "deer_vs_climate_correlations.csv"), row.names = FALSE)

  p_deer_corr_bar <- ggplot(deer_driver_corr, aes(reorder(driver, r), r, fill = method)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    coord_flip() +
    labs(
      title = ww_title("Deer Population vs Climate Inputs", "Pearson + Spearman Correlations"),
      subtitle = "Bars show correlation coefficients across available short/long drivers",
      x = "Driver",
      y = "Correlation with deer population",
      fill = "Method"
    ) +
    theme_pub(base_size = 15)

  ggsave(file.path(plot_dir, "13_deer_vs_climate_correlation_bars.png"), p_deer_corr_bar, width = 10.5, height = 6.4, dpi = save_dpi)

  deer_scatter_long <- bind_rows(lapply(possible_drivers, function(drv) {
    data.frame(
      year = deer_merge$year,
      driver = drv,
      climate_value = deer_merge[[drv]],
      deer_population = deer_merge$deer_population,
      stringsAsFactors = FALSE
    )
  }))

  deer_scatter_ann <- bind_rows(lapply(possible_drivers, function(drv) {
    cp <- safe_cor(deer_merge$deer_population, deer_merge[[drv]], method = "pearson")
    cs <- safe_cor(deer_merge$deer_population, deer_merge[[drv]], method = "spearman")
    data.frame(
      driver = drv,
      label = paste0(
        "Pearson r = ", format(round(cp$r, 3), nsmall = 3), " (", fmt_p(cp$p), ")",
        "\nSpearman rho = ", format(round(cs$r, 3), nsmall = 3), " (", fmt_p(cs$p), ")",
        "\nn = ", cp$n
      ),
      stringsAsFactors = FALSE
    )
  }))

  p_deer_clim_panels <- ggplot(deer_scatter_long, aes(climate_value, deer_population)) +
    geom_point(color = COL_DEER, alpha = 0.9, size = 2.0) +
    geom_smooth(method = "lm", se = TRUE, color = COL_FIT, linewidth = 1.0) +
    facet_wrap(~driver, scales = "free_x", ncol = 3) +
    geom_text(
      data = deer_scatter_ann,
      aes(x = Inf, y = Inf, label = label),
      hjust = 1.05,
      vjust = 1.1,
      inherit.aes = FALSE,
      size = 3.0
    ) +
    labs(
      title = ww_title("Deer Population vs Climate Inputs", "Scatter + Regression Panels"),
      subtitle = "Each panel shows deer response to one climate variable (with Pearson and Spearman stats)",
      x = "Climate value",
      y = "Deer population"
    ) +
    theme_pub(base_size = 15)

  ggsave(file.path(plot_dir, "13a_deer_vs_climate_scatter_panels.png"), p_deer_clim_panels, width = 13.4, height = 8.0, dpi = save_dpi)

  top_driver_tbl <- deer_driver_corr %>%
    filter(method == "pearson", !is.na(r)) %>%
    arrange(desc(abs(r))) %>%
    head(4)
  top_drivers <- unique(top_driver_tbl$driver)

  if (length(top_drivers) > 0) {
    deer_ts_long <- bind_rows(
      data.frame(year = deer_merge$year, z = zfun(deer_merge$deer_population), series = "Deer population", stringsAsFactors = FALSE),
      bind_rows(lapply(top_drivers, function(drv) {
        data.frame(
          year = deer_merge$year,
          z = zfun(deer_merge[[drv]]),
          series = paste0("Climate: ", drv),
          stringsAsFactors = FALSE
        )
      }))
    ) %>% filter(!is.na(z))

    p_deer_top_ts <- ggplot(deer_ts_long, aes(year, z, color = series)) +
      geom_line(linewidth = 0.9, alpha = 0.9) +
      labs(
        title = ww_title("Deer Population and Climate", "Top-Driver Time-Series Overlay"),
        subtitle = "Shows deer population with the top climate drivers by absolute Pearson correlation",
        x = "Year",
        y = "Standardized value (z)",
        color = "Line legend"
      ) +
      theme_pub(base_size = 15)

    ggsave(file.path(plot_dir, "13b_deer_vs_top_climate_timeseries.png"), p_deer_top_ts, width = 11.8, height = 6.0, dpi = save_dpi)
  }

  if (nrow(top_driver_tbl) > 0) {
    p_driver_rank <- top_driver_tbl %>%
      mutate(driver = reorder(driver, abs(r))) %>%
      ggplot(aes(abs(r), driver, color = r)) +
      geom_segment(aes(x = 0, xend = abs(r), y = driver, yend = driver), linewidth = 1.2, alpha = 0.8) +
      geom_point(size = 3.3) +
      scale_color_gradient2(low = COL_HEAT_NEG, mid = 'white', high = COL_HEAT_POS, midpoint = 0) +
      labs(
        title = ww_title("Deer Population vs Climate Inputs", "Top Driver Ranking"),
        subtitle = "Absolute Pearson correlation strength with full statistics in exported tables",
        x = "|Pearson r|",
        y = "Driver",
        color = "Signed r"
      ) +
      theme_pub(base_size = 15)

    ggsave(file.path(plot_dir, "20_deer_top_driver_rankings.png"), p_driver_rank, width = 10.8, height = 6.0, dpi = save_dpi)

    best_driver <- top_driver_tbl$driver[1]
    best_df <- deer_merge %>%
      transmute(year, deer_population = deer_population, driver_value = .data[[best_driver]]) %>%
      filter(complete.cases(.))

    p_best_driver <- ggplot(best_df, aes(driver_value, deer_population)) +
      geom_point(color = COL_DEER, size = 2.8, alpha = 0.9) +
      geom_smooth(method = "lm", se = TRUE, color = COL_FIT, linewidth = 1.1) +
      labs(
        title = ww_title("Deer Population vs Climate Inputs", "Best Driver Scatter"),
        subtitle = paste0("Deer residual vs top Pearson driver (", best_driver, ")"),
        x = best_driver,
        y = "Deer population"
      ) +
      theme_pub(base_size = 15)

    ggsave(file.path(plot_dir, "21_deer_vs_top_driver_scatter.png"), p_best_driver, width = 10.4, height = 6.2, dpi = save_dpi)
  }

  deer_lag_clim_pear <- make_lag_corr(deer_merge, "year", "deer_population", possible_drivers, max_lag = max_lag, method = "pearson")
  deer_lag_clim_spear <- make_lag_corr(deer_merge, "year", "deer_population", possible_drivers, max_lag = max_lag, method = "spearman")
  deer_lag_clim_all <- bind_rows(deer_lag_clim_pear, deer_lag_clim_spear)
  deer_lag_clim_export <- deer_lag_clim_all %>%
    rename(predictor = variable) %>%
    mutate(response = "deer_population", .before = 1)
  write.csv(deer_lag_clim_export, file.path(out_dir, "deer_lag_vs_climate_correlations.csv"), row.names = FALSE)

  plot_lag_heatmap(
    deer_lag_clim_pear,
    ww_title("Deer Population vs Climate Inputs", "Lag Correlation Heatmap (Pearson)"),
    "Lag 1 means prior-year climate value",
    file.path(plot_dir, "13c_deer_vs_climate_lag_heatmap_pearson.png")
  )

  plot_lag_heatmap(
    deer_lag_clim_spear,
    ww_title("Deer Population vs Climate Inputs", "Lag Correlation Heatmap (Spearman)"),
    "Lag 1 means prior-year climate value",
    file.path(plot_dir, "13d_deer_vs_climate_lag_heatmap_spearman.png")
  )

  if ("RWI" %in% names(deer_merge)) {
    deer_lag_pear <- make_lag_corr(deer_merge, "year", "RWI", c("deer_population"), max_lag = max_lag, method = "pearson")
    deer_lag_spear <- make_lag_corr(deer_merge, "year", "RWI", c("deer_population"), max_lag = max_lag, method = "spearman")
    deer_lag_all <- bind_rows(deer_lag_pear, deer_lag_spear)
    deer_lag_export <- deer_lag_all %>%
      rename(predictor = variable) %>%
      mutate(response = "RWI", .before = 1)
    write.csv(deer_lag_export, file.path(out_dir, "deer_lag_vs_rwi_correlations.csv"), row.names = FALSE)

    deer_lag_plot <- deer_lag_all %>%
      mutate(lag_f = paste0("Lag ", lag), label = heatmap_sig_label(p))

    p_deer_lag <- ggplot(deer_lag_plot, aes(lag_f, method, fill = r)) +
      geom_tile(color = "white") +
      geom_text(aes(label = label), size = 4.0, fontface = "bold") +
      scale_fill_gradient2(low = COL_HEAT_NEG, mid = 'white', high = COL_HEAT_POS, midpoint = 0, na.value = 'gray90') +
      labs(
        title = ww_title("RWI vs Deer Population", "Lag Correlation Summary"),
        subtitle = "Lag 1 means prior-year deer population",
        x = "Lag",
        y = "Correlation method",
        fill = "Correlation"
      ) +
      theme_pub(base_size = 15)

    ggsave(file.path(plot_dir, "14_deer_lag_vs_rwi_heatmap.png"), p_deer_lag, width = 7.8, height = 4.8, dpi = save_dpi)
  }

  deer_short_predictors <- intersect(c("EVT", "PPT", "TMEAN", "PDSI", "SPEI", "NDVI"), names(deer_merge))
  deer_long_predictors <- intersect(c("PPT_mm", "TMIN_C", "TMEAN_C", "TMAX_C", "VPDMIN", "VPDMAX", "SPEI6_gs", "SPEI6_min", "PDSI"), names(deer_merge))

  if (length(deer_short_predictors) > 0) {
    build_all_lag_regression_panels(
      df = deer_merge,
      year_col = "year",
      response = "deer_population",
      predictors = deer_short_predictors,
      max_lag = max_lag,
      title_txt = ww_title("Deer Population vs MODIS-Era Predictors", "All Lag Regression Panels"),
      subtitle_txt = "Every panel shows one MODIS-era predictor and lag (0 to 2 years) against deer population with Pearson fit and sample size.",
      file_out = file.path(plot_dir, "13e_deer_modis_all_lag_regression_panels.png"),
      stats_out = file.path(out_dir, "deer_modis_all_lag_regression_stats.csv"),
      ncol = 4
    )
  }

  if (length(deer_long_predictors) > 0) {
    build_all_lag_regression_panels(
      df = deer_merge,
      year_col = "year",
      response = "deer_population",
      predictors = deer_long_predictors,
      max_lag = max_lag,
      title_txt = ww_title("Deer Population vs Climatology Climate", "All Lag Regression Panels"),
      subtitle_txt = "Every panel shows one climatology predictor and lag (0 to 2 years) against deer population with Pearson fit and sample size.",
      file_out = file.path(plot_dir, "13f_deer_climatology_all_lag_regression_panels.png"),
      stats_out = file.path(out_dir, "deer_climatology_all_lag_regression_stats.csv"),
      ncol = 4
    )
  }
}

# Species-level support if workbook includes additional numeric response columns
species_cols <- setdiff(names(deer_df), c("year", "deer_population", "deer_sheet_rwi"))
if (length(species_cols) > 0 && length(possible_drivers) > 0) {
  sp_rows <- list()
  k <- 1
  for (sp in species_cols) {
    for (drv in possible_drivers) {
      for (meth in c("pearson", "spearman")) {
        ct <- safe_cor(deer_merge[[sp]], deer_merge[[drv]], method = meth)
        sp_rows[[k]] <- data.frame(species = sp, driver = drv, method = meth, r = ct$r, r2 = ct$r2, p = ct$p, n = ct$n, stringsAsFactors = FALSE)
        k <- k + 1
      }
    }
  }
  species_corr <- bind_rows(sp_rows)
  write.csv(species_corr, file.path(out_dir, "species_driver_correlations.csv"), row.names = FALSE)

  for (meth in c("pearson", "spearman")) {
    d <- species_corr %>% filter(method == meth) %>%
      mutate(
        species = factor(species, levels = rev(unique(species))),
        driver = factor(driver, levels = unique(driver)),
        label = heatmap_sig_label(p)
      )

    p_sp <- ggplot(d, aes(driver, species, fill = r)) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_text(aes(label = label), size = 3.2, fontface = "bold") +
      scale_fill_gradient2(low = COL_HEAT_NEG, mid = 'white', high = COL_HEAT_POS, midpoint = 0, na.value = 'gray90') +
      labs(
        title = ww_title("Species Responses", paste0("Driver Correlations (", tools::toTitleCase(meth), ")")),
        subtitle = "Stars show significance only; exported CSV tables contain r, R2, p, and n",
        x = "Driver",
        y = "Species",
        fill = "Correlation"
      ) +
      theme_pub(base_size = 13) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))

    ggsave(file.path(plot_dir, paste0("15_species_driver_heatmap_", meth, ".png")), p_sp, width = 12.8, height = 8.8, dpi = save_dpi)
  }
}

if (!file.exists(file.path(plot_dir, "15_species_driver_heatmap_pearson.png")) && exists("deer_driver_corr")) {
  for (meth in c("pearson", "spearman")) {
    d_fallback <- deer_driver_corr %>%
      filter(method == meth) %>%
      mutate(
        species = "Deer population",
        driver = factor(driver, levels = unique(driver)),
        species = factor(species),
        label = heatmap_sig_label(p)
      )

    p_sp_fallback <- ggplot(d_fallback, aes(driver, species, fill = r)) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_text(aes(label = label), size = 4.2, fontface = "bold") +
      scale_fill_gradient2(low = COL_HEAT_NEG, mid = 'white', high = COL_HEAT_POS, midpoint = 0, na.value = 'gray90') +
      labs(
        title = ww_title("Species Responses", paste0("Driver Correlations (", tools::toTitleCase(meth), ")")),
        subtitle = "Stars show significance only; exported CSV tables contain r, R2, p, and n",
        x = "Driver",
        y = "Response",
        fill = "Correlation"
      ) +
      theme_pub(base_size = 13) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))

    ggsave(file.path(plot_dir, paste0("15_species_driver_heatmap_", meth, ".png")), p_sp_fallback, width = 12.8, height = 4.8, dpi = save_dpi)
  }
}

# ==========================
# G) STAND AGE DYNAMICS
# ==========================

grow_bulk <- clean_rwl_for_detrend(safe_read_rwl(rwl_bulk_file))
bulk_mat <- as.matrix(grow_bulk)
bulk_years <- as.integer(rownames(grow_bulk))

est_year <- sapply(seq_len(ncol(bulk_mat)), function(i) {
  idx <- which(!is.na(bulk_mat[, i]))
  if (length(idx) == 0) return(NA_real_)
  bulk_years[min(idx)]
})

est_df <- data.frame(series = colnames(bulk_mat), est_year = as.numeric(est_year), stringsAsFactors = FALSE)
est_df <- est_df[!is.na(est_df$est_year), ]
write.csv(est_df, file.path(out_dir, "stand_establishment_years.csv"), row.names = FALSE)

bulk_long <- data.frame(
  year = rep(bulk_years, times = ncol(bulk_mat)),
  series = rep(colnames(bulk_mat), each = length(bulk_years)),
  rw = as.vector(bulk_mat),
  stringsAsFactors = FALSE
)
bulk_long <- bulk_long[!is.na(bulk_long$rw), ]

est_lookup <- setNames(est_df$est_year, est_df$series)
bulk_long$age <- bulk_long$year - est_lookup[bulk_long$series] + 1

annual_age <- bulk_long %>%
  group_by(year) %>%
  summarise(
    mean_age = mean(age, na.rm = TRUE),
    median_age = median(age, na.rm = TRUE),
    n_series = dplyr::n_distinct(series),
    mean_rw = mean(rw, na.rm = TRUE),
    .groups = "drop"
  )

annual_age <- annual_age %>%
  left_join(short_core %>% dplyr::select(all_of(c(year_var, short_core_vars))), by = c("year" = year_var)) %>%
  left_join(deer_df %>% dplyr::select(year, deer_population), by = "year")

write.csv(annual_age, file.path(out_dir, "annual_stand_age_summary.csv"), row.names = FALSE)

age_ts <- bind_rows(
  data.frame(year = annual_age$year, z = zfun(annual_age$mean_age), series = "Mean stand age"),
  data.frame(year = annual_age$year, z = zfun(annual_age$median_age), series = "Median stand age"),
  data.frame(year = annual_age$year, z = zfun(annual_age$RWI), series = "RWI"),
  data.frame(year = annual_age$year, z = zfun(annual_age$deer_population), series = "Deer population")
) %>% filter(!is.na(z))

p_age_ts <- ggplot(age_ts, aes(year, z, color = series)) +
  geom_line(linewidth = 1.1, alpha = 0.9) +
  scale_color_manual(values = c('Mean stand age' = '#6B7280', 'Median stand age' = '#9CA3AF', 'RWI' = COL_RWI, 'Deer population' = COL_DEER)) +
  labs(
    title = ww_title("Stand Age Dynamics", "Age + Growth + Deer Overlay"),
    subtitle = "All series standardized to z-scores",
    x = "Year",
    y = "Standardized value (z)",
    color = "Line legend"
  ) +
  theme_pub(base_size = 15)

ggsave(file.path(plot_dir, "16_stand_age_dynamics_overlay.png"), p_age_ts, width = 12.2, height = 6.0, dpi = save_dpi)

if ("RWI" %in% names(annual_age)) {
  age_rwi_p <- safe_cor(annual_age$mean_age, annual_age$RWI, method = "pearson")
  age_rwi_s <- safe_cor(annual_age$mean_age, annual_age$RWI, method = "spearman")

  p_age_scatter <- ggplot(annual_age, aes(mean_age, RWI)) +
    geom_point(color = COL_RWI, alpha = 0.9, size = 2.5) +
    geom_smooth(method = "lm", se = TRUE, color = COL_FIT, linewidth = 1.0) +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste0(
        "Pearson r = ", format(round(age_rwi_p$r, 3), nsmall = 3), " (", fmt_p(age_rwi_p$p), ")",
        "\nSpearman rho = ", format(round(age_rwi_s$r, 3), nsmall = 3), " (", fmt_p(age_rwi_s$p), ")"
      ),
      hjust = 1.05, vjust = 1.1, size = 3.4
    ) +
    labs(
      title = ww_title("Stand Mean Age vs Growth", "Correlation Regression"),
      x = "Mean stand age",
      y = "RWI"
    ) +
    theme_pub(base_size = 15)

  ggsave(file.path(plot_dir, "17_stand_mean_age_vs_rwi.png"), p_age_scatter, width = 9.0, height = 6.0, dpi = save_dpi)
}

p_est_hist <- ggplot(est_df, aes(est_year)) +
  geom_histogram(binwidth = 10, fill = COL_ETA, color = 'white') +
  labs(
    title = ww_title("Stand Establishment Structure", "Tree Establishment-Year Histogram"),
    x = "Establishment year",
    y = "Number of series"
  ) +
  theme_pub(base_size = 15)

ggsave(file.path(plot_dir, "18_stand_establishment_histogram.png"), p_est_hist, width = 10.6, height = 5.8, dpi = save_dpi)

# ==========================
# H) UPDATED COMMENTARY FILE
# ==========================

prism_year_min <- min(annual_prism$Year, na.rm = TRUE)
prism_year_max <- max(annual_prism$Year, na.rm = TRUE)
short_year_min <- min(short_core[[year_var]], na.rm = TRUE)
short_year_max <- max(short_core[[year_var]], na.rm = TRUE)

best_long <- long_lag_pear %>%
  filter(!is.na(r), lag == 0) %>%
  arrange(desc(abs(r))) %>%
  head(6)

best_short <- short_mat_pear %>%
  filter(var_x == response_var, var_y != response_var, !is.na(r)) %>%
  arrange(desc(abs(r))) %>%
  head(6)

best_deer <- NULL
if (exists("deer_driver_corr")) {
  best_deer <- deer_driver_corr %>%
    filter(method == "pearson", !is.na(r)) %>%
    arrange(desc(abs(r))) %>%
    head(6)
}

age_slope <- tryCatch(coef(lm(mean_age ~ year, data = annual_age))[2], error = function(e) NA_real_)

comment_lines <- c(
  "Wesselman Woods integrated analysis comments",
  paste0("Run timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("PRISM file used: ", basename(prism_file)),
  paste0("Short climate file used: ", short_file_used),
  paste0("PDSI source used: ", pdsi_source_file),
  paste0("NDVI source used: ", ndvi_source_file),
  paste0("ENSO source used: ", enso_source_file),
  paste0("PRISM annual coverage: ", prism_year_min, " to ", prism_year_max),
  paste0("MODIS-era coverage: ", short_year_min, " to ", short_year_max),
  paste0("Deer workbook columns detected: ", paste(names(deer_raw), collapse = ", ")),
  paste0("Incomplete PRISM years excluded from annual regressions: ", if (length(prism_incomplete_years) == 0) "none" else paste(prism_incomplete_years, collapse = ", ")),
  paste0("Climate-derived marker years used for shaded bands: ", if (length(marker_band_years) == 0) "none" else paste(marker_band_years, collapse = ", ")),
  "",
  "Top long-period PRISM vs RWI correlations (Pearson, lag 0):"
)

if (nrow(best_long) > 0) {
  for (i in seq_len(nrow(best_long))) {
    comment_lines <- c(comment_lines, paste0(
      "- ", best_long$variable[i], ": r = ", format(round(best_long$r[i], 3), nsmall = 3),
      ", R2 = ", format(round(best_long$r2[i], 3), nsmall = 3), ", ", fmt_p(best_long$p[i]), ", n = ", best_long$n[i]
    ))
  }
} else {
  comment_lines <- c(comment_lines, "- No valid long-period correlations were available.")
}

comment_lines <- c(comment_lines, "", "Top MODIS-era RWI relationships (Pearson matrix):")
if (nrow(best_short) > 0) {
  for (i in seq_len(nrow(best_short))) {
    comment_lines <- c(comment_lines, paste0(
      "- ", best_short$var_y[i], ": r = ", format(round(best_short$r[i], 3), nsmall = 3),
      ", R2 = ", format(round(best_short$r2[i], 3), nsmall = 3), ", ", fmt_p(best_short$p[i]), ", n = ", best_short$n[i]
    ))
  }
} else {
  comment_lines <- c(comment_lines, "- No valid short-period correlations were available.")
}

comment_lines <- c(comment_lines, "", "Top deer vs climate relationships (Pearson):")
if (!is.null(best_deer) && nrow(best_deer) > 0) {
  for (i in seq_len(nrow(best_deer))) {
    comment_lines <- c(comment_lines, paste0(
      "- ", best_deer$driver[i], ": r = ", format(round(best_deer$r[i], 3), nsmall = 3),
      ", R2 = ", format(round(best_deer$r2[i], 3), nsmall = 3), ", ", fmt_p(best_deer$p[i]), ", n = ", best_deer$n[i]
    ))
  }
} else {
  comment_lines <- c(comment_lines, "- Deer-vs-climate correlations were not available.")
}

comment_lines <- c(
  comment_lines,
  "",
  paste0("Stand mean-age trend slope (years of age per calendar year): ", format(round(age_slope, 4), nsmall = 4)),
  "",
  "Interpretation guardrails:",
  "- Pearson captures linear sensitivity; Spearman captures monotonic ranking behavior.",
  "- Long PRISM period improves climate signal stability relative to 25-year MODIS window.",
  "- Negative PDSI and SPEI values are valid drought-index values and were retained.",
  "- If species-level herbivory columns are added to the deer workbook, this script auto-builds species-driver heatmaps.",
  "- Stand-age dynamics should be interpreted alongside sample depth and chronology coverage."
)

writeLines(comment_lines, con = file.path(out_dir, "updated_analysis_comments.txt"))

cat("Done. Output folder:", out_dir, "\n")
cat("Primary plots folder:", plot_dir, "\n")
cat("PRISM file used:", prism_file, "\n")
cat("Short climate file used:", short_climate_active_file, "\n")
cat("PDSI source used:", pdsi_source_file, "\n")
cat("NDVI source used:", ndvi_source_file, "\n")
cat("ENSO source used:", enso_source_file, "\n")
cat("Created tables and plots for:\n")
cat("- Long PRISM Pearson/Spearman lag heatmaps + best-lag regression panels\n")
cat("- MODIS-era Pearson/Spearman matrices, lag heatmaps, and best-lag regression panels\n")
cat("- Faceted climate-variable time series with shaded marker years (climatology + MODIS-era)\n")
cat("- All-lag regression panel sets for RWI and deer population (0 to 2 years)\n")
cat("- Composite climate and predictor panels\n")
cat("- Deer vs growth/climate relationships\n")
cat("- Stand age dynamics and establishment structure\n")
cat("- Updated analysis comments text file\n")
cat("- climate_qc_flags.csv (incomplete PRISM years + impossible-value checks)\n")
cat("- RWI vs ENSO regression + standard wavelet transform outputs (when ENSO file is available)\n")



































