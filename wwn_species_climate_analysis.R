#!/usr/bin/env Rscript

# Species-by-species climate analysis built to mirror the multispecies workflow.
# Assumptions:
# - Short-period analysis uses species standard chronologies, matching the legacy short-period RWI.
# - Long-period analysis uses species residual chronologies, matching the prewhitened long chronology.

if (dir.exists(".Rlibs")) {
  .libPaths(c(normalizePath(".Rlibs"), .libPaths()))
}

req_pkgs <- c("ggplot2", "dplyr")
missing_pkgs <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them and re-run."
  )
}

suppressWarnings(suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
}))

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else "."
output_dir <- if (length(args) >= 2) args[[2]] else file.path(base_dir, "wwn_species_climate_outputs")
max_lag <- if (length(args) >= 3) as.integer(args[[3]]) else 2L
species_zoom_window_years <- 120L

save_dpi <- 420
COL_RWI <- "#7F0000"
COL_DEER <- "#A67C52"
COL_VPD <- "#4169E1"
COL_SPEI <- "#DAA520"
COL_PPT <- "#0B3C8A"
COL_ETA <- "#2E8B57"
COL_NDVI <- "#39D353"
COL_PDSI <- "#CC7722"
COL_TEMP <- "#7FA9D9"
COL_FIT <- "#2F3E56"
COL_HEAT_NEG <- "#0B3C8A"
COL_HEAT_POS <- "#CC7722"

ww_title <- function(subject, analysis) {
  paste0("Wesselman Woods - ", subject, " (", analysis, ")")
}

theme_pub <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.title.position = "plot",
      plot.margin = margin(14, 22, 18, 14),
      panel.spacing = grid::unit(0.9, "lines"),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10))
    )
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
  lab <- sig_code(p)
  lab[is.na(lab)] <- ""
  lab
}

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", format.pval(p, digits = 3, eps = 0.001))
}

pick_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

find_latest_table <- function(filename, preferred_subdir = NULL) {
  candidates <- character(0)

  if (!is.null(preferred_subdir)) {
    candidates <- c(candidates, file.path(base_dir, preferred_subdir, filename))
  }

  candidates <- c(
    candidates,
    file.path(base_dir, "wwn_handoff_outputs", filename),
    file.path(base_dir, "02_archive_data_exports", filename)
  )

  latest_dirs <- sort(list.files(base_dir, pattern = "^wwn_outputs_latest_", full.names = TRUE), decreasing = TRUE)
  if (length(latest_dirs) > 0) {
    candidates <- c(candidates, file.path(latest_dirs, "tables", filename))
  }

  pick_existing(unique(candidates))
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 4) {
    return(list(r = NA_real_, p = NA_real_, n = n, r2 = NA_real_))
  }

  x_ok <- x[ok]
  y_ok <- y[ok]
  if (sd(x_ok) == 0 || sd(y_ok) == 0) {
    return(list(r = NA_real_, p = NA_real_, n = n, r2 = NA_real_))
  }

  cor_args <- list(x = x_ok, y = y_ok, method = method)
  if (identical(method, "spearman")) cor_args$exact <- FALSE

  ct <- tryCatch(do.call(cor.test, cor_args), error = function(e) NULL)
  if (is.null(ct)) {
    return(list(r = NA_real_, p = NA_real_, n = n, r2 = NA_real_))
  }

  r <- as.numeric(ct$estimate)
  list(r = r, p = as.numeric(ct$p.value), n = n, r2 = r^2)
}

zfun <- function(x) {
  x <- as.numeric(x)
  if (sum(is.finite(x)) < 2) return(rep(NA_real_, length(x)))
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

lag_vec <- function(x, k) {
  if (k <= 0) return(x)
  c(rep(NA_real_, k), head(as.numeric(x), -k))
}

to_empty_lag_table <- function() {
  data.frame(
    species = character(0),
    window = character(0),
    response_series = character(0),
    driver = character(0),
    lag = integer(0),
    method = character(0),
    r = numeric(0),
    r2 = numeric(0),
    p = numeric(0),
    n = integer(0),
    year_start = integer(0),
    year_end = integer(0),
    stringsAsFactors = FALSE
  )
}

build_species_lag_table <- function(species_df, climate_df, year_col_species, year_col_climate, species_cols,
                                    predictors, max_lag = 2L, method = "pearson",
                                    window_label = "short", response_label = "standard") {
  if (length(species_cols) == 0 || length(predictors) == 0) return(to_empty_lag_table())

  out_rows <- vector("list", length(species_cols) * length(predictors) * (max_lag + 1L))
  k <- 1L

  for (sp in species_cols) {
    merged <- merge(
      species_df[, c(year_col_species, sp), drop = FALSE],
      climate_df[, c(year_col_climate, predictors), drop = FALSE],
      by.x = year_col_species,
      by.y = year_col_climate,
      all = FALSE,
      sort = TRUE
    )

    names(merged)[1:2] <- c("year", "response")
    merged <- merged[order(merged$year), , drop = FALSE]

    resp_years <- merged$year[is.finite(merged$response)]
    year_start <- if (length(resp_years) == 0) NA_integer_ else min(resp_years)
    year_end <- if (length(resp_years) == 0) NA_integer_ else max(resp_years)

    for (drv in predictors) {
      for (lg in 0:max_lag) {
        ct <- safe_cor(merged$response, lag_vec(merged[[drv]], lg), method = method)
        out_rows[[k]] <- data.frame(
          species = sp,
          window = window_label,
          response_series = response_label,
          driver = drv,
          lag = lg,
          method = method,
          r = ct$r,
          r2 = ct$r2,
          p = ct$p,
          n = ct$n,
          year_start = year_start,
          year_end = year_end,
          stringsAsFactors = FALSE
        )
        k <- k + 1L
      }
    }
  }

  bind_rows(out_rows[seq_len(k - 1L)])
}

extract_lag0_matrix <- function(lag_df) {
  lag_df %>%
    filter(lag == 0L) %>%
    arrange(species, driver)
}

extract_best_lag_matrix <- function(lag_df) {
  lag_df %>%
    filter(!is.na(r)) %>%
    group_by(species, driver) %>%
    arrange(desc(abs(r)), p, lag, .by_group = TRUE) %>%
    slice(1L) %>%
    ungroup() %>%
    arrange(species, driver)
}

extract_best_driver_summary <- function(lag_df) {
  lag_df %>%
    filter(!is.na(r)) %>%
    group_by(window, response_series, method, species) %>%
    arrange(desc(abs(r)), p, lag, .by_group = TRUE) %>%
    slice(1L) %>%
    ungroup() %>%
    mutate(direction = ifelse(r >= 0, "positive", "negative")) %>%
    arrange(window, method, desc(abs(r)), species)
}

make_matrix_corr <- function(df, vars, method = "pearson") {
  rows <- list()
  k <- 1L
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
      k <- k + 1L
    }
  }
  bind_rows(rows)
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
    scale_fill_gradient2(low = COL_HEAT_NEG, mid = "white", high = COL_HEAT_POS, midpoint = 0, na.value = "gray90") +
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

  invisible(d)
}

plot_species_heatmap <- function(df, title_txt, subtitle_txt, file_out) {
  if (nrow(df) == 0) return(FALSE)

  driver_order <- unique(df$driver)
  species_order <- rev(sort(unique(df$species)))
  show_lag_label <- "lag" %in% names(df) && any(is.finite(df$lag) & df$lag != 0)

  plot_df <- df %>%
    mutate(
      driver = factor(driver, levels = driver_order),
      species = factor(species, levels = species_order),
      sig_label = heatmap_sig_label(p),
      label = if (show_lag_label) {
        ifelse(
          is.na(lag),
          sig_label,
          ifelse(sig_label == "", paste0("L", lag), paste0("L", lag, "\n", sig_label))
        )
      } else {
        sig_label
      }
    )

  p <- ggplot(plot_df, aes(driver, species, fill = r)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = label), size = if (show_lag_label) 2.8 else 3.1, fontface = "bold", lineheight = 0.92) +
    scale_fill_gradient2(
      low = COL_HEAT_NEG,
      mid = "white",
      high = COL_HEAT_POS,
      midpoint = 0,
      na.value = "gray90"
    ) +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "Climate driver",
      y = "Species",
      fill = "Correlation"
    ) +
    theme_pub(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      axis.text.y = element_text(margin = margin(r = 6))
    )

  plot_width <- max(12.8, 0.72 * length(unique(plot_df$driver)) + 8.1)
  plot_height <- max(8.8, 0.5 * length(unique(plot_df$species)) + 4.0)
  ggsave(file_out, p, width = plot_width, height = plot_height, dpi = save_dpi)
  TRUE
}

climate_driver_colors <- function(vars) {
  base_cols <- c(
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
  out <- rep("NA", length(x))
  keep <- !is.na(x)
  out[keep] <- format(round(x[keep], digits), nsmall = digits)
  out
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
  k <- 1L

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
    k <- k + 1L
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
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = COL_FIT, linewidth = 0.95) +
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
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = COL_FIT, linewidth = 0.9) +
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

extract_top_driver_table <- function(lag_df, top_n = 5L) {
  lag_df %>%
    filter(!is.na(r)) %>%
    group_by(window, response_series, method, species) %>%
    arrange(desc(abs(r)), p, lag, .by_group = TRUE) %>%
    mutate(rank_within_species = row_number()) %>%
    filter(rank_within_species <= top_n) %>%
    ungroup() %>%
    mutate(direction = ifelse(r >= 0, "positive", "negative")) %>%
    arrange(window, method, species, rank_within_species)
}

make_species_window_df <- function(species_df, climate_df, species_name, predictors, climate_year_col) {
  out <- merge(
    species_df[, c("year", species_name), drop = FALSE],
    climate_df[, c(climate_year_col, predictors), drop = FALSE],
    by.x = "year",
    by.y = climate_year_col,
    all = FALSE,
    sort = TRUE
  )
  names(out)[2] <- "response"
  out[order(out$year), , drop = FALSE]
}

write_species_regression_panels <- function(
  species_name,
  merged_df,
  lag_df,
  predictors,
  response_label,
  title_prefix,
  plot_subdir,
  stats_subdir,
  file_prefix
) {
  if (nrow(merged_df) == 0 || nrow(lag_df) == 0) return(invisible(NULL))

  plot_df <- merged_df
  names(plot_df)[names(plot_df) == "response"] <- species_name
  lag_src <- lag_df %>% rename(variable = driver)
  build_best_lag_regression_panels(
    df = plot_df,
    year_col = "year",
    response = species_name,
    predictors = predictors,
    lag_table = lag_src,
    title_txt = ww_title(title_prefix, "Best-Lag Regression Panels"),
    subtitle_txt = "Panels use the strongest absolute Pearson lag per predictor for this species.",
    file_out = file.path(plot_subdir, paste0(file_prefix, "_bestlag_regression_panels.png")),
    stats_out = file.path(stats_subdir, paste0(file_prefix, "_bestlag_regression_stats.csv")),
    ncol = 3
  )

  build_all_lag_regression_panels(
    df = plot_df,
    year_col = "year",
    response = species_name,
    predictors = predictors,
    max_lag = max_lag,
    title_txt = ww_title(title_prefix, "All-Lag Regression Panels"),
    subtitle_txt = "Every panel shows one predictor-lag combination against the species chronology with Pearson fit and statistics.",
    file_out = file.path(plot_subdir, paste0(file_prefix, "_all_lag_regression_panels.png")),
    stats_out = file.path(stats_subdir, paste0(file_prefix, "_all_lag_regression_stats.csv")),
    ncol = 4
  )

  invisible(response_label)
}

write_species_subtables <- function(lag_df, species_dir, prefix_base) {
  if (nrow(lag_df) == 0) return(invisible(NULL))
  for (sp in unique(lag_df$species)) {
    sp_dir <- file.path(species_dir, sp)
    dir.create(sp_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(
      lag_df[lag_df$species == sp, , drop = FALSE],
      file.path(sp_dir, paste0(sp, "_", prefix_base, ".csv")),
      row.names = FALSE,
      na = ""
    )
  }
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

plot_species_timeseries_with_markers <- function(
  wide_df,
  species_cols,
  marker_years,
  title_txt,
  subtitle_txt,
  file_out,
  year_min = NULL,
  year_max = NULL,
  line_color = COL_RWI
) {
  rows <- bind_rows(lapply(species_cols, function(sp) {
    data.frame(
      year = wide_df$year,
      species = sp,
      value = as.numeric(wide_df[[sp]]),
      stringsAsFactors = FALSE
    )
  }))

  rows <- rows[is.finite(rows$year), , drop = FALSE]
  if (!is.null(year_min)) rows <- rows[rows$year >= year_min, , drop = FALSE]
  if (!is.null(year_max)) rows <- rows[rows$year <= year_max, , drop = FALSE]
  rows <- rows[!is.na(rows$value), , drop = FALSE]
  if (nrow(rows) == 0) return(invisible(NULL))

  band_df <- make_marker_bands(intersect(marker_years, rows$year))

  p <- ggplot(rows, aes(year, value)) +
    geom_rect(
      data = band_df,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "grey80",
      alpha = 0.35
    ) +
    geom_line(color = line_color, linewidth = 0.9, na.rm = TRUE) +
    geom_point(color = line_color, size = 0.8, alpha = 0.65, na.rm = TRUE) +
    facet_wrap(~species, scales = "free_y", ncol = 3) +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "Year",
      y = "Chronology value"
    ) +
    theme_pub(base_size = 15) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  plot_height <- max(8.0, 2.35 * ceiling(length(species_cols) / 3))
  ggsave(file_out, p, width = 13.8, height = plot_height, dpi = save_dpi)

  invisible(rows)
}

build_marker_lookup <- function(marker_df) {
  if (nrow(marker_df) == 0 || !("Year" %in% names(marker_df))) {
    return(data.frame(
      year = integer(0),
      is_flood = logical(0),
      is_drought = logical(0),
      event_group = character(0),
      stringsAsFactors = FALSE
    ))
  }

  year_vals <- as.integer(marker_df$Year)
  event_txt <- if ("event_type" %in% names(marker_df)) tolower(as.character(marker_df$event_type)) else rep("", nrow(marker_df))
  is_flood <- if ("is_flood" %in% names(marker_df)) as.logical(marker_df$is_flood) else grepl("flood", event_txt)
  is_drought <- if ("is_drought" %in% names(marker_df)) as.logical(marker_df$is_drought) else grepl("drought", event_txt)

  lookup <- data.frame(
    year = year_vals,
    is_flood = is_flood,
    is_drought = is_drought,
    stringsAsFactors = FALSE
  )
  lookup <- lookup[is.finite(lookup$year), , drop = FALSE]
  if (nrow(lookup) == 0) {
    return(data.frame(
      year = integer(0),
      is_flood = logical(0),
      is_drought = logical(0),
      event_group = character(0),
      stringsAsFactors = FALSE
    ))
  }

  lookup <- lookup %>%
    group_by(year) %>%
    summarise(
      is_flood = any(is_flood, na.rm = TRUE),
      is_drought = any(is_drought, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      event_group = case_when(
        is_flood & is_drought ~ "Flood + Drought marker",
        is_flood ~ "Flood marker",
        is_drought ~ "Drought marker",
        TRUE ~ "Non-marker"
      )
    )

  as.data.frame(lookup, stringsAsFactors = FALSE)
}

make_species_value_long <- function(wide_df, species_cols) {
  out <- bind_rows(lapply(species_cols, function(sp) {
    data.frame(
      year = wide_df$year,
      species = sp,
      value = as.numeric(wide_df[[sp]]),
      stringsAsFactors = FALSE
    )
  }))

  out <- out[is.finite(out$year), , drop = FALSE]
  out %>%
    group_by(species) %>%
    mutate(value_z = zfun(value)) %>%
    ungroup()
}

plot_species_event_boxplots <- function(event_df, title_txt, subtitle_txt, file_out) {
  plot_df <- event_df %>% filter(!is.na(value_z))
  if (nrow(plot_df) == 0) return(invisible(NULL))

  lvl <- c("Flood marker", "Drought marker", "Flood + Drought marker", "Non-marker")
  lvl <- lvl[lvl %in% unique(plot_df$event_group)]
  col_map <- c(
    "Flood marker" = COL_PPT,
    "Drought marker" = COL_PDSI,
    "Flood + Drought marker" = COL_DEER,
    "Non-marker" = "grey65"
  )

  plot_df$event_group <- factor(plot_df$event_group, levels = lvl)

  p <- ggplot(plot_df, aes(event_group, value_z, fill = event_group)) +
    geom_hline(yintercept = 0, color = "grey35", linewidth = 0.4) +
    geom_boxplot(outlier.shape = NA, width = 0.72, alpha = 0.88) +
    geom_point(
      position = position_jitter(width = 0.12, height = 0),
      size = 0.9,
      alpha = 0.38,
      color = "black"
    ) +
    facet_wrap(~species, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = col_map[lvl], drop = FALSE) +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "Marker-year class",
      y = "Species chronology (z within species)",
      fill = "Year class"
    ) +
    theme_pub(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      strip.text = element_text(face = "bold")
    )

  plot_height <- max(8.5, 2.45 * ceiling(length(unique(plot_df$species)) / 3))
  ggsave(file_out, p, width = 14.2, height = plot_height, dpi = save_dpi)

  invisible(plot_df)
}

extract_top_driver_choices <- function(best_lag_df, top_n = 2L) {
  best_lag_df %>%
    group_by(species) %>%
    arrange(desc(abs(r)), p, lag, .by_group = TRUE) %>%
    mutate(driver_rank = row_number()) %>%
    filter(driver_rank <= top_n) %>%
    ungroup()
}

build_species_overlay_df <- function(species_df, climate_df, driver_choices, climate_year_col = "Year") {
  if (nrow(driver_choices) == 0) {
    return(data.frame(
      year = integer(0),
      species = character(0),
      series_role = character(0),
      series_label = character(0),
      color_key = character(0),
      value = numeric(0),
      z = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  rows <- list()
  ann_rows <- list()
  k <- 1L
  j <- 1L

  for (sp in unique(driver_choices$species)) {
    sp_choices <- driver_choices[driver_choices$species == sp, , drop = FALSE]
    merge_vars <- unique(c(climate_year_col, sp_choices$driver))
    merged <- merge(
      species_df[, c("year", sp), drop = FALSE],
      climate_df[, merge_vars, drop = FALSE],
      by.x = "year",
      by.y = climate_year_col,
      all = FALSE,
      sort = TRUE
    )
    names(merged)[2] <- "chronology"
    merged <- merged[order(merged$year), , drop = FALSE]

    rows[[k]] <- data.frame(
      year = merged$year,
      species = sp,
      series_role = "Chronology",
      series_label = "Chronology",
      color_key = "Chronology",
      value = merged$chronology,
      stringsAsFactors = FALSE
    )
    k <- k + 1L

    ann_text <- character(0)
    for (ii in seq_len(nrow(sp_choices))) {
      drv <- sp_choices$driver[ii]
      lg <- as.integer(sp_choices$lag[ii])
      drv_vals <- if (lg == 0) merged[[drv]] else lag_vec(merged[[drv]], lg)
      role <- paste0("Top driver ", ii)
      label <- paste0(ii, ". ", drv, " (lag ", lg, ")")

      rows[[k]] <- data.frame(
        year = merged$year,
        species = sp,
        series_role = role,
        series_label = label,
        color_key = role,
        value = drv_vals,
        stringsAsFactors = FALSE
      )
      k <- k + 1L

      ann_text <- c(ann_text, paste0(label, ", r = ", fmt_num(sp_choices$r[ii])))
    }

    ann_rows[[j]] <- data.frame(
      species = sp,
      label = paste(ann_text, collapse = "\n"),
      stringsAsFactors = FALSE
    )
    j <- j + 1L
  }

  overlay_df <- bind_rows(rows) %>%
    group_by(species, series_role, series_label, color_key) %>%
    mutate(z = zfun(value)) %>%
    ungroup()

  list(
    series = overlay_df,
    annotations = bind_rows(ann_rows)
  )
}

plot_species_top_driver_overlays <- function(
  overlay_obj,
  marker_years,
  title_txt,
  subtitle_txt,
  file_out,
  year_min = NULL,
  year_max = NULL
) {
  plot_df <- overlay_obj$series
  ann_df <- overlay_obj$annotations
  if (nrow(plot_df) == 0) return(invisible(NULL))

  if (!is.null(year_min)) plot_df <- plot_df[plot_df$year >= year_min, , drop = FALSE]
  if (!is.null(year_max)) plot_df <- plot_df[plot_df$year <= year_max, , drop = FALSE]
  plot_df <- plot_df[!is.na(plot_df$z), , drop = FALSE]
  if (nrow(plot_df) == 0) return(invisible(NULL))

  ann_pos <- build_label_positions(plot_df, "species", "year", "z", x_frac = 0.03, y_frac = 0.08)
  names(ann_pos)[names(ann_pos) == "group_id"] <- "species"
  ann_df <- left_join(ann_df, ann_pos, by = "species")

  band_df <- make_marker_bands(intersect(marker_years, plot_df$year))
  col_map <- c(
    "Chronology" = COL_RWI,
    "Top driver 1" = COL_PPT,
    "Top driver 2" = COL_SPEI
  )
  col_map <- col_map[names(col_map) %in% unique(plot_df$color_key)]

  p <- ggplot(plot_df, aes(year, z, color = color_key)) +
    geom_rect(
      data = band_df,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "grey80",
      alpha = 0.30
    ) +
    geom_line(linewidth = 0.85, na.rm = TRUE) +
    facet_wrap(~species, ncol = 3, scales = "free_x") +
    scale_color_manual(values = col_map) +
    geom_text(
      data = ann_df,
      aes(x = label_x, y = label_y, label = label),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      color = COL_FIT,
      size = 2.5,
      lineheight = 0.94
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "Year",
      y = "Standardized value (z)",
      color = "Series"
    ) +
    theme_pub(base_size = 14) +
    theme(strip.text = element_text(face = "bold"))

  plot_height <- max(8.9, 2.55 * ceiling(length(unique(plot_df$species)) / 3))
  ggsave(file_out, p, width = 14.8, height = plot_height, dpi = save_dpi)

  invisible(plot_df)
}

plot_species_response_rankings <- function(rank_df, title_txt, subtitle_txt, file_out) {
  plot_df <- rank_df
  if (nrow(plot_df) == 0) return(invisible(NULL))

  plot_df$label <- paste0(plot_df$driver_rank, ". ", plot_df$driver, " (lag ", plot_df$lag, ")")
  plot_df$label <- factor(plot_df$label, levels = rev(unique(plot_df$label)))

  xmax <- max(abs(plot_df$r), na.rm = TRUE)
  if (!is.finite(xmax)) xmax <- 1

  p <- ggplot(plot_df, aes(abs(r), label, color = r)) +
    geom_segment(aes(x = 0, xend = abs(r), y = label, yend = label), linewidth = 1.0, alpha = 0.8) +
    geom_point(size = 3.0) +
    geom_text(aes(label = paste0("r = ", fmt_num(r))), hjust = -0.15, size = 2.7, color = COL_FIT) +
    facet_wrap(~species, ncol = 3, scales = "free_y") +
    scale_color_gradient2(low = COL_HEAT_NEG, mid = "white", high = COL_HEAT_POS, midpoint = 0) +
    coord_cartesian(xlim = c(0, xmax * 1.42), clip = "off") +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "|Pearson r|",
      y = "Driver ranking",
      color = "Signed r"
    ) +
    theme_pub(base_size = 14) +
    theme(strip.text = element_text(face = "bold"))

  plot_height <- max(8.9, 2.55 * ceiling(length(unique(plot_df$species)) / 3))
  ggsave(file_out, p, width = 14.8, height = plot_height, dpi = save_dpi)

  invisible(plot_df)
}

plot_species_event_year_heatmap <- function(event_df, title_txt, subtitle_txt, file_out) {
  plot_df <- event_df %>% filter(event_group != "Non-marker", !is.na(value_z))
  if (nrow(plot_df) == 0) return(invisible(NULL))

  plot_df$event_panel <- factor(
    plot_df$event_group,
    levels = c("Flood marker", "Drought marker", "Flood + Drought marker")
  )
  plot_df <- plot_df[!is.na(plot_df$event_panel), , drop = FALSE]
  plot_df$species <- factor(plot_df$species, levels = rev(sort(unique(plot_df$species))))
  plot_df$year_f <- factor(plot_df$year, levels = sort(unique(plot_df$year)))

  plot_width <- max(13.5, 0.22 * length(unique(plot_df$year)) + 6)

  p <- ggplot(plot_df, aes(year_f, species, fill = value_z)) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_grid(. ~ event_panel, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = COL_HEAT_NEG, mid = "white", high = COL_HEAT_POS, midpoint = 0, na.value = "gray90") +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "Marker year",
      y = "Species",
      fill = "Growth anomaly\n(z)"
    ) +
    theme_pub(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.text = element_text(face = "bold")
    )

  ggsave(file_out, p, width = plot_width, height = 6.8, dpi = save_dpi)

  invisible(plot_df)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(output_dir, "plots")
species_out_dir <- file.path(output_dir, "species")
species_plot_root <- file.path(plot_dir, "species")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(species_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(species_plot_root, recursive = TRUE, showWarnings = FALSE)

species_standard_path <- file.path(base_dir, "wwn_species_deer_outputs", "species", "species_standard_chronologies_wide.csv")
species_residual_path <- file.path(base_dir, "wwn_species_deer_outputs", "species", "species_residual_chronologies_wide.csv")

if (!file.exists(species_standard_path)) stop("Missing species standard chronology table: ", species_standard_path)
if (!file.exists(species_residual_path)) stop("Missing species residual chronology table: ", species_residual_path)

short_climate_path <- find_latest_table("short_modis_with_composite.csv")
long_climate_path <- find_latest_table("long_prism_rwi_annual_merge.csv")
marker_years_path <- find_latest_table("flood_drought_marker_years.csv")

if (is.na(short_climate_path)) stop("Could not find short_modis_with_composite.csv in expected output folders.")
if (is.na(long_climate_path)) stop("Could not find long_prism_rwi_annual_merge.csv in expected output folders.")

species_standard <- read.csv(species_standard_path, check.names = FALSE, stringsAsFactors = FALSE)
species_residual <- read.csv(species_residual_path, check.names = FALSE, stringsAsFactors = FALSE)
short_df <- read.csv(short_climate_path, check.names = FALSE, stringsAsFactors = FALSE)
long_df <- read.csv(long_climate_path, check.names = FALSE, stringsAsFactors = FALSE)
marker_years_df <- if (!is.na(marker_years_path)) read.csv(marker_years_path, check.names = FALSE, stringsAsFactors = FALSE) else data.frame()

names(species_standard)[1] <- "year"
names(species_residual)[1] <- "year"

short_df$year <- as.integer(short_df$year)
long_df$Year <- as.integer(long_df$Year)
species_standard$year <- as.integer(species_standard$year)
species_residual$year <- as.integer(species_residual$year)

species_cols <- intersect(names(species_standard)[-1], names(species_residual)[-1])
if (length(species_cols) == 0) stop("No overlapping species columns found between standard and residual chronology tables.")

marker_years <- integer(0)
if (nrow(marker_years_df) > 0 && "Year" %in% names(marker_years_df)) {
  marker_years <- sort(unique(as.integer(marker_years_df$Year[is.finite(as.integer(marker_years_df$Year))])))
}

short_drivers <- intersect(c("ClimateComposite", "ENSO_NINO34", "EVT", "PPT", "TMEAN", "PDSI", "SPEI", "NDVI"), names(short_df))
long_drivers <- intersect(c("PPT_mm", "TMIN_C", "TMEAN_C", "TMAX_C", "VPDMIN", "VPDMAX", "SPEI6_gs", "SPEI6_min", "PDSI", "ENSO_NINO34", "ClimateComposite"), names(long_df))

if (length(short_drivers) == 0) stop("No short-period climate drivers were found in: ", short_climate_path)
if (length(long_drivers) == 0) stop("No long-period climate drivers were found in: ", long_climate_path)

short_climate_core <- short_df[, c("year", short_drivers), drop = FALSE]
long_climate_core <- long_df[, c("Year", long_drivers), drop = FALSE]
species_standard_core <- species_standard[, c("year", species_cols), drop = FALSE]
species_residual_core <- species_residual[, c("year", species_cols), drop = FALSE]

short_clim_mat_pear <- make_matrix_corr(short_climate_core, short_drivers, method = "pearson")
short_clim_mat_spear <- make_matrix_corr(short_climate_core, short_drivers, method = "spearman")
long_clim_mat_pear <- make_matrix_corr(long_climate_core, long_drivers, method = "pearson")
long_clim_mat_spear <- make_matrix_corr(long_climate_core, long_drivers, method = "spearman")
species_std_mat_pear <- make_matrix_corr(species_standard_core, species_cols, method = "pearson")
species_std_mat_spear <- make_matrix_corr(species_standard_core, species_cols, method = "spearman")
species_res_mat_pear <- make_matrix_corr(species_residual_core, species_cols, method = "pearson")
species_res_mat_spear <- make_matrix_corr(species_residual_core, species_cols, method = "spearman")

write.csv(short_clim_mat_pear, file.path(output_dir, "climate_short_self_matrix_pearson.csv"), row.names = FALSE, na = "")
write.csv(short_clim_mat_spear, file.path(output_dir, "climate_short_self_matrix_spearman.csv"), row.names = FALSE, na = "")
write.csv(long_clim_mat_pear, file.path(output_dir, "climate_long_self_matrix_pearson.csv"), row.names = FALSE, na = "")
write.csv(long_clim_mat_spear, file.path(output_dir, "climate_long_self_matrix_spearman.csv"), row.names = FALSE, na = "")
write.csv(species_std_mat_pear, file.path(output_dir, "species_species_standard_matrix_pearson.csv"), row.names = FALSE, na = "")
write.csv(species_std_mat_spear, file.path(output_dir, "species_species_standard_matrix_spearman.csv"), row.names = FALSE, na = "")
write.csv(species_res_mat_pear, file.path(output_dir, "species_species_residual_matrix_pearson.csv"), row.names = FALSE, na = "")
write.csv(species_res_mat_spear, file.path(output_dir, "species_species_residual_matrix_spearman.csv"), row.names = FALSE, na = "")

invisible(plot_matrix_heatmap(
  short_clim_mat_pear,
  ww_title("Short-Period Climate Variables", "Climate-to-Climate Correlation Matrix (Pearson)"),
  "Correlation matrix among the MODIS-era climate and productivity variables used in the species analysis.",
  file.path(plot_dir, "16_climate_short_self_matrix_pearson.png")
))
invisible(plot_matrix_heatmap(
  short_clim_mat_spear,
  ww_title("Short-Period Climate Variables", "Climate-to-Climate Correlation Matrix (Spearman)"),
  "Correlation matrix among the MODIS-era climate and productivity variables used in the species analysis.",
  file.path(plot_dir, "17_climate_short_self_matrix_spearman.png")
))
invisible(plot_matrix_heatmap(
  long_clim_mat_pear,
  ww_title("Long-Period Climate Variables", "Climate-to-Climate Correlation Matrix (Pearson)"),
  "Correlation matrix among the long-term PRISM, drought, ENSO, and climate-composite variables.",
  file.path(plot_dir, "18_climate_long_self_matrix_pearson.png")
))
invisible(plot_matrix_heatmap(
  long_clim_mat_spear,
  ww_title("Long-Period Climate Variables", "Climate-to-Climate Correlation Matrix (Spearman)"),
  "Correlation matrix among the long-term PRISM, drought, ENSO, and climate-composite variables.",
  file.path(plot_dir, "19_climate_long_self_matrix_spearman.png")
))
invisible(plot_matrix_heatmap(
  species_std_mat_pear,
  ww_title("Species Standard Chronologies", "Species-to-Species Correlation Matrix (Pearson)"),
  "Correlation matrix among species standard chronologies from the .rwl-derived species runs.",
  file.path(plot_dir, "20_species_species_standard_matrix_pearson.png")
))
invisible(plot_matrix_heatmap(
  species_std_mat_spear,
  ww_title("Species Standard Chronologies", "Species-to-Species Correlation Matrix (Spearman)"),
  "Correlation matrix among species standard chronologies from the .rwl-derived species runs.",
  file.path(plot_dir, "21_species_species_standard_matrix_spearman.png")
))
invisible(plot_matrix_heatmap(
  species_res_mat_pear,
  ww_title("Species Residual Chronologies", "Species-to-Species Correlation Matrix (Pearson)"),
  "Correlation matrix among species residual chronologies from the .rwl-derived species runs.",
  file.path(plot_dir, "22_species_species_residual_matrix_pearson.png")
))
invisible(plot_matrix_heatmap(
  species_res_mat_spear,
  ww_title("Species Residual Chronologies", "Species-to-Species Correlation Matrix (Spearman)"),
  "Correlation matrix among species residual chronologies from the .rwl-derived species runs.",
  file.path(plot_dir, "23_species_species_residual_matrix_spearman.png")
))

short_pear <- build_species_lag_table(
  species_df = species_standard,
  climate_df = short_df,
  year_col_species = "year",
  year_col_climate = "year",
  species_cols = species_cols,
  predictors = short_drivers,
  max_lag = max_lag,
  method = "pearson",
  window_label = "short",
  response_label = "standard"
)
short_spear <- build_species_lag_table(
  species_df = species_standard,
  climate_df = short_df,
  year_col_species = "year",
  year_col_climate = "year",
  species_cols = species_cols,
  predictors = short_drivers,
  max_lag = max_lag,
  method = "spearman",
  window_label = "short",
  response_label = "standard"
)
long_pear <- build_species_lag_table(
  species_df = species_residual,
  climate_df = long_df,
  year_col_species = "year",
  year_col_climate = "Year",
  species_cols = species_cols,
  predictors = long_drivers,
  max_lag = max_lag,
  method = "pearson",
  window_label = "long",
  response_label = "residual"
)
long_spear <- build_species_lag_table(
  species_df = species_residual,
  climate_df = long_df,
  year_col_species = "year",
  year_col_climate = "Year",
  species_cols = species_cols,
  predictors = long_drivers,
  max_lag = max_lag,
  method = "spearman",
  window_label = "long",
  response_label = "residual"
)

write.csv(short_pear, file.path(output_dir, "species_short_lag_correlations_pearson.csv"), row.names = FALSE, na = "")
write.csv(short_spear, file.path(output_dir, "species_short_lag_correlations_spearman.csv"), row.names = FALSE, na = "")
write.csv(long_pear, file.path(output_dir, "species_long_lag_correlations_pearson.csv"), row.names = FALSE, na = "")
write.csv(long_spear, file.path(output_dir, "species_long_lag_correlations_spearman.csv"), row.names = FALSE, na = "")

write_species_subtables(short_pear, species_out_dir, "short_lag_correlations_pearson")
write_species_subtables(short_spear, species_out_dir, "short_lag_correlations_spearman")
write_species_subtables(long_pear, species_out_dir, "long_lag_correlations_pearson")
write_species_subtables(long_spear, species_out_dir, "long_lag_correlations_spearman")

short_lag0_pear <- extract_lag0_matrix(short_pear)
short_lag0_spear <- extract_lag0_matrix(short_spear)
long_lag0_pear <- extract_lag0_matrix(long_pear)
long_lag0_spear <- extract_lag0_matrix(long_spear)

short_best_pear <- extract_best_lag_matrix(short_pear)
short_best_spear <- extract_best_lag_matrix(short_spear)
long_best_pear <- extract_best_lag_matrix(long_pear)
long_best_spear <- extract_best_lag_matrix(long_spear)

write.csv(short_lag0_pear, file.path(output_dir, "species_short_climate_matrix_lag0_pearson.csv"), row.names = FALSE, na = "")
write.csv(short_lag0_spear, file.path(output_dir, "species_short_climate_matrix_lag0_spearman.csv"), row.names = FALSE, na = "")
write.csv(long_lag0_pear, file.path(output_dir, "species_long_climate_matrix_lag0_pearson.csv"), row.names = FALSE, na = "")
write.csv(long_lag0_spear, file.path(output_dir, "species_long_climate_matrix_lag0_spearman.csv"), row.names = FALSE, na = "")

write.csv(short_best_pear, file.path(output_dir, "species_short_climate_matrix_bestlag_pearson.csv"), row.names = FALSE, na = "")
write.csv(short_best_spear, file.path(output_dir, "species_short_climate_matrix_bestlag_spearman.csv"), row.names = FALSE, na = "")
write.csv(long_best_pear, file.path(output_dir, "species_long_climate_matrix_bestlag_pearson.csv"), row.names = FALSE, na = "")
write.csv(long_best_spear, file.path(output_dir, "species_long_climate_matrix_bestlag_spearman.csv"), row.names = FALSE, na = "")

full_year_max <- suppressWarnings(max(species_standard$year, na.rm = TRUE))
recent_year_min <- if (is.finite(full_year_max)) full_year_max - species_zoom_window_years + 1L else NULL

invisible(plot_species_timeseries_with_markers(
  wide_df = species_standard,
  species_cols = species_cols,
  marker_years = marker_years,
  title_txt = ww_title("Species Standard Chronologies", "Faceted Time Series with Flood/Drought Marker Years"),
  subtitle_txt = "Gray bands mark climate-derived flood and drought marker years from the long-period analysis.",
  file_out = file.path(plot_dir, "09_species_standard_chronologies_marker_bands_fullrange.png")
))

if (!is.null(recent_year_min)) {
  invisible(plot_species_timeseries_with_markers(
    wide_df = species_standard,
    species_cols = species_cols,
    marker_years = marker_years,
    title_txt = ww_title("Species Standard Chronologies", paste0("Recent ", species_zoom_window_years, "-Year Window with Marker Years")),
    subtitle_txt = "Gray bands mark climate-derived flood and drought marker years within the recent shared window.",
    file_out = file.path(plot_dir, "10_species_standard_chronologies_marker_bands_recent_window.png"),
    year_min = recent_year_min,
    year_max = full_year_max
  ))
}

best_summary <- bind_rows(
  extract_best_driver_summary(short_pear),
  extract_best_driver_summary(short_spear),
  extract_best_driver_summary(long_pear),
  extract_best_driver_summary(long_spear)
)
top_driver_table <- bind_rows(
  extract_top_driver_table(short_pear, top_n = 5L),
  extract_top_driver_table(short_spear, top_n = 5L),
  extract_top_driver_table(long_pear, top_n = 5L),
  extract_top_driver_table(long_spear, top_n = 5L)
)
write.csv(best_summary, file.path(output_dir, "species_best_climate_matches.csv"), row.names = FALSE, na = "")
write.csv(top_driver_table, file.path(output_dir, "species_top_climate_responses.csv"), row.names = FALSE, na = "")

marker_lookup <- build_marker_lookup(marker_years_df)
species_standard_long <- make_species_value_long(species_standard, species_cols) %>%
  left_join(marker_lookup, by = c("year" = "year"))
species_standard_long$event_group[is.na(species_standard_long$event_group)] <- "Non-marker"

event_summary <- species_standard_long %>%
  group_by(species, event_group) %>%
  summarise(
    n_years = sum(!is.na(value)),
    mean_growth = if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE),
    mean_growth_z = if (all(is.na(value_z))) NA_real_ else mean(value_z, na.rm = TRUE),
    median_growth_z = if (all(is.na(value_z))) NA_real_ else median(value_z, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(species_standard_long, file.path(output_dir, "species_marker_event_growth_values.csv"), row.names = FALSE, na = "")
write.csv(event_summary, file.path(output_dir, "species_marker_event_growth_summary.csv"), row.names = FALSE, na = "")

invisible(plot_species_event_boxplots(
  species_standard_long,
  ww_title("Species Growth by Marker-Year Class", "Flood vs Drought Boxplots"),
  "Each panel compares species chronology anomalies during climate-derived flood years, drought years, and non-marker years.",
  file.path(plot_dir, "11_species_marker_event_boxplots.png")
))

overlay_choices_long <- extract_top_driver_choices(long_best_pear, top_n = 2L)
ranking_choices_long <- extract_top_driver_choices(long_best_pear, top_n = 3L)
overlay_obj_long <- build_species_overlay_df(
  species_df = species_residual,
  climate_df = long_df,
  driver_choices = overlay_choices_long,
  climate_year_col = "Year"
)

write.csv(overlay_choices_long, file.path(output_dir, "species_long_top_driver_choices.csv"), row.names = FALSE, na = "")
write.csv(overlay_obj_long$series, file.path(output_dir, "species_long_top_driver_overlay_series.csv"), row.names = FALSE, na = "")
write.csv(ranking_choices_long, file.path(output_dir, "species_long_response_rankings_top3.csv"), row.names = FALSE, na = "")

invisible(plot_species_top_driver_overlays(
  overlay_obj = overlay_obj_long,
  marker_years = marker_years,
  title_txt = ww_title("Species Residual Chronologies", "Top Climate Driver Overlays"),
  subtitle_txt = "Each panel shows the species chronology with its top two long-window Pearson climate drivers (best lag for each driver).",
  file_out = file.path(plot_dir, "12_species_long_top_driver_overlays_fullrange.png")
))

if (!is.null(recent_year_min)) {
  invisible(plot_species_top_driver_overlays(
    overlay_obj = overlay_obj_long,
    marker_years = marker_years,
    title_txt = ww_title("Species Residual Chronologies", paste0("Top Climate Driver Overlays, Recent ", species_zoom_window_years, "-Year Window")),
    subtitle_txt = "Each panel shows the species chronology with its top two long-window Pearson climate drivers in the recent shared window.",
    file_out = file.path(plot_dir, "13_species_long_top_driver_overlays_recent_window.png"),
    year_min = recent_year_min,
    year_max = full_year_max
  ))
}

invisible(plot_species_response_rankings(
  ranking_choices_long,
  ww_title("Species Climate Responses", "Top Three Long-Window Pearson Rankings"),
  "Panels rank the strongest long-window climate drivers for each species by absolute Pearson correlation.",
  file.path(plot_dir, "14_species_long_response_rankings_top3.png")
))

event_heat_df <- species_standard_long %>% filter(event_group != "Non-marker")
write.csv(event_heat_df, file.path(output_dir, "species_marker_year_growth_heatmap_values.csv"), row.names = FALSE, na = "")

invisible(plot_species_event_year_heatmap(
  event_heat_df,
  ww_title("Species Growth During Marker Years", "Event-Year Heatmap"),
  "Tile color shows within-species growth anomaly during climate-derived flood and drought marker years.",
  file.path(plot_dir, "15_species_marker_year_growth_heatmap.png")
))

for (sp in species_cols) {
  sp_plot_dir <- file.path(species_plot_root, sp)
  sp_stats_dir <- file.path(species_out_dir, sp)
  dir.create(sp_plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(sp_stats_dir, recursive = TRUE, showWarnings = FALSE)

  short_merge_sp <- make_species_window_df(
    species_df = species_standard,
    climate_df = short_df,
    species_name = sp,
    predictors = short_drivers,
    climate_year_col = "year"
  )
  long_merge_sp <- make_species_window_df(
    species_df = species_residual,
    climate_df = long_df,
    species_name = sp,
    predictors = long_drivers,
    climate_year_col = "Year"
  )

  write_species_regression_panels(
    species_name = sp,
    merged_df = short_merge_sp,
    lag_df = short_pear %>% filter(species == sp),
    predictors = short_drivers,
    response_label = "standard chronology",
    title_prefix = paste0(sp, " vs MODIS-Era Climate"),
    plot_subdir = sp_plot_dir,
    stats_subdir = sp_stats_dir,
    file_prefix = paste0(sp, "_short")
  )

  write_species_regression_panels(
    species_name = sp,
    merged_df = long_merge_sp,
    lag_df = long_pear %>% filter(species == sp),
    predictors = long_drivers,
    response_label = "residual chronology",
    title_prefix = paste0(sp, " vs Long-Term Climate"),
    plot_subdir = sp_plot_dir,
    stats_subdir = sp_stats_dir,
    file_prefix = paste0(sp, "_long")
  )
}

invisible(plot_species_heatmap(
  short_lag0_pear,
  ww_title("Species vs MODIS-Era Climate", "Correlation Matrix, Lag 0 (Pearson)"),
  "Standard chronologies against the short-period climate variables used in the multispecies workflow.",
  file.path(plot_dir, "01_species_short_climate_matrix_lag0_pearson.png")
))
invisible(plot_species_heatmap(
  short_lag0_spear,
  ww_title("Species vs MODIS-Era Climate", "Correlation Matrix, Lag 0 (Spearman)"),
  "Standard chronologies against the short-period climate variables used in the multispecies workflow.",
  file.path(plot_dir, "02_species_short_climate_matrix_lag0_spearman.png")
))
invisible(plot_species_heatmap(
  short_best_pear,
  ww_title("Species vs MODIS-Era Climate", "Best-Lag Response Matrix (Pearson)"),
  "Each tile uses the strongest absolute correlation across lags 0 to 2 years. Tile labels show the winning lag (L0-L2); stars show significance.",
  file.path(plot_dir, "03_species_short_climate_matrix_bestlag_pearson.png")
))
invisible(plot_species_heatmap(
  short_best_spear,
  ww_title("Species vs MODIS-Era Climate", "Best-Lag Response Matrix (Spearman)"),
  "Each tile uses the strongest absolute rank correlation across lags 0 to 2 years. Tile labels show the winning lag (L0-L2); stars show significance.",
  file.path(plot_dir, "04_species_short_climate_matrix_bestlag_spearman.png")
))
invisible(plot_species_heatmap(
  long_lag0_pear,
  ww_title("Species vs Long-Term Climate", "Correlation Matrix, Lag 0 (Pearson)"),
  "Residual chronologies against PRISM, SPEI, PDSI, ENSO, and the long-period climate composite.",
  file.path(plot_dir, "05_species_long_climate_matrix_lag0_pearson.png")
))
invisible(plot_species_heatmap(
  long_lag0_spear,
  ww_title("Species vs Long-Term Climate", "Correlation Matrix, Lag 0 (Spearman)"),
  "Residual chronologies against PRISM, SPEI, PDSI, ENSO, and the long-period climate composite.",
  file.path(plot_dir, "06_species_long_climate_matrix_lag0_spearman.png")
))
invisible(plot_species_heatmap(
  long_best_pear,
  ww_title("Species vs Long-Term Climate", "Best-Lag Response Matrix (Pearson)"),
  "Each tile uses the strongest absolute correlation across lags 0 to 2 years. Tile labels show the winning lag (L0-L2); stars show significance.",
  file.path(plot_dir, "07_species_long_climate_matrix_bestlag_pearson.png")
))
invisible(plot_species_heatmap(
  long_best_spear,
  ww_title("Species vs Long-Term Climate", "Best-Lag Response Matrix (Spearman)"),
  "Each tile uses the strongest absolute rank correlation across lags 0 to 2 years. Tile labels show the winning lag (L0-L2); stars show significance.",
  file.path(plot_dir, "08_species_long_climate_matrix_bestlag_spearman.png")
))

summary_lines <- c(
  "Species climate analysis summary",
  paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Short climate source: ", normalizePath(short_climate_path, winslash = "/", mustWork = FALSE)),
  paste0("Long climate source: ", normalizePath(long_climate_path, winslash = "/", mustWork = FALSE)),
  paste0("Marker years source: ", if (is.na(marker_years_path)) "not found" else normalizePath(marker_years_path, winslash = "/", mustWork = FALSE)),
  paste0("Species standard chronology source: ", normalizePath(species_standard_path, winslash = "/", mustWork = FALSE)),
  paste0("Species residual chronology source: ", normalizePath(species_residual_path, winslash = "/", mustWork = FALSE)),
  paste0("Species analyzed: ", paste(species_cols, collapse = ", ")),
  paste0("Short drivers: ", paste(short_drivers, collapse = ", ")),
  paste0("Long drivers: ", paste(long_drivers, collapse = ", ")),
  paste0("Maximum lag tested: ", max_lag, " years"),
  paste0("Top-ranked driver table: ", normalizePath(file.path(output_dir, "species_top_climate_responses.csv"), winslash = "/", mustWork = FALSE)),
  paste0("Short climate self-correlation table: ", normalizePath(file.path(output_dir, "climate_short_self_matrix_pearson.csv"), winslash = "/", mustWork = FALSE)),
  paste0("Long climate self-correlation table: ", normalizePath(file.path(output_dir, "climate_long_self_matrix_pearson.csv"), winslash = "/", mustWork = FALSE)),
  paste0("Species-to-species standard table: ", normalizePath(file.path(output_dir, "species_species_standard_matrix_pearson.csv"), winslash = "/", mustWork = FALSE)),
  paste0("Species-to-species residual table: ", normalizePath(file.path(output_dir, "species_species_residual_matrix_pearson.csv"), winslash = "/", mustWork = FALSE)),
  paste0("Marker-year growth summary: ", normalizePath(file.path(output_dir, "species_marker_event_growth_summary.csv"), winslash = "/", mustWork = FALSE)),
  paste0("Top long-driver choices: ", normalizePath(file.path(output_dir, "species_long_top_driver_choices.csv"), winslash = "/", mustWork = FALSE)),
  paste0("Per-species regression plots: ", normalizePath(species_plot_root, winslash = "/", mustWork = FALSE)),
  paste0("Species marker-year time series (full range): ", normalizePath(file.path(plot_dir, "09_species_standard_chronologies_marker_bands_fullrange.png"), winslash = "/", mustWork = FALSE)),
  paste0("Species marker-year time series (recent window): ", normalizePath(file.path(plot_dir, "10_species_standard_chronologies_marker_bands_recent_window.png"), winslash = "/", mustWork = FALSE)),
  paste0("Species flood/drought boxplots: ", normalizePath(file.path(plot_dir, "11_species_marker_event_boxplots.png"), winslash = "/", mustWork = FALSE)),
  paste0("Species top-driver overlays (full range): ", normalizePath(file.path(plot_dir, "12_species_long_top_driver_overlays_fullrange.png"), winslash = "/", mustWork = FALSE)),
  paste0("Species top-driver overlays (recent window): ", normalizePath(file.path(plot_dir, "13_species_long_top_driver_overlays_recent_window.png"), winslash = "/", mustWork = FALSE)),
  paste0("Species response ranking plot: ", normalizePath(file.path(plot_dir, "14_species_long_response_rankings_top3.png"), winslash = "/", mustWork = FALSE)),
  paste0("Species marker-year heatmap: ", normalizePath(file.path(plot_dir, "15_species_marker_year_growth_heatmap.png"), winslash = "/", mustWork = FALSE)),
  paste0("Short climate self-correlation heatmap: ", normalizePath(file.path(plot_dir, "16_climate_short_self_matrix_pearson.png"), winslash = "/", mustWork = FALSE)),
  paste0("Long climate self-correlation heatmap: ", normalizePath(file.path(plot_dir, "18_climate_long_self_matrix_pearson.png"), winslash = "/", mustWork = FALSE)),
  paste0("Species-to-species standard heatmap: ", normalizePath(file.path(plot_dir, "20_species_species_standard_matrix_pearson.png"), winslash = "/", mustWork = FALSE)),
  paste0("Species-to-species residual heatmap: ", normalizePath(file.path(plot_dir, "22_species_species_residual_matrix_pearson.png"), winslash = "/", mustWork = FALSE)),
  "",
  "Top Pearson matches by window:"
)

top_short <- best_summary %>%
  filter(method == "pearson", window == "short") %>%
  transmute(line = paste0("- ", species, ": ", driver, " (lag ", lag, ", r = ", round(r, 3), ", p = ", fmt_p(p), ", n = ", n, ")"))
top_long <- best_summary %>%
  filter(method == "pearson", window == "long") %>%
  transmute(line = paste0("- ", species, ": ", driver, " (lag ", lag, ", r = ", round(r, 3), ", p = ", fmt_p(p), ", n = ", n, ")"))

summary_lines <- c(
  summary_lines,
  if (nrow(top_short) > 0) top_short$line else "- No short-window Pearson matches were available.",
  "",
  "Top long-window Pearson matches:",
  if (nrow(top_long) > 0) top_long$line else "- No long-window Pearson matches were available."
)

writeLines(summary_lines, con = file.path(output_dir, "species_climate_analysis_summary.txt"))

cat("Species climate analysis complete.\n")
cat("Outputs written to:", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n")
