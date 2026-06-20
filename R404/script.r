# clear environment
rm(list = ls())
if (!is.null(dev.list())) dev.off(dev.list()["RStudioGD"])

# set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load packages
library(tidyverse)
library(lubridate)
library(Cairo)
library(RColorBrewer)
library(locfit)
library(quantreg)

# read data and compute daily percent change as proportions (0.01 = 1%)
df <- read_csv("mbi10_values.csv") %>%
  mutate(Date = as.Date(Date, "%Y.%m.%d")) %>%
  filter(Date >= as.Date("2009-06-15")) %>%
  arrange(Date) %>%
  mutate(pct = if ("% change" %in% names(.))
    as.numeric(gsub("%", "", `% change`)) / 100
    else (value / lag(value) - 1)) %>%
  filter(!is.na(pct))

# histogram breaks (30 equal and FD) computed on same scale as df$pct
bins <- list(
  bins30 = seq(floor(min(df$pct)*100)/100, ceiling(max(df$pct)*100)/100, length.out = 31),
  fd = {
    nbins <- nclass.FD(df$pct)
    if (nbins < 1) nbins <- 10
    seq(floor(min(df$pct)*100)/100, ceiling(max(df$pct)*100)/100, length.out = nbins + 1)
  }
)

# kernels and bandwidths (for density overlays)
ks   <- c("rectangular","epanechnikov","biweight","triangular","gaussian")
bws  <- c("nrd0","SJ","ucv") 

# output folder
out_dir <- "hist_outputs"
dir.create(out_dir, FALSE, TRUE)

# safe bandwidth label
safe_bw_label <- function(bw) {
  if (is.character(bw)) paste0("bw_", bw)
  else paste0("bw_", gsub("\\.", "p", as.character(bw)))
}

# density helper with error handling
safe_density <- function(x, kernel, bw, ...) {
  bw_arg <- if (is.character(bw)) bw else as.numeric(bw)
  density(x, kernel = kernel, bw = bw_arg, ...)
}
density_df <- function(x, kernel, bw, ...) {
  d <- tryCatch(safe_density(x, kernel, bw, ...), error = function(e) NULL)
  if (is.null(d)) NULL else data.frame(x = d$x, y = d$y)
}

# x-axis limits (proportions)
x_lim <- c(-0.05, 0.05)
x_lim_month <- c(-0.3, 0.3)

# robust plot_and_save using exact breaks
plot_and_save <- function(data,
                          breaks,
                          dens_df = NULL,
                          title = "",
                          file,
                          w = 8,
                          h = 5,
                          hist_fill   = "grey70",
                          hist_colour = "grey30",
                          dens_colour = "grey10",
                          show_mean   = TRUE,
                          show_median = TRUE,
                          show_quant  = c(0.025, 0.975),
                          quant_col   = c(red = "red", green = "green")) {
  
  # -----------------------------------------------------------------
  # 1. Build the histogram (counts → relative %)
  # -----------------------------------------------------------------
  hst <- hist(data$pct,
              breaks = breaks,
              plot = FALSE,
              include.lowest = TRUE,
              right = FALSE)
  
  hist_df <- data.frame(
    left  = head(hst$breaks, -1),
    right = tail(hst$breaks, -1),
    rel   = (hst$counts / sum(hst$counts)) * 100
  ) %>% mutate(xmin = left, xmax = right)
  
  ymax_hist <- max(hist_df$rel, na.rm = TRUE) * 1.05
  
  # -----------------------------------------------------------------
  # 3. Base plot – histogram bars
  # -----------------------------------------------------------------
  p <- ggplot() +
    geom_rect(data = hist_df,
              aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = rel),
              fill   = hist_fill,
              colour = hist_colour) +
    labs(title = title,
         x = "Period return",
         y = "Relative Frequency (%)") +
    theme_minimal() +
    coord_cartesian(xlim = x_lim) +
    scale_y_continuous(
      limits = c(0, ymax_hist),
      expand = expansion(mult = c(0, 0.02)),
      # secondary axis for density (will be scaled later)
      sec.axis = dup_axis(name = "Density")
    )
  
  # -----------------------------------------------------------------
  # 4. Add kernel‑density overlay (scaled to secondary axis)
  # -----------------------------------------------------------------
  if (!is.null(dens_df)) {
    # Scale density to the secondary axis.
    # We map the density values onto the same range as the histogram
    # by using `scale_y_continuous`'s `sec.axis`.  The easiest way is
    # to transform the density values to the primary y‑scale.
    # Compute a simple linear scaling factor:
    dens_max   <- max(dens_df$y, na.rm = TRUE)
    scale_fact <- ymax_hist / dens_max   # histogram max / density max
    
    p <- p + geom_line(
      data = dens_df,
      aes(x = x, y = y * scale_fact),   # <-- rescale here
      colour = dens_colour,
      linewidth = 0.9
    )
  }
  
  # -----------------------------------------------------------------
  # 5. Optional reference lines (mean, median, quantiles)
  # -----------------------------------------------------------------
  if (show_mean) {
    p <- p + geom_vline(xintercept = mean(data$pct),
                        linetype = "dashed",
                        colour = "black")
  }
  if (show_median) {
    p <- p + geom_vline(xintercept = median(data$pct),
                        linetype = "dotted",
                        colour = "black")
  }
  if (!is.null(show_quant)) {
    probs <- quantile(data$pct, probs = show_quant, na.rm = TRUE)
    p <- p +
      geom_vline(xintercept = probs[1],
                 linetype = "longdash",
                 colour = quant_col[1],
                 linewidth = 0.7) +
      geom_vline(xintercept = probs[2],
                 linetype = "longdash",
                 colour = quant_col[2],
                 linewidth = 0.7)
  }
  
  # -----------------------------------------------------------------
  # 6. Save the plot
  # -----------------------------------------------------------------
  ggsave(filename = file, plot = p, width = w, height = h, device = cairo_pdf)
  invisible(p)
}


# main routine
overall_histograms_daily <- function(df, bins, ks, bws, out) {
  bw_dir   <- file.path(out, "bw_overlays")
  kern_dir <- file.path(out, "kernel_overlays")
  dir.create(bw_dir,   FALSE, TRUE)
  dir.create(kern_dir, FALSE, TRUE)
  
  for (bn in names(bins)) {
    brk <- bins[[bn]]
    
    # compute nbins for title when bn == "fd"
    nbins_for_title <- if (bn == "fd") {
      # nbins is breaks - 1
      length(brk) - 1
    } else NA_integer_
    
    # plain histograms
    for (k in ks) for (b in bws) {
      dens <- density_df(df$pct, k, b)
      ttl  <- if (is.na(nbins_for_title)) {
        sprintf("Hist: %s | %s | %s", bn, k, b)
      } else {
        sprintf("Hist: %s | %s | %s | nbins=%d", bn, k, b, nbins_for_title)
      }
      f <- file.path(out, sprintf("hist_%s_%s_%s.pdf", bn, k, safe_bw_label(b)))
      plot_and_save(df, brk, dens, ttl, f)
    }
    
    # bandwidth overlays (all kernels)
    for (b in bws) {
      dlist <- compact(lapply(ks, function(k) {
        d <- density_df(df$pct, k, b)
        if (!is.null(d)) d$kernel <- k
        d
      }))
      if (length(dlist) == 0) next
      ttl_bw <- if (is.na(nbins_for_title)) {
        sprintf("BW %s | %s", safe_bw_label(b), bn)
      } else {
        sprintf("BW %s | %s | nbins=%d", safe_bw_label(b), bn, nbins_for_title)
      }
      p <- ggplot(bind_rows(dlist),
                  aes(x = x, y = y, colour = kernel, linetype = kernel)) +
        geom_line(linewidth = 0.9) +
        labs(title = ttl_bw, x = "Period return", y = "Density") +
        theme_minimal() + theme(legend.position = "right")
      ggsave(
        filename = file.path(bw_dir,
                             sprintf("bw_%s_%s_overlay_kernels.pdf", safe_bw_label(b), bn)),
        plot = p,
        width = as.numeric(8),
        height = as.numeric(5),
        device = cairo_pdf
      )
    }
    
    # kernel overlays (all bandwidths)
    for (k in ks) {
      dlist <- compact(lapply(bws, function(b) {
        d <- density_df(df$pct, k, b)
        if (!is.null(d)) d$bw <- safe_bw_label(b)
        d
      }))
      if (length(dlist) == 0) next
      ttl_kern <- if (is.na(nbins_for_title)) {
        sprintf("Kernel %s | %s", k, bn)
      } else {
        sprintf("Kernel %s | %s | nbins=%d", k, bn, nbins_for_title)
      }
      p <- ggplot(bind_rows(dlist),
                  aes(x = x, y = y, colour = bw, linetype = bw)) +
        geom_line(linewidth = 0.9) +
        labs(title = ttl_kern, x = "Daily return (proportion)", y = "Density") +
        theme_minimal() + theme(legend.position = "right")
      ggsave(
        filename = file.path(kern_dir,
                             sprintf("kernel_%s_%s_overlay_bws.pdf", k, bn)),
        plot = p,
        width = as.numeric(8),
        height = as.numeric(5),
        device = cairo_pdf
      )
    }
  }
}
# run
overall_histograms_daily(df, bins, ks, bws, out_dir)










# safe compact (uses purrr if available, fallback)
compact <- function(x) Filter(Negate(is.null), x)

# compute density safely and return data frame on a common x grid
compute_density_df <- function(xvec, kernel, bw, x_grid = NULL) {
  d <- tryCatch({
    bw_arg <- if (is.character(bw)) bw else as.numeric(bw)
    dens <- density(xvec, kernel = kernel, bw = bw_arg)
    data.frame(x = dens$x, y = dens$y)
  }, error = function(e) NULL)
  if (!is.null(d) && !is.null(x_grid)) {
    # interpolate onto x_grid to permit averaging across groups
    y_int <- approx(d$x, d$y, xout = x_grid, rule = 2)$y
    data.frame(x = x_grid, y = y_int)
  } else d
}

# compute densities per weekday (returns data.frame with weekday, x, y)
compute_per_wd <- function(df_wd, kernel, bw, x_grid = NULL) {
  wd_levels <- unique(df_wd$weekday)
  res <- lapply(wd_levels, function(wd) {
    vals <- df_wd %>% filter(weekday == wd) %>% pull(pct)
    if (length(vals) < 2) return(NULL)
    d <- compute_density_df(vals, kernel, bw, x_grid)
    if (is.null(d)) return(NULL)
    d$weekday <- wd
    d
  }) %>% compact() %>% bind_rows()
  if (nrow(res) > 0) res$weekday <- factor(res$weekday, levels = wd_levels)
  res
}

# Ensure breaks span data range (inclusive)
expand_breaks_to_data <- function(breaks, x) {
  rmin <- min(x, na.rm = TRUE); rmax <- max(x, na.rm = TRUE)
  bmin <- min(breaks); bmax <- max(breaks)
  if (rmin < bmin) breaks[1] <- floor(rmin * 100) / 100
  if (rmax > bmax) breaks[length(breaks)] <- ceiling(rmax * 100) / 100
  sort(unique(breaks))
}

# plot_and_save_hist: same semantics as plot_and_save but accepts a density dataframe
plot_and_save_hist <- function(data,
                          breaks,
                          dens_df = NULL,
                          title = "",
                          file,
                          w = 8,
                          h = 5,
                          hist_fill   = "grey70",
                          hist_colour = "grey30",
                          dens_colour = "grey10",
                          show_mean   = TRUE,
                          show_median = TRUE,
                          show_quant  = c(0.025, 0.975),
                          quant_col   = c(red = "red", green = "green")) {
  
  # ---------------------------------------------------------------
  # 1. Build the histogram (relative %)
  # ---------------------------------------------------------------
  hst <- hist(data$pct,
              breaks = breaks,
              plot = FALSE,
              include.lowest = TRUE,
              right = FALSE)
  
  hist_df <- data.frame(
    left  = head(hst$breaks, -1),
    right = tail(hst$breaks, -1),
    rel   = (hst$counts / sum(hst$counts)) * 100
  ) %>% mutate(xmin = left, xmax = right)
  
  ymax_hist <- max(hist_df$rel, na.rm = TRUE) * 1.05
  
  # ---------------------------------------------------------------
  # 2. Base plot – histogram bars
  # ---------------------------------------------------------------
  p <- ggplot() +
    geom_rect(data = hist_df,
              aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = rel),
              fill   = hist_fill,
              colour = hist_colour) +
    labs(title = title,
         x = "Period return",
         y = "Relative Frequency (%)") +
    theme_minimal() +
    coord_cartesian(xlim = x_lim) +
    scale_y_continuous(
      limits = c(0, ymax_hist),
      expand = expansion(mult = c(0, 0.02)),
      sec.axis = dup_axis(name = "Density")
    )
  
  # ---------------------------------------------------------------
  # 3. Add kernel‑density overlay (scaled to secondary axis)
  # ---------------------------------------------------------------
  if (!is.null(dens_df)) {
    dens_max   <- max(dens_df$y, na.rm = TRUE)
    scale_fact <- ymax_hist / dens_max   # histogram max / density max
    
    p <- p + geom_line(
      data = dens_df,
      aes(x = x, y = y * scale_fact),
      colour = dens_colour,
      linewidth = 0.9
    )
  }
  
  # ---------------------------------------------------------------
  # 4. Optional reference lines (mean, median, quantiles)
  # ---------------------------------------------------------------
  if (show_mean) {
    p <- p + geom_vline(xintercept = mean(data$pct),
                        linetype = "dashed",
                        colour = "black")
  }
  if (show_median) {
    p <- p + geom_vline(xintercept = median(data$pct),
                        linetype = "dotted",
                        colour = "black")
  }
  if (!is.null(show_quant)) {
    probs <- quantile(data$pct, probs = show_quant, na.rm = TRUE)
    p <- p +
      geom_vline(xintercept = probs[1],
                 linetype = "longdash",
                 colour = quant_col[1],
                 linewidth = 0.7) +
      geom_vline(xintercept = probs[2],
                 linetype = "longdash",
                 colour = quant_col[2],
                 linewidth = 0.7)
  }
  
  # ---------------------------------------------------------------
  # 5. Save the plot
  # ---------------------------------------------------------------
  ggsave(filename = file, plot = p, width = w, height = h, device = cairo_pdf)
  invisible(p)
}



# -----------------------------
# Prepare weekday data & bins
# -----------------------------
df_wd <- df %>%
  mutate(weekday = wday(Date, label = TRUE, week_start = 1)) %>%
  filter(weekday %in% c("Mon","Tue","Wed","Thu","Fri")) %>%
  mutate(weekday = recode(as.character(weekday),
                          Mon = "Monday", Tue = "Tuesday",
                          Wed = "Wednesday", Thu = "Thursday", Fri = "Friday"))

if (!exists("x_grid")) x_grid <- seq(min(df_wd$pct, na.rm = TRUE), max(df_wd$pct, na.rm = TRUE), length.out = 512)

bin_breaks_wd <- list(
  bins30 = seq(x_lim[1], x_lim[2], length.out = 31),
  fd = {
    nb <- nclass.FD(df_wd$pct)
    if (nb < 1) nb <- 10
    seq(x_lim[1], x_lim[2], length.out = nb + 1)
  }
)

kernels <- if (exists("ks")) ks else c("rectangular","epanechnikov","biweight","triangular","gaussian")
bw_list  <- if (exists("bws")) bws else c("nrd0","SJ","ucv")

# New output folder for weekday-specific plots (outside hist_outputs)
out_dir_wd_root <- "hist_outputs_weekdays"
dir.create(out_dir_wd_root, FALSE, TRUE)

# subfolders for overlays
out_dir_wd_bw    <- file.path(out_dir_wd_root, "bw_overlays")
out_dir_wd_kern  <- file.path(out_dir_wd_root, "kernel_overlays")
dir.create(out_dir_wd_bw,    FALSE, TRUE)
dir.create(out_dir_wd_kern,  FALSE, TRUE)

# folder for individual weekday histogram PDFs
out_dir_wd <- out_dir_wd_root

# -----------------------------
# Weekday-only plotting routine
# -----------------------------
weekday_histograms <- function(df_wd, bin_breaks_wd, kernels, bw_list, out_dir_wd, out_dir_bw, out_dir_kern, x_grid) {
  for (bb in names(bin_breaks_wd)) {
    for (k in kernels) {
      for (bw in bw_list) {
        # ensure breaks span full data before use
        brk_wd <- expand_breaks_to_data(bin_breaks_wd[[bb]], df_wd$pct)
        
        # densities per weekday on common grid
        dens_wd <- compute_per_wd(df_wd, k, bw, x_grid)
        
        # 1) individual weekday histograms with density overlay
        if (!is.null(dens_wd) && nrow(dens_wd) > 0) {
          for (day in sort(unique(dens_wd$weekday))) {
            df_day   <- df_wd %>% filter(weekday == day)
            dens_day <- dens_wd %>% filter(weekday == day) %>% select(x,y)
            ttl <- paste(day, "|", bb, "| Kernel:", k, "| bw:", safe_bw_label(bw))
            fname <- file.path(out_dir_wd, paste0("weekday_", day, "_", bb,
                                                  "_kernel_", k, "_", safe_bw_label(bw), ".pdf"))
            plot_and_save_hist(df_day, brk_wd, dens_day, ttl, fname, w = 7, h = 5)
          }
        }
        
        # 2) overlay of kernels (avg over weekdays) for this bandwidth -> save under kernel_overlays
        dlist_kernels <- lapply(kernels, function(kern_i) {
          dd <- compute_per_wd(df_wd, kern_i, bw, x_grid)
          if (!is.null(dd) && nrow(dd) > 0) {
            dd_summary <- dd %>% group_by(x) %>% summarise(y = mean(y, na.rm = TRUE), .groups = "drop")
            dd_summary$kernel <- kern_i
            dd_summary
          } else NULL
        }) %>% compact() %>% bind_rows()
        
        if (nrow(dlist_kernels) > 0) {
          p_all_kerns <- ggplot(dlist_kernels, aes(x = x, y = y, colour = kernel, linetype = kernel)) +
            geom_line(size = 0.9) +
            labs(title = paste("Overlay kernels (avg over weekdays) |", bb, "| bw:", safe_bw_label(bw)),
                 x = "Daily return (proportion)", y = "Density") +
            theme_minimal() + theme(legend.position = "right")
          ggsave(file.path(out_dir_kern, paste0("overlay_kernels_avg_wd_", bb, "_bw_", safe_bw_label(bw), ".pdf")),
                 p_all_kerns, width = 8, height = 5, device = cairo_pdf)
        }
        
        # 3) overlay of bandwidths averaged across weekdays (for this kernel k) -> save under bw_overlays
        dlist_bws <- lapply(bw_list, function(bw_i) {
          dd <- compute_per_wd(df_wd, k, bw_i, x_grid)
          if (!is.null(dd) && nrow(dd) > 0) {
            dd_summary <- dd %>% group_by(x) %>% summarise(y = mean(y, na.rm = TRUE), .groups = "drop")
            dd_summary$bw <- safe_bw_label(bw_i)
            dd_summary
          } else NULL
        }) %>% compact() %>% bind_rows()
        
        if (nrow(dlist_bws) > 0) {
          p_bw_overlay <- ggplot(dlist_bws, aes(x = x, y = y, colour = bw, linetype = bw)) +
            geom_line(size = 0.9) +
            labs(title = paste("Overlay bandwidths (avg over weekdays) |", bb, "| Kernel:", k),
                 x = "Daily return (proportion)", y = "Density") +
            theme_minimal() + theme(legend.position = "right")
          ggsave(file.path(out_dir_bw, paste0("overlay_bws_avg_wd_", bb, "_kernel_", k, ".pdf")),
                 p_bw_overlay, width = 8, height = 5, device = cairo_pdf)
        }
      } # end bw loop
    } # end kernel loop
  } # end bins loop
}

# Run weekday routine (writes into hist_outputs_weekdays/*)
weekday_histograms(df_wd = df_wd,
                   bin_breaks_wd = bin_breaks_wd,
                   kernels = kernels,
                   bw_list = bw_list,
                   out_dir_wd = out_dir_wd,
                   out_dir_bw = out_dir_wd_bw,
                   out_dir_kern = out_dir_wd_kern,
                   x_grid = x_grid)

# Parameters: set preferred week start / end names (use full weekday names)
preferred_start_name <- "Monday"   # preferred start-of-week (will fall back to first available day)
preferred_end_name   <- "Friday"   # preferred end-of-week (will fall back to last available day)

# prepare df_wd
df_wd <- df_wd %>%
  mutate(Date = as.Date(Date),
         Value = as.numeric(Value),
         pct = if ("% change" %in% names(.))
           as.numeric(gsub("%", "", `% change`)) / 100
         else if ("pct" %in% names(.))
           as.numeric(pct)
         else
           (Value / lag(Value) - 1)) %>%
  arrange(Date)

# week id by ISO week
week_id <- function(date) paste0(isoyear(date), "-W", sprintf("%02d", isoweek(date)))

weekly_progress <- df_wd %>%
  arrange(Date) %>%
  mutate(wid = week_id(Date),
         wday_name = weekdays(Date)) %>%
  group_by(wid) %>%
  group_modify(~{
    dat <- .x
    # preferred start: choose row with preferred_start_name if present, else first row
    if (any(dat$wday_name == preferred_start_name)) {
      start_row <- dat %>% filter(wday_name == preferred_start_name) %>% slice_min(Date, n = 1)
    } else start_row <- dat %>% slice_min(Date, n = 1)
    # preferred end: choose row with preferred_end_name if present, else last row
    if (any(dat$wday_name == preferred_end_name)) {
      end_row <- dat %>% filter(wday_name == preferred_end_name) %>% slice_max(Date, n = 1)
    } else end_row <- dat %>% slice_max(Date, n = 1)
    tibble(
      week_start = start_row$Date[1],
      week_start_value = start_row$Value[1],
      week_end = end_row$Date[1],
      week_end_value = end_row$Value[1]
    )
  }, .keep = FALSE) %>%
  ungroup() %>%
  filter(!is.na(week_start_value), !is.na(week_end_value), week_end >= week_start) %>%
  mutate(weekly_change = (week_end_value / week_start_value) - 1) %>%
  arrange(week_end)

weekly_vals <- weekly_progress$weekly_change

# build bin breaks and run plotting loop (reuse earlier code)
bin_breaks_list <- list(
  bins30 = seq(floor(min(weekly_vals, na.rm = TRUE)*100)/100,
               ceiling(max(weekly_vals, na.rm = TRUE)*100)/100,
               length.out = 31),
  fd = {
    nbins <- nclass.FD(weekly_vals)
    if (is.na(nbins) || nbins < 1) nbins <- 10
    seq(floor(min(weekly_vals, na.rm = TRUE)*100)/100,
        ceiling(max(weekly_vals, na.rm = TRUE)*100)/100,
        length.out = nbins + 1)
  }
)

out_dir_weekly <- "hist_outputs_weekly"
dir.create(out_dir_weekly, FALSE, TRUE)
bw_dir_wk   <- file.path(out_dir_weekly, "bw_overlays"); dir.create(bw_dir_wk, FALSE, TRUE)
kern_dir_wk <- file.path(out_dir_weekly, "kernel_overlays"); dir.create(kern_dir_wk, FALSE, TRUE)

compute_density_df <- function(x, kernel, bw, ...) {
  d <- tryCatch(safe_density(x, kernel, bw, ...), error = function(e) NULL)
  if (is.null(d)) NULL else data.frame(x = d$x, y = d$y)
}

kernels <- if (exists("ks")) ks else c("rectangular","epanechnikov","biweight","triangular","gaussian")
bw_list  <- if (exists("bws")) bws else c("nrd0","SJ","ucv")

for (bb in names(bin_breaks_list)) {
  brk_week <- bin_breaks_list[[bb]]
  nbins_for_title <- if (bb == "fd") length(brk_week) - 1 else NA_integer_
  
  for (k in kernels) for (bw in bw_list) {
    dens_wk <- compute_density_df(weekly_vals, k, bw)
    ttl <- if (is.na(nbins_for_title)) sprintf("Weekly: %s | %s | %s", bb, k, bw) else sprintf("Weekly: %s | %s | %s | nbins=%d", bb, k, bw, nbins_for_title)
    fname <- file.path(out_dir_weekly, sprintf("weekly_hist_%s_%s_%s.pdf", bb, k, safe_bw_label(bw)))
    plot_and_save(data.frame(pct = weekly_vals), brk_week, dens_wk, ttl, fname)
  }
  
  for (bw in bw_list) {
    dlist <- compact(lapply(kernels, function(k) {
      d <- compute_density_df(weekly_vals, k, bw); if (!is.null(d)) d$kernel <- k; d
    }))
    if (length(dlist) == 0) next
    p <- ggplot(bind_rows(dlist), aes(x = x, y = y, colour = kernel, linetype = kernel)) + geom_line(linewidth = 0.9) + labs(title = sprintf("BW %s | %s", safe_bw_label(bw), bb), x = "Weekly change (proportion)", y = "Density") + theme_minimal() + theme(legend.position = "right")
    ggsave(filename = file.path(bw_dir_wk, sprintf("bw_%s_%s_overlay_kernels.pdf", safe_bw_label(bw), bb)), plot = p, width = 8, height = 5, device = cairo_pdf)
  }
  
  for (k in kernels) {
    dlist <- compact(lapply(bw_list, function(bw) {
      d <- compute_density_df(weekly_vals, k, bw); if (!is.null(d)) d$bw <- safe_bw_label(bw); d
    }))
    if (length(dlist) == 0) next
    p <- ggplot(bind_rows(dlist), aes(x = x, y = y, colour = bw, linetype = bw)) + geom_line(linewidth = 0.9) + labs(title = sprintf("Kernel %s | %s", k, bb), x = "Weekly change (proportion)", y = "Density") + theme_minimal() + theme(legend.position = "right")
    ggsave(filename = file.path(kern_dir_wk, sprintf("kernel_%s_%s_overlay_bws.pdf", k, bb)), plot = p, width = 8, height = 5, device = cairo_pdf)
  }
}





df_wd <- df_wd %>% mutate(Date = as.Date(Date), Value = as.numeric(Value)) %>% arrange(Date)

month_id <- function(date) paste0(year(date), "-", sprintf("%02d", month(date)))

monthly_progress <- df_wd %>%
  arrange(Date) %>%
  mutate(mid = month_id(Date),
         mday_name = weekdays(Date)) %>%
  group_by(mid) %>%
  group_modify(~{
    dat <- .x
    start_row <- dat %>% slice_min(Date, n = 1)
    end_row   <- dat %>% slice_max(Date, n = 1)
    tibble(
      month_start = start_row$Date[1],
      month_start_value = start_row$Value[1],
      month_end = end_row$Date[1],
      month_end_value = end_row$Value[1]
    )
  }, .keep = FALSE) %>%
  ungroup() %>%
  filter(!is.na(month_start_value), !is.na(month_end_value), month_end >= month_start) %>%
  mutate(monthly_change = (month_end_value / month_start_value) - 1) %>%
  arrange(month_end)

monthly_vals <- monthly_progress$monthly_change

# bin breaks (same style as daily/weekly)
bin_breaks_list_month <- list(
  bins30 = seq(floor(min(monthly_vals, na.rm = TRUE)*100)/100,
               ceiling(max(monthly_vals, na.rm = TRUE)*100)/100,
               length.out = 31),
  fd = {
    nbins <- nclass.FD(monthly_vals)
    if (is.na(nbins) || nbins < 1) nbins <- 10
    seq(floor(min(monthly_vals, na.rm = TRUE)*100)/100,
        ceiling(max(monthly_vals, na.rm = TRUE)*100)/100,
        length.out = nbins + 1)
  }
)

# output dirs
out_dir_monthly <- "hist_outputs_monthly"
dir.create(out_dir_monthly, FALSE, TRUE)
bw_dir_mo   <- file.path(out_dir_monthly, "bw_overlays"); dir.create(bw_dir_mo, FALSE, TRUE)
kern_dir_mo <- file.path(out_dir_monthly, "kernel_overlays"); dir.create(kern_dir_mo, FALSE, TRUE)

# ensure compute_density_df, kernels, bw_list exist (reuse from earlier)
compute_density_df <- function(x, kernel, bw, ...) {
  d <- tryCatch(safe_density(x, kernel, bw, ...), error = function(e) NULL)
  if (is.null(d)) NULL else data.frame(x = d$x, y = d$y)
}
kernels <- if (exists("ks")) ks else c("rectangular","epanechnikov","biweight","triangular","gaussian")
bw_list  <- if (exists("bws")) bws else c("nrd0","SJ","ucv")

# monthly x-axis limits
x_lim_month <- c(-0.3, 0.3)

# plotting loop (mirrors weekly section) with monthly x-limits applied
for (bb in names(bin_breaks_list_month)) {
  brk_month <- bin_breaks_list_month[[bb]]
  nbins_for_title <- if (bb == "fd") length(brk_month) - 1 else NA_integer_
  
  # 9.1 individual histograms (temporarily override global x_lim)
  for (k in kernels) for (bw in bw_list) {
    dens_mo <- compute_density_df(monthly_vals, k, bw)
    ttl <- if (is.na(nbins_for_title)) {
      sprintf("Monthly: %s | %s | %s", bb, k, bw)
    } else {
      sprintf("Monthly: %s | %s | %s | nbins=%d", bb, k, bw, nbins_for_title)
    }
    fname <- file.path(out_dir_monthly, sprintf("monthly_hist_%s_%s_%s.pdf", bb, k, safe_bw_label(bw)))
    old_x_lim <- x_lim
    x_lim <- x_lim_month
    plot_and_save(data.frame(pct = monthly_vals), brk_month, dens_mo, ttl, fname)
    x_lim <- old_x_lim
  }
  
  # 9.2 bandwidth overlays (all kernels) — add coord_cartesian with x_lim_month
  for (bw in bw_list) {
    dlist <- compact(lapply(kernels, function(k) {
      d <- compute_density_df(monthly_vals, k, bw)
      if (!is.null(d)) d$kernel <- k
      d
    }))
    if (length(dlist) == 0) next
    ttl_bw <- if (is.na(nbins_for_title)) sprintf("BW %s | %s", safe_bw_label(bw), bb) else sprintf("BW %s | %s | nbins=%d", safe_bw_label(bw), bb, nbins_for_title)
    p <- ggplot(bind_rows(dlist), aes(x = x, y = y, colour = kernel, linetype = kernel)) +
      geom_line(linewidth = 0.9) +
      labs(title = ttl_bw, x = "Monthly change (proportion)", y = "Density") +
      theme_minimal() + theme(legend.position = "right") +
      coord_cartesian(xlim = x_lim_month)
    ggsave(filename = file.path(bw_dir_mo, sprintf("bw_%s_%s_overlay_kernels.pdf", safe_bw_label(bw), bb)),
           plot = p, width = 8, height = 5, device = cairo_pdf)
  }
  
  # 9.3 kernel overlays (all bandwidths) — add coord_cartesian with x_lim_month
  for (k in kernels) {
    dlist <- compact(lapply(bw_list, function(bw) {
      d <- compute_density_df(monthly_vals, k, bw)
      if (!is.null(d)) d$bw <- safe_bw_label(bw)
      d
    }))
    if (length(dlist) == 0) next
    ttl_kern <- if (is.na(nbins_for_title)) sprintf("Kernel %s | %s", k, bb) else sprintf("Kernel %s | %s | nbins=%d", k, bb, nbins_for_title)
    p <- ggplot(bind_rows(dlist), aes(x = x, y = y, colour = bw, linetype = bw)) +
      geom_line(linewidth = 0.9) +
      labs(title = ttl_kern, x = "Monthly change (proportion)", y = "Density") +
      theme_minimal() + theme(legend.position = "right") +
      coord_cartesian(xlim = x_lim_month)
    ggsave(filename = file.path(kern_dir_mo, sprintf("kernel_%s_%s_overlay_bws.pdf", k, bb)),
           plot = p, width = 8, height = 5, device = cairo_pdf)
  }
}


# Function to create year identifier
year_id <- function(date) {
  year(date)
}

# Yearly Progress Calculation
yearly_progress <- df_wd %>%
  arrange(Date) %>%
  mutate(year = year_id(Date)) %>%
  group_by(year) %>%
  group_modify(~{
    dat <- .x
    start_row <- dat %>% slice_min(Date, n = 1)
    end_row   <- dat %>% slice_max(Date, n = 1)
    tibble(
      year_start = start_row$Date[1],
      year_start_value = start_row$Value[1],
      year_end = end_row$Date[1],
      year_end_value = end_row$Value[1]
    )
  }, .keep = FALSE) %>%
  ungroup() %>%
  filter(!is.na(year_start_value), !is.na(year_end_value), year_end >= year_start) %>%
  mutate(yearly_change = (year_end_value / year_start_value) - 1) %>%
  arrange(year_end)

yearly_vals <- yearly_progress$yearly_change

# Bin breaks for yearly values
bin_breaks_list_year <- list(
  bins10 = seq(floor(min(yearly_vals, na.rm = TRUE) * 100) / 100,
               ceiling(max(yearly_vals, na.rm = TRUE) * 100) / 100,
               length.out = 11),
  fd = {
    nbins <- nclass.FD(yearly_vals)
    if (is.na(nbins) || nbins < 1) nbins <- 10
    seq(floor(min(yearly_vals, na.rm = TRUE) * 100) / 100,
        ceiling(max(yearly_vals, na.rm = TRUE) * 100) / 100,
        length.out = nbins + 1)
  }
)

# Output directories for yearly results
out_dir_yearly <- "hist_outputs_yearly"
dir.create(out_dir_yearly, FALSE, TRUE)
bw_dir_yr <- file.path(out_dir_yearly, "bw_overlays"); dir.create(bw_dir_yr, FALSE, TRUE)
kern_dir_yr <- file.path(out_dir_yearly, "kernel_overlays"); dir.create(kern_dir_yr, FALSE, TRUE)

# Yearly x-axis limits
x_lim_year <- c(-0.5, 0.5)

# Plotting loop for yearly data
for (bb in names(bin_breaks_list_year)) {
  brk_year <- bin_breaks_list_year[[bb]]
  nbins_for_title <- if (bb == "fd") length(brk_year) - 1 else NA_integer_
  
  # 10.1 Individual histograms (temporarily override global x_lim)
  for (k in kernels) {
    for (bw in bw_list) {
      dens_yr <- compute_density_df(yearly_vals, k, bw)
      ttl <- if (is.na(nbins_for_title)) {
        sprintf("Yearly: %s | %s | %s", bb, k, bw)
      } else {
        sprintf("Yearly: %s | %s | %s | nbins=%d", bb, k, bw, nbins_for_title)
      }
      fname <- file.path(out_dir_yearly, sprintf("yearly_hist_%s_%s_%s.pdf", bb, k, safe_bw_label(bw)))
      old_x_lim <- x_lim
      x_lim <- x_lim_year
      plot_and_save(data.frame(pct = yearly_vals), brk_year, dens_yr, ttl, fname)
      x_lim <- old_x_lim
    }
  }
  
  # 10.2 Bandwidth overlays (all kernels)
  for (bw in bw_list) {
    dlist <- compact(lapply(kernels, function(k) {
      d <- compute_density_df(yearly_vals, k, bw)
      if (!is.null(d)) d$kernel <- k
      d
    }))
    if (length(dlist) == 0) next
    ttl_bw <- if (is.na(nbins_for_title)) sprintf("BW %s | %s", safe_bw_label(bw), bb) else sprintf("BW %s | %s | nbins=%d", safe_bw_label(bw), bb, nbins_for_title)
    
    p <- ggplot(bind_rows(dlist), aes(x = x, y = y, colour = kernel, linetype = kernel)) +
      geom_line(linewidth = 0.9) +
      labs(title = ttl_bw, x = "Yearly Change (Proportion)", y = "Density") +
      theme_minimal() + theme(legend.position = "right") +
      coord_cartesian(xlim = x_lim_year)
    
    ggsave(filename = file.path(bw_dir_yr, sprintf("bw_%s_%s_overlay_kernels.pdf", safe_bw_label(bw), bb)),
           plot = p, width = 8, height = 5, device = cairo_pdf)
  }
  
  # 10.3 Kernel overlays (all bandwidths)
  for (k in kernels) {
    dlist <- compact(lapply(bw_list, function(bw) {
      d <- compute_density_df(yearly_vals, k, bw)
      if (!is.null(d)) d$bw <- safe_bw_label(bw)
      d
    }))
    if (length(dlist) == 0) next
    ttl_kern <- if (is.na(nbins_for_title)) sprintf("Kernel %s | %s", k, bb) else sprintf("Kernel %s | %s | nbins=%d", k, bb, nbins_for_title)
    
    p <- ggplot(bind_rows(dlist), aes(x = x, y = y, colour = bw, linetype = bw)) +
      geom_line(linewidth = 0.9) +
      labs(title = ttl_kern, x = "Yearly Change (Proportion)", y = "Density") +
      theme_minimal() + theme(legend.position = "right") +
      coord_cartesian(xlim = x_lim_year)
    
    ggsave(filename = file.path(kern_dir_yr, sprintf("kernel_%s_%s_overlay_bws.pdf", k, bb)),
           plot = p, width = 8, height = 5, device = cairo_pdf)
  }
}






# Define kernel names and bandwidth selection methods
rek_k <- c("rectangular", "epanechnikov", "biweight", "triangular", "gaussian")
bw_rules <- c("nrd0", "SJ", "ucv")

# Function to get bandwidth based on chosen method
get_bw <- function(rule, x) {
  switch(rule,
         nrd0 = stats::bw.nrd0(x),
         SJ   = stats::bw.SJ(x),
         ucv  = stats::bw.ucv(x),
         stop("Unknown rule"))
}

# Safe bandwidth label
safe_bw_label <- function(bw_rule) { as.character(bw_rule) }

# Smooth function using locfit
smooth_df <- function(x, y, kernel_name_unused, bw_numeric) {
  tryCatch({
    x <- as.numeric(x)
    n <- sum(!is.na(x) & !is.na(y))
    if (n < 5) return(NULL)
    rng <- diff(range(x, na.rm = TRUE))
    if (is.na(rng) || rng == 0) return(NULL)
    
    # compute span from bandwidth and clamp to ensure enough neighbors
    span <- bw_numeric / rng
    k_min <- max(10, ceiling(0.02 * n))    # at least 10 neighbors or 2% of points
    span_min <- k_min / n
    span <- min(max(span, span_min), 1)
    
    model <- stats::loess(y ~ x, span = span, family = "gaussian",
                          control = loess.control(surface = "interpolate"))
    newx <- sort(unique(x))
    pred <- predict(model, newdata = data.frame(x = newx))
    data.frame(x = newx, y = pred, stringsAsFactors = FALSE)
  }, error = function(e) {
    message(sprintf("Error in loess smoothing with bw=%s: %s", bw_numeric, e$message))
    NULL
  })
}


# Plotting function
plot_overlay <- function(df, aescol, ttl, xlim, ylim, x_as_date=FALSE, date_origin="1970-01-01") {
  p <- ggplot(df, aes(x = if(x_as_date) as.Date(x, origin=date_origin) else x, 
                      y = y, colour = .data[[aescol]])) +
    geom_line(linewidth = .9) +
    labs(title = ttl, x = if(x_as_date) "Date" else "x", y = "Value") +
    theme_minimal() + theme(legend.position = "right") +
    coord_cartesian(xlim = if(x_as_date) as.Date(xlim, origin=date_origin) else xlim,
                    ylim = ylim)
}

ensure_dir <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

# create top-level period folders and weekday subfolders up-front
base_out <- "local_regressions"
ensure_dir(base_out)

periods <- c("daily","workdays","weekly","monthly","yearly")
for (pdir in periods) ensure_dir(file.path(base_out, pdir))

weekdays_list <- c("Monday","Tuesday","Wednesday","Thursday","Friday")
for (wd in weekdays_list) ensure_dir(file.path(base_out, "workdays", wd))

# create overlay subfolders for every period
for (pdir in periods) {
  ensure_dir(file.path(base_out, pdir, "bw_overlays"))
  ensure_dir(file.path(base_out, pdir, "kernel_overlays"))
}

# List of data frames based on period (ensure these data frames exist before running)
period_map <- list(
  daily = list(df = df, xcol = "Date", ycol = "Value", x_as_date = TRUE),
  workdays = list(df = df_wd, xcol = "Date", ycol = "Value", x_as_date = TRUE),
  weekly = list(df = weekly_progress, xcol = "week_end", ycol = "week_end_value", x_as_date = TRUE),
  monthly = list(df = monthly_progress, xcol = "month_end", ycol = "month_end_value", x_as_date = TRUE),
  yearly = list(df = yearly_progress, xcol = "year_end", ycol = "year_end_value", x_as_date = TRUE)
)


# Main processing loop
for (p in names(period_map)) {
  info <- period_map[[p]]
  dfp <- info$df
  xcol <- info$xcol
  ycol <- info$ycol
  x_as_date <- info$x_as_date
  
  if (p == "workdays") {
    dfp$weekday <- weekdays(as.Date(dfp[[xcol]]))
    weekdays_list <- c("Monday","Tuesday","Wednesday","Thursday","Friday")
    
    for (gname in weekdays_list) {
      sel <- dfp[dfp$weekday == gname, , drop = FALSE]
      if (nrow(sel) < 5) next
      
      x <- as.numeric(sel[[xcol]])
      y <- sel[[ycol]]
      raw_pts <- tibble(x = x, y = y)
      ylim <- range(y, na.rm = TRUE)
      xlim <- range(x, na.rm = TRUE)
      
      numeric_bws <- map_dbl(bw_rules, function(rule) {
        bwv <- tryCatch(get_bw(rule, x), error = function(e) NA_real_)
        if (is.na(bwv) || !is.finite(bwv) || bwv <= .Machine$double.eps * 100) bwv <- tryCatch(get_bw("nrd0", x), error = function(e) NA_real_)
        bwv
      })
      names(numeric_bws) <- bw_rules
      rng <- diff(range(x, na.rm = TRUE))
      if (!is.na(rng) && rng > 0) {
        min_bw <- rng * 1e-3
        max_bw <- rng * 0.5
        numeric_bws <- pmin(pmax(numeric_bws, min_bw), max_bw)
      }
      
      dir.create(file.path(base_out, p, gname), recursive = TRUE, showWarnings = FALSE)
      
      for (kernel_name in rek_k) {
        for (bw_rule in bw_rules) {
          bw_numeric <- numeric_bws[bw_rule]
          if (is.na(bw_numeric) || !is.finite(bw_numeric)) next
          d <- smooth_df(x, y, kernel_name, bw_numeric)
          if (is.null(d)) next
          d$kernel <- kernel_name; d$bw <- safe_bw_label(bw_rule)
          ttl <- sprintf("%s %s local regression: kernel=%s | bw=%s", toupper(p), gname, kernel_name, safe_bw_label(bw_rule))
          pplot <- ggplot() +
            geom_line(data = d, aes(x = if(x_as_date) as.Date(x, origin='1970-01-01') else x, y = y), colour = "steelblue", size = 0.9) +
            geom_point(data = raw_pts, aes(x = if(x_as_date) as.Date(x, origin='1970-01-01') else x, y = y), colour = "black", size = 1.2) +
            labs(title = ttl, x = if(x_as_date) "Date" else "x", y = "Value") +
            theme_minimal() + theme(legend.position = "none") +
            coord_cartesian(xlim = if(x_as_date) as.Date(xlim, origin='1970-01-01') else xlim, ylim = ylim)
          fname <- file.path(base_out, p, gname, sprintf("%s_locreg_%s_%s.pdf", gname, kernel_name, safe_bw_label(bw_rule)))
          ggsave(filename = fname, plot = pplot, width = 8, height = 5, device = cairo_pdf)
        }
      }
      
      for (bw_rule in bw_rules) {
        bw_numeric <- numeric_bws[bw_rule]
        if (is.na(bw_numeric) || !is.finite(bw_numeric)) next
        lst <- compact(lapply(rek_k, function(kernel_name) {
          d <- smooth_df(x, y, kernel_name, bw_numeric)
          if (is.null(d)) return(NULL)
          d$kernel <- kernel_name; d$bw <- safe_bw_label(bw_rule); d
        }))
        if (length(lst) == 0) next
        pplot <- plot_overlay(bind_rows(lst), "kernel",
                              sprintf("BW %s | %s local regressions", safe_bw_label(bw_rule), gname),
                              xlim, ylim, x_as_date = x_as_date)
        dir.create(file.path(base_out, p, gname, "bw_overlays"), recursive = TRUE, showWarnings = FALSE)
        ggsave(filename = file.path(base_out, p, gname, "bw_overlays",
                                    sprintf("bw_%s_%s_overlay_kernels.pdf", safe_bw_label(bw_rule), gname)),
               plot = pplot, width = 8, height = 5, device = cairo_pdf)
        
      }
      
      for (kernel_name in rek_k) {
        lst <- compact(lapply(bw_rules, function(bw_rule) {
          bw_numeric <- numeric_bws[bw_rule]
          if (is.na(bw_numeric) || !is.finite(bw_numeric)) return(NULL)
          d <- smooth_df(x, y, kernel_name, bw_numeric)
          if (is.null(d)) return(NULL)
          d$kernel <- kernel_name; d$bw <- safe_bw_label(bw_rule); d
        }))
        if (length(lst) == 0) next
        pplot <- plot_overlay(bind_rows(lst), "bw",
                              sprintf("Kernel %s | %s local regressions", kernel_name, gname),
                              xlim, ylim, x_as_date = x_as_date)
        dir.create(file.path(base_out, p, gname, "kernel_overlays"), recursive = TRUE, showWarnings = FALSE)
        ggsave(filename = file.path(base_out, p, gname, "kernel_overlays", sprintf("kernel_%s_overlay_bws_localreg.pdf", kernel_name)),
               plot = pplot, width = 8, height = 5, device = cairo_pdf)
        
      }
    } # end weekday loop
    
    next
  }
  
  if (!(xcol %in% names(dfp)) || !(ycol %in% names(dfp))) next
  x <- as.numeric(dfp[[xcol]])
  y <- dfp[[ycol]]
  raw_pts <- tibble(x = x, y = y)
  
  # Calculate bandwidths for all rules
  numeric_bws <- map_dbl(bw_rules, function(rule) {
    bwv <- get_bw(rule, x)
    if (is.na(bwv) || !is.finite(bwv) || bwv <= .Machine$double.eps * 100) {
      bwv <- get_bw("nrd0", x)  # Fallback if invalid
    }
    bwv
  })
  names(numeric_bws) <- bw_rules
  
  # Clamp bws to reasonable fraction of data range (prevents extreme bw.ucv results)
  rng <- diff(range(x, na.rm = TRUE))
  if (!is.na(rng) && rng > 0) {
    min_bw <- rng * 1e-3
    max_bw <- rng * 0.5
    numeric_bws <- pmin(pmax(numeric_bws, min_bw), max_bw)
  }
  
  names(numeric_bws) <- bw_rules
  ylim <- range(y, na.rm = TRUE)
  xlim <- range(x, na.rm = TRUE)
  
  # Iterate over each combination of kernel and bandwidth
  for (kernel_name in rek_k) {
    for (bw_rule in bw_rules) {
      bw_numeric <- numeric_bws[bw_rule]
      if (is.na(bw_numeric) || !is.finite(bw_numeric)) next
      
      d <- smooth_df(x, y, kernel_name, bw_numeric)
      
      if (is.null(d)) {
        message(sprintf("Failed to smooth with kernel '%s' and bandwidth '%s'", kernel_name, bw_numeric))
        next
      }
      
      d$kernel <- kernel_name
      d$bw <- safe_bw_label(bw_rule)
      ttl <- sprintf("%s local regression: kernel=%s | bw=%s", toupper(p), kernel_name, safe_bw_label(bw_rule))
      
      # Create the plot
      pplot <- ggplot() + 
        geom_line(data = d, aes(x = if(x_as_date) as.Date(x, origin="1970-01-01") else x, 
                                y = y), colour = "steelblue", size = 0.9) +
        geom_point(data = raw_pts, aes(x = if(x_as_date) as.Date(x, origin="1970-01-01") else x, 
                                       y = y), colour = "black", size = 1.2) +
        labs(title = ttl, x = if(x_as_date) "Date" else "x", y = "Value") +
        theme_minimal() + theme(legend.position = "none") +
        coord_cartesian(xlim = if(x_as_date) as.Date(xlim, origin="1970-01-01") else xlim, ylim = ylim)
      
      fname <- file.path(base_out, p, sprintf("%s_locreg_%s_%s.pdf", p, kernel_name, safe_bw_label(bw_rule)))
      ggsave(fname, pplot, width = 8, height = 5, device = cairo_pdf)
    }
  }
  
  
  # Overlay plots for bandwidths and kernels
  for (bw_rule in bw_rules) {
    bw_numeric <- numeric_bws[bw_rule]
    if (is.na(bw_numeric) || !is.finite(bw_numeric)) next
    
    lst <- compact(lapply(rek_k, function(kernel_name) {
      d <- smooth_df(x, y, kernel_name, bw_numeric)
      if (is.null(d)) return(NULL)
      d$kernel <- kernel_name
      d$bw <- safe_bw_label(bw_rule)
      d
    }))
    
    if (length(lst) == 0) next
    
    pplot <- plot_overlay(bind_rows(lst), "kernel",
                          sprintf("BW %s | %s local regressions", safe_bw_label(bw_rule), toupper(p)),
                          xlim, ylim, x_as_date = x_as_date)
    
    ggsave(file.path(base_out, p, "bw_overlays", sprintf("bw_%s_overlay_kernels_localreg.pdf", safe_bw_label(bw_rule))),
           pplot, width = 8, height = 5, device = cairo_pdf)
  }
  
  for (kernel_name in rek_k) {
    lst <- compact(lapply(bw_rules, function(bw_rule) {
      bw_numeric <- numeric_bws[bw_rule]
      if (is.na(bw_numeric) || !is.finite(bw_numeric)) return(NULL)
      d <- smooth_df(x, y, kernel_name, bw_numeric)
      if (is.null(d)) return(NULL)
      d$kernel <- kernel_name
      d$bw <- safe_bw_label(bw_rule)
      d
    }))
    
    if (length(lst) == 0) next
    
    pplot <- plot_overlay(bind_rows(lst), "bw",
                          sprintf("Kernel %s | %s local regressions", kernel_name, toupper(p)),
                          xlim, ylim, x_as_date = x_as_date)
    
    ggsave(filename = file.path(base_out, p, "kernel_overlays", sprintf("kernel_%s_overlay_bws_localreg.pdf", kernel_name)),
           plot = pplot, width = 8, height = 5, device = "pdf")
  }
}





















# Define kernel names and bandwidth selection methods
rek_k <- c("rectangular", "epanechnikov", "biweight", "triangular", "gaussian")
bw_rules <- c("nrd0", "SJ", "ucv")

# Function to get bandwidth based on chosen method
get_bw <- function(rule, x) {
  fallback_bw <- stats::bw.nrd0(x)
  
  bw <- switch(rule,
               nrd0 = fallback_bw,
               SJ   = stats::bw.SJ(x),
               ucv  = stats::bw.ucv(x))
  
  if (is.na(bw) || bw <= .Machine$double.eps) {
    return(NA)
  }
  
  if (rule == "ucv") {
    range_x <- diff(range(x, na.rm = TRUE))
    if (bw >= range_x * 0.5) {
      message("Warning: UCV bandwidth too large; adjusting to a smaller value based on range.")
      bw <- range_x * 0.8
      if (bw <= 0) {
        bw <- fallback_bw
      }
    }
  }
  
  return(bw)
}

# Safe bandwidth label
safe_bw_label <- function(bw_rule) { as.character(bw_rule) }

smooth_quantile_df <- function(x, y, tau, degree = 3) {
  tryCatch({
    model <- rq(y ~ poly(x, degree), tau = tau) 
    new_x <- sort(unique(x))
    pred <- predict(model, newdata = data.frame(x = new_x))
    data.frame(x = new_x, y = pred)
  }, error = function(e) {
    message(sprintf("Error in quantile regression for tau=%s: %s", tau, e$message))
    NULL
  })
}

# Plotting function for quantile regression
plot_with_quantile_regression_smooth <- function(raw_data, smooth_data_list, title, xlim, ylim, x_as_date=FALSE, date_origin="1970-01-01") {
  p <- ggplot() +
    geom_point(data = raw_data, aes(x = if(x_as_date) as.Date(x, origin=date_origin) else x, 
                                    y = y), colour = "black", size = 1.2) +
    geom_line(data = smooth_data_list[["0.025"]], aes(x = x, y = y), colour = "red", size = 0.9) +
    geom_line(data = smooth_data_list[["0.25"]], aes(x = x, y = y), colour = "brown", size = 0.9) +
    geom_line(data = smooth_data_list[["0.5"]], aes(x = x, y = y), colour = "blue", size = 0.9) +
    geom_line(data = smooth_data_list[["0.75"]], aes(x = x, y = y), colour = "cyan", size = 0.9) +
    geom_line(data = smooth_data_list[["0.975"]], aes(x = x, y = y), colour = "green", size = 0.9) +
    labs(title = title, x = if(x_as_date) "Date" else "x", y = "Value") +
    theme_minimal() + 
    theme(legend.position = "none") +
    coord_cartesian(xlim = if(x_as_date) as.Date(xlim, origin=date_origin) else xlim, ylim = ylim)
  return(p)
}

# Overlay plotting function for quantile regression without raw points
plot_overlay_quantile <- function(smooth_data, aescol, ttl, xlim, ylim, x_as_date=FALSE, date_origin="1970-01-01") {
  if (!is.data.frame(smooth_data)) return(NULL)
  
  p <- ggplot(smooth_data, aes(x = if(x_as_date) as.Date(x, origin=date_origin) else x, 
                               y = y, colour = .data[[aescol]])) +
    geom_line(linewidth = 0.9) +
    labs(title = ttl, x = if(x_as_date) "Date" else "x", y = "Value") +
    theme_minimal() +
    theme(legend.position = "right") +
    coord_cartesian(xlim = if(x_as_date) as.Date(xlim, origin=date_origin) else xlim, 
                    ylim = ylim)
  return(p)
}

# Ensure directory exists function
ensure_dir <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

# Directory structure setup
base_out <- "quantile_regressions"
ensure_dir(base_out)

# Example Data Frame Setup
# Assuming df, df_wd, weekly_progress, monthly_progress, yearly_progress are already defined

period_map <- list(
  daily = list(df = df, xcol = "Date", ycol = "Value", x_as_date = TRUE),
  workdays = list(df = df_wd, xcol = "Date", ycol = "Value", x_as_date = TRUE),
  weekly = list(df = weekly_progress, xcol = "week_end", ycol = "week_end_value", x_as_date = TRUE),
  monthly = list(df = monthly_progress, xcol = "month_end", ycol = "month_end_value", x_as_date = TRUE),
  yearly = list(df = yearly_progress, xcol = "year_end", ycol = "year_end_value", x_as_date = TRUE)
)

# Main processing loop for quantile regression with overlay
for (p in names(period_map)) {
  info <- period_map[[p]]
  dfp <- info$df
  xcol <- info$xcol
  ycol <- info$ycol
  x_as_date <- info$x_as_date
  
  if (!(xcol %in% names(dfp)) || !(ycol %in% names(dfp))) next
  x <- as.numeric(dfp[[xcol]])
  y <- dfp[[ycol]]
  
  raw_pts <- tibble(x = x, y = y)
  ylim <- range(y, na.rm = TRUE)
  xlim <- range(x, na.rm = TRUE)
  
  # Create a period-specific output directory
  period_dir <- file.path(base_out, p)
  ensure_dir(period_dir)
  
  # Choose quantiles to process
  quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)
  smooth_data_list <- list()
  
  for (tau in quantiles) {
    d <- smooth_quantile_df(x, y, tau)
    if (is.null(d)) next
    smooth_data_list[[as.character(tau)]] <- d
  }
  
  # Create the plot with actual data and quantile regression lines
  pplot <- plot_with_quantile_regression_smooth(raw_pts, smooth_data_list, 
                                                sprintf("%s Quantile Regression", toupper(p)),
                                                xlim, ylim, x_as_date = x_as_date)
  
  # Save the plot in the period-specific directory
  ggsave(file.path(period_dir, sprintf("%s_quantile_regression.pdf", p)), 
         plot = pplot, width = 8, height = 5, device = "pdf")
  
  # Overlay plots for each kernel and bandwidth
  for (bw_rule in bw_rules) {
    bw_numeric <- get_bw(bw_rule, x)
    if (is.na(bw_numeric) || !is.finite(bw_numeric)) next
    
    lst <- compact(lapply(rek_k, function(kernel_name) {
      d_tau_025 <- smooth_quantile_df(x, y, 0.025)
      d_tau_975 <- smooth_quantile_df(x, y, 0.975)
      
      if (is.null(d_tau_025) || is.null(d_tau_975)) return(NULL)
      
      d_tau_025$kernel <- kernel_name
      d_tau_975$kernel <- kernel_name
      
      return(bind_rows(d_tau_025, d_tau_975))
    }))
    
    if (length(lst) == 0) next
    
    overlay_data <- bind_rows(lst)
    pplot_overlay <- plot_overlay_quantile(overlay_data, "kernel",
                                           sprintf("Quantile Overlays | BW: %s", safe_bw_label(bw_rule)),
                                           xlim, ylim, x_as_date = x_as_date)
    
    ggsave(file.path(period_dir, sprintf("%s_quantile_bw_overlay_%s.pdf", p, safe_bw_label(bw_rule))),
           plot = pplot_overlay, width = 8, height = 5, device = "pdf")
  }
  
  # Overlay for quantile regression based on kernels
  for (kernel_name in rek_k) {
    overlay_data <- compact(lapply(bw_rules, function(bw_rule) {
      bw_numeric <- get_bw(bw_rule, x)
      if (is.na(bw_numeric) || !is.finite(bw_numeric)) return(NULL)
      
      d_tau_025 <- smooth_quantile_df(x, y, 0.025)
      d_tau_975 <- smooth_quantile_df(x, y, 0.975)
      if (is.null(d_tau_025) || is.null(d_tau_975)) return(NULL)
      
      d_tau_025$kernel <- kernel_name
      d_tau_975$kernel <- kernel_name
      
      return(bind_rows(d_tau_025, d_tau_975))
    }))
    
    if (length(overlay_data) == 0) next
    
    pplot_kernel_overlay <- plot_overlay_quantile(overlay_data, "kernel",
                                                  sprintf("Kernel %s | Quantile Overlays", kernel_name),
                                                  xlim, ylim, x_as_date = x_as_date)
    
    ggsave(file.path(period_dir, sprintf("%s_quantile_kernel_overlay_%s.pdf", p, kernel_name)),
           plot = pplot_kernel_overlay, width = 8, height = 5, device = "pdf")
  }
}