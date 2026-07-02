############################################################
# 0. CLEAN SESSION
############################################################
rm(list = ls())
graphics.off()

############################################################
# 1. PACKAGES
############################################################
library(tidyverse)
library(lubridate)
library(Cairo)
library(quantreg)

############################################################
# 2. GENERIC HELPERS
############################################################
compact <- function(x) Filter(Negate(is.null), x)

safe_bw_label <- function(bw) {
  if (is.character(bw)) paste0("bw_", bw)
  else paste0("bw_", gsub("\\.", "p", as.character(bw)))
}

############################################################
# 3. SAFE KERNEL DENSITY ENGINE
############################################################
compute_density_df <- function(x, kernel, bw) {
  tryCatch({
    d <- density(x, kernel = kernel, bw = bw)
    tibble(x = d$x, y = d$y)
  }, warning = function(w) {
    if (bw == "ucv") NULL else invokeRestart("muffleWarning")
  }, error = function(e) NULL)
}

############################################################
# 4. HISTOGRAM + DENSITY PLOT
############################################################
plot_hist_density <- function(
    data, breaks, dens_df, title, file, xlim,
    width = 8, height = 5
) {
  
  h <- hist(data$pct, breaks = breaks, plot = FALSE)
  hist_df <- tibble(
    xmin = head(h$breaks, -1),
    xmax = tail(h$breaks, -1),
    rel  = h$counts / sum(h$counts) * 100
  )
  
  ymax <- max(hist_df$rel) * 1.05
  
  p <- ggplot() +
    geom_rect(
      data = hist_df,
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = rel),
      fill = "grey75", colour = "grey30"
    ) +
    scale_y_continuous(
      limits = c(0, ymax),
      sec.axis = dup_axis(name = "Density")
    ) +
    coord_cartesian(xlim = xlim) +
    labs(title = title, x = "Return", y = "Relative frequency (%)") +
    theme_minimal()
  
  if (!is.null(dens_df)) {
    scale_factor <- ymax / max(dens_df$y)
    p <- p + geom_line(
      data = dens_df,
      aes(x = x, y = y * scale_factor),
      colour = "black", linewidth = 0.9
    )
  }
  
  ggsave(file, p, width = width, height = height, device = cairo_pdf)
}

############################################################
# 5. PERIOD-AGNOSTIC HISTOGRAM PIPELINE
############################################################
process_period <- function(
    values, period_name, out_dir, xlim,
    kernels = c("rectangular","epanechnikov","biweight","triangular","gaussian"),
    bws = c("nrd0","SJ","ucv")
) {
  
  dir.create(out_dir, FALSE, TRUE)
  dir.create(file.path(out_dir, "bw_overlays"), FALSE, TRUE)
  dir.create(file.path(out_dir, "kernel_overlays"), FALSE, TRUE)
  
  df <- tibble(pct = values)
  
  bins <- list(
    bins30 = seq(floor(min(values)*100)/100,
                 ceiling(max(values)*100)/100,
                 length.out = 31),
    fd = {
      n <- nclass.FD(values)
      if (n < 1) n <- 10
      seq(min(values), max(values), length.out = n + 1)
    }
  )
  
  for (bn in names(bins)) {
    brk <- bins[[bn]]
    
    for (k in kernels) {
      for (bw in bws) {
        
        dens <- compute_density_df(values, k, bw)
        fname <- file.path(
          out_dir,
          sprintf("%s_%s_%s_%s.pdf", period_name, bn, k, safe_bw_label(bw))
        )
        
        plot_hist_density(
          df, brk, dens,
          sprintf("%s | %s | %s | %s", period_name, bn, k, bw),
          fname, xlim
        )
      }
    }
    
    # BW overlays
    for (bw in bws) {
      dlist <- compact(lapply(kernels, function(k) {
        d <- compute_density_df(values, k, bw)
        if (!is.null(d)) d$kernel <- k
        d
      }))
      
      if (length(dlist)) {
        p <- ggplot(bind_rows(dlist),
                    aes(x = x, y = y, colour = kernel)) +
          geom_line() +
          theme_minimal() +
          coord_cartesian(xlim = xlim) +
          labs(title = paste(period_name, "| BW", bw))
        
        ggsave(
          file.path(out_dir, "bw_overlays",
                    sprintf("%s_bw_%s.pdf", period_name, bw)),
          p, 8, 5, device = cairo_pdf
        )
      }
    }
    
    # Kernel overlays
    for (k in kernels) {
      dlist <- compact(lapply(bws, function(bw) {
        d <- compute_density_df(values, k, bw)
        if (!is.null(d)) d$bw <- bw
        d
      }))
      
      if (length(dlist)) {
        p <- ggplot(bind_rows(dlist),
                    aes(x = x, y = y, colour = bw)) +
          geom_line() +
          theme_minimal() +
          coord_cartesian(xlim = xlim) +
          labs(title = paste(period_name, "| Kernel", k))
        
        ggsave(
          file.path(out_dir, "kernel_overlays",
                    sprintf("%s_kernel_%s.pdf", period_name, k)),
          p, 8, 5, device = cairo_pdf
        )
      }
    }
  }
}

############################################################
# 6. LOCAL REGRESSION (LOESS, BANDWIDTH-DRIVEN)
############################################################
bw_rules <- c("nrd0","SJ","ucv")
kernels  <- c("rectangular","epanechnikov","biweight","triangular","gaussian")

get_bw_safe <- function(rule, x) {
  bw <- tryCatch(
    switch(rule,
           nrd0 = bw.nrd0(x),
           SJ   = bw.SJ(x),
           ucv  = bw.ucv(x)),
    error = function(e) NA_real_
  )
  
  rng <- diff(range(x))
  if (!is.finite(bw) || bw <= 0) bw <- bw.nrd0(x)
  min(max(bw, rng * 1e-3), rng * 0.5)
}

smooth_local_df <- function(x, y, bw) {
  n <- sum(is.finite(x) & is.finite(y))
  if (n < 10) return(NULL)
  
  span <- bw / diff(range(x))
  span <- min(max(span, 10 / n), 1)
  
  fit <- tryCatch(
    loess(y ~ x, span = span,
          control = loess.control(surface = "interpolate")),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  
  newx <- sort(unique(x))
  tibble(x = newx, y = predict(fit, data.frame(x = newx)))
}

run_local_regressions <- function(df, xcol, ycol, out_dir, x_as_date = TRUE) {
  
  x <- as.numeric(df[[xcol]])
  y <- df[[ycol]]
  raw <- tibble(x = x, y = y)
  
  xlim <- range(x)
  ylim <- range(y)
  
  dir.create(out_dir, FALSE, TRUE)
  dir.create(file.path(out_dir, "bw_overlays"), FALSE, TRUE)
  
  bws <- map_dbl(bw_rules, get_bw_safe, x = x)
  names(bws) <- bw_rules
  
  for (k in kernels) {
    for (bw_name in bw_rules) {
      
      d <- smooth_local_df(x, y, bws[bw_name])
      if (is.null(d)) next
      
      p <- ggplot(raw,
                  aes(x = if (x_as_date) as.Date(x, origin="1970-01-01") else x, y = y)) +
        geom_point(size = 1.1) +
        geom_line(
          data = d,
          aes(x = if (x_as_date) as.Date(x, origin="1970-01-01") else x, y = y),
          colour = "steelblue", linewidth = 0.9
        ) +
        theme_minimal() +
        coord_cartesian(xlim = xlim, ylim = ylim) +
        labs(title = paste("Local regression |", k, "|", bw_name))
      
      ggsave(
        file.path(out_dir,
                  sprintf("locreg_%s_%s.pdf", k, bw_name)),
        p, 8, 5, device = cairo_pdf
      )
    }
  }
}

############################################################
# 7. QUANTILE / PERCENTILE REGRESSION
############################################################
smooth_quantile_df <- function(x, y, tau, degree = 2) {
  tryCatch({
    fit <- rq(y ~ poly(x, degree), tau = tau)
    newx <- sort(unique(x))
    tibble(x = newx, y = predict(fit, data.frame(x = newx)))
  }, error = function(e) NULL)
}

run_quantile_regressions <- function(df, xcol, ycol, out_dir, x_as_date = TRUE) {
  
  x <- as.numeric(df[[xcol]])
  y <- df[[ycol]]
  raw <- tibble(x = x, y = y)
  
  xlim <- range(x)
  ylim <- range(y)
  
  taus <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  
  qlist <- compact(lapply(taus, function(t) {
    d <- smooth_quantile_df(x, y, t)
    if (!is.null(d)) d
  }))
  names(qlist) <- taus
  
  dir.create(out_dir, FALSE, TRUE)
  
  p <- ggplot(raw,
              aes(x = if (x_as_date) as.Date(x, origin="1970-01-01") else x, y = y)) +
    geom_point(size = 1.1) +
    theme_minimal() +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "Quantile regression")
  
  cols <- c("red","orange","blue","cyan","green")
  
  i <- 1
  for (nm in names(qlist)) {
    p <- p + geom_line(
      data = qlist[[nm]],
      aes(x = if (x_as_date) as.Date(x, origin="1970-01-01") else x, y = y),
      colour = cols[i], linewidth = 0.9
    )
    i <- i + 1
  }
  
  ggsave(
    file.path(out_dir, "quantile_regression.pdf"),
    p, 8, 5, device = cairo_pdf
  )
}

############################################################
# 8. DATA PREPARATION
############################################################
df <- read_csv("mbi10_values.csv") %>%
  mutate(Date = as.Date(Date),
         Value = as.numeric(Value),
         pct = Value / lag(Value) - 1) %>%
  drop_na()

############################################################
# 9. PERIOD AGGREGATIONS
############################################################
daily_vals <- df$pct

weekly_vals <- df %>%
  group_by(isoyear(Date), isoweek(Date)) %>%
  summarise(val = last(Value) / first(Value) - 1, .groups = "drop") %>%
  pull(val)

monthly_vals <- df %>%
  group_by(year(Date), month(Date)) %>%
  summarise(val = last(Value) / first(Value) - 1, .groups = "drop") %>%
  pull(val)

yearly_vals <- df %>%
  group_by(year(Date)) %>%
  summarise(val = last(Value) / first(Value) - 1, .groups = "drop") %>%
  pull(val)

############################################################
# 10. RUN EVERYTHING
############################################################
process_period(daily_vals,   "daily",   "histograms/daily",   c(-0.05, 0.05))
process_period(weekly_vals,  "weekly",  "histograms/weekly",  c(-0.2, 0.2))
process_period(monthly_vals, "monthly", "histograms/monthly", c(-0.3, 0.3))
process_period(yearly_vals,  "yearly",  "histograms/yearly",  c(-0.5, 0.5))

run_local_regressions(df, "Date", "Value", "local_regressions/daily")
run_quantile_regressions(df, "Date", "Value", "quantile_regressions/daily")
