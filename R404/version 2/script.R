# ==============================================================================
# OPTIMIZED DISTRIBUTION & REGRESSION ANALYSIS
# Enhanced Visualization & Performance
# FIXED: Proper UCV with filtering out if it dosen't work
# ==============================================================================

# Clear environment
rm(list = ls())
if (!is.null(dev.list())) dev.off(dev.list()["RStudioGD"])

# Set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load packages
library(tidyverse)
library(lubridate)
library(Cairo)
library(locfit)
library(quantreg)
library(future) 
library(future.apply)
plan(multisession, workers = parallel::detectCores() - 1)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Analysis parameters
KERNELS <- c("rectangular", "epanechnikov", "biweight", "triangular", "gaussian")
BW_RULES <- c("nrd0", "SJ", "ucv")

# Axis limits
X_LIM_DAILY <- c(-0.05, 0.05)
X_LIM_WEEKLY <- c(-0.05, 0.05)
X_LIM_MONTHLY <- c(-0.3, 0.3)
X_LIM_YEARLY <- c(-0.5, 0.5)

# Output directories
DIRS <- list(
  hist = "hist_outputs",
  hist_wd = "hist_outputs_weekdays",
  hist_weekly = "hist_outputs_weekly",
  hist_monthly = "hist_outputs_monthly",
  hist_yearly = "hist_outputs_yearly",
  locreg = "local_regressions",
  quantreg = "quantile_regressions"
)

# Create all directories
lapply(DIRS, function(d) dir.create(d, recursive = TRUE, showWarnings = FALSE))

# Enhanced color palettes
PALETTE <- list(
  kernel = c(rectangular = "#E74C3C", epanechnikov = "#3498DB", 
             biweight = "#2ECC71", triangular = "#F39C12", gaussian = "#9B59B6"),
  bw = c(nrd0 = "#E74C3C", SJ = "#3498DB", ucv = "#2ECC71"),
  quantile = c("0.025" = "#E74C3C", "0.25" = "#E67E22", "0.5" = "#3498DB", 
               "0.75" = "#1ABC9C", "0.975" = "#2ECC71"),
  hist_fill = "#3498DB",
  hist_border = "#2C3E50",
  density_line = "#E74C3C",
  stats_mean = "#2C3E50",
  stats_median = "#34495E",
  stats_quant = c("#E74C3C", "#2ECC71")
)

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

compact <- function(x) Filter(Negate(is.null), x)

safe_bw_label <- function(bw) {
  if (is.character(bw)) paste0("bw_", bw)
  else paste0("bw_", gsub("\\.", "p", sprintf("%.4f", bw)))
}

ensure_dir <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

expand_breaks_to_data <- function(breaks, x) {
  rng <- range(x, na.rm = TRUE)
  brk_rng <- range(breaks)
  
  if (rng[1] < brk_rng[1]) breaks[1] <- floor(rng[1] * 100) / 100
  if (rng[2] > brk_rng[2]) breaks[length(breaks)] <- ceiling(rng[2] * 100) / 100
  
  sort(unique(breaks))
}

# Simple progress tracker
progress_bar <- function(current, total, width = 50, label = "") {
  pct <- current / total
  filled <- round(width * pct)
  bar <- paste0(
    "\r", label, " [",
    paste(rep("=", filled), collapse = ""),
    paste(rep(" ", width - filled), collapse = ""),
    "] ", sprintf("%3.0f%%", pct * 100),
    sprintf(" (%d/%d)", current, total)
  )
  cat(bar)
  if (current == total) cat("\n")
  flush.console()
}

# Enhanced theme for all plots
theme_enhanced <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 3, 
                                colour = "#2C3E50", margin = margin(b = 5)),
      plot.subtitle = element_text(colour = "#7F8C8D", size = base_size, 
                                   margin = margin(b = 10)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#ECF0F1", linewidth = 0.3),
      axis.title = element_text(face = "bold", size = base_size, colour = "#2C3E50"),
      axis.text = element_text(colour = "#34495E", size = base_size - 1),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text = element_text(size = base_size - 1),
      legend.background = element_rect(fill = "white", colour = "#BDC3C7", linewidth = 0.3),
      legend.key.size = unit(0.8, "cm"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "#FAFAFA", colour = NA),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ==============================================================================
# MANUAL UCV BANDWIDTH SELECTOR (FALLBACK IMPLEMENTATION)
# ==============================================================================

#' Manual UCV bandwidth selector using numerical optimization
#' 
#' This is a FALLBACK implementation used only when stats::bw.ucv() fails.
#' It implements the unbiased cross-validation criterion for Gaussian kernels.
#' 
#' @param x Numeric vector of data points
#' @param lower Lower bound for bandwidth search
#' @param upper Upper bound for bandwidth search
#' @param tol Tolerance for optimization
#' @return Optimal bandwidth or NA if optimization fails
manual_bw_ucv <- function(x, lower, upper, tol = 0.1 * lower) {
  
  # UCV criterion function for Gaussian kernel
  # Based on: Scott (1992), Silverman (1986)
  ucv_objective <- function(h, x) {
    n <- length(x)
    if (h <= 0) return(Inf)
    
    # Avoid numerical issues
    if (h < 1e-10 || !is.finite(h)) return(Inf)
    
    # Calculate all pairwise differences
    diff_matrix <- outer(x, x, "-")
    
    # UCV for Gaussian kernel involves:
    # 1. Leave-one-out kernel density estimates
    # 2. Integration of squared density estimate
    # 3. Subtraction of diagonal terms (leave-one-out)
    
    # Term 1: R(f̂) approximation - integrated squared density
    # Using the convolution formula for Gaussian kernels
    # K * K has variance 2h²
    sqrt2h <- sqrt(2) * h
    
    # Off-diagonal terms (i ≠ j)
    term1 <- sum(dnorm(diff_matrix / sqrt2h)) / (n^2 * sqrt2h)
    
    # Term 2: Leave-one-out correction
    # Remove diagonal contributions
    # Each diagonal element contributes K(0) = dnorm(0)
    term2 <- 2 * dnorm(0) / (n * h)
    
    # UCV criterion (minimize this)
    ucv <- term1 - term2
    
    # Add small penalty for extreme bandwidths to ensure stability
    if (h < lower * 1.1) ucv <- ucv + (lower * 1.1 - h)^2
    if (h > upper * 0.9) ucv <- ucv + (h - upper * 0.9)^2
    
    if (!is.finite(ucv)) return(Inf)
    
    return(ucv)
  }
  
  # Try optimization with multiple methods
  result <- tryCatch({
    
    # Method 1: optimize() - works well for 1D problems
    opt1 <- optimize(
      f = function(h) ucv_objective(h, x),
      interval = c(lower, upper),
      tol = tol
    )
    
    if (is.finite(opt1$objective) && opt1$minimum > lower && opt1$minimum < upper) {
      return(opt1$minimum)
    }
    
    # Method 2: Grid search with refinement (fallback)
    grid_size <- 50
    h_grid <- exp(seq(log(lower * 1.05), log(upper * 0.95), length.out = grid_size))
    ucv_values <- sapply(h_grid, function(h) ucv_objective(h, x))
    
    valid_idx <- is.finite(ucv_values)
    if (sum(valid_idx) == 0) return(NA_real_)
    
    best_h <- h_grid[valid_idx][which.min(ucv_values[valid_idx])]
    
    # Refine around best point
    h_refine <- seq(
      max(best_h * 0.8, lower),
      min(best_h * 1.2, upper),
      length.out = 20
    )
    ucv_refine <- sapply(h_refine, function(h) ucv_objective(h, x))
    valid_refine <- is.finite(ucv_refine)
    
    if (sum(valid_refine) > 0) {
      best_h <- h_refine[valid_refine][which.min(ucv_refine[valid_refine])]
    }
    
    return(best_h)
    
  }, error = function(e) {
    return(NA_real_)
  })
  
  # Validate result
  if (is.na(result) || !is.finite(result) || result <= 0) {
    return(NA_real_)
  }
  
  return(result)
}

# ==============================================================================
# IMPROVED BANDWIDTH CALCULATION WITH FALLBACK
# ==============================================================================

get_bw <- function(rule, x, fallback = TRUE, verbose = FALSE) {
  fallback_bw <- stats::bw.nrd0(x)
  
  bw <- tryCatch({
    switch(rule,
           nrd0 = fallback_bw,
           
           SJ = stats::bw.SJ(x),
           
           ucv = {
             suppressWarnings({
               range_x <- diff(range(x, na.rm = TRUE))
               n <- length(x)
               
               # ========== SMART CONSTRAINTS ==========
               min_reasonable_bw <- max(
                 range_x * 0.05,                    # At least 5% of range
                 range_x / sqrt(n),                 # Sample size scaling
                 fallback_bw * 0.5                  # Never < 50% of nrd0
               )
               
               lower <- min_reasonable_bw
               
               # Upper bound: prevent undersmoothing
               upper <- min(
                 range_x * 0.25,
                 fallback_bw * 3
               )
               
               if (verbose) {
                 message(sprintf("      UCV search range: [%.6f, %.6f] (range=%.2f, nrd0=%.6f)", 
                                 lower, upper, range_x, fallback_bw))
               }
               # =========================================
               
               # PRIMARY: Try R's built-in UCV
               ucv_result <- tryCatch({
                 stats::bw.ucv(x, lower = lower, upper = upper)
               }, error = function(e) {
                 if (verbose) message("      UCV (native): Failed - ", conditionMessage(e))
                 stop("UCV native method failed")  # Explicitly stop to trigger fallback
               }, warning = function(w) {
                 # Check if it's the "minimum at boundary" warning
                 if (grepl("occurred at one end", conditionMessage(w))) {
                   stop("UCV native method failed at boundary")  # Explicitly stop
                 } else {
                   invokeRestart("muffleWarning")
                 }
               })
               
               # This block will only be reached if no error was thrown
               # Add additional validation
               if (!is.finite(ucv_result) || ucv_result <= 0) {
                 stop("Invalid UCV result")
               }
               
               # ========== SAFETY CAP ==========
               if (ucv_result > fallback_bw * 3) {
                 if (verbose) {
                   message(sprintf("      UCV: Capped from %.6f to %.6f (3× nrd0)", 
                                   ucv_result, fallback_bw * 3))
                 }
                 ucv_result <- fallback_bw * 3
               }
               
               # Additional constraint: never below 50% of nrd0
               if (ucv_result < fallback_bw * 0.5) {
                 if (verbose) {
                   message(sprintf("      UCV: Raised from %.6f to %.6f (50%% nrd0 minimum)", 
                                   ucv_result, fallback_bw * 0.5))
                 }
                 ucv_result <- fallback_bw * 0.5
               }
               
               # Report final selection
               if (verbose) {
                 ratio <- ucv_result / fallback_bw
                 pct_of_range <- 100 * ucv_result / range_x
                 message(sprintf("      UCV: Selected h=%.6f (%.2f%% of range, %.2f× nrd0)", 
                                 ucv_result, pct_of_range, ratio))
               }
               # ================================
               
               return(ucv_result)
             })
           },
           
           # Default case if no rule matches
           fallback_bw)
  }, error = function(e) {
    if (verbose) message("      Error in ", rule, ": ", conditionMessage(e))
    # If fallback is FALSE, return NA instead of falling back
    if (!fallback) return(NA_real_)
    # If fallback is TRUE, move to the next method in the sequence
    NA_real_
  })
  
  # If the result is invalid and fallback is TRUE, return NA to trigger next method
  if (is.na(bw) || !is.finite(bw) || bw <= .Machine$double.eps * 100) {
    return(if (fallback) NA_real_ else NA_real_)
  }
  
  # Additional safety clamp
  range_x <- diff(range(x, na.rm = TRUE))
  if (!is.na(range_x) && range_x > 0) {
    min_bw <- range_x * 0.001
    max_bw <- range_x * 0.5
    bw <- pmax(pmin(bw, max_bw), min_bw)
  }
  
  bw
}

calc_bandwidths <- function(x, rules = BW_RULES, verbose = FALSE) {
  if (verbose) {
    n <- length(x)
    rng <- diff(range(x, na.rm = TRUE))
    message(sprintf("  Calculating bandwidths for n=%d points, range=%.6f:", n, rng))
  }
  
  bws <- setNames(sapply(rules, function(r) {
    if (verbose) message(sprintf("    Rule: %s", toupper(r)))
    get_bw(r, x, verbose = verbose)
  }), rules)
  
  # Final range check
  rng <- diff(range(x, na.rm = TRUE))
  if (!is.na(rng) && rng > 0) {
    bws <- pmin(pmax(bws, rng * 1e-3), rng * 0.5)
  }
  
  if (verbose) {
    message("  Final bandwidths:")
    for (rule in names(bws)) {
      message(sprintf("    %s = %.6f", toupper(rule), bws[rule]))
    }
    message("")
  }
  
  bws
}

calc_bandwidths <- function(x, rules = BW_RULES, verbose = FALSE) {
  if (verbose) {
    n <- length(x)
    rng <- diff(range(x, na.rm = TRUE))
    message(sprintf("  Calculating bandwidths for n=%d points, range=%.6f:", n, rng))
  }
  
  bws <- setNames(sapply(rules, function(r) {
    if (verbose) message(sprintf("    Rule: %s", toupper(r)))
    get_bw(r, x, verbose = verbose)
  }), rules)
  
  # Final range check
  rng <- diff(range(x, na.rm = TRUE))
  if (!is.na(rng) && rng > 0) {
    bws <- pmin(pmax(bws, rng * 1e-3), rng * 0.5)
  }
  
  if (verbose) {
    message("  Final bandwidths:")
    for (rule in names(bws)) {
      message(sprintf("    %s = %.6f", toupper(rule), bws[rule]))
    }
    message("")
  }
  
  bws
}

# ==============================================================================
# DENSITY ESTIMATION
# ==============================================================================
compute_density_df <- function(xvec, kernel, bw, x_grid = NULL) {
  if (is.na(bw) || !is.finite(bw) || bw <= 0) return(NULL)
  
  xvec <- na.omit(xvec)
  if (length(xvec) < 2) return(NULL)
  
  d <- tryCatch({
    dens <- density(xvec, kernel = kernel, bw = bw)
    data.frame(x = dens$x, y = dens$y)
  }, warning = function(w) {
    if (grepl("minimum occurred at one end", conditionMessage(w))) return(NULL)
    invokeRestart("muffleWarning")
  }, error = function(e) NULL)
  
  if (!is.null(d)) return(d)
  
  if (is.character(bw)) {
    for (alt_bw in setdiff(c("nrd0", "SJ"), bw)) {
      d <- tryCatch({
        dens <- density(xvec, kernel = kernel, bw = alt_bw)
        data.frame(x = dens$x, y = dens$y)
      }, error = function(e) NULL)
      if (!is.null(d)) return(d)
    }
  }
  NULL
}


compute_per_weekday <- function(df_wd, kernel, bw, x_grid = NULL) {
  if (is.na(bw) || !is.finite(bw) || bw <= 0) return(NULL)
  
  wd_levels <- unique(df_wd$weekday)
  
  bind_rows(lapply(wd_levels, function(wd) {
    vals <- df_wd %>% filter(weekday == wd) %>% pull(pct) %>% na.omit()
    if (length(vals) < 2) return(NULL)
    
    d <- compute_density_df(vals, kernel, bw, x_grid)
    if (is.null(d)) return(NULL)
    
    d$weekday <- wd
    d
  })) %>%
    mutate(weekday = factor(weekday, levels = wd_levels))
}


# ==============================================================================
# ENHANCED HISTOGRAM PLOTTING
# ==============================================================================

plot_histogram <- function(data, breaks, dens_df = NULL, title = "", subtitle = "",
                           file, w = 10, h = 6, x_lim = c(-0.05, 0.05)) {
  
  # Ensure breaks span the full data range
  breaks <- expand_breaks_to_data(breaks, data$pct)
  
  # Build histogram
  hst <- hist(data$pct, breaks = breaks, plot = FALSE, 
              include.lowest = TRUE, right = FALSE)
  
  hist_df <- tibble(
    left = head(hst$breaks, -1),
    right = tail(hst$breaks, -1),
    rel = (hst$counts / sum(hst$counts)) * 100
  ) %>% mutate(xmin = left, xmax = right)
  
  ymax_hist <- max(hist_df$rel, na.rm = TRUE) * 1.05
  
  # Statistics
  stats_data <- tibble(
    mean_val = mean(data$pct, na.rm = TRUE),
    median_val = median(data$pct, na.rm = TRUE),
    q025 = quantile(data$pct, 0.025, na.rm = TRUE),
    q975 = quantile(data$pct, 0.975, na.rm = TRUE)
  )
  
  # Base plot with enhanced styling
  p <- ggplot() +
    geom_rect(data = hist_df,
              aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = rel),
              fill = PALETTE$hist_fill, colour = PALETTE$hist_border, 
              alpha = 0.7, linewidth = 0.3) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Return",
      y = "Relative Frequency (%)"
    ) +
    theme_enhanced() +
    coord_cartesian(xlim = x_lim) +
    scale_y_continuous(
      limits = c(0, ymax_hist),
      expand = expansion(mult = c(0, 0.02)),
      sec.axis = dup_axis(name = "Density")
    ) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 0.1))
  
  # Add density overlay
  if (!is.null(dens_df)) {
    dens_max <- max(dens_df$y, na.rm = TRUE)
    scale_fact <- ymax_hist / dens_max
    p <- p + geom_line(data = dens_df, aes(x = x, y = y * scale_fact),
                       colour = PALETTE$density_line, linewidth = 1.2, alpha = 0.9)
  }
  
  # Add reference lines with enhanced styling
  p <- p +
    geom_vline(xintercept = stats_data$mean_val, 
               linetype = "dashed", colour = PALETTE$stats_mean, 
               linewidth = 0.8, alpha = 0.8) +
    geom_vline(xintercept = stats_data$median_val,
               linetype = "dotted", colour = PALETTE$stats_median, 
               linewidth = 0.8, alpha = 0.8) +
    geom_vline(xintercept = stats_data$q025, linetype = "longdash",
               colour = PALETTE$stats_quant[1], linewidth = 0.7, alpha = 0.7) +
    geom_vline(xintercept = stats_data$q975, linetype = "longdash",
               colour = PALETTE$stats_quant[2], linewidth = 0.7, alpha = 0.7)
  
  # Add legend for reference lines
  p <- p +
    annotate("text", x = x_lim[2] * 0.95, y = ymax_hist * 0.95,
             label = sprintf("Mean: %.2f%%\nMedian: %.2f%%\n2.5%%: %.2f%% | 97.5%%: %.2f%%",
                             stats_data$mean_val * 100, stats_data$median_val * 100,
                             stats_data$q025 * 100, stats_data$q975 * 100),
             hjust = 1, vjust = 1, size = 3.5, colour = "#2C3E50",
             family = "sans", fontface = "italic",
             fill = "white", alpha = 0.9, 
             label.padding = unit(0.5, "lines"),
             label.r = unit(0.15, "lines"))
  
  # Save
  ggsave(filename = file, plot = p, width = w, height = h, device = cairo_pdf)
  invisible(p)
}

# ==============================================================================
# LOCAL REGRESSION
# ==============================================================================

smooth_loess <- function(x, y, bw_numeric) {
  if (is.na(bw_numeric) || !is.finite(bw_numeric) || bw_numeric <= 0) return(NULL)
  
  tryCatch({
    x <- as.numeric(x)
    n <- sum(!is.na(x) & !is.na(y))
    if (n < 5) return(NULL)
    
    rng <- diff(range(x, na.rm = TRUE))
    if (is.na(rng) || rng == 0) return(NULL)
    
    span <- min(max(bw_numeric / rng, max(10, ceiling(0.02 * n)) / n), 1)
    
    model <- loess(y ~ x, span = span, family = "gaussian",
                   control = loess.control(surface = "interpolate"))
    newx <- sort(unique(x))
    
    data.frame(x = newx, y = predict(model, newdata = data.frame(x = newx)))
  }, error = function(e) NULL)
}


# ==============================================================================
# QUANTILE REGRESSION FRAMEWORK
# ==============================================================================

kernel_weight <- function(u, kernel = "epanechnikov") {
  w <- switch(kernel,
              rectangular = ifelse(abs(u) <= 1, 0.5, 0),
              epanechnikov = ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0),
              biweight = ifelse(abs(u) <= 1, (15/16) * (1 - u^2)^2, 0),
              triangular = ifelse(abs(u) <= 1, 1 - abs(u), 0),
              gaussian = dnorm(u),
              0)
  w / sum(w)
}

local_quantile_reg <- function(x, y, tau, kernel, bw, x_eval = NULL) {
  if (is.na(bw) || !is.finite(bw) || bw <= 0) return(NULL)
  
  tryCatch({
    if (is.null(x_eval)) x_eval <- sort(unique(x))
    
    predictions <- vapply(x_eval, function(x0) {
      u <- (x - x0) / bw
      w <- kernel_weight(u, kernel)
      
      if (sum(w > 1e-10) < 5) return(NA_real_)
      
      fit <- tryCatch(rq(y ~ x, tau = tau, weights = w), error = function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      
      predict(fit, newdata = data.frame(x = x0))
    }, numeric(1))
    
    # NEW FIX
    if (all(is.na(predictions))) return(NULL)
    
    valid <- !is.na(predictions)
    if (sum(valid) < 2) return(NULL)
    
    data.frame(x = x_eval[valid], y = predictions[valid])
  }, error = function(e) NULL)
}

parametric_quantile_reg <- function(x, y, tau, degree = 2) {
  tryCatch({
    model <- rq(y ~ poly(x, degree), tau = tau)
    new_x <- sort(unique(x))
    data.frame(x = new_x, y = predict(model, newdata = data.frame(x = new_x)))
  }, error = function(e) NULL)
}

bootstrap_quantile_ci <- function(x, y, tau, kernel, bw, n_boot = 100, alpha = 0.05) {
  if (is.na(bw) || !is.finite(bw) || bw <= 0) return(NULL)
  
  tryCatch({
    x_eval <- sort(unique(x))
    n <- length(x)
    
    boot_results <- future_replicate(n_boot, {
      idx <- sample(n, n, replace = TRUE)
      fit <- local_quantile_reg(x[idx], y[idx], tau, kernel, bw, x_eval)
      
      if (is.null(fit)) return(rep(NA, length(x_eval)))
      
      pred <- rep(NA, length(x_eval))
      pred[match(fit$x, x_eval)] <- fit$y
      pred
    }, future.seed = TRUE)
    
    data.frame(
      x = x_eval,
      lower = apply(boot_results, 1, quantile, alpha/2, na.rm = TRUE),
      upper = apply(boot_results, 1, quantile, 1 - alpha/2, na.rm = TRUE)
    )
  }, error = function(e) NULL)
}


check_quantile_crossing <- function(quantile_fits) {
  if (length(quantile_fits) < 2) return(NULL)
  
  taus_sorted <- sort(as.numeric(names(quantile_fits)))
  crossings <- list()
  
  for (i in seq_len(length(taus_sorted) - 1)) {
    tau_lower <- taus_sorted[i]
    tau_upper <- taus_sorted[i + 1]
    
    fit_lower <- quantile_fits[[as.character(tau_lower)]]
    fit_upper <- quantile_fits[[as.character(tau_upper)]]
    
    if (is.null(fit_lower) || is.null(fit_upper)) next
    
    merged <- merge(fit_lower, fit_upper, by = "x", suffixes = c("_lower", "_upper"))
    n_crossings <- sum(merged$y_upper < merged$y_lower, na.rm = TRUE)
    
    if (n_crossings > 0) {
      crossings[[paste(tau_lower, "vs", tau_upper)]] <- n_crossings
    }
  }
  
  crossings
}

calc_mise <- function(fit_nonparam, fit_param) {
  if (is.null(fit_nonparam) || is.null(fit_param)) return(NA_real_)
  
  merged <- merge(fit_nonparam, fit_param, by = "x", suffixes = c("_np", "_p"))
  if (nrow(merged) < 2) return(NA_real_)
  
  mean((merged$y_np - merged$y_p)^2, na.rm = TRUE)
}

# ==============================================================================
# DATA LOADING & PREPARATION
# ==============================================================================

message("Loading and preparing data...")

# Read daily data
df <- read_csv("mbi10_values.csv", show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date, "%Y.%m.%d")) %>%
  filter(Date >= as.Date("2009-06-15")) %>%
  arrange(Date) %>%
  mutate(pct = if ("% change" %in% names(.))
    as.numeric(gsub("%", "", `% change`)) / 100
    else (value / lag(value) - 1)) %>%
  filter(!is.na(pct))

# Weekday data
df_wd <- df %>%
  mutate(weekday = wday(Date, label = TRUE, week_start = 1)) %>%
  filter(weekday %in% c("Mon", "Tue", "Wed", "Thu", "Fri")) %>%
  mutate(weekday = recode(as.character(weekday),
                          Mon = "Monday", Tue = "Tuesday", Wed = "Wednesday",
                          Thu = "Thursday", Fri = "Friday"))

# Helper functions for period calculations
week_id <- function(date) paste0(isoyear(date), "-W", sprintf("%02d", isoweek(date)))
month_id <- function(date) paste0(year(date), "-", sprintf("%02d", month(date)))

# Calculate period returns
calc_period_returns <- function(df, period_func, id_col, start_col, end_col, value_col) {
  df %>%
    arrange(Date) %>%
    mutate(period_id = period_func(Date)) %>%
    group_by(period_id) %>%
    summarise(
      period_start = min(Date),
      start_value = first(Value),
      period_end = max(Date),
      end_value = last(Value),
      .groups = "drop"
    ) %>%
    filter(!is.na(start_value), !is.na(end_value), period_end >= period_start) %>%
    mutate(period_change = (end_value / start_value) - 1) %>%
    arrange(period_end)
}

weekly_progress <- calc_period_returns(df_wd, week_id, "wid", 
                                       "week_start", "week_end", "weekly_change")
monthly_progress <- calc_period_returns(df_wd, month_id, "mid",
                                        "month_start", "month_end", "monthly_change")
yearly_progress <- df_wd %>%
  arrange(Date) %>%
  mutate(year = year(Date)) %>%
  group_by(year) %>%
  summarise(
    year_start = min(Date),
    year_start_value = first(Value),
    year_end = max(Date),
    year_end_value = last(Value),
    .groups = "drop"
  ) %>%
  filter(!is.na(year_start_value), !is.na(year_end_value)) %>%
  mutate(yearly_change = (year_end_value / year_start_value) - 1)

# ==============================================================================
# HISTOGRAM BIN DEFINITIONS
# ==============================================================================

create_bins <- function(vals, x_lim = NULL) {
  range_vals <- if (is.null(x_lim)) {
    c(floor(min(vals, na.rm = TRUE) * 100) / 100,
      ceiling(max(vals, na.rm = TRUE) * 100) / 100)
  } else x_lim
  
  nbins_fd <- nclass.FD(vals)
  if (is.na(nbins_fd) || nbins_fd < 1) nbins_fd <- 10
  
  list(
    bins30 = seq(range_vals[1], range_vals[2], length.out = 31),
    fd = seq(range_vals[1], range_vals[2], length.out = nbins_fd + 1)
  )
}

bins_daily <- create_bins(df$pct, X_LIM_DAILY)
bins_weekday <- create_bins(df_wd$pct, X_LIM_DAILY)
bins_weekly <- create_bins(weekly_progress$period_change, X_LIM_WEEKLY)
bins_monthly <- create_bins(monthly_progress$period_change, X_LIM_MONTHLY)
bins_yearly <- create_bins(yearly_progress$yearly_change, X_LIM_YEARLY)

# ==============================================================================
# MAIN HISTOGRAM GENERATION FUNCTION
# ==============================================================================

generate_histograms <- function(data, bins, kernels, bw_rules, out_dir, 
                                x_lim, period_name = "Period") {
  
  # Compute numeric bandwidths once
  bws <- calc_bandwidths(data$pct, bw_rules, verbose = FALSE)
  bw_dir <- file.path(out_dir, "bw_overlays")
  kern_dir <- file.path(out_dir, "kernel_overlays")
  ensure_dir(bw_dir)
  ensure_dir(kern_dir)
  
  for (bn in names(bins)) {
    brk <- bins[[bn]]
    nbins <- if (bn == "fd") length(brk) - 1 else 30
    
    # Individual histograms
    # === PARALLEL: Individual histograms ===
    future_lapply(kernels, function(k) {
      lapply(bw_rules, function(b) {
        
        bw_val <- bws[[b]]
        if (is.na(bw_val)) return(NULL)
        
        dens <- compute_density_df(data$pct, k, bw_val)
        
        ttl <- sprintf("%s Distribution Analysis", period_name)
        subt <- sprintf("Bins: %d | Kernel: %s | Bandwidth: %s",
                        nbins, str_to_title(k), toupper(b))
        
        fname <- file.path(out_dir,
                           sprintf("hist_%s_%s_%s.pdf", bn, k, safe_bw_label(b)))
        
        plot_histogram(data, brk, dens, ttl, subt, fname, x_lim = x_lim)
      })
    })
    
    # Bandwidth overlays
    # === PARALLEL: Bandwidth overlays ===
    for (b in bw_rules) {
      bw_val <- bws[[b]]
      
      dlist <- compact(future_lapply(kernels, function(k) {
        if (is.na(bw_val)) return(NULL)
        d <- compute_density_df(data$pct, k, bw_val)
        if (!is.null(d)) d$kernel <- k
        d
      }))
      
      if (length(dlist) == 0) next
      
      p <- ggplot(bind_rows(dlist), aes(x = x, y = y, colour = kernel, linetype = kernel)) +
        geom_line(linewidth = 1.1) +
        scale_colour_manual(values = PALETTE$kernel, name = "Kernel") +
        scale_linetype_manual(
          values = c(rectangular = "solid", epanechnikov = "dashed",
                     biweight = "dotted", triangular = "dotdash", gaussian = "longdash"),
          name = "Kernel"
        ) +
        labs(
          title = sprintf("%s Density Comparison", period_name),
          subtitle = sprintf("Bandwidth: %s | Bins: %d", toupper(b), nbins),
          x = "Return", y = "Density"
        ) +
        theme_enhanced() +
        coord_cartesian(xlim = x_lim) +
        scale_x_continuous(labels = scales::percent_format(accuracy = 0.1))
      
      ggsave(file.path(bw_dir, sprintf("bw_%s_%s_overlay.pdf", safe_bw_label(b), bn)),
             p, width = 10, height = 6, device = cairo_pdf)
    }
    
    # === PARALLEL: Kernel overlays ===
    for (k in kernels) {
      
      dlist <- compact(future_lapply(bw_rules, function(b) {
        bw_val <- bws[[b]]
        if (is.na(bw_val)) return(NULL)
        d <- compute_density_df(data$pct, k, bw_val)
        if (!is.null(d)) d$bw <- toupper(b)
        d
      }))
      
      if (length(dlist) == 0) next
      
      p <- ggplot(bind_rows(dlist), aes(x = x, y = y, colour = bw, linetype = bw)) +
        geom_line(linewidth = 1.1) +
        scale_colour_manual(values = PALETTE$bw, name = "Bandwidth") +
        scale_linetype_manual(
          values = c(NRD0 = "solid", SJ = "dashed", UCV = "dotted"),
          name = "Bandwidth"
        ) +
        labs(
          title = sprintf("%s Density Comparison", period_name),
          subtitle = sprintf("Kernel: %s | Bins: %d", str_to_title(k), nbins),
          x = "Return", y = "Density"
        ) +
        theme_enhanced() +
        coord_cartesian(xlim = x_lim) +
        scale_x_continuous(labels = scales::percent_format(accuracy = 0.1))
      
      ggsave(file.path(kern_dir, sprintf("kernel_%s_%s_overlay.pdf", k, bn)),
             p, width = 10, height = 6, device = cairo_pdf)
    }
  }
}

# ==============================================================================
# GENERATE ALL HISTOGRAMS
# ==============================================================================

message("\n=== GENERATING HISTOGRAMS ===\n")

message("Daily distributions:")
calc_bandwidths(df$pct, BW_RULES, verbose = TRUE)
generate_histograms(df, bins_daily, KERNELS, BW_RULES, DIRS$hist, 
                    X_LIM_DAILY, "Daily Returns")

message("Weekly distributions:")
calc_bandwidths(weekly_progress$period_change, BW_RULES, verbose = TRUE)
generate_histograms(
  data.frame(pct = weekly_progress$period_change),
  bins_weekly, KERNELS, BW_RULES, DIRS$hist_weekly,
  X_LIM_WEEKLY, "Weekly Returns"
)

message("Monthly distributions:")
calc_bandwidths(monthly_progress$period_change, BW_RULES, verbose = TRUE)
generate_histograms(
  data.frame(pct = monthly_progress$period_change),
  bins_monthly, KERNELS, BW_RULES, DIRS$hist_monthly,
  X_LIM_MONTHLY, "Monthly Returns"
)

message("Yearly distributions:")
calc_bandwidths(yearly_progress$yearly_change, BW_RULES, verbose = TRUE)
generate_histograms(
  data.frame(pct = yearly_progress$yearly_change),
  bins_yearly, KERNELS, BW_RULES, DIRS$hist_yearly,
  X_LIM_YEARLY, "Yearly Returns"
)

# ==============================================================================
# WEEKDAY-SPECIFIC HISTOGRAMS
# ==============================================================================

message("\n=== GENERATING WEEKDAY-SPECIFIC HISTOGRAMS ===\n")

x_grid <- seq(min(df_wd$pct, na.rm = TRUE), max(df_wd$pct, na.rm = TRUE), length.out = 512)
bw_dir_wd <- file.path(DIRS$hist_wd, "bw_overlays")
kern_dir_wd <- file.path(DIRS$hist_wd, "kernel_overlays")
ensure_dir(bw_dir_wd)
ensure_dir(kern_dir_wd)

for (bn in names(bins_weekday)) {
  brk <- expand_breaks_to_data(bins_weekday[[bn]], df_wd$pct)
  nbins <- if (bn == "fd") length(brk) - 1 else 30
  
  for (k in KERNELS) {
    for (bw in BW_RULES) {
      dens_wd <- compute_per_weekday(df_wd, k, bw, x_grid)
      
      if (!is.null(dens_wd) && nrow(dens_wd) > 0) {
        for (day in sort(unique(dens_wd$weekday))) {
          df_day <- df_wd %>% filter(weekday == day)
          dens_day <- dens_wd %>% filter(weekday == day) %>% select(x, y)
          
          ttl <- sprintf("%s Returns Distribution", day)
          subt <- sprintf("Bins: %d | Kernel: %s | Bandwidth: %s", 
                          nbins, str_to_title(k), toupper(bw))
          
          fname <- file.path(DIRS$hist_wd, 
                             sprintf("weekday_%s_%s_%s_%s.pdf", day, bn, k, safe_bw_label(bw)))
          plot_histogram(df_day, brk, dens_day, ttl, subt, fname, x_lim = X_LIM_DAILY)
        }
      }
    }
  }
  
  # Weekday overlays (bandwidth comparison)
  for (bw in BW_RULES) {
    dlist <- compact(lapply(KERNELS, function(k) {
      dd <- compute_per_weekday(df_wd, k, bw, x_grid)
      if (!is.null(dd) && nrow(dd) > 0) {
        dd %>% 
          group_by(x) %>% 
          summarise(y = mean(y, na.rm = TRUE), .groups = "drop") %>%
          mutate(kernel = k)
      }
    }))
    
    if (length(dlist) > 0) {
      p <- ggplot(bind_rows(dlist), aes(x = x, y = y, colour = kernel, linetype = kernel)) +
        geom_line(linewidth = 1.1) +
        scale_colour_manual(values = PALETTE$kernel, name = "Kernel") +
        labs(
          title = "Weekday Average Density Comparison",
          subtitle = sprintf("Bandwidth: %s | Bins: %d", toupper(bw), nbins),
          x = "Daily Return", y = "Density"
        ) +
        theme_enhanced() +
        coord_cartesian(xlim = X_LIM_DAILY) +
        scale_x_continuous(labels = scales::percent_format(accuracy = 0.1))
      
      ggsave(file.path(kern_dir_wd, sprintf("avg_wd_%s_bw_%s.pdf", bn, safe_bw_label(bw))),
             p, width = 10, height = 6, device = cairo_pdf)
    }
  }
  
  # Weekday overlays (kernel comparison)
  for (k in KERNELS) {
    dlist <- compact(lapply(BW_RULES, function(bw) {
      dd <- compute_per_weekday(df_wd, k, bw, x_grid)
      if (!is.null(dd) && nrow(dd) > 0) {
        dd %>% 
          group_by(x) %>% 
          summarise(y = mean(y, na.rm = TRUE), .groups = "drop") %>%
          mutate(bw = toupper(bw))
      }
    }))
    
    if (length(dlist) > 0) {
      p <- ggplot(bind_rows(dlist), aes(x = x, y = y, colour = bw, linetype = bw)) +
        geom_line(linewidth = 1.1) +
        scale_colour_manual(values = PALETTE$bw, name = "Bandwidth") +
        labs(
          title = "Weekday Average Density Comparison",
          subtitle = sprintf("Kernel: %s | Bins: %d", str_to_title(k), nbins),
          x = "Daily Return", y = "Density"
        ) +
        theme_enhanced() +
        coord_cartesian(xlim = X_LIM_DAILY) +
        scale_x_continuous(labels = scales::percent_format(accuracy = 0.1))
      
      ggsave(file.path(bw_dir_wd, sprintf("avg_wd_%s_kernel_%s.pdf", bn, k)),
             p, width = 10, height = 6, device = cairo_pdf)
    }
  }
}

# ==============================================================================
# LOCAL REGRESSIONS
# ==============================================================================

message("\n=== GENERATING LOCAL REGRESSIONS ===\n")

period_config <- list(
  daily = list(df = df, xcol = "Date", ycol = "Value", x_as_date = TRUE),
  weekly = list(df = weekly_progress, xcol = "period_end", ycol = "end_value", x_as_date = TRUE),
  monthly = list(df = monthly_progress, xcol = "period_end", ycol = "end_value", x_as_date = TRUE),
  yearly = list(df = yearly_progress, xcol = "year_end", ycol = "year_end_value", x_as_date = TRUE)
)

for (period in names(period_config)) {
  message(sprintf("%s local regressions:", str_to_title(period)))
  
  cfg <- period_config[[period]]
  pdir <- file.path(DIRS$locreg, period)
  ensure_dir(pdir)
  ensure_dir(file.path(pdir, "bw_overlays"))
  ensure_dir(file.path(pdir, "kernel_overlays"))
  
  dfp <- cfg$df
  x <- as.numeric(dfp[[cfg$xcol]])
  y <- dfp[[cfg$ycol]]
  
  raw_pts <- tibble(x = x, y = y)
  xlim <- range(x, na.rm = TRUE)
  ylim <- range(y, na.rm = TRUE)
  
  bws <- calc_bandwidths(x, BW_RULES, verbose = TRUE)
  
  date_scale <- if (cfg$x_as_date) {
    scale_x_continuous(
      labels = function(x) format(as.Date(x, origin = "1970-01-01"), "%Y-%m"),
      breaks = scales::pretty_breaks(n = 8)
    )
  } else scale_x_continuous(breaks = scales::pretty_breaks(n = 8))
  
  # Individual regressions
  for (k in KERNELS) {
    for (bw_name in names(bws)) {
      if (is.na(bws[bw_name])) next
      
      d <- smooth_loess(x, y, bws[bw_name])
      if (is.null(d)) next
      
      p <- ggplot() +
        geom_point(data = raw_pts, aes(x = x, y = y),
                   colour = "#95A5A6", size = 2, alpha = 0.3) +
        geom_line(data = d, aes(x = x, y = y), 
                  colour = PALETTE$kernel[k], linewidth = 1.3) +
        labs(
          title = sprintf("%s Local Regression", str_to_title(period)),
          subtitle = sprintf("Kernel: %s | Bandwidth: %s", str_to_title(k), toupper(bw_name)),
          x = if (cfg$x_as_date) "Date" else "X", 
          y = "Value"
        ) +
        theme_enhanced() +
        coord_cartesian(xlim = xlim, ylim = ylim) +
        date_scale
      
      ggsave(file.path(pdir, sprintf("%s_locreg_%s_%s.pdf", period, k, bw_name)),
             p, width = 10, height = 6, device = cairo_pdf)
    }
  }
  
  # Bandwidth overlays
  for (bw_name in names(bws)) {
    if (is.na(bws[bw_name])) next
    
    dlist <- compact(lapply(KERNELS, function(k) {
      d <- smooth_loess(x, y, bws[bw_name])
      if (!is.null(d)) d$kernel <- k
      d
    }))
    
    if (length(dlist) > 0) {
      p <- ggplot(bind_rows(dlist), aes(x = x, y = y, colour = kernel, linetype = kernel)) +
        geom_line(linewidth = 1.2) +
        scale_colour_manual(values = PALETTE$kernel, name = "Kernel") +
        scale_linetype_manual(
          values = c(rectangular = "solid", epanechnikov = "dashed",
                     biweight = "dotted", triangular = "dotdash", gaussian = "longdash"),
          name = "Kernel"
        ) +
        labs(
          title = sprintf("%s Local Regression - Kernel Comparison", str_to_title(period)),
          subtitle = sprintf("Bandwidth: %s", toupper(bw_name)),
          x = if (cfg$x_as_date) "Date" else "X", 
          y = "Value"
        ) +
        theme_enhanced() +
        coord_cartesian(xlim = xlim, ylim = ylim) +
        date_scale
      
      ggsave(file.path(pdir, "bw_overlays", sprintf("bw_%s_overlay.pdf", bw_name)),
             p, width = 11, height = 6, device = cairo_pdf)
    }
  }
  
  # Kernel overlays
  for (k in KERNELS) {
    dlist <- compact(lapply(names(bws), function(bw_name) {
      if (is.na(bws[bw_name])) return(NULL)
      d <- smooth_loess(x, y, bws[bw_name])
      if (!is.null(d)) d$bw <- toupper(bw_name)
      d
    }))
    
    if (length(dlist) > 0) {
      p <- ggplot(bind_rows(dlist), aes(x = x, y = y, colour = bw, linetype = bw)) +
        geom_line(linewidth = 1.2) +
        scale_colour_manual(values = PALETTE$bw, name = "Bandwidth") +
        scale_linetype_manual(
          values = c(NRD0 = "solid", SJ = "dashed", UCV = "dotted"),
          name = "Bandwidth"
        ) +
        labs(
          title = sprintf("%s Local Regression - Bandwidth Comparison", str_to_title(period)),
          subtitle = sprintf("Kernel: %s", str_to_title(k)),
          x = if (cfg$x_as_date) "Date" else "X", 
          y = "Value"
        ) +
        theme_enhanced() +
        coord_cartesian(xlim = xlim, ylim = ylim) +
        date_scale
      
      ggsave(file.path(pdir, "kernel_overlays", sprintf("kernel_%s_overlay.pdf", k)),
             p, width = 11, height = 6, device = cairo_pdf)
    }
  }
}

# ==============================================================================
# QUANTILE REGRESSIONS - COMPREHENSIVE ANALYSIS
# ==============================================================================

message("\n=== GENERATING QUANTILE REGRESSIONS ===\n")

QUANTILES <- c(0.025, 0.25, 0.5, 0.75, 0.975)
quantile_summary <- list()

for (period in names(period_config)) {
  message(sprintf("%s quantile regressions:", str_to_title(period)))
  
  cfg <- period_config[[period]]
  qdir <- file.path(DIRS$quantreg, period)
  
  # Create subdirectories
  subdirs <- c("nonparametric", "parametric", "comparison", "bw_overlays",
               "kernel_overlays", "quantile_overlays", "confidence_bands")
  lapply(file.path(qdir, subdirs), ensure_dir)
  
  dfp <- cfg$df
  x <- as.numeric(dfp[[cfg$xcol]])
  y <- dfp[[cfg$ycol]]
  
  raw_pts <- tibble(x = x, y = y)
  xlim <- range(x, na.rm = TRUE)
  ylim <- range(y, na.rm = TRUE)
  
  bws <- calc_bandwidths(x, BW_RULES, verbose = TRUE)
  
  period_summary <- list(cv_scores = list(), mise_scores = list(), crossing_violations = list())
  
  # --- PARAMETRIC BASELINE ---
  message(sprintf("  [1/6] Parametric baseline (%d quantiles)...", length(QUANTILES)))
  
  parametric_fits <- setNames(
    lapply(seq_along(QUANTILES), function(i) {
      tau <- QUANTILES[i]
      progress_bar(i, length(QUANTILES), label = "    Fitting")
      parametric_quantile_reg(x, y, tau)
    }),
    as.character(QUANTILES)
  ) %>% compact()
  
  # Plot individual parametric
  plot_num <- 0
  for (tau_str in names(parametric_fits)) {
    plot_num <- plot_num + 1
    progress_bar(plot_num, length(parametric_fits), label = "    Plotting")
    
    fit_param <- parametric_fits[[tau_str]]
    tau <- as.numeric(tau_str)
    
    p <- ggplot() +
      geom_point(data = raw_pts, aes(x = x, y = y),
                 colour = "grey60", size = 1.5, alpha = 0.3) +
      geom_line(data = fit_param, aes(x = x, y = y),
                colour = PALETTE$quantile[tau_str], linewidth = 1.2) +
      labs(
        title = sprintf("%s Parametric Quantile Regression", str_to_title(period)),
        subtitle = sprintf("Quantile: %.3f | Method: Polynomial (degree 2)", tau),
        x = if (cfg$x_as_date) "Date" else "X",
        y = "Value"
      ) +
      theme_enhanced() +
      coord_cartesian(xlim = xlim, ylim = ylim)
    
    ggsave(file.path(qdir, "parametric", sprintf("param_tau%.3f.pdf", tau)),
           p, width = 10, height = 6, device = cairo_pdf)
  }
  
  # Parametric: All quantiles overlay
  if (length(parametric_fits) > 0) {
    param_all <- bind_rows(lapply(names(parametric_fits), function(tau) {
      parametric_fits[[tau]] %>% mutate(tau = tau)
    }))
    
    p <- ggplot() +
      geom_point(data = raw_pts, aes(x = x, y = y),
                 colour = "grey70", size = 1.5, alpha = 0.2) +
      geom_line(data = param_all, aes(x = x, y = y, colour = tau), linewidth = 1.1) +
      scale_colour_manual(values = PALETTE$quantile, name = "Quantile") +
      labs(
        title = sprintf("%s Parametric Quantile Regression - All Quantiles", str_to_title(period)),
        subtitle = "Method: Polynomial (degree 2)",
        x = if (cfg$x_as_date) "Date" else "X",
        y = "Value"
      ) +
      theme_enhanced() +
      coord_cartesian(xlim = xlim, ylim = ylim)
    
    ggsave(file.path(qdir, "parametric", "parametric_all_quantiles.pdf"),
           p, width = 10, height = 6, device = cairo_pdf)
  }
  
  # --- NONPARAMETRIC QUANTILE REGRESSION ---
  message(sprintf("  [2/6] Nonparametric quantile regressions (%d kernels × %d bandwidths × %d quantiles = %d fits)...", 
                  length(KERNELS), sum(!is.na(bws)), length(QUANTILES),
                  length(KERNELS) * sum(!is.na(bws)) * length(QUANTILES)))
  
  total_fits <- length(KERNELS) * sum(!is.na(bws)) * length(QUANTILES)
  current_fit <- 0
  
  for (k in KERNELS) {
    for (bw_name in names(bws)) {
      if (is.na(bws[bw_name])) next
      
      combo_fits <- setNames(
        lapply(QUANTILES, function(tau) {
          current_fit <<- current_fit + 1
          progress_bar(current_fit, total_fits, label = "    Fitting")
          local_quantile_reg(x, y, tau, k, bws[bw_name])
        }),
        as.character(QUANTILES)
      ) %>% compact()
      
      # Calculate MISE scores
      for (tau_str in names(combo_fits)) {
        if (tau_str %in% names(parametric_fits)) {
          mise <- calc_mise(combo_fits[[tau_str]], parametric_fits[[tau_str]])
          period_summary$mise_scores[[sprintf("%s_%s_%s", k, bw_name, tau_str)]] <- mise
        }
      }
      
      # Individual nonparametric plots
      for (tau_str in names(combo_fits)) {
        tau <- as.numeric(tau_str)
        
        p <- ggplot() +
          geom_point(data = raw_pts, aes(x = x, y = y),
                     colour = "grey60", size = 1.5, alpha = 0.3) +
          geom_line(data = combo_fits[[tau_str]], aes(x = x, y = y),
                    colour = PALETTE$quantile[tau_str], linewidth = 1.2) +
          labs(
            title = sprintf("%s Nonparametric Quantile Regression", str_to_title(period)),
            subtitle = sprintf("Quantile: %.3f | Kernel: %s | Bandwidth: %s", 
                               tau, str_to_title(k), toupper(bw_name)),
            x = if (cfg$x_as_date) "Date" else "X",
            y = "Value"
          ) +
          theme_enhanced() +
          coord_cartesian(xlim = xlim, ylim = ylim)
        
        ggsave(file.path(qdir, "nonparametric",
                         sprintf("nonparam_%s_%s_tau%.3f.pdf", k, bw_name, tau)),
               p, width = 10, height = 6, device = cairo_pdf)
      }
      
      # Check quantile crossings
      if (length(combo_fits) > 1) {
        crossings <- check_quantile_crossing(combo_fits)
        if (length(crossings) > 0) {
          period_summary$crossing_violations[[sprintf("%s_%s", k, bw_name)]] <- crossings
        }
      }
      
      # All quantiles overlay
      if (length(combo_fits) > 0) {
        combo_all <- bind_rows(lapply(names(combo_fits), function(tau) {
          combo_fits[[tau]] %>% mutate(tau = tau)
        }))
        
        p <- ggplot() +
          geom_point(data = raw_pts, aes(x = x, y = y),
                     colour = "grey70", size = 1.5, alpha = 0.2) +
          geom_line(data = combo_all, aes(x = x, y = y, colour = tau), linewidth = 1.1) +
          scale_colour_manual(values = PALETTE$quantile, name = "Quantile") +
          labs(
            title = sprintf("%s Nonparametric Quantile Regression", str_to_title(period)),
            subtitle = sprintf("Kernel: %s | Bandwidth: %s", str_to_title(k), toupper(bw_name)),
            x = if (cfg$x_as_date) "Date" else "X",
            y = "Value"
          ) +
          theme_enhanced() +
          coord_cartesian(xlim = xlim, ylim = ylim)
        
        ggsave(file.path(qdir, "quantile_overlays", sprintf("quantiles_%s_%s.pdf", k, bw_name)),
               p, width = 10, height = 6, device = cairo_pdf)
      }
    }
  }
  
  # --- COMPARISON PLOTS ---
  message(sprintf("  [3/6] Parametric vs nonparametric comparisons (%d plots)...", 
                  length(parametric_fits)))
  
  plot_num <- 0
  for (tau_str in names(parametric_fits)) {
    plot_num <- plot_num + 1
    progress_bar(plot_num, length(parametric_fits), label = "    Plotting")
    
    tau <- as.numeric(tau_str)
    
    if ("nrd0" %in% names(bws) && !is.na(bws["nrd0"])) {
      fit_np <- local_quantile_reg(x, y, tau, "gaussian", bws["nrd0"])
      
      if (!is.null(fit_np)) {
        comparison_data <- bind_rows(
          parametric_fits[[tau_str]] %>% mutate(method = "Parametric"),
          fit_np %>% mutate(method = "Nonparametric")
        )
        
        p <- ggplot() +
          geom_point(data = raw_pts, aes(x = x, y = y),
                     colour = "grey70", size = 1.5, alpha = 0.2) +
          geom_line(data = comparison_data, 
                    aes(x = x, y = y, colour = method, linetype = method), linewidth = 1.1) +
          scale_colour_manual(values = c(Parametric = "#E74C3C", Nonparametric = "#3498DB"),
                              name = "Method") +
          scale_linetype_manual(values = c(Parametric = "dashed", Nonparametric = "solid"),
                                name = "Method") +
          labs(
            title = sprintf("%s Quantile Regression Comparison", str_to_title(period)),
            subtitle = sprintf("Quantile: %.3f | Nonparametric: Gaussian kernel, NRD0 bandwidth", tau),
            x = if (cfg$x_as_date) "Date" else "X",
            y = "Value"
          ) +
          theme_enhanced() +
          coord_cartesian(xlim = xlim, ylim = ylim)
        
        ggsave(file.path(qdir, "comparison", sprintf("comparison_tau%.3f.pdf", tau)),
               p, width = 10, height = 6, device = cairo_pdf)
      }
    }
  }
  
  # --- CONFIDENCE BANDS (median only) ---
  ci_kernels <- c("gaussian", "epanechnikov")
  valid_bws <- names(bws)[!is.na(bws)]
  total_ci <- length(ci_kernels) * length(valid_bws)
  
  message(sprintf("  [4/6] Bootstrap confidence bands (%d combinations, n_boot=50)...", total_ci))
  
  ci_num <- 0
  for (k in ci_kernels) {
    for (bw_name in valid_bws) {
      ci_num <- ci_num + 1
      progress_bar(ci_num, total_ci, label = "    Computing")
      
      if (is.na(bws[bw_name])) next
      
      fit_np <- local_quantile_reg(x, y, 0.5, k, bws[bw_name])
      if (is.null(fit_np)) next
      
      ci_bands <- bootstrap_quantile_ci(x, y, 0.5, k, bws[bw_name], n_boot = 50)
      
      if (!is.null(ci_bands)) {
        p <- ggplot() +
          geom_ribbon(data = ci_bands, aes(x = x, ymin = lower, ymax = upper),
                      fill = "#3498DB", alpha = 0.2) +
          geom_point(data = raw_pts, aes(x = x, y = y),
                     colour = "grey60", size = 1.5, alpha = 0.3) +
          geom_line(data = fit_np, aes(x = x, y = y),
                    colour = "#3498DB", linewidth = 1.3) +
          labs(
            title = sprintf("%s Quantile Regression with 95%% Confidence Bands", str_to_title(period)),
            subtitle = sprintf("Median | Kernel: %s | Bandwidth: %s | Bootstrap n=50", 
                               str_to_title(k), toupper(bw_name)),
            x = if (cfg$x_as_date) "Date" else "X",
            y = "Value"
          ) +
          theme_enhanced() +
          coord_cartesian(xlim = xlim, ylim = ylim)
        
        ggsave(file.path(qdir, "confidence_bands", sprintf("ci_%s_%s.pdf", k, bw_name)),
               p, width = 10, height = 6, device = cairo_pdf)
      }
    }
  }
  
  # --- BANDWIDTH OVERLAYS ---
  total_bw_overlays <- length(KERNELS) * length(QUANTILES)
  message(sprintf("  [5/6] Bandwidth comparison overlays (%d plots)...", total_bw_overlays))
  
  overlay_num <- 0
  for (k in KERNELS) {
    for (tau in QUANTILES) {
      overlay_num <- overlay_num + 1
      progress_bar(overlay_num, total_bw_overlays, label = "    Generating")
      
      dlist <- compact(lapply(names(bws), function(bw_name) {
        if (is.na(bws[bw_name])) return(NULL)
        d <- local_quantile_reg(x, y, tau, k, bws[bw_name])
        if (!is.null(d)) d$bw <- toupper(bw_name)
        d
      }))
      
      if (length(dlist) > 0) {
        p <- ggplot(bind_rows(dlist), aes(x = x, y = y, colour = bw, linetype = bw)) +
          geom_line(linewidth = 1.1) +
          scale_colour_manual(values = PALETTE$bw, name = "Bandwidth") +
          labs(
            title = sprintf("%s Bandwidth Comparison", str_to_title(period)),
            subtitle = sprintf("Quantile: %.3f | Kernel: %s", tau, str_to_title(k)),
            x = if (cfg$x_as_date) "Date" else "X",
            y = "Value"
          ) +
          theme_enhanced() +
          coord_cartesian(xlim = xlim, ylim = ylim)
        
        ggsave(file.path(qdir, "bw_overlays", sprintf("bw_tau%.3f_%s.pdf", tau, k)),
               p, width = 10, height = 6, device = cairo_pdf)
      }
    }
  }
  
  # --- KERNEL OVERLAYS ---
  total_kernel_overlays <- length(valid_bws) * length(QUANTILES)
  message(sprintf("  [6/6] Kernel comparison overlays (%d plots)...", total_kernel_overlays))
  
  overlay_num <- 0
  for (bw_name in valid_bws) {
    if (is.na(bws[bw_name])) next
    
    for (tau in QUANTILES) {
      overlay_num <- overlay_num + 1
      progress_bar(overlay_num, total_kernel_overlays, label = "    Generating")
      
      dlist <- compact(lapply(KERNELS, function(k) {
        d <- local_quantile_reg(x, y, tau, k, bws[bw_name])
        if (!is.null(d)) d$kernel <- k
        d
      }))
      
      if (length(dlist) > 0) {
        p <- ggplot(bind_rows(dlist), aes(x = x, y = y, colour = kernel, linetype = kernel)) +
          geom_line(linewidth = 1.1) +
          scale_colour_manual(values = PALETTE$kernel, name = "Kernel") +
          labs(
            title = sprintf("%s Kernel Comparison", str_to_title(period)),
            subtitle = sprintf("Quantile: %.3f | Bandwidth: %s", tau, toupper(bw_name)),
            x = if (cfg$x_as_date) "Date" else "X",
            y = "Value"
          ) +
          theme_enhanced() +
          coord_cartesian(xlim = xlim, ylim = ylim)
        
        ggsave(file.path(qdir, "kernel_overlays", sprintf("kernel_tau%.3f_%s.pdf", tau, bw_name)),
               p, width = 10, height = 6, device = cairo_pdf)
      }
    }
  }
  
  quantile_summary[[period]] <- period_summary
}

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

message("\n=== GENERATING SUMMARY STATISTICS ===\n")

summary_dir <- file.path(DIRS$quantreg, "summary_statistics")
ensure_dir(summary_dir)

# Bandwidth summary with detailed diagnostics
bandwidth_summary <- bind_rows(lapply(names(period_config), function(period) {
  cfg <- period_config[[period]]
  dfp <- cfg$df
  x <- as.numeric(dfp[[cfg$xcol]])
  
  range_x <- diff(range(x, na.rm = TRUE))
  n <- length(x)
  
  # Calculate all three bandwidths
  nrd0_bw <- get_bw("nrd0", x, verbose = FALSE)
  sj_bw <- get_bw("SJ", x, verbose = FALSE)
  ucv_bw <- get_bw("ucv", x, verbose = FALSE)
  
  tibble(
    period = period,
    n_points = n,
    data_range = range_x,
    nrd0 = nrd0_bw,
    SJ = sj_bw,
    UCV = ucv_bw,
    ucv_vs_nrd0_ratio = ucv_bw / nrd0_bw,
    ucv_pct_of_range = 100 * ucv_bw / range_x,
    sj_vs_nrd0_ratio = sj_bw / nrd0_bw
  )
}))

write.csv(bandwidth_summary, file.path(summary_dir, "bandwidth_summary.csv"), row.names = FALSE)

# Print bandwidth summary to console
message("\n=== BANDWIDTH SUMMARY ===")
print(bandwidth_summary, width = Inf)
message("=========================\n")

# MISE Scores
mise_df <- bind_rows(lapply(names(quantile_summary), function(period) {
  mise_list <- quantile_summary[[period]]$mise_scores
  if (length(mise_list) == 0) return(NULL)
  data.frame(period = period, config = names(mise_list), 
             mise = unlist(mise_list), stringsAsFactors = FALSE)
}))

if (!is.null(mise_df) && nrow(mise_df) > 0) {
  write.csv(mise_df, file.path(summary_dir, "mise_scores.csv"), row.names = FALSE)
  
  mise_df %>%
    group_by(period) %>%
    slice_min(mise, n = 3) %>%
    write.csv(file.path(summary_dir, "best_mise_configs.csv"), row.names = FALSE)
}

# Crossing Violations
crossing_df <- bind_rows(lapply(names(quantile_summary), function(period) {
  cv_list <- quantile_summary[[period]]$crossing_violations
  if (length(cv_list) == 0) return(NULL)
  
  bind_rows(lapply(names(cv_list), function(config) {
    data.frame(
      period = period, config = config,
      comparison = names(cv_list[[config]]),
      n_crossings = unlist(cv_list[[config]]),
      stringsAsFactors = FALSE
    )
  }))
}))

if (!is.null(crossing_df) && nrow(crossing_df) > 0) {
  write.csv(crossing_df, file.path(summary_dir, "quantile_crossings.csv"), row.names = FALSE)
}

message("\n=== ANALYSIS COMPLETE ===")
message("All outputs generated with:")
message("  ✓ Corrected UCV: R's bw.ucv() (primary)")
message("  ✓ Manual UCV fallback (when native fails)")
message("  ✓ Global fallback across all analysis types")
message("  ✓ Smart constraints (0.5-3× nrd0, 5-25% range)")
message("  ✓ Enhanced diagnostics and verbose output")
message("\nCheck quantile_regressions/summary_statistics/ for detailed metrics.")
message("==========================\n")