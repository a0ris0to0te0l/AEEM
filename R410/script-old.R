################################################################################
# STATE-DEPENDENT GOVERNMENT SPENDING MULTIPLIERS
# A Regional Study for the Western Balkans / South-Eastern Europe
#
# Countries: Albania (AL), Bosnia & Herzegovina (BA), Bulgaria (BG),
#            Croatia (HR), Greece (EL), Kosovo (XK), Montenegro (ME),
#            North Macedonia (MK), Romania (RO), Serbia (RS), Slovenia (SI)
#
# Methodology: Owyang, Ramey & Zubairy (2013) — State-Dependent Local Projections
#              Ramey & Zubairy (2018)          — IV Cumulative Multiplier Scaling
#              Blanchard & Perotti (2002)      — Spending shock identification
#
# Data: Eurostat API ({eurostat} R package)
#       Countries with insufficient Eurostat coverage are excluded automatically.
# Author: [Your name]
# Date:   Sys.Date()
################################################################################

set.seed(20240601)   # Reproducibility

# ==============================================================================
# 00. WORKING DIRECTORY
# ==============================================================================

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  cat("Working directory:", getwd(), "\n\n")
} else {
  cat("Not in RStudio. Set your working directory manually if needed.\n\n")
}

# ==============================================================================
# 01. PACKAGES
# ==============================================================================

required_pkgs <- c(
  # Data retrieval
  "eurostat",
  # Wrangling
  "dplyr", "tidyr", "lubridate", "zoo",
  # Econometrics
  "mFilter", "AER", "car", "lmtest", "sandwich", "tseries",
  # Parallel / progress
  "pbapply", "parallel",
  # Visualisation
  "ggplot2", "patchwork", "scales", "ggrepel", "viridis"
)

new_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(new_pkgs)) {
  message("Installing missing packages: ", paste(new_pkgs, collapse = ", "))
  install.packages(new_pkgs, repos = "https://cloud.r-project.org", quiet = TRUE)
}

invisible(lapply(required_pkgs, library, character.only = TRUE))
options(eurostat_verbose = FALSE, scipen = 999)

# ==============================================================================
# 02. COUNTRY CONFIGURATION
# ==============================================================================

# Eurostat codes and display labels for the full Balkan region.
# Countries with fewer than MIN_OBS usable quarters are automatically excluded
# from estimation (see Section 06). No external fallback is used.

COUNTRIES <- c(
  "AL",  # Albania
  "BA",  # Bosnia & Herzegovina
  "BG",  # Bulgaria
  "EL",  # Greece
  "HR",  # Croatia
  "ME",  # Montenegro
  "MK",  # North Macedonia
  "RO",  # Romania
  "RS",  # Serbia
  "SI",  # Slovenia
  "XK"   # Kosovo
)

COUNTRY_LABELS <- c(
  AL = "Albania",
  BA = "Bosnia & Herz.",
  BG = "Bulgaria",
  EL = "Greece",
  HR = "Croatia",
  ME = "Montenegro",
  MK = "North Macedonia",
  RO = "Romania",
  RS = "Serbia",
  SI = "Slovenia",
  XK = "Kosovo"
)

# National currency multipliers differ; for this study we use
# chain-linked volumes (CLV10_MNAC) so cross-country comparisons
# are in real local-currency units — multipliers are unit-free ratios.

H          <- 12   # Maximum LP horizon (quarters)
HP_FREQ    <- 1600 # HP filter smoothing parameter (quarterly data)
CI90_Z     <- 1.645
CI68_Z     <- 1.000
NW_PREWHITE <- FALSE

# ==============================================================================
# 03. DATA DOWNLOAD
# ==============================================================================

safe_eurostat <- function(id, filters, ...) {
  tryCatch(
    get_eurostat(id, filters = filters, time_format = "date", cache = TRUE, ...),
    error = function(e) {
      warning(sprintf("Eurostat download failed for '%s': %s", id, e$message))
      NULL
    }
  )
}

MIN_OBS <- 24  # Need at least 6 years of quarters for reliable LP estimation

cat("Downloading Eurostat GDP and government expenditure data...\n")
gdp_raw <- safe_eurostat(
  id      = "namq_10_gdp",
  filters = list(
    geo     = COUNTRIES,
    na_item = c("B1GQ", "P3_S13"),  # GDP and Government final consumption
    unit    = "CLV10_MNAC",          # Chain-linked volumes, 2010 prices, national currency
    s_adj   = "SCA"                  # Seasonally & calendar adjusted
  )
)

cat("Downloading Eurostat unemployment data...\n")
unemp_raw <- safe_eurostat(
  id      = "lfsq_urgan",
  filters = list(
    geo     = COUNTRIES,
    sex     = "T",
    unit    = "PC",
    age     = "Y15-74",
    citizen = "TOTAL"
  )
)


# ==============================================================================
# 04. DATA CLEANING
# ==============================================================================

cat("\nCleaning Eurostat data...\n")

clean_eurostat <- function(gdp_raw, unemp_raw) {
  gdp_clean <- gdp_raw %>%
    select(geo, time, na_item, values) %>%
    filter(!is.na(values)) %>%
    pivot_wider(names_from = na_item, values_from = values) %>%
    rename(country = geo, date = time, gdp = B1GQ, gov_exp = P3_S13) %>%
    arrange(country, date)
  
  unemp_clean <- unemp_raw %>%
    select(geo, time, values) %>%
    filter(!is.na(values)) %>%
    rename(country = geo, date = time, unemp = values) %>%
    arrange(country, date)
  
  left_join(gdp_clean, unemp_clean, by = c("country", "date")) %>%
    filter(!is.na(gdp), !is.na(gov_exp))
}

panel_raw <- if (!is.null(gdp_raw) && !is.null(unemp_raw)) {
  clean_eurostat(gdp_raw, unemp_raw)
} else {
  stop("Eurostat download returned NULL for GDP or unemployment. Check your connection and try again.")
}

cat(sprintf("  Eurostat panel: %d observations across %d countries.\n",
            nrow(panel_raw), length(unique(panel_raw$country))))

# ==============================================================================
# 05. VARIABLE CONSTRUCTION
# ==============================================================================

cat("\nConstructing variables...\n")

safe_hp <- function(x, freq = HP_FREQ, min_obs = 8) {
  idx    <- which(!is.na(x))
  result <- rep(NA_real_, length(x))
  if (length(idx) >= min_obs) {
    hp           <- hpfilter(x[idx], freq = freq)
    result[idx]  <- hp$trend
  }
  result
}

ldiff <- function(x) c(NA, diff(log(x))) * 100

panel <- panel_raw %>%
  group_by(country) %>%
  arrange(date) %>%
  mutate(
    # Log levels
    ln_gdp = log(gdp),
    ln_g   = log(gov_exp),
    
    # Log differences (growth rates)
    dy = ldiff(gdp),
    dg = ldiff(gov_exp),
    
    # Ramey-Zubairy (2018) scaling: levels relative to lagged GDP
    gdp_lag1  = lag(gdp, 1),
    g_scaled  = (gov_exp - lag(gov_exp)) / gdp_lag1 * 100,
    y_scaled  = (gdp    - lag(gdp))     / gdp_lag1 * 100,
    
    # HP-filtered unemployment gap
    unemp_hp_trend = safe_hp(unemp),
    unemp_gap      = unemp - unemp_hp_trend,
    
    # State indicator: 1 = slack, 0 = expansion (lagged 1 quarter)
    slack      = as.integer(unemp_gap > 0),
    slack_lag1 = lag(slack, 1),
    
    # Control variable lags
    dy_lag1       = lag(dy, 1),
    dy_lag2       = lag(dy, 2),
    dg_lag1       = lag(dg, 1),
    dg_lag2       = lag(dg, 2),
    g_scaled_lag1 = lag(g_scaled, 1),
    g_scaled_lag2 = lag(g_scaled, 2),
    y_scaled_lag1 = lag(y_scaled, 1),
    y_scaled_lag2 = lag(y_scaled, 2)
  ) %>%
  ungroup()

cat(sprintf("  Variables constructed. Panel rows: %d\n", nrow(panel)))

# ==============================================================================
# 06. DATA QUALITY CHECK
# ==============================================================================

cat("\nData coverage by country:\n")
coverage <- panel %>%
  filter(!is.na(g_scaled), !is.na(y_scaled), !is.na(slack_lag1)) %>%
  group_by(country) %>%
  summarise(
    start_q   = min(date),
    end_q     = max(date),
    n_obs     = n(),
    n_slack   = sum(slack_lag1, na.rm = TRUE),
    pct_slack = round(mean(slack_lag1, na.rm = TRUE) * 100, 1),
    .groups   = "drop"
  ) %>%
  mutate(country_name = COUNTRY_LABELS[country])

print(coverage)

# Drop countries with too few observations to estimate LP reliably
USABLE_COUNTRIES <- coverage %>%
  filter(n_obs >= MIN_OBS, n_slack >= 6, n_obs - n_slack >= 6) %>%
  pull(country)

cat(sprintf("\nUsable countries (>= %d obs, >= 6 quarters in each state):\n  %s\n",
            MIN_OBS, paste(USABLE_COUNTRIES, collapse = ", ")))

panel_est <- panel %>%
  filter(
    country %in% USABLE_COUNTRIES,
    !is.na(g_scaled), !is.na(y_scaled),
    !is.na(slack_lag1),
    !is.na(g_scaled_lag1), !is.na(g_scaled_lag2),
    !is.na(y_scaled_lag1), !is.na(y_scaled_lag2)
  )

# ==============================================================================
# 07. UNIT ROOT TESTS (ADF)
# ==============================================================================

cat("\n--- Augmented Dickey-Fuller tests ---\n")

run_adf <- function(x, name) {
  x <- x[!is.na(x)]
  if (length(x) < 10) return(NULL)
  tryCatch({
    t <- adf.test(x)
    data.frame(
      variable   = name,
      adf_stat   = round(t$statistic, 3),
      p_value    = round(t$p.value, 3),
      stationary = t$p.value < 0.05
    )
  }, error = function(e) NULL)
}

adf_table <- panel_est %>%
  group_by(country) %>%
  group_modify(function(d, k) {
    bind_rows(
      run_adf(d$ln_gdp, "ln_gdp"),
      run_adf(d$ln_g,   "ln_g"),
      run_adf(d$dy,     "dy"),
      run_adf(d$dg,     "dg")
    )
  }) %>%
  ungroup()

print(adf_table)

# ==============================================================================
# 08. LOCAL PROJECTION IV ESTIMATOR
# ==============================================================================
#
# Identification: Blanchard-Perotti (2002) — government spending is assumed
# predetermined within the quarter, so the contemporaneous spending change
# (g_scaled) is a valid instrument for cumulative spending at horizon h.
#
# Estimator: 2SLS via ivreg() with:
#   - Dependent variable:  Σ_{s=0}^{h} y_{t+s} / GDP_{t-1}  (cum output)
#   - Endogenous variable: Σ_{s=0}^{h} g_{t+s} / GDP_{t-1}  (cum spending)
#   - Instrument:          g_t / GDP_{t-1}  (contemporaneous shock)
# Interacted with state dummies (Slack / Expansion) à la ORZ (2013).
#
# Standard errors: Newey-West HAC with bandwidth = h+1 (Ramey & Zubairy 2018).
# Robustness: clustered by country (for panel specs).

run_lp_iv <- function(data, h, cluster_se = FALSE) {
  
  CTRL <- c("y_scaled_lag1", "y_scaled_lag2", "g_scaled_lag1", "g_scaled_lag2")
  
  df <- data %>%
    group_by(country) %>%
    arrange(date) %>%
    mutate(
      cum_y_h = rollapply(y_scaled, width = h + 1, FUN = sum,
                          align = "left", fill = NA, partial = FALSE),
      cum_g_h = rollapply(g_scaled, width = h + 1, FUN = sum,
                          align = "left", fill = NA, partial = FALSE)
    ) %>%
    ungroup() %>%
    filter(!is.na(cum_y_h), !is.na(cum_g_h), !is.na(slack_lag1)) %>%
    mutate(
      state_slack = slack_lag1,
      state_exp   = 1L - slack_lag1,
      
      # Interacted endogenous regressors
      cum_g_slack = state_slack * cum_g_h,
      cum_g_exp   = state_exp   * cum_g_h,
      
      # Interacted instruments (Blanchard-Perotti: contemporaneous shock)
      shock_slack = state_slack * g_scaled,
      shock_exp   = state_exp   * g_scaled,
      
      country = factor(country)
    ) %>%
    filter(complete.cases(.[, CTRL]))
  
  if (nrow(df) < 20) return(NULL)
  
  n_countries <- length(unique(df$country))
  ctrl_str    <- paste(CTRL, collapse = " + ")
  fe_str      <- if (n_countries > 1) " + country - 1" else ""
  
  fml <- as.formula(paste0(
    "cum_y_h ~ cum_g_slack + cum_g_exp + ", ctrl_str, fe_str, " | ",
    "shock_slack + shock_exp + ", ctrl_str, fe_str
  ))
  
  tryCatch({
    fit <- ivreg(fml, data = df)
    
    # --- Standard errors ---
    if (cluster_se && n_countries > 1) {
      # Cluster by country (preferred for panel)
      vcov_use <- vcovCL(fit, cluster = ~country)
    } else {
      # Newey-West HAC (preferred for time series / single country)
      vcov_use <- NeweyWest(fit, lag = h + 1, prewhite = NW_PREWHITE)
    }
    
    ct   <- coeftest(fit, vcov = vcov_use)
    wald <- car::linearHypothesis(fit, "cum_g_slack = cum_g_exp",
                                  vcov. = vcov_use, test = "F")
    
    # --- First-stage diagnostics (weak instrument check) ---
    # F-stat from first stage for the slack interaction
    fs_slack <- tryCatch({
      fs_fit <- lm(as.formula(paste0(
        "cum_g_slack ~ shock_slack + shock_exp + ", ctrl_str, fe_str
      )), data = df)
      summary(fs_fit)$fstatistic[1]
    }, error = function(e) NA_real_)
    
    # coeftest column name is "Pr(>|t|)" for ivreg objects
    pval_col <- grep("Pr\\(", colnames(ct), value = TRUE)[1]
    
    list(
      h             = h,
      beta_slack    = ct["cum_g_slack", "Estimate"],
      beta_exp      = ct["cum_g_exp",   "Estimate"],
      se_slack      = ct["cum_g_slack", "Std. Error"],
      se_exp        = ct["cum_g_exp",   "Std. Error"],
      pval_slack    = ct["cum_g_slack", pval_col],
      pval_exp      = ct["cum_g_exp",   pval_col],
      pval_diff     = wald$`Pr(>F)`[2],
      fs_f_stat     = fs_slack,
      n             = nrow(df),
      n_countries   = n_countries,
      df_resid      = fit$df.residual
    )
  }, error = function(e) {
    warning(sprintf("LP failed at h=%d: %s", h, conditionMessage(e)))
    NULL
  })
}

# ==============================================================================
# 09. ESTIMATE POOLED PANEL AND COUNTRY-BY-COUNTRY LPs
# ==============================================================================

n_cores <- max(1L, detectCores() - 1L)
cat(sprintf("\nEstimating LP-IV across horizons 0:%d using %d cores...\n", H, n_cores))

# — Pooled panel (all usable countries, country FE, clustered SE)
cat("  [1/2] Pooled panel...\n")
lp_panel <- pblapply(0:H, function(h) run_lp_iv(panel_est, h, cluster_se = TRUE))
lp_panel <- Filter(Negate(is.null), lp_panel)

if (length(lp_panel) == 0) {
  stop(paste(
    "LP-IV produced no results for the pooled panel.",
    "Check that panel_est has enough complete observations (run View(panel_est))",
    "and that both slack states are represented.",
    sep = "\n"
  ))
}

# — Country-by-country (Newey-West SE)
cat("  [2/2] Country-by-country...\n")
lp_by_country <- lapply(USABLE_COUNTRIES, function(cc) {
  sub <- panel_est %>% filter(country == cc)
  res <- pblapply(0:H, function(h) {
    r <- run_lp_iv(sub, h, cluster_se = FALSE)
    if (!is.null(r)) r$country <- cc
    r
  })
  Filter(Negate(is.null), res)
})
names(lp_by_country) <- USABLE_COUNTRIES

# ==============================================================================
# 10. EXTRACT IRF DATA FRAMES
# ==============================================================================

extract_irf <- function(lp_list, label) {
  if (is.null(lp_list) || length(lp_list) == 0) {
    warning(sprintf("extract_irf: no LP results for '%s' — skipping.", label))
    return(NULL)
  }
  do.call(rbind, lapply(lp_list, function(r) {
    data.frame(
      group          = label,
      horizon        = r$h,
      beta_slack     = r$beta_slack,
      se_slack       = r$se_slack,
      beta_exp       = r$beta_exp,
      se_exp         = r$se_exp,
      pval_diff      = r$pval_diff,
      fs_f_stat      = r$fs_f_stat,
      n              = r$n,
      ci90_lo_slack  = r$beta_slack - CI90_Z * r$se_slack,
      ci90_hi_slack  = r$beta_slack + CI90_Z * r$se_slack,
      ci90_lo_exp    = r$beta_exp   - CI90_Z * r$se_exp,
      ci90_hi_exp    = r$beta_exp   + CI90_Z * r$se_exp,
      ci68_lo_slack  = r$beta_slack - CI68_Z * r$se_slack,
      ci68_hi_slack  = r$beta_slack + CI68_Z * r$se_slack,
      ci68_lo_exp    = r$beta_exp   - CI68_Z * r$se_exp,
      ci68_hi_exp    = r$beta_exp   + CI68_Z * r$se_exp
    )
  }))
}

irf_panel <- extract_irf(lp_panel, "Panel")
if (is.null(irf_panel)) stop("irf_panel is NULL — check LP estimation output above.")

irf_countries <- do.call(rbind, Filter(Negate(is.null), lapply(USABLE_COUNTRIES, function(cc) {
  extract_irf(lp_by_country[[cc]], cc)
})))

irf_all <- rbind(irf_panel, irf_countries)
if (is.null(irf_all) || nrow(irf_all) == 0) stop("irf_all is empty — no IRF results to work with.")

# Multiplier table: extract point estimates at chosen horizons
make_mult_table <- function(irf_df, horizons = c(4, 8, 12)) {
  if (is.null(irf_df) || nrow(irf_df) == 0) stop("make_mult_table: irf_df is NULL or empty.")
  do.call(rbind, lapply(unique(irf_df$group), function(grp) {
    sub <- irf_df[irf_df$group == grp, ]
    rows <- lapply(horizons, function(hh) {
      r <- sub[sub$horizon == hh, ]
      if (nrow(r) == 0) return(NULL)
      data.frame(
        group      = grp,
        horizon    = hh,
        mult_slack = round(r$beta_slack, 3),
        mult_exp   = round(r$beta_exp,   3),
        diff       = round(r$beta_slack - r$beta_exp, 3),
        pval_diff  = round(r$pval_diff,  3)
      )
    })
    do.call(rbind, Filter(Negate(is.null), rows))
  }))
}

multipliers <- make_mult_table(irf_all)

cat("\n=== CUMULATIVE MULTIPLIERS (Ramey-Zubairy Scaling) ===\n")
print(multipliers)

# ==============================================================================
# 11. FIRST STAGE DIAGNOSTICS
# ==============================================================================

cat("\n--- First-stage F-statistics (weak instrument check; rule of thumb: F > 10) ---\n")

fs_diag <- irf_panel %>%
  select(horizon, fs_f_stat) %>%
  mutate(
    fs_f_stat = round(fs_f_stat, 2),
    weak_iv   = fs_f_stat < 10
  )
print(fs_diag)

# ==============================================================================
# 12. HYPOTHESIS TESTS
# ==============================================================================

cat("\n--- Wald tests: H0: β_slack = β_exp at each horizon (panel) ---\n")

wald_table <- irf_panel %>%
  transmute(
    horizon    = horizon,
    beta_slack = round(beta_slack, 3),
    beta_exp   = round(beta_exp, 3),
    difference = round(beta_slack - beta_exp, 3),
    p_value    = round(pval_diff, 3),
    sig        = case_when(
      pval_diff < 0.01 ~ "***",
      pval_diff < 0.05 ~ "**",
      pval_diff < 0.10 ~ "*",
      TRUE             ~ ""
    )
  )

print(wald_table)

# ==============================================================================
# 13. ROBUSTNESS — ALTERNATIVE SLACK THRESHOLDS
# ==============================================================================

cat("\n--- Robustness: percentile-based slack thresholds (25th, 50th, 75th) ---\n")

panel_rob <- panel_est %>%
  group_by(country) %>%
  arrange(date) %>%
  mutate(
    slack_p25 = lag(as.integer(unemp > quantile(unemp, .25, na.rm = TRUE)), 1),
    slack_p50 = lag(as.integer(unemp > quantile(unemp, .50, na.rm = TRUE)), 1),
    slack_p75 = lag(as.integer(unemp > quantile(unemp, .75, na.rm = TRUE)), 1)
  ) %>%
  ungroup()

run_rob_spec <- function(data, slack_col, label) {
  df <- data %>% mutate(slack_lag1 = .data[[slack_col]])
  lp <- pblapply(0:H, function(h) run_lp_iv(df, h, cluster_se = TRUE))
  lp <- Filter(Negate(is.null), lp)
  extract_irf(lp, label)
}

irf_rob <- bind_rows(
  run_rob_spec(panel_rob, "slack_p25", "Rob: u > P25"),
  run_rob_spec(panel_rob, "slack_p50", "Rob: u > P50"),
  run_rob_spec(panel_rob, "slack_p75", "Rob: u > P75")
)

mult_rob <- make_mult_table(irf_rob)
cat("\n=== ROBUSTNESS MULTIPLIERS ===\n")
print(mult_rob)

# ==============================================================================
# 14. SUMMARY STATISTICS
# ==============================================================================

sumstats <- panel_est %>%
  filter(!is.na(slack_lag1)) %>%
  mutate(state = ifelse(slack_lag1 == 1, "Slack", "Expansion")) %>%
  group_by(country = COUNTRY_LABELS[country], state) %>%
  summarise(
    quarters      = n(),
    mean_gdp_gr   = round(mean(dy,    na.rm = TRUE), 2),
    mean_g_gr     = round(mean(dg,    na.rm = TRUE), 2),
    mean_unemp    = round(mean(unemp, na.rm = TRUE), 2),
    sd_gdp_gr     = round(sd(dy,      na.rm = TRUE), 2),
    sd_g_gr       = round(sd(dg,      na.rm = TRUE), 2),
    .groups = "drop"
  )

cat("\n=== SUMMARY STATISTICS BY COUNTRY AND STATE ===\n")
print(sumstats)

# ==============================================================================
# 15. PLOTS
# ==============================================================================

cat("\nGenerating figures...\n")

COL_SLACK <- "#C0392B"
COL_EXP   <- "#2980B9"
COL_BG    <- "#F8F9FA"

theme_paper <- theme_minimal(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = COL_BG,  color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#DEE2E6", linewidth = 0.35),
    legend.position  = "bottom",
    strip.text       = element_text(face = "bold", size = 9.5),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(color = "grey40", size = 9),
    axis.title       = element_text(size = 9.5)
  )

panel_labelled <- panel_est %>%
  mutate(country_label = COUNTRY_LABELS[country])

# --- Figure 1: Unemployment + HP trend (faceted) ----------------------------

p1 <- panel_labelled %>%
  filter(!is.na(unemp), !is.na(unemp_hp_trend)) %>%
  ggplot(aes(x = date)) +
  geom_rect(
    data = panel_labelled %>%
      filter(!is.na(slack_lag1), slack_lag1 == 1) %>%
      select(country_label, date),
    aes(xmin = date - 45, xmax = date + 45, ymin = -Inf, ymax = Inf),
    fill = COL_SLACK, alpha = 0.18, inherit.aes = FALSE
  ) +
  geom_line(aes(y = unemp_hp_trend, color = "HP Trend"),
            linewidth = 0.7, linetype = "dashed") +
  geom_line(aes(y = unemp, color = "Unemployment"),
            linewidth = 0.9) +
  facet_wrap(~ country_label, scales = "free_y", ncol = 3) +
  scale_color_manual(
    values = c("Unemployment" = "#2C3E50", "HP Trend" = COL_EXP),
    name   = NULL
  ) +
  labs(
    title    = "Figure 1: Unemployment Rate and HP Trend — Balkan Region",
    subtitle = "Shaded: quarters identified as 'slack' (unemp > HP trend, lagged 1Q). Source: Eurostat/World Bank",
    x = NULL, y = "Unemployment Rate (%)"
  ) + theme_paper

# --- Figure 2: Gov expenditure growth ----------------------------------------

country_colors <- setNames(
  viridis::viridis(length(USABLE_COUNTRIES), option = "D"),
  COUNTRY_LABELS[USABLE_COUNTRIES]
)

p2 <- panel_labelled %>%
  filter(!is.na(dg)) %>%
  ggplot(aes(x = date, y = dg, fill = country_label)) +
  geom_col(alpha = 0.8, width = 70) +
  geom_hline(yintercept = 0, linewidth = 0.45) +
  facet_wrap(~ country_label, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = country_colors, guide = "none") +
  labs(
    title    = "Figure 2: Real Government Expenditure Growth (QoQ, %)",
    subtitle = "Seasonally adjusted chain-linked volumes. Source: Eurostat namq_10_gdp",
    x = NULL, y = "Growth Rate (%)"
  ) + theme_paper

# --- Figure 3: Pooled panel IRF ----------------------------------------------

irf_panel_long <- irf_panel %>%
  select(horizon, starts_with("beta"), starts_with("ci")) %>%
  pivot_longer(
    cols      = -horizon,
    names_to  = c(".value", "state"),
    names_pattern = "(.+)_(slack|exp)"
  ) %>%
  mutate(
    State = ifelse(state == "slack", "Slack (Recession)", "Expansion"),
    State = factor(State, levels = c("Slack (Recession)", "Expansion"))
  )

p3 <- ggplot(irf_panel_long, aes(x = horizon, color = State, fill = State)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_ribbon(aes(ymin = ci90_lo, ymax = ci90_hi), alpha = 0.12, color = NA) +
  geom_ribbon(aes(ymin = ci68_lo, ymax = ci68_hi), alpha = 0.22, color = NA) +
  geom_line(aes(y = beta), linewidth = 1.2) +
  geom_point(aes(y = beta), size = 2) +
  scale_color_manual(values = c("Slack (Recession)" = COL_SLACK, "Expansion" = COL_EXP)) +
  scale_fill_manual( values = c("Slack (Recession)" = COL_SLACK, "Expansion" = COL_EXP)) +
  scale_x_continuous(breaks = seq(0, H, 2)) +
  labs(
    title    = "Figure 3: State-Dependent IRFs — Pooled Balkan Panel",
    subtitle = "Cumulative GDP response (% of lagged GDP) to 1pp spending shock\n90% & 68% CI. Clustered SE by country. Country FE included.",
    x = "Horizon (quarters)", y = "Cumulative Response (% of lagged GDP)",
    color = "State", fill = "State"
  ) + theme_paper

# --- Figure 4: Country-by-country IRFs ---------------------------------------

irf_ctry_long <- irf_countries %>%
  mutate(country_label = COUNTRY_LABELS[group]) %>%
  select(country_label, horizon,
         beta_slack, ci90_lo_slack, ci90_hi_slack,
         beta_exp,   ci90_lo_exp,   ci90_hi_exp) %>%
  pivot_longer(
    cols      = -c(country_label, horizon),
    names_to  = c(".value", "state"),
    names_pattern = "(.+)_(slack|exp)"
  ) %>%
  mutate(
    State = ifelse(state == "slack", "Slack", "Expansion"),
    State = factor(State, levels = c("Slack", "Expansion"))
  )

p4 <- ggplot(irf_ctry_long, aes(x = horizon, color = State, fill = State)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_ribbon(aes(ymin = ci90_lo, ymax = ci90_hi), alpha = 0.13, color = NA) +
  geom_line(aes(y = beta), linewidth = 1) +
  geom_point(aes(y = beta), size = 1.5) +
  scale_color_manual(values = c("Slack" = COL_SLACK, "Expansion" = COL_EXP)) +
  scale_fill_manual( values = c("Slack" = COL_SLACK, "Expansion" = COL_EXP)) +
  scale_x_continuous(breaks = seq(0, H, 4)) +
  facet_wrap(~ country_label, ncol = 3, scales = "free_y") +
  labs(
    title    = "Figure 4: State-Dependent IRFs by Country",
    subtitle = "Cumulative GDP response (% of lagged GDP) to 1pp spending shock. 90% CI (Newey-West HAC)",
    x = "Horizon (quarters)", y = "Cumulative Response (% of lagged GDP)",
    color = "State", fill = "State"
  ) + theme_paper

# --- Figure 5: Multiplier bar chart ------------------------------------------

mult_long <- multipliers %>%
  filter(group %in% c("Panel", USABLE_COUNTRIES)) %>%
  pivot_longer(
    cols      = c(mult_slack, mult_exp),
    names_to  = "state",
    values_to = "multiplier"
  ) %>%
  mutate(
    state         = ifelse(state == "mult_slack", "Slack", "Expansion"),
    state         = factor(state, levels = c("Slack", "Expansion")),
    horizon_label = paste0(horizon, "Q"),
    country_label = ifelse(group == "Panel", "Panel",
                           COUNTRY_LABELS[group]),
    country_label = factor(country_label,
                           levels = c("Panel", COUNTRY_LABELS[USABLE_COUNTRIES]))
  )

p5 <- ggplot(mult_long, aes(x = country_label, y = multiplier, fill = state)) +
  geom_col(position = position_dodge(0.75), width = 0.65) +
  geom_hline(yintercept = 0,  linewidth = 0.5, color = "black") +
  geom_hline(yintercept = 1,  linewidth = 0.7, color = "darkgreen", linetype = "dotted") +
  facet_wrap(~ horizon_label) +
  scale_fill_manual(values = c("Slack" = COL_SLACK, "Expansion" = COL_EXP)) +
  labs(
    title    = "Figure 5: Cumulative Government Spending Multipliers",
    subtitle = "Green dotted line = multiplier of 1.0 (spending pays for itself)\nRamey-Zubairy (2018) IV scaling",
    x = NULL, y = "Cumulative Multiplier", fill = "State"
  ) +
  theme_paper +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8))

# --- Figure 6: State-share timeline ------------------------------------------

state_share <- panel_est %>%
  filter(!is.na(slack_lag1)) %>%
  mutate(year = year(date)) %>%
  group_by(country, year) %>%
  summarise(pct_slack = mean(slack_lag1) * 100, .groups = "drop") %>%
  mutate(country_label = COUNTRY_LABELS[country])

p6 <- ggplot(state_share, aes(x = year, y = pct_slack, fill = country_label)) +
  geom_col(position = "dodge", alpha = 0.82, width = 0.7) +
  scale_fill_manual(values = country_colors, name = NULL) +
  labs(
    title    = "Figure 6: Share of Quarters in 'Slack' State by Year",
    subtitle = "% of quarters per year where lagged unemployment > HP trend",
    x = NULL, y = "% Quarters in Slack"
  ) + theme_paper +
  theme(legend.key.size = unit(0.4, "cm"))

# --- Figure 7: Multiplier heatmap across countries ---------------------------

heat_data <- multipliers %>%
  filter(group != "Panel") %>%
  mutate(country_label = COUNTRY_LABELS[group]) %>%
  filter(!is.na(country_label)) %>%
  select(country_label, horizon, mult_slack, mult_exp, diff) %>%
  pivot_longer(cols = c(mult_slack, mult_exp, diff),
               names_to = "series", values_to = "value") %>%
  mutate(
    series = case_match(series,
                        "mult_slack" ~ "Slack Multiplier",
                        "mult_exp"   ~ "Expansion Multiplier",
                        "diff"       ~ "Difference (Slack - Exp)"
    )
  )

p7 <- ggplot(
  heat_data %>% filter(series %in% c("Slack Multiplier", "Expansion Multiplier")),
  aes(x = factor(horizon), y = country_label, fill = value)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 2.6, color = "white",
            fontface = "bold") +
  scale_fill_gradient2(
    low  = "#B71C1C", mid = "grey90", high = "#1A237E",
    midpoint = 0, name = "Multiplier"
  ) +
  facet_wrap(~ series) +
  labs(
    title    = "Figure 7: Multiplier Heatmap — State and Country",
    subtitle = "Cumulative GDP response per 1pp spending shock at selected horizons",
    x = "Horizon (quarters)", y = NULL
  ) + theme_paper

# ==============================================================================
# 16. SAVE OUTPUTS
# ==============================================================================

out_dir <- file.path("outputs", format(Sys.Date(), "%Y%m%d"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("\nSaving outputs to '%s/'\n", out_dir))

fig_specs <- list(
  list(plot = p1, file = "fig1_unemployment_slack.pdf",  w = 11, h = 10),
  list(plot = p2, file = "fig2_gov_expenditure.pdf",     w = 11, h = 10),
  list(plot = p3, file = "fig3_irf_panel.pdf",           w =  9, h =  6),
  list(plot = p4, file = "fig4_irf_by_country.pdf",      w = 12, h = 12),
  list(plot = p5, file = "fig5_multipliers.pdf",         w = 10, h =  6),
  list(plot = p6, file = "fig6_state_share.pdf",         w = 11, h =  5),
  list(plot = p7, file = "fig7_multiplier_heatmap.pdf",  w = 11, h =  5)
)

invisible(lapply(fig_specs, function(s) {
  ggsave(file.path(out_dir, s$file), s$plot,
         width = s$w, height = s$h, dpi = 150)
  cat(sprintf("  Saved: %s\n", s$file))
}))

# CSV results
write.csv(irf_all,    file.path(out_dir, "irf_results.csv"),         row.names = FALSE)
write.csv(multipliers,file.path(out_dir, "multipliers.csv"),         row.names = FALSE)
write.csv(wald_table, file.path(out_dir, "hypothesis_tests.csv"),    row.names = FALSE)
write.csv(sumstats,   file.path(out_dir, "summary_statistics.csv"),  row.names = FALSE)
write.csv(coverage,   file.path(out_dir, "data_coverage.csv"),       row.names = FALSE)
write.csv(adf_table,  file.path(out_dir, "adf_unit_root_tests.csv"), row.names = FALSE)
write.csv(mult_rob,   file.path(out_dir, "robustness_multipliers.csv"), row.names = FALSE)
write.csv(fs_diag,    file.path(out_dir, "first_stage_fstats.csv"),  row.names = FALSE)

cat("\n=== Done. All results saved. ===\n")

# Session info for reproducibility
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "session_info.txt"))