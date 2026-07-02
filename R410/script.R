################################################################################
# STATE-DEPENDENT GOVERNMENT SPENDING MULTIPLIERS
# A Regional Study: Western Balkans / South-Eastern Europe
#   + G5 Advanced Economy Comparison Group
#   (USA, Canada, United Kingdom, Germany, France)
#
# Countries — Balkans:
#   Albania (AL), Bosnia & Herzegovina (BA), Bulgaria (BG),
#   Croatia (HR), Greece (EL), Kosovo (XK), Montenegro (ME),
#   North Macedonia (MK), Romania (RO), Serbia (RS), Slovenia (SI)
#
# Countries — G5 Comparison:
#   United States (US), Canada (CA), United Kingdom (GB),
#   Germany (DE), France (FR)
#
# Methodology:
#   Owyang, Ramey & Zubairy (2013) — State-Dependent Local Projections
#   Ramey & Zubairy (2018)         — IV Cumulative Multiplier Scaling
#   Blanchard & Perotti (2002)     — Spending shock identification
#
# Data Sources:
#   Balkan countries  → Eurostat API  ({eurostat} R package)
#   G5 countries      → OECD.Stat API ({OECD} R package)
#     GDP & Gov. Expenditure : OECD QNA (Quarterly National Accounts)
#     Unemployment           : OECD MEI (Main Economic Indicators)
#
# Sample restriction for comparability:
#   All groups are estimated over the SAME sample window that is available
#   for the Balkan panel — nominally 1998 Q1 to latest available quarter.
#   Countries with fewer than MIN_OBS usable quarters are dropped automatically.
#
# NOTE ON ORZ (2013) REPLICATION:
#   The original paper used US data back to 1890 (war-based instruments).
#   Here we use the Blanchard-Perotti (2002) identification (no war instruments)
#   and restrict ALL countries — including the G5 — to the post-1998 window
#   so that every multiplier estimate is directly comparable across groups.
#
# Author : [Your name]
# Date   : Sys.Date()
################################################################################

set.seed(20240601)

# ==============================================================================
# 00. WORKING DIRECTORY
# ==============================================================================

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  cat("Working directory:", getwd(), "\n\n")
} else {
  cat("Not in RStudio. Set working directory manually if needed.\n\n")
}

# ==============================================================================
# 01. PACKAGES
# ==============================================================================

required_pkgs <- c(
  # Data retrieval
  "eurostat", "rdbnomics", 
  # Wrangling
  "dplyr", "tidyr", "lubridate", "zoo", "purrr",
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

# ── Balkan countries (Eurostat) ──────────────────────────────────────────────

BALKAN_COUNTRIES <- c(
  "AL",   # Albania
  "BA",   # Bosnia & Herzegovina
  "BG",   # Bulgaria
  "EL",   # Greece
  "HR",   # Croatia
  "ME",   # Montenegro
  "MK",   # North Macedonia
  "RO",   # Romania
  "RS",   # Serbia
  "SI",   # Slovenia
  "XK"    # Kosovo
)

BALKAN_LABELS <- c(
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

# ── G5 countries (OECD) ──────────────────────────────────────────────────────
# OECD location codes used in the QNA and MEI datasets

G5_COUNTRIES_OECD <- c(
  "USA",  # United States
  "CAN",  # Canada
  "GBR",  # United Kingdom
  "DEU",  # Germany
  "FRA"   # France
)

# Harmonised 2-character codes used internally (match ISO-2 where possible)
G5_COUNTRIES <- c("US", "CA", "GB", "DE", "FR")

G5_LABELS <- c(
  US = "United States",
  CA = "Canada",
  GB = "United Kingdom",
  DE = "Germany",
  FR = "France"
)

# Lookup: internal 2-char code → OECD 3-char code
G5_OECD_MAP <- setNames(G5_COUNTRIES_OECD, G5_COUNTRIES)

# ── Combined lookup (used in plotting) ───────────────────────────────────────

ALL_LABELS <- c(BALKAN_LABELS, G5_LABELS)

# ── Global estimation parameters ─────────────────────────────────────────────

H            <- 12    # Maximum LP horizon (quarters)
HP_FREQ      <- 1600  # HP filter smoothing (quarterly)
CI90_Z       <- 1.645
CI68_Z       <- 1.000
NW_PREWHITE  <- FALSE
MIN_OBS      <- 24    # Minimum usable quarters per country

# Sample window: restrict G5 to the same era as Balkans for comparability.
# Data back to Q1 1995 is requested; the minimum-obs filter does the trimming.
SAMPLE_START <- as.Date("1995-01-01")   # earliest date to keep
SAMPLE_END   <- Sys.Date()              # up to latest available

# ==============================================================================
# 03A. EUROSTAT DOWNLOAD — BALKAN COUNTRIES
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

cat("=== [A] EUROSTAT: Balkan GDP & Government Expenditure ===\n")
gdp_raw_bk <- safe_eurostat(
  id      = "namq_10_gdp",
  filters = list(
    geo     = BALKAN_COUNTRIES,
    na_item = c("B1GQ", "P3_S13"),
    unit    = "CLV10_MNAC",
    s_adj   = "SCA"
  )
)

cat("Downloading Balkan unemployment (Eurostat lfsq_urgan)...\n")
unemp_raw_bk <- safe_eurostat(
  id      = "lfsq_urgan",
  filters = list(
    geo     = BALKAN_COUNTRIES,
    sex     = "T",
    unit    = "PC",
    age     = "Y15-74",
    citizen = "TOTAL"
  )
)

# ==============================================================================
# 03B. DBNOMICS DOWNLOAD — G5 COUNTRIES (Replacing OECD package)
# ==============================================================================

safe_rdb <- function(dataset, dimensions, ...) {
  tryCatch({
    rdbnomics::rdb(provider_code = "OECD", dataset_code = dataset, dimensions = dimensions)
  }, error = function(e) {
    warning(sprintf("DBnomics download failed for '%s': %s", dataset, e$message))
    NULL
  })
}

cat("\n=== [B] DBnomics: G5 Quarterly National Accounts ===\n")

# QNA filter string: LOCATION.SUBJECT.MEASURE.FREQUENCY
qna_filter <- list(
  LOCATION  = G5_COUNTRIES_OECD,
  SUBJECT   = c("B1_GE", "P3S13"),
  MEASURE   = c("LNBQRSA"),
  FREQUENCY = c("Q")
)

qna_raw <- safe_rdb("QNA", qna_filter)

cat("Downloading G5 unemployment (DBnomics / OECD MEI)...\n")
mei_filter <- list(
  LOCATION  = G5_COUNTRIES_OECD,
  SUBJECT   = c("LRUNTTTT"), 
  MEASURE   = c("ST"),       
  FREQUENCY = c("M")
)

unemp_raw_g5 <- safe_rdb("MEI", mei_filter)

# ==============================================================================
# 04. DATA CLEANING — BALKAN (Eurostat)
# ==============================================================================

cat("\n=== Cleaning Balkan data ===\n")

clean_eurostat_panel <- function(gdp_raw, unemp_raw) {
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
    filter(!is.na(gdp), !is.na(gov_exp)) %>%
    filter(date >= SAMPLE_START, date <= SAMPLE_END) %>%
    mutate(group_id = "Balkan")
}

if (is.null(gdp_raw_bk) || is.null(unemp_raw_bk)) {
  stop("Eurostat returned NULL for Balkan data. Check network connection.")
}
panel_bk <- clean_eurostat_panel(gdp_raw_bk, unemp_raw_bk)

cat(sprintf("  Balkan panel: %d obs, %d countries\n",
            nrow(panel_bk), length(unique(panel_bk$country))))

# ==============================================================================
# 05. DATA CLEANING — G5 (DBnomics)
# ==============================================================================

cat("\n=== Cleaning G5 data ===\n")

# ── 5a. DBnomics QNA → GDP and Gov Expenditure ───────────────────────────────

clean_oecd_qna <- function(qna_raw) {
  if (is.null(qna_raw) || nrow(qna_raw) == 0) return(NULL)
  
  qna_raw %>%
    filter(!is.na(value)) %>%
    select(LOCATION, SUBJECT, period, value) %>%
    rename(
      oecd_code = LOCATION,
      date      = period
    ) %>%
    pivot_wider(names_from = SUBJECT, values_from = value) %>%
    rename(
      gdp       = B1_GE,
      gov_exp   = P3S13
    ) %>%
    filter(!is.na(gdp), !is.na(gov_exp))
}

# ── 5b. DBnomics MEI → Unemployment (monthly → quarterly average) ─────────────

clean_oecd_mei <- function(unemp_raw_g5) {
  if (is.null(unemp_raw_g5) || nrow(unemp_raw_g5) == 0) return(NULL)
  
  unemp_raw_g5 %>%
    filter(!is.na(value)) %>%
    mutate(
      # The 'period' column is already a Date object (first of the month)
      date_q = as.Date(zoo::as.yearqtr(period))
    ) %>%
    group_by(LOCATION, date_q) %>%
    summarise(unemp = mean(value, na.rm = TRUE), .groups = "drop") %>%
    rename(oecd_code = LOCATION, date = date_q)
}

# ── 5c. Merge QNA + MEI and harmonise codes ──────────────────────────────────

build_g5_panel <- function(qna_raw, unemp_raw_g5) {
  qna_clean   <- clean_oecd_qna(qna_raw)
  unemp_clean <- clean_oecd_mei(unemp_raw_g5)
  
  if (is.null(qna_clean)) stop("QNA cleaning returned NULL.")
  if (is.null(unemp_clean)) {
    warning("MEI cleaning returned NULL — unemployment missing for G5.")
    unemp_clean <- tibble(oecd_code = character(), date = as.Date(character()),
                          unemp = numeric())
  }
  
  panel <- left_join(qna_clean, unemp_clean, by = c("oecd_code", "date")) %>%
    filter(!is.na(gdp), !is.na(gov_exp)) %>%
    # Map OECD 3-char codes back to internal 2-char codes
    mutate(
      country = names(G5_OECD_MAP)[match(oecd_code, G5_OECD_MAP)],
      group_id = "G5"
    ) %>%
    filter(!is.na(country)) %>%
    filter(date >= SAMPLE_START, date <= SAMPLE_END) %>%
    select(country, date, gdp, gov_exp, unemp, group_id) %>%
    arrange(country, date)
  
  panel
}

panel_g5 <- build_g5_panel(qna_raw, unemp_raw_g5)

cat(sprintf("  G5 panel: %d obs, %d countries\n",
            nrow(panel_g5), length(unique(panel_g5$country))))

# ==============================================================================
# 06. COMBINE PANELS
# ==============================================================================

cat("\n=== Combining panels ===\n")

panel_raw <- bind_rows(panel_bk, panel_g5) %>%
  arrange(group_id, country, date)

cat(sprintf("  Combined panel: %d obs, %d countries\n",
            nrow(panel_raw), length(unique(panel_raw$country))))

# ==============================================================================
# 07. VARIABLE CONSTRUCTION
# ==============================================================================

cat("\nConstructing variables...\n")

safe_hp <- function(x, freq = HP_FREQ, min_obs = 8) {
  idx    <- which(!is.na(x))
  result <- rep(NA_real_, length(x))
  if (length(idx) >= min_obs) {
    hp          <- hpfilter(x[idx], freq = freq)
    result[idx] <- hp$trend
  }
  result
}

ldiff <- function(x) c(NA, diff(log(x))) * 100

panel <- panel_raw %>%
  group_by(country) %>%
  arrange(date) %>%
  mutate(
    ln_gdp = log(gdp),
    ln_g   = log(gov_exp),
    
    dy = ldiff(gdp),
    dg = ldiff(gov_exp),
    
    gdp_lag1 = lag(gdp, 1),
    
    # Ramey-Zubairy (2018) scaling: change relative to lagged GDP
    g_scaled = (gov_exp - lag(gov_exp)) / gdp_lag1 * 100,
    y_scaled = (gdp    - lag(gdp))     / gdp_lag1 * 100,
    
    # Unemployment gap (HP-based)
    unemp_hp_trend = safe_hp(unemp),
    unemp_gap      = unemp - unemp_hp_trend,
    
    # State indicator: 1 = slack, 0 = expansion
    slack      = as.integer(unemp_gap > 0),
    slack_lag1 = lag(slack, 1),
    
    # Lags for controls
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
# 08. DATA QUALITY CHECK & COUNTRY SELECTION
# ==============================================================================

cat("\n=== Data Coverage by Country (both groups) ===\n")

coverage <- panel %>%
  filter(!is.na(g_scaled), !is.na(y_scaled), !is.na(slack_lag1)) %>%
  group_by(country, group_id) %>%
  summarise(
    start_q   = min(date),
    end_q     = max(date),
    n_obs     = n(),
    n_slack   = sum(slack_lag1, na.rm = TRUE),
    pct_slack = round(mean(slack_lag1, na.rm = TRUE) * 100, 1),
    .groups   = "drop"
  ) %>%
  mutate(country_name = ALL_LABELS[country]) %>%
  arrange(group_id, country)

print(coverage)
write.csv(coverage, "data_coverage_all.csv", row.names = FALSE)

# Countries passing the minimum-observations filter
USABLE_COUNTRIES <- coverage %>%
  filter(n_obs >= MIN_OBS, n_slack >= 6, n_obs - n_slack >= 6) %>%
  pull(country)

# BUG FIX #2: use intersect() with names(BALKAN_LABELS) directly —
# the original used a circular self-filter that could silently drop countries.
USABLE_BALKAN <- intersect(USABLE_COUNTRIES, names(BALKAN_LABELS))
USABLE_G5     <- intersect(USABLE_COUNTRIES, G5_COUNTRIES)

cat(sprintf("\nUsable Balkan countries : %s\n", paste(USABLE_BALKAN, collapse = ", ")))
cat(sprintf("Usable G5 countries     : %s\n", paste(USABLE_G5,     collapse = ", ")))

panel_est <- panel %>%
  filter(
    country %in% USABLE_COUNTRIES,
    !is.na(g_scaled), !is.na(y_scaled),
    !is.na(slack_lag1),
    !is.na(g_scaled_lag1), !is.na(g_scaled_lag2),
    !is.na(y_scaled_lag1), !is.na(y_scaled_lag2)
  )

# ==============================================================================
# 09. UNIT ROOT TESTS (ADF)
# ==============================================================================

cat("\n--- Augmented Dickey-Fuller tests (by country) ---\n")

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
# 10. LOCAL PROJECTION IV ESTIMATOR
# ==============================================================================
#
# Identification : Blanchard-Perotti (2002)
# Estimator      : 2SLS via ivreg()
# SE             : Newey-West HAC (bandwidth = h+1) or clustered by country
# State          : slack_lag1 based on HP-filtered unemployment gap

run_lp_iv <- function(data, h, cluster_se = FALSE) {
  
  CTRL <- c("y_scaled_lag1", "y_scaled_lag2", "g_scaled_lag1", "g_scaled_lag2")
  
  df <- data %>%
    group_by(country) %>%
    arrange(date) %>%
    mutate(
      # BUG FIX #3: pass na.rm = TRUE to the rollapply FUN so partial windows
      # with leading/trailing NAs do not silently propagate NA into cum_y_h /
      # cum_g_h for otherwise-valid observations.
      cum_y_h = rollapply(y_scaled, width = h + 1,
                          FUN = function(x) sum(x, na.rm = TRUE),
                          align = "left", fill = NA, partial = FALSE),
      cum_g_h = rollapply(g_scaled, width = h + 1,
                          FUN = function(x) sum(x, na.rm = TRUE),
                          align = "left", fill = NA, partial = FALSE)
    ) %>%
    ungroup() %>%
    filter(!is.na(cum_y_h), !is.na(cum_g_h), !is.na(slack_lag1)) %>%
    mutate(
      state_slack = slack_lag1,
      state_exp   = 1L - slack_lag1,
      cum_g_slack = state_slack * cum_g_h,
      cum_g_exp   = state_exp   * cum_g_h,
      shock_slack = state_slack * g_scaled,
      shock_exp   = state_exp   * g_scaled,
      country     = factor(country)
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
    fit      <- ivreg(fml, data = df)
    vcov_use <- if (cluster_se && n_countries > 1) {
      vcovCL(fit, cluster = ~country)
    } else {
      NeweyWest(fit, lag = h + 1, prewhite = NW_PREWHITE)
    }
    
    ct   <- coeftest(fit, vcov = vcov_use)
    wald <- car::linearHypothesis(fit, "cum_g_slack = cum_g_exp",
                                  vcov. = vcov_use, test = "F")
    
    fs_slack <- tryCatch({
      fs_fit <- lm(as.formula(paste0(
        "cum_g_slack ~ shock_slack + shock_exp + ", ctrl_str, fe_str
      )), data = df)
      summary(fs_fit)$fstatistic[1]
    }, error = function(e) NA_real_)
    
    pval_col <- grep("Pr\\(", colnames(ct), value = TRUE)[1]
    
    list(
      h           = h,
      beta_slack  = ct["cum_g_slack", "Estimate"],
      beta_exp    = ct["cum_g_exp",   "Estimate"],
      se_slack    = ct["cum_g_slack", "Std. Error"],
      se_exp      = ct["cum_g_exp",   "Std. Error"],
      pval_slack  = ct["cum_g_slack", pval_col],
      pval_exp    = ct["cum_g_exp",   pval_col],
      pval_diff   = wald$`Pr(>F)`[2],
      fs_f_stat   = fs_slack,
      n           = nrow(df),
      n_countries = n_countries,
      df_resid    = fit$df.residual
    )
  }, error = function(e) {
    warning(sprintf("LP failed at h=%d: %s", h, conditionMessage(e)))
    NULL
  })
}

# ==============================================================================
# 11. ESTIMATION — FOUR POOLED PANELS
#
#   (1) Balkan panel          — country FE, clustered SE
#   (2) G5 panel              — country FE, clustered SE
#   (3) Full panel (All)      — country FE, clustered SE
#   (4) Country-by-country    — Newey-West SE
# ==============================================================================

n_cores <- max(1L, detectCores() - 1L)
cat(sprintf("\nEstimating LP-IV across horizons 0:%d ...\n", H))

run_panel_lp <- function(data, label, cluster = TRUE) {
  cat(sprintf("  [%s] ...\n", label))
  res <- pblapply(0:H, function(h) run_lp_iv(data, h, cluster_se = cluster))
  res <- Filter(Negate(is.null), res)
  if (length(res) == 0) {
    warning(sprintf("No LP results for panel '%s'.", label))
    return(NULL)
  }
  res
}

# ─── (1) Balkan pooled ───────────────────────────────────────────────────────
lp_balkan <- run_panel_lp(
  panel_est %>% filter(country %in% USABLE_BALKAN), "Balkan panel"
)

# ─── (2) G5 pooled ───────────────────────────────────────────────────────────
lp_g5 <- run_panel_lp(
  panel_est %>% filter(country %in% USABLE_G5), "G5 panel"
)

# ─── (3) Full pooled (all countries) ─────────────────────────────────────────
lp_all <- run_panel_lp(panel_est, "Full panel (All countries)")

# ─── (4) Country-by-country ──────────────────────────────────────────────────
cat("  [Country-by-country] ...\n")
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
# 12. EXTRACT IRF DATA FRAMES
# ==============================================================================

extract_irf <- function(lp_list, label) {
  if (is.null(lp_list) || length(lp_list) == 0) {
    warning(sprintf("extract_irf: no results for '%s'.", label))
    return(NULL)
  }
  do.call(rbind, lapply(lp_list, function(r) {
    data.frame(
      group         = label,
      horizon       = r$h,
      beta_slack    = r$beta_slack,
      se_slack      = r$se_slack,
      beta_exp      = r$beta_exp,
      se_exp        = r$se_exp,
      pval_diff     = r$pval_diff,
      fs_f_stat     = r$fs_f_stat,
      n             = r$n,
      ci90_lo_slack = r$beta_slack - CI90_Z * r$se_slack,
      ci90_hi_slack = r$beta_slack + CI90_Z * r$se_slack,
      ci90_lo_exp   = r$beta_exp   - CI90_Z * r$se_exp,
      ci90_hi_exp   = r$beta_exp   + CI90_Z * r$se_exp,
      ci68_lo_slack = r$beta_slack - CI68_Z * r$se_slack,
      ci68_hi_slack = r$beta_slack + CI68_Z * r$se_slack,
      ci68_lo_exp   = r$beta_exp   - CI68_Z * r$se_exp,
      ci68_hi_exp   = r$beta_exp   + CI68_Z * r$se_exp
    )
  }))
}

irf_balkan  <- extract_irf(lp_balkan, "Balkan")
irf_g5      <- extract_irf(lp_g5,    "G5")
irf_all_pan <- extract_irf(lp_all,   "All Countries")

irf_countries <- do.call(rbind, Filter(Negate(is.null), lapply(USABLE_COUNTRIES, function(cc) {
  extract_irf(lp_by_country[[cc]], cc)
})))

# Master IRF table
irf_all <- bind_rows(irf_balkan, irf_g5, irf_all_pan, irf_countries)

# ==============================================================================
# 13. MULTIPLIER TABLES
# ==============================================================================

make_mult_table <- function(irf_df, horizons = c(4, 8, 12)) {
  if (is.null(irf_df) || nrow(irf_df) == 0) return(NULL)
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

cat("\n=== CUMULATIVE MULTIPLIERS — ALL GROUPS ===\n")
print(multipliers)

# Separate tables for readability
mult_panels <- multipliers %>%
  filter(group %in% c("Balkan", "G5", "All Countries"))

# BUG FIX #6: derive economy_type correctly from the group (2-char code) column
# so that country-level rows get "Balkan" or "G5" rather than NA.
mult_countries <- multipliers %>%
  filter(!group %in% c("Balkan", "G5", "All Countries")) %>%
  mutate(
    country_name = ALL_LABELS[group],
    economy_type = ifelse(group %in% G5_COUNTRIES, "G5", "Balkan")
  )

cat("\n=== PANEL MULTIPLIERS (Balkan vs G5 vs All) ===\n")
print(mult_panels)

cat("\n=== COUNTRY-LEVEL MULTIPLIERS ===\n")
print(mult_countries)

# ==============================================================================
# 14. FIRST-STAGE DIAGNOSTICS
# ==============================================================================

cat("\n--- First-stage F-statistics (rule of thumb: F > 10) ---\n")

fs_diag <- bind_rows(
  irf_balkan  %>% mutate(panel = "Balkan"),
  irf_g5      %>% mutate(panel = "G5"),
  irf_all_pan %>% mutate(panel = "All")
) %>%
  select(panel, horizon, fs_f_stat) %>%
  mutate(fs_f_stat = round(fs_f_stat, 2), weak_iv = fs_f_stat < 10)

print(fs_diag)

# ==============================================================================
# 15. HYPOTHESIS TESTS — WALD TESTS ON PANEL ESTIMATES
# ==============================================================================

cat("\n--- Wald tests: H0: β_slack = β_exp ---\n")

format_wald <- function(irf_df, label) {
  irf_df %>%
    transmute(
      panel      = label,
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
}

wald_table <- bind_rows(
  if (!is.null(irf_balkan))  format_wald(irf_balkan,  "Balkan"),
  if (!is.null(irf_g5))      format_wald(irf_g5,      "G5"),
  if (!is.null(irf_all_pan)) format_wald(irf_all_pan, "All")
)

print(wald_table)

# ==============================================================================
# 16. ROBUSTNESS — PERCENTILE-BASED SLACK THRESHOLDS
# ==============================================================================

cat("\n--- Robustness: percentile-based slack thresholds ---\n")

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

# Run robustness for each group separately at P25, P50, P75
run_rob_group <- function(countries_vec, prefix) {
  sub <- panel_rob %>% filter(country %in% countries_vec)
  bind_rows(
    run_rob_spec(sub, "slack_p25", paste0(prefix, ": u > P25")),
    run_rob_spec(sub, "slack_p50", paste0(prefix, ": u > P50")),
    run_rob_spec(sub, "slack_p75", paste0(prefix, ": u > P75"))
  )
}

irf_rob <- bind_rows(
  run_rob_group(USABLE_BALKAN, "Balkan"),
  run_rob_group(USABLE_G5,     "G5")
)

mult_rob <- make_mult_table(irf_rob)
cat("\n=== ROBUSTNESS MULTIPLIERS ===\n")
print(mult_rob)

# ==============================================================================
# 17. SUMMARY STATISTICS
# ==============================================================================

sumstats <- panel_est %>%
  filter(!is.na(slack_lag1)) %>%
  mutate(
    state        = ifelse(slack_lag1 == 1, "Slack", "Expansion"),
    country_name = ALL_LABELS[country],
    economy_type = ifelse(country %in% G5_COUNTRIES, "G5", "Balkan")
  ) %>%
  group_by(economy_type, country_name, state) %>%
  summarise(
    quarters    = n(),
    mean_gdp_gr = round(mean(dy,    na.rm = TRUE), 2),
    mean_g_gr   = round(mean(dg,    na.rm = TRUE), 2),
    mean_unemp  = round(mean(unemp, na.rm = TRUE), 2),
    sd_gdp_gr   = round(sd(dy,      na.rm = TRUE), 2),
    sd_g_gr     = round(sd(dg,      na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  arrange(economy_type, country_name, state)

cat("\n=== SUMMARY STATISTICS ===\n")
print(sumstats)

# ==============================================================================
# 18. PLOTS
# ==============================================================================

cat("\nGenerating figures...\n")

COL_SLACK  <- "#C0392B"
COL_EXP    <- "#2980B9"
COL_BG     <- "#F8F9FA"
COL_BALKAN <- "#7D3C98"   # Purple for Balkan panel
COL_G5     <- "#1E8449"   # Green for G5 panel
COL_ALL    <- "#CA6F1E"   # Orange for full panel

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
  mutate(
    country_label = ALL_LABELS[country],
    economy_type  = ifelse(country %in% G5_COUNTRIES, "G5", "Balkan")
  )

all_countries_ordered <- c(
  # BUG FIX #1: was COUNTRY_LABELS[USABLE_BALKAN] — undefined object.
  # Corrected to BALKAN_LABELS[USABLE_BALKAN].
  BALKAN_LABELS[USABLE_BALKAN],
  G5_LABELS[USABLE_G5]
)

country_colors <- setNames(
  c(viridis::viridis(length(USABLE_BALKAN), option = "D"),
    viridis::viridis(length(USABLE_G5),     option = "C")),
  all_countries_ordered
)

# ── Figure 1: Unemployment + HP trend (Balkan, faceted) ──────────────────────

bk_labelled <- panel_labelled %>% filter(country %in% USABLE_BALKAN)

p1 <- bk_labelled %>%
  filter(!is.na(unemp), !is.na(unemp_hp_trend)) %>%
  ggplot(aes(x = date)) +
  geom_rect(
    data = bk_labelled %>%
      filter(!is.na(slack_lag1), slack_lag1 == 1) %>%
      select(country_label, date),
    aes(xmin = date - 45, xmax = date + 45, ymin = -Inf, ymax = Inf),
    fill = COL_SLACK, alpha = 0.18, inherit.aes = FALSE
  ) +
  geom_line(aes(y = unemp_hp_trend, color = "HP Trend"),
            linewidth = 0.7, linetype = "dashed") +
  geom_line(aes(y = unemp, color = "Unemployment"), linewidth = 0.9) +
  facet_wrap(~ country_label, scales = "free_y", ncol = 3) +
  scale_color_manual(
    values = c("Unemployment" = "#2C3E50", "HP Trend" = COL_EXP), name = NULL
  ) +
  labs(
    title    = "Figure 1a: Unemployment Rate & HP Trend — Balkan Countries",
    subtitle = "Shaded: slack quarters (unemp > HP trend, lagged 1Q). Source: Eurostat",
    x = NULL, y = "Unemployment Rate (%)"
  ) + theme_paper

# ── Figure 1b: Unemployment + HP trend (G5, faceted) ─────────────────────────

g5_labelled <- panel_labelled %>% filter(country %in% USABLE_G5)

p1b <- g5_labelled %>%
  filter(!is.na(unemp), !is.na(unemp_hp_trend)) %>%
  ggplot(aes(x = date)) +
  geom_rect(
    data = g5_labelled %>%
      filter(!is.na(slack_lag1), slack_lag1 == 1) %>%
      select(country_label, date),
    aes(xmin = date - 45, xmax = date + 45, ymin = -Inf, ymax = Inf),
    fill = COL_SLACK, alpha = 0.18, inherit.aes = FALSE
  ) +
  geom_line(aes(y = unemp_hp_trend, color = "HP Trend"),
            linewidth = 0.7, linetype = "dashed") +
  geom_line(aes(y = unemp, color = "Unemployment"), linewidth = 0.9) +
  facet_wrap(~ country_label, scales = "free_y", ncol = 3) +
  scale_color_manual(
    values = c("Unemployment" = "#2C3E50", "HP Trend" = COL_EXP), name = NULL
  ) +
  labs(
    title    = "Figure 1b: Unemployment Rate & HP Trend — G5 Countries",
    subtitle = "Shaded: slack quarters (unemp > HP trend, lagged 1Q). Source: OECD MEI",
    x = NULL, y = "Unemployment Rate (%)"
  ) + theme_paper

# ── Figure 2: Gov expenditure growth (both groups, faceted) ───────────────────

p2 <- panel_labelled %>%
  filter(!is.na(dg)) %>%
  ggplot(aes(x = date, y = dg, fill = country_label)) +
  geom_col(alpha = 0.8, width = 70) +
  geom_hline(yintercept = 0, linewidth = 0.45) +
  facet_wrap(~ country_label, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = country_colors, guide = "none") +
  labs(
    title    = "Figure 2: Real Government Expenditure Growth (QoQ, %)",
    subtitle = "Balkan: Eurostat CLV10_MNAC  |  G5: OECD QNA LNBQRSA",
    x = NULL, y = "Growth Rate (%)"
  ) + theme_paper

# ── Figure 3: CORE COMPARISON — Pooled Balkan vs G5 IRFs ─────────────────────

prep_irf_long <- function(irf_df, label_col = "group") {
  irf_df %>%
    select(all_of(label_col), horizon,
           beta_slack, se_slack, ci90_lo_slack, ci90_hi_slack, ci68_lo_slack, ci68_hi_slack,
           beta_exp,   se_exp,   ci90_lo_exp,   ci90_hi_exp,   ci68_lo_exp,   ci68_hi_exp) %>%
    pivot_longer(
      cols      = -c(all_of(label_col), horizon),
      names_to  = c(".value", "state"),
      names_pattern = "(.+)_(slack|exp)"
    ) %>%
    mutate(
      State = ifelse(state == "slack", "Slack (Recession)", "Expansion"),
      State = factor(State, levels = c("Slack (Recession)", "Expansion"))
    )
}

# Combine Balkan and G5 panel IRFs
irf_compare <- bind_rows(irf_balkan, irf_g5)

irf_compare_long <- prep_irf_long(irf_compare) %>%
  mutate(panel_label = group)

p3 <- ggplot(irf_compare_long,
             aes(x = horizon, color = State, fill = State,
                 linetype = panel_label, shape = panel_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_ribbon(aes(ymin = ci90_lo, ymax = ci90_hi), alpha = 0.10, color = NA) +
  geom_ribbon(aes(ymin = ci68_lo, ymax = ci68_hi), alpha = 0.18, color = NA) +
  geom_line(aes(y = beta), linewidth = 1.2) +
  geom_point(aes(y = beta), size = 2.5) +
  scale_color_manual(values = c("Slack (Recession)" = COL_SLACK, "Expansion" = COL_EXP)) +
  scale_fill_manual( values = c("Slack (Recession)" = COL_SLACK, "Expansion" = COL_EXP)) +
  scale_linetype_manual(values = c("Balkan" = "solid", "G5" = "dashed"), name = "Panel") +
  scale_shape_manual(  values = c("Balkan" = 16,      "G5" = 17),        name = "Panel") +
  scale_x_continuous(breaks = seq(0, H, 2)) +
  facet_wrap(~ State, ncol = 2) +
  labs(
    title    = "Figure 3: State-Dependent IRFs — Balkan vs G5 (CORE COMPARISON)",
    subtitle = paste0("Cumulative GDP response (% of lagged GDP) to 1pp spending shock\n",
                      "90% & 68% CI. Clustered SE. Country FE. Same sample window."),
    x = "Horizon (quarters)", y = "Cumulative Response (% of lagged GDP)"
  ) + theme_paper

# ── Figure 4: Three-panel comparison (Balkan / G5 / All) ─────────────────────

irf_3panel <- bind_rows(
  irf_balkan  %>% mutate(group = "Balkan"),
  irf_g5      %>% mutate(group = "G5"),
  irf_all_pan %>% mutate(group = "All Countries")
) %>% mutate(group = factor(group, levels = c("Balkan", "G5", "All Countries")))

irf_3long <- irf_3panel %>%
  select(group, horizon,
         beta_slack, ci90_lo_slack, ci90_hi_slack,
         beta_exp,   ci90_lo_exp,   ci90_hi_exp) %>%
  pivot_longer(
    cols      = -c(group, horizon),
    names_to  = c(".value", "state"),
    names_pattern = "(.+)_(slack|exp)"
  ) %>%
  mutate(
    State = ifelse(state == "slack", "Slack", "Expansion"),
    State = factor(State, levels = c("Slack", "Expansion"))
  )

p4 <- ggplot(irf_3long, aes(x = horizon, color = State, fill = State)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_ribbon(aes(ymin = ci90_lo, ymax = ci90_hi), alpha = 0.12, color = NA) +
  geom_line(aes(y = beta), linewidth = 1.1) +
  geom_point(aes(y = beta), size = 1.8) +
  scale_color_manual(values = c("Slack" = COL_SLACK, "Expansion" = COL_EXP)) +
  scale_fill_manual( values = c("Slack" = COL_SLACK, "Expansion" = COL_EXP)) +
  scale_x_continuous(breaks = seq(0, H, 2)) +
  facet_wrap(~ group, ncol = 3) +
  labs(
    title    = "Figure 4: IRFs by Panel Group",
    subtitle = "Balkan | G5 | Full (all countries). 90% CI.",
    x = "Horizon (quarters)", y = "Cumulative Response (% of lagged GDP)",
    color = "State", fill = "State"
  ) + theme_paper

# ── Figure 5: Country-by-country IRFs — Balkan ───────────────────────────────

irf_ctry_bk_long <- irf_countries %>%
  filter(group %in% USABLE_BALKAN) %>%
  mutate(country_label = BALKAN_LABELS[group]) %>%
  filter(!is.na(country_label)) %>%
  select(country_label, horizon,
         beta_slack, ci90_lo_slack, ci90_hi_slack,
         beta_exp,   ci90_lo_exp,   ci90_hi_exp) %>%
  pivot_longer(
    cols = -c(country_label, horizon),
    names_to  = c(".value", "state"),
    names_pattern = "(.+)_(slack|exp)"
  ) %>%
  mutate(
    State = ifelse(state == "slack", "Slack", "Expansion"),
    State = factor(State, levels = c("Slack", "Expansion"))
  )

p5 <- ggplot(irf_ctry_bk_long, aes(x = horizon, color = State, fill = State)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_ribbon(aes(ymin = ci90_lo, ymax = ci90_hi), alpha = 0.13, color = NA) +
  geom_line(aes(y = beta), linewidth = 1) +
  geom_point(aes(y = beta), size = 1.5) +
  scale_color_manual(values = c("Slack" = COL_SLACK, "Expansion" = COL_EXP)) +
  scale_fill_manual( values = c("Slack" = COL_SLACK, "Expansion" = COL_EXP)) +
  scale_x_continuous(breaks = seq(0, H, 4)) +
  facet_wrap(~ country_label, ncol = 3, scales = "free_y") +
  labs(
    title    = "Figure 5: Country IRFs — Balkan Region",
    subtitle = "Cumulative GDP response (% of lagged GDP) to 1pp spending shock. 90% CI (NW HAC)",
    x = "Horizon (quarters)", y = "Cumulative Response",
    color = "State", fill = "State"
  ) + theme_paper

# ── Figure 6: Country-by-country IRFs — G5 ───────────────────────────────────

irf_ctry_g5_long <- irf_countries %>%
  filter(group %in% USABLE_G5) %>%
  mutate(country_label = G5_LABELS[group]) %>%
  filter(!is.na(country_label)) %>%
  select(country_label, horizon,
         beta_slack, ci90_lo_slack, ci90_hi_slack,
         beta_exp,   ci90_lo_exp,   ci90_hi_exp) %>%
  pivot_longer(
    cols = -c(country_label, horizon),
    names_to  = c(".value", "state"),
    names_pattern = "(.+)_(slack|exp)"
  ) %>%
  mutate(
    State = ifelse(state == "slack", "Slack", "Expansion"),
    State = factor(State, levels = c("Slack", "Expansion"))
  )

p6 <- ggplot(irf_ctry_g5_long, aes(x = horizon, color = State, fill = State)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_ribbon(aes(ymin = ci90_lo, ymax = ci90_hi), alpha = 0.13, color = NA) +
  geom_line(aes(y = beta), linewidth = 1) +
  geom_point(aes(y = beta), size = 1.8) +
  scale_color_manual(values = c("Slack" = COL_SLACK, "Expansion" = COL_EXP)) +
  scale_fill_manual( values = c("Slack" = COL_SLACK, "Expansion" = COL_EXP)) +
  scale_x_continuous(breaks = seq(0, H, 4)) +
  facet_wrap(~ country_label, ncol = 3, scales = "free_y") +
  labs(
    title    = "Figure 6: Country IRFs — G5 Economies",
    subtitle = "Cumulative GDP response (% of lagged GDP) to 1pp spending shock. 90% CI (NW HAC)",
    x = "Horizon (quarters)", y = "Cumulative Response",
    color = "State", fill = "State"
  ) + theme_paper

# ── Figure 7: Multiplier bar chart — Balkan vs G5 (panel + country) ──────────

mult_bar <- mult_countries %>%
  bind_rows(
    mult_panels %>%
      filter(group %in% c("Balkan", "G5")) %>%
      mutate(
        country_name = group,
        economy_type = group
      )
  ) %>%
  mutate(
    label_f = factor(
      country_name,
      levels = c("Balkan", "G5",
                 BALKAN_LABELS[USABLE_BALKAN],
                 G5_LABELS[USABLE_G5])
    )
  ) %>%
  pivot_longer(
    cols      = c(mult_slack, mult_exp),
    names_to  = "state",
    values_to = "multiplier"
  ) %>%
  mutate(
    state = ifelse(state == "mult_slack", "Slack", "Expansion"),
    state = factor(state, levels = c("Slack", "Expansion")),
    horizon_label = paste0(horizon, "Q")
  )

# BUG FIX #7: scale_alpha_manual was keyed by full country-name strings but the
# aesthetic maps to economy_type which only holds "Balkan" / "G5" / individual
# country codes.  Simplified to just the two panel-level keys so the scale
# resolves correctly for every row.
p7 <- ggplot(mult_bar, aes(x = label_f, y = multiplier, fill = state,
                           alpha = economy_type)) +
  geom_col(position = position_dodge(0.75), width = 0.65) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_hline(yintercept = 1, linewidth = 0.7, color = "darkgreen", linetype = "dotted") +
  facet_wrap(~ horizon_label, nrow = 1) +
  scale_fill_manual(values = c("Slack" = COL_SLACK, "Expansion" = COL_EXP)) +
  scale_alpha_manual(
    values = c("Balkan" = 1.00, "G5" = 1.00),
    guide  = "none"
  ) +
  labs(
    title    = "Figure 7: Cumulative Multipliers — Balkan vs G5",
    subtitle = "Green dotted = multiplier of 1.0 | Ramey-Zubairy (2018) IV scaling",
    x = NULL, y = "Cumulative Multiplier", fill = "State"
  ) + theme_paper +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7.5))

# ── Figure 8: Multiplier heatmap — all countries ──────────────────────────────

heat_data <- mult_countries %>%
  bind_rows(
    mult_panels %>%
      filter(group %in% c("Balkan", "G5")) %>%
      mutate(country_name = group, economy_type = group)
  ) %>%
  mutate(
    label_f = factor(
      country_name,
      levels = rev(c("Balkan", "G5",
                     BALKAN_LABELS[USABLE_BALKAN],
                     G5_LABELS[USABLE_G5]))
    )
  ) %>%
  pivot_longer(cols = c(mult_slack, mult_exp),
               names_to = "state", values_to = "value") %>%
  mutate(
    state = case_match(state,
                       "mult_slack" ~ "Slack Multiplier",
                       "mult_exp"   ~ "Expansion Multiplier"
    )
  )

p8 <- ggplot(heat_data,
             aes(x = factor(horizon), y = label_f, fill = value)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 2.6,
            color = "white", fontface = "bold") +
  scale_fill_gradient2(
    low = "#B71C1C", mid = "grey90", high = "#1A237E",
    midpoint = 0, name = "Multiplier"
  ) +
  facet_wrap(~ state) +
  labs(
    title    = "Figure 8: Multiplier Heatmap — Balkan & G5 (All Countries)",
    subtitle = "Cumulative GDP response per 1pp spending shock at selected horizons",
    x = "Horizon (quarters)", y = NULL
  ) + theme_paper

# ── Figure 9: State share timeline — Balkan vs G5 side-by-side ───────────────

state_share <- panel_est %>%
  filter(!is.na(slack_lag1)) %>%
  mutate(
    year         = year(date),
    country_label = ALL_LABELS[country],
    economy_type  = ifelse(country %in% G5_COUNTRIES, "G5", "Balkan")
  ) %>%
  group_by(economy_type, country_label, year) %>%
  summarise(pct_slack = mean(slack_lag1) * 100, .groups = "drop")

p9 <- ggplot(state_share, aes(x = year, y = pct_slack, fill = country_label)) +
  geom_col(position = "dodge", alpha = 0.82, width = 0.8) +
  facet_wrap(~ economy_type, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = country_colors, name = NULL) +
  labs(
    title    = "Figure 9: Share of Quarters in Slack State by Year",
    subtitle = "% of quarters per year where lagged unemployment > HP trend",
    x = NULL, y = "% Quarters in Slack"
  ) + theme_paper +
  theme(legend.key.size = unit(0.4, "cm"))

# ── Figure 10: Slack / Expansion multiplier scatter (Balkan vs G5) ────────────
# Scatter of (slack mult at H=8) vs (exp mult at H=8) per country

scatter_data <- mult_countries %>%
  filter(horizon == 8) %>%
  mutate(
    country_label = country_name,
    economy_type  = economy_type
  )

p10 <- ggplot(scatter_data,
              aes(x = mult_exp, y = mult_slack,
                  color = economy_type, label = country_label)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(size = 3.2, max.overlaps = 15) +
  scale_color_manual(
    values = c("Balkan" = COL_BALKAN, "G5" = COL_G5), name = "Economy Type"
  ) +
  labs(
    title    = "Figure 10: Slack vs Expansion Multiplier by Country (H = 8 Quarters)",
    subtitle = "Points above the 45° line: higher multiplier in slack than expansion",
    x = "Expansion Multiplier", y = "Slack Multiplier"
  ) + theme_paper

# ==============================================================================
# 19. SAVE OUTPUTS
# ==============================================================================

out_dir <- file.path("outputs_all", format(Sys.Date(), "%Y%m%d"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("\nSaving outputs to '%s/'\n", out_dir))

fig_specs <- list(
  list(plot = p1,   file = "fig1a_unemployment_balkan.pdf",   w = 11, h = 10),
  list(plot = p1b,  file = "fig1b_unemployment_g5.pdf",        w = 11, h =  7),
  list(plot = p2,   file = "fig2_gov_expenditure_all.pdf",     w = 14, h = 13),
  list(plot = p3,   file = "fig3_irf_balkan_vs_g5.pdf",        w = 11, h =  7),
  list(plot = p4,   file = "fig4_irf_three_panels.pdf",        w = 14, h =  6),
  list(plot = p5,   file = "fig5_irf_by_country_balkan.pdf",   w = 12, h = 12),
  list(plot = p6,   file = "fig6_irf_by_country_g5.pdf",       w = 12, h =  8),
  list(plot = p7,   file = "fig7_multipliers_balkan_g5.pdf",   w = 14, h =  7),
  list(plot = p8,   file = "fig8_multiplier_heatmap_all.pdf",  w = 12, h =  9),
  list(plot = p9,   file = "fig9_state_share_timeline.pdf",    w = 12, h =  8),
  list(plot = p10,  file = "fig10_slack_exp_scatter.pdf",      w =  9, h =  7)
)

invisible(lapply(fig_specs, function(s) {
  ggsave(file.path(out_dir, s$file), s$plot,
         width = s$w, height = s$h, dpi = 150)
  cat(sprintf("  Saved: %s\n", s$file))
}))

# ── CSV outputs ───────────────────────────────────────────────────────────────

write.csv(irf_all,      file.path(out_dir, "irf_results_all.csv"),              row.names = FALSE)
write.csv(irf_balkan,   file.path(out_dir, "irf_balkan_panel.csv"),             row.names = FALSE)
write.csv(irf_g5,       file.path(out_dir, "irf_g5_panel.csv"),                 row.names = FALSE)
write.csv(irf_all_pan,  file.path(out_dir, "irf_full_panel.csv"),               row.names = FALSE)
write.csv(multipliers,  file.path(out_dir, "multipliers_all.csv"),              row.names = FALSE)
write.csv(mult_panels,  file.path(out_dir, "multipliers_panels.csv"),           row.names = FALSE)
write.csv(mult_countries, file.path(out_dir, "multipliers_countries.csv"),      row.names = FALSE)
write.csv(mult_rob,     file.path(out_dir, "robustness_multipliers.csv"),       row.names = FALSE)
write.csv(wald_table,   file.path(out_dir, "hypothesis_tests.csv"),             row.names = FALSE)
write.csv(sumstats,     file.path(out_dir, "summary_statistics.csv"),           row.names = FALSE)
write.csv(coverage,     file.path(out_dir, "data_coverage_all.csv"),            row.names = FALSE)
write.csv(adf_table,    file.path(out_dir, "adf_unit_root_tests.csv"),          row.names = FALSE)
write.csv(fs_diag,      file.path(out_dir, "first_stage_fstats.csv"),           row.names = FALSE)

cat("\n=== Done. All results saved. ===\n")

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "session_info.txt"))