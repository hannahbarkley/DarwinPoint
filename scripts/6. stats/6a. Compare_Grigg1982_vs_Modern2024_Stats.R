# ==============================================================================
# Script:  6a. Compare_Grigg1982_vs_Modern2024_Stats.R
# Purpose: Statistical evaluation of decadal shifts (1982 Grigg vs. 2024 P90)
#          treating 1982 as fixed historical constants:
#            1. Island-by-Island Wald z-tests and 95% bootstrap CIs against the 
#               1982 constant.
#            2. Archipelago-wide stratified bootstrap test (B = 10,000) for the
#               pooled mean decadal offset (Delta) and empirical p-value vs. 0.
# Inputs:  - data/grigg data/HA_1982_Grigg.csv
#          - data/grigg data/HA_2024_Grigg_Equivalent_Sites.csv
# Outputs: - outputs/stats/HA_6a_Grigg1982_vs_Modern2024_Stats.csv
#          - outputs/stats/HA_6a_Grigg1982_vs_Modern2024_Stats_Report.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

out_dir <- here("outputs", "stats")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. HELPERS & INGESTION
# ------------------------------------------------------------------------------
island_lookup <- c(
  "kure"                         = "H\u014dlanik\u016b",
  "kure atoll"                   = "H\u014dlanik\u016b",
  "midway"                       = "Kuaihelani",
  "midway atoll"                 = "Kuaihelani",
  "pearl & hermes"               = "Manawai",
  "pearl and hermes"             = "Manawai",
  "pearl & hermes atoll"         = "Manawai",
  "pearl_hermes"                 = "Manawai",
  "lisianski"                    = "Kapou",
  "lisianski island"             = "Kapou",
  "laysan"                       = "Kamole",
  "laysan island"                = "Kamole",
  "maro"                         = "Kamokuokamohoali'i",
  "maro reef"                    = "Kamokuokamohoali'i",
  "gardner"                      = "P\u016bh\u0101honu",
  "gardner pinnacles"            = "P\u016bh\u0101honu",
  "gardner_pinnacles"            = "P\u016bh\u0101honu",
  "french frigate"               = "Lalo",
  "french frigate shoals"        = "Lalo",
  "french_frigate"               = "Lalo",
  "french_frigate_shoals"        = "Lalo",
  "ffs"                          = "Lalo",
  "necker"                       = "Mokumanamana",
  "necker island"                = "Mokumanamana",
  "nihoa"                        = "Nihoa",
  "nihoa island"                 = "Nihoa",
  "kauai"                        = "Kaua'i",
  "niihau"                       = "Ni'ihau",
  "oahu"                         = "O'ahu",
  "molokai"                      = "Moloka'i",
  "lanai"                        = "L\u0101na'i",
  "maui"                         = "Maui",
  "kahoolawe"                    = "Kaho'olawe",
  "hawaii"                       = "Hawai'i"
)
island_lookup <- island_lookup[!duplicated(names(island_lookup))]

format_islands <- function(x) {
  clean_key <- tolower(trimws(as.character(x)))
  clean_key <- gsub("_", " ", clean_key)
  mapped    <- island_lookup[clean_key]
  coalesce(unname(mapped), as.character(x))
}

saferead_csv <- function(candidates) {
  search_dirs <- c(
    here("data", "grigg data"),
    here("data", "strs data", "processed data"),
    here("data", "strs data", "field data"),
    here("data"),
    here()
  )
  for (f in candidates) {
    for (d in search_dirs) {
      target <- file.path(d, f)
      if (file.exists(target)) return(read_csv(target, show_col_types = FALSE))
    }
  }
  stop(paste("None of candidate files found:", paste(candidates, collapse = " | ")))
}

# ------------------------------------------------------------------------------
# 2. DATA INGESTION & HARMONIZATION
# ------------------------------------------------------------------------------
df_grigg <- saferead_csv("HA_1982_Grigg.csv") %>%
  mutate(
    PLOT_LABEL = format_islands(LOCATION),
    GRIGG_COV  = coalesce(!!!select(., any_of(c("CORAL_COVER_PERC", "HARD_CORAL_COVER_PCT")))),
    GRIGG_PROD = coalesce(!!!select(., any_of(c("GROSS_PROD", "GROSS_CARB_PROD_KG_M2_YR_MEAN", "PROD"))))
  ) %>%
  select(PLOT_LABEL, LATITUDE, GRIGG_COV, GRIGG_PROD) %>%
  distinct(PLOT_LABEL, .keep_all = TRUE)

df_strs_sw <- saferead_csv(c(
  "HA_2024_Grigg_Equivalent_Sites.csv", 
  "HA_2024_StRS_Net_SiteLevel_GriggSubset.csv"
)) %>%
  mutate(
    LOCATION            = coalesce(!!!select(., any_of(c("LOCATION", "ISLAND")))),
    PLOT_LABEL          = format_islands(LOCATION),
    CORAL_COVER_PERC    = coalesce(!!!select(., any_of(c("CORAL_COVER_PERC", "HARD_CORAL_COVER_PCT_MEAN", "HARD_CORAL_COVER_PCT")))),
    GROSS_PROD_KG_M2_YR = coalesce(!!!select(., any_of(c("GROSS_PROD_KG_M2_YR", "GROSS_CARB_PROD_KG_M2_YR_MEAN", "PROD"))))
  )

# Filter modern data to islands shared with Grigg
shared_islands <- intersect(unique(df_grigg$PLOT_LABEL), unique(df_strs_sw$PLOT_LABEL))
df_grigg_sub   <- df_grigg %>% filter(PLOT_LABEL %in% shared_islands)
df_sw_sub      <- df_strs_sw %>% filter(PLOT_LABEL %in% shared_islands)

n_islands <- length(shared_islands)

# ------------------------------------------------------------------------------
# 3. ISLAND-LEVEL BOOTSTRAP QUANTIFICATION
# ------------------------------------------------------------------------------
set.seed(42)
B <- 10000

# Compute modern 90th percentile, bootstrap SE, and 95% bootstrap percentile CIs
island_boot_summary <- df_sw_sub %>%
  group_by(PLOT_LABEL) %>%
  summarise(
    N_Sites = n(),
    
    # Coral Cover
    MODERN_COV_P90     = quantile(CORAL_COVER_PERC, 0.90, na.rm = TRUE),
    MODERN_COV_SE      = sd(replicate(B, quantile(sample(na.omit(CORAL_COVER_PERC), replace = TRUE), 0.90))),
    MODERN_COV_CI_LOW  = quantile(replicate(B, quantile(sample(na.omit(CORAL_COVER_PERC), replace = TRUE), 0.90)), 0.025),
    MODERN_COV_CI_HIGH = quantile(replicate(B, quantile(sample(na.omit(CORAL_COVER_PERC), replace = TRUE), 0.90)), 0.975),
    
    # Gross Production
    MODERN_PROD_P90     = quantile(GROSS_PROD_KG_M2_YR, 0.90, na.rm = TRUE),
    MODERN_PROD_SE      = sd(replicate(B, quantile(sample(na.omit(GROSS_PROD_KG_M2_YR), replace = TRUE), 0.90))),
    MODERN_PROD_CI_LOW  = quantile(replicate(B, quantile(sample(na.omit(GROSS_PROD_KG_M2_YR), replace = TRUE), 0.90)), 0.025),
    MODERN_PROD_CI_HIGH = quantile(replicate(B, quantile(sample(na.omit(GROSS_PROD_KG_M2_YR), replace = TRUE), 0.90)), 0.975),
    .groups = "drop"
  )

# Join with 1982 fixed constants & compute island Wald z-tests against the constant
island_stats <- df_grigg_sub %>%
  inner_join(island_boot_summary, by = "PLOT_LABEL") %>%
  mutate(
    # Coral Cover vs. Constant
    DELTA_COV       = MODERN_COV_P90 - GRIGG_COV,
    Z_COV           = DELTA_COV / MODERN_COV_SE,
    P_COV_ISLAND    = 2 * (1 - pnorm(abs(Z_COV))),
    COV_SIG_DIFF    = (GRIGG_COV < MODERN_COV_CI_LOW) | (GRIGG_COV > MODERN_COV_CI_HIGH),
    
    # Gross Production vs. Constant
    DELTA_PROD      = MODERN_PROD_P90 - GRIGG_PROD,
    Z_PROD          = DELTA_PROD / MODERN_PROD_SE,
    P_PROD_ISLAND   = 2 * (1 - pnorm(abs(Z_PROD))),
    PROD_SIG_DIFF   = (GRIGG_PROD < MODERN_PROD_CI_LOW) | (GRIGG_PROD > MODERN_PROD_CI_HIGH)
  ) %>%
  arrange(desc(LATITUDE))

# ------------------------------------------------------------------------------
# 4. ARCHIPELAGO-WIDE STRATIFIED BOOTSTRAP TEST (vs. FIXED CONSTANTS)
# ------------------------------------------------------------------------------
cat(sprintf("Running stratified bootstrap (B = %d) across %d islands...\n", B, n_islands))

# Extract split vectors for fast bootstrapping
sites_by_island_cov  <- split(df_sw_sub$CORAL_COVER_PERC, df_sw_sub$PLOT_LABEL)
sites_by_island_prod <- split(df_sw_sub$GROSS_PROD_KG_M2_YR, df_sw_sub$PLOT_LABEL)

# Match order of constants to split list names
constants_cov  <- setNames(island_stats$GRIGG_COV, island_stats$PLOT_LABEL)[names(sites_by_island_cov)]
constants_prod <- setNames(island_stats$GRIGG_PROD, island_stats$PLOT_LABEL)[names(sites_by_island_prod)]

boot_mean_delta_cov  <- numeric(B)
boot_mean_delta_prod <- numeric(B)

for (b in seq_len(B)) {
  # Resample within each island and calculate P90
  p90_cov_b  <- vapply(sites_by_island_cov,  function(x) quantile(sample(na.omit(x), replace = TRUE), 0.90), numeric(1))
  p90_prod_b <- vapply(sites_by_island_prod, function(x) quantile(sample(na.omit(x), replace = TRUE), 0.90), numeric(1))
  
  # Subtract fixed 1982 constants and compute archipelago mean offset
  boot_mean_delta_cov[b]  <- mean(p90_cov_b - constants_cov)
  boot_mean_delta_prod[b] <- mean(p90_prod_b - constants_prod)
}

# Empirical Bootstrap Metrics
mean_d_cov   <- mean(island_stats$DELTA_COV)
se_boot_cov  <- sd(boot_mean_delta_cov)
ci_boot_cov  <- quantile(boot_mean_delta_cov, probs = c(0.025, 0.975))
p_boot_cov   <- 2 * min(mean(boot_mean_delta_cov <= 0), mean(boot_mean_delta_cov >= 0))

mean_d_prod  <- mean(island_stats$DELTA_PROD)
se_boot_prod <- sd(boot_mean_delta_prod)
ci_boot_prod <- quantile(boot_mean_delta_prod, probs = c(0.025, 0.975))
p_boot_prod  <- 2 * min(mean(boot_mean_delta_prod <= 0), mean(boot_mean_delta_prod >= 0))

# ------------------------------------------------------------------------------
# 5. SUMMARY REPORT & EXPORT
# ------------------------------------------------------------------------------
report_lines <- c(
  "================================================================================",
  " SCRIPT 6a: BOOTSTRAP HYPOTHESIS TESTING (MODERN P90 vs. 1982 CONSTANTS)",
  "================================================================================",
  sprintf("Number of compared shared islands: %d", n_islands),
  sprintf("Bootstrap iterations:              %d", B),
  "",
  "--- 1. ARCHIPELAGO-WIDE STRATIFIED BOOTSTRAP RESULTS ---",
  "[ Hard Coral Cover (%) ]",
  sprintf("  • Mean Decadal Shift (Delta):        %+5.2f%%", mean_d_cov),
  sprintf("  • Stratified Bootstrap SE:           %.2f%%", se_boot_cov),
  sprintf("  • 95%% Bootstrap Percentile CI:       [%.2f%%, %.2f%%]", ci_boot_cov[1], ci_boot_cov[2]),
  sprintf("  • Empirical Bootstrap p-value:       p = %.4f %s\n", 
          p_boot_cov, if (p_boot_cov < 0.05) "(Significant)" else "(Not Significant)"),
  
  "[ Gross Carbonate Production (kg m^-2 yr^-1) ]",
  sprintf("  • Mean Decadal Shift (Delta):        %+5.2f kg m^-2 yr^-1", mean_d_prod),
  sprintf("  • Stratified Bootstrap SE:           %.2f kg m^-2 yr^-1", se_boot_prod),
  sprintf("  • 95%% Bootstrap Percentile CI:       [%.2f, %.2f] kg m^-2 yr^-1", ci_boot_prod[1], ci_boot_prod[2]),
  sprintf("  • Empirical Bootstrap p-value:       p = %.4f %s\n", 
          p_boot_prod, if (p_boot_prod < 0.05) "(Significant)" else "(Not Significant)"),
  
  "--- 2. ISLAND-BY-ISLAND TESTS AGAINST FIXED 1982 CONSTANTS ---",
  sprintf("Coral Cover: %d of %d islands significantly differ from 1982 constant (p < 0.05)",
          sum(island_stats$COV_SIG_DIFF), n_islands),
  sprintf("Gross Prod:  %d of %d islands significantly differ from 1982 constant (p < 0.05)",
          sum(island_stats$PROD_SIG_DIFF), n_islands),
  "================================================================================"
)

cat(paste(report_lines, collapse = "\n"), "\n\n")

# Display island-level results table
print(
  island_stats %>%
    select(PLOT_LABEL, GRIGG_COV, MODERN_COV_P90, MODERN_COV_SE, DELTA_COV, Z_COV, P_COV_ISLAND,
           GRIGG_PROD, MODERN_PROD_P90, MODERN_PROD_SE, DELTA_PROD, Z_PROD, P_PROD_ISLAND),
  n = Inf
)

# Export files
out_csv <- file.path(out_dir, "HA_6a_Grigg1982_vs_Modern2024_Stats.csv")
out_txt <- file.path(out_dir, "HA_6a_Grigg1982_vs_Modern2024_Stats_Report.txt")

write_csv(island_stats, out_csv)
writeLines(report_lines, out_txt)

cat(sprintf("\n✅ Island-level test data exported to: %s\n", out_csv))
cat(sprintf("✅ Statistical report saved to:       %s\n", out_txt))

# ------------------------------------------------------------------------------
# 6. MANUSCRIPT TEXT PHRASING
# ------------------------------------------------------------------------------
fmt_p <- function(p, B = 10000) {
  if (is.na(p) || p < 1 / B) sprintf("p < %.4f", 1 / B) else sprintf("p = %.4f", p)
}

p_txt_cov  <- fmt_p(p_boot_cov, B)
p_txt_prod <- fmt_p(p_boot_prod, B)

cat("\n=== MANUSCRIPT TEXT PHRASING ===\n")
cat(sprintf(
  "When evaluating modern 90th-percentile southwestern fore-reef estimates against historical 1982 island constants, stratified bootstrap hypothesis testing (B = %s) revealed an archipelago-wide shift of %+.1f ± %.1f%% in coral cover (mean ± bootstrap SE; 95%% CI: [%.1f%%, %.1f%%]; bootstrap %s) and %+.2f ± %.2f kg m⁻² yr⁻¹ in gross carbonate production (95%% CI: [%.2f, %.2f]; bootstrap %s).\n",
  format(B, big.mark = ","),
  mean_d_cov, se_boot_cov, ci_boot_cov[1], ci_boot_cov[2], p_txt_cov,
  mean_d_prod, se_boot_prod, ci_boot_prod[1], ci_boot_prod[2], p_txt_prod
))