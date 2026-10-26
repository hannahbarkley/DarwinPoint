# ==============================================================================
# Script:  6c. DarwinPoint_fixed_SfM_vs_CoralNet.R
# Purpose: Statistical Parity and Methodological Validation of SfM vs. CoralNet
# Inputs:  data/fixed data/processed data/HA_2024_Fixed_Net_SiteLevel_sfm.csv
#            (or HA_2024_Fixed_Net_SiteLevel_cnet.csv / HA_2024_Fixed_Net_SiteLevel.csv)
#          data/strs data/field data/HA_2024_StRS_Prod_SiteLevel.csv
# Outputs: outputs/tables/stats_sfm_vs_coralnet_parity.txt
#          outputs/tables/stats_sfm_vs_coralnet_parity.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(car) # for linearHypothesis()
})

out_dir <- here("outputs", "tables")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

saferead_csv <- function(candidates) {
  search_dirs <- c(
    here("data", "fixed data", "processed data"),
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
  stop(paste("None of the candidate files found:", paste(candidates, collapse = " | ")))
}

# ------------------------------------------------------------------------------
# 1. Load and Pair Datasets
# ------------------------------------------------------------------------------
df_fixed_net <- saferead_csv(c(
  "HA_2024_Fixed_Net_SiteLevel_sfm.csv",
  "HA_2024_Fixed_Net_SiteLevel_cnet.csv",
  "HA_2024_Fixed_Net_SiteLevel.csv"
)) %>%
  mutate(
    SITE_ID   = coalesce(if ("OCC_SITEID" %in% names(.)) OCC_SITEID else NULL,
                         if ("SITE" %in% names(.)) SITE else NULL),
    SFM_PROD  = as.numeric(coalesce(if ("GROSS_CARB_PROD_KG_M2_YR_MEAN" %in% names(.)) GROSS_CARB_PROD_KG_M2_YR_MEAN else NULL,
                                    if ("PROD" %in% names(.)) PROD else NULL,
                                    NA_real_))
  )

df_coralnet_prod <- saferead_csv("HA_2024_StRS_Prod_SiteLevel.csv") %>%
  transmute(
    SITE_ID       = SITE,
    CORALNET_PROD = as.numeric(PROD)
  )

matched_data <- df_fixed_net %>%
  select(SITE_ID, SFM_PROD) %>%
  inner_join(df_coralnet_prod, by = "SITE_ID") %>%
  drop_na(SFM_PROD, CORALNET_PROD)

n_sites <- nrow(matched_data)
cat(sprintf("Loaded and matched %d survey sites.\n\n", n_sites))

# ------------------------------------------------------------------------------
# 2. Paired t-Test (Mean Bias)
# ------------------------------------------------------------------------------
diff_prod <- matched_data$CORALNET_PROD - matched_data$SFM_PROD
mean_diff <- mean(diff_prod)
sd_diff   <- sd(diff_prod)

ttest_res <- t.test(matched_data$CORALNET_PROD, matched_data$SFM_PROD, paired = TRUE)

# ------------------------------------------------------------------------------
# 3. Wilcoxon Signed-Rank Test (Non-parametric Bias)
# ------------------------------------------------------------------------------
wilcox_res <- wilcox.test(matched_data$CORALNET_PROD, matched_data$SFM_PROD, paired = TRUE, exact = FALSE)

# ------------------------------------------------------------------------------
# 4. Ordinary Least Squares (OLS) Linear Regression & 95% CIs
# ------------------------------------------------------------------------------
fit_lm  <- lm(CORALNET_PROD ~ SFM_PROD, data = matched_data)
sum_lm  <- summary(fit_lm)
ci_lm   <- confint(fit_lm, level = 0.95)

intercept_est <- coef(fit_lm)[1]
slope_est     <- coef(fit_lm)[2]
r_squared     <- sum_lm$r.squared

# ------------------------------------------------------------------------------
# 5. Joint Hypothesis Test (Simultaneous Parity: Slope = 1 and Intercept = 0)
# ------------------------------------------------------------------------------
joint_test <- linearHypothesis(fit_lm, c("(Intercept) = 0", "SFM_PROD = 1"))

f_stat  <- joint_test$F[2]
df_num  <- joint_test$Df[2]
df_den  <- joint_test$Res.Df[2]; p_val_f <- joint_test$`Pr(>F)`[2]

# ------------------------------------------------------------------------------
# 6. Console Report & Result Export
# ------------------------------------------------------------------------------
report_text <- sprintf(
  "================================================================================
STATISTICAL VALIDATION REPORT: SfM vs. CoralNet Gross Production
================================================================================
Sample size: %d matched sites

1. Paired Differences (CoralNet - SfM):
   - Mean Difference: %+.3f +/- %.3f kg m^-2 yr^-1 (mean +/- SD)
   - Paired t-test:   t(%d) = %.3f, p = %.4f
   - 95%% CI:          [%.3f, %.3f]

2. Non-Parametric Signed-Rank Test:
   - Wilcoxon:        W = %.1f, p = %.4f

3. OLS Linear Model (CoralNet = beta_1 * SfM + beta_0):
   - Slope:           %.3f (95%% CI: [%.3f, %.3f])
   - Intercept:       %.3f (95%% CI: [%.3f, %.3f])
   - R^2:             %.3f

4. Joint Hypothesis Test of 1:1 Parity (H0: beta_0 = 0 and beta_1 = 1):
   - F(%d, %d) = %.3f, p = %.4f
================================================================================
",
  n_sites,
  mean_diff, sd_diff,
  ttest_res$parameter, ttest_res$statistic, ttest_res$p.value,
  ttest_res$conf.int[1], ttest_res$conf.int[2],
  wilcox_res$statistic, wilcox_res$p.value,
  slope_est, ci_lm["SFM_PROD", 1], ci_lm["SFM_PROD", 2],
  intercept_est, ci_lm["(Intercept)", 1], ci_lm["(Intercept)", 2],
  r_squared,
  df_num, df_den, f_stat, p_val_f
)

cat(report_text)
writeLines(report_text, file.path(out_dir, "stats_sfm_vs_coralnet_parity.txt"))

stats_summary_table <- tibble(
  metric = c("N", "Mean_Diff", "SD_Diff", "t_stat", "t_df", "t_pvalue", 
             "t_ci_lower", "t_ci_upper", "wilcox_W", "wilcox_pvalue",
             "OLS_slope", "slope_ci_lower", "slope_ci_upper",
             "OLS_intercept", "intercept_ci_lower", "intercept_ci_upper",
             "R2", "F_stat", "F_df_num", "F_df_den", "F_pvalue"),
  value = c(n_sites, mean_diff, sd_diff, ttest_res$statistic, ttest_res$parameter, ttest_res$p.value,
            ttest_res$conf.int[1], ttest_res$conf.int[2], wilcox_res$statistic, wilcox_res$p.value,
            slope_est, ci_lm["SFM_PROD", 1], ci_lm["SFM_PROD", 2],
            intercept_est, ci_lm["(Intercept)", 1], ci_lm["(Intercept)", 2],
            r_squared, f_stat, df_num, df_den, p_val_f)
)

write_csv(stats_summary_table, file.path(out_dir, "stats_sfm_vs_coralnet_parity.csv"))
cat("✅ Analysis complete. Summary table and log saved to outputs/tables/\n")