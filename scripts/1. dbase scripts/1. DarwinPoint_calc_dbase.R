# ==============================================================================
# Script:  1. DarwinPoint_calc_dbase.R
# Purpose: Master wrapper to sequentially execute the 3-stage reference database 
#          generation pipeline in 'scripts/1.dbase scripts':
#            Stage 1: Empirical aggregation, sister-taxa imputation & coordinate lookup
#                     via 1a. DarwinPoint_calc_dbase_ref.R
#            Stage 2: Island crossing, Grigg (1982) scaling & CCA rate integration
#                     via 1b. DarwinPoint_calc_dbase_scale.R
#            Stage 3: Growth curve regression coefficient and intercept calculation
#                     via 1c. DarwinPoint_calc_dbase_coeffs.R
# Outputs: - data/dbase/prod_dbase_ha_ref.csv          (Stage 1)
#          - data/dbase/prod_dbase_ha_fixed.csv        (Stage 2)
#          - data/dbase/prod_dbase_ha_strs.csv         (Stage 2)
#          - data/dbase/prod_dbase_ha_strs_coeffs.csv  (Stage 3)
# ==============================================================================

library(here)

# ------------------------------------------------------------------------------
# 1. DIRECT PATH RESOLUTION
# ------------------------------------------------------------------------------

scripts_dir <- if (dir.exists(here("scripts", "1.dbase scripts"))) {
  here("scripts", "1.dbase scripts")
} else {
  here("scripts", "1. dbase scripts")
}

script_1a <- file.path(scripts_dir, "1a. DarwinPoint_calc_dbase_ref.R")
script_1b <- file.path(scripts_dir, "1b. DarwinPoint_calc_dbase_scale.R")
script_1c <- file.path(scripts_dir, "1c. DarwinPoint_calc_dbase_coeffs.R")

stopifnot("Stage 1 script not found!" = file.exists(script_1a))
stopifnot("Stage 2 script not found!" = file.exists(script_1b))
stopifnot("Stage 3 script not found!" = file.exists(script_1c))

# Direct target spreadsheet outputs
dbase_dir         <- here("data", "dbase")
stage1_out        <- file.path(dbase_dir, "prod_dbase_ha_ref.csv")
stage2_fixed_out  <- file.path(dbase_dir, "prod_dbase_ha_fixed.csv")
stage2_strs_out   <- file.path(dbase_dir, "prod_dbase_ha_strs.csv")
stage3_coeffs_out <- file.path(dbase_dir, "prod_dbase_ha_strs_coeffs.csv")

cat("\n==============================================================================\n")
cat("🚀 STARTING DARWIN POINT MASTER PIPELINE (3 STAGES)\n")
cat("==============================================================================\n")
cat(sprintf("Scripts Directory : %s\n", scripts_dir))
cat(sprintf("Stage 1 Script    : %s\n", basename(script_1a)))
cat(sprintf("Stage 2 Script    : %s\n", basename(script_1b)))
cat(sprintf("Stage 3 Script    : %s\n", basename(script_1c)))

pipeline_start_time <- Sys.time()

# ------------------------------------------------------------------------------
# 2. STAGE 1: EMPIRICAL HARMONIZATION & PROXY IMPUTATION
# ------------------------------------------------------------------------------

cat("\n------------------------------------------------------------------------------\n")
cat(sprintf("⏳ STAGE 1: Running %s ...\n", basename(script_1a)))
cat("------------------------------------------------------------------------------\n")

stage1_start <- Sys.time()
stage1_env   <- new.env(parent = globalenv())

tryCatch(
  {
    source(script_1a, local = stage1_env)
    stage1_duration <- round(difftime(Sys.time(), stage1_start, units = "secs"), 2)
    cat(sprintf("✅ Stage 1 completed in %s seconds.\n", stage1_duration))
  },
  error = function(e) {
    stop(sprintf("\n❌ Pipeline aborted at Stage 1:\n%s", e$message), call. = FALSE)
  }
)

if (!file.exists(stage1_out) || file.info(stage1_out)$size == 0) {
  stop(sprintf("\n❌ Pipeline aborted: Target output '%s' was not created or is empty.", stage1_out), call. = FALSE)
}
cat(sprintf("📁 Stage 1 output verified: %s (%d bytes)\n", stage1_out, file.info(stage1_out)$size))

# ------------------------------------------------------------------------------
# 3. STAGE 2: TEMPLATE GENERATION, GRIGG SCALING & CCA INTEGRATION
# ------------------------------------------------------------------------------

cat("\n------------------------------------------------------------------------------\n")
cat(sprintf("⏳ STAGE 2: Running %s ...\n", basename(script_1b)))
cat("------------------------------------------------------------------------------\n")

stage2_start <- Sys.time()
stage2_env   <- new.env(parent = globalenv())

tryCatch(
  {
    source(script_1b, local = stage2_env)
    stage2_duration <- round(difftime(Sys.time(), stage2_start, units = "secs"), 2)
    cat(sprintf("✅ Stage 2 completed in %s seconds.\n", stage2_duration))
  },
  error = function(e) {
    stop(sprintf("\n❌ Pipeline aborted at Stage 2:\n%s", e$message), call. = FALSE)
  }
)

if (!file.exists(stage2_fixed_out) || !file.exists(stage2_strs_out)) {
  stop(sprintf("\n❌ Pipeline aborted: Stage 2 failed to produce '%s' and/or '%s'.", stage2_fixed_out, stage2_strs_out), call. = FALSE)
}
cat(sprintf("📁 Stage 2 fixed output verified: %s (%d bytes)\n", stage2_fixed_out, file.info(stage2_fixed_out)$size))
cat(sprintf("📁 Stage 2 strs output verified : %s (%d bytes)\n", stage2_strs_out, file.info(stage2_strs_out)$size))

# ------------------------------------------------------------------------------
# 4. STAGE 3: REGRESSION COEFFICIENT & INTERCEPT CALCULATION
# ------------------------------------------------------------------------------

cat("\n------------------------------------------------------------------------------\n")
cat(sprintf("⏳ STAGE 3: Running %s ...\n", basename(script_1c)))
cat("------------------------------------------------------------------------------\n")

stage3_start <- Sys.time()
stage3_env   <- new.env(parent = globalenv())

tryCatch(
  {
    source(script_1c, local = stage3_env)
    stage3_duration <- round(difftime(Sys.time(), stage3_start, units = "secs"), 2)
    cat(sprintf("✅ Stage 3 completed in %s seconds.\n", stage3_duration))
  },
  error = function(e) {
    stop(sprintf("\n❌ Pipeline aborted at Stage 3:\n%s", e$message), call. = FALSE)
  }
)

if (!file.exists(stage3_coeffs_out) || file.info(stage3_coeffs_out)$size == 0) {
  stop(sprintf("\n❌ Pipeline aborted: Stage 3 output '%s' was not created or is empty.", stage3_coeffs_out), call. = FALSE)
}
cat(sprintf("📁 Stage 3 coeffs output verified: %s (%d bytes)\n", stage3_coeffs_out, file.info(stage3_coeffs_out)$size))

# ------------------------------------------------------------------------------
# 5. PIPELINE SUMMARY
# ------------------------------------------------------------------------------

total_duration <- round(difftime(Sys.time(), pipeline_start_time, units = "secs"), 2)

cat("\n==============================================================================\n")
cat("🎉 FULL DARWIN POINT PIPELINE EXECUTION SUCCESSFUL\n")
cat("==============================================================================\n")
cat(sprintf("Total Elapsed Time : %s seconds\n", total_duration))
cat(sprintf("Stage 1 Master Ref : %s\n", stage1_out))
cat(sprintf("Stage 2 Fixed Dbase: %s\n", stage2_fixed_out))
cat(sprintf("Stage 2 STRS Dbase : %s\n", stage2_strs_out))
cat(sprintf("Stage 3 STRS Coeffs: %s\n", stage3_coeffs_out))
cat("==============================================================================\n\n")