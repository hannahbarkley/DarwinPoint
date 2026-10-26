# ==============================================================================
# Script:  0. DarwinPoint_run_all.R
# Purpose: Master End-to-End Orchestrator for the Hawaiian Archipelago Carbonate
#          Budget Pipeline ("Darwin Point"). Executes all 4 project stages:
#            Stage 1: Reference Database Generation (1a, 1b, 1c)
#            Stage 2: Fixed OCC Site Net Budgets (2a, 2b)
#            Stage 3: StRS Regional Budgets & Survey Weighting (3a, 3b, 3c)
#            Stage 4: Grigg Fore-Reef Subsetting & P90 Modeling (4a, 4b, 4c)
# Outputs: Complete suite of calibrated databases, site-level net budgets,
#          survey design RDS objects, island-level summaries, and figures.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(crayon)
  library(tibble)
  library(dplyr)
})

# ==============================================================================
# 1. GLOBAL PIPELINE CONFIGURATION
# ==============================================================================

master_start_time <- Sys.time()

region <- "HA"
year   <- "2024"

base_dir <- here()
scripts_base <- file.path(base_dir, "scripts")

# Resolve stage subdirectories across standard directory patterns
find_stage_dir <- function(pattern, fallback) {
  candidates <- list.dirs(scripts_base, recursive = FALSE, full.names = TRUE)
  match <- grep(pattern, candidates, value = TRUE, ignore.case = TRUE)
  if (length(match) > 0) match[1] else file.path(scripts_base, fallback)
}

stage_dirs <- list(
  stage1 = find_stage_dir("1.*dbase", "1. dbase scripts"),
  stage2 = find_stage_dir("2.*fixed", "2. fixed scripts"),
  stage3 = find_stage_dir("3.*strs",  "3. strs scripts"),
  stage4 = find_stage_dir("4.*grigg", "4. grigg scripts")
)

# Output directory targets
data_dir   <- file.path(base_dir, "data")
dbase_dir  <- file.path(data_dir, "dbase")
fixed_dir  <- file.path(data_dir, "fixed data", "processed data")
strs_dir   <- file.path(data_dir, "strs data", "processed data")
grigg_dir  <- file.path(data_dir, "grigg data")
plot_dir   <- file.path(base_dir, "outputs", "supporting figures")

# Resolve stage wrapper scripts
find_stage_wrapper <- function(dir_path, filenames) {
  for (f in filenames) {
    p <- file.path(dir_path, f)
    if (file.exists(p)) return(p)
    p_root <- file.path(base_dir, f)
    if (file.exists(p_root)) return(p_root)
  }
  return(file.path(dir_path, filenames[1]))
}

pipeline_stages <- list(
  stage1 = list(
    id          = "Stage 1",
    name        = "Reference Database Calibration",
    path        = find_stage_wrapper(stage_dirs$stage1, c("1. DarwinPoint_calc_dbase.R", "1. DarwinPoint_calc_dbase_2.R")),
    dir         = stage_dirs$stage1,
    outputs     = c(
      file.path(dbase_dir, "prod_dbase_ha_ref.csv"),
      file.path(dbase_dir, "prod_dbase_ha_fixed.csv"),
      file.path(dbase_dir, "prod_dbase_ha_strs.csv"),
      file.path(dbase_dir, "prod_dbase_ha_strs_coeffs.csv")
    )
  ),
  stage2 = list(
    id          = "Stage 2",
    name        = "Fixed-Site Carbonate Budgets",
    path        = find_stage_wrapper(stage_dirs$stage2, c("2. DarwinPoint_process_fixed.R")),
    dir         = stage_dirs$stage2,
    outputs     = c(
      file.path(fixed_dir, sprintf("%s_%s_Fixed_Urchins_SiteLevel.csv", region, year)),
      file.path(fixed_dir, sprintf("%s_%s_Fixed_Fish_SiteLevel.csv", region, year)),
      file.path(fixed_dir, sprintf("%s_%s_Fixed_Net_SiteLevel_sfm.csv", region, year)),
      file.path(fixed_dir, sprintf("%s_%s_Fixed_Net_SiteLevel_cnet.csv", region, year))
    )
  ),
  stage3 = list(
    id          = "Stage 3",
    name        = "StRS Surveys, Porosity & Island Scaling",
    path        = find_stage_wrapper(stage_dirs$stage3, c("3. DarwinPoint_process_strs.R")),
    dir         = stage_dirs$stage3,
    outputs     = c(
      file.path(strs_dir, sprintf("%s_%s_StRS_Fish_SiteLevel.csv", region, year)),
      file.path(strs_dir, sprintf("%s_%s_StRS_Net_SiteLevel_cnet.csv", region, year)),
      file.path(strs_dir, sprintf("%s_%s_StRS_Porosity.csv", region, year)),
      file.path(strs_dir, sprintf("%s_%s_StRS_SiteLevel_Weighted.csv", region, year)),
      file.path(strs_dir, sprintf("%s_%s_StRS_SurveyDesign.rds", region, year)),
      file.path(strs_dir, sprintf("%s_%s_StRS_IslandLevel.csv", region, year))
    )
  ),
  stage4 = list(
    id          = "Stage 4",
    name        = "Grigg Fore-Reef Subsetting & P90 Latitudinal Modeling",
    path        = find_stage_wrapper(stage_dirs$stage4, c("4. DarwinPoint_grigg_comp.R", "4. DarwinPoint_process_grigg.R", "4. DarwinPoint_grigg.R")),
    dir         = stage_dirs$stage4,
    outputs     = c(
      file.path(grigg_dir, sprintf("%s_%s_Grigg_Equivalent_Sites.csv", region, year)),
      file.path(grigg_dir, sprintf("%s_%s_Grigg_IslandLevel_P90.csv", region, year))
    )
  )
)

# ==============================================================================
# 2. PRE-FLIGHT VERIFICATION
# ==============================================================================

cat(bold(cyan("\n==============================================================================\n")))
cat(bold(cyan(sprintf(" 🌺 MASTER END-TO-END PIPELINE: HAWAIIAN ARCHIPELAGO CARBONATE BUDGET (%s)\n", year))))
cat(bold(cyan("==============================================================================\n\n")))

cat(blue("🔍 Running pre-flight system checks across all 4 stages...\n"))

# Ensure base destination folders exist
invisible(lapply(c(dbase_dir, fixed_dir, strs_dir, grigg_dir, plot_dir), function(d) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}))

missing_wrappers <- character(0)
for (st in pipeline_stages) {
  if (!file.exists(st$path)) {
    missing_wrappers <- c(missing_wrappers, sprintf("%s: %s (Expected at: %s)", st$id, basename(st$path), st$path))
  }
}

if (length(missing_wrappers) > 0) {
  cat(red("\n❌ PRE-FLIGHT ABORTED: Missing stage wrapper scripts:\n"))
  cat(red(paste0("  • ", missing_wrappers, collapse = "\n")), "\n\n")
  stop("Master execution halted. Verify script paths in 'scripts/'.")
}

cat(green("✅ All 4 stage wrapper scripts located and validated.\n\n"))

# ==============================================================================
# 3. STAGE EXECUTION ENGINE
# ==============================================================================

execute_stage <- function(stage_info) {
  cat(bold(magenta(sprintf("\n▶ EXECUTING %s: %s\n", toupper(stage_info$id), stage_info$name))))
  cat(silver(sprintf("  Script: %s\n", stage_info$path)))
  
  stage_t0 <- Sys.time()
  stage_env <- new.env(parent = globalenv())
  
  tryCatch(
    expr = {
      source(stage_info$path, local = stage_env, echo = FALSE)
    },
    error = function(e) {
      cat(red(sprintf("\n❌ CRITICAL CRASH IN %s [%s]:\n   %s\n", 
                      stage_info$id, basename(stage_info$path), e$message)))
      stop(e)
    }
  )
  
  duration_sec <- round(as.numeric(difftime(Sys.time(), stage_t0, units = "secs")), 1)
  
  # Audit outputs
  missing_deliverables <- stage_info$outputs[!file.exists(stage_info$outputs) | file.info(stage_info$outputs)$size == 0]
  if (length(missing_deliverables) > 0) {
    cat(red(sprintf("\n⚠️ Deliverable verification failed for %s. Missing or empty files:\n", stage_info$id)))
    cat(red(paste0("  • ", basename(missing_deliverables), collapse = "\n")), "\n")
    stop(sprintf("%s failed to generate required outputs.", stage_info$id))
  }
  
  cat(green(sprintf("✅ %s finished successfully in %s s. (%d deliverables verified)\n", 
                    stage_info$id, duration_sec, length(stage_info$outputs))))
  
  return(tibble(
    Stage        = stage_info$id,
    Name         = stage_info$name,
    Script       = basename(stage_info$path),
    Status       = "Success",
    Deliverables = sprintf("%d / %d verified", length(stage_info$outputs), length(stage_info$outputs)),
    Duration     = sprintf("%.1f s", duration_sec)
  ))
}

# ==============================================================================
# 4. MASTER SEQUENTIAL RUNNER
# ==============================================================================

master_audit_log <- list()

master_audit_log[[1]] <- execute_stage(pipeline_stages$stage1)
master_audit_log[[2]] <- execute_stage(pipeline_stages$stage2)
master_audit_log[[3]] <- execute_stage(pipeline_stages$stage3)
master_audit_log[[4]] <- execute_stage(pipeline_stages$stage4)

# ==============================================================================
# 5. PIPELINE AUDIT REPORT & SUMMARY
# ==============================================================================

total_min <- round(as.numeric(difftime(Sys.time(), master_start_time, units = "mins")), 2)
master_summary_df <- bind_rows(master_audit_log)

cat(bold(cyan("\n==============================================================================\n")))
cat(bold(cyan(" 📋 DARWIN POINT COMPLETE PIPELINE AUDIT REPORT\n")))
cat(bold(cyan("==============================================================================\n")))
print(as.data.frame(master_summary_df), row.names = FALSE)

cat(bold(green(sprintf("\n🎉 ALL 4 PIPELINE STAGES COMPLETED IN %.2f MINUTES!\n", total_min))))
cat(blue(sprintf("📁 Calibrated Databases  : %s\n", dbase_dir)))
cat(blue(sprintf("📁 Fixed Budgets         : %s\n", fixed_dir)))
cat(blue(sprintf("📁 StRS Regional Budgets : %s\n", strs_dir)))
cat(blue(sprintf("📁 Grigg P90 Modeling    : %s\n", grigg_dir)))
cat(blue(sprintf("📊 Publication Figures   : %s\n\n", plot_dir)))