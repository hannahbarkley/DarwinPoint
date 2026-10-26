# ==============================================================================
# Script:  3. DarwinPoint_process_strs.R
# Purpose: Master execution wrapper for the Hawaiian Archipelago StRS carbonate
#          budget pipeline. Orchestrates end-to-end execution, runtime profiling,
#          dependency checks, and output validation for:
#            1. 3a. DarwinPoint_processCB_strs_cnet.R (Net Budget & Figures)
#            2. 3b. DarwinPoint_calc_porosity.R       (Framework Porosity)
#            3. 3c. DarwinPoint_strs2island.R         (Stratum Weighting & Scaling)
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(crayon)
  library(tibble)
  library(dplyr)
})

# ==============================================================================
# 1. PIPELINE CONFIGURATION & PATH SETUP
# ==============================================================================

pipeline_start_time <- Sys.time()

region <- "HA"
year   <- "2024"

base_dir    <- here()
scripts_dir <- here("scripts","3. strs scripts")
data_dir    <- file.path(base_dir, "data")
proc_dir    <- file.path(data_dir, "strs data", "processed data")
plot_dir    <- file.path(base_dir, "outputs", "supporting figures")

# Define target script paths
scripts <- list(
  step_3a = list(
    name        = "3a. DarwinPoint_processCB_strs_cnet.R",
    path        = file.path(scripts_dir, "3a. DarwinPoint_processCB_strs_cnet.R"),
    description = "StRS Net Carbonate Budget Calculation & Plot Generation",
    outputs     = c(
      file.path(proc_dir, sprintf("%s_%s_StRS_Fish_SiteLevel.csv", region, year)),
      file.path(proc_dir, sprintf("%s_%s_StRS_Net_SiteLevel_cnet.csv", region, year)),
      file.path(plot_dir, sprintf("%s_%s_netprod_strs.png", region, year))
    )
  ),
  step_3b = list(
    name        = "3b. DarwinPoint_calc_porosity.R",
    path        = file.path(scripts_dir, "3b. DarwinPoint_calc_porosity.R"),
    description = "Structural Framework Porosity Calculation",
    outputs     = c(
      file.path(proc_dir, sprintf("%s_%s_StRS_Porosity.csv", region, year))
    )
  ),
  step_3c = list(
    name        = "3c. DarwinPoint_strs2island.R",
    path        = file.path(scripts_dir, "3c. DarwinPoint_strs2island.R"),
    description = "Strata Sampling Weighting, Survey Design & Island Scaling",
    outputs     = c(
      file.path(proc_dir, sprintf("%s_%s_StRS_SiteLevel_Weighted.csv", region, year)),
      file.path(proc_dir, sprintf("%s_%s_StRS_SurveyDesign.rds", region, year)),
      file.path(proc_dir, sprintf("%s_%s_StRS_IslandLevel.csv", region, year))
    )
  )
)

# ==============================================================================
# 2. PRE-FLIGHT VERIFICATION
# ==============================================================================

cat(bold(cyan("\n==============================================================================\n")))
cat(bold(cyan(sprintf(" 🌊 MASTER RUNNER: HAWAIIAN ARCHIPELAGO StRS CARBONATE BUDGET PIPELINE (%s)\n", year))))
cat(bold(cyan("==============================================================================\n\n")))

cat(blue("🔍 Running pre-flight script and directory checks...\n"))

# Ensure scripts exist before execution begins
missing_scripts <- character(0)
for (s in scripts) {
  if (!file.exists(s$path)) {
    missing_scripts <- c(missing_scripts, sprintf("%s (Expected at: %s)", s$name, s$path))
  }
}

if (length(missing_scripts) > 0) {
  cat(red("\n❌ ERROR: One or more pipeline scripts were not found:\n"))
  cat(red(paste0("  • ", missing_scripts, collapse = "\n")), "\n\n")
  stop("Execution aborted due to missing source scripts.")
}

cat(green("✅ All pipeline scripts located successfully.\n\n"))

# ==============================================================================
# 3. PIPELINE STEP EXECUTION ENGINE
# ==============================================================================

execute_step <- function(step_obj) {
  cat(bold(magenta(sprintf("▶ Executing: %s\n", step_obj$name))))
  cat(silver(sprintf("  Purpose: %s\n", step_obj$description)))
  
  step_start <- Sys.time()
  
  # Run the target script in a fresh child environment or global environment
  tryCatch(
    expr = {
      source(step_obj$path, local = FALSE, echo = FALSE)
    },
    error = function(e) {
      cat(red(sprintf("\n❌ FATAL CRASH in [%s]: %s\n", step_obj$name, e$message)))
      stop(e)
    }
  )
  
  step_duration <- round(as.numeric(difftime(Sys.time(), step_start, units = "secs")), 2)
  
  # Validate that expected deliverables were created
  missing_outs <- step_obj$outputs[!file.exists(step_obj$outputs)]
  if (length(missing_outs) > 0) {
    cat(red(sprintf("\n⚠️ Output validation failed for %s. Missing files:\n", step_obj$name)))
    cat(red(paste0("  • ", basename(missing_outs), collapse = "\n")), "\n")
    stop(sprintf("Step %s failed to produce expected deliverables.", step_obj$name))
  }
  
  cat(green(sprintf("✅ Completed [%s] in %s seconds.\n\n", step_obj$name, step_duration)))
  
  return(tibble(
    Script   = step_obj$name,
    Status   = "Success",
    Duration = sprintf("%.2f s", step_duration)
  ))
}

# ==============================================================================
# 4. MASTER PIPELINE ORCHESTRATION
# ==============================================================================

run_summary <- list()

# Step 1: Site-Level Net Production Calculation & Overview Figures (3a)
run_summary[[1]] <- execute_step(scripts$step_3a)

# Step 2: Site-Level Morphological Framework Porosity Calculation (3b)
run_summary[[2]] <- execute_step(scripts$step_3b)

# Step 3: Spatial Weighting, Survey Design & Island-Level Aggregation (3c)
run_summary[[3]] <- execute_step(scripts$step_3c)

# ==============================================================================
# 5. PIPELINE SUMMARY & AUDIT LOG
# ==============================================================================

total_duration <- round(as.numeric(difftime(Sys.time(), pipeline_start_time, units = "mins")), 2)
summary_df     <- bind_rows(run_summary)

cat(bold(cyan("==============================================================================\n")))
cat(bold(cyan(" 📋 PIPELINE EXECUTION SUMMARY REPORT\n")))
cat(bold(cyan("==============================================================================\n")))
print(as.data.frame(summary_df), row.names = FALSE)

cat(bold(green(sprintf("\n🎉 ALL PIPELINE MODULES COMPLETED SUCCESSFULLY in %.2f minutes!\n", total_duration))))
cat(blue(sprintf("📁 Primary Processed Data: %s\n", proc_dir)))
cat(blue(sprintf("📊 Primary Figures Output:  %s\n\n", plot_dir)))