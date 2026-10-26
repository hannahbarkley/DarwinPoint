# ==============================================================================
# Script:  2. DarwinPoint_process_fixed.R
# Purpose: Master execution wrapper for the Hawaiian Archipelago Fixed-Site 
#          carbonate budget pipelines. Orchestrates runtime profiling, 
#          dependency validation, error handling, and output verification for:
#            1. 2a. DarwinPoint_processCB_fixed_sfm.R  (SfM Photogrammetry Budget)
#            2. 2b. DarwinPoint_processCB_fixed_cnet.R (StRS-matched Fixed Budget)
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
scripts_dir <- here("scripts", "2. fixed scripts")
data_dir    <- file.path(base_dir, "data")
proc_dir    <- file.path(data_dir, "fixed data", "processed data")
plot_dir    <- file.path(base_dir, "outputs", "supporting figures")

# Define target script paths and expected primary deliverables
scripts <- list(
  step_2a = list(
    name        = "2a. DarwinPoint_processCB_fixed_sfm.R",
    path        = file.path(scripts_dir, "2a. DarwinPoint_processCB_fixed_sfm.R"),
    description = "Structure-from-Motion (SfM) Benthic Production & Net Carbonate Budget",
    outputs     = c(
      file.path(proc_dir, sprintf("%s_%s_Fixed_Urchins_SiteLevel.csv", region, year)),
      file.path(proc_dir, sprintf("%s_%s_Fixed_Fish_SiteLevel.csv", region, year)),
      file.path(proc_dir, sprintf("%s_%s_Fixed_Net_SiteLevel_sfm.csv", region, year)),
      file.path(plot_dir, sprintf("%s_%s_gross_overview_sfm.png", region, year)),
      file.path(plot_dir, sprintf("%s_%s_net_overview_sfm.png", region, year)),
      file.path(plot_dir, sprintf("%s_%s_netprod_fixed_sfm.png", region, year))
    )
  ),
  step_2b = list(
    name        = "2b. DarwinPoint_processCB_fixed_cnet.R",
    path        = file.path(scripts_dir, "2b. DarwinPoint_processCB_fixed_cnet.R"),
    description = "StRS Gross Production Match, Coral Cover & Fixed Site Net Budget",
    outputs     = c(
      file.path(proc_dir, sprintf("%s_%s_Fixed_Net_SiteLevel_cnet.csv", region, year)),
      file.path(plot_dir, sprintf("%s_%s_gross_overview.png", region, year)),
      file.path(plot_dir, sprintf("%s_%s_net_overview.png", region, year)),
      file.path(plot_dir, sprintf("%s_%s_netprod_fixed.png", region, year)),
      file.path(plot_dir, sprintf("%s_%s_Fixed_netprod_cnet.png", region, year)),
      file.path(plot_dir, sprintf("%s_%s_Fixed_netprod_vs_coralcover_cnet.png", region, year))
    )
  )
)

# ==============================================================================
# 2. PRE-FLIGHT VERIFICATION
# ==============================================================================

cat(bold(cyan("\n==============================================================================\n")))
cat(bold(cyan(sprintf(" 🪸 MASTER RUNNER: HAWAIIAN ARCHIPELAGO FIXED-SITE BUDGETS (%s)\n", year))))
cat(bold(cyan("==============================================================================\n\n")))

cat(blue("🔍 Running pre-flight script and directory checks...\n"))

# Ensure destination output directories exist
invisible(lapply(c(proc_dir, plot_dir), function(d) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}))

# Check that pipeline scripts exist before execution
missing_scripts <- character(0)
for (s in scripts) {
  if (!file.exists(s$path)) {
    missing_scripts <- c(missing_scripts, sprintf("%s (Expected at: %s)", s$name, s$path))
  }
}

if (length(missing_scripts) > 0) {
  cat(red("\n❌ ERROR: Pipeline scripts were not found:\n"))
  cat(red(paste0("  • ", missing_scripts, collapse = "\n")), "\n\n")
  stop("Execution halted due to missing source scripts.")
}

cat(green("✅ All fixed-site pipeline scripts verified.\n\n"))

# ==============================================================================
# 3. PIPELINE STEP EXECUTION ENGINE
# ==============================================================================

execute_step <- function(step_obj) {
  cat(bold(magenta(sprintf("▶ Executing: %s\n", step_obj$name))))
  cat(silver(sprintf("  Purpose: %s\n", step_obj$description)))
  
  step_start <- Sys.time()
  
  # Execute target script in the global workspace
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
    cat(red(sprintf("\n⚠️ Output validation failed for %s. Missing deliverables:\n", step_obj$name)))
    cat(red(paste0("  • ", basename(missing_outs), collapse = "\n")), "\n")
    stop(sprintf("Step %s failed to generate all required outputs.", step_obj$name))
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

# Step 1: Run SfM Photogrammetry Net Carbonate Budget & Visualizations (2a)[cite: 4]
run_summary[[1]] <- execute_step(scripts$step_2a)

# Step 2: Run StRS-matched Fixed Net Carbonate Budget & Figures (2b)[cite: 5]
run_summary[[2]] <- execute_step(scripts$step_2b)

# ==============================================================================
# 5. PIPELINE SUMMARY & AUDIT LOG
# ==============================================================================

total_duration <- round(as.numeric(difftime(Sys.time(), pipeline_start_time, units = "mins")), 2)
summary_df     <- bind_rows(run_summary)

cat(bold(cyan("==============================================================================\n")))
cat(bold(cyan(" 📋 FIXED-SITE PIPELINE EXECUTION SUMMARY REPORT\n")))
cat(bold(cyan("==============================================================================\n")))
print(as.data.frame(summary_df), row.names = FALSE)

cat(bold(green(sprintf("\n🎉 FIXED-SITE PIPELINES COMPLETED SUCCESSFULLY in %.2f minutes!\n", total_duration))))
cat(blue(sprintf("📁 Processed Data: %s\n", proc_dir)))
cat(blue(sprintf("📊 Supporting Figures: %s\n\n", plot_dir)))