# ==============================================================================
# Script:  DarwinPoint_run_all_figures.R
# Purpose: Master execution wrapper for all publication-ready main figures (1–4),
#          Table 1, and supplementary figures (S1–S7) for the Hawaiian Darwin 
#          Point carbonate budget study.
#          Orchestrates dependency checks, pre-flight file validation, runtime 
#          profiling, error isolation, and deliverable verification.
# Outputs: - outputs/main figures/ (Figures 1–4, Table 1)
#          - outputs/supplementary figures/ (Figures S1–S6)
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(crayon)
  library(tibble)
  library(dplyr)
})

# ==============================================================================
# 1. PIPELINE CONFIGURATION & DIRECTORY RESOLUTION
# ==============================================================================

pipeline_start_time <- Sys.time()

base_dir <- here()
main_dir <- file.path(base_dir, "outputs", "main figures")
supp_dir <- file.path(base_dir, "outputs", "supplementary figures")

# Ensure output directories exist
invisible(lapply(c(main_dir, supp_dir), function(d) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}))

# Candidate directories where individual figure scripts might be located
script_search_dirs <- c(
  file.path(base_dir, "scripts", "figures"),
  file.path(base_dir, "scripts", "5. figures"),
  file.path(base_dir, "scripts"),
  base_dir
)

find_script_path <- function(script_name) {
  for (d in script_search_dirs) {
    p <- file.path(d, script_name)
    if (file.exists(p)) return(p)
  }
  return(file.path(script_search_dirs[1], script_name))
}

# Define each figure task, its source path, description, and expected deliverables
figure_tasks <- list(
  # --- MAIN FIGURES & TABLE ---
  list(
    id          = "Figure 1",
    name        = "Figure1_HA_map.R",
    path        = find_script_path("Figure1_HA_map.R"),
    section     = "Main Figures",
    description = "Full Archipelago Bathymetric Map and Global Inset",
    outputs     = c(
      file.path(main_dir, "Figure1_HA_map.png"),
      file.path(main_dir, "Figure1_HA_map.tiff")
    )
  ),
  list(
    id          = "Figure 2",
    name        = "Figure2_Grigg_vs_NCRMP.R",
    path        = find_script_path("Figure2_Grigg_vs_NCRMP.R"),
    section     = "Main Figures",
    description = "Historical (1982 Grigg) vs. Modern (2024 NCRMP P90) Comparison",
    outputs     = c(
      file.path(main_dir, "Figure2_Grigg_vs_NCRMP.png"),
      file.path(main_dir, "Figure2_Grigg_vs_NCRMP.tiff")
    )
  ),
  list(
    id          = "Figure 3",
    name        = "Figure3_FullStRS.R",
    path        = find_script_path("Figure3_FullStRS.R"),
    section     = "Main Figures",
    description = "Archipelago-Wide StRS Survey Metric Scaling across Latitude",
    outputs     = c(
      file.path(main_dir, "Figure3_FullStRS.png"),
      file.path(main_dir, "Figure3_FullStRS.tiff")
    )
  ),
  list(
    id          = "Figure 4",
    name        = "Figure4_VARvsSLR.R",
    path        = find_script_path("Figure4_VARvsSLR.R"),
    section     = "Main Figures",
    description = "Vertical Accretion Rate (VAR) vs. 2050 Predicted Sea-Level Rise",
    outputs     = c(
      file.path(main_dir, "Figure4_VARvsSLR.png"),
      file.path(main_dir, "Figure4_VARvsSLR.tiff")
    ),
    alt_outputs = c(
      file.path(main_dir, "Figure4_VARvsSLR.pdf")
    )
  ),
  list(
    id          = "Table 1",
    name        = "Table1_IslandSummary.R",
    path        = find_script_path("Table1_IslandSummary.R"),
    section     = "Main Figures",
    description = "Journal Summary Table of Island Carbonate Budgets",
    outputs     = c(
      file.path(main_dir, "Table1_IslandSummary.docx")
    ),
    alt_outputs = c(
      file.path(main_dir, "Table1_IslandSummary.rtf")
    )
  ),
  
  # --- SUPPLEMENTARY FIGURES (S1–S7) ---
  list(
    id          = "Figure S1",
    name        = "FigureS1_Sites.R",
    path        = find_script_path("FigureS1_Sites.R"),
    section     = "Supplementary",
    description = "Spatial Sampling Distribution across MHI and NWHI Insets",
    outputs     = c(
      file.path(supp_dir, "FigureS1_Sites.png"),
      file.path(supp_dir, "FigureS1_Sites.tiff")
    )
  ),
  list(
    id          = "Figure S2",
    name        = "FigureS2_growthrates.R",
    path        = find_script_path("FigureS2_growthrates.R"),
    section     = "Supplementary",
    description = "Coral Linear Extension, Skeletal Density, and Colony Production",
    outputs     = c(
      file.path(supp_dir, "FigureS2_growthrates.png"),
      file.path(supp_dir, "FigureS2_growthrates.tiff")
    )
  ),
  list(
    id          = "Figure S3",
    name        = "FigureS3_Prod_Scaling_Comparison.R",
    path        = find_script_path("FigureS3_Prod_Scaling_Comparison.R"),
    section     = "Supplementary",
    description = "Island-Level Scaled (Darwin Point) vs. Unscaled Gross Production",
    outputs     = c(
      file.path(supp_dir, "FigureS3_Prod_Scaling_Comparison.png"),
      file.path(supp_dir, "FigureS3_Prod_Scaling_Comparison.tiff")
    )
  ),
  list(
    id          = "Figure S4",
    name        = "FigureS4_SfM_vs_CoralNet.R",
    path        = find_script_path("FigureS4_SfM_vs_CoralNet.R"),
    section     = "Supplementary",
    description = "SfM Photogrammetry vs. CoralNet Gross Production Validation",
    outputs     = c(
      file.path(supp_dir, "FigureS4_SfM_vs_CoralNet.png"),
      file.path(supp_dir, "FigureS4_SfM_vs_CoralNet.tiff")
    )
  ),
  list(
    id          = "Figure S5",
    name        = "FigureS5_Sensitivity.R",
    path        = find_script_path("FigureS5_Sensitivity.R"),
    section     = "Supplementary",
    description = "Sensitivity Analysis: All Habitats vs. SW Fore-Reef Subsets",
    outputs     = c(
      file.path(supp_dir, "FigureS5_Sensitivity.png"),
      file.path(supp_dir, "FigureS5_Sensitivity.tiff")
    )
  ),
  list(
    id          = "Figure S6",
    name        = "FigureS6_Bioerosion.R",
    path        = find_script_path("FigureS6_Bioerosion.R"),
    section     = "Supplementary",
    description = "Urchin and Parrotfish Bioerosion Rates across Latitude",
    outputs     = c(
      file.path(supp_dir, "FigureS6_Bioerosion.png"),
      file.path(supp_dir, "FigureS6_Bioerosion.tiff")
    )
  )
)

# ==============================================================================
# 2. PRE-FLIGHT VERIFICATION
# ==============================================================================

cat(bold(cyan("\n==============================================================================\n")))
cat(bold(cyan(" 🎨 MASTER RUNNER: HAWAIIAN DARWIN POINT FIGURES & TABLES\n")))
cat(bold(cyan("==============================================================================\n\n")))

cat(blue("🔍 Running pre-flight script discovery...\n"))

missing_scripts <- character(0)
for (task in figure_tasks) {
  if (!file.exists(task$path)) {
    missing_scripts <- c(missing_scripts, sprintf("%s: %s (Searched in: %s)", task$id, task$name, task$path))
  }
}

if (length(missing_scripts) > 0) {
  cat(red("\n❌ ERROR: One or more figure generation scripts were not found:\n"))
  cat(red(paste0("  • ", missing_scripts, collapse = "\n")), "\n\n")
  stop("Execution aborted due to missing source scripts. Verify file placement in 'scripts/figures/' or project root.")
}

cat(green(sprintf("✅ All %d figure generation scripts located successfully.\n\n", length(figure_tasks))))

# ==============================================================================
# 3. FIGURE EXECUTION ENGINE
# ==============================================================================

execute_figure <- function(task) {
  cat(bold(magenta(sprintf("▶ Executing [%s]: %s\n", task$id, task$name))))
  cat(silver(sprintf("  Purpose: %s\n", task$description)))
  
  step_start <- Sys.time()
  step_env   <- new.env(parent = globalenv())
  
  # Run in an isolated child environment to prevent variable cross-contamination
  tryCatch(
    expr = {
      source(task$path, local = step_env, echo = FALSE)
    },
    error = function(e) {
      cat(red(sprintf("\n❌ FATAL CRASH in [%s]: %s\n", task$name, e$message)))
      stop(e)
    }
  )
  
  step_duration <- round(as.numeric(difftime(Sys.time(), step_start, units = "secs")), 1)
  
  # Check primary outputs (or alternative outputs like RTF for Table 1 or PDF for Figure 4)
  primary_missing <- task$outputs[!file.exists(task$outputs) | file.info(task$outputs)$size == 0]
  
  verified_count <- sum(file.exists(task$outputs) & file.info(task$outputs)$size > 0)
  total_expected <- length(task$outputs)
  
  if (length(primary_missing) > 0) {
    if (!is.null(task$alt_outputs) && any(file.exists(task$alt_outputs) & file.info(task$alt_outputs)$size > 0)) {
      verified_count <- max(verified_count, 1)
    } else {
      cat(red(sprintf("\n⚠️ Deliverable validation failed for %s. Missing files:\n", task$name)))
      cat(red(paste0("  • ", basename(primary_missing), collapse = "\n")), "\n")
      stop(sprintf("Task %s failed to produce expected output files.", task$name))
    }
  }
  
  cat(green(sprintf("✅ [%s] completed in %s s. (%d deliverables verified)\n\n", 
                    task$id, step_duration, verified_count)))
  
  return(tibble(
    ID           = task$id,
    Section      = task$section,
    Script       = task$name,
    Status       = "Success",
    Deliverables = sprintf("%d / %d verified", verified_count, total_expected),
    Duration     = sprintf("%.1f s", step_duration)
  ))
}

# ==============================================================================
# 4. SEQUENTIAL EXECUTION
# ==============================================================================

execution_log <- list()

for (i in seq_along(figure_tasks)) {
  execution_log[[i]] <- execute_figure(figure_tasks[[i]])
}

# ==============================================================================
# 5. AUDIT REPORT & SUMMARY
# ==============================================================================

total_min  <- round(as.numeric(difftime(Sys.time(), pipeline_start_time, units = "mins")), 2)
summary_df <- bind_rows(execution_log)

cat(bold(cyan("==============================================================================\n")))
cat(bold(cyan(" 📋 FIGURE & TABLE GENERATION SUMMARY REPORT\n")))
cat(bold(cyan("==============================================================================\n")))
print(as.data.frame(summary_df), row.names = FALSE)

cat(bold(green(sprintf("\n🎉 ALL FIGURES & TABLES GENERATED SUCCESSFULLY in %.2f minutes!\n", total_min))))
cat(blue(sprintf("📁 Main Figures & Table 1 : %s\n", main_dir)))
cat(blue(sprintf("📁 Supplementary Figures   : %s\n\n", supp_dir)))