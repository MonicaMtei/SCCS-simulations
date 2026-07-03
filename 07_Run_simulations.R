# Run one simulation replicate for a given scenario variation, fit all methods
one_run_one_variation2 <- function(i, variation_id, RI, years, master_seed = 12) {
  m      <- .Machine$integer.max
  sseed  <- scenario_seed(master_seed, RI, years)
  
  # Base seed: varies by replicate i 
  seed_base <- as.integer(abs(((sseed*1103515245 + i*12345) %% m)))
  base_df   <- simulate_base(RI = RI, years = years, seed = seed_base)
  
  # Variation seed: independent per (scenario, i, variation_id)
  seed_var  <- variation_seed(sseed, i, variation_id)
  
  dat <- if (variation_id == 0) base_df else
    make_variation(
      base_df,
      variation_id = variation_id,
      RI = RI,
      years = years,
      seed = seed_var
    )
  
  purrr::imap_dfr(fitters, function(fitter, method_name) {
    fit <- tryCatch(fitter(dat), error = function(e) NULL)
    extract_ri_se(fit, approach = method_name) %>%
      dplyr::mutate(sim = i, variation = variation_id, RI_true = RI, years = years)
  })
}


# variations to run
variations <- 0:6
nsim       <- 5000
# Run all RI/year scenarios across all variations and simulation replicates
run_all_scenarios <- function(RI_values = c(1, 2),
                              year_values = c(1, 4),
                              variations = 0:6,
                              nsim = 5000) {
  scen <- tidyr::expand_grid(RI = RI_values, years = year_values)
  
  purrr::pmap_dfr(scen, function(RI, years) {
    message(sprintf("Running RI=%s, years=%s ...", RI, years))
    
    purrr::map_dfr(variations, function(v) {
      # loop sims for one (variation, RI, years)
      purrr::map_dfr(seq_len(nsim), function(i) {
        one_run_one_variation2(i = i, variation_id = v, RI = RI, years = years)
      })
    }) %>%
      mutate(scenario = paste0("RI=", RI, ", Years=", years))
  })
}



# ---- run everything ----
results_all <- run_all_scenarios(
  RI_values  = c(1, 2),
  year_values = c(1, 4),
  variations = variations,
  nsim = nsim
)

# a results folder 
dir.create("results", showWarnings = FALSE)

# timestamped filename
fn <- sprintf("results_all_%s.rds",
              format(Sys.time(), "%Y%m%d_%H%M%S"))

path <- file.path("results", fn)
saveRDS(results_all, file = path)
cat("Saved to:", path, "\n")
print(file.info(path)$mtime)

saveRDS(results_all, file = "results/results_all.rds")


