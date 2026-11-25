

#function to simulate replicates for each scenario from base case to variation
# Run one simulation replicate for a given scenario variation, fit all methods
one_run_one_variation2 <- function(i, variation_id, RI, years, master_seed = 12) {
  m      <- .Machine$integer.max
  sseed  <- scenario_seed(master_seed, RI, years) #seed for base simulation (varies by scenario only i.e., RI and years)
  
  # Base seed: varies by replicate i 
  seed_base <- as.integer(abs(((sseed*1103515245 + i*12345) %% m))) #seed for base case and varies by simulation replicate only and constant ffor scenario (years and RI)
  base_df   <- simulate_base(RI = RI, years = years, seed = seed_base) #base case simulated using seed base variation=0
  
  # Variation seed: independent per (scenario, i, variation_id)
  seed_var  <- variation_seed(sseed, i, variation_id) #seed for variation variation=0:6 and where 0 is base case, varies by simulation replicate and variation_id and constant for scenario (years and RI)

  #simulate data for a particular replicate
  dat <- if (variation_id == 0) base_df else #if variation id=0, it returns the base case for that replicate, otheriwise it makes a particular variation from the base case for that replicate i
    make_variation(
      base_df,
      variation_id = variation_id,
      RI = RI,
      years = years,
      seed = seed_var #when making variation it uses the seed var
    )

  #the simulated data for replicate i is then fitted for all the methods, estimates and standard errors from the model are extracted and results are stored along with their seeds that used to simulate the data to ensure reproducibility
  purrr::imap_dfr(fitters, function(fitter, method_name) {
    fit <- tryCatch(fitter(dat), error = function(e) NULL)
    extract_ri_se(fit, approach = method_name) %>%
      dplyr::mutate(sim = i, variation = variation_id, RI_true = RI, years = years, master_seed = master_seed,
    seed_base=seed_base, seed_var=seed_var)
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


