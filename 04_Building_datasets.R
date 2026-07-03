
## Pipeline: build 4 base cases + 6 variations from each base
simulate_with_base_and_variations <- function(
    RI_values    = c(1, 2),
    year_values  = c(1, 4),
    n_variations = 6,
    base_seed    = 12,
    include_base = TRUE
) {
  # grid of scenarios (the 4 bases)
  scenarios <- tidyr::expand_grid(RI = RI_values, years = year_values) %>%
    mutate(
      base_seed = scenario_seed(base_seed, RI, years),
      base_data = pmap(list(RI, years, base_seed), simulate_base)
    ) %>%
    mutate(
      data_list = pmap(list(base_data, RI, years, base_seed), function(bdf, RI, years, sseed) {
        var_list <- purrr::map(1:n_variations, function(vid) {
          vseed <- variation_seed(sseed, i = 1L, variation_id = vid)
          make_variation(bdf, variation_id = vid, RI = RI, years = years, seed = vseed)
        })
        if (include_base) c(list(bdf), var_list) else var_list
      })
    )
  
  
  # unnest to one row per (RI, years, dataset), assign variation 0..n_variations (or 1..n)
  results <- scenarios %>%
    select(RI, years, data_list) %>%
    tidyr::unnest_longer(data_list, values_to = "data") %>%
    group_by(RI, years) %>%
    mutate(variation = row_number() - if (include_base) 1L else 0L) %>%
    ungroup()
  
  # columns: RI, years, variation (0..n), data (data.frame)
  results
}


results <- simulate_with_base_and_variations(
  RI_values    = c(1, 2),
  year_values  = c(1, 4),
  n_variations = 6,
  base_seed    = 12,
  include_base = TRUE
)
