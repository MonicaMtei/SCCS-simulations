

#seed function for base case simulation (varies by scenario only i.e., RI and years), the large multiplied numbers ensures that the seed is different for each scenario to avoid simulating datasets that are correlated across scenarios (RI and years)
scenario_seed <- function(master_seed, RI, years) {
  m <- .Machine$integer.max
  as.integer(abs((master_seed*39461758 + RI*16492649 + years*19562750) %% m))
}


#seed function for variation simulation (varies by scenario, replicate, and variation_id), the large multiplied numbers ensures that the seed is different for each scenario to avoid simulating datasets that are correlated across scenarios (RI and years) and also different for each replicate and variation_id to avoid simulating datasets that are correlated across replicates and variations
variation_seed <- function(sseed, i, variation_id) {
  m <- .Machine$integer.max
  as.integer(abs((sseed*2654435761 + i*1013904223 + variation_id*1664525) %% m))
}

