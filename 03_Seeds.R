scenario_seed <- function(master_seed, RI, years) {
  m <- .Machine$integer.max
  as.integer(abs((master_seed*39461758 + RI*16492649 + years*19562750) %% m))
}

variation_seed <- function(sseed, i, variation_id) {
  m <- .Machine$integer.max
  as.integer(abs((sseed*2654435761 + i*1013904223 + variation_id*1664525) %% m))
}

