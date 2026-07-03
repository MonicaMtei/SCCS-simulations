# check distribution is as intended
analyze_one_dist <- function(sim, id = "dataset") {
  max_day <- max(sim$aend, na.rm = TRUE)
  xbreaks <- if (max_day <= 365) seq(0, max_day, 30) else seq(0, max_day, 180)
  
  counts <- tibble(
    id = id,
    over365_d1  = sum(sim$adrug1 > 365,  na.rm = TRUE),
    over365_d2  = sum(sim$adrug2 > 365,  na.rm = TRUE),
    over1460_d1 = sum(sim$adrug1 > 1460, na.rm = TRUE),
    over1460_d2 = sum(sim$adrug2 > 1460, na.rm = TRUE)
  )
  
  d1s <- summary(sim$adrug1); d2s <- summary(sim$adrug2)
  mins_maxs <- tibble(
    id = id,
    d1_min = unname(d1s["Min."]), d1_max = unname(d1s["Max."]),
    d2_min = unname(d2s["Min."]), d2_max = unname(d2s["Max."])
  )
  
  p_dose1 <- ggplot(sim, aes(x = adrug1)) +
    geom_density(fill = "skyblue", alpha = 0.4) +
    scale_x_continuous(limits = c(0, max_day), breaks = xbreaks) +
    labs(title = paste0("Dose 1 density — ", id),
         x = "Day of observation", y = "Density") +
    theme_minimal()
  
  p_dose2 <- ggplot(sim, aes(x = adrug2)) +
    geom_density(fill = "skyblue", alpha = 0.4) +
    scale_x_continuous(limits = c(0, max_day), breaks = xbreaks) +
    labs(title = paste0("Dose 2 density — ", id),
         x = "Day of observation", y = "Density") +
    theme_minimal()
  
  list(counts = counts, mins_maxs = mins_maxs, plots = list(dose1 = p_dose1, dose2 = p_dose2))
}

# apply to all datasets 2 RIs, 2 obs periods × (1 base + 6 vars) = 28 datasets)
analyze_all_dists <- function(results_tbl) {
  results_tbl %>%
    mutate(
      id    = paste0("RI", RI, "years", years, "variation", variation),
      out   = map2(data, id, analyze_one_dist),
      counts    = map(out, "counts"),
      mins_maxs = map(out, "mins_maxs"),
      plots     = map(out, "plots")
    ) %>%
    select(RI, years, variation, id, counts, mins_maxs, plots)
}


all_dist <- analyze_all_dists(results)
all_dist
