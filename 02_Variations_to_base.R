make_variation <- function(base_df, variation_id, RI, years, seed) {
  set.seed(seed)
  out <- base_df
  
  ## ---- selection helpers
  pick_fixed_fraction <- function(idx, fraction) {
    n_sel <- round(fraction * length(idx))
    if (n_sel > 0) sample(idx, size = n_sel, replace = FALSE) else integer(0)
  }
  pick_fixed_count <- function(idx, n_sel) {
    if (n_sel > 0) sample(idx, size = n_sel, replace = FALSE) else integer(0)
  }
  
  
  ## ---- eligibility helpers
  
  elig_14_before_d1 <- function(d) {
    gap <- d$adrug1 - d$aevent
    which(!is.na(d$adrug1) & !is.na(d$aevent) & gap >= 1L & gap <= 14L)
  }
  elig_14_before_d2 <- function(d) {
    gap <- d$adrug2 - d$aevent
    which(!is.na(d$adrug1) & !is.na(d$adrug2) & !is.na(d$aevent) &
            d$aevent > d$adrug1 & gap >= 1L & gap <= 14L)
  }
  
  
  #  1..89 days BEFORE Dose 1
  elig_89_before_d1 <- function(d) {
    gap <- d$adrug1 - d$aevent
    which(!is.na(d$adrug1) & !is.na(d$aevent) & gap >= 1L & gap <= 89L)
  }
  
  #  1..89 days BEFORE Dose 2 AND after Dose 1
  elig_89_before_d2 <- function(d) {
    gap <- d$adrug2 - d$aevent
    which(!is.na(d$adrug1) & !is.na(d$adrug2) & !is.na(d$aevent) &
            (d$aevent > d$adrug1) & gap >= 1L & gap <= 89L)
  }
  
  elig_any_before_d1 <- function(d) which(!is.na(d$aevent) & !is.na(d$adrug1) & (d$aevent < d$adrug1))
  elig_any_between_d1d2 <- function(d) which(!is.na(d$aevent) & !is.na(d$adrug1) & !is.na(d$adrug2) &
                                               d$aevent > d$adrug1 & d$aevent < d$adrug2)
  
  
  
  
  ## ---- shift helpers
  
  shift_d1_plus_d2_days <- function(d, idx, delay_days) {
    if (!length(idx)) return(d)
    d$adrug1[idx]  <- d$adrug1[idx] + delay_days
    d$aedrug1[idx] <- d$adrug1[idx] + 28L
    has_d2 <- !is.na(d$adrug2[idx])
    rows_d2 <- idx[has_d2]
    if (length(rows_d2)) {
      d$adrug2[rows_d2]  <- d$adrug2[rows_d2] + delay_days
      d$aedrug2[rows_d2] <- d$adrug2[rows_d2] + 28L
    }
    d
  }
  shift_d2_only_days <- function(d, idx, delay_days) {
    if (!length(idx)) return(d)
    d$adrug2[idx]  <- d$adrug2[idx] + delay_days
    d$aedrug2[idx] <- d$adrug2[idx] + 28L
    d
  }
  
  
  
  # set D1 to event+L and propagate the same delta to D2 (preserve spacing)
  shift_d1_to_event_plusL <- function(d, idx, L = 90L, rp_len = 28L) {
    if (!length(idx)) return(d)
    target <- d$aevent[idx] + L
    delta  <- target - d$adrug1[idx]
    move   <- which(delta > 0L)   # exact 90-day leads to no change
    if (!length(move)) return(d)
    rows <- idx[move]; dlt <- delta[move]
    
    d$adrug1[rows]  <- d$adrug1[rows] + dlt
    d$aedrug1[rows] <- d$adrug1[rows] + rp_len
    
    has_d2 <- !is.na(d$adrug2[rows])
    if (any(has_d2)) {
      r2  <- rows[has_d2]; d2d <- dlt[has_d2]
      d$adrug2[r2]  <- d$adrug2[r2] + d2d
      d$aedrug2[r2] <- d$adrug2[r2] + rp_len
    }
    d
  }
  
  # set D2 to event+L directly (only forward; exact 90-day lead -> no change)
  set_d2_to_event_plusL <- function(d, idx, L = 90L, rp_len = 28L) {
    if (!length(idx)) return(d)
    target <- d$aevent[idx] + L
    move   <- which(target > d$adrug2[idx])  # exact 90-day lead -> no change
    if (!length(move)) return(d)
    rows <- idx[move]; tgt <- target[move]
    d$adrug2[rows]  <- tgt
    d$aedrug2[rows] <- d$adrug2[rows] + rp_len
    d
  }
  
  
  
  ## ---- variations -------
  if (variation_id == 1) {
    # event within 14d before each dose, delay vaccine by 14 days, p=0.5
    sim <- base_df
    idxA <- elig_14_before_d1(sim); selA <- pick_fixed_fraction(idxA, 0.5)
    sim  <- shift_d1_plus_d2_days(sim, selA, delay_days = 14L)
    idxB <- elig_14_before_d2(sim); selB <- pick_fixed_fraction(idxB, 0.5)
    sim  <- shift_d2_only_days(sim, selB, delay_days = 14L)
    out <- sim
    
  } else if (variation_id == 2) {
    # event within 14d before each dose, delay vaccine by 14 days, p=1
    sim <- base_df
    idxA <- elig_14_before_d1(sim); sim <- shift_d1_plus_d2_days(sim, idxA, 14L)
    idxB <- elig_14_before_d2(sim); sim <- shift_d2_only_days(sim, idxB, 14L)
    out <- sim
    
    
  } else if (variation_id == 3) {
    # Vaccination happens 90 days post-event, p = 0.5 
    sim <- base_df
    idxA0 <- elig_89_before_d1(base_df)
    idxB0 <- elig_89_before_d2(base_df)
    stopifnot(length(intersect(idxA0, idxB0)) == 0L)
    
    selA <- pick_fixed_fraction(idxA0, 0.5)
    selB <- pick_fixed_fraction(idxB0, 0.5)
    
    sim <- shift_d1_to_event_plusL(sim, selA, L = 90L, rp_len = 28L)  # D1 + propagate to D2
    sim <- set_d2_to_event_plusL(sim, selB, L = 90L, rp_len = 28L)    # D2 direct to event+90
    
    out <- sim
    
  } else if (variation_id == 4) {
    # Vaccination happens 90 days post-event, p = 1 
    sim <- base_df
    idxA0 <- elig_89_before_d1(base_df)
    idxB0 <- elig_89_before_d2(base_df)
    stopifnot(length(intersect(idxA0, idxB0)) == 0L)
    
    sim <- shift_d1_to_event_plusL(sim, idxA0, L = 90L, rp_len = 28L)
    sim <- set_d2_to_event_plusL(sim, idxB0, L = 90L, rp_len = 28L)
    
    out <- sim
    
    
    
    
  } else if (variation_id == 5) {
    # event anytime before each dose, PERMANENT CANCEL, p=0.5
    sim <- base_df
    idxA <- elig_any_before_d1(sim); idxB <- elig_any_between_d1d2(sim)
    pick <- pick_fixed_fraction(c(idxA, idxB), 0.5)
    idx_cancel_d1 <- intersect(pick, idxA)
    idx_cancel_d2 <- intersect(pick, idxB)
    if (length(idx_cancel_d1)) {
      sim$adrug1[idx_cancel_d1]  <- NA_integer_
      sim$aedrug1[idx_cancel_d1] <- NA_integer_
      sim$adrug2[idx_cancel_d1]  <- NA_integer_
      sim$aedrug2[idx_cancel_d1] <- NA_integer_
    }
    if (length(idx_cancel_d2)) {
      sim$adrug2[idx_cancel_d2]  <- NA_integer_
      sim$aedrug2[idx_cancel_d2] <- NA_integer_
    }
    out <- sim
    
  } else if (variation_id == 6) {
    # event anytime before each dose, CONTRAINDICATION, p=1
    sim <- base_df
    idxA <- elig_any_before_d1(sim); idxB <- elig_any_between_d1d2(sim)
    if (length(idxA)) {
      sim$adrug1[idxA]  <- NA_integer_
      sim$aedrug1[idxA] <- NA_integer_
      sim$adrug2[idxA]  <- NA_integer_
      sim$aedrug2[idxA] <- NA_integer_
    }
    if (length(idxB)) {
      sim$adrug2[idxB]  <- NA_integer_
      sim$aedrug2[idxB] <- NA_integer_
    }
    out <- sim
    
  } else{
    stop("unknown variation_id", variation_id)
  }
  
  out$variation <- variation_id
  out
}
