#function to simulate base case under left-skewed distribution of vaccine uptake over the observation period
simulate_base <- function(RI, years, seed = seed) {
  set.seed(seed)
  
  # === simulate ===
  
  # --- Inputs ---
  n      <- if (years == 1) 500 else 2000 #number of individuals in the simulation, 500 for 1 year and 2000 for 4 years
  astart <- 0 #start of observation period
  aend   <- years * 425 #end of observation period, 425 days for 1 year and 1460 days for 4 years
  startdose1 <- 60 #start of dose 1 distribution, 60 days after start of observation period
  enddose1   <- if (years == 1) 425 else 1460 #end of dose 1 distribution, 425 days for 1 year and 1460 days for 4 years
  risklen <- 28L #length of risk window, 28 days after vaccination
  
  
  # left-skewed distribution function to simulate dose 1 distribution where dose 1 days are concentrated earlier on in the observation period
  rleft_skew_int <- function(n, lo, hi, shape1 = 1.5, shape2 = 4) {
    u <- stats::rbeta(n, shape1, shape2)  
    x <- lo + floor(u * (hi - lo + 1L))
    pmin(x, hi)
  }
  
  dose1 <- rleft_skew_int(n, startdose1, enddose1, shape1 = 1.5, shape2 = 4) #similate dose 1 day following the left-skewed distribution
  gap12 <- sample(60:90, n, replace = TRUE)  #gap between dose 1 and dose 2, sampled uniformly from 60 to 90 days to ensure no overlap with the risk window of 28 days           
  dose2 <- dose1 + gap12 #dose 2 day is dose 1 day plus the gap between dose 1 and dose 2
  
  
  r1s <- dose1; r1e <- dose1 + risklen #start and end of risk window for dose 1
  r2s <- dose2; r2e <- dose2 + risklen #start and end of risk window for dose 2
  
  # Person-Time accounting
  S1 <- pmax(astart, r1s); E1 <- pmin(aend, r1e); len1 <- pmax(0L, E1 - S1 + 1L) #S1 is start of risk window for dose 1, E1 is end of risk window for dose 1, len1 is length of risk window for dose 1
  S2 <- pmax(astart, r2s); E2 <- pmin(aend, r2e); len2 <- pmax(0L, E2 - S2 + 1L) #S2 is start of risk window for dose 2, E2 is end of risk window for dose 2, len2 is length of risk window for dose 2
  Lrisk <- len1 + len2 #total length of risk window for both doses
  
  Lpre <- pmax(0L, dose1 - astart) #length of pre-risk window, from start of observation period to start of risk window for dose 1
  gap_start <- pmax(astart, E1 + 1L) #start of gap between risk window for dose 1 and start of risk window for dose 2, it starts at the end of risk period of dose 1+1day
  gap_end   <- pmin(aend,   S2 - 1L) #end of gap between risk window for dose 1 and start of risk window for dose 2, it ends at the start of risk period for dose 2-1day
  Lgap      <- pmax(0L, gap_end - gap_start + 1L) #length of gap between risk window for dose 1 and start of risk window for dose 2
  post_start <- pmax(astart, E2 + 1L) #start of post-risk window (from end of risk window for dose 2 to end of observation period)
  post_end   <- aend #end of post-risk window (end of observation period)
  Lpost_tail <- pmax(0L, post_end - post_start + 1L) #length of post-risk window, from end of risk window for dose 2 to end of observation period
  Lpost      <- Lgap + Lpost_tail #total length of post-risk windows 
  
  den <- Lpre + RI * Lrisk + Lpost #denominator for multinomial probabilities. It has three blocks (Lpre time before first risk window, Lrisk risk windows time and Lpost time after last risk window) 
  p1 <- ifelse(den > 0, Lpre / den, 0) #probability of event occurring in the pre-risk window, calculated as the length of pre-risk window divided by the total length of all three blocks
  p2 <- ifelse(den > 0, (RI * Lrisk) / den, 0) #probability of event occurring in the risk windows, calculated as the length of risk windows multiplied by the relative incidence (RI) divided by the total length of all three blocks
  #p3 is implicit as 1-(p1+p2) and is the probability of event occurring in the post-risk window
  
  u <- runif(n)  #draw a uniform random number for each individual to determine which block the event will occur in
  block_int <- ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L)) #draw block (1=pre, 2=risk, 3=post) based on the probabilities calculated above. If the uniform random number is less than p1, the event occurs in the pre-risk window (block 1). If it is between p1 and p1+p2, the event occurs in the risk windows (block 2). Otherwise, it occurs in the post-risk window (block 3).
  block <- factor(block_int, levels = 1:3, labels = c("pre","risk","post")) #create a factor variable for the block, with levels 1, 2, and 3 corresponding to "pre", "risk", and "post" respectively.
  
  # draw event day uniformly within chosen block (1.pre/window before dose 1, 2.risk period, and 3.post/windows between end of risk period for dose 1 and start of risk 2 period and post/windows after risk period for dose 2 to the end of observation period)
  event_day <- integer(n) 
  # pre
  idx1 <- which(block_int == 1L & Lpre > 0L) #find the indices of individuals whose event occurs in the pre-risk window (block 1) and have a non-zero length of pre-risk window (Lpre > 0)
  if (length(idx1)) event_day[idx1] <- astart + floor(runif(length(idx1)) * Lpre[idx1])
  
  # risk (choose window 1 vs 2 weighted by lengths)
  idx2 <- which(block_int == 2L & Lrisk > 0L) #find the indices for those whose event occurs in the risk windows (block 2) and have a non-zero length of risk windows (Lrisk > 0)
  if (length(idx2)) {
    w1 <- len1[idx2]; w2 <- len2[idx2]; W <- w1 + w2 #calculate the weights for each risk window (len1 and len2) and the total weight (W)
    u2 <- runif(length(idx2)) #draw a uniform random number for each individual in the risk windows to determine which risk window the event will occur in
    pick1 <- (W > 0) & (u2 < w1 / W) #determine which individuals will have their event occur in risk window 1 (dose 1) based on the weights calculated above.  
    pick2 <- (W > 0) & (!pick1) #if did not pick risk window 1, then they will pick risk window 2 (dose 2)
    if (any(pick1)) {
      L1draw <- len1[idx2[pick1]] #length of risk window 1 for those who picked it
      event_day[idx2[pick1]] <- S1[idx2[pick1]] + floor(runif(sum(pick1)) * L1draw) #draw a uniform random number for each individual who picked risk window 1 to determine the day of the event within that window. The event day is calculated as the start of risk window 1 (S1) plus a random number between 0 and the length of risk window 1 (L1draw)
    }
    if (any(pick2)) {
      L2draw <- len2[idx2[pick2]] #length of risk window 2 for those who picked it
      event_day[idx2[pick2]] <- S2[idx2[pick2]] + floor(runif(sum(pick2)) * L2draw) #draw a uniform random number for each individual who picked risk window 2 to determine the day of the event within that window. The event day is calculated as the start of risk window 2 (S2) plus a random number between 0 and the length of risk window 2 (L2draw)
    }
  }
  # post (gap vs tail) #same process as above for the post-risk window, but now we have two sub-windows (gap and tail) to choose from. The gap is the period between the end of risk window 1 and the start of risk window 2, while the tail is the period after the end of risk window 2 to the end of observation period.
  idx3 <- which(block_int == 3L & Lpost > 0L)
  if (length(idx3)) {
    wgap <- Lgap[idx3]; wtail <- Lpost_tail[idx3]; Wp <- wgap + wtail
    u3 <- runif(length(idx3))
    choose_gap  <- (Wp > 0) & (u3 < wgap / Wp)
    choose_tail <- (Wp > 0) & (!choose_gap)
    if (any(choose_gap)) {
      Lg <- Lgap[idx3[choose_gap]]
      event_day[idx3[choose_gap]] <- gap_start[idx3[choose_gap]] + floor(runif(sum(choose_gap)) * Lg)
    }
    if (any(choose_tail)) {
      Lt <- Lpost_tail[idx3[choose_tail]]
      event_day[idx3[choose_tail]] <- post_start[idx3[choose_tail]] + floor(runif(sum(choose_tail)) * Lt)
    }
  }
  
  # SCCS-ready base
  sim <- data.frame(
    indiv = 1:n,
    astart = astart, aend = aend,
    adrug1 = dose1,  aedrug1 = r1e,
    adrug2 = dose2,  aedrug2 = r2e,
    aevent = event_day
  )
  
  return(sim)
}
