simulate_base <- function(RI, years, seed = 2536) {
  set.seed(seed)
  
  # === simulate ===
  
  # --- Inputs ---
  n      <- if (years == 1) 500 else 2000
  astart <- 0
  aend   <- years * 425
  startdose1 <- 60
  enddose1   <- if (years == 1) 425 else 1460
  risklen <- 28L
  
  
  # left-skewed integers on [lo, hi] via Beta(shape1 > shape2)
  rleft_skew_int <- function(n, lo, hi, shape1 = 1.5, shape2 = 4) {
    u <- stats::rbeta(n, shape1, shape2)   # mass near 0 -> earlier days
    x <- lo + floor(u * (hi - lo + 1L))
    pmin(x, hi)
  }
  
  dose1 <- rleft_skew_int(n, startdose1, enddose1, shape1 = 1.5, shape2 = 4)
  gap12 <- sample(60:90, n, replace = TRUE)             # preserves ≥30 gap > risklen=28
  dose2 <- dose1 + gap12
  
  # untruncated risk windows (may exceed aend)
  r1s <- dose1; r1e <- dose1 + risklen
  r2s <- dose2; r2e <- dose2 + risklen
  
  # intersections with [astart, aend] for Person-Time accounting
  S1 <- pmax(astart, r1s); E1 <- pmin(aend, r1e); len1 <- pmax(0L, E1 - S1 + 1L)
  S2 <- pmax(astart, r2s); E2 <- pmin(aend, r2e); len2 <- pmax(0L, E2 - S2 + 1L)
  Lrisk <- len1 + len2
  
  Lpre <- pmax(0L, dose1 - astart)
  gap_start <- pmax(astart, E1 + 1L)
  gap_end   <- pmin(aend,   S2 - 1L)
  Lgap      <- pmax(0L, gap_end - gap_start + 1L)
  post_start <- pmax(astart, E2 + 1L)
  post_end   <- aend
  Lpost_tail <- pmax(0L, post_end - post_start + 1L)
  Lpost      <- Lgap + Lpost_tail
  
  den <- Lpre + RI * Lrisk + Lpost
  p1 <- ifelse(den > 0, Lpre / den, 0)
  p2 <- ifelse(den > 0, (RI * Lrisk) / den, 0)
  
  u <- runif(n)
  block_int <- ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L))
  block <- factor(block_int, levels = 1:3, labels = c("pre","risk","post"))
  
  # draw event day uniformly within chosen block
  event_day <- integer(n)
  # pre
  idx1 <- which(block_int == 1L & Lpre > 0L)
  if (length(idx1)) event_day[idx1] <- astart + floor(runif(length(idx1)) * Lpre[idx1])
  # risk (choose window 1 vs 2 weighted by lengths)
  idx2 <- which(block_int == 2L & Lrisk > 0L)
  if (length(idx2)) {
    w1 <- len1[idx2]; w2 <- len2[idx2]; W <- w1 + w2
    u2 <- runif(length(idx2))
    pick1 <- (W > 0) & (u2 < w1 / W)
    pick2 <- (W > 0) & (!pick1)
    if (any(pick1)) {
      L1draw <- len1[idx2[pick1]]
      event_day[idx2[pick1]] <- S1[idx2[pick1]] + floor(runif(sum(pick1)) * L1draw)
    }
    if (any(pick2)) {
      L2draw <- len2[idx2[pick2]]
      event_day[idx2[pick2]] <- S2[idx2[pick2]] + floor(runif(sum(pick2)) * L2draw)
    }
  }
  # post (gap vs tail)
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
