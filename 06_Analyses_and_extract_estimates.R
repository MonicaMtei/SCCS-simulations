## --- define the 5 analysis methods ---
fitters <- list(
  standard = function(dat) standardsccs(
    event  ~ adrug1,
    indiv  = indiv,
    astart = astart,
    aend   = aend,
    aevent = aevent,
    adrug  = cbind(adrug1, adrug2),
    aedrug = cbind(aedrug1, aedrug2),
    dataformat = "multi",
    sameexpopar = FALSE,
    data   = dat
  ),
  pre14 = function(dat) standardsccs(
    event  ~ adrug1,
    indiv  = indiv,
    astart = astart,
    aend   = aend,
    aevent = aevent,
    adrug  = cbind(adrug1, adrug2),
    aedrug = cbind(aedrug1, aedrug2),
    expogrp = c(-14, 0),
    dataformat = "multi",
    sameexpopar = FALSE,
    data   = dat
  ),
  pre28 = function(dat) standardsccs(
    event  ~ adrug1,
    indiv  = indiv,
    astart = astart,
    aend   = aend,
    aevent = aevent,
    adrug  = cbind(adrug1, adrug2),
    aedrug = cbind(aedrug1, aedrug2),
    expogrp = c(-28, 0),
    dataformat = "multi",
    sameexpopar = FALSE,
    data   = dat
  ),
  eventdep = function(dat) eventdepenexp(
    indiv  = indiv,
    astart = astart,
    aend   = aend,
    aevent = aevent,
    adrug  = cbind(adrug1, adrug2),
    aedrug = cbind(aedrug1, aedrug2),
    dataformat = "multi",
    sameexpopar = FALSE,
    data   = dat
  ),
  
  pre60 = function(dat) {
    # build pre-windows safely: start = max(astart, candidate), end = start-1 of dose
    adrug1pre <- pmax(dat$adrug1 - 60L, dat$astart)
    adrug2pre_raw <- dat$adrug2 - 60L
    # if risk1 overlaps into that pre2 window, start pre2 the day after risk1 ends
    adrug2pre <- pmax(adrug2pre_raw, dat$aedrug1 + 1L)
    adrug2pre <- pmax(adrug2pre, dat$astart)
    
    standardsccs(
      event  ~ adrug1 + adrug1pre + adrug2 + adrug2pre,
      indiv  = indiv,
      astart = astart,
      aend   = aend,
      aevent = aevent,
      adrug  = cbind(adrug1, adrug1pre, adrug2, adrug2pre),
      aedrug = cbind(aedrug1, adrug1 - 1L, aedrug2, adrug2 - 1L),
      dataformat   = "stack",
      sameexpopar  = FALSE,
      data   = dat
    )
  }
)

#extract estimates

extract_ri_se <- function(fit, approach) {
  co <- tryCatch(fit[["coefficients"]], error = function(e) NULL)
  if (is.null(co)) {
    return(tibble::tibble(
      approach    = approach,
      dose        = c("dose1","dose2"),
      log_RI_est  = NA_real_,   # coef on log scale
      RI_est      = NA_real_,   # exp(coef)
      SE_est      = NA_real_    # SE on log scale
    ))
  }
  
  cn <- colnames(co)
  
  # log-scale coefficient
  logv <- if ("coef" %in% cn) {
    co[, "coef"]
  } else if ("exp(coef)" %in% cn) {
    log(co[, "exp(coef)"])
  } else {
    stop("Couldn't find coef or exp(coef) in coefficients table.")
  }
  
  # RI on original scale
  expv <- if ("exp(coef)" %in% cn) co[, "exp(coef)"] else exp(logv)
  
  # SE on log scale
  sev <- if ("se(coef)" %in% cn) co[, "se(coef)"] else co[, "Std. Error"]
  
  # Pick the rows that correspond to dose1 / dose2 for this approach
  idx <- switch(approach,
                standard = c(1, 2),
                eventdep = c(1, 2),
                pre14    = c(2, 4),
                pre28    = c(2, 4),
                pre60    = c(1, 3),
                c(1, 2)
  )
  idx <- idx[idx <= nrow(co)]
  
  tibble::tibble(
    approach   = approach,
    dose       = c("dose1","dose2")[seq_along(idx)],
    log_RI_est = as.numeric(logv[idx]),   # <-- log scale coef
    RI_est     = as.numeric(expv[idx]),   # <-- exp(coef)
    SE_est     = as.numeric(sev[idx])     # <-- SE on log scale
  )
}
