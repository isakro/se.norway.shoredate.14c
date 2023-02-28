# Reverse a single calendar date to a corresponding elevation
# Repurposed from Bchron::unCalibrate()
reverse_shoredate <- function(shoreline_date,
                              displacement_curve,
                              elev_reso){

  if(inherits(displacement_curve, "list")){
    displacement_curve <- rbind.data.frame(displacement_curve)
  }

  lwelev <- round(stats::approx(displacement_curve$bce,
                                displacement_curve$lowerelev,
                                xout = shoreline_date,
                                rule = 1
  )$y, 1)

  upelev <- round(stats::approx(displacement_curve$bce,
                                displacement_curve$upperelev,
                                xout = shoreline_date,
                                rule = 1
  )$y, 1)

  if(!(is.na(lwelev) | is.na(upelev))){
    # In case the displacement curves intersect
    minelev <- min(lwelev, upelev)
    maxelev <- max(lwelev, upelev)

    # Make sure elevation is never 0 (i.e. present day sea-level)
    if(minelev < elev_reso){
      minelev <- elev_reso
    }
    if(maxelev < elev_reso){
      maxelev <- elev_reso
    }

    return(sample(seq(minelev, maxelev, elev_reso), 1))
  } else {
    return(NA)
  }
}


simulation_summary <- function(spd, simulation_results, cut_off, nsim){

  # Reformat data frame of normalised simulated SPDs
  sim <- simulation_results %>%
    filter(bce <= cut_off) %>%
    group_by(bce) %>%
    dplyr::reframe(prob_sum_normalised, simn) %>%
    tidyr::pivot_wider(names_from = simn, values_from = prob_sum_normalised) %>%
    tibble::column_to_rownames("bce")

  # Format observed SPD and retrieve the probability
  fSPD <- as.data.frame(spd) %>%
    filter(sum.bce <= cut_off) %>%
    select(sum.probability) %>%
    unlist(use.names = FALSE)

  Zsim <- t(apply(sim, 1, scale))
  zLo <- apply(Zsim, 1, quantile, prob = 0.025, na.rm = TRUE)
  zHi <- apply(Zsim, 1, quantile, prob = 0.975, na.rm = TRUE)

  Zscore_empirical <- (fSPD - apply(sim, 1, mean))/apply(sim, 1, sd)
  busts <- which(Zscore_empirical < zLo)
  booms <- which(Zscore_empirical > zHi)

  observedStatistic <- sum(c(zLo[busts] - Zscore_empirical[busts]),
                           c(Zscore_empirical[booms] - zHi[booms]))
  expectedstatistic <- abs(apply(Zsim, 2, function(x, y){a = x - y; i = which(a < 0); return(sum(a[i]))}, y = zLo)) +
    apply(Zsim, 2, function(x, y){a = x - y; i = which(a > 0); return(sum(a[i]))}, y = zHi)


  pvalue <- c(length(expectedstatistic[expectedstatistic > observedStatistic]) + 1)/c(nsim + 1)

  return(list(pvalue = pvalue,
              busts = busts,
              booms = booms))
}
