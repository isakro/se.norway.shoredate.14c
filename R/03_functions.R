# This script defines functions used for Monte Carlo simulation and plotting of
# the results.

# Reverse a single calendar date to a corresponding elevation
# Repurposed from Bchron.
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

# Return summary of a Monte Carlo simulation run. Repurposed from ADMUR
simulation_summary <- function(spd, simulation_results,
                               cut_offs, mod, sample_size){

  # spd needs to be formated for ADMUR
  SPD.obs <- spd[,1]

  # Format simulation results for ADMUR
  SPD.sims <- simresults %>%
    dplyr::filter(bce <= cut_offs[1] & bce >= cut_offs[2]) %>%
    dplyr::mutate(bp = (bce * -1) + 1950) %>%
    arrange(bp) %>%
    dplyr::select(-bce, -prob_sum) %>%
    dplyr::group_by(bp) %>%
    dplyr::filter(!is.na(prob_sum_normalised)) %>%
    dplyr::group_by(simn) %>%
    tidyr::pivot_wider(names_from = simn, values_from = prob_sum_normalised) %>%
    tibble::column_to_rownames("bp")


  calBP <- as.numeric(rownames(spd))

  # Code below repurposed from ADMUR

  # expected simulation
  expected.sim <- rowMeans(SPD.sims)

  # local standard deviation
  SD <- apply(SPD.sims, 1, sd)

  CI <- t(apply(SPD.sims,1,quantile,prob=c(0.025, 0.125, 0.25,
                                           0.75, 0.875, 0.975)))
  upper.95 <- CI[, dimnames(CI)[[2]] == "97.5%"]
  lower.95 <- CI[, dimnames(CI)[[2]] == "2.5%"]
  index <- as.numeric(SPD.obs >= upper.95) - as.numeric(SPD.obs <= lower.95)
  SS.obs <- sum(SPD.obs > upper.95 | SPD.obs < lower.95) / length(SPD.obs)

  N <- ncol(SPD.sims)
  SS.sims <- numeric(N)
  for(n in 1:N){
    SPD.tmp <- SPD.sims[,n]
    SS.sims[n] <- sum(SPD.tmp > upper.95 | SPD.tmp < lower.95) / length(SPD.obs)
  }
  pvalue <- sum(SS.sims >= SS.obs)/N

  timeseries <- cbind(data.frame(calBP = calBP, expected.sim = expected.sim,
                                 local.sd = SD, model = mod, SPD = SPD.obs,
                                 index = index), CI)

  n.dates.all <- n.dates.effective <- n.phases.all <- n.phases.effective <-
    n.phases.internal <- sample_size

  return(list(timeseries=timeseries,
              pvalue=pvalue,
              observed.stat=SS.obs,
              simulated.stat=SS.sims,
              n.dates.all=n.dates.all,
              n.dates.effective=n.dates.effective,
              n.phases.all=n.phases.all,
              n.phases.effective=n.phases.effective,
              n.phases.internal=n.phases.internal))
}

# Custom function for plotting results of a Monte Carlo simulation.
# Based partially on rcarbon.
plot_mc <- function(mc_summary, xppos, yppos, title = NULL) {
  booms <- mc_summary$timeseries$calBP[mc_summary$timeseries$index == 1]
  busts <- mc_summary$timeseries$calBP[mc_summary$timeseries$index == -1]

  spdmax <- max(mc_summary$timeseries$SPD)
  spdmin <- min(mc_summary$timeseries$SPD)
  spdadj <- (spdmax - spdmin) * 1.1

  p <- round(mc_summary$pvalue, 3)
  plabel <- ifelse(p == 0, "p < 0.001", paste("p =", p))

  ggplot(data = mc_summary$timeseries, aes(x = calBP)) +
    geom_vline(xintercept = busts,
               col = "firebrick", alpha = 0.05) +
    geom_vline(xintercept = booms,
               col = "darkgreen", alpha = 0.05) +
    geom_ribbon(aes(ymin = `2.5%`,
                    ymax = `97.5%`),
                fill = "grey60", alpha = 0.8) +
    geom_line(aes(y = model), linewidth = 0.5, col = "red") +
    geom_line(aes(y = SPD)) +
    labs(x = "BCE", y = "Summed probability", title = title) +
    scale_x_reverse(limits = c(11950, 4450),
                    breaks = seq(11950, 4450, -1000),
                    expand = expansion(mult = c(0, 0)),
                    labels = function(x)(x-1950)*-1) +
    geom_text(aes(11000, spdadj), label = plabel) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                       labels = scales::label_comma()) +
    theme_bw()
}
