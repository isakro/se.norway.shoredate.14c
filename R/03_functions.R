# Reverse a single calendar date to a corresponding elevation
# Repurposed from Bchron::unCalibrate()
reverse_shoredate <- function(shoreline_date,
                              displacement_curve,
                              elev_reso){

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
    if(minelev <= 0){
      minelev <- elev_reso
    }

    return(sample(seq(minelev, maxelev, elev_reso), 1))
  } else{
    return(NA)
  }
}
