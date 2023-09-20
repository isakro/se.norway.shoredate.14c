# This script performs Monte Carlo simulation using the radiocarbon data.

library(dplyr)
library(rcarbon)
library(ADMUR)
library(patchwork)
library(ggplot2)
library(DEoptimR)
library(here)

set.seed(42)

source(here("R/03_functions.R"))

c14 <- read.csv(here::here("analysis/data/raw_data/radiocarbon.csv"))
c14 <- c14 %>% filter(context != "Food crust")

# Use the calibration procedure from rcarbon to identify min and max ages
caldates <- calibrate(x = c14$c14_bp, normalised = FALSE,
                      errors = c14$error, calCurves = "intcal20")
c14spd <- spd(caldates, timeRange = c(12000, 4500))

max(c14spd$grid[c14spd$grid$PrDens > 0,]$calBP)
min(c14spd$grid[c14spd$grid$PrDens > 0,]$calBP)

minage <- min(c14spd$grid[c14spd$grid$PrDens > 0,]$calBP)
maxage <- max(c14spd$grid[c14spd$grid$PrDens > 0,]$calBP)


# Switch to ADMUR
names(c14) <- c("site", "context", "material", "lab_ref", "age", "sd")

CalArray <- makeCalArray(intcal20, calrange = c(minage, maxage))
PD <- phaseCalibrator(c14, CalArray, remove.external = TRUE)

exp <- JDEoptim(lower = -0.01, upper = 0.01, fn = objectiveFunction,
                PDarray = PD, type = 'exp', NP = 20)
log <- JDEoptim(lower = c(0, 0000), upper = c(1, 10000), fn = objectiveFunction,
                PDarray = PD, type = 'logistic', NP = 40)
uniform <- objectiveFunction(NULL, PD, type = 'uniform')

save(exp, log, uniform,
     file = here("analysis/data/derived_data/rcarbon_models.RData"))
save(exp, log, uniform, PD,
     file = here("analysis/data/derived_data/rcarbon_models_pd.RData"))

c14$datingType <- '14C'

expsum <- SPDsimulationTest(c14, calcurve = intcal20,
                             calrange = c(minage, maxage), pars = exp$par,
                             type = 'exp', N = 10000)
logsum <- SPDsimulationTest(c14, calcurve = intcal20,
                             calrange = c(minage, maxage), pars = log$par,
                             type = 'logistic', N = 10000)
unisum <- SPDsimulationTest(c14, calcurve = intcal20,
                            calrange = c(minage, maxage), pars = NULL,
                            type = 'uniform', N = 10000)

# Plot for inspection using ADMUR
plotSimulationSummary(expsum)
plotSimulationSummary(logsum)
plotSimulationSummary(unisum)

# Plotting using custom plot function
expplotr <- plot_mc(expsum)
logplotr <- plot_mc(logsum)
uniplotr <- plot_mc(unisum)

# Load MC plots for shoreline dates, stored externally
load("../external_data/shorespd/expplts.rda")
load("../external_data/shorespd/logplts.rda")
load("../external_data/shorespd/uniplts.rda")

(expplts + logplts + uniplts)/
(expplotr + logplotr + uniplotr) +
  plot_annotation(tag_levels = "A")

ggsave(here::here("analysis/figures/mc.png"),
       units = "px", width = 4000, height = 1250*2)
