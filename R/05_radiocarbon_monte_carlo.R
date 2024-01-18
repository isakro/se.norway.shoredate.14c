# This script performs Monte Carlo simulation using the radiocarbon data and
# plots results from both simulations with both shoreline and radiocarbon dates.
# Based on the failed rejection of the logistic model, a version of the model
# fit when accounting for taphonomy is then created (figure included in
# as supplementary material)

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

minage <- min(c14spd$grid[c14spd$grid$PrDens > 0,]$calBP)
maxage <- max(c14spd$grid[c14spd$grid$PrDens > 0,]$calBP)

# c14spd_taph <- transformSPD(c14spd)

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

# Load MC plots for shoreline dates
load(here("analysis/data/derived_data/expplts.rda"))
load(here("analysis/data/derived_data/logplts.rda"))
load(here("analysis/data/derived_data/uniplts.rda"))

(expplts + logplts + uniplts)/
(expplotr + logplotr + uniplotr) +
  plot_annotation(tag_levels = "a")

ggsave(here::here("analysis/figures/mc.png"),
       units = "px", width = 6000, height = 3750, dpi = 400)


#### Logistic model with taphonomic loss ####

# Use the default parameter settings from ADMUR, as there are no available
# prior data on which to adjust these from the region (see the ADMUR vignette)
logtaph <- JDEoptim(lower = c(0, 0000, 0, -3), upper = c(1, 10000, 20000,0),
                    fn = objectiveFunction, PDarray = PD, type = 'logistic',
                    NP = 40, taphonomy = TRUE)

# Retrieve the logistic model fit without accounting for taphonomy
pop_log <- convertPars(pars=log$par,
                       years=minage:maxage,
                       type='logistic')

# Retrieve the population dynamics component of the logistic model fit when
# accounting for taphonomy
pop_log_taph <- convertPars(pars=logtaph$par[1:2],
                            years=minage:maxage,
                            type='logistic')

# Retrieve the taphonomic component
taph_log_taph <- convertPars(pars=logtaph$par[3:4],
                            years=minage:maxage,
                            type='power')

# Assemble the SPD
SPD <- as.data.frame(rowSums(PD))
SPD <- SPD/( sum(SPD) *5 )

popplt <- ggplot() +
  geom_line(data = SPD, aes(x = as.numeric(rownames(SPD)), SPD[,1])) +
  geom_line(data = pop_log, aes(x = year, pdf, col = "pop")) +
  geom_line(data = pop_log_taph, aes(x = year, y = pdf, col = "poptaph")) +
  labs(x = "BCE", y = "Summed probability",
       title = "Population dynamics component") +
  scale_colour_manual(name = "",
                      labels = c("Logistic model with taphonomy",
                                 "Logistic model"),
                      breaks = c("poptaph", "pop"),
                      values = c("red", "blue")) +
  scale_x_reverse(limits = c(11950, 4450),
                  breaks = seq(11950, 4450, -1000),
                  expand = expansion(mult = c(0, 0)),
                  labels = function(x)(x-1950)*-1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     labels = scales::label_comma()) +
  theme_bw() +
  theme(legend.position = "bottom")

taphplt <- ggplot() +
  geom_line(data = SPD, aes(x = as.numeric(rownames(SPD)), SPD[,1]), col = NA) +
  geom_line(data = taph_log_taph, aes(x = year, pdf), col = "green4") +
  labs(x = "BCE", y = "Probability distribution",
       title = "Taphonomic component") +
  scale_x_reverse(limits = c(11950, 4450),
                  breaks = seq(11950, 4450, -1000),
                  expand = expansion(mult = c(0, 0)),
                  labels = function(x)(x-1950)*-1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     labels = scales::label_comma()) +
  theme_bw()

popplt + taphplt +
  plot_annotation(tag_levels = 'a')

ggsave(here::here("analysis/figures/taphonomic.png"),
       units = "px", width = 4500, height = 2250, dpi = 400)
