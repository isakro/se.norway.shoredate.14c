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

names(c14) <- c("site", "context", "material", "lab_ref", "age", "sd")

minage <- 4500
maxage <- 10555

CalArray <- makeCalArray(intcal20, calrange = c(minage, maxage))
PD <- phaseCalibrator(c14, CalArray, remove.external = TRUE)

SPD <- as.data.frame(rowSums(PD))
# normalise
SPD <- SPD/(sum(SPD) * CalArray$inc)

CPL.1 <- JDEoptim(lower = 0, upper = 1, fn = objectiveFunction,
                  PDarray = PD, type = 'CPL', NP = 20)
CPL.2 <- JDEoptim(lower = rep(0, 3), upper = rep(1, 3),
                  fn = objectiveFunction, PDarray = PD, type = 'CPL', NP = 60)
CPL.3 <- JDEoptim(lower = rep(0, 5), upper = rep(1, 5),
                  fn = objectiveFunction, PDarray = PD, type = 'CPL', NP = 100)
CPL.4 <- JDEoptim(lower = rep(0, 7), upper = rep(1, 7),
                  fn = objectiveFunction, PDarray = PD, type = 'CPL', NP = 140)
exp <- JDEoptim(lower = -0.01, upper = 0.01, fn = objectiveFunction,
                PDarray = PD, type = 'exp', NP = 20)
log <- JDEoptim(lower = c(0, 0000), upper = c(1, 10000), fn = objectiveFunction,
                PDarray = PD, type = 'logistic', NP = 40)
uniform <- objectiveFunction(NULL, PD, type = 'uniform')

save(exp, log, uniform, CPL.1, CPL.2, CPL.3, CPL.4,
     file = here("analysis/data/derived_data/rcarbon_models.RData"))
load(file = here("analysis/data/derived_data/rcarbon_models.RData"))

CPL1 <- convertPars(pars=CPL.1$par, years=minage:maxage, type='CPL')
CPL2 <- convertPars(pars=CPL.2$par, years=minage:maxage, type='CPL')
CPL3 <- convertPars(pars=CPL.3$par, years=minage:maxage, type='CPL')
CPL4 <- convertPars(pars=CPL.4$par, years=minage:maxage, type='CPL')
EXP <- convertPars(pars=exp$par, years=minage:maxage, type='exp')
LOG <- convertPars(pars = log$par, years = minage:maxage, type = 'logistic')
UNI <- convertPars(pars = NULL, years = minage:maxage, type = "uniform")

data.frame(L1= -CPL.1$value,
           L2= -CPL.2$value,
           L3= -CPL.3$value,
           L4= -CPL.4$value,
           Lexp= -exp$value,
           Llog = -log$value,
           Lunif= -uniform)

BIC.1 <- 1*log(155) - 2*(-CPL.1$value)
BIC.2 <- 3*log(155) - 2*(-CPL.2$value)
BIC.3 <- 5*log(155) - 2*(-CPL.3$value)
BIC.4 <- 7*log(155) - 2*(-CPL.4$value)
BIC.exp <- 1*log(155) - 2*(-exp$value)
BIC.log <- 2*log(155) - 2*(-log$value)
BIC.uniform <- 0 - 2*(-uniform)

bics <- data.frame(model = c("1-CPL", "2-CPL", "3-CPL",
                             "4-CPL", "Exponential", "Logistic", "Uniform"),
                   bic = c(BIC.1, BIC.2, BIC.3, BIC.4, BIC.exp,
                           BIC.log, BIC.uniform))

bicplt <- ggplot(data = bics, aes(x = bic, y = model)) + geom_point(size = 2) +
  labs(x = "BIC", y = "Model") +
  theme_bw()

CPL1$name <- "1-CPL"
CPL2$name <- "2-CPL"
CPL3$name <- "3-CPL"
CPL4$name <- "4-CPL"
EXP$name <- "Exponential"
LOG$name <- "Logistic"
UNI$name <- "Uniform"

rcarbonmodels <- rbind(CPL1, CPL2, CPL3, CPL4,  LOG, EXP, UNI)

mplt <- ggplot() +
  geom_bar(aes(x = as.numeric(rownames(SPD)),
               SPD[,1]),
           stat = "identity", col = "grey") +
  geom_line(aes(x = as.numeric(rownames(SPD)),
                y = SPD[,1])) +
  geom_line(data = rcarbonmodels, aes(year, pdf, col = name),
            linewidth = 1) +
  labs(x = "BCE", y = "Summed probability") +
  scale_colour_manual(values = c(
    "Uniform" = "#e41a1c",
    "Logistic" = "#377eb8",
    "Exponential" = "#4daf4a",
    "1-CPL" = "#984ea3",
    "2-CPL" = "#ff7f00",
    "3-CPL" = "#e0d616",
    "4-CPL" = "#a65628",
    "5-CPL" = "black"
  )) +
  # scale_x_continuous(limit = c(-10000, -2750)) +
  # scale_y_continuous(limit = c(0, 0.001)) +
  scale_x_reverse(labels = function(x)(x-1950)*-1) +
  theme_bw() +
  theme(legend.title = element_blank())

mplt + bicplt

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

save(expsum, logsum, unisum,
     file = here("analysis/data/derived_data/admur_model_test.RData"))
load(file = here("analysis/data/derived_data/admur_model_test.RData"))

# Plot for inspection (setting up plots below that follow same style as
# that to be used with shoreline dates)
plotSimulationSummary(expsum)
plotSimulationSummary(logsum)
plotSimulationSummary(unisum)

expplotr <- plot_mc(expsum)
logplotr <- plot_mc(logsum)
uniplotr <- plot_mc(unisum)

load("../external_data/shorespd/expplts.rda")
load("../external_data/shorespd/logplts.rda")
load("../external_data/shorespd/uniplts.rda")

(expplts + logplts + uniplts)/
(expplotr + logplotr + uniplotr) +
  plot_annotation(tag_levels = "A")

ggsave(here::here("analysis/figures/mc.png"),
       units = "px", width = 4000, height = 1250*2)

# Save pd, logistic model and min and max age for mcmc
save(log, PD, minage, maxage,
     file = here("analysis/data/derived_data/rcarbon_mcmc_data.RData"))
