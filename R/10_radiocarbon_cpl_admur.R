library(dplyr)
library(rcarbon)
library(ADMUR)
library(patchwork)
library(ggplot2)
library(DEoptimR)

set.seed(42)

c14 <- read.csv(here::here("analysis/data/raw_data/radiocarbon.csv"))
c14 <- c14 %>% filter(context != "Food crust")

names(c14) <- c("site", "context", "material", "lab_ref", "age", "sd")

CalArray <- makeCalArray(intcal20, calrange = c(minage, maxage))
PD <- phaseCalibrator(c14, CalArray, remove.external = TRUE)

SPD <- as.data.frame( rowSums(PD) )
# normalise
SPD <- SPD/( sum(SPD) * CalArray$inc )


# CPL parameters must be between 0 and 1, and an odd length.
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
BIC.log <- 2*log(ssize) - 2*(-log$value)
BIC.uniform <- 0 - 2*(-uniform)

bics <- data.frame(model = c("BIC.1", "BIC.2", "BIC.3",
                             "BIC.4", "BIC.exp", "BIC.log", "BIC.uniform"),
                   bic = c(BIC.1, BIC.2, BIC.3, BIC.4, BIC.exp,
                           BIC.log, BIC.uniform))

ggplot(data = bics, aes(x = model, y = bic)) + geom_point()

plotPD(SPD)
cols <- c('firebrick','orchid2','coral2','steelblue','goldenrod3')
lines(CPL1$year, CPL1$pdf*-1, col=cols[1], lwd=2)
lines(CPL2$year, CPL2$pdf*-1, col=cols[2], lwd=2)
lines(CPL3$year, CPL3$pdf*-1, col=cols[3], lwd=2)
lines(CPL4$year, CPL4$pdf*-1, col=cols[4], lwd=2)
lines(EXP$year, EXP$pdf*-1, col=cols[5], lwd=2)


CPL1$name <- "1-CPL"
CPL2$name <- "2-CPL"
CPL3$name <- "3-CPL"
CPL4$name <- "4-CPL"
EXP$name <- "Exponential"
LOG$name <- "Logistic"
UNI$name <- "Uniform"

rcarbonmodels <- rbind(CPL1, CPL2, CPL3, CPL4, LOG, EXP, UNI)

ggplot() +
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
    "3-CPL" = "#ffff33",
    "4-CPL" = "#a65628",
    "5-CPL" = "black"
  )) +
  # scale_x_continuous(limit = c(-10000, -2750)) +
  # scale_y_continuous(limit = c(0, 0.001)) +
  scale_x_reverse(labels = function(x)(x-2000)*-1) +
  theme_bw() +
  theme(legend.title = element_blank())

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

# Plot for inspection (setting up plots below that follow same style as
# that to be used with shoreline dates)
plotSimulationSummary(expsum)
plotSimulationSummary(logsum)
plotSimulationSummary(unisum)

# logsum$n.dates.effective
# round(logsum$n.phases.effective)
expbooms <- expsum$timeseries$calBP[expsum$timeseries$index == 1]
expbusts <- expsum$timeseries$calBP[expsum$timeseries$index == -1]

expplt <- ggplot(data = expsum$timeseries, aes(x = calBP)) +
  geom_vline(xintercept = expbusts,
             col = "firebrick", alpha = 0.06) +
  geom_vline(xintercept = expbooms,
             col = "darkgreen", alpha = 0.06) +
  geom_ribbon(aes(ymin = expsum$timeseries$`2.5%`,
                  ymax = expsum$timeseries$`97.5%`),
              fill = "grey60", alpha = 0.8) +
  geom_line(aes(y = expsum$timeseries$SPD), linewidth = 0.5) +
  geom_line(aes(y = expsum$timeseries$model),
            linewidth = 0.5, col = "red") +
  labs(x = "BCE", y = "Summed probability") +
  scale_x_reverse(limits = c(12000, 4500),
                  breaks = seq(12000, 4500, -1000),
                  expand = expansion(mult = c(0, 0)),
                  labels = function(x)(x-2000)*-1) +
  geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
                label = paste("p =", expsum$pvalue))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_bw()

logbooms <- logsum$timeseries$calBP[logsum$timeseries$index == 1]
logbusts <- logsum$timeseries$calBP[logsum$timeseries$index == -1]

logplt <- ggplot(data = logsum$timeseries, aes(x = calBP)) +
  geom_vline(xintercept = logbusts,
             col = "firebrick", alpha = 0.06) +
  geom_vline(xintercept = logbooms,
             col = "darkgreen", alpha = 0.06) +
  geom_ribbon(aes(ymin = logsum$timeseries$`2.5%`,
                  ymax = logsum$timeseries$`97.5%`),
                fill = "grey60", alpha = 0.8) +
  geom_line(aes(y = logsum$timeseries$SPD), linewidth = 0.5) +
  geom_line(aes(y = logsum$timeseries$model),
            linewidth = 0.5, col = "red") +
  labs(x = "BCE", y = "Summed probability") +
  scale_x_reverse(limits = c(12000, 4500),
                  breaks = seq(12000, 4500, -1000),
                  expand = expansion(mult = c(0, 0)),
                  labels = function(x)(x-2000)*-1) +
  # geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
  #               label = paste("p =", round(lognull$pval, 3)))) +
  geom_text(aes(-Inf, Inf, hjust = 9.5, vjust = 1.5,
                label = paste("p =", logsum$pvalue))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_bw()

unibooms <- unisum$timeseries$calBP[unisum$timeseries$index == 1]
unibusts <- unisum$timeseries$calBP[unisum$timeseries$index == -1]

uniplt <- ggplot(data = unisum$timeseries, aes(x = calBP)) +
  geom_vline(xintercept = unibusts,
             col = "firebrick", alpha = 0.06) +
  geom_vline(xintercept = unibooms,
             col = "darkgreen", alpha = 0.06) +
  geom_ribbon(aes(ymin = unisum$timeseries$`2.5%`,
                  ymax = unisum$timeseries$`97.5%`),
              fill = "grey60", alpha = 0.8) +
  geom_line(aes(y = unisum$timeseries$SPD), linewidth = 0.5) +
  labs(x = "BCE", y = "Summed probability") +
  geom_line(aes(y = unisum$timeseries$model),
            linewidth = 0.5, col = "red") +
  scale_x_reverse(limits = c(12000, 4500),
                  breaks = seq(12000, 4500, -1000),
                  expand = expansion(mult = c(0, 0)),
                  labels = function(x)(x-2000)*-1) +
  # geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
  #               label = paste("p =", round(uninull$pval, 3)))) +
  geom_text(aes(-Inf, Inf, hjust = 9.5, vjust = 1.5,
                label = paste("p =", unisum$pvalue))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_bw()


expplt + logplt + uniplt

