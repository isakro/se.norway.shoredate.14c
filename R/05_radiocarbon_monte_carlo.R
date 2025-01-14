# This script performs Monte Carlo simulation using the radiocarbon data and
# plots results from simulations with both shoreline and radiocarbon dates.
# The impact of normalising and not normalising the dates is then explored
# (included as a supplementary figure), and, based on the failed rejection of
# the logistic model, a version of the model fit when accounting for taphonomy
# is then explored (included as supplementary figures).

library(dplyr)
library(rcarbon)
library(ADMUR)
library(patchwork)
library(ggplot2)
library(ggnewscale)
library(DEoptimR)
library(here)
library(ggplotify)

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


# Plot using custom plot function
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

##### Fitting models to non-normalised dates #####

# The code below using the rcarbon approach for calibrating and summing without
# normalising, and then fits the models using the ADMUR approach

# Read in 14C-data again
c14 <- read.csv(here::here("analysis/data/raw_data/radiocarbon.csv"))
c14 <- c14 %>% filter(context != "Food crust")

# Calibrate and bin using rcarbon
caldates <- calibrate(x = c14$c14_bp, normalised = FALSE,
                      errors = c14$error, calCurves = "intcal20")
subcal <- subset(caldates, BP <= 12000 & BP >= 4500, p = 0.05)
c14sub <- c14[which.CalDates(caldates, BP <= 12000 &
                                BP >= 4500, p = 0.05),]
subbins = binPrep(sites =  c14sub$site_name, ages = c14sub$c14_bp, h = 200)

# Retrieve dates (based on code from rcarbon::spd())
x <- subcal
binNames <- unique(subbins)

calyears <- data.frame(calBP = seq(12000, 4500, -1))
binnedMatrix <- matrix(NA, nrow= nrow(calyears), ncol = length(binNames))
colnames(binnedMatrix) <- binNames
rownames(binnedMatrix) <- rev(calyears[,1])

# Set up dates calibrated without normalisation using rcarbon for use with ADMUR
for(b in 1:length(binNames)){
  index <- which(subbins == binNames[b])
  slist <- x$grids[index]
  slist <- lapply(slist,FUN = function(x) merge(calyears, x, all.x = TRUE))
  slist <- rapply(slist, f = function(x) ifelse(is.na(x), 0, x),
                  how = "replace")
  slist <- lapply(slist, FUN = function(x) x[with(x, order(-calBP)), ])

  tmp <- do.call(rbind, slist)
  tmp <- aggregate(tmp$PrDens, by = list(Category = tmp$calBP), FUN = sum)
  spdtmp <- tmp[,-1]
  spdtmp <- as.data.frame(spdtmp)
  rownames(spdtmp) <- tmp[, 1]
  names(spdtmp) <- "PrDens"

  binnedMatrix[,b] <- spdtmp[,1]/length(index)
}

pd <- as.data.frame(binnedMatrix)

# Sanity check: Compare rcarbon SPD and SPD when formatted for ADMUR
c14spd <- rcarbon::spd(subcal, spdnormalised = TRUE,
                       bins = subbins, timeRange = c(12000, 4500))
SPD <- as.data.frame(rowSums(pd))
SPD <- SPD/sum(SPD)

par(mfrow = c(2,1))
plot(c14spd)
plotPD(SPD)

### Fit models using ADMUR
# Find min and max ages
maxage <- max(c14spd$grid[c14spd$grid$PrDens > 0,]$calBP)
minage <- min(c14spd$grid[c14spd$grid$PrDens > 0,]$calBP)

# Maximum likelihood search using DEoptimR
logi <- JDEoptim(lower = c(0, 0000), upper = c(1, 10000),
                 fn = objectiveFunction, PDarray = pd, type = 'logistic',
                 NP = 40, trace = TRUE)
expnonnorm <- JDEoptim(lower = -0.01, upper = 0.01, fn = objectiveFunction,
                    PDarray = pd, type = 'exp', NP = 20)
uninonnorm <- objectiveFunction(NULL, pd, type = 'uniform')

save(logi, expnonnorm, uninonnorm,
     file = here("analysis/data/derived_data/nonnorm_models.RData"))
load(here("analysis/data/derived_data/nonnorm_models.RData"))

# Retrieve distribution
log_nonnormalised <- convertPars(pars = logi$par, years = minage:maxage,
                                 type = 'logistic')
log_nonnormalised$model <- "Logistic, non-normalised"

# Retrieve distributions for models fit to SPD based on normalised dates
log_normalised <- convertPars(pars = log$par, years = minage:maxage,
                              type = 'logistic')
log_normalised$model <- "Logistic, normalised"

# Combine with models fit to SPD based on non-normalised dates for plotting
rmodels <- rbind(log_nonnormalised, log_normalised)

# Compile SPD
spd <- as.data.frame(rowSums(PD))
spd <- spd/( sum(spd) *5 )

ggplot() +
  geom_line(data = spd, aes(x = as.numeric(rownames(spd)), spd[,1],
                            colour = "Normalised")) +
  geom_line(data = SPD, aes(x = as.numeric(rownames(SPD)), SPD[,1],
                            colour = "Non-normalised")) +
  scale_colour_manual(name = "RSPD",
                      values = c("Normalised" = "darkgray",
                                 "Non-normalised" = "black")) +
  new_scale_color() +
  geom_line(data = rmodels, aes(year, pdf, col = model),
          linewidth = 1) +
  scale_color_brewer(name = "Models", palette = 'Set1') +
  scale_x_reverse(limits = c(11950, 4450),
                  breaks = seq(11950, 4450, -1000),
                  expand = expansion(mult = c(0, 0)),
                  labels = function(x)(x-1950)*-1) +
  labs(x = "BCE", y = "Summed probability") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     labels = scales::label_comma()) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(here::here("analysis/figures/normalisation.png"),
       units = "px", width = 4500, height = 2250, dpi = 400)

#### Monte Carlo simulation with non-normalised dates using rcarbon but
# models fit with ADMUR ####
exp_nonnormalised <- convertPars(pars = exp$par, years = minage:maxage,
                                 type = 'exp')
exp_nonnormalised$model <- "Exponential, non-normalised"
uni_nonnormalised <- convertPars(pars = NULL, years = minage:maxage,
                                 type = 'uniform')
uni_nonnormalised$model <- "Uniform, non-normalised"

# Format for rcarbon
log_nonnormalised <- log_nonnormalised[,c(1,2)]
names(log_nonnormalised) <- c("calBP", "PrDens")
exp_nonnormalised <- exp_nonnormalised[,c(1,2)]
names(exp_nonnormalised) <- c("calBP", "PrDens")
uni_nonnormalised <- uni_nonnormalised[,c(1,2)]
names(uni_nonnormalised) <- c("calBP", "PrDens")


# rcarbon Monte Carlo tests with non-normalised dates
log_nonnrom_test <- modelTest(x = subcal, errors = c14sub$error,
                         bins = subbins, model="custom",
                         predgrid = log_nonnormalised,
                         nsim = 10000, spdnormalised = TRUE,
                         timeRange = c(maxage, minage))

exp_nonnrom_test <- modelTest(x = subcal, errors = c14sub$error,
                              bins = subbins, model="custom",
                              predgrid = exp_nonnormalised,
                              nsim = 10000, spdnormalised = TRUE,
                              timeRange = c(maxage, minage))

uni_nonnrom_test <- modelTest(x = subcal, errors = c14sub$error,
                              bins = subbins, model="custom",
                              predgrid = uni_nonnormalised,
                              nsim = 10000, spdnormalised = TRUE,
                              timeRange = c(maxage, minage))

save(log_nonnrom_test, exp_nonnrom_test, uni_nonnrom_test,
     file = here("analysis/data/derived_data/modelTest_nonnorm.RData"))
load(here("analysis/data/derived_data/modelTest_nonnorm.RData"))

obs <- log_nonnrom_test$result[, 1:2]
envelope <- log_nonnrom_test$result[, 3:4]
booms <- which(obs$PrDens > envelope[, 2])
busts <- which(obs$PrDens < envelope[, 1])

spdmax <- max(obs$PrDens)
spdmin <- min(obs$PrDens)
spdadj <- (spdmax - spdmin) * 1.1

p <- round(log_nonnrom_test$pval, 3)
plabel <- ifelse(p == 0, "p < 0.001", paste("p =", p))

logi_nonnor_plt <- ggplot(log_nonnrom_test$result, aes(x = calBP)) +
  geom_vline(xintercept = obs$calBP[busts],
             col = "firebrick", alpha = 0.01) +
  geom_vline(xintercept = obs$calBP[booms],
             col = "darkgreen", alpha = 0.01) +
  geom_ribbon(aes(ymin = lo,
                  ymax = hi),
              fill = "grey60", alpha = 0.8) +
  geom_line(aes(y = PrDens)) +
  geom_line(data = log_nonnormalised, aes(y = PrDens),
            col = "red") +
  theme_bw() +
  labs(x = "BCE", y = "Summed probability") +
  scale_x_reverse(limits = c(11950, 4450),
                  breaks = seq(11950, 4450, -1000),
                  expand = expansion(mult = c(0, 0)),
                  labels = function(x)(x-1950)*-1) +
  geom_text(aes(11000, spdadj), label = plabel) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     labels = scales::label_comma())

obs <- exp_nonnrom_test$result[, 1:2]
envelope <- exp_nonnrom_test$result[, 3:4]
booms <- which(obs$PrDens > envelope[, 2])
busts <- which(obs$PrDens < envelope[, 1])

spdmax <- max(obs$PrDens)
spdmin <- min(obs$PrDens)
spdadj <- (spdmax - spdmin) * 1.1

p <- round(exp_nonnrom_test$pval, 3)
plabel <- ifelse(p == 0, "p < 0.001", paste("p =", p))

exp_nonnor_plt <- ggplot(exp_nonnrom_test$result, aes(x = calBP)) +
  geom_vline(xintercept = obs$calBP[busts],
             col = "firebrick", alpha = 0.01) +
  geom_vline(xintercept = obs$calBP[booms],
             col = "darkgreen", alpha = 0.01) +
  geom_ribbon(aes(ymin = lo,
                  ymax = hi),
              fill = "grey60", alpha = 0.8) +
  geom_line(aes(y = PrDens)) +
  geom_line(data = exp_nonnormalised, aes(y = PrDens),
            col = "red") +
  theme_bw() +
  labs(x = "BCE", y = "Summed probability") +
  scale_x_reverse(limits = c(11950, 4450),
                  breaks = seq(11950, 4450, -1000),
                  expand = expansion(mult = c(0, 0)),
                  labels = function(x)(x-1950)*-1) +
  geom_text(aes(11000, spdadj), label = plabel) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     labels = scales::label_comma())

obs <- uni_nonnrom_test$result[, 1:2]
envelope <- uni_nonnrom_test$result[, 3:4]
booms <- which(obs$PrDens > envelope[, 2])
busts <- which(obs$PrDens < envelope[, 1])

spdmax <- max(obs$PrDens)
spdmin <- min(obs$PrDens)
spdadj <- (spdmax - spdmin) * 1.1

p <- round(uni_nonnrom_test$pval, 3)
plabel <- ifelse(p == 0, "p < 0.001", paste("p =", p))

uni_nonnor_plt <- ggplot(uni_nonnrom_test$result, aes(x = calBP)) +
  geom_vline(xintercept = obs$calBP[busts],
             col = "firebrick", alpha = 0.01) +
  geom_vline(xintercept = obs$calBP[booms],
             col = "darkgreen", alpha = 0.01) +
  geom_ribbon(aes(ymin = lo,
                  ymax = hi),
              fill = "grey60", alpha = 0.8) +
  geom_line(aes(y = PrDens)) +
  geom_line(data = uni_nonnormalised, aes(y = PrDens),
            col = "red") +
  theme_bw() +
  labs(x = "BCE", y = "Summed probability") +
  scale_x_reverse(limits = c(11950, 4450),
                  breaks = seq(11950, 4450, -1000),
                  expand = expansion(mult = c(0, 0)),
                  labels = function(x)(x-1950)*-1) +
  geom_text(aes(11000, spdadj), label = plabel) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     labels = scales::label_comma())


(exp_nonnor_plt + logi_nonnor_plt + uni_nonnor_plt) +
  plot_annotation(tag_levels = "a")

ggsave(here::here("analysis/figures/mc_nonnorm.png"),
       units = "px", width = 6000, height = 3750/2, dpi = 400)

#### Logistic model with taphonomic loss ####
# First MCMC of population model without taphonomy

chainpop <- mcmc(PDarray = PD,
                   startPars = log$par,
                   type = 'logistic',
                   N = 100000,
                   burn = 2000,
                   thin = 5,
                   jumps = c(0.011, 200))

print(chainpop$acceptance.ratio)
par(mfrow=c(2,1))
col <- 'steelblue'
for(n in 1:2){
  plot(chainpop$all.pars[,n], type='l', col=col, xlab='', ylab='',
       main=paste('par',n))
}

save(chainpop,
     file = here("../external_data/chain.RData"))

log_pop <- convertPars(pars = chainpop$res[sample(1:nrow(chainpop$res),
                                                  size = 1000), ],
                       years = seq(minage, maxage, 50), type = 'logistic')

# Introduce taphonomy
# Use the default prior from ADMUR, as there are no available
# prior data on which to adjust these from the region
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
                            years=seq(minage, maxage, 5),
                            type='power')

# Assemble the SPD
SPD <- as.data.frame(rowSums(PD))
SPD <- SPD/(sum(SPD)*5)

popplt <- ggplot() +
  geom_line(data = SPD, aes(x = as.numeric(rownames(SPD)), SPD[,1])) +
  geom_line(data = pop_log, aes(x = year, pdf, col = "pop")) +
  geom_line(data = pop_log_taph, aes(x = year, y = pdf, col = "poptaph")) +
  labs(x = "BCE", y = "Summed probability",
       title = "MLE - Population dynamics component") +
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
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

taphplt <- ggplot() +
  geom_line(data = SPD, aes(x = as.numeric(rownames(SPD)), SPD[,1]), col = NA) +
  geom_line(data = taph_log_taph, aes(x = year, pdf), col = "green4") +
  labs(x = "BCE", y = "Probability distribution",
       title = "MLE - Taphonomic component") +
  scale_x_reverse(limits = c(11950, 4450),
                  breaks = seq(11950, 4450, -1000),
                  expand = expansion(mult = c(0, 0)),
                  labels = function(x)(x-1950)*-1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     labels = scales::label_comma()) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

# MCMC search accounting for taphonomy
chain.taph <- mcmc(PDarray = PD,
                   startPars = logtaph$par,
                   type = 'logistic', taphonomy = TRUE,
                   N = 100000,
                   burn = 5000,
                   thin = 5,
                   jumps = c(0.011, 125, 0.011, 0.011))

print(chain.taph$acceptance.ratio)
par(mfrow=c(2,2), mar=c(4,3,3,1))
col <- 'steelblue'
for(n in 1:4){
  plot(chain.taph$all.pars[,n], type='l', col=col, xlab='', ylab='',
       main=paste('par',n))
}

save(chain.taph,
     file = here("../external_data/n_100k_jump_011.RData"))
load(here("../external_data/n_100k_jump_011.RData"))


# Retrieve the demographic component
taph_pop <- convertPars(pars = chain.taph$res[5000:6000, 1:2],
                    years = seq(minage, maxage, 50), type = 'logistic')

# Retrieve the taphonomic component
taph <- convertPars(pars = chain.taph$res[sample(1:nrow(chain.taph$res),
                                                 size = 1000), 3:4],
                    years = seq(1000, 30000, by = 50), type = 'power')

dfpoptaph <- as.data.frame(t(taph_pop))
dfpoptaph$year <- as.numeric(rownames(dfpoptaph))
dfpoptaph <- reshape2::melt(dfpoptaph,  id.vars = 'year',
                            variable.name = 'series')

dflog <- as.data.frame(t(log_pop))
dflog$year <- as.numeric(rownames(dflog))
dflog <- reshape2::melt(dflog,  id.vars = 'year', variable.name = 'series')

popcomp <- ggplot() +
  geom_line(data = dflog, aes(year, value, group = series, col = "jpdpop"),
            alpha = 0.1) +
  geom_line(data = dfpoptaph, aes(year, value, group = series, col = "jpdtaph"),
            alpha = 0.1) +
  geom_line(data = pop_log_taph, aes(x = year, y = pdf, col = "pop"),
            linewidth = 1) +
  geom_line(data = pop_log, aes(x = year, y = pdf, col = "poptaph"),
            linewidth = 1) +
  scale_colour_manual(name = "",
                      labels = c("ML",
                                "ML with taphonomy",
                                "Joint posterior distributions",
                                "Joint posterior distributions with taphonomy"),
                      breaks = c("pop", "poptaph", "jpdpop", "jpdtaph"),
                      values = c("green3", "black", "red", "grey60")) +
  guides(col = guide_legend(override.aes = list(alpha = 1)))  +
  scale_x_reverse() +
  # scale_x_reverse(limits = c(11950, 4450),
  #                 breaks = seq(11950, 4450, -1000),
  #                 expand = expansion(mult = c(0, 0)),
  #                 labels = function(x)(x-1950)*-1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "cal BP", y = "Probability density",
       title = "Population dynamics") +
  theme_bw() +
  theme(legend.position = "bottom")

taphpar <- as.data.frame(chain.taph$res[,3:4])
names(taphpar) <- c("b", "c") # Name following ADMUR convention

taphparplot <- ggplot(taphpar) +
  geom_point(aes(x = b, y = c, col = "jpdp"), alpha = 0.3) +
  geom_point(aes(x = logtaph$par[3], y = logtaph$par[4], col = "mlp"),
             size = 2) +
  scale_colour_manual(name = "",
                      labels = c("Joint posterior distributions",
                                 "ML parameters"),
                      breaks = c("jpdp", "mlp"),
                      values = c("grey60", "black")) +
  labs(x = "Parameter b", y = "Parameter c",
       title = "Taphonomic parameters") +
  guides(col = guide_legend(override.aes = list(alpha = 1)))  +
  theme_bw() +
  theme(legend.position = "bottom")

popcomp / taphparplot +
  plot_annotation(tag_levels = 'a')

ggsave(here::here("analysis/figures/taphonomic.png"),
       units = "px", width = 4000, height = 4500, dpi = 400)

# ML parameters for the demographic component
ml_pars <- round(log$par, 3)
ml_taph_pars <- round(logtaph$par, 3)[1:2]

ci_pop_par1 <- round(quantile(chainpop$res[,1], prob=c(0.025,0.975)), 3)
ci_pop_par2 <- round(quantile(chainpop$res[,2], prob=c(0.025,0.975)))
ci_poptaph_par1 <- round(quantile(chain.taph$res[,1], prob=c(0.025,0.975)), 3)
ci_poptaph_par2 <- round(quantile(chain.taph$res[,2], prob=c(0.025,0.975)))

params <- data.frame(ml_par1 = c(ml_pars[1], ml_taph_pars[1]),
           ci_par1_2_5 = c(ci_pop_par1[1], ci_poptaph_par1[1]),
           ci_par1_97_5 = c(ci_pop_par1[2], ci_poptaph_par1[2]),
           ml_par2 = c(ml_pars[2], ml_taph_pars[2]),
           ci_par2_2_5 = c(ci_pop_par2[1], ci_poptaph_par2[1]),
           ci_par2_97_5 = c(ci_pop_par2[2], ci_poptaph_par2[2]))

rownames(params) <- c("Logistic model", "Logistic model with taphonomy")
write.csv(params,
          here::here("analysis/data/derived_data/logistic_parameters.csv"))
