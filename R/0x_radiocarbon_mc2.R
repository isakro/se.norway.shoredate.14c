library(dplyr)
library(rcarbon)
library(ADMUR)
library(patchwork)
library(ggplot2)
library(DEoptimR)

set.seed(42)

c14 <- read.csv(here::here("analysis/data/raw_data/radiocarbon.csv"))
c14 <- c14 %>% filter(context != "Food crust")

caldates <- calibrate(x = c14$c14_bp, normalised = FALSE,
                      errors = c14$error, calCurves = "intcal20")
c14bins <- binPrep(sites =  c14$site_name, ages = c14$c14_bp, h = 200)

# SPD of entire range in use for fitting models
c14spd <- spd(caldates, bins = c14bins, timeRange = c(12000, 3500))


# Using rcarbon, calibrate each date individually or sum dates that are binned
subcal2 <- subset(caldates, BP <= 12000 & BP >= 4500, p = 0.05)
c14sub2 <- c14[which.CalDates(caldates, BP <= 12000 &
                                BP >= 4500, p = 0.05),]
subbins2 = binPrep(sites =  c14sub2$site_name, ages = c14sub2$c14_bp, h = 200)

# Retrieve binned and unbinned dates (based on code from rcarbon::spd())
x <- subcal2
binNames <- unique(subbins2)

calyears <- data.frame(calBP = seq(12000, 4500, -1))
binnedMatrix <- matrix(NA, nrow= nrow(calyears), ncol = length(binNames))
colnames(binnedMatrix) <- binNames
rownames(binnedMatrix) <- rev(calyears[,1])

for(b in 1:length(binNames)){
  index <- which(subbins2 == binNames[b])
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

  # tmp <- lapply(slist,`[`,2)
  # spdtmp <- Reduce("+", tmp)
  binnedMatrix[,b] <- spdtmp[,1]/length(index)
}

pd <- as.data.frame(binnedMatrix)

# Sanity check: Compare rcarbon SPD and SPD adjusted for ADMUR
c14spd <- rcarbon::spd(subcal2, bins = subbins2, timeRange = c(12000, 4500))
SPD <- as.data.frame(rowSums(pd))

par(mfrow = c(2,1))
plot(c14spd)
plotPD(SPD)

maxage <- max(c14spd$grid[c14spd$grid$PrDens > 0,]$calBP)
minage <- min(c14spd$grid[c14spd$grid$PrDens > 0,]$calBP)

# Maximum likelihood search using DEoptimR
exp <- JDEoptim(lower = -0.01, upper = 0.01, fn = objectiveFunction,
                PDarray = pd, type = "exp", NP = 20)
logi <- JDEoptim(lower = c(0, 0000), upper = c(1, 10000),
                 fn = objectiveFunction, PDarray = pd, type = 'logistic',
                 NP = 40, trace = TRUE)
# No search required for uniform
unif <- -objectiveFunction(pars = NULL, PDarray = pd, type = 'uniform')

expsum <- SPDsimulationTest(c14, calcurve = intcal20,
                            calrange = c(minage, maxage), pars = exp$par,
                            type = 'exp', N = 10000)
logsum <- SPDsimulationTest(c14, calcurve = intcal20,
                            calrange = c(minage, maxage), pars = log$par,
                            type = 'logistic', N = 10000)
unisum <- SPDsimulationTest(c14, calcurve = intcal20,
                            calrange = c(minage, maxage), pars = NULL,
                            type = 'uniform', N = 10000)
plotSimulationSummary(expsum)
plotSimulationSummary(logsum)
plotSimulationSummary(unisum)


expp <- convertPars(pars = exp$par, years = minage:maxage, type = 'exp')
expp$model <- "Exponential"
logip <- convertPars(pars = logi$par, years = minage:maxage, type = 'logistic')
logip$model <- "Logistic"
unifp <- convertPars(pars = NULL, years = minage:maxage, type = "uniform")
unifp$model <- "Uniform"

save(pd, exp, logi, unif, maxage, minage,
     file = here("analysis/data/derived_data/rcarbon_pd_models.RData"))

rspdmodels <- rbind(expp, logip, unifp)

#  Inspect models
ggplot() +
  geom_bar(aes(x = c14spd$grid$calBP,
               y = c14spd$grid$PrDens / sum(c14spd$grid$PrDens)),
           stat = "identity", col = "grey") +
  geom_line(data = rspdmodels, aes(year, pdf, col = model),
            linewidth = 1.1) +
  scale_x_reverse()

# Rename for rcarbon
names(expp) <- c("calBP", "PrDens", "model")
expp <- expp[,1:2]
names(logip) <- c("calBP", "PrDens", "model")
logip <- logip[,1:2]
names(unifp) <- c("calBP", "PrDens", "model")
unifp  <- unifp[,1:2]

nsim <- 10000

# Monte carlo simulation using rcarbon
expnull <- modelTest(x = subcal2, errors = c14$error, bins = subbins2,
                    model = "custom", predgrid = expp, nsim = nsim,
                    spdnormalised = TRUE,
                    timeRange = c(12000, 3500))
lognull <- modelTest(x = subcal2, errors = c14$error, bins = subbins2,
                    model = "custom", predgrid = logip, nsim = nsim,
                    spdnormalised = TRUE,
                    timeRange = c(12000, 3500))
uninull <- modelTest(x = subcal2, errors = c14$error, bins = subbins2,
                    model = "custom", predgrid = unifp, nsim = nsim,
                    spdnormalised = TRUE,
                    timeRange = c(12000, 3500))

save(expnull, lognull, uninull,
     file = here("analysis/data/derived_data/rcarbon_model_test.RData"))
load(file = here("analysis/data/derived_data/rcarbon_model_test.RData"))


# Plot exponential
obs <- expnull$result[,1:2]
envelope <- expnull$result[,3:4]
booms <- which(obs$PrDens>envelope[,2])
busts <- which(obs$PrDens<envelope[,1])

expnull$result <- expnull$result %>% filter(PrDens > 0)
expnull$fit <- expnull$fit %>% filter(calBP < max(expnull$result$calBP))

exprplt <- ggplot(data = expnull$result, aes(x = (calBP-2000)*-1)) +
  geom_vline(xintercept = (obs$calBP[busts]-2000)*-1,
             col = "firebrick", alpha = 0.01) +
  geom_vline(xintercept = (obs$calBP[booms]-2000)*-1,
             col = "darkgreen", alpha = 0.01) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              fill = "grey60", alpha = 0.8) +
  # geom_line(data = expnull$fit, aes(x = (calBP-2000)*-1, y = PrDens),
  #           linetype = "dashed", colour = "black") +
  geom_line(aes(y = PrDens)) +
  geom_line(data = expp, aes((calBP-2000)*-1, PrDens),
            linewidth = 0.5, col = "red") +
  labs(x = "BCE", y = "Summed probability") +
  scale_x_continuous(breaks = seq(-10000, 2000, 1000),
                     limits = c(-10000, -2500),
                     expand = expansion(mult = c(0, 0))) +
  # geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
  #               label = paste("p =", round(expnull$pval, 3)))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme_bw()

# Plot logistic
obs <- lognull$result[,1:2]
envelope <- lognull$result[,3:4]
booms <- which(obs$PrDens>envelope[,2])
busts <- which(obs$PrDens<envelope[,1])

lognull$result <- lognull$result %>% filter(PrDens > 0)
lognull$fit <- lognull$fit %>% filter(calBP < max(lognull$result$calBP))


logrplt <- ggplot(data = lognull$result, aes(x = (calBP-2000)*-1)) +
  geom_vline(xintercept = (obs$calBP[busts]-2000)*-1,
             col = "firebrick", alpha = 0.01) +
  geom_vline(xintercept = (obs$calBP[booms]-2000)*-1,
             col = "darkgreen", alpha = 0.01) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              fill = "grey60", alpha = 0.8) +
  # geom_line(data = lognull$fit, aes(x = (calBP-2000)*-1, y = PrDens),
  #           linetype = "dashed", colour = "black") +
  geom_line(aes(y = PrDens)) +
  geom_line(data = logip, aes((calBP-2000)*-1, PrDens),
            linewidth = 0.5, col = "red") +
  labs(x = "BCE", y = "Summed probability") +
  scale_x_continuous(breaks = seq(-10000, 2000, 1000),
                     limits = c(-10000, -2500),
                     expand = expansion(mult = c(0, 0))) +
  # geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
  #               label = paste("p =", round(lognull$pval, 3)))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme_bw()

# Plot uniform
obs <- uninull$result[,1:2]
envelope <- uninull$result[,3:4]
booms <- which(obs$PrDens>envelope[,2])
busts <- which(obs$PrDens<envelope[,1])

uninull$result <- uninull$result %>% filter(PrDens > 0)
uninull$fit <- uninull$fit %>% filter(calBP < max(uninull$result$calBP))

unirplt <- ggplot(data = uninull$result, aes(x = (calBP-2000)*-1)) +
  geom_vline(xintercept = (obs$calBP[busts]-2000)*-1,
             col = "firebrick", alpha = 0.01) +
  geom_vline(xintercept = (obs$calBP[booms]-2000)*-1,
             col = "darkgreen", alpha = 0.01) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              fill = "grey60", alpha = 0.8) +
  # geom_line(data = uninull$fit, aes(x = (calBP-2000)*-1, y = PrDens),
  #           linetype = "dashed", colour = "black") +
  geom_line(aes(y = PrDens)) +
  geom_line(data = unifp, aes((calBP-2000)*-1, PrDens),
            linewidth = 0.5, col = "red") +
  labs(x = "BCE", y = "Summed probability") +
  scale_x_continuous(breaks = seq(-10000, 2000, 1000),
                     limits = c(-10000, -2500),
                     expand = expansion(mult = c(0, 0))) +
  # geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
  #               label = paste("p =", round(uninull$pval, 3)))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme_bw()

if(expnull$pval < 0.0001) {
  exprplt <- exprplt +
    geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
                  label = paste("p < 0.0001")))
} else {
  exprplt <- exprplt +
    geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
                  label = paste("p =", round(expnull$pval, 3))))
}

if(lognull$pval < 0.0001) {
  logrplt <- logrplt +
    geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
                  label = paste("p < 0.0001")))
} else {
  logrplt <- logrplt +
    geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
                  label = paste("p =", round(lognull$pval, 3))))
}

if(uninull$pval < 0.0001) {
  unirplt <- unirplt +
    geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
                  label = paste("p < 0.0001")))
} else {
  unirplt <- unirplt +
    geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
                  label = paste("p =", round(uninull$pval, 3))))
}

exprplt + logrplt + unirplt

ggsave(here::here("analysis/figures/rcarbon_mc.png"),
       units = "px", width = 4000, height = 1250)

