library(dplyr)
library(rcarbon)
library(ADMUR)
library(patchwork)
library(ggplot2)

c14 <- read.csv(here::here("analysis/data/raw_data/radiocarbon.csv"))
c14 <- c14 %>% filter(context != "Food crust")

caldates <- calibrate(x = c14$c14_bp, normalised = FALSE,
                      errors = c14$error, calCurves = "intcal20")
c14bins <- binPrep(sites =  c14$site_name, ages = c14$c14_bp, h = 200)

# SPD of entire range in use for fitting models
c14spd <- spd(caldates, bins = c14bins, timeRange = c(12000, 3500))

# Inspect the impact of binning
# binsense(x = caldates, y = c14$site_name, h = seq(0,500,100),
#         timeRange=c(12000,3700))

# Fit exponential model
expfit <- modelTest(x = caldates, errors = c14$error, bins = c14bins,
                    model = "exponential", timeRange = c(12000, 3500),
                    fitonly = TRUE)

# Fit logistic model
grd_c14 <- c14spd$grid %>% filter(PrDens > 0)
x <- grd_c14$calBP
y <- grd_c14$PrDens
logfit <- nls(y ~ SSlogis(x, Asym, xmid, scal),
              control = nls.control(warnOnly = TRUE, maxiter = 200),
              start = list(Asym = 0.2, xmid = 4500, scal = -100))
logdf <- data.frame(calBP = x, PrDens = predict(logfit), model = "Logistic")

# Fit uniform model
unifit <- modelTest(x = caldates, errors = c14$error, bins = c14bins,
                    model = "uniform", timeRange = c(12000, 3500),
                    fitonly = TRUE)

# Remove zero prob for visualisation to match shoreline date treatment
c14spd$grid <- c14spd$grid %>% filter(PrDens > 0)
expfit$fit <- expfit$fit %>% filter(calBP <= max(c14spd$grid$calBP))
expfit$fit$model <- "Exponential"
unifit$fit <- unifit$fit %>% filter(calBP <= max(c14spd$grid$calBP))
unifit$fit$model <- "Uniform"

models_fit <- rbind(logdf, expfit$fit, unifit$fit)

# timeRange argument in rcarbon::spd() appears to return sample sizes
# for the entire dataset, irrespective of the specified time range.
# Experimenting with the value for p found that the subset of dates where 0.0001
# of the cumulative prob falls within the employed date range finds the
# effective sample size.
subcal <- subset(caldates, BP <= 12000 & BP >= 3500, p = 0.000001)
c14sub <- c14[which.CalDates(caldates, BP <= 12000 & BP >= 3500, p = 0.000001),]
subbins = binPrep(sites =  c14sub$site_name, ages = c14sub$c14_bp, h = 200)

ggplot() +
  geom_bar(aes(x = (c14spd$grid$calBP-2000) * -1,
                y = c14spd$grid$PrDens), stat = "identity", col = "grey") +
  geom_line(data = models_fit, aes((calBP - 2000)*-1, PrDens, col = model),
            linewidth = 1.1) +
  geom_line(aes(x = (c14spd$grid$calBP-2000) * -1,
                y = c14spd$grid$PrDens)) +
  # geom_vline(xintercept = (3700-2000)*-1) +
  geom_text(aes(-8500, Inf, hjust = 0, vjust = 1.2,
                 label = paste0("Sites = ", length(unique(c14sub$site_name)),
                               "\nDates = ", length(subcal),
                               "\nBins = ", length(unique(subbins))))) +
  scale_x_continuous(breaks = seq(-9000, -2000, 1000),
                     expand = expansion(mult = c(0, 0), add = c(100, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0),  add = c(0, 0.005))) +
  labs(x = "BCE", y = "Summed probability",
     title = "Summed probability distribution of radiocarbon dates (unnormalised)") +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave(here::here("analysis/figures/rspd.png"),
       units = "px", width = 2250, height = 1750)

nsim <- 1000

subcal2 <- subset(caldates, BP <= 12000 & BP >= 4500, p = 0.000001)
c14sub2 <- c14[which.CalDates(caldates, BP <= 12000 & BP >= 4500, p = 0.000001),]
subbins2 = binPrep(sites =  c14sub2$site_name, ages = c14sub2$c14_bp, h = 200)

### Exponential
expnull <- modelTest(x = caldates, errors = c14$error,
                     timeRange = c(12000, 4500),
                     bins = c14bins, nsim = nsim, model = "exponential")

obs <- expnull$result[,1:2]
envelope <- expnull$result[,3:4]
booms <- which(obs$PrDens>envelope[,2])
busts <- which(obs$PrDens<envelope[,1])

expnull$result <- expnull$result %>% filter(PrDens > 0)
expnull$fit <- expnull$fit %>% filter(calBP < max(expnull$result$calBP))


plt4 <- ggplot(data = expnull$result, aes(x = (calBP-2000)*-1)) +
  geom_vline(xintercept = (obs$calBP[busts]-2000)*-1,
             col = "firebrick", alpha = 0.01) +
  geom_vline(xintercept = (obs$calBP[booms]-2000)*-1,
             col = "darkgreen", alpha = 0.01) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              fill = "grey60", alpha = 0.8) +
  geom_line(data = expnull$fit, aes(x = (calBP-2000)*-1, y = PrDens),
            linetype = "dashed", colour = "black") +
  labs(x = "BCE", y = "Summed probability") +
  scale_x_continuous(breaks = seq(-10000, 2000, 1000),
                     limits = c(-10000, -2500),
                     expand = expansion(mult = c(0, 0))) +
  # geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
  #               label = paste("p =", round(expnull$pval, 3)))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0),  add = c(0, 0.005))) +
  geom_line(aes(y = PrDens)) +
  theme_bw()

### Logistic
lognull <- modelTest(caldates, errors = c14$error, bins = c14bins, nsim = nsim,
                     timeRange = c(12000, 4500), predgrid = logdf,
                     model = "custom")

obs <- lognull$result[, 1:2]
envelope <- lognull$result[, 3:4]
booms <- which(obs$PrDens > envelope[, 2])
busts <- which(obs$PrDens < envelope[, 1])

lognull$result <- lognull$result %>% filter(PrDens > 0)
lognull$fit<- lognull$fit %>% filter(calBP < max(lognull$result$calBP))

plt5 <- ggplot(data = lognull$result, aes(x = (calBP - 2000)* -1)) +
  geom_vline(xintercept = (obs$calBP[busts]-2000)*-1,
             col = "firebrick", alpha = 0.01) +
  geom_vline(xintercept = (obs$calBP[booms]-2000)*-1,
             col = "darkgreen", alpha = 0.01) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              fill = "grey60", alpha = 0.8) +
  geom_line(data = lognull$fit, aes(x = (calBP-2000)*-1, y = PrDens),
            linetype = "dashed", colour = "black") +
  labs(x = "BCE", y = "Summed probability") +
  scale_x_continuous(breaks = seq(-10000, 2000, 1000),
                     limits = c(-10000, -2500),
                     expand = expansion(mult = c(0, 0))) +
  # geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
  #               label = paste("p =", round(lognull$pval, 3)))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0),  add = c(0, 0.005))) +
  geom_line(aes(y = PrDens)) +
  theme_bw()

### Uniform

uninull <- modelTest(x = caldates, errors = c14$error,
                     bins = c14bins, nsim = nsim,
                     timeRange = c(12000, 4500), model = "uniform")

obs <- uninull$result[,1:2]
envelope <- uninull$result[,3:4]
booms <- which(obs$PrDens>envelope[,2])
busts <- which(obs$PrDens<envelope[,1])

uninull$result <- uninull$result %>% filter(PrDens > 0)
uninull$fit <- uninull$fit %>% filter(calBP < max(uninull$result$calBP))

plt6 <- ggplot(data = uninull$result, aes(x = (calBP-2000)*-1)) +
  geom_vline(xintercept = (obs$calBP[busts]-2000)*-1,
             col = "firebrick", alpha = 0.01) +
  geom_vline(xintercept = (obs$calBP[booms]-2000)*-1,
             col = "darkgreen", alpha = 0.01) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              fill = "grey60", alpha = 0.8) +
  geom_line(data = uninull$fit, aes(x = (calBP-2000)*-1, y = PrDens),
            linetype = "dashed", colour = "black") +
  labs(x = "BCE", y = "Summed probability") +
  scale_x_continuous(breaks = seq(-10000, 2000, 1000),
                     limits = c(-10000, -2500),
                     expand = expansion(mult = c(0, 0))) +
  # geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
  #               label = paste("p =", round(expnull$pval, 3)))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0),  add = c(0, 0.005))) +
  geom_line(aes(y = PrDens)) +
  theme_bw()

(plt1 + plt2 + plt3) /
(plt4 + plt5 + plt6)

ggsave(here::here("analysis/figures/mc.png"),
       units = "px", width = 4000, height = 1250*2)

## Model comparison using ADMUR


# Using rcarbon, calibrate each date individually or sum dates that are binned
subcal2 <- subset(caldates, BP <= 12000 & BP >= 4500, p = 0.000001)
c14sub2 <- c14[which.CalDates(caldates, BP <= 12000 & BP >= 4500, p = 0.000001),]
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
  slist <- rapply(slist, f = function(x) ifelse(is.na(x), 0, x), how = "replace")
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

# # Rename model columns to match ADMUR syntax
# expmodel <- expnull$fit %>%
#   rename("year" = "calBP", "pdf" = "PrDens") %>%
#   filter(year < 10554) %>%
#   select(year, pdf)
#
# logmodel <- lognull$fit %>%
#   rename("year" = "calBP", "pdf" = "PrDens") %>%
#   select(year, pdf)
#
# unimodel <- uninull$fit %>%
#   rename("year" = "calBP", "pdf" = "PrDens") %>%
#   filter(year < 10554) %>%
#   select(year, pdf)

# Make year range in pd match that of the models
pd <- pd[as.numeric(rownames(pd)) <= max(expmodel$year),]
plotPD(pd)
# Check that SPD and models matches plt4-plt6 , above
SPD <- as.data.frame(rowSums(pd))
plotPD(SPD/(sum(SPD) * 5))

library(DEoptimR)
# Fit CPL-models to radiocarbon data
cpl_1 <- JDEoptim(lower = 0, upper = 1, fn = objectiveFunction, PDarray = pd,
                  type = 'CPL', NP = 20, trace = TRUE)
cpl_2 <- JDEoptim(lower = rep(0, 3), upper = rep(1, 3), fn = objectiveFunction,
                  PDarray = pd, type = 'CPL', NP = 60, trace = TRUE)
# Increased maxiter for the remaining searches to achieve convergence
cpl_3 <- JDEoptim(lower = rep(0, 5), upper = rep(1,5), fn = objectiveFunction,
                  PDarray = pd, type = 'CPL', NP = 100, maxiter = 400 * 5,
                  trace = TRUE)
cpl_4 <- JDEoptim(lower = rep(0, 7), upper = rep(1,7), fn = objectiveFunction,
                  PDarray = pd, type = 'CPL', NP = 140, maxiter = 400 * 7,
                  trace = TRUE)

# Exponential
exp <- JDEoptim(lower = -0.01, upper = 0.01, fn = objectiveFunction, PDarray = pd, type = 'exp', NP = 20)
# Uniform
uniform <- objectiveFunction(NULL, pd, type='uniform')

maxyr <- max(as.numeric(rownames(SPD)))
minyr <- min(as.numeric(rownames(SPD)))

cpl1 <- convertPars(pars = cpl_1$par, years = maxyr:minyr, type='CPL')
cpl1$model <- "CPL-1"
cpl2 <- convertPars(pars = cpl_2$par, years = -9445:-2500, type='CPL')
cpl2$model <- "CPL-2"
EXP <- convertPars(pars=exp$par, years = maxyr:minyr, type = 'exp')


plotPD(SPD)
lines(cpl1$year, cpl1$pdf, lwd = 2, col = cols[1])

ggplot() +
  geom_bar(aes(x = c14spd$grid$calBP,
               y = c14spd$grid$PrDens/sum(c14spd$grid$PrDens)), stat = "identity", col = "grey") +
  geom_line(aes(EXP$year, EXP$pdf))
  geom_hline(yintercept = uniform)

  geom_line(aes((cpl1$year-1950) *-1, cpl1$pdf))

cols <- c('steelblue','firebrick','orange')
plotPD(SPD)
lines(expmodel$year, expmodel$pdf, lwd = 2, col = cols[1])
lines(logmodel$year, logmodel$pdf, lwd = 2, col = cols[2])
lines(unimodel$year, unimodel$pdf, lwd = 2, col = cols[3])

# Find log-likelihood
explik <- loglik(PD = pd, model = expmodel)
loglik <- loglik(PD = pd, model = logmodel)
unilik <- loglik(PD = pd, model = unimodel)

# Find BIC
bic_exp <- 1*log(ncol(pd)) - 2*(explik)
bic_log <- 3*log(ncol(pd)) - 2*(loglik)
bic_uni <- 0 - 2*(unilik)

# Save results to be used in paper
c14_lik <- c("Exponential" = explik, "Logistic" = loglik, "Uniform" = unilik)
c14_bic <- c("Exponential" = bic_exp, "Logistic" = bic_log, "Uniform" = bic_uni)

c14_pvals <- c("Exponential" = expnull$pval, "Logistic" = lognull$pval,
               "Uniform" = uninull$pval)

c14_sample <- c("Sites" = length(unique(c14sub2$site_name)),
                "Dates" = length(subcal2), "Bins" = length(unique(subbins2)))



