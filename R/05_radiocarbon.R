library(dplyr)
library(rcarbon)


c14 <- read.csv(here::here("analysis/data/raw_data/radiocarbon.csv"))
c14 <- c14 %>% filter(context != "Food crust")

caldates <- calibrate(x = c14$c14_bp, normalised = FALSE,
                      errors = c14$error, calCurves = "intcal20")

c14bins = binPrep(sites = c14$site_name, ages = c14$c14_bp, h = 200)
c14spd <- spd(caldates, bins = c14bins, timeRange = c(11000, 2000))

nsim <- 10

### Exponential
expnull <- modelTest(x = caldates, errors = c14$error, bins = c14bins, nsim = nsim,
                     timeRange = c(12000, 4500), model = "exponential")

plot(expnull, calendar = "BCAD")

ggplot(data = expnull$result, aes(x = (calBP-2000)*-1)) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              fill = "grey", alpha = 0.9) +
  geom_line(data = expnull$fit, aes(x = (calBP-2000)*-1, y = PrDens),
            linetype = "dashed", colour = "grey40") +
  labs(x = "BCE", y = "Summed probability") +
  scale_x_continuous(breaks = seq(-10000, 2000, 1000),
                     limits = c(-10000, -2500)) +
  geom_line(aes(y = PrDens)) +
  theme_bw()

### Logistic

grd_c14 <- c14spd$grid %>% filter(PrDens > 0)
x <- grd_c14$calBP
y <- grd_c14$PrDens
logfit <- nls(y ~ SSlogis(x, Asym, xmid, scal),
             control = nls.control(warnOnly = TRUE, maxiter = 200),
             start = list(Asym = 0.2, xmid = 4500, scal = -100))
logdf <- data.frame(calBP = x, PrDens = predict(logfit))

plot(grd_c14$calBP, grd_c14$PrDens, type = "l")
lines(grd_c14$calBP, predict(logfit))

lognull <- modelTest(caldates, errors = c14$error, bins = c14bins, nsim = nsim,
                     timeRange = c(12000, 4500), predgrid = logdf,
                     model = "custom", runm = 100)

plot(lognull, calendar = "BCAD")

ggplot(data = lognull$result, aes(x = (calBP-2000)*-1)) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              fill = "grey", alpha = 0.9) +
  geom_line(data = logdf, aes(x = (calBP-2000)*-1, y = PrDens),
            linetype = "dashed", colour = "grey40") +
  labs(x = "BCE", y = "Summed probability") +
  scale_x_continuous(breaks = seq(-10000, 2000, 1000),
                     limits = c(-10000, -2500)) +
  scale_y_continuous(limits = c(0, 0.1)) +
  geom_line(aes(y = PrDens)) +
  theme_bw()

### Uniform

library(rcarbon)
data(euroval)

set.seed(1)

dk <- subset(euroevol, Country == "Denmark")
dk_cal <- calibrate(x = dk$C14Age, errors = dk$C14SD, calCurves = "intcal20")
dk_bins <- binPrep(sites = dk$SiteID, ages = dk$C14Age, h = 100)


swe <- subset(euroevol, Country == "Sweden")
swe_cal <- calibrate(x = swe$C14Age, errors = swe$C14SD, calCurves = "intcal20")
swe_bins <- binPrep(sites = swe$SiteID, ages = swe$C14Age, h = 100)


dk_exp <- modelTest(dk_cal, errors = dk$C14SD, bins = dk_bins, nsim = 100,
                    timeRange = c(8000, 4000), model = "exponential")

dk_uni <- modelTest(dk_cal, errors = dk$C14SD, bins = dk_bins, nsim = 100,
                    timeRange = c(8000, 4000), model = "uniform")

swe_exp <- modelTest(swe_cal, errors = swe$C14SD, bins = swe_bins, nsim = 100,
                     timeRange = c(8000, 4000), model = "exponential")

swe_uni <- modelTest(swe_cal, errors = swe$C14SD, bins = swe_bins, nsim = 100,
                     timeRange = c(8000, 4000), model = "uniform")

unique(c(dk_exp$pval, dk_uni$pval, swe_exp$pval, swe_uni$pval))

plot(dk_exp)
plot(dk_uni)
plot(swe_exp)
plot(swe_uni)
