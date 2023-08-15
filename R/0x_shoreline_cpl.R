library(ADMUR)
library(DEoptimR)
library(ggplot2)
library(patchwork)

set.seed(42)

# Load data from shoreline_monte_carlo.R: models fit for MC, the pd and the
# min and max ages
load(file = here::here("analysis/data/derived_data/shore_models_pd.RData"))

# Fit CPL-models to shoreline dates
cpl_1 <- JDEoptim(lower = 0, upper = 1, fn = objectiveFunction, PDarray = pd,
                  type = 'CPL', NP = 20, trace = TRUE)
cpl_2 <- JDEoptim(lower = rep(0, 3), upper = rep(1, 3), fn = objectiveFunction,
                  PDarray = pd, type = 'CPL', NP = 60, trace = TRUE)
cpl_3 <- JDEoptim(lower = rep(0, 5), upper = rep(1, 5), fn = objectiveFunction,
                  PDarray = pd, type = 'CPL', NP = 100, trace = TRUE)
cpl_4 <- JDEoptim(lower = rep(0, 7), upper = rep(1, 7), fn = objectiveFunction,
                  PDarray = pd, type = 'CPL', NP = 140, trace = TRUE)
cpl_5 <- JDEoptim(lower = rep(0, 9), upper = rep(1, 9), fn = objectiveFunction,
                  PDarray = pd, type='CPL', NP = 180, trace = TRUE)
cpl_6 <- JDEoptim(lower = rep(0, 11), upper = rep(1,11), fn = objectiveFunction,
                  PDarray = pd, type = 'CPL', NP = 220, trace = TRUE)

save(cpl_1, cpl_2, cpl_3, cpl_4, cpl_5, cpl_6,
     file = here::here("analysis/data/derived_data/shore_models.RData"))
load(file = here::here("analysis/data/derived_data/shore_models.RData"))

cpl1 <- convertPars(pars = cpl_1$par, years = minage:maxage, type = 'CPL')
cpl1$model <- "1-CPL"
cpl2 <- convertPars(pars = cpl_2$par, years = minage:maxage, type = 'CPL')
cpl2$model <- "2-CPL"
cpl3 <- convertPars(pars = cpl_3$par, years = minage:maxage, type = 'CPL')
cpl3$model <- "3-CPL"
cpl4 <- convertPars(pars = cpl_4$par, years = minage:maxage, type = 'CPL')
cpl4$model <- "4-CPL"

shmodels <- rbind(expp, logip, unifp, cpl1, cpl2, cpl3, cpl4, cpl5, cpl6)

sh_lik <- c("Exponential" = -exp$value,
                 "Logistic" = -logi$value,
                 "Uniform" = unif,
                 "1-CPL" = -cpl_1$value,
                 "2-CPL" = -cpl_2$value,
                 "3-CPL" = -cpl_3$value,
                 "4-CPL" = -cpl_4$value,
                 "5-CPL" = -cpl_5$value,
                 "6-CPL" = -cpl_6$value)

ssize <- ncol(pd)
sh_bic <- c(log(ssize)*1 - 2*sh_lik[1], # Exponential
                 log(ssize)*2 - 2*sh_lik[2], # Logistic
                 log(ssize)*0  - 2*sh_lik[3], # Uniform
                 log(ssize)*1 - 2*sh_lik[4], # 1-CPL
                 log(ssize)*3 - 2*sh_lik[5], # 2-CPL
                 log(ssize)*5 - 2*sh_lik[6], # 3-CPL
                 log(ssize)*7 - 2*sh_lik[7], # 4-CPL
                 log(ssize)*9 - 2*sh_lik[8], # 5-CPL
                 log(ssize)*11 - 2*sh_lik[9]) # 6-CPL

SPD <- as.data.frame(rowSums(pd))
SPD <- SPD/(sum(SPD) * 5)

cplplt <- ggplot() +
  geom_bar(aes(x = as.numeric(rownames(SPD)),
               SPD[,1]),
           stat = "identity", col = "grey") +
  geom_line(aes(x = as.numeric(rownames(SPD)),
                y = SPD[,1])) +
  geom_line(data = shmodels, aes(year, pdf, col = model),
            linewidth = 1.1) +
  labs(x = "BCE", y = "Summed probability") +
  scale_x_reverse(labels = function(x)(x-1950)*-1) +
  theme_bw() +
  theme(legend.title = element_blank())

bicplt_shore <- ggplot() +
  geom_point(aes(sh_bic, names(sh_bic)), size = 2) +
  labs(x = "BIC", y = "", subtitle = "Models fit to SSPD") +
  theme_bw() +
  theme(panel.border = element_rect(linewidth = 0.5))

cplplt + bicplt_shore

# Save pd, 4-CPL model and min and max age for mcmc
save(cpl_4, pd, minage, maxage,
     file = here("analysis/data/derived_data/shore_mcmc_data.RData"))
