library(ADMUR)
library(DEoptimR)

set.seed(42)

# Load data from radiocarbon_mc.R: models fit for MC, the pd and the min and
# max ages
load(file = here("analysis/data/derived_data/rcarbon_pd_models.RData"))

# Fit CPL-models to radiocarbon data
cpl_1 <- JDEoptim(lower = 0, upper = 1, fn = objectiveFunction, PDarray = pd,
                  type = 'CPL', NP = 20, trace = TRUE)
cpl_2 <- JDEoptim(lower = rep(0, 3), upper = rep(1, 3), fn = objectiveFunction,
                  PDarray = pd, type = 'CPL', NP = 60,  maxiter = 400 * 3,
                  trace = TRUE)
cpl_3 <- JDEoptim(lower = rep(0, 5), upper = rep(1, 5), fn = objectiveFunction,
                  PDarray = pd, type = 'CPL', NP = 100, maxiter = 400 * 5,
                  trace = TRUE)
cpl_4 <- JDEoptim(lower = rep(0, 7), upper = rep(1, 7), fn = objectiveFunction,
                  PDarray = pd, type = 'CPL', NP = 140, maxiter = 400 * 7,
                  trace = TRUE)
cpl_5 <- JDEoptim(lower = rep(0, 9), upper = rep(1, 9), fn = objectiveFunction,
                  PDarray = pd, type='CPL', NP = 180, maxiter = 400 * 9,
                  trace = TRUE)
cpl_6 <- JDEoptim(lower = rep(0, 11), upper = rep(1,11), fn = objectiveFunction,
                  PDarray = pd, type = 'CPL', NP = 220, maxiter =  400 * 11,
                  trace = TRUE)

save(exp, logi, cpl_1, cpl_2, cpl_3, cpl_4, cpl_5,
          file = here("analysis/data/derived_data/rcarbon_models.RData"))

load(file = here("analysis/data/derived_data/rcarbon_models.RData"))

expp <- convertPars(pars = exp$par, years = minage:maxage, type = 'exp')
expp$model <- "Exponential"
logip <- convertPars(pars = logi$par, years = minage:maxage, type = 'logistic')
logip$model <- "Logistic"
unifp <- convertPars(pars = NULL, years = minage:maxage, type = "uniform")
unifp$model <- "Uniform"
cpl1 <- convertPars(pars = cpl_1$par, years = minage:maxage, type = 'CPL')
cpl1$model <- "CPL-1"
cpl2 <- convertPars(pars = cpl_2$par, years = minage:maxage, type = 'CPL')
cpl2$model <- "CPL-2"
cpl3 <- convertPars(pars = cpl_3$par, years = minage:maxage, type = 'CPL')
cpl3$model <- "CPL-3"
cpl4 <- convertPars(pars = cpl_4$par, years = minage:maxage, type = 'CPL')
cpl4$model <- "CPL-4"
cpl5 <- convertPars(pars = cpl_5$par, years = minage:maxage, type = 'CPL')
cpl5$model <- "CPL-5"

rcarbonmodels <- rbind(expp, logip, unifp, cpl1, cpl2, cpl3, cpl4, cpl5)

rcarbon_lik <- c("Exponential" = -exp$value,
               "Logistic" = -logi$value,
               "Uniform" = unif,
               "CPL-1" = -cpl_1$value,
               "CPL-2" = -cpl_2$value,
               "CPL-3" = -cpl_3$value,
               "CPL-4" = -cpl_4$value,
               "CPL-5" = -cpl_5$value)

# "CPL-3" = -cpl_3$value,
# "CPL-4" = -cpl_4$value

# Effective sample size
ssize <- ncol(pd)

rcarbon_bic <- c(log(ssize)*1 - 2*rcarbon_lik[1], # Exponential
               log(ssize)*2 - 2*rcarbon_lik[2], # Logistic
               log(ssize)*0  - 2*rcarbon_lik[3], # Uniform
               log(ssize)*1 - 2*rcarbon_lik[4], # CPL-1
               log(ssize)*3 - 2*rcarbon_lik[5], # CPL-2
               log(ssize)*5 - 2*rcarbon_lik[6], # CPL-3
               log(ssize)*7 - 2*rcarbon_lik[7], # CPL-4
               log(ssize)*9 - 2*rcarbon_lik[8] # CPL-5
               ) #log(ssize)*11 - 2*rcarbon_lik[9] # CPL-6

# log(ssize)*5 - 2*shore_lik[6],
# log(ssize)*7 - 2*shore_lik[6],
# log(ssize)*9 - 2*shore_lik[7],
# log(ssize)*11 - 2*shore_lik[8],
# log(ssize)*13 - 2*shore_lik[9],
# log(ssize)*15 - 2*shore_lik[10]

SPD <- as.data.frame(rowSums(pd))

cplplt <- ggplot() +
  geom_bar(aes(x = (as.numeric(rownames(SPD)) - 2000) * -1,
               SPD[,1] / sum(SPD[,1])),
           stat = "identity", col = "grey") +
  geom_line(aes(x = (as.numeric(rownames(SPD)) - 2000) * -1,
                y = SPD[,1] / sum(SPD[,1]))) +
  geom_line(data = rcarbonmodels, aes((year - 2000) * -1, pdf, col = model),
            linewidth = 1.1) +
  labs(x = "BCE", y = "Summed probability") +
  # scale_x_continuous(limit = c(-10000, -2750)) +
  # scale_y_continuous(limit = c(0, 0.001)) +
  theme_bw() +
  theme(legend.title = element_blank())


bicplt <- ggplot() +
  geom_point(aes(rcarbon_bic, names(rcarbon_bic)), size = 3) +
  labs(x = "BIC", y = "Model") +
  theme_bw()

cplplt + bicplt
