library(ADMUR)

# Load pd, logistic model, and in and max age for RSPD
load(file = here::here("analysis/data/derived_data/rcarbon_mcmc_data.RData"))

# Rename
minage_rcarbon <- minage
maxage_rcarbon <- maxage
pd_rcarbon <- PD
log_rcarbon <- log

# Load pd, 4-CPL model, and in and max age for SSPD
load(file = here::here("analysis/data/derived_data/shore_mcmc_data.RData"))

# Rename
minage_shore <- minage
maxage_shore <- maxage
pd_shore <- pd
cpl4_shore <- cpl_4

chain_rcarbon <- mcmc(PDarray = pd_rcarbon,
                      startPars = log_rcarbon$par,
                      type = 'logistic',
                      N = 30000,
                      burn = 2000,
                      thin = 5,
                      jumps = 0.025)

chain_shore <- mcmc(PDarray = pd_shore,
                      startPars = cpl4_shore$par,
                      type = 'CPL',
                      N = 30000,
                      burn = 2000,
                      thin = 5,
                      jumps = 0.025)

save(chain_rcarbon, chain_shore,
     file = here("analysis/data/derived_data/mcmc_chains.RData"))

chain_rcarbon$acceptance.ratio
chain_shore$acceptance.ratio

chain_rcarbon$all.pars

