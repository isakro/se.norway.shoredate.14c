library(ggplot2)
library(shoredate)
library(dplyr)
library(sf)
library(here)
library(patchwork)
library(ADMUR)
library(DEoptimR)

set.seed(42)

# Source functions
source(here("R/03_functions.R"))

# Load combined sites
sites <- st_read(here("analysis/data/derived_data/combined_sites.gpkg"))

# Load polygons with displacement curves
load(here("analysis/data/derived_data/displacement_polygons.RData"))

# Load shoreline dates
load(here("analysis/data/derived_data/sdates.RData"))

# Sum the probability distributions for the shoreline dates
sumdates <- sum_shoredates(shorelinedates)

# Make the results a data frame
sumdatesdf <- as.data.frame(sumdates) %>%
  filter(sum.probability != 0)

# Find density of sites across polygons
incpolys$dens <- lengths(st_intersects(incpolys, sites)) /
  sum(lengths(st_intersects(incpolys, sites)))

# Sample size
ssize = sumdates$dates_n

# Number of simulations
nsim <- 1000

# Below follow sections that fit a uniform and a exponential model to
# the shoredate SPD, and from this simulates a 95% critical envelope with which
# to compare the observed data.

# Fit exponential model (re-purposed from modelTest() from rcarbon)
fit <- nls(y ~ (exp(a + b * x)),
           data = data.frame(x = sumdatesdf$sum.bce,
                             y = sumdatesdf$sum.probability),
           start = list(a = 0, b = 0))
est <- predict(fit, list(x = sumdatesdf$sum.bce))
expgrid <- data.frame(bce = sumdatesdf$sum.bce, prob_dens = est)

sumdatesdf <- filter(sumdatesdf, sum.bce <= -2500)
expgrid <- filter(expgrid, bce <= -2500)

# Inspect SPD and fitted model
ggplot() +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y = sum.probability)) +
  geom_line(data = expgrid, aes(x = bce, y = prob_dens),
            linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000)) +
  theme_bw()

# Data frame for simulated dates
random_dates <- data.frame(matrix(nrow = ssize*nsim, ncol = 3))
names(random_dates) <- c("sample", "simn", "displacement_curve")


random_dates$sample = sample(expgrid$bce, replace = TRUE,
                             size = ssize*nsim, prob = expgrid$prob_dens)
random_dates$simn <- rep(1:nsim, each = nrow(random_dates)/nsim)
random_dates$displacement_curve <- incpolys[sample(incpolys$id, replace = TRUE,
                                              size = ssize*nsim,
                                              prob = incpolys$dens),]$disp[[1]]

# Reverse shoreline date using reverse_shoredate from 03_functions.R
random_dates$reverse_elevation <- mapply(reverse_shoredate,
                                         random_dates$sample,
                                         random_dates$displacement_curve,
                                         step_reso)

simdates <- vector("list", ssize)
start_time <- Sys.time()
for (i in 1:nrow(random_dates)){
  print(paste(i, "/", nrow(random_dates)))

  simdates[[i]] <- as.data.frame(shoredate::shoreline_date(
    site = random_dates$simn[i],
    elev_reso = step_reso,
    cal_reso = 5,
    elevation = random_dates$reverse_elevation[i],
    interpolated_curve = random_dates$displacement_curve[i],
    sparse = TRUE))
  simdates[[i]]$simn <- random_dates$simn[i]

  if(i %% ssize == 0){
    # Save data externally and overwrite simdates for memory purposes
      tmp <- bind_rows(simdates) %>% group_by(bce, simn) %>%
        filter(!is.na(probability)) %>%
        summarise(prob_sum = sum(probability))

      save(tmp,
           file = paste0(here::here("../external_data/shorespd/exp/"),
                         i, "shore.rds"))

      simdates <- vector("list", ssize)
  }
}
end_time <- Sys.time()
end_time - start_time


resfiles <- list.files(here::here("../external_data/shorespd/exp/"),
                       full.names = TRUE)

resfiles <- gtools::mixedsort(resfiles)

results <- list()
for(i in seq_along(resfiles)){
  load(resfiles[i])
  results[[i]] <- bind_rows(tmp)
  rm(tmp)
}

# Retrieve results and normlise each simulated SPD
simresults <- do.call(rbind.data.frame, results) %>%
  group_by(simn) %>%
  mutate(prob_sum_normalised = prob_sum / sum(prob_sum, na.rm = TRUE))

expsum <- simulation_summary(sumdates, simresults, cut_off = -2500, nsim)

simresults <- simresults %>%
  group_by(bce) %>%
  mutate(low = quantile(prob_sum_normalised, prob = 0.025, na.rm = TRUE),
         high = quantile(prob_sum_normalised, prob = 0.975, na.rm = TRUE),
         mean = mean(prob_sum_normalised)) %>%
  distinct(bce, .keep_all = TRUE)

plt1 <- ggplot() +
  geom_vline(xintercept = as.numeric(names(expsum$busts)),
             col = "firebrick", alpha = 0.04) +
  geom_vline(xintercept = as.numeric(names(expsum$booms)),
             col = "darkgreen", alpha = 0.04) +
  ggplot2::geom_ribbon(data = simresults, aes(x = bce, ymin = low, ymax = high),
                       fill = "grey60", alpha = 0.8) +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y = sum.probability)) +
  geom_line(data = expgrid, aes(x = bce, y = prob_dens),
            linetype = "dashed", colour = "black") +
  # geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
  #               label = paste("p =", round(expsum$pvalue, 3)))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0),  add = c(0, 0.0001))) +
  labs(title = "Exponential", x = "BCE", y = "Summed probability") +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000),
                     limits = c(-10000, -2500),
                     expand = expansion(mult = c(0, 0))) +
  theme_bw()

###### Logistic model ######
x = sumdatesdf$sum.bce
y = sumdatesdf$sum.probability

logfit <- nls(y ~ SSlogis(x, Asym, xmid, scal))

loggrid <- data.frame(bce = x, prob_dens = predict(logfit))
loggrid <- filter(loggrid, bce <= -2500)

# Inspect SPD and fitted model
ggplot() +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y = sum.probability)) +
  geom_line(data = loggrid, aes(x = bce, y = prob_dens),
            linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000)) +
  ggtitle("Logistic model") +
  theme_bw()

# Data frame for simulated dates
random_dates <- data.frame(matrix(nrow = ssize*nsim, ncol = 3))
names(random_dates) <- c("sample", "simn", "displacement_curve")

random_dates$sample = sample(loggrid$bce, replace = TRUE,
                             size = ssize*nsim, prob = loggrid$prob_dens)
random_dates$simn <- rep(1:nsim, each = nrow(random_dates)/nsim)
random_dates$displacement_curve <- incpolys[sample(incpolys$id, replace = TRUE,
                                                   size = ssize*nsim,
                                                   prob = incpolys$dens),]$disp[[1]]

# Reverse shoreline date using reverse_shoredate from 03_functions.R
random_dates$reverse_elevation <- mapply(reverse_shoredate,
                                         random_dates$sample,
                                         random_dates$displacement_curve,
                                         step_reso)

simdates <- vector("list", ssize)
start_time <- Sys.time()
for (i in 1:nrow(random_dates)){
  print(paste(i, "/", nrow(random_dates)))

  simdates[[i]] <- as.data.frame(shoredate::shoreline_date(
    site = random_dates$simn[i],
    elev_reso = step_reso,
    cal_reso = 5,
    elevation = random_dates$reverse_elevation[i],
    interpolated_curve = random_dates$displacement_curve[i],
    sparse = TRUE))
  simdates[[i]]$simn <- random_dates$simn[i]

  if(i %% ssize == 0){
    # Save data externally and overwrite simdates for memory purposes
    tmp <- bind_rows(simdates) %>% dplyr::group_by(bce, simn) %>%
      dplyr::filter(!is.na(probability)) %>%
      dplyr::summarise(prob_sum = sum(probability))

    save(tmp,
         file = paste0(here::here("../external_data/shorespd/logi/"),
                       i, "shore.rds"))

    simdates <- vector("list", ssize)
  }
}
end_time <- Sys.time()
end_time - start_time

resfiles <- list.files(here::here("../external_data/shorespd/logi/"),
                       full.names = TRUE)

resfiles <- gtools::mixedsort(resfiles)

results <- list()
for(i in seq_along(resfiles)){
  load(resfiles[i])
  results[[i]] <- bind_rows(tmp)
  rm(tmp)
}

simresults <- do.call(rbind.data.frame, results) %>%
  dplyr::group_by(simn) %>%
  dplyr::mutate(prob_sum_normalised = prob_sum / sum(prob_sum, na.rm = TRUE))

logisum <- simulation_summary(spd = sumdates,
                              simulation_results = simresults,
                              cut_off = -2500, nsim)

simresults <- simresults %>%
  dplyr::group_by(bce) %>%
  dplyr::mutate(low = quantile(prob_sum_normalised, prob = 0.025, na.rm = TRUE),
                high = quantile(prob_sum_normalised, prob = 0.975, na.rm = TRUE),
                mean = mean(prob_sum_normalised)) %>%
  dplyr::distinct(bce, .keep_all = TRUE)

plt2 <- ggplot() +
  geom_vline(xintercept = as.numeric(names(logisum$busts)),
             col = "firebrick", alpha = 0.04) +
  geom_vline(xintercept = as.numeric(names(logisum$booms)),
             col = "darkgreen", alpha = 0.04) +
  ggplot2::geom_ribbon(data = simresults, aes(x = bce, ymin = low, ymax = high),
                       fill = "grey60", alpha = 0.8) +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y = sum.probability)) +
  # geom_line(data = simresults, aes(x = bce, y = mean), col = "red", lwd = 0.5) +
  geom_line(data = loggrid, aes(x = bce, y = prob_dens),
            linetype = "dashed", colour = "black") +
  labs(title = "Logistic", x = "BCE", y = "Summed probability") +
  # geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
  #               label = paste("p =", round(logisum$pvalue, 3)))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0),  add = c(0, 0.0001))) +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000),
                     limits = c(-10000, -2500),
                     expand = expansion(mult = c(0, 0))) +
  theme_bw()


###### Uniform model ######

# Fit uniform model (re-purposed from modelTest() from rcarbon)
unigrid <- data.frame(bce = sumdatesdf$sum.bce,
                       prob_dens = mean(sumdatesdf$sum.probability))
unigrid <- filter(predgrid, bce <= -2500)

# Inspect SPD and fitted model
ggplot() +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y = sum.probability)) +
  geom_line(data = unigrid, aes(x = bce, y = prob_dens),
            linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000)) +
  theme_bw()

# Data frame for simulated dates
random_dates <- data.frame(matrix(nrow = ssize*nsim, ncol = 3))
names(random_dates) <- c("sample", "simn", "displacement_curve")


random_dates$sample = sample(unigrid$bce, replace = TRUE,
                             size = ssize*nsim, prob = unigrid$prob_dens)
random_dates$simn <- rep(1:nsim, each = nrow(random_dates)/nsim)
random_dates$displacement_curve <- incpolys[sample(incpolys$id, replace = TRUE,
                                            size = ssize*nsim,
                                            prob = incpolys$dens),]$disp[[1]]

# Reverse shoreline date using reverse_shoredate from 03_functions.R
random_dates$reverse_elevation <- mapply(reverse_shoredate,
                                         random_dates$sample,
                                         random_dates$displacement_curve,
                                         step_reso)

simdates <- vector("list", ssize)
start_time <- Sys.time()
for (i in 1:nrow(random_dates)){
  print(paste(i, "/", nrow(random_dates)))

  simdates[[i]] <- as.data.frame(shoredate::shoreline_date(
    site = random_dates$simn[i],
    elev_reso = step_reso,
    cal_reso = 5,
    elevation = random_dates$reverse_elevation[i],
    interpolated_curve = random_dates$displacement_curve[i],
    sparse = TRUE))
  simdates[[i]]$simn <- random_dates$simn[i]

  if(i %% ssize == 0){
    # Save data externally and overwrite simdates for memory purposes
    tmp <- bind_rows(simdates) %>% dplyr::group_by(bce, simn) %>%
      dplyr::filter(!is.na(probability)) %>%
      dplyr::summarise(prob_sum = sum(probability))

    save(tmp,
         file = paste0(here::here("../external_data/shorespd/uni/"),
                       i, "shore.rds"))

    simdates <- vector("list", ssize)
  }
}
end_time <- Sys.time()
end_time - start_time

resfiles <- list.files(here::here("../external_data/shorespd/uni/"),
                       full.names = TRUE)

resfiles <- gtools::mixedsort(resfiles)

results <- list()
for(i in seq_along(resfiles)){
  load(resfiles[i])
  results[[i]] <- bind_rows(tmp)
  rm(tmp)
}

simresults <- do.call(rbind.data.frame, results) %>%
  dplyr::group_by(simn) %>%
  dplyr::mutate(prob_sum_normalised = prob_sum / sum(prob_sum, na.rm = TRUE))

unisum <- simulation_summary(sumdates, simresults, cut_off = -2500, nsim)

simresults <- simresults %>%
  dplyr::group_by(bce) %>%
  dplyr::mutate(low = quantile(prob_sum_normalised, prob = 0.025, na.rm = TRUE),
         high = quantile(prob_sum_normalised, prob = 0.975, na.rm = TRUE),
         mean = mean(prob_sum_normalised)) %>%
  dplyr::distinct(bce, .keep_all = TRUE)

plt3 <- ggplot() +
  geom_vline(xintercept = as.numeric(names(unisum$busts)),
             col = "firebrick", alpha = 0.04) +
  geom_vline(xintercept = as.numeric(names(unisum$booms)),
             col = "darkgreen", alpha = 0.04) +
  ggplot2::geom_ribbon(data = simresults, aes(x = bce, ymin = low, ymax = high),
                       fill = "grey60", alpha = 0.8) +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y = sum.probability)) +
  geom_line(data = unigrid, aes(x = bce, y = prob_dens),
            linetype = "dashed", colour = "black") + #"grey40"
  labs(title = "Uniform", x = "BCE", y = "Summed probability") +
  # geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
  #               label = paste("p =", round(unisum$pvalue, 3)))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0),  add = c(0, 0.0001))) +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000),
                     limits = c(-10000, -2500),
                     expand = expansion(mult = c(0, 0))) +
  theme_bw()


plt1 + plt2 + plt3

ggsave(here::here("analysis/figures/shoreline_mc.png"),
       units = "px", width = 4000, height = 1250)

save(plt1, plt2, plt3,
     file = here("analysis/data/derived_data/shoredate_plots.RData"))

shore_pvals <- c("Exponential" = expsum$pvalue,
                 "Logistic" = logisum$pvalue,
                 "Uniform" = unisum$pvalue)

# Finding log likelihood using ADMUR

names(expgrid) <- c("year", "pdf")
expgrid <- purrr::map_df(expgrid, rev)
names(loggrid) <- c("year", "pdf")
loggrid <- purrr::map_df(loggrid, rev)
names(unigrid) <- c("year", "pdf")
unigrid <- purrr::map_df(unigrid, rev)

# Retrieve original dates before summing and make these a list of data frames
dates_dfs <- lapply(shorelinedates, as.data.frame)

# Combine these into a single data frame
sdates <- do.call(rbind, dates_dfs)

# Assign simulation results for each site a unique ID
sdates$site <- rep(1:length(dates_dfs) , each = nrow(sdates)/length(dates_dfs))

pd <- as.data.frame(sdates) %>%
  dplyr::filter(bce <= -2500 & bce >= min(expgrid$year)) %>%
  dplyr::group_by(bce) %>%
  dplyr::filter(!is.na(probability)) %>%
  dplyr::group_by(site) %>%
  tidyr::pivot_wider(names_from = site, values_from = probability) %>%
  tibble::column_to_rownames("bce")

pd <- sweep(pd, 2, colSums(pd),"/")

SPD <- as.data.frame(rowSums(pd))
SPD <- SPD/( sum(SPD) *5 )

# Check that this looks good with ADMUR
cols <- c('steelblue','firebrick','orange')
plotPD(SPD)
lines(expgrid$year, expgrid$pdf, lwd = 2, col = cols[1])
lines(loggrid$year, loggrid$pdf, lwd = 2, col = cols[2])
lines(unigrid$year, unigrid$pdf, lwd = 2, col = cols[3])
lines(expp$year, expp$pdf, lwd = 2)
lines(cpl1$year, cpl1$pdf, lwd = 2)
lines(cpl2$year, cpl2$pdf, lwd = 2)
lines(cpl3$year, cpl3$pdf, lwd = 2)

exploglik <- ADMUR::loglik(PD = pd, model = expgrid)
logloglik <- ADMUR::loglik(PD = pd, model = loggrid)
uniloglik <- ADMUR::loglik(PD = pd, model = unigrid)

exp <- JDEoptim(lower = -0.01, upper = 0.01, fn = objectiveFunction,
                PDarray = pd, type = "exp", NP = 20)
logi <- JDEoptim(lower = c(0, 0000), upper = c(1, 10000),
                 fn = objectiveFunction, PDarray = pd, type = 'logistic',
                 NP = 40, trace = TRUE)
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
cpl_5 <- JDEoptim(lower = rep(0, 9), upper = rep(1,9), fn = objectiveFunction,
                 PDarray = pd, type='CPL',
                 NP = 180, maxiter = 800 * 9,
                 trace = TRUE)
cpl_6 <- JDEoptim(lower = rep(0, 11), upper = rep(1,11), fn = objectiveFunction,
                  PDarray = pd, type = 'CPL', NP = 220, maxiter = 400 * 11,
                  trace = TRUE)
cpl_7 <- JDEoptim(lower = rep(0, 13), upper = rep(1, 13), fn = objectiveFunction,
                  PDarray = pd, type = 'CPL', NP = 260, maxiter = 400 * 13,
                  trace = TRUE)


save(cpl_1, cpl_2, cpl_3, cpl_4, cpl_5, cpl_6, exp, logi,
     file = here("analysis/data/derived_data/shoremodels.RData"))
load(here("analysis/data/derived_data/shoremodels.RData"))

expbic <- 2*log(ssize) - 2*(exploglik)
logbic <- 3*log(ssize) - 2*(logloglik)
unibic <- 0 - 2*(uniloglik)

shore_lik <- c("Exponential" = -exp$value,
               "Logistic" = -logi$value,
             "CPL-1" = -cpl_1$value,
             "CPL-2" = -cpl_2$value,
             "CPL-3" = -cpl_3$value,
             "CPL-4" = -cpl_4$value,
             "CPL-5" = -cpl_5$value,
             "CPL-6" = -cpl_6$value)

shore_bic <- c(log(ssize)*1 - 2*shore_lik[1],
               log(ssize)*2 - 2*shore_lik[2],
               log(ssize)*1 - 2*shore_lik[3],
               log(ssize)*3 - 2*shore_lik[4],
               log(ssize)*5 - 2*shore_lik[5],
               log(ssize)*7 - 2*shore_lik[6],
               log(ssize)*9 - 2*shore_lik[7],
               log(ssize)*11 - 2*shore_lik[8])

shore_aic <- c(2*1 - 2*shore_lik[1],
               2*2 - 2*shore_lik[2],
               2*3 - 2*shore_lik[4],
               2*5 - 2*shore_lik[5],
               2*7 - 2*shore_lik[6],
               2*9 - 2*shore_lik[7],
               2*11 - 2*shore_lik[8])

save(shore_lik,
     file = here("analysis/data/derived_data/shore_logliks.rds"))

expp <- convertPars(pars = exp$par, years = -9445:-2500, type = 'exp')
expp$model <- "Exponential"
logip <- convertPars(pars = logi$par, years = -9445:-2500, type = 'logistic')
logip$model <- "Logistic"
cpl1 <- convertPars(pars = cpl_1$par, years = -9445:-2500, type='CPL')
cpl1$model <- "CPL-1"
cpl2 <- convertPars(pars = cpl_2$par, years = -9445:-2500, type='CPL')
cpl2$model <- "CPL-2"
cpl3 <- convertPars(pars = cpl_3$par, years = -9445:-2500, type='CPL')
cpl3$model <- "CPL-3"
cpl4 <- convertPars(pars = cpl_4$par, years = -9445:-2500, type='CPL')
cpl4$model <- "CPL-4"
cpl5 <- convertPars(pars = cpl_5$par, years = -9445:-2500, type='CPL')
cpl5$model <- "CPL-5"
cpl6 <- convertPars(pars = cpl_6$par, years = -9445:-2500, type='CPL')
cpl6$model <- "CPL-6"
cpl7 <- convertPars(pars = cpl_7$par, years = -9445:-2500, type='CPL')
cpl7$model <- "CPL-7"

shoremodels <- rbind(expp, logip, cpl1, cpl2, cpl3, cpl4, cpl5, cpl6)

mplt <- ggplot() +
  geom_bar(aes(x = as.numeric(rownames(SPD)), SPD[,1]),
           stat = "identity", col = "grey") +
  geom_line(aes(x = as.numeric(rownames(SPD)),
                y = SPD[,1])) +
  geom_line(data = shoremodels, aes(year, pdf, col = model),
            linewidth = 1.1) +
  # geom_line(data = cpl6, aes(year, pdf),
  #           linewidth = 1.1) +
  # geom_vline(xintercept = (3700-2000)*-1) +
  # scale_x_continuous(breaks = seq(-9000, -2000, 1000),
  #                    expand = expansion(mult = c(0, 0), add = c(100, 0))) +
  # scale_y_continuous(expand = expansion(mult = c(0, 0),  add = c(0, 0.005))) +
  labs(x = "BCE", y = "Summed probability (unnormalised)") +
  theme_bw() +
  theme(legend.title = element_blank())

bicplt <- ggplot() +
  geom_point(aes(names(shore_bic), shore_bic), size = 3) +
  labs(x = "Model", y = "BIC") +
  theme_bw()

aicplt <- ggplot() +
  geom_point(aes(names(shore_aic), shore_aic), size = 3) +
  labs(x = "Model", y = "BIC") +
  theme_bw()

mplt + bicplt + aicplt
