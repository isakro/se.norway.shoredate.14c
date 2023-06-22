library(ggplot2)
library(shoredate)
library(dplyr)
library(sf)
library(here)
library(patchwork)
library(ADMUR)
library(DEoptimR)

set.seed(1)

# Source functions
source(here("R/03_functions.R"))

# Load combined sites
sites <- st_read(here("analysis/data/derived_data/combined_sites.gpkg"))

# Load polygons with displacement curves
load(here("analysis/data/derived_data/displacement_polygons.RData"))

# Load shoreline dates
load(here("analysis/data/derived_data/sdates.RData"))

dates_dfs <- lapply(shorelinedates, as.data.frame)

# Combine these into a single data frame
sdates <- do.call(rbind, dates_dfs)

# Assign results for each site a unique ID
sdates$site <- rep(1:length(dates_dfs) , each = nrow(sdates)/length(dates_dfs))

# Format for ADMUR, changing BCE to BP
pd <- as.data.frame(sdates) %>%
  dplyr::filter(bce <= -2500 & bce >= -9445) %>%
  dplyr::mutate(bp = (bce * -1) + 1950) %>%
  dplyr::select(-bce) %>%
  dplyr::group_by(bp) %>%
  dplyr::filter(!is.na(probability)) %>%
  dplyr::group_by(site) %>%
  tidyr::pivot_wider(names_from = site, values_from = probability) %>%
  tibble::column_to_rownames("bp")
pd <- sweep(pd, 2, colSums(pd),"/")

# Create SPD using the ADMUR approach
SPD <- as.data.frame(rowSums(pd))
SPD <- SPD/( sum(SPD) *5 )

# Find mininum and maximum ages (BP)
minage <- min(as.numeric(row.names(pd)))
maxage <- max(as.numeric(row.names(pd)))

# Maximum likelihood search using DEoptimR
exp <- JDEoptim(lower = -0.01, upper = 0.01, fn = objectiveFunction,
                PDarray = pd, type = "exp", NP = 20)
# No search required for uniform
unif <- -objectiveFunction(pars = NULL, PDarray = pd, type = 'uniform')

expp <- convertPars(pars = exp$par, years = minage:maxage, type = 'exp')
expp$model <- "Exponential"
unifp <- convertPars(pars = NULL, years =  minage:maxage, type = "uniform")
unifp$model <- "Uniform"

# The logistic model JDEoptim search appears to only succeed with another
# transformation of the BCE dates
pd2 <- as.data.frame(sdates) %>%
  dplyr::filter(bce <= -2500 & bce >= -9445) %>%
  dplyr::mutate(bp = (bce + 1950) * -1) %>%
  dplyr::select(-bce) %>%
  dplyr::group_by(bp) %>%
  dplyr::filter(!is.na(probability)) %>%
  dplyr::group_by(site) %>%
  tidyr::pivot_wider(names_from = site, values_from = probability) %>%
  tibble::column_to_rownames("bp")
pd2 <- sweep(pd2, 2, colSums(pd2),"/")
minage2 <- min(as.numeric(row.names(pd2)))
maxage2 <- max(as.numeric(row.names(pd2)))

logi <- JDEoptim(lower = c(0, 0000), upper = c(1, 10000),
                 fn = objectiveFunction, PDarray = pd2, type = 'logistic',
                 NP = 40, trace = TRUE)

logip <- convertPars(pars = logi$par, years = minage2:maxage2, type = 'logistic')
logip$model <- "Logistic"
logip$year <- logip$year + 3900 # Transforming back to BP


# Save
save(pd, expp, logip, unifp, maxage, minage,
     file = here("analysis/data/derived_data/shore_models_pd.RData"))
load(here("analysis/data/derived_data/shore_models_pd.RData"))

# Combine models for plotting
shoremodels <- rbind(expp, logip, unifp)

# Call to plot
ggplot() +
  geom_bar(aes(x = as.numeric(rownames(SPD)), SPD[,1]),
           stat = "identity", col = "grey") +
  geom_line(aes(x = as.numeric(rownames(SPD)),
                y = SPD[,1])) +
  geom_line(data = shoremodels, aes(year, pdf, col = model),
            linewidth = 1.1) +
  scale_x_reverse()

# Change from BP to BCE to use with shoredate

expp$year <- (expp$year - 1950) * -1
logip$year <- (logip$year - 1950) * -1
unifp$year <- (unifp$year - 1950) * -1

#### Set up Monte Carlo simulation ####

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
nsim <- 10000

names(expp) <- c("bce", "prob_dens", "model")

# Restrict SSPD and normalise probability
sumdatesdf <- filter(sumdatesdf, sum.bce <= max(expp$bce))
sumdatesdf$probability <- sumdatesdf$sum.probability /
                          (sum(sumdatesdf$sum.probability) * 5)

##### Monte Carlo simulation - Exponential function #####

# Create data frame to hold simulated dates
random_dates <- data.frame(matrix(nrow = ssize*nsim, ncol = 3))
names(random_dates) <- c("sample", "simn", "displacement_curve")

random_dates$sample = sample(expp$bce, replace = TRUE,
                             size = ssize*nsim, prob = expp$prob_dens)
random_dates$simn <- rep(1:nsim, each = nrow(random_dates)/nsim)
random_dates$displacement_curve <- incpolys[sample(incpolys$id,
                                              replace = TRUE,
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

results <- list()
for(i in 1:1000){
  load(resfiles[i])
  results[[i]] <- bind_rows(tmp)
  rm(tmp)
}

# Retrieve results and normalise each simulated SPD
simresults <- do.call(rbind.data.frame, results) %>%
  group_by(simn) %>%
  mutate(prob_sum_normalised = prob_sum /(sum(prob_sum, na.rm = TRUE) * 5))

# Model summary (can also be plotted with ADMUR)
expp <- convertPars(pars = exp$par,
                    years = as.numeric(rownames(SPD)), type = 'exp')
mod <- approx(x = expp$year, y = expp$pdf, xout = as.numeric(rownames(SPD)),
              ties = 'ordered', rule = 2)$y
exp_summary <- simulation_summary(SPD, simulation_results,
                                  c(-2500, -9445), mod, ncol(pd))

# ADMUR plot
plotSimulationSummary(exp_summary)

# Custom plot
plot_mc(exp_summary)

save(expplt, file = "../external_data/shorespd/expplt.rda")

##### Monte Carlo simulation - Logistic #####

names(logip) <- c("bce", "prob_dens", "model")

# Data frame for simulated dates
random_dates <- data.frame(matrix(nrow = ssize*nsim, ncol = 3))
names(random_dates) <- c("sample", "simn", "displacement_curve")

random_dates$sample = sample(logip$bce, replace = TRUE,
                             size = ssize*nsim, prob = logip$prob_dens)
random_dates$simn <- rep(1:nsim, each = nrow(random_dates)/nsim)
random_dates$displacement_curve <- incpolys[sample(incpolys$id,
                                              replace = TRUE,
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
  dplyr::mutate(
    prob_sum_normalised = prob_sum / (sum(prob_sum, na.rm = TRUE) * 5))

logisum <- simulation_summary(spd = sumdates,
                              simulation_results = simresults,
                              cut_off = -2500, nsim)

simresults <- simresults %>%
  dplyr::group_by(bce) %>%
  dplyr::mutate(low = quantile(prob_sum_normalised, prob = 0.025,
                               na.rm = TRUE),
                high = quantile(prob_sum_normalised, prob = 0.975,
                                na.rm = TRUE),
                mean = mean(prob_sum_normalised)) %>%
  dplyr::distinct(bce, .keep_all = TRUE)

logiplt <- ggplot() +
  geom_vline(xintercept = as.numeric(names(logisum$busts)),
             col = "firebrick", alpha = 0.04) +
  geom_vline(xintercept = as.numeric(names(logisum$booms)),
             col = "darkgreen", alpha = 0.04) +
  ggplot2::geom_ribbon(data = simresults, aes(x = bce, ymin = low, ymax = high),
                       fill = "grey60", alpha = 0.8) +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y = probability)) +
  # geom_line(data = simresults, aes(x = bce, y = mean), col = "red", lwd = 0.5) +
  geom_line(data = logip, aes(x = bce, y = prob_dens),
            linetype = "dashed", colour = "black") +
  labs(title = "Logistic", x = "BCE", y = "Summed probability") +
  # geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
  #               label = paste("p =", round(logisum$pvalue, 3)))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0),  add = c(0, 0.0001))) +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000),
                     limits = c(-10000, -2500),
                     expand = expansion(mult = c(0, 0))) +
  theme_bw()

##### Monte Carlo simulation - Uniform #####

names(unifp) <-  c("bce", "prob_dens", "model")

# Inspect SPD and fitted model
ggplot() +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y = probability)) +
  geom_line(data = unifp, aes(x = bce, y = prob_dens),
            linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000)) +
  theme_bw()

# Data frame for simulated dates
random_dates <- data.frame(matrix(nrow = ssize*nsim, ncol = 3))
names(random_dates) <- c("sample", "simn", "displacement_curve")


random_dates$sample = sample(unifp$bce, replace = TRUE,
                             size = ssize*nsim, prob = unifp$prob_dens)
random_dates$simn <- rep(1:nsim, each = nrow(random_dates)/nsim)
random_dates$displacement_curve <- incpolys[sample(incpolys$id,
                                              replace = TRUE,
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
  dplyr::mutate(
    prob_sum_normalised = prob_sum / (sum(prob_sum, na.rm = TRUE) * 5))

unisum <- simulation_summary(sumdates, simresults, cut_off = -2500, nsim)

simresults <- simresults %>%
  dplyr::group_by(bce) %>%
  dplyr::mutate(low = quantile(prob_sum_normalised, prob = 0.025,
                               na.rm = TRUE),
                high = quantile(prob_sum_normalised, prob = 0.975,
                                na.rm = TRUE),
                mean = mean(prob_sum_normalised)) %>%
  dplyr::distinct(bce, .keep_all = TRUE)


ggplot() +
  geom_vline(xintercept = as.numeric(names(expsum$busts)),
             col = "firebrick", alpha = 0.06) +
  geom_vline(xintercept = as.numeric(names(expsum$booms)),
             col = "darkgreen", alpha = 0.06) +
  ggplot2::geom_ribbon(data = simresults, aes(x = bce, ymin = low, ymax = high),
                       fill = "grey60", alpha = 0.8) +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y = probability)) +
  geom_line(data = expp, aes(x = bce, y = prob_dens),
            linewidth = 0.5, col = "red") +
  # geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
  #               label = paste("p =", round(expsum$pvalue, 3)))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0),  add = c(0, 0.0001))) +
  labs(x = "BCE", y = "Summed probability") +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000),
                     limits = c(-10000, -2500),
                     expand = expansion(mult = c(0, 0))) +
  theme_bw()

uniplt <- ggplot() +
  geom_vline(xintercept = as.numeric(names(unisum$busts)),
             col = "firebrick", alpha = 0.06) +
  geom_vline(xintercept = as.numeric(names(unisum$booms)),
             col = "darkgreen", alpha = 0.06) +
  ggplot2::geom_ribbon(data = simresults,
                       aes(x = bce, ymin = low, ymax = high),
                       fill = "grey60", alpha = 0.8) +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y = probability)) +
  geom_line(data = unifp, aes(x = year, y = pdf),
            linewidth = 0.5, colour = "red") + #"grey40"
  labs(title = "Uniform", x = "BCE", y = "Summed probability") +
  # geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
  #               label = paste("p =", round(unisum$pvalue, 3)))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0),  add = c(0, 0.0001))) +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000),
                     limits = c(-10000, -2500),
                     expand = expansion(mult = c(0, 0))) +
  theme_bw()

if(unisum$pvalue < 0.0001) {
  uniplt <- uniplt +
    geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
                  label = paste("p < 0.0001")))
} else {
  uniplt <- uniplt +
    geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
                  label = paste("p =", round(unisum$pvalue, 3))))
}

save(uniplt, file = "../external_data/shorespd/uniplt.rda")

expplt + logiplt + uniplt

ggsave(here::here("analysis/figures/shoreline_mc.png"),
       units = "px", width = 4000, height = 1250)
