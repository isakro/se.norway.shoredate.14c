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

dates_dfs <- lapply(shorelinedates, as.data.frame)

# Combine these into a single data frame
sdates <- do.call(rbind, dates_dfs)

# Assign simulation results for each site a unique ID
sdates$site <- rep(1:length(dates_dfs) , each = nrow(sdates)/length(dates_dfs))

pd <- as.data.frame(sdates) %>%
  dplyr::filter(bce <= -2500 & bce >= -9445) %>%
  dplyr::group_by(bce) %>%
  dplyr::filter(!is.na(probability)) %>%
  dplyr::group_by(site) %>%
  tidyr::pivot_wider(names_from = site, values_from = probability) %>%
  tibble::column_to_rownames("bce")

pd <- sweep(pd, 2, colSums(pd),"/")

SPD <- as.data.frame(rowSums(pd))
SPD <- SPD/( sum(SPD) *5 )

# Maximum likelihood search using DEoptimR
exp <- JDEoptim(lower = -0.01, upper = 0.01, fn = objectiveFunction,
                PDarray = pd, type = "exp", NP = 20)
logi <- JDEoptim(lower = c(0, 0000), upper = c(1, 10000),
                 fn = objectiveFunction, PDarray = pd, type = 'logistic',
                 NP = 40, trace = TRUE)
# No search required for uniform
unif <- -objectiveFunction(pars = NULL, PDarray = pd, type = 'uniform')

expp <- convertPars(pars = exp$par, years = -9445:-2500, type = 'exp')
expp$model <- "Exponential"
logip <- convertPars(pars = logi$par, years = -9445:-2500, type = 'logistic')
logip$model <- "Logistic"
unifp <- convertPars(pars = NULL, years =  -9445:-2500, type = "uniform")
unifp$model <- "Uniform"


# Combine models for plotting
shoremodels <- rbind(expp, logip, unifp)

# Call to plot
ggplot() +
  geom_bar(aes(x = as.numeric(rownames(SPD)), SPD[,1]),
           stat = "identity", col = "grey") +
  geom_line(aes(x = as.numeric(rownames(SPD)),
                y = SPD[,1])) +
  geom_line(data = shoremodels, aes(year, pdf, col = model),
            linewidth = 1.1)


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
nsim <- 10

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

# Retrieve results and normalise each simulated SPD
simresults <- do.call(rbind.data.frame, results) %>%
  group_by(simn) %>%
  mutate(prob_sum_normalised = prob_sum /(sum(prob_sum, na.rm = TRUE) * 5))

expsum <- simulation_summary(sumdates, simresults, cut_off = -2500, nsim)

simresults <- simresults %>%
  group_by(bce) %>%
  mutate(low = quantile(prob_sum_normalised, prob = 0.025, na.rm = TRUE),
         high = quantile(prob_sum_normalised, prob = 0.975, na.rm = TRUE),
         mean = mean(prob_sum_normalised)) %>%
  distinct(bce, .keep_all = TRUE)

expplt <- ggplot() +
  geom_vline(xintercept = as.numeric(names(expsum$busts)),
             col = "firebrick", alpha = 0.04) +
  geom_vline(xintercept = as.numeric(names(expsum$booms)),
             col = "darkgreen", alpha = 0.04) +
  ggplot2::geom_ribbon(data = simresults, aes(x = bce, ymin = low, ymax = high),
                       fill = "grey60", alpha = 0.8) +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y = probability)) +
  geom_line(data = expp, aes(x = bce, y = prob_dens),
            linetype = "dashed", colour = "black") +
  # geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
  #               label = paste("p =", round(expsum$pvalue, 3)))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0),  add = c(0, 0.0001))) +
  labs(title = "Exponential", x = "BCE", y = "Summed probability") +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000),
                     limits = c(-10000, -2500),
                     expand = expansion(mult = c(0, 0))) +
  theme_bw()

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

uniplt <- ggplot() +
  geom_vline(xintercept = as.numeric(names(unisum$busts)),
             col = "firebrick", alpha = 0.04) +
  geom_vline(xintercept = as.numeric(names(unisum$booms)),
             col = "darkgreen", alpha = 0.04) +
  ggplot2::geom_ribbon(data = simresults,
                       aes(x = bce, ymin = low, ymax = high),
                       fill = "grey60", alpha = 0.8) +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y = probability)) +
  geom_line(data = unifp, aes(x = bce, y = prob_dens),
            linetype = "dashed", colour = "black") + #"grey40"
  labs(title = "Uniform", x = "BCE", y = "Summed probability") +
  # geom_text(aes(-Inf, Inf, hjust = -0.25, vjust = 2,
  #               label = paste("p =", round(unisum$pvalue, 3)))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0),  add = c(0, 0.0001))) +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000),
                     limits = c(-10000, -2500),
                     expand = expansion(mult = c(0, 0))) +
  theme_bw()
uniplt
