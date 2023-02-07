library(ggplot2)
library(shoredate)
library(dplyr)
library(sf)
library(here)

set.seed(1)

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
predgrid <- data.frame(bce = sumdatesdf$sum.bce, prob_dens = est)

sumdatesdf <- filter(sumdatesdf, sum.bce <= -2500)
predgrid <- filter(predgrid, bce <= -2500)

# Inspect SPD and fitted model
ggplot() +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y = sum.probability)) +
  geom_line(data = predgrid, aes(x = bce, y = prob_dens),
            linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000)) +
  theme_bw()

# Data frame for simulated dates
random_dates <- data.frame(matrix(nrow = ssize*nsim, ncol = 3))
names(random_dates) <- c("sample", "simn", "displacement_curve")


random_dates$sample = sample(predgrid$bce, replace = TRUE,
                             size = ssize*nsim, prob = predgrid$prob_dens)
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


modelsum <- simulation_summary(sumdates, simresults, cut_off = -2500, nsim)

simresults <- simresults %>%
  group_by(bce) %>%
  mutate(low = quantile(prob_sum_normalised, prob = 0.025, na.rm = TRUE),
         high = quantile(prob_sum_normalised, prob = 0.975, na.rm = TRUE),
         mean = mean(prob_sum_normalised)) %>%
  distinct(bce, .keep_all = TRUE)

ggplot() +
  ggplot2::geom_ribbon(data = simresults, aes(x = bce, ymin = low, ymax = high),
                       fill = "grey", alpha = 0.9) +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y = sum.probability)) +
  # geom_line(data = simdates, aes(x = bce, y = mean), col = "red", lwd = 0.5) +
  geom_line(data = predgrid, aes(x = bce, y = prob_dens),
            linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000),
                     limits = c(-10000, -2500)) +
  theme_bw()

###### Logistic model ######
x = sumdatesdf$sum.bce
y = sumdatesdf$sum.probability

logfit <- nls(y ~ SSlogis(x, Asym, xmid, scal))

predgrid <- data.frame(bce = x, prob_dens = predict(logfit))
predgrid <- filter(predgrid, bce <= -2500)

# Inspect SPD and fitted model
ggplot() +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y = sum.probability)) +
  geom_line(data = predgrid, aes(x = bce, y = prob_dens),
            linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000)) +
  ggtitle("Logistic model") +
  theme_bw()


nsim <- 1000 #1000

# Data frame for simulated dates
random_dates <- data.frame(matrix(nrow = ssize*nsim, ncol = 3))
names(random_dates) <- c("sample", "simn", "displacement_curve")


random_dates$sample = sample(predgrid$bce, replace = TRUE,
                             size = ssize*nsim, prob = predgrid$prob_dens)
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
for (i in 92101:nrow(random_dates)){
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

simresults <- do.call(rbind.data.frame, results)

pval <- simulation_summary(sumdates, simresults)

simresults <- simresults %>%
  dplyr::group_by(simn) %>%
  dplyr::mutate(prob_sum_normalised = prob_sum / sum(prob_sum, na.rm = TRUE)) %>%
  dplyr::group_by(bce) %>%
  dplyr::mutate(low = quantile(prob_sum_normalised, prob = 0.025, na.rm = TRUE),
                high = quantile(prob_sum_normalised, prob = 0.975, na.rm = TRUE),
                mean = mean(prob_sum_normalised)) %>%
  dplyr::distinct(bce, .keep_all = TRUE)

ggplot() +
  geom_vline(xintercept = as.numeric(names(modelsum$busts)),
             col = "firebrick", alpha = 0.05) +
  geom_vline(xintercept = as.numeric(names(modelsum$booms)),
             col = "darkgreen", alpha = 0.05) +
  ggplot2::geom_ribbon(data = simresults, aes(x = bce, ymin = low, ymax = high),
                       fill = "darkgrey", alpha = 0.9) +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y = sum.probability)) +
  # geom_line(data = simresults, aes(x = bce, y = mean), col = "red", lwd = 0.5) +
  geom_line(data = predgrid, aes(x = bce, y = prob_dens),
            linetype = "dashed", colour = "black") +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000),
                     limits = c(-10000, -2500)) +
  theme_bw()


###### Uniform model ######

# Fit uniform model (re-purposed from modelTest() from rcarbon)
predgrid <- data.frame(bce = sumdatesdf$sum.bce,
                       prob_dens = mean(sumdatesdf$sum.probability))
predgrid <- filter(predgrid, bce <= -2500)

# Inspect SPD and fitted model
ggplot() +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y = sum.probability)) +
  geom_line(data = predgrid, aes(x = bce, y = prob_dens),
            linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000)) +
  theme_bw()

# Data frame for simulated dates
random_dates <- data.frame(matrix(nrow = ssize*nsim, ncol = 3))
names(random_dates) <- c("sample", "simn", "displacement_curve")


random_dates$sample = sample(predgrid$bce, replace = TRUE,
                             size = ssize*nsim, prob = predgrid$prob_dens)
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
for (i in 184201:nrow(random_dates)){
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

modelsum <- simulation_summary(sumdates, simresults, cut_off = -2500, nsim)

simresults <- simresults %>%
  dplyr::group_by(bce) %>%
  dplyr::mutate(low = quantile(prob_sum_normalised, prob = 0.025, na.rm = TRUE),
         high = quantile(prob_sum_normalised, prob = 0.975, na.rm = TRUE),
         mean = mean(prob_sum_normalised)) %>%
  dplyr::distinct(bce, .keep_all = TRUE)

ggplot() +
  geom_vline(xintercept = as.numeric(names(modelsum$busts)),
             col = "firebrick", alpha = 0.05) +
  geom_vline(xintercept = as.numeric(names(modelsum$booms)),
             col = "darkgreen", alpha = 0.05) +
  ggplot2::geom_ribbon(data = simresults, aes(x = bce, ymin = low, ymax = high),
                       fill = "darkgrey", alpha = 0.9) +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y = sum.probability)) +
  # geom_line(data = simresults, aes(x = bce, y = mean), col = "red", lwd = 0.5) +
  geom_line(data = predgrid, aes(x = bce, y = prob_dens),
            linetype = "dashed", colour = "black") + #"grey40"
  labs(x = "BCE", y = "Summed probability") +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000),
                     limits = c(-10000, -2500)) +
  theme_bw()


