library(shoredate)
library(here)
library(ADMUR)


set.seed(1)

source(here("R/03_functions.R"))

# Load polygons with displacement curves
load(here("analysis/data/derived_data/displacement_polygons.RData"))

# Load shoreline dates
load(here("analysis/data/derived_data/sdates.RData"))

# Load CPL models
load(file = here::here("analysis/data/derived_data/shore_models.RData"))

cpl4 <- convertPars(pars = cpl_4$par, years = minage:maxage, type = 'CPL')
cpl4$model <- "4-CPL"

cpl4$year <- (cpl4$year - 1950) * -1
names(cpl4) <- c("bce", "prob_dens", "model")

# Sum the probability distributions for the shoreline dates
sumdates <- sum_shoredates(shorelinedates)
sumdatesdf <- as.data.frame(sumdates) %>%
  filter(sum.probability != 0)

# Restrict SSPD and normalise probability
sumdatesdf <- filter(sumdatesdf, sum.bce <= max(cpl4$bce))
sumdatesdf$probability <- sumdatesdf$sum.probability /
  (sum(sumdatesdf$sum.probability) * 5)

ggplot() +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y= probability)) +
  geom_line(data = cpl4, aes(x = bce, y = prob_dens))

# Find density of sites across polygons
incpolys$dens <- lengths(st_intersects(incpolys, sites)) /
  sum(lengths(st_intersects(incpolys, sites)))

# Sample size
ssize = sumdates$dates_n

# Number of simulations
nsim <- 10000

##### Monte Carlo simulation - 4-CPL #####

# Create data frame to hold simulated dates
random_dates <- data.frame(matrix(nrow = ssize*nsim, ncol = 3))
names(random_dates) <- c("sample", "simn", "displacement_curve")

random_dates$sample = sample(cpl4$bce, replace = TRUE,
                             size = ssize*nsim, prob = cpl4$prob_dens)
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
         file = paste0(here::here("../external_data/shorespd/4cpl/"),
                       i, "shore.rds"))

    simdates <- vector("list", ssize)
  }
}
end_time <- Sys.time()
end_time - start_time


resfiles <- list.files(here::here("../external_data/shorespd/4cpl/"),
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

mod <- approx(x = cpl4$bce, y = cpl4$prob_dens,
              xout = as.numeric(rownames(SPD)),
              ties = 'ordered', rule = 2)$y

# Model summary (can also be plotted with ADMUR)
cpl4_summary <- simulation_summary(SPD, simulation_results,
                                  c(-2500, -9445), mod, ncol(pd))

# ADMUR plot
plotSimulationSummary(cpl4_summary)

# Custom plot
cpl4plts <- plot_mc(cpl4_summary)

save(cpl4plts, file = "../external_data/shorespd/cpl4plts.rda")
