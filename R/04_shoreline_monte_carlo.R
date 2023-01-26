library(ggplot2)
library(shoredate)
library(here)

# Source functions
source(here("R/03_functions.R"))

# Load polygons with displacement curves
load(here("analysis/data/derived_data/displacement_polygons.RData"))

# Load shoreline dates
load(here("analysis/data/derived_data/sdates.RData"))

# Sum the probability distributions for the shoreline dates
sumdates <- sum_shoredates(shorelinedates)

# Make the results a data frame
sumdatesdf <- as.data.frame(sumdates)


# Code below is repurposed from modelTest() from rcarbon

# Fit exponential model
fit <- nls(y ~ (exp(a + b * x)),
           data = data.frame(x = sumdatesdf$sum.bce,
                             y = sumdatesdf$sum.probability),
           start = list(a = 0, b = 0))
est <- predict(fit, list(x = sumdatesdf$sum.bce))
predgrid <- data.frame(bce = sumdatesdf$sum.bce, prob_dens = est)

# Inspect SPD and fitted model
shoredate_sumplot(sumdates)  +
  geom_line(data = predgrid, aes(x = bce, y = prob_dens)) +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000)) +
  theme_classic()

# Find density of site across polygons
incpolys$dens <- lengths(st_intersects(incpolys, sites)) /
  sum(lengths(st_intersects(incpolys, sites)))

# Sample size
ssize = sumdates$dates_n

# Number of simulations
nsim <- 10

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


simdates <- list()
start_time <- Sys.time()
for (i in 1:nrow(random_dates)){
  print(paste(i, "/", nrow(random_dates)))

  simdates[[i]] <- shoredate::shoreline_date(
    site = random_dates$simn[i],
    elev_reso = step_reso,
    cal_reso = 5,
    elevation = random_dates$reverse_elevation[i],
    interpolated_curve = random_dates$displacement_curve[i],
    sparse = TRUE)

  if(i %% ssize == 0){
    # Save data externally and overwrite simdates for memory purposes
      tmp <- bind_rows(simdates) %>% group_by(bce, site_name) %>%
        filter(!is.na(probability)) %>%
        summarise(prob_sum = sum(probability))

      save(tmp,
           file = paste0(here::here("../external_data/shorespd/exp/"),
                         i, "shore.rds"))

      simdates <- list()
  }
}
end_time <- Sys.time()
end_time - start_time
