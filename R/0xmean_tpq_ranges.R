# Find the average range of TPQ ranges for the geological displacement curves,
# to be reported in the main text.

library(shoredate)

# Load centre points for the displacement curves
centrepoints <- sf::st_read(
  system.file("extdata/isobase_centrepts.gpkg",
              package = "shoredate",
              mustWork = TRUE), quiet = TRUE)

# Horten

# Retreieve the point on the Horten isobase
horten <- centrepoints[1,]

# Create a sequence of eleveation values from well above the marine limit for
# the curves and down to 20 masl
elevations <- seq(190, 1, -1)

# Make these a series of "sites" to be dated
sites <- sf::st_sf(rep(sf::st_geometry(horten), length(elevations)))
sites$name <- paste(elevations, "masl")

# Perform dating procedure, setting model to "none" and HDR to 1 to account for
# the full range.
hsdates <- shoreline_date(sites,
                          elevation = elevations,
                          hdr_prob = 1,
                          model = "none")

# Plot for inspection
shoredate_plot(hsdates, date_probability = FALSE, multiplot = TRUE)

# Porsgrunn curve
porsgrunn <- centrepoints[2,]
sites <- sf::st_sf(rep(sf::st_geometry(porsgrunn), length(elevations)))
sites$name <- paste(elevations, "masl")
psdates <- shoreline_date(sites,
                          elevation = elevations,
                          hdr_prob = 1,
                          model = "none")

shoredate_plot(psdates, date_probability = FALSE, multiplot = TRUE)

# Arendal curve
arendal <- centrepoints[3,]
sites <- sf::st_sf(rep(sf::st_geometry(arendal), length(elevations)))
sites$name <- paste(elevations, "masl")
asdates <- shoreline_date(sites,
                          elevation = elevations,
                          hdr_prob = 1,
                          model = "none")

shoredate_plot(asdates, date_probability = FALSE, multiplot = TRUE)

# Tvedestrand curve
tvedestrand <- centrepoints[4,]
sites <- sf::st_sf(rep(sf::st_geometry(tvedestrand), length(elevations)))
sites$name <- paste(elevations, "masl")
tsdates <- shoreline_date(sites,
                          elevation = elevations,
                          hdr_prob = 1,
                          model = "none")

shoredate_plot(tsdates, date_probability = FALSE, multiplot = TRUE)

# Find the date ranges
prange <- c()
for(i in 1:length(psdates)){
  tmp <- psdates[[i]][[1]]$hdr_end - psdates[[i]][[1]]$hdr_start
  prange <- c(prange, tmp)
}

hrange <- c()
for(i in 1:length(hsdates)){
  tmp <- hsdates[[i]][[1]]$hdr_end - hsdates[[i]][[1]]$hdr_start
  hrange <- c(hrange, tmp)
}

arange <- c()
for(i in 1:length(asdates)){
  tmp <- asdates[[i]][[1]]$hdr_end - asdates[[i]][[1]]$hdr_start
  arange <- c(arange, tmp)
}

trange <- c()
for(i in 1:length(tsdates)){
  tmp <- tsdates[[i]][[1]]$hdr_end - tsdates[[i]][[1]]$hdr_start
  trange <- c(trange, tmp)
}

# Find the mean and standard deviation across all date ranges
tpqmean <- mean(c(prange, hrange, arange, trange), na.rm = TRUE)
tpqsd <- sd(c(prange, hrange, arange, trange), na.rm = TRUE)

# save results
save(tpqmean, tpqsd,
     file = here("analysis/data/derived_data/disp_tpq_ranges.RData"))

library(dplyr)

dispcurves <- get(load(system.file("extdata/displacement_curves.rda",
              package = "shoredate")))

rocdisp <- dispcurves %>%
  group_by(name) %>%
  arrange(bce, .by_group = TRUE) %>%
  mutate(rate_upper = upperelev - lag(upperelev),
         rate_lower = lowerelev - lag(lowerelev)) %>%
  rowwise() %>%
  mutate(roc_mean = mean(c(rate_upper, rate_lower), na.rm = TRUE))

rocdisp <- rocdisp %>%
  group_by(name) %>%
  mutate( roc_mean = rowMeans(select(rocdisp,
                  c(rate_upper, rate_lower)), na.rm = TRUE))

ggplot(data = rocdisp) +
  # geom_line(aes(x = bce, y = rate_upper, col = name)) +
  # geom_line(aes(x = bce, y = rate_lower, col = name)) +
  geom_line(aes(x = bce, y = roc_mean, col = name)) +
  scale_x_continuous(breaks = seq(-10000, 0, 1000))
