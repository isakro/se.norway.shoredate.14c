# The first part of this script prepares data for shoreline dating and the
# creates the SPD. The second part creates plots exemplifying different ways of
# summing shoreline dates.

library(tidyverse)
library(sf)
library(terra)
library(here)
library(shoredate) # Note that shoredate has to be v.1.0.2
library(patchwork)

# For reproducibility
set.seed(1)

# Load data
surveyed <- st_read(here("analysis/data/raw_data/surveyed_sites.gpkg"))
excavated <- st_zm(read_sf(here("analysis/data/raw_data/excavated_sites.gpkg")))
dtm <- rast("/home/isak/phd/eaa_presentation/dtm10/dtm10.tif")
site_dat <- read.csv(here("analysis/data/raw_data/sites.csv"))

###### Prepare site data for shoreline dating ######

# Exclude sites that haven't been excavated and which haven't been shoreline
# dated
shore_sites <- site_dat %>% dplyr::filter(investigated == "t") %>%
  dplyr::filter(dating_method == "shore_typo" & !is.na(reported_earliest_bce)
                & !is.na(reported_min_elev))

# Join shoreline dated excavated sites with the spatial data for the sites
excshore <- st_as_sf(
  dplyr::left_join(shore_sites, excavated,
                    by = c("name" = "site_name", "ask_id"))
  ) %>%
  dplyr::filter(!(is.na(st_dimension(geom))))

# Find mean reported elevation for the excavated sites
excshore$elev <- rowMeans(subset(st_drop_geometry(excshore),
                                 select = c(reported_min_elev,
                                            reported_max_elev)))

# Find elevation of surveyed sites
surveyed$elev <- terra::extract(dtm, vect(surveyed), fun = mean)[, -1]

# Exclude surveyed sites with a quality score worse than 3
survq <- dplyr::filter(surveyed, quality < 4)

# Select and rename columns needed to perform shoreline dating
rcolexcshore <- excshore %>% dplyr::select(ask_id, elev)
rcolsurvq <- survq  %>% dplyr::select(askeladden_id, elev) %>%
  dplyr::rename(ask_id = askeladden_id)

# Combine the data
sites <- rbind(rcolexcshore, rcolsurvq)

###### Shoreline dating ######
step_reso = 0.1

# Find weighted means dates for Figure 3 in the main text
wmean_dates <- shoreline_date(sites, elevation = dtm,
                              elev_reso = step_reso,
                              cal_reso = 5,
                              model = "none",
                              verbose = TRUE)

wmeans <- c()
for(i in 1:length(wmean_dates)){
  wmeans <- c(wmeans, wmean_dates[[i]][[1]]$weighted_mean)
}
wmeans <- wmeans[wmeans  <= -2500]

# Perform shoreline dating using the probabilistic method
shorelinedates <- shoreline_date(sites,
                                 elevation = dtm,
                                 elev_reso = step_reso,
                                 cal_reso = 5,
                                 sparse = TRUE,
                                 verbose = TRUE)


save(step_reso, wmeans, shorelinedates,
     file = here("analysis/data/derived_data/sdates.RData"))

# Sum the probability distributions for the shoreline dates
sumdates <- sum_shoredates(shorelinedates)
sumdatesdf <- as.data.frame(sumdates) %>%
  filter(sum.probability != 0)
sumdatesdf <- filter(sumdatesdf, sum.bce <= -2500)

# Exclude sites that were given a NA date (elevation equates to a date older
# than 9465 BCE or younger than 2500 BCE).

                  # Later start date than 2500 BCE
out_of_range <- c("159969", "94301-1", "144111-1", "144113-1", "144114-1",
                  "138457-1", "230584-0", "265781-0", "265784-0",
                  # Earlier start date than 9465 BCE
                  "158183-1", "171065-1", "171069-1", "220991-1", "220992-1",
                  "244498-0", "250359-0", "262662-0", "270912-0", "270916-0",
                  "274933-0", "276689-0", "287819-0")

combined_sites <- filter(sites, !ask_id %in% out_of_range)

st_write(combined_sites,
         here("analysis/data/derived_data/combined_sites.gpkg"),
         append = FALSE)

# Exclude out of bounds sites from the two data sets for plotting
surveyed <- filter(surveyed, !askeladden_id %in% out_of_range)
excshore <- filter(excshore, !ask_id %in% out_of_range)

save(sites, surveyed, excshore,
     file = here("analysis/data/derived_data/prepared_sitedata.RData"))

# Create plot of site elevations, mean TPQ dates and SPD of shoreline dates

wplot <- ggplot() +
  geom_histogram(aes(x = wmeans, y = ..density..),
                 binwidth = 200,
                 col = "grey", fill = "grey") +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000),
                     limits = c(-9500, -2500)) +
  labs(title = "Mean TPQ dates", x = "BCE", y = "Density") +
  theme_bw()

splot <- ggplot() +
  # geom_line(data = sumdatesdf, aes(x = sum.bce, y = sum.probability)) +
  geom_bar(data = sumdatesdf, aes(x = sum.bce, y = sum.probability),
                 stat = "identity",
                 col = "grey", fill = "grey") +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000)) +
  labs(title = "Summed probability distribution\nof shoreline dates",
       x = "BCE", y = "Summed probability") +
  theme_bw()

eplot <- ggplot() +
  geom_histogram(aes(x = combined_sites$elev, y = ..density..),
                 binwidth = 4,
                 col = "grey", fill = "grey") +
  # geom_density(aes(x = combined_sites$elev),
  #            col = "black", alpha = 0.7) +
  labs(title = "Site elevations", x = "Elevation (m)", y = "Density") +
  scale_x_reverse() +
  theme_bw()


eplot + wplot + splot + plot_annotation(tag_levels = 'A')
