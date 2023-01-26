library(tidyverse)
library(sf)
library(terra)
library(here)
library(shoredate)

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

#
# surveyed <- surveyed %>% dplyr::mutate(lab = ifelse(quality < 4, "t", "f"))
# excavated$lab <- "y"

# Select and rename columns needed to perform shoreline dating
rcolexcshore <- excshore %>% dplyr::select(ask_id, elev)
rcolsurvq <- survq  %>% dplyr::select(askeladden_id, elev) %>%
  dplyr::rename(ask_id = askeladden_id)

# Combine the data
sites <- rbind(rcolexcshore, rcolsurvq)

###### Shoreline dating ######
step_reso = 0.1

shorelinedates <- shoredate::shoreline_date(sites,
                                            elevation = dtm,
                                            elev_reso = step_reso,
                                            cal_reso = 5,
                                            sparse = TRUE,
                                            verbose = TRUE)

# sdates <- do.call(rbind, lapply(shorelinedates, as.data.frame))

save(step_reso, shorelinedates,
     file = here("analysis/data/derived_data/sdates.RData"))


st_write(sites, here("analysis/data/derived_data/combined_sites.gpkg"))
