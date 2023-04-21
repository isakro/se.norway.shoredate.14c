library(tidyverse)
library(sf)
library(terra)
library(here)
library(shoredate)
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

# Find weighted means dates for figure X in text
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
                  # Younger start date than 2500 BCE
out_of_range <- c("159969", "94301-1", "144111-1", "144113-1", "144114-1",
                  "138457-1", "230584-0", "265781-0", "265784-0",
                  # Earlier start date than 9465 BCE
                  "158183-1", "171065-1", "171069-1", "220991-1", "220992-1",
                  "244498-0", "250359-0", "262662-0", "270912-0", "270916-0",
                  "274933-0", "276689-0", "287819-0")

combined_sites <- filter(sites, !ask_id %in% out_of_range)

wplot <- ggplot() +
  geom_histogram(aes(x = wmeans, y = ..density..),
                 binwidth = 200,
                 col = "black", fill = "grey", alpha = 0.7) +
  geom_density(aes(x = wmeans),
               col = "black") +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000)) +
  labs(title = "Mean shoreline date", x = "BCE", y = "Density") +
  theme_bw()

splot <- ggplot() +
  geom_line(data = sumdatesdf, aes(x = sum.bce, y = sum.probability)) +
  scale_x_continuous(breaks = seq(-10000, 2500, 1000)) +
  labs(title = "SSPD", x = "BCE", y = "Summed probability") +
  theme_bw()

eplot <- ggplot() +
  geom_histogram(aes(x = combined_sites$elev, y = ..density..),
                 binwidth = 4,
                 col = "black", fill = "grey", alpha = 0.7) +
  geom_density(aes(x = combined_sites$elev),
             col = "black", alpha = 0.7) +
  labs(title = "Site elevations", x = "Elevation (m)", y = "Density") +
  scale_x_reverse() +
  theme_bw()


eplot + wplot + splot + plot_annotation(tag_levels = 'A')

st_write(sites,
         here("analysis/data/derived_data/combined_sites.gpkg"),
         append = FALSE)

# Create map of sites

# Exclude out of bounds sites from the two data sets
surveyed <- filter(surveyed, !askeladden_id %in% out_of_range)
excshore <- filter(excshore, !ask_id %in% out_of_range)

# Create label column for plotting
surveyed <- surveyed %>%  dplyr::mutate(lab = ifelse(quality < 4, "t", "f"))
excavated$lab <- "y"

worldmap <- st_read(here("analysis/data/raw_data/naturalearth_countries.gpkg"))
norway <- st_read(here("analysis/data/raw_data/naturalearth_norway.gpkg"))
muncipalities <- st_read(here("analysis/data/raw_data/municipalities.gpkg"))

# bboxpolys <- st_bbox(sites)
# bboxpolys[1] <- bboxpolys[1] - 900
# bboxpolys[3] <- bboxpolys[3] + 900
# bboxpolys[2] <- bboxpolys[2] - 900
# bboxpolys[4] <- bboxpolys[4] + 900
# bboxpolyspoly <- st_as_sf(st_as_sfc(bboxpolys))

sitbbox <- st_bbox(sites)
sitbbox[1:2] <- sitbbox[1:2] - 1000000
sitbbox[3:4] <- sitbbox[3:4] + 1000000
boundingpoly <- st_as_sf(st_as_sfc(sitbbox))

studyareabox <- st_bbox(sites)
studyareabox[1] <- studyareabox[1] - 15000
studyareabox[3] <- studyareabox[3] + 15000
studyareabox[2] <- studyareabox[2] - 5000
studyareabox[4] <- studyareabox[4] + 5000

# Reproject the bounding box to match world map, and crop the world map
# with the bounding box polygon.
bound_reproj <- st_transform(boundingpoly, st_crs(worldmap))
mapcountries <- worldmap %>%
  dplyr::filter(st_intersects(., bound_reproj, sparse = FALSE))
count_reproj <- st_transform(mapcountries, st_crs(sites))

overview <-
  ggplot() +
  geom_sf(data = count_reproj, fill = "grey", colour = NA) +
  geom_sf(data = bboxpolyspoly,
          fill = NA, colour = "black", size = 0.5) +
  coord_sf(xlim = c(sitbbox[1], sitbbox[3]), ylim = c(sitbbox[2], sitbbox[4]),
           expand = FALSE) +
  theme_bw() + theme(axis.title=element_blank(),
                     axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     rect = element_rect(),
                     axis.ticks = element_blank(),
                     panel.grid.major = element_blank())

norw_reproj <- st_transform(norway, st_crs(surveyed))

anc <- as.numeric(c(sitbbox$ymin, sitbbox$xmax))
sa <- ggplot() +
  geom_sf(data = norw_reproj, fill = "grey", colour = NA) +
  geom_sf(data = st_transform(muncipalities, st_crs(sites)), fill = NA,
          colour = "black", lwd = 0.25) +
  geom_sf(data = st_centroid(surveyed), aes(fill = lab),
          size = 2, shape = 21,
          colour = "black", show.legend = "point") +
  geom_sf(data = st_centroid(dplyr::filter(surveyed, quality < 4)),
          aes(fill = lab),
          size = 2, shape = 21,
          colour = "black", show.legend = "point") +
  geom_sf(data = st_centroid(excshore), aes(fill = lab),
          size = 2, shape = 21,
          colour = "black", show.legend = "point") +
  scale_fill_manual(labels = c(paste0("Surveyed sites,\nincluded (n = ",
                                      nrow(dplyr::filter(surveyed, quality < 4)) ,")"),
                               paste0("Surveyed sites,\nexcluded (n = ",
                                      nrow(dplyr::filter(surveyed, quality > 3)) ,")"),
                               paste0("Excavated sites (n = ", nrow(excshore) ,
                                      ")")),
                    values = c("t" = "darkgoldenrod1", "f" = "black", "y" = "white"),
                    name = "") +
  ggsn::scalebar(data = sites, dist = 20, dist_unit = "km",
                 transform = FALSE, st.size = 4, height = 0.02,
                 border.size = 0.1, st.dist = 0.03,
                 anchor = c(x = anc[2] - 15500, y = anc[1]) + 8000) +
  coord_sf(xlim = c(studyareabox[1], studyareabox[3]),
           ylim = c(studyareabox[2], studyareabox[4]),
           expand = FALSE) +
  theme_bw() + theme(axis.title=element_blank(),
                     axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     rect = element_rect(),
                     axis.ticks = element_blank(),
                     panel.grid.major = element_blank(),
                     legend.position = "bottom",
                     legend.text=element_text(size = 11))

cowplot::ggdraw() +
  cowplot::draw_plot(sa) +
  cowplot::draw_plot(overview, x = 0.126,
                     y = 0.65, width = 0.35, height = 0.35)


# overview + sa +   plot_layout(widths = c(1, 2.5))
ggsave(here::here("analysis/figures/map.png"),
       units = "px", width = 2194*1.7, height = 2380)
