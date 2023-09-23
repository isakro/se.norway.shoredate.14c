#### Create maps and plots providing example of shoreline dating #####

library(ggplot2)
library(patchwork)
library(sf)
library(here)
library(shoredate)

load(here("analysis/data/derived_data/prepared_sitedata.RData"))
load(here("analysis/data/derived_data/rcarbon_models_pd.RData"))
surveyed <- st_read(here("analysis/data/raw_data/surveyed_sites.gpkg"))
excavated <- st_zm(read_sf(here("analysis/data/raw_data/excavated_sites.gpkg")))

# Create label column for plotting
surveyed <- surveyed %>%  dplyr::mutate(lab = ifelse(quality < 4, "t", "f"))
excavated$lab <- "y"

worldmap <- st_read(here("analysis/data/raw_data/naturalearth_countries.gpkg"))
norway <- st_read(here("analysis/data/raw_data/naturalearth_norway.gpkg"))
muncipalities <- st_read(here("analysis/data/raw_data/municipalities.gpkg"))

sitbbox <- st_bbox(sites)
sitbbox[1:2] <- sitbbox[1:2] - 1000000
sitbbox[3:4] <- sitbbox[3:4] + 1000000
boundingpoly <- st_as_sf(st_as_sfc(sitbbox))

studyareabox <- st_bbox(sites)
studyareabox[1] <- studyareabox[1] - 15000
studyareabox[3] <- studyareabox[3] + 15000
studyareabox[2] <- studyareabox[2] - 5000
studyareabox[4] <- studyareabox[4] + 5000
bboxpolyspoly <- st_as_sf(st_as_sfc(studyareabox))

# Reproject the bounding box to match world map, and crop the world map
# with the bounding box polygon.
bound_reproj <- st_transform(boundingpoly, st_crs(worldmap))
mapcountries <- worldmap %>%
  dplyr::filter(st_intersects(., bound_reproj, sparse = FALSE))
count_reproj <- st_transform(mapcountries, st_crs(surveyed))

overview <-
  ggplot() +
  geom_sf(data = count_reproj, fill = "grey", colour = NA) +
  geom_sf(data = bboxpolyspoly,
          fill = NA, colour = "black", size = 1) +
  coord_sf(xlim = c(sitbbox[1], sitbbox[3]), ylim = c(sitbbox[2], sitbbox[4]),
           expand = FALSE) +
  theme_bw() + theme(axis.title = element_blank(),
                     axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     rect = element_rect(),
                     axis.ticks = element_blank(),
                     panel.grid.major = element_blank())

norw_reproj <- st_transform(norway, st_crs(surveyed))

sa <- ggplot() +
  geom_sf(data = norw_reproj, fill = "grey", colour = NA) +
  geom_sf(data = st_transform(muncipalities, st_crs(surveyed)), fill = NA,
          colour = "black", lwd = 0.25) +
  geom_sf(data = st_centroid(surveyed), aes(fill = lab),
          size = 1.5, shape = 21,
          colour = "black", show.legend = "point") +
  geom_sf(data = st_centroid(dplyr::filter(surveyed, quality < 4)),
          aes(fill = lab),
          size = 1.5, shape = 21,
          colour = "black", show.legend = "point") +
  geom_sf(data = st_centroid(excshore), aes(fill = "y"),
          size = 1.5, shape = 21,
          colour = "black", show.legend = "point") +
  scale_fill_manual(labels = c(paste0("Surveyed sites,\nincluded (n = ",
                                      nrow(dplyr::filter(surveyed,
                                                         quality < 4)) ,")"),
                               paste0("Surveyed sites,\nexcluded (n = ",
                                      nrow(dplyr::filter(surveyed,
                                                         quality > 3)) ,")"),
                               paste0("Excavated, shoreline\ndated sites (n = ",
                                      nrow(excshore) ,
                                      ")")),
                    values = c("t" = "darkgoldenrod1",
                               "f" = "black",
                               "y" = "white"),
                    name = "") +
  coord_sf(xlim = c(studyareabox[1], studyareabox[3]),
           ylim = c(studyareabox[2], studyareabox[4]),
           expand = FALSE) +
  ggspatial::annotation_scale(
    location = "br",
    width_hint = 0.3,
    style = "ticks") +
  theme_bw() + theme(axis.title=element_blank(),
                     axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     rect = element_rect(),
                     axis.ticks = element_blank(),
                     panel.grid.major = element_blank(),
                     legend.position = "bottom",
                     legend.text=element_text(size = 11))

# Retrieve site names for plotting of radiocarbon data
c14_names <- tools::file_path_sans_ext(colnames(PD))

c14_sites <- excavated %>%
  dplyr::filter(site_name %in% c14_names)

pltc14 <- ggplot() +
  geom_sf(data = norw_reproj, fill = "grey", colour = NA) +
  geom_sf(data = st_transform(muncipalities, st_crs(surveyed)), fill = NA,
          colour = "black", lwd = 0.25) +
  geom_sf(data = st_centroid(c14_sites), aes(fill = lab),
          size = 1.5, shape = 24,
          colour = "black", show.legend = "point") +
  scale_fill_manual(values = c("y" = "red"),
                    labels = c(paste0("Sites with radiocarbon dates (n = ",
                                      nrow(c14_sites), ")")),
                    name = "") +
  coord_sf(xlim = c(studyareabox[1], studyareabox[3]),
           ylim = c(studyareabox[2], studyareabox[4]),
           expand = FALSE) +
  ggspatial::annotation_scale(
    location = "br",
    width_hint = 0.3,
    style = "ticks") +
  theme_bw() + theme(axis.title=element_blank(),
                     axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     rect = element_rect(),
                     axis.ticks = element_blank(),
                     panel.grid.major = element_blank(),
                     legend.position = "bottom",
                     legend.text=element_text(size = 11))

plt1 <- cowplot::ggdraw() +
  cowplot::draw_plot(sa) +
  cowplot::draw_plot(overview, x = 0.179,
                     y = 0.7, width = 0.3, height = 0.3)

(plt1)/
  pltc14 + plot_annotation(tag_levels = "A")

# Example date plot

target_site <- sites[535,]
target_curve <- interpolate_curve(target_site)
tpq_date <- shoreline_date(target_site,
                           elevation = target_site$elev,
                           model = "none", hdr_prob = 1)
shdate <-  shoreline_date(target_site,
                          elevation = target_site$elev)

tplt <- target_plot(target_site, basemap = norw_reproj) +
  coord_sf(xlim = c(studyareabox[1], studyareabox[3]),
           ylim = c(studyareabox[2], studyareabox[4]),
           expand = FALSE)


dplt <- displacement_plot(target_curve,
                          displacement_alpha = 0.4)
dplt <- dplt + theme(legend.position = "left", legend.direction = "vertical")
shoredate_plot(tpq_date,
               date_probability_scale = 6000,
               greyscale = TRUE,
               hdr_label_yadj = 0.5)
tpqplt <- last_plot()

shoredate_plot(shdate, greyscale = TRUE)
splt <- last_plot()

(tplt + dplt) /
  (tpqplt + splt) + plot_annotation(tag_levels = "A")
