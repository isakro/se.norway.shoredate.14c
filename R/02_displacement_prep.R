# This script creates 2km wide line segments that run perpendicular to the
# shoreline gradient, interpolates a displacement curve to centre of these,
# and finds the distribution of sites across them.

library(tidyverse)
library(sf)
library(terra)
library(here)
library(patchwork)

sites <- st_read(here("analysis/data/derived_data/combined_sites.gpkg"))
bmap <- st_read(here("analysis/data/raw_data/naturalearth_norway.gpkg")) %>%
  st_transform(st_crs(sites))

# Specify distance of increments for which to interpolate displacement curve
increment <- 2000
# Degree direction of isobases
deg <- 327
# Arbritrary long distance for lines
linedist <- 100000

# Find bounding box of sites
outcoords <- st_bbox(sites)
# Find bottom left and top right points
startpt <- st_point(outcoords[1:2]) - 1000
endpt <- st_point(outcoords[3:4]) + 1000

xx <- startpt[1] + (linedist/3) * (cos(deg * pi / 180))
yy <- startpt[2] + (linedist/3) * (sin(deg * pi / 180))
xx2 <- startpt[1] - (linedist/3) * (cos(deg * pi / 180))
yy2 <- startpt[2] - (linedist/3) * (sin(deg * pi / 180))
startline <- st_cast(c(st_point(c(xx,yy)),  st_point(c(xx2,yy2))), "LINESTRING")

xx <- endpt[1] + (linedist/3) * (cos(deg * pi / 180))
yy <- endpt[2] + (linedist/3) * (sin(deg * pi / 180))
xx2 <- endpt[1] - (linedist/3) * (cos(deg * pi / 180))
yy2 <- endpt[2] - (linedist/3) * (sin(deg * pi / 180))
endline <- st_cast(c(st_point(c(xx,yy)),  st_point(c(xx2,yy2))), "LINESTRING")

# Find line between closest points on the two lines
perpline <- st_nearest_points(startline, endline)
# Find number of increments given the distance between the lines and a set
# increment (2km)
ninc <- round(st_length(perpline)/increment)
# Create ninc number of points at regular intervals along this line
incpts <- st_line_sample(perpline, ninc) %>%
  st_sfc() %>%
  st_set_crs(32632) %>%
  st_cast("POINT")

# Object to hold lines
inclines <- st_sfc() %>%  st_set_crs(32632)

# Loop over points
for (i in 1:length(incpts)){
  # Find x and y coords
  x <- st_coordinates(incpts[i])[1]
  y <- st_coordinates(incpts[i])[2]

  # Find coords at the specified distance from the point at deg degree angle

  xx <- x + linedist * (cos(deg * pi / 180))
  yy <- y + linedist * (sin(deg * pi / 180))
  xx2 <- x - linedist * (cos(deg * pi / 180))
  yy2 <- y - linedist * (sin(deg * pi / 180))

  pts <- st_sfc(st_multipoint(rbind(c(xx, yy), c(xx2, yy2)))) %>%
    st_set_crs(32632) %>% # WGS 84 / UTM 32N
    st_sf()

  inclines <- rbind(st_cast(pts, "LINESTRING"), inclines)

  # assign(paste0("increment", i), st_cast(pts, to = 'LINESTRING'))
}

plot(inclines)
# plot(c(startpt, endpt), add = TRUE)
plot(st_centroid(sites), pch = 19, cex = 0.5, colour = "black", add = TRUE)

incpolys <- st_sfc() %>%  st_set_crs(32632)

for(i in 1:(nrow(inclines) - 1)){
  s1 <- st_cast(inclines[i,], "POINT")[1,]
  e1 <- st_cast(inclines[i,], "POINT")[2,]

  s2 <- st_cast(inclines[i + 1,], "POINT")[1,]
  e2 <- st_cast(inclines[i + 1,], "POINT")[2,]

  st_cast(rbind(inclines[i,], inclines[i+1,]), "LINESTRING")

  pol <- rbind(s1, s2, e2, e1) %>%
    summarise(do_union = FALSE) %>%
    st_cast("POLYGON")

  incpolys <- rbind(pol, incpolys)
}

#  Find centroids of polygons
inccents <- st_centroid(incpolys)
inccents$id <- seq(1, nrow(inccents))

incpolys$disp <- NA
incpolys$id <- seq(1, nrow(incpolys))

# Interpolate displacement curves to centroids, assigning the curve
# to the polygon feature
for(i in 1:nrow(incpolys)){
  incpolys[i,]$disp <- list(shoredate::interpolate_curve(target = inccents[i,],
                                                         cal_reso = 5))
}

# Assign weight to the polygon features based on the distribution of sites
incpolys$dens <- lengths(st_intersects(incpolys, sites)) /
  sum(lengths(st_intersects(incpolys, sites)))

save(incpolys,
     file = here("analysis/data/derived_data/displacement_polygons.RData"))

# Retrieve displacement curves and weights for plotting
displist <- list()
for(i in 1:nrow(incpolys)){
  tmp_curve <- incpolys[i,]$disp[[1]][[1]]
  tmp_curve$dens <- incpolys[i,]$dens
  tmp_curve$id <- incpolys[i,]$id
  displist[[i]] <- tmp_curve
}

dispcurves <- do.call(rbind, displist)

displt <- ggplot(as.data.frame(dispcurves)) +
  ggplot2::geom_ribbon(ggplot2::aes(x = bce,
                                    ymin = lowerelev,
                                    ymax = upperelev,
                                    group = id,
                                    alpha = dens)) +
  labs(title = "B", x = "BCE/CE", y = "Meters above present sea-level") +
  scale_x_continuous(expand = c(0,0), limits = c(-9469, 2000)) +
  theme_bw() +
  theme(legend.position = "none")

bboxpolys <- st_bbox(sites)
bboxpolys[1] <- bboxpolys[1] - 900
bboxpolys[3] <- bboxpolys[3] + 900
bboxpolys[2] <- bboxpolys[2] - 900
bboxpolys[4] <- bboxpolys[4] + 900
bboxpolyspoly <- st_as_sf(st_as_sfc(bboxpolys))


mplt <- ggplot() +
  geom_sf(data = incpolys, aes(fill = dens), colour = NA, alpha = 0.7) +
  geom_sf(data = bmap, fill = NA, colour = "black") +
  geom_sf(data = st_centroid(sites), size = 0.5, colour = "black") +
  # ggsn::scalebar(data = surveyed, dist = 20, dist_unit = "km",
  #                transform = FALSE, st.size = 4, height = 0.02,
  #                border.size = 0.1, st.dist = 0.03,
  #                anchor = c(x = anc[2] - 15500, y = anc[1]) + 8000) +
  scale_fill_gradient(low = "white", high = "black", name = "Site density") +
  ggtitle("A") +
  coord_sf(xlim = c(bboxpolys[1], bboxpolys[3]),
           ylim = c(bboxpolys[2], bboxpolys[4]),
           expand = FALSE) +
  theme_bw() + theme(axis.title=element_blank(),
                     axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     rect = element_rect(),
                     axis.ticks = element_blank(),
                     panel.grid.major = element_blank(),
                     legend.position = "left")

mplt + displt

ggsave(here::here("analysis/figures/incpolys.png"),
       units = "px", width = 2650, height = 1600)
