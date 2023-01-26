library(tidyverse)
library(sf)
library(terra)
library(here)
# library(spatstat)

sites <- st_read(here("analysis/data/derived_data/combined_sites.gpkg"))
bmap <- st_read(here("analysis/data/raw_data/naturalearth_countries.gpkg")) %>%
  st_transform(st_crs(sites))

# Specify distance of increments for which to interpolate displacement curve
increment <- 2000
# Degree direction of isobases
deg <- 327
# Arbritrary long distance for lines
linedist <- 100000

# Find bounding box of sites adding 1km
outcoords <- st_bbox(sites) + 1000
# Find bottom left and top right points
startpt <- st_point(outcoords[1:2])
endpt <- st_point(outcoords[3:4])

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
# increment (5km)
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
save(incpolys,
         file = here("analysis/data/derived_data/displacement_polygons.RData"))

ggplot() +
  geom_sf(data = incpolys, aes(col = id)) +
  # geom_sf(data = inccents, aes(col = id)) +
  geom_sf(data = st_centroid(sites), pch = 21, fill = "red") +
  theme_classic()

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
cols <- sample(color, nrow(incpolys))

plot(x = incpolys[1,]$disp[[1]][[1]]$bce,
     y = incpolys[1,]$disp[[1]][[1]]$lowerelev,
     type = "l", col = cols[1])
lines(incpolys[1,]$disp[[1]][[1]]$bce,
      incpolys[1,]$disp[[1]][[1]]$upperelev, col = cols[1])

for(i in 2:nrow(incpolys)){
  lines(x = incpolys[i,]$disp[[1]][[1]]$bce,
        y = incpolys[i,]$disp[[1]][[1]]$lowerelev,
        type = "l", col = cols[i])
  lines(incpolys[i,]$disp[[1]][[1]]$bce,
        incpolys[i,]$disp[[1]][[1]]$upperelev, col = cols[i])
}
