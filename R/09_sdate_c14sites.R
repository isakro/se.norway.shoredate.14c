library(here)
library(shoredate)
library(terra)
library(sf)
library(ggplot2)
library(dplyr)
library(ADMUR)

dtm <- rast(here("../external_data/dtm/dtm10.tif")) # see 00_dtm_prep.R
# Load RSPD models and probability distribution
load(here("analysis/data/derived_data/rcarbon_models_pd.RData"))
# Load excavated sites data
excavated <- st_zm(read_sf(here("analysis/data/raw_data/excavated_sites.gpkg")))
# Background map for plotting
bmap <- st_read(here("analysis/data/raw_data/naturalearth_norway.gpkg")) %>%
  st_transform(st_crs(excavated))

# Find names of sites with radiocarbon dates
c14_names <- tools::file_path_sans_ext(colnames(PD))

# Filter for sites with radiocarbon dates
c14_sites <- excavated %>%
  dplyr::filter(site_name %in% c14_names)

# Find the elevations of the sites with radiocarbon dates
c14_sites$elev <- terra::extract(dtm, vect(c14_sites), fun = mean)[, -1]

# Define here for use in Monte Carlo simulation as well
step_reso = 0.1

# Perform shoreline dating
c14sites_shoredates <- shoreline_date(c14_sites,
                                 elevation = c14_sites$elev,
                                 elev_reso = step_reso,
                                 cal_reso = 5,
                                 sparse = TRUE,
                                 verbose = FALSE)

# Save results
save(c14sites_shoredates,
     file = here("analysis/data/derived_data/sdates_c14.RData"))
load(here("analysis/data/derived_data/sdates_c14.RData"))

# Create ADMUR compatible PD and SPD
dates_dfs <- lapply(c14sites_shoredates, as.data.frame)
sdates <- do.call(rbind, dates_dfs)
sdates$site <- rep(1:length(dates_dfs),
                   each = nrow(sdates)/length(dates_dfs))

pd <- as.data.frame(sdates) %>%
  dplyr::filter(bce <= -2550 & bce >= -9445) %>%
  dplyr::mutate(bp = (bce * -1) + 1950) %>%
  dplyr::select(-bce) %>%
  dplyr::group_by(bp) %>%
  dplyr::filter(!is.na(probability)) %>%
  dplyr::group_by(site) %>%
  tidyr::pivot_wider(names_from = site, values_from = probability) %>%
  tibble::column_to_rownames("bp")
# pd <- sweep(pd, 2, colSums(pd),"/")

SPD <- as.data.frame(rowSums(pd))
SPD <- SPD/( sum(SPD) *5 )

c14_shoresum <- sum_shoredates(c14sites_shoredates)

# Exclude probabilities of 0 and use same cut-off as before.
c14_shoresumdf <- as.data.frame(c14_shoresum[[1]]) %>%
  filter(probability != 0)
c14_shoresumdf <- filter(c14_shoresumdf, bce <= -2550 & bce >= -9445)
# Normalise after cut-off
c14_shoresumdf$probability <- c14_shoresumdf$probability /
                              (sum(c14_shoresumdf$probability)* 5)

# Retrieve logistic model used for RSPD
log_c14 <- convertPars(pars = log$par,
                       years = minage:maxage,
                       type = 'logistic')

# Transform to BCE/CE
log_c14$year <- (log_c14$year - 1950) * -1


# Plot the new SSPD with the logistic model from the RSPD for inspection
ggplot() +
  geom_bar(data = c14_shoresumdf, aes(x = bce, y = probability),
           stat = "identity",
           col = "grey", fill = "grey") +
  scale_x_continuous(breaks = seq(-10000, -2550, 1000)) +
  geom_line(data = log_c14, aes(x = year, y = pdf)) +
  labs(title = "Summed probability distribution of \nshoreline dated sites with radiocarbon dates",
       x = "BCE", y = "Summed probability") +
  theme_bw()

#### Spatial distribution of radiocarbon dated sites ####

# Identify the spatial density distribution of sites with radiocarbon dates
# in the landscape for Monte Carlo simulation

# Specify distance of increments for which to interpolate displacement curve
increment <- 2000
# Degree direction of isobases
deg <- 327
# Arbritrary long distance for lines
linedist <- 100000

# Find bounding box of sites
outcoords <- st_bbox(c14_sites)
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
plot(st_centroid(c14_sites), pch = 19, cex = 0.5, colour = "black", add = TRUE)

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

incpolys$dens <- lengths(st_intersects(incpolys, c14_sites)) /
  sum(lengths(st_intersects(incpolys, c14_sites)))

# Store the object with a new name to separate it from the one
# saved in 02_displacement_prep.R
incpolys_c14 <- incpolys

save(incpolys_c14,
     file = here("analysis/data/derived_data/displacement_polygons_c14.RData"))
load(here("analysis/data/derived_data/displacement_polygons_c14.RData"))

# # Create plot
# displist <- list()
# for(i in 1:nrow(incpolys_c14)){
#   tmp_curve <- incpolys_c14[i,]$disp[[1]]
#   tmp_curve$dens <- incpolys_c14[i,]$dens
#   tmp_curve$id <- incpolys_c14[i,]$id
#   displist[[i]] <- tmp_curve
# }
#
# dispcurves <- do.call(rbind, displist)
#
# ggplot(as.data.frame(dispcurves)) +
#   ggplot2::geom_ribbon(ggplot2::aes(x = bce,
#                                     ymin = lowerelev,
#                                     ymax = upperelev,
#                                     group = id,
#                                     alpha = dens)) +
#   labs(title = "B", x = "BCE/CE", y = "Meters above present sea-level") +
#   scale_x_continuous(expand = c(0,0), limits = c(-9469, 2000)) +
#   theme_bw() +
#   theme(legend.position = "none")


bboxpolys <- st_bbox(c14_sites)
bboxpolys[1] <- bboxpolys[1] - 900
bboxpolys[3] <- bboxpolys[3] + 900
bboxpolys[2] <- bboxpolys[2] - 900
bboxpolys[4] <- bboxpolys[4] + 900
bboxpolyspoly <- st_as_sf(st_as_sfc(bboxpolys))


ggplot() +
  geom_sf(data = incpolys_c14, aes(fill = dens), colour = NA, alpha = 0.7) +
  geom_sf(data = bmap, fill = NA, colour = "black") +
  geom_sf(data = st_centroid(c14_sites), size = 0.5, colour = "black") +
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


#### Monte Carlo simulation ####

# The logistic model does not cover the entire date range for the SSPD.
# Insert year with a probability of 0 for these years to the model
missing <- as.data.frame(cbind(setdiff(c14_shoresumdf$bce, log_c14$year),
                 rep(0, length(setdiff(c14_shoresumdf$bce, log_c14$year)))))
names(missing) <- names(log_c14)

log_c14 <- rbind(log_c14, missing)

log_c14$model <- "Logistic"
names(log_c14) <- c("bce", "prob_dens", "model")

ssize <- c14_shoresum[[2]]

# Set number of simulations
nsim <- 10

# Create data frame to hold simulated dates
random_dates <- data.frame(matrix(nrow = ssize*nsim, ncol = 3))
names(random_dates) <- c("sample", "simn", "displacement_curve")

# Draw random samples from the model and assign displacement curve
# based on the distribution of sites.
random_dates$sample = sample(log_c14$bce, replace = TRUE,
                             size = ssize*nsim, prob = log_c14$prob_dens)
random_dates$simn <- rep(1:nsim, each = nrow(random_dates)/nsim)
random_dates$displacement_curve <- incpolys_c14[sample(incpolys_c14$id,
                                                   replace = TRUE,
                                                   size = ssize*nsim,
                                                   prob = incpolys_c14$dens),]$disp

# Reverse sampled date using reverse_shoredate from 03_functions.R
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
    target_curve = random_dates$displacement_curve[[i]],
    sparse = TRUE))
  simdates[[i]]$simn <- random_dates$simn[i]

  if(i %% ssize == 0){
    # Save data externally and overwrite simdates for memory purposes
    tmp <- bind_rows(simdates) %>% group_by(bce, simn) %>%
      filter(!is.na(probability)) %>%
      summarise(prob_sum = sum(probability))

    save(tmp,
         file = paste0(here::here("../external_data/shorespd/log_c14_shore/"),
                       i, "shore.rds"))

    simdates <- vector("list", ssize)
  }
}
end_time <- Sys.time()
end_time - start_time


resfiles <- list.files(here::here("../external_data/shorespd/log_c14_shore/"),
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

mod <- approx(x = (log_c14$bce * -1) + 1950, y = log_c14$prob_dens,
              xout = as.numeric(rownames(SPD)),
              ties = 'ordered', rule = 2)$y

# Model summary (can also be plotted with ADMUR)
log_c14_summary <- simulation_summary(SPD, simresults,
                                  c(-2550, -9445), mod, ssize)

plot_mc(log_c14_summary)

save(log_c14_summary,
     file = here("analysis/data/derived_data/log_mc_shorec14_summary.rda"))

