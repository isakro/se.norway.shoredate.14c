library(here)
library(supportR)

# Due to its size (4.1 GB), the DTM used for shoreline dating in
# 01_shoreline_spd.R is not provided directly with the repository of this paper,
# but is instead downloaded in tiled form from a GitHub repository associated
# with another paper at https://github.com/isakro/assessing.sealevel.dating.
# As an alternative it is also possible to download the raster data directly
# from the Norwegian Mapping Authority (https://hoydedata.no/LaserInnsyn2/),
# making sure the downloaded raster covers the sites provided in
# surveyed_sites.gpkg and excavated_sites.gpkg in the raw_data directory
# of this repository.

# The DTM tiles and file are stored in a local directory external from the R
# project. This has to be created, but can also be changed to another
# directory of choice. However, this will require editing of the file path to
# the dtm specified in 01_shoreline_spd.R

# Uncomment to create external directory if does not already exist
# dir.create(here("../external_data/data/"))

# Define destination of DTM
destdir <- here("../external_data/data/")

# List files in the GitHub repository folder
dtm_tiles <- github_ls("https://github.com/isakro/assessing.sealevel.dating",
                       folder = "analysis/data/raw_data/tiled_dtm")

# Download the files to the destination directory
for(i in 1:nrow(dtm_tiles)) { #
  print(paste0(i, "/", nrow(dtm_tiles)))
  download.file(url = paste0("https://raw.github.com/isakro/assessing.sealevel.dating/master/analysis/data/raw_data/tiled_dtm/",
                       dtm_tiles$name[i]), destfile = paste0(destdir, dtm_tiles$name[i]))
}

# Retrieve path to DTM tiles
dtmtiles <- list.files(destdir, pattern= ".*.tif$", full.names = TRUE)

# Load the tiles and merge
tiles <- sapply(dtmtiles, terra::rast)
names(tiles)[1:2] <- c('x', 'y')
dtm <- do.call(terra::merge, tiles)

# Delete the DTM tiles
files <- list.files(destdir, full.names = TRUE)
sapply(files, file.remove)

# Save merged raster in the external directory
terra::writeRaster(dtm, file.path(destdir, "dtm10.tif"), overwrite = TRUE)
