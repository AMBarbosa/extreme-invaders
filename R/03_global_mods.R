setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("Rsc00_packages.R")


# IMPORT WORLDCLIM VARIABLES ####

download.file("https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_10m_bio.zip", destfile = "wclim.zip")
unzip("wclim.zip", exdir = "wclim")
wclim_stack <- raster::stack(list.files("wclim", pattern = "\\.tif$", full.names = TRUE))
plot(wclim_stack[[1]])
unlink("wclim.zip")


# IMPORT AND CLEAN GLOBAL OCCURRENCE DATA ####

data_dir <- "../data/GBIF_global"
temp_dir <- "tmp"

spp <- tools::file_path_sans_ext(list.files(data_dir, pattern = ".zip"))
spp

spatial_uncert_limit <- (15000 * sqrt(2)) / 2  # distance from centroid to corner (half the diagonal) of a square with 15-km side

rast <- wclim_stack[[1]]

dir.create("../data/GBIF_global/thinned")


for (spc in spp) {
  cat(spc, "\n")

  data_zip <- paste0(data_dir, "/", spc, ".zip")
  data_csv <- unzip(data_zip, list = TRUE)[1, 1]
  unzip(data_zip, exdir = temp_dir)
  gbif <- data.table::fread(paste(temp_dir, data_csv, sep = "/"), header = TRUE, sep = "\t", quote = "", select = c("species", "occurrenceStatus", "decimalLongitude", "decimalLatitude", "coordinateUncertaintyInMeters"))

  # remove absence records:
  gbif_pres <- subset(gbif, occurrenceStatus == "PRESENT")

  # remove wrong coordinates:
  gbif_clean <- coord_incomplete(coord_imprecise(coord_impossible(coord_unlikely(gbif_pres))))
  cat("Number of GBIF data points:                  ", nrow(gbif), "\n")
  cat("After coordinate cleaning with 'scrubr':     ", nrow(gbif_clean), "\n")

  # remove records exceding spatial uncertainty limit:
  gbif_smalluncert <- scrubr::coord_uncertain(gbif_clean, coorduncertainityLimit = spatial_uncert_limit)
  cat("After coordinate uncertainty filtering:      ", nrow(gbif_smalluncert), "\n")

  # extract to template raster resolution:
  gbif_sf <- st_as_sf(gbif_smalluncert, coords = c("decimalLongitude", "decimalLatitude"))
  st_crs(gbif_sf) <- 4326
  gbif_extract <- raster::extract(rast, gbif_sf, cellnumbers = TRUE, df = TRUE)

  # remove points on NA pixels:
  gbif_extract_noNA <- na.omit(gbif_extract)
  cat("After removing points on NA raster pixels:   ", nrow(gbif_extract_noNA), "\n")

  # get centroid of pixel where each point falls:
  gbif_pixcentroids <- raster::xyFromCell(rast[[1]], gbif_extract_noNA$cells)
  gbif_pixcentroids <- data.frame(cell = gbif_extract_noNA$cells, gbif_pixcentroids)

  # remove duplicated pixel centroids:
  gbif_pixcentroids_noduplicates <- unique(gbif_pixcentroids)
  cat("After removing duplicates within same pixel: ", nrow(gbif_pixcentroids_noduplicates), "\n\n")

  # export cleaned occurrence data to .csv:
  write.csv(gbif_pixcentroids_noduplicates, paste0("../data/GBIF_global/thinned/", spc, ".csv"), row.names = FALSE)

  # map occurrence data (takes long for species with many presences):
  maps::map("world", col = "grey")
  title(spc)
  points(gbif[ , c("decimalLongitude", "decimalLatitude")], pch = 20, col = "red", cex = 0.2)
  points(gbif_pixcentroids_noduplicates[ , c("x", "y")], pch = 20, col = "darkblue", cex = 0.2)
  legend("bottom", legend = c("raw", "clean"), pch = 20, col = c("red", "darkblue"), bty = "n")
}


# GENERATE GLOBAL PSEUDO-ABSENCES ####

dir.create("../data/GBIF_global/pseudoabs")

seed <- 0
for (spc in spp) {
  message(spc)
  pres <- read.csv(paste0("../data/GBIF_global/thinned/", spc, ".csv"))
  pres_sf <- st_as_sf(pres, coords = c("x", "y"))
  buff <- st_union(st_buffer(pres_sf, 2))
  buff_rast <- mask(crop(rast, as_Spatial(buff)), as_Spatial(buff))
  buff_rast_np <- mask(buff_rast, pres_sf, inverse = TRUE)
  plot(buff_rast_np, main = spc)
  seed <- seed + 1
  set.seed(seed)
  pseudoabs <- sampleRandom(buff_rast_np, size = nrow(pres), cells = TRUE, xy = TRUE)
  points(pseudoabs[ , c("x", "y")], pch = ".")
  write.csv(pseudoabs[ , 1:3], paste0("../data/GBIF_global/pseudoabs/", spc, ".csv"), row.names = FALSE)
}


# COMPUTE GLOBAL MODELS ####

mods <- prob <- fav <- vector("list", length(spp))
names(mods) <- names(prob) <- names(fav) <- spp


# create raster with variables only on study area for prediction:
utm10 <- st_read("../maps/utm10.shp")
plot(utm10[ , 1])
utm10 <- st_transform(utm10, st_crs(wclim_stack))
wclim_ib <- mask(crop(wclim_stack, utm10, snap = "out"), utm10)
plot(wclim_ib[[1]])

for (spc in spp) {
  message(spc)
  pres <- read.csv(paste0("../data/GBIF_global/thinned/", spc, ".csv"))
  pseudoabs <- read.csv(paste0("../data/GBIF_global/pseudoabs/", spc, ".csv"))
  pres_extr <- raster::extract(wclim_stack, pres[ , c("x", "y")], cellnumbers = TRUE)
  pseudoabs_extr <- raster::extract(wclim_stack, pseudoabs[ , c("x", "y")], cellnumbers = TRUE)
  dat <- base::rbind(pres_extr, pseudoabs_extr)
  dat <- data.frame(dat, xyFromCell(wclim_stack, dat[ , "cells"]))
  dat$presence <- c(rep(1, nrow(pres)), rep(0, nrow(pseudoabs)))
  seed <- seed + 1
  set.seed(seed)
  mods[[spc]] <- bart(x.train = dat[ , grep("wc", names(dat))], y.train = dat$presence, keeptrees = TRUE)
  prob[[spc]] <- predict(mods[[spc]], wclim_ib)
  fav[[spc]] <- Fav(pred = prob[[spc]], sample.preval = prevalence(dat$presence))
}


fav_stack <- raster::stack()
for (f in 1:length(fav))  fav_stack <- raster::stack(fav_stack, fav[[f]])
names(fav_stack) <- names(fav)

spectral0 <- rev(brewer.pal(10, "Spectral"))  # max was 11
spectral <- colorRampPalette(spectral0)(100)

plot(fav_stack, col = spectral, axes = FALSE, box = FALSE, zlim = c(0, 1))


