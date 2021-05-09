# R scripts for the article:
# Baquero R.A., Barbosa A.M., Ayllón D., Guerra C., Sánchez E., Araújo M.B. & Nicola G.G. Potential distributions of invasive vertebrates in the Iberian Peninsula under projected changes in extreme climate events.


source("00_packages.R")


# IMPORT MAP OF THE POLYGON GRID TO WHICH THE VARIABLES REFER ####

utm10 <- readOGR("maps/utm10_ib.gpkg", "utm10_ib")
plot(utm10)


# IMPORT TABLE WITH THE VALUES OF THE VARIABLES ON THIS GRID ####

dat <- read.csv("data/variables_Baquero_et_al.csv")
head(dat)
names(dat)
names(dat) <- gsub(pattern = "\\.", replacement = "_", x = names(dat))
names(dat)

var_cols <- 3:ncol(dat)

names(dat)[var_cols]

species_names <- c("Amandava amandava", "Myiopsitta monachus", "Psittacula krameri", "Estrilda astrild", "Neovison vison", "Trachemys scripta")
spp <- spCodes(species_names, 1, 1)
species <- data.frame(spp = spp, species = species_names)
species


# add centroid coordinates to map and dat:

utm10@data$lonUTM <- coordinates(utm10)[ , 1]
utm10@data$latUTM <- coordinates(utm10)[ , 2]

head(utm10@data)
names(dat)
nrow(dat)  # 5665
names(dat)
nrow(utm10@data)  # 6215
dat <- merge(dat, utm10@data, by.x = "UTM_10X10", by.y = "utm10", all = TRUE)

nrow(dat)  # 6215
nrow(unique(dat))  # now 6215

names(dat)[1] <- "utm10"

names(dat)
coord_cols <- 93:94

choroLayer(spdf = utm10, df = dat, spdfid = "utm10", dfid = "utm10", var = "lonUTM", border = NA)  # OK
choroLayer(spdf = utm10, df = dat, spdfid = "utm10", dfid = "utm10", var = "latUTM", border = NA)  # OK


# ADD SPECIES OCCURRENCE DATA ####

# added data from national (Portuguese and Spanish) atlases
# of distributions of birds, mammals, and amphibians & reptiles
# (available from their sources upon acceptance of authors' conditions)
# added more data from a range of additional sources
# (see Appendix S1 of the article)
# with this script, we will add data from two sources:
# GBIF (https://www.gbif.org)
# and Ascensao et al. (https://neobiota.pensoft.net/article/55597/)


# ADD OCCURRENCE DATA FROM GBIF ####

gbif <- read.csv("data/GBIFdownload_2020-12-01.csv", sep = "\t", as.is = TRUE)  # citation for these data: GBIF.org (2020) GBIF Occurrence Download, 1 December 2020. DOI: https://doi.org/10.15468/dl.kw2dcy

head(gbif)
names(gbif)
sort(unique(gbif$species))

sort(unique(gbif$occurrenceStatus), na.last = TRUE)
sort(unique(gbif$individualCount), na.last = TRUE)
sort(unique(gbif$coordinateUncertaintyInMeters), na.last = TRUE)  # up to 895041
sort(unique(gbif$coordinatePrecision), na.last = TRUE)  # 1e-05  NA
nrow(gbif)  # 60994
gbif <- subset(gbif, occurrenceStatus != "ABSENT")
nrow(gbif)  # 60994
gbif <- subset(gbif, is.na(individualCount) | individualCount > 0)
nrow(gbif)  # 60994

sort(unique(gbif$coordinateUncertaintyInMeters))

names(gbif)
gbif <- gbif[ , c("species", "decimalLongitude", "decimalLatitude", "occurrenceStatus", "coordinateUncertaintyInMeters")]
head(gbif)
unique(gbif$occurrenceStatus)  # "PRESENT" ""


# CLEAN GBIF DATA WITH 'scrubr' PACKAGE:

nrow(gbif)  # 60994
gbif <- coord_incomplete(coord_imprecise(coord_impossible(coord_unlikely(gbif, lat = "decimalLatitude", lon = "decimalLongitude", drop = TRUE))))
nrow(gbif)  # 60801

gbif <- coord_uncertain(gbif, coorduncertainityLimit = 7000)
nrow(gbif)  # 59382


# map GBIF data and extract to utm10:

utm10_wgs <- spTransform(utm10, CRS = "+proj=longlat +ellps=WGS84")

plot(utm10_wgs, border = "grey")
points(gbif[ , c("decimalLongitude", "decimalLatitude")], pch = 20, col = as.factor(gbif$species))

gbif_pts <- gbif
coordinates(gbif_pts) <- ~ decimalLongitude + decimalLatitude
crs(gbif_pts) <- "+proj=longlat +ellps=WGS84"

species <- unique(gbif$species)
species

gbif_utm10 <- vector("list", length(species))
names(gbif_utm10) <- species
gbif_utm10

for (s in species) {
  gbif_utm10[[s]] <- over(gbif_pts[gbif_pts$species == s, ], utm10_wgs)
}

lapply(gbif_utm10, head)
sapply(gbif_utm10, nrow)
sapply(gbif_utm10, function(x) nrow(unique(x)))

gbif_utm10 <- lapply(gbif_utm10, function(x) unique(x$utm10))
sapply(gbif_utm10, length)
gbif_utm10 <- lapply(gbif_utm10, function(x) x[!is.na(x)])
sapply(gbif_utm10, length)

names(gbif_utm10)
match(names(gbif_utm10), species)  # I changed "species" object at some point
match(names(gbif_utm10), species_names)
gbif_utm10 <- gbif_utm10[match(names(gbif_utm10), species_names)]
names(gbif_utm10)
names(gbif_utm10) <- spCodes(names(gbif_utm10), 1, 1)

nrow(dat)  # 6215
dat_gbif <- data.frame(utm10 = dat$utm10)
for (s in spp) {
  dat_gbif <- data.frame(dat_gbif, rep(0, nrow(dat_gbif)))
  names(dat_gbif)[ncol(dat_gbif)] <- s
  dat_gbif[ , s][dat_gbif[ , "utm10"] %in% gbif_utm10[[s]]] <- 1
}

head(dat_gbif)
names(utm10)
names(dat_gbif)

par(mfrow = arrangePlots(length(spp)), mar = c(0, 0, 1, 0))
for(s in spp) {
  choroLayer(spdf = utm10, df = dat_gbif, spdfid = "utm10", dfid = "utm10", var = s, border = NA)
}

names(dat)
sapply(dat_gbif[ , -1], sum, na.rm = TRUE)


# ADD SPECIES OCCURRENCE DATA TO THE MAIN DATA TABLE:

names(dat)
names(dat_gbif)
dat <- data.frame(dat_gbif, dat[ , -1])

species <- data.frame(spp = spp, species = species_names)
species

#pdf("images/Fig_1_occurrence.pdf", width = 7, height = 8)
par(mfrow = arrangePlots(length(spp)), mar = c(0, 0, 2, 0))
for(s in spp) {
  choroLayer(spdf = utm10, df = dat, spdfid = "utm10", dfid = "utm10", var = s, border = NA, legend.pos = NA, col = c("grey86", "black"))
  title(species[species$spp == s, "species"], font.main = 3, cex.main = 2.5)
}
#dev.off()


# ADD OCCURRENCE DATA FROM ASCENSAO ET AL. (2021) ####

asc <- read.csv("https://zenodo.org/record/4018706/files/Data_Ascensa%CC%83oEtAl_Neobiota.csv?download=1")
# citation for these data:
# Ascensao F., D’Amico M., Martins R.C., Rebelo R., Barbosa A.M., Bencatel J., Barrientos R., Abellan P., Tella J.L., Cardador L., Anadon J.D., Carrete M., Murgui E., Fernandes P., Santos S.M., Mira A., Mathias M.L., Tiago P., Casabella E., Reino L., Paulo O.S., Pereira H.M. & Capinha C. (2021) Distribution of alien tetrapods in the Iberian Peninsula. NeoBiota, 64: 1-21
str(asc)
head(asc)
nrow(asc)  # 159676
names(asc)
sort(unique(asc$Species))
species
asc <- subset(asc, Species %in% species$species)
nrow(asc)  # 102324

asc$utm10 <- substr(asc$UTM, 4, 7)

head(asc)
head(dat)
length(unique(dat$utm10)) == nrow(dat)  # TRUE, so all UTM cells are in the table

# create list with Ascensao et al.'s presences for each species:

asc_pres <- vector("list", nrow(species))
names(asc_pres) <- unique(species$species)
for (spc in names(asc_pres)) {
  asc_pres[[spc]] <- subset(asc, Species == spc)$utm10
}
asc_pres

names(asc_pres)
names(asc_pres) <- spCodes(names(asc_pres), 1, 1)

sapply(asc_pres, length)
sapply(asc_pres, function(x) length(unique(x)))  # different, so there were repeats (probably same cell different sources)
sapply(dat[species$spp], sum)
# mostly smaller numbers, so Ascensao et al. provided several additional presences

for (i in species$spp) {
  dat[ , i][as.character(dat$utm10) %in% asc_pres[[i]]] <- 1
}

sapply(dat[as.character(species$spp)], sum)  # ok - numbers visibly increased

for (i in spp) {
  choroLayer(spdf = utm10, df = dat, spdfid = "utm10", dfid = "utm10", var = i, border = NA, legend.pos = NA); title(i)
}


# CLIMATE MODELS ####

names(dat)
var_cols <- 9:98
names(dat)[var_cols]

# simplify composed variable names:
names(dat)[var_cols] <- gsub("Hum_relat", "Humrel", names(dat)[var_cols])
names(dat)[var_cols]

climods <- unique(sapply(strsplit(names(dat)[var_cols], "_"), `[`, 2))
climods
# "REM" "WRA" "WRB" "PRO" "MM5"

vars <- names(dat)[grep("REM", names(dat))]
vars
vars <- unique(sapply(strsplit(vars, "_"), `[`, 1))
vars

vars_present <- names(dat)[endsWith(names(dat), "1971_2000")]
vars_future <- names(dat)[endsWith(names(dat), "2021_2050")]
vars_present
vars_future


# MAP DATA ####

par(mfrow = c(2, 2))
choroLayer(spdf = utm10, df = dat, spdfid = "utm10", dfid = "utm10", var = "lonUTM", border = NA)
choroLayer(spdf = utm10, df = dat, spdfid = "utm10", dfid = "utm10", var = "latUTM", border = NA)
choroLayer(spdf = utm10, df = dat, spdfid = "utm10", dfid = "utm10", var = names(dat)[var_cols][1], border = NA)
choroLayer(spdf = utm10, df = dat, spdfid = "utm10", dfid = "utm10", var = names(dat)[var_cols][90], border = NA)

# map species:
par(mfrow = arrangePlots(length(spp)))
for (i in spp) {
  choroLayer(spdf = utm10, df = dat, spdfid = "utm10", dfid = "utm10", var = i, border = NA, legend.pos = NA)
  title(i)
}


# make species presence NA where variables are NA:

unique(dat$COUNTRY)  # NA "SPAIN" "PORTUGAL" "SPAIN-PORTUGAL"
dat$COUNTRY[is.na(dat$COUNTRY)] <- "SPAIN"

for (s in spp) {
  dat[which(dat$COUNTRY != "PORTUGAL" & is.na(dat$Tmax30_REM_1971_2000)), s] <- NA
  choroLayer(spdf = utm10, df = dat, spdfid = "utm10", dfid = "utm10", var = s, border = NA, legend.pos = NA, col = c("grey90", "black"))
  title(species[species$spp == s, "species"], font = 4)
}


nrow(dat)  # 6215
nrow(unique(dat))  # 6215
length(unique(dat$utm10))  # 6215
length(utm10)  # 6215
nrow(na.omit(dat))  # 5665
unique(sapply(dat[ , var_cols], function(x) sum(is.finite(x))))  # 5665


# PUT CLIMODS AT BEGINNING OF VAR NAMES INSTEAD OF END (as they were in the first version of the dataset when I wrote the script)
splits <- strsplit(names(dat)[var_cols], "_")
splits
sapply(splits, `[`, 1)
sapply(splits, `[`, 2)
names(dat)[var_cols] <- paste(sapply(splits, `[`, 2), sapply(splits, `[`, 1), sapply(splits, `[`, 3), sapply(splits, `[`, 4), sep = "_")
names(dat)

