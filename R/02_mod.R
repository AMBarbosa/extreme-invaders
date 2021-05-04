# R scripts for the article:
# Baquero R.A., Barbosa A.M., Ayllón D., Guerra C., Sánchez E., Araújo M.B. & Nicola G.G. Potential distributions of invasive vertebrates in the Iberian Peninsula under projected changes in extreme climate events.


source("00_packages.R")
source("01_data.R")


# ORGANIZE DATA ####

dat$rownum <- 1:nrow(dat)

# get UTM10 x and y numbers to fill UTM zone unions later:
dat$utm10x <- as.integer(substr(dat$utm10, 3, 3))
dat$utm10y <- as.integer(substr(dat$utm10, 4, 4))

head(dat)
nrow(dat)  # 6215
names(dat)
ncol(dat)  # 103
reordered_cols <- c(101, 1, 102:103, 99:100, 8, 2:7, 9:98)
length(reordered_cols) == ncol(dat)  # TRUE - ok
dat <- dat[ , reordered_cols]
head(dat)
names(dat)

coord_cols <- 5:6
sp_cols <- 8:13
var_cols <- 14:103
names(dat)[coord_cols]
names(dat)[sp_cols]
names(dat)[var_cols]


# CALCULATE SPATIAL TRENDS (TSA) ####

tsa <- multTSA(data = dat, sp.cols = sp_cols, coord.cols = coord_cols, id.col = "utm10", type = "F")
head(tsa)
names(tsa)[-1] <- paste0(names(tsa)[-1], "_fav")
head(tsa)


# plot the spatial trends:

spectral0 <- rev(brewer.pal(10, "Spectral"))  # max was 11
spectral <- colorRampPalette(spectral0)(100)

greys <- brewer.pal(9, "Greys")
greys100 <- colorRampPalette(greys[2:length(greys)])(100)

brks <- seq(0, 1, 0.01)

brks_leg <- seq(0, 1, 0.1)
spectral_leg <- colorRampPalette(spectral0)(length(brks_leg))

#jpeg("images/Fig_S2_spatial_trends.jpg", width = 800, height = 900)
par(mar = c(0, 0, 1.5, 0), mfrow = modEvA::arrangePlots(length(spp)))
for (i in 2:ncol(tsa)) {
  spc <- substr(names(tsa)[i], 1, 2)
  choroLayer(spdf = utm10, df = tsa, spdfid = "utm10", dfid = "utm10", var = names(tsa)[i], border = NA, legend.pos = NA, col = spectral, breaks = brks)
  #points(dat[dat[ , spc] == 1, c("lonUTM", "latUTM")], pch = 20, cex = 0.2)
  title(species[species$spp == spc, "species"], font.main = 3, cex.main = 2.5)
}
legendChoro(pos = "bottomright", border = NA, breaks = brks_leg, col = spectral_leg, symbol = "box", title.txt = "", nodata = FALSE, values.cex = 1.3)
#dev.off()

names(dat)
nrow(dat)  # 6215
dat <- merge(dat, tsa)
nrow(dat)  # 6215
nrow(unique(dat))  # 6215
tail(dat)

par(mfrow = c(1, 1))
choroLayer(spdf = utm10, df = dat, spdfid = "utm10", dfid = "utm10", var = names(dat)[var_cols][1], border = NA)


# SEPARATE DATA INTO PRESENT AND FUTURE ####
# (so that variables don't have different names when projecting)

vars_present
vars_future

names(dat)

dat_pres <- dat[ , c(1:13, grep("1971_2000", names(dat)), 104:ncol(dat))]
dat_fut <- dat[ , c(1:13, grep("2021_2050", names(dat)), 104:ncol(dat))]

names(dat_pres)
names(dat_fut)
ncol(dat)  # 109
ncol(dat_pres)  # 64
ncol(dat_fut)  # 64

names(dat_pres)
var_cols <- 14:58  # now with only vars from the present
names(dat_pres)[var_cols]
names(dat_fut)[var_cols]


# remove the years from the column (var) names:
names(dat_pres)[var_cols] <- substr(names(dat_pres)[var_cols], 1, nchar(names(dat_pres)[var_cols]) - 10)
names(dat_fut)[var_cols] <- substr(names(dat_fut)[var_cols], 1, nchar(names(dat_fut)[var_cols]) - 10)
all.equal(names(dat_pres), names(dat_fut))  # TRUE - ok


# DIVIDE STUDY AREA FOR 10-FOLD CROSS-VALIDATION ####

# restrict 'dat' to rows with variables:

finite_rows <- is.finite(dat_pres$REM_Tmax30)
finite_rows2 <- is.finite(dat_fut$REM_Tmax30)
all.equal(finite_rows, finite_rows2)  # TRUE
rm(finite_rows2)

dat_pres <- dat_pres[finite_rows, ]
dat_fut <- dat_fut[finite_rows, ]

nrow(dat_pres)  # 5665
nrow(dat_fut)  # 5665
head(dat_pres)

nrow(dat)  # 6215


# spatial block folds using 'blockCV' package:

names(dat_pres)
utm_centroids <- dat_pres[ , c("lonUTM", "latUTM")]
coordinates(utm_centroids) <- utm_centroids
crs(utm_centroids) <- crs(utm10)
plot(utm_centroids)

set.seed(2020)
blocks <- spatialBlock(speciesData = utm_centroids, species = NULL, rasterLayer = NULL, theRange = 100000, k = 10L, selection = "systematic", iteration = 100L, numLimit = 0L, border = utm10, showBlocks = TRUE, biomod2Format = TRUE, xOffset = 0, yOffset = 0, progress = TRUE, verbose = TRUE)
# previously working fine, now Error in st_geos_binop("intersects", x, y, sparse = sparse, prepared = prepared, : st_crs(x) == st_crs(y) is not TRUE

#saveRDS(blocks, "cvblocks.rds")  # saved from previous session where this had worked
blocks <- readRDS("R/cvblocks.rds")

blocks$folds
blocks$foldID

#pdf("images/Fig_S1_CVblocks.pdf")
#plot(blocks)
#dev.off()

# add folds to data table:
dat_pres$cvfolds <- blocks$foldID

cvfolds <- vector("list", length(blocks$folds))
for (i in 1:length(cvfolds)) {
  cvfolds[[i]] <- which(dat_pres$cvfolds == i)
}
cvfolds
sapply(cvfolds, length)  # 542 728 728 704 648 615 543 400 348 409

folds <- 1:length(cvfolds)


# COMPUTE THE MODELS WITH 4 ALGORITHMS: GLM, GAM, RF, BART ####

#source("predict_bart_df.R")  # I edited predict2.bart to work on data frames rather than raster layers
source("https://raw.githubusercontent.com/AMBarbosa/unpackaged/master/predict_bart_df")

nrow(dat_pres)  # 5665

preds_cv <- data.frame(dat_pres[ , c("utm10", "rownum", "lonUTM", "latUTM", spp, "cvfolds")])
nrow(preds_cv)  # 5665
head(preds_cv)

# scale the variables to get standardized (comparable) coefficients (variable importance):
names(dat_pres)
head(dat_pres)
names(dat_pres)[var_cols]
names(dat_fut)[var_cols]
dat_pres[ , var_cols] <- scale(dat_pres[ , var_cols])
dat_fut[ , var_cols] <- scale(dat_fut[ , var_cols])
head(dat_pres)
head(dat_fut)


# list to store variable importance:

meths <- c("GLM", "GAM", "RF", "BART")
vars_chosen_general <- c("Tmax30", "TmaxAvemin", "PrAvemax5d", "Pr1mm5d")

n_spp <- length(spp)
n_meths <- length(meths)
n_vars <- length(vars_chosen_general)

var_contribs <- data.frame(species = rep(spp, each = n_meths * n_vars),
                           method = rep(rep(meths, each = n_vars), n_spp),
                           variable = rep(vars_chosen_general, n_spp * n_meths))
var_contribs
for (climod in climods) {
  var_contribs[ , climod] <- NA
}
names(var_contribs)[4:ncol(var_contribs)] <- climods
head(var_contribs)

nrow(var_contribs)  # 96
n_vars * n_meths * n_spp  # 96 - ok

gc()


# loop to compute the models and calculate variable importance (MAY TAKE A FEW HOURS!):
start_time_init <- Sys.time()
n <- 0
for (spc in spp)  for (climod in climods) {

  # report progress:
  n <- n + 1
  start_time <- Sys.time()
  cat("Modelling set", n, "of", length(spp) * length(climods), ":", spc, climod, "...")

  vars_chosen <- paste(climod, vars_chosen_general, sep = "_")

  spc_fact <- as.factor(dat_pres[ , spc])

  form_glm <- as.formula(paste0(spc, "~", paste(vars_chosen, collapse = "+")))
  form_gam <- as.formula(paste0(spc, "~", paste0("s(", vars_chosen, ")", collapse = "+")))  # GAM with smoothing splines
  form_rf <- as.formula(paste("spc_fact", "~", paste(vars_chosen, collapse = "+")))


  for (fold in folds) {

    cat(" fold", fold, "...")

    dat_train <- dat_pres[-cvfolds[[fold]], ]
    dat_test <- dat_pres[cvfolds[[fold]], ]

    spc_fact <- as.factor(dat_train[ , spc])  # for RandomForest

    mod_glm <- glm(form_glm, family = binomial, data = dat_train)
    mod_gam <- gam(form_gam, family = binomial, data = dat_train)
    mod_rf <- randomForest(form_rf, data = dat_train, na.action = na.exclude)
    # try also weighting presences and absences equally (as random forests don't deal well with imbalanced data), although ultimately we need to use the originally imbalanced data so that predicted probabilities are commensurable:
    #mod_rf <- randomForest(form_rf, data = dat_train, na.action = na.exclude, sampsize = rep(sum(dat_train[ , spc], na.rm = TRUE), 2), ntree = 1000)  # 'sampsize' to avoid randomforest problems with imbalanced data: https://stats.stackexchange.com/questions/163251/creating-a-test-set-with-imbalanced-data/163567#163567
    mod_bart <- bart(y.train = dat_train[ , spc], x.train = dat_train[ , vars_chosen], keeptrees = TRUE, verbose = FALSE)

    foldname <- ifelse(fold < 10, paste0("fold0", fold), paste0("fold", fold))
    preds_cv[ , paste0(spc, "_", climod, "_GLM_", foldname)] <- predict(mod_glm, dat_pres, type = "response")
    preds_cv[ , paste0(spc, "_", climod, "_GAM_", foldname)] <- predict(mod_gam, dat_pres, type = "response")
    preds_cv[ , paste0(spc, "_", climod, "_RF_", foldname)] <- predict(mod_rf, dat_pres, type = "prob")[ , "1"]
    preds_cv[ , paste0(spc, "_", climod, "_BART_", foldname)] <- predict_bart_df(mod_bart, dat_pres)
  }  # end for fold


  # do also the models with all data (no folding):

  cat(" ALL ...")

  spc_fact <- as.factor(dat_pres[ , spc])

  mod_glm <- glm(form_glm, family = binomial, data = dat_pres)
  mod_gam <- gam(form_gam, family = binomial, data = dat_pres)
  mod_rf <- randomForest(form_rf, data = dat_train, na.action = na.exclude)
  # try also weighting presences and absences equally (as random forests don't deal well with imbalanced data), although ultimately we need to use the originally imbalanced data so that predicted probabilities are commensurable:
  #mod_rf <- randomForest(form_rf, data = dat_pres, na.action = na.exclude, sampsize = rep(sum(dat_train[ , spc], na.rm = TRUE), 2), ntree = 1000)  # 'sampsize' to avoid randomforest problems with imbalanced data: https://stats.stackexchange.com/questions/163251/creating-a-test-set-with-imbalanced-data/163567#163567
  mod_bart <- bart(y.train = dat_pres[ , spc], x.train = dat_pres[ , vars_chosen], keeptrees = TRUE, verbose = FALSE)

  var_contribs[var_contribs$species == spc & var_contribs$method == "GLM", climod] <- coefficients(mod_glm)[-1]  # excludes the intercept
  var_contribs[var_contribs$species == spc & var_contribs$method == "GAM", climod] <- coefficients(mod_gam)[-1]  # excludes the intercept
  var_contribs[var_contribs$species == spc & var_contribs$method == "RF", climod] <- importance(mod_rf)[, "MeanDecreaseGini"]
  var_contribs[var_contribs$species == spc & var_contribs$method == "BART", climod] <- varimp(mod_bart)[, "varimps"]

  preds_cv[ , paste0(spc, "_", climod, "_GLM_alldat")] <- predict(mod_glm, dat_pres, type = "response")
  preds_cv[ , paste0(spc, "_", climod, "_GAM_alldat")] <- predict(mod_gam, dat_pres, type = "response")
  preds_cv[ , paste0(spc, "_", climod, "_RF_alldat")] <- predict(mod_rf, dat_pres, type = "prob")[ , "1"]
  preds_cv[ , paste0(spc, "_", climod, "_BART_alldat")] <- predict_bart_df(mod_bart, dat_pres)


  # predict with models of all data to the future:

  cat(" FUTURE ...\n")
  preds_cv[ , paste0(spc, "_", climod, "_GLM_future")] <- predict(mod_glm, dat_fut, type = "response")
  preds_cv[ , paste0(spc, "_", climod, "_GAM_future")] <- predict(mod_gam, dat_fut, type = "response")
  preds_cv[ , paste0(spc, "_", climod, "_RF_future")] <- predict(mod_rf, dat_fut, type = "prob")[ , "1"]
  preds_cv[ , paste0(spc, "_", climod, "_BART_future")] <- predict_bart_df(mod_bart, dat_fut)

  print(head(preds_cv[ , (ncol(preds_cv) - 10) : ncol(preds_cv)]))
  gc()
  fuzzySim::timer(start_time)
  message("From the beginning:")
  fuzzySim::timer(start_time_init)
}  # end models spp x climods x folds


head(preds_cv[  , 1:30])
ncol(preds_cv)  # 1451
head(names(preds_cv), 20)

nrow(preds_cv)  # 5665
names(preds_cv)[1:20]
tail(names(preds_cv))  # "_alldat" and "_future"
min(sapply(preds_cv[ , 12:ncol(preds_cv)], min, na.rm = TRUE))  # 0
max(sapply(preds_cv[ , 12:ncol(preds_cv)], max, na.rm = TRUE))  # 1

names(dat_pres)
nrow(dat_pres)  # 5665
nrow(preds_cv)  # 5665

write.csv(preds_cv, "preds_cv.csv", row.names = FALSE)

par(mfrow = c(2, 2))
choroLayer(spdf = utm10, df = preds_cv, spdfid = "utm10", dfid = "utm10", var = "Mm_REM_BART_future", border = NA)
choroLayer(spdf = utm10, df = preds_cv, spdfid = "utm10", dfid = "utm10", var = "Mm_REM_GAM_future", border = NA)
choroLayer(spdf = utm10, df = preds_cv, spdfid = "utm10", dfid = "utm10", var = "Mm_REM_GLM_future", border = NA)
choroLayer(spdf = utm10, df = preds_cv, spdfid = "utm10", dfid = "utm10", var = "Mm_REM_RF_future", border = NA)


# get correlations among selected variables:
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
for (v in vars_chosen) { hist(dat_pres[ , v], main = v) }
sel_var_corrs <- cor(dat_pres[ , vars_chosen], method = "spearman")
sel_var_corrs
diag(sel_var_corrs) <- NA
max(abs(sel_var_corrs), na.rm = TRUE)


# REORDER PREDS:

names(preds_cv)[1:20]
preds_cv <- preds_cv[ , c(names(preds_cv)[1:11], sort(names(preds_cv)[12:ncol(preds_cv)]))]
names(preds_cv)[1:20]


# MAP SOME PREDICTIONS ####

par(mfrow = c(4, 3), mar = c(0, 0, 1, 0))

for (i in 12:23) {
  choroLayer(spdf = utm10, df = preds_cv, spdfid = "utm10", dfid = "utm10", var = names(preds_cv)[i], border = NA, col = spectral, breaks = brks, legend.pos = "n")
  title(names(preds_cv)[i])
}

# clean up:
rm(i, n, mod_bart, mod_gam, mod_glm, mod_rf, form_gam, form_glm, form_rf, dat_train, dat_test, start_time, start_time_init)
gc()


# EVALUATE FOLDS ####

names(preds_cv)
names(dat_pres)
nrow(preds_cv)  # 5665
nrow(na.omit(preds_cv)) # 5665

head(preds_cv[ , 1:20])

fold_names <- c(paste0(0, 1:9), "10")

spp
climods
meths
fold_names

# first, create an empty list which will receive the AUC and Miller slope values for each species, climod, method, and fold:

crossval <- vector("list", length(spp) * length(climods) * length(meths))

par(mfrow = c(5, 2), mar = c(2, 2, 1, 1), oma = c(0, 0, 1.5, 0))

element <- 0; for (spc in spp) for (climod in climods) for (meth in meths) {
  element <- element + 1
  crossval[[element]] <- as.data.frame(matrix(NA, nrow = length(fold_names), ncol = 3))
  names(crossval)[element] <- paste(spc, climod, meth, sep = "_")
  names(crossval[[element]]) <- c("AUC", "TSS", "Miller")
  fold_pred_names <- paste0(paste(spc, climod, meth, "fold", sep = "_"), fold_names)

  fold_preds <- preds_cv[ , fold_pred_names]
  presabs <- preds_cv[ , spc]

#  jpeg(paste0("images/AUC_Miller_folds/", names(crossval)[element], ".jpg"), width = 400, height = 600)
  par(mfrow = c(5, 4), mar = c(2.3, 2.1, 1, 1), oma = c(0, 0, 1.5, 0))
  for(f in 1:ncol(fold_preds)) {
    foldname <- sapply(strsplit(names(fold_preds)[f], "_"), `[`, 4)
    crossval[[element]][f, "AUC"] <- AUC(obs = presabs, pred = fold_preds[,f], simplif = TRUE, plot = TRUE, main = foldname)
    crossval[[element]][f, "TSS"] <- threshMeasures(obs = presabs, pred = fold_preds[,f], measures = "TSS", simplif = TRUE, thresh = "preval", plot = FALSE, standardize = FALSE)
    crossval[[element]][f, "Miller"] <- MillerCalib(obs = presabs, pred = fold_preds[,f], main = foldname)$slope
    mtext(names(crossval)[element], side = 3, outer = TRUE)
  }  # end for f
#  dev.off()
}  # end triple for

crossval


# BAR PLOTS OF EVALUATION MEASURES:

Millers <- sort(sapply(sapply(crossval, `[`, "Miller"), mean), decreasing = TRUE)
AUCs <- sort(sapply(sapply(crossval, `[`, "AUC"), mean))
TSSs <- sort(sapply(sapply(crossval, `[`, "TSS"), mean))
names(Millers) <- sapply(strsplit(names(Millers), "\\."), `[`, 1)
names(AUCs) <- sapply(strsplit(names(AUCs), "\\."), `[`, 1)
names(TSSs) <- sapply(strsplit(names(TSSs), "\\."), `[`, 1)

par(mfrow = c(1, 1), mar = c(2.5, 6.2, 1.5, 1), oma = c(0,0,0,0))
barplot(Millers, las = 2, cex.names = 0.4, cex.ax = 0.8, space = 1, border = NA, main = "Miller slope", horiz = TRUE)
abline(v = 1, lwd = 0.5, col = "blue")
abline(v = 1.5, lwd = 0.5, col = "red")
barplot(AUCs, las = 2, cex.names = 0.4, cex.ax = 0.8, space = 1, border = NA, main = "AUC", horiz = TRUE, xlim = c(0, 1))
abline(v = 1, lwd = 0.5, col = "blue")
abline(v = 0.7, lwd = 0.5, col = "red")
barplot(TSSs, las = 2, cex.names = 0.4, cex.ax = 0.8, space = 1, border = NA, main = "TSS", horiz = TRUE, xlim = c(0, 1))
abline(v = 1, lwd = 0.5, col = "blue")
abline(v = 0.4, lwd = 0.5, col = "red")


circular_barplots <- function(metrics, pass_val, ideal_val = 1, line_labels) {
  # code adapted from https://www.r-graph-gallery.com/296-add-labels-to-circular-barplot.html

  data <- data.frame(model = names(metrics), value = metrics)
  empty_bar <- 2  # Set a number of 'empty bars'
  to_add <- matrix(NA, empty_bar, ncol(data))  # Add lines to the initial dataset
  colnames(to_add) <- colnames(data)
  data <- rbind(data, to_add)
  data$id <- seq(1, nrow(data))

  # Get the name and the y position of each label
  label_data <- data
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id - 0.5) / number_of_bar  # I subtract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle + 180, angle)

  line_labels <- c(0, pass_val, ideal_val)  # euze

  ggplot(data, aes(x = as.factor(id), y = value)) +   # Note that id is a factor. If x is numeric, there is some space between the first bar
    geom_bar(stat = "identity", fill = "lightcyan3") +
    geom_hline(yintercept = ideal_val, col = "darkblue", lty = 2, lwd = 0.2) +
    geom_hline(yintercept = pass_val, col = "red", lty = 2, lwd = 0.2) +
    ylim(-1, 4) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1, 4), "cm")
    ) +
    coord_polar(start = 0) +
    geom_text(
      data = label_data,
      aes(x = id, y = value + 0.5, label = model, hjust = hjust),
      color = "black",
      #fontface = "bold",
      size = 2,
      angle = label_data$angle,
      inherit.aes = FALSE
    ) +
    annotate(geom = "text", x = -1, y = line_labels, label = line_labels, size = 1.8)  # euze
}

#pdf("images/Fig_S4_AUCbarplot.pdf")
circular_barplots(AUCs, pass_val = c(0.7))
#dev.off()

#pdf("images/Fig_S5_TSSbarplot.pdf")
circular_barplots(TSSs, pass_val = c(0.4))
#dev.off()

#pdf("images/Fig_S6_Millerbarplot.pdf")
circular_barplots(Millers, pass_val = c(0.5, 1.5))
#dev.off()


# calculate mean cross-validation performance:
crossval_means <- lapply(crossval, function(x) sapply(x, mean))
crossval_means <- as.data.frame(t(as.data.frame(crossval_means)))
head(crossval_means)

# calculate difference between each metric and the ideal value 1:
crossval_diff <- 1 - crossval_means
crossval_diff$Miller <- abs(1 - abs(crossval_means$Miller))
head(crossval_diff)
sapply(crossval_diff, range)

# rank models according to their crossval differences from 1:
rownames(crossval_diff)[order(crossval_diff$AUC)]
rownames(crossval_diff)[order(crossval_diff$TSS)]
rownames(crossval_diff)[order(crossval_diff$Miller)]

# see which models pass with acceptable metrics:
head(crossval_diff)
crossval_diff$pass <- FALSE
crossval_diff$pass[crossval_diff$AUC <= 1-0.7 & crossval_diff$TSS <= 1-0.4 & crossval_diff$Miller <= 0.5] <- TRUE
head(crossval_diff)

selected_mods <- rownames(crossval_diff)[crossval_diff$pass == TRUE]
selected_mods
length(selected_mods)  # 59
setdiff(names(crossval), selected_mods)

selected_mods_method <- sapply(strsplit(selected_mods, "_"), `[`, 3)
table(selected_mods_method)
# BART  GAM  GLM
#   30   15   14
# all BART models selected, all RF models excluded, ~half GLMs and GAMs selected


# BOXPLOTS OF VARIABLE CONTRIBUTIONS ####

var_contribs
head(var_contribs)
var_contribs_sel <- var_contribs

selected_mods

for(i in 1:nrow(var_contribs_sel)) for (climod in climods) {
  spc <- as.character(var_contribs_sel[i, "species"])
  meth <- as.character(var_contribs_sel[i, "method"])
  mod <- paste(spc, climod, meth, sep = "_")
  if (!(mod %in% selected_mods))
    var_contribs_sel[var_contribs_sel$species == spc & var_contribs_sel$method == meth, climod] <- NA
}

var_contribs_sel
nrow(na.omit(var_contribs_sel))  # 44
sum(!is.na(var_contribs_sel[ , climods])) / 4  # 59
length(selected_mods)  # 59 - ok


#pdf("images/Fig_S3_varimp_boxplots.pdf", width = 6, height = 7.5)
par(mfrow = c(6, 3))
for (spc in spp) for(m in meths[-grep("RF", meths)]) {
  bp <- var_contribs_sel[var_contribs_sel$species == spc & var_contribs_sel$method == m, c("variable", climods)]
  tbp <- as.data.frame(t(bp[ , -1]))
  colnames(tbp) <- bp[, "variable"]
  upper_margin <- ifelse(spc == "Aa", 2, 0)
  left_margin <- ifelse(m == "GLM", 4.2, 3.2)
  par(mar = c(3, left_margin, upper_margin, 0.5))
  #empty <- nrow(na.omit(bp)) == 0
  empty <- all(is.na(bp[ , -1]))
  if (empty) plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  else boxplot(tbp, names = NA, las = 2)
  if (m == "GLM") mtext(spc, side = 2, line = 2.8)
  if (spc == "Aa") mtext(m, side = 3, line = 0.5)
  if (!empty) mtext(bp[1, "variable"], side = 1, at = 1, line = 0.5, cex = 0.65)
  if (!empty) mtext(bp[2, "variable"], side = 1, at = 2, line = 1.5, cex = 0.65)
  if (!empty) mtext(bp[3, "variable"], side = 1, at = 3, line = 0.5, cex = 0.65)
  if (!empty) mtext(bp[4, "variable"], side = 1, at = 4, line = 1.5, cex = 0.65)
  n_mods <- sum(!is.na(bp[1, -1]))
  if (!empty) text(x = 0.5, y = max(bp[ , -1], na.rm = TRUE), adj = c(0, 1), labels = paste("N =", n_mods))  # , col = "grey20"
}
#dev.off()


# CALCULATE FAVOURABILITY ####

curr_pred_names <- names(preds_cv)[grep(pattern = "_alldat", x = names(preds_cv))]

future_pred_names <- names(preds_cv)[grep(pattern = "_future", x = names(preds_cv))]
future_pred_names

# add favourability to the preds table:
for (p in curr_pred_names) {
  preds_cv[ , paste0(substr(p, 1, nchar(p) - 7), "_fav")] <- fuzzySim::Fav(obs = dat_pres[ , substr(p, 1, 2)], pred = preds_cv[ , p])
}
rm(p)
tail(names(preds_cv))

for (p in future_pred_names) {
  preds_cv[ , paste0(substr(p, 1, nchar(p) - 3), "_fav")] <- fuzzySim::Fav(obs = dat_pres[ , substr(p, 1, 2)], pred = preds_cv[ , p])
}
names(preds_cv)[1:30]
tail(names(preds_cv))

# get a table with only fav predictions:
fav_preds <- preds_cv[ , grep("_fav", names(preds_cv))]
head(fav_preds)
names(fav_preds)
fav_preds_curr <- fav_preds[1:120]
fav_preds_fut <- fav_preds[grep("fut", names(fav_preds))]
head(fav_preds_curr)
head(fav_preds_fut)


# GET FAVOURABILITY ALSO FOR UTM ZONE UNIONS ####

addPolygons(addProviderTiles(leaflet(spTransform(utm10, "+proj=longlat")), provider = "OpenStreetMap"), label = ~ utm10, opacity = 0.4, col = "blue")  # mouse over the polygons to see the label of each one

utm10char <- as.character(dat$utm10)
dat$zone_unions <- as.integer((startsWith(utm10char, "Q") & dat$utm10x %in% 4:6) | (startsWith(utm10char, "T") & dat$utm10x %in% 3:5) | (startsWith(utm10char, "Y") & dat$utm10x %in% 4:6) | (startsWith(utm10char, "B") & dat$utm10x %in% 3:5))
par(mfrow = c(1, 1))
choroLayer(spdf = utm10, df = dat, spdfid = "utm10", dfid = "utm10", var = "zone_unions", border = NA)  # check OK
dat$zone_unions <- as.integer(dat$zone_unions & is.na(dat$REM_Tmax30_1971_2000))
choroLayer(spdf = utm10, df = dat, spdfid = "utm10", dfid = "utm10", var = "zone_unions", border = NA)  # check OK

# remove some cells that also get selected and are not exactly on UTM union:
zone_unions <- subset(utm10, utm10 %in% dat[dat$zone_unions == 1, "utm10"])
addPolygons(addProviderTiles(leaflet(spTransform(zone_unions, "+proj=longlat")), provider = "OpenStreetMap"), label = ~ utm10, opacity = 0.4, col = "blue")  # mouse over the polygons to see the label of each one
dat$zone_unions[dat$utm10 %in% c("QA44", "QA43", "QA42", "TE59", "YK40", "YQ49", "BE53", "YJ44", "YJ42", "YJ43", "YJ49", "YH46", "BD50", "BC58")] <- 0
choroLayer(spdf = utm10, df = dat, spdfid = "utm10", dfid = "utm10", var = "zone_unions", border = NA)  # check OK
zone_unions <- subset(utm10, utm10 %in% dat[dat$zone_unions == 1, "utm10"])
addPolygons(addProviderTiles(leaflet(spTransform(zone_unions, "+proj=longlat")), provider = "OpenStreetMap"), label = ~ utm10, opacity = 0.4, col = "blue")  # # check OK

length(zone_unions)  # 318 UTM zone union cells
nrow(fav_preds)  # 5665
names(fav_preds)
names(preds_cv)
names(zone_unions)
names(dat)

names(fav_preds_curr)

fav_preds_withZoneUnions <- merge(dat[ , c("utm10", "rownum", "utm10x", "utm10y", "lonUTM", "latUTM", "zone_unions")], preds_cv[ , c("utm10", names(fav_preds_curr), names(fav_preds_fut))], all = TRUE)
ncol(fav_preds_withZoneUnions)  # 247
nrow(fav_preds_withZoneUnions)  # 6215
head(fav_preds_withZoneUnions)

sum(is.na(fav_preds_withZoneUnions$Aa_MM5_BART_fav))  # 550
sum(fav_preds_withZoneUnions$zone_unions)  # 318

fav_preds_withZoneUnions <- fav_preds_withZoneUnions[!(is.na(fav_preds_withZoneUnions$Aa_MM5_BART_fav) & fav_preds_withZoneUnions$zone_unions == 0), ]
nrow(fav_preds_withZoneUnions)  # 5983
6215 - 550 + 318  # 5983
5665 + 318  # 5983
head(fav_preds_withZoneUnions[ , 1:20])


# interpolate favourability predictions:

names(utm10)
names(preds_cv)
fav_interp <- slot(utm10, "data")[ , "utm10", drop = FALSE]

head(fav_interp)
nrow(fav_interp)  # 6215

utm10$zone_unions <- 0
utm10$zone_unions[utm10$utm10 %in% dat[dat$zone_unions == 1, "utm10"]] <- 1
head(utm10)

names(preds_cv)
favs <- names(preds_cv)[grep("_fav", names(preds_cv))]
favs
length(favs)

slot(utm10, "data") <- data.frame(slot(utm10, "data"), preds_cv[match(as.character(slot(utm10, "data")$utm10), as.character(preds_cv$utm10)), favs])
na_rows <- is.na(utm10$Aa_MM5_BART_fav)

# loop to interpolate each favourability prediction (MAY TAKE A FEW HOURS TO COMPUTE!):
gc()
n <- 0
for (f in favs) {
  n <- n + 1
  message("interpolating fav ", n, " of ", length(favs), " - ", f, "...\n")
  fav_idw <- gstat::idw(as.formula(paste(f, "~ 1")), utm10[!na_rows, ], newdata = utm10, maxdist = 19500)
  fav_interp[ , paste0(f, "_idw")] <- fav_idw@data[ , "var1.pred"]
  rm(fav_idw)
  gc()
}  # end for v


choroLayer(spdf = utm10, df = fav_interp, spdfid = "utm10", dfid = "utm10", var = "Aa_MM5_BART_fav_idw", border = NA, col = spectral, breaks = brks, legend.pos = "n")

choroLayer(spdf = utm10, df = fav_interp, spdfid = "utm10", dfid = "utm10", var = "Aa_MM5_RF_fav_idw", border = NA, col = spectral, breaks = brks, legend.pos = "n")


head(fav_interp)
ncol(fav_interp)  # 241
nrow(fav_interp)  # 6215
nrow(utm10)  # 6215


# correlations between original and interpolated fav:

interp_favs <- names(fav_interp)[grep("_idw", names(fav_interp))]
interp_favs
interp_corrs <- rep(NA, length(interp_favs))
names(interp_corrs) <- interp_favs
for (f in interp_favs) {
  interp <- fav_interp[ , f]
  orig <- utm10@data[ , substr(f, 1, nchar(f) - 4)]
  interp_corrs[f] <- cor(orig, interp, use = "pairwise.complete.obs")
}  # end for v

summary(interp_corrs)
boxplot(interp_corrs)


# MAP FAVOURABILITY PREDICTIONS ####

names(fav_interp)

fav_interp_pres <- fav_interp[ , 1:121]
fav_interp_fut <- fav_interp[ , c(1, grep("_fut", names(fav_interp)))]
names(fav_interp_pres)
names(fav_interp_fut)

for (spc in spp) {
  #jpeg(paste0("images/pred_maps/", spc, "_fav_present.jpg"), width = 800, height = 900)
  par(mar = c(0, 0, 1.2, 0), mfrow = c(5, 4))
  for (f in names(fav_interp_pres)[grep(spc, names(fav_interp_pres))]) {
    mod <- substr(f, 1, nchar(f) - 8)
    choroLayer(spdf = utm10, df = fav_interp_pres, spdfid = "utm10", dfid = "utm10", var = f, border = NA, col = spectral, breaks = brks, legend.pos = "n")
    if (mod %in% selected_mods)  title(mod, cex.main = 2, font.main = 2)
    else  title(paste0("(", mod, ")"), cex.main = 2)
  }
  legendChoro(pos = "bottomright", border = NA, breaks = brks_leg, col = spectral_leg, symbol = "box", title.txt = "", nodata = FALSE, cex = 1.45, values.cex = 1)
  #dev.off()
}

for (spc in spp) {
  #jpeg(paste0("images/pred_maps/", spc, "_fav_future.jpg"), width = 800, height = 900)
  par(mar = c(0, 0, 1.2, 0), mfrow = c(5, 4))
  for (v in names(fav_interp_fut)[grep(spc, names(fav_interp_fut))]) {
    choroLayer(spdf = utm10, df = fav_interp_fut, spdfid = "utm10", dfid = "utm10", var = v, border = NA, col = spectral, breaks = brks, legend.pos = "n")
    mod <- substr(v, 1, nchar(v) - 12)
    if (mod %in% selected_mods)  title(mod, cex.main = 2, font.main = 2)
    else  title(paste0("(", mod, ")"), cex.main = 2)
  }
  legendChoro(pos = "bottomright", border = NA, breaks = brks_leg, col = spectral_leg, symbol = "box", title.txt = "", nodata = FALSE, cex = 1.45, values.cex = 1)
  #dev.off()
}


# INTERSECT PREDICTIONS WITH SPATIAL TREND ####

# don't use interp here because it's for analysis (PCA)
names(dat_pres)
names(fav_preds_curr)
nrow(dat_pres)
nrow(fav_preds_curr)


# present:

fav_preds_curr_withTSA <- fav_preds_curr
fav_preds_curr_withTSA[] <- NA
fav_preds_curr_withTSA
names(fav_preds_curr_withTSA)

for(i in names(fav_preds_curr_withTSA)) {
  spc <- substr(i, 1, 2)
  fav_preds_curr_withTSA[ , i] <- fuzzyOverlay(data = data.frame(dat_pres, fav_preds_curr), overlay.cols = c(i, paste0(spc, "_TS_fav")))
}
head(fav_preds_curr_withTSA)


# future:

fav_preds_fut_withTSA <- fav_preds_fut
fav_preds_fut_withTSA[] <- NA
fav_preds_fut_withTSA
names(fav_preds_fut_withTSA)

for(i in names(fav_preds_fut_withTSA)) {
  spc <- substr(i, 1, 2)
  fav_preds_fut_withTSA[ , i] <- fuzzyOverlay(data = data.frame(dat_pres, fav_preds_fut), overlay.cols = c(i, paste0(spc, "_TS_fav")))
}
head(fav_preds_fut_withTSA)


# PCA OF FAV PREDS OF DIFFERENT ALGORITHMS ####

names(fav_preds_curr_withTSA)
selected_mods

PCAs <- vector("list", length(spp))
names(PCAs) <- spp

for (spc in spp) {
  selmods <- selected_mods[grep(spc, selected_mods)]
  PCAs[[spc]] <- prcomp(fav_preds_curr_withTSA[ , names(fav_preds_curr_withTSA) %in% paste0(selmods, "_fav")], center = FALSE, scale. = FALSE)
}

PCAs

PCAs_future <- vector("list", length(spp))
names(PCAs_future) <- spp

for (spc in spp) {
  selmods <- selected_mods[grep(spc, selected_mods)]
  PCAs_future[[spc]] <- prcomp(fav_preds_fut_withTSA[ , names(fav_preds_fut_withTSA) %in% paste0(selmods, "_fut_fav")], center = FALSE, scale. = FALSE)
}

length(PCAs)  # 6
lapply(PCAs, summary)
par(mfrow = modEvA::arrangePlots(length(PCAs)), mar = c(2, 2, 1, 1))
lapply(PCAs, plot)
lapply(PCAs, biplot)
sort(sapply(PCAs, function(x) summary(x)$importance["Proportion of Variance", "PC1"]))
sort(sapply(PCAs_future, function(x) summary(x)$importance["Proportion of Variance", "PC1"]))


range(sapply(PCAs, function(x) range(as.data.frame(x$x)[ , "PC1"])))
range(sapply(PCAs_future, function(x) range(as.data.frame(x$x)[ , "PC1"])))
range(sapply(PCAs, function(x) range(as.data.frame(x$rotation)[ , "PC1"])))
range(sapply(PCAs_future, function(x) range(as.data.frame(x$rotation)[ , "PC1"])))

dim(PCAs[["Aa"]]$x)  # 5665 14
dim(PCAs[["Aa"]]$rotation)  # 14 14
dim(PCAs[["Mm"]]$x)  # now 5665 5
dim(PCAs[["Mm"]]$rotation)  # 5 5

PCAs[["Mm"]]$rotation[ , "PC1"]

sapply(spp, function(x) nrow(PCAs[[x]]$rotation))
sapply(spp, function(x) range(PCAs[[x]]$rotation[ , "PC1"]))


hist(sapply(PCAs, function(x) range(as.data.frame(x$x)[ , "PC1"])))
hist(sapply(PCAs_future, function(x) range(as.data.frame(x$x)[ , "PC1"])))

PC1s <- data.frame(utm10 = dat_pres$utm10, sapply(PCAs, function(x) x$x[ , "PC1"]))
nrow(PC1s)  # 5665
nrow(dat_pres)  # 5665
head(PC1s)

PC1s_fut <- data.frame(utm10 = dat_pres$utm10, sapply(PCAs_future, function(x) x$x[ , "PC1"]))
nrow(PC1s_fut)  # 5665
head(PC1s_fut)


# PLOT CONSENSUS PCA MAPS ####

spectral2 <- rev(colorRampPalette(spectral)(13))
greys <- brewer.pal(9, "Greys")
greys <- colorRampPalette(greys[2:length(greys)])(13)

#jpeg("images/consensus_present.jpg", width = 800, height = 900)
par(mfrow = modEvA::arrangePlots(length(spp)), mar = c(0, 0, 2, 0))
for (i in 2:length(PC1s)) {
  choroLayer(spdf = utm10, df = PC1s, spdfid = "utm10", dfid = "utm10", var = names(PC1s)[i], border = NA, legend.pos = "n", method = "equal", col = spectral2)  # rev(greys)
  title(species[species$spp == names(PC1s)[i], "species"], font.main = 3, cex.main = 2.5)
}
#dev.off()

#jpeg("images/consensus_future.jpg", width = 800, height = 900)
par(mfrow = modEvA::arrangePlots(length(spp)), mar = c(0, 0, 2, 0))
for (i in 2:length(PC1s)) {
  choroLayer(spdf = utm10, df = PC1s_fut, spdfid = "utm10", dfid = "utm10", var = names(PC1s)[i], border = NA, legend.pos = "n", method = "equal", col = spectral2)  # rev(greys)
  title(species[species$spp == names(PC1s)[i], "species"], font.main = 3, cex.main = 2.5)
}
#dev.off()


# CALCULATE MINIMUM AND PCA-WEIGHTED MEAN FAVOURABILITY ####

lapply(PC1s[ , -1], summary)

par(mfrow = modEvA::arrangePlots(length(spp)), mar = c(2, 2, 1, 1))
lapply(PC1s[ , -1], function(x) plot(sort(x), pch = 20))

preds_PCAweighted <- preds_PCAweighted_fut <- preds_min <- preds_min_fut <- data.frame(matrix(NA, nrow(dat_pres), ncol = length(spp)))
names(preds_PCAweighted) <- names(preds_PCAweighted_fut) <- names(preds_min) <- names(preds_min_fut) <- spp
head(preds_PCAweighted)
head(preds_PCAweighted_fut)
head(preds_min)
head(preds_min_fut)

for (spc in spp) {
  preds_selmods_names <- selected_mods[grep(spc, selected_mods)]
  preds_selmods <- fav_preds_curr_withTSA[ , paste0(preds_selmods_names, "_fav")]
  preds_selmods_fut <- fav_preds_fut_withTSA[ , paste0(preds_selmods_names, "_fut_fav")]
  preds_PCAweighted[ , spc] <- apply(preds_selmods, 1, weighted.mean, w = abs(PCAs[[spc]]$rotation[ , "PC1"]))
  preds_PCAweighted_fut[ , spc] <- apply(preds_selmods_fut, 1, weighted.mean, w = abs(PCAs_future[[spc]]$rotation[ , "PC1"]))
  preds_min[ , spc] <- apply(preds_selmods, 1, min)
  preds_min_fut[ , spc] <- apply(preds_selmods_fut, 1, min)
}

head(preds_selmods)
head(preds_PCAweighted)
head(preds_PCAweighted_fut)
head(preds_min)
head(preds_min_fut)

preds_PCAweighted <- data.frame(utm10 = dat_pres$utm10, preds_PCAweighted)
preds_PCAweighted_fut <- data.frame(utm10 = dat_pres$utm10, preds_PCAweighted_fut)
preds_min <- data.frame(utm10 = dat_pres$utm10, preds_min)
preds_min_fut <- data.frame(utm10 = dat_pres$utm10, preds_min_fut)

nrow(preds_PCAweighted)  # 5665
nrow(preds_min)  # 5665


# interpolate consensus-weighted and minimum favourability maps:

names(preds_PCAweighted)[-1] <- paste0(names(preds_PCAweighted)[-1], "_consens_pres")
names(preds_PCAweighted_fut)[-1] <- paste0(names(preds_PCAweighted_fut)[-1], "_consens_fut")
names(preds_PCAweighted)
names(preds_PCAweighted_fut)

slot(utm10, "data") <- data.frame(slot(utm10, "data"), preds_PCAweighted[match(as.character(slot(utm10, "data")$utm10), as.character(preds_PCAweighted$utm10)), -1])
slot(utm10, "data") <- data.frame(slot(utm10, "data"), preds_PCAweighted_fut[match(as.character(slot(utm10, "data")$utm10), as.character(preds_PCAweighted_fut$utm10)), -1])
na_rows2 <- is.na(utm10$Ts_consens_fut)
all.equal(na_rows, na_rows2)  # TRUE


names(preds_min)
names(preds_min_fut)

names(preds_min)[-1] <- paste0(names(preds_min)[-1], "_Fmin_pres")
names(preds_min_fut)[-1] <- paste0(names(preds_min_fut)[-1], "_Fmin_fut")
names(preds_min)
names(preds_min_fut)

slot(utm10, "data") <- data.frame(slot(utm10, "data"), preds_min[match(as.character(slot(utm10, "data")$utm10), as.character(preds_min$utm10)), -1])
slot(utm10, "data") <- data.frame(slot(utm10, "data"), preds_min_fut[match(as.character(slot(utm10, "data")$utm10), as.character(preds_min_fut$utm10)), -1])
na_rows2 <- is.na(utm10$Ts_Fmin_fut)
all.equal(na_rows, na_rows2)  # TRUE

names(utm10)

consens_interp <- slot(utm10, "data")[ , "utm10", drop = FALSE]

par(mfrow = c(3, 2))
gc()
cons <- c(grep("_consens", names(utm10)), grep("_Fmin", names(utm10)))
n <- 0
for (f in cons) {
  n <- n + 1
  name <- names(utm10)[f]
  message("interpolating column ", n, " of ", length(cons), " - ", name, "...\n")
  fav_idw <- idw(as.formula(paste(name, "~ 1")), utm10[!na_rows, ], newdata = utm10, maxdist = 19500)
  consens_interp[ , paste0(name, "_idw")] <- fav_idw@data[ , "var1.pred"]
  rm(fav_idw)
  choroLayer(spdf = utm10, df = consens_interp, spdfid = "utm10", dfid = "utm10", var = paste0(name, "_idw"), border = NA, col = spectral, breaks = brks, legend.pos = "n"); title(name)
  gc()
}


sapply(preds_PCAweighted[ , -1], max)
sapply(preds_PCAweighted_fut[ , -1], max)

sort(sapply(preds_min[ , -1], max))
sort(sapply(preds_min_fut[ , -1], max))

names(consens_interp)

#pdf("images/Fig_2_fav_wmean_present.pdf", width = 7, height = 8)
par(mfrow = modEvA::arrangePlots(length(spp)), mar = c(0, 0, 2, 0))
for (i in grep("_consens_pres", names(consens_interp))) {
  choroLayer(spdf = utm10, df = consens_interp, spdfid = "utm10", dfid = "utm10", var = names(consens_interp)[i], border = NA, legend.pos = "n", breaks = brks, col = spectral)  # greys100
  title(species[species$spp == substr(names(consens_interp)[i], 1, 2), "species"], font.main = 3, cex.main = 2.5)
}
legendChoro(pos = "bottomright", border = NA, breaks = brks_leg, col = spectral_leg, symbol = "box", title.txt = "", nodata = FALSE)
#dev.off()

#pdf("images/Fig_3_fav_wmean_future.pdf", width = 7, height = 8)
par(mfrow = modEvA::arrangePlots(length(spp)), mar = c(0, 0, 2, 0))
for (i in grep("_consens_fut", names(consens_interp))) {
  choroLayer(spdf = utm10, df = consens_interp, spdfid = "utm10", dfid = "utm10", var = names(consens_interp)[i], border = NA, legend.pos = "n", breaks = brks, col = spectral)  # greys100
  #title(names(PCAs)[i])
  title(species[species$spp == substr(names(consens_interp)[i], 1, 2), "species"], font.main = 3, cex.main = 2.5)
}
legendChoro(pos = "bottomright", border = NA, breaks = brks_leg, col = spectral_leg, symbol = "box", title.txt = "", nodata = FALSE)
#dev.off()


for (spc in spp) {
  print(paste(spc, cor(abs(PC1s[ , spc]), preds_PCAweighted[ , grep(spc, names(preds_PCAweighted))])))
}  # all > 0.9999


names(consens_interp)

#jpeg("images/Fig_S19_favmin_pres.jpg", width = 800, height = 900)
par(mfrow = modEvA::arrangePlots(length(spp)), mar = c(0, 0, 2, 0))
for (i in grep("Fmin_pres", names(consens_interp))) {
  choroLayer(spdf = utm10, df = consens_interp, spdfid = "utm10", dfid = "utm10", var = names(consens_interp)[i], border = NA, legend.pos = "n", breaks = brks, col = spectral)  # greys100
  title(species[species$spp == substr(names(consens_interp)[i], 1, 2), "species"], font.main = 3, cex.main = 2.5)
}
legendChoro(pos = "bottomright", border = NA, breaks = brks_leg, col = spectral_leg, symbol = "box", title.txt = "", nodata = FALSE, values.cex = 1.3)
#dev.off()

#jpeg("images/Fig_S20_favmin_fut.jpg", width = 800, height = 900)
par(mfrow = modEvA::arrangePlots(length(spp)), mar = c(0, 0, 2, 0))
for (i in grep("Fmin_fut", names(consens_interp))) {
  choroLayer(spdf = utm10, df = consens_interp, spdfid = "utm10", dfid = "utm10", var = names(consens_interp)[i], border = NA, legend.pos = "n", breaks = brks, col = spectral)  # greys100
  #title(names(PCAs)[i])
  title(species[species$spp == substr(names(consens_interp)[i], 1, 2), "species"], font.main = 3, cex.main = 2.5)
}
legendChoro(pos = "bottomright", border = NA, breaks = brks_leg, col = spectral_leg, symbol = "box", title.txt = "", nodata = FALSE, values.cex = 1.3)
#dev.off()


# QUANTIFY CHANGE ####

head(preds_PCAweighted)
head(preds_PCAweighted_fut)

preds_PCAweighted_change <- data.frame(utm10 = preds_PCAweighted[ , "utm10"], preds_PCAweighted_fut[ , -1] - preds_PCAweighted[ , -1])
head(preds_PCAweighted_change)

sapply(preds_PCAweighted_change[ , -1], range)
max(abs(preds_PCAweighted_change[ , -1]))
max(preds_PCAweighted_change[ , -1])
min(preds_PCAweighted_change[ , -1])
brks2 <- seq(-0.3, 0.3, by = 0.05)

ramp_posneg <- colorRampPalette(c('blue', 'gray90', 'red')) (length(brks2))

par(mfrow = modEvA::arrangePlots(length(spp)), mar = c(0, 0, 2, 0))
for (i in 2:ncol(preds_PCAweighted_change)) {
  choroLayer(spdf = utm10, df = preds_PCAweighted_change, spdfid = "utm10", dfid = "utm10", var = names(preds_PCAweighted_change)[i], border = NA, breaks = brks2, legend.values.rnd = 2, legend.title.txt = "", col = ramp_posneg)
  title(species[species$spp == names(preds_PCAweighted_change)[i], "species"], font.main = 3, cex.main = 2.5)
}

for (i in 2:ncol(preds_PCAweighted)) {
  print(names(preds_PCAweighted)[i])
  print(t.test(preds_PCAweighted[ , i], preds_PCAweighted_fut[ , i]), paired = TRUE)
}

for (i in 2:ncol(preds_PCAweighted)) {
  print(names(preds_PCAweighted)[i])
  print(wilcox.test(preds_PCAweighted[ , i], preds_PCAweighted_fut[ , i]), paired = TRUE)
}


# INTERPOLATE CHANGE ####

head(preds_PCAweighted_change)
names(preds_PCAweighted_change)[-1] <- paste0(names(preds_PCAweighted_change)[-1], "_change")

slot(utm10, "data") <- data.frame(slot(utm10, "data"), preds_PCAweighted_change[match(as.character(slot(utm10, "data")$utm10), as.character(preds_PCAweighted_change$utm10)), -1])

names(utm10)

change_interp <- slot(utm10, "data")[ , "utm10", drop = FALSE]

change_columns <- grep("_change", names(utm10))
par(mfrow = modEvA::arrangePlots(length(spp)))
gc()
n <- 0
for (f in change_columns) {
  n <- n + 1
  name <- names(utm10)[f]
  message("interpolating column ", n, " of ", length(change_columns), " - ", name, "...\n")
  fav_idw <- idw(as.formula(paste(name, "~ 1")), utm10[!na_rows, ], newdata = utm10, maxdist = 19500)
  change_interp[ , paste0(name, "_idw")] <- fav_idw@data[ , "var1.pred"]
  rm(fav_idw)
  choroLayer(spdf = utm10, df = change_interp, spdfid = "utm10", dfid = "utm10", var = paste0(name, "_idw"), border = NA, col = ramp_posneg, breaks = brks2, legend.pos = "n"); title(name)
  gc()
}  # end for f

head(change_interp)
summary(change_interp[ , -1])
range(change_interp[ , -1], na.rm = TRUE)
brks3 <- seq(-0.25, 0.25, by = 0.025)
ramp_posneg2 <- colorRampPalette(c('blue', 'gray86', 'red')) (length(brks3))

brks3_leg <- seq(brks3[1], brks3[length(brks3)], by = 0.05)
ramp_posneg2_leg <- colorRampPalette(c('blue', 'gray86', 'red')) (length(brks3_leg))

bbox(utm10)

#pdf("images/Fig_5_fav_change_maps.pdf", width = 7, height = 8)
par(mfrow = modEvA::arrangePlots(length(spp)), mar = c(0, 0, 2, 0))
for (i in 2:ncol(preds_PCAweighted_change)) {
  choroLayer(spdf = utm10, df = change_interp, spdfid = "utm10", dfid = "utm10", var = names(change_interp)[i], border = NA, breaks = brks3, legend.values.rnd = 2, legend.title.txt = "", col = ramp_posneg2, legend.pos = "n")
  title(species[species$spp == substr(names(change_interp)[i], 1, 2), "species"], font.main = 3, cex.main = 2.5)
}
legendChoro(pos = "bottomright", border = NA, breaks = brks3_leg, col = ramp_posneg2_leg, symbol = "box", title.txt = "", nodata = FALSE, values.cex = 1.3)
#dev.off()


# BARPLOT OF FUZZY RANGE CHANGE ####

names(preds_PCAweighted)
names(preds_PCAweighted_fut)

#pdf("images/Fig_4_range_change_barplots.pdf", width = 4, height = 6)
par(mfrow = modEvA::arrangePlots(length(spp)), mar = c(6, 3, 2, 1))
for (i in 2:ncol(preds_PCAweighted)) {
  fuzzyRangeChange(preds_PCAweighted[ , i], preds_PCAweighted_fut[ , i], ylim = c(-0.65, 0.65), las = 2, col = "grey50", border = NA)
  title(species[species$spp == substr(names(change_interp)[i], 1, 2), "species"], font.main = 3, cex.main = 1.3)
}
#dev.off()
