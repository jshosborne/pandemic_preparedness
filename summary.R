in_dir <- commandArgs()[4]
out_dir <- commandArgs()[5]
yr_max <- commandArgs()[6]
monthly <- as.logical(commandArgs()[7])
mo <- commandArgs()[8]
f_name <- commandArgs()[9]

package_list <- c('seegSDM', 'raster', 'maptools', 'sp', 'rgdal', 'rgeos')
lapply(package_list, library, character.only = TRUE)
source('helper_functions.R')
source('spillover_functions.R')

# read in covariates and data
cov_name_list <- read.csv(paste0(in_dir, '/cov_name_list.csv'), as.is=T)
cov_names <- cov_name_list$cov_name

covs <- stack(paste0(in_dir, '/bricks/pred_stack_', yr_max, '.grd'))
names(covs) <- cov_names

## SUMMARY CALCULATIONS

# make lists of all files in output directories
data_files <- list.files(path = paste0(out_dir, '/data'), full.names = TRUE)
dat_all <- c()
stats_files <- list.files(path = paste0(out_dir, '/stats_output'), full.names = TRUE)
stats_list <- c()
model_files <- list.files(path = paste0(out_dir, '/model_output'), full.names = TRUE)
model_list <- c()
for (i in 1:length(data_files)){
  model_list[[i]] <- get(load(model_files[[i]]))
  stats_list[[i]] <- get(load(stats_files[[i]]))
  dat_all <- rbind(dat_all, read.csv(data_files[[i]]))
}

# write data files - all bootstrap data
if (!monthly) {
  dat_all.bg <- dat_all[dat_all$PA==0,]
  write.csv(dat_all.bg, paste0(out_dir, '/dat_all_bg.csv'), row.names=F)
  dat_all.occ <- dat_all[dat_all$PA==1,]
  write.csv(dat_all.occ, paste0(out_dir, '/dat_all_occ.csv'), row.names=F)
}

# get prediction raster and save it as GeoTiff
if (monthly) {
  preds_files <- list.files(path = paste0(out_dir, '/preds_output_', mo), full.names = TRUE)
} else {
  preds_files <- list.files(path = paste0(out_dir, '/preds_output'), full.names = TRUE)
}
preds_list <- lapply(preds_files, raster)
preds_list <- stack(preds_list)

preds_sry <- combinePreds(preds_list)
names(preds_sry) <- c('mean', 'median', 'lowerCI', 'upperCI')

# get synoptic fit statistics for summarized model...
# ...mean
syn_stats <- get_fit_stats(dat_all, preds_sry$mean, get_ROC=T)
write.csv(syn_stats$fit_stats, paste0(out_dir, '/synoptic_stats_mean.csv'))
write.csv(cbind(syn_stats$fpr, syn_stats$sensitivity), paste0(out_dir, '/synoptic_stats_mean_ROC.csv'))
# ...upper CI
syn_stats_uppCI <- get_fit_stats(dat_all, preds_sry$upperCI)
write.csv(syn_stats_uppCI, paste0(out_dir, '/synoptic_stats_uppCI.csv'))
#... lower CI
syn_stats_lowCI <- get_fit_stats(dat_all, preds_sry$lowerCI)
write.csv(syn_stats_lowCI, paste0(out_dir, '/synoptic_stats_lowCI.csv'))

# calculate uncertainty in the predictions 
preds_sry$uncertainty <- preds_sry[[4]] - preds_sry[[3]]

# classify predictions based on optimal threshold of mean predictions
preds_sry$binary <- preds_sry$mean
values(preds_sry$binary)[values(preds_sry$binary) >= opt_thresh] <- 1
values(preds_sry$binary)[values(preds_sry$binary) < opt_thresh] <- 0

# write mean prediction, binary prediciton, and uncertainty rasters as GeoTiffs
preds_sry <- combinePreds(preds_list)
names(preds_sry) <- c('mean', 'median', 'lowerCI', 'upperCI')

# save matrix of summary statistics across all bootstraps
stats <- do.call('rbind', stats_list)

# get and save relative influence scores of covariates; make a box plot
relinf <- getRelInf(model_list, plot=F) #for gbm fitting, use summary(gbm.fit$finalModel)

# get MESS raster for preds_sry
MESS <- mess(preds_sry, dat_all.occ[,c('longitude', 'latitude')])

# write outputs - indexed by month if applicable
write.csv(names(preds_sry), paste0(out_dir, '/preds_layer_names.csv'))
if (monthly) {
  writeRaster(preds_sry,
              file = paste0(out_dir, '/preds', '_', mo,  '.tif'),
              format = 'GTiff',
              overwrite=T)
  
  write.csv(relinf, file = paste0(out_dir, '/relinf', '_', mo,  '.csv'))
  write.csv(stats, file = paste0(out_dir, '/stats', '_', mo,  '.csv'))
  writeRaster(x=MESS, filename=paste0(out_dir, '/MESS_calc', '_', mo), format="GTiff", overwrite=T)
} else {
  writeRaster(preds_sry,
              file = paste0(out_dir, '/preds.tif'),
              format = 'GTiff',
              overwrite=T)
  
  write.csv(relinf, file = paste0(out_dir, '/relinf.csv'))
  write.csv(stats, file = paste0(out_dir, '/stats.csv'))
  writeRaster(x=MESS, filename=paste0(out_dir, '/MESS_calc'), format="GTiff", overwrite=T)
}

## SPILLOVER POTENTIAL CALCULATION

index_list <- read.csv('config.csv')$index

for (expmt_idx in index_list){
  print(expmt_idx)
  out_dir <- 
  
  #load in list of base rasters to contextualize risk map in
  rast_list <- read.csv(paste0(in_dir, '/rast_list.csv'), as.is=T)
  rast_stack <- stack()
  for (k in 1:ncol(rast_list)){
    rast <- raster(rast_list[,k])
    rast_stack <- addLayer(rast_stack, rast)
  }
  names(rast_stack) <- names(rast_list)
  
  #make comorbidity summary raster from outputs downloaded from vizhub
  comorb_fname <- paste0(' ')
  admin2_shp <- shapefile(" ")
  if (!file.exists(comorb_fname)){
    dbkd <- read.csv(paste0( ), as.is=T)
    card <- read.csv(paste0( ), as.is=T)
    admin2_df <- admin2_shp@data
    admin2_df$B8_prev <- dbkd$Value[match(admin2_df$ADM0_NAME, dbkd$Location)]
    admin2_df$B2_prev <- card$Value[match(admin2_df$ADM0_NAME, card$Location)]
    comorb_shp <- admin2_shp
    comorb_shp@data <- admin2_df
    writeSpatialShape(comorb_shp, comorb_fname)
  } else {
    comorb_shp <- shapefile(comorb_fname)
  }
  admin2_raster <- raster(' ')
  
  #load in the raster of interest
  full_measure_raster <- stack(paste0(out_dir, '/preds.tif'))
  # what summarization functions are meaningful for each layer? -> 1:1 w/ rast_stack layers
  meas_func <- c('sum', 'sum', 'sum')
  ctxt_func <- c('sum', 'sum', 'sum')
  
  if (remove_disputed == T) {
    disputed_territories <- shapefile(' ')
    rast_stack$population <- mask(rast_stack$population, disputed_territories, updateValue = 0, inverse = TRUE)
  }
  
  for (i in c(1,3,4)) {
    admin2_out <- spillover_map(full_measure_raster[[i]], 
                                rast_stack,
                                mask=NULL,
                                rast_stack_idx=c(1:length(rast_list)),
                                meas_func, 
                                ctxt_func, 
                                out_dir, 
                                poly_dir,
                                admin2_raster,
                                admin2_shp)
    
    pop_at_risk_spill <- get_pop_quint(rr_sf=admin2_out, 
                                       id_col='ADM2_CODE', 
                                       fill_col='geom_hum_cam')
    
    pop_at_risk_comrb <- get_pop_quint(rr_sf=admin2_out, 
                                       id_col='ADM2_CODE', 
                                       fill_col='geom_all')
    
    if (i==1) {
      write.csv(pop_at_risk_spill, paste0(out_dir, '/PaRisk_geom_hum_cam.csv'), row.names=T) 
      write.csv(pop_at_risk_comrb, paste0(out_dir, '/PaRisk_geom_all.csv'), row.names=T) 
      #save output shapefile
      writeSpatialShape(admin2_out, paste0(out_dir, '/admin2_mers_cov_spillover.shp'))
    } else if (i==3) {
      write.csv(pop_at_risk_spill, paste0(out_dir, '/PaRisk_geom_hum_cam_5p.csv'), row.names=T) 
      write.csv(pop_at_risk_comrb, paste0(out_dir, '/PaRisk_geom_all_5p.csv'), row.names=T) 
      #save output shapefile
      writeSpatialShape(admin2_out, paste0(out_dir, '/admin2_mers_cov_spillover_5p.shp'))
    } else if (i==4) {
      write.csv(pop_at_risk_spill, paste0(out_dir, '/PaRisk_geom_hum_cam_95p.csv'), row.names=T) 
      write.csv(pop_at_risk_comrb, paste0(out_dir, '/PaRisk_geom_all_95p.csv'), row.names=T) 
      #save output shapefile
      writeSpatialShape(admin2_out, paste0(out_dir, '/admin2_mers_cov_spillover_95p.shp'))
    }
  }
}
