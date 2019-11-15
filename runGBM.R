jobnum <- commandArgs()[4]
opt_type <- commandArgs()[5]
in_dir <- commandArgs()[6]
out_dir <- commandArgs()[7]
data_file <- commandArgs()[8]
stack_out_dir <- commandArgs()[9]
dat_dir <- commandArgs()[10]
index <- commandArgs()[11]
cov_dir <- commandArgs()[12]
yr_max <- commandArgs()[13]
cov_pctg <- commandArgs()[14]
f_name <- commandArgs()[15]

source('helper_functions.R')
source('model_functions.R')
package_list <- c('seegSDM', 'gbm', 'rgdal', 'dismo', 'maptools')
lapply(package_list, library, character.only = TRUE)

# read in data
dat_orig <- read.csv(data_file)
pop_bias <- raster(paste0(in_dir, '/bias_grid/bias_mask.tif'))
pop_bias <- 1 - pop_bias # invert pop mask values for biasing
covs <- brick('covs.grd')

# sample from w/in time bins: mo_start, yr_start -> mo_end, yr_end
tbins <- calc_tbin_prob(dat_orig)

## CONSTRUCT DATA FOR POINT RECORDS
pt.dat <- point.data(dat_orig, tbins$yr_freqs, tbins$mo_freqs)
  
## SAMPLE BACKGROUND FROM EXPERIMENT-SPECIFIC BUFFER FILE -- 1:1 W/ NO. OF OCCURRENCES
bg.dat <- bg.data(dat_orig, buffer, bias=pop_bias, yr_freqs, mo_freqs)

## CONSTRUCT DATA FOR BUFFER RECORDS
if (f_name == 'none') {
  buffer_data <- shapefile(paste0(dat_dir, '/sampling_shapefiles_', cov_pctg, '/', index, '_buffer_data.shp'))
} else {
  buffer_data <- shapefile(paste0(dat_dir, '/sampling_shapefiles_', cov_pctg, '_', f_name, '/', index, '_buffer_data.shp'))
}
buf.dat <- buffer.data(dat_orig, buffer_data, bias, yr_freqs, mo_freqs)
    
## CONSTRUCT DATA FOR POLYGON RECORDS
if (f_name == 'none') {
  polygon_data <- shapefile(paste0(dat_dir, '/sampling_shapefiles_', cov_pctg, '/', index, '_polygon_data.shp'))
} else {
  polygon_data <- shapefile(paste0(dat_dir, '/sampling_shapefiles_', cov_pctg, '_', f_name, '/', index, '_polygon_data.shp'))
}
smp.dat <- polygon.data(dat_orig, polygon_data, bias, yr_freqs, mo_freqs)
    
## EXTRACT COVARIATE VALUES FOR OCCURRENCE RECORDS
cov_name_list <- read.csv(paste0(in_dir, '/cov_name_list.csv'), as.is=T)
cov_names <- cov_name_list$cov_name
dat.all(smp.dat, buf.dat, pt.dat, bg.dat, cov_name_list, cov_dir, stack_out_dir)

# save data going into bootstrap model
write.csv(dat_all, paste0(out_dir, '/data/dat_all_', jobnum, '.csv'), row.names=F)

## MAKE .CSV OF OPTIMIZED HPARS

# split to make OoS hold-out data
indx <- sample(nrow(dat_all), 0.80*nrow(dat_all))
data_train <- dat_all[indx,]
data_test <- dat_all[-indx,]
write.csv(data_train, file = paste0(out_dir, '/data_train_', jobnum, '.csv'), row.names=F)

# run optimization in python via system
run_optimizerPy(python <- '',
                funcs.file_path <- '',
                funcs.file <- '/optimizers.py',
                bounds.file_path <- in_dir,
                bounds.file <- '/space_bounds.csv',
                data.file_path <- out_dir,
                data.loc <- paste0('/data_train_', jobnum, '.csv'),
                optimizer <- paste0(opt_type),
                learner <- 'brtR',
                cv_folds <- '10',
                n_calls <- '150',
                jobnum <- jobnum,
                col_start <- toString(which(names(data_train) == cov_names[1])))

# delete extraneous file
file.remove(paste0(out_dir, '/data_train_', jobnum, '.csv'))

## RUN BRT MODEL
# read in optimized hyperparameters
par <- read.csv(paste0(out_dir, '/best_pars/best_pars_', jobnum, '.csv'))

# get performance for selected hyperparameter values on training set
model_train <- run_brt_model(data_train, par, covs, cov_names, final=F)
save(model_train$stats, file = paste0(out_dir,'/stats_train_output/stats_train_', jobnum,'.Rdata'))

# get performance for selected hyperparameter values on hold out set
model_test <- run_brt_model(data_test, par, covs, cov_names, final=F)
save(model_test$stats, file = paste0(out_dir,'/stats_test_output/stats_test_', jobnum,'.Rdata'))

# get model outputs and performance for selected hyperparameter values on whole set
model <- run_brt_model(dat_all, par, covs, cov_names, final=T)
save(model$model.list, file = paste0(out_dir,'/model_output/model_', jobnum,'.Rdata'))
save(model$stats, file = paste0(out_dir,'/stats_output/stats_', jobnum,'.Rdata'))
writeRaster(model$pred.raster,
            paste0(out_dir, '/preds_output/preds_', jobnum, '.tif'),
            format = 'GTiff', overwrite=T)

