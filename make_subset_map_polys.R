# script to parse experiment subsets into mappable shapefiles at each admin scale

# read in script variables from qsub
in_dir <- commandArgs()[4]
expmt_idx <- as.numeric(commandArgs()[5])
f_name <- commandArgs()[6]
expmt_dir <- commandArgs()[7]
cov_pctg <- commandArgs()[8]

# load in libraries
source('data_prep_functions.R')
package_list <- c('seegSDM', 'rgdal', 'maptools', 'sp')
lapply(package_list, library, character.only = TRUE)

# make directories based on f_name
if (f_name!='none'){
  out_dir <- paste0(in_dir, '/experiment_mapping_', f_name)
  shp_dir <- paste0(in_dir, '/sampling_shapefiles_', cov_pctg, '_', f_name)
} else {
  out_dir <- paste0(in_dir, '/experiment_mapping')
  shp_dir <- paste0(in_dir, '/sampling_shapefiles_', cov_pctg)
}

# make sure path to new shapefiles exists
if (!dir.exists(out_dir)) dir.create(out_dir)

# select the shapefiles to use from sampling_shapefiles
shp_file_list <- list.files(shp_dir, full.names = T)
poly_file <- shp_file_list[grep(paste0('experiment', expmt_idx, '_polygon_data.shp'), shp_file_list)]
buff_file <- shp_file_list[grep(paste0('experiment', expmt_idx, '_buffer_data.shp'), shp_file_list)]

# load in experiment shapefile for polygon and buffer data
poly_shp <- shapefile(poly_file)
buff_shp <- shapefile(buff_file)

# combine buffer and polygon data
all_shp <- bind(poly_shp, buff_shp)

# select data from list of subsets
expmt_file_list <- list.files(expmt_dir, full.names = T)
expmt_file <- expmt_file_list[expmt_idx]

# load in experiment data
expmt_dat <- read.csv(expmt_file, as.is=T)

# replace missing values with organism_type values
for (i in 1:nrow(expmt_dat)){
  if (is.na(expmt_dat$patient_type[i])) expmt_dat$patient_type[i] <- expmt_dat$organism_type[i]
}

# define standard coordinate reference system
crs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# what values should be retained? -- must match input data
dat_cols <- c('occ_id', 'patient_type', 'pathogen')

# create experiment subset shapefiles
write_subsets(poly_shp, buff_shp, expmt_dat, crs, dat_cols, out_dir, shp_dir)


