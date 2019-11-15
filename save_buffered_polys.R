in_dir <- commandArgs()[4]
out_dir <- commandArgs()[5]
poly_idx <- as.numeric(commandArgs()[6])
expmt_idx <- as.numeric(commandArgs()[7])
f_name <- commandArgs()[8]
cov_pctg <- commandArgs()[9]

source('helper_functions.R')
source('data_prep_functions.R')
package_list <- c('seegSDM', 'rgdal', 'maptools')
lapply(package_list, library, character.only = TRUE)

land <- shapefile("land.shp")
crs <- CRS(paste0(proj4string(land))) # define common crs to be mbg standard crs

if (f_name!='none') {
  temp_dir <- paste0('temp_shapefiles_', f_name, '_', expmt_idx)
} else {
  temp_dir <- paste0('temp_shapefiles_', expmt_idx)
}
if(!dir.exists(temp_dir)) dir.create(temp_dir)
setwd(temp_dir)

dat_files <- list.files(in_dir)
dat_files <- dat_files[grep('experiment', dat_files)]
expmt_file <- dat_files[expmt_idx]

if (f_name!='none'){
  distances <- read.csv(paste0(out_dir, '/buffer_distances_', f_name, '.csv'))
} else {
  distances <- read.csv(paste0(out_dir, '/buffer_distances.csv'))
}
buffer_dist <- distances[grep(substr(expmt_file, 1, 11), names(distances))][[1]]

# subset data
expmt <- read.csv(paste0(in_dir, '/', expmt_file))
expmt.polys <- expmt[expmt$shape_type=='polygon',]
expmt.polys <- expmt.polys[expmt.polys$poly_type!='buffer',]
expmt.polys <- expmt.polys[expmt.polys$poly_type != "",]
expmt.polys <- expmt.polys[!is.na(expmt.polys$poly_reference),]

# read individual data polygons and buffer them
pbuff <- polygons.data_buffer(expmt.polys, buffer_dist, crs)
poly_final.buffer <- pbuff$buffer
poly_final <- pbuff$data

# write out polygons and their buffers
writeOGR(as(poly_final, 'SpatialPolygonsDataFrame'), 
         ".", 
         paste0("poly_", poly_idx),
         "ESRI Shapefile", 
         overwrite_layer=T)

writeOGR(as(poly_final.buffer, 'SpatialPolygonsDataFrame'), 
         ".", 
         paste0("buff_", poly_idx),
         "ESRI Shapefile", 
         overwrite_layer=T)