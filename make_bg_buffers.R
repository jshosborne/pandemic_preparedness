# code for making objects to buffer for each modelling experiment

out_dir <- commandArgs()[4]
in_dir <- commandArgs()[5]
expmt_idx <- as.numeric(commandArgs()[6])
area_max <- as.numeric(commandArgs()[7]) # in km^2, 1 pixel = 25 km^2
cov_pctg <- commandArgs()[8]
f_name <- commandArgs()[9]

source('helper_functions.R')
package_list <- c('seegSDM', 'rgdal', 'maptools', 'ggplot2', 'sf', 'raster')
lapply(package_list, library, character.only = TRUE)

land <- shapefile("land.shp")
crs <- CRS(paste0(proj4string(land))) # define common crs to be mbg standard crs

# location of modelling experiment datasets
dat_files <- list.files(in_dir)
dat_files <- dat_files[grep('experiment', dat_files)]
if (f_name!='none'){
  distances <- read.csv(paste0(out_dir, '/buffer_distances_', f_name, '.csv'))
} else {
  distances <- read.csv(paste0(out_dir, '/buffer_distances.csv'))
}
expmt_file <- dat_files[expmt_idx]

# subset data
expmt <- read.csv(paste0(in_dir, '/', expmt_file))
expmt.pts <- expmt[expmt$shape_type=='point',]

buffer_dist <- distances[grep(substr(expmt_file, 1, 11), names(distances))][[1]]

expmt.pts_sp <- SpatialPoints(expmt.pts[,c('long', 'lat')], proj4string = crs)

pts.buffer <- buffer(expmt.pts_sp, width=buffer_dist)

expmt.buffs <- expmt[expmt$poly_type=='buffer',]
expmt.polys <- expmt[expmt$shape_type=='polygon',]
expmt.polys <- expmt.polys[expmt.polys$poly_type!='buffer',]
expmt.polys <- expmt.polys[expmt.polys$poly_type != "",]
expmt.polys <- expmt.polys[!is.na(expmt.polys$poly_reference),]

dat.buff <- buffers.data_buffer(expmt.buffs, buffer_dist, area_max, land, crs)
expmt.buff_shp <- dat.buff$data
buffs.buffer <- dat.buff$buffer

# write to file for occurrence sampling w/in bootstraps
if (f_name!='none'){
  setwd(paste0(out_dir, '/sampling_shapefiles_', cov_pctg, '_', f_name))
} else {
  setwd(paste0(out_dir, '/sampling_shapefiles_', cov_pctg))
}
writeOGR(as(expmt.buff_shp, 'SpatialPolygonsDataFrame'), 
	       ".", 
	       paste0(substr(expmt_file, 1, 11), "_buffer_data"),
	       "ESRI Shapefile", 
	       overwrite_layer=T)
  

# parallelize   
poly.run <- 'save_buffered_polys.R'
occ_id_list <- c()
i_list <- c()
for (i in 1:nrow(expmt.polys)){
  occ_id_list <- append(occ_id_list, expmt.polys$occ_id[i])
  i_list <- append(i_list, i)
##
# parallelize save_buffered_polys here
##
}
  
if (f_name!='none') {
  temp_dir <- paste0('', f_name, '_', expmt_idx)
} else {
  temp_dir <- paste0('', expmt_idx)
}
  
com <- combine.polys_buffs(expmt.polys, area_max, temp_dir)
polys.buffer <- com$buffer
expmt.poly_shp <- com$data

# combine all buffers and dissolve boundaries between them
all <- combine.all(pts.buffer, buffs.buffer, polys.buffer, expmt.poly_shp, expmt.buff_shp, expmt.pts_sp)
all.buffers <- all$final
smooth_buffer <- all$smooth
  
# write to file for background sampling w/in bootstraps
writeOGR(as(all.buffers, 'SpatialPolygonsDataFrame'), 
     	   ".",
	       paste0(substr(expmt_file, 1, 11), "_final_buffer"),
	       "ESRI Shapefile",
	       overwrite_layer=T)

writeOGR(as(smooth_buffer, 'SpatialPolygonsDataFrame'), 
         ".",
         paste0(substr(expmt_file, 1, 11), "_smooth_buffer"),
         "ESRI Shapefile",
         overwrite_layer=T)


# clear out temp_dir once individual files are not needed
do.call(file.remove, list(list.files(temp_dir, full.names = TRUE)))
