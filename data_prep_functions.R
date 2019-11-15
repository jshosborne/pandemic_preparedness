package_list <- c('seegSDM', 'rgdal', 'maptools', 'ggplot2', 'sf', 'sp', 'raster', 'rgeos', 'dismo')
lapply(package_list, library, character.only = TRUE)

## write_subsets(): function for partitioning model subset data by map layers
## returns: NULL, writes spatial shape objects to GeoTIFF
## params:
## poly_shp  <- experimental subset polygon data shapefile
## buff_shp  <- experimental subset buffer data shapefile
## expmt_dat <- experimental subset data.frame
## crs       <- a CRS() object, should match LBD standard crs
## dat_cols  <- list of strings; column names to retain in output shapefiles for symbolization
## out_dir   <- output directory to save shapefiles to
## shp_dir   <- directory to shapefiles used in model run of experimental subset

write_subsets <- function(poly_shp, buff_shp, expmt_dat, crs, dat_cols, out_dir, shp_dir){
  # POLYGONS: check for polygon data and parse out by admin level; combine with experiment data and save
  if (nrow(poly_shp)!=0){
    # get matching rows
    poly_rows <- which(expmt_dat$occ_id %in% poly_shp$occ_id)
    
    # subset data
    expmt_dat.poly <- expmt_dat[poly_rows, ]
    
    # combine subset polygon shapefile and data columns
    row.names(expmt_dat.poly) <- row.names(poly_shp)
    poly_dat <- spCbind(poly_shp, expmt_dat.poly)
    
    # subset to each admin level layer of the occurrence map and save separately for...
    
    # ... admin 0
    adm0_idx <- which(poly_dat$poly_field == 'ADM0_CODE')
    if (length(adm0_idx)!=0){
      poly_dat.adm0 <- spCbind(poly_dat[adm0_idx, ], expmt_dat.poly[adm0_idx, dat_cols])
      writeSpatialShape(poly_dat.adm0, paste0(out_dir, '/experiment', expmt_idx, '_adm0.shp'))
    }
    # ... admin 1
    adm1_idx <- which(poly_dat$poly_field == 'ADM1_CODE')
    if (length(adm1_idx)!=0){
      poly_dat.adm1 <- spCbind(poly_dat[adm1_idx, ], expmt_dat.poly[adm1_idx, dat_cols])
      writeSpatialShape(poly_dat.adm1, paste0(out_dir, '/experiment', expmt_idx, '_adm1.shp'))
    }
    # ... admin 2 
    adm2_idx <- which(poly_dat$poly_field == 'ADM2_CODE')
    if (length(adm2_idx)!=0){
      poly_dat.adm2 <- spCbind(poly_dat[adm2_idx, ], expmt_dat.poly[adm2_idx, dat_cols])
    }
    # ... custom polygons
    cust_idx <- which(poly_dat$poly_field == 'GAUL_CODE')
    if (length(cust_idx)!=0){
      poly_dat.custom <- spCbind(poly_dat[cust_idx, ], expmt_dat.poly[cust_idx, dat_cols])
      if (length(adm2_idx)!=0) poly_dat.adm2_custom <- bind(poly_dat.adm2, poly_dat.custom)
      writeSpatialShape(poly_dat.adm2_custom, paste0(out_dir, '/experiment', expmt_idx, '_adm2_custom.shp'))
    } else if (length(adm2_idx)!=0) {
      writeSpatialShape(poly_dat.adm2, paste0(out_dir, '/experiment', expmt_idx, '_adm2.shp'))
    }
  } 
  # BUFFERS: check for buffer data polygons; combine with experiment data and save
  if (nrow(buff_shp)!=0){
    # get matching rows
    buff_rows <- which(expmt_dat$occ_id %in% buff_shp$occ_id)
    
    # subset data
    expmt_dat.buff <- expmt_dat[buff_rows, append(dat_cols, c('buffer_radius', 'long', 'lat'))]
    expmt_dat.buff$buffer_radius <- 1000 * expmt_dat.buff$buffer_radius
    
    # combine subset buffer shapefile and data columns
    row.names(expmt_dat.buff) <- row.names(buff_shp)
    buff_dat <- spCbind(buff_shp, expmt_dat.buff)
    writeSpatialShape(buff_dat, paste0(out_dir, '/experiment', expmt_idx, '_buffers.shp'))
  }
  # POINTS: check for point data; retain experiment data and save as spatial object
  point_idx <- which(expmt_dat$shape_type == 'point')
  if (length(point_idx)!=0){
    expmt_dat.pt <- expmt_dat[point_idx,]
    # convert point data to spatial object with same CRS as polygons
    expmt_dat.pt <- SpatialPointsDataFrame(coords=expmt_dat[point_idx, c('long', 'lat')],
                                           data=data.frame(expmt_dat[point_idx, dat_cols]),
                                           proj4string=crs)
    writeSpatialShape(expmt_dat.pt, paste0(out_dir, '/experiment', expmt_idx, '_points.shp'))
  }
}

## buffers.data_buffer(): function to draw buffered data regions and buffer them by buffer_dist
## returns: tuple of buffered region spatial shapes and buffers of buffered region spatial shapes
## params:
## expmt.buff   -> data.frame of centroids of buffered regions of occurrence
## buffer_dist  -> distance to buffer the buffered regions by
## area_max     -> calculated maximum area tolerance (for sampling/boostrap efficiency)
## land         -> shapefile defining coastlines
## crs          -> coordinate reference system, standardize to convert to/from for buffering in units of meters

buffers.data_buffer <- function(expmt.buffs, buffer_dist, area_max, land, crs){
  buffs_list <- c()
  buffs_areas <- c()
  for (i in 1:nrow(expmt.buffs)) {
    loc <- data.frame(expmt.buffs[i, c('long', 'lat')])
    rad <- as.numeric(expmt.buffs$buffer_radius[i])* 1000 # standard radius unit is km -> convert to meters
    # create buffer
    loc.sp <- SpatialPoints(loc, proj4string = crs)
    buffer <- buffer(loc.sp, width=rad)
    buffer_area <- gArea(spTransform(buffer, CRS=CRS("+proj=merc +ellps=GRS80"))) * 1e-6 # convert to meters^2
    buffs_list <- append(buffs_list, buffer)
    buffs_areas <- append(buffs_areas, buffer_area)
  }
  buffs_list <- buffs_list[which(buffs_areas <= area_max)]
  
  # merge buffered regions
  expmt.buff_shp <- do.call(bind, buffs_list)
  expmt.buff_shp <- as(expmt.buff_shp, 'SpatialPolygonsDataFrame')
  expmt.buff_shp <- spCbind(expmt.buff_shp, expmt.buffs[which(buffs_areas <= area_max), 'occ_id'])
  names(expmt.buff_shp)[ncol(expmt.buff_shp)] <- 'occ_id'
  
  expmt.buff_shp <- intersect(expmt.buff_shp, land)
  
  # buffer the buffer data by buffer_dist
  buffs.buffer <- spTransform(buffer(spTransform(expmt.buff_shp, CRS=CRS("+proj=merc +ellps=GRS80")), buffer_dist), CRS=crs)
  
  return(list(data=expmt.buff_shp, buffer=buffs.buffer))
}

## combine.polys_buffs(): function to read in polygon data from shapefile and buffer it by buffer_dist
##   -- runs best in parallel
## returns: tuple of polygonal region spatial shapes and buffers of polygonal region spatial shapes
## params:
## expmt.polys  -> spatial shape of polygonal regions of occurrence
## buffer_dist  -> distance to buffer the buffered regions by
## crs          -> coordinate reference system, standardize to convert to/from for buffering in units of meters

polygons.data_buffer <- function(expmt.polys, buffer_dist, crs){
  # convert filepath from windows to unix convention (for cluster use)
  poly_path <- path_converter(expmt.polys$poly_reference[poly_idx], "J")
  # read shapefile
  shapeDF <- readOGR(poly_path, stringsAsFactors=F)
  # get polygon id (not necessarily Gaul), make sure it's a number
  currGaul <- as.numeric(expmt.polys$poly_id[poly_idx])
  # get field name where poly id is to be found and makes ure it's a string
  poly_field <- toString(expmt.polys$poly_field[poly_idx])
  # subset shapefile (a SpatialPolygonsDataFrame object) to specified polygon
  poly <- shapeDF[shapeDF[[poly_field]] == currGaul,]
  # calculate area of equal-area-projected subset and convert from m^2 to km^2
  poly_area <- gArea(spTransform(poly, CRS=CRS("+proj=merc +ellps=GRS80"))) * 1e-6
  # save occ_id for troubleshooting
  occ_id <- expmt.polys$occ_id[poly_idx]
  
  poly_final <- spCbind(spCbind(poly, poly_area), occ_id)
  
  poly_final.buffer <- spTransform(buffer(spTransform(poly_final, CRS=CRS("+proj=merc +ellps=GRS80")), buffer_dist), CRS=crs)
  
  return(list(data=poly_final, buffer=poly_final.buffer))
}

## combine.polys_buffs(): function to combine combine.polys_buffs() outputs
## returns: tuple of polygonal region spatial shapes and buffers of polygonal region spatial shapes
## params:
## expmt.polys  -> spatial shape of polygonal regions of occurrence
## area_max     -> calculated maximum area tolerance (for sampling/boostrap efficiency)
## temp_dir     -> directory to polygons.data_buffer parallelized outputs

combine.polys_buffs <- function(expmt.polys, area_max, temp_dir){
  polys_list <- c()
  polys_areas <- c()
  poly_buffs_list <- c()
  for (i in 1:nrow(expmt.polys)){#as.numeric(jobnum)){
    poly_file <- paste0(temp_dir, '/poly_', i, '.shp')
    poly <- shapefile(poly_file)
    polys_list <- append(polys_list, poly)
    polys_areas <- append(polys_areas, poly$poly_area)
    buff_file <- paste0(temp_dir, '/buff_', i, '.shp')
    poly_buffs_list <- append(poly_buffs_list, shapefile(buff_file))
  }
  
  # merge shapefiles together
  expmt.poly_shp <- do.call(bind, polys_list)
  expmt.poly_shp <- expmt.poly_shp[which(polys_areas <= area_max),]
  poly_buffs_list <- poly_buffs_list[which(polys_areas <= area_max)]
  polys.buffer <- do.call(bind, poly_buffs_list)
  polys.buffer <- gUnaryUnion(polys.buffer)
  
  return(list(data=expmt.poly_shp, buffer=polys.buffer))
}

## combine.all(): function to merge all background sampling shapefile regions
## returns: tuple of background sampling region with and without spatial shape occurrence sampling regions removed
## params:
## pts.buffer     -> spatial shape of point occurrences background sampling region
## buffs.buffer   -> spatial shape of buffer occurrences background sampling region
## polys.buffer   -> spatial shape of polygon occurrences background sampling region
## expmt.poly_shp -> spatial shape of polygon occurrences
## expmt.buff_shp -> spatial shape of buffer occurrences
## expmt.pts_sp   -> spatial shape of point occurrences

combine.all <- function(pts.buffer, buffs.buffer, polys.buffer, expmt.poly_shp, expmt.buff_shp, expmt.pts_sp){
  # combine all buffers and dissolve boundaries between them
  all.buffers <- do.call(bind, c(polys.buffer, pts.buffer, buffs.buffer))
  all.buffers <- gUnaryUnion(all.buffers)
  
  # retain intersection with land (drop buffered regions over water)
  all.buffers <- intersect(all.buffers, land)
  all.buffers_smooth <- all.buffers
  
  # erase features occurrence samples are being drawn from
  for (i in 1:nrow(expmt.poly_shp)){
    all.buffers <- erase(all.buffers, expmt.poly_shp[i,])
  }
  all.buffers <- erase(all.buffers, expmt.buff_shp)
  pts.intr_buffer <- circles(expmt.pts_sp, lonlat=TRUE, d=2500)@polygons
  all.buffers <- erase(all.buffers, pts.intr_buffer)
  
  return(list(smooth=all.buffers_smooth, final=all.buffers))
}
