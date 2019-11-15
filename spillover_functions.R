package_list <- c('raster', 'maptools', 'sp', 'rgdal', 'rgeos')
lapply(package_list, library, character.only = TRUE)

## standardize_field(): function for producting 1-10 standardized scores for an arbitrary field
## returns: shapefile@data with new column of values
## params:
## paly_val   -> base shapefile@data object
## field      -> shapefile@data column to standardize
## field_name -> string, name of field in returned object

standardize_field <- function(poly_val,
                              field,
                              field_name)
{
  poly_val$ctxt <- field
  
  for (i in 1:nrow(poly_val)) {
    if (is.na(poly_val$ctxt[i])) poly_val$ctxt[i]<-0
    if (poly_val$ctxt[i]==0) {
      poly_val$log_ctxt[i]<-0
    } else {
      poly_val$log_ctxt[i]<-log10(poly_val$ctxt[i])
    }
  }
  
  max_ctxt<-max(poly_val$log_ctxt)
  min_ctxt<-min(poly_val$log_ctxt[which(poly_val$log_ctxt>0)])
  
  for (i in 1:nrow(poly_val)) {
    if (poly_val$ctxt[i] == 0) {
      poly_val$standard_ctxt[i] <- 0
    } else {
      poly_val$standard_ctxt[i]<-((poly_val$log_ctxt[i] - min_ctxt)/(max_ctxt-min_ctxt))*10
    }
  }
  
  for (i in 1:nrow(poly_val)) {
    if (poly_val$standard_ctxt[i] < 0) {
      poly_val$standard_ctxt[i] <- 0
    }
  }
  names(poly_val) <- gsub('ctxt', field_name, names(poly_val))
  return(poly_val)
}

## zone_summary(): function for summarizing (risk) rasters over polygon areas with contextual, rasterized data
##                  -- in addition, rescales values to 1-10 using log transform
## returns: shapefile@data with new column of values
## params:
## measure_raster   -> base raster to do zonal summarization with
## poly_raster      -> rasterized shapefiles defining extent of zones
## poly_data        -> polygon metadata taken from shapefile
## context_raster   -> raster data to compound with measure_raster -- e.g. -- human/host species' populations
## sum_func_meas    -> zonal() parameter value (string) for summarizing measure_raster
## sum_func_ctxt    -> zonal() parameter value (string) for summarizing context_raster*measure_raster

zone_summary <- function(poly_data,
                         measure_raster,
                         poly_raster,
                         context_raster,
                         sum_func_meas='sum',
                         sum_func_ctxt='sum')
{
  #create population living in areas of risk
  context_risk <- context_raster*measure_raster
  #crop to area
  context_risk <- crop(context_risk, poly_raster)
  
  #perform a zonal summary
  zonal_summary <- zonal(context_risk,
                         poly_raster,
                         fun=sum_func_meas)
  zonal_summary <- data.frame(zonal_summary)
  
  #append a null row to dataframe
  poly_val<-data.frame(poly_data, ctxt_abs=rep(0,nrow(poly_data)))
  #merge
  for (i in 1:nrow(zonal_summary)){
    id<-zonal_summary$zone[i]
    poly_val$ctxt_abs[which(poly_val$ADM2_CODE==id)]<-zonal_summary$sum[i]
  }
  
  #rescale all the values where max value = 10, min value = 0, after log transforming
  
  for (i in 1:nrow(poly_val)) {
    if (poly_val$ctxt_abs[i]==0) {
      poly_val$log_ctxt[i]<-0
    } else {
      poly_val$log_ctxt[i]<-log10(poly_val$ctxt_abs[i])
    }
  }
  
  max_ctxt<-max(poly_val$log_ctxt)
  min_ctxt<-min(poly_val$log_ctxt[which(poly_val$log_ctxt>0)])
  
  for (i in 1:nrow(poly_val)) {
    if (poly_val$ctxt_abs[i] == 0) {
      poly_val$standard_ctxt[i] <- 0
    } else {
      poly_val$standard_ctxt[i]<-((poly_val$log_ctxt[i] - min_ctxt)/(max_ctxt-min_ctxt))*10
    }
  }
  
  for (i in 1:nrow(poly_val)) {
    if (poly_val$standard_ctxt[i] < 0) {
      poly_val$standard_ctxt[i] <- 0
    }
  }
  
  #calculate total populations - have to go beyond buffer as some districts are partially in buffer
  cropped_context <- crop(context_raster, poly_raster)
  zonal_context <- zonal(cropped_context,
                         poly_raster,
                         fun=sum_func_ctxt)
  zonal_context <- data.frame(zonal_context)
  
  #append total pops to the end of the row NOTE that some locations are smaller than 5km x 5km so are currently missing
  poly_val$total_ctxt <- rep(0, nrow(poly_val))
  for(i in 1:nrow(poly_val)){
    id <- zonal_context$zone[i]
    poly_val$total_ctxt[which(poly_val$ADM2_CODE == id)] <- zonal_context$sum[i]
  }
  
  #evaluate the proportion of the total population exposed in the district
  poly_val$prop_ctxt<-rep(0, nrow(poly_val))
  for(i in 1:nrow(poly_val)){
    if(poly_val$total_ctxt[i] == 0){
      poly_val$prop_ctxt[i] <- 0
    } else {
      poly_val$prop_ctxt[i] <- poly_val$ctxt_abs[i]/poly_val$total_ctxt[i]
    }
  }
  
  #standardize the proportion
  for(i in 1:nrow(poly_val)){
    if (poly_val$prop_ctxt[i] == 0){
      poly_val$standard_prop_ctxt[i] <- 0
    } else {
      poly_val$standard_prop_ctxt[i] <- ((poly_val$prop_ctxt[i] - 0)/(1-0))*10
    }
  }
  
  #take the geometric mean of proportion population and standardized total population
  for(i in 1:nrow(poly_val)){
    poly_val$geom_ctxt[i]<-sqrt((poly_val$standard_ctxt[i]*poly_val$standard_prop_ctxt[i]))
  }
  
  names(poly_val) <- gsub('ctxt', names(context_raster), names(poly_val))
  return(poly_val)
}

#get population total per quintile in spillover maps
get_pop_quint <- function(rr_sf, id_col, fill_col){
  # identical quintile creation from get_rel_risk_map()
  incl <- which(rr_sf@data[,fill_col] != 0)
  rr_sf_na_excl <- rr_sf[incl,]
  spill_over_df <- fortify(rr_sf_na_excl, region = id_col)
  metad <- rr_sf_na_excl@data
  spill_over_df <- merge(spill_over_df, metad, by.x = 'id', by.y = id_col)
  quint.int <- quantile(spill_over_df[,fill_col], probs=seq(0, 1, by = 1/5))
  
  # get quintile population totals
  population_at_risk <- c()
  for (i in 1:5){
    quint_idx <- which(rr_sf[[fill_col]] <= quint.int[[i+1]] & rr_sf[[fill_col]] > quint.int[[i]])
    quint_pop <- sum(rr_sf[quint_idx,]$population_abs)
    population_at_risk <- append(population_at_risk, quint_pop)
  }
  
  pop_tot <- data.frame(population_at_risk)
  row.names(pop_tot) <- c('quintile1','quintile2','quintile3','quintile4','quintile5')
  pop_tot <- rbind(pop_tot, total = sum(pop_tot$population_at_risk))
  
  return(pop_tot)
}

# geometric mean function,
gm_mean = function(a){prod(a)^(1/length(a))}

## spillover_map(): function for computing spillover/detection potential of administrative units
## returns: shapefile with new column of values
## params:
## full_measure_raster   -> base raster 
## mask                  -> mask to apply to zonal summary
## rast_stack            -> stack of context rasters to perform zonal summary
## rast_stack_idx        -> stack element indices of context rasters
## meas_func             -> zonal() parameter value (string) for summarizing measure_raster
## ctxt_func             -> zonal() parameter value (string) for summarizing context_raster*measure_raster
## threshold             -> ROC optimized threshold for env. suitabilit predictions, defines pixels to use in zonal_summary
## admin2_raster         -> raster of 2nd level administrative units
## admin2_shp            -> shapefile of 2nd level administrative units
## all_shp               -> shapefile of regions to exclude from assessment

spillover_map <- function(full_measure_raster, 
                          mask=NULL, 
                          rast_stack,
                          rast_stack_idx,
                          meas_func, 
                          ctxt_func, 
                          threshold,
                          admin2_raster,
                          admin2_shp,
                          all_shp=NULL){
  #load in relevant mask
  if (!is.null(mask)) {
    #rasterize it
    if(toString(class(mask))!="RasterLayer") mask_raster <- rasterize(mask, full_measure_raster)
    measure_raster <- mask_raster*full_measure_raster
  } else {
    measure_raster <- full_measure_raster
  }
  
  #threshold value for binary distinctions
  f.list <- list.files(out_dir, full.names=T)
  threshold<-read.csv(f.list[grep("thresh_", f.list)])$opt_thresh
  
  #convert suitability raster to binary raster using threshold
  measure_raster[measure_raster >= threshold] <- 1
  measure_raster[measure_raster < threshold] <- 0
  
  #create known region mask
  if (!is.null(all_shp)){
    poly_rast <- rasterize(all_shp, full_measure_raster)
    poly_rast[poly_rast > 0] <- 0
    poly_rast[is.na(poly_rast)] <- 1
  }
  
  admin_val <- admin2_shp@data
  poly_val <- admin_val
  
  for (j in rast_stack_idx){
    if (!is.null(all_shp)) {
      poly_val <- zone_summary(poly_val,
                               measure_raster*poly_rast, 
                               admin2_raster,
                               rast_stack[[j]],
                               meas_func[j],
                               ctxt_func[j])
    } else {
      poly_val <- zone_summary(poly_val,
                               measure_raster,
                               admin2_raster,
                               rast_stack[[j]],
                               meas_func[j],
                               ctxt_func[j])
    }
  }
  
  # standardize diabetes and kidney disease prevalence
  poly_val <- standardize_field(poly_val, comorb_shp$B8_prev, 'dbkd')
  
  # standardize cardiovascular disease prevalence
  poly_val <- standardize_field(poly_val, comorb_shp$B2_prev, 'card')
  
  # calculate geometric mean of two factors
  poly_val$geom_hum_cam <- gm_mean(c(poly_val$geom_camel_density,poly_val$geom_population))
  
  # calculate geometric mean of four factors
  poly_val$geom_all <- gm_mean(c(poly_val$geom_camel_density, poly_val$geom_population, poly_val$card, poly_val$dbkd))
  
  #update shapefile vals
  admin2_out <- admin2_shp 
  admin2_out@data <- poly_val
  
  return(admin2_out)  
}
