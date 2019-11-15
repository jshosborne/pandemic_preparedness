# code for computing pair-wise distance of point long-lats for each modelling experiment
# - each result is the value to buffer each poly_type with

# convert to radians
deg2rad <- function(deg) return(deg*pi/180)

# calculate the geodesic distance between two points (coords in rads) using Haversine formula
gcd.hf <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1, sqrt(a)))
  d = (R * c) * 1000
  return(d) # Distance in meters
}

# calculate matrix of distances between each two sites
# longlats - data.frame w/ columns c("long","lat")
CalcDists <- function(longlats) {
	name <- list(rownames(longlats), rownames(longlats))
	n <- nrow(longlats)
	# a matrix containing linear distances between all pairwise sites (in meters)
	z <- matrix(0, n, n, dimnames = name)
    for (i in 1:n) {
        for (j in 1:n) {
			z[i, j] <- gcd.hf(long1 = deg2rad(longlats[i, 2]), 
							  lat1 = deg2rad(longlats[i, 1]), long2 = deg2rad(longlats[j, 2]), 
							  lat2 = deg2rad(longlats[j, 1]))
		}
    }
    return(z)
}

# find the maximum minimum pair-wise distance (in m) between all point data
getMaxMin <- function(longlats) {
	dists <- CalcDists(longlats)
	mins <- c()
	n <- ncol(dists)
	for (i in 1:n) {
		colm <- dists[,i]
		colm[colm == 0] <- NA
		colm <- na.omit(colm)
		mins[i] <- min(colm)
	}
	return(max(mins))
}

# location of modelling experiment datasets
# out_dir <- '/snfs1/WORK/11_geospatial/zoonotic_disease/02_henipavirus/01_data/09_modelling_datasets/'
out_dir <- commandArgs()[4]
in_dir <- commandArgs()[5]
f_name <- commandArgs()[6]
dat_files <- list.files(in_dir)
dat_files <- dat_files[grep('experiment', dat_files)]

# calculate and write m.m.p-w. distances to file
buffer_dists <- data.frame(matrix(NA, nrow = 1, ncol = length(dat_files)))
for (i in 1:length(dat_files)) {
	expmt <- read.csv(paste0(in_dir, '/', dat_files[i]))
	expmt.pts <- expmt[expmt$shape_type=='point',]
	buffer_dists[1,i]  <- getMaxMin(na.omit(expmt.pts[,c("long","lat")]))
	names(buffer_dists)[i] <- substr(dat_files[i], 1, 11)
}

if (f_name!='none'){
  write.csv(buffer_dists, paste0(out_dir, '/buffer_distances_', f_name, '.csv'), row.names=F)
} else {
  write.csv(buffer_dists, paste0(out_dir, '/buffer_distances.csv'), row.names=F)
}
