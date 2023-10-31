# running the CA from only user input of ignition lon lat OR existing fire (VIIRS)
# assume that the ignition is within the domain
# read the entire domain. know its extent in lat lon
# sample a square around the ignition point OR existing fire (how big of square?)
# generate initial state raster for CA
# write rasters to CA input folder
# run CA and save images to output folder
# compile images into GIF somehow

###################
# retrieve user input variables from Database.
library("RPostgreSQL", lib.loc="/usr/lib64/R/library")
library(rpostgis)
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname="fire_activity",user="topofire",password="topofire")
# Just checking the connection, don't need this really..
pgPostGIS(con)
print("Searching Job Queue...")
# Grab the huc polygon
# Grab the record info
query<-paste("select id,lon,lat from firespread_jobs where status = 'Pending' limit 1;")
rowset<-dbSendQuery(con,query)
# the query results
results<-fetch(rowset,n=0)

row_id<-results[1,1]
huc<-results[1,2]
huc_name<-results[1,3]
user_email<-results[1,4]
status<-results[1,5]

# the query results: If no pending jobs found, close program. else switch to running
if (is.null(row_id)){
  # handle empty result set
  print("No jobs found..")
  quit()
}else{
  print("Job found..")
  query<-paste0("update firespread_jobs set status = 'Running' where id =",row_id,";")
  rowset<-dbSendQuery(con,query)
}
# close database connection here
dbDisconnect(con)

ign_loc <- cbind(results[2],results[3])

require(terra)

# IMPORTANT: if you're using terra and large rasters, careful of memory
# mak a temporary directory for R:terra like this:
terraOptions(tempdir="/home/adam.viray/tmprast") # you need to make this dir. yourself

## read in data for extent
print("reading data")
# data locations
spat.dir <- "/mnt/DataDrive1/data/firespread_CA/spatial_inputs/"
wind.dir <- "/mnt/cphrstor2DataDrive1/NDFD/hourly_wind/"
pred.dir <- "/mnt/cphrstor2DataDrive1/R1_fire/daily_predict/"

# topography
slope <- terra::rast(paste0(spat.dir,"r1_slope_8s.tif"))
aspect <- terra::rast(paste0(spat.dir,"r1_aspect_8s.tif"))
dem <- terra::rast(paste0(spat.dir,"r1_dem_8s.tif"))

# vegetation
noveg <- paste0(spat.dir,"rap_2_noveg_2022_8s.tif")
forest <- paste0(spat.dir,"rap_6_forest_2022_8s.tif")
shrub <- paste0(spat.dir,"rap_5_shrub_2022_8s.tif")
veg <- terra::rast(c(noveg,shrub,forest))


# wind
a = nchar(paste0(in.dir,"hourly_wind/"))
wd.list <- list.files(wind.dir, pattern="WindDir", full.names = T)
wd.dates <- as.POSIXct(substr(wd.list, a+9, a+20), format="%Y%m%d%H%M", tz="MST")
wdir <- terra::project(terra::rast(wd.list), slope, method="near")

ws.list <- list.files(wind.dir, pattern="WindSpd", full.names=T)
ws.dates <- as.POSIXct(substr(ws.list, a+9, a+20), format="%Y%m%d%H%M", tz="MST")
wspd <- terra::project(terra::rast(ws.list), slope, method="bilinear")

# ignition probability
p.list <- list.files(paste0(in.dir,"daily_ignit"), full.names = T)
prob <- terra::rast(p.list)

# rate of growth
r.list <- list.files(paste0(in.dir,"daily_growth"), full.names = T)
RoS <- terra::rast(r.list)

## get smallest number of hours
k = min(c(dim(RoS)[3]*24,dim(prob)[3]*24,dim(wspd)[3],dim(wdir)[3]))
l = ceiling(k/24)

## get domain extent
domain <- terra::ext(slope)

## sample a square around ignition point
w <- max(c(domain[1], ign_loc[1]-0.1422))
e <- min(c(domain[2], ign_loc[1]+0.1422))
s <- max(c(domain[3], ign_loc[2]-0.1422))
n <- min(c(domain[4], ign_loc[2]+0.1422))

extent <- c(w,e,s,n)

##create a temp. directory for running each sim
sim.dir <- "/home/adam.viray/Documents/CA/data/temp/"
if(!dir.exists(sim.dir)) { dir.create(sim.dir) }

## generate initial state raster
state <- terra::crop(terra::rast(slope),extent)
terra::values(state) <- 0
m = round(dim(state)[[1]]/2)+1
state[m,m] <- 2

## write cropped rasters
print("writing data")
terra::writeRaster(terra::crop(prob[[1:l]],extent),paste0(sim.dir,"prob.tif"), overwrite=T)
terra::writeRaster(terra::crop(RoS[[1:l]],extent),paste0(sim.dir,"RoS.tif"), overwrite=T)
terra::writeRaster(terra::crop(slope,extent),paste0(sim.dir,"slope.tif"), overwrite=T)
terra::writeRaster(terra::crop(aspect,extent),paste0(sim.dir,"aspect.tif"), overwrite=T)
terra::writeRaster(terra::crop(dem,extent),paste0(sim.dir,"dem.tif"), overwrite=T)
terra::writeRaster(terra::crop(wdir[[1:k]],extent),paste0(sim.dir,"wdir.tif"), overwrite=T)
terra::writeRaster(terra::crop(wspd[[1:k]],extent),paste0(sim.dir,"wspd.tif"), overwrite=T)
terra::writeRaster(terra::crop(veg,extent),paste0(sim.dir,"veg.tif"), overwrite=T)
terra::writeRaster(state,paste0(sim.dir,"state.tif"), overwrite=T)

## run python CA
gif.dir <- "/home/adam.viray/Documents/CA/output/GIF/"
system(paste("python", "/home/adam.viray/Documents/CA/code/CA.py", sim.dir, gif.dir))

## turn images into gif
fname <- paste0(gsub(", ", "_", toString(ign_loc)),".gif")
system(paste("convert -delay 20 -loop 0", paste0(gif.dir, "*.png"), paste0(gif.dir, fname)))
