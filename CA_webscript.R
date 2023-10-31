# running the CA from only user input of ignition lon lat OR existing fire (VIIRS)
# assume that the ignition is within the domain
# read the entire domain. know its extent in lat lon
# sample a square around the ignition point OR existing fire (how big of square?)
# generate initial state raster for CA
# write rasters to CA input folder
# run CA and save images to output folder
# compile images into GIF somehow

## constants
# data extent
WEST <- -120.2
EAST <- -105.000000007141
SOUTH <- 43.5000000022993
NORTH <- 49.0000000003689
SPATIAL.RES <- 60
base.radius <- 0.008

###################
# retrieve user input variables from Database.quit()
#library("RPostgreSQL")#, lib.loc="/usr/lib64/R/library")
#library(rpostgis)
#drv <- dbDriver("PostgreSQL")
#con <- dbConnect(drv, dbname="fire_activity",user="topofire",password="topofire")
# Just checking the connection, don't need this really..
#pgPostGIS(con)
#print("Searching Job Queue...")
# Grab the ignition lon lat
#query<-paste("select id,lon,lat,user_email,firename, start_year from firespread_jobs where status = 'Pending' limit 1;")
#rowset<-dbSendQuery(con,query)
# the query results
#results<-fetch(rowset,n=0)
#print(results)
#
#row_id <- results[1,1]
#user_email <- results[1,4]
#firename=chartr(" ","_",results[1,5])
#start_date=results[1,6]
#print(paste("job number", row_id))
#lon <- results[1,2]
#lat <- results[1,3]
#print(lon); print(lat)
# check ignition point in extent
#if(lon < WEST || lon > EAST || lat < SOUTH || lat > NORTH) {
#  print("coordinate outside extent")
#  quit()
#}

row_id = 0
user_email = "adam.viray@gmail.com"
firename = "test"
start_date = "2023-09-11"
lon = -113.5862
lat = 47.91082

## create log file
log.dir <- "/mnt/DataDrive1/data/firespread_CA/logfiles/"
today <- format(as.Date(Sys.time()),"%Y%m%d")
logname <- paste0(log.dir, firename, "_", today, "_log.txt")
zz <- file(logname,open="wt")
sink(zz,split=T);sink(zz,type="message")
on.exit({sink(type="message");sink()})

print("log file created")
print(paste("spatial res. ", SPATIAL.RES))
# define start date as 1 day (86400 seconds) behind when the sim. gets called. 
if(!exists("start_date")) {
start_date=substr(Sys.time()-86400,1,10)
}

# format the date to match the file name dates. 
start_date2=format(as.Date(start_date), "%Y%m%d")
print(paste("sim. start date ", start_date2))

# check if user input name of fire
#if(length(firename)!=1) {firename <- paste0(round(lat,3),"_",round(lon,3))}
if(!exists("firename")) {   firename <- paste0(round(lat,3),"_",round(lon,3)) }
print(firename)

# the query results: If no pending jobs found, close program. else switch to running
#if (is.null(row_id)){
#  # handle empty result set
#  print("No jobs found..")
#  quit()
#}else{
#  print("Job found..")
#  query<-paste0("update firespread_jobs set status = 'Running' where id =",row_id,";")
#  rowset<-dbSendQuery(con,query)
#}
# close database connection here
#dbDisconnect(con)

## read in data for extent
print("building data")
# data locations
spat.dir <- "/mnt/DataDrive1/data/firespread_CA/spatial_inputs/"
wind.dir <- "/mnt/cphrstor2DataDrive1/NDFD/hourly_wind/"
pred.dir <- "/mnt/cphrstor2DataDrive1/R1_fire/daily_predict/"
out.dir <- "/mnt/DataDrive1/data/firespread_CA/output/"
##create a temp. directory for running each sim
usr.dir <- paste0(out.dir, user_email,"/")
if(!dir.exists(usr.dir)) { dir.create(usr.dir) }
sim.dir <- paste0(usr.dir,row_id,"_",firename,"_",start_date2,"-",lon,"_",lat,"_",SPATIAL.RES,"m/")
if(!dir.exists(sim.dir)) { dir.create(sim.dir) }


## make simulation domain
## IMPORTANT: if you're using terra and large rasters, careful of memory
# mak a temporary directory for R:terra like this:

require(terra)
terraOptions(tempdir="/home/adam.viray/tmprast") # you need to make this dir. yourself


## sample a square around ignition point
# some rules of thumb by trial and error. 
# .001 = 111 meters; .1 = 11 km; 1=111,000 meters 
# .01 = 1.1 km 
# start with base of .05; 
# make adjustment factor based on log(RoS)*(prob)
radius=base.radius
w <- max(c(WEST, lon[[1]]-radius))
e <- min(c(EAST, lon[[1]]+radius))
s <- max(c(SOUTH, lat[[1]]-radius))
n <- min(c(NORTH, lat[[1]]+radius))
base.domain <- ext(w,e,s,n)

# topography-30m conus DEM
dem <- terra::rast(paste0(spat.dir,"dem_epsg4326.compressed.tiled.oviews.tif"))


## dem_crop will be template for all other data intputs. 
# aggregate it based on desired spatial resolution. 30m is very slow. 
ag.fact=round(SPATIAL.RES/30)
print("cropping initial dem")
dem_crop <- terra::crop(dem, base.domain)
dem_crop <- terra::aggregate(dem_crop, fact=ag.fact, fun="mean")

# get RoS and prob. raster values for initial domain. If it is dry, enlarge extent. 

# ignition probability
p.list <- list.files(pred.dir, pattern="R1_fire_growth_day", full.names = T)
#print(p.list)
p.start <- grep(start_date2, p.list)
p.list.sub <- p.list[p.start:length(p.list)]
prob <- rast()
for(i in 1:length(p.list.sub)) {
	p1 <- terra::rast(p.list.sub[i])
	p1p <- terra::crop(p1, dem_crop)
	prob <- c(prob, p1p)
}

p_mean=mean(prob)
p_mean_val=mean(values(p_mean), na.rm=T)
print(paste(round(p_mean_val, digits=1), "mean prob"))
if(p_mean_val < .01) { p_mean_val=.01 }

# rate of growth
r.list <- list.files(pred.dir, pattern="R1_fire_growth_acres_........_noveg", full.names = T)
#print(r.list)
r.start <- grep(start_date2, r.list)
r.list.sub <- r.list[r.start:length(r.list)]
RoS <- rast()
for(i in 1:length(r.list.sub)) {
	r1 <- terra::rast(r.list.sub[i])
	r1p <- terra::crop(r1, dem_crop)
	RoS <- c(RoS, r1p)
}

# testing RoS effect here. 
#RoS=RoS*10

# get info on predicted RoS. Use this to scale size of sim. domain
ros_mean=mean(RoS)
ros_mean_val=mean(values(ros_mean), na.rm=T)
print(paste(round(ros_mean_val, digits=1), "mean RoS"))


# modify extent here based on initial RoS and prob. conditions. 
radius.adj=max(c(log(ros_mean_val)*p_mean_val/100,0))
radius=base.radius+radius.adj
w <- max(c(WEST, lon[[1]]-radius))
e <- min(c(EAST, lon[[1]]+radius))
s <- max(c(SOUTH, lat[[1]]-radius))
n <- min(c(NORTH, lat[[1]]+radius))
print(c(w,e,s,n))
domain <- ext(w,e,s,n)
sim.ext=round(2*radius*111111, digits=0)/1000
print(paste("sim. extent ", sim.ext, "km"))

print("cropping adjusted dem")
dem_crop <- terra::crop(dem, domain)
dem_crop <- terra::aggregate(dem_crop, fact=ag.fact, fun="mean")

print(paste("sim domain has ", ncell(dem_crop), " cells"))

# crop prob. and growth rasters with final sim extent
# ignition probability
p.list <- list.files(pred.dir, pattern="R1_fire_growth_day", full.names = T)
#print(p.list)
p.start <- grep(start_date2, p.list)
p.list.sub <- p.list[p.start:length(p.list)]
prob <- rast()
for(i in 1:length(p.list.sub)) {
	p1 <- terra::rast(p.list.sub[i])
	p1p <- terra::crop(p1, dem_crop)
	prob <- c(prob, p1p)
}

# rate of growth
r.list <- list.files(pred.dir, pattern="R1_fire_growth_acres_........_noveg", full.names = T)
#print(r.list)
r.start <- grep(start_date2, r.list)
r.list.sub <- r.list[r.start:length(r.list)]
RoS <- rast()
for(i in 1:length(r.list.sub)) {
	r1 <- terra::rast(r.list.sub[i])
	r1p <- terra::crop(r1, dem_crop)
	RoS <- c(RoS, r1p)
}



###########################
# crop N.Rockies 30m veg. grids. Currently shrub,forest and nonveg.
 
# no vegetation
rcm=c(0,10,10, 100,0, 1)
rcm=matrix(rcm, 2)
noveg <- terra::rast(paste0(spat.dir,"rap_noveg_2022_30m.tif"))
noveg_crop <- terra::classify(terra::crop(noveg, dem_crop), rcm)

# vegetation
forest <- paste0(spat.dir,"rap_forest_2022_30m.tif")
shrub <- paste0(spat.dir,"rap_shrub_2022_30m.tif")
veg <- terra::rast(c(shrub,forest))
veg_crop <- terra::crop(veg, dem_crop)


# wind forecast files are in different projection. 
# Project our small wgs extent into that projection and use it to crop
# wind speed and projection files


# wind direction
wd.list <- list.files(wind.dir, pattern="WindDir_............_mst.tif$", full.names = T)
wind1=rast(wd.list[1])
wind_extent=terra::project(dem_crop, crs(wind1))
rm(wind1)

wd.start <- grep(start_date2, wd.list)
wd.list.sub <- wd.list[wd.start:length(wd.list)]
wdir <- terra::rast(wd.list.sub)
wdir_crop <- terra::crop(wdir, wind_extent)

# wind speed
ws.list <- list.files(wind.dir, pattern="WindSpd_............_mst.tif$", full.names=T)
ws.start <- grep(start_date2, ws.list)
ws.list.sub <- ws.list[ws.start:length(ws.list)]
wspd <- terra::rast(ws.list.sub)
wspd_crop <- terra::crop(wspd, wind_extent)



# make terrain indices 
aspect <- terra::terrain(dem_crop, "aspect", unit="degrees")
slope <- terra::terrain(dem_crop, "slope", unit="degrees")
hillshade <- terra::shade(slope*pi/180, aspect*pi/180)

# quick plots here, check spatial inputs
pdf(paste0(sim.dir, "plots.pdf"), height=6, width=6)
plot(dem_crop, main="elev")
plot(aspect, main="aspect")
plot(slope, main="slope")
plot(hillshade, main="shade")
plot(RoS, 1, main="ros")
plot(prob, 1, main="prob")
plot(wspd_crop, 1, main="wspd")
plot(wdir_crop, 1, main="wdir")
plot(veg_crop, main="veg")
dev.off()

## generate initial state raster
#state <- terra::rast(dem_crop)
state=dem_crop
terra::values(state) <- 0

# set center cell alight
state[round(dim(state)[[1]]/2)+1,round(dim(state)[[2]]/2)+1] <- 2


## get smallest number of hours
#k = min(c(dim(RoS)[3]*24,dim(prob)[3]*24,dim(wspd)[3],dim(wdir)[3]))
print(dim(RoS))
k = min(c(dim(RoS)[3]*24,dim(prob)[3]*24,dim(wspd)[3],dim(wdir)[3]))
l = ceiling(k/24)
print(paste("K = ", k))
print(paste("L = ", l))

## projecting rasters
prob_proj = terra::project(prob[[1:l]],dem_crop,method='bilinear')
RoS_proj = terra::project(RoS[[1:l]],dem_crop,method='bilinear')
wdir_proj = terra::project(wdir_crop[[1:k]],dem_crop,method='near')
wspd_proj = terra::project(wspd_crop[[1:k]],dem_crop,method='bilinear')
noveg_proj = terra::project(noveg_crop,dem_crop,method='near')
veg_proj = terra::project(veg_crop,dem_crop,method='bilinear')


## fill NAs
wdir_fill <- focal(wdir_proj, 13, "modal", na.policy="only")
wspd_fill <- focal(wspd_proj, 13, "modal", na.policy="only")


# quick plots here, check spatial inputs
pdf(paste0(sim.dir, "plots2.pdf"), height=6, width=6)
plot(dem_crop, main="elev")
plot(RoS_proj, 1, main="ros")
plot(prob_proj, 1, main="prob")
plot(wspd_fill, 1, main="wspd")
plot(wdir_fill, 1, main="wdir")
plot(veg_proj, main="veg")
dev.off()

## write projected rasters
print("writing data")
terra::writeRaster(prob_proj,paste0(sim.dir,"prob.tif"), overwrite=T)
terra::writeRaster(RoS_proj,paste0(sim.dir,"RoS.tif"), overwrite=T)
terra::writeRaster(dem_crop,paste0(sim.dir,"dem.tif"), overwrite=T)
terra::writeRaster(slope,paste0(sim.dir,"slope.tif"), overwrite=T)
terra::writeRaster(aspect,paste0(sim.dir,"aspect.tif"), overwrite=T)
terra::writeRaster(hillshade,paste0(sim.dir,"hillshade.tif"), overwrite=T)
terra::writeRaster(wdir_fill,paste0(sim.dir,"wdir.tif"), overwrite=T)
terra::writeRaster(wspd_fill,paste0(sim.dir,"wspd.tif"), overwrite=T)
terra::writeRaster(noveg_proj,paste0(sim.dir,"noveg.tif"), overwrite=T)
terra::writeRaster(veg_proj,paste0(sim.dir,"veg.tif"), overwrite=T)
terra::writeRaster(state,paste0(sim.dir,"state.tif"), overwrite=T)


## run python CA
# a should be:
# 1. path to data
# 2. path to output
# 3. path to parameters
# 4. start date time
# 5. spatial resolution
print("running CA")
a <- paste(sim.dir, sim.dir, "/home/adam.viray/Documents/CA/code/settings/default.txt", start_date, SPATIAL.RES, row_id, firename)
system(paste("python", "/mnt/DataDrive1/data/firespread_CA/scripts/CA.py", a))


## copy gif to upload dir
print("copy gif")
gif.dir <- "/var/www/html/firespread_jobs/firespread_modeling/finished_jobs/"
system(paste("cp", paste0(sim.dir,gifname), paste0(gif.dir,gifname)))
 

## clean up simulation directory
#file.remove(list.files(sim.dir, '\\d\\d_\\d\\d:\\d\\d:\\d\\d.png', full.names=T))


## finish script
#drv <- dbDriver("PostgreSQL")
#con <- dbConnect(drv, dbname="fire_activity",user="topofire",password="topofire")
# row_id is the record ID
# Grab the record info
# pass a download link back to database.
#url.link <- paste0("https://topofire.dbs.umt.edu/firespread_jobs/firespread_modeling/finished_jobs/", gifname)
#print("Updating job status")
#query<-paste0("Update firespread_jobs set status = 'Finished', download_link = '", url.link,  "' where id = ",row_id,";")
#rowset<-dbSendQuery(con,query)
#dbDisconnect(con)
