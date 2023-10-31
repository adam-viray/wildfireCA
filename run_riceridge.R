require(terra)

# parameter
spat.res = 30

for(spat.res in c(60,120,240)) {
agg.fact = spat.res/30
disagg.fact = 240/spat.res

# IMPORTANT: if you're using terra and large rasters, careful of memory
# mak a temporary directory for R:terra like this:
terraOptions(tempdir="/home/adam.viray/tmprast") # you need to make this dir. yourself

## read in data for extent
print("reading data")
# data locations
data.dir <- "/home/adam.viray/Documents/CA/data/riceridge/"

# 30m rate of spread from 240m
if (disagg.fact > 1) {
    RoS <- terra::disagg(terra::rast(paste0(data.dir,"MT4726811348520170724_daily_lme_log10growth_weatheronly_gwqr_preds_240m.nc")), disagg.fact, method='bilinear')
} else {
    RoS <- terra::rast(paste0(data.dir,"MT4726811348520170724_daily_lme_log10growth_weatheronly_gwqr_preds_240m.nc"))
}

# 30m dem
if (agg.fact > 1) {
    dem <- terra::aggregate(terra::rast(paste0(data.dir,"dem_30.tif")), agg.fact, cores=8)
} else {
    dem <- terra::rast(paste0(data.dir,"dem_30.tif"))
}

# extent of RoS (which is smallest) in EPGS:4326
extent <- terra::project(terra::ext(RoS), terra::crs(RoS), terra::crs(dem))

#crop dem
dem <- terra::crop(dem, extent)

# change RoS projection
RoS <- terra::project(RoS, dem, threads=T)

# topography
aspect <- terra::terrain(dem, "aspect", unit="degrees")
slope <- terra::terrain(dem, "slope", unit="degrees")
hillshade <- terra::shade(slope*pi/180, aspect*pi/180)

# VIIRS burnday
burnday <- terra::project(terra::crop(terra::rast(paste0(data.dir,"riceridge_burnday_250m_wgs.tif")), extent), dem, method="near")

# RAP vegetation 2016
veg <- terra::rast(paste0(data.dir,"vegetation_riceridge_2016_250m.tif"))

# no vegetation
rcm=c(0,10,10, 100,0, 1)
rcm=matrix(rcm, 2)
noveg <- terra::project(terra::classify(terra::crop(veg[[2]], extent), rcm), dem, method='bilinear')

# vegetation
veg <- veg[[5:6]]
veg <- terra::project(terra::crop(veg, extent), dem, method="bilinear")

# wind
wdir <- terra::project(terra::crop(terra::rast(paste0(data.dir,"rtma_hourly_wdir_wgs_250m.tif")), extent), dem, method="near", threads=T)
wspd <- terra::project(terra::crop(terra::rast(paste0(data.dir,"rtma_hourly_wspd_wgs_250m.tif")), extent), dem, method="bilinear", threads=T)

# ignition probability
#prob <- terra::crop(terra::project(terra::rast(paste0(data.dir,"fire_prob_201706-201709_250m_v2.nc"))[[-1:-53]], dem_base, method="bilinear"), extent)

# rate of growth
#RoS <- terra::crop(terra::project(terra::rast(paste0(data.dir,"fire_growth_meters_201706-201709_250m.nc"))[[-1:-53]], dem_base, method="bilinear"), extent)

## get smallest number of hours
k = min(c(dim(RoS)[3]*24,dim(prob)[3]*24,dim(wspd)[3],dim(wdir)[3]))
l = ceiling(k/24)

##create a temp. directory for running each sim
sim.dir <- paste0("/home/adam.viray/Documents/CA/data/riceridge/",spat.res,"m/")
if(!dir.exists(sim.dir)) { dir.create(sim.dir) }

## generate initial state raster
state <- terra::init(dem,0)
state[terra::rowFromY(state,47.2699),terra::colFromX(state,-113.4814)] <- 2

## write cropped rasters
print("writing data")
#terra::writeRaster(terra::crop(prob[[1:l]],dem),paste0(sim.dir,"prob.tif"), overwrite=T)
terra::writeRaster(RoS[[1:l]],paste0(sim.dir,"RoS.tif"), overwrite=T)
terra::writeRaster(dem,paste0(sim.dir,"dem.tif"), overwrite=T)
terra::writeRaster(slope,paste0(sim.dir,"slope.tif"), overwrite=T)
terra::writeRaster(aspect,paste0(sim.dir,"aspect.tif"), overwrite=T)
terra::writeRaster(hillshade,paste0(sim.dir,"hillshade.tif"), overwrite=T)
terra::writeRaster(wdir[[1:k]],paste0(sim.dir,"wdir.tif"), overwrite=T)
terra::writeRaster(wspd[[1:k]],paste0(sim.dir,"wspd.tif"), overwrite=T)
terra::writeRaster(veg,paste0(sim.dir,"veg.tif"), overwrite=T)
terra::writeRaster(noveg,paste0(sim.dir,"noveg.tif"), overwrite=T)
terra::writeRaster(burnday,paste0(sim.dir,"burnday.tif"), overwrite=T)
terra::writeRaster(state,paste0(sim.dir,"state.tif"), overwrite=T)
}


#dem_250 <- terra::rast('/home/adam.viray/Documents/CA/data/riceridge/elevation_250m_wgs.tif')
#aspect_250 <- terra::terrain(dem_250, "aspect", unit="degrees")
#slope_250 <- terra::terrain(dem_250, "slope", unit="degrees")
#hillshade_250 <- terra::shade(slope_250*pi/180, aspect_250*pi/180)
#terra::writeRaster(hillshade_250, '/home/adam.viray/Documents/CA/data/riceridge/hillshade_250m_wgs.tif',overwrite=T)
#terra::writeRaster(slope_250, '/home/adam.viray/Documents/CA/data/riceridge/slope_250m_wgs.tif',overwrite=T)
#terra::writeRaster(aspect_250, '/home/adam.viray/Documents/CA/data/riceridge/riceridge_aspect_wgs_250m.tif',overwrite=T)
