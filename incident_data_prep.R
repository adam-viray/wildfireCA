require(terra)

# IMPORTANT: if you're using terra and large rasters, careful of memory
# mak a temporary directory for R:terra like this:
terraOptions(tempdir="/home/adam.viray/tmprast") # you need to make this dir. yourself

# data locations
wind.dir <- "/mnt/DataDrive1/data/rtma_wind/"
growth.dir <- "/mnt/DataDrive1/data/VIIRS/ncdf_preds/"
obs.dir <- "/mnt/DataDrive1/data/VIIRS/ncdf_v2/"
fname.dir <- "/mnt/DataDrive1/data/VIIRS/test_fires/"
sim.dir <- "/home/adam.viray/Documents/CA/data/test_fires/"

# get a list of incident codes from a folder. could be RoS, GIFs, wind, whatever
inc.list <- lapply(list.files(fname.dir), {function(x) sub("\\-.*", "", x)})

# read in full-extent data
dem <- terra::rast("/mnt/DataDrive1/data/elevation/dem_west_30m_crop_mask.tif") # EPSG 5070 Albers
veg <- terra::rast("/mnt/DataDrive1/data/RAP_cover/vegetation-cover-v3-2016.tif") # EPSG 4326 

# iterate over incidents
for (inc in inc.list){
    print(inc)
    
    # create a directory for running each incident
    inc.dir <- paste0(sim.dir,inc,"/")
    if(!dir.exists(inc.dir)) { dir.create(inc.dir) }
    
    for (spat.res in c(120,240)){
        print(spat.res)
        agg.fact = spat.res/30
        
        # read in per-incident data
        inc.wdir <- terra::rast(paste0(wind.dir,inc,"_winddir.nc")) # ESRI 102039 - same as EPSG 5070?
        inc.wspd <- terra::rast(paste0(wind.dir,inc,"_windspeed.nc")) # ESRI 102039 - same as EPSG 5070?
        inc.growth <- terra::rast(list.files(growth.dir,paste0(inc,".*"),full.names=T)) # ESRI 102039 - same as EPSG 5070?
        inc.burnday <- terra::rast(list.files(obs.dir,paste0(inc,"_static.*"),full.names=T),subds="krige_jday") # ESRI 102039 - same as EPSG 5070?
        
        # discard data before day 0
        jday <- where.min(inc.burnday)[3]
        start_date <- as.Date(jday, origin=as.Date('2017-01-01'))
        inc.wdir <- inc.wdir[[time(inc.wdir) > start_date-1]]
        inc.wspd <- inc.wspd[[time(inc.wspd) > start_date-1]]
        inc.growth <- inc.growth[[time(inc.growth) > start_date-1]]
        
        # get incident extent
        inc.ext <- terra::ext(inc.burnday)
        
        # crop full-extent data to incident
        inc.dem <- terra::crop(dem,inc.ext)
        inc.veg <- terra::project(terra::crop(veg, terra::project(inc.ext,terra::crs(inc.burnday),terra::crs(veg))), inc.dem)
        
        # aggregate dem and project all others
        inc.dem <- terra::aggregate(inc.dem, agg.fact, cores=8)
        inc.wdir <- terra::project(inc.wdir, inc.dem, method="near", threads=T)
        inc.wspd <- terra::project(inc.wspd, inc.dem, method="bilinear", threads=T)
        inc.growth <- terra::project(inc.growth, inc.dem, method="bilinear", threads=T)
        inc.burnday <- terra::project(inc.burnday, inc.dem, method="bilinear", threads=T)
        inc.veg <- terra::project(inc.veg, inc.dem, method="bilinear", threads=T)
        
        # produce topo derivatives
        inc.aspect <- terra::terrain(inc.dem, "aspect", unit="degrees")
        inc.slope <- terra::terrain(inc.dem, "slope", unit="degrees")
        inc.hillshade <- terra::shade(inc.slope*pi/180, inc.aspect*pi/180)
        
        # produce veg derivatives
        rcm <- c(0,10,10, 100,0, 1)
        rcm <- matrix(rcm, 2)
        inc.noveg <- terra::classify(inc.veg[[2]], rcm)
        inc.veg <- inc.veg[[5:6]]
        
        # generate initial state raster
        inc.state <- terra::init(inc.dem,0)
        inc.state[inc.burnday == where.min(inc.burnday)[3]] <- 2
        
        # get time span
        k = min(c(dim(inc.growth)[3]*24,dim(inc.wspd)[3],dim(inc.wdir)[3]))
        l = ceiling(k/24)
        
        # create a directory for running each resolution
        res.dir <- paste0(inc.dir,spat.res,"m/")
        if(!dir.exists(res.dir)) { dir.create(res.dir) }
        
        # write data
        terra::writeRaster(inc.growth[[1:l]],paste0(res.dir,"RoS.tif"), overwrite=T)
        terra::writeRaster(inc.dem,paste0(res.dir,"dem.tif"), overwrite=T)
        terra::writeRaster(inc.slope,paste0(res.dir,"slope.tif"), overwrite=T)
        terra::writeRaster(inc.aspect,paste0(res.dir,"aspect.tif"), overwrite=T)
        terra::writeRaster(inc.hillshade,paste0(res.dir,"hillshade.tif"), overwrite=T)
        terra::writeRaster(inc.wdir[[1:k]],paste0(res.dir,"wdir.tif"), overwrite=T)
        terra::writeRaster(inc.wspd[[1:k]],paste0(res.dir,"wspd.tif"), overwrite=T)
        terra::writeRaster(inc.veg,paste0(res.dir,"veg.tif"), overwrite=T)
        terra::writeRaster(inc.noveg,paste0(res.dir,"noveg.tif"), overwrite=T)
        terra::writeRaster(inc.burnday,paste0(res.dir,"burnday.tif"), overwrite=T)
        terra::writeRaster(inc.state,paste0(res.dir,"state.tif"), overwrite=T)
    }
}
