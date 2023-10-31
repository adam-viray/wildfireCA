require(terra)
require(doParallel)
require(foreach)

# IMPORTANT: if you're using terra and large rasters, careful of memory
# mak a temporary directory for R:terra like this:
terraOptions(tempdir="/home/adam.viray/tmprast") # you need to make this dir. yourself
