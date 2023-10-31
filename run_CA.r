require(doParallel)
require(foreach)

in.dir = "/home/adam.viray/Documents/CA/data/riceridge/"
out.dir = "/home/adam.viray/Documents/CA/output/"

# arguments should be:
# 1. path to VIIRS folder. should look like /mnt/mordor3/data/VIIRS/riceridge/model_inputs
# 2. path to output file. should look like /home/adam.viray/Documents/CA/output/1.csv
# 3. parameter 'weight'. default 0.5, could be between 0 and 10
# 4. parameter 'weight_d'. default 0.8, could be between 0 and 10
# 5. parameter 'weight_phi'. default 0.8, could be between 0 and 1
# 6. parameter 'spread_threshold'. default 0.95, could be between 0 and 1
# 7. parameter 'burn_speed'. default 0.06, could be between 0 and 10
# 8. parameter 'spot_sigma'.
# 9. parameter 'brand_sigma'.
r <- c(0.325,1.6,0.7,0.925,0.07,0.3,0.3)
a <- 0.025
b <- 0.2
c <- 0.5
d <- 0.00625
e <- 0.0025
f <- 0.1
g <- 0.1
pars_mtx = rbind(r,r+c(a,0,0,0,0,0,0),r-c(a,0,0,0,0,0,0),r+c(0,b,0,0,0,0,0),r-c(0,b,0,0,0,0,0)
	,r+c(0,0,c,0,0,0,0),r-c(0,0,c,0,0,0,0),r+c(0,0,0,d,0,0,0),r-c(0,0,0,d,0,0,0)
	,r+c(0,0,0,0,e,0,0),r-c(0,0,0,0,e,0,0),r+c(0,0,0,0,0,f,0),r-c(0,0,0,0,0,f,0)
	,r+c(0,0,0,0,0,0,g),r-c(0,0,0,0,0,0,g))

cl = makeCluster(11)
registerDoParallel(cl)

foreach(i=1:11) %dopar% {
	print(i)
	p = pars_mtx[i,]
	fname = paste0(out.dir,i,'.csv')
	system(paste('python /home/adam.viray/Documents/CA/code/CA.py', in.dir, fname, p[1], p[2], p[3], p[4], p[5], p[6], p[7]))
}

stopCluster(cl)


