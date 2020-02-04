# Script to calculate Moran's I to detect spatial autocorrelation
# (Supplementary figure)

# clear workspace
rm(list=ls())

# load packages
require(spdep)

# read bird data
# NOTE: Available upon request.

# read spatial coordinates of agricultural subdivisions (centroids)
coor <- as.matrix(read.csv("data/centroids.csv"))

# create a k=4 nearest neighbor set (prepare for Moran test function)
us.nb4 <- knearneigh(coor, k=4)
us.nb4 <- knn2nb(us.nb4)
us.nb4 <- make.sym.nb(us.nb4)
us.wt4 <- nb2listw(us.nb4, style="W")

# Moran's test: test bird count distributions for spatial autocorrelation
spautocorr = apply(bird, 2, moran.test, listw = us.wt4)

# format as a table for easier browsing
spautocorr_res = matrix(ncol = 2, nrow = ncol(bird), 
                        dimnames = list(colnames(bird), c("Moran's I", "p")))
for(i in 1:nrow(spautocorr_res)){
  spautocorr_res[i, 1] = signif(spautocorr[[i]]$statistic, digits = 3)
  spautocorr_res[i, 2] = signif(spautocorr[[i]]$p.value, digits = 1)
}
# save table
write.csv(spautocorr_res, "outputs/test_moran.csv")

# get names of species whose distributions are spatially autocorrelated (p<.05)
names(which(spautocorr_res[,2] < .05))
