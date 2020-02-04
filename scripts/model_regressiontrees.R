# Script to compute and plot univariate regression tree for each bird species
# (Supplementary figures)

# clear workspace
rm(list=ls())

# load packages
library(rpart)
library(rpart.plot)

# read bird and agricultural land cover data
# NOTE: Available upon request.

# run regression trees (and save plots and models)
tree = list()
for(i in 1:ncol(bird)){
  # regression tree: bird abundance as a function of all landcover variables
  tree[[i]] <- rpart(bird[,i] ~ ., data = land, 
                     method = "poisson",
                     maxdepth = 2, model = TRUE) 
  # visualize the tree (and save)
  png(paste0("figures/regressiontrees/",colnames(bird)[i],".png"), 
      width = 600, height = 500)
  rpart.plot(tree[[i]], box.palette = "#92C5DE",
             shadow.col = "gray", clip.right.labs = FALSE, nn = TRUE)
  dev.off()
}

#### rerun regression tree for species -----------------------------------------
#### with spatially autocorrelated abundances (i.e., "RTHA")

# read agricultural subdivisions' centroid coordinates
coor <- read.csv("data/centroids.csv")
# bind to land cover dataset
land <- cbind(land, coor)

# overwrite tree for RTHA with one including spatial coordinates
tree[[which(colnames(bird) %in% "RTHA")]] <- rpart(bird[,"RTHA"] ~ ., data = land, 
                                                    method = "poisson",
                                                    maxdepth = 2, model = TRUE) 
# visualize the tree (and save)
png("figures/regressiontrees/RTHA.png", 
    width = 600, height = 500)
rpart.plot(tree[[which(colnames(bird) %in% "RTHA")]], box.palette = "#92C5DE",
           shadow.col = "gray", clip.right.labs = FALSE, nn = TRUE)
dev.off()

# save list of models
saveRDS(tree, "outputs/regressiontrees.RDS")

