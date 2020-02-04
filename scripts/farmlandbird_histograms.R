# Create and save a histogram per farmland bird to see its distribution

# clear workspace
rm(list = ls())

# set working directory
setwd("~/Documents/GitHub/FarmlandBirdIndicators/plots/birddistribution")

# run script that loads data & extracts spatial coordinates of samples
source('~/Documents/GitHub/FarmlandBirdIndicators/scr/extract_coordinates.R', echo = FALSE)

# check distribution of counts
for(i in 1:ncol(bird)){
  png(filename = paste0(colnames(bird)[i],".png"))
  hist(bird[,i], main = colnames(bird)[i])
  dev.off()
}

