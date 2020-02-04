# Script to make correlation plot of agricultural land cover variables.
# (Supplementary figure)

# clear workspace
rm(list=ls())

# load packages
require(Hmisc)
require(corrplot)

# load agricultural land cover variables 
# (previously standardized using decostand)
land <- read.csv("data/landcover.csv")

# calculate correlation matrix (with p values)
res = rcorr(as.matrix(land), type = "pearson")

# visualize correlation matrix 
png("figures/correlationmatrix.png", width = 500, height = 500)
corrplot(res$r, type = "upper", order = "hclust", 
         tl.col = "black", sig.level = 0.05, tl.srt = 45)
dev.off()
