# univariate regression tree per bird species (compute & plot)

# clear workspace
rm(list=ls())

# set working directory
setwd("~/Documents/GitHub/FarmlandBirdIndicators")

# load packages
library(rpart)
library(rpart.plot)
library(adespatial)
library(spdep)
library(ggplot2)
library(magrittr)
library(tidyverse)

# run script that loads data & extracts spatial coordinates of samples
source('~/Documents/GitHub/FarmlandBirdIndicators/scr/extract_coordinates.R', echo = FALSE)

# set up coordinate matrix
coor <- matrix(0, nrow(spdat), 2, dimnames = list(NULL, c("xcoord", "ycoord")))
coor[,1] <- spdat$xcoord/100000
coor[,2] <- spdat$ycoord/100000

# Create a k=4 nearest neighbor set
us.nb4<-knearneigh(coor, k=4)
us.nb4<-knn2nb(us.nb4)
us.nb4<-make.sym.nb(us.nb4)
us.wt4<-nb2listw(us.nb4, style="W")

# test bird count distributions for spatial autocorrelations
spautocorr = apply(bird, 2, moran.test, listw = us.wt4)
spautocorr_res = matrix(ncol = 2, nrow = ncol(bird), dimnames = list(colnames(bird), c("Moran's I", "p")))
for(i in 1:nrow(spautocorr_res)){
  spautocorr_res[i, 1] = signif(spautocorr[[i]]$statistic, digits = 3)
  spautocorr_res[i, 2] = signif(spautocorr[[i]]$p.value, digits = 1)
}
write.csv(spautocorr_res, "~/Documents/GitHub/FarmlandBirdIndicators/MoransI.csv")

# bind coordinates to land data
land = cbind(land, coor)

# subset land data to required variables
land <- subset(land, select = -c(SDA_HA, CSDNAME, geometry))
colnames(land) = gsub("xcoord", "LONG", colnames(land))
colnames(land) = gsub("ycoord", "LAT", colnames(land))
land = subset(land, select = -c(LONG.1, LAT.1, Bpas))

# species that that showed positive spatial autocorr: AMKE, EUST, HOSP, SAVS
spasp = "RTHA"

## compute, plot, and save regression trees ----

# create list of trees
tree <- list()
# set directory to save .png of trees
setwd("~/Documents/GitHub/FarmlandBirdIndicators/plots/urt_perbird")

# print(tree[[1]])
# printcp(tree[[1]]) # display the results 
# plotcp(tree[[1]]) # visualize cross-validation results 
# summary(tree[[1]]) # detailed summary of splits

# run non-spatial RTs
#bird.nospa = c("HOSP","SAVS","EUST", "AMKE")

tree = list()
# for each bird species without spatial info
for(i in 1:ncol(bird)){
  # regression tree on bird abundances
  tree[[i]] <- rpart(bird[,i] ~ ., 
                     data = subset(land, select = -c(LAT, LONG)), 
                     method = "poisson",
                     maxdepth = 2,
                     model = TRUE) 
  # Visualize the decision tree with rpart.plot
  png(paste0(i,"_urt.png"), width = 600, height = 500)
  rpart.plot(tree[[i]], box.palette = "#92C5DE",
             shadow.col = "gray", 
             clip.right.labs = FALSE,
             nn = TRUE)
  dev.off()
}
# for each bird species
for(i in 1){
  # regression tree on bird abundances
  tree[[47]] <- rpart(bird[,spasp[i]] ~ ., data = land, 
                     method = "poisson",
                     maxdepth = 2,
                     model = TRUE) 
  # Visualize the decision tree with rpart.plot
  png(paste0(spasp[i],"_urt.png"), width = 600, height = 500)
  rpart.plot(tree[[47]], box.palette = "#92C5DE",
             shadow.col = "gray", 
             clip.right.labs = FALSE,
             nn = TRUE)
  dev.off()
}

r2 = as.data.frame(matrix(NA, nrow = length(tree), ncol = 2))
colnames(r2) = c("SpeciesCode", "R2")
r2$SpeciesCode= colnames(bird)
for(i in 1:length(tree)){
  r2$R2[i] = tree[[i]]$cptable[1,1] # extract complexity paramater at node 0, which is R2 of whole tree
}
write.csv(r2, "~/Documents/GitHub/FarmlandBirdIndicators/urt_R2.csv")

# extract important splits
labels(tree[[1]])
temp <- unlist(strsplit(labels(tree[[1]]), split = c("<")))
temp <- unlist(strsplit(temp, split = ">="))
chrtemp <- parse_character(temp)
xtemp <- parse_number(temp)
chrtemp <- chrtemp[which(is.na(xtemp))]
vars <- unique(chrtemp)
tree[[1]]$variable.importance
summary(tree[[1]])[1]


# prune the tree by selecting the complexity parameter associated 
ptree <- list()
for(i in 1:length(tree)){
    # with minimum error,
  ptree[[i]] <- prune(tree[[i]], 
                      cp = tree[[i]]$cptable[which.min(tree[[i]]$cptable[,"xerror"]),"CP"])
  
  if (nrow(ptree[[i]]$cptable) > 1) {
  # plot tree 
    png(paste0(colnames(bird.nospa)[i],"_prunedurt.png"), width = 600, height = 500)
    rpart.plot(ptree[[i]], box.palette = "#92C5DE", 
             shadow.col = "gray", 
             nn=TRUE)
    dev.off()
  }
}


## variable importance & inclusion in the trees ----

# pull out variables used in tree construction
impvar = list()
for(i in 1:length(tree)){
  impvar[[i]] = levels(droplevels(tree[[i]]$frame$var, "<leaf>"))  
}
# create table to store resilts
df = data.frame(matrix(NA, nrow = ncol(bird), ncol = (ncol(land))),
                row.names = r2$SpeciesCode)
colnames(df) = colnames(land)
for(i in 1:length(impvar)){
  df[i, impvar[[i]]] = 1
}
# write file
write.csv(df, "~/Documents/GitHub/FarmlandBirdIndicators/urt_importantvariables.csv")
df <- read.csv("~/Documents/GitHub/FarmlandBirdIndicators/urt_importantvariables.csv", row.names = 1)

# pull out variable importance
impvar.values = list()
for(i in 1:length(tree)){
  impvar.values[[i]] = tree[[i]]$variable.importance
}
# create table to store results
df2 = data.frame(matrix(NA, nrow = ncol(bird), ncol = (ncol(land))),
                row.names = r2$SpeciesCode)
colnames(df2) = colnames(land)
for(i in 1:length(impvar.values)){
  vars = names(impvar.values[[i]])
  for(n in 1:length(impvar.values[[i]])){
    df2[i, vars[n]] = unname(impvar.values[[i]][n])
  }
}
write.csv(df2, "~/Documents/GitHub/FarmlandBirdIndicators/urt_importantvariablesvalues.csv")
df2 <- read.csv("~/Documents/GitHub/FarmlandBirdIndicators/urt_importantvariablesvalues.csv", row.names = 1)

# rank "importance values" per row (i.e. per species)
#df2 = subset(df2, select = -c(LAT.1, LONG.1))
df2.ranks = apply(df2, 1, rank, na.last = "keep")

# prepare df for plotting
df = as.matrix(df)
df[which(is.na(df))] = 0
df = apply(df, 1:2, as.numeric)
#df = subset(df, select = -c(LAT.1, LONG.1))

# order according to functional guild
setwd("~/Documents/GitHub/FarmlandBirdIndicators")
# read in final bird classification list
# guilds = read.csv("~/Documents/GitHub/FarmlandBirdIndicators/data/finallist.csv",
#                   colClasses = "character")
# guilds$guild[which(guilds$species %in% "Gray Partridge")] <- "FG"
# guilds$guild[which(guilds$species %in% "Common Nighthawk")] <- "FG"
# guilds$guild[which(guilds$species %in% "Northern Mockingbird")] <- "FG"
# guilds$guild[which(guilds$species %in% "Red-headed Woodpecker")] <- "FG"
# 
# # correct species names to match NABBS spelling/capitalization
# guilds$species = gsub("Henslow\xd5s Sparrow", "Henslow's Sparrow", guilds$species)
# guilds$species = gsub("Red-tailed Hawk ", "Red-tailed Hawk", guilds$species)
# guilds$species = gsub(" Eastern Wood-pewee", "Eastern Wood-Pewee", guilds$species)

# read in expert guild classification
cats = read.csv("data/finalclassification.csv", colClasses = "character")
# correct abbreviation(s)
rownames(df)[which(rownames(df) %in% "CAGO")] <- "CANG"
rownames(df2)[which(rownames(df2) %in% "CAGO")] <- "CANG"
rownames(df)[which(rownames(df) %in% "GRPA")] <- "GRAP"
cats$abb[which(cats$abb %in% "GRPA")] <- "GRAP"
# bind bird guild categories to tables
cats.birds = cats[which(cats$abb %in% rownames(df)),]
# sort into functional guilds
cats.birds = cats.birds[order(cats.birds$final),]

# plot variable inclusion matrix
df.l = reshape2::melt(df)
df.l$Var1 = factor(df.l$Var1, levels = cats.birds$abb)
df.l$value = as.character(df.l$value)
# attach final classification here
df.l$finalclass = NA
for(i in 1:nrow(df.l)){
  df.l$finalclass[i] = cats.birds$final[which(cats.birds$abb %in% df.l$Var1[i])]
}
df.l$finalclass[which(df.l$value == 0)] <- NA
df.l$finalclass = as.factor(df.l$finalclass)
df.l$value = as.integer(df.l$value)


(gg.incl = ggplot(filter(df.l, !is.na(df.l$Var2)), aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = finalclass)) +
  labs(x = "predictor variables", y = "species") +
  ggtitle("a)") +
  ggsci::scale_fill_locuszoom(name="Final\nclassification", 
                              na.value = "white",
                              labels = c(levels(df.l$finalclass), "not selected")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, size = 10, vjust=0.9),
          axis.title=element_text(size=14,face="bold"),
          legend.title = element_text(size=14,face="bold"),
          legend.text = element_text(size=11),
          plot.title = element_text(size = 18, face = "bold"))
  )

# plot variable importance matrix
df2.ranks.l = reshape2::melt(df2.ranks)
df2.ranks.l$Var2 = factor(df2.ranks.l$Var2, levels = cats.birds$abb)
df2.ranks.l$value = as.integer(round(df2.ranks.l$value, digits = 1))
(gg.imp = ggplot(filter(df2.ranks.l, !is.na(df2.ranks.l$Var2)), aes(x = Var1, y = Var2)) +
                  geom_raster(aes(fill = value)) +
                  labs(x = "predictor variables", y = "species") +
                  ggtitle("b)")+ 
  scale_fill_viridis_c("Variable\nimportance", option = "D", na.value = "white") +
                  # scale_fill_gradient2("rank", 
                  #                      low = "#deebf7", mid = "#6baed6", high = "#0c2c84",
                  #                      na.value = "white", midpoint = 8) +
                  theme(axis.text.x = element_text(angle = 90, size = 10, vjust=0.9),
                        axis.title=element_text(size=14,face="bold"),
                        legend.title = element_text(size=14,face="bold"),
                        legend.text = element_text(size=11),
                        plot.title = element_text(size = 18, face = "bold"))
  )
# save plots side-by-side
gridExtra::grid.arrange(gg.incl, gg.imp, ncol=2)
ggsave("RT_variables.png", gridExtra::grid.arrange(gg.incl, gg.imp, ncol=2))
