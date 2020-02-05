# Script to create grid plots of variable selection and variable importance in
# the regression tree models.

# load packages
require(tidyverse)
require(reshape2)
require(ggsci)

# load formatted regression tree results
var <- read.csv("outputs/rt_variables.csv", row.names = 1)
imp <- read.csv("outputs/rt_importance.csv", row.names = 1)
# remove spatial coordinate columns
var <- subset(var, select = -c(xcoord, ycoord))
imp <- subset(imp, select = -c(xcoord, ycoord))

# rank importance values per row (i.e. per species)
imp <- as.data.frame(t(apply(imp, 1, rank, na.last = "keep")))

# order tables according to functional guild
guild <- read.csv("data/birdguilds.csv", colClasses = "character")
# bind guild categories to tables
var$guild = guild$final[which(guild$abb %in% gsub("CAGO", "CANG", rownames(var)))]
imp$guild = guild$final[which(guild$abb %in% gsub("CAGO", "CANG", rownames(imp)))]

# order according to groups
var = var[order(var$guild, decreasing = TRUE),]
imp = imp[order(imp$guild, decreasing = TRUE),]
# make species abbreviation its own column (factor to retain order when plotted)
var$sp = factor(rownames(var), levels = rownames(var))
imp$sp = factor(rownames(imp), levels = rownames(imp))

# plot variable inclusion matrix
df = melt(var)
df$guild = factor(df$guild, levels = c("RCS", "PSS", "FSS", "FEG", "FG"))
df$value = as.character(df$value)
# create plot
(gg.var = ggplot(filter(df, value == "1"), aes(x = variable, y = sp)) +
    geom_raster(aes(fill = guild)) +
    labs(x = "predictor variables", y = "species") +
    ggtitle("a)") +
    scale_fill_locuszoom(name="Final\nclassification", 
                          na.value = "white",
                          labels = c(levels(df.l$guild), "not selected")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, size = 10, vjust=0.9),
          axis.title=element_text(size=14,face="bold"),
          legend.title = element_text(size=14,face="bold"),
          legend.text = element_text(size=11),
          plot.title = element_text(size = 18, face = "bold"))
)

# plot variable importance matrix
df2 = melt(imp)
df2$guild = factor(df2$guild, levels = c("RCS", "PSS", "FSS", "FEG", "FG"))
df2$value = as.integer(round(df2$value, digits = 1))
(gg.imp = ggplot(filter(df2, !is.na(df2$value)), aes(x = variable, y = sp)) +
    geom_raster(aes(fill = value)) +
    labs(x = "predictor variables", y = "species") +
    ggtitle("b)")+ 
    scale_fill_viridis_c("Variable\nimportance", option = "D", na.value = "white") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, size = 10, vjust=0.9),
          axis.title=element_text(size=14,face="bold"),
          legend.title = element_text(size=14,face="bold"),
          legend.text = element_text(size=11),
          plot.title = element_text(size = 18, face = "bold"))
)
# save plots side-by-side
gridExtra::grid.arrange(gg.var, gg.imp, ncol=2)
ggsave("figures/rtvariables.png", gridExtra::grid.arrange(gg.var, gg.imp, ncol=2))