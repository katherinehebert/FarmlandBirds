# Script to make dotchart of annual population trends for bird species
# from the Ontario Breeding Bird Atlas (OBBA) between 2001 and 2005.
# (Supplementary figure)

# clear workspace
rm(list=ls())

# load packages
require(ggplot2)
require(tidyverse)
require(ggpubr)
require(ggsci)

# import OBBA trends (subsetted to final list of farmland birds for Ontario)
obba <- read.csv("data/obbatrends.csv")

# plot dotchart of specialists' trends 
S <- ggdotchart(obba[grep("S", obba$guild),], 
   x = "commonName", y = "trend",
   color = "guild",                              # Color by groups
   palette = c("#5CB85C", "#46B8DA", "#357EBD", "#9632B8"),
   sorting = "descending",                       # Sort value in descending order
   add = "segments",                             # Add segments from y = 0 to dots
   add.params = list(color = "lightgray", size = 1.5), # Change segment color and size
   group = "guild",                                # Order by groups
   dot.size = 5,  # Large dot size
   rotate = TRUE, 
   ggtheme = theme_pubr()) +                      # ggplot2 theme
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  xlab("Species") + ylab("OBBA trend (2001-2005)") + ylim(c(-2,2)) +
  theme(axis.text = element_text(size = 11), 
        axis.title = element_text(size = 12,face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        legend.position = "top",
        title = element_text(size = 12, face = "bold")) +
  ggtitle("a) Specialists")

# plot dotchart of generalists' trends 
G = ggdotchart(obba[grep("G", obba$guild),], x = "commonName", y = "trend",
                 color = "guild",
                 palette = c("#D43F3A", "#eea236"), 
                 # Color by groups
                 sorting = "descending", # Sort value in descending order
                 add = "segments", # Add segments from y = 0 to dots
                 add.params = list(color = "lightgray", size = 1.5), # Change segment color and size
                 group = "guild",  # Order by groups
                 dot.size = 5,  # Large dot size
                 rotate = TRUE, 
                 ggtheme = theme_pubr())+ # ggplot2 theme
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  ylab("OBBA trend (2001-2005)") + xlab(NULL) + ylim(c(-2,2)) +
  theme(axis.text = element_text(size = 11), 
        axis.title = element_text(size = 12,face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        legend.position = "top",
        title = element_text(size = 12, face = "bold")) +
  ggtitle("b) Generalists")

# generate 2-panelled plot and save figure
png("figures/obbatrend.png", width = 700, height = 600)
gridExtra::grid.arrange(S, G, ncol = 2)
dev.off()
