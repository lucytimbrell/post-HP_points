# Timbrell, L. et al. (2022). The post-Howiesons Poort points from Border Cave, KwaZulu Natal, South Africa. QSR. 
# This script performs the geometric morphometric analysis of Border Cave

# Data available to download at: https://github.com/lucytimbrell/post-HP_points

# To reproduce the analysis, download the data as a zip file, extract it to a folder, and set the
# working directory below to the folder where you extracted the data.

# Ensure to use the .rds files in GitHub to ensure replicability (as landmarking can be subjective and affect results)

# Clear R environment and set working directory
rm(list=ls())

setwd("...")

## Install and load packages
if(!require("Momocs")) install.packages('Momocs', repos='http://cran.us.r-project.org')  
if(!require("tidyverse")) install.packages('tidyverse', repos='http://cran.us.r-project.org')
if(!require("rio")) install.packages('rio', repos='http://cran.us.r-project.org')
if(!require("ggplot2")) install.packages('ggplot2', repos='http://cran.us.r-project.org')
if(!require("ggpubr")) install.packages('ggpubr', repos='http://cran.us.r-project.org')
if(!require("patchwork")) install.packages('patchwork', repos='http://cran.us.r-project.org')
if(!require("MASS")) install.packages('MASS', repos='http://cran.us.r-project.org')

library(Momocs)
library(rio)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(patchwork)
library(MASS)

remotes::install_github("yesdavid/outlineR", dependencies = TRUE, force = TRUE)
library(outlineR)

#### Data import and preparation ####

tpslines <- import_tps("BC_outlinedata.tps")

database <- read.csv("BC_data.csv") # Database

database$Tip.damage <- as.factor(database$Tip.damage) 
database$Site <- as.factor(database$Site) 
database$Layer <- as.factor(database$Layer)
database$Member <- as.factor(database$Member)
database$Relative.stratigraphy <- as.factor(database$Relative.stratigraphy)
database$Raw.materials <- as.factor(database$Raw.materials)
database$Raw.materials.2 <- as.factor(database$Raw.materials.2)
database$Knapping.modality <- as.factor(database$Knapping.modality)
database$New.strat <- as.factor(database$New.strat)
database$Proximal.dorsal.preparation <- as.factor(database$Proximal.dorsal.preparation)
database$Ridges.and.Nodes <- as.factor(database$Ridges.and.Nodes)
database$Dorsal.scar.pattern <- as.factor(database$Dorsal.scar.pattern)
database$Angle <- as.factor(database$Angle)
database$Orientation.of.the.tip <- as.factor(database$Orientation.of.the.tip)

saveRDS(tpslines, file =  "BC_tpslines.rds")
saveRDS(database, file =  "BC_database.rds")

#### Import data and create outline object ####
tpsdata <- import("BC_tpslines.rds") 
database <- import("BC_database.rds")   

shape <- Out(tpsdata$coo, fac = database)
names(shape) <- database$Number

panel(shape, main = "Outline data") # Visualization of points in their original orientation

#### Outline normalization, landmarking and artefact-specific processing ####

shapenorm <- shape %>% 
  coo_centre() %>% 
  coo_scale() %>% 
  coo_slidedirection("right") %>% 
  coo_untiltx() %>% 
  coo_close() %>% 
  coo_aligncalliper() 

panel(shapenorm, main = "Outline data") # Visualization of points in their original orientation

shapenorm1 <- def_ldk(shapenorm, 1) # add new landmarks to tip and base for additional orientation
shapenorm2 <- coo_untiltx(shapenorm1, ldk = 1) # option 2 - reduce tilt according to tip landmark - WE USE THIS
shapenorm3 <- coo_rotate(shapenorm2, theta = (pi/2)*3) # rotate so points are facing up
shapenorm4 <- coo_smooth(shapenorm3, 10)

# Artefact-specific processing
coo_plot(shapenorm3[14], col = "grey", centroid = TRUE, main = "Artefact 2297") # this outline needs further smoothing
shapenorm4[14] <- coo_smooth(shapenorm4[14], 1000)  
coo_plot(shapenorm4[14], col = "grey", centroid = TRUE, main = "Artefact 2297")

coo_plot(shapenorm3[51], col = "grey", centroid = TRUE, main = "Artefact P121") # this outline needs futher smoothing
shapenorm4[51] <- coo_smooth(shapenorm4[51], 500)  
coo_plot(shapenorm4[51], col = "grey", centroid = TRUE, main = "Artefact P121")

panel(shapenorm4, main = "Scaled outline data", names = TRUE)

shapenorm4 <- shapenorm4 %>% coo_centre() %>% coo_slidedirection("right") 

saveRDS(shapenorm4, file = "Normalised_outlinesBC.rds") # save normalised and transformed landmarks

#### Reload and subset the data  RUN FROM HERE FROM NOW ON  ####

BC_outlines <- import("Normalised_outlinesBC.rds")
database <- import("BC_database.rds")  

panel(BC_outlines, main = "Border Cave", names = TRUE)
stack(BC_outlines)

## Calculate mean number of points representing outlines
meanpoints <- matrix(0, nrow = length(BC_outlines), ncol = 1)
for(i in 1:length(BC_outlines)){
  artefact <- BC_outlines[i]
  nlandmarks <- length(unlist(artefact))/2
  meanpoints[i,] <- nlandmarks
}

mean(meanpoints)

# Data tidying
BC_outlines$fac$Raw.materials.2[BC_outlines$fac$Raw.materials.2=="other "] <- "other"
table(BC_outlines$fac$Raw.materials.2) 
BC_outlines$fac$Raw.materials.2 <- droplevels(BC_outlines$fac$Raw.materials.2)
table(BC_outlines$fac$Raw.materials.2) 

BC_outlines$fac$Proximal.dorsal.preparation[BC_outlines$fac$Proximal.dorsal.preparation==""] <- "N"
table(BC_outlines$fac$Proximal.dorsal.preparation) 
BC_outlines$fac$Proximal.dorsal.preparation <- droplevels(BC_outlines$fac$Proximal.dorsal.preparation)
table(BC_outlines$fac$Proximal.dorsal.preparation) 

BC_outlines$fac$Dorsal.scar.pattern[BC_outlines$fac$Dorsal.scar.pattern=="Unidirectional corvengent direct"] <- "Unidirectional convergent direct"
BC_outlines$fac$Dorsal.scar.pattern[BC_outlines$fac$Dorsal.scar.pattern=="Unidirectional direct convergent"] <- "Unidirectional convergent direct"
table(BC_outlines$fac$Dorsal.scar.pattern) 
BC_outlines$fac$Dorsal.scar.pattern <- droplevels(BC_outlines$fac$Dorsal.scar.pattern)
table(BC_outlines$fac$Dorsal.scar.pattern) 

#### EFA ####

calibrate_harmonicpower_efourier(BC_outlines, nb.h = 20, plot = FALSE) # 11 harmonics
calibrate_reconstructions_efourier(BC_outlines, range = 1:11)

efashape <- efourier(BC_outlines, nb.h = 11, norm = FALSE)

####  PCA  ####

pcashape <- PCA(efashape) 

# Plot 

plot.new()
gg <- PCcontrib(pcashape, nax = 1:3, plot = FALSE)
gg$gg + 
  geom_polygon(fill="gray", col="black") 


plot.new()
par(mfrow =c(2,2))
plot_PCA(pcashape, axes = c(1,2), morphospace_position = "range", zoom = 1, chull = FALSE, eigen = FALSE, legend = FALSE, palette=pal_manual(c("slategrey"))) %>% layer_points(cex = 1)  
plot_PCA(pcashape, axes = c(1,3),  morphospace_position = "range", zoom = 1, chull = FALSE, eigen = FALSE, legend = FALSE, palette=pal_manual(c("slategrey"))) %>% layer_points(cex = 1)  
plot_PCA(pcashape, axes = c(2,3), morphospace_position = "range", zoom = 1, chull = FALSE,  eigen = FALSE, legend = FALSE, palette=pal_manual(c("slategrey"))) %>% layer_points(cex = 1) 

vp.BottomRight <- viewport(height=unit(.5, "npc"), width=unit(0.5, "npc"), 
                           just=c("left","top"), 
                           y=0.5, x=0.5)

p1 <- scree_plot(pcashape, nax =1:10) # PC1-3 gives over 95% of cum variance in the data, PC1 = 79% of variance
p1 <- p1  + theme_minimal()
print(p1, vp = vp.BottomRight) 

### Identify and remove outliers #### 

outliers <- which_out(efashape, conf = 0.01) # 0 outliers

par(mfrow=c(1,2))
#  plot outliers
for(i in 1:length(outliers)){
  outlier <- outliers[i]
  coo_plot(BC_outlines[outliers[i]], col = "grey", centroid = TRUE, main = database[outliers[i], 3])
}


# efashape <- Momocs::slice(efashape, -outliers) # to remove them  
# pcashape <- PCA(efashape) 


#### Size analysis ####

## Build new database with PCs and centroid size
centroidsize <- as_tibble(coo_centsize(BC_outlines))
centroidsize <- rename(centroidsize, CS = "value")

pcascores <- as_tibble(pcashape$x)
databasedata <- cbind(BC_outlines$fac, centroidsize, pcascores) # new database with PCs and centroid size

## Length
p1 <- ggscatter(databasedata, x = "Length", y = "CS", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Length (mm)", ylab = "Centroid size")


cor.test(databasedata$CS, databasedata$Length)
summary(lm(Length~CS, data = databasedata))

## Width
p2 <- ggscatter(databasedata, x = "Width", y = "CS", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Width (mm)", ylab = "Centroid size")
cor.test(databasedata$CS, databasedata$Width)
summary(lm(Width~CS, data = databasedata))

## PC1
p3 <- ggscatter(databasedata, x = "PC1", y = "CS", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "PC1", ylab = "Centroid size")

cor.test(databasedata$PC1, databasedata$CS)
cor(databasedata$PC1, databasedata$CS)

summary(lm(CS~PC1, data = databasedata))

## PC2
p4 <- ggscatter(databasedata, x = "PC2", y = "CS", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "PC2", ylab = "Centroid size")

cor.test(databasedata$PC2, databasedata$CS)
cor(databasedata$PC2, databasedata$CS)

summary(lm(CS~PC2, data = databasedata))

## PC3
p5 <- ggscatter(databasedata, x = "PC3", y = "CS", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "PC3", ylab = "Centroid size")

cor.test(databasedata$PC3, databasedata$CS)
cor(databasedata$PC3, databasedata$CS)

summary(lm(CS~PC3, data = databasedata))

(p3| p4| p5) / (p1| p2 ) +
  plot_annotation(title = "Allometry", tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 24))) 


#### PCA plots ####
tidy.table <- cbind(as.tibble(pcashape$x[,1:3]), as.tibble(database))

tidy.table$Member[tidy.table$Member=="2BS"] <- "2 BS"
tidy.table$Member[tidy.table$Member=="2WA"] <- "2 WA"
tidy.table$Member[tidy.table$Member=="3BS"] <- "3 BS"

tidy.table$New.Strat[tidy.table$New.strat=="2BS.UP"] <- "2 BS.UP"
tidy.table$New.Stra[tidy.table$New.strat=="2BS.LR"] <- "2 BS.LR"
tidy.table$New.Stra[tidy.table$New.strat=="2WA.UP"] <- "2 WA.UP"
tidy.table$New.Stra[tidy.table$New.strat=="2WA.MD"] <- "2 WA.MD"
tidy.table$New.Stra[tidy.table$New.strat=="2WA.LR"] <- "2 WA.LR"
tidy.table$New.Stra[tidy.table$New.strat=="3BS"] <- "3 BS"
tidy.table$New.Stra[tidy.table$New.strat=="1RGBS"] <- "1 RGBS "

## Member

a <- ggplot(tidy.table, aes(PC1, PC2, colour = Member)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99", "#6699CC")) +
  theme_minimal() +
  theme(legend.position = "none")

b <- ggplot(tidy.table, aes(PC1, PC3, colour = Member)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99", "#6699CC")) +
  theme_minimal() +
  theme(legend.position = "none")

c <- ggplot(tidy.table, aes(PC2, PC3, colour = Member)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99", "#6699CC")) +
  theme_minimal() +
  theme(legend.position = "none")

d <- tidy.table %>% 
  dplyr::select(PC1, PC2, PC3, Member) %>%
  pivot_longer(cols = contains("PC")) %>%
  mutate(name = as.factor(name)) %>%
  mutate(Member = as.factor(Member)) %>%
  ggplot(aes(name, value, fill = Member)) +
  scale_fill_manual(values = c("#AA4499", "#44AA99", "#6699CC")) +
  labs(y = "Score") +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())

(a | b | c) / d +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 30))) 

## New chronology
# These are not in the right order (atm alphabetical but should be by stratrigraphy above) and the legend needs to be wider

new_strat <- tidy.table %>% 
  dplyr::select(PC1, PC2, PC3, New.strat) %>%
  pivot_longer(cols = contains("PC")) %>%
  mutate(name = as.factor(name)) %>%
  mutate(New.strat = as.factor(New.strat)) %>%
  mutate(New.strat = fct_relevel(New.strat, "2 BS.UP", "2 BS.LR", "2 WA.UP", "2 WA.MD", "2 WA.LR", "3 BS", "1 RGBS")) %>%
  ggplot(aes(name, value, fill = New.strat)) +
  scale_fill_manual(guide = guide_legend(nrow = 1, byrow = TRUE), limits = c("2 BS.UP", "2 BS.LR", "2 WA.UP", "2 WA.MD", "2 WA.LR", "3 BS", "1 RGBS"),values = c("#44AA99","#AA4499", "#6699CC", "#999933", "#888888", "#000000", "#000000")) +
  labs(y = "Score") +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())

rel_strat <- tidy.table %>% 
  dplyr::select(PC1, PC2, PC3, Relative.stratigraphy) %>%
  pivot_longer(cols = contains("PC")) %>%
  mutate(name = as.factor(name)) %>%
  mutate(Relative.stratigraphy = as.factor(Relative.stratigraphy)) %>%
  ggplot(aes(name, value, fill = Relative.stratigraphy)) +
  scale_fill_manual(guide = guide_legend(nrow = 1, byrow = TRUE), values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888", "#E69F00", "#F0E442", "#CC79A7", "#009E73", "#000000", "#000000", "0072B2",  "#56B4E9", "#D55E00", "#000000", "#000000")) +
  labs(y = "Score") +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())

new_strat / rel_strat +
  plot_annotation(title = "Chronology", tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 24))) 

## Reduction method

a <- ggplot(tidy.table, aes(PC1, PC2, colour = Knapping.modality)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  theme_minimal() +
  theme(legend.position = "none")

b <- ggplot(tidy.table, aes(PC1, PC3, colour = Knapping.modality)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  theme_minimal() +
  theme(legend.position = "none")

c <- ggplot(tidy.table, aes(PC2, PC3, colour = Knapping.modality)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  theme_minimal() +
  theme(legend.position = "none")

d <- tidy.table %>% 
  dplyr::select(PC1, PC2, PC3, Knapping.modality) %>%
  pivot_longer(cols = contains("PC")) %>%
  mutate(name = as.factor(name)) %>%
  mutate(Knapping.modality = as.factor(Knapping.modality)) %>%
  ggplot(aes(name, value, fill = Knapping.modality)) +
  scale_fill_manual(values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  labs(y = "Score") +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())

(a | b | c) / d +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 30))) 


## Number of scars 

a <- ggplot(tidy.table, aes(PC1, PC2, colour = Number.of.scars)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(limits = c("2S", "3S-Low", "3S-Mesial", "3S-High", "M"),values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  theme_minimal() +
  theme(legend.position = "none")

b <- ggplot(tidy.table, aes(PC1, PC3, colour = Number.of.scars)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(limits = c("2S", "3S-Low", "3S-Mesial", "3S-High", "M"),values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  theme_minimal() +
  theme(legend.position = "none")

c <- ggplot(tidy.table, aes(PC2, PC3, colour = Number.of.scars)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(limits = c("2S", "3S-Low", "3S-Mesial", "3S-High", "M"),values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  theme_minimal() +
  theme(legend.position = "none")

d <- tidy.table %>% 
  dplyr::select(PC1, PC2, PC3, Number.of.scars) %>%
  pivot_longer(cols = contains("PC")) %>%
  mutate(name = as.factor(name)) %>%
  mutate(Number.of.scars = as.factor(Number.of.scars)) %>%
  mutate(Number.of.scars = fct_relevel(Number.of.scars, "2S", "3S-Low", "3S-Mesial", "3S-High", "M")) %>%
  ggplot(aes(name, value, fill = Number.of.scars)) +
  scale_fill_manual(limits = c("2S", "3S-Low", "3S-Mesial", "3S-High", "M"),values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  labs(y = "Score") +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())

(a | b | c) / d +
  plot_annotation(title = "Number of scars", tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 24))) 


## Retouch
a <- ggplot(tidy.table, aes(PC1, PC2, colour = Retouch)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99", "#6699CC")) +
  theme_minimal() +
  theme(legend.position = "none")

b <- ggplot(tidy.table, aes(PC1, PC3, colour = Retouch)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99", "#6699CC")) +
  theme_minimal() +
  theme(legend.position = "none")

c <- ggplot(tidy.table, aes(PC2, PC3, colour = Retouch)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99", "#6699CC")) +
  theme_minimal() +
  theme(legend.position = "none")

d <- tidy.table %>% 
  dplyr::select(PC1, PC2, PC3, Retouch) %>%
  pivot_longer(cols = contains("PC")) %>%
  mutate(name = as.factor(name)) %>%
  mutate(Retouch = as.factor(Retouch)) %>%
  ggplot(aes(name, value, fill = Retouch)) +
  scale_fill_manual(values = c("#AA4499", "#44AA99", "#6699CC")) +
  labs(y = "Score") +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())

(a | b | c) / d +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 30))) 

## Orientation of tip
a <- ggplot(tidy.table, aes(PC1, PC2, colour = Orientation.of.the.tip)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99", "#6699CC")) +
  theme_minimal() +
  theme(legend.position = "none")

b <- ggplot(tidy.table, aes(PC1, PC3, colour = Orientation.of.the.tip)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99", "#6699CC")) +
  theme_minimal() +
  theme(legend.position = "none")

c <- ggplot(tidy.table, aes(PC2, PC3, colour = Orientation.of.the.tip)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99", "#6699CC")) +
  theme_minimal() +
  theme(legend.position = "none")

d <- tidy.table %>% 
  dplyr::select(PC1, PC2, PC3, Orientation.of.the.tip) %>%
  pivot_longer(cols = contains("PC")) %>%
  mutate(name = as.factor(name)) %>%
  mutate(Orientation.of.the.tip = as.factor(Orientation.of.the.tip)) %>%
  ggplot(aes(name, value, fill = Orientation.of.the.tip)) +
  scale_fill_manual(values = c("#AA4499", "#44AA99", "#6699CC")) +
  labs(y = "Score") +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())

(a | b | c) / d +
  plot_annotation( tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 30))) 

## PDP

a <- ggplot(tidy.table, aes(PC1, PC2, colour = Proximal.dorsal.preparation)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99")) +
  theme_minimal() +
  theme(legend.position = "none")

b <- ggplot(tidy.table, aes(PC1, PC3, colour = Proximal.dorsal.preparation)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99")) +
  theme_minimal() +
  theme(legend.position = "none")

c <- ggplot(tidy.table, aes(PC2, PC3, colour = Proximal.dorsal.preparation)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(values = c("#AA4499", "#44AA99")) +
  theme_minimal() +
  theme(legend.position = "none")

d <- tidy.table %>% 
  dplyr::select(PC1, PC2, PC3, Proximal.dorsal.preparation) %>%
  pivot_longer(cols = contains("PC")) %>%
  mutate(name = as.factor(name)) %>%
  mutate(Proximal.dorsal.preparation = as.factor(Proximal.dorsal.preparation)) %>%
  ggplot(aes(name, value, fill = Proximal.dorsal.preparation)) +
  scale_fill_manual(values = c("#AA4499", "#44AA99")) +
  labs(y = "Score") +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())

(a | b | c) / d +
  plot_annotation(title = "Proximal dorsal preparation", tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 24))) 

## Ridges and nodes/dorsal scar pattern

rn <- tidy.table %>% 
  dplyr::select(PC1, PC2, PC3, Ridges.and.Nodes) %>%
  pivot_longer(cols = contains("PC")) %>%
  mutate(name = as.factor(name)) %>%
  mutate(Ridges.and.Nodes = as.factor(Ridges.and.Nodes)) %>%
  mutate(Ridges.and.Nodes = fct_relevel(Ridges.and.Nodes, "1N", "1N+1R", "1R","2N", "2R", "3N", "4N", "8N", "Ind")) %>%
  ggplot(aes(name, value, fill = Ridges.and.Nodes)) +
  scale_fill_manual(guide = guide_legend(nrow =1, byrow = TRUE), limits = c("1N", "1N+1R", "1R","2N", "2R", "3N", "4N", "8N", "Ind"), values = c("#44AA99","#000000","#AA4499", "#6699CC", "#000000", "#999933", "#44AA99","#000000","#888888")) +
  labs(y = "Score") +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())+
  ggtitle("Ridges and nodes")


dsp <- tidy.table %>% 
  dplyr::select(PC1, PC2, PC3, Dorsal.scar.pattern) %>%
  pivot_longer(cols = contains("PC")) %>%
  mutate(name = as.factor(name)) %>%
  mutate(Dorsal.scar.pattern = as.factor(Dorsal.scar.pattern)) %>%
  mutate(Dorsal.scar.pattern = fct_relevel(Dorsal.scar.pattern, "Unidirectional direct","Unidirectional convergent direct","Bidirectional", "Transversal right", "Transversal Unidirectional direct", "Multidirectional", "Ind")) %>%
  ggplot(aes(name, value, fill = Dorsal.scar.pattern)) +
  scale_fill_manual(guide = guide_legend(nrow = 2, byrow = TRUE), limits = c("Unidirectional direct","Unidirectional convergent direct","Bidirectional", "Transversal right", "Transversal Unidirectional direct", "Multidirectional", "Ind"), values = c("#999933", "#AA4499", "#888888","#000000","#000000","#6699CC", "#44AA99")) +
  labs(y = "Score") +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())+
  ggtitle("Dorsal scar pattern")

rn / dsp +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 24))) 


## Angle

a <- ggplot(tidy.table, aes(PC1, PC2, colour =Angle)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(limits = c("<60", "60-80", "80", ">80", "Ind"), values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  theme_minimal() +
  theme(legend.position = "none")

b <- ggplot(tidy.table, aes(PC1, PC3, colour = Angle)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(limits = c("<60", "60-80", "80", ">80", "Ind"), values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  theme_minimal() +
  theme(legend.position = "none")

c <- ggplot(tidy.table, aes(PC2, PC3, colour = Angle)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(limits = c("<60", "60-80", "80", ">80", "Ind"), values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  theme_minimal() +
  theme(legend.position = "none")

d <- tidy.table %>% 
  dplyr::select(PC1, PC2, PC3, Angle) %>%
  pivot_longer(cols = contains("PC")) %>%
  mutate(name = as.factor(name)) %>%
  mutate(Angle = as.factor(Angle)) %>%
  mutate(Angle = fct_relevel(Angle, "<60", "60-80", "80", ">80", "Ind")) %>%
  ggplot(aes(name, value, fill = Angle)) +
  scale_fill_manual(limits = c("<60", "60-80", "80", ">80", "Ind"), values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  labs(y = "Score") +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())

(a | b | c) / d +
  plot_annotation(title = "Platform angle", tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 24))) 

#### Discriminant analysis ####

## Member
table(databasedata$Member) # priors
dashapefc <- LDA(efashape, ~Member, prior = c(36, 16, 2)/nrow(databasedata, cv = TRUE)) # Fourier coefficient (raw data)
dashapefc$CV.correct
dashapefc$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 2, chull = FALSE) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

dashape95 <- LDA(pcashape, ~Member, prior = c(36, 16, 2)/nrow(databasedata), retain = 0.95, cv = TRUE) # 95% cum var PC scores
dashape95$CV.correct
dashape95$CV.ce

## New strat
table(databasedata$New.strat) # priors
dashapefc <- LDA(efashape, ~New.strat, prior = c(1,34, 2, 7, 3, 6, 1)/nrow(databasedata, cv = TRUE)) # Fourier coefficient (raw data)
dashapefc$CV.correct
dashapefc$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 1, chull = FALSE, eigen = FALSE, legend = TRUE, palette = col_summer2) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

dashape95 <- LDA(pcashape, ~New.strat, prior =c(1,34, 2, 7, 3, 6, 1)/nrow(databasedata), retain = 0.95, cv = TRUE) # 95% cum var PC scores
dashape95$CV.correct
dashape95$CV.ce

## Raw material
table(databasedata$Raw.materials) # priors

dashapefc <- LDA(efashape, ~Raw.materials, prior = c(7, 7, 4, 2, 1, 1, 7, 1, 6, 3, 14,1)/nrow(databasedata), cv = TRUE) # Fourier coefficient (raw data)
dashapefc$CV.correct
dashapefc$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 1, chull = FALSE, eigen = FALSE, legend = TRUE, palette = col_summer2) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

dashape95 <- LDA(pcashape, ~Raw.materials,  prior = c(7, 7, 4, 2, 1, 1, 7, 1, 6, 3, 14,1)/nrow(databasedata), retain = 0.95, cv = TRUE) # 95% cum var PC scores
dashape95$CV.correct
dashape95$CV.ce

## Retouch
table(databasedata$Retouch) # priors

dashapefc <- LDA(efashape, ~Retouch, prior = c(5,6,43)/nrow(databasedata), cv = TRUE) # Fourier coefficient (raw data)
dashapefc$CV.correct
dashapefc$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 1, chull = FALSE, eigen = FALSE, legend = TRUE, palette = col_summer2) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

dashape95 <- LDA(pcashape, ~Retouch,  prior =c(5,6,43)/nrow(databasedata), retain = 0.95, cv = TRUE) # 95% cum var PC scores
dashape95$CV.correct
dashape95$CV.ce

## Reduction 
table(databasedata$Knapping.modality) # priors

dashapefc <- LDA(efashape, ~Knapping.modality, prior = c(18,8,11,15,2)/nrow(databasedata), cv = TRUE) # Fourier coefficient (raw data)
dashapefc$CV.correct
dashapefc$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 1, chull = FALSE, eigen = FALSE, legend = TRUE, palette = col_summer2) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

dashape95 <- LDA(pcashape, ~Knapping.modality,  prior = c(18,8,11,15,2)/nrow(databasedata), retain = 0.95, cv = TRUE) # 95% cum var PC scores
dashape95$CV.correct
dashape95$CV.ce

## Number of scars
table(databasedata$Number.of.scars) # priors

dashapefc <- LDA(efashape, ~Number.of.scars, prior = c(11,5,2,5,31)/nrow(databasedata), cv = TRUE) # Fourier coefficient (raw data)
dashapefc$CV.correct
dashapefc$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 2, chull = FALSE) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

dashape95 <- LDA(pcashape, ~Number.of.scars, prior = c(11,5,2,5,31)/nrow(databasedata), cv = TRUE, retain = 0.95) # 95% cum var PC scores
dashape95$CV.correct
dashape95$CV.ce

## Orientation of tip
table(databasedata$Orientation.of.the.tip) # priors

dashapefc <- LDA(efashape, ~Orientation.of.the.tip, prior = c(34, 12, 8)/nrow(databasedata), cv = TRUE) # Fourier coefficient (raw data)
dashapefc$CV.correct
dashapefc$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 2, chull = FALSE) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

dashape95 <- LDA(pcashape, ~Orientation.of.the.tip, prior = c(34, 12, 8)/nrow(databasedata), cv = TRUE, retain = 0.95) # 95% cum var PC scores
dashape95$CV.correct
dashape95$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 1, chull = FALSE) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

## Proximal dorsal preparation
table(databasedata$Proximal.dorsal.preparation) # priors

dashapefc <- LDA(efashape, ~Proximal.dorsal.preparation, prior = c(42,12)/nrow(databasedata), cv = TRUE) # Fourier coefficient (raw data)
dashapefc$CV.correct
dashapefc$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 2, chull = FALSE) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

dashape95 <- LDA(pcashape, ~Proximal.dorsal.preparation, prior = c(42,12)/nrow(databasedata), cv = TRUE, retain = 0.95) # 95% cum var PC scores
dashape95$CV.correct
dashape95$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 1, chull = FALSE) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

## Ridges and Nodes
table(databasedata$Ridges.and.Nodes) # priors

dashapefc <- LDA(efashape, ~Ridges.and.Nodes, prior = c(13, 1, 6, 11, 1, 12, 2, 1, 7)/nrow(databasedata), cv = TRUE) # Fourier coefficient (raw data)
dashapefc$CV.correct
dashapefc$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 2, chull = FALSE) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

dashape95 <- LDA(pcashape, ~Ridges.and.Nodes, prior = c(13, 1, 6, 11, 1, 12, 2, 1, 7)/nrow(databasedata), cv = TRUE, retain = 0.95) # 95% cum var PC scores
dashape95$CV.correct
dashape95$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 1, chull = FALSE) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

## Dorsal scar pattern
table(databasedata$Dorsal.scar.pattern) # priors

dashapefc <- LDA(efashape, ~Dorsal.scar.pattern, prior = c(1, 13, 25, 1, 1, 7, 6)/nrow(databasedata), cv = TRUE) # Fourier coefficient (raw data)
dashapefc$CV.correct
dashapefc$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 2, chull = FALSE) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

dashape95 <- LDA(pcashape, ~Dorsal.scar.pattern, prior = c(1, 13, 25, 1, 1, 7, 6)/nrow(databasedata), cv = TRUE, retain = 0.95) # 95% cum var PC scores
dashape95$CV.correct
dashape95$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 1, chull = FALSE) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

## Angle
table(databasedata$Angle) # priors

dashapefc <- LDA(efashape, ~Angle, prior = c(4,5,20,1,24)/nrow(databasedata), cv = TRUE) # Fourier coefficient (raw data)
dashapefc$CV.correct
dashapefc$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 2, chull = FALSE) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

dashape95 <- LDA(pcashape, ~Angle, prior = c(4,5,20,1,24)/nrow(databasedata), cv = TRUE, retain = 0.95) # 95% cum var PC scores
dashape95$CV.correct
dashape95$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 1, chull = FALSE) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

#### Regression of shape against time #### 
databasedata$Relative.stratigraphy <- as.numeric(databasedata$Relative.stratigraphy) # turn into numeric

# PC1 
ggplot(databasedata, aes(PC1, Relative.stratigraphy)) + 
  geom_smooth(method = lm, col = "red", se = FALSE) + 
  geom_point(size = 2, pch = 16, alpha = 0.4) + 
  theme(text = element_text(size = 8), axis.text = element_text(size = 8)) + 
  ylab("Time")

cor.test(databasedata$PC1, databasedata$Relative.stratigraphy)
lmPC1 <- lm(PC1 ~Relative.stratigraphy, data = databasedata)
summary(lmPC1)


ggscatter(databasedata, x = "PC1", y = "Relative.stratigraphy", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "PC1", ylab = "Chronology")

# PC2
ggplot(databasedata, aes(PC2, Relative.stratigraphy)) + 
  geom_point(size = 2, pch = 16, alpha = 0.4) + 
  geom_smooth(method = lm, se = FALSE) + 
  theme(text = element_text(size = 8), axis.text = element_text(size = 8)) + 
  ylab("Time")

cor.test(databasedata$PC2, databasedata$Relative.stratigraphy)
lmPC2 <- lm(PC2 ~Relative.stratigraphy, data = databasedata)
summary(lmPC2)

# PC3
ggplot(databasedata, aes(PC3, Relative.stratigraphy)) + 
  geom_point(size = 2, pch = 16, alpha = 0.4) + 
  geom_smooth(method = lm, se = FALSE) + 
  theme(text = element_text(size = 8), axis.text = element_text(size = 8)) + 
  ylab("Time")

cor.test(databasedata$PC3, databasedata$Relative.stratigraphy)
lmPC3 <- lm(PC3 ~Relative.stratigraphy, data = databasedata)
summary(lmPC3)

#### MANOVA and ANOVA on raw data and PCs ####

## Member
efashape %>% MANOVA(~Member)
pcashape %>% MANOVA(~Member, retain = 0.95)

summary(aov(PC1~Member, databasedata))
summary(aov(PC2~Member, databasedata))
summary(aov(PC3~Member, databasedata))
res <- aov(PC3~Member, databasedata) ## extra Tukey HSD as is significant 
TukeyHSD(res)

## New layer groupings
efashape %>% MANOVA(~New.strat)
pcashape %>% MANOVA(~New.strat, retain = 0.95)

summary(aov(PC1~New.strat, databasedata))
summary(aov(PC2~New.strat, databasedata))
summary(aov(PC3~New.strat, databasedata))

res <- aov(PC1~New.strat, databasedata) ## extra Tukey HSD as is significant 
TukeyHSD(res)

## Raw material
efashape %>% MANOVA(~Raw.materials)
pcashape %>% MANOVA(~Raw.materials, retain = 0.95)
summary(aov(PC1~Raw.materials, databasedata))
summary(aov(PC2~Raw.materials, databasedata))
summary(aov(PC3~Raw.materials, databasedata))

## Retouch
efashape %>% MANOVA(~Retouch)
pcashape %>% MANOVA(~Retouch, retain = 0.95)

summary(aov(PC1~Retouch, databasedata))
summary(aov(PC2~Retouch, databasedata))
summary(aov(PC3~Retouch, databasedata))

res <- aov(PC1~Retouch, databasedata) ## extra Tukey HSD as is significant 
TukeyHSD(res)

## Reduction method
efashape %>% MANOVA(~Knapping.modality)
pcashape %>% MANOVA(~Knapping.modality, retain = 0.95)

summary(aov(PC1~Knapping.modality, tidy.table))
summary(aov(PC2~Knapping.modality, tidy.table))
summary(aov(PC3~Knapping.modality, tidy.table))

res <- aov(PC1~Knapping.modality, tidy.table) ## extra Tukey HSD as is significant 
TukeyHSD(res)

res <- aov(PC3~Knapping.modality, tidy.table) ## extra Tukey HSD as is significant 
TukeyHSD(res)

## Number of scars
efashape %>% MANOVA(~Number.of.scars)
pcashape %>% MANOVA(~Number.of.scars, retain = 0.95)

## Tip orientation
efashape %>% MANOVA(~Orientation.of.the.tip)
pcashape %>% MANOVA(~Orientation.of.the.tip, retain = 0.95)

summary(aov(PC1~Orientation.of.the.tip, databasedata))
summary(aov(PC2~Orientation.of.the.tip, databasedata))
summary(aov(PC3~Orientation.of.the.tip, databasedata))

res <- aov(PC3~Orientation.of.the.tip, databasedata) ## Tukey HSD
TukeyHSD(res)

## Proximal dorsal preparation
efashape %>% MANOVA(~Proximal.dorsal.preparation)
pcashape %>% MANOVA(~Proximal.dorsal.preparation, retain = 0.95)

## Ridges and Nodes
efashape %>% MANOVA(~Ridges.and.Nodes)
pcashape %>% MANOVA(~Ridges.and.Nodes, retain = 0.95)
summary(aov(PC1~Ridges.and.Nodes, databasedata))
summary(aov(PC2~Ridges.and.Nodes, databasedata))
summary(aov(PC3~Ridges.and.Nodes, databasedata))

## Dorsal scar pattern
efashape %>% MANOVA(~Dorsal.scar.pattern)
pcashape %>% MANOVA(~Dorsal.scar.pattern, retain = 0.95)

## Angle
efashape %>% MANOVA(~Angle)
pcashape %>% MANOVA(~Angle, retain = 0.95)

#### Test effect of tip breakage on results ####
complete_outlines <- Momocs:: filter(BC_outlines, Tip.damage=="N")

efashape2 <- efourier(complete_outlines, nb.h = 9, norm = FALSE)
pcashape2 <- PCA(efashape2)

pcascores2 <- as_tibble(pcashape2$x)
databasedata2 <- cbind(complete_outlines$fac, pcascores2) # new database with PCs and centroid size

plot.new()
gg <- PCcontrib(pcashape2, nax = 1:3, plot = FALSE)
gg$gg + 
  geom_polygon(fill="gray", col="black") 

### MANOVA
pcashape %>% MANOVA(~Member, retain = 0.95)
pcashape2 %>% MANOVA(~Member, retain = 0.95)

## New stratigraphy
pcashape %>% MANOVA(~New.strat, retain = 0.95)
pcashape2 %>% MANOVA(~New.strat, retain = 0.95)

## Raw material
pcashape %>% MANOVA(~Raw.materials)
pcashape2 %>% MANOVA(~Raw.materials, retain = 0.95)

## Reduction
pcashape %>% MANOVA(~Knapping.modality, retain = 0.95)
pcashape2 %>% MANOVA(~Knapping.modality, retain = 0.95)

summary(aov(PC1~Knapping.modality, databasedata))
summary(aov(PC1~Knapping.modality, databasedata2))

summary(aov(PC2~Knapping.modality, databasedata))
summary(aov(PC2~Knapping.modality, databasedata2))

summary(aov(PC3~Knapping.modality, databasedata))
summary(aov(PC3~Knapping.modality, databasedata2))

## Retouch
pcashape %>% MANOVA(~Retouch)
pcashape2 %>% MANOVA(~Retouch, retain = 0.95)

summary(aov(PC1~Retouch, databasedata))
summary(aov(PC1~Retouch, databasedata2))

summary(aov(PC2~Retouch, databasedata))
summary(aov(PC2~Retouch, databasedata2))

summary(aov(PC3~Retouch, databasedata))
summary(aov(PC3~Retouch, databasedata2))

## Orientation of tip
pcashape %>% MANOVA(~Orientation.of.the.tip)
pcashape2 %>% MANOVA(~Orientation.of.the.tip, retain = 0.95)

summary(aov(PC1~Orientation.of.the.tip, databasedata))
summary(aov(PC1~Orientation.of.the.tip, databasedata2))

summary(aov(PC2~Orientation.of.the.tip, databasedata))
summary(aov(PC2~Orientation.of.the.tip, databasedata2))

summary(aov(PC3~Orientation.of.the.tip, databasedata))
summary(aov(PC3~Orientation.of.the.tip, databasedata2))




