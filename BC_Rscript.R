# Timbrell, L. et al. (2022). The post-Howiesons Poort points from Border Cave, KwaZulu Natal, South Africa. QSR. 
# This script performs the geometric morphometric analysis 

# Data available to download at: 

# To reproduce the analysis, download the entire folder on OSF as a zip file, extract it to a folder, and set the
# working directory below to the folder where you extracted the data.

# Clear R environment and set working directory
rm(list=ls())

setwd("/Volumes/Lucy HD/BorderCave_Sibudu")
#setwd("...")

## Install and load packages
if(!require("Momocs")) install.packages('Momocs', repos='http://cran.us.r-project.org')  
if(!require("tidyverse")) install.packages('tidyverse', repos='http://cran.us.r-project.org')
if(!require("rio")) install.packages('rio', repos='http://cran.us.r-project.org')
if(!require("ggplot2")) install.packages('ggplot2', repos='http://cran.us.r-project.org')
if(!require("ggpubr")) install.packages('ggpubr', repos='http://cran.us.r-project.org')

library("Momocs")
library("rio")
library("tidyverse")
library("ggplot2")
library("ggpubr")
library("gridExtra")
library("grid")

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

saveRDS(shapenorm4, file = "Normalised_outlinesBC.rds") # save normalised and transformed landmarks

#### Reload and subset the data  RUN FROM HERE FROM NOW ON ####

BC_outlines <- import("Normalised_outlinesBC.rds")
database <- import("BC_database.rds")  

panel(BC_outlines, main = "Border Cave", names = TRUE)

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

#### EFA ####

calibrate_harmonicpower_efourier(BC_outlines, nb.h = 20, plot = FALSE) # 11 harmonics
calibrate_reconstructions_efourier(BC_outlines, range = 1:11)

efashape <- efourier(BC_outlines, nb.h = 11, norm = FALSE)

####  PCA  ####

pcashape <- PCA(efashape) 

# Plot 
scree_plot(pcashape, nax =1:10) # PC1-3 gives over 95% of cum variance in the data, PC1 = 79% of variance

plot.new()
gg <- PCcontrib(pcashape, nax = 1:3, plot = FALSE)
par(mfrow = c(2,2))
plot_PCA(pcashape, axes = c(1,2), morphospace_position = "range", zoom = 1, chull = FALSE, eigen = FALSE, legend = FALSE, palette=pal_manual(c("slategrey"))) %>% layer_points(cex = 1)  
plot_PCA(pcashape, axes = c(1,3),  morphospace_position = "range", zoom = 1, chull = FALSE, eigen = FALSE, legend = FALSE, palette=pal_manual(c("slategrey"))) %>% layer_points(cex = 1)  
plot_PCA(pcashape, axes = c(2,3), morphospace_position = "range", zoom = 1, chull = FALSE,  eigen = FALSE, legend = FALSE, palette=pal_manual(c("slategrey"))) %>% layer_points(cex = 1) 

vp.BottomRight <- viewport(height=unit(.5, "npc"), width=unit(0.5, "npc"), 
                           just=c("left","top"), 
                           y=0.5, x=0.5)
GG <- gg$gg + 
  geom_polygon(fill="white", col="black") + 
  theme_classic() +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),  #remove y axis ticks
        axis.line.x = element_blank(),
        axis.line.y =element_blank())
print(GG, vp = vp.BottomRight)  


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


ggarrange(p1, p2, p3, p4, p5, ncol = 3, nrow = 2)

#### PCA plots ####
## Member
plot.new()
par(mfrow = c(2,2))
plot_PCA(pcashape, axes = c(1,2), ~ Member, morphospace_position = "range_axes", zoom = 1, chull = FALSE, palette=pal_manual(c("indianred3", "forestgreen", "royalblue2")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1)  
plot_PCA(pcashape, axes = c(1,3), ~ Member, morphospace_position = "range_axes", zoom = 1, chull = FALSE,palette=pal_manual(c("indianred3", "forestgreen", "royalblue2")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1)  
plot_PCA(pcashape, axes = c(2,3), ~ Member, morphospace_position = "range_axes", zoom = 1, chull = FALSE, palette=pal_manual(c("indianred3", "forestgreen", "royalblue2")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1) 

vp.BottomRight <- viewport(height=unit(.5, "npc"), width=unit(0.5, "npc"), 
                           just=c("left","top"), 
                           y=0.5, x=0.5)
p1 <- boxplot(pcashape, ~Member, nax = 1:3, col = c("indianred3", "forestgreen", "royalblue2"))
p1 <- p1  + theme_classic2()
print(p1, vp = vp.BottomRight)  

## New chronology
plot.new()

pcashape$fac$New.strat <- factor(pcashape$fac$New.strat, levels = c("2BS.UP", "2BS.LR", "2WA.UP", "2WA.MD", "2WA.LR", "3BS", "1RGBS"))

p1 <- boxplot(pcashape, ~ New.strat, nax = 1:3, lex.order = TRUE)
p1 <- p1 + ggtitle("Allostratrigraphic units") + theme_classic2()

p2 <- boxplot(pcashape, ~ Relative.stratigraphy, nax = 1:3, lex.order = FALSE)
p2 <- p2 + ggtitle("Relative stratigraphy")+ theme_classic2()

ggarrange(p1, p2, nrow = 2)

## Reduction method
plot.new()
par(mfrow = c(2,2))
plot_PCA(pcashape, axes = c(1,2), ~ Knapping.modality, morphospace_position = "range_axes", zoom = 1, chull = FALSE, palette=pal_manual(c("indianred3", "yellow4", "springgreen2", "deepskyblue2", "mediumorchid1")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1)  
plot_PCA(pcashape, axes = c(1,3), ~ Knapping.modality, morphospace_position = "range_axes", zoom = 1, chull = FALSE,palette=pal_manual(c("indianred3", "yellow4", "springgreen2", "deepskyblue2", "mediumorchid1")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1)  
plot_PCA(pcashape, axes = c(2,3), ~ Knapping.modality, morphospace_position = "range_axes", zoom = 1, chull = FALSE, palette=pal_manual(c("indianred3", "yellow4", "springgreen2", "deepskyblue2", "mediumorchid1")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1) 

vp.BottomRight <- viewport(height=unit(.5, "npc"), width=unit(0.5, "npc"), 
                           just=c("left","top"), 
                           y=0.5, x=0.5)
p1 <- boxplot(pcashape, ~Knapping.modality, nax = 1:3)
p1 <- p1  + theme_classic2()
print(p1, vp = vp.BottomRight) 

## Raw.materials 
plot.new()
par(mfrow = c(2,2))
plot_PCA(pcashape, axes = c(1,2), ~ Raw.materials, morphospace_position = "range_axes", zoom = 1, chull = FALSE, palette=pal_manual(c("indianred3", "yellow4", "springgreen2", "deepskyblue2", "mediumorchid1")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1)  
plot_PCA(pcashape, axes = c(1,3), ~ Raw.materials, morphospace_position = "range_axes", zoom = 1, chull = FALSE,palette=pal_manual(c("indianred3", "yellow4", "springgreen2", "deepskyblue2", "mediumorchid1")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1)  
plot_PCA(pcashape, axes = c(2,3), ~ Raw.materials, morphospace_position = "range_axes", zoom = 1, chull = FALSE, palette=pal_manual(c("indianred3", "yellow4", "springgreen2", "deepskyblue2", "mediumorchid1")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1) 

## Number of scars 
plot.new()
par(mfrow = c(2,2))
plot_PCA(pcashape, axes = c(1,2), ~ Number.of.scars, morphospace_position = "range_axes", zoom = 1 , chull = FALSE, palette=pal_manual(c("indianred3", "yellow4", "springgreen2", "deepskyblue2", "mediumorchid1")), eigen = FALSE) %>% layer_points(cex = 1.2, pch = 20) 
plot_PCA(pcashape, axes = c(1,3), ~ Number.of.scars, morphospace_position = "range_axes", zoom = 1, chull = FALSE, palette=pal_manual(c("indianred3", "yellow4", "springgreen2", "deepskyblue2", "mediumorchid1")), eigen = FALSE) %>% layer_points(cex = 1.2, pch = 20) 
plot_PCA(pcashape, axes = c(2,3), ~ Number.of.scars, morphospace_position = "range_axes", zoom = 1, chull = FALSE, palette=pal_manual(c("indianred3", "yellow4", "springgreen2", "deepskyblue2", "mediumorchid1")), eigen = FALSE) %>% layer_points(cex = 1.2, pch = 20) 
p1 <- boxplot(pcashape, ~Number.of.scars, nax = 1:3)
p1 <- p1  + theme_classic2()
print(p1, vp = vp.BottomRight) 

## Retouch
plot.new()
par(mfrow = c(2,2))
plot_PCA(pcashape, axes = c(1,2), ~ Retouch, morphospace_position = "range_axes", zoom = 1 , chull = FALSE, palette=pal_manual(c("indianred3", "forestgreen", "royalblue2")), eigen = FALSE) %>% layer_points(cex = 1.2, pch = 20) 
plot_PCA(pcashape, axes = c(1,3), ~ Retouch, morphospace_position = "range_axes", zoom = 1, chull = FALSE, palette=pal_manual(c("indianred3", "forestgreen", "royalblue2")), eigen = FALSE) %>% layer_points(cex = 1.2, pch = 20) 
plot_PCA(pcashape, axes = c(2,3), ~ Retouch, morphospace_position = "range_axes", zoom = 1, chull = FALSE, palette=pal_manual(c("indianred3", "forestgreen", "royalblue2")), eigen = FALSE) %>% layer_points(cex = 1.2, pch = 20) 
p1 <- boxplot(pcashape, ~Retouch, nax = 1:3)
p1 <- p1  + theme_classic2()
print(p1, vp = vp.BottomRight) 

## Orientation of tip
plot.new()
par(mfrow = c(2,2))
plot_PCA(pcashape, axes = c(1,2), ~ Orientation.of.the.tip, morphospace_position = "range_axes", zoom = 1, chull = FALSE, palette=pal_manual(c("indianred3", "forestgreen", "royalblue2")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1)  
plot_PCA(pcashape, axes = c(1,3), ~ Orientation.of.the.tip, morphospace_position = "range_axes", zoom = 1, chull = FALSE,palette=pal_manual(c("indianred3", "forestgreen", "royalblue2")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1)  
plot_PCA(pcashape, axes = c(2,3), ~ Orientation.of.the.tip, morphospace_position = "range_axes", zoom = 1, chull = FALSE, palette=pal_manual(c("indianred3", "forestgreen", "royalblue2")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1) 

vp.BottomRight <- viewport(height=unit(.5, "npc"), width=unit(0.5, "npc"), 
                           just=c("left","top"), 
                           y=0.5, x=0.5)
p1 <- boxplot(pcashape, ~Orientation.of.the.tip, nax = 1:3)
p1 <- p1  + theme_classic2()
print(p1, vp = vp.BottomRight) 

## PDP
plot.new()
par(mfrow = c(2,2))
plot_PCA(pcashape, axes = c(1,2), ~ Proximal.dorsal.preparation, morphospace_position = "range_axes", zoom = 1 , chull = FALSE, palette=pal_manual(c("indianred3", "royalblue2")), eigen = FALSE) %>% layer_points(cex = 1.2, pch = 20) 
plot_PCA(pcashape, axes = c(1,3), ~ Proximal.dorsal.preparation, morphospace_position = "range_axes", zoom = 1 , chull = FALSE, palette=pal_manual(c("indianred3", "royalblue2")), eigen = FALSE) %>% layer_points(cex = 1.2, pch = 20) 
plot_PCA(pcashape, axes = c(2,3), ~ Proximal.dorsal.preparation, morphospace_position = "range_axes", zoom = 1 , chull = FALSE, palette=pal_manual(c("indianred3", "royalblue2")), eigen = FALSE) %>% layer_points(cex = 1.2, pch = 20) 
p1 <- boxplot(pcashape, ~Proximal.dorsal.preparation, nax = 1:3)
p1 <- p1  + theme_classic2()
print(p1, vp = vp.BottomRight) 

## Ridges and nodes/dorsal scar pattern
plot.new()
par(mfrow = c(1,2))
p1 <- boxplot(pcashape, ~Ridges.and.Nodes, nax = 1:3)
p1 <- p1  + theme_classic2()

p2 <- boxplot(pcashape, ~Dorsal.scar.pattern, nax = 1:3)
p2 <- p2  + theme_classic2()

ggarrange(p1, p2, ncol=1)

## Angle
plot.new()
par(mfrow = c(2,2))
plot_PCA(pcashape, axes = c(1,2), ~ Angle, morphospace_position = "range_axes", zoom = 1 , chull = FALSE, palette=pal_manual(c("indianred3", "yellow4", "springgreen2", "deepskyblue2", "mediumorchid1")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1.2, pch = 20) 
plot_PCA(pcashape, axes = c(1,3), ~ Angle, morphospace_position = "range_axes", zoom = 1 , chull = FALSE, palette=pal_manual(c("indianred3", "yellow4", "springgreen2", "deepskyblue2", "mediumorchid1")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1.2, pch = 20) 
plot_PCA(pcashape, axes = c(2,3), ~ Angle, morphospace_position = "range_axes", zoom = 1 , chull = FALSE, palette=pal_manual(c("indianred3", "yellow4", "springgreen2", "deepskyblue2", "mediumorchid1")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1.2, pch = 20) 
p1 <- boxplot(pcashape, ~Angle, nax = 1:3)
p1 <- p1  + theme_classic2()
print(p1, vp = vp.BottomRight)

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

dashapefc <- LDA(efashape, ~Dorsal.scar.pattern, prior = c(1, 13, 25, 1, 1, 3, 1, 6, 3)/nrow(databasedata), cv = TRUE) # Fourier coefficient (raw data)
dashapefc$CV.correct
dashapefc$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 2, chull = FALSE) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

dashape95 <- LDA(pcashape, ~Dorsal.scar.pattern, prior = c(1, 13, 25, 1, 1, 3, 1, 6, 3)/nrow(databasedata), cv = TRUE, retain = 0.95) # 95% cum var PC scores
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

lmPC1 <- loess(abs(PC1) ~Relative.stratigraphy, data = databasedata)
summary(lmPC1)

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

## Retouch
efashape %>% MANOVA(~Retouch)
pcashape %>% MANOVA(~Retouch, retain = 0.95)

summary(aov(PC1~Retouch, databasedata))
summary(aov(PC2~Retouch, databasedata))
summary(aov(PC3~Retouch, databasedata))

res <- aov(PC3~Retouch, databasedata) ## extra Tukey HSD as is significant 
TukeyHSD(res)

## Reduction method
efashape %>% MANOVA(~Knapping.modality)
pcashape %>% MANOVA(~Knapping.modality, retain = 0.95)

summary(aov(PC1~Knapping.modality, databasedata))
summary(aov(PC2~Knapping.modality, databasedata))
summary(aov(PC3~Knapping.modality, databasedata))

res <- aov(PC3~Knapping.modality, databasedata) ## extra Tukey HSD as is significant 
TukeyHSD(res)

## Number of scars
efashape %>% MANOVA(~Number.of.scars)
pcashape %>% MANOVA(~Number.of.scars, retain = 0.95)

## Tip orientation
efashape %>% MANOVA(~Orientation.of.the.tip)
pcashape %>% MANOVA(~Orientation.of.the.tip, retain = 0.95)

res <- aov(PC3~Orientation.of.the.tip, databasedata) ## Tukey HSD
TukeyHSD(res)

## Proximal dorsal preparation
efashape %>% MANOVA(~Proximal.dorsal.preparation)
pcashape %>% MANOVA(~Proximal.dorsal.preparation, retain = 0.95)

## Ridges and Nodes
efashape %>% MANOVA(~Ridges.and.Nodes)
pcashape %>% MANOVA(~Ridges.and.Nodes, retain = 0.95)

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

gg <- PCcontrib(pcashape2, nax = 1:4)
p1 <- scree_plot(pcashape2, nax =1:10) # PC1-4 gives over 95% of cum variance in the data, PC1 = 77% of variance
p1 <- p1 + theme_bw() 
p2 <- gg$gg + 
  geom_polygon(fill="slategrey", col="black") + 
  ggtitle("Shapes along PC1-3") # PC1 = asymmetry with base heavy, PC2 = wide to elongated, PC3 - wide to elongated with base heavy
ggarrange(p1, p2, ncol =2)

### MANOVA

efashape2 %>% MANOVA(~Member)
pcashape2 %>% MANOVA(~Member, retain = 0.95)

## New stratigraphy
efashape2 %>% MANOVA(~New.strat)
pcashape2 %>% MANOVA(~New.strat, retain = 0.95)

## Raw material
efashape2 %>% MANOVA(~Raw.materials)
pcashape2 %>% MANOVA(~Raw.materials, retain = 0.95)

## Reduction
efashape2 %>% MANOVA(~Knapping.modality)
pcashape2 %>% MANOVA(~Knapping.modality, retain = 0.95)

## Retouch
efashape2 %>% MANOVA(~Retouch)
pcashape2 %>% MANOVA(~Retouch, retain = 0.95)

## Number of scars
efashape2 %>% MANOVA(~Number.of.scars)
pcashape2 %>% MANOVA(~Number.of.scars, retain = 0.95)

## Orientation of tip
efashape2 %>% MANOVA(~Orientation.of.the.tip)
pcashape2 %>% MANOVA(~Orientation.of.the.tip, retain = 0.95)

## Proximal dorsal preparation
efashape2 %>% MANOVA(~Proximal.dorsal.preparation)
pcashape2 %>% MANOVA(~Proximal.dorsal.preparation, retain = 0.95)

## Ridges and nodes
efashape2 %>% MANOVA(~Ridges.and.Nodes)
pcashape2 %>% MANOVA(~Ridges.and.Nodes, retain = 0.95)

## Dorsal scar pattern
efashape2 %>% MANOVA(~Dorsal.scar.pattern)
pcashape2 %>% MANOVA(~Dorsal.scar.pattern, retain = 0.95)

## Angle
efashape2 %>% MANOVA(~Angle)
pcashape2 %>% MANOVA(~Angle, retain = 0.95)


############## COMPARISON W/SIBUDU ############

remotes::install_github("yesdavid/outlineR",
                        dependencies = TRUE, force = TRUE)
library(outlineR)
setwd("/Volumes/Lucy HD/BorderCave_Sibudu/project_folder")

### set paths
inpath <- file.path(".", "prep")
outpath <- file.path(".", "prep", "single")
tpspath <- file.path(".", "raw","scaling.tps")

### seperate images
outlineR::separate_single_artefacts(inpath = inpath, 
                                    outpath = outpath)
### import scaling factor
tps_df <- outlineR::tps_to_df(tpspath)

### extract outlines according to scaling factor
single_outlines_list <- outlineR::get_outlines(outpath = outpath, 
                                               tps_file_rescale = tps_df)

### combine outlies into single Out object
SC_outlines_combined <- outlineR::combine_outlines(single_outlines_list = single_outlines_list)
saveRDS(SC_outlines_combined , file = "OutlinesSC.rds") # save normalised and transformed landmarks

#### Reload the data  RUN FROM HERE FROM NOW ON ####
setwd("/Volumes/Lucy HD/BorderCave_Sibudu/project_folder")

SC_outlines <- import("OutlinesSC.rds")
Morphotype <- c("Ndwndwe", "Ndwndwe","Ndwndwe","Ndwndwe","tongati", "tongati", "tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati","tongati")
SC_outlines  <- Out(SC_outlines$coo, fac = Morphotype)

### normalisation
SC_outlines <- SC_outlines %>% 
  coo_centre() %>% 
  coo_scale() %>% 
  coo_slidedirection("right") %>% 
  coo_close() 

panel(SC_outlines, main = "Outline data") # Visualization of points in their original orientation

### efa
SC_efa <- efourier(SC_outlines, nb.h = 11, norm = FALSE)

### predict PC scores from previous PCA
newPCA <- rePCA(pcashape, SC_efa)

### plot 
plot.new()
gg <- PCcontrib(pcashape, nax = 1:3, sd.r = c( -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3), plot = FALSE)

par(mfrow =c(2,2))
plot(pcashape$x[, c(1,2)], pch = 20, cex = 1, xlim = c(min(newPCA$x[,1]), max(newPCA$x[,1])), ylim = c(min(newPCA$x[,2]),max(newPCA$x[,2])))
points(newPCA$x[, c(1,2)], col= ifelse(newPCA$fac$value == "Ndwndwe", "red", "blue"), pch=20, cex = 1)
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

plot(pcashape$x[, c(1,3)], pch = 20, cex = 1, xlim = c(min(newPCA$x[,1]), max(newPCA$x[,1])), ylim = c(min(newPCA$x[,3]),max(newPCA$x[,3])))
points(newPCA$x[, c(1,3)], col= ifelse(newPCA$fac$value == "Ndwndwe", "red", "blue"), pch=20, cex = 1)
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

plot(pcashape$x[, c(2,3)], pch = 20, cex = 1, xlim = c(min(newPCA$x[,2]), max(newPCA$x[,2])), ylim = c(min(newPCA$x[,3]),max(newPCA$x[,3])))
points(newPCA$x[, c(2,3)], col= ifelse(newPCA$fac$value == "Ndwndwe", "red", "blue"), pch=20, cex = 1)
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

GG <- gg$gg + 
  geom_polygon(fill="white", col="black") + 
  theme_classic() +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),  #remove y axis ticks
        axis.line.x = element_blank(),
        axis.line.y =element_blank())
print(GG, vp = vp.BottomRight)  


centroidsize <- as_tibble(coo_centsize(SC_outlines))
centroidsize <- rename(centroidsize, CS = "value")

SC_pcascores <- as_tibble(newPCA$x)
SC_databasedata <- cbind(SC_outlines$fac, centroidsize, SC_pcascores) # new database with PCs and centroid size

setwd("/Volumes/Lucy HD/BorderCave_Sibudu")
write.csv(SC_databasedata, "Sibudu_pcscores.csv")
write.csv(databasedata, "BorderCave_pcscores.csv")
