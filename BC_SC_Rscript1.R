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

#### Data import and preparation ####

tpslines <- import_tps("BC_SC_outlinedata1.tps")

database <- read.csv("BC_Sibudu_data.csv") # Database

database$Tip.damage <- as.factor(database$Tip.damage) 
database$Site <- as.factor(database$Site) 
database$Layer <- as.factor(database$Layer)
database$Member <- as.factor(database$Member)
database$Relative.stratigraphy <- as.factor(database$Relative.stratigraphy)
database$Raw.materials <- as.factor(database$Raw.materials)
database$Raw.materials.2 <- as.factor(database$Raw.materials.2)
database$Knapping.modality <- as.factor(database$Knapping.modality)
database$New.strat <- as.factor(database$New.strat)

saveRDS(tpslines, file =  "tpslines.rds")
saveRDS(database, file =  "database.rds")

#### Import data and create outline object ####
tpsdata <- import("tpslines.rds") 
database <- import("database.rds")   

shape <- Out(tpsdata$coo, fac = database)
names(shape) <- database$Number

panel(shape, main = "Outline data", fac = "Site") # Visualization of points in their original orientation

#### Outline normalization, landmarking and artefact-specific processing ####

shapenorm <- shape %>% 
  coo_centre() %>% 
  coo_scale() %>% 
  coo_slidedirection("right") %>% 
  coo_untiltx() %>% 
  coo_close() %>% 
  coo_aligncalliper()

shapenorm1 <- def_ldk(shapenorm, 1) # add new landmarks to tip and base for additional orientation
shapenorm2 <- coo_untiltx(shapenorm1, ldk = 1) # option 2 - reduce tilt according to tip landmark - WE USE THIS

shapenorm3 <- coo_rotate(shapenorm2, theta = (pi/2)*3) # rotate so point is facing up
shapenorm4 <- coo_smooth(shapenorm3, 10)

# Artefact-specific processing
coo_plot(shapenorm3[14], col = "grey", centroid = TRUE, main = "Artefact 2297") # this outline needs smoothing
shapenorm4[14] <- coo_smooth(shapenorm4[14], 1000)  
coo_plot(shapenorm4[14], col = "grey", centroid = TRUE, main = "Artefact 2297")

coo_plot(shapenorm3[51], col = "grey", centroid = TRUE, main = "Artefact P121") # this outline needs smoothing
shapenorm4[51] <- coo_smooth(shapenorm4[51], 500)  
coo_plot(shapenorm4[51], col = "grey", centroid = TRUE, main = "Artefact P121")

stack(shapenorm4, title = "Stack:: Normalised Outlines")
panel(shapenorm4, main = "Scaled outline data", names = TRUE)

saveRDS(shapenorm4, file = "Normalised_outlinesBCSC.rds") # save normalised and transformed landmarks

#### Reload and subset the data  RUN FROM HERE FROM NOW ON ####

outlinedata <- import("Normalised_outlinesBCSC.rds")
database <- import("database.rds")  

## Calculate mean number of points representing outlines
meanpoints <- matrix(0, nrow = length(outlinedata), ncol = 1)
for(i in 1:length(outlinedata)){
  artefact <- outlinedata[i]
  nlandmarks <- length(unlist(artefact))/2
  meanpoints[i,] <- nlandmarks
}

mean(meanpoints)

## Seperate Border Cave and Sibudu 
BC_outlines <- Momocs:: filter(outlinedata, Site=="Border Cave")
SC_outlines <- Momocs:: filter(outlinedata, Site=="Sibudu")

BC_outlines$fac <- droplevels(BC_outlines$fac)
SC_outlines$fac <- droplevels(SC_outlines$fac)

############## BORDER CAVE ONLY ############
panel(BC_outlines, main = "Border Cave", names = TRUE)

# Data tidying
BC_outlines$fac$Raw.materials.2[BC_outlines$fac$Raw.materials.2=="other "] <- "other"
table(BC_outlines$fac$Raw.materials.2) # priors
BC_outlines$fac$Raw.materials.2 <- droplevels(BC_outlines$fac$Raw.materials.2)

#### EFA ####

calibrate_harmonicpower_efourier(BC_outlines, nb.h = 20, plot = FALSE) # 11 harmonics
calibrate_reconstructions_efourier(BC_outlines, range = 1:20)

efashape <- efourier(BC_outlines, nb.h = 11, norm = FALSE)

####  PCA  ####

pcashape <- PCA(efashape) 
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


# efashape1 <- Momocs::slice(efashape, -outliers) # to remove them  
# pcashape1 <- PCA(efashape1) # replace with efashape1 if want removed outliers

# BC_outlines <- Momocs::slice(BC_outlines, -outliers) # to remove them  
# scree_plot(pcashape1, nax =1:10) # PC1-3 gives over 95% of cum variance in the data, PC1 = 76% of variance

# gg <- PCcontrib(pcashape1, nax = 1:4)
# gg$gg + 
#  geom_polygon(fill="slategrey", col="black") + 
#  ggtitle("Shapes along PC1-4") # PC1 = asymmetry with base heavy, PC2 = wide to elongated, PC3 - wide to elongated with base heavy

#### Size analysis ####
centroidsize <- as_tibble(coo_centsize(BC_outlines))
centroidsize <- rename(centroidsize, CS = "value")

pcascores <- as_tibble(pcashape$x)
databasedata <- cbind(BC_outlines$fac, centroidsize, pcascores) # new database with PCs and centroid size

## Length
p1 <- ggplot(databasedata, aes(Length, CS)) + 
  geom_point(size = 2, pch = 16, alpha = 0.4) + 
  geom_smooth(method = lm, se = FALSE) + 
  theme(text = element_text(size = 8), axis.text = element_text(size = 8)) + 
  ylab("Centroid size")

cor.test(databasedata$CS, databasedata$Length)

summary(lm(Length~CS, data = databasedata))

## Width
p2 <- ggplot(databasedata, aes(Width, CS)) + 
  geom_point(size = 2, pch = 16, alpha = 0.4) + 
  geom_smooth(method = lm, se = FALSE) + 
  theme(text = element_text(size = 8), axis.text = element_text(size = 8)) + 
  ylab("Centroid size")

cor.test(databasedata$CS, databasedata$Width)

summary(lm(Width~CS, data = databasedata))

## PC1
p3 <- ggplot(databasedata, aes(PC1, CS)) + 
  geom_point(size = 2, pch = 16, alpha = 0.4) + 
  geom_smooth(method = lm, se = FALSE) + 
  theme(text = element_text(size = 8), axis.text = element_text(size = 8)) + 
  ylab("Centroid size")

cor.test(databasedata$PC1, databasedata$CS)
cor(databasedata$PC1, databasedata$CS)

summary(lm(CS~PC1, data = databasedata))

## PC2
p4 <- ggplot(databasedata, aes(PC2, CS)) + 
  geom_point(size = 2, pch = 16, alpha = 0.4) + 
  geom_smooth(method = lm, se = FALSE) + 
  theme(text = element_text(size = 8), axis.text = element_text(size = 8)) + 
  ylab("Centroid size")

cor.test(databasedata$PC2, databasedata$CS)
cor(databasedata$PC2, databasedata$CS)

summary(lm(CS~PC2, data = databasedata))

## PC3
p5 <- ggplot(databasedata, aes(PC3, CS)) + 
  geom_point(size = 2, pch = 16, alpha = 0.4) + 
  geom_smooth(method = lm, se = FALSE) + 
  theme(text = element_text(size = 8), axis.text = element_text(size = 8)) + 
  ylab("Centroid size")

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
plot_PCA(pcashape, axes = c(2,3), ~ Member, morphospace_position = "range", zoom = 1, chull = FALSE, palette=pal_manual(c("indianred3", "forestgreen", "royalblue2")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1) 

vp.BottomRight <- viewport(height=unit(.5, "npc"), width=unit(0.5, "npc"), 
                           just=c("left","top"), 
                           y=0.5, x=0.5)
p1 <- boxplot(pcashape, ~Member, nax = 1:3, col = c("indianred3", "forestgreen", "royalblue2"))
print(p1, vp = vp.BottomRight)  

## New chronology
plot.new()
p1 <- boxplot(pcashape, ~ Relative.stratigraphy, nax = 1:3, lex.order = FALSE)
p1 <- p1 + ggtitle("Relative stratigraphy")

p2 <- boxplot(pcashape, ~ New.strat, nax = 1:3, lex.order = FALSE)
p2 <- p2 + ggtitle("New layer groups")
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
p1 <- theme_bw()
print(p1, vp = vp.BottomRight) 

## Reduction method
plot.new()
par(mfrow = c(2,2))
plot_PCA(pcashape, axes = c(1,2), ~ Raw.materials, morphospace_position = "range_axes", zoom = 1, chull = FALSE, palette=pal_manual(c("indianred3", "yellow4", "springgreen2", "deepskyblue2", "mediumorchid1")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1)  
plot_PCA(pcashape, axes = c(1,3), ~ Knapping.modality, morphospace_position = "range_axes", zoom = 1, chull = FALSE,palette=pal_manual(c("indianred3", "yellow4", "springgreen2", "deepskyblue2", "mediumorchid1")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1)  
plot_PCA(pcashape, axes = c(2,3), ~ Knapping.modality, morphospace_position = "range_axes", zoom = 1, chull = FALSE, palette=pal_manual(c("indianred3", "yellow4", "springgreen2", "deepskyblue2", "mediumorchid1")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1) 


## Number of scars
boxplot(pcashape, ~ Number.of.scars, nax = 1:3)  #  No difference in mean along PC1 between the two post-HP members, 3BS is different 

plot.new()
plot_PCA(pcashape, axes = c(1,2), ~ Number.of.scars, morphospace_position = "range", zoom = 2 , chull = FALSE, palette=pal_manual(c("indianred3", "forestgreen", "royalblue2"))) %>% layer_points(cex = 1.2, pch = 20) 
plot_PCA(pcashape, axes = c(1,3), ~ Number.of.scars, morphospace_position = "range", zoom = 2, chull = FALSE, palette=pal_manual(c("indianred3", "forestgreen", "royalblue2"))) %>% layer_points(cex = 1.2, pch = 20) 

## Orientation of tip
boxplot(pcashape, ~ Orientation.of.the.tip, nax = 1:3)  #  No difference in mean along PC1 between the two post-HP members, 3BS is different 

plot.new()
plot_PCA(pcashape, axes = c(1,2), ~ Orientation.of.the.tip, morphospace_position = "range", zoom = 2 , chull = FALSE, palette=pal_manual(c("indianred3", "forestgreen", "royalblue2"))) %>% layer_points(cex = 1.2, pch = 20) 
plot_PCA(pcashape1, axes = c(1,3), ~ Orientation.of.the.tip, morphospace_position = "range_axes", zoom = 1, chull = FALSE, palette=pal_manual(c("indianred3", "forestgreen", "royalblue2"))) %>% layer_points(cex = 1.2, pch = 20) %>% layer_ellipses(conf = 0.95)

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
dashapefc <- LDA(efashape, ~New.strat, prior = c(33, 2, 6, 2, 7, 4)/nrow(databasedata, cv = TRUE)) # Fourier coefficient (raw data)
dashapefc$CV.correct
dashapefc$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 1, chull = FALSE, eigen = FALSE, legend = TRUE, palette = col_summer2) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

dashape95 <- LDA(pcashape, ~New.strat, prior = c(33, 2, 6, 2, 7, 4)/nrow(databasedata), retain = 0.95, cv = TRUE) # 95% cum var PC scores
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

## Reduction 
table(databasedata$Knapping.modality) # priors

dashapefc <- LDA(efashape, ~Knapping.modality, prior = c(29,10,6,1,8)/nrow(databasedata), cv = TRUE) # Fourier coefficient (raw data)
dashapefc$CV.correct
dashapefc$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 1, chull = FALSE, eigen = FALSE, legend = TRUE, palette = col_summer2) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

dashape95 <- LDA(pcashape, ~Knapping.modality,  prior = c(29,10,6,1,8)/nrow(databasedata), retain = 0.95, cv = TRUE) # 95% cum var PC scores
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

#### Regression of shape against time #### 
databasedata$Relative.stratigraphy <- as.numeric(databasedata$Relative.stratigraphy) # turn into numeric

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

ggplot(databasedata, aes(PC2, Relative.stratigraphy)) + 
  geom_point(size = 2, pch = 16, alpha = 0.4) + 
  geom_smooth(method = lm, se = FALSE) + 
  theme(text = element_text(size = 8), axis.text = element_text(size = 8)) + 
  ylab("Time")

cor.test(databasedata$PC2, databasedata$Relative.stratigraphy)
lmPC2 <- lm(PC2 ~Relative.stratigraphy, data = databasedata)
summary(lmPC2)

ggplot(databasedata, aes(PC3, Relative.stratigraphy)) + 
  geom_point(size = 2, pch = 16, alpha = 0.4) + 
  geom_smooth(method = lm, se = FALSE) + 
  theme(text = element_text(size = 8), axis.text = element_text(size = 8)) + 
  ylab("Time")

cor.test(databasedata$PC3, databasedata$Relative.stratigraphy)
lmPC3 <- lm(PC3 ~Relative.stratigraphy, data = databasedata)
summary(lmPC3)

#### MANOVA and ANOVA on PCs ####

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

## Raw material
efashape %>% MANOVA(~Raw.materials)
pcashape %>% MANOVA(~Raw.materials, retain = 0.95)
pcashape %>% MANOVA(~Raw.materials, retain = 0.99)

summary(aov(PC1~Raw.materials, databasedata))
summary(aov(PC2~Raw.materials, databasedata))
summary(aov(PC3~Raw.materials, databasedata))

res <- aov(PC1~Raw.materials, databasedata)
TukeyHSD(res)

res <- aov(PC2~Raw.materials, databasedata)
TukeyHSD(res)

res <- aov(PC3~Raw.materials, databasedata)
TukeyHSD(res)

## Reduction method
efashape %>% MANOVA(~Knapping.modality)
pcashape %>% MANOVA(~Knapping.modality, retain = 0.95)
pcashape %>% MANOVA(~Knapping.modality, retain = 0.99)

res <- aov(PC1~Knapping.modality, databasedata)
TukeyHSD(res)

summary(aov(PC2~Knapping.modality, databasedata))

res <- aov(PC3~Knapping.modality, databasedata)
TukeyHSD(res)

## Number of scars
efashape %>% MANOVA(~Number.of.scars)
pcashape %>% MANOVA(~Number.of.scars, retain = 0.95)

#### Cluster analyses #### 

## Member
CLUST(pcashape, ~Member, dist_method = "euclidean", type = "horizontal", cex = 0.6, hclust_method = "complete", palette=pal_manual(c("indianred3", "forestgreen", "royalblue2")))
CLUST(pcashape, ~New.strat, dist_method = "euclidean", type = "horizontal", cex = 0.6, hclust_method = "complete", palette = col_summer2)
CLUST(pcashape, ~Raw.materials, dist_method = "euclidean", type = "horizontal", cex = 0.6, hclust_method = "complete", palette = col_summer2)
CLUST(pcashape, ~Knapping.modality, dist_method = "euclidean", type = "horizontal", cex = 0.6, hclust_method = "complete", palette = col_summer2)
CLUST(pcashape, ~Number.of.scars, dist_method = "euclidean", type = "horizontal", cex = 0.6, hclust_method = "complete", palette = col_summer2)



#### Constructing mean shapes ####
## Member
meanshape <- MSHAPES(efashape, ~Member)
plot_MSHAPES(meanshape, size = 1, palette=pal_manual(c( "blue", "forestgreen", "red")))

## Relative stratigraphy
meanshape <- MSHAPES(efashape, ~Relative.stratigraphy)
plot_MSHAPES(meanshape, size = 0.75)

## Raw material
meanshape <- MSHAPES(efashape, ~Raw.materials)
plot_MSHAPES(meanshape, size = 0.75)

#### Test effect of tip breakage on results ####
complete_outlines <- Momocs:: filter(BC_outlines, Tip.damage=="N")

efashape2 <- efourier(complete_outlines, nb.h = 9, norm = FALSE)
pcashape2 <- PCA(efashape2)

scree_plot(pcashape2, nax =1:10) # PC1-4 gives over 95% of cum variance in the data, PC1 = 77% of variance

gg <- PCcontrib(pcashape2, nax = 1:4)
gg$gg + 
  geom_polygon(fill="slategrey", col="black") + 
  ggtitle("Shapes along PC1-3") # PC1 = asymmetry with base heavy, PC2 = wide to elongated, PC3 - wide to elongated with base heavy

## Member
boxplot(pcashape2, ~Member, nax = 1:4)  #  No difference in mean along PC1 between the two post-HP members, 3BS is different 

plot.new()
plot_PCA(pcashape2, axes = c(1,2), ~ Member, morphospace_position = "range_axes", zoom = 1, chull = TRUE, palette = col_summer) %>% layer_points(cex = 1.2, pch = 20) 
plot_PCA(pcashape2, axes = c(1,2), ~ Member, morphospace_position = "range", zoom = 1, chull = FALSE, palette = col_summer) %>% layer_points(cex = 1.2, pch = 20) %>% layer_ellipses(conf = 0.95)

plot_PCA(pcashape2, axes = c(2,3), ~ Member, morphospace_position = "range_axes", zoom = 1, chull = TRUE, legend = T, palette = col_summer) %>% layer_points(cex = 1.2, pch = 20) 
plot_PCA(pcashape2, axes = c(2,3), ~ Member, morphospace_position = "range", zoom = 1, chull = FALSE, legend = T, palette = col_summer) %>% layer_points(cex = 1.2, pch = 20) %>% layer_ellipses(conf = 0.95)

plot_PCA(pcashape2, axes = c(1,3), ~ Member, morphospace_position = "range_axes", zoom = 1, chull = TRUE, legend = T, palette = col_summer) %>% layer_points(cex = 1.2, pch = 20) 
plot_PCA(pcashape2, axes = c(1,3), ~ Member, morphospace_position = "range", zoom = 1, chull = FALSE, legend = T, palette = col_summer) %>% layer_points(cex = 1.2, pch = 20) %>% layer_ellipses(conf = 0.95)

## Relative stratigraphy
boxplot(pcashape2, ~ Relative.stratigraphy, nax = 1:4, lex.order = FALSE)

plot_PCA(pcashape2, axes = c(1,2), ~ Relative.stratigraphy, morphospace_position = "range_axes", zoom = 1, chull = FALSE, legend = T) %>% layer_points(cex = 1) %>% layer_ellipses(conf = 0.95)
plot_PCA(pcashape2, axes = c(1,2), ~ Relative.stratigraphy, morphospace_position = "range", zoom = 1, chull = FALSE, legend = T) %>% layer_points(cex = 1) %>% layer_ellipses(conf = 0.95)

plot_PCA(pcashape2, axes = c(2,3), ~ Relative.stratigraphy, morphospace_position = "range_axes", zoom = 1, chull = FALSE, legend = T) %>% layer_points(cex = 1) %>% layer_ellipses(conf = 0.95)
plot_PCA(pcashape2, axes = c(2,3), ~ Relative.stratigraphy, morphospace_position = "range", zoom = 1, chull = FALSE, legend = T) %>% layer_points(cex = 1) %>% layer_ellipses(conf = 0.95)

plot_PCA(pcashape2, axes = c(1,3), ~ Relative.stratigraphy, morphospace_position = "range_axes", zoom = 1, chull = FALSE, legend = T) %>% layer_points(cex = 1) %>% layer_ellipses(conf = 0.95)
plot_PCA(pcashape2, axes = c(1,3), ~ Relative.stratigraphy, morphospace_position = "range", zoom = 1, chull = FALSE, legend = T) %>% layer_points(cex = 1) %>% layer_ellipses(conf = 0.95)

## Raw materials
boxplot(pcashape2, ~Raw.materials, nax = 1:5)

plot_PCA(pcashape2, axes = c(1,2), ~ Raw.materials, morphospace_position = "range_axes", zoom = 1, chull = FALSE, legend = T) %>% layer_points(cex = 1) %>% layer_ellipses(conf = 0.95)
plot_PCA(pcashape2, axes = c(1,2), ~ Raw.materials, morphospace_position = "range", zoom = 1, chull = FALSE, legend = T) %>% layer_points(cex = 1) %>% layer_ellipses(conf = 0.95)

plot_PCA(pcashape2, axes = c(2,3), ~ Raw.materials, morphospace_position = "range_axes", zoom = 1, chull = FALSE, legend = T) %>% layer_points(cex = 1) %>% layer_ellipses(conf = 0.95)
plot_PCA(pcashape2, axes = c(2,3), ~ Raw.materials, morphospace_position = "range", zoom = 1, chull = FALSE, legend = T) %>% layer_points(cex = 1) %>% layer_ellipses(conf = 0.95)

plot_PCA(pcashape2, axes = c(1,3), ~ Raw.materials, morphospace_position = "range_axes", zoom = 1, chull = FALSE, legend = T) %>% layer_points(cex = 1) %>% layer_ellipses(conf = 0.95)
plot_PCA(pcashape2, axes = c(1,3), ~ Raw.materials, morphospace_position = "range", zoom = 1, chull = FALSE, legend = T) %>% layer_points(cex = 1) %>% layer_ellipses(conf = 0.95)

### Discriminant analysis 

## Member
dashapefc <- LDA(efashape2, ~Member) # Fourier coefficient (raw data)
dashapefc$CV.correct
dashapefc$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 2, chull = FALSE) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

dashape95 <- LDA(pcashape2, ~Member, retain = 0.95) # 95% cum var PC scores
dashape95$CV.correct
dashape95$CV.ce

dashape99 <- LDA(pcashape2, ~Member, retain = 0.99) # 99% cum var PC scores
dashape99$CV.correct
dashape99$CV.ce

## Relative stratigraphy
dashapefc <- LDA(efashape2, ~Relative.stratigraphy) # Fourier coefficient (raw data)
dashapefc$CV.correct
dashapefc$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 2, chull = FALSE) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

dashape95 <- LDA(pcashape2, ~Relative.stratigraphy, retain = 0.95) # 95% cum var PC scores
dashape95$CV.correct
dashape95$CV.ce

dashape99 <- LDA(pcashape2, ~Relative.stratigraphy, retain = 0.99) # 99% cum var PC scores
dashape99$CV.correct
dashape99$CV.ce

## Raw material
dashapefc <- LDA(efashape2, ~Raw.materials) # Fourier coefficient (raw data)
dashapefc$CV.correct
dashapefc$CV.ce

plot_LDA(dashapefc, axes = c(1,2) , zoom = 2, chull = FALSE) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range")

dashape95 <- LDA(pcashape2, ~Raw.materials, retain = 0.95) # 95% cum var PC scores
dashape95$CV.correct
dashape95$CV.ce

dashape99 <- LDA(pcashape2, ~Raw.materials, retain = 0.99) # 99% cum var PC scores
dashape99$CV.correct
dashape99$CV.ce

## MANOVA

## Member
efashape2 %>% MANOVA(~Member)
pcashape2 %>% MANOVA(~Member, retain = 0.95)
pcashape2 %>% MANOVA(~Member, retain = 0.99)

## Raw material
efashape2 %>% MANOVA(~Raw.materials)
pcashape2 %>% MANOVA(~Raw.materials, retain = 0.95)
pcashape2 %>% MANOVA(~Raw.materials, retain = 0.99)

## Reduction
efashape2 %>% MANOVA(~Knapping.modality)
pcashape2 %>% MANOVA(~Knapping.modality, retain = 0.95)
pcashape2 %>% MANOVA(~Raw.materials, retain = 0.99)



############## COMPARISON W/SIBUDU ############
# Full sample 
# PLot full sample
panel(outlinedata, main = "Border Cave and Sibudu", names = TRUE, fac = outlinedata$Site)

outlinedatasorted <- outlinedata %>% Momocs::arrange(Site, .by_group= TRUE)
#outlinedata <- Momocs::slice(outlinedata, -outliers) # to remove outliers

panel(outlinedatasorted, main = "Border Cave and Sibudu", names = TRUE, fac = outlinedatasorted$Site)
outlinedata <- Momocs::slice(outlinedata, -outliers) # to remove outlier 

# EFA
calibrate_harmonicpower_efourier(outlinedata, nb.h = 20, plot = FALSE) # 10 harmonics
calibrate_reconstructions_efourier(outlinedata, range = 1:20)
efashape <- efourier(outlinedata, nb.h = 10, norm = FALSE)

# PCA
pcashape <- PCA(efashape) 

# Size
centroidsize <- as_tibble(coo_centsize(outlinedata))
centroidsize <- rename(centroidsize, CS = "value")

pcascores <- as_tibble(pcashape$x)
databasedata <- cbind(outlinedata$fac, centroidsize, pcascores)

# Add in member/site variable
databasedata$SiteMember <- ifelse(databasedata$Member == "2BS", "2BS",
                                  ifelse(databasedata$Member == "2WA", "2WA",
                                         ifelse(databasedata$Member == "3BS", "3BS","Sibudu")))
newvar <- databasedata$SiteMember
newvar[is.na(newvar)] <- "Sibudu"

outlinedata$fac$SiteMember <- newvar 
efashape <- efourier(outlinedata, nb.h = 10, norm = FALSE)

pcashape <- PCA(efashape) 

# Plot 
scree_plot(pcashape, nax =1:10) # PC1-3 gives over 95% of cum variance in the data, PC1 = 70% of variance

gg <- PCcontrib(pcashape, nax = 1:4)
gg$gg + 
  geom_polygon(fill="slategrey", col="black") + 
  ggtitle("Shapes along PC1-4") # PC1 = asymmetry with base heavy, PC2 = wide to elongated, PC3 - wide to elongated with base heavy


pcascores <- as_tibble(pcashape$x)
databasedata <- cbind(outlinedata$fac, centroidsize, pcascores)

table(databasedata$Site) # priors

dashapefc <- LDA(efashape, ~Site, prior = c(54, 88)/nrow(databasedata), cv = TRUE) # Fourier coefficient (raw data)
dashapefc$CV.correct
dashapefc$CV.ce

plot.new()
par(mfrow = c(2,2))
plot_PCA(pcashape, axes = c(1,2), ~ Site, morphospace_position = "range_axes", zoom = 1, chull = FALSE, palette=pal_manual(c("indianred3",  "deepskyblue2")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1)  
plot_PCA(pcashape, axes = c(1,3), ~ Site, morphospace_position = "range_axes", zoom = 1, chull = FALSE,palette=pal_manual(c("indianred3", "deepskyblue2")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1)  
plot_PCA(pcashape, axes = c(2,3), ~ Site, morphospace_position = "range_axes", zoom = 1, chull = FALSE, palette=pal_manual(c("indianred3", "deepskyblue2")), eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1) 
vp.BottomRight <- viewport(height=unit(.5, "npc"), width=unit(0.5, "npc"), 
                           just=c("left","top"), 
                           y=0.5, x=0.5)
p1 <- boxplot(pcashape, ~Site, nax = 1:3, col = c("indianred3", "royalblue2"))
print(p1, vp = vp.BottomRight)  


plot_LDA(dashapefc)

dashapepc95 <- LDA(pcashape, ~Site, prior = c(54, 88)/nrow(databasedata), retain =0.95, cv = TRUE) # Fourier coefficient (raw data)
dashapepc95$CV.correct
dashapepc95$CV.ce

table(databasedata$SiteMember) # priors

dashapefc <- LDA(efashape, ~SiteMember, prior = c(36,16,2, 88)/nrow(databasedata), cv = TRUE) # Fourier coefficient (raw data)
dashapefc$CV.correct
dashapefc$CV.ce

plot.new()
plot_LDA(dashapefc , axes = c(1,2), zoom = 1, chull = TRUE, palette = col_summer2, morphospace_position = "range_axes", eigen = FALSE) %>% layer_points(cex = 1) %>% layer_morphospace_LDA(position = "range_axes")

dashapepc95 <- LDA(pcashape, ~SiteMember, prior = c(36,16,2, 88)/nrow(databasedata), retain =0.95, cv = TRUE) # Fourier coefficient (raw data)
dashapepc95$CV.correct
dashapepc95$CV.ce


efashape %>% MANOVA(~Site)
pcashape %>% MANOVA(~Site, retain = 0.95)
pcashape %>% MANOVA(~Site, retain = 0.99)

summary(aov(PC1~Site, databasedata))
summary(aov(PC2~Site, databasedata))
summary(aov(PC3~Site, databasedata))
summary(aov(PC4~Site, databasedata))

res <- aov(PC1~SiteMember, databasedata)
TukeyHSD(res)

res <- aov(PC2~SiteMember, databasedata)
TukeyHSD(res)

res <- aov(PC3~SiteMember, databasedata)
TukeyHSD(res)

twoBS<- subset(databasedata, SiteMember == "2BS")
twoWA <-  subset(databasedata, SiteMember == "2WA")
Sib <- subset(databasedata, SiteMember == "Sibudu")

t.test(twoBS$PC1, Sib$PC1)
t.test(twoWA$PC1, Sib$PC1)

t.test(twoBS$PC2, Sib$PC2)
t.test(twoWA$PC2, Sib$PC2)

t.test(twoBS$PC3, Sib$PC3)
t.test(twoWA$PC3, Sib$PC3)


## Seperate Sibudu 
Sibudu_pca <- Momocs:: filter(pcashape, Site=="Sibudu")
Sibudu_pca$fac <- droplevels(Sibudu_pca$fac)

## Relative stratigraphy
library(grid)
plot.new()
par(mfrow = c(2,2))
plot_PCA(Sibudu_pca, axes = c(1,2), ~ Relative.stratigraphy, morphospace_position = "range_axes", zoom = 0.5, chull = TRUE, eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1)  
plot_PCA(Sibudu_pca, axes = c(1,3), ~ Relative.stratigraphy, morphospace_position = "range_axes", zoom = 0.5, chull = TRUE, eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1)  
plot_PCA(Sibudu_pca, axes = c(2,3), ~ Relative.stratigraphy, morphospace_position = "range_axes", zoom = 0.5, chull = TRUE, eigen = FALSE, legend = FALSE) %>% layer_points(cex = 1) 
vp.BottomRight <- viewport(height=unit(.5, "npc"), width=unit(0.5, "npc"), 
                           just=c("left","top"), 
                           y=0.5, x=0.5)
p1 <- boxplot(Sibudu_pca, ~Relative.stratigraphy, nax = 1:3)
print(p1, vp = vp.BottomRight)  

databasedata2 <- subset(databasedata, Site == "Sibudu")

databasedata2$Relative.stratigraphy <- as.numeric(databasedata2$Relative.stratigraphy) # turn into numeric

ggplot(databasedata2, aes(PC1, Relative.stratigraphy)) + 
  geom_smooth(method = lm, col = "red", se = FALSE) + 
  geom_point(size = 2, pch = 16, alpha = 0.4) + 
  theme(text = element_text(size = 8), axis.text = element_text(size = 8)) + 
  ylab("Time")

cor.test(databasedata2$PC1, databasedata2$Relative.stratigraphy)
lmPC1 <- lm(PC1 ~Relative.stratigraphy, data = databasedata2)
summary(lmPC1)


ggplot(databasedata2, aes(PC2, Relative.stratigraphy)) + 
  geom_point(size = 2, pch = 16, alpha = 0.4) + 
  geom_smooth(method = lm, se = FALSE) + 
  theme(text = element_text(size = 8), axis.text = element_text(size = 8)) + 
  ylab("Time")

cor.test(databasedata2$PC2, databasedata2$Relative.stratigraphy)
lmPC2 <- lm(PC2 ~Relative.stratigraphy, data = databasedata2)
summary(lmPC2)

ggplot(databasedata2, aes(PC3, Relative.stratigraphy)) + 
  geom_point(size = 2, pch = 16, alpha = 0.4) + 
  geom_smooth(method = lm, se = FALSE) + 
  theme(text = element_text(size = 8), axis.text = element_text(size = 8)) + 
  ylab("Time")

cor.test(databasedata2$PC3, databasedata2$Relative.stratigraphy)
lmPC3 <- lm(PC3 ~Relative.stratigraphy, data = databasedata2)
summary(lmPC3)

############## SIBUDU ONLY ############

calibrate_harmonicpower_efourier(SC_outlines, nb.h = 20, plot = FALSE) # 11 harmonics
calibrate_reconstructions_efourier(SC_outlines, range = 1:20)

efashape <- efourier(SC_outlines, nb.h = 10, norm = FALSE)

