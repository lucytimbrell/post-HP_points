# Timbrell, L. et al. (2022). The post-Howiesons Poort points from Border Cave, KwaZulu Natal, South Africa. QSR. 
# This script performs the geometric morphometric analysis of the Sibudan morphotypes and Border Cave

# Data available to download at: https://github.com/lucytimbrell/post-HP_points

# This script follows protocols from Matzig et al (2021) - please note that the outline extraction is random and therefore to 
# fully reproduce the analysis, please use the "Normalised_outlinesSC.rds" file to ensure that the morphotype labels are 
# correctly assigned

# You also need to download the rds files "Normalised_outlinesBC.rds" and "BC_database.rds" - this is the Border Cave data.
 
############## TWO SITES TOGETHER ############

library(Momocs)
library(rio)
library(tidyverse)
library(outlineR)

############## Load in Border Cave ############## 
setwd("...")

BC_outlines <- import("Normalised_outlinesBC.rds")
database <- import("BC_database.rds")  

panel(BC_outlines, main = "Border Cave", names = TRUE)
stack(BC_outlines)

BC_efashape <- efourier(BC_outlines, nb.h = 11, norm = FALSE)
BC_pcashape <- PCA(BC_efashape) 

############## Get Sibudu outlines ############## 
setwd("/Volumes/Lucy HD/BorderCave_Sibudu/project_folder2")

### set paths
inpath <- file.path(".", "prep")
outpath <- file.path(".", "prep", "single")
tpspath <- file.path(".", "raw","WC_SC_59_TPS.TPS")

### seperate images
outlineR::separate_single_artefacts(inpath = inpath, 
                                    outpath = outpath)
### import scaling factor
tps_df <-outlineR::tps_to_df(tpspath)

### extract outlines according to scaling factor
single_outlines_list <- outlineR::get_outlines(outpath = outpath, 
                                               tps_file_rescale = tps_df) ### need to remove number 69

### combine outlines into single Out object
SC_outlines_combined <- outlineR::combine_outlines(single_outlines_list = single_outlines_list)

length(SC_outlines_combined) #how many outlines do you have?

Momocs::panel(SC_outlines_combined) # shows all outlines next to each other - note when combine does not do it in same order as original image

Morphotype <- c(rep("Tongati", 19), rep("Ndwedwe", 2), "Tongati", rep("Ndwedwe", 5), rep("ACT", 5), "Tongati", 
                    rep("ACT", 5), "Tongati", rep("Unretouched", 5), "Tongati", rep("Unretouched", 10), rep("Tongati", 4))


SC_outlines  <- Out(SC_outlines_combined$coo, fac = Morphotype)

### normalisation
SC_outlines <- SC_outlines %>% 
  coo_centre() %>% 
  coo_scale() %>% 
  coo_slidedirection("right") %>% 
  coo_close() 

colnames(SC_outlines$fac) <- "Morphotype"
SC_outlines$fac$Morphotype <- as.factor(SC_outlines$fac$Morphotype)
panel(SC_outlines, main = "Sibudu Cave outlines", fac = "Morphotype", palette =pal_manual(c("#888888","#999933","#6699CC","#44AA99" ))) # Visualization of points in their original orientation

setwd("/Volumes/Lucy HD/BorderCave_Sibudu/github_files")
saveRDS(SC_outlines, file = "Normalised_outlinesSC.rds") # save normalised and transformed landmarks

######  ANALYSIS ###### 
#### make sure to load in this specific rds file as the outline generation is random and therefore the morphotype assignment will be different
setwd("/Volumes/Lucy HD/BorderCave_Sibudu/github_files") 
SC_outlines <- import("Normalised_outlinesSC.rds")

SC_efashape <- efourier(SC_outlines, nb.h = 11, norm = FALSE)
SC_pcashape <- PCA(SC_efashape) 

#### Establish morphotypes ####

plot.new()
gg <- PCcontrib(BC_pcashape, nax = 1:3, plot = FALSE)
gg$gg + 
  geom_polygon(fill="white", col="black") 

tidy.table.SC <- cbind(as.tibble(SC_pcashape$x[,1:3]))
tidy.table.SC$Morphotype <- Morphotype # create new factor for morphotypes 

a <- ggplot(tidy.table.SC, aes(PC1, PC2, colour = Morphotype)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(limits = c("Unretouched", "Tongati", "Ndwedwe", "ACT"), values = c("#44AA99", "#6699CC", "#999933", "#888888")) +
  theme_minimal() +
  theme(legend.position = "none")


b <- ggplot(tidy.table.SC, aes(PC1, PC3, colour = Morphotype)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(limits = c("Unretouched", "Tongati", "Ndwedwe", "ACT"), values = c("#44AA99", "#6699CC", "#999933", "#888888")) +
  theme_minimal() +
  theme(legend.position = "none")

c <- ggplot(tidy.table.SC, aes(PC2, PC3, colour = Morphotype)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(limits = c("Unretouched", "Tongati", "Ndwedwe", "ACT"), values = c("#44AA99", "#6699CC", "#999933", "#888888")) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        axis.title.x = element_blank())


d <- tidy.table.SC %>% 
  dplyr::select(PC1, PC2, PC3, Morphotype) %>%
  pivot_longer(cols = contains("PC")) %>%
  mutate(name = as.factor(name)) %>%
  mutate(Morphotype = as.factor(Morphotype)) %>%
  mutate(Morphotype = fct_relevel(Morphotype, "Unretouched", "Tongati", "Ndwedwe", "ACT")) %>%
  ggplot(aes(name, value, colour = Morphotype)) +
  scale_colour_manual(limits =  c("Unretouched", "Tongati", "Ndwedwe", "ACT"),values = c("#44AA99", "#6699CC", "#999933", "#888888")) +
  labs(y = "Score") +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title.x = element_blank())

(a | b | c) / d +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 30))) 

### Analysis 
table(tidy.table.SC$Morphotype)

fit <- lda(Morphotype~., data = tidy.table.SC, prior = c(10,7,27,15)/nrow(tidy.table.SC), CV = TRUE) # 95% cum var PC scores

ct <- table(tidy.table.SC$Morphotype, fit$class) # Ndwedwe + Tongati are very distinguishable
diag(prop.table(ct, 1))
sum(diag(prop.table(ct)))

summary(manova(cbind(PC1, PC2, PC3) ~Morphotype, tidy.table.SC))

summary(aov(PC1~Morphotype, tidy.table.SC))
summary(aov(PC2~Morphotype, tidy.table.SC))
summary(aov(PC3~Morphotype, tidy.table.SC))

res <- aov(PC1~Morphotype, tidy.table.SC) ## extra Tukey HSD as is significant 
TukeyHSD(res)

res <- aov(PC2~Morphotype, tidy.table.SC) ## extra Tukey HSD as is significant 
TukeyHSD(res)

res <- aov(PC3~Morphotype, tidy.table.SC) ## extra Tukey HSD as is significant 
TukeyHSD(res)

######  predict PC scores from BC PCA ###### 
newPCA <- rePCA(BC_pcashape, SC_efashape)

tidy.table <- cbind(as.tibble(BC_pcashape$x[,1:3]))

### add tibbles together
# look at order of file names of SC_outlines and use that to create a vector of names
Morphotype <- c(rep("Tongati", 27), rep("Ndwedwe", 7), rep("ACT", 10), rep("Unretouched", 15))

tidy.table2 <- as.tibble(newPCA$x[,1:3]) # pc scores for morphotypes
tidy.table3 <- rbind(tidy.table2[,1:3], tidy.table[,1:3]) # join with pc scores for BC
tidy.table3$Morphotype <- c(Morphotype, rep("Border Cave", 54)) # create new factor for morphotypes vs BC

x <- ggplot(tidy.table3, aes(PC1, PC2, colour = Morphotype)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(limits = c("Border Cave", "Unretouched", "Tongati", "Ndwedwe", "ACT"), values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  theme_minimal() +
  theme(legend.position = "none")

y <- ggplot(tidy.table3, aes(PC1, PC3, colour = Morphotype)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(limits = c("Border Cave", "Unretouched", "Tongati", "Ndwedwe", "ACT"), values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  theme_minimal() +
  theme(legend.position = "none")

z <- ggplot(tidy.table3, aes(PC2, PC3, colour = Morphotype)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(limits = c("Border Cave", "Unretouched", "Tongati", "Ndwedwe", "ACT"), values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  theme_minimal() +
  theme(legend.position = "none")

d <- tidy.table3 %>% 
  dplyr::select(PC1, PC2, PC3, Morphotype) %>%
  pivot_longer(cols = contains("PC")) %>%
  mutate(name = as.factor(name)) %>%
  mutate(Morphotype = as.factor(Morphotype)) %>%
  mutate(Morphotype = fct_relevel(Morphotype, "Border Cave", "Unretouched", "Tongati", "Ndwedwe", "ACT")) %>%
  ggplot(aes(name, value, colour = Morphotype)) +
  scale_colour_manual(limits = c("Border Cave", "Unretouched", "Tongati", "Ndwedwe", "ACT"), values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  labs(y = "Score") +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())

(x | y | z) / d +
  plot_annotation(title = "'Sibudan' points on Border Cave PCs", tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 24))) 

setwd("/Volumes/Lucy HD/BorderCave_Sibudu")
write.csv(as.data.frame(tidy.table3), "morphotype_pcscores.csv")

### Analysis 
library(MASS)
table(tidy.table3$Morphotype)

fit <- lda(Morphotype~., data = tidy.table3, prior = c(10,54,7,27,15)/nrow(tidy.table3), CV = TRUE) # 95% cum var PC scores

ct <- table(tidy.table3$Morphotype, fit$class)
diag(prop.table(ct, 1))
sum(diag(prop.table(ct)))

summary(manova(cbind(PC1, PC2, PC3) ~Morphotype, tidy.table3))

summary(aov(PC1~Morphotype, tidy.table3))
summary(aov(PC2~Morphotype, tidy.table3))
summary(aov(PC3~Morphotype, tidy.table3))

res <- aov(PC1~Morphotype, tidy.table3) ## extra Tukey HSD as is significant 
TukeyHSD(res)

res <- aov(PC2~Morphotype, tidy.table3) ## extra Tukey HSD as is significant 
TukeyHSD(res)

res <- aov(PC3~Morphotype, tidy.table3) ## extra Tukey HSD as is significant 
TukeyHSD(res)

#### Reverse projection ####

newPCA2 <- rePCA(SC_pcashape, BC_efashape)

tidy.table2 <- as.tibble(newPCA2$x[,1:3]) # pc scores for BC
tidy.table3 <- rbind(tidy.table.SC[,1:3], tidy.table2[,1:3]) # join with pc scores for morphotypes
tidy.table3$Morphotype <- c(Morphotype, rep("Border Cave", 54)) # create new factor for morphotypes vs BC

l <- ggplot(tidy.table3, aes(PC1, PC2, colour = Morphotype)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(limits = c("Border Cave", "Unretouched", "Tongati", "Ndwedwe", "ACT"), values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  theme_minimal() +
  theme(legend.position = "none")


m <- ggplot(tidy.table3, aes(PC1, PC3, colour = Morphotype)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(limits = c("Border Cave", "Unretouched", "Tongati", "Ndwedwe", "ACT"), values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  theme_minimal() +
  theme(legend.position = "bottom")

n <- ggplot(tidy.table3, aes(PC2, PC3, colour = Morphotype)) + 
  geom_point(size = 3) +
  stat_ellipse() +
  scale_colour_manual(limits = c("Border Cave", "Unretouched", "Tongati", "Ndwedwe", "ACT"), values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  theme_minimal() +
  theme(legend.position = "none")

d <- tidy.table3 %>% 
  dplyr::select(PC1, PC2, PC3, Morphotype) %>%
  pivot_longer(cols = contains("PC")) %>%
  mutate(name = as.factor(name)) %>%
  mutate(Morphotype = as.factor(Morphotype)) %>%
  mutate(Morphotype = fct_relevel(Morphotype, "Border Cave", "Unretouched", "Tongati", "Ndwedwe", "ACT")) %>%
  ggplot(aes(name, value, colour = Morphotype)) +
  scale_colour_manual(limits =  c("Border Cave", "Unretouched", "Tongati", "Ndwedwe", "ACT"),values = c("#AA4499", "#44AA99", "#6699CC", "#999933", "#888888")) +
  labs(y = "Score") +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())

(l | m | n) / d +
  plot_annotation(title = "Border Cave points on 'Sibudan' PCs", tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 24))) 

### Analysis 
table(tidy.table3$Morphotype)

fit <- lda(Morphotype~., data = tidy.table3, prior = c(10,54,7,27,15)/nrow(tidy.table3), CV = TRUE) # 95% cum var PC scores

ct <- table(tidy.table3$Morphotype, fit$class)
diag(prop.table(ct, 1))
sum(diag(prop.table(ct)))

summary(manova(cbind(PC1, PC2, PC3) ~Morphotype, tidy.table3))

summary(aov(PC1~Morphotype, tidy.table3))
summary(aov(PC2~Morphotype, tidy.table3))
summary(aov(PC3~Morphotype, tidy.table3))

res <- aov(PC1~Morphotype, tidy.table3) ## extra Tukey HSD as is significant 
TukeyHSD(res)

res <- aov(PC2~Morphotype, tidy.table3) ## extra Tukey HSD as is significant 
TukeyHSD(res)

res <- aov(PC3~Morphotype, tidy.table3) ## extra Tukey HSD as is significant 
TukeyHSD(res)

############## PLOT ############## 
(x | y | z) / (l | m | n)  +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 30))) 

