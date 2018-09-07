

##############################################################################################
######################## SET WORKING DIRECTORY ###############################################
##############################################################################################

# MAC
setwd("~/Desktop/MetaAnalysis")

##############################################################################################
######################## INSTALL PACKAGES & LIBRARIES ########################################
##############################################################################################

install.packages("MASS")
library(MASS)
install.packages("meta")
library(meta)
install.packages("metafor")
library(metafor)
install.packages("dplyr")
library(dplyr)
install.packages("tidyr")
library(tidyr)
install.packages("ggfortify")
library(ggfortify)
install.packages("ggplot2")
library(ggplot2)
# package to for set function, used to replace NA's in metabolic rate data with 0's
install.packages("sets")
library(sets)
# to combine multiple graphs in one frame
library(gridExtra)
library(grid)
## mixed models
install.packages("lme4")
library(lme4)
install.packages("lmerTest")
library(lmerTest)


##############################################################################################
######################## IMPORT DATA & ASSIGN CLASSES ########################################
##############################################################################################

##import raw data
raw.data <- read.csv("madata_7_10_18_2.csv", header = TRUE, stringsAsFactors = TRUE)

## remove studies that we are not using right now (no effect size calculated)
data <- filter(raw.data, use == "yes") 

## assign classess to variables

data$study <- as.character(data$study)
data$study <- as.factor(data$study)

data$temp <- as.character(data$temp)
data$temp <- as.numeric(data$temp)

data$severity.score <- as.character(data$severity.score)
data$severity.score <- as.factor(data$severity.score)

data$longev <- as.character(data$longev)
data$longev <- as.numeric(data$longev)

data$log.longev <- as.character(data$log.longev)
data$log.longev <- as.numeric(data$log.longev)

data$f.age.mat <- as.character(data$f.age.mat)
data$f.age.mat <- as.numeric(data$f.age.mat)

data$c.size <- as.character(data$c.size)
data$c.size <- as.numeric(data$c.size)

data$c.year <- as.character(data$c.year)
data$c.year <- as.numeric(data$c.year)

data$latitude <- as.character(data$latitude)
data$latitude <- as.numeric(data$latitude)

##############################################################################################
############# LIFE HISTORY DATA PCA (MANUSCRIPT FIGURE) #######################################
##############################################################################################

## create a life history data set that only includes records with values for all 4 lh variables

#lh.data <- filter(data, log.longev > -100 & log.f.age.mat > -100 & log.c.size > -100 & log.c.year > -100)

lh.data <- filter(data, longev > 0 & f.age.mat > 0 & c.size > 0 & c.year > 0)

# I can't get select function to work so I'm using longer code below to select columns
#lh.data <- select(lh.data, class, order, species, mass, longev, f.age.mat, c.size, c.year)
#lh.data %>%
  #select(class)

#lh.data <- data.frame(lh.data$class, lh.data$order, lh.data$species, lh.data$log.mass, 
                   #lh.data$log.longev, lh.data$log.f.age.mat, lh.data$log.c.size, lh.data$log.c.year)

lh.data <- data.frame(lh.data$class, lh.data$order, lh.data$species, 
                      lh.data$longev, lh.data$f.age.mat, lh.data$c.size, lh.data$c.year)

## Generate headers without the numbers - we will use below when we bind all together
list.heads <- c("class", "order", "species", "longev", "fam", "c.size", "c.year")

names(lh.data) <- list.heads


## manual scaling of data
lh.data.m.sc <-  mutate(lh.data, s.longev = (longev - mean(longev)) / sd(longev), 
                   s.fam = (fam - mean(fam)) / sd(fam), 
                   s.c.size = (c.size - mean(c.size))/sd(c.size), 
                   s.c.year = (c.year - mean(c.year))/sd(c.year))


## Remove species duplicates
lh.data.m.sc <- unique(lh.data.m.sc[,1:11])

## we have 38 unique species for which we have all 4 lh variables to conduct pca

## Rounding each variable to 2 decimal points and dropping any rows that contain missing values (NA)
#lh.data %>% mutate(longev = round(longev,2), fam = round(fam,2), c.size = round(c.size,2), c.year = round(c.year,2))

## Removing columns that are not LH variables so that we can run PCA on this data frame
lh.data.m.sc.clean <- lh.data.m.sc[,8:11]



## auto scaling of data with r scaling function

## Remove species duplicates
#lh.data.a.sc <- unique(lh.data[,1:7])

## we have 37 unique species for which we have all 4 lh variables to conduct pca

## Removing columns that are not LH variables so that we can run PCA on this data frame
#lh.data.a.sc.clean <- lh.data.a.sc[,4:7]

#test <- scale(lh.data.a.sc.clean, center = TRUE, scale = TRUE)

#head(lh.data.m.sc.clean)

pca.manual <- prcomp(lh.data.m.sc.clean)

summary(pca.manual)

pca.data <- cbind(lh.data.m.sc, pca.manual$x[,1:2])

ggplot(pca.data, aes(x = PC1, y = PC2))+
  geom_point(aes(col = class))+
  theme_classic()


## plotting pc1 by pc2
pca.graph <- autoplot(prcomp(lh.data.m.sc.clean), data = lh.data.m.sc, 
         loadings = TRUE, loadings.colour = "black", loadings.label = TRUE,
         loadings.label.colour = "black", 
         loadings.label.hjust = -0.3,loadings.label.vjust = 0.08)+
  geom_point(aes(col = class), size = 4)+
  scale_color_manual(values = c("black", "grey50"))+
  theme_classic()+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))+
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 12))+
  guides(col=guide_legend(title="Taxonomic Class"))

ggsave("PCA_plot.jpg", plot = pca.graph, dpi = 500, height = 5, width = 6, units = "in")




## Let's see for how many species, we have incomplete lh info
## This tells us how many species are missing from pca b/c of incomplete info

all.lh.data <- filter(data, longev > 0 | f.age.mat > 0 | c.size > 0 | c.year > 0)

all.lh.data <- data.frame(all.lh.data$class, all.lh.data$order, all.lh.data$species, all.lh.data$mass, 
                      all.lh.data$longev, all.lh.data$f.age.mat, all.lh.data$c.size, all.lh.data$c.year)

## Generate headers without the numbers - we will use below when we bind all together
list.heads <- c("class", "order", "species", "mass", "longev", "fam", "c.size", "c.year")

names(all.lh.data) <- list.heads

## Remove species duplicates
all.lh.data <- unique(all.lh.data[,1:8])

## Rounding each variable to 2 decimal points and dropping any rows that contain missing values (NA)
lh.data %>% mutate(longev = round(longev,2), fam = round(fam,2), c.size = round(c.size,2), c.year = round(c.year,2))

## we have 51 species when we include incomplete records of lh data. 
# so by using only those species with all 4 variables, we exclude only 14 species from the pca.


##############################################################################################
################### MASS X METABOLISM (YES) (MANUSCRIPT FIGURE) ##############################
##############################################################################################



### Mass / Metabolism Graph

## create dataset with only rows that have metabolism
met.data <- filter(data, mr > 0)

## select only certain columns so that we can then select unique combos of met and mass
## we have to specifiy to use select from dplyr, otherwise, r tries to use it from the library MASS and 
# we get an error
met.data <- dplyr::select(met.data, class, order, species, mass, log.mass, mr, log.mr, ms.mr, log.ms.mr)

# headers without the numbers - we will use below when we bind all together
list.heads <- c("Class", "order", "species", "mass", "log.mass", "mr", "log.mr", "ms.mr", "log.ms.mr")

names(met.data) <- list.heads

met.data.unique <- unique(met.data[,1:9])

#### graph of mass x whole-body metabolic rate

wb.mr.graph <- ggplot(met.data.unique, aes(x = log.mass, y = log.mr, col = Class, shape = Class))+
  geom_point(size = 3)+
  theme_classic()+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))+
  scale_shape_discrete(solid = F)+
  scale_colour_manual(values = c("grey20", "black", "grey40"))+
  geom_smooth(method = "lm", col = "black", alpha = 0.3, size = 0.5)+
  ylab("Log Metabolic Rate")+
  xlab("Log Mass")+
  theme(legend.position=c(.7,0.3), legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

ggsave(wb.mr.graph, file = "mass.mr.jpg", width = 5, height = 5, dpi = 600, units = "in")

## graph of mass x mass-specific metabolic rate

#ms.mr.graph <- ggplot(met.data.unique, aes(x = log.mass, y = log.ms.mr, col = Class, shape = Class))+
#geom_point(size = 3)+
# theme_classic()+
# theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))+
# scale_shape_discrete(solid = F)+
# scale_colour_manual(values = c("grey20", "black", "grey40"))+
# geom_smooth(method = "lm", col = "black", alpha = 0.3, size = 0.5)+
# ylab("Log Mass-Specific MR")+
# xlab("Log Species Mass")+
# theme(legend.position=c(.8,0.85), legend.text = element_text(size = 12),
#       legend.title = element_text(size = 14))+
# theme(legend.position = "none")


#mass.mr <- grid.arrange(wb.mr.graph, ms.mr.graph, ncol = 1)

#ggsave(mass.mr, file = "mass.mr.jpg", width = 5, height = 10, units = "in", dpi = 600)


## models assessing relationship between mass and mr across taxa

## whole body mr ~ mass for all animals
mass.mr <- lm(log.mr ~ log.mass, data = met.data.unique)
summary(mass.mr)
anova(mass.mr)

## whole body mr ~ mass * taxa
mass.mr.taxa <- lm(log.mr ~ log.mass * Class, data = met.data.unique)
summary(mass.mr.taxa)

## whole body mr ~ mass * taxa
mass.mr.taxa <- lm(log.mr ~ log.mass + Class, data = met.data.unique)
summary(mass.mr.taxa)



## whole body mr ~ mass for birds
mass.mr.birds <- lm(log.mr ~ log.mass, data = subset(met.data.unique, Class == "Aves"))
summary(mass.mr.birds)


## mass-specific mr ~ mass for fish
mass.ms.mr.fish <- lm(log.mr ~ log.mass, data = subset(met.data.unique, Class == "Actinopterygii"))
summary(mass.ms.mr.fish)


## whole body mr ~ mass for mammals
mass.mr.mamm <- lm(log.mr ~ log.mass, data = subset(met.data.unique, Class == "Mammalia"))
summary(mass.mr.mamm)



## mass-specific mr ~ mass for mammals
mass.ms.mr.mamm <- lm(log.ms.mr ~ log.mass, data = subset(met.data.unique, Class = "Mammalia"))
summary(mass.ms.mr.mamm)


## whole body mr ~ mass for fish
mass.mr.fish <- lm(log.mr ~ log.mass, data = subset(met.data.unique, Class = "Actinopterygii"))
summary(mass.mr.fish)

## mass-specific mr ~ mass for fish
mass.ms.mr.fish <- lm(log.ms.mr ~ log.mass, data = subset(met.data.unique, Class = "Actinopterygii"))
summary(mass.ms.mr.fish)



##############################################################################################
################ ANIMAL TYPE (NO) (MANUSCRIPT FIGURE) ########################################
##############################################################################################

at.es <- ggplot(data, aes(x = animal.type, y = hedges.d))+
  geom_jitter(width = 0.15, alpha = 0.15, aes(size=1/variance))+
  theme_classic()+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 20))+
  #geom_abline(intercept = -1.6, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = -0.7, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = -0.22, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = 0.26, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #theme(legend.position = "none")+
  ylab("Effect Size")+
  xlab("Animal Taxanomic Group")+
  geom_point(data = at.mod1.results, aes(x = animal.type, y = hedges.d), pch = 15, size=4, col="black")+
  geom_errorbar(data=at.mod1.results, aes(x = animal.type, ymin = hedges.d-se, ymax = hedges.d+se), 
                col= "black", width = 0.15)+
  theme(legend.position = "none")
  #theme(legend.position = "right", legend.text = element_text(size = 12),
        #legend.title = element_text(size = 14))
  
ggsave(at.es, file = "taxa.hd.jpg", width = 5, height = 5, units = "in", dpi = 600)

data$animal.type <- relevel(data$animal.type, ref = "Mammal")


at.mod1 <- rma.mv(hedges.d, variance, mods = ~factor(animal.type),
                  random = list(~1|study),
                  W = 1/variance, method = "REML", data = data)
summary(at.mod1)

animal.type <- c("Bird", "Fish", "Herp", "Invert", "Mammal", "Shellfish")
hedges.d <- c(0.40, 0.85, 0.85, 0.64, 0.60, 0.67)
se <- c(0.46, 0.64, 0.59, 0.59, 0.54, 0.71)

at.mod1.results <- cbind.data.frame(animal.type, hedges.d, se)

as.data.frame(at.mod1.results, header = TRUE)


##############################################################################################
########################## HOMEOTHERMS VS POIKILOTHERMS (NO) #################################
##############################################################################################


### Homeotherm vs. Poikilotherm

ggplot(data, aes(x = h.p, y = hedges.d))+
  geom_point(aes(size = 1/variance))

mod.hp <- rma.mv(hedges.d, variance, mods = ~factor(h.p),
                 random = ~1|study,
                 W = 1/variance, method = "REML", data = data)

summary(mod.hp)


##############################################################################################
############################# RODENTS VS. OTHER MAMMALS (NO) #################################
##############################################################################################

mod.rodents1 <- rma.mv(hedges.d, ln.variance, mods = ~factor(Rodent),
                       random = ~1|study,
                       W = 1/ln.variance, method = "REML", data = subset(data, class = "Mammalia"))

summary(mod.rodents1)

ggplot(data, aes(x = Rodent, y = hedges.d))+
  geom_point(aes(size = 1/variance))


##############################################################################################
################################ FIELD VS. LAB (NO) ##########################################
##############################################################################################

ggplot(data, aes(x=field.lab, y=hedges.d))+
  geom_jitter(width = 0.2, alpha = 0.3, aes(size=(1/variance), col = animal.type))+
  theme_classic()+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))+
  #geom_abline(intercept = -1.6, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = -0.7, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = -0.22, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = 0.26, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  geom_point(data = fl.mod1.results, aes(x = field.lab, y = hedges.d), size = 4, col = "red", pch = 15)+
  geom_errorbar(data=fl.mod1.results, aes(x = field.lab, ymin = hedges.d-se, ymax = hedges.d+se), 
                col= "red", width = 0.1)+
  ylab("Hedges d")+
  xlab("Experimental Design")+
  scale_x_discrete(labels = c("Field", "Field Enclosure", "Lab"))+
  theme(legend.position = "right", legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))


data$field.lab <- relevel(data$field.lab, ref = "field")


fl.mod1 <- rma.mv(hedges.d, variance, mods = ~factor(field.lab), 
                  random = list(~1|study),
                  W = 1/variance, method = "REML", data = data)
summary(fl.mod1)

# generate dataframe with results for fl.mod1
field.lab <- c("field", "field.enclosure", "lab")
hedges.d <- c(0.50, 0.68, 0.89)
se <- c(0.28, 0.44, 0.48)

fl.mod1.results <- cbind.data.frame(field.lab, hedges.d, se)

as.data.frame(fl.mod1.results, header = TRUE)


## let's look to see if a given taxa responded to field.lab conditions

## birds - no
ggplot(data = subset(data, animal.type == "Bird"), aes(x=field.lab, y=hedges.d))+
  geom_jitter(width = 0.2, alpha = 0.3, aes(size=(1/variance)))+
  theme_classic()

fl.mod.birds <- rma.mv(hedges.d, variance, mods = ~factor(field.lab), 
                  random = list(~1|study),
                  W = 1/variance, method = "REML", data = subset(data, animal.type == "Bird"))
summary(fl.mod.birds)

## fish - no
ggplot(data = subset(data, animal.type == "Fish"), aes(x=field.lab, y=hedges.d))+
  geom_jitter(width = 0.2, alpha = 0.3, aes(size=(1/variance)))+
  theme_classic()

fl.mod.fish <- rma.mv(hedges.d, variance, mods = ~factor(field.lab), 
                       random = list(~1|study),
                       W = 1/variance, method = "REML", data = subset(data, animal.type == "Fish"))
summary(fl.mod.fish)


## herp - no
ggplot(data = subset(data, animal.type == "Herp"), aes(x=field.lab, y=hedges.d))+
  geom_jitter(width = 0.2, alpha = 0.3, aes(size=(1/variance)))+
  theme_classic()

fl.mod.herp <- rma.mv(hedges.d, variance, mods = ~factor(field.lab), 
                      random = list(~1|study),
                      W = 1/variance, method = "REML", data = subset(data, animal.type == "Herp"))
summary(fl.mod.herp)


data$field.lab <- relevel(data$field.lab, ref = "field")

## inverts - FIELD ENCLOSURE HAS LOWER EFFECT SIZE THAN FIELD AND LAB
# this doesn't seem ecologically relevant or logical - also small ss of 
#field enclosure studies with inverts
ggplot(data = subset(data, animal.type == "Invert"), aes(x=field.lab, y=hedges.d))+
  geom_jitter(width = 0.2, alpha = 0.3, aes(size=(1/variance), col = study))+
  theme_classic()

fl.mod.invert <- rma.mv(hedges.d, variance, mods = ~factor(field.lab), 
                      random = list(~1|study),
                      W = 1/variance, method = "REML", data = subset(data, animal.type == "Invert"))
summary(fl.mod.invert)

# mammals - no
ggplot(data = subset(data, animal.type == "Mammal"), aes(x=field.lab, y=hedges.d))+
  geom_jitter(width = 0.2, alpha = 0.3, aes(size=(1/variance)))+
  theme_classic()

fl.mod.mamm <- rma.mv(hedges.d, variance, mods = ~factor(field.lab), 
                        random = list(~1|study),
                        W = 1/variance, method = "REML", data = subset(data, animal.type == "Mammal"))
summary(fl.mod.mamm)

# shellfish - no
ggplot(data = subset(data, animal.type == "Shellfish"), aes(x=field.lab, y=hedges.d))+
  geom_jitter(width = 0.2, alpha = 0.3, aes(size=(1/variance)))+
  theme_classic()

fl.mod.sf <- rma.mv(hedges.d, variance, mods = ~factor(field.lab), 
                      random = list(~1|study),
                      W = 1/variance, method = "REML", data = subset(data, animal.type == "Shellfish"))
summary(fl.mod.sf)



##############################################################################################
################################ RESPONSE TYPE (NO) ##########################################
##############################################################################################

ggplot(data, aes (x = response.type, y = hedges.d))+
  geom_jitter(width = 0.15, aes(size=1/variance, col = animal.type), alpha = 0.5)+
  theme_classic()+
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16))+
  #geom_abline(intercept = -1.6, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = -0.7, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = -0.22, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = 0.26, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  geom_point(data = rt.mod1.results, aes(x = response.type, y=hedges.d), size = 4, col = "red", pch=15)+
  geom_errorbar(data=rt.mod1.results, aes(x = response.type, ymin = hedges.d-se, ymax = hedges.d+se), 
               col= "red", width = 0.15)+
  ylab("Hedges d")+
  xlab("Measured Foraging Response to Risk")
  #scale_x_discrete(labels = c("Intake", "Time", "Visits"))+
theme(legend.position = "right", legend.text = element_text(size = 12),
      legend.title = element_text(size = 14))

data$response.type <- relevel(data$response.type, ref = "time")

# model assessing differences in effect size among response types
rt.mod1 <- rma.mv(hedges.d, variance, mods = ~factor(response.type),
                  random = list(~1|study),
                  W = 1/variance, method = "REML", data = data)
summary(rt.mod1)

# generating a dataset with means, se, upper, lower from rt.mod.1

response.type <- c("intake", "other", "time", "visits")
hedges.d <- c(0.68, 0.81, 0.44, 0.43)
se <- c(0.18, 0.42, 0.36, 0.45)


rt.mod1.results <- cbind.data.frame(response.type, hedges.d, se)

as.data.frame(rt.mod1.results, header = TRUE)


## let's look at reponse type variation within taxa

## birds - no
ggplot(data = subset(data, animal.type == "Bird"), aes (x = response.type, y = hedges.d))+
  geom_jitter(width = 0.15, aes(size=1/variance), alpha = 0.5)+
  theme_classic()

## fish - no
ggplot(data = subset(data, animal.type == "Fish"), aes (x = response.type, y = hedges.d))+
  geom_jitter(width = 0.15, aes(size=1/variance), alpha = 0.5)+
  theme_classic()

## herp - no
ggplot(data = subset(data, animal.type == "Herp"), aes (x = response.type, y = hedges.d))+
  geom_jitter(width = 0.15, aes(size=1/variance, col = study), alpha = 0.5)+
  theme_classic()

## invert - no
ggplot(data = subset(data, animal.type == "Invert"), aes (x = response.type, y = hedges.d))+
  geom_jitter(width = 0.15, aes(size=1/variance, col = study), alpha = 0.5)+
  theme_classic()

## mammals - no
ggplot(data = subset(data, animal.type == "Mammal"), aes (x = response.type, y = hedges.d))+
  geom_jitter(width = 0.15, aes(size=1/variance), alpha = 0.5)+
  theme_classic()

## shellfish - no
ggplot(data = subset(data, animal.type == "Shellfish"), aes (x = response.type, y = hedges.d))+
  geom_jitter(width = 0.15, aes(size=1/variance, col = study), alpha = 0.5)+
  theme_classic()


##############################################################################################
################################ SEVERITY SCORE (NO) #########################################
##############################################################################################

ggplot(data, aes (x = severity.score, y = hedges.d))+
  geom_jitter(width = 0.15, aes(size=1/variance), alpha = 0.5)+
  theme_classic()+
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16))+
  #geom_abline(intercept = -1.6, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = -0.7, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = -0.22, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = 0.26, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  geom_point(data = ss.mod.results, aes(x = severity.score, y=hedges.d), size = 4, col = "red", pch=15)+
  geom_errorbar(data=ss.mod.results, aes(x = severity.score, ymin = hedges.d-se, ymax = hedges.d+se), 
  col= "red", width = 0.15)+
  ylab("Hedges d")+
  xlab("Severity Score")
  scale_x_discrete(labels = c("0", "1", "2"))+
  theme(legend.position = "right", legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))


data$severity.score <- relevel(data$severity.score, ref = "1")
  
## model
ss.mod <- rma.mv(hedges.d, variance, mods = ~factor(severity.score),
                  random = list(~1|study),
                  W = 1/variance, method = "REML", data = data)
summary(ss.mod)


### generate data to plot
severity.score <- c("0", "1", "2")
hedges.d <- c(0.56, 0.76, 0.74)
se <- c(0.27, 0.43, 0.54)


ss.mod.results <- cbind.data.frame(severity.score, hedges.d, se)

as.data.frame(ss.mod.results, header = TRUE)

## just look at whether the animal was food deprived - nope

data$food.restrict <- as.character(data$food.restrict)
data$food.restrict <- as.factor(data$food.restrict)

ggplot(data, aes (x = food.restrict, y = hedges.d))+
  geom_jitter(width = 0.15, aes(size=1/variance, col = animal.type), alpha = 0.5)+
  theme_classic()

fr.mod <- rma.mv(hedges.d, variance, mods = ~factor(food.restrict),
                 random = list(~1|study),
                 W = 1/variance, method = "REML", data = data)
summary(fr.mod)

## just look at whether the animal had a safe food option - nope

data$safe.food.option <- as.character(data$safe.food.option)
data$safe.food.option <- as.factor(data$safe.food.option)

ggplot(data, aes (x = safe.food.option, y = hedges.d))+
  geom_jitter(width = 0.15, aes(size=1/variance, col = animal.type), alpha = 0.5)+
  theme_classic()

so.mod <- rma.mv(hedges.d, variance, mods = ~factor(safe.food.option),
                 random = list(~1|study),
                 W = 1/variance, method = "REML", data = data)
summary(so.mod)

##############################################################################################
############################## RISK CUE (NO) #################################################
##############################################################################################

#data.one$risk.cue <- relevel(data.one$risk.cue, ref = "non.lethal")
#data.trunc$risk.cue <- relevel(data.trunc$risk.cue, ref = "enviro")

ggplot(data, aes (x = risk.cue, y = hedges.d))+
  geom_jitter(width = 0.15, aes(size=1/variance), alpha = 0.5)+
  theme_classic()+
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16))+
  #geom_abline(intercept = -1.6, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = -0.7, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = -0.22, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = 0.26, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  theme(legend.position = "none")+
  geom_point(data = rc.mod.results, aes(x = risk.cue, y=hedges.d), size = 4, col = "red", pch=15)+
  geom_errorbar(data=rc.mod.results, aes(x = risk.cue, ymin = hedges.d-se, ymax = hedges.d+se), 
                col= "red", width = 0.15)+
  ylab("Hedges d")+
  xlab("Risk Cue")+
  scale_x_discrete(labels = c("Enviromental", "Lethal Pred.", "Non-Lethal Pred."))+
  theme(legend.position = "right", legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

## model
data$risk.cue <- relevel(data$risk.cue, ref = "enviro")

rc.mod <- rma.mv(hedges.d, variance, mods = ~factor(risk.cue),
                 random = list(~1|study),
                 W = 1/variance, method = "REML", data = data)
summary(rc.mod)

### generate data to plot
risk.cue <- c("enviro", "lethal", "non.lethal")
hedges.d <- c(0.54, 0.70, 0.78)
se <- c(0.29, 0.36, 0.44)


rc.mod.results <- cbind.data.frame(risk.cue, hedges.d, se)

as.data.frame(rc.mod.results, header = TRUE)



### risk.enviro.pred cue

ep.mod <- rma.mv(hedges.d, variance, mods = ~factor(risk.enviro.pred),
                 random = list(~1|study),
                 W = 1/variance, method = "REML", data = data)
summary(ep.mod)


ggplot(data, aes(x = risk.enviro.pred, y = hedges.d))+
  geom_point(aes(size = 1/variance))+
  theme_classic()


## let's see if there's an effect of risk cue on ES within taxa groups

## birds - no
ggplot(data = subset(data, animal.type == "Bird"), aes (x = risk.cue, y = hedges.d))+
  geom_jitter(width = 0.15, aes(size=1/variance), alpha = 0.5)+
  theme_classic()

## fish - no
ggplot(data = subset(data, animal.type == "Fish"), aes (x = risk.cue, y = hedges.d))+
  geom_jitter(width = 0.15, aes(size=1/variance), alpha = 0.5)+
  theme_classic()

rc.mod.fish <- rma.mv(hedges.d, variance, mods = ~factor(risk.cue),
                 random = list(~1|study),
                 W = 1/variance, method = "REML", data = subset(data, animal.type == "Fish"))
summary(rc.mod.fish)

## herps - no
ggplot(data = subset(data, animal.type == "Herp"), aes (x = risk.cue, y = hedges.d))+
  geom_jitter(width = 0.15, aes(size=1/variance), alpha = 0.5)+
  theme_classic()

## inverts - NON-LETHAL CUES HAVE LARGER EFFECT THAN LETHAL CUES
ggplot(data = subset(data, animal.type == "Invert"), aes (x = risk.cue, y = hedges.d))+
  geom_jitter(width = 0.15, aes(size=1/variance, col = study), alpha = 0.5)+
  theme_classic()

rc.mod.invert <- rma.mv(hedges.d, variance, mods = ~factor(risk.cue),
                      random = list(~1|study),
                      W = 1/variance, method = "REML", data = subset(data, animal.type == "Invert"))
summary(rc.mod.invert)

## mammals - no
ggplot(data = subset(data, animal.type == "Mammal"), aes (x = risk.cue, y = hedges.d))+
  geom_jitter(width = 0.15, aes(size=1/variance), alpha = 0.5)+
  theme_classic()

## shellfish - no - all non-lethal
ggplot(data = subset(data, animal.type == "Shellfish"), aes (x = risk.cue, y = hedges.d))+
  geom_jitter(width = 0.15, aes(size=1/variance), alpha = 0.5)+
  theme_classic()




##############################################################################################
######################## STUDY DESIGN (MANUSCRIPT FIGURE) ###################################
##############################################################################################

## FIELD.LAB

fl.mod1 <- rma.mv(hedges.d, variance, mods = ~factor(field.lab), 
                  random = list(~1|study),
                  W = 1/variance, method = "REML", data = data)
summary(fl.mod1)

# generate dataframe with results for fl.mod1
field.lab <- c("field", "field.enclosure", "lab")
hedges.d <- c(0.50, 0.68, 0.89)
se <- c(0.28, 0.44, 0.48)

fl.mod1.results <- cbind.data.frame(field.lab, hedges.d, se)

as.data.frame(fl.mod1.results, header = TRUE)

f.l <- ggplot(data, aes(x=field.lab, y=hedges.d))+
  geom_jitter(width = 0.2, alpha = 0.15, aes(size=(1/variance)))+
  theme_classic()+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20),
        axis.title.y = element_blank(), legend.position = "none")+
  #geom_abline(intercept = -1.6, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = -0.7, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = -0.22, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = 0.26, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  geom_point(data = fl.mod1.results, aes(x = field.lab, y = hedges.d), size = 4, col = "black", pch = 15)+
  geom_errorbar(data=fl.mod1.results, aes(x = field.lab, ymin = hedges.d-se, ymax = hedges.d+se), 
                col= "black", width = 0.1)+
  ylab("Hedges d")+
  xlab("Experimental Design")+
  scale_x_discrete(labels = c("Field", "Field Enclosure", "Lab"))

## RISK CUE

## model
data$risk.cue <- relevel(data$risk.cue, ref = "lethal")

rc.mod <- rma.mv(hedges.d, variance, mods = ~factor(risk.cue),
                 random = list(~1|study),
                 W = 1/variance, method = "REML", data = data)
summary(rc.mod)

### generate data to plot
risk.cue <- c("enviro", "lethal", "non.lethal")
hedges.d <- c(0.54, 0.70, 0.78)
se <- c(0.29, 0.36, 0.44)


rc.mod.results <- cbind.data.frame(risk.cue, hedges.d, se)

as.data.frame(rc.mod.results, header = TRUE)

r.c <- ggplot(data, aes (x = risk.cue, y = hedges.d))+
  geom_jitter(width = 0.2, aes(size=1/variance), alpha = 0.15)+
  theme_classic()+
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16),
        axis.title.y = element_blank(), legend.position = "none")+
  #geom_abline(intercept = -1.6, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = -0.7, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = -0.22, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = 0.26, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  theme(legend.position = "none")+
  geom_point(data = rc.mod.results, aes(x = risk.cue, y=hedges.d), size = 4, col = "black", pch=15)+
  geom_errorbar(data=rc.mod.results, aes(x = risk.cue, ymin = hedges.d-se, ymax = hedges.d+se), 
                col= "black", width = 0.15)+
  ylab("Hedges d")+
  xlab("Risk Cue")+
  scale_x_discrete(labels = c("Enviro.", "Lethal", "Non-Lethal"))

## SEVERITY SCORE

## model
ss.mod <- rma.mv(hedges.d, variance, mods = ~factor(severity.score),
                 random = list(~1|study),
                 W = 1/variance, method = "REML", data = data)
summary(ss.mod)


### generate data to plot
severity.score <- c("0", "1", "2")
hedges.d <- c(0.56, 0.76, 0.74)
se <- c(0.27, 0.43, 0.54)


ss.mod.results <- cbind.data.frame(severity.score, hedges.d, se)

as.data.frame(ss.mod.results, header = TRUE)

s.s <- ggplot(data, aes (x = severity.score, y = hedges.d))+
  geom_jitter(width = 0.15, aes(size=1/variance), alpha = 0.15)+
  theme_classic()+
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16),
        axis.title.y = element_blank(), legend.position = "none")+
  #geom_abline(intercept = -1.6, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = -0.7, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = -0.22, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = 0.26, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  geom_point(data = ss.mod.results, aes(x = severity.score, y=hedges.d), size = 4, col = "black", pch=15)+
  geom_errorbar(data=ss.mod.results, aes(x = severity.score, ymin = hedges.d-se, ymax = hedges.d+se), 
                col= "black", width = 0.15)+
  ylab("Hedges d")+
  xlab("Severity Score")+
scale_x_discrete(labels = c("0", "1", "2"))

## FORAGING METRIC

# model assessing differences in effect size among response types
rt.mod1 <- rma.mv(hedges.d, variance, mods = ~factor(response.type),
                  random = list(~1|study),
                  W = 1/variance, method = "REML", data = data)
summary(rt.mod1)

# generating a dataset with means, se, upper, lower from rt.mod.1

response.type <- c("intake", "other", "time", "visits")
hedges.d <- c(0.68, 0.81, 0.44, 0.43)
se <- c(0.18, 0.42, 0.36, 0.45)


rt.mod1.results <- cbind.data.frame(response.type, hedges.d, se)

as.data.frame(rt.mod1.results, header = TRUE)

f.m <- ggplot(data, aes (x = response.type, y = hedges.d))+
  geom_jitter(width = 0.15, aes(size=1/variance), alpha = 0.15)+
  theme_classic()+
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16),
        axis.title.y = element_blank(), legend.position = "none")+
  #geom_abline(intercept = -1.6, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = -0.7, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = -0.22, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  #geom_abline(intercept = 0.26, slope = 0, col = "grey", linetype = "dashed", size = 1)+
  geom_point(data = rt.mod1.results, aes(x = response.type, y=hedges.d), size = 4, col = "black", pch=15)+
  geom_errorbar(data=rt.mod1.results, aes(x = response.type, ymin = hedges.d-se, ymax = hedges.d+se), 
                col= "black", width = 0.15)+
  ylab("Hedges d")+
  xlab("Foraging Metric")+
scale_x_discrete(labels = c("Consump.", "Other", "Time", "Visits"))

study.design.figure <- grid.arrange(arrangeGrob(f.l, r.c, s.s, f.m, ncol = 2, 
                          left = textGrob("Effect Size", rot = 90, gp=gpar(fontsize = 20))))

ggsave(study.design.figure, file = "study.design.jpg", width = 9, height = 9, dpi = 600, units = "in")

##############################################################################################
####################################### MASS - ALL (NO) ######################################
##############################################################################################

#data.one$animal.type <- relevel(data.one$animal.type, ref = "Mammal")

## histograms show we need to use log mass so that data looks normal
hist(data$mass)
hist(data$log.mass)

## graph of mass by hedges d
ggplot(data, aes(x = log.mass, y = hedges.d))+
  geom_point(aes(size = 1/variance, col = animal.type), alpha = 0.7)+
  theme_classic()

data$animal.type <- relevel(data$animal.type, ref = "Mammal")



### model mass, all taxa - no relationship
## these are stats for no relationship between mass and hd across taxa
mass.mod <- rma.mv(hedges.d, variance,
                      mods = ~log.mass + animal.type,
                      random = ~1|study,
                      W = 1/variance, method = "REML", data = data)

summary(mass.mod)

## ranges of mass for discussion

## 76 studies, 264 effects that include mass for the prey species
mass.data <- filter(data, mass > 0)
table(droplevels(mass.data$study))
nlevels(factor(mass.data$study))
nrow(mass.data)

## 39 studies, 151 effects that include mass between 10 and 100

# between vole and a rat?

# over 50% of studies and almost 60% of effects
## from a small mouse to a chipmunk
mass.data.1 <- filter(data, mass > 10 & mass < 100)
table(droplevels(mass.data.1$study))
nlevels(factor(mass.data.1$study))
nrow(mass.data.1)

## 60 studies, 233 effects that include mass under 100g

## chipmunk or smaller?

mass.data.2 <- filter(data, mass < 100)
table(droplevels(mass.data.2$study))
nlevels(factor(mass.data.2$study))
nrow(mass.data.2)


##############################################################################################
####################################### MASS - TAXA (NO) #####################################
##############################################################################################

### model: mass and animal.type (all taxa) - no patterns

#data.one.at$animal.type <- relevel(data.one.at$animal.type, ref = "Fish")
#data.one$animal.type <- relevel(data.one$animal.type, ref = "Fish")

mass.mod.at <- rma.mv(hedges.d, variance,
                      mods = ~log.mass * animal.type,
                      random = list(~1|study),
                      W = 1/variance, method = "REML", data = data)

summary(mass.mod.at)


## graph all taxa


## all 
all <- ggplot(data, aes(x = log.mass, y = hedges.d))+
  geom_point(aes(size = 1/variance, col = animal.type))+
  geom_smooth(method = "lm")+
  theme_update(plot.title = element_text(hjust = 0.5, size = 20))+
  #ggtitle("ALL TAXA")+
  theme_classic()+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18))
  #coord_cartesian(xlim = c(-1 , 4.8), ylim = c(-2.8, 2.8))+
  #theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        #axis.line=element_line(color='black'),panel.background = element_blank())
#theme(legend.position = "nonw")

## by animal.type
all+facet_wrap(~animal.type, ncol = 2)

## by homeo vs. poik
all + facet_grid(.~h.p)


##############################################################################################
####################################### MASS - HP (NO) #######################################
##############################################################################################

## model mass x h.p - no interaction
mass.mod.hp <- rma.mv(hedges.d, variance,
                   mods = ~log.mass * h.p,
                   random = list(~1|study),
                   W = 1/variance, method = "REML", data = data)

summary(mass.mod.hp)

## test homeotherms and poikilotherms on their own. No relationship
mass.mod.homeo <- rma.mv(hedges.d, variance,
                   mods = ~log.mass,
                   random = list(~1|study),
                   W = 1/variance, method = "REML", data = subset(data, h.p == "homeo"))

summary(mass.mod.homeo)


mass.mod.poik <- rma.mv(hedges.d, variance,
                        mods = ~log.mass,
                        random = list(~1|study),
                        W = 1/variance, method = "REML", data = subset(data, h.p == "poik"))

summary(mass.mod.poik)










##############################################################################################
##################################### MASS - BIRDS (NO) ######################################
##############################################################################################

## model: mass (birds only)  

mass.mod.birds <- rma.mv(hedges.d, variance,
                        mods = ~log.mass,
                        random = list(~1|study),
                        W = 1/variance, method = "REML", data = subset(data, animal.type == "Bird"))

summary(mass.mod.birds)



## generate data from model to plot on graph

new.mass.bird <- seq(0.95, 2.6, by = 0.02)  
pred.mass.bird <- predict(mass.mod.birds, newmods = new.mass.bird, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred

newdat.mass.bird <- cbind(new.mass.bird, pred.mass.bird$pred, pred.mass.bird$ci.lb, pred.mass.bird$ci.ub)

# convert newdat from matrix to dataframe
newdat.mass.bird <- as.data.frame(newdat.mass.bird)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.mass.bird) <- c("log.mass", "hedges.d", "ci.lb", "ci.ub")


## graph mass by effect size for birds with modelled regression

bird.graph <- ggplot(data = subset(data, animal.type == "Bird"), aes(x = log.mass, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.mass.bird, size = 1)+
  geom_ribbon(data = newdat.mass.bird, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  ggtitle("Birds", subtitle = "p = 0.57")+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16), 
        plot.title = element_text(size = 16, hjust = 0.5), 
        plot.subtitle = element_text(size = 12, hjust = 0.5))+
  coord_cartesian(xlim = c(-1,5), ylim = c(-3.0,10))


##############################################################################################
#################################### MASS MAMMALS (NO) #######################################
##############################################################################################

## model: mass (mammals only)  

mass.mod.mamm <- rma.mv(hedges.d, variance,
                         mods = ~log.mass,
                         random = list(~1|study),
                         W = 1/variance, method = "REML", data= subset(data, animal.type == "Mammal"))

summary(mass.mod.mamm)



## generate data from model to plot on graph

new.mass.mamm <- seq(0.70, 4.8, by = 0.1)  
pred.mass.mamm <- predict(mass.mod.mamm, newmods = new.mass.mamm, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred

newdat.mass.mamm <- cbind(new.mass.mamm, pred.mass.mamm$pred, pred.mass.mamm$ci.lb, pred.mass.mamm$ci.ub)

# convert newdat from matrix to dataframe
newdat.mass.mamm <- as.data.frame(newdat.mass.mamm)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.mass.mamm) <- c("log.mass", "hedges.d", "ci.lb", "ci.ub")


## graph mass by effect size for birds with modelled regression

mamm.graph <- ggplot(data = subset(data, animal.type == "Mammal"), aes(x = log.mass, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.mass.mamm, size = 1)+
  geom_ribbon(data = newdat.mass.mamm, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  ggtitle("Mammals", subtitle = "p = 0.97")+
  theme(axis.text = element_text(size = 16), axis.title = element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5), 
        plot.subtitle = element_text(size = 12, hjust = 0.5))+
  coord_cartesian(xlim = c(-1,5), ylim = c(-3.0,10))


##############################################################################################
######################################## MASS FISH (NO) #####################################
##############################################################################################

## model: mass (fish only)  

mass.mod.fish <- rma.mv(hedges.d, variance,
                        mods = ~log.mass,
                        random = list(~1|study),
                        W = 1/variance, method = "REML", data = subset(data, animal.type == "Fish"))

summary(mass.mod.fish)



## generate data from model to plot on graph

new.mass.fish <- seq(-0.85, 1.0, by = 0.05)  
pred.mass.fish <- predict(mass.mod.fish, newmods = new.mass.fish, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred

newdat.mass.fish <- cbind(new.mass.fish, pred.mass.fish$pred, pred.mass.fish$ci.lb, pred.mass.fish$ci.ub)

# convert newdat from matrix to dataframe
newdat.mass.fish <- as.data.frame(newdat.mass.fish)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.mass.fish) <- c("log.mass", "hedges.d", "ci.lb", "ci.ub")


## graph mass by effect size for birds with modelled regression

fish.graph <- ggplot(data = subset(data, animal.type == "Fish"), aes(x = log.mass, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.mass.fish, size = 1)+
  geom_ribbon(data = newdat.mass.fish, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  ggtitle("Fish", subtitle = "p = 0.96")+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16), 
        plot.title = element_text(size = 16, hjust = 0.5), 
        plot.subtitle = element_text(size = 12, hjust = 0.5))+
  coord_cartesian(xlim = c(-1,5), ylim = c(-3.0,10))



##############################################################################################
######################################## HERPS (NO) ##########################################
##############################################################################################
## model: mass (herps only)  

mass.mod.herp <- rma.mv(hedges.d, variance,
                        mods = ~log.mass,
                        random = list(~1|study),
                        W = 1/variance, method = "REML", data = subset(data, animal.type == "Herp"))

summary(mass.mod.herp)



## generate data from model to plot on graph

new.mass.herp <- seq(-0.80, 3.1, by = 0.05)  
pred.mass.herp <- predict(mass.mod.herp, newmods = new.mass.herp, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred

newdat.mass.herp <- cbind(new.mass.herp, pred.mass.herp$pred, pred.mass.herp$ci.lb, pred.mass.herp$ci.ub)

# convert newdat from matrix to dataframe
newdat.mass.herp <- as.data.frame(newdat.mass.herp)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.mass.herp) <- c("log.mass", "hedges.d", "ci.lb", "ci.ub")


## graph mass by effect size for birds with modelled regression

herp.graph <- ggplot(data = subset(data, animal.type == "Herp"), aes(x = log.mass, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.mass.herp, size = 1)+
  geom_ribbon(data = newdat.mass.herp, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  ggtitle("Herpetofauna", subtitle = "p = 0.31")+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16), 
        plot.title = element_text(size = 16, hjust = 0.5), 
        plot.subtitle = element_text(size = 12, hjust = 0.5))+
 coord_cartesian(xlim = c(-1,5), ylim = c(-3.0,10))

##############################################################################################
############ GRAPH - HEDGES X MASS FOR EA. TAXA (MANUSCRIPT FIGURE) ##########################
##############################################################################################


## Combine graphs in one frame
mass.taxa.figure <- grid.arrange(arrangeGrob(bird.graph, fish.graph, herp.graph, mamm.graph, ncol = 2,
                         left = textGrob("Effect Size", rot = 90, gp=gpar(fontsize = 20)),
                         bottom = textGrob("Log Mass (g)", gp=gpar(fontsize = 20))))

ggsave(mass.taxa.figure, file = "mass.taxa.jpg", width = 9, height = 9, dpi = 600, units = "in")


##############################################################################################
################################ METABOLIC RATE - ALL (NO) ####################################
##############################################################################################

## histograms show we have to use log of mr
hist(data$mr)
hist(data$log.mr)

## all 
all <- ggplot(data, aes(x = log.mr, y = hedges.d))+
  geom_point(aes(size = 1/variance, col = animal.type))+
  #scale_color_manual(values = my.colors)+
  theme_update(plot.title = element_text(hjust = 0.7, size = 20))+
  #ggtitle("ALL TAXA")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18))+
  xlab("Log Metabolic Rate")+
  ylab("Hedges d")+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        axis.line=element_line(color='black'),panel.background = element_blank())
  #theme(legend.position = "nonw")

## by animal.type
#all+facet_wrap(~animal.type)

## by homeo vs. poik
#all + facet_grid(.~h.p)


## model: mr 

mr.mod <- rma.mv(hedges.d, variance,
                    mods = ~log.mr + animal.type,
                    random = list(~1|study),
                    W = 1/variance, method = "REML", data = data)

summary(mr.mod)

## model: mr and animal type

mr.mod.hp <- rma.mv(hedges.d, variance,
                    mods = ~log.mr * h.p,
                    random = list(~1|study),
                    W = 1/variance, method = "REML", data = data)

summary(mr.mod.hp)

## model: mr and animal type

mr.mod.at <- rma.mv(hedges.d, variance,
                  mods = ~log.mr * animal.type,
                  random = list(~1|study),
                  W = 1/variance, method = "REML", data = data)

summary(mr.mod.at)

###############################################################################################
################################# METABOLIC RATE - HP (NO) ###################################
###############################################################################################


## model mass x h.p - no interaction
mr.mod.hp <- rma.mv(hedges.d, variance,
                      mods = ~log.mr * h.p,
                      random = list(~1|study),
                      W = 1/variance, method = "REML", data = data)

summary(mr.mod.hp)

## test homeotherms and poikilotherms on their own. No relationship
mr.mod.homeo <- rma.mv(hedges.d, variance,
                         mods = ~log.mr,
                         random = list(~1|study),
                         W = 1/variance, method = "REML", data = subset(data, h.p == "homeo"))

summary(mr.mod.homeo)


mr.mod.poik <- rma.mv(hedges.d, variance,
                        mods = ~log.mr,
                        random = list(~1|study),
                        W = 1/variance, method = "REML", data = subset(data, h.p == "poik"))

summary(mr.mod.poik)








##############################################################################################
################################ METABOLIC RATE - BIRDS (NO) #################################
##############################################################################################


### just birds - nope
mr.mod.bird <- rma.mv(hedges.d, variance,
                    mods = ~log.mr,
                    random = list(~1|study),
                    W = 1/variance, method = "REML", data = subset(data, animal.type == "Bird"))

summary(mr.mod.bird)


## generate data from model to plot on graph

new.mr.bird <- seq(-0.7, 0.21, by = 0.01)  
pred.mr.bird <- predict(mr.mod.bird, newmods = new.mr.bird, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred

newdat.mr.bird <- cbind(new.mr.bird, pred.mr.bird$pred, pred.mr.bird$ci.lb, pred.mr.bird$ci.ub)

# convert newdat from matrix to dataframe
newdat.mr.bird <- as.data.frame(newdat.mr.bird)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.mr.bird) <- c("log.mr", "hedges.d", "ci.lb", "ci.ub")


## graph mass by effect size for birds with modelled regression

bird.graph.mr <- ggplot(data = subset(data, animal.type == "Bird"), aes(x = log.mr, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.mr.bird, size = 1)+
  geom_ribbon(data = newdat.mr.bird, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  ggtitle("Birds", subtitle = "p = 0.67")+
  theme(axis.title = element_blank(),axis.text = element_text(size = 16), 
        plot.title = element_text(size = 16, hjust = 0.5), 
        plot.subtitle = element_text(size = 12, hjust = 0.5))+
  coord_cartesian(xlim = c(-4.5,2.3), ylim = c(-3.0,10))+
  scale_x_continuous(breaks = c(-4.0, -2.0, 0, 2.0))




##############################################################################################
################################ METABOLIC RATE - MAMMALS (NO) ###############################
##############################################################################################

### just mammals - nope
mr.mod.mamm <- rma.mv(hedges.d, variance,
                      mods = ~log.mr,
                      random = list(~1|study),
                      W = 1/variance, method = "REML", data = subset(data, animal.type == "Mammal"))

summary(mr.mod.mamm)


## generate data from model to plot on graph

new.mr.mamm <- seq(-0.95, 2.2, by = 0.02)  
pred.mr.mamm <- predict(mr.mod.mamm, newmods = new.mr.mamm, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred

newdat.mr.mamm <- cbind(new.mr.mamm, pred.mr.mamm$pred, pred.mr.mamm$ci.lb, pred.mr.mamm$ci.ub)

# convert newdat from matrix to dataframe
newdat.mr.mamm <- as.data.frame(newdat.mr.mamm)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.mr.mamm) <- c("log.mr", "hedges.d", "ci.lb", "ci.ub")


## graph mass by effect size for birds with modelled regression

mamm.graph.mr <- ggplot(data = subset(data, animal.type == "Mammal"), aes(x = log.mr, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.mr.mamm, size = 1)+
  geom_ribbon(data = newdat.mr.mamm, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  ggtitle("Mammals", subtitle = "p = 0.72")+
  theme(axis.title = element_blank(),axis.text = element_text(size = 16), 
        plot.title = element_text(size = 16, hjust = 0.5), 
        plot.subtitle = element_text(size = 12, hjust = 0.5))+
  coord_cartesian(xlim = c(-4.5,2.3), ylim = c(-3.0,10))+
  scale_x_continuous(breaks = c(-4.0, -2.0, 0, 2.0))



##############################################################################################
################################ METABOLIC RATE - FISH (NO) ##################################
##############################################################################################


## just fish- nope

mr.mod.fish <- rma.mv(hedges.d, variance,
                      mods = ~log.mr,
                      random = list(~1|study),
                      W = 1/variance, method = "REML", data = subset(data, animal.type == "Fish"))

summary(mr.mod.fish)


## generate data from model to plot on graph

new.mr.fish <- seq(-4.3, -2.6, by = 0.02)  
pred.mr.fish <- predict(mr.mod.fish, newmods = new.mr.fish, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred

newdat.mr.fish <- cbind(new.mr.fish, pred.mr.fish$pred, pred.mr.fish$ci.lb, pred.mr.fish$ci.ub)

# convert newdat from matrix to dataframe
newdat.mr.fish <- as.data.frame(newdat.mr.fish)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.mr.fish) <- c("log.mr", "hedges.d", "ci.lb", "ci.ub")


## graph mass by effect size for birds with modelled regression

fish.graph.mr <- ggplot(data = subset(data, animal.type == "Fish"), aes(x = log.mr, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.mr.fish, size = 1)+
  geom_ribbon(data = newdat.mr.fish, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  ggtitle("Fish", subtitle = "p = 0.97")+
  theme(axis.title = element_blank(), axis.text = element_text(size = 16), 
        plot.title = element_text(size = 16, hjust = 0.5), 
        plot.subtitle = element_text(size = 12, hjust = 0.5))+
  coord_cartesian(xlim = c(-4.5,2.3), ylim = c(-3.0,10))+
  scale_x_continuous(breaks = c(-4.0, -2.0, 0, 2.0))

##############################################################################################
####################### GRAPH - HEDGES X MR FOR EA.TAXA (MANUSCRIPT FIGURE) ###################
##############################################################################################

mr.taxa.figure <- grid.arrange(arrangeGrob(bird.graph.mr, fish.graph.mr, mamm.graph.mr, ncol = 1,
                         left = textGrob("Effect Size", rot = 90, gp=gpar(fontsize = 18)),
                         bottom = textGrob("Log Metabolic Rate", 
                                           gp=gpar(fontsize = 18))))

ggsave(mr.taxa.figure, file = "mr.taxa.jpg", width = 3, height = 9, dpi = 600, units = "in")

##############################################################################################
####################### MASS SPECIFIC METABOLIC RATE - ALL (NO) ############################
##############################################################################################
## no patterns here.

# histograms show we should use log of ms.mr
hist(data$ms.mr)
hist(data$log.ms.mr)

all <- ggplot(data, aes(x = log.ms.mr, y = hedges.d, col=animal.type))+
  geom_point(aes(size = 1/variance))+
  #scale_color_manual(values = my.colors)+
  #geom_smooth(method = "lm")+
  #theme_update(plot.title = element_text(hjust = 0.7, size = 20))+
  #ggtitle("ALL TAXA")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18))+
  #xlab("log average metabolic rate")+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        axis.line=element_line(color='black'),panel.background = element_blank())+
  theme(legend.position = "nonw")

#### mass specific metabolic rate only - nope

msmr.mod <- rma.mv(hedges.d, variance,
                    mods = ~log.ms.mr + animal.type,
                    random = ~1|study,
                    W = 1/variance, method = "REML", data = data)

summary(msmr.mod)


#### mass specific metabolic rate by hp - nope
msmr.mod.at <- rma.mv(hedges.d, variance,
                      mods = ~log.ms.mr * h.p,
                      random = list(~1|study),
                      W = 1/variance, method = "REML", data = data)

summary(msmr.mod.at)


#### mass specific metabolic rate by animal type - nope
msmr.mod.at <- rma.mv(hedges.d, variance,
                   mods = ~log.ms.mr * animal.type,
                   random = list(~1|study),
                   W = 1/variance, method = "REML", data = data)

summary(msmr.mod.at)




##############################################################################################
####################### MASS SPECIFIC METABOLIC RATE - BIRDS (NO) ############################
##############################################################################################
### birds only

msmr.mod.bird <- rma.mv(hedges.d, variance,
                     mods = ~log.ms.mr,
                     random = list(~1|study),
                     W = 1/variance, method = "REML", data = subset(data, animal.type == "Bird"))

summary(msmr.mod.bird)

## generate data from model to plot on graph

new.msmr.bird <- seq(-2.4, -1.6, by = 0.01)  
pred.msmr.bird <- predict(msmr.mod.bird, newmods = new.msmr.bird, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred

newdat.msmr.bird <- cbind(new.msmr.bird, pred.msmr.bird$pred, pred.msmr.bird$ci.lb, pred.msmr.bird$ci.ub)

# convert newdat from matrix to dataframe
newdat.msmr.bird <- as.data.frame(newdat.msmr.bird)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.msmr.bird) <- c("log.ms.mr", "hedges.d", "ci.lb", "ci.ub")


## graph mass by effect size for birds with modelled regression

bird.graph.msmr <- ggplot(data = subset(data, animal.type == "Bird"), aes(x = log.ms.mr, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.msmr.bird, size = 1)+
  geom_ribbon(data = newdat.msmr.bird, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  xlab("Log Mass SpecificMetabolic Rate")+
  ylab("Hedges d")+
  ggtitle("Birds", subtitle = "p = 0.45")+
  theme(axis.title = element_text(size = 20),axis.text = element_text(size = 16), 
        plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5))
#coord_cartesian(xlim = c(-0.9,4.8), ylim = c(-3.0,2.2))



##############################################################################################
##################### MASS SPECIFIC METABOLIC RATE - MAMMALS (NO) ############################
##############################################################################################

### mammals only

msmr.mod.mamm <- rma.mv(hedges.d, variance,
                        mods = ~log.ms.mr,
                        random = list(~1|study),
                        W = 1/variance, method = "REML", data = subset(data, animal.type == "Mammal"))

summary(msmr.mod.mamm)

## generate data from model to plot on graph

new.msmr.mamm <- seq(-2.7, -1.38, by = 0.02)  
pred.msmr.mamm <- predict(msmr.mod.mamm, newmods = new.msmr.mamm, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred

newdat.msmr.mamm <- cbind(new.msmr.mamm, pred.msmr.mamm$pred, pred.msmr.mamm$ci.lb, pred.msmr.mamm$ci.ub)

# convert newdat from matrix to dataframe
newdat.msmr.mamm <- as.data.frame(newdat.msmr.mamm)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.msmr.mamm) <- c("log.ms.mr", "hedges.d", "ci.lb", "ci.ub")


## graph mass by effect size for birds with modelled regression

mamm.graph.msmr <- ggplot(data = subset(data, animal.type == "Mammal"), aes(x = log.ms.mr, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.msmr.mamm, size = 1)+
  geom_ribbon(data = newdat.msmr.mamm, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  xlab("Log Mass Specific Metabolic Rate")+
  ylab("Hedges d")+
  ggtitle("Mammals", subtitle = "p = 0.42")+
  theme(axis.title = element_text(size = 20),axis.text = element_text(size = 16), 
        plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5))
#coord_cartesian(xlim = c(-0.9,4.8), ylim = c(-3.0,2.2))


##############################################################################################
######################## MASS SPECIFIC METABOLIC RATE - FISH (NO) ############################
##############################################################################################


### fish only
msmr.mod.fish <- rma.mv(hedges.d, variance,
                        mods = ~log.ms.mr,
                        random = list(~1|study),
                        W = 1/variance, method = "REML", data = subset(data, animal.type == "Fish"))

summary(msmr.mod.fish)

## generate data from model to plot on graph

new.msmr.fish <- seq(-3.8, -3.25, by = 0.02)  
pred.msmr.fish <- predict(msmr.mod.fish, newmods = new.msmr.fish, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred

newdat.msmr.fish <- cbind(new.msmr.fish, pred.msmr.fish$pred, pred.msmr.fish$ci.lb, pred.msmr.fish$ci.ub)

# convert newdat from matrix to dataframe
newdat.msmr.fish <- as.data.frame(newdat.msmr.fish)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.msmr.fish) <- c("log.ms.mr", "hedges.d", "ci.lb", "ci.ub")


## graph mass by effect size for birds with modelled regression

fish.graph.msmr <- ggplot(data = subset(data, animal.type == "Fish"), aes(x = log.ms.mr, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.msmr.fish, size = 1)+
  geom_ribbon(data = newdat.msmr.fish, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  xlab("Log Mass Specific Metabolic Rate")+
  ylab("Hedges d")+
  ggtitle("Fish", subtitle = "p = 0.98")+
  theme(axis.title = element_text(size = 20),axis.text = element_text(size = 16), 
        plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5))
#coord_cartesian(xlim = c(-0.9,4.8), ylim = c(-3.0,2.2))







##############################################################################################
####################### GRAPH - HEDGES X MR FOR EA.TAXA ######################################
##############################################################################################

grid.arrange(bird.graph.msmr, mamm.graph.msmr, fish.graph.msmr, ncol = 2)


##############################################################################################
######################## LATITUDE - ALL & EXPLORATION (NO) ###################################
##############################################################################################

## no need to transform latitude
hist(data$latitude)

#data.at$animal.type <- relevel(data.at$animal.type, ref = "Mammal")

all <- ggplot(data, aes(x = latitude, y = hedges.d))+
  geom_point(aes(size = 1/variance))+
  #geom_text(aes(label=e.num), hjust=-2)+
  #scale_color_manual(values = my.colors)+
  #theme_update(plot.title = element_text(hjust = 0.7, size = 20))+
  #ggtitle("ALL TAXA")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18))+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        axis.line=element_line(color='black'),panel.background = element_blank())
  theme(legend.position = "nonw")


## by animal.type
#all+facet_wrap(~animal.type, ncol = 2)

## by homeo vs. poik
#all + facet_grid(.~h.p)

# risk.enviro.pred
#all + facet_grid(~risk.enviro.pred)

## model: latitude - nope
lat.mod.rt <- rma.mv(hedges.d, variance,
                     mods = ~latitude,
                     random = list(~1|study),
                     W = 1/variance, method = "REML", data = data)

summary(lat.mod.rt)


## model: latitude and animal type  - nope  

data$animal.type <- relevel(data$animal.type, ref = "Mammal")

lat.mod.at <- rma.mv(hedges.d, variance,
                     mods = ~latitude + animal.type,
                     random = list(~1|study),
                     W = 1/variance, method = "REML", 
                     data = data)

summary(lat.mod.at)



## model: latitude and enviro / pred cue - nope

lat.mod.ep <- rma.mv(hedges.d, variance,
                     mods = ~latitude * risk.enviro.pred,
                     random = list(~1|study),
                     W = 1/variance, method = "REML", 
                     data = subset(data, risk.enviro.pred == "enviro" | risk.enviro.pred == "pred"))

summary(lat.mod.ep)


## range of latitude - for discussion

min(lat.data$latitude)
max(lat.data$latitude)
mean(lat.data$latitude)

## 86 studies, 311 effects for which we recorded latitude
lat.data <- filter(data, latitude > -100)
table(droplevels(lat.data$study))
nlevels(factor(lat.data$study))
nrow(lat.data)



lat.1 <- filter(lat.data, latitude >= 25 & latitude <= 45)

## 62 studies, 277 effects between 30 - 50 degrees from the equator
# over 70% of studies and almost 90% of all effects
## 30 is Houston Texas and northern Florida
## 50 is just north of the northern boundary of the united states
lat.data.2 <- filter(data, latitude >=30 & latitude <= 50)
table(droplevels(lat.data.2$study))
nlevels(factor(lat.data.2$study))
nrow(lat.data.2)

##############################################################################################
################################ LATITUDE - H.P. (NO) ####################################
############################################################################################

## model: latitude and homeo vs. poik - nope

lat.mod.hp <- rma.mv(hedges.d, variance,
                     mods = ~latitude * h.p,
                     random = list(~1|study),
                     W = 1/variance, method = "REML", data = data)

summary(lat.mod.hp)

## effect of latitude on poikilotherms - nope

mod.poik.lat<- rma.mv(hedges.d, variance, mods = ~latitude,
                      random = ~1|study,
                      W = 1/variance, method = "REML", data = subset(data, h.p == "poik"))

summary(mod.poik.lat)


### just homeotherms - nope
lat.mod.homeo <- rma.mv(hedges.d, variance,
                        mods = ~latitude,
                        random = list(~1|study),
                        W = 1/variance, method = "REML", data = subset(data, h.p == "homeo"))

summary(lat.mod.homeo)




##############################################################################################
################################ LATITUDE - MAMMALS (NO) #####################################
##############################################################################################

## mammals - nope
lat.mod.mamm <- rma.mv(hedges.d, variance,
                       mods = ~latitude,
                       random = list(~1|study),
                       W = 1/variance, method = "REML", data = subset(data, animal.type == "Mammal"))

summary(lat.mod.mamm)

## generate data from model to plot on graph

new.lat.mamm <- seq(6.5, 64, by = 0.5)  
pred.lat.mamm <- predict(lat.mod.mamm, newmods = new.lat.mamm, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred

newdat.lat.mamm <- cbind(new.lat.mamm, pred.lat.mamm$pred, pred.lat.mamm$ci.lb, pred.lat.mamm$ci.ub)

# convert newdat from matrix to dataframe
newdat.lat.mamm <- as.data.frame(newdat.lat.mamm)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.lat.mamm) <- c("latitude", "hedges.d", "ci.lb", "ci.ub")


## graph latitude by effect size for birds with modelled regression

mamm.graph <- ggplot(data = subset(data, animal.type = "Mammal"), aes(x = latitude, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.lat.mamm, size = 1)+
  geom_ribbon(data = newdat.lat.mamm, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  xlab("Distance from Equator (degrees)")+
  ylab("Hedges d")+
  ggtitle("Mammals", subtitle = "p = 0.53")+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16), 
        plot.title = element_text(size = 16, hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5))+
  coord_cartesian(xlim = c(8,62), ylim = c(-2.7,9))

## limited range for disucssion

                  
#############################################################################################
################################ LATITUDE - BIRDS (NO) #######################################
##############################################################################################

## birds - nope
lat.mod.bird <- rma.mv(hedges.d, variance,
                       mods = ~latitude,
                       random = list(~1|study),
                       W = 1/variance, method = "REML", data = subset(data, animal.type == "Bird"))

summary(lat.mod.bird)


## generate data from model to plot on graph

new.lat.bird <- seq(19, 60, by = 0.2)  
pred.lat.bird <- predict(lat.mod.bird, newmods = new.lat.bird, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred

newdat.lat.bird <- cbind(new.lat.bird, pred.lat.bird$pred, pred.lat.bird$ci.lb, pred.lat.bird$ci.ub)

# convert newdat from matrix to dataframe
newdat.lat.bird <- as.data.frame(newdat.lat.bird)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.lat.bird) <- c("latitude", "hedges.d", "ci.lb", "ci.ub")


## graph latitude by effect size for birds with modelled regression

bird.graph <- ggplot(data = subset(data, animal.type == "Bird"), aes(x = latitude, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.lat.bird, size = 1)+
  geom_ribbon(data = newdat.lat.bird, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  xlab("Distance from Equator (degrees)")+
  ylab("Hedges d")+
  ggtitle("Birds", subtitle = "p = 0.94")+
  theme(axis.text = element_text(size = 16),
        axis.title = element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5))+
  coord_cartesian(xlim = c(8,62), ylim = c(-2.7,9))


##############################################################################################
################################ LATITUDE - FISH (NO) ########################################
##############################################################################################

## fish - nope
lat.mod.fish <- rma.mv(hedges.d, variance,
                       mods = ~latitude,
                       random = list(~1|study),
                       W = 1/variance, method = "REML", data = subset(data, animal.type == "Fish"))

summary(lat.mod.fish)


## generate data from model to plot on graph

new.lat.fish <- seq(10, 57, by = 0.5)  
pred.lat.fish <- predict(lat.mod.fish, newmods = new.lat.fish, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred

newdat.lat.fish <- cbind(new.lat.fish, pred.lat.fish$pred, pred.lat.fish$ci.lb, pred.lat.fish$ci.ub)

# convert newdat from matrix to dataframe
newdat.lat.fish <- as.data.frame(newdat.lat.fish)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.lat.fish) <- c("latitude", "hedges.d", "ci.lb", "ci.ub")


## graph latitude by effect size for birds with modelled regression

fish.graph <- ggplot(data = subset(data, animal.type == "Fish"), aes(x = latitude, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.lat.fish, size = 1)+
  geom_ribbon(data = newdat.lat.fish, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  xlab("Distance from Equator (degrees)")+
  ylab("Hedges d")+
  ggtitle("Fish", subtitle = "p = 0.91")+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16), 
        plot.title = element_text(size = 16, hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5))+
  coord_cartesian(xlim = c(8,62), ylim = c(-2.7,9))


##############################################################################################
################################ LATITUDE - INVERTS (YES) ####################################
##############################################################################################

## invert 
lat.mod.invert <- rma.mv(hedges.d, variance,
                       mods = ~latitude,
                       random = list(~1|study),
                       W = 1/variance, method = "REML", data = subset(data, animal.type == "Invert"))

summary(lat.mod.invert)


## generate data from model to plot on graph

new.lat.invert <- seq(22, 42, by = 0.5)  
pred.lat.invert <- predict(lat.mod.invert, newmods = new.lat.invert, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred

newdat.lat.invert <- cbind(new.lat.invert, pred.lat.invert$pred, pred.lat.invert$ci.lb, pred.lat.invert$ci.ub)

# convert newdat from matrix to dataframe
newdat.lat.invert <- as.data.frame(newdat.lat.invert)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.lat.invert) <- c("latitude", "hedges.d", "ci.lb", "ci.ub")


## graph latitude by effect size for invertebrates with modelled regression

invert.graph <- ggplot(data = subset(data, animal.type == "Invert"), aes(x = latitude, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.lat.invert, size = 1)+
  geom_ribbon(data = newdat.lat.invert, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  xlab("Distance from Equator (degrees)")+
  #ylab("Hedges d")+
  ggtitle("Invertebrates", subtitle = "p = 0.0002")+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 16, hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5))+
  coord_cartesian(xlim = c(8,62), ylim = c(-2.7,9))


##############################################################################################
################################ LATITUDE - HERPS (NO) ####################################
##############################################################################################

### herps - nope
lat.mod.herp <- rma.mv(hedges.d, variance,
                         mods = ~latitude,
                         random = list(~1|study),
                         W = 1/variance, method = "REML", data = subset(data, animal.type == "Herp"))

summary(lat.mod.herp)


## generate data from model to plot on graph

new.lat.herp <- seq(10, 43, by = 0.5)  
pred.lat.herp <- predict(lat.mod.herp, newmods = new.lat.herp, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred

newdat.lat.herp <- cbind(new.lat.herp, pred.lat.herp$pred, pred.lat.herp$ci.lb, pred.lat.herp$ci.ub)

# convert newdat from matrix to dataframe
newdat.lat.herp <- as.data.frame(newdat.lat.herp)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.lat.herp) <- c("latitude", "hedges.d", "ci.lb", "ci.ub")



## graph latitude by effect size for herps with modelled regression

herp.graph <- ggplot(data = subset(data, animal.type == "Herp"), aes(x = latitude, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.lat.herp, size = 1)+
  geom_ribbon(data = newdat.lat.herp, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  xlab("Distance from Equator (degrees)")+
  #ylab("Hedges d")+
  ggtitle("Herpetefauna", subtitle = "p = 0.42")+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 16, hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5))+
  coord_cartesian(xlim = c(8,62), ylim = c(-2.7,9))

##############################################################################################
################################ LATITUDE - SHELLFISH (NO) ####################################
##############################################################################################

### shellfish - nope
lat.mod.shellfish <- rma.mv(hedges.d, variance,
                       mods = ~latitude,
                       random = list(~1|study),
                       W = 1/variance, method = "REML", data = subset(data, animal.type == "Shellfish"))

summary(lat.mod.shellfish)



## generate data from model to plot on graph

new.lat.sf <- seq(33, 52, by = 0.5)  
pred.lat.sf <- predict(lat.mod.shellfish, newmods = new.lat.sf, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred

newdat.lat.sf <- cbind(new.lat.sf, pred.lat.sf$pred, pred.lat.sf$ci.lb, pred.lat.sf$ci.ub)

# convert newdat from matrix to dataframe
newdat.lat.sf <- as.data.frame(newdat.lat.sf)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.lat.sf) <- c("latitude", "hedges.d", "ci.lb", "ci.ub")



## graph latitude by effect size for shellfish with modelled regression

shellfish.graph <- ggplot(data = subset(data, animal.type == "Shellfish"), aes(x = latitude, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.lat.sf, size = 1)+
  geom_ribbon(data = newdat.lat.sf, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  xlab("Distance from Equator (degrees)")+
  #ylab("Hedges d")+
  ggtitle("Shellfish", subtitle = "p = 0.62")+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 16, hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5))+
  coord_cartesian(xlim = c(8,62), ylim = c(-2.7,9))

##############################################################################################
##################### LATITUDE - GRAPH - EA. TAXA (MANUSCRIPT FIGURES) #######################
##############################################################################################

lat.taxa.graph <- grid.arrange((arrangeGrob(bird.graph, fish.graph, herp.graph,
                          invert.graph, mamm.graph, shellfish.graph,
                          ncol = 2, 
                          left = textGrob("Hedges d", 
                                          rot = 90, vjust = 1, 
                                          gp=gpar(fontsize=20)),
                          bottom = textGrob("Distance to Equator (degrees)", rot = 0, hjust = 0.5,
                                            gp=gpar(fontsize=20)))))
ggsave(lat.taxa.graph, file = "lat.taxa.jpg", height = 10, width = 6, dpi = 600, units = "in")

##############################################################################################
####################### LATITUDE - MORE EXPLORATION - RISK.CUE (NO) ##########################
##############################################################################################


### Environmental Cues - Not significant alone
lat.mod.enviro <- rma.mv(hedges.d, variance,
                     mods = ~latitude,
                     random = list(~1|study),
                     W = 1/variance, method = "REML", data = subset(data, risk.enviro.pred == "enviro"))


## Predator Cues - Not significant alone
lat.mod.pred <- rma.mv(hedges.d, variance,
                         mods = ~latitude,
                         random = list(~1|study),
                         W = 1/variance, method = "REML", data = subset(data, risk.enviro.pred == "pred"))

summary(lat.mod.pred)







##############################################################################################
############################### LIFE HISTORY - ALL TAXA (NO) #################################
##############################################################################################

## use log of all 4 LH variables
hist(data$longev)
hist(data$log.longev)
hist(data$f.age.mat)
hist(data$log.f.age.mat)
hist(data$c.size)
hist(data$log.c.size)
hist(data$c.year)
hist(data$log.c.year)

## longevity - no
all.longev <- ggplot(data, aes(x = log.longev, y = hedges.d))+
  geom_point(aes(size = 1/variance, col = animal.type))+
  #scale_color_manual(values = my.colors)+
  theme_update(plot.title = element_text(hjust = 0.7, size = 20))+
  #ggtitle("ALL TAXA")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18))+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        axis.line=element_line(color='black'),panel.background = element_blank())
  #theme(legend.position = "nonw")


## model: longevity - all - no
longev.mod <- rma.mv(hedges.d, variance,
                     mods = ~log.longev,
                     random = list(~1|study, ~1|animal.type),
                     W = 1/variance, method = "REML", data = data)

summary(longev.mod)


## female age at first maturity - no
all.fam <- ggplot(data, aes(x = log.f.age.mat, y = hedges.d))+
  geom_point(aes(size = 1/variance, col = animal.type))+
  #scale_color_manual(values = my.colors)+
  theme_update(plot.title = element_text(hjust = 0.7, size = 20))+
  #ggtitle("ALL TAXA")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18))+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        axis.line=element_line(color='black'),panel.background = element_blank())
#theme(legend.position = "nonw")


## model: female age at first maturity - all - no
fam.mod <- rma.mv(hedges.d, variance,
                     mods = ~log.f.age.mat + animal.type,
                     random = list(~1|study),
                     W = 1/variance, method = "REML", 
                  data = subset(data, longev > 0 & c.size > 0 & c.year > 0))

summary(fam.mod)


## clutch size - no
all.cs <- ggplot(data, aes(x = log.c.size, y = hedges.d))+
  geom_point(aes(size = 1/variance, col = animal.type))+
  #scale_color_manual(values = my.colors)+
  theme_update(plot.title = element_text(hjust = 0.7, size = 20))+
  #ggtitle("ALL TAXA")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18))+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        axis.line=element_line(color='black'),panel.background = element_blank())
#theme(legend.position = "nonw")


## model: clutch size - all - no
cs.mod <- rma.mv(hedges.d, variance,
                  mods = ~log.c.size,
                  random = list(~1|study, ~1|animal.type),
                  W = 1/variance, method = "REML", data = data)

summary(cs.mod)


## clutches per year - no
all.cy <- ggplot(data, aes(x = log.c.year, y = hedges.d))+
  geom_point(aes(size = 1/variance, col = animal.type))+
  #scale_color_manual(values = my.colors)+
  theme_update(plot.title = element_text(hjust = 0.7, size = 20))+
  #ggtitle("ALL TAXA")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18))+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        axis.line=element_line(color='black'),panel.background = element_blank())
#theme(legend.position = "nonw")


## model: clutch per year - all - no
cy.mod <- rma.mv(hedges.d, variance,
                 mods = ~log.c.year,
                 random = list(~1|study, ~1|animal.type),
                 W = 1/variance, method = "REML", data = data)

summary(cy.mod)


## put 4 graphs together
grid.arrange(all.longev, all.fam, all.cs, all.cy, ncol = 2)



##############################################################################################
############################### LIFE HISTORY * TAXA - MODELS (NO) ############################
##############################################################################################

## model: longevity and taxa - no
longev.mod.taxa <- rma.mv(hedges.d, variance,
                        mods = ~log.longev * animal.type,
                        random = list(~1|study),
                        W = 1/variance, method = "REML", data = data)
summary(longev.mod.taxa)


## model: fam and taxa - no
fam.mod.taxa <- rma.mv(hedges.d, variance,
                          mods = ~log.f.age.mat * animal.type,
                          random = list(~1|study),
                          W = 1/variance, method = "REML", data = data)
summary(fam.mod.taxa)


## model: c.size and taxa - no
cs.mod.taxa <- rma.mv(hedges.d, variance,
                       mods = ~log.c.size * animal.type,
                       random = list(~1|study),
                       W = 1/variance, method = "REML", data = data)
summary(cs.mod.taxa)


## model: c. per year and taxa - no
cy.mod.taxa <- rma.mv(hedges.d, variance,
                      mods = ~log.c.year * animal.type,
                      random = list(~1|study),
                      W = 1/variance, method = "REML", data = data)
summary(cy.mod.taxa)

##############################################################################################
############################### LIFE HISTORY - BIRDS - MODELS (NO) ###########################
##############################################################################################

## longev: only birds - no
longev.mod.birds <- rma.mv(hedges.d, variance,
                          mods = ~log.longev,
                          random = list(~1|study),
                          W = 1/variance, method = "REML", 
                          data = subset(data, animal.type == "Bird"))
summary(longev.mod.birds)

## fam: only birds - no
fam.mod.birds <- rma.mv(hedges.d, variance,
                           mods = ~log.f.age.mat,
                           random = list(~1|study),
                           W = 1/variance, method = "REML", 
                        data = subset(data, longev > 0,  animal.type == "Bird"))
summary(fam.mod.birds)


## cs: only birds - no
cs.mod.birds <- rma.mv(hedges.d, variance,
                        mods = ~log.c.size,
                        random = list(~1|study),
                        W = 1/variance, method = "REML", data = subset(data, animal.type == "Bird"))
summary(cs.mod.birds)

## cy: only birds - no
cy.mod.birds <- rma.mv(hedges.d, variance,
                       mods = ~log.c.year,
                       random = list(~1|study),
                       W = 1/variance, method = "REML", data = subset(data, animal.type == "Bird"))
summary(cy.mod.birds)


##############################################################################################
############################ LIFE HISTORY - MAMMALS - MODELS (NO) ############################
##############################################################################################

## longev: only mammals - no
longev.mod.mamm <- rma.mv(hedges.d, variance,
                           mods = ~log.longev,
                           random = list(~1|study),
                           W = 1/variance, method = "REML", data = subset(data, animal.type == "Mammal"))
summary(longev.mod.mamm)

## fam: only mammals - no
fam.mod.mamm <- rma.mv(hedges.d, variance,
                        mods = ~log.f.age.mat,
                        random = list(~1|study),
                        W = 1/variance, method = "REML", data = subset(data, animal.type == "Mammal"))
summary(fam.mod.mamm)


## cs: only mammals - no
cs.mod.mamm <- rma.mv(hedges.d, variance,
                       mods = ~log.c.size,
                       random = list(~1|study),
                       W = 1/variance, method = "REML", data = subset(data, animal.type == "Mammal"))
summary(cs.mod.mamm)

## cy: only mammals - no
cy.mod.mamm <- rma.mv(hedges.d, variance,
                       mods = ~log.c.year,
                       random = list(~1|study),
                       W = 1/variance, method = "REML", data = subset(data, animal.type == "Mammal"))
summary(cy.mod.mamm)



##############################################################################################
############################ LIFE HISTORY - LATITUDE - MODELS (NO) ###########################
##############################################################################################


## plot & model: longevity and latitude - no
ggplot(data, aes(x = latitude, y = log.longev))+
  geom_point(aes(col = animal.type))+
  theme_classic()

longev.mod.lat <- rma.mv(hedges.d, variance,
                          mods = ~log.longev + latitude,
                          random = list(~1|study),
                          W = 1/variance, method = "REML", data = data)
summary(longev.mod.lat)


## plot & model: fam and latitude - no
ggplot(data, aes(x = latitude, y = log.f.age.mat))+
  geom_point(aes(col = animal.type))+
  theme_classic()

fam.mod.lat <- rma.mv(hedges.d, variance,
                         mods = ~log.f.age.mat + latitude,
                         random = list(~1|study),
                         W = 1/variance, method = "REML", data = data)
summary(fam.mod.lat)


## plot & model: c.size and latitude - no
ggplot(data, aes(x = latitude, y = log.c.size))+
  geom_point(aes(col = animal.type))+
  theme_classic()

cs.mod.lat <- rma.mv(hedges.d, variance,
                      mods = ~log.c.size + latitude,
                      random = list(~1|study),
                      W = 1/variance, method = "REML", data = data)
summary(cs.mod.lat)

## plot & model: c.year and latitude - no
ggplot(data, aes(x = latitude, y = log.c.year))+
  geom_point(aes(col = animal.type))+
  theme_classic()

cy.mod.lat <- rma.mv(hedges.d, variance,
                     mods = ~log.c.year + latitude,
                     random = list(~1|study),
                     W = 1/variance, method = "REML", data = data)
summary(cy.mod.lat)



##############################################################################################
###################################### GUDs ONLY - (NO) ######################################
##############################################################################################

## guds are 104 of 358 observations
data.gud <- filter(data, response.gud == "yes")

## taxa - no
ggplot(data.gud, aes(x = animal.type, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

gud.at.mod <- rma.mv(hedges.d, variance, mods = animal.type,
                      random = ~1|study,
                      method = "REML", W = 1/variance, data = data.gud)
summary(gud.at.mod)


## field.lab - no
ggplot(data.gud, aes(x = field.lab, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

## risk.cue - no
ggplot(data.gud, aes(x = risk.cue, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

## severity.score - no
ggplot(data.gud, aes(x = severity.score, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

## mass - no
ggplot(data.gud, aes(x = log.mass, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

gud.mass.mod <- rma.mv(hedges.d, variance, mods = log.mass,
                       random = ~1|study,
                       method = "REML", W = 1/variance, data = data.gud)
summary(gud.mass.mod)


## wb metabolism - no
ggplot(data.gud, aes(x = log.mr, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

gud.mr.mod <- rma.mv(hedges.d, variance, mods = log.mr,
                       random = ~1|study,
                       method = "REML", W = 1/variance, data = data.gud)
summary(gud.mr.mod)


## ms metabolism - no
ggplot(data.gud, aes(x = log.ms.mr, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

gud.msmr.mod <- rma.mv(hedges.d, variance, mods = log.ms.mr,
                     random = ~1|study,
                     method = "REML", W = 1/variance, data = data.gud)
summary(gud.msmr.mod)

## latitude - no
ggplot(data.gud, aes(x = latitude, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

gud.lat.mod <- rma.mv(hedges.d, variance, mods = latitude,
                       random = ~1|study,
                       method = "REML", W = 1/variance, data = data.gud)
summary(gud.lat.mod)

## longevity - no
ggplot(data.gud, aes(x = log.longev, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

## fam - no
ggplot(data.gud, aes(x = log.f.age.mat, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

## cs - no
ggplot(data.gud, aes(x = log.c.size, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

## cy - no
ggplot(data.gud, aes(x = log.c.year, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()



##############################################################################################
###################################### INTAKE ONLY - (NO) ######################################
##############################################################################################

# intake is 173 out of 358 effects
intake <- filter(data, response.type == "intake")
#nrow(intake)
#nrow(data)

## taxa - no
ggplot(intake, aes(x = animal.type, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

#intake.at.mod <- rma.mv(hedges.d, variance, mods = animal.type,
                     #random = ~1|study,
                     #method = "REML", W = 1/variance, data = intake)
#summary(intake.at.mod)


## field.lab - no
ggplot(intake, aes(x = field.lab, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

## risk.cue - no
ggplot(intake, aes(x = risk.cue, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

## severity.score - no
ggplot(intake, aes(x = severity.score, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

## mass - no
ggplot(intake, aes(x = log.mass, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

intake.mass.mod <- rma.mv(hedges.d, variance, mods = log.mass,
                       random = list(~1|study, ~1|animal.type),
                       method = "REML", W = 1/variance, data = intake)
summary(intake.mass.mod)


## wb metabolism - no
ggplot(data.gud, aes(x = log.mr, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

intake.mr.mod <- rma.mv(hedges.d, variance, mods = log.mr,
                     random = list(~1|study, ~1|animal.type),
                     method = "REML", W = 1/variance, data = intake)
summary(intake.mr.mod)


## ms metabolism - no
ggplot(data.gud, aes(x = log.ms.mr, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

intake.msmr.mod <- rma.mv(hedges.d, variance, mods = log.ms.mr,
                       random = list(~1|study, ~1|animal.type),
                       method = "REML", W = 1/variance, data = intake)
summary(intake.msmr.mod)

## latitude - no
ggplot(data.gud, aes(x = latitude, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

gud.lat.mod <- rma.mv(hedges.d, variance, mods = latitude,
                      random = ~1|study,
                      method = "REML", W = 1/variance, data = data.gud)
summary(gud.lat.mod)

## longevity - no
ggplot(intake, aes(x = log.longev, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

## fam - no
ggplot(intake, aes(x = log.f.age.mat, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

## cs - no
ggplot(intake, aes(x = log.c.size, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()

## cy - no
ggplot(intake, aes(x = log.c.year, y = hedges.d))+
  geom_point(aes(col = animal.type, size = 1/variance))+
  theme_classic()





##############################################################################################
########################### STUDY DURATION & MEASUREMENT FREQUENCY ###########################
##############################################################################################

## we need to use log of both variables
hist(data$study.duration.days)
hist(data$log.sdd)

hist(data$measurement.frequ.min)
hist(data$log.mfm)


## significant positive relationship between animal mass and duration of foraging measurement
### shorter measurements for smaller animals
mass.mf <- ggplot(data, aes(x = log.mass, y = log.mfm))+
  geom_point(aes(shape = animal.type), size = 3)+
  geom_smooth(method = "lm", col = "black", size = 0.5)+
  theme_classic()+
  xlab("Log Mass (g)")+
  ylab("Log Measurement Frequency (min)")+
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14),
        axis.title.y = element_blank())+
  coord_cartesian(xlim = c(-1,4))+
  theme(legend.text = element_text(size = 14))+
  guides(shape = guide_legend(""))
    

mr.mf <- ggplot(data, aes(x = log.mr, y = log.mfm))+
  geom_point(aes(shape = animal.type), size = 3)+
  geom_smooth(method = "lm", col = "black", size = 0.5)+
  theme_classic()+
  xlab("Log Metabolic Rate (joules/sec)")+
  ylab("Log Measurement Frequency (min)")+
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14),
        axis.title.y = element_blank())+
  theme(legend.text = element_text(size = 14))+
  guides(shape = guide_legend(""))
        
frequ.graphs <- grid.arrange(arrangeGrob(mass.mf, mr.mf, ncol = 1,
                         left = textGrob("Log Measruement Frequency (min)", 
                                         rot = 90, gp=gpar(fontsize = 16))))
                          
ggsave(frequ.graphs, file = "FrequGraphs.jpg", width = 6, height = 8, dpi = 600)


## how many studies do we have mfm??
mfm.subset <- filter(data, log.mfm > -100)
table(droplevels(mfm.subset$study))
nlevels(factor(mfm.subset$study))
nrow(mfm.subset)

## model - measurement frequency x mass
mod.mfm.mass <- lm(log.mfm ~ log.mass, data = data)
summary(mod.mfm.mass)

## model - measurement frequency x metabolism
mod.mfm.met <- lm(log.mfm ~ log.mr, data = data)
summary(mod.mfm.met)

## model - effect size x measurement frequency
mod.mf <- rma.mv(hedges.d, variance, mods = ~log.mfm,
                 random = ~1|study,
                 W = 1/variance, method = "REML", data = data)
summary(mod.mf)

##############################################################################################
###################### FUNNEL PLOT ETC. (MANUSCRIPT FIGURES) #################################
##############################################################################################

mod <- rma.mv(hedges.d, variance, method = "FE", data = data)

funnel.plot <- funnel(mod, xlab = "Effect Size", ylab = "Variance", cex.lab = 2, cex.axis = 1.5)

ss.hd.graph <- ggplot(data, aes(x = hedges.d, y = total.ss))+
  geom_point(aes(size = 1/variance))+
  theme_classic()+
  theme(legend.position = "none")+
  ylab("Total Study Sample Size")+
  xlab("Effect Size")+
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16))
  

#theme(axis.text = element_text(size = 10), axis.title = element_text(size=16))+
  #scale_x_continuous(breaks = round(seq(min(data$hedges.d, max(data$hedges.d, by = 1),1))))+
  
  
  ggsave(ss.hd.graph, file = "Hd_SS.jpg", dpi = 600, width = 6, height = 6, units = "in")


##########################################################################################  
############################# EXPLORATION ############################################### 
##########################################################################################  
  
## no relationship between year and es  
ggplot(data, aes(x = year, y = hedges.d))+
    geom_point(aes(size = 1/variance))+
    theme_classic()

## no relationship between year and ss
ggplot(data, aes(x = year, y = total.ss))+
    geom_point(aes(size = 1/variance))+
    theme_classic()  

## no diff. in es in brown group
ggplot(data, aes(x = brown.group, y = hedges.d))+
  geom_point(aes(size = 1/variance))+
  theme_classic() 
  

##########################################################################################  
############################# OLD CODE - MODEL SELECTION #################################
##########################################################################################
## generate a dataset in which every effect size has a mass, latitude and taxa

data.ms <- subset(data.one, ln.avg.mass > -1)
data.ms <- subset(data.ms, latitude > 0)




mod1 <- rma.mv(ln.hedges.d, ln.variance, mods = ~ln.avg.mass,
                 random = ~1|study,
                 method = "REML", W = 1/ln.variance, data=data.ms)

mod2 <- rma.mv(ln.hedges.d, ln.variance, mods = ~animal.type,
               random = ~1|study,
               method = "REML", W = 1/ln.variance, data=data.ms)

mod3 <- rma.mv(ln.hedges.d, ln.variance, mods = ~latitude,
               random = ~1|study,
               method = "REML", W = 1/ln.variance, data=data.ms)

## mass & taxa

mod4 <- rma.mv(ln.hedges.d, ln.variance, mods = ~ln.avg.mass + animal.type,
               random = ~1|study,
               method = "REML", W = 1/ln.variance, data=data.ms)

mod5 <- rma.mv(ln.hedges.d, ln.variance, mods = ~ln.avg.mass : animal.type,
               random = ~1|study,
               method = "REML", W = 1/ln.variance, data=data.ms)

mod6 <- rma.mv(ln.hedges.d, ln.variance, mods = ~ln.avg.mass * animal.type,
               random = ~1|study,
               method = "REML", W = 1/ln.variance, data=data.ms)


## mass & latitude

mod7 <- rma.mv(ln.hedges.d, ln.variance, mods = ~ln.avg.mass + latitude,
               random = ~1|study,
               method = "REML", W = 1/ln.variance, data=data.ms)

mod8 <- rma.mv(ln.hedges.d, ln.variance, mods = ~ln.avg.mass : latitude,
               random = ~1|study,
               method = "REML", W = 1/ln.variance, data=data.ms)

mod9 <- rma.mv(ln.hedges.d, ln.variance, mods = ~ln.avg.mass * latitude,
                random = ~1|study,
                method = "REML", W = 1/ln.variance, data=data.ms)


## taxa & latitude

mod10 <- rma.mv(ln.hedges.d, ln.variance, mods = ~animal.type + latitude,
                random = ~1|study,
                method = "REML", W = 1/ln.variance, data=data.ms)

mod11 <- rma.mv(ln.hedges.d, ln.variance, mods = ~animal.type : latitude,
                random = ~1|study,
                method = "REML", W = 1/ln.variance, data=data.ms)

mod12 <- rma.mv(ln.hedges.d, ln.variance, mods = ~animal.type * latitude,
                random = ~1|study,
                method = "REML", W = 1/ln.variance, data=data.ms)



## mass, taxa, latitude


mod13 <- rma.mv(ln.hedges.d, ln.variance, mods = ~ln.avg.mass + animal.type + latitude,
                random = ~1|study,
                method = "REML", W = 1/ln.variance, data=data.ms)

mod14 <- rma.mv(ln.hedges.d, ln.variance, mods = ~ln.avg.mass : animal.type + latitude,
                random = ~1|study,
                method = "REML", W = 1/ln.variance, data=data.ms)

mod15 <- rma.mv(ln.hedges.d, ln.variance, mods = ~ln.avg.mass * animal.type + latitude,
                random = ~1|study,
                method = "REML", W = 1/ln.variance, data=data.ms)

mod16 <- rma.mv(ln.hedges.d, ln.variance, mods = ~ln.avg.mass + animal.type : latitude,
                random = ~1|study,
                method = "REML", W = 1/ln.variance, data=data.ms)

mod17 <- rma.mv(ln.hedges.d, ln.variance, mods = ~ln.avg.mass + animal.type * latitude,
                random = ~1|study,
                method = "REML", W = 1/ln.variance, data=data.ms)

mod18 <- rma.mv(ln.hedges.d, ln.variance, mods = ~ln.avg.mass : latitude + animal.type,
                random = ~1|study,
                method = "REML", W = 1/ln.variance, data=data.ms)

mod19 <- rma.mv(ln.hedges.d, ln.variance, mods = ~ln.avg.mass * latitude + animal.type,
                random = ~1|study,
                method = "REML", W = 1/ln.variance, data=data.ms)


## null models

mod.null1 <- rma.mv(ln.hedges.d, ln.variance, method = "FE", W = 1/ln.variance, data=data.ms)

mod.null2 <- rma.mv(ln.hedges.d, ln.variance, random = ~1|study, method = "REML", W = 1/ln.variance, data=data.ms)


m.s.results <- AIC.rma(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10,
             mod11, mod12, mod13, mod14, mod15, mod16, mod17, mod18, mod19, 
             mod.null1, mod.null2)

ggplot(data.one, aes(x = h.p, y = ln.hedges.d))+
  geom_point(aes(size = ln.variance))

### fish and temperature effect

ggplot(data.fish, aes(x = temp, y = ln.avg.mass))+
  geom_point(aes(size = 1/ln.variance))


## no correlation between temp and mass in fish
temp.mass.mod <- lmer(temp~ln.avg.mass + (1|study), data=data.fish)
summary(temp.mass.mod)

## no effect of temp on hedge's d
mod.temp <- rma.mv(ln.hedges.d, ln.variance, mods = ~temp*ln.avg.mass,
                random = ~1|study,
                method = "REML", W = 1/ln.variance, data=data.fish)
summary(mod.temp)














##########################################################################################  
#################################### COUNTING ############################################
##########################################################################################  

## all
# 100 studies, 358 effects
# 30 orders, 88 species
nlevels(factor(data$study))
nrow(data)
nlevels(factor(data$order))
nlevels(factor(data$species))
nlevels(factor(data$class))
levels(factor(data$class))


## latitude range
north.lat <- subset(data, latitude.hem == "north")
south.lat <- subset(data, latitude.hem == "south")

lat <- subset(data, latitude > 0)

max(north.lat$latitude)
max(south.lat$latitude)

nlevels(factor(lat$study))
nrow(lat)
levels(factor(lat$animal.type))

## No. of studies and effects with latitude recorded

## HP

## homeotherms
# 53 studies, 183 effects
# 7 orders, 46 species
homeo <- filter(data, h.p == "homeo") 
table(droplevels(homeo$study))
nlevels(factor(homeo$study))
nrow(homeo)
droplevels(homeo$order)
nlevels(factor(homeo$order))
droplevels(homeo$species)
nlevels(factor(homeo$species))

## poikilotherms
# 47 studies, 175 effects
# 23 orders, 42 species
poik <- filter(data, h.p == "poik") 
table(droplevels(poik$study))
nlevels(factor(poik$study))
nrow(poik)
droplevels(poik$order)
nlevels(factor(poik$order))
droplevels(poik$species)
nlevels(factor(poik$species))



## ANIMAL TYPE

## birds
# 8 studies, 26 effects
# 2 orders, 13 species
birds <- filter(data, animal.type == "Bird")
droplevels(birds$study)
nlevels(factor(birds$study))
nrow(birds)
droplevels(birds$order)
nlevels(factor(birds$order))
droplevels(birds$species)
nlevels(factor(birds$species))

## fish
# 21 studies, 98 effects
# 8 orders, 17 species
fish <- filter(data, animal.type == "Fish")
droplevels(fish$study)
nlevels(factor(fish$study))
nrow(fish)
droplevels(fish$order)
nlevels(factor(fish$order))
droplevels(fish$species)
nlevels(factor(fish$species))

## herps
# 9 studies, 21 effects
# 5 orders, 9 species
herp <- filter(data, animal.type == "Herp")
droplevels(herp$study)
nlevels(factor(herp$study))
nrow(herp)
droplevels(herp$order)
nlevels(factor(herp$order))
droplevels(herp$species)
nlevels(factor(herp$species))

## inverts
# 11 studies, 37 effects
# 7 orders, 11 species
invert <- filter(data, animal.type == "Invert")
droplevels(invert$study)
nlevels(factor(invert$study))
nrow(invert)
droplevels(invert$order)
nlevels(factor(invert$order))
droplevels(invert$species)
nlevels(factor(invert$species))

## mammals
# 45 studies, 157 effects
# 5 orders, 34 species
mamm <- filter(data, animal.type == "Mammal")
droplevels(mamm$study)
nlevels(factor(mamm$study))
nrow(mamm)
droplevels(mamm$order)
nlevels(factor(mamm$order))
droplevels(mamm$species)
nlevels(factor(mamm$species))

## shellfish
# 6 studies, 19 effects
# 3 orders, 5 species
sf <- filter(data, animal.type == "Shellfish")
droplevels(sf$study)
nlevels(factor(sf$study))
nrow(sf)
droplevels(sf$order)
nlevels(factor(sf$order))
droplevels(sf$species)
nlevels(factor(sf$species))


##  ANIMAL CLASSES

# aves
aves <- filter(data, class == "Aves")
droplevels(aves$study)
nlevels(factor(aves$study))
nrow(aves)
droplevels(aves$order)
nlevels(factor(aves$order))
droplevels(aves$species)
nlevels(factor(aves$species))

# reptilia
rep <- filter(data, class == "Reptilia")
droplevels(rep$study)
nlevels(factor(rep$study))
nrow(rep)
droplevels(rep$order)
nlevels(factor(rep$order))
droplevels(rep$species)
nlevels(factor(rep$species))

# amphibia
amph <- filter(data, class == "Amphibia")
droplevels(amph$study)
nlevels(factor(amph$study))
nrow(amph)
droplevels(amph$order)
nlevels(factor(amph$order))
droplevels(amph$species)
nlevels(factor(amph$species))

# insecta
insecta <- filter(data, class == "Insecta")
droplevels(insecta$study)
nlevels(factor(insecta$study))
nrow(insecta)
droplevels(insecta$order)
nlevels(factor(insecta$order))
droplevels(insecta$species)
nlevels(factor(insecta$species))

# actinopterygii
actin <- filter(data, class == "Actinopterygii")
droplevels(actin$study)
nlevels(factor(actin$study))
nrow(actin)
droplevels(actin$order)
nlevels(factor(actin$order))
droplevels(actin$species)
nlevels(factor(actin$species))

# gastropoda
gastro <- filter(data, class == "Gastropoda")
droplevels(gastro$study)
nlevels(factor(gastro$study))
nrow(gastro)
droplevels(gastro$order)
nlevels(factor(gastro$order))
droplevels(gastro$species)
nlevels(factor(gastro$species))

# mollusca
mollusc <- filter(data, class == "Mollusca")
droplevels(mollusc$study)
nlevels(factor(mollusc$study))
nrow(mollusc)
droplevels(mollusc$order)
nlevels(factor(mollusc$order))
droplevels(mollusc$species)
nlevels(factor(mollusc$species))

# bivalvia
bivalve <- filter(data, class == "Bivalvia")
droplevels(bivalve$study)
nlevels(factor(bivalve$study))
nrow(bivalve)
droplevels(bivalve$order)
nlevels(factor(bivalve$order))
droplevels(bivalve$species)
nlevels(factor(bivalve$species))

# malacostraca
mala <- filter(data, class == "Malacostraca")
droplevels(mala$study)
nlevels(factor(mala$study))
nrow(mala)
droplevels(mala$order)
nlevels(factor(mala$order))
droplevels(mala$species)
nlevels(factor(mala$species))

# mammalia
mamm <- filter(data, class == "Mammalia")
droplevels(mamm$study)
nlevels(factor(mamm$study))
nrow(mamm)
droplevels(mamm$order)
nlevels(factor(mamm$order))
droplevels(mamm$species)
nlevels(factor(mamm$species))




### FIELD LAB

#field
field <- filter(data, field.lab == "field")
droplevels(field$study)
nlevels(factor(field$study))
nrow(field)
droplevels(field$order)
nlevels(factor(field$order))
droplevels(field$species)
nlevels(factor(field$species))
levels(factor(field$animal.type))

#field enclosure
fe <- filter(data, field.lab == "field.enclosure")
droplevels(fe$study)
nlevels(factor(fe$study))
nrow(fe)
droplevels(fe$order)
nlevels(factor(fe$order))
droplevels(fe$species)
nlevels(factor(fe$species))
levels(factor(fe$animal.type))

#lab
lab <- filter(data, field.lab == "lab")
droplevels(lab$study)
nlevels(factor(lab$study))
nrow(lab)
droplevels(lab$order)
nlevels(factor(lab$order))
droplevels(lab$species)
nlevels(factor(lab$species))
levels(factor(lab$animal.type))

## RISK CUE

# enviro
enviro <- filter(data, risk.cue == "enviro")
droplevels(enviro$study)
nlevels(factor(enviro$study))
nrow(enviro)
droplevels(enviro$order)
nlevels(factor(enviro$order))
droplevels(enviro$species)
nlevels(factor(enviro$species))
levels(factor(enviro$animal.type))

# non.lethal
non.lethal <- filter(data, risk.cue == "non.lethal")
droplevels(non.lethal$study)
nlevels(factor(non.lethal$study))
nrow(non.lethal)
droplevels(non.lethal$order)
nlevels(factor(non.lethal$order))
droplevels(non.lethal$species)
nlevels(factor(non.lethal$species))
levels(factor(non.lethal$animal.type))

# lethal
lethal <- filter(data, risk.cue == "lethal")
droplevels(lethal$study)
nlevels(factor(lethal$study))
nrow(lethal)
droplevels(lethal$order)
nlevels(factor(lethal$order))
droplevels(lethal$species)
nlevels(factor(lethal$species))
levels(factor(lethal$animal.type))

## RESPONSE TYPE

# intake
intake <- filter(data, response.type == "intake")
droplevels(intake$study)
nlevels(factor(intake$study))
nrow(intake)
droplevels(intake$order)
nlevels(factor(intake$order))
droplevels(intake$species)
nlevels(factor(intake$species))
levels(factor(intake$animal.type))

# time
time <- filter(data, response.type == "time")
droplevels(time$study)
nlevels(factor(time$study))
nrow(time)
droplevels(time$order)
nlevels(factor(time$order))
droplevels(time$species)
nlevels(factor(time$species))
levels(factor(time$animal.type))

# visits
visits <- filter(data, response.type == "visits")
droplevels(visits$study)
nlevels(factor(visits$study))
nrow(visits)
droplevels(visits$order)
nlevels(factor(visits$order))
droplevels(visits$species)
nlevels(factor(visits$species))
levels(factor(visits$animal.type))

# other
other <- filter(data, response.type == "other")
droplevels(other$study)
nlevels(factor(other$study))
nrow(other)
droplevels(other$order)
nlevels(factor(other$order))
droplevels(other$species)
nlevels(factor(other$species))
levels(factor(other$animal.type))



# severity score 0
sszero <- filter(data, severity.score == "0")
droplevels(sszero$study)
nlevels(factor(sszero$study))
nrow(sszero)
droplevels(sszero$order)
nlevels(factor(sszero$order))
droplevels(sszero$species)
nlevels(factor(sszero$species))
levels(factor(sszero$animal.type))

# severity score 1
ssone <- filter(data, severity.score == "1")
droplevels(ssone$study)
nlevels(factor(ssone$study))
nrow(ssone)
droplevels(ssone$order)
nlevels(factor(ssone$order))
droplevels(ssone$species)
nlevels(factor(ssone$species))
levels(factor(ssone$animal.type))

# severity score 2
sstwo <- filter(data, severity.score == "2")
droplevels(sstwo$study)
nlevels(factor(sstwo$study))
nrow(sstwo)
droplevels(sstwo$order)
nlevels(factor(sstwo$order))
droplevels(sstwo$species)
nlevels(factor(sstwo$species))
levels(factor(sstwo$animal.type))

## counting # of species and studies for which we have mass
mass.data <- filter(data, mass > 0)
droplevels(factor(mass.data$study))
nlevels(factor(mass.data$study))
nrow(mass.data)
droplevels(factor(mass.data$animal.type))
levels(factor(mass.data$animal.type))

mamm.mass.data <- filter(mass.data, animal.type == "Mammal")
droplevels(factor(mamm.mass.data$order))
nlevels(factor(mamm.mass.data$order))
droplevels(factor(mamm.mass.data$species))
nlevels(factor(mamm.mass.data$species))

herp.mass.data <- filter(mass.data, animal.type == "Herp")
droplevels(factor(herp.mass.data$order))
nlevels(factor(herp.mass.data$order))
droplevels(factor(herp.mass.data$species))
nlevels(factor(herp.mass.data$species))

fish.mass.data <- filter(mass.data, animal.type == "Fish")
droplevels(factor(fish.mass.data$order))
nlevels(factor(fish.mass.data$order))
droplevels(factor(fish.mass.data$species))
nlevels(factor(fish.mass.data$species))

bird.mass.data <- filter(mass.data, animal.type == "Bird")
droplevels(factor(bird.mass.data$order))
nlevels(factor(bird.mass.data$order))
droplevels(factor(bird.mass.data$species))
nlevels(factor(bird.mass.data$species))


## counting # of species and studies for which we have metabolism
mr.data <- filter(data, mr > 0)
droplevels(factor(mr.data$study))
nlevels(factor(mr.data$study))
nrow(mr.data)
droplevels(factor(mr.data$animal.type))
levels(factor(mr.data$animal.type))

bird.mr.data <- filter(mr.data, animal.type == "Bird")
droplevels(factor(bird.mr.data$order))
nlevels(factor(bird.mr.data$order))
droplevels(factor(bird.mr.data$species))
nlevels(factor(bird.mr.data$species))

fish.mr.data <- filter(mr.data, animal.type == "Fish")
droplevels(factor(fish.mr.data$order))
nlevels(factor(fish.mr.data$order))
droplevels(factor(fish.mr.data$species))
nlevels(factor(fish.mr.data$species))

mamm.mr.data <- filter(mr.data, animal.type == "Mammal")
droplevels(factor(mamm.mr.data$order))
nlevels(factor(mamm.mr.data$order))
droplevels(factor(mamm.mr.data$species))
nlevels(factor(mamm.mr.data$species))

## female age at first maturity
b.fam.data <- filter(data, f.age.mat > 0 & animal.type == "Bird")
droplevels(factor(b.fam.data$order))
nlevels(factor(b.fam.data$order))
droplevels(factor(b.fam.data$species))
nlevels(factor(b.fam.data$species))

m.fam.data <- filter(data, f.age.mat > 0 & animal.type == "Mammal")
droplevels(factor(m.fam.data$order))
nlevels(factor(m.fam.data$order))
droplevels(factor(m.fam.data$species))
nlevels(factor(m.fam.data$species))

## number of studies and effects of inverts with latitude

lat.invert <- filter(data, latitude > -100 & animal.type == "Invert")

##########################################################################################  
################# LIFE HISTORY PCA x ES (MANUSCRIPT FIGURES) #############################  
##########################################################################################  

## PC1 - all taxa
ggplot(data, aes(y = hedges.d, x = pc1))+
  geom_point(aes(size = 1/variance))+
  theme_classic()

pc1.mod <- rma.mv(hedges.d, variance, mods = ~pc1,
                  random = list(~1|study),
                  W = 1/variance, method = "REML", data = data)
summary(pc1.mod)

## PC2 - all taxa
ggplot(data, aes(y = hedges.d, x = pc2))+
  geom_point(aes(size = 1/variance))+
  theme_classic()

pc2.mod <- rma.mv(hedges.d, variance, mods = ~pc2,
                  random = list(~1|study),
                  W = 1/variance, method = "REML", data = data)
summary(pc2.mod)





## PC1 - Birds
pc1.mod.birds <- rma.mv(hedges.d, variance, mods = ~pc1,
                  random = list(~1|study),
                  W = 1/variance, method = "REML", data = subset(data, animal.type == "Bird"))
summary(pc1.mod.birds)

## find range of pc1 values for birds
bird <- subset(data, animal.type == "Bird")
pc1.b <- bird$pc1
pc1.b <- pc1.b[!is.na(pc1.b)]
min(pc1.b)
max(pc1.b)

## generate data from model to plot on graph
new.pc1.bird <- seq(0.6, 2.6, by = 0.02)  
pred.pc1.bird <- predict(pc1.mod.birds, newmods = new.pc1.bird, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred
newdat.pc1.bird <- cbind(new.pc1.bird, pred.pc1.bird$pred, pred.pc1.bird$ci.lb, pred.pc1.bird$ci.ub)

# convert newdat from matrix to dataframe
newdat.pc1.bird <- as.data.frame(newdat.pc1.bird)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.pc1.bird) <- c("pc1", "hedges.d", "ci.lb", "ci.ub")


## graph pc1 by effect size for birds with modelled regression
pc1.bird.graph <- ggplot(data = subset(data, animal.type == "Bird"), aes(x = pc1, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.pc1.bird, size = 1)+
  geom_ribbon(data = newdat.pc1.bird, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  ggtitle("Birds", subtitle = "p = 0.27")+
  #xlab("PC 1")+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16), 
        plot.title = element_text(size = 16, hjust = 0.5), 
        plot.subtitle = element_text(size = 12, hjust = 0.5))+
  coord_cartesian(xlim = c(-3,3.5), ylim = c(-3,8))



## PC1 - Mammals
pc1.mod.mamm <- rma.mv(hedges.d, variance, mods = ~pc1,
                    random = list(~1|study),
                    W = 1/variance, method = "REML", data = subset(data, animal.type == "Mammal"))
summary(pc1.mod.mamm)

## find range of pc1 values for mammals
mamm <- subset(data, animal.type == "Mammal")
pc1.m <- mamm$pc1
pc1.m <- pc1.m[!is.na(pc1.m)]
min(pc1.m)
max(pc1.m)

## generate data from model to plot on graph
new.pc1.mamm <- seq(-2.7, 3.1, by = 0.02)  
pred.pc1.mamm <- predict(pc1.mod.mamm, newmods = new.pc1.mamm, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred
newdat.pc1.mamm <- cbind(new.pc1.mamm, pred.pc1.mamm$pred, pred.pc1.mamm$ci.lb, pred.pc1.mamm$ci.ub)

# convert newdat from matrix to dataframe
newdat.pc1.mamm <- as.data.frame(newdat.pc1.mamm)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.pc1.mamm) <- c("pc1", "hedges.d", "ci.lb", "ci.ub")

## graph pc1 by effect size for mammals with modelled regression
pc1.mamm.graph <- ggplot(data = subset(data, animal.type == "Mammal"), aes(x = pc1, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.pc1.mamm, size = 1)+
  geom_ribbon(data = newdat.pc1.mamm, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  ggtitle("Mammals", subtitle = "p = 0.82")+
  xlab("Principal Component 1")+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 20),
        axis.text = element_text(size = 16), 
        plot.title = element_text(size = 16, hjust = 0.5), 
        plot.subtitle = element_text(size = 12, hjust = 0.5))+
  coord_cartesian(xlim = c(-3,3.5), ylim = c(-3,8))






## PC2 - Birds
pc2.mod.birds <- rma.mv(hedges.d, variance, mods = ~pc2,
                    random = list(~1|study),
                    W = 1/variance, method = "REML", data = subset(data, animal.type == "Bird"))
summary(pc2.mod.birds)

## find range of pc2 values for birds
bird <- subset(data, animal.type == "Bird")
pc2.b <- bird$pc2
pc2.b <- pc2.b[!is.na(pc2.b)]
min(pc2.b)
max(pc2.b)

## generate data from model to plot on graph
new.pc2.bird <- seq(-0.8, 1.56, by = 0.02)  
pred.pc2.bird <- predict(pc2.mod.birds, newmods = new.pc2.bird, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred
newdat.pc2.bird <- cbind(new.pc2.bird, pred.pc2.bird$pred, pred.pc2.bird$ci.lb, pred.pc2.bird$ci.ub)

# convert newdat from matrix to dataframe
newdat.pc2.bird <- as.data.frame(newdat.pc2.bird)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.pc2.bird) <- c("pc2", "hedges.d", "ci.lb", "ci.ub")

## graph pc2 by effect size for birds with modelled regression
pc2.bird.graph <- ggplot(data = subset(data, animal.type == "Bird"), aes(x = pc2, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.pc2.bird, size = 1)+
  geom_ribbon(data = newdat.pc2.bird, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  ggtitle("Birds", subtitle = "p = 0.32")+
  #xlab("PC 2")+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16), 
        plot.title = element_text(size = 16, hjust = 0.5), 
        plot.subtitle = element_text(size = 12, hjust = 0.5))+
  coord_cartesian(xlim = c(-2,3.2), ylim = c(-3,8))



## PC2 - Mammals
pc2.mod.mamm <- rma.mv(hedges.d, variance, mods = ~pc2,
                    random = list(~1|study),
                    W = 1/variance, method = "REML", data = subset(data, animal.type == "Mammal"))
summary(pc2.mod.mamm)

## find range of pc1 values for mammals
mamm <- subset(data, animal.type == "Mammal")
pc2.m <- mamm$pc2
pc2.m <- pc2.m[!is.na(pc2.m)]
min(pc2.m)
max(pc2.m)


## generate data from model to plot on graph
new.pc2.mamm <- seq(-1.8, 3.12, by = 0.02)  
pred.pc2.mamm <- predict(pc2.mod.mamm, newmods = new.pc2.mamm, se.fit = TRUE)

# we need to combine new.mass with pred variable from data named pred
newdat.pc2.mamm <- cbind(new.pc2.mamm, pred.pc2.mamm$pred, pred.pc2.mamm$ci.lb, pred.pc2.mamm$ci.ub)

# convert newdat from matrix to dataframe
newdat.pc2.mamm <- as.data.frame(newdat.pc2.mamm)

# name the columns so that our ggplot commands recognizes the variables
colnames(newdat.pc2.mamm) <- c("pc2", "hedges.d", "ci.lb", "ci.ub")

## graph pc2 by effect size for mammals with modelled regression
pc2.mamm.graph <- ggplot(data = subset(data, animal.type == "Mammal"), aes(x = pc2, y = hedges.d))+
  geom_point(aes(size = 1/variance), alpha = 0.4)+
  geom_line(data = newdat.pc2.mamm, size = 1)+
  geom_ribbon(data = newdat.pc2.mamm, aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15)+
  theme_classic()+
  theme(legend.position = "nonw")+
  ggtitle("Mammals", subtitle = "p = 0.64")+
  xlab("Principal Component 2")+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 20),
        axis.text = element_text(size = 16), 
        plot.title = element_text(size = 16, hjust = 0.5), 
        plot.subtitle = element_text(size = 12, hjust = 0.5))+
  coord_cartesian(xlim = c(-2,3.2), ylim = c(-3,8))

## combine 4 graphs from above into figure for manuscript

pc.hd.graph <- grid.arrange(arrangeGrob(pc1.bird.graph, pc2.bird.graph, pc1.mamm.graph, pc2.mamm.graph, 
                                                ncol = 2, 
                                                left = textGrob("Effect Size", rot = 90, gp=gpar(fontsize = 20))))

ggsave(pc.hd.graph, file = "pc.hd.jpg", width = 9, height = 9, dpi = 600, units = "in")

## Do PCA1 and PCA1 correlate with mass?

ggplot(aes(x = log.mass, y = pc1), data=data)+
  geom_point()+
  theme_classic()

ggplot(aes(x = log.mass, y = pc2), data=data)+
  geom_point()+
  theme_classic()

## sig: pc1 x mass
pc1.mass <- lm(pc1 ~ log.mass, data=data)
summary(pc1.mass)

## sig: pc2 x mass
pc2.mass <- lm(pc2 ~ log.mass, data=data)
summary(pc2.mass)

## Do PCA1 and PCA1 correlate with measurement frequency?

ggplot(aes(x = log.mfm, y = pc1), data=data)+
  geom_point()+
  theme_classic()

ggplot(aes(x = log.mfm, y = pc2), data=data)+
  geom_point()+
  theme_classic()

## sig: pc1 x mfm
pc1.mfm <- lm(pc1 ~ log.mfm, data=data)
summary(pc1.mfm)

## not sig: pc2 x mfm
pc2.mfm <- lm(pc2 ~ log.mfm, data=data)
summary(pc2.mfm)

