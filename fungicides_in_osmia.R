# this is a copy of data_formatting_hq scrit
# but quickly adjusted to look for good osmia sites
# for spring 2022 experiments



#Tobias Mueller
# 2/2022




# load packages
library(ggplot2)
library(reshape2)
library(gridExtra)
library(fitdistrplus)
library(lme4)
library(emmeans)
library(tidyverse)
library(vegan)

#a few packages for graph colors
library(viridis)
library(RColorBrewer)
library(scales)

rm(list = ls())


### read in data --------------------------------------------------

alldf<- read.csv("input/all_pesticide2019.csv")
ld50<- read.csv("input/ld50.csv")
pesticidenames<-read.csv("input/pesticide_type.csv")
gis <- read.csv("input/appleorchard_gis_data2015_3000m.csv")
osmia_mary_fun <- read_csv("input/omsia_mary_heavyedited.csv")
  
#then fix some things
alldf$Sample.type <- as.factor(alldf$Sample.type)
pesticidenames$moa<-as.factor(pesticidenames$moa)
pesticidenames$pesticide<-as.factor(pesticidenames$pesticide)
pesticidenames$type<-as.factor(pesticidenames$type)
pesticidenames$class<-as.factor(pesticidenames$class)


# first lets only look at bees (may change this later)

osmiadf <- alldf %>% 
  filter(Sample.type == "Mason bee pollen provision") %>%
  select(Site:(ncol(.)-2)) %>%
  droplevels(.)



osmialong <- reshape2::melt(osmiadf, id.vars=("Site"))

colnames(osmialong)[2] <- "pesticide"

osmialong <- merge(pesticidenames,osmialong, by = "pesticide", all.y=TRUE)

osmiafun <- osmialong %>%
  filter(type == "Fungicide") %>%
  select(-c(class,moa))

summ<- osmiafun %>%
  group_by(Site)%>%
summarise(sum=sum(value)) %>%
  arrange(-sum)


ggplot() +
  geom_boxplot(aes(x=Site,))


colnames(osmia_mary_fun)[1] <- "pesticide"
colnames(osmia_mary_fun)[2] <- "type"

mary_long <- reshape2::melt(osmia_mary_fun)


summ.mary<- mary_long %>%
  group_by(variable)%>%
  summarise(sum=sum(value, na.rm=TRUE)) %>%
  arrange(-sum)



