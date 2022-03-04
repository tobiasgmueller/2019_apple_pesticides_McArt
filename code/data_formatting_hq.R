# this is a draft r scrip for figuring out HQ of pesticides



#Tobias Mueller
#9/2021


#setwd("C:/Users/obiew/Desktop/scott_pesticide_data")

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

# DATA PREP -------------------
### read in data --------------------------------------------------

alldf<- read.csv("input/all_pesticide2019.csv")
ld50<- read.csv("input/ld50.csv")
pesticidenames<-read.csv("input/pesticide_type.csv")
gis <- read.csv("input/appleorchard_gis_data2015_3000m.csv")

#then fix some things
alldf$Sample.type <- as.factor(alldf$Sample.type)
pesticidenames$moa<-as.factor(pesticidenames$moa)
pesticidenames$pesticide<-as.factor(pesticidenames$pesticide)
pesticidenames$type<-as.factor(pesticidenames$type)
pesticidenames$class<-as.factor(pesticidenames$class)

# first lets only look at bees (may change this later)
alldf1 <- subset(alldf, Sample.type!="Apple flower" &
                           Sample.type!="Bumble bee wax"&
                           Sample.type!="Dandelion flower"&
                           Sample.type!="Nectar"&
                           Sample.type!="Honey bee pollen"&
                           Sample.type!="Mason bee pollen provision")


alldf1<-droplevels(alldf1) # drop unused levels
rownames(alldf1) = seq(length=nrow(alldf1)) # now fix row # after subsetting




### adjusting to HQ --------------------------------------------
# now we'll try to make everything in alldf (ppb) as a percent of ld50 of honey bees
# which would entail dividing alldf by corresponding values of ld50

#drop identifier colummns
justpesticides <- alldf1[7:98]


ld50contact <- ld50[1,2:93]
ld50oral <- ld50[2,2:93]

# okay probably not the cleanest way but we can add the ld50 row to the pesticide df

# then matrix is and devide by the ld 50 row
oral<-rbind(justpesticides,ld50oral)

#row numbers got wonky again
rownames(oral) = seq(length=nrow(oral)) # now fix row # after subsetting


temp = as.matrix(oral)

oralhq<- t(t(temp)/temp[104,])

oralhq = as.data.frame(oralhq)

#now remove the ld 50 row
oralhq <- oralhq[1:103,]

# then we'll set all the inf and NA to 0 (where the ld50 was either 0 or NA)
oralhq[] <- lapply(oralhq, function(x) replace(x, !is.finite(x), 0))


#now do the same for contact
contact<-rbind(justpesticides,ld50contact)

#row numbers got wonky again
rownames(contact) = seq(length=nrow(contact)) # now fix row # after subsetting

temp = as.matrix(contact)

contacthq<- t(t(temp)/temp[104,])

contacthq = as.data.frame(contacthq)


contacthq <- contacthq[1:103,]

contacthq[] <- lapply(contacthq, function(x) replace(x, !is.finite(x), 0))



#now add back descriptors
oraldf<- cbind(alldf1[,1:6],oralhq)
contactdf<- cbind(alldf1[,1:6],contacthq)

#then on second thought lets drop some of those descriptor columns
oraldf<- cbind(alldf1[,c(3,4,6)],oralhq)
contactdf<- cbind(alldf1[,c(3,4,6)],contacthq)


# okay now just quickly double check the ug/g to ppb is correct
# first lets figure out the average weight of a hb in the sample
mean((oraldf$Mass..g.[oraldf$Sample.type == "Apis mellifera worker"]))

# okay .1 seems quite spot on





# now lets add a column to see total additive ld50 
#this assumes pesticides are purely additive
oraldf$chemsum <- rowSums(oraldf[4:95])
contactdf$chemsum <- rowSums(contactdf[4:95])


# and now add a column to see if any pesticides were found or not
oraldf$non_zero <- ifelse(oraldf$chemsum > 0, 1, 0)
contactdf$non_zero<- ifelse(contactdf$chemsum > 0, 1, 0)




#then quickly just check the frequency samples with no pesticides
table(oraldf$non_zero)
table(contactdf$non_zero)
# wow okay incredibly low



### remove undetected pesticides ------------------------------

# I think we should remove pesticides that were never detected
# or that have no ld50
# as theyre not informative in our ld50 data table
# and then just report which pesticides were removed in a table or something

# we can do this by summing columns and removing zeros 

i <- (colSums((oraldf[,4:95]), na.rm=T) != 0) # T if colSum is not 0, F otherwise
cols_to_drop_oral <- c(rep(TRUE, 3), i) # then add 3 T to keep meta data columns


oraldf1 <- oraldf[, cols_to_drop_oral] # all the non-zero columns plus meta data
notfoundpesticides.oral <- oraldf[, !cols_to_drop_oral]  # all the zero columns


# and the same for contact
i <- (colSums((contactdf[,4:95]), na.rm=T) != 0) # T if colSum is not 0, F otherwise
cols_to_drop_contact <- c(rep(TRUE, 3), i) # then add 3 T to keep meta data columns


contactdf1 <- contactdf[, cols_to_drop_contact] # all the non-zero columns plus meta data
notfoundpesticides.contact <- contactdf[, !cols_to_drop_contact]  # all the zero columns

# the below should sum to 97 if everything went well
# ncol(oraldf1) + ncol(notfoundpesticides.oral)
# ncol(contactdf1) + ncol(notfoundpesticides.contact)

# here is the list of pesticides that were not found or have no oral ld50
noimpactoral <- colnames(notfoundpesticides.oral)
noimpactcontact <- colnames(notfoundpesticides.contact)






### melting for graphing ----------------------------------
# now lets melt that data for graphing and playing with
oralmelt <- reshape2::melt(oraldf1[,1:((ncol(oraldf1))-2)], id.vars=c("Mass..g.", "Sample.type", "Site"))

colnames(oralmelt)[4] <- "pesticide"
colnames(oralmelt)[5] <- "oralhq"


contactmelt <- reshape2::melt(contactdf1[1:((ncol(contactdf1))-2)], id.vars=c("Mass..g.", "Sample.type", "Site"))

colnames(contactmelt)[4] <- "pesticide"
colnames(contactmelt)[5] <- "contacthq"


#now I would love to add pesticide type as a category

oralmelt <- merge(pesticidenames,oralmelt, by = "pesticide", all.y=TRUE)
contactmelt <- merge(pesticidenames,contactmelt, by = "pesticide", all.y=TRUE)



### adding GIS information ####

# we have 2015 ag data for all 2019 orchards except 1 (chip bailey) 
gis_sum <- gis[,55:61]

gis_sum$Site <- as.factor(gis_sum$Site)

oralmelt$Site <- as.factor(oralmelt$Site)
contactmelt$Site <- as.factor(contactmelt$Site)
oraldf1$Site <- as.factor(oraldf1$Site)
contactdf1$Site <- as.factor(contactdf1$Site)


oralmelt <- merge(oralmelt,gis_sum, by = "Site", all.x = TRUE)
oraldf1 <- merge(oraldf1,gis_sum, by = "Site", all.x = TRUE)

contactmelt <- merge(contactmelt,gis_sum, by = "Site", all.x = TRUE)
contactdf1 <- merge(contactdf1,gis_sum, by = "Site", all.x = TRUE)

# then lets make a melt of only pesticide findings
oralmelt.findings <- subset(oralmelt, oralmelt$oralhq >0)
contactmelt.findings <- subset(contactmelt, contactmelt$contacthq >0)

# lets drop levels from oralmelt.finding
oralmelt.findings$pesticide<-factor(oralmelt.findings$pesticide)
contactmelt.findings$pesticide<-factor(contactmelt.findings$pesticide)

# okay lets clean things up
# and remove df I dont need

rm(list = c("alldf",
            "alldf1",
            "contact",
            "contacthq",
            "oralhq",
            "oral",
            "ld50",
            "ld50contact",
            "ld50oral",
            "notfoundpesticides.contact",
            "notfoundpesticides.oral",
            "i",
            "pesticidenames",
            "temp",
            "cols_to_drop_oral",
            "cols_to_drop_contact",
            "noimpactcontact",
            "noimpactoral",
            "justpesticides"
            ))





# END OF DATA PREP ####



# some summary tables to looks at the data ####

# counts of samples
table(oraldf$Sample.type)

# count of samples for each sites
table(oraldf$Sample.type,oraldf$Site)

table(oralmelt$Sample.type,oralmelt$class) 
table(contactmelt$Sample.type,contactmelt$class)


# and now avg findings per bee by species
table(oralmelt.findings$Sample.type)/table(oraldf$Sample.type)

plot(table(oralmelt.findings$Sample.type)/table(oraldf$Sample.type))

# which pesticide MoA are showing up
table(oralmelt.findings$Sample.type,oralmelt.findings$moa)




beecountvectororal <- c(19,20,16,16,20,17,12)

# this is now average findings of that class per sample
table(oralmelt.findings$Sample.type,oralmelt.findings$moa)/beecountvectororal

#table(oralmelt.findings$Sample.type,oralmelt.findings$type)/beecountvectororal


aggregate(oraldf1[,35:36], list(oraldf1$Sample.type), mean)


oraldf1 %>%
  group_by(Sample.type) %>% 
  summarise(mean_hq = mean(chemsum), max_hq = max(chemsum),mid_hq = median(chemsum), perc_findings= mean(non_zero))%>%
  ungroup()


pesticide.by.bee<-oralmelt %>%
  group_by(Sample.type,pesticide) %>% 
  summarise(mean_hq = mean(oralhq), max_hq = max(oralhq))%>%
  ungroup()

# now lets make a table of just pesticides  over 3% of LD50
temp<-oralmelt[,0:((ncol(oralmelt)-6))]
ld50.3percoral<-subset(temp, oralhq >= .03)



temp<-ld50.3percoral %>%
  group_by(Sample.type,pesticide) %>% 
  summarise(mean_hq = mean(oralhq), max_hq = max(oralhq))%>%
  ungroup()



# table of which bees have a combined risk of  over epa limit 
# for oral (.3)
chemsumld50.4percoral<-subset(oraldf1, chemsum >= .04)
chemsumld50.4perccont<-subset(contactdf1, chemsum >= .04)

table(chemsumld50.4percoral$Sample.type)
table(chemsumld50.4percoral$Site)
table(chemsumld50.4percoral$Site, chemsumld50.4percoral$Sample.type)


table(chemsumld50.4perccont$Sample.type)
table(chemsumld50.4perccont$Site)

# and for contact





# now lets try to create a table showing all pesticides and what percent of their findings are in each bee group


oralfinding_counts<-as.data.frame.matrix(table(oralmelt.findings$pesticide,oralmelt.findings$Sample.type))



percfindingoralbypesticide<-(t(apply(oralfinding_counts,1, function(x) x/sum(x))))*100
percfindingoralbypesticide # this shows the frequency each pesticide is found in each bee
head(rowSums(percfindingoralbypesticide)) # this should all sum to 100
percfindingoralbypesticide <- as.data.frame(percfindingoralbypesticide) 
percfindingoralbypesticide$pesticide <- row.names(percfindingoralbypesticide) # add row names to df
percfindingoralbypesticidemelt <- reshape2::melt(percfindingoralbypesticide, id.vars=c("pesticide")) # melt for graphing

# okay the below is showing what percent of each pesticide findings is in each bee - that is, eveness of exposure
ggplot(percfindingoralbypesticidemelt, aes(variable, value, color=variable))+
  geom_point()+
  facet_wrap(~pesticide)+
  ylab("percent of pesticide findings")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



percfindingoralbybee <- (oralfinding_counts / rep(colSums(oralfinding_counts), each = nrow(oralfinding_counts)))*100
percfindingoralbybee # this shows the percent of a bees findings that is a certain pesticide
head(colSums(percfindingoralbee)) # this should all sum to 100

percfindingoralbybee <- as.data.frame(percfindingoralbybee)
percfindingoralbybee$pesticide <- row.names(percfindingoralbybee) # add row names to df
percfindingoralbybeemelt <- reshape2::melt(percfindingoralbybee, id.vars=c("pesticide")) # melt for graphing


ggplot(percfindingoralbybeemelt, aes(pesticide, value, color=pesticide))+
  geom_point()+
  facet_wrap(~variable)+
  ylab("frequency of pesticide by bee")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



#attempts at stacked barplot

ggplot(oralmelt.findings, aes(fill=Sample.type, x=pesticide))+
  geom_bar(position="fill")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(oralmelt.findings, aes(fill=pesticide, x=Sample.type))+
  geom_bar(position="fill")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





ggplot(oralmelt.findings, aes(fill=Sample.type, x=pesticide))+
  geom_bar(position="fill")




# not suoer useful but maybe later
#with(oralmelt.findings, table(oralmelt.findings$Sample.type,oralmelt.findings$moa)) %>% prop.table(margin=1)



ggplot(oralmelt.findings, aes(x=Sample.type, y=oralhq, color=moa))+
  geom_boxplot(alpha=.5)+
  facet_wrap(~moa, scales="free_y")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(contactmelt.findings, aes(x=Sample.type, y=contacthq, color=moa))+
  geom_boxplot(alpha=.5)+
  facet_wrap(~moa, scales="free_y")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# looking at hq by site ####

#first create a new df copy and drop chip bailey (no ag data)
oraldf1.by.ag <- subset(oraldf1, Site != "ChipBailey")


# first reorder levels high to low ag
oraldf1.by.ag$Site <- reorder(oraldf1.by.ag$Site, -oraldf1.by.ag$perc_ag)




#then add a rank column 
oraldf1.by.ag <- 
  transform(oraldf1.by.ag,agrank=as.numeric(factor(Site)))

# now plot chemsum by site
ggplot(oraldf1.by.ag, aes(x=Site,y=chemsum, color=Sample.type))+
  geom_point(size=8, shape=21, stroke=2)+
  geom_hline(yintercept = .2)

# and also with a log scaled Y axis
ggplot(oraldf1.by.ag, aes(x=agrank,y=chemsum, color=Sample.type))+
  geom_point(size=6, shape=21, stroke=2) + 
  geom_hline(yintercept = .2)+
  geom_smooth(method="lm")+
  scale_y_continuous(trans="log10")+
  annotate("text", x = 2.5, y = 100,
           label = "kendall correlation= -.0031, p= .96")



# then lets try a quick spearmans or kendall rank analysis
cor.test(oraldf1.by.ag$agrank,oraldf1.by.ag$chemsum, method = "kendall",exact=FALSE)







# plotting gis by pesticides ####
# on second thought I think we want this as some additive measure
# hmmmm.... 

ggplot(oralmelt, aes(x=perc_ag, y= oralhq, color=type))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~type)+
  scale_y_continuous(trans='log10')


ggplot(oralmelt, aes(x=perc_cornsoy, y= oralhq, color=type))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~type)+
  scale_y_continuous(trans='log10')

      
ggplot(oralmelt, aes(x=perc_apple, y= oralhq, color=type))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~type)+
  scale_y_continuous(trans='log10')



ggplot(oraldf1, aes(x=perc_ag, y= chemsum, color=Sample.type))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~Sample.type, scales = "free_y")+
  scale_y_continuous(trans='log10')

ggplot(oraldf1, aes(x=perc_cornsoy, y= chemsum, color=Sample.type))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~Sample.type, scales = "free_y")+
  scale_y_continuous(trans='log10')



ggplot(oraldf1, aes(x=perc_apple, y= chemsum, color=Sample.type))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~Sample.type, scales = "free_y")+
  scale_y_continuous(trans='log10')




ggplot(oraldf1, aes(x=perc_ag, y= chemsum))+
  geom_point()+
  geom_smooth()


# heat map test ####

test<- oraldf1 %>%
  group_by(Sample.type) %>% 
  summarise(mean_hq = mean(chemsum), max_hq = max(chemsum), perc_findings= mean(non_zero))%>%
  ungroup()

test2<- oraldf1 %>%
  group_by(Site) %>% 
  summarise(mean_hq = mean(chemsum), max_hq = max(chemsum), perc_findings= mean(non_zero))%>%
  ungroup()



bypesticide<- oralmelt.findings %>%
  group_by(Sample.type, pesticide) %>% 
  summarise(mean_hq = mean(oralhq), max_hq = max(oralhq), min_hq=min(oralhq))%>%
  ungroup()

bytype<- oralmelt.findings %>%
  group_by(Sample.type, type) %>% 
  summarise(mean_hq = mean(oralhq), max_hq = max(oralhq), min_hq=min(oralhq))%>%
  ungroup()



ggplot(bypesticide, aes(pesticide, Sample.type, fill= log(mean_hq))) + 
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_continuous(high = "#132B43", low = "#56B1F7")


ggplot(bypesticide, aes(pesticide, Sample.type, fill= log(max_hq))) + 
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_continuous(high = "#132B43", low = "#56B1F7")


ggplot(bytype, aes(type, Sample.type, fill= log(mean_hq))) + 
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_continuous(high = "#132B43", low = "#56B1F7")


ggplot(bytype, aes(type, Sample.type, fill= log(max_hq))) + 
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_continuous(high = "#132B43", low = "#56B1F7")



ggplot(contactdf1, aes(Site, Sample.type, fill= log(chemsum))) + 
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



# okay lets try to make a df showing percent of bees in a group (across sites) 
# that are over 3% oral ld50 for each pesticide



efsald50 <- oralmelt %>%
  group_by(Sample.type, pesticide) %>% 
  summarise(prop = (length(which(oralhq >= .03))/length(oralhq) ))%>%
  ungroup()

# this is showing the proportion of bees (aggregated across sites) that are over 3% of ld50 (oral)
ggplot(efsald50, aes(pesticide, prop, color=Sample.type))+
  geom_point()+
  facet_wrap(~Sample.type)+
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#heatmap of bees over 3% oral ld50
# which sites have which bees over 3% of ld50
ggplot(ld50.3percoral, aes(Site, Sample.type)) + 
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_continuous(high = "#132B43", low = "#56B1F7")






# Graphing/analysis ####
# if we log transform we can drop all zeros and actually see our data 

#first oralhq levels for bees across pesticide types
ggplot(oralmelt, aes(x=Sample.type, y=oralhq, color=type))+
  geom_boxplot(alpha=.5)+
  facet_wrap(~type)+
  scale_y_continuous(trans='log10')+
  geom_hline(yintercept = 1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(contactmelt, aes(x=Sample.type, y=contacthq, color=type))+
  geom_boxplot(alpha=.5)+
  facet_wrap(~type)+
  scale_y_continuous(trans='log10')+
  geom_hline(yintercept = 1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



# here additive oralhq across bees
# HB comes out significantly higher. everything else similiar
ggplot(oraldf1, aes(x=Sample.type, y=chemsum, color=Sample.type))+
  geom_boxplot(alpha=.5)+
  scale_y_continuous(trans='log10')+
  geom_hline(yintercept = .4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# breaks = scales::log_breaks(n = 8),labels = scales::comma   - changes scales to be non scientific notation, and sets number of breaks

ggplot(oraldf1, aes(x=Sample.type, y=chemsum, color=Sample.type))+
  geom_boxplot(lwd=1.2)+
  geom_point()+
  geom_hline(yintercept = .4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




oralaov <- aov(data=oraldf1, chemsum~Sample.type)
TukeyHSD(oralaov)

# and for contact. here HB and mason pollen provisions
# come out significantly higher than the rest
ggplot(contactdf1, aes(x=Sample.type, y=chemsum))+
  geom_boxplot(alpha=.5)+
  scale_y_continuous(trans='log10')+
  geom_hline(yintercept = .4)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

contactaov <- aov(data=contactdf1, chemsum~Sample.type)
TukeyHSD(contactaov)




#additive oral hq by percent corn/soy in 2015
ggplot(oraldf1, aes(x=perc_cornsoy, y=chemsum))+
  geom_point(alpha=.5)+
  geom_smooth(aes(color=Sample.type))+
  geom_hline(yintercept = 1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# and then by percent ag
ggplot(oraldf1, aes(x=perc_ag, y=chemsum))+
  geom_point(alpha=.5)+
  geom_smooth(aes(color=Sample.type))+
  geom_hline(yintercept = 1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(oralmelt, aes(x=Sample.type, y=oralhq, color=Sample.type))+
  geom_boxplot(alpha=.5)+
  facet_wrap(~pesticide, scales = "free")+
  scale_y_continuous(trans='log10')+
  geom_hline(yintercept = 1)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())




# then across just different main insecticide classes
ggplot(oralmelt, aes(x=Sample.type, y=oralhq, color=class))+
  geom_boxplot(alpha=.5)+
  facet_wrap(~class)+
  scale_y_continuous(trans='log10')+
  geom_hline(yintercept = 1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(oralmelt, aes(x=Sample.type, y=oralhq, color=class))+
  geom_boxplot(alpha=.5)+
  facet_wrap(~class)+
  geom_hline(yintercept = 1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# then aggregated
ggplot(oralmelt, aes(x=Sample.type, y=oralhq))+
  geom_boxplot(alpha=.5)+
  scale_y_continuous(trans='log10')+
  geom_hline(yintercept = 1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





ggplot(contactmelt, aes(x=Sample.type, y=contacthq, color=type))+
  geom_boxplot(alpha=.5)+
  facet_wrap(~type)+
  scale_y_continuous(trans='log10')+
  geom_hline(yintercept = 1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))






#### graphs for scott ####
# lets quickly for nice graphing add a shorter name for the bees
shortname <- read.csv("input/short_bee_name.csv")
shortname.contactdf1 <- merge(contactdf1,shortname, by="Sample.type")

dev.new()

ggplot(shortname.contactdf1, aes(x=bee.name, y=chemsum, color=Sample.type))+
  geom_violin(lwd=1.2)+
  geom_point()+
  scale_y_continuous(trans='log10')+
  geom_hline(yintercept = .4, color="darkred", size=1.2)+
  theme_bw()+
  ylab("log % Contact LD50")+
  annotate("text", x=5.7,y=1, label="EPA level of concern", size=8)+
  theme(text=element_text(size=35), 
        axis.text.x = element_text(angle = 70, vjust = 1, hjust=1,face="italic"),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin=unit(c(1,1,1.5,1.2),"cm"))



ggplot(shortname.contactdf1, aes(x=bee.name, y=chemsum, fill=Sample.type, alpha=.5))+
  geom_boxplot(lwd=1.2)+
  geom_point(size=5, pch=21, color="black", aes(fill=Sample.type))+
  geom_hline(yintercept = .4, color="darkred", size=1.2)+
  theme_bw()+
  ylab("% Contact LD50")+
  annotate("text", x=5.7,y=.43, label="EPA level of concern", size=8)+
  theme(text=element_text(size=35), 
        axis.text.x = element_text(angle = 70, vjust = 1, hjust=1,face="italic"),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin=unit(c(1,1,1.5,1.2),"cm"))











#NMDS  attempt  -------------------------------------
# now lets try to make an NMDS plot of pesticides by species



# okay below is code to give a unique name to each Sample (andrea 1, andrena2 etc)

df2<-contactdf # using all pesticides screened for. not sure if more informative?


df3 <- df2[order(df2[,3]),] # reorder df by sample type


df4<-df3 %>% # make a unique id column
  group_by(Sample.type) %>%
  mutate(id=row_number()) %>%
  ungroup()

df4$Sample.type<- paste(df4$Sample.type,df4$id,sep="_")

df4<-as.data.frame(df4) #make my tibble a df 


beedf<-df3[rowSums(df3[,4:(ncol(df3)-8)])>0,] #drop bees that have no pesticides to match nmds requirements
bee <- (beedf[,"Sample.type"]) # make vector of bees for later grouping

# now make new column and fill with nesting type of bee
beedf$beetype <- beedf$Sample.type

levels(beedf$beetype)[match(c("Andrena sp.","Andrena (Melandrena)","Bombus impatiens queen" ),levels(beedf$beetype))] <- "groundnesting"
levels(beedf$beetype)[match(c("Bombus impatiens worker","Apis mellifera worker"),levels(beedf$beetype))] <- "boxnesting"
levels(beedf$beetype)[match(c("Xylocopa virginica"),levels(beedf$beetype))] <- "woodnesting"

beetype <- (beedf[,"beetype"]) # make vector of nesting type for later grouping


bee <- as.vector(bee) # make a vector of the list of bees (useful later)
beetype <- as.vector(beetype) # and for beetype


df5<- df4[,3:(ncol(df3)-8)] # drop meta data

df6 <- df5[,-1]  # now rename rows to be my samples. df should have no descriptor data
rownames(df6) <- df5[,1]

df6<-df6[rowSums(df6)>0,] # drop 6 rows with no pesticides because it makes NMDS sad


# now format the data for a few different NMDS

nopesticides <- df5[rowSums(df5[,-1])<=0,]

df6sqrt<- sqrt(df6) # some people always square root their data to limit influence of extremes

# lets also look at presence absence data
presenceabsenceoral <-df6
presenceabsenceoral[presenceabsenceoral[,0:ncol(presenceabsenceoral)] >0] <- 1 # set all numbers over 0 to 1



nmds_results <- metaMDS(df6,  
                        distance = "bray",       
                        k = 3,
                        maxit = 999,
                        trymax = 500,
                        wascores = TRUE)


nmds_presence <- metaMDS(presenceabsenceoral,  
                        distance = "bray",       
                        k = 3,
                        maxit = 999,
                        trymax = 500,
                        wascores = TRUE)



nmds_sqrt <- metaMDS(df6sqrt,  
                        distance = "bray",       
                        k = 3,
                        maxit = 999,
                        trymax = 500,
                        wascores = TRUE)



beecolor= c("#F8766D","#B79F00", "#00BA38", "#00BFC4", "#619CFF","#F564E3")
beecolor_loop = c("#F8766D","#B79F00", "#00BA38", "#00BFC4", "#619CFF","#F564E3")[beedf$Sample.type]


# okay i think this is as pretty as my NMDS plot is going to get...

#first decide which nmds to plot
nmds<-nmds_results
#nmds <- nmds_presence
#nmds<-nmds_sqrt

dev.new() #function will open a new plot window, which then becomes the target for all plots.

# dev.off()

par(mar =  c(5, 6, 5, 5)) # set margins around the graph

plot(nmds,
     type="n",
     cex.lab=2.7,
     cex.axis=1.5) #create graph and make labels big

points(nmds, display = "sites", 
       pch=21,
       cex = 3,
      bg=beecolor_loop)


ordiellipse(nmds,groups= bee,
            draw="polygon",
            col = beecolor,
            alpha=70)


# you can also change how your ellipses are being drawn
# ordiellipse(nmds,groups= bee,
#             draw="polygon",
#             col = beecolor,
#             kind = "se",
#             conf = 0.95,
#             alpha=70)



legend(x="bottomleft",    # this adds a legend
       legend=levels(beedf$Sample.type),
       col=beecolor, 
       pch=21,
       pt.bg = beecolor,
       cex = 1.8)



# # then actually I think I want to only label certain drivers
# # first figure out top ten pesticides
# temp <- colSums(df6)
# sort(temp)
# topten <- temp[1:10]
# 
# temp<-df6[c("Dinotefuran","Thiabendazole","Thiamethoxam","Clothianidin","Imidacloprid",
#      "Acetamiprid","Tebuthiuron","Thiophanate.methyl","Pyrimethanil","Carbaryl")]
# 
# temp<-df6[c("Thiamethoxam")]
# 


# ARROWS
env<-envfit(nmds, df6, permutations=999)# make arrows using df your NMDS was created from

# add in arrows!
plot(env,
     col = "black", # arrow color
     cex = 2, # change thickness
     arrow.mul=4,# extends arrows
     p.max = 0.05 # only plot significant arrows
) 




# __________________________________________________





# ordispider(nmds,  # this adds spiders. ooo spooooky !!!
#            groups=bee,
#            col=c("red", "blue", "gold", "green","cyan", "black"))
# 


# # 
# ordihull(nmds,     # this adds polygons connecting points
#          groups=bee,
#          draw="polygon",
#          col=c("red", "blue", "gold", "green","cyan", "black"),
#          label=F)
# # 
text(nmds, display = "species", cex = 0.8, col = "darkred")     # adds pesticide locations and labels
# 
#      
#          
#    
text(nmds, display = "sites", cex = 0.8, col = "black")     # adds bee sample location and labels
# 


## NMDS nesting ####


plot(nmds)

ell1<-ordiellipse(nmds,groups= beetype, kind = "se", conf = 0.95, label=F,
                  col=c("red", "blue", "gold", "green","cyan", "black"))



legend(x="bottomleft",    # this adds a legend
       legend=levels(beedf$beetype),
       col=c("red", "blue", "gold", "green","cyan", "black"), 
       pch=21,
       pt.bg = c("red", "blue", "gold", "green","cyan", "black"))

env<-envfit(nmds, oraldf6, permutations=999)# make arrows
plot(env) # add in arrows!






## ggplot NMDS ####
#lets try to put this into ggplot for a prettier graph
data.scores = as.data.frame(scores(nmds)) # put nmds into dataframe
data.scores$bee = bee #add bees to datafram

en_coord_cont = as.data.frame(scores(env, "vectors")) * ordiArrowMul(env) # then extract arrows

# now plot everything
ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = bee), size = 3, alpha = 0.5) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) +
  labs(colour = "bee")+ 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30")+ 
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(en_coord_cont)) 

# oh wait ellipses! need those -- I have given up on this. I have no clue how vegancovellipse is working
# cant seem to pull out the data to plot it in ggplot


# NMDS.ell = data.frame(MDS1 = nmds$points[,1], MDS2 = nmds$points[,2],group=bee)
# 
# veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
# {
#   theta <- (0:npoints) * 2 * pi/npoints
#   Circle <- cbind(cos(theta), sin(theta))
#   t(center + scale * t(Circle %*% chol(cov)))
# }
# 
# df_ell <- data.frame()
# 
# for(g in levels(NMDS.ell$group)){
#   df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS.ell[NMDS.ell$group==g,],
#                                                    veganCovEllipse(ell1[[g]]$cov,ell1[[g]]$center,ell1[[g]]$scale)))
#                                 ,group=g))
# }




##no hb nmds ####
#first run the NMDS
nmdsnohb <- metaMDS(nohb,  
                    distance = "bray",       
                    k = 3,
                    maxit = 999,
                    trymax = 500,
                    wascores = TRUE)




# then create list of bees
nohbsim <- beesim[-c(38:53),]

beenohb <- (nohbsim[,"bee"])
beenohb <- as.vector(beenohb)
nohbsim$bee<-as.factor(nohbsim$bee)

plot(nmdsnohb)


ordiellipse(nmdsnohb,  # this adds ellipses
            groups=beenohb,
            col=c("red", "blue", "green","cyan", "tomato2"))


legend(x="bottomleft", 
       legend=levels(nohbsim$bee),
       col=c("red", "blue",  "green","cyan", "tomato2"), 
       pch=21,
       pt.bg = c("red", "blue",  "green","cyan", "tomato2"))



ordihull(nmdsnohb,
         groups=beenohb,
         draw="polygon",
         col=c("red", "blue",  "green","cyan", "tomato2"),
         label=F)
text(nmdsnohb, display = "species", cex = 0.8, col = "darkred")





### now I will also try a PCoA below ####

dist <- vegdist(oraldf6,  method = "bray") # make distance matrix

library(ape)
PCOA <- pcoa(dist) # make pcoa



biplot.pcoa(PCOA)



library(network)
beecolor <-as.color(bee, opacity = .9)





### NMDS TRASH BELOW ####

by.sampleoral <- oraldf1[order(oraldf1[,3]),] # reorder df by sample type

rownames(by.sampleoral) = seq(length=nrow(by.sampleoral)) # now fix row # after subsetting

test<- by.sampleoral[,4:29] # drop meta data



test1<-test[rowSums(test)>0,] # drop rows with no pesticides


bee <- as.vector(by.sampleoral[,3])



dis.matrix.oral <- vegdist(test1, method = "bray")
dis.matrix.oral <- as.matrix(dis.matrix.oral, labels=T)
  


nmds_results <- metaMDS(dis.matrix.oral,  
                        distance = "bray",       
                        k = 3,
                        maxit = 999,
                        trymax = 500,
                        wascores = TRUE)




plot(nmds_results, "sites")
orditorp(nmds_results, "sites")
ordihull(nmds_results, groups = bee, draw = "polygon")

test<-aggregate(. ~ Sample.type, oraldf1, mean)
test1 <- test[,-1]
rownames(test1) <- test[,1]
test2<- test1[,3:28] # drop meta data


nmds_results <- metaMDS(test2,  
                        distance = "bray",       
                        k = 2,
                        maxit = 999,
                        trymax = 500,
                        wascores = TRUE)



plot(nmds_results)
orditorp(nmds_results, "sites")
ordihull(nmds_results, groups = bee, draw = "polygon")




example_NMDS=metaMDS(community_matrix, # Our community-by-species matrix
                     k=2)


plot(example_NMDS)
orditorp(example_NMDS, "sites")
ordihull(example_NMDS, draw = "polygon")



# MANOVA attempt ####

# IDONT THINK WE shOULD TRUST ThE beLOW

# # and now I think the below is an MANOVA (if I didnt it correctly) run for each pesticide across bees
# responses <- paste( colnames(oraldf1)[4:29], collapse=",") # get the string of pesticide names
# myformula <- as.formula( paste0( "cbind(", responses , ")~ Sample.type" ) )
# beeman <- manova(myformula, data=oraldf1)
# 
# summary(beeman)
# summary.aov(beeman)
# 


# if I did this correctly ( i think i did..)
# then we see that there are significant differences of pesticide community between bee type
# but when looking at which pesticides are significantly different across bees
# we find just Thiamethoxam (***), and chlorpyrifos (*) come out significant

# HOWEVER manova is highly sensitive to alot of zeros. which my data is full of


# PERMANOVA attempt ####
## NOTE - currently the NMDS section must be run before this
# because... well.. bad coding.. and using some of those dfs I create there
# oof. laziness

beematrix <- as.matrix(df5[,2:(ncol(df5))])
bee.dis <- vegdist(beematrix, method = 'bray')

beesimple<- df3[3:29]

# and now lets try the above after removing bees with no pesticides. it makes things sad
beematrix2 <- as.matrix(df6)
bee.dis.nozero <- vegdist(beematrix2, method = 'bray')


beesim <- cbind(df6,bee)

beeperman<- adonis2(bee.dis.nozero~bee, data=beesimple, permutations=999, method="bray")

beeperman
# this shows that there is a strong impact of bee species on the pesticides that are in each species



# now something additional that some people do but im not sure about
# is taking the square root of data to minimize the influence of super high numbers

# lets give it a quick test and see if it changes things drastically
testmat <- sqrt(beematrix2)

bee.dis.nozero <- vegdist(testmat, method = 'bray')

beeperman<- adonis2(bee.dis.nozero~bee, data=beesimple, permutations=999, method="bray")

# still significant.
# now im curious if I remove HB, what will happen
# assume no significance between the rest(?)

nohb <- df6[-c(38:53),]
nohbsim <- beesim[-c(38:53),]

nohb.matrix <- as.matrix(nohb)
nohb.dis.nozero <- vegdist(nohb.matrix, method = 'bray')
nohbperma<- adonis2(nohb.dis.nozero~bee, data=nohbsim, permutations=999, method="bray")

nohbperma
# still significant but much less so. my guess was wrong. fascinating


# I think the only way I can do pairsise tests is if i do new dfs
# with just two species and run permanovas
















# attempts at hurdle model #### ----------------------------------

# now I think we want to try a gamma hurdle model
# which first looks at if there are pesticides (0 or 1)
# and then the amounts after excluding the zero portion of the data set

# for the moment im mostly following 
# https://seananderson.ca/2014/05/18/gamma-hurdle/



# first make a new column 1 if pesticide, 0 if no pesticides found
oralmelt$non_zero <- ifelse(oralmelt$oralhq > 0, 1, 0)


#now lets graph frequency of finding pesticides by sample type
table<- with(oralmelt, table(non_zero, Sample.type))

ggplot(as.data.frame(table), aes(factor(Sample.type), Freq, fill = non_zero)) +     
  geom_col(position="fill")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# then make a new df with no zeros
oralmelt1<-subset(oralmelt, non_zero == 1)


# now produce two glms, one looking at the change of zeros
# the other subsetting the data to only 1 and using gamma distribution
hist(log(oralmelt1$oralhq),breaks=20)
hist(oralmelt1$oralhq, breaks=15)


qqnorm(oralmelt1$oralhq)
qqnorm(log(oralmelt1$oralhq))

# a log transformation doesnt quite normalize it 
# still right skewed

m1 <- glm(non_zero ~ Sample.type, data = oralmelt, family = binomial(link = logit))
#m2 <- glm(oralhq ~ type, data = oralmelt1, family = Gamma(link = log))


m2 <- glm(oralhq ~ Sample.type, data=oralmelt1, family=inverse.gaussian())

emmeans(m1, pairwise ~ Sample.type, adjust="tukey")
emmeans(m2, pairwise ~ Sample.type, adjust="tukey")



# okay if I use gamma(link="log) 
# im getting an Error in glm.fit(x = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  : 
# NA/NaN/Inf in 'x'

#but both say theres no NA or INF... hmmmm...
any(is.na(oralmelt$non_zero))

any(is.infinite(oralmelt$non_zero))

d = rstandard(m2)

qqnorm(d)
qqline(d) 





#trash ####

# so these graphs look way different than christinas
# and after some poking around i figured out why
# shes graphing cumulative risk, that is, all hq added for each sample


# now looking at additive risk

#oral additive
ggplot(oral.add, aes(x=Sample.type, y=add.hq))+
  geom_boxplot(alpha=.5)+
  ylab("sum (oral hq)")+
  geom_hline(yintercept = 1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# contact additive
ggplot(contact.add, aes(x=Sample.type, y=add.hq))+
  geom_boxplot(alpha=.5)+
  ylab("sum(contact hq)")+
  geom_hline(yintercept = 1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# comparison residue levels in bees subset by chemical type
fungicide.aov <- aov(data=subset(bees.oral.melt, type=="Fungicide"),
                     oralhq ~ Sample.type)

TukeyHSD(fungicide.aov)


insecticide.aov <- aov(data=subset(bees.oral.melt, type=="Insecticide"),
                       oralhq ~ Sample.type)

TukeyHSD(insecticide.aov)


insecticide.aov <- aov(data=subset(bees.oral.melt, type=="Herbivore"),
                       oralhq ~ Sample.type)

TukeyHSD(insecticide.aov)

# andrena seem to have slightly elevated fungicide levels compared to a few other bees



# lets do a quick comparisons of bee species levels by class

neonic.lm <- lm(data=subset(bees.oral.melt, class=="neonicotinoid"),
                oralhq ~ Sample.type)

emmeans(neonic.lm, pairwise ~ Sample.type, adjust="tukey")

neonic.aov <- aov(data=subset(bees.oral.melt, class=="neonicotinoid"),
                  oralhq ~ Sample.type)

TukeyHSD(neonic.aov)


op.aov <- aov(data=subset(bees.oral.melt, class=="organophosphate"),
              oralhq ~ Sample.type)

TukeyHSD(op.aov)


# A. mellifera is coming out significantly higher (by alot!)
# than all other bees. 
# all other bee comparisons are non significant
# !! very interesting



