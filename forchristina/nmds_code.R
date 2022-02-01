# cleaned script for data prep, NMDS, and permanova
# put together for Christine
# by tobias -- 10/2021



setwd("C:/Users/obiew/Desktop/scott_pesticide_data") # set your working directory

# load packages 
library(ggplot2)
library(reshape2)
library(dplyr)
library(vegan)


rm(list = ls()) # clean your environment

# DATA PREP ---- Run from here to "END OF DATA PREP" to get your df created
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
rownames(alldf1) = seq(length=nrow(alldf1)) # now fix row numbers after subsetting




#### adjusting to HQ --------------------------------------------
# now we'll try to make everything in alldf (ppb) as a percent of ld50 of honey bees
# which would entail dividing alldf by corresponding values of ld50

#drop identifier colummns
justpesticides <- alldf1[7:98]


ld50contact <- ld50[1,2:93]
ld50oral <- ld50[2,2:93]

# okay probably not the cleanest way but we can add the ld50 row to the pesticide df

# add ld50s to pesticide df
oral<-rbind(justpesticides,ld50oral)

#fix row numbers (got wonky again)
rownames(oral) = seq(length=nrow(oral)) # now fix row # after subsetting

# then matrix it and divide by the ld 50 row
temp = as.matrix(oral)

oralhq<- t(t(temp)/temp[104,])

oralhq = as.data.frame(oralhq)

#now remove the ld 50 row
oralhq <- oralhq[1:103,]

# then we'll set all the inf and NA to 0 (where the ld50 was either 0 or NA)
oralhq[] <- lapply(oralhq, function(x) replace(x, !is.finite(x), 0))



#now do the same for contact (same code below, see above for annotations)
contact<-rbind(justpesticides,ld50contact)

rownames(contact) = seq(length=nrow(contact)) 

temp = as.matrix(contact)

contacthq<- t(t(temp)/temp[104,])

contacthq = as.data.frame(contacthq)


contacthq <- contacthq[1:103,]

contacthq[] <- lapply(contacthq, function(x) replace(x, !is.finite(x), 0))




#now add back descriptor columns
oraldf<- cbind(alldf1[,1:6],oralhq)
contactdf<- cbind(alldf1[,1:6],contacthq)

#then on second thought lets drop some of those descriptor columns
oraldf<- cbind(alldf1[,c(3,4,6)],oralhq)
contactdf<- cbind(alldf1[,c(3,4,6)],contacthq)



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



#### remove undetected pesticides ------------------------------

# I think we should remove pesticides that were never detected
# or that have no ld50
# as theyre not informative in our ld50 data table
# and then just report which pesticides were removed in a table or something

# we can do this by summing columns and removing zeros 

i <- (colSums((oraldf[,4:95]), na.rm=T) != 0) # T if colSum is not 0, F otherwise
cols_to_keep_oral <- c(rep(TRUE, 3), i) # then add 3 T to keep meta data columns


oraldf1 <- oraldf[, cols_to_keep_oral] # all the non-zero columns plus meta data
notfoundpesticides.oral <- oraldf[, !cols_to_keep_oral]  # all the zero columns


# and the same for contact
i <- (colSums((contactdf[,4:95]), na.rm=T) != 0) # T if colSum is not 0, F otherwise
cols_to_keep_contact <- c(rep(TRUE, 3), i) # then add 3 T to keep meta data columns


contactdf1 <- contactdf[, cols_to_keep_contact] # all the non-zero columns plus meta data
notfoundpesticides.contact <- contactdf[, !cols_to_keep_contact]  # all the zero columns

# the below should sum to 97 if everything went well
# ncol(oraldf1) + ncol(notfoundpesticides.oral)
# ncol(contactdf1) + ncol(notfoundpesticides.contact)

# here is the list of pesticides that were not found or have no oral ld50
noimpactoral <- colnames(notfoundpesticides.oral)
noimpactcontact <- colnames(notfoundpesticides.contact)






#### melting for graphing ----------------------------------
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



#### adding GIS information ####

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
#### clean environment ####

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
            "cols_to_keep_oral",
            "cols_to_keep_contact",
            "noimpactoral",
            "justpesticides"
))






# END OF DATA PREP ####




#NMDS  attempt  -------------------------------------
# now lets try to make an NMDS plot of pesticides by species



# okay below is code to give a unique name to each Sample (andrea1, andrena2 etc)
oraldf2<-oraldf1 # make new df

oraldf3 <- oraldf2[order(oraldf2[,3]),] # reorder df by sample type

oraldf3 <- oraldf3[rowSums(oraldf3[,4:29])>0,] # drop all rows with no pesticides because it makes NMDS sad

oraldf4<-oraldf3 %>% # make a unique id column
  group_by(Sample.type) %>%
  mutate(id=row_number()) %>%
  ungroup()

oraldf4$Sample.type<- paste(oraldf4$Sample.type,oraldf4$id,sep="_")

oraldf4<-as.data.frame(oraldf4) #make my tibble back into a df 

# make a vector of the list of bees (useful later)
# for specifying groups in graphing
bee <- (oraldf3[,"Sample.type"])
bee <- as.vector(bee) 

# drop meta data so the df is only pesticide levels and samples names
oraldf5<- oraldf4[,3:29] # drop meta data

oraldf6 <- oraldf5[,-1]  # now rename rows to be my samples
rownames(oraldf6) <- oraldf5[,1]

# if you want to reduce impact of outliers you can sqrt your data
oraldf6sqrt<- sqrt(oraldf6)

# lets also look at presence absence data
presenceabsenceoral <-oraldf6
presenceabsenceoral[presenceabsenceoral[,0:ncol(presenceabsenceoral)] >0] <- 1 # set all numbers over 0 to 1


# new i'll create the 3 diferent NMDS

nmds_results <- metaMDS(oraldf6,  
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



nmds_sqrt <- metaMDS(oraldf6sqrt,  
                     distance = "bray",       
                     k = 3,
                     maxit = 999,
                     trymax = 500,
                     wascores = TRUE)




#first decide which nmds to plot
nmds<-nmds_results
#nmds <- nmds_presence
#nmds<-nmds_sqrt


# the way NMDS gets plotted (that I could figure out)
# is you plot the base NMDS and then add layers ontop of it
# for example, ellipses, legend, text etc


plot(nmds)

# now create ellipses. I dont know whats best. paper I read 
# used SE 95% conf so ill do that too. << read more on this

ordiellipse(nmds,groups= bee, kind = "se", conf = 0.95, label=F,
                  col=c("red", "blue", "gold", "green","cyan", "black"))



legend(x="bottomleft",    # this adds a legend
       legend=levels(oraldf3$Sample.type),
       col=c("red", "blue", "gold", "green","cyan", "black"), 
       pch=21,
       pt.bg = c("red", "blue", "gold", "green","cyan", "black"))

env<-envfit(nmds, oraldf6, permutations=999)# make arrows
plot(env) # add in arrows!

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
# 
#text(nmds, display = "species", cex = 0.8, col = "darkred")     # adds pesticide labels

#text(nmds, display = "sites", cex = 0.8, col = "black")     # adds bee sample labels






#### ggplot NMDS ####
# below is my attempt to plot the nmds in ggplot2
# I couldnt quite figure out the ellipses so I gave up

#lets try to put this into ggplot for a prettier graph
data.scores = as.data.frame(scores(nmds)) # put nmds into dataframe
data.scores$bee = bee #add bees to dataframe

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

# oh wait ellipses! need those -- I have given up on this. I have no cluehow vegancovellipse is working
# cant seem to pull out the data to plot it in ggplot. below is my attempt


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





# Permanova code ####


# make df of only sample groups and pesticides
beesimple<- oraldf3[3:29]

# and same with just pesticides
beesimple.justpest <- beesimple[,-1]

# turn it into a matrix
beematrix <- as.matrix(beesimple.justpest)

#calculate dis.matrix
bee.dis.nozero <- vegdist(beematrix, method = 'bray') 

#run adonis2 (permanova)
beeperman<- adonis2(bee.dis.nozero~Sample.type, data=beesimple, permutations=999, method="bray")

beeperman









