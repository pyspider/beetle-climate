###################
#Fig.S2: NMDS for OTU table
###################
library(phyloseq) #phyloseq v1.30.0
library(ggplot2) #ggplot2 v3.3.4
library(plyr) #plyr v1.8.6
library(vegan) #vegan v2.5-6
library(anytime) #anytime v0.3.9
library(dplyr) #dplyr v1.0.2
library(ggordiplots) #ggordiplots v0.3.0
library(reshape) #reshape v0.8.8
library(patchwork) #patchwork v1.0.0
library(ggfortify) #ggfortify v0.4.11
library(cowplot) #cowplot v1.1.1
library(stringr) #stringr v1.4.0

setwd(HEREISYOURWORKINGDIRECTORY)

comm <- read.table("Data/GLGS_Community_Data_fixed_three_cases_20201019.txt",header=TRUE,sep = "\t",stringsAsFactors = FALSE)
data_biweek<-cbind(comm$Channel_ID,comm$Collection_Date_formatted_1,comm$Sumaclust_97_OTU,sub(" ","",paste0(ceiling(lubridate::week(as.Date(anydate(comm$Collection_Date_formatted_1))) / 2),format(as.Date(anydate(comm$Collection_Date_formatted_1)), "_%Y"))))
colnames(data_biweek)<-c("Channel_ID","Collection_Date_formatted_1","Sumaclust_97_OTU","Biweek_Number")
data_biweek<-as.data.frame(data_biweek)
biweek <- data_biweek[c('Biweek_Number','Channel_ID','Sumaclust_97_OTU','Collection_Date_formatted_1')]
#need to ignore the blank cell for the sample of failure finally.
biweek <- filter(biweek,biweek$Sumaclust_97_OTU!='')
biweek<-cbind(biweek,paste0(biweek$Channel_ID,"_",biweek$Biweek_Number))
colnames(biweek)<-c("Biweek_Number","Channel_ID","Sumaclust_97_OTU","Collection_Date_formatted_1","Combination_Number")

biweek_1 <- biweek[c('Combination_Number','Sumaclust_97_OTU')]
OTU.data_biweek<-table(biweek_1)

##create the table of OTU size
biweek_1<-as.data.frame(OTU.data_biweek)
OTU.data_biweek.pivot<-cast(biweek_1, Sumaclust_97_OTU ~ Combination_Number)
OTU.data_biweek.pivot<-apply(OTU.data_biweek.pivot,2,sum)
OTU.data_biweek.pivot<-as.data.frame(OTU.data_biweek.pivot)
OTU.data_biweek.pivot<-cbind(OTU.data_biweek.pivot,rownames(OTU.data_biweek.pivot))
biweek_size<-OTU.data_biweek.pivot
colnames(biweek_size)<-c("Freq","Combination_Number")

biweek$Collection_Date_formatted_1 <- lubridate::dmy(biweek$Collection_Date_formatted_1)
x<-dplyr::arrange(biweek, Collection_Date_formatted_1)

date.min<-data.frame(Biweek_Number="",min_date="")
date.min<-date.min[-1,]
date.list<-c()
for (i in unique(x$Biweek_Number)){ 
  
  date.list<-filter(x,x$Biweek_Number==i)
  date.min_0<-data.frame(Biweek_Number=i,min_date=min(as.Date(date.list$Collection_Date_formatted_1)))
  date.min<-rbind(date.min,date.min_0)
  date.list<-c() 
}
month <- date.min
data_month<-cbind(date.min,month=format(as.Date(anydate(date.min$min_date)), "%Y_%m"))

Biweek_Month_merge<-dplyr::inner_join(biweek,data_month,by=c("Biweek_Number"="Biweek_Number"))

#NMDS & PCA for OTUs of biweek
OTU.data_biweek.hel <- decostand(OTU.data_biweek, method = 'hellinger')
OTU.data_biweek.pivot<-cast(as.data.frame(OTU.data_biweek.hel), Combination_Number ~ Sumaclust_97_OTU)
test_11<-as.data.frame(OTU.data_biweek.hel)
#deleted the first column including biweek numbers.
row.names(OTU.data_biweek.pivot)<-OTU.data_biweek.pivot[,1]
OTU.data_biweek.pivot_fixed<-as.data.frame(OTU.data_biweek.pivot[,-1:-2])

OTU_biweek.nmds2 <- metaMDS(OTU.data_biweek.pivot_fixed, dist="bray", k=4)

#testing using Sørensen dissimilarity, binary=TRUE
# OTU.data_biweek_table<-as.data.frame(OTU.data_biweek)
# OTU.data_biweek_table[which(OTU.data_biweek_table$Freq>1),]$Freq<-1
# OTU.data_biweek_table_1<-cast(OTU.data_biweek_table, Combination_Number ~ Sumaclust_97_OTU)
# row.names(OTU.data_biweek_table_1)<-OTU.data_biweek_table_1[,1]
# OTU.data_biweek_table_1<-as.data.frame(OTU.data_biweek_table_1[,-1:-2])
# 
# OTU_biweek.nmds2 <- metaMDS(OTU.data_biweek_table_1, dist="bray", k=4)

#####################
#NMDS after the removal of the outilers
#####################
OTU_biweek.NMDS = data.frame(NMDS1 = OTU_biweek.nmds2$points[,1], NMDS2 = OTU_biweek.nmds2$points[,2])
OTU_biweek.NMDS$sample <- rownames(OTU_biweek.NMDS)

OTU_biweek.NMDS_merge<-dplyr::inner_join(OTU_biweek.NMDS,biweek_size,by=c("sample"="Combination_Number"))

OTU_biweek.NMDS_merge<-dplyr::inner_join(OTU_biweek.NMDS_merge,unique(Biweek_Month_merge[,-1:-4]),by=c("sample"="Combination_Number"))

OTU_biweek.NMDS$group <- as.factor(OTU_biweek.NMDS_merge$month)

#NMDS subgroup
OTU_biweek.NMDS_for_doug_Plot<-OTU_biweek.NMDS
OTU_biweek.NMDS_for_doug_Plot$month<-gsub("^[0-9]+_","",OTU_biweek.NMDS_for_doug_Plot$group,perl=TRUE)
OTU_biweek.NMDS_for_doug_Plot$elevation<-gsub("[1-4]2M?$","",gsub("_[0-9]+_[0-9]+$","",gsub("DLJ_","",OTU_biweek.NMDS_for_doug_Plot$sample,perl=TRUE),perl=TRUE),perl=TRUE)
OTU_biweek.NMDS_for_doug_Plot$Season<-gsub("^[0-9]+$","Summer",gsub("^11$|^12$|^01$|^02$|^03$","Winter",gsub("[0-9]+_","",OTU_biweek.NMDS_for_doug_Plot$group,perl=TRUE),perl=TRUE),perl=TRUE)
OTU_biweek.NMDS_for_doug_Plot$Location<-gsub("^201[45]$","South",gsub("^2016$","North",gsub("^\\S+[0-9]+_","",gsub("^\\S+17_2015$|^\\S+18_2015$|^\\S+19_2015$|^\\S+20_2015$|^\\S+21_2015$|^\\S+22_2015$|^\\S+23_2015$|^\\S+24_2015$|^\\S+25_2015$|^\\S+26_2015$","North",OTU_biweek.NMDS_for_doug_Plot$sample,perl=TRUE),perl=TRUE),perl=TRUE),perl=TRUE)

##from Zhengyang Fig4.R##
OTU_NMDS <- OTU_biweek.NMDS_for_doug_Plot

# definition of season needs to be adjusted
OTU_NMDS$Season[OTU_NMDS$month %in% c("03","04","05")] <- "Spring"
OTU_NMDS$Season[OTU_NMDS$month %in% c("06","07","08")] <- "Summer"
OTU_NMDS$Season[OTU_NMDS$month %in% c("09","10","11")] <- "Autumn"
OTU_NMDS$Season[OTU_NMDS$month %in% c("12","01","02")] <- "Winter"
OTU_NMDS$Location_and_elevation <- paste(OTU_NMDS$Location,OTU_NMDS$elevation)

outlier_0<-c("DLJ_140012_23_2015","DLJ_180022_8_2016","DLJ_180022_24_2015","DLJ_140012_5_2016","DLJ_140012_7_2016")
outlier_list_1400<-c("140032_2_2015","DLJ_140012_14_2016","DLJ_140022_26_2015","DLJ_140012M_26_2015")
outlier_list_1800<-c("DLJ_180022_25_2015","180012_3_2015")

outlier_list<-c(outlier_0,outlier_list_1400,outlier_list_1800)

outlier_list_merge<-dplyr::inner_join(Biweek_Month_merge,OTU_NMDS[which(OTU_NMDS$sample %in% outlier_list),],by=c("Combination_Number"="sample"))
##These outlier points are rare species.Then you need to remove them.

#################
#from NMDS results after the removal of all outlier points. 
#################
#remove outlier points

`%!in%` = Negate(`%in%`)
OTU.data_biweek.pivot_fixed_filtered<-OTU.data_biweek.pivot_fixed[which(row.names(OTU.data_biweek.pivot_fixed) %!in% outlier_list),]

#for Sørensen dissimilarity
# `%!in%` = Negate(`%in%`)
# OTU.data_biweek_table_1<-OTU.data_biweek_table_1[which(row.names(OTU.data_biweek_table_1) %!in% outlier_list),]
# dist<-vegdist(OTU.data_biweek_table_1,method="bray", binary=TRUE,diag=TRUE)
# OTU_biweek.nmds2<-metaMDS(dist)

OTU_biweek.nmds2 <- metaMDS(OTU.data_biweek.pivot_fixed_filtered, dist="bray", k=4)

OTU_biweek.NMDS = data.frame(NMDS1 = OTU_biweek.nmds2$points[,1], NMDS2 = OTU_biweek.nmds2$points[,2])
OTU_biweek.NMDS$sample <- rownames(OTU_biweek.NMDS)

OTU_biweek.NMDS_merge<-dplyr::inner_join(OTU_biweek.NMDS,biweek_size,by=c("sample"="Combination_Number"))

OTU_biweek.NMDS_merge<-dplyr::inner_join(OTU_biweek.NMDS_merge,unique(Biweek_Month_merge[,-1:-4]),by=c("sample"="Combination_Number"))

OTU_biweek.NMDS$group <- as.factor(OTU_biweek.NMDS_merge$month)

#NMDS subgroup
OTU_biweek.NMDS_for_doug_Plot<-OTU_biweek.NMDS
OTU_biweek.NMDS_for_doug_Plot$month<-gsub("^[0-9]+_","",OTU_biweek.NMDS_for_doug_Plot$group,perl=TRUE)
OTU_biweek.NMDS_for_doug_Plot$elevation<-gsub("[1-4]2M?$","",gsub("_[0-9]+_[0-9]+$","",gsub("DLJ_","",OTU_biweek.NMDS_for_doug_Plot$sample,perl=TRUE),perl=TRUE),perl=TRUE)
OTU_biweek.NMDS_for_doug_Plot$Season<-gsub("^[0-9]+$","Summer",gsub("^11$|^12$|^01$|^02$|^03$","Winter",gsub("[0-9]+_","",OTU_biweek.NMDS_for_doug_Plot$group,perl=TRUE),perl=TRUE),perl=TRUE)
OTU_biweek.NMDS_for_doug_Plot$Location<-gsub("^201[45]$","South",gsub("^2016$","North",gsub("^\\S+[0-9]+_","",gsub("^\\S+17_2015$|^\\S+18_2015$|^\\S+19_2015$|^\\S+20_2015$|^\\S+21_2015$|^\\S+22_2015$|^\\S+23_2015$|^\\S+24_2015$|^\\S+25_2015$|^\\S+26_2015$","North",OTU_biweek.NMDS_for_doug_Plot$sample,perl=TRUE),perl=TRUE),perl=TRUE),perl=TRUE)

OTU_NMDS <- OTU_biweek.NMDS_for_doug_Plot

# definition of season needs to be adjusted
OTU_NMDS$Season[OTU_NMDS$month %in% c("03","04","05")] <- "Spring"
OTU_NMDS$Season[OTU_NMDS$month %in% c("06","07","08")] <- "Summer"
OTU_NMDS$Season[OTU_NMDS$month %in% c("09","10","11")] <- "Autumn"
OTU_NMDS$Season[OTU_NMDS$month %in% c("12","01","02")] <- "Winter"
OTU_NMDS$Location_and_elevation <- paste(OTU_NMDS$Location,OTU_NMDS$elevation)

####################### Now start plotting

#my.colors <- c("gray","gray","#8a5656","#618dc6")
my.colors <- c("gray","gray","#1d953f","black")

#####################################
# Now let's plot different elevation seperately
#####################################

OTU.season.nmds.hul <- gg_ordiplot(OTU_biweek.nmds2, groups =OTU_NMDS$Season, show.groups = c("Summer","Winter") , hull = TRUE, spiders = TRUE,label = TRUE, ellipse = FALSE, plot = FALSE,conf=0.95,kind = "se")

OTU.location.nmds.spider <- gg_ordiplot(OTU_biweek.nmds2, groups =OTU_NMDS$Location, hull = TRUE, spiders = TRUE,label = TRUE, ellipse = FALSE, plot = FALSE,conf=0.95,kind = "se")

OTU_NMDS_1400m <- ggplot(OTU_NMDS,aes(x=NMDS1, y=NMDS2))+
  geom_point(data = OTU_NMDS[(OTU_NMDS$Season=="Summer"|OTU_NMDS$Season=="Winter")&OTU_NMDS$elevation==1400,],aes(color=Season,shape=Location_and_elevation), size =2)+
  scale_shape_manual(values=c(19,1,15,0))+
  scale_color_manual(values = alpha(my.colors,0.5))+ geom_segment(data=OTU.location.nmds.spider$df_spiders[which(rownames(OTU.location.nmds.spider$df_spiders) %in% rownames(OTU_NMDS[(OTU_NMDS$Season=="Summer"|OTU_NMDS$Season=="Winter")&OTU_NMDS$elevation==1400,])),][which(OTU.location.nmds.spider$df_spiders[which(rownames(OTU.location.nmds.spider$df_spiders) %in% rownames(OTU_NMDS[(OTU_NMDS$Season=="Summer"|OTU_NMDS$Season=="Winter")&OTU_NMDS$elevation==1400,])),]$Group=="North"),], aes(x=cntr.x, xend=x, y=cntr.y, yend=y,color=Group,linetype=1),alpha=0.4,show.legend=FALSE) +geom_segment(data=OTU.location.nmds.spider$df_spiders[which(rownames(OTU.location.nmds.spider$df_spiders) %in% rownames(OTU_NMDS[(OTU_NMDS$Season=="Summer"|OTU_NMDS$Season=="Winter")&OTU_NMDS$elevation==1400,])),][which(OTU.location.nmds.spider$df_spiders[which(rownames(OTU.location.nmds.spider$df_spiders) %in% rownames(OTU_NMDS[(OTU_NMDS$Season=="Summer"|OTU_NMDS$Season=="Winter")&OTU_NMDS$elevation==1400,])),]$Group=="South"),], aes(x=cntr.x, xend=x, y=cntr.y, yend=y,color=Group,linetype=2),alpha=0.4,show.legend=FALSE)+ scale_linetype_identity()+labs(x = paste("NMDS1"), y = paste("NMDS2")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  labs(title = paste("OTU no phylogenetic information (1400m)"))+annotate("text",x=-0.9,y=-1.3,label=paste('Stress =', round(OTU_biweek.nmds2$stress, 4)),size=3)+ theme(legend.position = "none")

OTU_NMDS_1400m


OTU_NMDS_1800m <- ggplot(OTU_NMDS,aes(x=NMDS1, y=NMDS2))+
  geom_point(data = OTU_NMDS[(OTU_NMDS$Season=="Summer"|OTU_NMDS$Season=="Winter")&OTU_NMDS$elevation==1800,],aes(color=Season,shape=Location), size =2)+
  scale_shape_manual(values=c(19,1,15,0))+
  scale_color_manual(values = alpha(my.colors,0.5))+ geom_segment(data=OTU.location.nmds.spider$df_spiders[which(rownames(OTU.location.nmds.spider$df_spiders) %in% rownames(OTU_NMDS[(OTU_NMDS$Season=="Summer"|OTU_NMDS$Season=="Winter")&OTU_NMDS$elevation==1800,])),][which(OTU.location.nmds.spider$df_spiders[which(rownames(OTU.location.nmds.spider$df_spiders) %in% rownames(OTU_NMDS[(OTU_NMDS$Season=="Summer"|OTU_NMDS$Season=="Winter")&OTU_NMDS$elevation==1800,])),]$Group=="North"),], aes(x=cntr.x, xend=x, y=cntr.y, yend=y,color=Group,linetype=1),alpha=0.4,show.legend=FALSE) +geom_segment(data=OTU.location.nmds.spider$df_spiders[which(rownames(OTU.location.nmds.spider$df_spiders) %in% rownames(OTU_NMDS[(OTU_NMDS$Season=="Summer"|OTU_NMDS$Season=="Winter")&OTU_NMDS$elevation==1800,])),][which(OTU.location.nmds.spider$df_spiders[which(rownames(OTU.location.nmds.spider$df_spiders) %in% rownames(OTU_NMDS[(OTU_NMDS$Season=="Summer"|OTU_NMDS$Season=="Winter")&OTU_NMDS$elevation==1800,])),]$Group=="South"),], aes(x=cntr.x, xend=x, y=cntr.y, yend=y,color=Group,linetype=2),alpha=0.4,show.legend=FALSE)+ scale_linetype_identity()+ labs(x = paste("NMDS1"), y = paste("NMDS2")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  labs(title = paste("OTU no phylogenetic information (1800m)"))+annotate("text",x=-1.1,y=-1.6,label=paste('Stress =', round(OTU_biweek.nmds2$stress, 4)),size=3)+
  theme(legend.position = "none")# +  coord_cartesian(xlim=c(-0.25,-0.175),ylim=c(0.125,0.2)) 

OTU_NMDS_1800m

#########################
##using Unifrac on ultrametric denovo tree
#########################
library(phyloseq) #phyloseq v1.30.0
library(ggplot2) #ggplot2 v3.3.4
library(plyr) #plyr v1.8.6
library(vegan) #vegan v2.5-6
library(anytime) #anytime v0.3.9
library(dplyr) #dplyr v1.0.2
library(ggordiplots) #ggordiplots v0.3.0
library(reshape) #reshape v0.8.8
library(patchwork) #patchwork v1.0.0
library(ggtree) #ggtree v2.0.4
library(stringr) #stringr v1.4.0

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

##create size group by biweek and malaise trap number.
comm <- read.table("Data/GLGS_Community_Data_fixed_three_cases_20201019.txt",header=TRUE,sep = "\t",stringsAsFactors = FALSE)
biweek<-cbind(comm$OTU,comm$Channel_ID,comm$Collection_Date_formatted_1,sub(" ","",paste0(ceiling(lubridate::week(as.Date(anydate(comm$Collection_Date_formatted_1))) / 2),format(as.Date(anydate(comm$Collection_Date_formatted_1)), "_%Y"))))
colnames(biweek)<-c("OTU","Channel_ID","Collection_Date_formatted_1","Biweek_Number")
biweek<-as.data.frame(biweek)
biweek<-cbind(biweek,paste0(biweek$Channel_ID,"_",biweek$Biweek_Number))
colnames(biweek)<-c("OTU","Channel_ID","Collection_Date_formatted_1","Biweek_Number","Combination_Number")
biweek<-as.data.frame(biweek)
biweek$Collection_Date_formatted_1 <- lubridate::dmy(biweek$Collection_Date_formatted_1)
x<-dplyr::arrange(biweek, Collection_Date_formatted_1)

date.min<-data.frame(Biweek_Number="",min_date="")
date.min<-date.min[-1,]
date.list<-c()
for (i in unique(x$Biweek_Number)){ 
  
  date.list<-filter(x,x$Biweek_Number==i)
  date.min_0<-data.frame(Biweek_Number=i,min_date=min(as.Date(date.list$Collection_Date_formatted_1)))
  date.min<-rbind(date.min,date.min_0)
  date.list<-c() 
}
month <- date.min
data_month<-cbind(date.min,month=format(as.Date(anydate(date.min$min_date)), "%Y_%m"))

Biweek_Month_merge<-dplyr::inner_join(biweek,data_month,by=c("Biweek_Number"="Biweek_Number"))

#need to ignore the blank cell for the sample of failure finally.
data_biweek <- filter(Biweek_Month_merge,Biweek_Month_merge$OTU!='')
data_biweek_1<-data_biweek[,-2:-4]
data_biweek_1<-data_biweek_1[,-3:-4]
data_biweek_1<-as.data.frame(table(data_biweek_1))
OTU.data_biweek.pivot<-cast(data_biweek_1, OTU ~ Combination_Number)
biweek_OTU<-OTU.data_biweek.pivot[-1,]
biweek_OTU<-biweek_OTU[,-2]
rownames(biweek_OTU)<-biweek_OTU[,1]
biweek_OTU<-biweek_OTU[,-1]
biweek_OTU<-data.frame(biweek_OTU)
colnames(biweek_OTU)<-gsub("^X","",colnames(biweek_OTU),perl=TRUE)
rownames(biweek_OTU)<-gsub("\\|","_",rownames(biweek_OTU),perl=TRUE)
rownames(biweek_OTU)<-gsub("\\:","___",rownames(biweek_OTU),perl=TRUE)
#ignore blank cell
biweek_OTU <- biweek_OTU[, colSums(biweek_OTU) != 0]

OTU.data_biweek.size<-apply(OTU.data_biweek.pivot,2,sum)
OTU.data_biweek.size<-as.data.frame(OTU.data_biweek.size)
OTU.data_biweek.size<-cbind(OTU.data_biweek.size,rownames(OTU.data_biweek.size))
biweek_size<-OTU.data_biweek.size[-1,]
colnames(biweek_size)<-c("Freq","Combination_Number")

iqtree<-read.tree(file="Data/iqtree_treePL_fixed_v4.newick")

iqtree.otu_table<-otu_table(biweek_OTU,taxa_are_rows = TRUE)

iqtree.phyloseq<-phyloseq(iqtree.otu_table,iqtree)

UniFrac.dist<-data.frame(as.matrix(UniFrac(iqtree.phyloseq,FALSE)))

colnames(UniFrac.dist)<-gsub("^X","",colnames(UniFrac.dist),perl=TRUE)

Unifrac_biweek.nmds2 <- metaMDS(UniFrac.dist, dist="jaccard",k=4)

#################
#NMDS for Doug
#################
Unifrac_biweek.NMDS = data.frame(NMDS1 = Unifrac_biweek.nmds2$points[,1], NMDS2 = Unifrac_biweek.nmds2$points[,2])
Unifrac_biweek.NMDS$sample <- rownames(Unifrac_biweek.NMDS)

Unifrac_biweek.NMDS_merge<-dplyr::inner_join(Unifrac_biweek.NMDS,biweek_size,by=c("sample"="Combination_Number"))

Unifrac_biweek.NMDS_merge<-dplyr::inner_join(Unifrac_biweek.NMDS_merge,unique(Biweek_Month_merge[,-1:-4]),by=c("sample"="Combination_Number"))

Unifrac_biweek.NMDS$group <- as.factor(Unifrac_biweek.NMDS_merge$month)

#NMDS subgroup
Unifrac_biweek.NMDS_for_doug_Plot<-Unifrac_biweek.NMDS
Unifrac_biweek.NMDS_for_doug_Plot$month<-gsub("^[0-9]+_","",Unifrac_biweek.NMDS_for_doug_Plot$group,perl=TRUE)
Unifrac_biweek.NMDS_for_doug_Plot$elevation<-gsub("[1-4]2M?$","",gsub("_[0-9]+_[0-9]+$","",gsub("DLJ_","",Unifrac_biweek.NMDS_for_doug_Plot$sample,perl=TRUE),perl=TRUE),perl=TRUE)
Unifrac_biweek.NMDS_for_doug_Plot$Season<-gsub("^[0-9]+$","Summer",gsub("^11$|^12$|^01$|^02$|^03$","Winter",gsub("[0-9]+_","",Unifrac_biweek.NMDS_for_doug_Plot$group,perl=TRUE),perl=TRUE),perl=TRUE)
Unifrac_biweek.NMDS_for_doug_Plot$Location<-gsub("^201[45]$","South",gsub("^2016$","North",gsub("^\\S+[0-9]+_","",gsub("^\\S+17_2015$|^\\S+18_2015$|^\\S+19_2015$|^\\S+20_2015$|^\\S+21_2015$|^\\S+22_2015$|^\\S+23_2015$|^\\S+24_2015$|^\\S+25_2015$|^\\S+26_2015$","North",Unifrac_biweek.NMDS_for_doug_Plot$sample,perl=TRUE),perl=TRUE),perl=TRUE),perl=TRUE)

Unifrac_NMDS <- Unifrac_biweek.NMDS_for_doug_Plot

# definition of season needs to be adjusted
Unifrac_NMDS$Season[Unifrac_NMDS$month %in% c("03","04","05")] <- "Spring"
Unifrac_NMDS$Season[Unifrac_NMDS$month %in% c("06","07","08")] <- "Summer"
Unifrac_NMDS$Season[Unifrac_NMDS$month %in% c("09","10","11")] <- "Autumn"
Unifrac_NMDS$Season[Unifrac_NMDS$month %in% c("12","01","02")] <- "Winter"
Unifrac_NMDS$Location_and_elevation <- paste(Unifrac_NMDS$Location,Unifrac_NMDS$elevation)

####################### Now start plotting

#my.colors <- c("gray","gray","#8a5656","#618dc6")
my.colors <- c("gray","gray","#1d953f","black")

#####################################
# Now let's plot different elevation seperately
#####################################

Unifrac.season.nmds.hul <- gg_ordiplot(Unifrac_biweek.nmds2, groups =Unifrac_NMDS$Season, show.groups = c("Summer","Winter") , hull = TRUE, spiders = TRUE,label = TRUE, ellipse = FALSE, plot = FALSE,conf=0.95,kind = "se")

Unifrac.location.nmds.spider <- gg_ordiplot(Unifrac_biweek.nmds2, groups =Unifrac_NMDS$Location, hull = TRUE, spiders = TRUE,label = TRUE, ellipse = FALSE, plot = FALSE,conf=0.95,kind = "se")

Unifrac_NMDS_1400m <- ggplot(Unifrac_NMDS,aes(x=NMDS1, y=NMDS2))+
  geom_point(data = Unifrac_NMDS[(Unifrac_NMDS$Season=="Summer"|Unifrac_NMDS$Season=="Winter")&Unifrac_NMDS$elevation==1400,],aes(color=Season,shape=Location_and_elevation), size =2)+
  scale_shape_manual(values=c(19,1,15,0))+
  scale_color_manual(values = alpha(my.colors,0.5))+ geom_segment(data=Unifrac.location.nmds.spider$df_spiders[which(rownames(Unifrac.location.nmds.spider$df_spiders) %in% rownames(Unifrac_NMDS[(Unifrac_NMDS$Season=="Summer"|Unifrac_NMDS$Season=="Winter")&Unifrac_NMDS$elevation==1400,])),][which(Unifrac.location.nmds.spider$df_spiders[which(rownames(Unifrac.location.nmds.spider$df_spiders) %in% rownames(Unifrac_NMDS[(Unifrac_NMDS$Season=="Summer"|Unifrac_NMDS$Season=="Winter")&Unifrac_NMDS$elevation==1400,])),]$Group=="North"),], aes(x=cntr.x, xend=x, y=cntr.y, yend=y,color=Group,linetype=1),alpha=0.4,show.legend=FALSE) +geom_segment(data=Unifrac.location.nmds.spider$df_spiders[which(rownames(Unifrac.location.nmds.spider$df_spiders) %in% rownames(Unifrac_NMDS[(Unifrac_NMDS$Season=="Summer"|Unifrac_NMDS$Season=="Winter")&Unifrac_NMDS$elevation==1400,])),][which(Unifrac.location.nmds.spider$df_spiders[which(rownames(Unifrac.location.nmds.spider$df_spiders) %in% rownames(Unifrac_NMDS[(Unifrac_NMDS$Season=="Summer"|Unifrac_NMDS$Season=="Winter")&Unifrac_NMDS$elevation==1400,])),]$Group=="South"),], aes(x=cntr.x, xend=x, y=cntr.y, yend=y,color=Group,linetype=2),alpha=0.4,show.legend=FALSE)+ scale_linetype_identity()+labs(x = paste("NMDS1"), y = paste("NMDS2")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  labs(title = paste("Unifrac on denovo tree (1400m)")) +annotate("text",x=-0.06,y=-0.30,label=paste('Stress =', round(Unifrac_biweek.nmds2$stress, 4)),size=3)+theme(legend.position = "none")
Unifrac_NMDS_1400m

Unifrac_NMDS_1800m <- ggplot(Unifrac_NMDS,aes(x=NMDS1, y=NMDS2))+
  geom_point(data = Unifrac_NMDS[(Unifrac_NMDS$Season=="Summer"|Unifrac_NMDS$Season=="Winter")&Unifrac_NMDS$elevation==1800,],aes(color=Season,shape=Location), size =2)+
  scale_shape_manual(values=c(19,1,15,0))+
  scale_color_manual(values = alpha(my.colors,0.5))+ geom_segment(data=Unifrac.location.nmds.spider$df_spiders[which(rownames(Unifrac.location.nmds.spider$df_spiders) %in% rownames(Unifrac_NMDS[(Unifrac_NMDS$Season=="Summer"|Unifrac_NMDS$Season=="Winter")&Unifrac_NMDS$elevation==1800,])),][which(Unifrac.location.nmds.spider$df_spiders[which(rownames(Unifrac.location.nmds.spider$df_spiders) %in% rownames(Unifrac_NMDS[(Unifrac_NMDS$Season=="Summer"|Unifrac_NMDS$Season=="Winter")&Unifrac_NMDS$elevation==1800,])),]$Group=="North"),], aes(x=cntr.x, xend=x, y=cntr.y, yend=y,color=Group,linetype=1),alpha=0.4,show.legend=FALSE) +geom_segment(data=Unifrac.location.nmds.spider$df_spiders[which(rownames(Unifrac.location.nmds.spider$df_spiders) %in% rownames(Unifrac_NMDS[(Unifrac_NMDS$Season=="Summer"|Unifrac_NMDS$Season=="Winter")&Unifrac_NMDS$elevation==1800,])),][which(Unifrac.location.nmds.spider$df_spiders[which(rownames(Unifrac.location.nmds.spider$df_spiders) %in% rownames(Unifrac_NMDS[(Unifrac_NMDS$Season=="Summer"|Unifrac_NMDS$Season=="Winter")&Unifrac_NMDS$elevation==1800,])),]$Group=="South"),], aes(x=cntr.x, xend=x, y=cntr.y, yend=y,color=Group,linetype=2),alpha=0.4,show.legend=FALSE)+ scale_linetype_identity()+ labs(x = paste("NMDS1"), y = paste("NMDS2")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  labs(title = paste("Unifrac on denovo tree (1800m)")) +annotate("text",x=-0.03,y=-0.19,label=paste('Stress =', round(Unifrac_biweek.nmds2$stress, 4)),size=3)+
  theme(legend.position = "none")

Unifrac_NMDS_1800m

###################
##using Unifrac on mitogenome tree
###################
library(phyloseq) #phyloseq v1.30.0
library(ggplot2) #ggplot2 v3.3.4
library(plyr) #plyr v1.8.6
library(vegan) #vegan v2.5-6
library(anytime) #anytime v0.3.9
library(dplyr) #dplyr v1.0.2
library(ggordiplots) #ggordiplots v0.3.0
library(reshape) #reshape v0.8.8
library(patchwork) #patchwork v1.0.0
library(ggtree) #ggtree v2.0.4
library(stringr) #stringr v1.4.0

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
##create size group by biweek and malaise trap number.
comm <- read.table("Data/GLGS_Community_Data_fixed_three_cases_20201019.txt",header=TRUE,sep = "\t",stringsAsFactors = FALSE)
biweek<-cbind(comm$OTU,comm$Channel_ID,comm$Collection_Date_formatted_1,sub(" ","",paste0(ceiling(lubridate::week(as.Date(anydate(comm$Collection_Date_formatted_1))) / 2),format(as.Date(anydate(comm$Collection_Date_formatted_1)), "_%Y"))))
colnames(biweek)<-c("OTU","Channel_ID","Collection_Date_formatted_1","Biweek_Number")
biweek<-as.data.frame(biweek)
biweek<-cbind(biweek,paste0(biweek$Channel_ID,"_",biweek$Biweek_Number))
colnames(biweek)<-c("OTU","Channel_ID","Collection_Date_formatted_1","Biweek_Number","Combination_Number")
biweek<-as.data.frame(biweek)
biweek$Collection_Date_formatted_1 <- lubridate::dmy(biweek$Collection_Date_formatted_1)
x<-dplyr::arrange(biweek, Collection_Date_formatted_1)

date.min<-data.frame(Biweek_Number="",min_date="")
date.min<-date.min[-1,]
date.list<-c()
for (i in unique(x$Biweek_Number)){ 
  
  date.list<-filter(x,x$Biweek_Number==i)
  date.min_0<-data.frame(Biweek_Number=i,min_date=min(as.Date(date.list$Collection_Date_formatted_1)))
  date.min<-rbind(date.min,date.min_0)
  date.list<-c() 
}
month <- date.min
data_month<-cbind(date.min,month=format(as.Date(anydate(date.min$min_date)), "%Y_%m"))

Biweek_Month_merge<-dplyr::inner_join(biweek,data_month,by=c("Biweek_Number"="Biweek_Number"))

#need to ignore the blank cell for the sample of failure finally.
data_biweek <- filter(Biweek_Month_merge,Biweek_Month_merge$OTU!='')
data_biweek_1<-data_biweek[,-2:-4]
data_biweek_1<-data_biweek_1[,-3:-4]
data_biweek_1 <-as.data.frame(table(data_biweek_1))
OTU.data_biweek.pivot<-cast(data_biweek_1, OTU ~ Combination_Number)
## Using Freq as value column.  Use the value argument to cast to override this choice
biweek_OTU<-OTU.data_biweek.pivot[-1,]
rownames(biweek_OTU)<-biweek_OTU[,1]
biweek_OTU<-biweek_OTU[,-1:-2]
biweek_OTU<-data.frame(biweek_OTU)
colnames(biweek_OTU)<-gsub("^X","",colnames(biweek_OTU),perl=TRUE)
#ignore blank cell
biweek_OTU <- biweek_OTU[, colSums(biweek_OTU) != 0]
#write.csv(biweek_OTU,file="biweek_OTU.csv")
OTU.data_biweek.pivot<-apply(OTU.data_biweek.pivot,2,sum)
OTU.data_biweek.pivot<-as.data.frame(OTU.data_biweek.pivot)
OTU.data_biweek.pivot<-cbind(OTU.data_biweek.pivot,rownames(OTU.data_biweek.pivot))
biweek_size<-OTU.data_biweek.pivot[-1,]
colnames(biweek_size)<-c("Freq","Combination_Number")

placement_tree<-read.tree(file="Data/mitogenome_treePL_fixed_v4.newick")

rownames(biweek_OTU)<-gsub(":","-",rownames(biweek_OTU),perl=TRUE)
placement.otu_table<-otu_table(biweek_OTU,taxa_are_rows = TRUE)

placement.phyloseq<-phyloseq(placement.otu_table,placement_tree)
#I failed to use unifrac on non-fully_resolve mitogenome tree, because many pairwise distance equal to 1 or 0. The ordination results are less pronounceable and don't make sense.
UniFrac.dist<-data.frame(as.matrix(UniFrac(placement.phyloseq,FALSE)))
colnames(UniFrac.dist)<-gsub("^X","",colnames(UniFrac.dist),perl=TRUE)

EPA_no_heur_whole_biweek.nmds2 <- metaMDS(UniFrac.dist,dist="jaccard",k=4)

#################
#NMDS for Doug
#################
Placement_biweek.NMDS = data.frame(NMDS1 = EPA_no_heur_whole_biweek.nmds2$points[,1], NMDS2 = EPA_no_heur_whole_biweek.nmds2$points[,2])
Placement_biweek.NMDS$sample <- rownames(Placement_biweek.NMDS)

Placement_biweek.NMDS_merge<-dplyr::inner_join(Placement_biweek.NMDS,biweek_size,by=c("sample"="Combination_Number"))

Placement_biweek.NMDS_merge<-dplyr::inner_join(Placement_biweek.NMDS_merge,unique(Biweek_Month_merge[,-1:-4]),by=c("sample"="Combination_Number"))

Placement_biweek.NMDS$group <- as.factor(Placement_biweek.NMDS_merge$month)

#NMDS subgroup
Placement_biweek.NMDS_for_doug_Plot<-Placement_biweek.NMDS
Placement_biweek.NMDS_for_doug_Plot$month<-gsub("^[0-9]+_","",Placement_biweek.NMDS_for_doug_Plot$group,perl=TRUE)
Placement_biweek.NMDS_for_doug_Plot$elevation<-gsub("[1-4]2M?$","",gsub("_[0-9]+_[0-9]+$","",gsub("DLJ_","",Placement_biweek.NMDS_for_doug_Plot$sample,perl=TRUE),perl=TRUE),perl=TRUE)
Placement_biweek.NMDS_for_doug_Plot$Season<-gsub("^[0-9]+$","Summer",gsub("^11$|^12$|^01$|^02$|^03$","Winter",gsub("[0-9]+_","",Placement_biweek.NMDS_for_doug_Plot$group,perl=TRUE),perl=TRUE),perl=TRUE)
Placement_biweek.NMDS_for_doug_Plot$Location<-gsub("^201[45]$","South",gsub("^2016$","North",gsub("^\\S+[0-9]+_","",gsub("^\\S+17_2015$|^\\S+18_2015$|^\\S+19_2015$|^\\S+20_2015$|^\\S+21_2015$|^\\S+22_2015$|^\\S+23_2015$|^\\S+24_2015$|^\\S+25_2015$|^\\S+26_2015$","North",Placement_biweek.NMDS_for_doug_Plot$sample,perl=TRUE),perl=TRUE),perl=TRUE),perl=TRUE)

Placement_NMDS <- Placement_biweek.NMDS_for_doug_Plot

# definition of season needs to be adjusted
Placement_NMDS$Season[Placement_NMDS$month %in% c("03","04","05")] <- "Spring"
Placement_NMDS$Season[Placement_NMDS$month %in% c("06","07","08")] <- "Summer"
Placement_NMDS$Season[Placement_NMDS$month %in% c("09","10","11")] <- "Autumn"
Placement_NMDS$Season[Placement_NMDS$month %in% c("12","01","02")] <- "Winter"
Placement_NMDS$Location_and_elevation <- paste(Placement_NMDS$Location,Placement_NMDS$elevation)


####################### Now start plotting

#my.colors <- c("gray","gray","#8a5656","#618dc6")
my.colors <- c("gray","gray","#1d953f","black")

#####################################
# Now let's plot different elevation seperately
#####################################
Placement.season.nmds.hul <- gg_ordiplot(EPA_no_heur_whole_biweek.nmds2, groups =Placement_NMDS$Season, show.groups = c("Summer","Winter") , hull = TRUE, spiders = TRUE,label = TRUE, ellipse = FALSE, plot = FALSE,conf=0.95,kind = "se")

Placement.location.nmds.spider <- gg_ordiplot(EPA_no_heur_whole_biweek.nmds2, groups =Placement_NMDS$Location, hull = TRUE, spiders = TRUE,label = TRUE, ellipse = FALSE, plot = FALSE,conf=0.95,kind = "se")

Placement_NMDS_1400m <- ggplot(Placement_NMDS,aes(x=NMDS1, y=NMDS2))+
  geom_point(data = Placement_NMDS[(Placement_NMDS$Season=="Summer"|Placement_NMDS$Season=="Winter")&Placement_NMDS$elevation==1400,],aes(color=Season,shape=Location_and_elevation), size =2)+
  scale_shape_manual(values=c(19,1,15,0))+
  scale_color_manual(values = alpha(my.colors,0.5))+ geom_segment(data=Placement.location.nmds.spider$df_spiders[which(rownames(Placement.location.nmds.spider$df_spiders) %in% rownames(Placement_NMDS[(Placement_NMDS$Season=="Summer"|Placement_NMDS$Season=="Winter")&Placement_NMDS$elevation==1400,])),][which(Placement.location.nmds.spider$df_spiders[which(rownames(Placement.location.nmds.spider$df_spiders) %in% rownames(Placement_NMDS[(Placement_NMDS$Season=="Summer"|Placement_NMDS$Season=="Winter")&Placement_NMDS$elevation==1400,])),]$Group=="North"),], aes(x=cntr.x, xend=x, y=cntr.y, yend=y,color=Group,linetype=1),alpha=0.4,show.legend=FALSE) +geom_segment(data=Placement.location.nmds.spider$df_spiders[which(rownames(Placement.location.nmds.spider$df_spiders) %in% rownames(Placement_NMDS[(Placement_NMDS$Season=="Summer"|Placement_NMDS$Season=="Winter")&Placement_NMDS$elevation==1400,])),][which(Placement.location.nmds.spider$df_spiders[which(rownames(Placement.location.nmds.spider$df_spiders) %in% rownames(Placement_NMDS[(Placement_NMDS$Season=="Summer"|Placement_NMDS$Season=="Winter")&Placement_NMDS$elevation==1400,])),]$Group=="South"),], aes(x=cntr.x, xend=x, y=cntr.y, yend=y,color=Group,linetype=2),alpha=0.4,show.legend=FALSE)+ scale_linetype_identity()+labs(x = paste("NMDS1"), y = paste("NMDS2")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  labs(title = paste("Unifrac on mitogenome tree (1400m)"))+annotate("text",x=-0.226,y=-0.45,label=paste('Stress =', round(EPA_no_heur_whole_biweek.nmds2$stress, 4)),size=3)+theme(legend.position = "none")
Placement_NMDS_1400m

Placement_NMDS_1800m <- ggplot(Placement_NMDS,aes(x=NMDS1, y=NMDS2))+
  geom_point(data = Placement_NMDS[(Placement_NMDS$Season=="Summer"|Placement_NMDS$Season=="Winter")&Placement_NMDS$elevation==1800,],aes(color=Season,shape=Location), size =2)+
  scale_shape_manual(values=c(19,1,15,0))+
  scale_color_manual(values = alpha(my.colors,0.5))+ geom_segment(data=Placement.location.nmds.spider$df_spiders[which(rownames(Placement.location.nmds.spider$df_spiders) %in% rownames(Placement_NMDS[(Placement_NMDS$Season=="Summer"|Placement_NMDS$Season=="Winter")&Placement_NMDS$elevation==1800,])),][which(Placement.location.nmds.spider$df_spiders[which(rownames(Placement.location.nmds.spider$df_spiders) %in% rownames(Placement_NMDS[(Placement_NMDS$Season=="Summer"|Placement_NMDS$Season=="Winter")&Placement_NMDS$elevation==1800,])),]$Group=="North"),], aes(x=cntr.x, xend=x, y=cntr.y, yend=y,color=Group,linetype=1),alpha=0.4,show.legend=FALSE) +geom_segment(data=Placement.location.nmds.spider$df_spiders[which(rownames(Placement.location.nmds.spider$df_spiders) %in% rownames(Placement_NMDS[(Placement_NMDS$Season=="Summer"|Placement_NMDS$Season=="Winter")&Placement_NMDS$elevation==1800,])),][which(Placement.location.nmds.spider$df_spiders[which(rownames(Placement.location.nmds.spider$df_spiders) %in% rownames(Placement_NMDS[(Placement_NMDS$Season=="Summer"|Placement_NMDS$Season=="Winter")&Placement_NMDS$elevation==1800,])),]$Group=="South"),], aes(x=cntr.x, xend=x, y=cntr.y, yend=y,color=Group,linetype=2),alpha=0.4,show.legend=FALSE)+ scale_linetype_identity()+ labs(x = paste("NMDS1"), y = paste("NMDS2")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  labs(title = paste("Unifrac on mitogenome tree (1800m)"))+annotate("text",x=-0.27,y=-0.33,label=paste('Stress =', round(EPA_no_heur_whole_biweek.nmds2$stress, 4)),size=3)+
  theme(legend.position = "none")
Placement_NMDS_1800m 

dev.new()
ggdraw() + draw_plot(OTU_NMDS_1400m, x = 0, y = .66, width = 0.5, height = .33) +
  draw_plot(Unifrac_NMDS_1400m, x = 0, y = .33, width = 0.5, height = .33) +
  draw_plot(Placement_NMDS_1400m, x = 0, y = .0, width = 0.5, height = .33)+
  draw_plot(OTU_NMDS_1800m, x = 0.5, y = .66, width = 0.5, height = .33) +
  draw_plot(Unifrac_NMDS_1800m, x = 0.5, y = .33, width = 0.5, height = .33) +
  draw_plot(Placement_NMDS_1800m, x = 0.5, y = .0, width = 0.5, height = .33)
ggsave("figS2_NMDS_non_ultrametric_20220124.pdf", width = 22, height = 20, units = "cm")
dev.off()