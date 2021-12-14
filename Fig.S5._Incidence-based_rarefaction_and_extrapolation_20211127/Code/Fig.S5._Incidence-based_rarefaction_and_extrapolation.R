###############################################################
###############################alpha diversity###################
###############################################################
#This is taking soooo long for case B and C and need to run on hpc.
#######################
#Case A OTU Table
#######################
library(ggplot2) #ggplot2 v3.3.4
library(plyr) #plyr v1.8.6
library(vegan) #vegan v2.5-6
library(ape) #ape v5.4
library(anytime) #anytime v0.3.9
library(dplyr) #dplyr v1.0.2
library(ggtree) #ggtree v2.0.4
library(reshape) #reshape v0.8.8
library(patchwork) #patchwork v1.0.0
library(iNEXT3D) #iNEXT3D v0.0.1
library(phytools) #phytools v0.7-80
library(phangorn) #phangorn v2.5.5
library(rlist) #rlist v0.4.6.1
library(phyclust) #phyclust v0.1-30

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
## Using Freq as value column.  Use the value argument to cast to override this choice
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

OTU.data_biweek.pivot<-cast(as.data.frame(OTU.data_biweek), Combination_Number ~ Sumaclust_97_OTU)
## Using Freq as value column.  Use the value argument to cast to override this choice

#deleted the first column including biweek numbers.
row.names(OTU.data_biweek.pivot)<-OTU.data_biweek.pivot[,1]
OTU.data_biweek.pivot_fixed<-as.data.frame(OTU.data_biweek.pivot[,-1:-2])#fixed to remove the blank column(e.g.,DATA[,-2]).

#create broad group 
OTU_biweek.PCA_merge.broad<-dplyr::inner_join(unique(Biweek_Month_merge[,-1:-4]),biweek_size,by=c("Combination_Number"="Combination_Number"))
OTU_biweek.PCA_merge.broad$month<-gsub("^[0-9]+_","",OTU_biweek.PCA_merge.broad$month,perl=TRUE)
OTU_biweek.PCA_merge.broad$season<-gsub("^09$|^10$|^11$","Autumn",gsub("^03$|^04$|^05$","Spring",gsub("^06$|^07$|^08$","Summer",gsub("^12$|^01$|^02$","Winter",gsub("[0-9]+_","",OTU_biweek.PCA_merge.broad$month,perl=TRUE),perl=TRUE),perl=TRUE),perl=TRUE),perl=TRUE)
OTU_biweek.PCA_merge.broad$location<-gsub("^201[45]$","South",gsub("^2016$","North",gsub("^\\S+[0-9]+_","",gsub("^\\S+17_2015$|^\\S+18_2015$|^\\S+19_2015$|^\\S+20_2015$|^\\S+21_2015$|^\\S+22_2015$|^\\S+23_2015$|^\\S+24_2015$|^\\S+25_2015$|^\\S+26_2015$","North",OTU_biweek.PCA_merge.broad$Combination_Number,perl=TRUE),perl=TRUE),perl=TRUE),perl=TRUE)
OTU_biweek.PCA_merge.broad$elevation<-gsub("^DLJ_","",gsub("[1-4]2M?_\\S+_\\S+$","",OTU_biweek.PCA_merge.broad$Combination_Number,perl=TRUE),perl=TRUE)
row.names(OTU_biweek.PCA_merge.broad)<-OTU_biweek.PCA_merge.broad[,1]
OTU_biweek.PCA_merge.broad <- OTU_biweek.PCA_merge.broad[which(rowSums(OTU_biweek.PCA_merge.broad==0)==0),]
#for Season&location
group = data.frame(sample = rownames(OTU_biweek.PCA_merge.broad),
                   group = paste0(OTU_biweek.PCA_merge.broad$location,"_",OTU_biweek.PCA_merge.broad$season))

group_name = unique(group$group)

loss.abun=list()
loss.inci=list()
list.loss.inci_V2=list()

for (i in 1:length(group_name)){
  group_i = subset(group, group %in% group_name[i])
  otu_i = t(OTU.data_biweek.pivot_fixed[group_i$sample, ])
  temp_i<-t(apply(otu_i,1,sum))
  loss.abun[i]<- unlist(list(list(temp_i),loss.abun), recursive=FALSE)
  names(loss.abun)[i] <- paste0(group_name[i])
  otu_i[otu_i>1] <- 1
  loss.inci[i]<- unlist(list(list(otu_i),loss.inci), recursive=FALSE)
  names(loss.inci)[i] <- paste0(group_name[i])
  
}

#for whole dataset
loss.abun["South_Spring"]<-NULL
loss.abun["North_Spring"]<-NULL
loss.abun["North_Autumn"]<-NULL
loss.abun["South_Autumn"]<-NULL

loss.inci["South_Spring"]<-NULL
loss.inci["North_Spring"]<-NULL
loss.inci["North_Autumn"]<-NULL
loss.inci["South_Autumn"]<-NULL

#0 set:whole dataset
entire.loss.inci<-list.cbind(loss.inci)
list.loss.inci_V2[1]<-unlist(list(list(entire.loss.inci),list.loss.inci_V2), recursive=FALSE)
names(list.loss.inci_V2)[1] <- paste0("Whole Dataset")

#1 loss set:South_Summer
minus_ss.loss.inci<-loss.inci
minus_ss.loss.inci["South_Summer"]<-NULL
minus_ss.loss.inci<-list.cbind(minus_ss.loss.inci)
list.loss.inci_V2[2]<-unlist(list(list(minus_ss.loss.inci),list.loss.inci_V2), recursive=FALSE)
names(list.loss.inci_V2)[2] <- paste0("South Summer Loss")

#2 loss set:North_Winter
minus_nw.loss.inci<-loss.inci
minus_nw.loss.inci["North_Winter"]<-NULL
minus_nw.loss.inci<-list.cbind(minus_nw.loss.inci)
list.loss.inci_V2[3]<-unlist(list(list(minus_nw.loss.inci),list.loss.inci_V2), recursive=FALSE)
names(list.loss.inci_V2)[3] <- paste0("North Winter Loss")

#3 loss set:North_Summer
minus_ns.loss.inci<-loss.inci
minus_ns.loss.inci["North_Summer"]<-NULL
minus_ns.loss.inci<-list.cbind(minus_ns.loss.inci)
list.loss.inci_V2[4]<-unlist(list(list(minus_ns.loss.inci),list.loss.inci_V2), recursive=FALSE)
names(list.loss.inci_V2)[4] <- paste0("North Summer Loss")

#4 loss set:South_Winter
minus_sw.loss.inci<-loss.inci
minus_sw.loss.inci["South_Winter"]<-NULL
minus_sw.loss.inci<-list.cbind(minus_sw.loss.inci)
list.loss.inci_V2[5]<-unlist(list(list(minus_sw.loss.inci),list.loss.inci_V2), recursive=FALSE)
names(list.loss.inci_V2)[5] <- paste0("South Winter Loss")

#5 loss set:Summer
minus_s.loss.inci<-loss.inci
minus_s.loss.inci["South_Summer"]<-NULL
minus_s.loss.inci["North_Summer"]<-NULL
minus_s.loss.inci<-list.cbind(minus_s.loss.inci)
list.loss.inci_V2[6]<-unlist(list(list(minus_s.loss.inci),list.loss.inci_V2), recursive=FALSE)
names(list.loss.inci_V2)[6] <- paste0("Summer Loss")

#6 loss set:Winter
minus_w.loss.inci<-loss.inci
minus_w.loss.inci["South_Winter"]<-NULL
minus_w.loss.inci["North_Winter"]<-NULL
minus_w.loss.inci<-list.cbind(minus_w.loss.inci)
list.loss.inci_V2[7]<-unlist(list(list(minus_w.loss.inci),list.loss.inci_V2), recursive=FALSE)
names(list.loss.inci_V2)[7] <- paste0("Winter Loss")

#7 loss set:1400
otu_loss_1400 = t(OTU.data_biweek.pivot_fixed[group[which(grepl("^1800|^DLJ_1800",group$sample)),][which(grepl("_Winter$|Summer$",group[which(grepl("^1800|^DLJ_1800",group$sample)),]$group)),]$sample, ])
otu_loss_1400[otu_loss_1400>1] <- 1
list.loss.inci_V2[8]<-unlist(list(list(otu_loss_1400),list.loss.inci_V2), recursive=FALSE)
names(list.loss.inci_V2)[8] <- paste0("1400 Loss")

#8 loss set:1800
otu_loss_1800 = t(OTU.data_biweek.pivot_fixed[group[which(grepl("^1400|^DLJ_1400",group$sample)),][which(grepl("_Winter$|Summer$",group[which(grepl("^1400|^DLJ_1400",group$sample)),]$group)),]$sample, ])
otu_loss_1800[otu_loss_1800>1] <- 1
list.loss.inci_V2[9]<-unlist(list(list(otu_loss_1800),list.loss.inci_V2), recursive=FALSE)
names(list.loss.inci_V2)[9] <- paste0("1800 Loss")

#9 loss set:south
otu_loss_south = t(OTU.data_biweek.pivot_fixed[group[which(grepl("^North_",group$group)),][which(grepl("_Winter$|Summer$",group[which(grepl("^North_",group$group)),]$group)),]$sample, ])
otu_loss_south[otu_loss_south>1] <- 1
list.loss.inci_V2[10]<-unlist(list(list(otu_loss_south),list.loss.inci_V2), recursive=FALSE)
names(list.loss.inci_V2)[10] <- paste0("South Loss")

#10 loss set:north
otu_loss_north = t(OTU.data_biweek.pivot_fixed[group[which(grepl("^South_",group$group)),][which(grepl("_Winter$|Summer$",group[which(grepl("^South_",group$group)),]$group)),]$sample, ])
otu_loss_north[otu_loss_north>1] <- 1
list.loss.inci_V2[11]<-unlist(list(list(otu_loss_north),list.loss.inci_V2), recursive=FALSE)
names(list.loss.inci_V2)[11] <- paste0("North Loss")

memory.limit(100000)
# Hill numbers (q):  0 = sp richness, 1 = Shannon, 2 = inverse Simpson
loss.out.inci <- iNEXT3D(list.loss.inci_V2, class = 'TD', q = c(0,1,2), datatype = "incidence_raw")
save(loss.out.inci,file="p.loss.inci_OTU_20211214")

###################
#CASE B for denovo tree
###################
library(ggplot2) #ggplot2 v3.3.4
library(plyr) #plyr v1.8.6
library(vegan) #vegan v2.5-6
library(ape) #ape v5.4
library(anytime) #anytime v0.3.9
library(dplyr) #dplyr v1.0.2
library(ggtree) #ggtree v2.0.4
library(reshape) #reshape v0.8.8
library(patchwork) #patchwork v1.0.0
library(iNEXT3D) #iNEXT3D v0.0.1
library(phytools) #phytools v0.7-80
library(phangorn) #phangorn v2.5.5
library(rlist) #rlist v0.4.6.1
library(phyclust) #phyclust v0.1-30

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
## Using Freq as value column.  Use the value argument to cast to override this choice
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

iqtree.otu_table<-biweek_OTU
#create subgroupings 
unifrac_biweek.PCA_merge<-dplyr::inner_join(unique(Biweek_Month_merge[,-1:-4]),biweek_size,by=c("Combination_Number"="Combination_Number"))
unifrac_biweek.PCA_merge$month<-gsub("^[0-9]+_","",unifrac_biweek.PCA_merge$month,perl=TRUE)
unifrac_biweek.PCA_merge$season<-gsub("^09$|^10$|^11$","Autumn",gsub("^03$|^04$|^05$","Spring",gsub("^06$|^07$|^08$","Summer",gsub("^12$|^01$|^02$","Winter",gsub("[0-9]+_","",unifrac_biweek.PCA_merge$month,perl=TRUE),perl=TRUE),perl=TRUE),perl=TRUE),perl=TRUE)
unifrac_biweek.PCA_merge$location<-gsub("^201[45]$","South",gsub("^2016$","North",gsub("^\\S+[0-9]+_","",gsub("^\\S+17_2015$|^\\S+18_2015$|^\\S+19_2015$|^\\S+20_2015$|^\\S+21_2015$|^\\S+22_2015$|^\\S+23_2015$|^\\S+24_2015$|^\\S+25_2015$|^\\S+26_2015$","North",unifrac_biweek.PCA_merge$Combination_Number,perl=TRUE),perl=TRUE),perl=TRUE),perl=TRUE)
unifrac_biweek.PCA_merge$elevation<-gsub("^DLJ_","",gsub("[1-4]2M?_\\S+_\\S+$","",unifrac_biweek.PCA_merge$Combination_Number,perl=TRUE),perl=TRUE)

row.names(unifrac_biweek.PCA_merge)<-unifrac_biweek.PCA_merge[,1]
unifrac_biweek.PCA_merge$group<-paste0(unifrac_biweek.PCA_merge$location,"_", unifrac_biweek.PCA_merge$season,"_",unifrac_biweek.PCA_merge$elevation)
unifrac_biweek.PCA_merge <- unifrac_biweek.PCA_merge[which(rowSums(unifrac_biweek.PCA_merge==0)==0),]
group = data.frame(sample = rownames(unifrac_biweek.PCA_merge),
                   group = paste0(unifrac_biweek.PCA_merge$location,"_",unifrac_biweek.PCA_merge$season))

group_name = unique(group$group)
#test for loss set
loss.unifrac.inci = list()
loss.unifrac.abun = list()
for (i in 1:length(group_name)){
  group_i = subset(group, group %in% group_name[i])
  otu_i = iqtree.otu_table[,group_i$sample]
  temp_i<-list(apply(otu_i,1,sum))
  loss.unifrac.abun[i]<- unlist(list(temp_i,loss.unifrac.abun), recursive=FALSE)
  names(loss.unifrac.abun)[i] <- paste0(group_name[i])
  otu_i[otu_i>1] <- 1
  loss.unifrac.inci[i]<- unlist(list(list(otu_i),loss.unifrac.inci), recursive=FALSE)
  names(loss.unifrac.inci)[i] <- paste0(group_name[i])
}

loss.unifrac.abun["South_Spring"]<-NULL
loss.unifrac.abun["North_Spring"]<-NULL
loss.unifrac.abun["North_Autumn"]<-NULL
loss.unifrac.abun["South_Autumn"]<-NULL

loss.unifrac.inci["South_Spring"]<-NULL
loss.unifrac.inci["North_Spring"]<-NULL
loss.unifrac.inci["North_Autumn"]<-NULL
loss.unifrac.inci["South_Autumn"]<-NULL

iqtree.treepl<-read.tree(file="Data/denovo.treePL")
sum(iqtree.treepl$edge.length == 0)
#[1] 0
is.ultrametric(iqtree.treepl)
#iqtree.treepl<-force.ultrametric(iqtree.treepl)#default method: "nnls"
iqtree.treepl<-force.ultrametric(iqtree.treepl,method=c("extend"))
is.ultrametric(iqtree.treepl)

##Assemble process ultrametric tree
iqtree.treepl<-drop.tip(iqtree.treepl,c("XXXXXXXXXXXXXXX----_NNy4Ny7myr_Nymphes_myrmeleonoides_NC_024825", "XXXXXXXXXXXXXXX----_NCh4Chr7nip_Chrysoperla_nipponensis_NC_015093", "XXXXXXXXXXXXXXX----_NR4Ra7xiz_Rapisma_xizangense_NC_023362", "XXXXXXXXXXXXXXX----_NO4Th7lan_Thyridosmylus_langii_NC_021415", "XXXXXXXXXXXXXXX----_NMy4My7imm_Myrmeleon_immanis_NC_024826", "XXXXXXXXXXXXXXX----_NA4As7app_Ascaloptynx_appendiculatus_NC_011277", "XXXXXXXXXXXXXXX----_NPo4Po7pun_Polystoechotes_punctatus_NC_011278", "XXXXXXXXXXXXXXX----_NMa4Di7bis_Ditaxis_biseriata_NC_013257"))

#non-ultrametric tree into a ultrametric tree
##Caution: MPL return negative branch lengths, meaning that it should not be used.
##unifrca.Ult.MPL <- chronoMPL(iqtree)
##computationally expensive
##unifrca.Ult.iqtree <- chronos(iqtree)
#phylo tree contains duplicated node labels, remove them
iqtree.treepl$node.label<-NULL

list.loss.unifrac.inci=list()
#0 set:whole dataset
entire.loss.unifrac.inci<-list.cbind(loss.unifrac.inci)
list.loss.unifrac.inci[1]<-unlist(list(list(entire.loss.unifrac.inci),list.loss.unifrac.inci), recursive=FALSE)
names(list.loss.unifrac.inci)[1] <- paste0("Whole Dataset")

#1 loss set:South_Summer
minus_ss.loss.unifrac.inci<-loss.unifrac.inci
minus_ss.loss.unifrac.inci["South_Summer"]<-NULL
minus_ss.loss.unifrac.inci<-list.cbind(minus_ss.loss.unifrac.inci)
list.loss.unifrac.inci[2]<-unlist(list(list(minus_ss.loss.unifrac.inci),list.loss.unifrac.inci), recursive=FALSE)
names(list.loss.unifrac.inci)[2] <- paste0("South Summer Loss")

#2 loss set:North_Winter
minus_nw.loss.unifrac.inci<-loss.unifrac.inci
minus_nw.loss.unifrac.inci["North_Winter"]<-NULL
minus_nw.loss.unifrac.inci<-list.cbind(minus_nw.loss.unifrac.inci)
list.loss.unifrac.inci[3]<-unlist(list(list(minus_nw.loss.unifrac.inci),list.loss.unifrac.inci), recursive=FALSE)
names(list.loss.unifrac.inci)[3] <- paste0("North Winter Loss")

#3 loss set:North_Summer
minus_ns.loss.unifrac.inci<-loss.unifrac.inci
minus_ns.loss.unifrac.inci["North_Summer"]<-NULL
minus_ns.loss.unifrac.inci<-list.cbind(minus_ns.loss.unifrac.inci)
list.loss.unifrac.inci[4]<-unlist(list(list(minus_ns.loss.unifrac.inci),list.loss.unifrac.inci), recursive=FALSE)
names(list.loss.unifrac.inci)[4] <- paste0("North Summer Loss")

#4 loss set:South_Winter
minus_sw.loss.unifrac.inci<-loss.unifrac.inci
minus_sw.loss.unifrac.inci["South_Winter"]<-NULL
minus_sw.loss.unifrac.inci<-list.cbind(minus_sw.loss.unifrac.inci)
list.loss.unifrac.inci[5]<-unlist(list(list(minus_sw.loss.unifrac.inci),list.loss.unifrac.inci), recursive=FALSE)
names(list.loss.unifrac.inci)[5] <- paste0("South Winter Loss")

#5 loss set:Summer
minus_s.loss.unifrac.inci<-loss.unifrac.inci
minus_s.loss.unifrac.inci["South_Summer"]<-NULL
minus_s.loss.unifrac.inci["North_Summer"]<-NULL
minus_s.loss.unifrac.inci<-list.cbind(minus_s.loss.unifrac.inci)
list.loss.unifrac.inci[6]<-unlist(list(list(minus_s.loss.unifrac.inci),list.loss.unifrac.inci), recursive=FALSE)
names(list.loss.unifrac.inci)[6] <- paste0("Summer Loss")

#6 loss set:Winter
minus_w.loss.unifrac.inci<-loss.unifrac.inci
minus_w.loss.unifrac.inci["South_Winter"]<-NULL
minus_w.loss.unifrac.inci["North_Winter"]<-NULL
minus_w.loss.unifrac.inci<-list.cbind(minus_w.loss.unifrac.inci)
list.loss.unifrac.inci[7]<-unlist(list(list(minus_w.loss.unifrac.inci),list.loss.unifrac.inci), recursive=FALSE)
names(list.loss.unifrac.inci)[7] <- paste0("Winter Loss")

#7 loss set:1400
iqtree_loss_1400 = iqtree.otu_table[,group[which(grepl("^1800|^DLJ_1800",group$sample)),][which(grepl("_Winter$|Summer$",group[which(grepl("^1800|^DLJ_1800",group$sample)),]$group)),]$sample]
iqtree_loss_1400[iqtree_loss_1400>1] <- 1
list.loss.unifrac.inci[8]<-unlist(list(list(iqtree_loss_1400),list.loss.unifrac.inci), recursive=FALSE)
names(list.loss.unifrac.inci)[8] <- paste0("1400 Loss")

#8 loss set:1800
iqtree_loss_1800 = iqtree.otu_table[,group[which(grepl("^1400|^DLJ_1400",group$sample)),][which(grepl("_Winter$|Summer$",group[which(grepl("^1400|^DLJ_1400",group$sample)),]$group)),]$sample]
iqtree_loss_1800[iqtree_loss_1800>1] <- 1
list.loss.unifrac.inci[9]<-unlist(list(list(iqtree_loss_1800),list.loss.unifrac.inci), recursive=FALSE)
names(list.loss.unifrac.inci)[9] <- paste0("1800 Loss")

#9 loss set:South
iqtree_loss_south = iqtree.otu_table[,group[which(grepl("^North_",group$group)),][which(grepl("_Winter$|Summer$",group[which(grepl("^North_",group$group)),]$group)),]$sample]
iqtree_loss_south[iqtree_loss_south>1] <- 1
list.loss.unifrac.inci[10]<-unlist(list(list(iqtree_loss_south),list.loss.unifrac.inci), recursive=FALSE)
names(list.loss.unifrac.inci)[10] <- paste0("South Loss")

#10 loss set:north
iqtree_loss_north = iqtree.otu_table[,group[which(grepl("^South_",group$group)),][which(grepl("_Winter$|Summer$",group[which(grepl("^South_",group$group)),]$group)),]$sample]
iqtree_loss_north[iqtree_loss_north>1] <- 1
list.loss.unifrac.inci[11]<-unlist(list(list(iqtree_loss_north),list.loss.unifrac.inci), recursive=FALSE)
names(list.loss.unifrac.inci)[11] <- paste0("North Loss")

unifrac.inci.loss.set<-data.frame(cbind(list.loss.unifrac.inci$`Whole Dataset`,list.loss.unifrac.inci$`South Summer Loss`,list.loss.unifrac.inci$`North Winter Loss`,list.loss.unifrac.inci$`North Summer Loss`,list.loss.unifrac.inci$`South Winter Loss`,list.loss.unifrac.inci$`Summer Loss`,list.loss.unifrac.inci$`Winter Loss`,list.loss.unifrac.inci$`1400 Loss`,list.loss.unifrac.inci$`1800 Loss`,list.loss.unifrac.inci$`South Loss`,list.loss.unifrac.inci$`North Loss`))
nt<-c("Whole Dataset"=ncol(list.loss.unifrac.inci$`Whole Dataset`),"South Summer Loss"=ncol(list.loss.unifrac.inci$`South Summer Loss`),"North Winter Loss"=ncol(list.loss.unifrac.inci$`North Winter Loss`),"North Summer Loss"=ncol(list.loss.unifrac.inci$`North Summer Loss`),"South Winter Loss"=ncol(list.loss.unifrac.inci$`South Winter Loss`),"Summer Loss"=ncol(list.loss.unifrac.inci$`Summer Loss`),"Winter Loss"=ncol(list.loss.unifrac.inci$`Winter Loss`),"1400 Loss"=ncol(list.loss.unifrac.inci$`1400 Loss`),"1800 Loss"=ncol(list.loss.unifrac.inci$`1800 Loss`),"South Loss"=ncol(list.loss.unifrac.inci$`South Loss`),"North Loss"=ncol(list.loss.unifrac.inci$`North Loss`))

# Hill numbers (q):  0 = sp richness, 1 = Shannon, 2 = inverse Simpson
pd.unifrac.inci <- iNEXT3D(unifrac.inci.loss.set, nT=nt, class = 'PD', tree = iqtree.treepl, q=c(0,1,2), datatype="incidence_raw")
save(pd.unifrac.inci,file="pd.inci_denovo_treepl_20211125") 

###################
#CASE C for placement tree
###################
library(ggplot2) #ggplot2 v3.3.4
library(plyr) #plyr v1.8.6
library(vegan) #vegan v2.5-6
library(ape) #ape v5.4
library(anytime) #anytime v0.3.9
library(dplyr) #dplyr v1.0.2
library(ggtree) #ggtree v2.0.4
library(reshape) #reshape v0.8.8
library(patchwork) #patchwork v1.0.0
library(iNEXT3D) #iNEXT3D v0.0.1
library(phytools) #phytools v0.7-80
library(phangorn) #phangorn v2.5.5
library(rlist) #rlist v0.4.6.1
library(phyclust) #phyclust v0.1-30

##create size group by biweek and malaise trap number.
comm <- read.table("Data//GLGS_Community_Data_fixed_three_cases_20201019.txt",header=TRUE,sep = "\t",stringsAsFactors = FALSE)
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
## Using Freq as value column.  Use the value argument to cast to override this choice
biweek_OTU<-OTU.data_biweek.pivot[-1,]
rownames(biweek_OTU)<-biweek_OTU[,1]
biweek_OTU<-biweek_OTU[,-1:-2]
biweek_OTU<-data.frame(biweek_OTU)
colnames(biweek_OTU)<-gsub("^X","",colnames(biweek_OTU),perl=TRUE)
#ignore blank cell
biweek_OTU <- biweek_OTU[, colSums(biweek_OTU) != 0]

OTU.data_biweek.size<-apply(OTU.data_biweek.pivot,2,sum)
OTU.data_biweek.size<-as.data.frame(OTU.data_biweek.size)
OTU.data_biweek.size<-cbind(OTU.data_biweek.size,rownames(OTU.data_biweek.size))
biweek_size<-OTU.data_biweek.size[-1,]
colnames(biweek_size)<-c("Freq","Combination_Number")

mito.tree.otu_table<-biweek_OTU
#create subgroupings 
unifrac_biweek.PCA_merge<-dplyr::inner_join(unique(Biweek_Month_merge[,-1:-4]),biweek_size,by=c("Combination_Number"="Combination_Number"))
unifrac_biweek.PCA_merge$month<-gsub("^[0-9]+_","",unifrac_biweek.PCA_merge$month,perl=TRUE)
unifrac_biweek.PCA_merge$season<-gsub("^09$|^10$|^11$","Autumn",gsub("^03$|^04$|^05$","Spring",gsub("^06$|^07$|^08$","Summer",gsub("^12$|^01$|^02$","Winter",gsub("[0-9]+_","",unifrac_biweek.PCA_merge$month,perl=TRUE),perl=TRUE),perl=TRUE),perl=TRUE),perl=TRUE)
unifrac_biweek.PCA_merge$location<-gsub("^201[45]$","South",gsub("^2016$","North",gsub("^\\S+[0-9]+_","",gsub("^\\S+17_2015$|^\\S+18_2015$|^\\S+19_2015$|^\\S+20_2015$|^\\S+21_2015$|^\\S+22_2015$|^\\S+23_2015$|^\\S+24_2015$|^\\S+25_2015$|^\\S+26_2015$","North",unifrac_biweek.PCA_merge$Combination_Number,perl=TRUE),perl=TRUE),perl=TRUE),perl=TRUE)
unifrac_biweek.PCA_merge$elevation<-gsub("^DLJ_","",gsub("[1-4]2M?_\\S+_\\S+$","",unifrac_biweek.PCA_merge$Combination_Number,perl=TRUE),perl=TRUE)

row.names(unifrac_biweek.PCA_merge)<-unifrac_biweek.PCA_merge[,1]
unifrac_biweek.PCA_merge$group<-paste0(unifrac_biweek.PCA_merge$location,"_", unifrac_biweek.PCA_merge$season,"_",unifrac_biweek.PCA_merge$elevation)
unifrac_biweek.PCA_merge <- unifrac_biweek.PCA_merge[which(rowSums(unifrac_biweek.PCA_merge==0)==0),]
group = data.frame(sample = rownames(unifrac_biweek.PCA_merge),
                   group = paste0(unifrac_biweek.PCA_merge$location,"_",unifrac_biweek.PCA_merge$season))

group_name = unique(group$group)
#test for loss set
loss.mito.inci = list()
loss.mito.abun = list()
for (i in 1:length(group_name)){
  group_i = subset(group, group %in% group_name[i])
  otu_i = mito.tree.otu_table[,group_i$sample]
  temp_i<-list(apply(otu_i,1,sum))
  loss.mito.abun[i]<- unlist(list(temp_i,loss.mito.abun), recursive=FALSE)
  names(loss.mito.abun)[i] <- paste0(group_name[i])
  otu_i[otu_i>1] <- 1
  loss.mito.inci[i]<- unlist(list(list(otu_i),loss.mito.inci), recursive=FALSE)
  names(loss.mito.inci)[i] <- paste0(group_name[i])
}


loss.mito.abun["South_Spring"]<-NULL
loss.mito.abun["North_Spring"]<-NULL
loss.mito.abun["North_Autumn"]<-NULL
loss.mito.abun["South_Autumn"]<-NULL

loss.mito.inci["South_Spring"]<-NULL
loss.mito.inci["North_Spring"]<-NULL
loss.mito.inci["North_Autumn"]<-NULL
loss.mito.inci["South_Autumn"]<-NULL

#read tree
mito.tree<-read.tree("Data/fully_resolved_queries_10524_whole_13995_mitoalignment_V3_RAPPAS_Reduced_no_heur.newick")
mito.treepl<-read.tree("Data/fully_resolved_queries_10524.treePL")
#fixed the treepl tree tip labels
mito.treepl$tip.label<-mito.tree$tip.label
# count zero length branches
sum(mito.treepl$edge.length == 0)
#[1] 2988
#Then set zero-length branches to be 1/1000000 total tree length
mito.treepl$edge.length[mito.treepl$edge.length==0]<-max(nodeHeights(mito.treepl))*1e-6
is.ultrametric(mito.treepl)
#nnls method have a negative branch, thus do not use this.
#mito.treepl<-force.ultrametric(mito.treepl)#default method: "nnls"
#However, extend method may convert these negative branches to zero. 
mito.treepl<-force.ultrametric(mito.treepl,method=c("extend"))
# count negative branches
sum(mito.treepl$edge.length < 0)
# count zero length branches
sum(mito.treepl$edge.length == 0)
is.ultrametric(mito.treepl)
mito.treepl<-drop.tip(mito.treepl,c("XXXXXXXXXXXXXXX----_NNy4Ny7myr_Nymphes_myrmeleonoides_NC_024825", "XXXXXXXXXXXXXXX----_NCh4Chr7nip_Chrysoperla_nipponensis_NC_015093", "XXXXXXXXXXXXXXX----_NR4Ra7xiz_Rapisma_xizangense_NC_023362", "XXXXXXXXXXXXXXX----_NO4Th7lan_Thyridosmylus_langii_NC_021415", "XXXXXXXXXXXXXXX----_NMy4My7imm_Myrmeleon_immanis_NC_024826", "XXXXXXXXXXXXXXX----_NA4As7app_Ascaloptynx_appendiculatus_NC_011277", "XXXXXXXXXXXXXXX----_NPo4Po7pun_Polystoechotes_punctatus_NC_011278", "XXXXXXXXXXXXXXX----_NMa4Di7bis_Ditaxis_biseriata_NC_013257"))

#for incidence
list.loss.mito.inci=list()
#0 set:whole dataset
entire.loss.mito.inci<-list.cbind(loss.mito.inci)
list.loss.mito.inci[1]<-unlist(list(list(entire.loss.mito.inci),list.loss.mito.inci), recursive=FALSE)
names(list.loss.mito.inci)[1] <- paste0("Whole Dataset")

#1 loss set:South_Summer
minus_ss.loss.mito.inci<-loss.mito.inci
minus_ss.loss.mito.inci["South_Summer"]<-NULL
minus_ss.loss.mito.inci<-list.cbind(minus_ss.loss.mito.inci)
list.loss.mito.inci[2]<-unlist(list(list(minus_ss.loss.mito.inci),list.loss.mito.inci), recursive=FALSE)
names(list.loss.mito.inci)[2] <- paste0("South Summer Loss")

#2 loss set:North_Winter
minus_nw.loss.mito.inci<-loss.mito.inci
minus_nw.loss.mito.inci["North_Winter"]<-NULL
minus_nw.loss.mito.inci<-list.cbind(minus_nw.loss.mito.inci)
list.loss.mito.inci[3]<-unlist(list(list(minus_nw.loss.mito.inci),list.loss.mito.inci), recursive=FALSE)
names(list.loss.mito.inci)[3] <- paste0("North Winter Loss")

#3 loss set:North_Summer
minus_ns.loss.mito.inci<-loss.mito.inci
minus_ns.loss.mito.inci["North_Summer"]<-NULL
minus_ns.loss.mito.inci<-list.cbind(minus_ns.loss.mito.inci)
list.loss.mito.inci[4]<-unlist(list(list(minus_ns.loss.mito.inci),list.loss.mito.inci), recursive=FALSE)
names(list.loss.mito.inci)[4] <- paste0("North Summer Loss")

#4 loss set:South_Winter
minus_sw.loss.mito.inci<-loss.mito.inci
minus_sw.loss.mito.inci["South_Winter"]<-NULL
minus_sw.loss.mito.inci<-list.cbind(minus_sw.loss.mito.inci)
list.loss.mito.inci[5]<-unlist(list(list(minus_sw.loss.mito.inci),list.loss.mito.inci), recursive=FALSE)
names(list.loss.mito.inci)[5] <- paste0("South Winter Loss")

#5 loss set:Summer
minus_s.loss.mito.inci<-loss.mito.inci
minus_s.loss.mito.inci["South_Summer"]<-NULL
minus_s.loss.mito.inci["North_Summer"]<-NULL
minus_s.loss.mito.inci<-list.cbind(minus_s.loss.mito.inci)
list.loss.mito.inci[6]<-unlist(list(list(minus_s.loss.mito.inci),list.loss.mito.inci), recursive=FALSE)
names(list.loss.mito.inci)[6] <- paste0("Summer Loss")

#6 loss set:Winter
minus_w.loss.mito.inci<-loss.mito.inci
minus_w.loss.mito.inci["South_Winter"]<-NULL
minus_w.loss.mito.inci["North_Winter"]<-NULL
minus_w.loss.mito.inci<-list.cbind(minus_w.loss.mito.inci)
list.loss.mito.inci[7]<-unlist(list(list(minus_w.loss.mito.inci),list.loss.mito.inci), recursive=FALSE)
names(list.loss.mito.inci)[7] <- paste0("Winter Loss")

#7 loss set:1400
mito_loss_1400 = mito.tree.otu_table[,group[which(grepl("^1800|^DLJ_1800",group$sample)),][which(grepl("_Winter$|Summer$",group[which(grepl("^1800|^DLJ_1800",group$sample)),]$group)),]$sample]
mito_loss_1400[mito_loss_1400>1] <- 1
list.loss.mito.inci[8]<-unlist(list(list(mito_loss_1400),list.loss.mito.inci), recursive=FALSE)
names(list.loss.mito.inci)[8] <- paste0("1400 Loss")

#8 loss set:1800
mito_loss_1800 = mito.tree.otu_table[,group[which(grepl("^1400|^DLJ_1400",group$sample)),][which(grepl("_Winter$|Summer$",group[which(grepl("^1400|^DLJ_1400",group$sample)),]$group)),]$sample]
mito_loss_1800[mito_loss_1800>1] <- 1
list.loss.mito.inci[9]<-unlist(list(list(mito_loss_1800),list.loss.mito.inci), recursive=FALSE)
names(list.loss.mito.inci)[9] <- paste0("1800 Loss")

#9 loss set:South
mito_loss_south = mito.tree.otu_table[,group[which(grepl("^North_",group$group)),][which(grepl("_Winter$|Summer$",group[which(grepl("^North_",group$group)),]$group)),]$sample]
mito_loss_south[mito_loss_south>1] <- 1
list.loss.mito.inci[10]<-unlist(list(list(mito_loss_south),list.loss.mito.inci), recursive=FALSE)
names(list.loss.mito.inci)[10] <- paste0("South Loss")

#10 loss set:north
mito_loss_north = mito.tree.otu_table[,group[which(grepl("^South_",group$group)),][which(grepl("_Winter$|Summer$",group[which(grepl("^South_",group$group)),]$group)),]$sample]
mito_loss_north[mito_loss_north>1] <- 1
list.loss.mito.inci[11]<-unlist(list(list(mito_loss_north),list.loss.mito.inci), recursive=FALSE)
names(list.loss.mito.inci)[11] <- paste0("North Loss")

mito.inci.loss.set<-data.frame(cbind(list.loss.mito.inci$`Whole Dataset`,list.loss.mito.inci$`South Summer Loss`,list.loss.mito.inci$`North Winter Loss`,list.loss.mito.inci$`North Summer Loss`,list.loss.mito.inci$`South Winter Loss`,list.loss.mito.inci$`Summer Loss`,list.loss.mito.inci$`Winter Loss`,list.loss.mito.inci$`1400 Loss`,list.loss.mito.inci$`1800 Loss`,list.loss.mito.inci$`South Loss`,list.loss.mito.inci$`North Loss`))
nt<-c("Whole Dataset"=ncol(list.loss.mito.inci$`Whole Dataset`),"South Summer Loss"=ncol(list.loss.mito.inci$`South Summer Loss`),"North Winter Loss"=ncol(list.loss.mito.inci$`North Winter Loss`),"North Summer Loss"=ncol(list.loss.mito.inci$`North Summer Loss`),"South Winter Loss"=ncol(list.loss.mito.inci$`South Winter Loss`),"Summer Loss"=ncol(list.loss.mito.inci$`Summer Loss`),"Winter Loss"=ncol(list.loss.mito.inci$`Winter Loss`),"1400 Loss"=ncol(list.loss.mito.inci$`1400 Loss`),"1800 Loss"=ncol(list.loss.mito.inci$`1800 Loss`),"South Loss"=ncol(list.loss.mito.inci$`South Loss`),"North Loss"=ncol(list.loss.mito.inci$`North Loss`))

# Hill numbers (q):  0 = sp richness, 1 = Shannon, 2 = inverse Simpson
pd.mito.inci <- iNEXT3D(mito.inci.loss.set, nT=nt, class = 'PD', tree = mito.treepl, q=c(0,1,2), datatype="incidence_raw")
save(pd.mito.inci,file="pd.inci_mito_treepl_20211125") 

##################
#plot the results of incidence-based rarefaction and extrapolation 
##################

cols <- c("#1d953f","#1d953f","#1d953f","black","black","black","#8a5656","#8a5656","#8a5656","#D55E00","#D55E00")
load("Data/p.loss.inci_OTU_20211214")
p.loss.inci<-ggiNEXT3D(loss.out.inci, type=1, facet.var = "Order.q", color.var ="None")[[1]]+scale_shape_manual(name="Experimental",values=c(1,19,0,0,1,19,2,10,13,7,12),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+scale_colour_manual(name="Experimental",values = alpha(cols,0.5), aesthetics = c("colour", "fill"),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+theme(plot.title = element_text(hjust = 0.5))+theme_bw()
p.loss.inci

load("Data/pd.inci_denovo_treepl_20211125")
pd.unifrac.inci<-ggiNEXT3D(pd.unifrac.inci, type=1, facet.var = "Order.q", color.var ="None")[[1]]+scale_shape_manual(name="Experimental",values=c(1,19,0,0,1,19,2,10,13,7,12),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+scale_colour_manual(name="Experimental",values = alpha(cols,0.5), aesthetics = c("colour", "fill"),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+theme(plot.title = element_text(hjust = 0.5))+theme_bw()
pd.unifrac.inci

load("Data/pd.inci_mito_treepl_20211125")
pd.mito.inci<-ggiNEXT3D(pd.mito.inci, type=1, facet.var = "Order.q", color.var ="None")[[1]]+scale_shape_manual(name="Experimental",values=c(1,19,0,0,1,19,2,10,13,7,12),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+scale_colour_manual(name="Experimental",values = alpha(cols,0.5), aesthetics = c("colour", "fill"),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+theme(plot.title = element_text(hjust = 0.5))+theme_bw()
pd.mito.inci

dev.new()
pdf("Fig.S5.inext3d_inci_20211127.pdf",width=14,height=18)
p.loss.inci/pd.unifrac.inci/pd.mito.inci
dev.off()