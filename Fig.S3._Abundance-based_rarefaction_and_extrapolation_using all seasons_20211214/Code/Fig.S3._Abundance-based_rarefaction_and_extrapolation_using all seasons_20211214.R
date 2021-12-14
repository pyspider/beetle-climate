###############################################################
##############alpha diversity using all seasons################
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
loss_V2=list()
loss.inci=list()

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
# loss.abun["South_Spring"]<-NULL
# loss.abun["North_Spring"]<-NULL
# loss.abun["North_Autumn"]<-NULL
# loss.abun["South_Autumn"]<-NULL
# 
# loss.inci["South_Spring"]<-NULL
# loss.inci["North_Spring"]<-NULL
# loss.inci["North_Autumn"]<-NULL
# loss.inci["South_Autumn"]<-NULL


#0 set:whole dataset
entire.loss.abun<-apply(list.rbind(loss.abun),2,sum)
sum(entire.loss.abun)
loss_V2[1]<-unlist(list(list(entire.loss.abun),loss_V2), recursive=FALSE)
names(loss_V2)[1] <- paste0("Whole Dataset")

#1 loss set:South_Summer
minus_ss.loss.abun<-loss.abun
minus_ss.loss.abun["South_Summer"]<-NULL
minus_ss.loss.abun<-apply(list.rbind(minus_ss.loss.abun),2,sum)
#sum(minus_ss.loss.abun)
loss_V2[2]<-unlist(list(list(minus_ss.loss.abun),loss_V2), recursive=FALSE)
names(loss_V2)[2] <- paste0("South Summer Loss")

#2 loss set:North_Winter
minus_nw.loss.abun<-loss.abun
minus_nw.loss.abun["North_Winter"]<-NULL
minus_nw.loss.abun<-apply(list.rbind(minus_nw.loss.abun),2,sum)
loss_V2[3]<-unlist(list(list(minus_nw.loss.abun),loss_V2), recursive=FALSE)
names(loss_V2)[3] <- paste0("North Winter Loss")

#3 loss set:North_Summer
minus_ns.loss.abun<-loss.abun
minus_ns.loss.abun["North_Summer"]<-NULL
minus_ns.loss.abun<-apply(list.rbind(minus_ns.loss.abun),2,sum)
loss_V2[4]<-unlist(list(list(minus_ns.loss.abun),loss_V2), recursive=FALSE)
names(loss_V2)[4] <- paste0("North Summer Loss")

#4 loss set:South_Winter
minus_sw.loss.abun<-loss.abun
minus_sw.loss.abun["South_Winter"]<-NULL
minus_sw.loss.abun<-apply(list.rbind(minus_sw.loss.abun),2,sum)
loss_V2[5]<-unlist(list(list(minus_sw.loss.abun),loss_V2), recursive=FALSE)
names(loss_V2)[5] <- paste0("South Winter Loss")

#5 loss set:Summer
minus_s.loss.abun<-loss.abun
minus_s.loss.abun["South_Summer"]<-NULL
minus_s.loss.abun["North_Summer"]<-NULL
minus_s.loss.abun<-apply(list.rbind(minus_s.loss.abun),2,sum)
#sum(minus_s.loss.abun)
loss_V2[6]<-unlist(list(list(minus_s.loss.abun),loss_V2), recursive=FALSE)
names(loss_V2)[6] <- paste0("Summer Loss")

#6 loss set:Winter
minus_w.loss.abun<-loss.abun
minus_w.loss.abun["South_Winter"]<-NULL
minus_w.loss.abun["North_Winter"]<-NULL
minus_w.loss.abun<-apply(list.rbind(minus_w.loss.abun),2,sum)
#sum(minus_w.loss.abun)
loss_V2[7]<-unlist(list(list(minus_w.loss.abun),loss_V2), recursive=FALSE)
names(loss_V2)[7] <- paste0("Winter Loss")

#7 loss set:1400
otu_loss_1400 = t(OTU.data_biweek.pivot_fixed[group[which(grepl("^1800|^DLJ_1800",group$sample)),]$sample, ])
temp_loss_1400<-apply(otu_loss_1400,1,sum)
loss_V2[8]<-unlist(list(list(temp_loss_1400),loss_V2), recursive=FALSE)
names(loss_V2)[8] <- paste0("1400 Loss")

#8 loss set:1800
otu_loss_1800 = t(OTU.data_biweek.pivot_fixed[group[which(grepl("^1400|^DLJ_1400",group$sample)),]$sample, ])
temp_loss_1800<-apply(otu_loss_1800,1,sum)
loss_V2[9]<-unlist(list(list(temp_loss_1800),loss_V2), recursive=FALSE)
names(loss_V2)[9] <- paste0("1800 Loss")

#9 loss set:South'
otu_loss_south = t(OTU.data_biweek.pivot_fixed[group[which(grepl("^North_",group$group)),]$sample, ])
temp_loss_south<-apply(otu_loss_south,1,sum)
loss_V2[10]<-unlist(list(list(temp_loss_south),loss_V2), recursive=FALSE)
names(loss_V2)[10] <- paste0("south Loss")

#10 loss set:North
otu_loss_north = t(OTU.data_biweek.pivot_fixed[group[which(grepl("^South_",group$group)),]$sample, ])
temp_loss_north<-apply(otu_loss_north,1,sum)
loss_V2[11]<-unlist(list(list(temp_loss_north),loss_V2), recursive=FALSE)
names(loss_V2)[11] <- paste0("north Loss")

loss.abun.data<-data.frame(cbind(loss_V2$`Whole Dataset`,loss_V2$`South Summer Loss`,loss_V2$`North Winter Loss`,loss_V2$`North Summer Loss`,loss_V2$`South Winter Loss`,loss_V2$`Summer Loss`,loss_V2$`Winter Loss`,loss_V2$`1400 Loss`,loss_V2$`1800 Loss`,loss_V2$`south Loss`,loss_V2$`north Loss`))
colnames(loss.abun.data)<-c("Whole Dataset","South Summer Loss","North Winter Loss","North Summer Loss","South Winter Loss","Summer Loss","Winter Loss","1400 Loss","1800 Loss","South Loss","North Loss")
#loss.abun.data.filter<-loss.abun.data[which(rowSums(loss.abun.data)>0),]
apply(loss.abun.data,2,sum)
# Whole Dataset South Summer Loss North Winter Loss North Summer Loss South Winter Loss       Summer Loss       Winter Loss         1400 Loss 
# 10258              5965             10094              8790              9950              4497              9786              2879 
# 1800 Loss        South Loss        North Loss 
# 7379              3182              7076 

# Hill numbers (q):  0 = sp richness, 1 = Shannon, 2 = inverse Simpson
otu.abun.whole <- iNEXT3D(loss.abun.data, class = 'TD', q = c(0,1,2), datatype = "abundance")
save(otu.abun.whole,file="pd.otu_whole_20211214")


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
#ultrametric tree
iqtree.treepl<-read.tree(file="Data/denovo.treePL")
# count zero length branches
sum(iqtree.treepl$edge.length == 0)
#[1] 0
is.ultrametric(iqtree.treepl)
#iqtree.treepl<-force.ultrametric(iqtree.treepl)#default method: "nnls"
iqtree.treepl<-force.ultrametric(iqtree.treepl,method=c("extend"))
is.ultrametric(iqtree.treepl)

##Assemble process ultrametric tree
iqtree.treepl<-drop.tip(iqtree.treepl,c("XXXXXXXXXXXXXXX----_NNy4Ny7myr_Nymphes_myrmeleonoides_NC_024825", "XXXXXXXXXXXXXXX----_NCh4Chr7nip_Chrysoperla_nipponensis_NC_015093", "XXXXXXXXXXXXXXX----_NR4Ra7xiz_Rapisma_xizangense_NC_023362", "XXXXXXXXXXXXXXX----_NO4Th7lan_Thyridosmylus_langii_NC_021415", "XXXXXXXXXXXXXXX----_NMy4My7imm_Myrmeleon_immanis_NC_024826", "XXXXXXXXXXXXXXX----_NA4As7app_Ascaloptynx_appendiculatus_NC_011277", "XXXXXXXXXXXXXXX----_NPo4Po7pun_Polystoechotes_punctatus_NC_011278", "XXXXXXXXXXXXXXX----_NMa4Di7bis_Ditaxis_biseriata_NC_013257"))

#phylo tree contains duplicated node labels, remove them
iqtree.treepl$node.label<-NULL

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

# loss.unifrac.abun["South_Spring"]<-NULL
# loss.unifrac.abun["North_Spring"]<-NULL
# loss.unifrac.abun["North_Autumn"]<-NULL
# loss.unifrac.abun["South_Autumn"]<-NULL
# 
# loss.unifrac.inci["South_Spring"]<-NULL
# loss.unifrac.inci["North_Spring"]<-NULL
# loss.unifrac.inci["North_Autumn"]<-NULL
# loss.unifrac.inci["South_Autumn"]<-NULL


list.loss=list()
#0 set:whole dataset
entire.loss.abun<-list.cbind(loss.unifrac.abun)
list.loss[1]<-unlist(list(list(apply(entire.loss.abun,1,sum)),list.loss), recursive=FALSE)
names(list.loss)[1] <- paste0("Whole Dataset")

#1 loss set:South_Summer
minus_ss.loss.abun<-loss.unifrac.abun
minus_ss.loss.abun["South_Summer"]<-NULL
minus_ss.loss.abun<-list.cbind(minus_ss.loss.abun)
list.loss[2]<-unlist(list(list(apply(minus_ss.loss.abun,1,sum)),list.loss), recursive=FALSE)
names(list.loss)[2] <- paste0("South Summer Loss")

#2 loss set:North_Winter
minus_nw.loss.abun<-loss.unifrac.abun
minus_nw.loss.abun["North_Winter"]<-NULL
minus_nw.loss.abun<-list.cbind(minus_nw.loss.abun)
list.loss[3]<-unlist(list(list(apply(minus_nw.loss.abun,1,sum)),list.loss), recursive=FALSE)
names(list.loss)[3] <- paste0("North Winter Loss")

#3 loss set:North_Summer
minus_ns.loss.abun<-loss.unifrac.abun
minus_ns.loss.abun["North_Summer"]<-NULL
minus_ns.loss.abun<-list.cbind(minus_ns.loss.abun)
list.loss[4]<-unlist(list(list(apply(minus_ns.loss.abun,1,sum)),list.loss), recursive=FALSE)
names(list.loss)[4] <- paste0("North Summer Loss")

#4 loss set:South_Winter
minus_sw.loss.abun<-loss.unifrac.abun
minus_sw.loss.abun["South_Winter"]<-NULL
minus_sw.loss.abun<-list.cbind(minus_sw.loss.abun)
list.loss[5]<-unlist(list(list(apply(minus_sw.loss.abun,1,sum)),list.loss), recursive=FALSE)
names(list.loss)[5] <- paste0("South Winter Loss")

#5 loss set:Summer
minus_s.loss.abun<-loss.unifrac.abun
minus_s.loss.abun["South_Summer"]<-NULL
minus_s.loss.abun["North_Summer"]<-NULL
minus_s.loss.abun<-list.cbind(minus_s.loss.abun)
list.loss[6]<-unlist(list(list(apply(minus_s.loss.abun,1,sum)),list.loss), recursive=FALSE)
names(list.loss)[6] <- paste0("Summer Loss")

#6 loss set:Winter
minus_w.loss.abun<-loss.unifrac.abun
minus_w.loss.abun["South_Winter"]<-NULL
minus_w.loss.abun["North_Winter"]<-NULL
minus_w.loss.abun<-list.cbind(minus_w.loss.abun)
list.loss[7]<-unlist(list(list(apply(minus_w.loss.abun,1,sum)),list.loss), recursive=FALSE)
names(list.loss)[7] <- paste0("Winter Loss")

#7 loss set:1400
iqtree_loss_1400 = iqtree.otu_table[,group[which(grepl("^1800|^DLJ_1800",group$sample)),]$sample]
temp_loss_1400<-apply(iqtree_loss_1400,1,sum)
list.loss[8]<-unlist(list(list(temp_loss_1400),list.loss), recursive=FALSE)
names(list.loss)[8] <- paste0("1400 Loss")

#8 loss set:1800
iqtree_loss_1800 = iqtree.otu_table[,group[which(grepl("^1400|^DLJ_1400",group$sample)),]$sample]
temp_loss_1800<-apply(iqtree_loss_1800,1,sum)
list.loss[9]<-unlist(list(list(temp_loss_1800),list.loss), recursive=FALSE)
names(list.loss)[9] <- paste0("1800 Loss")

#9 loss set:South
iqtree_loss_south = iqtree.otu_table[,group[which(grepl("^North_",group$group)),]$sample]
temp_loss_south<-apply(iqtree_loss_south,1,sum)
list.loss[10]<-unlist(list(list(temp_loss_south),list.loss), recursive=FALSE)
names(list.loss)[10] <- paste0("South Loss")

#10 loss set:north
iqtree_loss_north = iqtree.otu_table[,group[which(grepl("^South_",group$group)),]$sample]
temp_loss_north<-apply(iqtree_loss_north,1,sum)
list.loss[11]<-unlist(list(list(temp_loss_north),list.loss), recursive=FALSE)
names(list.loss)[11] <- paste0("North Loss")

unifrac.abun.loss.set<-data.frame(cbind(list.loss$`Whole Dataset`,list.loss$`South Summer Loss`,list.loss$`North Winter Loss`,list.loss$`North Summer Loss`,list.loss$`South Winter Loss`,list.loss$`Summer Loss`,list.loss$`Winter Loss`,list.loss$`1400 Loss`,list.loss$`1800 Loss`,list.loss$`South Loss`,list.loss$`North Loss`))
colnames(unifrac.abun.loss.set)<-c("Whole Dataset","South Summer Loss","North Winter Loss","North Summer Loss","South Winter Loss","Summer Loss","Winter Loss","1400 Loss","1800 Loss","South Loss","North Loss")
apply(unifrac.abun.loss.set,2,sum)
# Whole Dataset South Summer Loss North Winter Loss North Summer Loss South Winter Loss       Summer Loss       Winter Loss         1400 Loss 
# 10524              6164             10360              9055             10173              4695             10009              2949 
# 1800 Loss        South Loss        North Loss 
# 7575              3183              7341 

# Hill numbers (q):  0 = sp richness, 1 = Shannon, 2 = inverse Simpson
pd.abun <- iNEXT3D(unifrac.abun.loss.set, class = 'PD', tree = iqtree.treepl, datatype = "abundance", q = c(0, 1, 2))
save(pd.abun,file="pd.abun_denovo_treepl_whole_20211125")

####################################
#CASE C for placement tree
####################################
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


# loss.mito.abun["South_Spring"]<-NULL
# loss.mito.abun["North_Spring"]<-NULL
# loss.mito.abun["North_Autumn"]<-NULL
# loss.mito.abun["South_Autumn"]<-NULL
# 
# loss.mito.inci["South_Spring"]<-NULL
# loss.mito.inci["North_Spring"]<-NULL
# loss.mito.inci["North_Autumn"]<-NULL
# loss.mito.inci["South_Autumn"]<-NULL

list.loss=list()
#0 set:whole dataset
entire.loss.abun<-list.cbind(loss.mito.abun)
list.loss[1]<-unlist(list(list(apply(entire.loss.abun,1,sum)),list.loss), recursive=FALSE)
names(list.loss)[1] <- paste0("Whole Dataset")

#1 loss set:South_Summer
minus_ss.loss.abun<-loss.mito.abun
minus_ss.loss.abun["South_Summer"]<-NULL
minus_ss.loss.abun<-list.cbind(minus_ss.loss.abun)
list.loss[2]<-unlist(list(list(apply(minus_ss.loss.abun,1,sum)),list.loss), recursive=FALSE)
names(list.loss)[2] <- paste0("South Summer Loss")

#2 loss set:North_Winter
minus_nw.loss.abun<-loss.mito.abun
minus_nw.loss.abun["North_Winter"]<-NULL
minus_nw.loss.abun<-list.cbind(minus_nw.loss.abun)
list.loss[3]<-unlist(list(list(apply(minus_nw.loss.abun,1,sum)),list.loss), recursive=FALSE)
names(list.loss)[3] <- paste0("North Winter Loss")

#3 loss set:North_Summer
minus_ns.loss.abun<-loss.mito.abun
minus_ns.loss.abun["North_Summer"]<-NULL
minus_ns.loss.abun<-list.cbind(minus_ns.loss.abun)
list.loss[4]<-unlist(list(list(apply(minus_ns.loss.abun,1,sum)),list.loss), recursive=FALSE)
names(list.loss)[4] <- paste0("North Summer Loss")

#4 loss set:South_Winter
minus_sw.loss.abun<-loss.mito.abun
minus_sw.loss.abun["South_Winter"]<-NULL
minus_sw.loss.abun<-list.cbind(minus_sw.loss.abun)
list.loss[5]<-unlist(list(list(apply(minus_sw.loss.abun,1,sum)),list.loss), recursive=FALSE)
names(list.loss)[5] <- paste0("South Winter Loss")

#5 loss set:Summer
minus_s.loss.abun<-loss.mito.abun
minus_s.loss.abun["South_Summer"]<-NULL
minus_s.loss.abun["North_Summer"]<-NULL
minus_s.loss.abun<-list.cbind(minus_s.loss.abun)
list.loss[6]<-unlist(list(list(apply(minus_s.loss.abun,1,sum)),list.loss), recursive=FALSE)
names(list.loss)[6] <- paste0("Summer Loss")

#6 loss set:Winter
minus_w.loss.abun<-loss.mito.abun
minus_w.loss.abun["South_Winter"]<-NULL
minus_w.loss.abun["North_Winter"]<-NULL
minus_w.loss.abun<-list.cbind(minus_w.loss.abun)
list.loss[7]<-unlist(list(list(apply(minus_w.loss.abun,1,sum)),list.loss), recursive=FALSE)
names(list.loss)[7] <- paste0("Winter Loss")

#7 loss set:1400
mito_loss_1400 = mito.tree.otu_table[,group[which(grepl("^1800|^DLJ_1800",group$sample)),]$sample]
temp_loss_1400<-apply(mito_loss_1400,1,sum)
list.loss[8]<-unlist(list(list(temp_loss_1400),list.loss), recursive=FALSE)
names(list.loss)[8] <- paste0("1400 Loss")

#8 loss set:1800
mito_loss_1800 = mito.tree.otu_table[,group[which(grepl("^1400|^DLJ_1400",group$sample)),]$sample]
temp_loss_1800<-apply(mito_loss_1800,1,sum)
list.loss[9]<-unlist(list(list(temp_loss_1800),list.loss), recursive=FALSE)
names(list.loss)[9] <- paste0("1800 Loss")

#9 loss set:South
mito_loss_south = mito.tree.otu_table[,group[which(grepl("^North_",group$group)),]$sample]
temp_loss_south<-apply(mito_loss_south,1,sum)
list.loss[10]<-unlist(list(list(temp_loss_south),list.loss), recursive=FALSE)
names(list.loss)[10] <- paste0("South Loss")

#10 loss set:north
mito_loss_north = mito.tree.otu_table[,group[which(grepl("^South_",group$group)),]$sample]
temp_loss_north<-apply(mito_loss_north,1,sum)
list.loss[11]<-unlist(list(list(temp_loss_north),list.loss), recursive=FALSE)
names(list.loss)[11] <- paste0("North Loss")

mito.abun.loss.set<-data.frame(cbind(list.loss$`Whole Dataset`,list.loss$`South Summer Loss`,list.loss$`North Winter Loss`,list.loss$`North Summer Loss`,list.loss$`South Winter Loss`,list.loss$`Summer Loss`,list.loss$`Winter Loss`,list.loss$`1400 Loss`,list.loss$`1800 Loss`,list.loss$`South Loss`,list.loss$`North Loss`))
colnames(mito.abun.loss.set)<-c("Whole Dataset","South Summer Loss","North Winter Loss","North Summer Loss","South Winter Loss","Summer Loss","Winter Loss","1400 Loss","1800 Loss","South Loss","North Loss")

apply(mito.abun.loss.set,2,sum)
# Whole Dataset South Summer Loss North Winter Loss North Summer Loss South Winter Loss       Summer Loss       Winter Loss         1400 Loss 
# 10524              6164             10360              9055             10173              4695             10009              2949 
# 1800 Loss        South Loss        North Loss 
# 7575              3183              7341 

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
if (sort(rownames(mito.abun.loss.set))==sort(mito.treepl$tip.label)){paste0("tip label:passed")}

# Hill numbers (q):  0 = sp richness, 1 = Shannon, 2 = inverse Simpson
pd.abun <- iNEXT3D(mito.abun.loss.set, class = 'PD', tree = mito.treepl, datatype = "abundance", q = c(0, 1, 2))
save(pd.abun,file="pd.abun_mito_whole_20211125")

##################
#plot the results of Abundance-based rarefaction and extrapolation using all seasons 
##################

#cols <- c("#1d953f","#1d953f","#1d953f","black","black","black","#8a5656","#8a5656","#8a5656","#896a45","#896a45")

load("Data/pd.otu_whole_20211214")
#p.loss.abun.whole<-ggiNEXT3D(otu.abun.whole, type=1, facet.var = "Order.q", color.var ="None")[[1]]+scale_shape_manual(name="Experimental",values=c(1,19,0,0,1,19,2,10,13,7,12),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+scale_colour_manual(name="Experimental",values = alpha(cols,0.5), aesthetics = c("colour", "fill"),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+theme(plot.title = element_text(hjust = 0.5))+theme_bw()
#p.loss.abun.whole

otu.abun.whole<-otu.abun.whole$iNextEst$size_based
df.point <- otu.abun.whole[which(otu.abun.whole$Order.q=="0"),]
df.point <- df.point[which(df.point$Method=="Observed"),]
df.line <- otu.abun.whole[which(otu.abun.whole$Order.q=="0"),]
df.line.Extrapolation <- otu.abun.whole[which(otu.abun.whole$Order.q=="0"&otu.abun.whole$Method=="Extrapolation"),]
cols <- c("#1d953f","#1d953f","#1d953f","black","black","black","#8a5656","#8a5656","#8a5656","#D55E00","#D55E00")

p0.otu.whole<-ggplot(filter(otu.abun.whole,otu.abun.whole$Method!="Extrapolation"), aes(x=m, y=qD, colour=Assemblage))  +
  geom_point(aes(shape=Assemblage), size=5, data=df.point) +
  geom_line(aes(linetype=c("1")), lwd=1.5,data =filter(df.line,df.line$Method!="Extrapolation")) +
  geom_line(aes(linetype=c("2")), lwd=1.5,data =df.line.Extrapolation) +
  geom_ribbon(aes(ymin=qD.LCL, ymax=qD.UCL,
                  fill=Assemblage, colour=NULL), data = df.line,alpha=0.2)+
  theme(legend.position = "bottom", 
        legend.title=element_blank(),
        text=element_text(size=18),
        legend.box = "vertical")+scale_shape_manual(name="Experimental",values=c(1,19,0,0,1,19,2,10,13,7,12),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+scale_colour_manual(name="Experimental",values = alpha(cols,0.5), aesthetics = c("colour", "fill"),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+theme(plot.title = element_text(hjust = 0.5))+theme_bw()+theme(legend.position = "none")

df.point <- otu.abun.whole[which(otu.abun.whole$Order.q=="1"),]
df.point <- df.point[which(df.point$Method=="Observed"),]
df.line <- otu.abun.whole[which(otu.abun.whole$Order.q=="1"),]
df.line.Extrapolation <- otu.abun.whole[which(otu.abun.whole$Order.q=="1"&otu.abun.whole$Method=="Extrapolation"),]

p1.otu.whole<-ggplot(filter(otu.abun.whole,otu.abun.whole$Method!="Extrapolation"), aes(x=m, y=qD, colour=Assemblage))  +
  geom_point(aes(shape=Assemblage), size=5, data=df.point) +
  geom_line(aes(linetype=c("1")), lwd=1.5,data =filter(df.line,df.line$Method!="Extrapolation")) +
  geom_line(aes(linetype=c("2")), lwd=1.5,data =df.line.Extrapolation) +
  geom_ribbon(aes(ymin=qD.LCL, ymax=qD.UCL,
                  fill=Assemblage, colour=NULL), data = df.line,alpha=0.2)+
  theme(legend.position = "bottom", 
        legend.title=element_blank(),
        text=element_text(size=18),
        legend.box = "vertical")+scale_shape_manual(name="Experimental",values=c(1,19,0,0,1,19,2,10,13,7,12),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+scale_colour_manual(name="Experimental",values = alpha(cols,0.5), aesthetics = c("colour", "fill"),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+theme(plot.title = element_text(hjust = 0.5))+theme_bw()+theme(legend.position = "none")

df.point <- otu.abun.whole[which(otu.abun.whole$Order.q=="2"),]
df.point <- df.point[which(df.point$Method=="Observed"),]
df.line <- otu.abun.whole[which(otu.abun.whole$Order.q=="2"),]
df.line.Extrapolation <- otu.abun.whole[which(otu.abun.whole$Order.q=="2"&otu.abun.whole$Method=="Extrapolation"),]

p2.otu.whole<-ggplot(filter(otu.abun.whole,otu.abun.whole$Method!="Extrapolation"), aes(x=m, y=qD, colour=Assemblage))  +
  geom_point(aes(shape=Assemblage), size=5, data=df.point) +
  geom_line(aes(linetype=c("1")), lwd=1.5,data =filter(df.line,df.line$Method!="Extrapolation")) +
  geom_line(aes(linetype=c("2")), lwd=1.5,data =df.line.Extrapolation) +
  geom_ribbon(aes(ymin=qD.LCL, ymax=qD.UCL,
                  fill=Assemblage, colour=NULL), data = df.line,alpha=0.2)+
  theme(legend.position = "bottom", 
        legend.title=element_blank(),
        text=element_text(size=18),
        legend.box = "vertical")+scale_shape_manual(name="Experimental",values=c(1,19,0,0,1,19,2,10,13,7,12),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+scale_colour_manual(name="Experimental",values = alpha(cols,0.5), aesthetics = c("colour", "fill"),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+theme(plot.title = element_text(hjust = 0.5))+theme_bw()

p.otu.whole <- p0.otu.whole+p1.otu.whole+p2.otu.whole

load("Data/pd.abun_denovo_treepl_whole_20211125")
#pd.unifrac.abun.whole<-ggiNEXT3D(pd.abun, type=1, facet.var = "Order.q", color.var ="None")[[1]]+scale_shape_manual(name="Experimental",values=c(1,19,0,0,1,19,2,10,13,7,12),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+scale_colour_manual(name="Experimental",values = alpha(cols,0.5), aesthetics = c("colour", "fill"),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+theme(plot.title = element_text(hjust = 0.5))+theme_bw()
#pd.unifrac.abun.whole

pd.abun.denovo.whole<-pd.abun$PDiNextEst$size_based
df.point <- pd.abun.denovo.whole[which(pd.abun.denovo.whole$Order.q=="0"),]
df.point <- df.point[which(df.point$Method=="Observed"),]
df.line <- pd.abun.denovo.whole[which(pd.abun.denovo.whole$Order.q=="0"),]
df.line.Extrapolation <- pd.abun.denovo.whole[which(pd.abun.denovo.whole$Order.q=="0"&pd.abun.denovo.whole$Method=="Extrapolation"),]
cols <- c("#1d953f","#1d953f","#1d953f","black","black","black","#8a5656","#8a5656","#8a5656","#D55E00","#D55E00")

p0.denovo.whole<-ggplot(filter(pd.abun.denovo.whole,pd.abun.denovo.whole$Method!="Extrapolation"), aes(x=m, y=qPD, colour=Assemblage))  +
  geom_point(aes(shape=Assemblage), size=5, data=df.point) +
  geom_line(aes(linetype=c("1")), lwd=1.5,data =filter(df.line,df.line$Method!="Extrapolation")) +
  geom_line(aes(linetype=c("2")), lwd=1.5,data =df.line.Extrapolation) +
  geom_ribbon(aes(ymin=qPD.LCL, ymax=qPD.UCL,
                  fill=Assemblage, colour=NULL), data = df.line,alpha=0.2)+
  theme(legend.position = "bottom", 
        legend.title=element_blank(),
        text=element_text(size=18),
        legend.box = "vertical")+scale_shape_manual(name="Experimental",values=c(1,19,0,0,1,19,2,10,13,7,12),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+scale_colour_manual(name="Experimental",values = alpha(cols,0.5), aesthetics = c("colour", "fill"),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+theme(plot.title = element_text(hjust = 0.5))+theme_bw()+theme(legend.position = "none")

df.point <- pd.abun.denovo.whole[which(pd.abun.denovo.whole$Order.q=="1"),]
df.point <- df.point[which(df.point$Method=="Observed"),]
df.line <- pd.abun.denovo.whole[which(pd.abun.denovo.whole$Order.q=="1"),]
df.line.Extrapolation <- pd.abun.denovo.whole[which(pd.abun.denovo.whole$Order.q=="1"&pd.abun.denovo.whole$Method=="Extrapolation"),]
cols <- c("#1d953f","#1d953f","#1d953f","black","black","black","#8a5656","#8a5656","#8a5656","#D55E00","#D55E00")

p1.denovo.whole<-ggplot(filter(pd.abun.denovo.whole,pd.abun.denovo.whole$Method!="Extrapolation"), aes(x=m, y=qPD, colour=Assemblage))  +
  geom_point(aes(shape=Assemblage), size=5, data=df.point) +
  geom_line(aes(linetype=c("1")), lwd=1.5,data =filter(df.line,df.line$Method!="Extrapolation")) +
  geom_line(aes(linetype=c("2")), lwd=1.5,data =df.line.Extrapolation) +
  geom_ribbon(aes(ymin=qPD.LCL, ymax=qPD.UCL,
                  fill=Assemblage, colour=NULL), data = df.line,alpha=0.2)+
  theme(legend.position = "bottom", 
        legend.title=element_blank(),
        text=element_text(size=18),
        legend.box = "vertical")+scale_shape_manual(name="Experimental",values=c(1,19,0,0,1,19,2,10,13,7,12),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+scale_colour_manual(name="Experimental",values = alpha(cols,0.5), aesthetics = c("colour", "fill"),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+theme(plot.title = element_text(hjust = 0.5))+theme_bw()+theme(legend.position = "none")

df.point <- pd.abun.denovo.whole[which(pd.abun.denovo.whole$Order.q=="2"),]
df.point <- df.point[which(df.point$Method=="Observed"),]
df.line <- pd.abun.denovo.whole[which(pd.abun.denovo.whole$Order.q=="2"),]
df.line.Extrapolation <- pd.abun.denovo.whole[which(pd.abun.denovo.whole$Order.q=="2"&pd.abun.denovo.whole$Method=="Extrapolation"),]
cols <- c("#1d953f","#1d953f","#1d953f","black","black","black","#8a5656","#8a5656","#8a5656","#D55E00","#D55E00")

p2.denovo.whole<-ggplot(filter(pd.abun.denovo.whole,pd.abun.denovo.whole$Method!="Extrapolation"), aes(x=m, y=qPD, colour=Assemblage))  +
  geom_point(aes(shape=Assemblage), size=5, data=df.point) +
  geom_line(aes(linetype=c("1")), lwd=1.5,data =filter(df.line,df.line$Method!="Extrapolation")) +
  geom_line(aes(linetype=c("2")), lwd=1.5,data =df.line.Extrapolation) +
  geom_ribbon(aes(ymin=qPD.LCL, ymax=qPD.UCL,
                  fill=Assemblage, colour=NULL), data = df.line,alpha=0.2)+
  theme(legend.position = "bottom", 
        legend.title=element_blank(),
        text=element_text(size=18),
        legend.box = "vertical")+scale_shape_manual(name="Experimental",values=c(1,19,0,0,1,19,2,10,13,7,12),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+scale_colour_manual(name="Experimental",values = alpha(cols,0.5), aesthetics = c("colour", "fill"),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+theme(plot.title = element_text(hjust = 0.5))+theme_bw()

p.denovo.whole <- p0.denovo.whole+p1.denovo.whole+p2.denovo.whole

load("Data/pd.abun_mito_whole_20211125")
#pd.mito.abun.whole<-ggiNEXT3D(pd.abun, type=1, facet.var = "Order.q", color.var ="None")[[1]]+scale_shape_manual(name="Experimental",values=c(1,19,0,0,1,19,2,10,13,7,12),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+scale_colour_manual(name="Experimental",values = alpha(cols,0.5), aesthetics = c("colour", "fill"),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+theme(plot.title = element_text(hjust = 0.5))+theme_bw()
#pd.mito.abun.whole

pd.abun.mito.whole<-pd.abun$PDiNextEst$size_based
df.point <- pd.abun.mito.whole[which(pd.abun.mito.whole$Order.q=="0"),]
df.point <- df.point[which(df.point$Method=="Observed"),]
df.line <- pd.abun.mito.whole[which(pd.abun.mito.whole$Order.q=="0"),]
df.line.Extrapolation <- pd.abun.mito.whole[which(pd.abun.mito.whole$Order.q=="0"&pd.abun.mito.whole$Method=="Extrapolation"),]
cols <- c("#1d953f","#1d953f","#1d953f","black","black","black","#8a5656","#8a5656","#8a5656","#D55E00","#D55E00")

p0.mito.whole<-ggplot(filter(pd.abun.mito.whole,pd.abun.mito.whole$Method!="Extrapolation"), aes(x=m, y=qPD, colour=Assemblage))  +
  geom_point(aes(shape=Assemblage), size=5, data=df.point) +
  geom_line(aes(linetype=c("1")), lwd=1.5,data =filter(df.line,df.line$Method!="Extrapolation")) +
  geom_line(aes(linetype=c("2")), lwd=1.5,data =df.line.Extrapolation) +
  geom_ribbon(aes(ymin=qPD.LCL, ymax=qPD.UCL,
                  fill=Assemblage, colour=NULL), data = df.line,alpha=0.2)+
  theme(legend.position = "bottom", 
        legend.title=element_blank(),
        text=element_text(size=18),
        legend.box = "vertical")+scale_shape_manual(name="Experimental",values=c(1,19,0,0,1,19,2,10,13,7,12),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+scale_colour_manual(name="Experimental",values = alpha(cols,0.5), aesthetics = c("colour", "fill"),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+theme(plot.title = element_text(hjust = 0.5))+theme_bw()+theme(legend.position = "none")

pd.abun.mito.whole<-pd.abun$PDiNextEst$size_based
df.point <- pd.abun.mito.whole[which(pd.abun.mito.whole$Order.q=="1"),]
df.point <- df.point[which(df.point$Method=="Observed"),]
df.line <- pd.abun.mito.whole[which(pd.abun.mito.whole$Order.q=="1"),]
df.line.Extrapolation <- pd.abun.mito.whole[which(pd.abun.mito.whole$Order.q=="1"&pd.abun.mito.whole$Method=="Extrapolation"),]
cols <- c("#1d953f","#1d953f","#1d953f","black","black","black","#8a5656","#8a5656","#8a5656","#D55E00","#D55E00")

p1.mito.whole<-ggplot(filter(pd.abun.mito.whole,pd.abun.mito.whole$Method!="Extrapolation"), aes(x=m, y=qPD, colour=Assemblage))  +
  geom_point(aes(shape=Assemblage), size=5, data=df.point) +
  geom_line(aes(linetype=c("1")), lwd=1.5,data =filter(df.line,df.line$Method!="Extrapolation")) +
  geom_line(aes(linetype=c("2")), lwd=1.5,data =df.line.Extrapolation) +
  geom_ribbon(aes(ymin=qPD.LCL, ymax=qPD.UCL,
                  fill=Assemblage, colour=NULL), data = df.line,alpha=0.2)+
  theme(legend.position = "bottom", 
        legend.title=element_blank(),
        text=element_text(size=18),
        legend.box = "vertical")+scale_shape_manual(name="Experimental",values=c(1,19,0,0,1,19,2,10,13,7,12),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+scale_colour_manual(name="Experimental",values = alpha(cols,0.5), aesthetics = c("colour", "fill"),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+theme(plot.title = element_text(hjust = 0.5))+theme_bw()+theme(legend.position = "none")

pd.abun.mito.whole<-pd.abun$PDiNextEst$size_based
df.point <- pd.abun.mito.whole[which(pd.abun.mito.whole$Order.q=="2"),]
df.point <- df.point[which(df.point$Method=="Observed"),]
df.line <- pd.abun.mito.whole[which(pd.abun.mito.whole$Order.q=="2"),]
df.line.Extrapolation <- pd.abun.mito.whole[which(pd.abun.mito.whole$Order.q=="2"&pd.abun.mito.whole$Method=="Extrapolation"),]
cols <- c("#1d953f","#1d953f","#1d953f","black","black","black","#8a5656","#8a5656","#8a5656","#D55E00","#D55E00")

p2.mito.whole<-ggplot(filter(pd.abun.mito.whole,pd.abun.mito.whole$Method!="Extrapolation"), aes(x=m, y=qPD, colour=Assemblage))  +
  geom_point(aes(shape=Assemblage), size=5, data=df.point) +
  geom_line(aes(linetype=c("1")), lwd=1.5,data =filter(df.line,df.line$Method!="Extrapolation")) +
  geom_line(aes(linetype=c("2")), lwd=1.5,data =df.line.Extrapolation) +
  geom_ribbon(aes(ymin=qPD.LCL, ymax=qPD.UCL,
                  fill=Assemblage, colour=NULL), data = df.line,alpha=0.2)+
  theme(legend.position = "bottom", 
        legend.title=element_blank(),
        text=element_text(size=18),
        legend.box = "vertical")+scale_shape_manual(name="Experimental",values=c(1,19,0,0,1,19,2,10,13,7,12),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+scale_colour_manual(name="Experimental",values = alpha(cols,0.5), aesthetics = c("colour", "fill"),breaks = c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"))+theme(plot.title = element_text(hjust = 0.5))+theme_bw()+theme(legend.position = "none")

p.mito.whole <- p0.mito.whole+p1.mito.whole+p2.mito.whole

dev.new()
pdf("Fig.S3._Abundance-based_inext3d_uisng_all_seasons_20211215.pdf",width=14,height=18)
p.otu.whole/p.denovo.whole/p.mito.whole
dev.off()
