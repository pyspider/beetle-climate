##################
#Ratio of loss sets, specifying the diversity metric
#################
library(ggplot2) #‘ggplot2’ version 3.3.5
library(dplyr) #‘dplyr’ version 1.0.2
library(patchwork) #‘patchwork’ version 1.0.0

#calculate the proportional decreases in Species richness/Shannon/Simpson and Faith's PD/Phylogenetic entropy/Rao’s quadratic entropy for OTU table, de novo tree, and placement tree
##################
#1. OTU table
##################
#import the inext3d result for Otu table 
load("Data/pd.otu_20211218")
otu.rate<-as.data.frame(otu.abun$iNextEst$size_based)
#using the observed diversity that just at the beginning of the extrapolation
otu.rate_0<-filter(otu.rate,otu.rate$Method=="Observed"&otu.rate$Order.q=="0")
#The ratio of diversity loss was calculated by subtracting the ratio of observed diversity of loss set to the whole dataset from one.
otu.rate_0$ratio<-otu.rate_0$qD/otu.rate_0[which(otu.rate_0$Assemblage=="Whole Dataset"),]$qD
otu.rate_1<-filter(otu.rate,otu.rate$Method=="Observed"&otu.rate$Order.q=="1")
otu.rate_1$ratio<-otu.rate_1$qD/otu.rate_1[which(otu.rate_1$Assemblage=="Whole Dataset"),]$qD
otu.rate_2<-filter(otu.rate,otu.rate$Method=="Observed"&otu.rate$Order.q=="2")
otu.rate_2$ratio<-otu.rate_2$qD/otu.rate_2[which(otu.rate_2$Assemblage=="Whole Dataset"),]$qD
otu.rate<-rbind(otu.rate_0,otu.rate_1,otu.rate_2)
otu.rate$Method<-c("OTU Table")

##################
#2. de novo tree
##################
#import the inext3d result for de novo tree 
load("Data/pd.abun_denovo_treepl_pd_20211219")
denovo.rate<-as.data.frame(pd.abun$PDiNextEst$size_based)
denovo.rate_0<-filter(denovo.rate,denovo.rate$Method=="Observed"&denovo.rate$Order.q=="0")
denovo.rate_0$ratio<-denovo.rate_0$qPD/denovo.rate_0[which(denovo.rate_0$Assemblage=="Whole Dataset"),]$qPD
denovo.rate_1<-filter(denovo.rate,denovo.rate$Method=="Observed"&denovo.rate$Order.q=="1")
denovo.rate_1$ratio<-denovo.rate_1$qPD/denovo.rate_1[which(denovo.rate_1$Assemblage=="Whole Dataset"),]$qPD
denovo.rate_2<-filter(denovo.rate,denovo.rate$Method=="Observed"&denovo.rate$Order.q=="2")
denovo.rate_2$ratio<-denovo.rate_2$qPD/denovo.rate_2[which(denovo.rate_2$Assemblage=="Whole Dataset"),]$qPD
denovo.rate<-rbind(denovo.rate_0,denovo.rate_1,denovo.rate_2)
denovo.rate<-denovo.rate[,-11:-12]
denovo.rate$Method<-c("De novo")

##################
#3. placement tree
##################
#import the inext3d result for placement tree 
load("Data/pd.abun_mito_pd_20211219")
mito.rate<-as.data.frame(pd.abun$PDiNextEst$size_based)
mito.rate_0<-filter(mito.rate,mito.rate$Method=="Observed"&mito.rate$Order.q=="0")
mito.rate_0$ratio<-mito.rate_0$qPD/mito.rate_0[which(mito.rate_0$Assemblage=="Whole Dataset"),]$qPD
mito.rate_1<-filter(mito.rate,mito.rate$Method=="Observed"&mito.rate$Order.q=="1")
mito.rate_1$ratio<-mito.rate_1$qPD/mito.rate_1[which(mito.rate_1$Assemblage=="Whole Dataset"),]$qPD
mito.rate_2<-filter(mito.rate,mito.rate$Method=="Observed"&mito.rate$Order.q=="2")
mito.rate_2$ratio<-mito.rate_2$qPD/mito.rate_2[which(mito.rate_2$Assemblage=="Whole Dataset"),]$qPD
mito.rate<-rbind(mito.rate_0,mito.rate_1,mito.rate_2)
mito.rate<-mito.rate[,-11:-12]
mito.rate$Method<-c("Placement")

rate<-bind_rows(otu.rate,denovo.rate,mito.rate)
rate$number<-paste0(rate$Order.q,"_",rate$Assemblage)
#write.csv(rate,file="rate_of_loss_set_20211211.csv")

#rate$number <- factor(rate$number, levels=c("0_North Summer Loss","0_South Summer Loss","0_Summer Loss","0_Winter Loss","0_North Winter Loss","0_South Winter Loss","0_Whole Dataset","0_1400 Loss","0_1800 Loss","0_South Loss","0_North Loss","1_North Summer Loss","1_South Summer Loss","1_Summer Loss","1_Winter Loss","1_North Winter Loss","1_South Winter Loss","1_Whole Dataset","1_1400 Loss","1_1800 Loss","1_South Loss","1_North Loss","2_North Summer Loss","2_South Summer Loss","2_Summer Loss","2_Winter Loss","2_North Winter Loss","2_South Winter Loss","2_Whole Dataset","2_1400 Loss","2_1800 Loss","2_South Loss","2_North Loss"), ordered=TRUE)

p=ggplot(filter(rate,rate$Method=="OTU Table"),mapping = aes(x=factor(as.character(sub("^[0-2]_","",number)),levels=c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"),ordered=TRUE),y=ratio,group=as.factor(Order.q), color=as.factor(Order.q)))+geom_line()+scale_y_continuous(limits = c(0,1.26))+theme_bw()+ theme(axis.text.x=element_text(angle=90, hjust=1))+labs(x = "")+theme(legend.position = "none")
p1=ggplot(filter(rate,rate$Method=="De novo"),mapping = aes(x=factor(as.character(sub("^[0-2]_","",number)),levels=c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"),ordered=TRUE),y=ratio,group=as.factor(Order.q), color=as.factor(Order.q)))+geom_line()+scale_y_continuous(limits = c(0,1.26))+theme_bw()+ theme(axis.text.x=element_text(angle=90, hjust=1))+labs(x = "")+theme(legend.position = "none")
p2=ggplot(filter(rate,rate$Method=="Placement"),mapping = aes(x=factor(as.character(sub("^[0-2]_","",number)),levels=c("North Summer Loss","South Summer Loss","Summer Loss","Winter Loss","North Winter Loss","South Winter Loss","Whole Dataset","1400 Loss","1800 Loss","South Loss","North Loss"),ordered=TRUE),y=ratio,group=as.factor(Order.q), color=as.factor(Order.q)))+geom_line()+scale_y_continuous(limits = c(0,1.26))+theme_bw()+ theme(axis.text.x=element_text(angle=90, hjust=1))+labs(x = "")+scale_color_discrete(name="Diversity facets",breaks=c("0","1","2"),labels=c("Richness/Faith'PD","Shannon/Phylogenetic entropy","Simpson/Rao’s quadratic entropy"))#+theme(legend.position = "none")

p.loss.ratio<-(p+labs(y="Ratio of diversity loss",title = "Diversity loss on OTU table")+scale_x_discrete("", labels = c("North Summer Loss"="– Summer_North","South Summer Loss"="– Summer_South","Summer Loss"="– Summer","Winter Loss"="– Winter","North Winter Loss"="– Winter_North","South Winter Loss"="– Winter_South","Whole Dataset"="    Whole dataset","1400 Loss"="– 1400 m","1800 Loss"="– 1800 m","South Loss"="– South","North Loss"="– North"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+theme(plot.title = element_text(size=16)))+(p1+labs(y="Ratio of diversity loss",title = "Diversity loss on de novo tree")+scale_x_discrete("", labels = c("North Summer Loss"="– Summer_North","South Summer Loss"="– Summer_South","Summer Loss"="– Summer","Winter Loss"="– Winter","North Winter Loss"="– Winter_North","South Winter Loss"="– Winter_South","Whole Dataset"="    Whole dataset","1400 Loss"="– 1400 m","1800 Loss"="– 1800 m","South Loss"="– South","North Loss"="– North"))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+theme(plot.title = element_text(size=16)))+(p2+labs(y="Ratio of diversity loss",title = "Diversity loss on placement tree")+scale_x_discrete("", labels = c("North Summer Loss"="– Summer_North","South Summer Loss"="– Summer_South","Summer Loss"="– Summer","Winter Loss"="– Winter","North Winter Loss"="– Winter_North","South Winter Loss"="– Winter_South","Whole Dataset"="    Whole dataset","1400 Loss"="– 1400 m","1800 Loss"="– 1800 m","South Loss"="– South","North Loss"="– North"))+
                                                                                            theme(axis.line = element_line(colour = "black"),
                                                                                                  panel.grid.major = element_blank(),
                                                                                                  panel.grid.minor = element_blank(),
                                                                                                  panel.border = element_blank(),
                                                                                                  panel.background = element_blank())+theme(plot.title = element_text(size=16)))
p.loss.ratio
pdf("inextpd_ratio_reduction_20211221.pdf",width=17,height=6)
p.loss.ratio
dev.off()
