library(ggplot2) #ggplot2 v3.3.4
library(reshape2) #reshape2 v1.4.4
library(patchwork) #patchwork v1.0.0
setwd(HEREISYOURWORKINGDIRECTORY)
      
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

###################
#lwr_histogram
###################
#extract the likelihood weight ratios (LWRs) by gappa v0.5.0. The command line:
#`gappa examine lwr --jplace-path ./placement_10524_whole_13995_mitoalignment_V3_RAPPAS_Reduced_no_heur.jplace --out-dir Whole_10524_EPA_no_heur_lwr`
#import the lwr
EPA_no_heur.whole.lwr_histogram <- read.table("Data/lwr_histogram.csv",header=TRUE,sep = ",",stringsAsFactors = FALSE)
EPA_no_heur.whole.lwr_histogram = melt(EPA_no_heur.whole.lwr_histogram,id.vars="Range",measure.vars=c("Value.1","Value.2","Value.3"),variable.name="LWRs",value.name = "percentage")
EPA_no_heur.whole.lwr<-ggplot(EPA_no_heur.whole.lwr_histogram,aes(x=Range,y=percentage,fill=LWRs))+  geom_bar(stat="identity",position="dodge",width=0.8,colour="black")+ scale_fill_brewer(palette="Reds",direction=-1,labels=c("LWR 1","LWR 2","LWR 3")) +theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))+theme+theme(legend.title=element_text(colour="black"))+xlab("Likelihood Weight Ratios")+ylab("Number of Placements")+ ggtitle("Histogram of the first three most probable Likelihood Weight Ratios:\n Whole_10524_EPA_no_heur.")+theme(plot.title = element_text(hjust = 0.5, vjust=1,size=10))

###################
#edpl_histogram
###################
#calculate the Expected Distance between Placement Locations (EDPL) by gappa v0.5.0. The command line:
#`gappa examine edpl --jplace-path ./placement_10524_whole_13995_mitoalignment_V3_RAPPAS_Reduced_no_heur.jplace --out-dir Whole_10524_EPA_no_heur_edpl`
#import the edpl
EPA_no_heur.whole.edpl_histogram <- read.table("Data//edpl_histogram.csv",header=TRUE,sep = ",",stringsAsFactors = FALSE)
EPA_no_heur.whole.edpl<-ggplot(EPA_no_heur.whole.edpl_histogram,aes(x=Range,y=Value))+  geom_bar(stat="identity",position="dodge",width=0.8,colour="black")+ scale_fill_brewer(palette="Reds",direction=-1,labels=c("LWR 1","LWR 2","LWR 3")) +theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))+theme+theme(legend.title=element_text(colour="black"))+xlab("Expected Distance between Placement Locations (EDPL) ")+ylab("Number of Placements")+ ggtitle("Histogram of The Expected Distance between Placement Locations (EDPL) :\n Whole_10524_EPA_no_heur.")+theme(plot.title = element_text(hjust = 0.5, vjust=1,size=10))


pdf("FigS9_lwr_and_edpl_histogram_20211214.pdf",width=14,height=12)
EPA_no_heur.whole.lwr / EPA_no_heur.whole.edpl
dev.off()

