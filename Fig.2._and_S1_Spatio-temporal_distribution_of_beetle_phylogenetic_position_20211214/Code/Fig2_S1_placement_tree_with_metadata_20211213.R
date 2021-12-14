############################
#Spatio-temporal distribution of beetle on the placement tree 
############################
library(ggtree) #ggtree v2.0.4
library("treeio") #treeio v1.10.0
library(ggplot2) #ggplot2 v3.3.4
library(ape) #ape v5.4
library(patchwork) #patchwork v1.0.0
library(phylotate) #phylotate v1.3
library(treeman) #treeman v1.1.4
library(ips) #ips v0.0.11
library(dplyr) #dplyr v1.0.2
library(anytime) #anytime v0.3.9
library(ggnewscale) #ggnewscale v0.4.5

options(encoding="UTF-8")
setwd(HEREISYOURWORKINGDIRECTORY)
whole.newick.tree<-read.tree("Data/queries_10524_whole_13995_mitoalignment_V3_RAPPAS_Reduced_no_heur.newick")
whole.labeltree <- whole.newick.tree

whole.labeltree_p1<-ggtree(whole.newick.tree, layout="rectangular",branch.length='none')

whole.labeltree.pc1<-whole.labeltree_p1 + xlim(1, 20000) + theme(legend.position="none")

whole.labeltree.comm <- read.table("Data/GLGS_Community_Data_fixed_three_cases_20201019.txt",header=TRUE,sep = "\t",stringsAsFactors = FALSE)
whole.labeltree.m<-as.data.frame(whole.labeltree.comm[c('OTU','Collection_Date_formatted_1')])
#need to ignore the blank cell for the sample of failure finally.
whole.labeltree.data_m <- filter(whole.labeltree.m,whole.labeltree.m$OTU!='')

whole.labeltree.site <- as.data.frame(whole.labeltree.comm[c('OTU','Site_ID')])
whole.labeltree.site <- filter(whole.labeltree.site,whole.labeltree.site$OTU!='')

whole.labeltree.month <- whole.labeltree.data_m
whole.labeltree.data_month<-cbind(whole.labeltree.month,format(as.Date(anydate(whole.labeltree.data_m$Collection_Date_formatted_1)), "%Y_%m"))

whole.labeltree.metadata_Month_merge<-dplyr::inner_join(data.frame(whole.newick.tree$tip.label),whole.labeltree.data_month,by=c("whole.newick.tree.tip.label"="OTU"))

whole.labeltree.metadata_Month_Site_merge<-dplyr::inner_join(whole.labeltree.metadata_Month_merge,whole.labeltree.site,by=c("whole.newick.tree.tip.label"="OTU"))

whole.labeltree.metadata_Month_Site_merge <- cbind(as.data.frame(whole.labeltree.metadata_Month_Site_merge), Season =gsub("^[0-9]+$","Summer",gsub("^11$|^12$|^01$|^02$|^03$","Winter",gsub("[0-9]+_","",whole.labeltree.metadata_Month_Site_merge$`format(as.Date(anydate(whole.labeltree.data_m$Collection_Date_formatted_1)), `,perl=TRUE),perl=TRUE),perl=TRUE))

whole.labeltree.metadata_Month_Site_merge <- cbind(as.data.frame(whole.labeltree.metadata_Month_Site_merge), Location = gsub("^2015$","North",gsub("^2016$","North",gsub("_[0-9]+","",gsub("^2014_[0-9]+$|^2015_01$|^2015_02$|^2015_03$|^2015_04$|^2015_05$","South",whole.labeltree.metadata_Month_Site_merge$`format(as.Date(anydate(whole.labeltree.data_m$Collection_Date_formatted_1)), `,perl=TRUE),perl=TRUE),perl=TRUE),perl=TRUE))

whole.labeltree.metadata_Month_Site_merge <- cbind(as.data.frame(whole.labeltree.metadata_Month_Site_merge), Altitude = gsub("[1-4]$","",gsub("M","", gsub(".+_","", whole.labeltree.metadata_Month_Site_merge$Site_ID,perl=TRUE)),perl=TRUE))

whole.labeltree.tips.label<-data.frame(tips=whole.labeltree$tip.label)

whole.labeltree.label_df_merge<-dplyr::left_join(whole.labeltree.tips.label,whole.labeltree.metadata_Month_Site_merge,by=c("tips"="whole.newick.tree.tip.label"))

whole.labeltree.df <- data.frame(whole.labeltree.label_df_merge[,-1:-4])
rownames(whole.labeltree.df)<- whole.labeltree.label_df_merge$tips

whole.labeltree.df_test <- data.frame(whole.labeltree.df$Season,whole.labeltree.df$Season,whole.labeltree.df$Location,whole.labeltree.df$Location,whole.labeltree.df$Altitude,whole.labeltree.df$Altitude)
rownames(whole.labeltree.df_test)<-rownames(whole.labeltree.df)
whole.labeltree.df_test$whole.labeltree.df.Location<-gsub("North","Background",whole.labeltree.df_test$whole.labeltree.df.Location)
whole.labeltree.df_test$whole.labeltree.df.Location.1<-gsub("South","Background",whole.labeltree.df_test$whole.labeltree.df.Location.1)
whole.labeltree.df_test$whole.labeltree.df.Altitude<-gsub("1800","Background",whole.labeltree.df_test$whole.labeltree.df.Altitude)
whole.labeltree.df_test$whole.labeltree.df.Altitude.1<-gsub("1400","Background",whole.labeltree.df_test$whole.labeltree.df.Altitude.1)
whole.labeltree.df_test[is.na(whole.labeltree.df_test)]<-0

colnames(whole.labeltree.df_test)<- c("South","North","1400","1800")

spring.cols <- c("1400" = "Steelblue", "1800" = "Steelblue", "0"="white", "Background"="white")
summer.cols <- c("1400" = "#1d953f", "1800" = "#1d953f", "0"="white", "Background"="white")
autumn.cols <- c("1400" = "#896a45","1800" = "#896a45", "0"="white", "Background"="white")
winter.cols <- c("1400" = "black","1800" = "black", "0"="white", "Background"="white")

month.list<-unlist(rep(1:12,by=1))

month.list<-gsub("(^[1-9]$)","0\\1",month.list,perl=TRUE)

for (i in month.list){
  
  #Location discretion
  south.color<-data.frame(whole.labeltree.df_test[,3])
  rownames(south.color)<-rownames(whole.labeltree.df_test)
  colnames(south.color)<-"South"
  
  north.color<-data.frame(whole.labeltree.df_test[,4])
  rownames(north.color)<-rownames(whole.labeltree.df_test)
  colnames(north.color)<-"North"
  
  #rectangular plot for Altitude
  m1400.color<-data.frame(whole.labeltree.df_test[,5])
  rownames(m1400.color)<-rownames(whole.labeltree.df_test)
  colnames(m1400.color)<-"1400"
  m1400.color_filter<-filter(m1400.color,grepl("^1400$",m1400.color$`1400`,perl=TRUE))
  m1400.color_filter<-cbind(m1400.color_filter,name=row.names(m1400.color_filter))
  
  whole.labeltree.label_df_merge_03<-filter(whole.labeltree.label_df_merge, grepl(paste0("_",i,"$"),whole.labeltree.label_df_merge$`format(as.Date(anydate(whole.labeltree.data_m$Collection_Date_formatted_1)), `,perl=TRUE))
  whole.labeltree.label_df_merge_03<-whole.labeltree.label_df_merge_03[,1:3]
  
  south.color_filter<-filter(south.color,grepl("^South$",south.color$South,perl=TRUE))
  south.color_filter<-cbind(south.color_filter,name=row.names(south.color_filter))
  
  south<-dplyr::inner_join(m1400.color_filter,south.color_filter)
  south<-dplyr::inner_join(south,whole.labeltree.label_df_merge_03,by=c("name"="tips"))
  
  
  north.color_filter<-filter(north.color,grepl("^North$",north.color$North,perl=TRUE))
  north.color_filter<-cbind(north.color_filter,name=row.names(north.color_filter))
  
  north<-dplyr::inner_join(m1400.color_filter,north.color_filter)
  north<-dplyr::inner_join(north,whole.labeltree.label_df_merge_03,by=c("name"="tips"))
  
  
  south.merge<-dplyr::left_join(whole.labeltree.label_df_merge,south,by=c("tips"="name"))
  south.merge.inner<-dplyr::inner_join(whole.labeltree.label_df_merge,south,by=c("tips"="name"))
  south.merge_1 <- filter(south.merge,south.merge$tips %in% south.merge[is.na(south.merge$Altitude),]$tips == F)
  south.merge_1 <- filter(south.merge_1,south.merge_1$tips %in% south$name == F)
  south.merge_1$`1400`<-0
  
  south.merge_0 <- filter(south.merge, south.merge$tips %in% south.merge[is.na(south.merge$Altitude),]$tips == T)
  #south.merge_0$`1400`<-"background"
  
  south.merge_3<-rbind(south.merge_1,south.merge_0)
  south.merge.1400<-rbind(south.merge_3,south.merge.inner)
  
  south.merge.1400<-south.merge.1400[,-2:-7]
  south.merge.1400<-south.merge.1400[,-3:-5]
  
  #North discretion
  north.merge<-dplyr::left_join(whole.labeltree.label_df_merge,north,by=c("tips"="name"))
  north.merge.inner<-dplyr::inner_join(whole.labeltree.label_df_merge,north,by=c("tips"="name"))
  north.merge_1 <- filter(north.merge,north.merge$tips %in% north.merge[is.na(north.merge$Altitude),]$tips == F)
  north.merge_1 <- filter(north.merge_1,north.merge_1$tips %in% north$name == F)
  north.merge_1$`1400`<-0
  
  north.merge_0 <- filter(north.merge, north.merge$tips %in% north.merge[is.na(north.merge$Altitude),]$tips == T)
  #north.merge_0$`1400`<-"background"
  
  north.merge_3<-rbind(north.merge_1,north.merge_0)
  north.merge.1400<-rbind(north.merge_3,north.merge.inner)
  
  north.merge.1400<-north.merge.1400[,-2:-7]
  north.merge.1400<-north.merge.1400[,-3:-5]
  
  assign(paste0("month.",i),inner_join(south.merge.1400,north.merge.1400,by = c("tips"="tips")))
  
  assign(paste0("north.merge.1400.SN",i),data.frame("name"=eval(parse(text=paste0("month.",i)))[,1],South=eval(parse(text=paste0("month.",i)))[,2],North=eval(parse(text=paste0("month.",i)))[,3]))
  
}

rownames(north.merge.1400.SN01)<-north.merge.1400.SN01[,1]
north.merge.1400.SN01<-north.merge.1400.SN01[,-1]
colnames(north.merge.1400.SN01)<-c("South_01","North_01")

rownames(north.merge.1400.SN02)<-north.merge.1400.SN02[,1]
north.merge.1400.SN02<-north.merge.1400.SN02[,-1]
colnames(north.merge.1400.SN02)<-c("South_02","North_02")

rownames(north.merge.1400.SN03)<-north.merge.1400.SN03[,1]
north.merge.1400.SN03<-north.merge.1400.SN03[,-1]
colnames(north.merge.1400.SN03)<-c("South_03","North_03")

rownames(north.merge.1400.SN04)<-north.merge.1400.SN04[,1]
north.merge.1400.SN04<-north.merge.1400.SN04[,-1]
colnames(north.merge.1400.SN04)<-c("South_04","North_04")

rownames(north.merge.1400.SN05)<-north.merge.1400.SN05[,1]
north.merge.1400.SN05<-north.merge.1400.SN05[,-1]
colnames(north.merge.1400.SN05)<-c("South_05","North_05")

rownames(north.merge.1400.SN06)<-north.merge.1400.SN06[,1]
north.merge.1400.SN06<-north.merge.1400.SN06[,-1]
colnames(north.merge.1400.SN06)<-c("South_06","North_06")

rownames(north.merge.1400.SN07)<-north.merge.1400.SN07[,1]
north.merge.1400.SN07<-north.merge.1400.SN07[,-1]
colnames(north.merge.1400.SN07)<-c("South_07","North_07")

rownames(north.merge.1400.SN08)<-north.merge.1400.SN08[,1]
north.merge.1400.SN08<-north.merge.1400.SN08[,-1]
colnames(north.merge.1400.SN08)<-c("South_08","North_08")

rownames(north.merge.1400.SN09)<-north.merge.1400.SN09[,1]
north.merge.1400.SN09<-north.merge.1400.SN09[,-1]
colnames(north.merge.1400.SN09)<-c("South_09","North_09")

rownames(north.merge.1400.SN10)<-north.merge.1400.SN10[,1]
north.merge.1400.SN10<-north.merge.1400.SN10[,-1]
colnames(north.merge.1400.SN10)<-c("South_10","North_10")

rownames(north.merge.1400.SN11)<-north.merge.1400.SN11[,1]
north.merge.1400.SN11<-north.merge.1400.SN11[,-1]
colnames(north.merge.1400.SN11)<-c("South_11","North_11")

rownames(north.merge.1400.SN12)<-north.merge.1400.SN12[,1]
north.merge.1400.SN12<-north.merge.1400.SN12[,-1]
colnames(north.merge.1400.SN12)<-c("South_12","North_12")

whole.labeltree.1400.p1<-gheatmap(whole.labeltree.pc1, north.merge.1400.SN03, offset=-300, width=3,
                                  colnames = TRUE,color = FALSE, font.size=2, 
                                  colnames_angle=-45, hjust=0)  + scale_fill_manual( values = spring.cols,
                                                                                     breaks = c("1400", "0","Background"),
                                                                                     labels = c("1400", "0","Background"),
                                                                                     name="South") + theme(legend.position="bottom") 
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1400.p2 <- whole.labeltree.1400.p1 + new_scale_fill()

whole.labeltree.1400.p3<-gheatmap(whole.labeltree.1400.p2,  north.merge.1400.SN04, offset=1244.4, width=3,
                                  colnames = TRUE,color = FALSE, font.size=2, 
                                  colnames_angle=-45, hjust=0)     + scale_fill_manual( values = spring.cols,
                                                                                        breaks = c("1400","0","Background"),
                                                                                        labels = c("1400","0","Background"),
                                                                                        name="North") + theme(legend.position="bottom") 
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1400.p4 <- whole.labeltree.1400.p3 + new_scale_fill()

whole.labeltree.1400.p5<-gheatmap(whole.labeltree.1400.p4,  north.merge.1400.SN05, offset=2788.8, width=3,
                                  colnames = TRUE,color = FALSE, font.size=2, 
                                  colnames_angle=-45, hjust=0)     + scale_fill_manual( values = spring.cols,
                                                                                        breaks = c("1400","0","Background"),
                                                                                        labels = c("1400","0","Background"),
                                                                                        name="North") + theme(legend.position="bottom") 


## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1400.p6 <- whole.labeltree.1400.p5 + new_scale_fill()

whole.labeltree.1400.p7<-gheatmap(whole.labeltree.1400.p6,  north.merge.1400.SN06, offset=4333.2, width=3,
                                  colnames = TRUE,color = FALSE, font.size=2, 
                                  colnames_angle=-45, hjust=0)     + scale_fill_manual( values = summer.cols,
                                                                                        breaks = c("1400","0","Background"),
                                                                                        labels = c("1400","0","Background"),
                                                                                        name="North") + theme(legend.position="bottom")

## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1400.p8 <- whole.labeltree.1400.p7 + new_scale_fill()

whole.labeltree.1400.p9<-gheatmap(whole.labeltree.1400.p8,  north.merge.1400.SN07, offset=5877.6, width=3,
                                  colnames = TRUE,color = FALSE, font.size=2, 
                                  colnames_angle=-45, hjust=0)     + scale_fill_manual( values = summer.cols,
                                                                                        breaks = c("1400","0","Background"),
                                                                                        labels = c("1400","0","Background"),
                                                                                        name="North") + theme(legend.position="bottom")

## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1400.p10 <- whole.labeltree.1400.p9 + new_scale_fill()

whole.labeltree.1400.p11<-gheatmap(whole.labeltree.1400.p10,  north.merge.1400.SN08, offset=7422, width=3,
                                   colnames = TRUE,color = FALSE, font.size=2, 
                                   colnames_angle=-45, hjust=0)     + scale_fill_manual( values = summer.cols,
                                                                                         breaks = c("1400","0","Background"),
                                                                                         labels = c("1400","0","Background"),
                                                                                         name="North") + theme(legend.position="bottom")
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1400.p12 <- whole.labeltree.1400.p11 + new_scale_fill()

whole.labeltree.1400.p13<-gheatmap(whole.labeltree.1400.p12,  north.merge.1400.SN09, offset=8966.4, width=3,
                                   colnames = TRUE,color = FALSE, font.size=2, 
                                   colnames_angle=-45, hjust=0)     + scale_fill_manual( values = autumn.cols,
                                                                                         breaks = c("1400","0","Background"),
                                                                                         labels = c("1400","0","Background"),
                                                                                         name="North") + theme(legend.position="bottom")
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1400.p14 <- whole.labeltree.1400.p13 + new_scale_fill()

whole.labeltree.1400.p15<-gheatmap(whole.labeltree.1400.p14,  north.merge.1400.SN10, offset=10510.8, width=3,
                                   colnames = TRUE,color = FALSE, font.size=2, 
                                   colnames_angle=-45, hjust=0)     + scale_fill_manual( values = autumn.cols,
                                                                                         breaks = c("1400","0","Background"),
                                                                                         labels = c("1400","0","Background"),
                                                                                         name="North") + theme(legend.position="bottom")
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1400.p16 <- whole.labeltree.1400.p15 + new_scale_fill()

whole.labeltree.1400.p17<-gheatmap(whole.labeltree.1400.p16,  north.merge.1400.SN11, offset=12055.2, width=3,
                                   colnames = TRUE,color = FALSE, font.size=2, 
                                   colnames_angle=-45, hjust=0)     + scale_fill_manual( values = autumn.cols,
                                                                                         breaks = c("1400","0","Background"),
                                                                                         labels = c("1400","0","Background"),
                                                                                         name="North") + theme(legend.position="bottom")
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1400.p18 <- whole.labeltree.1400.p17 + new_scale_fill()

whole.labeltree.1400.p19<-gheatmap(whole.labeltree.1400.p18,  north.merge.1400.SN12, offset=13599.6, width=3,
                                   colnames = TRUE,color = FALSE, font.size=2, 
                                   colnames_angle=-45, hjust=0)     + scale_fill_manual( values = winter.cols,
                                                                                         breaks = c("1400","0","Background"),
                                                                                         labels = c("1400","0","Background"),
                                                                                         name="North") + theme(legend.position="bottom")
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1400.p20 <- whole.labeltree.1400.p19 + new_scale_fill()

whole.labeltree.1400.p21<-gheatmap(whole.labeltree.1400.p20,  north.merge.1400.SN01, offset=15144, width=3,
                                   colnames = TRUE,color = FALSE, font.size=2, 
                                   colnames_angle=-45, hjust=0)     + scale_fill_manual( values = winter.cols,
                                                                                         breaks = c("1400","0","Background"),
                                                                                         labels = c("1400","0","Background"),
                                                                                         name="North") + theme(legend.position="bottom")
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.


whole.labeltree.1400.p22 <- whole.labeltree.1400.p21 + new_scale_fill()

whole.labeltree.1400.p23<-gheatmap(whole.labeltree.1400.p22,  north.merge.1400.SN02, offset=16688.4, width=3,
                                   colnames = TRUE,color = FALSE, font.size=2, 
                                   colnames_angle=-45, hjust=0)     + scale_fill_manual( values = winter.cols,
                                                                                         breaks = c("1400","0","Background"),
                                                                                         labels = c("1400","0","Background"),
                                                                                         name="North") + theme(legend.position="bottom")
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
whole.labeltree.1400.p23
# dev.new()
# svg(file="./fig2_1400_20211214.svg",width=80, height=80)
# whole.labeltree.1400.p23
# dev.off()

month.list<-unlist(rep(1:12,by=1))

month.list<-gsub("(^[1-9]$)","0\\1",month.list,perl=TRUE)

for (i in month.list){
  
  #Location discretion
  south.color<-data.frame(whole.labeltree.df_test[,3])
  rownames(south.color)<-rownames(whole.labeltree.df_test)
  colnames(south.color)<-"South"
  
  north.color<-data.frame(whole.labeltree.df_test[,4])
  rownames(north.color)<-rownames(whole.labeltree.df_test)
  colnames(north.color)<-"North"
  
  #rectangular plot for Altitude
  m1800.color<-data.frame(whole.labeltree.df_test[,6])
  rownames(m1800.color)<-rownames(whole.labeltree.df_test)
  colnames(m1800.color)<-"1800"
  m1800.color_filter<-filter(m1800.color,grepl("^1800$",m1800.color$`1800`,perl=TRUE))
  m1800.color_filter<-cbind(m1800.color_filter,name=row.names(m1800.color_filter))
  
  whole.labeltree.label_df_merge_03<-filter(whole.labeltree.label_df_merge, grepl(paste0("_",i,"$"),whole.labeltree.label_df_merge$`format(as.Date(anydate(whole.labeltree.data_m$Collection_Date_formatted_1)), `,perl=TRUE))
  whole.labeltree.label_df_merge_03<-whole.labeltree.label_df_merge_03[,1:3]
  
  south.color_filter<-filter(south.color,grepl("^South$",south.color$South,perl=TRUE))
  south.color_filter<-cbind(south.color_filter,name=row.names(south.color_filter))
  
  south<-dplyr::inner_join(m1800.color_filter,south.color_filter)
  south<-dplyr::inner_join(south,whole.labeltree.label_df_merge_03,by=c("name"="tips"))
  
  
  north.color_filter<-filter(north.color,grepl("^North$",north.color$North,perl=TRUE))
  north.color_filter<-cbind(north.color_filter,name=row.names(north.color_filter))
  
  north<-dplyr::inner_join(m1800.color_filter,north.color_filter)
  north<-dplyr::inner_join(north,whole.labeltree.label_df_merge_03,by=c("name"="tips"))
  
  
  south.merge<-dplyr::left_join(whole.labeltree.label_df_merge,south,by=c("tips"="name"))
  south.merge.inner<-dplyr::inner_join(whole.labeltree.label_df_merge,south,by=c("tips"="name"))
  south.merge_1 <- filter(south.merge,south.merge$tips %in% south.merge[is.na(south.merge$Altitude),]$tips == F)
  south.merge_1 <- filter(south.merge_1,south.merge_1$tips %in% south$name == F)
  south.merge_1$`1800`<-0
  
  south.merge_0 <- filter(south.merge, south.merge$tips %in% south.merge[is.na(south.merge$Altitude),]$tips == T)
  #south.merge_0$`1800`<-"background"
  
  south.merge_3<-rbind(south.merge_1,south.merge_0)
  south.merge.1800<-rbind(south.merge_3,south.merge.inner)
  
  south.merge.1800<-south.merge.1800[,-2:-7]
  south.merge.1800<-south.merge.1800[,-3:-5]
  
  #North discretion
  north.merge<-dplyr::left_join(whole.labeltree.label_df_merge,north,by=c("tips"="name"))
  north.merge.inner<-dplyr::inner_join(whole.labeltree.label_df_merge,north,by=c("tips"="name"))
  north.merge_1 <- filter(north.merge,north.merge$tips %in% north.merge[is.na(north.merge$Altitude),]$tips == F)
  north.merge_1 <- filter(north.merge_1,north.merge_1$tips %in% north$name == F)
  north.merge_1$`1800`<-0
  
  north.merge_0 <- filter(north.merge, north.merge$tips %in% north.merge[is.na(north.merge$Altitude),]$tips == T)
  #north.merge_0$`1800`<-"background"
  
  north.merge_3<-rbind(north.merge_1,north.merge_0)
  north.merge.1800<-rbind(north.merge_3,north.merge.inner)
  
  north.merge.1800<-north.merge.1800[,-2:-7]
  north.merge.1800<-north.merge.1800[,-3:-5]
  
  assign(paste0("month.",i),inner_join(south.merge.1800,north.merge.1800,by = c("tips"="tips")))
  
  assign(paste0("north.merge.1800.SN",i),data.frame("name"=eval(parse(text=paste0("month.",i)))[,1],South=eval(parse(text=paste0("month.",i)))[,2],North=eval(parse(text=paste0("month.",i)))[,3]))
  
}


rownames(north.merge.1800.SN01)<-north.merge.1800.SN01[,1]
north.merge.1800.SN01<-north.merge.1800.SN01[,-1]
colnames(north.merge.1800.SN01)<-c("South_01","North_01")

rownames(north.merge.1800.SN02)<-north.merge.1800.SN02[,1]
north.merge.1800.SN02<-north.merge.1800.SN02[,-1]
colnames(north.merge.1800.SN02)<-c("South_02","North_02")

rownames(north.merge.1800.SN03)<-north.merge.1800.SN03[,1]
north.merge.1800.SN03<-north.merge.1800.SN03[,-1]
colnames(north.merge.1800.SN03)<-c("South_03","North_03")

rownames(north.merge.1800.SN04)<-north.merge.1800.SN04[,1]
north.merge.1800.SN04<-north.merge.1800.SN04[,-1]
colnames(north.merge.1800.SN04)<-c("South_04","North_04")

rownames(north.merge.1800.SN05)<-north.merge.1800.SN05[,1]
north.merge.1800.SN05<-north.merge.1800.SN05[,-1]
colnames(north.merge.1800.SN05)<-c("South_05","North_05")

rownames(north.merge.1800.SN06)<-north.merge.1800.SN06[,1]
north.merge.1800.SN06<-north.merge.1800.SN06[,-1]
colnames(north.merge.1800.SN06)<-c("South_06","North_06")

rownames(north.merge.1800.SN07)<-north.merge.1800.SN07[,1]
north.merge.1800.SN07<-north.merge.1800.SN07[,-1]
colnames(north.merge.1800.SN07)<-c("South_07","North_07")

rownames(north.merge.1800.SN08)<-north.merge.1800.SN08[,1]
north.merge.1800.SN08<-north.merge.1800.SN08[,-1]
colnames(north.merge.1800.SN08)<-c("South_08","North_08")

rownames(north.merge.1800.SN09)<-north.merge.1800.SN09[,1]
north.merge.1800.SN09<-north.merge.1800.SN09[,-1]
colnames(north.merge.1800.SN09)<-c("South_09","North_09")

rownames(north.merge.1800.SN10)<-north.merge.1800.SN10[,1]
north.merge.1800.SN10<-north.merge.1800.SN10[,-1]
colnames(north.merge.1800.SN10)<-c("South_10","North_10")

rownames(north.merge.1800.SN11)<-north.merge.1800.SN11[,1]
north.merge.1800.SN11<-north.merge.1800.SN11[,-1]
colnames(north.merge.1800.SN11)<-c("South_11","North_11")

rownames(north.merge.1800.SN12)<-north.merge.1800.SN12[,1]
north.merge.1800.SN12<-north.merge.1800.SN12[,-1]
colnames(north.merge.1800.SN12)<-c("South_12","North_12")

whole.labeltree.1800.p1<-gheatmap(whole.labeltree.pc1, north.merge.1800.SN03, offset=-300, width=3,
                                  colnames = TRUE,color = FALSE, font.size=2, 
                                  colnames_angle=-45, hjust=0)  + scale_fill_manual( values = spring.cols,
                                                                                     breaks = c("1800", "0","Background"),
                                                                                     labels = c("1800", "0","Background"),
                                                                                     name="South") + theme(legend.position="bottom") 
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1800.p2 <- whole.labeltree.1800.p1 + new_scale_fill()

whole.labeltree.1800.p3<-gheatmap(whole.labeltree.1800.p2,  north.merge.1800.SN04, offset=1244.4, width=3,
                                  colnames = TRUE,color = FALSE, font.size=2, 
                                  colnames_angle=-45, hjust=0)     + scale_fill_manual( values = spring.cols,
                                                                                        breaks = c("1800","0","Background"),
                                                                                        labels = c("1800","0","Background"),
                                                                                        name="North") + theme(legend.position="bottom") 
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1800.p4 <- whole.labeltree.1800.p3 + new_scale_fill()

whole.labeltree.1800.p5<-gheatmap(whole.labeltree.1800.p4,  north.merge.1800.SN05, offset=2788.8, width=3,
                                  colnames = TRUE,color = FALSE, font.size=2, 
                                  colnames_angle=-45, hjust=0)     + scale_fill_manual( values = spring.cols,
                                                                                        breaks = c("1800","0","Background"),
                                                                                        labels = c("1800","0","Background"),
                                                                                        name="North") + theme(legend.position="bottom") 


## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1800.p6 <- whole.labeltree.1800.p5 + new_scale_fill()

whole.labeltree.1800.p7<-gheatmap(whole.labeltree.1800.p6,  north.merge.1800.SN06, offset=4333.2, width=3,
                                  colnames = TRUE,color = FALSE, font.size=2, 
                                  colnames_angle=-45, hjust=0)     + scale_fill_manual( values = summer.cols,
                                                                                        breaks = c("1800","0","Background"),
                                                                                        labels = c("1800","0","Background"),
                                                                                        name="North") + theme(legend.position="bottom")

## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1800.p8 <- whole.labeltree.1800.p7 + new_scale_fill()

whole.labeltree.1800.p9<-gheatmap(whole.labeltree.1800.p8,  north.merge.1800.SN07, offset=5877.6, width=3,
                                  colnames = TRUE,color = FALSE, font.size=2, 
                                  colnames_angle=-45, hjust=0)     + scale_fill_manual( values = summer.cols,
                                                                                        breaks = c("1800","0","Background"),
                                                                                        labels = c("1800","0","Background"),
                                                                                        name="North") + theme(legend.position="bottom")

## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1800.p10 <- whole.labeltree.1800.p9 + new_scale_fill()

whole.labeltree.1800.p11<-gheatmap(whole.labeltree.1800.p10,  north.merge.1800.SN08, offset=7422, width=3,
                                   colnames = TRUE,color = FALSE, font.size=2, 
                                   colnames_angle=-45, hjust=0)     + scale_fill_manual( values = summer.cols,
                                                                                         breaks = c("1800","0","Background"),
                                                                                         labels = c("1800","0","Background"),
                                                                                         name="North") + theme(legend.position="bottom")
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1800.p12 <- whole.labeltree.1800.p11 + new_scale_fill()

whole.labeltree.1800.p13<-gheatmap(whole.labeltree.1800.p12,  north.merge.1800.SN09, offset=8966.4, width=3,
                                   colnames = TRUE,color = FALSE, font.size=2, 
                                   colnames_angle=-45, hjust=0)     + scale_fill_manual( values = autumn.cols,
                                                                                         breaks = c("1800","0","Background"),
                                                                                         labels = c("1800","0","Background"),
                                                                                         name="North") + theme(legend.position="bottom")
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1800.p14 <- whole.labeltree.1800.p13 + new_scale_fill()

whole.labeltree.1800.p15<-gheatmap(whole.labeltree.1800.p14,  north.merge.1800.SN10, offset=10510.8, width=3,
                                   colnames = TRUE,color = FALSE, font.size=2, 
                                   colnames_angle=-45, hjust=0)     + scale_fill_manual( values = autumn.cols,
                                                                                         breaks = c("1800","0","Background"),
                                                                                         labels = c("1800","0","Background"),
                                                                                         name="North") + theme(legend.position="bottom")
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1800.p16 <- whole.labeltree.1800.p15 + new_scale_fill()

whole.labeltree.1800.p17<-gheatmap(whole.labeltree.1800.p16,  north.merge.1800.SN11, offset=12055.2, width=3,
                                   colnames = TRUE,color = FALSE, font.size=2, 
                                   colnames_angle=-45, hjust=0)     + scale_fill_manual( values = autumn.cols,
                                                                                         breaks = c("1800","0","Background"),
                                                                                         labels = c("1800","0","Background"),
                                                                                         name="North") + theme(legend.position="bottom")
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1800.p18 <- whole.labeltree.1800.p17 + new_scale_fill()

whole.labeltree.1800.p19<-gheatmap(whole.labeltree.1800.p18,  north.merge.1800.SN12, offset=13599.6, width=3,
                                   colnames = TRUE,color = FALSE, font.size=2, 
                                   colnames_angle=-45, hjust=0)     + scale_fill_manual( values = winter.cols,
                                                                                         breaks = c("1800","0","Background"),
                                                                                         labels = c("1800","0","Background"),
                                                                                         name="North") + theme(legend.position="bottom")
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.

whole.labeltree.1800.p20 <- whole.labeltree.1800.p19 + new_scale_fill()

whole.labeltree.1800.p21<-gheatmap(whole.labeltree.1800.p20,  north.merge.1800.SN01, offset=15144, width=3,
                                   colnames = TRUE,color = FALSE, font.size=2, 
                                   colnames_angle=-45, hjust=0)     + scale_fill_manual( values = winter.cols,
                                                                                         breaks = c("1800","0","Background"),
                                                                                         labels = c("1800","0","Background"),
                                                                                         name="North") + theme(legend.position="bottom")
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.


whole.labeltree.1800.p22 <- whole.labeltree.1800.p21 + new_scale_fill()

whole.labeltree.1800.p23<-gheatmap(whole.labeltree.1800.p22,  north.merge.1800.SN02, offset=16688.4, width=3,
                                   colnames = TRUE,color = FALSE, font.size=2, 
                                   colnames_angle=-45, hjust=0)     + scale_fill_manual( values = winter.cols,
                                                                                         breaks = c("1800","0","Background"),
                                                                                         labels = c("1800","0","Background"),
                                                                                         name="North") + theme(legend.position="bottom")

whole.labeltree.1800.p23
# dev.new()
# svg(file="./figS1_1800_20211214.svg",width=80, height=80)
# whole.labeltree.1800.p23
# dev.off()

############################
#add climate variable
############################
library(vegan) #vegan v2.5-6
library(ggplot2) #ggplot2 v3.3.4
library(patchwork) #patchwork v1.0.0
library(ggordiplots) #ggordiplots v0.3.0
library(anytime) #anytime v0.3.9
library(dplyr) #dplyr v1.0.2
library(reshape) #reshape v0.8.8
library(ggfortify) #ggfortify v0.4.11
library(cowplot) #cowplot v1.1.1
library(stringr) #stringr v1.4.0

setwd(HEREISYOURWORKINGDIRECTORY)

#daily climate variables
#month precipitation
tem.precipitation <- read.table("Data/Climate_daily_precipitation_20210402.csv",header=TRUE,sep = ",",stringsAsFactors = FALSE)
tem.precipitation$Number<-paste0(tem.precipitation$x,"_",tem.precipitation$y,"_",tem.precipitation$date)
tem.tempertrue <- read.table("Data/Climate_daily_tempertrue_mean_20210402.csv",header=TRUE,sep = ",",stringsAsFactors = FALSE)
tem.tempertrue$Number<-paste0(tem.tempertrue$x,"_",tem.tempertrue$y,"_",tem.tempertrue$date)
summary.tempertrue.precipitation<-dplyr::inner_join(tem.precipitation,tem.tempertrue,by=c("Number"="Number"))

summary.tempertrue.precipitation$Date<-paste0(lubridate::ymd(as.Date(anydate(summary.tempertrue.precipitation$date.x))))
dlj_list<-c("DLJ_14001","DLJ_14001M","DLJ_14002","DLJ_14003","DLJ_14004","DLJ_18001","DLJ_18002")
tem.precipitation.dlj<-summary.tempertrue.precipitation[(summary.tempertrue.precipitation$Site_ID.x %in% dlj_list),]
`%!in%` = Negate(`%in%`)
tem.precipitation.bhl<-summary.tempertrue.precipitation[(summary.tempertrue.precipitation$Site_ID.x %!in% dlj_list),]
tem.precipitation.bhl<-filter(tem.precipitation.bhl,tem.precipitation.bhl$precipitation>0)
tem.precipitation.dlj<-filter(tem.precipitation.dlj,tem.precipitation.dlj$precipitation>0)

#for BHL 1400
tem.precipitation.bhl.1400<-filter(tem.precipitation.bhl,tem.precipitation.bhl$glgs_DEM_30m_area_4326.x<1800)
tem.precipitation.bhl.1400.year<-filter(tem.precipitation.bhl.1400,tem.precipitation.bhl.1400$Date>as.Date("2014-06-06"))
tem.precipitation.bhl.1400.year<-filter(tem.precipitation.bhl.1400,tem.precipitation.bhl.1400$Date<as.Date("2015-05-29"))

tem.precipitation.bhl.1400.year$month <- factor(lubridate::month(tem.precipitation.bhl.1400.year$Date), levels=c("3", "4", "5","6","7","8","9","10","11","12","1","2"), ordered=TRUE)
tem.precipitation.bhl.1400.year.order<-tem.precipitation.bhl.1400.year[order(tem.precipitation.bhl.1400.year$month, tem.precipitation.bhl.1400.year$Date),]
tem.precipitation.bhl.1400.year.order$order<-c("NULL")
tem.precipitation.bhl.1400.year.order[which(lubridate::month(tem.precipitation.bhl.1400.year.order$Date)>=3),]$order <- paste0("2014/",lubridate::month(tem.precipitation.bhl.1400.year.order[which(lubridate::month(tem.precipitation.bhl.1400.year.order$Date)>=3),]$Date),"/",lubridate::day(tem.precipitation.bhl.1400.year.order[which(lubridate::month(tem.precipitation.bhl.1400.year.order$Date)>=3),]$Date))
tem.precipitation.bhl.1400.year.order[which(lubridate::month(tem.precipitation.bhl.1400.year.order$Date)<3),]$order <- paste0("2015/",lubridate::month(tem.precipitation.bhl.1400.year.order[which(lubridate::month(tem.precipitation.bhl.1400.year.order$Date)<3),]$Date),"/",lubridate::day(tem.precipitation.bhl.1400.year.order[which(lubridate::month(tem.precipitation.bhl.1400.year.order$Date)<3),]$Date))


# ggplot(tem.precipitation.bhl.1400.year.order, mapping = aes(x = factor(as.Date(order)), y = precipitation, group = 1)) + geom_point(col = "grey")+geom_smooth(col = "blue",linetype=1,span = 0.5, method = "loess", size = 0.5) + xlab('Year')+
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())+theme(axis.text.x = element_blank())+theme(legend.position = "none")+theme(legend.title=element_text(colour="black"))

#for DLJ 1400
tem.precipitation.dlj.1400<-filter(tem.precipitation.dlj,tem.precipitation.dlj$glgs_DEM_30m_area_4326.x<1800)
tem.precipitation.dlj.1400.year<-filter(tem.precipitation.dlj.1400,tem.precipitation.dlj.1400$Date>as.Date("2015-08-18"))
tem.precipitation.dlj.1400.year<-filter(tem.precipitation.dlj.1400.year,tem.precipitation.dlj.1400.year$Date<as.Date("2016-08-02"))

tem.precipitation.dlj.1400.year$month <- factor(lubridate::month(tem.precipitation.dlj.1400.year$Date), levels=c("3", "4", "5","6","7","8","9","10","11","12","1","2"), ordered=TRUE)
tem.precipitation.dlj.1400.year.order<-tem.precipitation.dlj.1400.year[order(tem.precipitation.dlj.1400.year$month, tem.precipitation.dlj.1400.year$Date),]
tem.precipitation.dlj.1400.year.order$order<-c("NULL")
tem.precipitation.dlj.1400.year.order[which(lubridate::month(tem.precipitation.dlj.1400.year.order$Date)>=3),]$order <- paste0("2014/",lubridate::month(tem.precipitation.dlj.1400.year.order[which(lubridate::month(tem.precipitation.dlj.1400.year.order$Date)>=3),]$Date),"/",lubridate::day(tem.precipitation.dlj.1400.year.order[which(lubridate::month(tem.precipitation.dlj.1400.year.order$Date)>=3),]$Date))
tem.precipitation.dlj.1400.year.order[which(lubridate::month(tem.precipitation.dlj.1400.year.order$Date)<3),]$order <- paste0("2015/",lubridate::month(tem.precipitation.dlj.1400.year.order[which(lubridate::month(tem.precipitation.dlj.1400.year.order$Date)<3),]$Date),"/",lubridate::day(tem.precipitation.dlj.1400.year.order[which(lubridate::month(tem.precipitation.dlj.1400.year.order$Date)<3),]$Date))

# ggplot(tem.precipitation.dlj.1400.year.order, mapping = aes(x = factor(as.Date(order)), y = precipitation, group = 1)) + geom_point(col = "grey")+geom_smooth(col = "blue",linetype=1,span = 0.5, method = "loess", size = 0.5) + xlab('Year')+
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())+theme(axis.text.x = element_blank())+theme(legend.position = "left")+theme(legend.title=element_text(colour="black"))

tem.precipitation.dlj.1400.year.order$group<-c("DLJ")
tem.precipitation.bhl.1400.year.order$group<-c("BHL")
tem.precipitation<-rbind(tem.precipitation.bhl.1400.year.order,tem.precipitation.dlj.1400.year.order)

climate.plot<-ggplot(tem.precipitation, mapping = aes(x = factor(as.Date(order)), y = precipitation, group = group)) + geom_point(tem.precipitation[which(tem.precipitation$group=="BHL"),],mapping = aes(x = factor(as.Date(order)), y = precipitation),col = "#D55E00",shape = 0,alpha = 1/20)+
  geom_point(tem.precipitation[which(tem.precipitation$group=="DLJ"),],mapping = aes(x = factor(as.Date(order)), y = precipitation),col = "#009E73",shape = 6,alpha = 1/20)+
  geom_smooth(tem.precipitation[which(tem.precipitation$group=="BHL"),],mapping = aes(x = factor(as.Date(order)), y = precipitation),col = "#D55E00",linetype=3, method = "loess", size = 0.5, fill="#D55E00") +
  geom_smooth(tem.precipitation[which(tem.precipitation$group=="DLJ"),],mapping = aes(x = factor(as.Date(order)), y = precipitation),col = "#009E73",linetype=3, method = "loess", size = 0.5,fill = "#009E73") + xlab('Day')+
  geom_point(tem.precipitation[which(tem.precipitation$group=="BHL"),],mapping = aes(x = factor(as.Date(order)), y = temperture/0.6),col = "#D55E00",shape = 4,alpha = 1/20)+
  geom_point(tem.precipitation[which(tem.precipitation$group=="DLJ"),],mapping = aes(x = factor(as.Date(order)), y = temperture/0.6),col = "#009E73",shape = 3,alpha = 1/20)+
  geom_smooth(tem.precipitation[which(tem.precipitation$group=="BHL"),],mapping = aes(x = factor(as.Date(order)), y = temperture/0.6),col = "#D55E00",linetype=1, method = "loess", size = 0.5, fill="#D55E00") +
  geom_smooth(tem.precipitation[which(tem.precipitation$group=="DLJ"),],mapping = aes(x = factor(as.Date(order)), y = temperture/0.6),col = "#009E73",linetype=1, method = "loess", size = 0.5,fill = "#009E73") + xlab('Day')+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+theme(axis.text.x = element_blank())+theme(legend.position = "none")+theme(legend.title=element_text(colour="black"))+ylab("Daily average precipitation")+scale_y_continuous(sec.axis = sec_axis(~.*0.6, name = "Daily average temperture"))

climate.plot.ref<-climate.plot + geom_vline(xintercept = factor(as.Date(c("2014-05-31","2014-08-31","2014-11-30"))))
climate.plot.ref

# ggplot(tem.precipitation, mapping = aes(x = factor(as.Date(order)), y = precipitation, group = group)) + geom_point(tem.precipitation[which(tem.precipitation$group=="BHL"),],mapping = aes(x = factor(as.Date(order)), y = precipitation),col = "#D55E00",shape = 0,alpha = 1/20)+
#   geom_point(tem.precipitation[which(tem.precipitation$group=="DLJ"),],mapping = aes(x = factor(as.Date(order)), y = precipitation),col = "#009E73",shape = 6,alpha = 1/20)+
#   geom_smooth(tem.precipitation[which(tem.precipitation$group=="BHL"),],mapping = aes(x = factor(as.Date(order)), y = precipitation),col = "#D55E00",linetype=3, method = "loess", size = 0.5, fill="#D55E00") +
#   geom_smooth(tem.precipitation[which(tem.precipitation$group=="DLJ"),],mapping = aes(x = factor(as.Date(order)), y = precipitation),col = "#009E73",linetype=3, method = "loess", size = 0.5,fill = "#009E73") + xlab('Day')+
#   geom_point(tem.precipitation[which(tem.precipitation$group=="BHL"),],mapping = aes(x = factor(as.Date(order)), y = temperture/0.6),col = "#D55E00",shape = 4,alpha = 1/20)+
#   geom_point(tem.precipitation[which(tem.precipitation$group=="DLJ"),],mapping = aes(x = factor(as.Date(order)), y = temperture/0.6),col = "#009E73",shape = 3,alpha = 1/20)+
#   geom_smooth(tem.precipitation[which(tem.precipitation$group=="BHL"),],mapping = aes(x = factor(as.Date(order)), y = temperture/0.6),col = "#D55E00",linetype=1, method = "loess", size = 0.5, fill="#D55E00") +
#   geom_smooth(tem.precipitation[which(tem.precipitation$group=="DLJ"),],mapping = aes(x = factor(as.Date(order)), y = temperture/0.6),col = "#009E73",linetype=1, method = "loess", size = 0.5,fill = "#009E73") + xlab('Day')+theme(legend.position = "left")+theme(legend.title=element_text(colour="black"))

dev.new()
ggdraw() + draw_plot(whole.labeltree.1400.p23, x = 0, y = 0, width = 1, height = 1) +
  draw_plot(climate.plot.ref, x = 0.057, y = 0.06, width = 0.873, height = .09) 
ggsave("fig2_climate_20211213.pdf", width = 80, height = 80, units = "cm")
dev.off()

# dev.new()
# svg(file="./fig2_20211213.svg",width=80, height=80)
# ggdraw() + draw_plot(whole.labeltree.1400.p23, x = 0, y = 0, width = 1, height = 1) +
#   draw_plot(climate.plot.ref, x = 0.057, y = 0.06, width = 0.873, height = .09)
# dev.off()

#for 1800m
#for BHL 1800
tem.precipitation.bhl.1800<-filter(tem.precipitation.bhl,tem.precipitation.bhl$glgs_DEM_30m_area_4326.x>1800)
tem.precipitation.bhl.1800.year<-filter(tem.precipitation.bhl.1800,tem.precipitation.bhl.1800$Date>as.Date("2014-06-06"))
tem.precipitation.bhl.1800.year<-filter(tem.precipitation.bhl.1800,tem.precipitation.bhl.1800$Date<as.Date("2015-05-29"))

tem.precipitation.bhl.1800.year$month <- factor(lubridate::month(tem.precipitation.bhl.1800.year$Date), levels=c("3", "4", "5","6","7","8","9","10","11","12","1","2"), ordered=TRUE)
tem.precipitation.bhl.1800.year.order<-tem.precipitation.bhl.1800.year[order(tem.precipitation.bhl.1800.year$month, tem.precipitation.bhl.1800.year$Date),]
tem.precipitation.bhl.1800.year.order$order<-c("NULL")
tem.precipitation.bhl.1800.year.order[which(lubridate::month(tem.precipitation.bhl.1800.year.order$Date)>=3),]$order <- paste0("2014/",lubridate::month(tem.precipitation.bhl.1800.year.order[which(lubridate::month(tem.precipitation.bhl.1800.year.order$Date)>=3),]$Date),"/",lubridate::day(tem.precipitation.bhl.1800.year.order[which(lubridate::month(tem.precipitation.bhl.1800.year.order$Date)>=3),]$Date))
tem.precipitation.bhl.1800.year.order[which(lubridate::month(tem.precipitation.bhl.1800.year.order$Date)<3),]$order <- paste0("2015/",lubridate::month(tem.precipitation.bhl.1800.year.order[which(lubridate::month(tem.precipitation.bhl.1800.year.order$Date)<3),]$Date),"/",lubridate::day(tem.precipitation.bhl.1800.year.order[which(lubridate::month(tem.precipitation.bhl.1800.year.order$Date)<3),]$Date))


# ggplot(tem.precipitation.bhl.1800.year.order, mapping = aes(x = factor(as.Date(order)), y = precipitation, group = 1)) + geom_point(col = "grey")+geom_smooth(col = "blue",linetype=1,span = 0.5, method = "loess", size = 0.5) + xlab('Year')+
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())+theme(axis.text.x = element_blank())+theme(legend.position = "none")+theme(legend.title=element_text(colour="black"))

#for DLJ 1800
tem.precipitation.dlj.1800<-filter(tem.precipitation.dlj,tem.precipitation.dlj$glgs_DEM_30m_area_4326.x>1800)
tem.precipitation.dlj.1800.year<-filter(tem.precipitation.dlj.1800,tem.precipitation.dlj.1800$Date>as.Date("2015-08-18"))
tem.precipitation.dlj.1800.year<-filter(tem.precipitation.dlj.1800.year,tem.precipitation.dlj.1800.year$Date<as.Date("2016-08-02"))

tem.precipitation.dlj.1800.year$month <- factor(lubridate::month(tem.precipitation.dlj.1800.year$Date), levels=c("3", "4", "5","6","7","8","9","10","11","12","1","2"), ordered=TRUE)
tem.precipitation.dlj.1800.year.order<-tem.precipitation.dlj.1800.year[order(tem.precipitation.dlj.1800.year$month, tem.precipitation.dlj.1800.year$Date),]
tem.precipitation.dlj.1800.year.order$order<-c("NULL")
tem.precipitation.dlj.1800.year.order[which(lubridate::month(tem.precipitation.dlj.1800.year.order$Date)>=3),]$order <- paste0("2014/",lubridate::month(tem.precipitation.dlj.1800.year.order[which(lubridate::month(tem.precipitation.dlj.1800.year.order$Date)>=3),]$Date),"/",lubridate::day(tem.precipitation.dlj.1800.year.order[which(lubridate::month(tem.precipitation.dlj.1800.year.order$Date)>=3),]$Date))
tem.precipitation.dlj.1800.year.order[which(lubridate::month(tem.precipitation.dlj.1800.year.order$Date)<3),]$order <- paste0("2015/",lubridate::month(tem.precipitation.dlj.1800.year.order[which(lubridate::month(tem.precipitation.dlj.1800.year.order$Date)<3),]$Date),"/",lubridate::day(tem.precipitation.dlj.1800.year.order[which(lubridate::month(tem.precipitation.dlj.1800.year.order$Date)<3),]$Date))

# ggplot(tem.precipitation.dlj.1800.year.order, mapping = aes(x = factor(as.Date(order)), y = precipitation, group = 1)) + geom_point(col = "grey")+geom_smooth(col = "blue",linetype=1,span = 0.5, method = "loess", size = 0.5) + xlab('Year')+
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())+theme(axis.text.x = element_blank())+theme(legend.position = "left")+theme(legend.title=element_text(colour="black"))

tem.precipitation.dlj.1800.year.order$group<-c("DLJ")
tem.precipitation.bhl.1800.year.order$group<-c("BHL")
tem.precipitation<-rbind(tem.precipitation.bhl.1800.year.order,tem.precipitation.dlj.1800.year.order)

climate.plot<-ggplot(tem.precipitation, mapping = aes(x = factor(as.Date(order)), y = precipitation, group = group)) + geom_point(tem.precipitation[which(tem.precipitation$group=="BHL"),],mapping = aes(x = factor(as.Date(order)), y = precipitation),col = "#D55E00",shape = 0,alpha = 1/20)+
  geom_point(tem.precipitation[which(tem.precipitation$group=="DLJ"),],mapping = aes(x = factor(as.Date(order)), y = precipitation),col = "#009E73",shape = 6,alpha = 1/20)+
  geom_smooth(tem.precipitation[which(tem.precipitation$group=="BHL"),],mapping = aes(x = factor(as.Date(order)), y = precipitation),col = "#D55E00",linetype=3, method = "loess", size = 0.5, fill="#D55E00") +
  geom_smooth(tem.precipitation[which(tem.precipitation$group=="DLJ"),],mapping = aes(x = factor(as.Date(order)), y = precipitation),col = "#009E73",linetype=3, method = "loess", size = 0.5,fill = "#009E73") + xlab('Day')+
  geom_point(tem.precipitation[which(tem.precipitation$group=="BHL"),],mapping = aes(x = factor(as.Date(order)), y = temperture/0.4),col = "#D55E00",shape = 4,alpha = 1/20)+
  geom_point(tem.precipitation[which(tem.precipitation$group=="DLJ"),],mapping = aes(x = factor(as.Date(order)), y = temperture/0.4),col = "#009E73",shape = 3,alpha = 1/20)+
  geom_smooth(tem.precipitation[which(tem.precipitation$group=="BHL"),],mapping = aes(x = factor(as.Date(order)), y = temperture/0.4),col = "#D55E00",linetype=1, method = "loess", size = 0.5, fill="#D55E00") +
  geom_smooth(tem.precipitation[which(tem.precipitation$group=="DLJ"),],mapping = aes(x = factor(as.Date(order)), y = temperture/0.4),col = "#009E73",linetype=1, method = "loess", size = 0.5,fill = "#009E73") + xlab('Day')+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+theme(axis.text.x = element_blank())+theme(legend.position = "none")+theme(legend.title=element_text(colour="black"))+ylab("Daily average precipitation")+scale_y_continuous(sec.axis = sec_axis(~.*0.4, name = "Daily average temperture"))

#climate.plot

climate.plot.ref<-climate.plot + geom_vline(xintercept = factor(as.Date(c("2014-05-31","2014-08-31","2014-11-30"))))
climate.plot.ref

# ggplot(tem.precipitation, mapping = aes(x = factor(as.Date(order)), y = precipitation, group = group)) + geom_point(tem.precipitation[which(tem.precipitation$group=="BHL"),],mapping = aes(x = factor(as.Date(order)), y = precipitation),col = "#D55E00",shape = 0,alpha = 1/20)+
#   geom_point(tem.precipitation[which(tem.precipitation$group=="DLJ"),],mapping = aes(x = factor(as.Date(order)), y = precipitation),col = "#009E73",shape = 6,alpha = 1/20)+
#   geom_smooth(tem.precipitation[which(tem.precipitation$group=="BHL"),],mapping = aes(x = factor(as.Date(order)), y = precipitation),col = "#D55E00",linetype=3, method = "loess", size = 0.5, fill="#D55E00") +
#   geom_smooth(tem.precipitation[which(tem.precipitation$group=="DLJ"),],mapping = aes(x = factor(as.Date(order)), y = precipitation),col = "#009E73",linetype=3, method = "loess", size = 0.5,fill = "#009E73") + xlab('Day')+
#   geom_point(tem.precipitation[which(tem.precipitation$group=="BHL"),],mapping = aes(x = factor(as.Date(order)), y = temperture/0.6),col = "#D55E00",shape = 4,alpha = 1/20)+
#   geom_point(tem.precipitation[which(tem.precipitation$group=="DLJ"),],mapping = aes(x = factor(as.Date(order)), y = temperture/0.6),col = "#009E73",shape = 3,alpha = 1/20)+
#   geom_smooth(tem.precipitation[which(tem.precipitation$group=="BHL"),],mapping = aes(x = factor(as.Date(order)), y = temperture/0.6),col = "#D55E00",linetype=1, method = "loess", size = 0.5, fill="#D55E00") +
#   geom_smooth(tem.precipitation[which(tem.precipitation$group=="DLJ"),],mapping = aes(x = factor(as.Date(order)), y = temperture/0.6),col = "#009E73",linetype=1, method = "loess", size = 0.5,fill = "#009E73") + xlab('Day')+theme(legend.position = "left")+theme(legend.title=element_text(colour="black"))

dev.new()
ggdraw() + draw_plot(whole.labeltree.1800.p23, x = 0, y = 0, width = 1, height = 1) +
  draw_plot(climate.plot.ref, x = 0.057, y = 0.06, width = 0.873, height = .09) 
ggsave("figS1_climate_1800_20211213.pdf", width = 80, height = 80, units = "cm")
dev.off()

# dev.new()
# svg(file="./fig2_20211213.svg",width=80, height=80)
# ggdraw() + draw_plot(whole.labeltree.1800.p23, x = 0, y = 0, width = 1, height = 1) +
#   draw_plot(climate.plot.ref, x = 0.057, y = 0.06, width = 0.873, height = .09)
# dev.off()
