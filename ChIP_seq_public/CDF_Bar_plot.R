library(ggplot2)
setwd("/share/home/zhangkx/uploadFromHuicao/Figure3/")
# read data
known_withPromoter<-read.table("./known_withPromoter_Intersect_S1ENDSeqHDNA.txt",header = F,sep = "\t")
known_withPromoter_shared<-read.table("./known_withPromoter_Intersect_S1ENDSeqHDNA_shared.txt",header = F,sep = "\t")
known_withoutPromoter<-read.table("./known_withoutPromoter_Intersect_S1ENDSeqHDNA_shared.txt",header = F,sep = "\t")
ms_withPromoter<-read.table("./MS_withPromoter_Intersect_S1ENDSeqHDNA.txt",header = F,sep = "\t")
ms_withPromoter_shared<-read.table("./ms_withPromoter_Intersect_S1ENDSeqHDNA_shared.txt",header = F,sep = "\t")
ms_withoutPromoter<-read.table("./MS_withoutPromoter_Intersect_S1ENDSeqHDNA_shared.txt",header = F,sep = "\t")
k562_withPromoter<-read.table("./TF_K562_withPromoter_Intersect_S1ENDSeq.txt",header = F,sep = "\t")
k562_withPromoter_shared<-read.table("./TF_K562_withPromoter_Intersect_S1ENDSeq_shared.txt",header = F,sep = "\t")
k562_withoutPromoter<-read.table("./TF_K562_withoutPromoter_Intersect_S1ENDSeq_shared.txt",header = F,sep = "\t")
hepg2_withPromoter<-read.table("./TF_HepG2_withPromoter_Intersect_S1ENDSeq.txt",header = F,sep = "\t")
hepg2_withPromoter_shared<-read.table("./TF_HepG2_withPromoter_Intersect_S1ENDSeq_shared.txt",header = F,sep = "\t")
hepg2_withoutPromoter<-read.table("TF_HepG2_withoutPromoter_Intersect_S1ENDSeq_shared.txt",header = F,sep = "\t")


k562_withPromoter$V1 <- gsub("-human", "", k562_withPromoter$V1)
k562_withPromoter_shared$V1 <- gsub("-human", "", k562_withPromoter_shared$V1)
k562_withoutPromoter$V1 <- gsub("-human", "", k562_withoutPromoter$V1)
hepg2_withPromoter$V1 <- gsub("-human", "", hepg2_withPromoter$V1)
hepg2_withPromoter_shared$V1 <- gsub("-human", "", hepg2_withPromoter_shared$V1)
hepg2_withoutPromoter$V1 <- gsub("-human", "", hepg2_withoutPromoter$V1)


tf_withPromoter<-merge(k562_withPromoter, hepg2_withPromoter, by = "V1", all = TRUE)
tf_withPromoter[is.na(tf_withPromoter)] <- 0
tf_withPromoter_shared<-merge(k562_withPromoter_shared, hepg2_withPromoter_shared, by = "V1", all = TRUE)
tf_withPromoter_shared[is.na(tf_withPromoter_shared)] <- 0
tf_withoutPromoter<-merge(k562_withoutPromoter, hepg2_withoutPromoter, by = "V1", all = TRUE)
tf_withoutPromoter[is.na(tf_withoutPromoter)] <- 0

# calculate count and total (need total number)
tf_withPromoter$total <- ifelse(tf_withPromoter$V3.x > tf_withPromoter$V3.y, 
                                tf_withPromoter$V2.x, tf_withPromoter$V2.y) # total
tf_withPromoter_shared$total <- ifelse(tf_withPromoter_shared$V3.x > tf_withPromoter_shared$V3.y, 
                                       tf_withPromoter_shared$V2.x, tf_withPromoter_shared$V2.y) # total
tf_withoutPromoter$total <- ifelse(tf_withoutPromoter$V3.x > tf_withoutPromoter$V3.y, 
                                   tf_withoutPromoter$V2.x, tf_withoutPromoter$V2.y) # total
tf_withPromoter$count <- pmax(tf_withPromoter$V3.x, tf_withPromoter$V3.y, na.rm = TRUE) # max
tf_withPromoter_shared$count <- pmax(tf_withPromoter_shared$V3.x, tf_withPromoter_shared$V3.y, na.rm = TRUE) # max
tf_withoutPromoter$count <- pmax(tf_withoutPromoter$V3.x, tf_withoutPromoter$V3.y, na.rm = TRUE) # max


colnames(tf_withPromoter)[1]<-"protein"
colnames(tf_withPromoter_shared)[1]<-"protein"
colnames(tf_withoutPromoter)[1]<-"protein"

TF_withPromoter<-tf_withPromoter[,c(1,6,7)]
TF_withPromoter_shared<-tf_withPromoter_shared[,c(1,6,7)]
TF_withoutPromoter<-tf_withoutPromoter[,c(1,6,7)]

protein_filter<-c("CHD4","DDX17","DDX21","DDX5","DNMT1","HNRNPC","HNRNPK","MATR3","NUP93"
                  ,"PCBP1","PCBP2","PDS5A","PHB2","PSIP1", "SFPQ","SMARCA5","TIF1B","TOP2A","TOP2B","XRCC5",
                  "DDX11","BRIP1","WRN","BLM","DHX9","p53","RPA1","HNRNPL","HNRNPA2B1","PTBP1","ORC4",
                  "NONO","U2AF2","Vim", "GFAP","DES")

TF_withPromoter<-TF_withPromoter[!(TF_withPromoter$protein %in% protein_filter),]
TF_withPromoter_shared<-TF_withPromoter_shared[!(TF_withPromoter_shared$protein %in% protein_filter),]
TF_withoutPromoter<-TF_withoutPromoter[!(TF_withoutPromoter$protein %in% protein_filter),]

TF_withPromoter$type<-"TF_withPromoter"
TF_withPromoter_shared$type<-"TF_withPromoter_shared"
TF_withoutPromoter$type<-"TF_withoutPromoter"

TF_withPromoter$ratio<-TF_withPromoter$count/TF_withPromoter$total
TF_withPromoter_shared$ratio<-TF_withPromoter_shared$count/TF_withPromoter_shared$total
TF_withoutPromoter$ratio<-TF_withoutPromoter$count/TF_withoutPromoter$total

ms_withPromoter$type<-"ms_withPromoter"
ms_withPromoter_shared$type<-"ms_withPromoter_shared"
ms_withoutPromoter$type<-"ms_withoutPromoter"

ms_withPromoter <- ms_withPromoter[order(ms_withPromoter[, 3], decreasing = TRUE), ]
ms_withPromoter_shared <- ms_withPromoter_shared[order(ms_withPromoter_shared[, 3], decreasing = TRUE), ]
ms_withoutPromoter <- ms_withoutPromoter[order(ms_withoutPromoter[, 3], decreasing = TRUE), ]

ms_withPromoter$ratio<-ms_withPromoter$V3/ms_withPromoter$V2
ms_withPromoter_shared$ratio<-ms_withPromoter_shared$V3/ms_withPromoter_shared$V2
ms_withoutPromoter$ratio<-ms_withoutPromoter$V3/ms_withoutPromoter$V2

known_withPromoter$type<-"known_withPromoter"
known_withPromoter_shared$type<-"known_withPromoter_shared"
known_withoutPromoter$type<-"known_withoutPromoter"

known_withPromoter$ratio<-known_withPromoter$V3/known_withPromoter$V2
known_withPromoter_shared$ratio<-known_withPromoter_shared$V3/known_withPromoter_shared$V2
known_withoutPromoter$ratio<-known_withoutPromoter$V3/known_withoutPromoter$V2

colnames(ms_withPromoter)<-c("protein","total","count","type")
colnames(ms_withPromoter_shared)<-c("protein","total","count","type")
colnames(ms_withoutPromoter)<-c("protein","total","count","type")
colnames(known_withPromoter)<-c("protein","total","count","type")
colnames(known_withPromoter_shared)<-c("protein","total","count","type")
colnames(known_withoutPromoter)<-c("protein","total","count","type")
colnames(TF_withPromoter)<-c("protein","total","count","type")
colnames(TF_withPromoter_shared)<-c("protein","total","count","type")
colnames(TF_withoutPromoter)<-c("protein","total","count","type")


MS_withPromoter<-ms_withPromoter[order(-ms_withPromoter$count),]
MS_withPromoter_shared<-ms_withPromoter_shared[order(-ms_withPromoter_shared$count),]
MS_withoutPromoter<-ms_withoutPromoter[order(-ms_withoutPromoter$count),]

df1_withPromoter<-rbind(ms_withPromoter,known_withPromoter)
df1_withPromoter_shared<-rbind(ms_withPromoter_shared,known_withPromoter_shared)
df1_withoutPromoter<-rbind(ms_withoutPromoter,known_withoutPromoter)

df2_withPromoter<-rbind(df1_withPromoter,TF_withPromoter)
df2_withPromoter_shared<-rbind(df1_withPromoter_shared,TF_withPromoter_shared)
df2_withoutPromoter<-rbind(df1_withoutPromoter,TF_withoutPromoter)

#CDF
p1_withPromoter<-ggplot(df2_withPromoter,aes(x=count,color=type))+
  stat_ecdf()+
  labs(x="Number of ChIP-seq peaks \n intersected with merged HDNA",
       y="cumulative distribution funcition (CDF)")+
  scale_color_manual(values =c("#6A3D9A","#1F78B4" ,"#33A02C"))+
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.title.y = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.text = element_text(color="black",size=11),
        panel.border = element_rect(linetype = "solid",colour = "black", fill = "NA",size=1),
        panel.spacing.x = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(color="black",size=10),
        legend.spacing.x=unit(0.1,'cm'),
        legend.background=element_blank(),
        legend.position = c(1,0), legend.justification = c(1,0))

p1_withPromoter + coord_cartesian(xlim = c(0, 2000)) + scale_y_continuous(breaks = c(0, 0.2, 0.5, 0.8, 1))


p1_withPromoter_shared<-ggplot(df2_withPromoter_shared,aes(x=count,color=type))+
  stat_ecdf()+
  labs(x="Number of ChIP-seq peaks \n intersected with shared HDNA",
       y="cumulative distribution funcition (CDF)")+
  scale_color_manual(values =c("#6A3D9A","#1F78B4" ,"#33A02C"))+
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.title.y = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.text = element_text(color="black",size=11),
        panel.border = element_rect(linetype = "solid",colour = "black", fill = "NA",size=1),
        panel.spacing.x = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(color="black",size=10),
        legend.spacing.x=unit(0.1,'cm'),
        legend.background=element_blank(),
        legend.position = c(1,0), legend.justification = c(1,0))

p1_withPromoter_shared + coord_cartesian(xlim = c(0, 200)) + scale_y_continuous(breaks = c(0, 0.2, 0.5, 0.8, 1))

p1_withoutPromoter<-ggplot(df2_withoutPromoter,aes(x=count,color=type))+
  stat_ecdf()+
  labs(x="Number of ChIP-seq peaks \n intersected with shared HDNA",
       y="cumulative distribution funcition (CDF)")+
  scale_color_manual(values =c("#6A3D9A","#1F78B4" ,"#33A02C"))+
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.title.y = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.text = element_text(color="black",size=11),
        panel.border = element_rect(linetype = "solid",colour = "black", fill = "NA",size=1),
        panel.spacing.x = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(color="black",size=10),
        legend.spacing.x=unit(0.1,'cm'),
        legend.background=element_blank(),
        legend.position = c(1,0), legend.justification = c(1,0))

p1_withoutPromoter + coord_cartesian(xlim = c(0, 200)) + scale_y_continuous(breaks = c(0, 0.2, 0.5, 0.8, 1))

# ks.test
ks.test(ms_withPromoter$count,tf_withPromoter$count) #p-value = 8.166e-07
ks.test(known_withPromoter$count,tf_withPromoter$count) #p-value = 0.002757
ks.test(ms_withPromoter$count,known_withPromoter$count) #p-value = 0.7473

ks.test(ms_withoutPromoter$count,tf_withoutPromoter$count) #p-value = 5.685e-06
ks.test(known_withoutPromoter$count,tf_withoutPromoter$count) #p-value = 0.007002
ks.test(ms_withoutPromoter$count,known_withoutPromoter$count) #p-value = 0.5259

ks.test(ms_withPromoter_shared$count,tf_withPromoter_shared$count) #p-value = 1.822e-06
ks.test(known_withPromoter_shared$count,tf_withPromoter_shared$count) #p-value = 0.009219
ks.test(ms_withPromoter_shared$count,known_withPromoter_shared$count) #p-value = 0.5259


#bar
colnames(ms_withPromoter)<-c("protein","total","count","type","ratio")
colnames(ms_withPromoter_shared)<-c("protein","total","count","type","ratio")
colnames(ms_withoutPromoter)<-c("protein","total","count","type","ratio")

p2_withPromoter<-ggplot(ms_withPromoter,aes(x=protein,y=count,fill=type))+
  geom_bar(stat = 'identity',position = "stack")+scale_x_discrete(limits = ms_withPromoter$protein)+
  labs(x=NULL)+
  scale_fill_manual(values = c("#6A3D9A"))+
  theme(panel.background = element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.title.y = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.text = element_text(color="black",size=11),axis.text.x = element_text(angle = 90,
                                                                                   hjust = 1),
        panel.border = element_rect(linetype = "solid",colour = "black", fill = "NA",size=1),
        panel.spacing.x = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(color="black",size=10),
        legend.spacing.x=unit(0.1,'cm'),
        legend.background=element_blank(),
        legend.position = c(1,1), legend.justification = c(1,1))+
  geom_text(aes(x= protein,y=count,label=count),vjust=-0.5,size=3.5,fontface='bold')+
  #geom_line(aes(x= protein,y=ratio*100*60,group=1),linetype=3,cex=1)+
  #geom_point(aes(x= protein,y=ratio*100*60),color='#33A02C',size=3.5)+
  scale_y_continuous(expand = c(0,0),limits = c(0,1200),sec.axis = sec_axis(~ . / 60,
                                                                            name = "Ratio"),name = "Count")
p2_withPromoter

p2_withPromoter_shared<-ggplot(ms_withPromoter_shared,aes(x=protein,y=count,fill=type))+
  geom_bar(stat = 'identity',position = "stack")+scale_x_discrete(limits = ms_withPromoter_shared$protein)+
  labs(x=NULL)+
  scale_fill_manual(values = c("#6A3D9A"))+
  theme(panel.background = element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.title.y = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.text = element_text(color="black",size=11),axis.text.x = element_text(angle = 90,
                                                                                   hjust = 1),
        panel.border = element_rect(linetype = "solid",colour = "black", fill = "NA",size=1),
        panel.spacing.x = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(color="black",size=10),
        legend.spacing.x=unit(0.1,'cm'),
        legend.background=element_blank(),
        legend.position = c(1,1), legend.justification = c(1,1))+
  geom_text(aes(x= protein,y=count,label=count),vjust=-0.5,size=3.5,fontface='bold')+
  #geom_line(aes(x= protein,y=ratio*100*60,group=1),linetype=3,cex=1)+
  #geom_point(aes(x= protein,y=ratio*100*60),color='#33A02C',size=3.5)+
  scale_y_continuous(expand = c(0,0),limits = c(0,65),sec.axis = sec_axis(~ . / 60,
                                                                            name = "Ratio"),name = "Count")
p2_withPromoter_shared


p2_withoutPromoter<-ggplot(ms_withoutPromoter,aes(x=protein,y=count,fill=type))+
  geom_bar(stat = 'identity',position = "stack")+scale_x_discrete(limits = ms_withoutPromoter$protein)+
  labs(x=NULL)+
  scale_fill_manual(values = c("#6A3D9A"))+
  theme(panel.background = element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.title.y = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.text = element_text(color="black",size=11),axis.text.x = element_text(angle = 90,
                                                                                   hjust = 1),
        panel.border = element_rect(linetype = "solid",colour = "black", fill = "NA",size=1),
        panel.spacing.x = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(color="black",size=10),
        legend.spacing.x=unit(0.1,'cm'),
        legend.background=element_blank(),
        legend.position = c(1,1), legend.justification = c(1,1))+
  geom_text(aes(x= protein,y=count,label=count),vjust=-0.5,size=3.5,fontface='bold')+
  #geom_line(aes(x= protein,y=ratio*100*60,group=1),linetype=3,cex=1)+
  #geom_point(aes(x= protein,y=ratio*100*60),color='#33A02C',size=3.5)+
  scale_y_continuous(expand = c(0,0),limits = c(0,65),sec.axis = sec_axis(~ . / 60,
                                                                            name = "Ratio"),name = "Count")
p2_withoutPromoter


ggsave("TF_MS_known_withPromoter_merged_peak.pdf")
ggsave("TF_MS_known_withPromoter_shared_peak.pdf")
ggsave("TF_MS_known_withoutPromoter_shared_peak.pdf")  # 3.99*3.28
ggsave("MS_withPromoter_merged_peak.pdf")
ggsave("MS_withPromoter_shared_peak.pdf")
ggsave("MS_withoutPromoter_shared_peak.pdf")
ggsave("TF_MS_known_withoutPromoter_shared_peak_x200.pdf",width =3.44, height =3.25) 

ggsave("MS_withPromoter_merged_peak_notdot.pdf")
ggsave("MS_withPromoter_shared_peak_notdot.pdf")
ggsave("MS_withoutPromoter_shared_peak_notdot.pdf")
