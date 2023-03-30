library("tidyr")
library("ggplot2")
library(data.table)
library(stringr)
library(plyr)
library(dplyr)

########################################################################################################################################################################################   
######################################################################################   Fig.4A   ######################################################################################   
########################################################################################################################################################################################

IDsignatures=fread("ID83_S4_Signatures.txt")
ID_mutation_classes=IDsignatures$MutationType
IDsignatures=IDsignatures[,-1]
ID_sigminer_signatures=as.data.frame(IDsignatures)
colnames(ID_sigminer_signatures)=c("Sig 1","Sig 2","Sig 3","Sig 4")
rownames(ID_sigminer_signatures)=ID_mutation_classes
library(sigminer)
show_sig_profile(as.matrix(ID_sigminer_signatures),mode="ID",style="cosmic",sig_names=colnames(IDsignatures),sig_orders = c(1,3,2,4))+
  annotate("text", x=0.95, y=0.952, hjust=1, label= "Stability 0.99") +
  annotate("text", x=0.95, y=0.932, hjust=1, label= "Sig. Mutations 20,946/61.9%")+
  annotate("text", x=0.95, y=0.722, hjust=1, label= "Stability 0.63") +
  annotate("text", x=0.95, y=0.702, hjust=1, label= "Sig. Mutations 3,405/10.1%")+
  annotate("text", x=0.95, y=0.492, hjust=1, label= "Stability 1.0") +
  annotate("text", x=0.95, y=0.472, hjust=1, label= "Sig. Mutations 6,745/19.9%")+
  annotate("text", x=0.95, y=0.262, hjust=1, label= "Stability 0.97") +
  annotate("text", x=0.95, y=0.242, hjust=1, label= "Sig. Mutations 2,767/8.2%")
ggsave("Fig4A.pdf",plot = last_plot(),device="pdf", width = 5000,height = 3000,units = "px")
########################################################################################################################################################################################
########################################################################################################################################################################################


########################################################################################################################################################################################   
######################################################################################   Fig.4B   ######################################################################################   
########################################################################################################################################################################################

ID_activities=fread("ID83_S4_NMF_Activities.txt")
colnames(ID_activities)=c("Samples","ID83A","ID83C","ID83B","ID83D") #change names to match the original set
ID_activities %>% select(sort(names(.))) -> ID_activities
ID_activities=ID_activities[1:21,]
ID_activities=cbind(Cell=substr(ID_activities$Samples,1,4),ID_activities)
ID_activities$Samples=sapply(str_split(ID_activities$Samples,"_"), "[[", 2)
ID_activities$Samples=mapvalues(ID_activities$Samples,from=c("CA","CC","HE1A","HE2B","PR1B","PR2B","CB","HE2A","PR1A","PR2A","C","HE1B","HE1C","PR2C"),
                                    to=c("Control_A","Control_C","Helium_1A","Helium_2B","Protons_1B","Protons_2B","Control_B","Helium_2A","Protons_1A","Protons_2A","Control_A","Helium_1B","Helium_1C","Protons_2C"))
ID_activities=cbind(Clone=paste(ID_activities$Cell,"  -  ",ID_activities$Samples,sep=""),ID_activities)

order_df=data.frame(Mutations_total=rowSums(ID_activities[,3:6]),Clone=ID_activities$Clone)
clone_order_decreasing=order_df[order(order_df$Mutations_total,decreasing=T),"Clone"]
clone_order_normal=rev(ID_activities$Clone)
pivot_longer(ID_activities, cols=3:6, names_to = "Signature", values_to = "Mutations") -> ID_activities
#ID_activities$Clone <- factor(ID_activities$Clone, levels = clone_order_decreasing) #this orders by total mutations
ID_activities$Cell <- factor(ID_activities$Cell, levels = c("HAP1", "A549", "MCF7"))

###ID_activities=ID_activities[-which(ID_activities$Clone=="MCF7  -  Helium_2B"),]  #when removing the HE2B clone
ggplot(ID_activities, aes(fill=Signature, y=Mutations, x=Samples)) +
  theme_bw()+
  geom_bar(position="stack", stat="identity",colour="black")+
  scale_y_continuous(breaks = seq(0, 700, by = 100))+
  facet_wrap(~Cell, nrow=1,scales = "free_x")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust=1),
        text = element_text(size = 20),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"))+
  ylab("Number of attributed mutations\n")+
  xlab("\nClone")+
  scale_fill_manual(values=c("#E69F00", "#CC79A7","#0072B2","#009E73")) #change to this when you change the names
ggsave("Fig4B.pdf",plot = last_plot(),device="pdf", width = 5000,height = 2000,units = "px")
########################################################################################################################################################################################
########################################################################################################################################################################################


########################################################################################################################################################################################   
####################################################################################   ID decomp.   ####################################################################################   
########################################################################################################################################################################################

ID_decomposition=data.frame(
  Signature=c("ID83A","ID83B","ID83B","ID83B","ID83B","ID83C","ID83D","ID83D"),
  Component=c("ID2","ID12","ID8","ID4","ID2","ID1","ID3","ID9"),
  Contribution=c(100.00,37.86,27.58,23.76,10.8,100.00,79.78,20.22),
  stringsAsFactors = F)

ID_decomposition$Component <- factor(ID_decomposition$Component, levels = unique(mixedsort(ID_decomposition$Component)))
colours=c("#E68FAC","#0067A5","#8DB600","#008856","#A1CAF1","#B3446C","#848482")
ggplot(ID_decomposition, aes(fill=Component, y=Contribution, x=Signature)) +
  geom_bar(position="stack", stat="identity",colour="black")+
  scale_fill_manual(values=colours)+
  scale_y_continuous(limits = c(0, 105), breaks = c(0,25,50,75,100))+
  geom_text(aes(label = Component), position = position_stack(vjust = 0.5),size=5.5) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 25),
        plot.title = element_text(hjust = 0.5,size=20,vjust=0),
        panel.background=element_rect(fill="white"),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"),
        legend.position = "bottom",
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))+
  ylab("Contribution (%)")+
  labs(fill="COSMIC signatures")+
  annotate("text", x=1, y=105, label= "0.999",size=6) + 
  annotate("text", x = 2, y=105, label = "0.781",size=6)+
  annotate("text", x = 3, y=105, label = "0.998",size=6)+
  annotate("text", x = 4, y=105, label = "0.956",size=6)+
  ggtitle("Cosine similarity with reconstructed signature")+
  guides(fill = guide_legend(nrow = 1,title.position="top", title.hjust = 0.5))
ggsave("Fig4E.pdf",plot = last_plot(),device="pdf", width = 2700,height = 3000,units = "px")

########################################################################################################################################################################################
########################################################################################################################################################################################


########################################################################################################################################################################################   
####################################################################################   ID heatmap   ####################################################################################   
########################################################################################################################################################################################

library(MutationalPatterns)
library(data.table)
COSMIC_ID=fread("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/COSMIC_signatures/COSMIC_v3.3_ID_GRCh37.txt")
COSMIC_ID=COSMIC_ID[,-1]
IDsignatures=fread("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/indels/Signatures/All_samples/ID83/All_Solutions/ID83_4_Signatures/Signatures/ID83_S4_Signatures.txt")
colnames(IDsignatures)=c("MutationType","ID83A","ID83C","ID83B","ID83D") #switch SBSB with SBSC
IDsignatures=cbind(IDsignatures[,1],IDsignatures[,2],IDsignatures[,4],IDsignatures[,3],IDsignatures[,5]) #reorder
IDsignatures=IDsignatures[,-1]

csm = cos_sim_matrix(COSMIC_ID, IDsignatures) 
#cluster cosmic signatures
hc.sample <- hclust(dist(csm), method = "complete")
col_order <- rownames(csm)[hc.sample$order]
#make dendrogram
dhc <- as.dendrogram(hc.sample)
ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
dendrogram_rows <- ggplot(ggdendro::segment(ddata)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  ggdendro::theme_dendro()
dendrogram_rows
ggsave("ID_dendrogram.pdf",plot = last_plot(),device="pdf", width = 2300,height = 300,units = "px")

#create cosine similarity df
csm=as.data.frame(csm)
csm$COSMIC_sig=rownames(csm)
pivot_longer(csm, cols=1:4, names_to = "Signature", values_to = "Similarity") -> csm
csm$COSMIC_sig <- factor(csm$COSMIC_sig, levels = col_order)
#plot
ggplot(data = csm, aes(x =Signature, y =COSMIC_sig)) +
  geom_raster(aes(fill = Similarity))+
  scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Cosine similarity", limits = c(0, 1.000000001)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=1, size = 15),
        axis.text.y = element_text(angle=0,hjust = 0, vjust = 0.5, size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.ticks = element_blank(),
        #legend.title = element_text(angle=270),
        legend.position = "bottom",
        axis.title = element_blank(),
        legend.key.width = unit(1.5, "cm"),
        text = element_text(size = 20))+
  scale_x_discrete(position = "top") +
  guides(fill = guide_colourbar(raster=T,title.position="top", title.hjust = 0.5,nbin=100))+
  labs(x = NULL, y = NULL)
ggsave("ID_heatmap.pdf",plot = last_plot(),device="pdf", width = 1500,height = 2300,units = "px")
########################################################################################################################################################################################
########################################################################################################################################################################################


########################################################################################################################################################################################   
####################################################################################   Supp Fig.1B   ###################################################################################   
########################################################################################################################################################################################

library(xlsx)
library(stringr)
naming=c("CA"="Control_A","CC"="Control_C","HE1A"="Helium_1A","HE2B"="Helium_2B","PR1B"="Protons_1B","PR2B"="Protons_2B","CB"="Control_B","HE2A"="Helium_2A","PR1A"="Protons_1A",
         "PR2A"="Protons_2A","C"="Control_A","HE1B"="Helium_1B","HE1C"="Helium_1C","PR2C"="Protons_2C")

Zou = read.table("denovo_indels_43genes_Zou2021.txt", h=T, stringsAsFactors = F, sep="\t")
gene_KOs_with_sig=c("OGG1","UNG","EXO1","RNF168","MLH1","MSH2","MSH6","PMS2","PMS1") #genes marked as having signatures in Zou 2021
Zou = Zou[which(Zou$Ko_gene %in% gene_KOs_with_sig),]
Kucab_supp = read.xlsx("suppTable2_Kucab2019.xlsx", sheetIndex = 1)

ID_activities=fread("ID83_S4_NMF_Activities.txt")
colnames(ID_activities)=c("Samples","ID83A","ID83C","ID83B","ID83D") #change names to match the original set
ID_activities=cbind(ID_activities[,1],ID_activities[,2],ID_activities[,4],ID_activities[,3],ID_activities[,5]) #reorder
ID_activities$Samples=gsub("MSM0_","MSM0.",ID_activities$Samples); ID_activities$Samples=gsub("MSK0_","MSK0.",ID_activities$Samples)
ID_activities_row_norm=ID_activities
ID_activities_row_norm[, 2:5]=sweep(ID_activities[, 2:5],MARGIN=1,FUN="/",STATS=rowSums(ID_activities[, 2:5])) #normalise by row to get % contribution per sample
ID83A=ID_activities_row_norm[which(ID_activities_row_norm$ID83A>=0.65),c("Samples","ID83A")]
colnames(ID83A)=c("Samples","Proportion"); ID83A$Signature=rep("ID83A",nrow(ID83A))
ID83B=ID_activities_row_norm[which(ID_activities_row_norm$ID83B>=0.65),c("Samples","ID83B")]
colnames(ID83B)=c("Samples","Proportion"); ID83B$Signature=rep("ID83B",nrow(ID83B))
ID83C=ID_activities_row_norm[which(ID_activities_row_norm$ID83C>=0.65),c("Samples","ID83C")]
colnames(ID83C)=c("Samples","Proportion"); ID83C$Signature=rep("ID83C",nrow(ID83C))
ID83D=ID_activities_row_norm[which(ID_activities_row_norm$ID83D>=0.65),c("Samples","ID83D")]
colnames(ID83D)=c("Samples","Proportion"); ID83D$Signature=rep("ID83D",nrow(ID83D))
ID_activities_by_treatment=rbind(ID83A,ID83B,ID83C,ID83D)
ID_activities_by_treatment$Treatment=NA
ID_activities_by_treatment$Study=NA

for (i in 1:nrow(ID_activities_by_treatment)){
    if (sapply(str_split(ID_activities_by_treatment$Samples,"_"), "[[", 1)[i] %in% c("A549","HAP1","MCF7")){
      ID_activities_by_treatment$Treatment[i]=paste(sapply(str_split(ID_activities_by_treatment$Samples,"_"), "[[", 1)[i]," - ",naming[sapply(str_split(ID_activities_by_treatment$Samples,"_"), "[[", 2)[i]],sep="")
      ID_activities_by_treatment$Study[i]="Our study"}
    if (ID_activities_by_treatment$Samples[i] %in% unique(Zou$Sample)){
      ID_activities_by_treatment$Treatment[i]=Zou[which(Zou$Sample==ID_activities_by_treatment$Samples[i]),"Ko_gene"][1]
      ID_activities_by_treatment$Study[i]="Zou"}
    if (sapply(str_split(ID_activities_by_treatment$Samples,"_"), "[[", 1)[i] %in% Kucab_supp$Sample.Name){
      ID_activities_by_treatment$Study[i]="Kucab"
      if(is.na(as.character(Kucab_supp[which(Kucab_supp$Sample.Name==sapply(str_split(ID_activities_by_treatment$Samples,"_"), "[[", 1)[i]),"Compound.Abbreviation."]))){
      ID_activities_by_treatment$Treatment[i]=as.character(Kucab_supp[which(Kucab_supp$Sample.Name==sapply(str_split(ID_activities_by_treatment$Samples,"_"), "[[", 1)[i]),"Compound"])}
      else{ID_activities_by_treatment$Treatment[i]=as.character(Kucab_supp[which(Kucab_supp$Sample.Name==sapply(str_split(ID_activities_by_treatment$Samples,"_"), "[[", 1)[i]),"Compound.Abbreviation."])}}}
ID_activities_by_treatment$Study <- factor(ID_activities_by_treatment$Study, levels = c("Kucab","Zou","Our study"))

ggplot(data = ID_activities_by_treatment, aes(x=Proportion, y=Treatment, color=Study)) + 
  geom_point(size=6,shape=18)+
  scale_color_manual(values=c("#AA4499","#6699CC","#332288"))+
  facet_wrap(~Signature,ncol=4,scales = "free_y")+
  xlab("\nProportion of total mutations contributed to by given signature")+
  ylab("Sample treatment\n")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5, size = 17),
        axis.text.y = element_text(vjust = 0.5, hjust=1, size = 17),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(size=20),
        legend.title = element_text(size=20),
        legend.position = "bottom",
        strip.text = element_text(size = 20),
        text = element_text(size = 20),
        panel.spacing = unit(2, "lines"))
ggsave("SFig3B.pdf",plot = last_plot(),device="pdf",width = 7000,height = 3000,units = "px")
########################################################################################################################################################################################
########################################################################################################################################################################################

