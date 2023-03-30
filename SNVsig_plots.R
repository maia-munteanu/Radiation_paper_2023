
########################################################################################################################################################################################   
####################################################################################   Supp Fig.3A   ###################################################################################   
########################################################################################################################################################################################

naming=c("CA"="Control_A","CC"="Control_C","HE1A"="Helium_1A","HE2B"="Helium_2B","PR1B"="Protons_1B","PR2B"="Protons_2B","CB"="Control_B","HE2A"="Helium_2A","PR1A"="Protons_1A",
         "PR2A"="Protons_2A","C"="Control_A","HE1B"="Helium_1B","HE1C"="Helium_1C","PR2C"="Protons_2C")

Zou = read.table("denovo_indels_43genes_Zou2021.txt", h=T, stringsAsFactors = F, sep="\t") #From Zou 2021
gene_KOs_with_sig=c("OGG1","UNG","EXO1","RNF168","MLH1","MSH2","MSH6","PMS2","PMS1") #genes marked as having signatures in Zou 2021
Zou = Zou[which(Zou$Ko_gene %in% gene_KOs_with_sig),]
Kucab_supp = read.xlsx("suppTable2_Kucab2019.xlsx", sheetIndex = 1) #from Kucab 2019

SNV_activities=fread("SBS96_S10_NMF_Activities.txt")
colnames(SNV_activities)=c("Samples","SBS96B","SBS96A","SBS96F","SBS96C","SBS96H","SBS96I","SBS96E","SBS96J","SBS96G","SBS96D")
SNV_activities %>% dplyr::select(sort(names(.))) -> SNV_activities
SNV_activities_row_norm=as.data.frame(SNV_activities)
SNV_activities_row_norm[, 2:11]=sweep(SNV_activities[, 2:11],MARGIN=1,FUN="/",STATS=rowSums(SNV_activities[, 2:11]))

for (i in colnames(SNV_activities)[colnames(SNV_activities)!="Samples"]){
    df=SNV_activities_row_norm[which(SNV_activities_row_norm[,i]>=0.45),c("Samples",i)]
    colnames(df)=c("Samples","Proportion"); df$Signature=rep(i,nrow(df))
    assign(i,df)}

SNV_activities_by_treatment=rbind(SBS96A,SBS96B,SBS96C,SBS96D,SBS96E,SBS96F,SBS96G,SBS96H,SBS96I,SBS96J)
SNV_activities_by_treatment$Treatment=NA
SNV_activities_by_treatment$Study=NA

for (i in 1:nrow(SNV_activities_by_treatment)){
  if (sapply(str_split(SNV_activities_by_treatment$Samples,"_"), "[[", 1)[i] %in% c("A549","HAP1","MCF7")){
    SNV_activities_by_treatment$Treatment[i]=paste(sapply(str_split(SNV_activities_by_treatment$Samples,"_"), "[[", 1)[i]," - ",naming[sapply(str_split(SNV_activities_by_treatment$Samples,"_"), "[[", 2)[i]],sep="")
    SNV_activities_by_treatment$Study[i]="Our study"}
  if (SNV_activities_by_treatment$Samples[i] %in% unique(Zou$Sample)){
    SNV_activities_by_treatment$Treatment[i]=Zou[which(Zou$Sample==SNV_activities_by_treatment$Samples[i]),"Ko_gene"][1]
    SNV_activities_by_treatment$Study[i]="Zou"}
  if (sapply(str_split(SNV_activities_by_treatment$Samples,"_"), "[[", 1)[i] %in% Kucab_supp$Sample.Name){
    SNV_activities_by_treatment$Study[i]="Kucab"
    if(is.na(as.character(Kucab_supp[which(Kucab_supp$Sample.Name==sapply(str_split(SNV_activities_by_treatment$Samples,"_"), "[[", 1)[i]),"Compound.Abbreviation."]))){
      SNV_activities_by_treatment$Treatment[i]=as.character(Kucab_supp[which(Kucab_supp$Sample.Name==sapply(str_split(SNV_activities_by_treatment$Samples,"_"), "[[", 1)[i]),"Compound"])}
    else{SNV_activities_by_treatment$Treatment[i]=as.character(Kucab_supp[which(Kucab_supp$Sample.Name==sapply(str_split(SNV_activities_by_treatment$Samples,"_"), "[[", 1)[i]),"Compound.Abbreviation."])}}}
SNV_activities_by_treatment$Study <- factor(SNV_activities_by_treatment$Study, levels = c("Kucab","Zou","Our study"))

ggplot(data = SNV_activities_by_treatment, aes(x=Proportion, y=Treatment, color=Study)) + 
  geom_point(size=6,shape=18)+
  scale_x_continuous(breaks=c(0.5,0.6,0.7,0.8,0.9,1))+
  scale_color_manual(values=c("#AA4499","#6699CC","#332288"))+
  facet_wrap(~Signature,ncol=5,scales = "free_y")+
  xlab("\nProportion of total mutations contributed to by given signature")+
  ylab("Sample treatment\n")+
  theme_bw()+
  theme(#panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5, size = 17),
        axis.text.y = element_text(vjust = 0.5, hjust=1, size = 17),
        axis.ticks = element_blank(),
        axis.title = element_text(size=20),
        axis.title.y = element_blank(),
        legend.title = element_text(size=20),
        legend.position = "bottom",
        strip.text = element_text(size = 20),
        text = element_text(size = 20),
        panel.spacing = unit(2, "lines"))
ggsave("SFig3A.pdf",plot = last_plot(),device="pdf", width = 8000,height = 5000,units = "px")

########################################################################################################################################################################################
########################################################################################################################################################################################


########################################################################################################################################################################################   
####################################################################################   Supp Fig.2A   ####################################################################################   
########################################################################################################################################################################################

library(sigminer)
SNVsignatures=fread("SBS96_S10_Signatures.txt")
SNV_mutation_classes=SNVsignatures$MutationType
colnames(SNVsignatures)=c("MutationType","SBS96B","SBS96A","SBS96F","SBS96C","SBS96H","SBS96I","SBS96E","SBS96J","SBS96G","SBS96D")
SNVsignatures %>% dplyr::select(sort(names(.))) -> SNVsignatures
SNVsignatures=SNVsignatures[,-1]
SNV_sigminer_signatures=as.data.frame(SNVsignatures)

colnames(SNV_sigminer_signatures)=paste0("Sig",1:10)
rownames(SNV_sigminer_signatures)=SNV_mutation_classes

show_sig_profile(as.matrix(SNV_sigminer_signatures),mode="SBS",style="cosmic",sig_names=colnames(SNVsignatures))+
  annotate("text", x=0.93, y=0.98, hjust=1, label= "Stability 0.96") +
  annotate("text", x=0.93, y=0.972, hjust=1, label= "Sig. Mutations 34,323/14.0%")+
  annotate("text", x=0.93, y=0.883, hjust=1, label= "Stability 1.0") +
  annotate("text", x=0.93, y=0.875, hjust=1, label= "Sig. Mutations 35,523/14.5%")+
  annotate("text", x=0.93, y=0.786, hjust=1, label= "Stability 0.94") +
  annotate("text", x=0.93, y=0.778, hjust=1, label= "Sig. Mutations 25,173/10.3%")+
  annotate("text", x=0.93, y=0.689, hjust=1, label= "Stability 0.85") +
  annotate("text", x=0.93, y=0.681, hjust=1, label= "Sig. Mutations 14,814/6.1%")+
  annotate("text", x=0.93, y=0.592, hjust=1, label= "Stability 0.98") +
  annotate("text", x=0.93, y=0.584, hjust=1, label= "Sig. Mutations 21,444/8.8%")+
  annotate("text", x=0.93, y=0.495, hjust=1, label= "Stability 0.97") +
  annotate("text", x=0.93, y=0.487, hjust=1, label= "Sig. Mutations 31,814/13.0%")+
  annotate("text", x=0.93, y=0.398, hjust=1, label= "Stability 0.99") +
  annotate("text", x=0.93, y=0.390, hjust=1, label= "Sig. Mutations 16,052/6.6%")+
  annotate("text", x=0.93, y=0.301, hjust=1, label= "Stability 0.75") +
  annotate("text", x=0.93, y=0.293, hjust=1, label= "Sig. Mutations 23,621/9.7%")+
  annotate("text", x=0.93, y=0.204, hjust=1, label= "Stability 0.79") +
  annotate("text", x=0.93, y=0.196, hjust=1, label= "Sig. Mutations 22,859/9.4%")+
  annotate("text", x=0.93, y=0.107, hjust=1, label= "Stability 0.74") +
  annotate("text", x=0.93, y=0.099, hjust=1, label= "Sig. Mutations 18,669/7.6%")
ggsave("SFig2A.pdf",plot = last_plot(),device="pdf", width = 4500,height = 8000,units = "px")

########################################################################################################################################################################################
########################################################################################################################################################################################


########################################################################################################################################################################################   
##################################################################################   Fig.3B   ##########################################################################################   
########################################################################################################################################################################################

SNV_activities=fread("SBS96_S10_NMF_Activities.txt")
colnames(SNV_activities)=c("Samples","SBS96B","SBS96A","SBS96F","SBS96C","SBS96H","SBS96I","SBS96E","SBS96J","SBS96G","SBS96D")
SNV_activities %>% dplyr::select(sort(names(.))) -> SNV_activities

SNV_activities=SNV_activities[1:21,]
SNV_activities=cbind(Cell=substr(SNV_activities$Samples,1,4),SNV_activities)
SNV_activities$Samples=sapply(str_split(SNV_activities$Samples,"_"), "[[", 2)
SNV_activities$Samples=mapvalues(SNV_activities$Samples,from=c("CA","CC","HE1A","HE2B","PR1B","PR2B","CB","HE2A","PR1A","PR2A","C","HE1B","HE1C","PR2C"),
                                to=c("Control_A","Control_C","Helium_1A","Helium_2B","Protons_1B","Protons_2B","Control_B","Helium_2A","Protons_1A","Protons_2A","Control_A","Helium_1B","Helium_1C","Protons_2C"))
SNV_activities=cbind(Clone=paste(SNV_activities$Cell,"  -  ",SNV_activities$Samples,sep=""),SNV_activities)

#get by row and by column normalised activities
SNV_activities_row_norm=SNV_activities
SNV_activities_col_norm=SNV_activities
SNV_activities_row_norm[, 4:13]=sweep(SNV_activities[, 4:13],MARGIN=1,FUN="/",STATS=rowSums(SNV_activities[, 4:13]))
SNV_activities_col_norm[, 4:13]=sweep(SNV_activities[, 4:13],MARGIN=2,FUN="/",STATS=colSums(SNV_activities[, 4:13]))
pivot_longer(SNV_activities_row_norm, cols=4:13, names_to = "Signature", values_to = "Mutations") -> SNV_activities_row_norm
pivot_longer(SNV_activities_col_norm, cols=4:13, names_to = "Signature", values_to = "Mutations") -> SNV_activities_col_norm

pivot_longer(SNV_activities, cols=4:13, names_to = "Signature", values_to = "Mutations") -> SNV_activities
SNV_activities$Cell <- factor(SNV_activities$Cell, levels = c("HAP1", "A549", "MCF7"))

###SNV_activities=SNV_activities[-which(SNV_activities$Clone=="MCF7  -  Helium_2B"),]  #when removing the HE2B clone
ggplot(SNV_activities, aes(fill=Signature, y=Mutations, x=Samples)) +
  theme_bw()+
  geom_bar(position="stack", stat="identity",colour="black")+
  scale_y_continuous(breaks = seq(0, 5000, by = 1000),limits=c(0,5000))+
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
  scale_fill_manual(values=c("#E69F00","#CC79A7","#0072B2","#009E73","#F0E442","#D55E00","#56B4E9","#882255","#999999","#000000")) 
ggsave("Fig3B.pdf",plot = last_plot(),device="pdf", width = 5000,height = 2000,units = "px")

########################################################################################################################################################################################
########################################################################################################################################################################################


########################################################################################################################################################################################   
##################################################################################   Fig.3C   ##########################################################################################   
########################################################################################################################################################################################

#Fig 3C
SNV_activities_row_norm$Mutations[SNV_activities_row_norm$Mutations == 0] <- NA
###SNV_activities_row_norm=SNV_activities_row_norm[-which(SNV_activities_row_norm$Clone=="MCF7  -  Helium_2B"),]  #when removing the HE2B clone

ggplot(SNV_activities_row_norm, aes(x = Signature, y = Clone, color = Signature)) +
  geom_count(aes(size = Mutations))+
  scale_color_manual(values=c("#E69F00","#CC79A7","#0072B2","#009E73","#F0E442","#D55E00","#56B4E9","#882255","#999999","#000000")) +
  scale_size(breaks = c(0.005,0.01,0.025,0.05,0.1,0.2,0.3,0.4,0.5),range = c(1.5,12))+
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_text(angle=270, hjust=0.5,size=12),
        axis.text.y = element_text(hjust=0,size=12),
        text = element_text(size = 15),
        panel.background=element_rect(fill="white"),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(3, "lines"),
        legend.key = element_rect(fill = "white"),
        #legend.position = "bottom",
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1))+
  geom_vline(xintercept = seq(1.5, 9.5, by = 1),color="grey40")+
  geom_hline(yintercept = seq(1.5, 20.5, by = 1),color="grey40")+
  ###geom_hline(yintercept = seq(1.5, 19.5, by = 1),color="grey40")+
  labs(size = "Attibuted mutation\nfraction per clone")+
  guides(colour = guide_legend(order = 1,override.aes = list(size=5)),
         size = guide_legend(order = 2))

ggsave("Fig3C.pdf",plot = last_plot(),device="pdf", width = 2300,height = 2700,units = "px")

########################################################################################################################################################################################
########################################################################################################################################################################################


########################################################################################################################################################################################   
##############################################################################   Fig.3C col-norm   #####################################################################################   
########################################################################################################################################################################################

#Fig 3C
SNV_activities_col_norm$Mutations[SNV_activities_col_norm$Mutations == 0] <- NA
###SNV_activities_col_norm=SNV_activities_col_norm[-which(SNV_activities_col_norm$Clone=="MCF7  -  Helium_2B"),]  #when removing the HE2B clone

ggplot(SNV_activities_col_norm, aes(x = Signature, y = Clone, color = Signature)) +
  geom_count(aes(size = Mutations))+
  scale_color_manual(values=c("#E69F00","#CC79A7","#0072B2","#009E73","#F0E442","#D55E00","#56B4E9","#882255","#999999","#000000")) +
  scale_size(breaks = c(0.001,0.005,0.01,0.025,0.05,0.1,0.15,0.2,0.25,0.3),range = c(0.5,10))+
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_text(angle=270, hjust=0.5,size=12),
        axis.text.y = element_text(hjust=0,size=12),
        text = element_text(size = 15),
        panel.background=element_rect(fill="white"),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(3, "lines"),
        legend.key = element_rect(fill = "white"),
        #legend.position = "bottom",
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1))+
  geom_vline(xintercept = seq(1.5, 9.5, by = 1),color="grey40")+
  geom_hline(yintercept = seq(1.5, 20.5, by = 1),color="grey40")+
  labs(size = "Assigned mutation\nfraction per signature")+
  guides(colour = guide_legend(order = 1,override.aes = list(size=5)),
         size = guide_legend(order = 2))

ggsave("Fig3C_by_col.pdf",plot = last_plot(),device="pdf", width = 2300,height = 2700,units = "px")

########################################################################################################################################################################################
########################################################################################################################################################################################


########################################################################################################################################################################################   
###############################################################################   SFig2C  ##############################################################################################   
########################################################################################################################################################################################

SNV_activities=fread("SBS96_S10_NMF_Activities.txt")
colnames(SNV_activities)=c("Samples","SBS96B","SBS96A","SBS96F","SBS96C","SBS96H","SBS96I","SBS96E","SBS96J","SBS96G","SBS96D")
SNV_activities %>% dplyr::select(sort(names(.))) -> SNV_activities

SNV_activities=SNV_activities[1:21,]
SNV_activities=cbind(Cell=substr(SNV_activities$Samples,1,4),SNV_activities)
SNV_activities$Samples=sapply(str_split(SNV_activities$Samples,"_"), "[[", 2)
SNV_activities$Samples=mapvalues(SNV_activities$Samples,from=c("CA","CC","HE1A","HE2B","PR1B","PR2B","CB","HE2A","PR1A","PR2A","C","HE1B","HE1C","PR2C"),
                                 to=c("Control_A","Control_C","Helium_1A","Helium_2B","Protons_1B","Protons_2B","Control_B","Helium_2A","Protons_1A","Protons_2A","Control_A","Helium_1B","Helium_1C","Protons_2C"))
SNV_activities=cbind(Clone=paste(SNV_activities$Cell,"  -  ",SNV_activities$Samples,sep=""),SNV_activities)

pivot_longer(SNV_activities, cols=4:13, names_to = "Signature", values_to = "Mutations") -> SNV_activities
SNV_activities$Cell <- factor(SNV_activities$Cell, levels = c("HAP1", "A549", "MCF7"))
SNV_activities$Treatment=sapply(str_split(SNV_activities$Samples,"_"), "[[", 1)
###SNV_activities=SNV_activities[-which(SNV_activities$Clone=="MCF7  -  Helium_2B"),]  #when removing the HE2B clone

ggplot(SNV_activities, aes(x = Treatment, y = Mutations, fill=Signature)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  facet_wrap(~Signature, nrow=5,scales = "free_y") +
  scale_fill_manual(values=c("#E69F00","#CC79A7","#0072B2","#009E73","#F0E442","#D55E00","#56B4E9","#882255","#999999","#000000")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 20),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(size=15),
        legend.position = "bottom",
        panel.spacing = unit(1, "lines"))+
  ylab("Number of attributed single base substitutions\n")
ggsave("SFig2C.pdf",plot = last_plot(),device="pdf",path="/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/plots/sig_plots",  width = 3000,height = 3000,units = "px")

########################################################################################################################################################################################
########################################################################################################################################################################################


########################################################################################################################################################################################   
##################################################################################   Fig.3D   ##########################################################################################  
########################################################################################################################################################################################

library(MutationalPatterns)
library(data.table)
library(dplyr)
library(tidyr)

COSMIC_SBS=fread("COSMIC_v3.3.1_SBS_GRCh38.txt") #From the COSMIC website
COSMIC_SBS=COSMIC_SBS[,-1]

SNVsignatures=as.data.frame(fread("SBS96_S10_Signatures.txt"))
colnames(SNVsignatures)=c("MutationType","SBS96B","SBS96A","SBS96F","SBS96C","SBS96H","SBS96I","SBS96E","SBS96J","SBS96G","SBS96D")
SNVsignatures %>% dplyr::select(sort(names(.))) -> SNVsignatures
SNVsignatures=SNVsignatures[,-1]

csm = cos_sim_matrix(COSMIC_SBS, SNVsignatures) 
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
ggsave("SNV_dendrogram.pdf",plot = last_plot(),device="pdf", width = 5450,height = 350,units = "px")

#create cosine similarity df
csm=as.data.frame(csm)
csm$COSMIC_sig=rownames(csm)
pivot_longer(csm, cols=1:10, names_to = "Signature", values_to = "Similarity") -> csm
csm$COSMIC_sig <- factor(csm$COSMIC_sig, levels = col_order)
#plot
ggplot(data = csm, aes(x = COSMIC_sig, y = Signature)) +
  geom_raster(aes(fill = Similarity))+
  scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Cosine \nsimilarity\n", limits = c(0, 1.000000001)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust=1, size = 15),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.ticks = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 20))+
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL)
ggsave("SNV_heatmap.pdf",plot = last_plot(),device="pdf", width = 6000,height = 1500,units = "px")
########################################################################################################################################################################################
########################################################################################################################################################################################


########################################################################################################################################################################################   
####################################################################################   SNV decomp.   ###################################################################################   
########################################################################################################################################################################################
library(gtools)

SBS_decomposition=data.frame(
  Signature=c("SBS96A","SBS96A","SBS96A","SBS96B","SBS96C","SBS96C","SBS96C","SBS96D","SBS96D","SBS96E","SBS96E","SBS96F","SBS96F","SBS96G","SBS96G","SBS96G","SBS96H","SBS96H","SBS96I","SBS96I","SBS96I","SBS96J","SBS96J"),
  Component=c("SBS45","SBS36","SBS5","SBS18","SBS46","SBS57","SBS5","SBS5","SBS21","SBS22","SBS5","SBS44","SBS1","SBS7b","SBS31","SBS7a","SBS4","SBS24","SBS40","SBS5","SBS1","SBS92","SBS26"),
  Contribution=c(64.26,28.32,7.42,100.00,54.96,34.36,10.68,69.78,30.22,88.62,11.22,94.44,5.56,38.28,38.16,23.56,60.06,39.94,78.84,19.84,1.32,52.76,47.24),
  stringsAsFactors = F)
SBS_decomposition$Component <- factor(SBS_decomposition$Component, levels = unique(mixedsort(SBS_decomposition$Component)))
colours=c("#F3C300","#875692","#F38400","#A1CAF1","#BE0032","#848482","#008856","#E68FAC","#0067A5","#F99379","#B3446C","#DCD300","#8DB600","#654522","#99ddff","#bbbbbb","#44bb99","#E25822")

ggplot(SBS_decomposition, aes(fill=Component, y=Contribution, x=Signature)) +
  geom_bar(position="stack", stat="identity",colour="black")+
  scale_fill_manual(values=colours)+
  scale_y_continuous(limits = c(0, 105), breaks = c(0,25,50,75,100))+
  geom_text(aes(label = ifelse(Contribution != 1.32, as.character(Component), "")), position = position_stack(vjust = 0.5),size=5.5) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 25),
        plot.title = element_text(hjust = 0.5,size=20,vjust=0),
        panel.background=element_rect(fill="white"),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))+
  ylab("Contribution (%)")+
  labs(fill="COSMIC\nsignatures")+
  ggtitle("Cosine similarity with reconstructed signature")+
  annotate("text", x=1, y=105, label= "0.968",size=6) +
  annotate("text", x = 2, y=105, label = "0.92",size=6)+
  annotate("text", x = 3, y=105, label = "0.854",size=6)+
  annotate("text", x = 4, y=105, label = "0.834",size=6)+
  annotate("text", x = 5, y=105, label = "0.949",size=6)+
  annotate("text", x = 6, y=105, label = "0.954",size=6)+
  annotate("text", x = 7, y=105, label = "0.978",size=6)+
  annotate("text", x = 8, y=105, label = "0.87",size=6)+
  annotate("text", x = 9, y=105, label = "0.854",size=6)+
  annotate("text", x = 10, y=105, label = "0.909",size=6)
ggsave("SNV_decomp.pdf",plot = last_plot(),device="pdf", width = 3300,height = 2800,units = "px")

########################################################################################################################################################################################
########################################################################################################################################################################################


########################################################################################################################################################################################   
####################################################################################   Subsampling   ###################################################################################   
########################################################################################################################################################################################

library(lsa)
library(ggcorrplot)
library(patchwork)
library(data.table)
library(dplyr)

SNVsignatures=as.data.frame(fread("SBS96_S10_Signatures.txt"))
SNV_mutation_classes=SNVsignatures$MutationType
colnames(SNVsignatures)=c("MutationType","SBS96B","SBS96A","SBS96F","SBS96C","SBS96H","SBS96I","SBS96E","SBS96J","SBS96G","SBS96D")
SNVsignatures %>% dplyr::select(sort(names(.))) -> SNVsignatures

subsignatures_1=as.data.frame(fread("SBS96_S10_Signatures_Sub1.txt"))
subsignatures_2=as.data.frame(fread("SBS96_S10_Signatures_Sub2.txt"))
subsignatures_3=as.data.frame(fread("SBS96_S10_Signatures_Sub3.txt"))

sigs=c("SBS96A","SBS96B","SBS96C","SBS96D","SBS96E","SBS96F","SBS96G","SBS96H","SBS96I","SBS96J")
sub_corr_df1=data.frame(Original_sig=c(),Sub_sig=c(),Cosine=c(),stringsAsFactors = F)
sub_corr_df2=data.frame(Original_sig=c(),Sub_sig=c(),Cosine=c(),stringsAsFactors = F)
sub_corr_df3=data.frame(Original_sig=c(),Sub_sig=c(),Cosine=c(),stringsAsFactors = F)

for (i in sigs){for (j in sigs){
  row1=data.frame(Original_sig=i,Sub_sig=j,Cosine=round(cosine(SNVsignatures[,i], subsignatures_1[,j]),2),stringsAsFactors = F)
  sub_corr_df1=rbind(sub_corr_df1,row1)
  row2=data.frame(Original_sig=i,Sub_sig=j,Cosine=round(cosine(SNVsignatures[,i], subsignatures_2[,j]),2),stringsAsFactors = F)
  sub_corr_df2=rbind(sub_corr_df2,row2)
  row3=data.frame(Original_sig=i,Sub_sig=j,Cosine=round(cosine(SNVsignatures[,i], subsignatures_3[,j]),2),stringsAsFactors = F)
  sub_corr_df3=rbind(sub_corr_df3,row3)}}

sub_corr_df1$Run=rep("Run1",100)
sub_corr_df2$Run=rep("Run2",100)
sub_corr_df3$Run=rep("Run3",100)
sub_corr_df=rbind(sub_corr_df1,sub_corr_df2,sub_corr_df3)

ggplot(data = sub_corr_df, aes(x=Sub_sig, y=Original_sig, fill=Cosine)) + 
  geom_tile()+
  facet_wrap(~Run,ncol=1)+
  geom_text(aes(Sub_sig, Original_sig, label = Cosine), color = "grey50", size = 6) +
  scale_fill_distiller(palette = "YlGnBu",limits=c(0,1),direction = 1,breaks=c(0,0.5,1))+
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust=0.5, size = 15),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=20),
        legend.title = element_text(size=20),
        legend.position = "bottom",
        panel.spacing = unit(2, "lines"),
        panel.background = element_rect(fill="grey85"),
        text = element_text(size = 25))+
  scale_x_discrete(position = "top")+
  labs(fill="Cosine")+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))+
  xlab("Subsampled SBS96 signatures (117 samples)\n")+
  ylab("Current SBS96 signatures (214 samples)\n") -> p1


SNV_activities=as.data.frame(fread("SBS96_S10_NMF_Activities.txt"))
colnames(SNV_activities)=c("Samples","SBS96B","SBS96A","SBS96F","SBS96C","SBS96H","SBS96I","SBS96E","SBS96J","SBS96G","SBS96D")
SNV_activities %>% dplyr::select(sort(names(.))) -> SNV_activities

#get subsampled activities, rename columns and reorder so we're using the original signature names
subactivities_1=as.data.frame(fread("SBS96_S10_NMF_Activities_Sub1.txt"))
subactivities_2=as.data.frame(fread("SBS96_S10_NMF_Activities_Sub2.txt"))
subactivities_3=as.data.frame(fread("SBS96_S10_NMF_Activities_Sub3.txt"))

SNV_activities=SNV_activities[,-1]
subactivities_1=subactivities_1[,-1]
subactivities_2=subactivities_2[,-1]
subactivities_3=subactivities_3[,-1]
sub_corract_df1=data.frame(Original_sig=c(),Sub_sig=c(),Cosine=c(),stringsAsFactors = F)
sub_corract_df2=data.frame(Original_sig=c(),Sub_sig=c(),Cosine=c(),stringsAsFactors = F)
sub_corract_df3=data.frame(Original_sig=c(),Sub_sig=c(),Cosine=c(),stringsAsFactors = F)

for (i in 1:10){for (j in 1:10){
  row1=data.frame(Original_sig=sigs[i],Sub_sig=sigs[j],Cosine=round(cosine(SNV_activities[1:21,i], subactivities_1[1:21,j]),2),stringsAsFactors = F)
  sub_corract_df1=rbind(sub_corract_df1,row1)
  row2=data.frame(Original_sig=sigs[i],Sub_sig=sigs[j],Cosine=round(cosine(SNV_activities[1:21,i], subactivities_2[1:21,j]),2),stringsAsFactors = F)
  sub_corract_df2=rbind(sub_corract_df2,row2)
  row3=data.frame(Original_sig=sigs[i],Sub_sig=sigs[j],Cosine=round(cosine(SNV_activities[1:21,i], subactivities_3[1:21,j]),2),stringsAsFactors = F)
  sub_corract_df3=rbind(sub_corract_df3,row3)}}
sub_corract_df1$Run=rep("Run1",100)
sub_corract_df2$Run=rep("Run2",100)
sub_corract_df3$Run=rep("Run3",100)
sub_corract_df=rbind(sub_corract_df1,sub_corract_df2,sub_corract_df3)

ggplot(data = sub_corract_df, aes(x=Sub_sig, y=Original_sig, fill=Cosine)) + 
  geom_tile()+
  facet_wrap(~Run,ncol=1)+
  geom_text(aes(Sub_sig, Original_sig, label = Cosine), color = "grey50", size = 6) +
  scale_fill_distiller(palette = "YlGnBu",limits=c(0,1),direction = 1,breaks=c(0,0.5,1))+
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust=0.5, size = 15),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=20),
        legend.title = element_text(size=20),
        legend.position = "bottom",
        panel.background = element_rect(fill="grey85"),
        panel.spacing = unit(2, "lines"),
        text = element_text(size = 25))+
  scale_x_discrete(position = "top")+
  xlab("Subsampled SBS96 signature exposures\nover all 21 clones (117 samples)\n")+
  ylab("Current SBS96 signature exposures\nover all 21 clones (214 samples)\n")+
  labs(fill="Cosine")+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) -> p2

p1 + plot_spacer() + p2 + plot_layout(widths = c(10, 1 ,10))
ggsave("SFig5A.pdf",plot = last_plot(),device="pdf", width = 6200,height = 6200,units = "px")

sub_corr_df_barplot=data.frame(Original_sig=c(),Sub_sig=c(),Run=c(),Max=c(),Class=c(),stringsAsFactors = F)
for (sig in sigs){for (run in c("Run1","Run2","Run3")){
         df=sub_corr_df[which(sub_corr_df$Original_sig==sig & sub_corr_df$Run==run),]
         row=df[which.max(df$Cosine),]
         entry=data.frame(Original_sig=row$Original_sig,Sub_sig=row$Sub_sig,Run=row$Run,Max=row$Cosine,Class="Signature spectra",stringsAsFactors = F)
         sub_corr_df_barplot=rbind(sub_corr_df_barplot,entry)}}

#we take the exposure of the sub signature with the highest similarity to the original from sub_corr_df_barplot (instead of just taking the highest value for each original signature)
sub_corract_df_barplot=data.frame(Original_sig=c(),Sub_sig=c(),Run=c(),Max=c(),Class=c(),stringsAsFactors = F)
for (sig in sigs){for (run in c("Run1","Run2","Run3")){
        sub=sub_corr_df_barplot[which(sub_corr_df_barplot$Original_sig==sig & sub_corr_df_barplot$Run==run),"Sub_sig"]
        row=sub_corract_df[which(sub_corract_df$Original_sig==sig & sub_corract_df$Run==run & sub_corract_df$Sub_sig==sub),]
        entry=data.frame(Original_sig=row$Original_sig,Sub_sig=row$Sub_sig,Run=row$Run,Max=row$Cosine,Class="Exposures",stringsAsFactors = F)
        sub_corract_df_barplot=rbind(sub_corract_df_barplot,entry)}}

sub_corr_barplot=rbind(sub_corr_df_barplot,sub_corract_df_barplot)
sub_corr_barplot$Class=factor(sub_corr_barplot$Class,levels=c("Signature spectra","Exposures"))

ggplot(data = sub_corr_barplot, aes(x=Original_sig, y=Max,fill=Original_sig)) + geom_boxplot() + geom_jitter(width=0.2,size=5,shape=18,color="grey25") + facet_wrap(~ Class)+
  scale_fill_manual(values=c("#E69F00","#CC79A7","#0072B2","#009E73","#F0E442","#D55E00","#56B4E9","#882255","#999999","#000000")) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0, size = 20),
        axis.text.y = element_text(hjust = 0.5, vjust = 0, size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=25),
        legend.position = "none",
        panel.spacing = unit(2, "lines"),
        text = element_text(size = 25))+
  xlab("\nSignatures")+
  ylab("Cosine similarity\n")
ggsave("SFig5B.pdf",plot = last_plot(),device="pdf", width = 5000,height = 3000,units = "px")
########################################################################################################################################################################################
########################################################################################################################################################################################

########################################################################################################################################################################################   
#####################################################################################   SNV v ID   #####################################################################################  
########################################################################################################################################################################################

library(data.table)
library(dplyr)
library(ggplot2)
SNV_activities=fread("SBS96_S10_NMF_Activities.txt")
colnames(SNV_activities)=c("Samples","SBS96B","SBS96A","SBS96F","SBS96C","SBS96H","SBS96I","SBS96E","SBS96J","SBS96G","SBS96D")
SNV_activities %>% dplyr::select(sort(names(.))) -> SNV_activities
SNV_activities_row_norm=as.data.frame(SNV_activities)
SNV_activities_row_norm[, 2:11]=sweep(SNV_activities[, 2:11],MARGIN=1,FUN="/",STATS=rowSums(SNV_activities[, 2:11]))
SNV_activities_row_norm=SNV_activities_row_norm[1:21,-1]

ID_activities=fread("ID83_S4_NMF_Activities.txt")
colnames(ID_activities)=c("Samples","ID83A","ID83C","ID83B","ID83D") #change names to match the original set
ID_activities=cbind(ID_activities[,1],ID_activities[,2],ID_activities[,4],ID_activities[,3],ID_activities[,5]) #reorder
ID_activities_row_norm=as.data.frame(ID_activities)
ID_activities_row_norm[, 2:5]=sweep(ID_activities[, 2:5],MARGIN=1,FUN="/",STATS=rowSums(ID_activities[, 2:5])) #normalise by row to get % contribution per sample
ID_activities_row_norm=ID_activities_row_norm[1:21,-1]

SNV_ID_corr=data.frame(SNV_sig=c(),ID_sig=c(),corr=c(),p=c(),significance=c())
for (i in 1:10){for (j in 1:4){
       x=SNV_activities_row_norm[,i]
       y=ID_activities_row_norm[,j]
       z=cor.test(x, y, method = c("pearson")) 
       if(round(unname(z$p.value),4)<0.001){sig="***"}
       if(round(unname(z$p.value),4)>=0.001 && round(unname(z$p.value),4)<0.01){sig="**"}
       if(round(unname(z$p.value),4)>=0.01 && round(unname(z$p.value),4)<0.05){sig="*"}
       if(round(unname(z$p.value),4)>=0.05){sig=""}
       df=data.frame(SNV_sig=colnames(SNV_activities_row_norm)[i],ID_sig=colnames(ID_activities_row_norm)[j], corr=round(unname(z$estimate),2),p=round(unname(z$p.value),4),significance=sig)
       SNV_ID_corr=rbind(SNV_ID_corr,df)}}

ggplot(data = SNV_ID_corr, aes(x=SNV_sig, y=ID_sig, fill=corr)) + 
  geom_raster(aes(x=SNV_sig, y=ID_sig, fill=corr)) +
  geom_text(aes(SNV_sig, ID_sig, label = paste(corr,significance,sep="")), color = "grey20", size = 5) +
  scale_fill_gradient2( low = "red", high = "blue",limits=c(-1,1),breaks=c(-1,-0.75,-0.50,-0.25,0,0.25,0.50,0.75,1)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust=1, size = 15),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        legend.text = element_text(angle = 90,hjust=1),
        axis.title = element_blank(),
        panel.background = element_rect(fill="grey95"),
        text = element_text(size = 20))+
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,fill="Pearson correlation coefficient")+
  guides(fill = guide_colourbar(raster=T,title.position="top", title.hjust = 0.5,nbin=100))
ggsave("SFig6.pdf",plot = last_plot(),device="pdf", width = 3000,height = 1800,units = "px")

  
