#This script outputs plots relating to the SNV/ID content of our samples. It doesn't contain any signature or cluster related info.
library(stringr)
library(plyr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggridges)
library(data.table)

########################################################################################################################################################################################   
######################################################################################   Fig.2A   ######################################################################################   
########################################################################################################################################################################################

SNVs=read.table("Our_Zou_Kucab_SNVs.txt",header = T)[,1:22]
#make counts dataframe
Counts=data.frame(Cell=rep(sapply(str_split(colnames(SNVs[,2:22]),"_"), "[[", 1),3), Samples=rep(sapply(str_split(colnames(SNVs[,2:22]),"_"), "[[", 2),3), Variant=c(rep("SNV",21),rep("INS",21),rep("DEL",21)), Count=rep(NA,63))
Counts$Samples=mapvalues(Counts$Samples,from=c("CA","CC","HE1A","HE2B","PR1B","PR2B","CB","HE2A","PR1A","PR2A","C","HE1B","HE1C","PR2C"),
                                    to=c("Control_A","Control_C","Helium_1A","Helium_2B","Protons_1B","Protons_2B","Control_B","Helium_2A","Protons_1A","Protons_2A","Control_A","Helium_1B","Helium_1C","Protons_2C"))              
Counts=cbind(Clone=paste(Counts$Cell,"  -  ",Counts$Samples,sep=""),Counts)
Counts$Count[1:21]=unname(colSums(SNVs[,2:22])) #add number of SNVs in our dataframe
IDs=read.table("Our_Zou_Kucab_indels.txt",header = T)[,1:22]
IDs$MutationType=sapply(str_split(IDs$MutationType,":"), "[[", 2)
pivot_longer(IDs, cols=2:22, names_to = "Sample", values_to = "Mutations") %>% group_by(Sample, MutationType) %>% summarise(across(Mutations, sum)) %>% filter(MutationType == "Ins") %>% pull(Mutations) -> Counts$Count[22:42]
pivot_longer(IDs, cols=2:22, names_to = "Sample", values_to = "Mutations") %>% group_by(Sample, MutationType) %>% summarise(across(Mutations, sum)) %>% filter(MutationType == "Del") %>% pull(Mutations) -> Counts$Count[43:63]
Counts$Cell <- factor(Counts$Cell, levels = c("HAP1", "A549", "MCF7"))
Counts$Variant <- factor(Counts$Variant, levels = c("SNV","INS","DEL"))

###Counts=Counts[-which(Counts$Clone=="MCF7  -  Helium_2B"),]  #when removing the HE2B clone
ggplot(Counts, aes(fill=Variant, y=Count, x=Samples)) + facet_grid(Variant~Cell, scale="free") +
  geom_bar(position="stack", stat="identity",colour="black")+
  geom_text(aes(label=Count), position = position_stack(vjust = 0.5),size=5.5)+
  scale_fill_discrete(name = "Variant type")+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        text = element_text(size = 25),
        panel.background=element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"),
        legend.position="bottom",
        strip.text.x = element_text(size = 25,face = "bold"))+
  ylab("Number of mutations\n")+
  xlab("\nClone")
ggsave("Fig2A.pdf",plot = last_plot(),device="pdf", width = 6000,height = 4500,units = "px")

########################################################################################################################################################################################
########################################################################################################################################################################################


########################################################################################################################################################################################   
######################################################################################   Fig.2B   ######################################################################################   
########################################################################################################################################################################################

SNV_VAFs=data.frame(Samples=c(),VAF=c())
for (sample in list.files("Filtered_SNVs.zip")){
     file=read.table(paste0("Filtered_SNVs.zip",sample),header=T)
     SNV_VAFs=rbind(SNV_VAFs,data.frame(Samples=rep(sample,nrow(file)),VAF=file$ALTVAF))}

SNV_VAFs$Cell=sapply(str_split(SNV_VAFs$Samples,"_"), "[[", 1)
SNV_VAFs$Samples=gsub(".txt","",sapply(str_split(SNV_VAFs$Samples,"_"), "[[", 2))
SNV_VAFs$Samples=mapvalues(SNV_VAFs$Samples,from=c("CA","CC","HE1A","HE2B","PR1B","PR2B","CB","HE2A","PR1A","PR2A","C","HE1B","HE1C","PR2C"),
                         to=c("Control_A","Control_C","Helium_1A","Helium_2B","Protons_1B","Protons_2B","Control_B","Helium_2A","Protons_1A","Protons_2A","Control_A","Helium_1B","Helium_1C","Protons_2C")) 
SNV_VAFs=cbind(Clone=paste(SNV_VAFs$Cell,"  -  ",SNV_VAFs$Samples,sep=""),SNV_VAFs)
SNV_VAFs$Treatment=sapply(str_split(SNV_VAFs$Samples,"_"), "[[", 1); SNV_VAFs$Treatment=factor(SNV_VAFs$Treatment,levels=c("Protons","Helium" ,"Control"))
SNV_VAFs$Cell <- factor(SNV_VAFs$Cell, levels = c("HAP1", "A549", "MCF7"))

###SNV_VAFs=SNV_VAFs[-which(SNV_VAFs$Clone=="MCF7  -  Helium_2B"),]  #when removing the HE2B clone
ggplot(SNV_VAFs, aes(x = VAF, y = Samples, fill=Treatment)) +
  geom_density_ridges(scale = 1)+
  facet_wrap(~Cell, nrow=1,scales = "free_y")+
  scale_x_continuous(breaks = seq(0, 1, by = 0.25))+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        strip.background = element_blank(),
        text = element_text(size = 25),
        panel.background=element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"),
        legend.position="none",
        strip.text.x = element_text(size = 25,face = "bold",hjust=0.5,vjust=1,margin = margin(0,0,2,0, "cm")))+
  xlab("\nSNV Variant allele fraction")+
  ylab("Clone\n")
ggsave("Fig2B.pdf",plot = last_plot(),device="pdf", width = 6500,height = 4000,units = "px")

########################################################################################################################################################################################
########################################################################################################################################################################################


########################################################################################################################################################################################   
######################################################################################   Fig.2C   ######################################################################################   
########################################################################################################################################################################################

ID_VAFs=data.frame(Samples=c(),VAF=c())
for (sample in list.files("Filtered_INDELs.zip")){
  file=read.table(paste0("Filtered_INDELs.zip",sample),header=T)
  ID_VAFs=rbind(ID_VAFs,data.frame(Samples=rep(sample,nrow(file)),VAF=file$ALTVAF))}

ID_VAFs$Cell=sapply(str_split(ID_VAFs$Samples,"_"), "[[", 1)
ID_VAFs$Samples=gsub(".txt","",sapply(str_split(ID_VAFs$Samples,"_"), "[[", 2))
ID_VAFs$Samples=mapvalues(ID_VAFs$Samples,from=c("CA","CC","HE1A","HE2B","PR1B","PR2B","CB","HE2A","PR1A","PR2A","C","HE1B","HE1C","PR2C"),
                           to=c("Control_A","Control_C","Helium_1A","Helium_2B","Protons_1B","Protons_2B","Control_B","Helium_2A","Protons_1A","Protons_2A","Control_A","Helium_1B","Helium_1C","Protons_2C")) 
ID_VAFs=cbind(Clone=paste(ID_VAFs$Cell,"  -  ",ID_VAFs$Samples,sep=""),ID_VAFs)
ID_VAFs$Treatment=sapply(str_split(ID_VAFs$Samples,"_"), "[[", 1); ID_VAFs$Treatment=factor(ID_VAFs$Treatment,levels=c("Protons","Helium" ,"Control"))
ID_VAFs$Cell <- factor(ID_VAFs$Cell, levels = c("HAP1", "A549", "MCF7"))

###ID_VAFs=ID_VAFs[-which(ID_VAFs$Clone=="MCF7  -  Helium_2B"),]  #when removing the HE2B clone
ggplot(ID_VAFs, aes(x = VAF, y = Samples, fill=Treatment)) +
  geom_density_ridges(scale = 1)+
  facet_wrap(~Cell, nrow=1,scales = "free_y")+
  scale_x_continuous(breaks = seq(0, 1, by = 0.25))+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        strip.background = element_blank(),
        text = element_text(size = 25),
        panel.background=element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"),
        legend.position="none",
        strip.text.x = element_text(size = 25,face = "bold",hjust=0.5,vjust=1,margin = margin(0,0,2,0, "cm")))+
  xlab("\nID Variant allele fraction")+
  ylab("Clone\n")
ggsave("Fig2C.pdf",plot = last_plot(),device="pdf", width = 6500,height = 4000,units = "px")

########################################################################################################################################################################################
########################################################################################################################################################################################


########################################################################################################################################################################################   
######################################################################################   Fig.3A   ######################################################################################   
########################################################################################################################################################################################

pivot_longer(SNVs, cols=2:22, names_to = "Clone", values_to = "Mutations") %>% mutate(SBS6=substr(MutationType,3,5)) %>% 
group_by(Clone, SBS6) %>% summarise(across(Mutations, sum)) %>% mutate(Cell=substr(Clone,1,4))-> SNVs_SBS6
SNVs_SBS6$Samples=sapply(str_split(SNVs_SBS6$Clone,"_"), "[[", 2)
SNVs_SBS6$Samples=mapvalues(SNVs_SBS6$Samples,from=c("CA","CC","HE1A","HE2B","PR1B","PR2B","CB","HE2A","PR1A","PR2A","C","HE1B","HE1C","PR2C"),
                          to=c("Control_A","Control_C","Helium_1A","Helium_2B","Protons_1B","Protons_2B","Control_B","Helium_2A","Protons_1A","Protons_2A","Control_A","Helium_1B","Helium_1C","Protons_2C")) 
SNVs_SBS6$Clone=paste(SNVs_SBS6$Cell,"  -  ",SNVs_SBS6$Samples,sep="")
SNVs_SBS6$SBS6 <- factor(SNVs_SBS6$SBS6, levels = c("T>G","T>C","T>A","C>T","C>G","C>A"))
SNVs_SBS6$Cell <- factor(SNVs_SBS6$Cell, levels = c("HAP1", "A549", "MCF7"))
###SNVs_SBS6=SNVs_SBS6[-which(SNVs_SBS6$Clone=="MCF7  -  Helium_2B"),]  #when removing the HE2B clone
###data_vline <- data.frame(Cell = c(rep("A549",5),rep("HAP1",5),rep("MCF7",7)), vline = c(seq(from=1.5, to=5.5, by = 1),seq(from=1.5, to=5.5, by = 1),seq(from=1.5, to=7.5, by = 1)))
data_vline <- data.frame(Cell = c(rep("A549",5),rep("HAP1",5),rep("MCF7",8)), vline = c(seq(from=1.5, to=5.5, by = 1),seq(from=1.5, to=5.5, by = 1),seq(from=1.5, to=8.5, by = 1)))
data_hline <- data.frame(Cell = c(rep("A549",5),rep("HAP1",5),rep("MCF7",5)),  hline = c(seq(from=1.5, to=5.5, by = 1),seq(from=1.5, to=5.5, by = 1),seq(from=1.5, to=5.5, by = 1)))

ggplot(SNVs_SBS6, aes(x = Samples, y = SBS6)) +
  geom_count(aes(size = Mutations,color=SBS6))+
  scale_size(breaks = c(55,75,150,300,600,900,1200,1500,1800),range = c(2,24)) +
  ###scale_size(breaks = c(90,150,300,600,900,1200,1500,1800),range = c(2,24)) +
  facet_grid(~Cell,space="free_x",scales = "free_x")+
  geom_vline(data = data_vline,aes(xintercept = vline),size=1)+
  geom_hline(data = data_hline,aes(yintercept = hline),size=1)+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        strip.background = element_blank(),
        text = element_text(size = 25),
        panel.background=element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(3, "lines"),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        strip.text.x = element_text(size = 25,face = "bold",hjust=0.5,vjust=1,margin = margin(0,0,1,0, "cm")))+
  scale_colour_manual(values = c("#E7C9C6FF","#ABCD72FF","#CBCACBFF","#D33C32FF","#050708FF","#5ABCEBFF"))+
  guides(size = guide_legend(nrow = 1),color="none")+
  xlab("\nClone")+
  ylab("Mutation type\n")

ggsave("Fig3A.pdf",plot = last_plot(),device="pdf", width = 6500,height = 3000,units = "px")
write.table(SNVs_SBS6, file = "all_SNVs_SBS6.txt", sep = "\t",col.names = TRUE,quote=F,row.names = F)

########################################################################################################################################################################################
########################################################################################################################################################################################


########################################################################################################################################################################################   
####################################################################################   Fig.3A norm #####################################################################################   
########################################################################################################################################################################################

SNVs_SBS6_norm=SNVs_SBS6
setDT(SNVs_SBS6_norm)[, Mutations_norm := Mutations / sum(Mutations), by = Clone]
###SNVs_SBS6_norm=SNVs_SBS6_norm[-which(SNVs_SBS6_norm$Clone=="MCF7  -  Helium_2B"),]  #when removing the HE2B clone
###data_vline <- data.frame(Cell = c(rep("A549",5),rep("HAP1",5),rep("MCF7",7)), vline = c(seq(from=1.5, to=5.5, by = 1),seq(from=1.5, to=5.5, by = 1),seq(from=1.5, to=7.5, by = 1)))

ggplot(SNVs_SBS6_norm, aes(x = Samples, y = SBS6)) +
  geom_count(aes(size = Mutations_norm,color=SBS6))+
  scale_size(breaks = c(0.025,0.05,0.1,0.2,0.3,0.4,0.5),range = c(2,24)) +
  facet_grid(~Cell,space="free_x",scales = "free_x")+
  geom_vline(data = data_vline,aes(xintercept = vline),size=1)+
  geom_hline(data = data_hline,aes(yintercept = hline),size=1)+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        strip.background = element_blank(),
        text = element_text(size = 25),
        panel.background=element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(3, "lines"),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        strip.text.x = element_text(size = 25,face = "bold",hjust=0.5,vjust=1,margin = margin(0,0,1,0, "cm")))+
  scale_colour_manual(values = c("#E7C9C6FF","#ABCD72FF","#CBCACBFF","#D33C32FF","#050708FF","#5ABCEBFF"))+
  guides(size = guide_legend(nrow = 1),color="none")+
  xlab("\nClone")+
  ylab("Mutation type\n")+
  guides(size = guide_legend(title = "Fraction",nrow=1))
ggsave("SFig2B.pdf",plot = last_plot(),device="pdf", width = 6500,height = 3000,units = "px")

########################################################################################################################################################################################
########################################################################################################################################################################################

########################################################################################################################################################################################   
######################################################################################   Fig.4C   ######################################################################################   
########################################################################################################################################################################################

IDs=read.table("Our_Zou_Kucab_indels.txt",header = T)[,1:22]
IDs$MutationType=substr(IDs$MutationType,3,5)
pivot_longer(IDs, cols=2:22, names_to = "Clone", values_to = "Mutations") %>% mutate(Cell=substr(Clone,1,4)) -> IDs_ratio
IDs_ratio$Samples=sapply(str_split(IDs_ratio$Clone,"_"), "[[", 2)
IDs_ratio$Samples=mapvalues(IDs_ratio$Samples,from=c("CA","CC","HE1A","HE2B","PR1B","PR2B","CB","HE2A","PR1A","PR2A","C","HE1B","HE1C","PR2C"),
                            to=c("Control_A","Control_C","Helium_1A","Helium_2B","Protons_1B","Protons_2B","Control_B","Helium_2A","Protons_1A","Protons_2A","Control_A","Helium_1B","Helium_1C","Protons_2C"))
aggregate(Mutations~Clone+MutationType+Cell+Samples, data=IDs_ratio, FUN=sum) -> IDs_ratio
Indel_ratios=data.frame(Cell=c(), Samples=c(), Treatment=c(),Ratio=c(), stringsAsFactors = F)
for (cell in c("A549","HAP1","MCF7")){
  df=IDs_ratio[which(IDs_ratio$Cell==cell),]
  samples=unique(df$Samples)
  for (sample in samples){
    ins=df[which(df$Samples==sample & df$MutationType=="Ins"),"Mutations"]
    del=df[which(df$Samples==sample & df$MutationType=="Del"),"Mutations"]
    row=data.frame(Cell=cell,Samples=sample,Treatment=sapply(str_split(sample,"_"), "[[", 1),Ratio=del/ins)
    Indel_ratios=rbind(Indel_ratios,row) }}
Indel_ratios$Cell <- factor(Indel_ratios$Cell, levels = c("HAP1", "A549", "MCF7"))
###Indel_ratios=Indel_ratios[-which(Indel_ratios$Samples=="Helium_2B" & Indel_ratios$Cell=="MCF7"),]

ggplot(Indel_ratios, aes(x=Treatment, y=Ratio)) +
  geom_boxplot(fill=NA, colour="grey20",outlier.shape = NA)+
  geom_jitter(size=3, shape=21, colour="grey20",fill="grey20",width = 0.2)+
  facet_grid(~Cell)+
  theme(axis.text.x = element_text(angle = 90, hjust=0),
        text = element_text(size = 25),
        strip.background = element_rect(fill="white"),
        strip.text.x = element_text(size=30,face="bold"),
        panel.background=element_rect(fill="white"),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"))+
  ylab("Deletions / Insertions ratio\n")+
  xlab("\nTreatment")

ggsave("Fig4C.pdf",plot = last_plot(),device="pdf", width = 5500,height = 2500,units = "px")

########################################################################################################################################################################################
########################################################################################################################################################################################