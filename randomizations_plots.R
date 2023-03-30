library(ggplot2)

p_values_SNVs=read.table("all_comparisons_SNVs_sigs.txt",header=T,sep=" ")
all_signatures_SNVs=read.table("all_signatures_SNVs.txt",header=T,sep=" ")
all_segments_SNVs=read.table("all_segments_SNVs.txt",header=T,sep=" ")

p_values_IDs=read.table("all_comparisons_IDs_sigs.txt",header=T,sep=" ")
all_signatures_IDs=read.table("all_signatures_IDs.txt",header=T,sep=" ")
all_segments_IDs=read.table("all_segments_IDs.txt",header=T,sep=" ")

all_signatures_SNVs$Samples=factor(all_signatures_SNVs$Samples,unique(p_values_SNVs$Sample))
colnames(all_segments_SNVs)=c("x","xend","y","yend","Signature")

all_signatures_IDs$Samples=factor(all_signatures_IDs$Samples,unique(p_values_IDs$Sample))
colnames(all_segments_IDs)=c("x","xend","y","yend","Signature")

p_values=p_values_SNVs
all_signatures=all_signatures_SNVs
all_segments=all_segments_SNVs

p_values=p_values_IDs
all_signatures=all_signatures_IDs
all_segments=all_segments_IDs

p_values[which(p_values$P_value<=0.1),] -> significant_p_values
significant_p_values$FDR=round(significant_p_values$P_adjusted,3)
significant_p_values$Stars=NA
significant_p_values$x=NA
significant_p_values$y=NA

for (i in 1:nrow(significant_p_values)){
  significant_p_values$x[i]=which(unique(p_values$Sample) %in% significant_p_values$Sample[i])
  if(significant_p_values$H0[i]=="Higher"){significant_p_values$y[i]=1.35*min(all_signatures[which(all_signatures$Signature==significant_p_values$Signature[i] & all_signatures$Samples==significant_p_values$Sample[i] & all_signatures$Iteration=="Random"),"Mean_diff"])
  }else{significant_p_values$y[i]=1.35*max(all_signatures[which(all_signatures$Signature==significant_p_values$Signature[i] & all_signatures$Samples==significant_p_values$Sample[i] & all_signatures$Iteration=="Random"),"Mean_diff"])}
  if(significant_p_values$P_value[i]<=0.01){sig="***"}
  if(significant_p_values$P_value[i]>0.01 && significant_p_values$P_value[i]<=0.05){sig="**"}
  if(significant_p_values$P_value[i]>0.05 && significant_p_values$P_value[i]<=0.1){sig="*"}
  significant_p_values$Stars[i]=sig}

ggplot(all_signatures[all_signatures$Iteration=="Random",], aes(x=Samples, y=Mean_diff)) +
  geom_violin()+
  xlab("\n Pairwise comparisons")+
  ylab("Pairwise difference in group real and randomized exposure means\n")+
  facet_wrap(~Signature,ncol=2,scales="free_y")+
  geom_segment(data = all_segments, aes(x = x, y = y, xend = xend, yend = yend),size=1,color="red")+
  geom_text(data = significant_p_values, aes(x = x, y = y, label = FDR),size=4,fontface="bold")+
  geom_text(data = significant_p_values, aes(x = x, y = rep(0,times=nrow(significant_p_values)), label = Stars),size=6,fontface="bold")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 15,margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 15),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.switch.pad.wrap = unit(4, "cm"),
        axis.ticks = element_blank(),
        axis.title = element_text(size=20),
        legend.title = element_text(size=20),
        text = element_text(size = 20)) -> plot
ggsave("SFig4A.pdf",plot = plot,device="pdf", width = 8500,height = 5000,units = "px")
ggsave("SFig4B.pdf",plot = plot,device="pdf", width = 8500,height = 2500,units = "px")

################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################

p_values_TMBs=read.table("all_comparisons_TMBs_sigs.txt",header=T,sep=" ")
all_signatures_TMBs=read.table("all_signatures_TMBs.txt",header=T,sep=" ")
all_segments_TMBs=read.table("all_segments_TMBs.txt",header=T,sep=" ")

colnames(all_segments_TMBs)=c("x","xend","y","yend","Signature")
all_signatures_TMBs$Samples=factor(all_signatures_TMBs$Samples,unique(p_values_TMBs$Sample))

p_values=p_values_TMBs
all_signatures=all_signatures_TMBs
all_segments=all_segments_TMBs
all_signatures$Signature=factor(all_signatures$Signature,levels=c("SNVs","IDs","IDs_del_ins_ratio","Del1_Dels_ratio_no_reps","Ins1_Inss_ratio_no_reps","MH_Dels_Del5_ratio","MH1_Del_MH_Dels_ratio","DEL","DUP","INS","INV","Del_InsDup_ratio"))
                                  
p_values[which(p_values$P_value<=0.1),] -> significant_p_values
significant_p_values$FDR=round(significant_p_values$P_adjusted,3)
significant_p_values$Stars=NA
significant_p_values$x=NA
significant_p_values$y=NA

for (i in 1:nrow(significant_p_values)){
  significant_p_values$x[i]=which(unique(p_values$Sample) %in% significant_p_values$Sample[i])
  if(significant_p_values$H0[i]=="Higher"){significant_p_values$y[i]=1.35*min(all_signatures[which(all_signatures$Signature==significant_p_values$Signature[i] & all_signatures$Samples==significant_p_values$Sample[i] & all_signatures$Iteration=="Random"),"Mean_diff"])
  }else{significant_p_values$y[i]=1.35*max(all_signatures[which(all_signatures$Signature==significant_p_values$Signature[i] & all_signatures$Samples==significant_p_values$Sample[i] & all_signatures$Iteration=="Random"),"Mean_diff"])}
  if(significant_p_values$P_value[i]<=0.01){sig="***"}
  if(significant_p_values$P_value[i]>0.01 && significant_p_values$P_value[i]<=0.05){sig="**"}
  if(significant_p_values$P_value[i]>0.05 && significant_p_values$P_value[i]<=0.1){sig="*"}
  significant_p_values$Stars[i]=sig}

ggplot(all_signatures[all_signatures$Iteration=="Random",], aes(x=Samples, y=Mean_diff)) +
  geom_violin()+
  xlab("\n Pairwise comparisons")+
  ylab("Pairwise difference in group real and randomized mutation count means \n")+
  facet_wrap(~Signature,ncol=2,scales="free_y")+
  geom_segment(data = all_segments, aes(x = x, y = y, xend = xend, yend = yend),size=1,color="red")+
  geom_text(data = significant_p_values, aes(x = x, y = y, label = FDR),size=4,fontface="bold")+
  geom_text(data = significant_p_values, aes(x = x, y = rep(0,times=nrow(significant_p_values)), label = Stars),size=6,fontface="bold")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 15,margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 15),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.switch.pad.wrap = unit(4, "cm"),
        axis.ticks = element_blank(),
        axis.title = element_text(size=20),
        legend.title = element_text(size=20),
        text = element_text(size = 20)) -> plot
ggsave("SFig1.pdf",plot = plot,device="pdf", width = 8500,height = 8500,units = "px")


