library(data.table)
library(plyr)
library(magrittr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(qvalue)
library(stringr)
options(scipen=999) 

iter=100000
iter1=iter+1
iter2=iter+2
norm=FALSE

treatments=c("Control","Helium","Protons")
cells=c("A549","HAP1","MCF7")
SNV_sigs=c("SBS96A","SBS96B","SBS96C","SBS96D","SBS96E","SBS96F","SBS96G","SBS96H","SBS96I","SBS96J")
ID_sigs=c("ID83A","ID83B","ID83C","ID83D")

SNV_activities=as.data.frame(fread("SBS96_S10_NMF_Activities.txt")) #SNV signature activities file (from SigProfiler)
colnames(SNV_activities)=c("Samples","SBS96B","SBS96A","SBS96F","SBS96C","SBS96H","SBS96I","SBS96E","SBS96J","SBS96G","SBS96D") #reorder sigs
SNV_activities %>% dplyr::select(sort(names(.))) -> SNV_activities
SNV_activities=SNV_activities[1:21,]
SNV_activities$Samples=mapvalues(SNV_activities$Samples, from=c("A549_CA", "A549_CC", "A549_HE1A","A549_HE2B","A549_PR1B","A549_PR2B",
                                                                "HAP1_CA", "HAP1_CB", "HAP1_HE1A","HAP1_HE2A","HAP1_PR1A","HAP1_PR2A",
                                                                "MCF7_C", "MCF7_HE1B", "MCF7_HE1C","MCF7_HE2A","MCF7_HE2B","MCF7_PR1A","MCF7_PR1B","MCF7_PR2A","MCF7_PR2C"),
                                 to=c("A549_Control", "A549_Control", "A549_Helium","A549_Helium","A549_Protons","A549_Protons",
                                      "HAP1_Control", "HAP1_Control", "HAP1_Helium","HAP1_Helium","HAP1_Protons","HAP1_Protons",
                                      "MCF7_Control", "MCF7_Helium","MCF7_Helium","MCF7_Helium","MCF7_Helium","MCF7_Protons","MCF7_Protons","MCF7_Protons","MCF7_Protons"))
if (norm==TRUE){SNV_activities[, 2:11]=sweep(SNV_activities[, 2:11],MARGIN=1,FUN="/",STATS=rowSums(SNV_activities[, 2:11]))}

#get randomised exposure matrices for each signature
randomised_exposures=vector(mode='list', length=length(SNV_sigs))
set.seed(1)
print("Randomising labels")
for (i in 1:length(SNV_sigs)){
    print(i)
    sig=SNV_sigs[i]
    activities=data.frame(Samples=SNV_activities$Samples,Signature=SNV_activities[,sig])
    for (j in 1:iter){
          column=sample(activities[,2],size=21,replace = F)
          activities=cbind(activities, column)}
    randomised_exposures[[i]]=activities
    colnames(randomised_exposures[[i]])=c("Samples",sig,paste0("i",1:iter))}

group_means=vector(mode='list', length=length(SNV_sigs))
print("\n \n Calculating means")
for(i in 1:10){
        print(i)
        for (cell in cells){
                for (treatment in treatments){
                        df=randomised_exposures[[i]][which(randomised_exposures[[i]]$Samples==paste0(cell,"_",treatment)),]
                        row=c(paste0(cell,"_",treatment),colMeans(df[,2:iter2]))
                        group_means[[i]]=rbind(group_means[[i]],row)}}
        for (treatment in treatments){
              df=randomised_exposures[[i]][grep(treatment,randomised_exposures[[i]][,1]),]
              row=c(treatment,colMeans(df[,2:iter2]))
              group_means[[i]]=rbind(group_means[[i]],row)}
        for (cell in cells){
            df=randomised_exposures[[i]][grep(cell,randomised_exposures[[i]][,1]),]
            row=c(cell,colMeans(df[,2:iter2]))
            group_means[[i]]=rbind(group_means[[i]],row)}}

mean_differences=vector(mode='list', length=length(SNV_sigs))
print("\n \n Calculating differences")
for (i in 1:length(SNV_sigs)){
    print(i)
    df=as.data.frame(group_means[[i]],stringsAsFactors=F)
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Helium_All",as.numeric(df[which(df$V1=="Control"),2:iter2])- as.numeric(df[which(df$V1=="Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Protons_All",as.numeric(df[which(df$V1=="Control"),2:iter2])- as.numeric(df[which(df$V1=="Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Helium-Protons_All",as.numeric(df[which(df$V1=="Helium"),2:iter2])- as.numeric(df[which(df$V1=="Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Helium_A549",as.numeric(df[which(df$V1=="A549_Control"),2:iter2])- as.numeric(df[which(df$V1=="A549_Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Protons_A549",as.numeric(df[which(df$V1=="A549_Control"),2:iter2])- as.numeric(df[which(df$V1=="A549_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Helium-Protons_A549",as.numeric(df[which(df$V1=="A549_Helium"),2:iter2])- as.numeric(df[which(df$V1=="A549_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Helium_HAP1",as.numeric(df[which(df$V1=="HAP1_Control"),2:iter2])- as.numeric(df[which(df$V1=="HAP1_Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Protons_HAP1",as.numeric(df[which(df$V1=="HAP1_Control"),2:iter2])- as.numeric(df[which(df$V1=="HAP1_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Helium-Protons_HAP1",as.numeric(df[which(df$V1=="HAP1_Helium"),2:iter2])- as.numeric(df[which(df$V1=="HAP1_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Helium_MCF7",as.numeric(df[which(df$V1=="MCF7_Control"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Protons_MCF7",as.numeric(df[which(df$V1=="MCF7_Control"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Helium-Protons_MCF7",as.numeric(df[which(df$V1=="MCF7_Helium"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-HAP1_All",as.numeric(df[which(df$V1=="A549"),2:iter2])- as.numeric(df[which(df$V1=="HAP1"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-MCF7_All",as.numeric(df[which(df$V1=="A549"),2:iter2])- as.numeric(df[which(df$V1=="MCF7"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("HAP1-MCF7_All",as.numeric(df[which(df$V1=="HAP1"),2:iter2])- as.numeric(df[which(df$V1=="MCF7"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-HAP1_Control",as.numeric(df[which(df$V1=="A549_Control"),2:iter2])- as.numeric(df[which(df$V1=="HAP1_Control"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-MCF7_Control",as.numeric(df[which(df$V1=="A549_Control"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Control"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("HAP1-MCF7_Control",as.numeric(df[which(df$V1=="HAP1_Control"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Control"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-HAP1_Helium",as.numeric(df[which(df$V1=="A549_Helium"),2:iter2])- as.numeric(df[which(df$V1=="HAP1_Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-MCF7_Helium",as.numeric(df[which(df$V1=="A549_Helium"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("HAP1-MCF7_Helium",as.numeric(df[which(df$V1=="HAP1_Helium"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-HAP1_Protons",as.numeric(df[which(df$V1=="A549_Protons"),2:iter2])- as.numeric(df[which(df$V1=="HAP1_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-MCF7_Protons",as.numeric(df[which(df$V1=="A549_Protons"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("HAP1-MCF7_Protons",as.numeric(df[which(df$V1=="HAP1_Protons"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Protons"),2:iter2])))
    colnames(mean_differences[[i]])=c("Samples",SNV_sigs[i],paste0("i",1:iter))}

print("\n \n Creating data frames")
for (i in 1:length(SNV_sigs)){
    print(i)
    pivot_longer(as.data.frame(mean_differences[[i]]),cols=2:iter2, names_to = "Iteration", values_to = "Mean_diff") %>% mutate(Mean_diff=as.numeric(as.character(Mean_diff))) -> df
    df$Iteration=rep(c("Real",rep("Random",times=iter)),times=24)
    df$Iteration=factor(df$Iteration,levels=c("Random","Real"))
    df$Samples=factor(df$Samples,levels=mean_differences[[1]][,1])
    df$Signature=SNV_sigs[i]
    segment_data = data.frame(x = seq(0.7,23.7,by=1), xend = seq(1.3,24.3,by=1), y = unname(df[df$Iteration=="Real", "Mean_diff"]), yend = unname(df[df$Iteration=="Real", "Mean_diff"]))
    segment_data$Signature=SNV_sigs[i]
    assign(paste0(SNV_sigs[i],"_df"),df)
    assign(paste0(SNV_sigs[i],"_segments"),segment_data)}

all_signatures=rbind(SBS96A_df,SBS96B_df,SBS96C_df,SBS96D_df,SBS96E_df,SBS96F_df,SBS96G_df,SBS96H_df,SBS96I_df,SBS96J_df)
all_segments=rbind(SBS96A_segments,SBS96B_segments,SBS96C_segments,SBS96D_segments,SBS96E_segments,SBS96F_segments,SBS96G_segments,SBS96H_segments,SBS96I_segments,SBS96J_segments)

p_values=data.frame(Signature=c(),Sample=c(),P_value=c(),H0=c(),stringsAsFactors = F)
print("\n \n Calculating empirical p values")
for (sig in SNV_sigs){
            for (sample in unique(all_signatures$Samples)){
                    df=all_signatures[which(all_signatures$Signature==sig & all_signatures$Samples==sample),]
                    real=deframe(df[which(df$Iteration=="Real"),"Mean_diff"])
                    random=deframe(df[which(df$Iteration=="Random"),"Mean_diff"])
                    is.higher=(length(which(real>=random))+1)/iter1
                    is.lower=(length(which(real<=random))+1)/iter1
                    data.frame(Signature=sig,Sample=sample,P_value=is.higher,H0="Higher")
                    p_values=rbind(p_values,data.frame(Signature=sig,Sample=sample,P_value=is.higher,H0="Higher"))
                    p_values=rbind(p_values,data.frame(Signature=sig,Sample=sample,P_value=is.lower,H0="Lower"))}}

p_values$P_adjusted=p.adjust(p_values$P_value,method="fdr") #get FDR correction
p_values_treatments=p_values[which(p_values$Sample %in% unique(p_values$Sample)[1:12]),] #take only treatment v treatment comparisons
p_values_treatments$P_adjusted=p.adjust(p_values_treatments$P_value,method="fdr") #get FDR correction for treatments

print("\n \n Exporting files")
if (norm==FALSE){
    write.table(all_signatures, file = "all_signatures_SNVs.txt", sep = " ",col.names = TRUE,quote=F)
    write.table(all_segments, file = "all_segments_SNVs.txt", sep = " ",col.names = TRUE,quote=F)
    write.table(p_values, file = "all_comparisons_SNVs_sigs.txt", sep = " ",col.names = TRUE,quote=F)
    write.table(p_values_treatments, file = "treatments_SNVs_sigs.txt", sep = " ",col.names = TRUE,quote=F)
    }else{
    write.table(all_signatures, file = "all_signatures_SNVs_norm.txt", sep = " ",col.names = TRUE,quote=F)
    write.table(all_segments, file = "all_segments_SNVs_norm.txt", sep = " ",col.names = TRUE,quote=F)
    write.table(p_values, file = "all_comparisons_SNVs_sigs_norm.txt", sep = " ",col.names = TRUE,quote=F)
    write.table(p_values_treatments, file = "treatments_SNVs_sigs_norm.txt", sep = " ",col.names = TRUE,quote=F)}

##############################################################################################################################################################################################################
##############################################################################################################################################################################################################
##############################################################################################################################################################################################################

ID_activities=as.data.frame(fread("ID83_S4_NMF_Activities.txt")) #Indel signature activities file (from SigProfiler)
colnames(ID_activities)=c("Samples","ID83A","ID83C","ID83B","ID83D") #reorder
ID_activities[,2:5] %>% select(sort(names(.))) -> ID_activities[,2:5]
ID_activities=ID_activities[1:21,]
colnames(ID_activities)=c("Samples","ID83A","ID83B","ID83C","ID83D")

ID_activities$Samples=mapvalues(ID_activities$Samples, from=c("A549_CA", "A549_CC", "A549_HE1A","A549_HE2B","A549_PR1B","A549_PR2B",
                                                                "HAP1_CA", "HAP1_CB", "HAP1_HE1A","HAP1_HE2A","HAP1_PR1A","HAP1_PR2A",
                                                                "MCF7_C", "MCF7_HE1B", "MCF7_HE1C","MCF7_HE2A","MCF7_HE2B","MCF7_PR1A","MCF7_PR1B","MCF7_PR2A","MCF7_PR2C"),
                                 to=c("A549_Control", "A549_Control", "A549_Helium","A549_Helium","A549_Protons","A549_Protons",
                                      "HAP1_Control", "HAP1_Control", "HAP1_Helium","HAP1_Helium","HAP1_Protons","HAP1_Protons",
                                      "MCF7_Control", "MCF7_Helium","MCF7_Helium","MCF7_Helium","MCF7_Helium","MCF7_Protons","MCF7_Protons","MCF7_Protons","MCF7_Protons"))
if (norm==TRUE){ID_activities[, 2:5]=sweep(ID_activities[, 2:5],MARGIN=1,FUN="/",STATS=rowSums(ID_activities[, 2:5]))}

#get randomised exposure matrices for each signature
randomised_exposures=vector(mode='list', length=length(ID_sigs))
set.seed(1)
print("Randomising labels")
for (i in 1:length(ID_sigs)){
    print(i)
    sig=ID_sigs[i]
    activities=data.frame(Samples=ID_activities$Samples,Signature=ID_activities[,sig])
    for (j in 1:iter){
        column=sample(activities[,2],size=21,replace = F)
        activities=cbind(activities, column)}
    randomised_exposures[[i]]=activities
    colnames(randomised_exposures[[i]])=c("Samples",sig,paste0("i",1:iter))}

group_means=vector(mode='list', length=length(ID_sigs))
print("\n \n Calculating means")
for(i in 1:length(ID_sigs)){
    print(i)
    for (cell in cells){
        for (treatment in treatments){
            df=randomised_exposures[[i]][which(randomised_exposures[[i]]$Samples==paste0(cell,"_",treatment)),]
            row=c(paste0(cell,"_",treatment),colMeans(df[,2:iter2]))
            group_means[[i]]=rbind(group_means[[i]],row)}}
    for (treatment in treatments){
        df=randomised_exposures[[i]][grep(treatment,randomised_exposures[[i]][,1]),]
        row=c(treatment,colMeans(df[,2:iter2]))
        group_means[[i]]=rbind(group_means[[i]],row)}
    for (cell in cells){
        df=randomised_exposures[[i]][grep(cell,randomised_exposures[[i]][,1]),]
        row=c(cell,colMeans(df[,2:iter2]))
        group_means[[i]]=rbind(group_means[[i]],row)}}

mean_differences=vector(mode='list', length=length(ID_sigs))
print("\n \n Calculating differences")
for (i in 1:length(ID_sigs)){
    print(i)
    df=as.data.frame(group_means[[i]],stringsAsFactors=F)
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Helium_All",as.numeric(df[which(df$V1=="Control"),2:iter2])- as.numeric(df[which(df$V1=="Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Protons_All",as.numeric(df[which(df$V1=="Control"),2:iter2])- as.numeric(df[which(df$V1=="Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Helium-Protons_All",as.numeric(df[which(df$V1=="Helium"),2:iter2])- as.numeric(df[which(df$V1=="Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Helium_A549",as.numeric(df[which(df$V1=="A549_Control"),2:iter2])- as.numeric(df[which(df$V1=="A549_Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Protons_A549",as.numeric(df[which(df$V1=="A549_Control"),2:iter2])- as.numeric(df[which(df$V1=="A549_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Helium-Protons_A549",as.numeric(df[which(df$V1=="A549_Helium"),2:iter2])- as.numeric(df[which(df$V1=="A549_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Helium_HAP1",as.numeric(df[which(df$V1=="HAP1_Control"),2:iter2])- as.numeric(df[which(df$V1=="HAP1_Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Protons_HAP1",as.numeric(df[which(df$V1=="HAP1_Control"),2:iter2])- as.numeric(df[which(df$V1=="HAP1_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Helium-Protons_HAP1",as.numeric(df[which(df$V1=="HAP1_Helium"),2:iter2])- as.numeric(df[which(df$V1=="HAP1_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Helium_MCF7",as.numeric(df[which(df$V1=="MCF7_Control"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Protons_MCF7",as.numeric(df[which(df$V1=="MCF7_Control"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Helium-Protons_MCF7",as.numeric(df[which(df$V1=="MCF7_Helium"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-HAP1_All",as.numeric(df[which(df$V1=="A549"),2:iter2])- as.numeric(df[which(df$V1=="HAP1"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-MCF7_All",as.numeric(df[which(df$V1=="A549"),2:iter2])- as.numeric(df[which(df$V1=="MCF7"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("HAP1-MCF7_All",as.numeric(df[which(df$V1=="HAP1"),2:iter2])- as.numeric(df[which(df$V1=="MCF7"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-HAP1_Control",as.numeric(df[which(df$V1=="A549_Control"),2:iter2])- as.numeric(df[which(df$V1=="HAP1_Control"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-MCF7_Control",as.numeric(df[which(df$V1=="A549_Control"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Control"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("HAP1-MCF7_Control",as.numeric(df[which(df$V1=="HAP1_Control"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Control"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-HAP1_Helium",as.numeric(df[which(df$V1=="A549_Helium"),2:iter2])- as.numeric(df[which(df$V1=="HAP1_Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-MCF7_Helium",as.numeric(df[which(df$V1=="A549_Helium"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("HAP1-MCF7_Helium",as.numeric(df[which(df$V1=="HAP1_Helium"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-HAP1_Protons",as.numeric(df[which(df$V1=="A549_Protons"),2:iter2])- as.numeric(df[which(df$V1=="HAP1_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-MCF7_Protons",as.numeric(df[which(df$V1=="A549_Protons"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("HAP1-MCF7_Protons",as.numeric(df[which(df$V1=="HAP1_Protons"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Protons"),2:iter2])))
    colnames(mean_differences[[i]])=c("Samples",ID_sigs[i],paste0("i",1:iter))}

print("\n \n Creating data frames")
for (i in 1:length(ID_sigs)){
    print(i)
    pivot_longer(as.data.frame(mean_differences[[i]]),cols=2:iter2, names_to = "Iteration", values_to = "Mean_diff") %>% mutate(Mean_diff=as.numeric(as.character(Mean_diff))) -> df
    df$Iteration=rep(c("Real",rep("Random",times=iter)),times=24)
    df$Iteration=factor(df$Iteration,levels=c("Random","Real"))
    df$Samples=factor(df$Samples,levels=mean_differences[[1]][,1])
    df$Signature=ID_sigs[i]
    segment_data = data.frame(x = seq(0.7,23.7,by=1), xend = seq(1.3,24.3,by=1), y = unname(df[df$Iteration=="Real", "Mean_diff"]), yend = unname(df[df$Iteration=="Real", "Mean_diff"]))
    segment_data$Signature=ID_sigs[i]
    assign(paste0(ID_sigs[i],"_df"),df)
    assign(paste0(ID_sigs[i],"_segments"),segment_data)}

all_signatures=rbind(ID83A_df,ID83B_df,ID83C_df,ID83D_df)
all_segments=rbind(ID83A_segments,ID83B_segments,ID83C_segments,ID83D_segments)

p_values=data.frame(Signature=c(),Sample=c(),P_value=c(),H0=c(),stringsAsFactors = F)
print("\n \n Calculating empirical p values")
for (sig in ID_sigs){
    for (sample in unique(all_signatures$Samples)){
        df=all_signatures[which(all_signatures$Signature==sig & all_signatures$Samples==sample),]
        real=deframe(df[which(df$Iteration=="Real"),"Mean_diff"])
        random=deframe(df[which(df$Iteration=="Random"),"Mean_diff"])
        is.higher=(length(which(real>=random))+1)/iter1
        is.lower=(length(which(real<=random))+1)/iter1
        data.frame(Signature=sig,Sample=sample,P_value=is.higher,H0="Higher")
        p_values=rbind(p_values,data.frame(Signature=sig,Sample=sample,P_value=is.higher,H0="Higher"))
        p_values=rbind(p_values,data.frame(Signature=sig,Sample=sample,P_value=is.lower,H0="Lower"))}}

p_values$P_adjusted=p.adjust(p_values$P_value,method="fdr") #get FDR correction
p_values_treatments=p_values[which(p_values$Sample %in% unique(p_values$Sample)[1:12]),] #take only treatment v treatment comparisons
p_values_treatments$P_adjusted=p.adjust(p_values_treatments$P_value,method="fdr") #get FDR correction for treatments

print("\n \n Exporting files")
if (norm==FALSE){
    write.table(all_signatures, file = "all_signatures_IDs.txt", sep = " ",col.names = TRUE,quote=F)
    write.table(all_segments, file = "all_segments_IDs.txt", sep = " ",col.names = TRUE,quote=F)
    write.table(p_values, file = "all_comparisons_IDs_sigs.txt", sep = " ",col.names = TRUE,quote=F)
    write.table(p_values_treatments, file = "treatments_IDs_sigs.txt", sep = " ",col.names = TRUE,quote=F)
}else{
    write.table(all_signatures, file = "all_signatures_IDs_norm.txt", sep = " ",col.names = TRUE,quote=F)
    write.table(all_segments, file = "all_segments_IDs_norm.txt", sep = " ",col.names = TRUE,quote=F)
    write.table(p_values, file = "all_comparisons_IDs_sigs_norm.txt", sep = " ",col.names = TRUE,quote=F)
    write.table(p_values_treatments, file = "treatments_IDs_sigs_norm.txt", sep = " ",col.names = TRUE,quote=F)}

##############################################################################################################################################################################################################
##############################################################################################################################################################################################################
##############################################################################################################################################################################################################

SNVs=read.table("Our_Zou_Kucab_SNVs.txt",header = T)[,1:22]
IDs=read.table("Our_Zou_Kucab_indels.txt",header = T)[,1:22]
mutation_counts=data.frame(Samples=colnames(SNVs)[-1],SNVs=unname(colSums(SNVs[,2:22])),IDs=unname(colSums(IDs[,2:22])))
SVs=read.table("SV_counts.txt",header=T,sep=" ")
sv_types=c("DEL","DUP","INS","INV")
sv_counts=data.frame(Samples=colnames(SNVs)[-1],DEL=rep(NA,times=21),DUP=rep(NA,times=21),INS=rep(NA,times=21),INV=rep(NA,times=21))


mutation_counts$Del1_Dels_ratio_no_reps=unname(colSums(IDs[c(1,7),2:22]))/unname(colSums(IDs[c(25,31,37,43),2:22]))
mutation_counts$Ins1_Inss_ratio_no_reps=unname(colSums(IDs[c(13,19),2:22]))/unname(colSums(IDs[c(49,55,60,67),2:22]))
mutation_counts$MH_Dels_Del5_ratio=unname(colSums(IDs[73:78,2:22]))/unname(colSums(IDs[79:83,2:22]))
mutation_counts$MH1_Del_MH_Dels_ratio=unname(colSums(IDs[c(73,74,76,79),2:22]))/unname(colSums(IDs[c(75,77,78,80,81,82,83),2:22]))

IDs$MutationType=sapply(str_split(IDs$MutationType,":"), "[[", 2)
pivot_longer(IDs, cols=2:22, names_to = "Sample", values_to = "Mutations") %>% group_by(Sample, MutationType) %>% summarise(across(Mutations, sum)) %>% filter(MutationType == "Ins") %>% pull(Mutations) -> ins
pivot_longer(IDs, cols=2:22, names_to = "Sample", values_to = "Mutations") %>% group_by(Sample, MutationType) %>% summarise(across(Mutations, sum)) %>% filter(MutationType == "Del") %>% pull(Mutations) -> del
mutation_counts$IDs_del_ins_ratio=del/ins

for (sample in sv_counts$Samples){
    df=SVs[which(SVs$Clone==sample),]
    i=which(sv_counts$Samples %in% sample)
    for (sv in sv_types){
                 if (sv %in% df$SV){
                        temp=df[which(df$SV==sv),]
                        sv_counts[i,sv]=temp$Count}else{sv_counts[i,sv]=0}}}
mutation_counts=cbind(mutation_counts,sv_counts[,2:5])

mutation_counts$Samples=mapvalues(mutation_counts$Samples, from=c("A549_CA", "A549_CC", "A549_HE1A","A549_HE2B","A549_PR1B","A549_PR2B",
                                                              "HAP1_CA", "HAP1_CB", "HAP1_HE1A","HAP1_HE2A","HAP1_PR1A","HAP1_PR2A",
                                                              "MCF7_C", "MCF7_HE1B", "MCF7_HE1C","MCF7_HE2A","MCF7_HE2B","MCF7_PR1A","MCF7_PR1B","MCF7_PR2A","MCF7_PR2C"),
                                to=c("A549_Control", "A549_Control", "A549_Helium","A549_Helium","A549_Protons","A549_Protons",
                                     "HAP1_Control", "HAP1_Control", "HAP1_Helium","HAP1_Helium","HAP1_Protons","HAP1_Protons",
                                     "MCF7_Control", "MCF7_Helium","MCF7_Helium","MCF7_Helium","MCF7_Helium","MCF7_Protons","MCF7_Protons","MCF7_Protons","MCF7_Protons"))
mutation_counts$Del_InsDup_ratio=mutation_counts$DEL/(mutation_counts$INS+mutation_counts$DUP)
mut_types=c("SNVs","IDs","IDs_del_ins_ratio","Del1_Dels_ratio_no_reps","Ins1_Inss_ratio_no_reps","MH_Dels_Del5_ratio","MH1_Del_MH_Dels_ratio",sv_types,"Del_InsDup_ratio")


#get randomised exposure matrices for each signature
randomised_exposures=vector(mode='list', length=length(mut_types))
set.seed(1)
print("Randomising labels")
for (i in 1:length(mut_types)){
    print(i)
    sig=mut_types[i]
    activities=data.frame(Samples=mutation_counts$Samples,Signature=mutation_counts[,sig])
    for (j in 1:iter){
        column=sample(activities[,2],size=21,replace = F)
        activities=cbind(activities, column)}
    randomised_exposures[[i]]=activities
    colnames(randomised_exposures[[i]])=c("Samples",sig,paste0("i",1:iter))}


group_means=vector(mode='list', length=length(mut_types))
print("\n \n Calculating means")
for(i in 1:length(mut_types)){
    print(i)
    for (cell in cells){
        for (treatment in treatments){
            df=randomised_exposures[[i]][which(randomised_exposures[[i]]$Samples==paste0(cell,"_",treatment)),]
            row=c(paste0(cell,"_",treatment),colMeans(df[,2:iter2]))
            group_means[[i]]=rbind(group_means[[i]],row)}}
    for (treatment in treatments){
        df=randomised_exposures[[i]][grep(treatment,randomised_exposures[[i]][,1]),]
        row=c(treatment,colMeans(df[,2:iter2]))
        group_means[[i]]=rbind(group_means[[i]],row)}
    for (cell in cells){
        df=randomised_exposures[[i]][grep(cell,randomised_exposures[[i]][,1]),]
        row=c(cell,colMeans(df[,2:iter2]))
        group_means[[i]]=rbind(group_means[[i]],row)}}

mean_differences=vector(mode='list', length=length(mut_types))
print("\n \n Calculating differences")
for (i in 1:length(mut_types)){
    print(i)
    df=as.data.frame(group_means[[i]],stringsAsFactors=F)
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Helium_All",as.numeric(df[which(df$V1=="Control"),2:iter2])- as.numeric(df[which(df$V1=="Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Protons_All",as.numeric(df[which(df$V1=="Control"),2:iter2])- as.numeric(df[which(df$V1=="Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Helium-Protons_All",as.numeric(df[which(df$V1=="Helium"),2:iter2])- as.numeric(df[which(df$V1=="Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Helium_A549",as.numeric(df[which(df$V1=="A549_Control"),2:iter2])- as.numeric(df[which(df$V1=="A549_Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Protons_A549",as.numeric(df[which(df$V1=="A549_Control"),2:iter2])- as.numeric(df[which(df$V1=="A549_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Helium-Protons_A549",as.numeric(df[which(df$V1=="A549_Helium"),2:iter2])- as.numeric(df[which(df$V1=="A549_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Helium_HAP1",as.numeric(df[which(df$V1=="HAP1_Control"),2:iter2])- as.numeric(df[which(df$V1=="HAP1_Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Protons_HAP1",as.numeric(df[which(df$V1=="HAP1_Control"),2:iter2])- as.numeric(df[which(df$V1=="HAP1_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Helium-Protons_HAP1",as.numeric(df[which(df$V1=="HAP1_Helium"),2:iter2])- as.numeric(df[which(df$V1=="HAP1_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Helium_MCF7",as.numeric(df[which(df$V1=="MCF7_Control"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Control-Protons_MCF7",as.numeric(df[which(df$V1=="MCF7_Control"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("Helium-Protons_MCF7",as.numeric(df[which(df$V1=="MCF7_Helium"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-HAP1_All",as.numeric(df[which(df$V1=="A549"),2:iter2])- as.numeric(df[which(df$V1=="HAP1"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-MCF7_All",as.numeric(df[which(df$V1=="A549"),2:iter2])- as.numeric(df[which(df$V1=="MCF7"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("HAP1-MCF7_All",as.numeric(df[which(df$V1=="HAP1"),2:iter2])- as.numeric(df[which(df$V1=="MCF7"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-HAP1_Control",as.numeric(df[which(df$V1=="A549_Control"),2:iter2])- as.numeric(df[which(df$V1=="HAP1_Control"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-MCF7_Control",as.numeric(df[which(df$V1=="A549_Control"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Control"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("HAP1-MCF7_Control",as.numeric(df[which(df$V1=="HAP1_Control"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Control"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-HAP1_Helium",as.numeric(df[which(df$V1=="A549_Helium"),2:iter2])- as.numeric(df[which(df$V1=="HAP1_Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-MCF7_Helium",as.numeric(df[which(df$V1=="A549_Helium"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("HAP1-MCF7_Helium",as.numeric(df[which(df$V1=="HAP1_Helium"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Helium"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-HAP1_Protons",as.numeric(df[which(df$V1=="A549_Protons"),2:iter2])- as.numeric(df[which(df$V1=="HAP1_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("A549-MCF7_Protons",as.numeric(df[which(df$V1=="A549_Protons"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Protons"),2:iter2])))
    mean_differences[[i]]=rbind(mean_differences[[i]], c("HAP1-MCF7_Protons",as.numeric(df[which(df$V1=="HAP1_Protons"),2:iter2])- as.numeric(df[which(df$V1=="MCF7_Protons"),2:iter2])))
    colnames(mean_differences[[i]])=c("Samples",mut_types[i],paste0("i",1:iter))}

print("\n \n Creating data frames")
for (i in 1:length(mut_types)){
    print(i)
    pivot_longer(as.data.frame(mean_differences[[i]]),cols=2:iter2, names_to = "Iteration", values_to = "Mean_diff") %>% mutate(Mean_diff=as.numeric(as.character(Mean_diff))) -> df
    df$Iteration=rep(c("Real",rep("Random",times=iter)),times=24)
    df$Iteration=factor(df$Iteration,levels=c("Random","Real"))
    df$Samples=factor(df$Samples,levels=mean_differences[[1]][,1])
    df$Signature=mut_types[i]
    segment_data = data.frame(x = seq(0.7,23.7,by=1), xend = seq(1.3,24.3,by=1), y = unname(df[df$Iteration=="Real", "Mean_diff"]), yend = unname(df[df$Iteration=="Real", "Mean_diff"]))
    segment_data$Signature=mut_types[i]
    assign(paste0(mut_types[i],"_df"),df)
    assign(paste0(mut_types[i],"_segments"),segment_data)}   

all_signatures=rbind(SNVs_df,IDs_df,IDs_del_ins_ratio_df,Del1_Dels_ratio_no_reps_df,Ins1_Inss_ratio_no_reps_df,MH_Dels_Del5_ratio_df,MH1_Del_MH_Dels_ratio_df,DEL_df,DUP_df,INS_df,INV_df,Del_InsDup_ratio_df)
all_segments=rbind(SNVs_segments,IDs_segments,IDs_del_ins_ratio_segments,Del1_Dels_ratio_no_reps_segments,Ins1_Inss_ratio_no_reps_segments,MH_Dels_Del5_ratio_segments,MH1_Del_MH_Dels_ratio_segments,DEL_segments,DUP_segments,INS_segments,INV_segments,Del_InsDup_ratio_segments)

p_values=data.frame(Signature=c(),Sample=c(),P_value=c(),H0=c(),stringsAsFactors = F)
print("\n \n Calculating empirical p values")
for (sig in mut_types){
    for (sample in unique(all_signatures$Samples)){
        df=all_signatures[which(all_signatures$Signature==sig & all_signatures$Samples==sample),]
        real=deframe(df[which(df$Iteration=="Real"),"Mean_diff"])
        random=deframe(df[which(df$Iteration=="Random"),"Mean_diff"])
        is.higher=(length(which(real>=random))+1)/iter1
        is.lower=(length(which(real<=random))+1)/iter1
        data.frame(Signature=sig,Sample=sample,P_value=is.higher,H0="Higher")
        p_values=rbind(p_values,data.frame(Signature=sig,Sample=sample,P_value=is.higher,H0="Higher"))
        p_values=rbind(p_values,data.frame(Signature=sig,Sample=sample,P_value=is.lower,H0="Lower"))}}

p_values$P_adjusted=p.adjust(p_values$P_value,method="fdr") #get FDR correction
p_values_treatments=p_values[which(p_values$Sample %in% unique(p_values$Sample)[1:12]),] #take only treatment v treatment comparisons 
p_values_treatments$P_adjusted=p.adjust(p_values_treatments$P_value,method="fdr") #get FDR correction for treatments

print("\n \n Exporting files")
    write.table(all_signatures, file = "all_signatures_TMBs.txt", sep = " ",col.names = TRUE,quote=F)
    write.table(all_segments, file = "all_segments_TMBs.txt", sep = " ",col.names = TRUE,quote=F)
    write.table(p_values, file = "all_comparisons_TMBs_sigs.txt", sep = " ",col.names = TRUE,quote=F)
    write.table(p_values_treatments, file = "treatments_TMBs_sigs.txt", sep = " ",col.names = TRUE,quote=F)

