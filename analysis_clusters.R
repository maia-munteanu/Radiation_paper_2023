library(tidygraph)
library(ggpubr)
library(gridExtra)
library(plyr)
library(ggraph)
library(igraph)
library(MutationalPatterns)
library(dplyr)
library(tibble)
library(tidyr)

#read files from Clustered_SNVs.zip
all_SNVs_HAP1_CA=read.table("HAP1_CA.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_HAP1_CB=read.table("HAP1_CB.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_HAP1_HE1A=read.table("HAP1_HE1A.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_HAP1_HE2A=read.table("HAP1_HE2A.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_HAP1_PR1A=read.table("HAP1_PR1A.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_HAP1_PR2A=read.table("HAP1_PR2A.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_A549_CA=read.table("A549_CA.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_A549_CC=read.table("A549_CC.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_A549_HE1A=read.table("A549_HE1A.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_A549_HE2B=read.table("A549_HE2B.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_A549_PR1B=read.table("A549_PR1B.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_A549_PR2B=read.table("A549_PR2B.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_MCF7_untreated_F8_REP2=read.table("MCF7_C.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_MCF7_HE_1B=read.table("MCF7_HE1B.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_MCF7_HE_1C=read.table("MCF7_HE1C.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_MCF7_HE_2A=read.table("MCF7_HE2A.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_MCF7_HE_2B=read.table("MCF7_HE2B.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_MCF7_PR_1A=read.table("MCF7_PR1A.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_MCF7_PR_1B=read.table("MCF7_PR1B.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_MCF7_PR_2A=read.table("MCF7_PR2A.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)
all_SNVs_MCF7_PR_2C=read.table("MCF7_PR2C.txt",header=T) %>% mutate(Variants=paste(CHROM,":",POS,"_",REF,"/",ALT,sep="")) %>% pull(Variants)

graph_clusters <- function(list_mutations, max_len = 1000){ # CAUTION should be the SAME chromosome
  allchr = unlist(lapply(list_mutations, function(mut) unlist(strsplit(mut, ":"))[1]))
  if(length(table(allchr)) > 1) stop("input mutations should be in the same chromosome")
  allpos = as.numeric(unlist(lapply(list_mutations, function(mut) unlist(strsplit(unlist(strsplit(mut, ":"))[2], "_"))[1])))
  allmuttype =  unlist(lapply(list_mutations, function(mut) unlist(strsplit(unlist(strsplit(mut, ":"))[2], "_"))[2]))
  if(exists("edges0")) rm(list="edges0")
  if(exists("nodes0")) rm(list="nodes0")
  if(exists("alldistances")) rm(list="alldistances")
  for(i in 1:length(list_mutations)){
    for(j in 1:length(list_mutations)){
      cur_distance = max(allpos[i], allpos[j]) - min(allpos[i], allpos[j])
      if(cur_distance < max_len & (i != j)) if(exists("edges0")) { 
        edges0 = rbind(edges0, t(data.frame(c(i, j)))) ; nodes = unique(c(nodes, i ,j)) ; alldistances = c(alldistances, cur_distance)
      } else { edges0 = t(data.frame(c(i, j))) ; nodes = unique(i, j) ; alldistances = cur_distance} }}
  if(exists("edges0")){ 
    nodes0 = data.frame(unlist(nodes))
    muts = list_mutations[unique(unlist(nodes))] ; names(muts) = unique(unlist(nodes)) #list_mutations[unique(edges0[,1], edges0[,2])]
    return(list(edges0=edges0, nodes0=nodes0, allmuttype=allmuttype, alldistances = alldistances, muts = muts)) } else { return(NA) }}

plot_graph <- function(nodes, edges, allmuttype, main="", mycols = null){
  routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
  routes_igraph2 <- simplify(routes_igraph)
  V(routes_igraph2)$color <- alpha(mycols[allmuttype[nodes[,1]]], 0.75)
  plot.igraph(routes_igraph2, edge.arrow.size = 0.2, main=main, vertex.size= 40)}

for(cellline in c("A549", "HAP1", "MCF7")){
  print(paste("Starting cell line: ", cellline, sep=""))
  alldat = list(ls(pattern = paste("all_SNVs_", cellline, sep="")))[[1]]
  allsm = unlist(lapply(alldat, function(x) unlist(strsplit(x, paste(cellline,"_",sep="")))[2]))
  if(exists("allmuts_largeclust")) rm("allmuts_largeclust")
  for(sm in allsm){
    print(sm)
    allmut0 = get(paste("all_SNVs", cellline, sm, sep = "_"))
    if(exists(paste("dist_", sm, sep=""))) rm(list=paste("dist_", sm, sep=""))
    for(chr in paste("chr", 1:22, sep="")){
      print(chr)
      allmutchr = allmut0[which(grepl(paste(chr, ":", sep=""), allmut0))]
      res = graph_clusters(allmutchr)
      if(!is.na(res)){
        compres = components(graph_from_data_frame(d = res$edges0, vertices = res$nodes0, directed = FALSE)); comps = compres$csize
        if(sum(comps > 5)>=1){ # we have at least one large cluster: save the mutations
          save_clust_mut = T
          for(largeclust in which(comps > 5)){
            id_muts_largeclust = names(compres$membership)[which(compres$membership == largeclust)]
            muts_largeclust0 = as.character(res$muts[id_muts_largeclust])
            if(largeclust == which(comps > 5)[1]) {assign("muts_largeclust", muts_largeclust0)} else {
              assign("muts_largeclust", c(muts_largeclust, muts_largeclust0))
            }
          }
        } else { save_clust_mut = F }
      }
      if(!is.na(res)) { if(!exists(paste("dist_", sm, sep=""))) { # if exists, we already build the vectors for previous chrs
        assign(paste("dist_", sm, sep=""), res$alldistances)
        assign(paste("muts_", sm, sep=""), res$muts)
        assign(paste("comps_", sm, sep=""), comps)
      } else { 
        assign(paste("dist_", sm, sep=""), c( get(paste("dist_", sm, sep="")), res$alldistances))
        assign(paste("muts_", sm, sep=""), c( get(paste("muts_", sm, sep="")), res$muts))
        assign(paste("comps_", sm, sep=""), c( get(paste("comps_", sm, sep="")), comps))
      }
      }
      if(save_clust_mut == T){
        if(exists(paste("muts_largeclust_", sm, cellline, sep=""))) {
          assign(paste("muts_largeclust_", sm, cellline, sep=""), c( get(paste("muts_largeclust_", sm, cellline, sep="")), muts_largeclust))
        } else { assign(paste("muts_largeclust_", sm, cellline, sep=""), muts_largeclust) }
      } 
    }
    if(match(sm, allsm) == 1) {
      alldist=data.frame(d=get(paste("dist_", sm, sep="")), sm=sm)
      allmuts=data.frame(d=get(paste("muts_", sm, sep="")), sm=sm)
      allcomps=data.frame(d=get(paste("comps_", sm, sep="")), sm=sm)
      if(exists(paste("muts_largeclust_", sm, cellline, sep=""))) allmuts_largeclust=data.frame(d=get(paste("muts_largeclust_", sm, cellline, sep="")), sm=sm)
      } else {
      alldist = rbind(alldist, data.frame(d=get(paste("dist_", sm, sep="")), sm=sm) )
      allmuts = rbind(allmuts, data.frame(d=get(paste("muts_", sm, sep="")), sm=sm) )
      allcomps = rbind(allcomps, data.frame(d=get(paste("comps_", sm, sep="")), sm=sm) )
      if(exists(paste("muts_largeclust_", sm, cellline, sep=""))){
        if(!exists("allmuts_largeclust")){ allmuts_largeclust=data.frame(d=get(paste("muts_largeclust_", sm, cellline, sep="")), sm=sm) }
        if(exists("allmuts_largeclust")){ allmuts_largeclust = rbind(allmuts_largeclust, data.frame(d=get(paste("muts_largeclust_", sm, cellline, sep="")), sm=sm))}
      } 
      }
  }
  assign(paste("alldist_", cellline, sep=""), alldist)
  assign(paste("allmuts_", cellline, sep=""), allmuts)
  assign(paste("allcomps_", cellline, sep=""), allcomps)
  assign(paste("allmuts_largeclust_", cellline, sep=""), allmuts_largeclust)
}

for (cell in c("A549", "HAP1", "MCF7")){
   df = get(paste("alldist_", cell, sep=""))
   df %>% group_by(sm,d) %>% summarize(count=n()) %>% mutate(count=count/2) %>% uncount(weights = count) -> df   #this removes duplicate distances caused by counting the same edge twice
   df$Class = rep("Control", nrow(df))
   df[which(grepl("HE", df$sm)), "Class"] = "Helium"
   df[which(grepl("PR", df$sm)), "Class"] = "Protons"
   df$d = log10(df$d)
   assign(paste("alldist2_", cell, sep=""), df)}

alldist2_A549$Cell="A549";alldist2_HAP1$Cell="HAP1";alldist2_MCF7$Cell="MCF7"
###alldist2_MCF7=alldist2_MCF7[-which(alldist2_MCF7$sm=="HE_2B"),] #removes the MCF7 HE2B clone
alldist2=rbind(alldist2_A549,alldist2_HAP1,alldist2_MCF7)
alldist2$Cell=factor(alldist2$Cell,levels=c("HAP1","A549","MCF7"))

mycomps=list(c("Control", "Helium"), c("Helium", "Protons"), c("Control", "Protons"))
give.n <- function(x){return(c(y = mean(fivenum(x)[3:4]), label = length(x)))}
ggviolin(alldist2, x = "Class", y = "d", fill = "Class", add = "boxplot", add.params = list(fill = "white"), palette="jco",trim=TRUE) + 
  stat_compare_means(comparisons = mycomps,step.increase=c(1,0.1,0.1),bracket.size = 0.55, size=4.2) +
  facet_wrap(~Cell)+
  ylab("Inter-mutational distance (log10)\n")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(vjust = 0.5, hjust=1, size = 15),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size=25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        strip.text.x = element_text(size=25,face="bold",margin = margin( b = 25, t = 0)),
        panel.background=element_rect(fill="white"),
        strip.background = element_rect(fill="white",color="white"),
        legend.position = "bottom",
        text = element_text(size = 20))+
  guides(fill = guide_legend(title = "Treatment"))+
  stat_summary(fun.data = give.n, geom = "text", size=4.5)
ggsave("Fig6A.pdf",plot = last_plot(),device="pdf", width = 5500,height = 2500,units = "px")


### plot numbers ###
for(cellline in c("A549", "HAP1", "MCF7")){
  dd = get(paste("allcomps_", cellline, sep=""))
  cl_sm = gsub(paste(cellline, "_", sep=""), "", unique(dd$sm))
  for(sm in cl_sm){
    curdat = dd[which(grepl(sm, dd$sm)), ]
    if(grepl("HE", sm)) { cond = "HE" } else { if(grepl("PR", sm)) { cond = "PR"} else { cond = "CT"} }
    if(cellline == "A549" & sm == cl_sm[1]) { alldat = data.frame(clone=sm, cellline=cellline, comp_size=curdat[,1], condition=cond) } else {
      alldat = rbind(alldat, data.frame(clone=sm, cellline=cellline, comp_size=curdat[,1], condition=cond))}}}
alldat$condition=mapvalues(alldat$condition, from=c("CT","PR","HE"),to=c("Control","Protons","Helium"))
alldat$cellline=factor(alldat$cellline,levels=c("HAP1","A549","MCF7"))
alldat$sample=paste0(alldat$cellline,"_",alldat$clone)
###alldat=alldat[-which(alldat$clone=="HE_2B" & alldat$cellline=="MCF7"),] #removes the MCF7 HE2B clone

ggplot(alldat, aes(x=condition, y=comp_size,fill=condition)) +
  scale_fill_manual(values=c("#0073C2","#EFC000","#868686"))+
  geom_jitter(size=3.5, shape=21,width = 0.35,height=0.3,alpha=0.8)+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~cellline,scales="free")+
  scale_y_continuous(breaks=c(2,4,6,8,10,12,14,16))+
  guides(fill = guide_legend(title = "Treatment"))+
  ylab("Number of mutations per cluster\n")+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(vjust = 0.5, hjust=1, size = 15),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size=25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.position = "bottom",
        strip.text.x = element_text(size=25,face="bold",margin = margin( b = 25, t = 0)),
        panel.background=element_rect(fill="white"),
        strip.background = element_rect(fill="white",color="white"),
        text = element_text(size = 20)) 
ggsave("Fig6B.pdf",plot = last_plot(),device="pdf", width = 5500,height = 2500,units = "px")


allmuts_A549$Mutation=gsub("/",">",sub("^.*_", "", allmuts_A549$d))
allmuts_A549$Cellline="A549"
allmuts_HAP1$Mutation=gsub("/",">",sub("^.*_", "", allmuts_HAP1$d))
allmuts_HAP1$Cellline="HAP1"
allmuts_MCF7$Mutation=gsub("/",">",sub("^.*_", "", allmuts_MCF7$d))
allmuts_MCF7$Cellline="MCF7"
allmuts=rbind(allmuts_A549,allmuts_HAP1,allmuts_MCF7)
allmuts$Clone=paste(allmuts$Cellline,"_",allmuts$sm,sep="")
allmuts$Mutation=mapvalues(allmuts$Mutation,from = c("A>C","A>G","A>T","G>A","G>C","G>T"),to=c("T>G","T>C","T>A","C>T","C>G","C>A"))
allmuts$Mutation <- factor(allmuts$Mutation, levels = c("T>G","T>C","T>A","C>T","C>G","C>A"))

allmuts_df=data.frame(Name=allmuts$Clone,Mutation=allmuts$Mutation)
allmuts_df %>% group_by(Name,Mutation) %>% dplyr::summarise(Count = n()) %>% ungroup() %>% add_row(Name="HAP1_CB",Mutation="T>G",Count=0) %>% mutate(Cell=substr(Name,1,4)) -> allmuts_df
allmuts_df$Clone=mapvalues(allmuts_df$Name, from=c("A549_CA", "A549_CC", "A549_HE1A","A549_HE2B","A549_PR1B","A549_PR2B",
                                                        "HAP1_CA", "HAP1_CB", "HAP1_HE1A","HAP1_HE2A","HAP1_PR1A","HAP1_PR2A",
                                                        "MCF7_untreated_F8_REP2", "MCF7_HE_1B", "MCF7_HE_1C","MCF7_HE_2A","MCF7_HE_2B","MCF7_PR_1A","MCF7_PR_1B","MCF7_PR_2A","MCF7_PR_2C"), 
                             to=c("Control_A","Control_C","Helium_1A","Helium_2B","Protons_1B","Protons_2B",
                                  "Control_A", "Control_B", "Helium_1A","Helium_2A","Protons_1A","Protons_2A",
                                  "Control_A", "Helium_1B","Helium_1C","Helium_2A","Helium_2B","Protons_1A","Protons_1B","Protons_2A","Protons_2C"))

allmuts_df$Cell <- factor(allmuts_df$Cell, levels = c("HAP1", "A549", "MCF7"))
allmuts_df$Mutation <- factor(allmuts_df$Mutation, levels = c("T>G","T>C","T>A","C>T","C>G","C>A"))
allmuts_df=with(allmuts_df, allmuts_df[order(Cell, Clone, Mutation),])
#get mutation count normalization by sample
allmuts_df %>% group_by(Name) %>% mutate(Mutation_norm = Count / sum(Count)) -> allmuts_df


counts_df=data.frame(Name=allmuts$Clone)
counts_df %>% group_by(Name) %>% dplyr::summarise(Count = n()) -> counts_df
counts_df$Type=rep("Mutations",times=21)

alldat %>% group_by(sample) %>% dplyr::summarise(Count = n()) -> comps
colnames(comps)=c("Name","Count")
comps$Type=rep("Clusters",times=21)
counts_df=rbind(counts_df,comps)
counts_df$Cell=rep(c("A549","A549","A549","A549","A549","A549","HAP1","HAP1","HAP1","HAP1","HAP1","HAP1","MCF7","MCF7","MCF7","MCF7","MCF7","MCF7","MCF7","MCF7","MCF7"),times=2)
counts_df$Name=mapvalues(counts_df$Name, from=c("A549_CA", "A549_CC", "A549_HE1A","A549_HE2B","A549_PR1B","A549_PR2B",
                                                   "HAP1_CA", "HAP1_CB", "HAP1_HE1A","HAP1_HE2A","HAP1_PR1A","HAP1_PR2A",
                                                   "MCF7_untreated_F8_REP2", "MCF7_HE_1B", "MCF7_HE_1C","MCF7_HE_2A","MCF7_HE_2B","MCF7_PR_1A","MCF7_PR_1B","MCF7_PR_2A","MCF7_PR_2C"), 
                           to=c("Control_A","Control_C","Helium_1A","Helium_2B","Protons_1B","Protons_2B",
                                "Control_A", "Control_B", "Helium_1A","Helium_2A","Protons_1A","Protons_2A",
                                "Control_A", "Helium_1B","Helium_1C","Helium_2A","Helium_2B","Protons_1A","Protons_1B","Protons_2A","Protons_2C"))
counts_df$Cell <- factor(counts_df$Cell, levels = c("HAP1", "A549", "MCF7"))
counts_df$Name <- factor(counts_df$Name,levels=sort(as.character(unique(counts_df$Name))))
###counts_df=counts_df[-which(counts_df$Cell=="MCF7" & counts_df$Name=="Helium_2B"),] #removes the MCF7 HE2B clone

ggplot(counts_df, aes(fill=Type, y=Count, x=Name)) + facet_grid(Type~Cell, scale="free") +
  geom_bar(position="stack", stat="identity",colour="black")+
  geom_text(aes(label=Count), position = position_stack(vjust = 0.5),size=5.5)+
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
  xlab("\nClone")
ggsave("SFig9A.pdf",plot = last_plot(),device="pdf",  width = 6000,height = 3000,units = "px")

all_SNVs_SBS6=read.table("all_SNVs_SBS6.txt",header=T,sep="\t") #file which contains the counts of all SNVs - by SBS6 classes
all_SNVs_SBS6$Cell <- factor(all_SNVs_SBS6$Cell, levels = c("HAP1", "A549", "MCF7"))
all_SNVs_SBS6$SBS6 <- factor(all_SNVs_SBS6$SBS6, levels = c("T>G","T>C","T>A","C>T","C>G","C>A"))
all_SNVs_SBS6=with(all_SNVs_SBS6, all_SNVs_SBS6[order(Cell, Samples, SBS6),])
allmuts_df$Total_muts=all_SNVs_SBS6$Mutations
allmuts_df %>% group_by(Name) %>% mutate(Total_muts_norm = Total_muts / sum(Total_muts)) -> allmuts_df
allmuts_df[allmuts_df == 0] <- NA
allmuts_df$log2Ratio=log2(allmuts_df$Count/allmuts_df$Total_muts)
allmuts_df$RatioNorm=allmuts_df$Mutation_norm/allmuts_df$Total_muts_norm
allmuts_df$Difference=allmuts_df$Mutation_norm - allmuts_df$Total_muts_norm


###allmuts_df=allmuts_df[-which(allmuts_df$Name=="MCF7_HE_2B"),]  #when removing the HE2B clone
###data_vline <- data.frame(Cell = c(rep("A549",5),rep("HAP1",5),rep("MCF7",7)), vline = c(seq(from=1.5, to=5.5, by = 1),seq(from=1.5, to=5.5, by = 1),seq(from=1.5, to=7.5, by = 1)))
data_vline <- data.frame(Cell = c(rep("A549",5),rep("HAP1",5),rep("MCF7",8)), vline = c(seq(from=1.5, to=5.5, by = 1),seq(from=1.5, to=5.5, by = 1),seq(from=1.5, to=8.5, by = 1)))
data_hline <- data.frame(Cell = c(rep("A549",5),rep("HAP1",5),rep("MCF7",5)),  hline = c(seq(from=1.5, to=5.5, by = 1),seq(from=1.5, to=5.5, by = 1),seq(from=1.5, to=5.5, by = 1)))

ggplot(allmuts_df, aes(x = Clone, y = Mutation)) +
  geom_count(aes(size = Mutation_norm,color=Mutation))+
  scale_size(breaks = c(0.01,0.025,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4),range = c(3,18)) +
  facet_grid(~Cell,space="free_x",scales = "free_x")+
  geom_vline(data = data_vline,aes(xintercept = vline))+
  geom_hline(data = data_hline,aes(yintercept = hline))+
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
  guides(size = guide_legend(title = "\n\nContribution to clone clustered mutation burden",nrow=1,title.position="top", title.hjust = 0.5))
ggsave("Fig6C.pdf",plot = last_plot(),device="pdf", width = 6500,height = 3000,units = "px")

ggplot(allmuts_df, aes(x = Clone, y = Mutation,fill=Difference)) +
  geom_raster(aes(x=Clone, y=Mutation, fill=Difference)) +
  facet_grid(~Cell,space="free_x",scales = "free_x")+
  geom_text(aes(Clone, Mutation, label = round(Difference,2)), color = "grey20", size = 6) +
  scale_fill_gradient2(low = "red",mid = "white",high = "blue",breaks=c(-0.3,-0.2,-0.1,0,0.1,0.2),labels=c(-0.3,-0.2,-0.1,0,0.1,0.2))+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        strip.background = element_blank(),
        text = element_text(size = 25),
        panel.background=element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(3, "lines"),
        legend.position = "right",
        legend.text = element_text(hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        strip.text.x = element_text(size = 25,face = "bold",hjust=0.5,vjust=1,margin = margin(0,0,1,0, "cm")))+
  xlab("\nClone")+
  ylab("Mutation type\n")
ggsave("SFig9B.pdf",plot = last_plot(),device="pdf", width = 7500,height = 3000,units = "px")



