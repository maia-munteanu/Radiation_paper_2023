##### MANTA STRUCTURAL VARIANTS #####
library(VariantAnnotation)
library(gridExtra)
library(ggplot2)
library(ggpubr)

A549 = "A549_diploidSV_PASS.vcf.gz"
HAP1 = "HAP1_diploidSV_PASS.vcf.gz"
MCF7 = "MCF7_diploidSV_PASS.vcf.gz"

A549_inv = "A549_diploidSV_PASS_convertInvertion.vcf.gz"
HAP1_inv = "HAP1_diploidSV_PASS_convertInvertion.vcf.gz"
MCF7_inv = "MCF7_diploidSV_PASS_convertInvertion_with_control.vcf.gz"

nb_positive_samples = function(id, matrix_geno, matrix_filter){
  genos = matrix_geno[id, ]
  filts = matrix_filter[id, ]
  genos_kept = genos #[which(filts == "PASS")] #add this will keep a lot of germline with VAF~1
  sum ( genos_kept!="0/0" & genos_kept!="0|0" ) }

for(cellline in c("A549", "HAP1", "MCF7")){
  print(paste(date(), " INFO: working on cell line: ", cellline, sep=""))
  vcf = open(VcfFile(get(cellline),  yieldSize=100000))
  vcf_chunk0 = readVcf(vcf, "hg38")
  vcf_chunk = vcf_chunk0[which(!grepl("BND", rownames(vcf_chunk0))) , ] # remove BND because we run the convertInversions.py
  vcf_chunk = vcf_chunk[which(is.na(info(vcf_chunk)$EVENT) | !duplicated(info(vcf_chunk)$EVENT)) , ] # same event for non-inv == one SV

  vcf_inv = open(VcfFile(get(paste(cellline, "_inv", sep="")),  yieldSize=100000))
  vcf_chunk_inv0 = readVcf(vcf_inv, "hg38")
  vcf_chunk_inv = vcf_chunk_inv0[which(grepl("INV", rownames(vcf_chunk_inv0))) , ] # add INV after removing BND
  events_ids = info(vcf_chunk_inv)$EVENT # compute the events IDs-- unique or NA are single inverted junctions
  vcf_chunk_inv = vcf_chunk_inv[ which( is.na(events_ids) | !events_ids %in% events_ids[duplicated(events_ids)]), ] # remove reciprocal inversions (should be artefacts)

  # merge SV without BND and the inversions
  gg0 = rbind(geno(vcf_chunk, "GT"), geno(vcf_chunk_inv, "GT"))
  ft0 = rbind(geno(vcf_chunk, "FT"), geno(vcf_chunk_inv, "FT"))
  
  pr0 = rbind(geno(vcf_chunk, "PR"), geno(vcf_chunk_inv, "PR")) 
  sr0 = rbind(geno(vcf_chunk, "SR"), geno(vcf_chunk_inv, "SR"))
  
  type0 = c(info(vcf_chunk)[,"SVTYPE"], info(vcf_chunk_inv)[,"SVTYPE"])
  len0 = c(info(vcf_chunk)[,"SVLEN"], info(vcf_chunk_inv)[,"SVLEN"])
  chr0=c(as.character(seqnames(rowRanges(vcf_chunk,"seqnames"))), as.character(seqnames(rowRanges(vcf_chunk_inv,"seqnames"))))
  loc0=c(start(ranges(rowRanges(vcf_chunk,"seqnames"))), start(ranges(rowRanges(vcf_chunk_inv,"seqnames"))))
  ref0=rep("A", length(loc0))
  alts0=ref0 
  
  ids_keep = which( unlist(lapply(1:nrow(gg0), function(id) nb_positive_samples(id, gg0, ft0))) < 2 )
  gg = gg0[ids_keep, ]
  ft = ft0[ids_keep, ]
  pr = pr0[ids_keep, ] 
  sr = sr0[ids_keep, ] 
  len = len0[ids_keep, ]
  muts_id = paste(chr0[ids_keep], ":", loc0[ids_keep], "_", ref0[ids_keep], "/", alts0[ids_keep], sep="")
  is.na(len) <- lengths(len) == 0 # replace integer(0) by NAs, otherwise the unlist() removes them
  len = abs(unlist(len))
  len[which(is.na(len))] = 1000 # large insertions more than 1000 are length NA in the VCF (CAUTION: BND are also NA)
  
  allsm = colnames(gg0)
  allsm = unlist(lapply(allsm, function(s) gsub("_REP1","",s)))
  for(sm in allsm){
    ggtmp = gg[, grepl(sm, colnames(gg))] ; fttmp = ft[, grepl(sm, colnames(ft))]
    prtmp = pr[, grepl(sm, colnames(pr))]; prtmp = sapply(prtmp, "[[", 2) 
    srtmp = sr[, grepl(sm, colnames(sr))]; srtmp_names=c(); srtmp_values=c()
    for (i in 1:length(srtmp)){
         element=srtmp[[i]]
         if(identical(element, integer(0))){value=0}else{value=element[2]}
         srtmp_names=c(srtmp_names,names(srtmp)[i])
         srtmp_values=c(srtmp_values,value)}
    srtmp=srtmp_values; names(srtmp)=srtmp_names 

    muts_id_sm = muts_id[which(ggtmp!="0/0" & ggtmp!="." & ggtmp!="0|0" & fttmp == "PASS" & srtmp>=2)] #new filters where we also require at least 2 split reads
    mutSV = names(ggtmp)[which(ggtmp!="0/0" & ggtmp!="." & ggtmp!="0|0" & fttmp == "PASS" & srtmp>=2)] 

    kept_ids_sv = 1:length(mutSV) # because we do not look at SV we take all
    lenSV = len[which(ggtmp!="0/0" & ggtmp!="." & ggtmp!="0|0" & fttmp == "PASS")][kept_ids_sv] # length of SVs where the sample is mutated and PASS and remove multiple BND
    mutSV = mutSV[kept_ids_sv] # we remove multiple BND
    muts_id_sm = muts_id_sm[kept_ids_sv] # we remove multiple BND
    assign(paste("all_SVs",sm,sep="_"), muts_id_sm) ; save(list=paste("all_SVs",sm,sep = "_"), file=paste("Rdata/all_SVs_", sm, ".Rdata", sep = ""))
    
    
    mutSVtype = gsub("Manta", "", unlist(lapply(mutSV, function(x) unlist(strsplit(x, ":"))[1]))) # type of SVs
    d = data.frame(clone = sm, nb_SV = as.numeric(table(mutSVtype)), type = names(table(mutSVtype)))
    dlen = data.frame(clone = sm, SV_length = log10(lenSV), type = mutSVtype)
    if(sm == allsm[1]) {finald = d} else {finald = rbind(finald, d)}
    if(sm == allsm[1]) {finaldlen = dlen} else {finaldlen = rbind(finaldlen, dlen)}
  }
  # consider BND as INV because we took them unique
  finald$type = as.character(finald$type) # do not consider it as factor
  finald[which(finald$type == "BND"), "type"] = "INV"
  # assign data frames
  assign(paste("finald_", cellline, sep=""), finald)
  assign(paste("finaldlen_", cellline, sep=""), finaldlen)
  # assign plot of number of mutation per type to each sample
  assign(paste("p_", cellline, sep=""),
         ggbarplot(finald, x="clone", y="nb_SV", fill="type", label=T, palette = "jco") + ggtitle(cellline) +  rotate_x_text(45)
  )
  # assign plot of length of mutation per type to each sample
  assign(paste("p_length_", cellline, sep=""),
         ggviolin(finaldlen, x="type", y="SV_length", fill="clone", add = "boxplot") + ggtitle(cellline))
}

SV_vector=c("all_SVs_A549_CA","all_SVs_A549_CC","all_SVs_A549_HE1A","all_SVs_A549_HE2B","all_SVs_A549_PR1B","all_SVs_A549_PR2B",
            "all_SVs_HAP1_CA","all_SVs_HAP1_CB","all_SVs_HAP1_HE1A","all_SVs_HAP1_HE2A","all_SVs_HAP1_PR1A","all_SVs_HAP1_PR2A",
            "all_SVs_MCF7_untreated_F8_REP2","all_SVs_MCF7_HE_1B","all_SVs_MCF7_HE_1C","all_SVs_MCF7_HE_2A","all_SVs_MCF7_HE_2B","all_SVs_MCF7_PR_1A","all_SVs_MCF7_PR_1B","all_SVs_MCF7_PR_2A","all_SVs_MCF7_PR_2C" )

for (SV in SV_vector){
  muts = get(SV)
  chr = unlist(lapply(muts, function(x) unlist(strsplit(x, ":"))[1]))
  start = as.numeric(unlist(lapply(muts, function(x) unlist(strsplit(unlist(strsplit(x, ":"))[2], "_"))[1])))
  end = as.numeric(unlist(lapply(muts, function(x) unlist(strsplit(unlist(strsplit(x, ":"))[2], "_"))[1])))
  ref = unlist(lapply(muts, function(x) unlist(strsplit(unlist(strsplit(unlist(strsplit(x, ":"))[2], "_"))[2], "/"))[1] ))
  alt = unlist(lapply(muts, function(x) unlist(strsplit(unlist(strsplit(unlist(strsplit(x, ":"))[2], "_"))[2], "/"))[2] ))
  dd = data.frame(sample = SV, chr = chr, start = start, end = end, ref = ref, alt = alt)
  assign(paste0(SV,"_txt"),dd)}

all_SVs_A549_CA_txt$sample=rep("A549_CA",times=nrow(all_SVs_A549_CA_txt))
all_SVs_A549_CC_txt$sample=rep("A549_CC",times=nrow(all_SVs_A549_CC_txt))
all_SVs_A549_HE1A_txt$sample=rep("A549_HE1A",times=nrow(all_SVs_A549_HE1A_txt))
all_SVs_A549_HE2B_txt$sample=rep("A549_HE2B",times=nrow(all_SVs_A549_HE2B_txt))
all_SVs_A549_PR1B_txt$sample=rep("A549_PR1B",times=nrow(all_SVs_A549_PR1B_txt))
all_SVs_A549_PR2B_txt$sample=rep("A549_PR2B",times=nrow(all_SVs_A549_PR2B_txt))
all_SVs_HAP1_CA_txt$sample=rep("HAP1_CA",times=nrow(all_SVs_HAP1_CA_txt))
all_SVs_HAP1_CB_txt$sample=rep("HAP1_CB",times=nrow(all_SVs_HAP1_CB_txt))
all_SVs_HAP1_HE1A_txt$sample=rep("HAP1_HE1A",times=nrow(all_SVs_HAP1_HE1A_txt))
all_SVs_HAP1_HE2A_txt$sample=rep("HAP1_HE2A",times=nrow(all_SVs_HAP1_HE2A_txt))
all_SVs_HAP1_PR1A_txt$sample=rep("HAP1_PR1A",times=nrow(all_SVs_HAP1_PR1A_txt))
all_SVs_HAP1_PR2A_txt$sample=rep("HAP1_PR2A",times=nrow(all_SVs_HAP1_PR2A_txt))
all_SVs_MCF7_untreated_F8_REP2_txt$sample=rep("MCF7_C",times=nrow(all_SVs_MCF7_untreated_F8_REP2_txt))
all_SVs_MCF7_HE_1B_txt$sample=rep("MCF7_HE1B",times=nrow(all_SVs_MCF7_HE_1B_txt))
all_SVs_MCF7_HE_1C_txt$sample=rep("MCF7_HE1C",times=nrow(all_SVs_MCF7_HE_1C_txt))
all_SVs_MCF7_HE_2A_txt$sample=rep("MCF7_HE2A",times=nrow(all_SVs_MCF7_HE_2A_txt))
all_SVs_MCF7_HE_2B_txt$sample=rep("MCF7_HE2B",times=nrow(all_SVs_MCF7_HE_2B_txt))
all_SVs_MCF7_PR_1A_txt$sample=rep("MCF7_PR1A",times=nrow(all_SVs_MCF7_PR_1A_txt))
all_SVs_MCF7_PR_1B_txt$sample=rep("MCF7_PR1B",times=nrow(all_SVs_MCF7_PR_1B_txt))
all_SVs_MCF7_PR_2A_txt$sample=rep("MCF7_PR2A",times=nrow(all_SVs_MCF7_PR_2A_txt))
all_SVs_MCF7_PR_2C_txt$sample=rep("MCF7_PR2C",times=nrow(all_SVs_MCF7_PR_2C_txt))

finald_MCF7$clone=mapvalues(finald_MCF7$clone,from=c("MCF7_untreated_F8_REP2","MCF7_HE_1B","MCF7_HE_1C","MCF7_HE_2A","MCF7_HE_2B","MCF7_PR_1A","MCF7_PR_1B","MCF7_PR_2A","MCF7_PR_2C"),to=c("MCF7_C","MCF7_HE1B","MCF7_HE1C","MCF7_HE2A","MCF7_HE2B","MCF7_PR1A","MCF7_PR1B","MCF7_PR2A","MCF7_PR2C"))
finald=rbind(finald_A549,finald_HAP1,finald_MCF7)
colnames(finald)=c("Clone","Count","SV")
finald$Cell=sapply(str_split(finald$Clone,"_"), "[[", 1)
finald$Samples=sapply(str_split(finald$Clone,"_"), "[[", 2)
finald$Samples=mapvalues(finald$Samples,from=c("CA","CC","HE1A","HE2B","PR1B","PR2B","CB","HE2A","PR1A","PR2A","C","HE1B","HE1C","PR2C"),to=c("Control_A","Control_C","Helium_1A","Helium_2B","Protons_1B","Protons_2B","Control_B","Helium_2A","Protons_1A","Protons_2A","Control_A","Helium_1B","Helium_1C","Protons_2C"))
finald$Treatment=sapply(str_split(finald$Samples,"_"), "[[", 1)
finald$Cell=factor(finald$Cell,levels=c("HAP1","A549","MCF7"))
finald$SV=factor(finald$SV,levels=c("DEL","DUP","INS","INV"))
###finald=finald[which(finald$Clone!="MCF7_HE2B"),] #remove the HE2B clone 

write.table(finald, file = "SV_counts.txt", sep = " ",col.names = TRUE,quote=F)

ggplot(finald, aes(fill=SV, y=Count, x=Samples)) +
  geom_bar(position="stack", stat="identity",colour="black")+
  geom_text(aes(label=Count), position = position_stack(vjust = 0.5),size=5.5)+
  facet_wrap(~Cell, nrow=1,scales = "free_x")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "bottom",
        text = element_text(size = 20),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"))+
  ylab("Number of structural variants\n")+
  xlab("\nClone")+
  scale_fill_manual(values=c("#0073C2","#EFC000","#868686","#CD534C")) +
  labs(fill="Structural variant type")
ggsave("Fig5A.pdf",plot = last_plot(),device="pdf", width = 5000,height = 2800,units = "px")

ratios=data.frame(Treatment=c(),Sample=c(),Ratio=c(),Cell=c())
for (cell in c("HAP1","A549","MCF7")){
    for (sample in unique(finald[which(finald$Cell==cell),"Samples"])){
         del=finald[which(finald$Cell==cell & finald$Samples==sample & finald$SV=="DEL"),"Count"]
         ins_dup=sum(finald[which(finald$Cell==cell & finald$Samples==sample & (finald$SV=="INS" | finald$SV=="DUP")),"Count"])
         df=data.frame(Treatment=finald[which(finald$Cell==cell & finald$Samples==sample & finald$SV=="DEL"),"Treatment"],Sample=sample, Ratio=del/ins_dup,Cell=cell)
         ratios=rbind(ratios,df)}}

ggplot(ratios, aes(x=Treatment, y=Ratio)) +
  geom_boxplot(fill=NA, colour="grey20",outlier.shape = NA)+
  geom_jitter(size=3, shape=21, colour="grey20",fill="grey20",width = 0.2)+
  facet_grid(~Cell)+
  scale_y_continuous(breaks = seq(0, 5, by = 1),limits=c(0,5))+
  theme(axis.text.x = element_text(angle = 90, hjust=0),
        text = element_text(size = 25),
        strip.background = element_rect(fill="white"),
        strip.text.x = element_text(size=30,face="bold"),
        panel.background=element_rect(fill="white"),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"))+
  ylab("Deletions / (Insertions + Duplications) ratio\n")+
  xlab("\nTreatment")
ggsave("Fig5B.pdf",plot = last_plot(),device="pdf", width = 5500,height = 2500,units = "px")

finaldlen_MCF7$clone=mapvalues(finaldlen_MCF7$clone,from=c("MCF7_untreated_F8_REP2","MCF7_HE_1B","MCF7_HE_1C","MCF7_HE_2A","MCF7_HE_2B","MCF7_PR_1A","MCF7_PR_1B","MCF7_PR_2A","MCF7_PR_2C"),to=c("MCF7_C","MCF7_HE1B","MCF7_HE1C","MCF7_HE2A","MCF7_HE2B","MCF7_PR1A","MCF7_PR1B","MCF7_PR2A","MCF7_PR2C"))
finaldlen=rbind(finaldlen_A549,finaldlen_HAP1,finaldlen_MCF7)
colnames(finaldlen)=c("Clone","Length","SV")
finaldlen$Cell=sapply(str_split(finaldlen$Clone,"_"), "[[", 1)
finaldlen$Samples=sapply(str_split(finaldlen$Clone,"_"), "[[", 2)
finaldlen$Samples=mapvalues(finaldlen$Samples,from=c("CA","CC","HE1A","HE2B","PR1B","PR2B","CB","HE2A","PR1A","PR2A","C","HE1B","HE1C","PR2C"),to=c("Control_A","Control_C","Helium_1A","Helium_2B","Protons_1B","Protons_2B","Control_B","Helium_2A","Protons_1A","Protons_2A","Control_A","Helium_1B","Helium_1C","Protons_2C"))
finaldlen$Treatment=sapply(str_split(finaldlen$Samples,"_"), "[[", 1)
finaldlen$Cell=factor(finaldlen$Cell,levels=c("HAP1","A549","MCF7"))
finaldlen$SV=factor(finaldlen$SV,levels=c("DEL","INS","DUP","INV"))
###finaldlen=finaldlen[which(finaldlen$Clone!="MCF7_HE2B"),] #remove the HE2B clone

mycomps=list(c("Control", "Helium"),c("Helium", "Protons"),c("Control", "Protons"))
ggplot(finaldlen, aes(y=Length, x=Treatment)) +
  geom_boxplot(aes(color = Treatment))+
  geom_point(aes(color = Treatment))+
  facet_grid(Cell~SV)+
  theme_bw()+
  stat_compare_means(comparisons = mycomps,step.increase = 0.2,bracket.size = 0.55, size=4.2) +
  ylim(0,12)+
  scale_fill_manual(values=c("#0073C2","#EFC000","#868686"))+
  scale_color_manual(values=c("#0073C2","#EFC000","#868686"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 20),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"))+
  ylab("Structural variant length (log10) \n")

ggsave("SFig7.pdf",plot = last_plot(),device="pdf", width = 4600,height = 3700,units = "px")


