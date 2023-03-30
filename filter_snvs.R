library(stringr)
library (dplyr)
library(data.table)
setwd("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2")
A549_CA=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/A549_CA_snvs.txt",sep="\t")
A549_CC=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/A549_CC_snvs.txt",sep="\t")
A549_HE1A=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/A549_HE1A_snvs.txt",sep="\t")
A549_HE2B=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/A549_HE2B_snvs.txt",sep="\t")
A549_PR1B=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/A549_PR1B_snvs.txt",sep="\t")
A549_PR2B=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/A549_PR2B_snvs.txt",sep="\t")
MCF7_C=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/MCF7_C_snvs.txt",sep="\t")
MCF7_HE1B=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/MCF7_HE1B_snvs.txt",sep="\t")
MCF7_HE1C=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/MCF7_HE1C_snvs.txt",sep="\t")
MCF7_HE2A=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/MCF7_HE2A_snvs.txt",sep="\t")
MCF7_HE2B=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/MCF7_HE2B_snvs.txt",sep="\t")
MCF7_PR1A=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/MCF7_PR1A_snvs.txt",sep="\t")
MCF7_PR1B=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/MCF7_PR1B_snvs.txt",sep="\t")
MCF7_PR2A=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/MCF7_PR2A_snvs.txt",sep="\t")
MCF7_PR2C=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/MCF7_PR2C_snvs.txt",sep="\t")
HAP1_CA=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/HAP1_CA_snvs.txt",sep="\t")
HAP1_CB=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/HAP1_CB_snvs.txt",sep="\t")
HAP1_HE1A=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/HAP1_HE1A_snvs.txt",sep="\t")
HAP1_HE2A=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/HAP1_HE2A_snvs.txt",sep="\t")
HAP1_PR1A=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/HAP1_PR1A_snvs.txt",sep="\t")
HAP1_PR2A=read.table("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/txts/HAP1_PR2A_snvs.txt",sep="\t")

samples=c("MCF7_C","MCF7_HE1B","MCF7_HE1C","MCF7_HE2A","MCF7_HE2B",
          "MCF7_PR1A","MCF7_PR1B","MCF7_PR2A","MCF7_PR2C",
          "A549_CA","A549_CC","A549_HE1A","A549_HE2B","A549_PR1B","A549_PR2B",
          "HAP1_CA","HAP1_CB","HAP1_HE1A","HAP1_HE2A","HAP1_PR1A","HAP1_PR2A")

for (sample in samples){
    cat("Working on sample: ",sample,"\n")
    cellline=str_split(sample,"_")[[1]][1] #get the cell line
    unique_snvs=get(sample) # get the appropriate txt file
    
    unique_snvs$V9=lengths(str_split(unique_snvs$V7,",")) #get the number of alleles at each position (we want to filter out positions with 2 or more alt alleles)
    unique_snvs$V10=sapply(str_split(unique_snvs$V7,","), "[[", 2) #get the number of supporting reads for the alt allele (we want to filter out variants with only 1 supporting read)
    unique_snvs$V11=as.numeric(unique_snvs$V10)/as.numeric(unique_snvs$V8) #calculate VAF for each snv

    unique_snvs=unique_snvs[which(unique_snvs$V9==2),]  #only keep variants with 2 alleles (1 alt allele)
    unique_snvs=unique_snvs[which(unique_snvs$V10>1),]  #only keep variants with more than 1 supporting read
    
    vcf=data.frame(CHROM = unique_snvs$V1, POS = unique_snvs$V2, REF = unique_snvs$V3, ALT = unique_snvs$V4, ID="rsX", QUAL=".", FILTER="PASS")  
    vcf = vcf[,c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")] # correct order for sigprofiler
    colnames(vcf)[which(grepl("CHROM", colnames(vcf)))] = "#CHROM"
    setwd("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/filtered_vcfs")
    write.table(c("##fileformat=VCFv4.2"), paste0(sample,".vcf"), row.names = F,col.names = F,quote = F)
    fwrite(vcf, paste0(sample,".vcf"), row.names = F,col.names = T,quote = F,append = T,sep="\t")
    
    setwd("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/Samples_strelka2/snvs/filtered_txts")
    txt=data.frame(CHROM = unique_snvs$V1, POS = unique_snvs$V2, REF = unique_snvs$V3, ALT = unique_snvs$V4,GT = unique_snvs$V5, GQ = unique_snvs$V6, ALTAD = unique_snvs$V10, ALTVAF = unique_snvs$V11,stringsAsFactors = F)
    write.table(txt, paste0(sample,".txt"), row.names = F,col.names = T,quote = F,sep="\t")
    
    txt2=data.frame(sample=rep(sample,times=nrow(txt)),chr=txt$CHROM,start=txt$POS,end=txt$POS,ref=txt$REF,alt=txt$ALT,stringsAsFactors = F)
    assign(paste0(sample,"_txt"),txt2)
}


setwd("/g/strcombio/fsupek_cancer1/radiation_project/radiation_project/Maia_analysis/REA_redone/TABLES_MERGED")
A549_C=rbind(A549_CA_txt,A549_CC_txt)
A549_H=rbind(A549_HE1A_txt,A549_HE2B_txt)
A549_P=rbind(A549_PR1B_txt,A549_PR2B_txt)
HAP1_C=rbind(HAP1_CA_txt,HAP1_CB_txt)
HAP1_H=rbind(HAP1_HE1A_txt,HAP1_HE2A_txt)
HAP1_P=rbind(HAP1_PR1A_txt,HAP1_PR2A_txt)
MCF7_C=rbind(MCF7_C_txt)
MCF7_H=rbind(MCF7_HE1B_txt,MCF7_HE1C_txt,MCF7_HE2A_txt,MCF7_HE2B_txt)
MCF7_P=rbind(MCF7_PR1A_txt,MCF7_PR1B_txt,MCF7_PR2A_txt,MCF7_PR2C_txt)

write.table(A549_C, file="A549_C.txt", sep="\t", quote=F, row.names = F)
write.table(A549_H, file="A549_H.txt", sep="\t", quote=F, row.names = F)
write.table(A549_P, file="A549_P.txt", sep="\t", quote=F, row.names = F)
write.table(HAP1_C, file="HAP1_C.txt", sep="\t", quote=F, row.names = F)
write.table(HAP1_H, file="HAP1_H.txt", sep="\t", quote=F, row.names = F)
write.table(HAP1_P, file="HAP1_P.txt", sep="\t", quote=F, row.names = F)
write.table(MCF7_C, file="MCF7_C.txt", sep="\t", quote=F, row.names = F)
write.table(MCF7_H, file="MCF7_H.txt", sep="\t", quote=F, row.names = F)
write.table(MCF7_P, file="MCF7_P.txt", sep="\t", quote=F, row.names = F)


