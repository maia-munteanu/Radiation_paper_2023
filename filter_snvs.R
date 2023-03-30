library(stringr)
library (dplyr)
library(data.table)
A549_CA=read.table("A549_CA_snvs.txt",sep="\t")
A549_CC=read.table("A549_CC_snvs.txt",sep="\t")
A549_HE1A=read.table("A549_HE1A_snvs.txt",sep="\t")
A549_HE2B=read.table("A549_HE2B_snvs.txt",sep="\t")
A549_PR1B=read.table("A549_PR1B_snvs.txt",sep="\t")
A549_PR2B=read.table("A549_PR2B_snvs.txt",sep="\t")
MCF7_C=read.table("MCF7_C_snvs.txt",sep="\t")
MCF7_HE1B=read.table("MCF7_HE1B_snvs.txt",sep="\t")
MCF7_HE1C=read.table("MCF7_HE1C_snvs.txt",sep="\t")
MCF7_HE2A=read.table("MCF7_HE2A_snvs.txt",sep="\t")
MCF7_HE2B=read.table("MCF7_HE2B_snvs.txt",sep="\t")
MCF7_PR1A=read.table("MCF7_PR1A_snvs.txt",sep="\t")
MCF7_PR1B=read.table("MCF7_PR1B_snvs.txt",sep="\t")
MCF7_PR2A=read.table("MCF7_PR2A_snvs.txt",sep="\t")
MCF7_PR2C=read.table("MCF7_PR2C_snvs.txt",sep="\t")
HAP1_CA=read.table("HAP1_CA_snvs.txt",sep="\t")
HAP1_CB=read.table("HAP1_CB_snvs.txt",sep="\t")
HAP1_HE1A=read.table("HAP1_HE1A_snvs.txt",sep="\t")
HAP1_HE2A=read.table("HAP1_HE2A_snvs.txt",sep="\t")
HAP1_PR1A=read.table("HAP1_PR1A_snvs.txt",sep="\t")
HAP1_PR2A=read.table("HAP1_PR2A_snvs.txt",sep="\t")

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
    setwd("Filtered_SNVs")
    write.table(c("##fileformat=VCFv4.2"), paste0(sample,".vcf"), row.names = F,col.names = F,quote = F)
    fwrite(vcf, paste0(sample,".vcf"), row.names = F,col.names = T,quote = F,append = T,sep="\t")
}




