library(stringr)
library(dplyr)
library(data.table)

A549_CA=read.table("A549_CA_indels.txt",sep="\t")
A549_CC=read.table("A549_CC_indels.txt",sep="\t")
A549_HE1A=read.table("A549_HE1A_indels.txt",sep="\t")
A549_HE2B=read.table("A549_HE2B_indels.txt",sep="\t")
A549_PR1B=read.table("A549_PR1B_indels.txt",sep="\t")
A549_PR2B=read.table("A549_PR2B_indels.txt",sep="\t")
MCF7_C=read.table("MCF7_C_indels.txt",sep="\t")
MCF7_HE1B=read.table("MCF7_HE1B_indels.txt",sep="\t")
MCF7_HE1C=read.table("MCF7_HE1C_indels.txt",sep="\t")
MCF7_HE2A=read.table("MCF7_HE2A_indels.txt",sep="\t")
MCF7_HE2B=read.table("MCF7_HE2B_indels.txt",sep="\t")
MCF7_PR1A=read.table("MCF7_PR1A_indels.txt",sep="\t")
MCF7_PR1B=read.table("MCF7_PR1B_indels.txt",sep="\t")
MCF7_PR2A=read.table("MCF7_PR2A_indels.txt",sep="\t")
MCF7_PR2C=read.table("MCF7_PR2C_indels.txt",sep="\t")
HAP1_CA=read.table("HAP1_CA_indels.txt",sep="\t")
HAP1_CB=read.table("HAP1_CB_indels.txt",sep="\t")
HAP1_HE1A=read.table("HAP1_HE1A_indels.txt",sep="\t")
HAP1_HE2A=read.table("HAP1_HE2A_indels.txt",sep="\t")
HAP1_PR1A=read.table("HAP1_PR1A_indels.txt",sep="\t")
HAP1_PR2A=read.table("HAP1_PR2A_indels.txt",sep="\t")

samples=c("A549_CA","A549_CC","A549_HE1A","A549_HE2B","A549_PR1B","A549_PR2B",
          "MCF7_C","MCF7_HE1B","MCF7_HE1C","MCF7_HE2A","MCF7_HE2B",
          "MCF7_PR1A","MCF7_PR1B","MCF7_PR2A","MCF7_PR2C",
          "HAP1_CA","HAP1_CB","HAP1_HE1A","HAP1_HE2A","HAP1_PR1A","HAP1_PR2A")

samples_nr=c("A549"=6,"HAP1"=6,"MCF7"=9)


for (sample in samples){
    cat("Working on sample: ",sample,"\n")
    cellline=str_split(sample,"_")[[1]][1] #get the cell line
    unique_indels=get(sample) # get the appropriate txt file
    
    unique_indels$V9=lengths(str_split(unique_indels$V7,",")) #get the number of alleles at each position (we want to filter out positions with 2 or more alt alleles)
    unique_indels$V10=sapply(str_split(unique_indels$V7,","), "[[", 2) #get the number of supporting reads for the alt allele (we want to filter out variants with only 1 supporting read)
    unique_indels$V11=as.numeric(unique_indels$V10)/as.numeric(unique_indels$V8) #calculate VAF for each indel

    unique_indels=unique_indels[which(unique_indels$V6>=10),] #filter by GQ, minimum value of 10
    unique_indels=unique_indels[which(unique_indels$V9==2),]  #only keep variants with 2 alleles (1 alt allele)
    unique_indels=unique_indels[which(unique_indels$V10>1),]  #only keep variants with more than 1 supporting read
    

    #now we loop through each variant to check it's read support in other samples
    all_samples=get(paste0(cellline,"_all")) #get the correct file for each cell line
    all_samples=all_samples[which(all_samples$V2 %in% unique_indels$V2),] #reduced the number of rows to speed this
    AD_end=4+2*as.numeric(samples_nr[cellline])
    AD_start=AD_end-as.numeric(samples_nr[cellline])+1
    excluded_vector=c()
    for (i in 1:nrow(unique_indels)){
      print(i)
      my_variant=filter(all_samples, as.character(V1) == as.character(unique_indels$V1[i]) & as.character(V2) == as.character(unique_indels$V2[i]) & as.character(V3) == as.character(unique_indels$V3[i]) & as.character(V4) == as.character(unique_indels$V4[i])) #get our variant from the main data frame so we have information about genotype and AD for the other samples as well
      variant_altAD=c(); for (j in AD_start:AD_end){variant_altAD=c(variant_altAD,sapply(str_split(my_variant[,j],","), "[[", 2))} #get the number of samples for which have 2 or more supporting reads for the alt allele (we only have 1 alt because we filtered out multiallelic loci)
      if(sum(variant_altAD>=2, na.rm=TRUE)>2 | sum(variant_altAD>=1, na.rm=TRUE)>3){excluded_vector=c(excluded_vector,i)}} #now add variants which have supporting reads that pass our filter in other samples
    
    unique_indels=unique_indels[-excluded_vector,] #remove variants with supporting reads
    vcf=data.frame(CHROM = unique_indels$V1, POS = unique_indels$V2, REF = unique_indels$V3, ALT = unique_indels$V4, ID="rsX", QUAL=".", FILTER="PASS")  
    vcf = vcf[,c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")] # correct order for sigprofiler
    colnames(vcf)[which(grepl("CHROM", colnames(vcf)))] = "#CHROM"
    setwd("Filtered_INDELs")
    write.table(c("##fileformat=VCFv4.2"), paste0(sample,".vcf"), row.names = F,col.names = F,quote = F)
    fwrite(vcf, paste0(sample,".vcf"), row.names = F,col.names = T,quote = F,append = T,sep="\t")
}
