### This script outputs indels from the Zou and Kucab studies in vcf form, so they can be later classified and added to our signature extraction pipeline

library("readxl")
library("MutationalPatterns")
library("xlsx")
library("plyr")
library("gtools")
library("stringi")
library("data.table")
library("dplyr")
library("stringr")

human_chrs = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22","chrX", "chrY")


############################################################################# ZOU #############################################################################
Zou = read.table("denovo_indels_43genes_Zou2021.txt", h=T, stringsAsFactors = F, sep="\t")
gene_KOs_with_sig=c("OGG1","UNG","EXO1","RNF168","MLH1","MSH2","MSH6","PMS2","PMS1") #genes marked as having signatures in Zou 2021

Zou = Zou[which(Zou$Ko_gene %in% gene_KOs_with_sig),] #now we remove all variants from other gene KOs not on our list
Zou_samples=unique(Zou$Sample)
Zou_samples=str_replace(mixedsort(str_replace(Zou_samples,"[.]","_")),"_",".")
Zou$Sample=factor(Zou$Sample,Zou_samples)
Zou=Zou[order(Zou$Sample),]
Zou$Sample = gsub("\\.", "_", Zou$Sample)
Zou_samples = gsub("\\.", "_", Zou_samples)

#change chromosome names and factorise
Zou$Chrom=mapvalues(Zou$Chrom, from=c(23,24), to=c("X","Y"))
Zou$Chrom=paste("chr",Zou$Chrom,sep="")
Zou$Chrom=factor(Zou$Chrom, levels = human_chrs)

setwd("Zou_vcfs")
for(sample in Zou_samples){ #loop through each sample and output an indel VCF
  df = Zou[which(Zou$Sample == sample), ]
  vcf = data.frame(CHROM = df$Chrom, POS = df$Pos, REF = df$Ref, ALT = df$Alt, ID="rsX", QUAL=".", FILTER="PASS")
  vcf = vcf[,c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")] # correct order for sigprofiler
  vcf$CHROM=factor(vcf$CHROM,levels=human_chrs)
  vcf=vcf[order(vcf$CHROM, vcf$POS),]
  colnames(vcf)[which(grepl("CHROM", colnames(vcf)))] = "#CHROM"
  write.table(c("##fileformat=VCFv4.2"), paste0(sample,".vcf"), row.names = F,col.names = F,quote = F)
  fwrite(vcf, paste0(sample,".vcf"), row.names = F,col.names = T,quote = F,append = T,sep="\t")}


############################################################################ KUCAB ############################################################################
Kucab = read.table("denovo_subclone_indels_final_Kucab2019.txt", h=T, stringsAsFactors = F, sep="\t")
Kucab_supp = read.xlsx("suppTable2_Kucab2019.xlsx", sheetIndex = 1)

#keep the tretaments with signature and the radiation treatment 
Kucab_treatments=c(as.character(Kucab_supp[which(Kucab_supp$Has.Signature=="YES"),"Sample.Name"]),"MSM0.71")  
Kucab_treatments=str_replace(mixedsort(str_replace(Kucab_treatments,"[.]","_")),"_",".")
Kucab=Kucab[which(Kucab$Sample.Name %in% Kucab_treatments),]

#some treatments have multiple clones, reorder the names 
Kucab_samples=unique(Kucab$Sample)
Kucab_samples=str_replace(mixedsort(str_replace(Kucab_samples,"[.]","_")),"_",".")
Kucab$Sample=factor(Kucab$Sample,Kucab_samples)
Kucab=Kucab[order(Kucab$Sample),]
Kucab$Sample = gsub("\\.", "_", Kucab$Sample)
Kucab_samples = gsub("\\.", "_", Kucab_samples)

#change chromosome names and factorise
#Kucab$Chrom=mapvalues(Kucab$Chrom, from=c(23,24), to=c("X","Y")) #no chromosomes 23 and 24 in the indel Kucab dataset
Kucab$Chrom=paste("chr",Kucab$Chrom,sep="")
Kucab$Chrom=factor(Kucab$Chrom, levels = human_chrs)

setwd("Kucab_vcfs")
for(sample in Kucab_samples){ #loop through each sample and output an indel VCF
  df = Kucab[which(Kucab$Sample == sample), ]
  vcf = data.frame(CHROM = df$Chrom, POS = df$Pos, REF = df$Ref, ALT = df$Alt, ID="rsX", QUAL=".", FILTER="PASS")
  vcf = vcf[,c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")] # correct order for sigprofiler
  vcf$CHROM=factor(vcf$CHROM,levels=human_chrs)
  vcf=vcf[order(vcf$CHROM, vcf$POS),]
  colnames(vcf)[which(grepl("CHROM", colnames(vcf)))] = "#CHROM"
  write.table(c("##fileformat=VCFv4.2"), paste0(sample,".vcf"), row.names = F,col.names = F,quote = F)
  fwrite(vcf, paste0(sample,".vcf"), row.names = F,col.names = T,quote = F,append = T,sep="\t")}
