### This script outputs snvs from the Zou and Kucab studies in mutation matrix form, so they can be later added to our signature extraction pipeline

library("readxl")
library("MutationalPatterns")
library("xlsx")
library("plyr")
library("gtools")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg19")
library("stringi")
library("data.table")
library("dplyr")
library("stringr")

human_chrs = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22","chrX", "chrY")
COSMIC_SBS=fread("COSMIC_v3.3.1_SBS_GRCh38.txt")
dnaRevCompl <- function(nucSeq){return(stri_reverse(chartr("acgtACGT", "tgcaTGCA", nucSeq)))}

############################################################################# OUR #############################################################################
Our=fread("Filtered_SNVs/output/SBS/snvs.SBS96.all")
rownames(Our)=COSMIC_SBS$Type
Our=Our[,-1]


############################################################################ KUCAB ############################################################################
#reading the mutation data and sample information
Kucab = read.table("denovo_subclone_subs_final_Kucab2019.txt", h=T, stringsAsFactors = F, sep="\t")
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

#change chromosome names and factorise
Kucab$Chrom=mapvalues(Kucab$Chrom, from=c(23,24), to=c("X","Y"))
Kucab$Chrom=paste("chr",Kucab$Chrom,sep="")
Kucab$Chrom=factor(Kucab$Chrom, levels = human_chrs)

#Get upstream and downstream bases for each variants
Kucab$Upstream=as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19,Kucab$Chrom,Kucab$Pos-1,Kucab$Pos-1))  #get upstream base
Kucab$Downstream=as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19,Kucab$Chrom,Kucab$Pos+1,Kucab$Pos+1)) 
Kucab=Kucab[which(Kucab$Upstream %in% c("A","C","G","T") & Kucab$Downstream %in% c("A","C","G","T")),]  #check all upstream and downstream bases are ok

#now classify the mutations based on context (always expressed as a pyrimidine mutation)
Kucab$Mutation=ifelse(Kucab$Ref %in% c("C","T"), paste(Kucab$Upstream,"[",Kucab$Ref,">",Kucab$Alt,"]",Kucab$Downstream,sep=""), paste(dnaRevCompl(Kucab$Downstream),"[",dnaRevCompl(Kucab$Ref),">",dnaRevCompl(Kucab$Alt),"]",dnaRevCompl(Kucab$Upstream),sep=""))

#create an empty data frame for the Kucab mutations
Kucab_df=data.frame(matrix(0, ncol=length(Kucab_samples), nrow=96))
rownames(Kucab_df)=COSMIC_SBS$Type
colnames(Kucab_df)=Kucab_samples

#look through every sample in Kucab and add mutations to the df
for (i in 1:length(Kucab_samples)){ 
  sample=Kucab_samples[i]
  df=Kucab[which(Kucab$Sample==sample),]
  df %>% group_by(Mutation) %>% dplyr::summarise(Count = n()) -> counts
  Kucab_df[counts$Mutation,i]=counts$Count}


############################################################################# ZOU #############################################################################
#reading the mutation data
Zou = read.table("denovo_subs_43genes_Zou2021.txt", h=T, stringsAsFactors = F, sep="\t")
gene_KOs_with_sig=c("OGG1","UNG","EXO1","RNF168","MLH1","MSH2","MSH6","PMS2","PMS1") #genes marked as having signatures in Zou 2021

Zou = Zou[which(Zou$Ko_gene %in% gene_KOs_with_sig),] #now we remove all variants from other gene KOs not on our list
Zou_samples=unique(Zou$Sample)
Zou_samples=str_replace(mixedsort(str_replace(Zou_samples,"[.]","_")),"_",".")
Zou$Sample=factor(Zou$Sample,Zou_samples)
Zou=Zou[order(Zou$Sample),]

#change chromosome names and factorise
Zou$Chrom=mapvalues(Zou$Chrom, from=c(23,24), to=c("X","Y"))
Zou$Chrom=paste("chr",Zou$Chrom,sep="")
Zou$Chrom=factor(Zou$Chrom, levels = human_chrs)

#Get upstream and downstream bases for each variants
Zou$Upstream=as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19,Zou$Chrom,Zou$Pos-1,Zou$Pos-1))  #get upstream base
Zou$Downstream=as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19,Zou$Chrom,Zou$Pos+1,Zou$Pos+1))  #get downstream base
Zou=Zou[which(Zou$Upstream %in% c("A","C","G","T") & Zou$Downstream %in% c("A","C","G","T")),]  #check all upstream and downstream bases are ok

#get the mutation with context making sure to express it in terms of pyrimidine 
Zou$Mutation=ifelse(Zou$Ref %in% c("C","T"), paste(Zou$Upstream,"[",Zou$Ref,">",Zou$Alt,"]",Zou$Downstream,sep=""), paste(dnaRevCompl(Zou$Downstream),"[",dnaRevCompl(Zou$Ref),">",dnaRevCompl(Zou$Alt),"]",dnaRevCompl(Zou$Upstream),sep=""))

#create an empty data frame for the Zou mutations
Zou_df=data.frame(matrix(0, ncol=length(Zou_samples), nrow=96))
rownames(Zou_df)=COSMIC_SBS$Type
colnames(Zou_df)=Zou_samples

#look through every sample in Zou and add mutations to the df
for (i in 1:length(Zou_samples)){ 
  sample=Zou_samples[i]
  df=Zou[which(Zou$Sample==sample),]
  df %>% group_by(Mutation) %>% dplyr::summarise(Count = n()) -> counts
  Zou_df[counts$Mutation,i]=counts$Count}


############################################################################ MERGE ############################################################################
nmf_samples = cbind(Our, Zou_df, Kucab_df)
nmf_samples = cbind(data.frame(MutationType = COSMIC_SBS$Type), nmf_samples)
write.table(nmf_samples, "Our_Zou_Kucab_SNVs.txt", append=F, sep="\t", row.names = FALSE, col.names = TRUE, quote=F)


#now randomly removing half of all clones from both Zou and Kucab, so we can run NMF on a smaller set of samples
set.seed(1)
Kucab_sub1=t(sample_n(as.data.frame(t(Kucab_df)),77)); Kucab_sub2=t(sample_n(as.data.frame(t(Kucab_df)),77)); Kucab_sub3=t(sample_n(as.data.frame(t(Kucab_df)),77))
Kucab_sub1=Kucab_sub1[,Kucab_samples[Kucab_samples %in% colnames(Kucab_sub1)]]
Kucab_sub2=Kucab_sub2[,Kucab_samples[Kucab_samples %in% colnames(Kucab_sub2)]]
Kucab_sub3=Kucab_sub3[,Kucab_samples[Kucab_samples %in% colnames(Kucab_sub3)]]

set.seed(1)
Zou_sub1=t(sample_n(as.data.frame(t(Zou_df)),19)); Zou_sub2=t(sample_n(as.data.frame(t(Zou_df)),19)); Zou_sub3=t(sample_n(as.data.frame(t(Zou_df)),19))
Zou_sub1=Zou_sub1[,Zou_samples[Zou_samples %in% colnames(Zou_sub1)]]
Zou_sub2=Zou_sub2[,Zou_samples[Zou_samples %in% colnames(Zou_sub2)]]
Zou_sub3=Zou_sub3[,Zou_samples[Zou_samples %in% colnames(Zou_sub3)]]

#join all
nmf_samples_sub1 = cbind(Our, Zou_sub1, Kucab_sub1); nmf_samples_sub2 = cbind(Our, Zou_sub2, Kucab_sub2); nmf_samples_sub3 = cbind(Our, Zou_sub3, Kucab_sub3)
nmf_samples_sub1 = cbind(data.frame(MutationType = COSMIC_SBS$Type), nmf_samples_sub1)
nmf_samples_sub2 = cbind(data.frame(MutationType = COSMIC_SBS$Type), nmf_samples_sub2)
nmf_samples_sub3 = cbind(data.frame(MutationType = COSMIC_SBS$Type), nmf_samples_sub3)
write.table(nmf_samples_sub1, "Our_Zou_Kucab_SNVs_sub1.txt", append=F, sep="\t", row.names = FALSE, col.names = TRUE, quote=F)
write.table(nmf_samples_sub2, "Our_Zou_Kucab_SNVs_sub2.txt", append=F, sep="\t", row.names = FALSE, col.names = TRUE, quote=F)
write.table(nmf_samples_sub3, "Our_Zou_Kucab_SNVs_sub3.txt", append=F, sep="\t", row.names = FALSE, col.names = TRUE, quote=F)
