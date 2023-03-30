library(reticulate)
library(SigProfilerMatrixGeneratorR)

inputfolder = "Filtered_SNVs"
matrices <- SigProfilerMatrixGeneratorR("snvs", "GRCh38", inputfolder, plot=T, exome=F,
                                       bed_file=NULL, chrom_based=F, tsb_stat=F, seqInfo=F, cushion=100)

inputfolder = "Clustered_SNVs"
matrices <- SigProfilerMatrixGeneratorR("snvs", "GRCh38", inputfolder, plot=T, exome=F,
                                        bed_file=NULL, chrom_based=F, tsb_stat=F, seqInfo=F, cushion=100)

inputfolder = "Filtered_INDELs"
matrices <- SigProfilerMatrixGeneratorR("indels", "GRCh38", inputfolder, plot=T, exome=F,
                                        bed_file=NULL, chrom_based=F, tsb_stat=F, seqInfo=F, cushion=100)

inputfolder = "Zou_vcfs"
matrices <- SigProfilerMatrixGeneratorR("indels", "GRCh37", inputfolder, plot=T, exome=F,
                                        bed_file=NULL, chrom_based=F, tsb_stat=F, seqInfo=F, cushion=100)

inputfolder = "Kucab_vcfs"
matrices <- SigProfilerMatrixGeneratorR("indels", "GRCh37", inputfolder, plot=T, exome=F,
                                        bed_file=NULL, chrom_based=F, tsb_stat=F, seqInfo=F, cushion=100)

#now merging all indel calls
library(gtools)
out1 = read.table("Filtered_INDELs/output/ID/indels.ID83.all",h=T, stringsAsFactors = F)
rownames(out1) = out1$MutationType
outK = read.table("Kucab_vcfs/output/ID/indels.ID83.all",h=T, stringsAsFactors = F)
rownames(outK) = outK$MutationType ; outK = outK[rownames(out1),mixedsort(colnames(outK[,-1]))]
outZ = read.table("Zou_vcfs/output/ID/indels.ID83.all",h=T, stringsAsFactors = F)
rownames(outZ) = outZ$MutationType ; outZ = outZ[rownames(out1),mixedsort(colnames(outZ[,-1]))]

allout = cbind(out1, outZ, outK)
write.table(allout, "Our_Zou_Kucab_indels.txt", append=F, sep="\t", row.names = FALSE, col.names = TRUE, quote=F)
