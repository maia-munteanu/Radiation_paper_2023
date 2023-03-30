####################### SNVs ####################### 


#first, we need to filter the vcfs so that we only keep variants with one non-ref genotype
bcftools view -i 'count(GT!="ref")==1' ../MCF7_PASS.hg38_multianno_k50.vcf.gz | bcftools sort -Oz > MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz
bcftools view -i 'count(GT!="ref")==1' ../HAP1_PASS.hg38_multianno_k50.vcf.gz | bcftools sort -Oz > HAP1_PASS.hg38_multianno_k50_discrete_variants.vcf.gz
bcftools view -i 'count(GT!="ref")==1' ../A549_PASS.hg38_multianno_k50.vcf.gz | bcftools sort -Oz > A549_PASS.hg38_multianno_k50_discrete_variants.vcf.gz

#now getting indels and snvs for each cell line
bcftools query -l A549_PASS.hg38_multianno_k50_discrete_variants.vcf.gz
bcftools view -s A549_CA_REP1 A549_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > A549_CA_snvs.vcf.gz
bcftools view -s A549_CC_REP1 A549_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > A549_CC_snvs.vcf.gz
bcftools view -s A549_HE1A_REP1 A549_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > A549_HE1A_snvs.vcf.gz
bcftools view -s A549_HE2B_REP1 A549_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > A549_HE2B_snvs.vcf.gz
bcftools view -s A549_PR1B_REP1 A549_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > A549_PR1B_snvs.vcf.gz
bcftools view -s A549_PR2B_REP1 A549_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > A549_PR2B_snvs.vcf.gz

bcftools query -l HAP1_PASS.hg38_multianno_k50_discrete_variants.vcf.gz
bcftools view -s HAP1_CA_REP1 HAP1_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > HAP1_CA_snvs.vcf.gz
bcftools view -s HAP1_CB_REP1 HAP1_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > HAP1_CB_snvs.vcf.gz
bcftools view -s HAP1_HE1A_REP1 HAP1_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > HAP1_HE1A_snvs.vcf.gz
bcftools view -s HAP1_HE2A_REP1 HAP1_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > HAP1_HE2A_snvs.vcf.gz
bcftools view -s HAP1_PR1A_REP1 HAP1_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > HAP1_PR1A_snvs.vcf.gz
bcftools view -s HAP1_PR2A_REP1 HAP1_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > HAP1_PR2A_snvs.vcf.gz

bcftools query -l MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz
bcftools view -s MCF7_untreated_F8_REP2 MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > MCF7_C_snvs.vcf.gz
bcftools view -s MCF7_HE_1B_REP1 MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > MCF7_HE1B_snvs.vcf.gz
bcftools view -s MCF7_HE_1C_REP1 MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > MCF7_HE1C_snvs.vcf.gz
bcftools view -s MCF7_HE_2A_REP1 MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > MCF7_HE2A_snvs.vcf.gz
bcftools view -s MCF7_HE_2B_REP1 MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > MCF7_HE2B_snvs.vcf.gz
bcftools view -s MCF7_PR_1A_REP1 MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > MCF7_PR1A_snvs.vcf.gz
bcftools view -s MCF7_PR_1B_REP1 MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > MCF7_PR1B_snvs.vcf.gz
bcftools view -s MCF7_PR_2A_REP1 MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > MCF7_PR2A_snvs.vcf.gz
bcftools view -s MCF7_PR_2C_REP1 MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types snps | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > MCF7_PR2C_snvs.vcf.gz

gunzip *.vcf.gz
mkdir vcfs
mv *vcf ./vcfs




####################### INDELS ####################### 

#first, we need to filter the vcfs so that we only keep variants with one non-ref genotype
bcftools view -i 'count(GT!="ref")==1' ../MCF7_PASS.hg38_multianno_k50.vcf.gz | bcftools sort -Oz > MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz
bcftools view -i 'count(GT!="ref")==1' ../HAP1_PASS.hg38_multianno_k50.vcf.gz | bcftools sort -Oz > HAP1_PASS.hg38_multianno_k50_discrete_variants.vcf.gz
bcftools view -i 'count(GT!="ref")==1' ../A549_PASS.hg38_multianno_k50.vcf.gz | bcftools sort -Oz > A549_PASS.hg38_multianno_k50_discrete_variants.vcf.gz

#now getting indels and snvs for each cell line
bcftools query -l A549_PASS.hg38_multianno_k50_discrete_variants.vcf.gz
bcftools view -s A549_CA_REP1 A549_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > A5
49_CA_indels.vcf.gz
bcftools view -s A549_CC_REP1 A549_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > A5
49_CC_indels.vcf.gz
bcftools view -s A549_HE1A_REP1 A549_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > 
A549_HE1A_indels.vcf.gz
bcftools view -s A549_HE2B_REP1 A549_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > 
A549_HE2B_indels.vcf.gz
bcftools view -s A549_PR1B_REP1 A549_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > 
A549_PR1B_indels.vcf.gz
bcftools view -s A549_PR2B_REP1 A549_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > A549_PR2B_indels.vcf.gz

bcftools query -l HAP1_PASS.hg38_multianno_k50_discrete_variants.vcf.gz
bcftools view -s HAP1_CA_REP1 HAP1_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > HAP1_CA_indels.vcf.gz
bcftools view -s HAP1_CB_REP1 HAP1_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > HAP1_CB_indels.vcf.gz
bcftools view -s HAP1_HE1A_REP1 HAP1_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > HAP1_HE1A_indels.vcf.gz
bcftools view -s HAP1_HE2A_REP1 HAP1_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > HAP1_HE2A_indels.vcf.gz
bcftools view -s HAP1_PR1A_REP1 HAP1_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > HAP1_PR1A_indels.vcf.gz
bcftools view -s HAP1_PR2A_REP1 HAP1_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > HAP1_PR2A_indels.vcf.gz

bcftools query -l MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz
bcftools view -s MCF7_untreated_F8_REP2 MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > MCF7_C_indels.vcf.gz
bcftools view -s MCF7_HE_1B_REP1 MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > MCF7_HE1B_indels.vcf.gz
bcftools view -s MCF7_HE_1C_REP1 MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > MCF7_HE1C_indels.vcf.gz
bcftools view -s MCF7_HE_2A_REP1 MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > MCF7_HE2A_indels.vcf.gz
bcftools view -s MCF7_HE_2B_REP1 MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > MCF7_HE2B_indels.vcf.gz
bcftools view -s MCF7_PR_1A_REP1 MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > MCF7_PR1A_indels.vcf.gz
bcftools view -s MCF7_PR_1B_REP1 MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > MCF7_PR1B_indels.vcf.gz
bcftools view -s MCF7_PR_2A_REP1 MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > MCF7_PR2A_indels.vcf.gz
bcftools view -s MCF7_PR_2C_REP1 MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz | bcftools view -i 'GT="alt"' | bcftools view --types indels | bcftools view -e 'gnomAD_genome_ALL>=0.001' | bcftools view -i 'FT="PASS"' | bcftools sort -Oz > MCF7_PR2C_indels.vcf.gz

rm MCF7_PASS.hg38_multianno_k50_discrete_variants.vcf.gz
rm A549_PASS.hg38_multianno_k50_discrete_variants.vcf.gz
rm HAP1_PASS.hg38_multianno_k50_discrete_variants.vcf.gz

gunzip *vcf.gz
mdir vcfs
mv *vcf ./vcfs
