#For SNVs
mkdir txts
for sample in ./vcfs/*.vcf
   do
          name=${sample%.vcf}
          echo $name
          bcftools query -f "%CHROM\t%POS\t%REF\t%ALT[\t%GT][\t%GQ][\t%AD][\t%DP]\n" $sample > "$name.txt"
          mv "$name.txt" ./txts
   done

#For INDELs
mkdir txts
for sample in ./vcfs/*.vcf
   do
          name=${sample%.vcf}
          echo $name
          bcftools query -f "%CHROM\t%POS\t%REF\t%ALT[\t%GT][\t%GQ][\t%AD][\t%DPI]\n" $sample > "$name.txt"
          mv "$name.txt" ./txts
   done
