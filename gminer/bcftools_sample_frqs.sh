#!/bin/bash -l 
#SBATCH -J clin_annot
#SBATCH -e errors_%j.txt 
#SBATCH -c 2
#SBATCH --mem=8096 

#bcftools concat -f file_names.txt | bcftools view -S adme_samples.list > clin_annot.vcf.gz  

#Recalculate AC, AN, AF

for i in G1 G2 G3 G4 G5; do 
    bcftools +fill-tags clin_annot_${i}.vcf.gz -Oz | bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%AN\t%AC\t%MAF\n' | bgzip > clin_annot_${i}maf.tab.gz
done

paste gene_ids_unsorted.txt clin_annot_${i}maf.tab > clin_annot_${i}_maf.tab

#rm clin_annot_G3maf.tab

for i in G1 G2 G3 G4 G5; do 
    awk '{print $1,$2,$9}' clin_annot_${i}_maf.tab > clin_annot_${i}maf.tab
    sort clin_annot_${i}maf.tab > clin_annot_${i}_maf.tab
    rm clin_annot_${i}maf.tab
done
