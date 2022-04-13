#!/bin/bash

fileItemString=$(cat /home/parke/ssdb2/UNCseq_oropharynx_fastqs/oropharynx_fastq_links/oropharynx_tumor_manifest.txt |tr "\n" " ")
arr=($fileItemString)


#save the last run jut incase
#cp all_tumor_hpv_snps_no_E4.vcf ./all_tumor_hpv_snps_no_E4.vcf.archive
#make new empty .vcf to pipe output to
#rm all_Tumor_hpv_snps_no_E4.vcf
grep -i  '#' header_snpeff.vcf > all_tumors_hpv_snps_all_genes.vcf
mkdir /home/parke/Desktop/UNCseq_cervix_op/oropahrynx_cohort/6_grouped_snpeff/all_genes

for fn in "${arr[@]}";
do

mkdir /home/parke/Desktop/UNCseq_cervix_op/oropahrynx_cohort/6_grouped_snpeff/all_genes/${fn}

snpEff -v NC_001526.4 /home/parke/ssdb2/UNCseq_oropharynx_fastqs/ViFi_Run/${fn}/output.viral.cs.bam.varscan.rename.vcf > /home/parke/Desktop/UNCseq_cervix_op/oropahrynx_cohort/6_grouped_snpeff/all_genes/${fn}/hpv16_snp.snpeff_all_genes.vcf


 
grep -i "NC_001526" /home/parke/Desktop/UNCseq_cervix_op/oropahrynx_cohort/6_grouped_snpeff/all_genes/${fn}/hpv16_snp.snpeff_all_genes.vcf | sed "s/$/ ${fn}/" >>  all_tumors_hpv_snps_all_genes.vcf
 
 
done
