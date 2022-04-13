#!/bin/bash
# run with sudo


fileItemString=$(cat /home/parke/ssdb2/UNCseq_oropharynx_fastqs/oropharynx_fastq_links/oropharynx_tumor_manifest.txt |tr "\n" " ")
arr=($fileItemString)


for fn in "${arr[@]}";
do


rm /home/parke/ssdb2/UNCseq_oropharynx_fastqs/ViFi_Run/${fn}/output.viral.cs.bam.varscan.rename.vcf

sed 's/hpv16ref_1/NC_001526/g' /home/parke/ssdb2/UNCseq_oropharynx_fastqs/ViFi_Run/${fn}/output.viral.cs.bam.varscan.vcf > /home/parke/ssdb2/UNCseq_oropharynx_fastqs/ViFi_Run/${fn}/output.viral.cs.bam.varscan.rename.vcf


 
done
