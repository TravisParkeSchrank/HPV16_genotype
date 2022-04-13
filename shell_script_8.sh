#!/bin/bash
# run with sudo


fileItemString=$(cat /home/parke/ssdb2/UNCseq_oropharynx_fastqs/oropharynx_fastq_links/oropharynx_tumor_manifest.txt |tr "\n" " ")
arr=($fileItemString)


for fn in "${arr[@]}";
do

#samtools mpileup -f /home/parke/ViFi/data/hpv/hg19_hpv.fas /home/parke/ssdb2/UNCseq_oropharynx_fastqs/ViFi_Run/${fn}/output.viral.cs.bam -o /home/parke/ssdb2/UNCseq_oropharynx_fastqs/ViFi_Run/${fn}/output.viral.cs.bam.mplieup

java -jar ~/VarScan.v2.3.9.jar mpileup2snp /home/parke/ssdb2/UNCseq_oropharynx_fastqs/ViFi_Run/${fn}/output.viral.cs.bam.mplieup --min-coverage 10 --min-reads2 5 --min-var-freq 0.001 --p-value 0.005 --output-vcf 1 > /home/parke/ssdb2/UNCseq_oropharynx_fastqs/ViFi_Run/${fn}/output.viral.cs.bam.varscan.vcf

sed 's/hpv16ref_1/NC_001526.4/g' /home/parke/ssdb2/UNCseq_oropharynx_fastqs/ViFi_Run/${fn}/output.viral.cs.bam.varscan.vcf > /home/parke/ssdb2/UNCseq_oropharynx_fastqs/ViFi_Run/${fn}/output.viral.cs.bam.varscan.rename.vcf


#must activate conda within shell script
source /home/parke/miniconda3/etc/profile.d/conda.sh
conda activate base
# this runds the conda version of snpEff
# the libraries are in /home/parke/miniconda3/share/snpeff???/data/
snpEff -v NC_001526.4 /home/parke/ssdb2/UNCseq_oropharynx_fastqs/ViFi_Run/${fn}/output.viral.cs.bam.varscan.rename.vcf > /home/parke/ssdb2/UNCseq_oropharynx_fastqs/ViFi_Run/${fn}/hpv16_snp.snpeff.vcf
 
#grep -i "NC_001526\t" /media/parke/ssdd/TCGA_HPVpos_ViFi/${fn}/hpv16_snp.snpeff.vcf >>  ./all_tumor_hpv_snps_7_23_21.vcf
 
 
done
