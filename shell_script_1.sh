#!/bin/bash

fileItemString=$(cat /home/parke/ssdb2/UNCseq_oropharynx_fastqs/oropharynx_fastq_links/oropharynx_tumor_manifest.txt |tr "\n" " ")
arr=($fileItemString)

for fn in "${arr[@]}";
do


samtools idxstats /home/parke/ssdb2/UNCseq_oropharynx_fastqs/ViFi_Run/${fn}/output.sorted.bam > /home/parke/ssdb2/UNCseq_oropharynx_fastqs/ViFi_Run/${fn}/output.sorted.bam.idxstat.txt



done
