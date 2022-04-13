#!/bin/bash

#cd
#cd ViFi
#docker pull namphuon/vifi
#source ~/.bashrc
#conda deactivate

fileItemString=$(cat /home/parke/ssdb2/UNCseq_oropharynx_fastqs/oropharynx_fastq_links/oropharynx_tumor_manifest.txt |tr "\n" " ")
arr=($fileItemString)
/
mkdir /home/parke/ssdb2/UNCseq_oropharynx_fastqs/ViFi_Run

for fn in "${arr[@]}";
do 
conda deactivate
docker pull namphuon/vifi
conda init bash
source ~/.bashrc
source ~/.bashrc
mkdir /home/parke/ssdb2/UNCseq_oropharynx_fastqs/ViFi_Run/${fn}

python $VIFI_DIR/scripts/run_vifi.py --cpus 42 -f /home/parke/ssdb2/UNCseq_oropharynx_fastqs/oropharynx_fastq_links/${fn}_T_R1.fastq.gz -r /home/parke/ssdb2/UNCseq_oropharynx_fastqs/oropharynx_fastq_links/${fn}_T_R2.fastq.gz -o /home/parke/ssdb2/UNCseq_oropharynx_fastqs/ViFi_Run/${fn}/ --docker

done
