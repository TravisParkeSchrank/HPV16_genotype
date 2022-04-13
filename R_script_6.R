library(here)
library(stringr)

hla_data <- read.table("out_hla_data.csv", header = TRUE,sep = ",")
hla_data  <- subset(hla_data, hla_data$sample_TN %in% c("T"))
hla_data$hla_count <- apply(hla_data, 1, function(x){length(unique(x))}) - 2
hpv_pos <- read.csv("104_hpv16_pos_cases.csv")
hla_data$hpv_status <- hla_data$study_id %in% hpv_pos$sample
hla_data <- subset(hla_data, hla_data$hpv_status == TRUE)


all_HLAs <- names(table(c(hla_data$A1, hla_data$A2, hla_data$B1, hla_data$B2, hla_data$C1, hla_data$C2)))
all_HLAs <- paste("HLA-", all_HLAs, sep = "")

reference_hlas <- data.frame(tumor = rep("HPV16_reference", times = length(all_HLAs)), hla = all_HLAs)

reference_hlas$hla <- str_replace_all(reference_hlas$hla, "[*]", "")

individ_hla_data <- data.frame(tumor = rep(hla_data$study_id, times = 6), hla = c(hla_data$A1, hla_data$A2, hla_data$B1, hla_data$B2, hla_data$C1, hla_data$C2))
individ_hla_data$hla <- paste("HLA-", individ_hla_data$hla, sep = "")
individ_hla_data$hla <- str_replace_all(individ_hla_data$hla, "[*]", "")


fsa_files <- read.table("unc_seq_fsa_manifest.txt", header = FALSE)
tmp_sample_ID <- c()


for(i in fsa_files$V1){
  holder <- str_split(i, pattern = "_")[[1]][1]
  tmp_sample_ID <- c(tmp_sample_ID, holder)
}

for(i in 1:8){
  tmp_sample_ID[i] <- "HPV16_reference" 
}

fsa_files$sample_ID <- tmp_sample_ID
names(fsa_files) <- c("file", "sample_ID")


my_lines <- c("#!/bin/bash")
reference_fsa <- subset(fsa_files, fsa_files$sample_ID %in% c("HPV16_reference"))
for(i in 1:length(reference_fsa$file)){
  for(ii in 1:length(reference_hlas$hla)){
    holder <- paste("/home/parke/netMHCpan-4.1/netMHCpan -BA -a ", as.character(reference_hlas$hla[ii]), " /home/parke/Desktop/HPV16_genotype_paper_revision/30_netMHCpan_analysis/Gen_Protein_Sequences/individ_fastas/", as.character(reference_fsa$file[i]), " > /home/parke/Desktop/HPV16_genotype_paper_revision/30_netMHCpan_analysis/outfiles/", as.character(reference_fsa$file[i]), ".",as.character(reference_hlas$hla[ii]), ".netMHCpan.txt" ,  sep = "" )
    my_lines <- c(my_lines, holder)
  }
}
writeLines(my_lines, "reference_hpv_script.sh")


my_lines <- c()
sample_fsa <- subset(fsa_files, !(fsa_files$sample_ID %in% c("HPV16_reference")))
for(i in 1:length(sample_fsa$file)){
  holder_hla <- subset(individ_hla_data, individ_hla_data$tumor == sample_fsa$sample_ID[i])
  for(ii in 1:6){
    holder <- paste("/home/parke/netMHCpan-4.1/netMHCpan -BA -a ", as.character(holder_hla$hla[ii]), " /home/parke/Desktop/HPV16_genotype_paper_revision/30_netMHCpan_analysis/Gen_Protein_Sequences/individ_fastas/", as.character(sample_fsa$file[i]), " > /home/parke/Desktop/HPV16_genotype_paper_revision/30_netMHCpan_analysis/outfiles/", as.character(sample_fsa$file[i]), ".",as.character(holder_hla$hla[ii]), ".netMHCpan.txt" ,  sep = "" )
    my_lines <- c(my_lines, holder)
  }
}
n_real_lines <- length(my_lines)
my_lines <- c(my_lines, rep("", times = 5501 - n_real_lines))

for(i in 0:10){
  write_lines <- c("#!/bin/bash", my_lines[500*i+1:500])
  writeLines(write_lines, paste("sample_script_", as.character(i), ".sh", sep = ""))
}


