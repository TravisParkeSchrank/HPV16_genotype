
require(DT)
require(dplyr)
require(SummarizedExperiment)
library(here)
library(stringr)
library(seqinr)
library(ape)
library(phangorn)
library(treeman)
library(gridExtra)
library(mclust)
##IMPORT DATA FOR DIST MATRIX
my_data <- read.phyDat(here("DNA_RNA_sample_genotypes.fasta"), format = "fasta")
dm  <- dist.ml(my_data)

## GET NAME INGO
aligned_seqs <- read.dna("DNA_RNA_sample_genotypes.fasta", format = "fasta")
names <- row.names(aligned_seqs)
rna_test <- str_sub(names, start = 1L, end = 3L) != "UNC"
names_data <- data.frame(names = names, rnaseq = rna_test)
rna_names_data <- subset(names_data, names_data$rnaseq == TRUE)
ref_names_data <- subset(names_data, names_data$rnaseq == FALSE)
ref_names <- ref_names_data$names
test_names <- rna_names_data$names


### Get distance data construct
distance_matrix <- as.matrix(dm)[, test_names]
seq_dist <- data.frame(distance_matrix)
load("dna_pars_phylogeny_data.RData")

test_names_dot <- names(seq_dist)


my_df_list <- list()

for(i in test_names_dot){
  my_df <- data.frame(distance = seq_dist[[i]], study_id = row.names(seq_dist))
  my_df <- merge(my_df, dna_pars_phylogeny_data, by = "study_id")
  my_df <- my_df[order(my_df$distance), ]
  my_df_list[[i]] <- my_df
  
}



predictions <- c()
nns <- c()
for(i in 1:length(test_names_dot)){
  holder <- my_df_list[[test_names_dot[i]]]
  print(test_names_dot[i])
  holder <- holder[!(holder$study_id %in% c(test_names_dot[i])), ]
  #print(majorityVote(holder$far_clade[1:3])$majority)
  print(head(holder, n = 3))
  predictions <- c(predictions, majorityVote(holder$far_clade[1:3])$majority)
  nns <- c(nns, holder$study_id[1])
}


prediction_data <- data.frame(study_id = test_names, study_id_dot = test_names_dot, prediction = predictions, nn = nns)
load("RNA_NOT_UNCseq_Sample_Data_hpv16_POS.RData")
prediction_data <- merge(prediction_data, NOT_UNCseq_Sample_Data_hpv16_POS, by= "study_id")

far_prediction_data <- subset(prediction_data, prediction_data$prediction ==  TRUE)
near_prediction_data <- subset(prediction_data, prediction_data$prediction == FALSE)
wilcox.test(far_prediction_data$E6E7toE2E5, near_prediction_data$E6E7toE2E5)

load("hpv_pos_ssGSEA_data.RData")
merged_genomics_genotype_prediction <- merge(prediction_data, hpv_pos_ssGSEA_data, by = "study_id")

save(merged_genomics_genotype_prediction, file = "merged_genomics_genotpye_prediction.RData")



