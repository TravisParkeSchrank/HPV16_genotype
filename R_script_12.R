
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
rna_test <- str_sub(names, start = -3L, end = -1L) == "RNA"
names_data <- data.frame(names = names, rnaseq = rna_test)
rna_names_data <- subset(names_data, names_data$rnaseq == TRUE)
ref_names_data <- subset(names_data, names_data$rnaseq == FALSE)
ref_names <- ref_names_data$names
test_names <- rna_names_data$names


### Get distance data construct
distance_matrix <- as.matrix(dm)[, test_names]
seq_dist <- data.frame(distance_matrix)
load("dna_pars_phylogeny_data.RData")

my_df_list <- list()

for(i in test_names){
  my_df <- data.frame(distance = seq_dist[[i]], study_id = row.names(seq_dist))
  my_df <- merge(my_df, dna_pars_phylogeny_data, by = "study_id")
  my_df <- my_df[order(my_df$distance), ]
  my_df_list[[i]] <- my_df
  
}

test_names_short <- str_sub(test_names, start = 1L, end = 10L)
test_names_short
test_phylogeny_data <- subset(dna_pars_phylogeny_data, dna_pars_phylogeny_data$study_id %in% test_names_short)

predictions <- c()
for(i in 1:length(test_names)){
  holder <- my_df_list[[test_names[i]]]
  print(test_names[i])
  holder <- holder[!(holder$study_id %in% c(test_names_short[i])), ]
  #print(majorityVote(holder$far_clade[1:3])$majority)
  print(head(holder, n = 3))
  predictions <- c(predictions, majorityVote(holder$far_clade[1:3])$majority)
}


prediction_data <- data.frame(study_id = test_names_short, prediction = predictions)
load("RNA_UNCseq_Sample_Data_hpv16_POS.RData")
prediction_data <- merge(prediction_data, dna_pars_phylogeny_data, by = "study_id")
prediction_data <- merge(prediction_data, UNCseq_Sample_Data_hpv16_POS, by= "study_id")
prediction_data$log_E67_E25 <- log(prediction_data$E6E7toE2E5, base = 2)


far_prediction_data <- subset(prediction_data, prediction_data$far_clade == TRUE)
near_prediction_data <- subset(prediction_data, prediction_data$far_clade == FALSE)
wilcox.test(far_prediction_data$E6E7toE2E5, near_prediction_data$E6E7toE2E5)


