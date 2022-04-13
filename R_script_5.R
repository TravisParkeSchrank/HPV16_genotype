
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

library(muscle)
library(Biostrings)
library(msa)
library(bios2mds)



aligned_seqs <- read.dna("aligned_reference_and_HPV16_sequences.fasta", format = "fasta")

names <- row.names(aligned_seqs)
unc_test <- str_sub(names, start = 1L, end = 2L) == "UN"
names_data <- data.frame(names = names, uncseq = unc_test)
unc_names_data <- subset(names_data, names_data$uncseq == TRUE)
ref_names_data <- subset(names_data, names_data$uncseq == FALSE)



ref_names <- ref_names_data$names
test_names <- unc_names_data$names

sublineage_assigns <- c()
sublineage_assigns_2 <- c()
for(i in test_names){
  test_align <- aligned_seqs[c(ref_names, i), ]
  hpv_phyDat <- phyDat(test_align, type = "DNA", levels = NULL)
  mt <- modelTest(hpv_phyDat)
  dna_dist <- dist.ml(hpv_phyDat, model="JC69")
  dna_dist2 <- dist.ml(hpv_phyDat, model="F81")
  #test_UPGMA <- upgma(dna_dist)
  #test_UPGMA2 <- upgma(dna_dist2)
  
  sublineage_assigns <- c(sublineage_assigns, which.min(as.matrix(dna_dist)[,17][1:16]))
  sublineage_assigns_2 <- c(sublineage_assigns_2, which.min(as.matrix(dna_dist2)[,17][1:16]))
}

out_sublineage_data <- data.frame(study_id = str_sub(test_names, start = 1L, end = 10L), sublineage_JC69 = str_sub(names(sublineage_assigns), start = 1L, end = 2L), sublineage_F81 = str_sub(names(sublineage_assigns_2), start = 1L, end = 2L))

save(out_sublineage_data, file = "out_sublineage_data.RData")

