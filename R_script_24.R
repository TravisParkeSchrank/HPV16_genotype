
### see https://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html#4_overview_of_the_package

#require(TCGAbiolinksGUI)
require(DT)
require(dplyr)
require(SummarizedExperiment)
library(here)
library(stringr)
library(seqinr)

hpv16_ref1 <- read.fasta("reference_hpv16_snpeff.fasta")
hpv16_ref<- toupper(hpv16_ref1$NC_001526.4)
hpv16_ref <- as.character(hpv16_ref)






##### truncate a sample name to make it in normal TCGA format
truncate <- function(x){
  return(substr(x, 1, 12))
}
##### truncate a sample name to make it in normal TCGA format
truncate_new <- function(x){
  return(str_sub(x, start = -12L, end = -3L))
}
##### import vars
in_vars <- read.table("all_tumors_hpv_snps_all_genes.vcf", header = FALSE, comment.char = "#")
names(in_vars) <- c("chromosome", "position", "blank_A", "REF", "ALT", "blank_B", "filter", "effect_raw", "scheme", "stats", "study_id")
##var effect extrations
gene_altered <- c()
var_effect_class <- c()
aa_change <- c()
for(i in as.character(in_vars$effect_raw)){
  var_effect_class <- c(var_effect_class, str_split(i, "[|]")[[1]][2])
  gene_altered <- c(gene_altered, str_split(i, "[|]")[[1]][4])
  aa_change <- c(aa_change, str_split(i, "[|]")[[1]][11])
}
#### depth and VAF extraction
depth <- c()
vaf <- c()
for(ii in as.character(in_vars$stats)){
  depth <- c(depth, str_split(ii, "[:]")[[1]][4])
  holder <- str_split(ii, "[:]")[[1]][7]
  vaf<- c(vaf, as.numeric(str_split(holder, "[%]")[[1]][1]))
}

new_vars <- data.frame(chromosome = in_vars$chromosome, POS = in_vars$position, REF = in_vars$REF, ALT = in_vars$ALT, gene = gene_altered, effect_class = var_effect_class, aa_change = aa_change, depth = depth, vaf = vaf, study_id = in_vars$study_id)
new_vars$tumor_id <- new_vars$study_id

clonal_vars <- subset(new_vars, new_vars$vaf >= 90)

vcf_variants <- clonal_vars[, names(clonal_vars) %in% c("ALT", "POS", "tumor_id")]


#####
#vcf_data <- read.table("all_tumors_hpv_snps_all_genes.vcf", header = FALSE, comment.char = "#")
#names(vcf_data) <- c("chromosome", "position", "blank_A", "REF", "ALT", "blank_B", "filter", "effect_raw", "scheme", "stats", "study_id")



##### VCF normal #######
vcf_normal <- data.frame("POS" = 1:7906, "REF" = hpv16_ref )
vcf_normal <- vcf_normal[order(vcf_normal$POS), ]

##### VCF Variants #######


vcf_data <- vcf_variants


############
positions_altered <- names(table(vcf_data$POS))
all_altered_tumors <- names(table(vcf_data$tumor_id))
### BLANK DF
blank_df_initial <- data.frame(vcf_normal)
row.names(blank_df_initial) <- blank_df_initial$POS
for(i in all_altered_tumors){
  blank_df_initial[[i]] <- vcf_normal$REF
}


#blank_df_initial <- blank_df_initial[, !(names(blank_df_initial) %in% c("REF"))]
##########
##########
genotype_df <- blank_df_initial

for(i in 1:length(vcf_variants$tumor_id)){
  genotype_df[as.character(vcf_variants$POS[i]), vcf_variants$tumor_id[i]] <- vcf_variants$ALT[i]
}


test <- genotype_df[, !(names(genotype_df) %in% c("POS", "REF"))]

genotopye_rows <- data.frame(t(test))
unique_genotypes<- unique(genotopye_rows)



sequences <- c()
for(i in names(test)){
  sequences <- c(sequences, str_c(test[[i]], collapse = ""))
}

#sequences <- as.list(sequences)

library(seqinr)


A1 <- read.fasta("A1.fasta")
A1 <- toupper(A1$A1)
A1 <- as.character(A1)


A2 <- read.fasta("A2.fasta")
A2 <- toupper(A2$A2)
A2 <- as.character(A2)

A3 <- read.fasta("A3.fasta")
A3 <- toupper(A3$A3)
A3 <- as.character(A3)

A4 <- read.fasta("A4.fasta")
A4 <- toupper(A4$A4)
A4 <- as.character(A4)

B1 <- read.fasta("B1.fasta")
B1 <- toupper(B1$B1)
B1 <- as.character(B1)

B2 <- read.fasta("B2.fasta")
B2 <- toupper(B2$B2)
B2 <- as.character(B2)

B3 <- read.fasta("B3.fasta")
B3 <- toupper(B3$B3)
B3 <- as.character(B3)

B4 <- read.fasta("B4.fasta")
B4 <- toupper(B4$B4)
B4 <- as.character(B4)

C1 <- read.fasta("C1.fasta")
C1 <- toupper(C1$C1)
C1 <- as.character(C1)

C2 <- read.fasta("C2.fasta")
C2 <- toupper(C2$C2)
C2 <- as.character(C2)

C3 <- read.fasta("C3.fasta")
C3 <- toupper(C3$C3)
C3 <- as.character(C3)

C4 <- read.fasta("C4.fasta")
C4 <- toupper(C4$C4)
C4 <- as.character(C4)

D1 <- read.fasta("D1.fasta")
D1 <- toupper(D1$D1)
D1 <- as.character(D1)

D2 <- read.fasta("D2.fasta")
D2 <- toupper(D2$D2)
D2 <- as.character(D2)

D3 <- read.fasta("D3.fasta")
D3 <- toupper(D3$D3)
D3 <- as.character(D3)

D4 <- read.fasta("D4.fasta")
D4 <- toupper(D4$D4)
D4 <- as.character(D4)

#str_c(test[[i]], collapse = "")
sequences <- c(sequences, str_c(A1, collapse = ""), str_c(A2, collapse = ""), str_c(A3, collapse = ""), str_c(A4, collapse = ""), str_c(B1, collapse = ""),str_c(B2, collapse = ""), str_c(B3, collapse = ""), str_c(B4, collapse = ""), str_c(C1, collapse = ""), str_c(C2, collapse = ""), str_c(C3, collapse = ""), str_c(C4, collapse = ""), str_c(D1, collapse = ""), str_c(D2, collapse = ""), str_c(D3, collapse = ""), str_c(D4, collapse = ""))

out_seq <- as.list(sequences)
seq_names <- c(names(test), "A1", "A2", "A3", "A4", "B1","B2", "B3", "B4", "C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4")


write.fasta(out_seq, seq_names, file.out = "reference_and_sample_HPV16_genotypes.fasta", open = "w", nbchar = 60, as.string = FALSE)








#


