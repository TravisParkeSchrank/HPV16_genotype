
### see https://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html#4_overview_of_the_package

#require(TCGAbiolinksGUI)
require(DT)
require(dplyr)
require(SummarizedExperiment)
library(here)
library(stringr)
library(seqinr)
library(genbankr)

hpv16_genome <- readGenBank("sequence.gb")
protein_genes_in <- data.frame(cds(hpv16_genome))

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

clonal_missense_vars <- subset(clonal_vars, clonal_vars$effect_class %in% c("missense_variant"))
clonal_missense_vars$AA_POS <- as.numeric(str_sub(clonal_missense_vars$aa_change, start = 6L, end = -4L))
clonal_missense_vars$AA_ALT <- a(str_sub(clonal_missense_vars$aa_change, start = -3L, end = -1L))
clonal_missense_vars$AA_REF <- a(str_sub(clonal_missense_vars$aa_change, start = 3L, end = 5L))

vcf_variants <- clonal_missense_vars[, names(clonal_missense_vars) %in% c("AA_ALT", "AA_POS", "AA_REF", "tumor_id", "gene")]

vcf_data <- vcf_variants


############
all_altered_tumors <- names(table(vcf_data$tumor_id))
### BLANK DF
blank_list_initial <- list("E1^E4" = str_split(protein_genes_in$translation[1], pattern = ""))
blank_list_initial[["E1"]] <- str_split(protein_genes_in$translation[3], pattern = "")
blank_list_initial[["E2"]] <- str_split(protein_genes_in$translation[4], pattern = "")
blank_list_initial[["E5"]] <- str_split(protein_genes_in$translation[5], pattern = "")
blank_list_initial[["L2"]] <- str_split(protein_genes_in$translation[6], pattern = "")
blank_list_initial[["L1"]] <- str_split(protein_genes_in$translation[7], pattern = "")
blank_list_initial[["E6"]] <- str_split(protein_genes_in$translation[8], pattern = "")
blank_list_initial[["E7"]] <- str_split(protein_genes_in$translation[11], pattern = "")

blank_nested_list <- list("HPV16_A1" = blank_list_initial)
for(i in all_altered_tumors){
  blank_nested_list[[i]] <- blank_list_initial
}

## mutate proteotypes
for(i in 1:length(vcf_data$tumor_id)){
  blank_nested_list[[vcf_data$tumor_id[i]]][vcf_data$gene[i]][[1]][[1]][vcf_data$AA_POS[i]] <- vcf_data$AA_ALT[i]
}

proteotype_data <- blank_nested_list
proteotype_names <- names(proteotype_data)
gene_names <- c("E1^E4", "E1", "E2", "E5",  "L2", "L1", "E6", "E7")

reference_type <- ""
for(i in gene_names){
  reference_type <- paste(reference_type, paste(proteotype_data$HPV16_A1[i][[1]][[1]], collapse = ""), collapse = "")
  reference_type <- gsub(" ", "", reference_type, fixed = TRUE)
}

proteotypes <- list("HPV16_A1" = reference_type)

for(ii in proteotype_names){
  type <- ""
  for(i in gene_names){
    type <- paste(type, paste(proteotype_data[[ii]][i][[1]][[1]], collapse = ""), collapse = "")
    type <- gsub(" ", "", type, fixed = TRUE)
  }
  proteotypes[[ii]] <- type
}


write.fasta(proteotypes, names(proteotypes), file.out = "reference_and_sample_HPV16_proteotypes.fasta", open = "w", nbchar = 60, as.string = FALSE)




