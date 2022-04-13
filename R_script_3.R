#library(Gviz)
#library(GenomicRanges)
library(Biostrings)
library(stringr)
#library(GenomicRanges)
library(janitor)
#library(rtracklayer)
#library(trackViewer)

fasta_data <- readDNAStringSet("hg19_hpv.fa")
 data <- read.table("hpv_annotation_data.csv", header = TRUE, sep = ",")

in_vars <- read.table("all_tumors_hpv_snps_all_genes.vcf", header = FALSE, comment.char = "#")

names(in_vars) <- c("chromosome", "position", "blank_A", "REF", "ALT", "blank_B", "filter", "effect_raw", "scheme", "stats", "study_id")


gene_altered <- c()
var_effect_class <- c()
aa_change <- c()


for(i in as.character(in_vars$effect_raw)){
	var_effect_class <- c(var_effect_class, str_split(i, "[|]")[[1]][2])
	gene_altered <- c(gene_altered, str_split(i, "[|]")[[1]][4])
	aa_change <- c(aa_change, str_split(i, "[|]")[[1]][11])
}


depth <- c()
vaf <- c()

for(ii in as.character(in_vars$stats)){
	depth <- c(depth, str_split(ii, "[:]")[[1]][4])
	holder <- str_split(ii, "[:]")[[1]][7]
	vaf<- c(vaf, as.numeric(str_split(holder, "[%]")[[1]][1]))
}


new_vars <- data.frame(chromosome = in_vars$chromosome, position = in_vars$position, REF = in_vars$REF, ALT = in_vars$ALT, gene = gene_altered, effect_class = var_effect_class, aa_change = aa_change, depth = depth, vaf = vaf, study_id = in_vars$study_id)
new_vars$study_id <- str_sub(new_vars$study_id, start = 1L, end = 10L)

subclonal_vars <- subset(new_vars, new_vars$vaf <= 90)
subclonal_coding_vars <- subset(subclonal_vars, subclonal_vars$effect_class != "synonymous_variant")

#subclonal_vars <- subset(subclonal_vars, subclonal_vars$vaf >= 5)
save(subclonal_vars, file = "subclonal_vars.RData")
#subclonal_vars<- subset(subclonal_vars, subclonal_vars$effect_class %in% c("missense_variant", "initiator_codon_variant", "stop_gained"))

new_vars <- subset(new_vars, new_vars$vaf > 90)


missense_vars <- subset(new_vars, new_vars$effect_class %in% c("missense_variant"))

stop_gain_vars <- subset(new_vars, new_vars$effect_class %in% c("initiator_codon_variant", "stop_gained"))

silent_vars <- subset(new_vars, new_vars$effect_class %in% c("synonymous_variant"))

unique_missense_vars <- unique(missense_vars)


my_hist <- hist(missense_vars$position, n=7906, plot = FALSE)

hist_data <- data.frame(position = my_hist$breaks[1:length(my_hist$counts)], counts = my_hist$counts)


n_clonal_missense_vars <- data.frame(table(missense_vars$study_id))
names(n_clonal_missense_vars) <- c("study_id", "n_clonal_missense")
save(n_clonal_missense_vars, file = "n_clonal_missense.RData")

L2 <- subset(missense_vars, missense_vars$gene %in% c("L2"))
L2$aa_pos <- str_sub(L2$aa_change, start = 6L, end = -4L)
L2_pos_269 <- subset(L2, L2$aa_pos == "269")

save(L2_pos_269, file = "L2_pos_269.RData")


pos_data  <- read.csv("104_hpv16_pos_cases.csv")

new_vars$var_id <- paste(new_vars$gene, new_vars$aa_change, sep  = "_")
clonal_missense <- subset(new_vars, new_vars$effect_class %in% c("missense_variant"))     
clonal_missense <- subset(clonal_missense, clonal_missense$study_id %in% pos_data$sample)
common_poly_data_op <- data.frame(table(clonal_missense$var_id)[rev(order(table(clonal_missense$var_id)))])
names(common_poly_data_op) <- c("var_id", "Freq")
n_tumors <- 104

prop <- c()
low_ci <- c()
high_ci <- c()

for(i in 1:length(common_poly_data_op$var_id)){
  prop <- c(prop, common_poly_data_op$Freq[i]/n_tumors)
  low_ci <- c(low_ci, prop.test(common_poly_data_op$Freq[i],n_tumors)[[6]][1])
  high_ci <- c(high_ci, prop.test(common_poly_data_op$Freq[i],n_tumors)[[6]][2])
}

common_poly_data_op$prop <- prop
common_poly_data_op$low_ci <- low_ci
common_poly_data_op$high_ci  <- high_ci

save(common_poly_data_op, file = "common_poly_data_oropharynx.RData")



     