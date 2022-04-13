library(Gviz)
library(GenomicRanges)
library(Biostrings)
library(stringr)
library(GenomicRanges)
library(janitor)
library(rtracklayer)
library(trackViewer)
library(stringr)
library(ggsignif)
library(ggplot2)
library(wesanderson)

 
fasta_data <- readDNAStringSet("hg19_hpv.fa")
data <- read.table("hpv_annotation_data.csv", header = TRUE, sep = ",")
#data <- data[!(data$symbol %in% c("E1E4")), ]
data$chromosome <- rep("hpv16ref_1", times = length(data$chromosome))
grtrack <- GeneRegionTrack(data, genome = "hpv16", chromosome = "hpv16ref_1", name = "HPV16", options(ucscChromosomeNames=FALSE), geneSymbol= TRUE, showId=TRUE, symbol = data$symbol, stacking="squish", fill = "grey")
gtrack <- GenomeAxisTrack()
#quartz(10, 1.5, title = "test")
#plotTracks(c(gtrack, grtrack), from = 1, to = 7906)
#alTrack <- AlignmentsTrack("output.viral.cs.bam", isPaired = TRUE)
sTrack <- SequenceTrack(fasta_data, chromosome = "hpv16ref_1")
#plotTracks(c(alTrack, sTrack), chromosome = "hpv16ref_1", from = 1, to = 7906)f
####
#### Read in exaustive var possibilities
####
in_vars_all_poss <- read.table("all_pos_vars.snpeff.vcf", header = FALSE)
names(in_vars_all_poss) <- c("chromosome", "position", "blank_A", "REF", "ALT", "blank_B", "filter", "effect_raw", "scheme", "stats")
#in_vars_all_poss$study_id <- str_sub(in_vars_all_poss$study_id, start = 1L, end = 10L)
gene_altered <- c()
var_effect_class <- c()
aa_change <- c()
for(i in as.character(in_vars_all_poss$effect_raw)){
var_effect_class <- c(var_effect_class, str_split(i, "[|]")[[1]][2])
gene_altered <- c(gene_altered, str_split(i, "[|]")[[1]][4])
aa_change <- c(aa_change, str_split(i, "[|]")[[1]][11])
}
depth <- c()
vaf <- c()
for(ii in as.character(in_vars_all_poss$stats)){
depth <- c(depth, str_split(ii, "[:]")[[1]][4])
holder <- str_split(ii, "[:]")[[1]][7]
vaf<- c(vaf, as.numeric(str_split(holder, "[%]")[[1]][1]))
}
new_vars_all_poss <- data.frame(chromosome = in_vars_all_poss$chromosome, position = in_vars_all_poss$position, REF = in_vars_all_poss$REF, ALT = in_vars_all_poss$ALT, gene = gene_altered, effect_class = var_effect_class, aa_change = aa_change, depth = depth, vaf = vaf)
all_poss_vars <- new_vars_all_poss
####
#### Read in Reference Vars
####
in_vars <- read.table("all_tumors_hpv_snps_all_genes.vcf", header = FALSE)
names(in_vars) <- c("chromosome", "position", "blank_A", "REF", "ALT", "blank_B", "filter", "effect_raw", "scheme", "stats", "study_id")
in_vars$study_id <- str_sub(in_vars$study_id, start = 1L, end = 10L)
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
subclonal_vars <- subset(new_vars, new_vars$vaf <= 90)
subclonal_vars <- subset(subclonal_vars, subclonal_vars$vaf >= 5)
subclonal_vars <- subset(subclonal_vars, subclonal_vars$depth >= 100)
subclonal_vars<- subset(subclonal_vars, subclonal_vars$effect_class %in% c("missense_variant", "initiator_codon_variant", "stop_gained"))
#new_vars <- subset(new_vars, new_vars$vaf > 90)
new_vars$subclonal <- new_vars$vaf <= 90
new_vars <- subset(new_vars, new_vars$vaf >= 5)
new_vars <- subset(new_vars, new_vars$depth >= 100)
new_vars$var_id <- paste(new_vars$REF, new_vars$ALT, new_vars$position, sep = "")

var_id_df<- data.frame(table(new_vars$var_id))
var_id_freq_lib <- var_id_df$Freq
names(var_id_freq_lib) <- var_id_df$Var1

new_vars$freq <- var_id_freq_lib[new_vars$var_id]

#### Subset down to common variants in more than 25% of tumors 
#new_vars <- subset(new_vars, new_vars$freq < 20)

missense_vars <- subset(new_vars, new_vars$effect_class %in% c("missense_variant"))







new_vars_effects <- subset(new_vars, new_vars$effect_class %in% c("missense_variant", "stop_gained", "synonymous_variant", "upstream_gene_variant", "initiator_codon_variant"))
#new_vars_effects <- subset(new_vars, new_vars$effect_class %in% c("missense_variant", "synonymous_variant"))
test_effect <-ggplot(new_vars_effects, aes(x=subclonal, fill=effect_class))
mybar_effect <- test_effect+ geom_bar(position = "fill") + theme_classic() + scale_fill_manual(values = c(wes_palette("Moonrise2"),wes_palette("Moonrise3")[1:3], wes_palette("Moonrise3")[5]))
pdf(file = "subclonal_vs_effect.pdf", width = 4, height = 4)
mybar_effect
dev.off()

test <-ggplot(new_vars_effects, aes(x=subclonal, fill=gene))
mybar <- test+ geom_bar(position = "fill") + theme_classic() + scale_fill_manual(values = c(wes_palette("Moonrise2"),wes_palette("Moonrise3")[1:3], wes_palette("Moonrise3")[5]))
pdf(file = "subclonal_vs_genes.pdf", width = 3, height = 4)
mybar
dev.off()
new_vars_clonal <- subset(new_vars_effects, new_vars$subclonal == FALSE)
new_vars_subclonal <- subset(new_vars_effects, new_vars$subclonal == TRUE)
gene_table <- t(table(new_vars_effects$subclonal, new_vars_effects$gene))
print(chisq.test(gene_table))





new_vars_effects_clonal <- subset(new_vars_effects, new_vars_effects$subclonal == FALSE)
new_vars_effects_subclonal <- subset(new_vars_effects, new_vars_effects$subclonal == TRUE)
var_class_table <- t(table(new_vars_effects$subclonal, new_vars_effects$effect_class))
print(chisq.test(var_class_table))
chisq_effect <- chisq.test(var_class_table)




     