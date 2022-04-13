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
subclonal_vars <- subset(subclonal_vars, subclonal_vars$vaf >= 1)
subclonal_vars<- subset(subclonal_vars, subclonal_vars$effect_class %in% c("missense_variant", "initiator_codon_variant", "stop_gained"))
new_vars <- subset(new_vars, new_vars$vaf > 90)
new_vars$var_id <- paste(new_vars$REF, new_vars$ALT, new_vars$position, sep = "")

var_id_df<- data.frame(table(new_vars$var_id))
var_id_freq_lib <- var_id_df$Freq
names(var_id_freq_lib) <- var_id_df$Var1

new_vars$freq <- var_id_freq_lib[new_vars$var_id]

#### Subset down to common variants in more than 25% of tumors 
#new_vars <- subset(new_vars, new_vars$freq < 20)

missense_vars <- subset(new_vars, new_vars$effect_class %in% c("missense_variant"))
rare_miss_vars <- missense_vars
stop_gain_vars <- subset(new_vars, new_vars$effect_class %in% c("initiator_codon_variant", "stop_gained"))
silent_vars <- subset(new_vars, new_vars$effect_class %in% c("synonymous_variant"))
unique_missense_vars <- unique(missense_vars)
my_hist <- hist(missense_vars$position, n=7906, plot = FALSE)
hist_data <- data.frame(position = my_hist$breaks[1:length(my_hist$counts)], counts = my_hist$counts)
hist_data <- subset(hist_data, hist_data$counts > 0 )
my_gr_var_counts <- GRanges(
seqnames = rep("hpv16ref_1", times = length(hist_data$counts)),
ranges = IRanges(hist_data$position, end = hist_data$position),
strand = rep("+", times = length(hist_data$counts)),
var_counts = hist_data$counts
)
dTrack <- DataTrack(my_gr_var_counts, name="missense", type = "histogram")
my_hist <- hist(unique_missense_vars$position, n=7906, plot = FALSE)
hist_data <- data.frame(position = my_hist$breaks[1:length(my_hist$counts)], counts = my_hist$counts)
hist_data <- subset(hist_data, hist_data$counts > 0 )
my_gr_var_counts_unique <- GRanges(
seqnames = rep("hpv16ref_1", times = length(hist_data$counts)),
ranges = IRanges(hist_data$position, end = hist_data$position),
strand = rep("+", times = length(hist_data$counts)),
var_counts = hist_data$counts
)
dTrack_unique <- DataTrack(my_gr_var_counts_unique, name="unique_varaints", type = "histogram")
my_hist <- hist(silent_vars$position, n=7906, plot = FALSE)
hist_data <- data.frame(position = my_hist$breaks[1:length(my_hist$counts)], counts = my_hist$counts)
hist_data <- subset(hist_data, hist_data$counts > 0 )
my_gr_var_counts <- GRanges(
seqnames = rep("hpv16ref_1", times = length(hist_data$counts)),
ranges = IRanges(hist_data$position, end = hist_data$position),
strand = rep("+", times = length(hist_data$counts)),
var_counts = hist_data$counts
)
dTrack_silent <- DataTrack(my_gr_var_counts, name="silent", type = "histogram")
my_hist <- hist(stop_gain_vars$position, n=7906, plot = FALSE)
hist_data <- data.frame(position = my_hist$breaks[1:length(my_hist$counts)], counts = my_hist$counts)
hist_data <- subset(hist_data, hist_data$counts > 0 )
my_gr_var_counts <- GRanges(
seqnames = rep("hpv16ref_1", times = length(hist_data$counts)),
ranges = IRanges(hist_data$position, end = hist_data$position),
strand = rep("+", times = length(hist_data$counts)),
var_counts = hist_data$counts
)
dTrack_stop <- DataTrack(my_gr_var_counts, name="stop gained", type = "histogram")
my_hist <- hist(subclonal_vars$position, n=7906, plot = FALSE)
hist_data <- data.frame(position = my_hist$breaks[1:length(my_hist$counts)], counts = my_hist$counts)
hist_data <- subset(hist_data, hist_data$counts > 0 )
my_gr_var_counts <- GRanges(
seqnames = rep("hpv16ref_1", times = length(hist_data$counts)),
ranges = IRanges(hist_data$position, end = hist_data$position),
strand = rep("+", times = length(hist_data$counts)),
var_counts = hist_data$counts
)
dTrack_subclonal <- DataTrack(my_gr_var_counts, name="subclonal", type = "histogram")
plotTracks(c(gtrack, grtrack, dTrack,  dTrack_stop,dTrack_silent,  dTrack_subclonal), chromosome = "hpv16ref_1", from = 1, to = 7906)
old_new_vars <- new_vars
#test <- get_dupes(new_vars, position, REF, ALT)
#dupe_vars <- data.frame(test)
#rare_vars <- subset(dupe_vars, dupe_vars$dupe_count < 20)
#rare_miss_vars <- subset(rare_vars, rare_vars$effect_class %in% c("missense_variant"))
#my_new_vars <- subset(new_vars,  new_vars$position %in% rare_vars$position)
#my_new_vars_E2 <- my_new_vars
#my_new_vars <- subset(my_new_vars, !(my_new_vars$gene %in% c("E2", "E1^E4")))
my_new_vars <- new_vars
data$length <- data$end - data$start
stat_gene <- c()
miss_misssynon_ratio_test <- c()
miss_misssynon_ratio_control <- c()
miss_misssynon_ratio_test_CI_low <- c()
miss_misssynon_ratio_control_CI_low <- c()
miss_misssynon_ratio_test_CI_high <- c()
miss_misssynon_ratio_control_CI_high <- c()
pval_missynon_ratio <- c()
miss_per_length <- c()
pval_miss_per_length <- c()
#### SET UP E7
E7_vars <- subset(my_new_vars, my_new_vars$gene %in% c("E7"))
not_E7_vars <- subset(my_new_vars, !(my_new_vars$gene %in% c("E7")))
all_poss_E7_vars <- subset(all_poss_vars, all_poss_vars$gene %in% c("E7"))
#print("binomial test of conservation of E7 COMPARED TO RANDOM")
#possible_vars_n = table(all_poss_E7_vars$effect_class)[["missense_variant"]] + table(all_poss_E7_vars$effect_class)[["synonymous_variant"]]
#possible_missense_n = table(all_poss_E7_vars$effect_class)[["missense_variant"]]
#vars_n = table(E7_vars$effect_class)[["missense_variant"]] + table(E7_vars$effect_class)[["synonymous_variant"]]
#missense_n = table(E7_vars$effect_class)[["missense_variant"]]
#print(binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n ))
save_E7_vars <- E7_vars
print("binomial test of conservation of E7 COMPARED TO HPV")
print("binomial test of conservation of E7 COMPARED TO HPV")
print("binomial test of conservation of E7 COMPARED TO HPV")
print("binomial test of conservation of E7 COMPARED TO HPV")
possible_vars_n = table(not_E7_vars$effect_class)[["missense_variant"]] + table(not_E7_vars$effect_class)[["synonymous_variant"]]
possible_missense_n = table(not_E7_vars$effect_class)[["missense_variant"]]
vars_n = table(E7_vars$effect_class)[["missense_variant"]] + table(E7_vars$effect_class)[["synonymous_variant"]]
missense_n = table(E7_vars$effect_class)[["missense_variant"]]
print(binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n ))
print(prop.test(c(missense_n, possible_missense_n), c(vars_n, possible_vars_n)))
print("number of missense")
print(missense_n)
print(binom.test(x = missense_n , n = data$length[9]  , p=dim(rare_miss_vars)[1]/sum(data$length)))
stat_gene <- c(stat_gene, "E7")
miss_misssynon_ratio_test <- c(miss_misssynon_ratio_test, missense_n/vars_n)
miss_misssynon_ratio_control <- c(miss_misssynon_ratio_control, possible_missense_n/possible_vars_n)
pval_missynon_ratio <- c(pval_missynon_ratio, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )$p.value)
miss_per_length <- c(miss_per_length, missense_n/data$length[9])
pval_miss_per_length <- c(pval_miss_per_length, binom.test(x = missense_n , n = data$length[9]  , p=dim(rare_miss_vars)[1]/sum(data$length))$p.value)
miss_misssynon_ratio_test_CI_low <- c(miss_misssynon_ratio_test_CI_low, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )[[4]][1])
miss_misssynon_ratio_test_CI_high <- c(miss_misssynon_ratio_test_CI_high, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )[[4]][2])
miss_misssynon_ratio_control_CI_low <- c(miss_misssynon_ratio_control_CI_low, binom.test(x = possible_missense_n , n = possible_vars_n , p=possible_missense_n/possible_vars_n )[[4]][1])
miss_misssynon_ratio_control_CI_high <- c(miss_misssynon_ratio_control_CI_high, binom.test(x = possible_missense_n , n = possible_vars_n , p=possible_missense_n/possible_vars_n )[[4]][2])
#### SET UP E2
E2_vars <- subset(my_new_vars, my_new_vars$gene %in% c("E2"))
not_E2_vars <- subset(my_new_vars, !(my_new_vars$gene %in% c("E2")))
all_poss_E2_vars <- subset(all_poss_vars, all_poss_vars$gene %in% c("E2"))
#print("binomial test of conservation of E2 COMPARED TO RANDOM")
possible_vars_n = table(all_poss_E2_vars$effect_class)[["missense_variant"]] + table(all_poss_E2_vars$effect_class)[["synonymous_variant"]] #### grep - synonomous out of all vars because some are hidden in E1E4 as missense
possible_missense_n = table(all_poss_E2_vars$effect_class)[["missense_variant"]]
vars_n = table(E2_vars$effect_class)[["missense_variant"]] + table(E2_vars$effect_class)[["synonymous_variant"]] #### grep - synonomous out of all vars because some are hidden in E1E4 as missense
missense_n = table(E2_vars$effect_class)[["missense_variant"]]
#print(binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n ))
print("binomial test of conservation of E2 COMPARED TO HPV")
print("binomial test of conservation of E2 COMPARED TO HPV")
print("binomial test of conservation of E2 COMPARED TO HPV")
print("binomial test of conservation of E2 COMPARED TO HPV")
print("binomial test of conservation of E2 COMPARED TO HPV")
print("binomial test of conservation of E2 COMPARED TO HPV")
possible_vars_n = table(not_E2_vars$effect_class)[["missense_variant"]] + table(not_E2_vars$effect_class)[["synonymous_variant"]]
possible_missense_n = table(not_E2_vars$effect_class)[["missense_variant"]]
vars_n = table(E2_vars$effect_class)[["missense_variant"]] +  table(E2_vars$effect_class)[["synonymous_variant"]] #### grep - synonomous out of all vars because some are hidden in E1E4 as missense
missense_n = table(E2_vars$effect_class)[["missense_variant"]]
print(binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n ))
print(prop.test(c(missense_n, possible_missense_n), c(vars_n, possible_vars_n)))
print("number of missense")
print(missense_n)
print(binom.test(x = missense_n , n = data$length[4]  , p=dim(missense_vars)[1]/sum(data$length)))
save_E2_vars <- E2_vars
stat_gene <- c(stat_gene, "E2")
miss_misssynon_ratio_test <- c(miss_misssynon_ratio_test, missense_n/vars_n)
miss_misssynon_ratio_control <- c(miss_misssynon_ratio_control, possible_missense_n/possible_vars_n)
pval_missynon_ratio <- c(pval_missynon_ratio, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )$p.value)
miss_per_length <- c(miss_per_length, missense_n/data$length[4])
pval_miss_per_length <- c(pval_miss_per_length, binom.test(x = missense_n , n = data$length[4]  , p=dim(missense_vars)[1]/sum(data$length))$p.value)
miss_misssynon_ratio_test_CI_low <- c(miss_misssynon_ratio_test_CI_low, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )[[4]][1])
miss_misssynon_ratio_test_CI_high <- c(miss_misssynon_ratio_test_CI_high, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )[[4]][2])
miss_misssynon_ratio_control_CI_low <- c(miss_misssynon_ratio_control_CI_low, binom.test(x = possible_missense_n , n = possible_vars_n , p=possible_missense_n/possible_vars_n )[[4]][1])
miss_misssynon_ratio_control_CI_high <- c(miss_misssynon_ratio_control_CI_high, binom.test(x = possible_missense_n , n = possible_vars_n , p=possible_missense_n/possible_vars_n )[[4]][2])
##### E14 Section
#### SET UP E14
E1E4_vars <- subset(my_new_vars, my_new_vars$gene %in% c("E1^E4"))
not_E1E4_vars <- subset(my_new_vars, !(my_new_vars$gene %in% c("E1^E4")))
all_poss_E1E4_vars <- subset(all_poss_vars, all_poss_vars$gene %in% c("E1^E4"))
#print("binomial test of conservation of E2 COMPARED TO RANDOM")
possible_vars_n = table(all_poss_E1E4_vars$effect_class)[["missense_variant"]] + table(all_poss_E1E4_vars$effect_class)[["synonymous_variant"]] #### grep - synonomous out of all vars because some are hidden in E1E4 as missense
possible_missense_n = table(all_poss_E1E4_vars$effect_class)[["missense_variant"]]
vars_n = table(E1E4_vars$effect_class)[["missense_variant"]] + table(E1E4_vars$effect_class)[["synonymous_variant"]] #### grep - synonomous out of all vars because some are hidden in E1E4 as missense
missense_n = table(E1E4_vars$effect_class)[["missense_variant"]]
print(binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n ))
print("binomial test of conservation of E1^E4 COMPARED TO HPV")
print("binomial test of conservation of E1^E4 COMPARED TO HPV")
print("binomial test of conservation of E1^E4 COMPARED TO HPV")
print("binomial test of conservation of E1^E4 COMPARED TO HPV")
print("binomial test of conservation of E1^E4 COMPARED TO HPV")
print("binomial test of conservation of E1^E4 COMPARED TO HPV")
possible_vars_n = table(not_E1E4_vars$effect_class)[["missense_variant"]] + table(not_E1E4_vars$effect_class)[["synonymous_variant"]]
possible_missense_n = table(not_E1E4_vars$effect_class)[["missense_variant"]]
vars_n = table(E1E4_vars$effect_class)[["missense_variant"]] +  table(E1E4_vars$effect_class)[["synonymous_variant"]] #### grep - synonomous out of all vars because some are hidden in E1E4 as missense
missense_n = table(E1E4_vars$effect_class)[["missense_variant"]]
print(binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n ))
print(prop.test(c(missense_n, possible_missense_n), c(vars_n, possible_vars_n)))
print("number of missense")
print(missense_n)
print(binom.test(x = missense_n , n = data$length[4]  , p=dim(missense_vars)[1]/sum(data$length)))
save_E1E4_vars <- E1E4_vars
stat_gene <- c(stat_gene, "E1E4")
miss_misssynon_ratio_test <- c(miss_misssynon_ratio_test, missense_n/vars_n)
miss_misssynon_ratio_control <- c(miss_misssynon_ratio_control, possible_missense_n/possible_vars_n)
pval_missynon_ratio <- c(pval_missynon_ratio, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )$p.value)
miss_per_length <- c(miss_per_length, missense_n/data$length[4])
pval_miss_per_length <- c(pval_miss_per_length, binom.test(x = missense_n , n = data$length[4]  , p=dim(rare_miss_vars)[1]/sum(data$length))$p.value)
miss_misssynon_ratio_test_CI_low <- c(miss_misssynon_ratio_test_CI_low, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )[[4]][1])
miss_misssynon_ratio_test_CI_high <- c(miss_misssynon_ratio_test_CI_high, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )[[4]][2])
miss_misssynon_ratio_control_CI_low <- c(miss_misssynon_ratio_control_CI_low, binom.test(x = possible_missense_n , n = possible_vars_n , p=possible_missense_n/possible_vars_n )[[4]][1])
miss_misssynon_ratio_control_CI_high <- c(miss_misssynon_ratio_control_CI_high, binom.test(x = possible_missense_n , n = possible_vars_n , p=possible_missense_n/possible_vars_n )[[4]][2])
#### SET UP E6
E6_vars <- subset(my_new_vars, my_new_vars$gene %in% c("E6"))
not_E6_vars <- subset(my_new_vars, !(my_new_vars$gene %in% c("E6")))
all_poss_E6_vars <- subset(all_poss_vars, all_poss_vars$gene %in% c("E6"))
print("binomial test of conservation of E6 COMPARED TO RANDOM")
possible_vars_n = table(all_poss_E6_vars$effect_class)[["missense_variant"]] + table(all_poss_E6_vars$effect_class)[["synonymous_variant"]]
possible_missense_n = table(all_poss_E6_vars$effect_class)[["missense_variant"]]
vars_n = table(E6_vars$effect_class)[["missense_variant"]] + table(E6_vars$effect_class)[["synonymous_variant"]]
missense_n = table(E6_vars$effect_class)[["missense_variant"]]
#print(binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n ))
print("binomial test of conservation of E6 COMPARED TO HPV")
print("binomial test of conservation of E6 COMPARED TO HPV")
print("binomial test of conservation of E6 COMPARED TO HPV")
print("binomial test of conservation of E6 COMPARED TO HPV")
possible_vars_n = table(not_E6_vars$effect_class)[["missense_variant"]] + table(not_E6_vars$effect_class)[["synonymous_variant"]]
possible_missense_n = table(not_E6_vars$effect_class)[["missense_variant"]]
vars_n = table(E6_vars$effect_class)[["missense_variant"]] + table(E6_vars$effect_class)[["synonymous_variant"]]
missense_n = table(E6_vars$effect_class)[["missense_variant"]]
print(binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n ))
print(prop.test(c(missense_n, possible_missense_n), c(vars_n, possible_vars_n)))
print("number of missense")
print(missense_n)
print("number of missense")
print(missense_n)
print(binom.test(x = missense_n , n = data$length[8]  , p=dim(rare_miss_vars)[1]/sum(data$length)))
stat_gene <- c(stat_gene, "E6")
miss_misssynon_ratio_test <- c(miss_misssynon_ratio_test, missense_n/vars_n)
miss_misssynon_ratio_control <- c(miss_misssynon_ratio_control, possible_missense_n/possible_vars_n)
pval_missynon_ratio <- c(pval_missynon_ratio, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )$p.value)
miss_per_length <- c(miss_per_length, missense_n/data$length[8])
pval_miss_per_length <- c(pval_miss_per_length, binom.test(x = missense_n , n = data$length[8]  , p=dim(missense_vars)[1]/sum(data$length))$p.value)
miss_misssynon_ratio_test_CI_low <- c(miss_misssynon_ratio_test_CI_low, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )[[4]][1])
miss_misssynon_ratio_test_CI_high <- c(miss_misssynon_ratio_test_CI_high, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )[[4]][2])
miss_misssynon_ratio_control_CI_low <- c(miss_misssynon_ratio_control_CI_low, binom.test(x = possible_missense_n , n = possible_vars_n , p=possible_missense_n/possible_vars_n )[[4]][1])
miss_misssynon_ratio_control_CI_high <- c(miss_misssynon_ratio_control_CI_high, binom.test(x = possible_missense_n , n = possible_vars_n , p=possible_missense_n/possible_vars_n )[[4]][2])
#### SET UP E5
E5_vars <- subset(my_new_vars, my_new_vars$gene %in% c("E5"))
not_E5_vars <- subset(my_new_vars, !(my_new_vars$gene %in% c("E5")))
all_poss_E5_vars <- subset(all_poss_vars, all_poss_vars$gene %in% c("E5"))
print("binomial test of conservation of E5 COMPARED TO RANDOM")
possible_vars_n = table(all_poss_E5_vars$effect_class)[["missense_variant"]] + table(all_poss_E5_vars$effect_class)[["synonymous_variant"]]
possible_missense_n = table(all_poss_E5_vars$effect_class)[["missense_variant"]]
vars_n = table(E5_vars$effect_class)[["missense_variant"]] + table(E5_vars$effect_class)[["synonymous_variant"]]
missense_n = table(E5_vars$effect_class)[["missense_variant"]]
#print(binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n ))
print("binomial test of conservation of E5 COMPARED TO HPV")
print("binomial test of conservation of E5 COMPARED TO HPV")
print("binomial test of conservation of E5 COMPARED TO HPV")
print("binomial test of conservation of E5 COMPARED TO HPV")
print("binomial test of conservation of E5 COMPARED TO HPV")
possible_vars_n = table(not_E5_vars$effect_class)[["missense_variant"]] + table(not_E5_vars$effect_class)[["synonymous_variant"]]
possible_missense_n = table(not_E5_vars$effect_class)[["missense_variant"]]
vars_n = table(E5_vars$effect_class)[["missense_variant"]] + table(E5_vars$effect_class)[["synonymous_variant"]]
missense_n = table(E5_vars$effect_class)[["missense_variant"]]
print(binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n ))
print(prop.test(c(missense_n, possible_missense_n), c(vars_n, possible_vars_n)))
print("number of missense")
print(missense_n)
print(binom.test(x = missense_n , n = data$length[5]  , p=dim(rare_miss_vars)[1]/sum(data$length)))
stat_gene <- c(stat_gene, "E5")
miss_misssynon_ratio_test <- c(miss_misssynon_ratio_test, missense_n/vars_n)
miss_misssynon_ratio_control <- c(miss_misssynon_ratio_control, possible_missense_n/possible_vars_n)
pval_missynon_ratio <- c(pval_missynon_ratio, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )$p.value)
miss_per_length <- c(miss_per_length, missense_n/data$length[5])
pval_miss_per_length <- c(pval_miss_per_length, binom.test(x = missense_n , n = data$length[5]  , p=dim(missense_vars)[1]/sum(data$length))$p.value)
miss_misssynon_ratio_test_CI_low <- c(miss_misssynon_ratio_test_CI_low, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )[[4]][1])
miss_misssynon_ratio_test_CI_high <- c(miss_misssynon_ratio_test_CI_high, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )[[4]][2])
miss_misssynon_ratio_control_CI_low <- c(miss_misssynon_ratio_control_CI_low, binom.test(x = possible_missense_n , n = possible_vars_n , p=possible_missense_n/possible_vars_n )[[4]][1])
miss_misssynon_ratio_control_CI_high <- c(miss_misssynon_ratio_control_CI_high, binom.test(x = possible_missense_n , n = possible_vars_n , p=possible_missense_n/possible_vars_n )[[4]][2])
#### SET UP E1
E1_vars <- subset(my_new_vars, my_new_vars$gene %in% c("E1"))
not_E1_vars <- subset(my_new_vars, !(my_new_vars$gene %in% c("E1")))
all_poss_E1_vars <- subset(all_poss_vars, all_poss_vars$gene %in% c("E1"))
print("binomial test of conservation of E1 COMPARED TO RANDOM")
possible_vars_n = table(all_poss_E1_vars$effect_class)[["missense_variant"]] + table(all_poss_E1_vars$effect_class)[["synonymous_variant"]]
possible_missense_n = table(all_poss_E1_vars$effect_class)[["missense_variant"]]
vars_n = table(E1_vars$effect_class)[["missense_variant"]] + table(E1_vars$effect_class)[["synonymous_variant"]]
missense_n = table(E1_vars$effect_class)[["missense_variant"]]
#print(binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n ))
print("binomial test of conservation of E1 COMPARED TO HPV")
print("binomial test of conservation of E1 COMPARED TO HPV")
print("binomial test of conservation of E1 COMPARED TO HPV")
print("binomial test of conservation of E1 COMPARED TO HPV")
print("binomial test of conservation of E1 COMPARED TO HPV")
print("binomial test of conservation of E1 COMPARED TO HPV")
possible_vars_n = table(not_E1_vars$effect_class)[["missense_variant"]] + table(not_E1_vars$effect_class)[["synonymous_variant"]]
possible_missense_n = table(not_E1_vars$effect_class)[["missense_variant"]]
vars_n = table(E1_vars$effect_class)[["missense_variant"]] + table(E1_vars$effect_class)[["synonymous_variant"]]
missense_n = table(E1_vars$effect_class)[["missense_variant"]]
print(binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n ))
print(prop.test(c(missense_n, possible_missense_n), c(vars_n, possible_vars_n)))
print("number of missense")
print(missense_n)
print(binom.test(x = missense_n , n = data$length[3]  , p=dim(rare_miss_vars)[1]/sum(data$length)))
stat_gene <- c(stat_gene, "E1")
miss_misssynon_ratio_test <- c(miss_misssynon_ratio_test, missense_n/vars_n)
miss_misssynon_ratio_control <- c(miss_misssynon_ratio_control, possible_missense_n/possible_vars_n)
pval_missynon_ratio <- c(pval_missynon_ratio, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )$p.value)
miss_per_length <- c(miss_per_length, missense_n/data$length[3])
pval_miss_per_length <- c(pval_miss_per_length, binom.test(x = missense_n , n = data$length[3]  , p=dim(missense_vars)[1]/sum(data$length))$p.value)
miss_misssynon_ratio_test_CI_low <- c(miss_misssynon_ratio_test_CI_low, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )[[4]][1])
miss_misssynon_ratio_test_CI_high <- c(miss_misssynon_ratio_test_CI_high, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )[[4]][2])
miss_misssynon_ratio_control_CI_low <- c(miss_misssynon_ratio_control_CI_low, binom.test(x = possible_missense_n , n = possible_vars_n , p=possible_missense_n/possible_vars_n )[[4]][1])
miss_misssynon_ratio_control_CI_high <- c(miss_misssynon_ratio_control_CI_high, binom.test(x = possible_missense_n , n = possible_vars_n , p=possible_missense_n/possible_vars_n )[[4]][2])
#### SET UP L2
L2_vars <- subset(my_new_vars, my_new_vars$gene %in% c("L2"))
not_L2_vars <- subset(my_new_vars, !(my_new_vars$gene %in% c("L2")))
all_poss_L2_vars <- subset(all_poss_vars, all_poss_vars$gene %in% c("L2"))
print("binomial test of conservation of L2 COMPARED TO RANDOM")
possible_vars_n = table(all_poss_L2_vars$effect_class)[["missense_variant"]] + table(all_poss_L2_vars$effect_class)[["synonymous_variant"]]
possible_missense_n = table(all_poss_L2_vars$effect_class)[["missense_variant"]]
vars_n = table(L2_vars$effect_class)[["missense_variant"]] + table(L2_vars$effect_class)[["synonymous_variant"]]
missense_n = table(L2_vars$effect_class)[["missense_variant"]]
#print(binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n ))
print("binomial test of conservation of L2 COMPARED TO HPV")
print("binomial test of conservation of L2 COMPARED TO HPV")
print("binomial test of conservation of L2 COMPARED TO HPV")
print("binomial test of conservation of L2 COMPARED TO HPV")
print("binomial test of conservation of L2 COMPARED TO HPV")
possible_vars_n = table(not_L2_vars$effect_class)[["missense_variant"]] + table(not_L2_vars$effect_class)[["synonymous_variant"]]
possible_missense_n = table(not_L2_vars$effect_class)[["missense_variant"]]
vars_n = table(L2_vars$effect_class)[["missense_variant"]] + table(L2_vars$effect_class)[["synonymous_variant"]]
missense_n = table(L2_vars$effect_class)[["missense_variant"]]
print(binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n ))
print(prop.test(c(missense_n, possible_missense_n), c(vars_n, possible_vars_n)))
print("number of missense")
print(missense_n)
print(binom.test(x = missense_n , n = data$length[7]  , p=dim(rare_miss_vars)[1]/sum(data$length)))
stat_gene <- c(stat_gene, "L2")
miss_misssynon_ratio_test <- c(miss_misssynon_ratio_test, missense_n/vars_n)
miss_misssynon_ratio_control <- c(miss_misssynon_ratio_control, possible_missense_n/possible_vars_n)
pval_missynon_ratio <- c(pval_missynon_ratio, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )$p.value)
miss_per_length <- c(miss_per_length, missense_n/data$length[7])
pval_miss_per_length <- c(pval_miss_per_length, binom.test(x = missense_n , n = data$length[7]  , p=dim(missense_vars)[1]/sum(data$length))$p.value)
miss_misssynon_ratio_test_CI_low <- c(miss_misssynon_ratio_test_CI_low, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )[[4]][1])
miss_misssynon_ratio_test_CI_high <- c(miss_misssynon_ratio_test_CI_high, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )[[4]][2])
miss_misssynon_ratio_control_CI_low <- c(miss_misssynon_ratio_control_CI_low, binom.test(x = possible_missense_n , n = possible_vars_n , p=possible_missense_n/possible_vars_n )[[4]][1])
miss_misssynon_ratio_control_CI_high <- c(miss_misssynon_ratio_control_CI_high, binom.test(x = possible_missense_n , n = possible_vars_n , p=possible_missense_n/possible_vars_n )[[4]][2])
#### SET UP L1
L1_vars <- subset(my_new_vars, my_new_vars$gene %in% c("L1"))
not_L1_vars <- subset(my_new_vars, !(my_new_vars$gene %in% c("L1")))
all_poss_L1_vars <- subset(all_poss_vars, all_poss_vars$gene %in% c("L1"))
print("binomial test of conservation of L1 COMPARED TO RANDOM")
possible_vars_n = table(all_poss_L1_vars$effect_class)[["missense_variant"]] + table(all_poss_L1_vars$effect_class)[["synonymous_variant"]]
possible_missense_n = table(all_poss_L1_vars$effect_class)[["missense_variant"]]
vars_n = table(L1_vars$effect_class)[["missense_variant"]] + table(L1_vars$effect_class)[["synonymous_variant"]]
missense_n = table(L1_vars$effect_class)[["missense_variant"]]
#print(binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n ))
print("binomial test of conservation of L1 COMPARED TO HPV")
print("binomial test of conservation of L1 COMPARED TO HPV")
print("binomial test of conservation of L1 COMPARED TO HPV")
print("binomial test of conservation of L1 COMPARED TO HPV")
print("binomial test of conservation of L1 COMPARED TO HPV")
print("binomial test of conservation of L1 COMPARED TO HPV")
possible_vars_n = table(not_L1_vars$effect_class)[["missense_variant"]] + table(not_L1_vars$effect_class)[["synonymous_variant"]]
possible_missense_n = table(not_L1_vars$effect_class)[["missense_variant"]]
vars_n = table(L1_vars$effect_class)[["missense_variant"]] + table(L1_vars$effect_class)[["synonymous_variant"]]
missense_n = table(L1_vars$effect_class)[["missense_variant"]]
print(binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n ))
print(prop.test(c(missense_n, possible_missense_n), c(vars_n, possible_vars_n)))
print("number of missense")
print(missense_n)
print(binom.test(x = missense_n , n = data$length[6]  , p=dim(missense_vars)[1]/sum(data$length)))
stat_gene <- c(stat_gene, "L1")
miss_misssynon_ratio_test <- c(miss_misssynon_ratio_test, missense_n/vars_n)
miss_misssynon_ratio_control <- c(miss_misssynon_ratio_control, possible_missense_n/possible_vars_n)
pval_missynon_ratio <- c(pval_missynon_ratio, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )$p.value)
miss_per_length <- c(miss_per_length, missense_n/data$length[6])
pval_miss_per_length <- c(pval_miss_per_length, binom.test(x = missense_n , n = data$length[6]  , p=dim(rare_miss_vars)[1]/sum(data$length))$p.value)
miss_misssynon_ratio_test_CI_low <- c(miss_misssynon_ratio_test_CI_low, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )[[4]][1])
miss_misssynon_ratio_test_CI_high <- c(miss_misssynon_ratio_test_CI_high, binom.test(x = missense_n , n = vars_n , p=possible_missense_n/possible_vars_n )[[4]][2])
miss_misssynon_ratio_control_CI_low <- c(miss_misssynon_ratio_control_CI_low, binom.test(x = possible_missense_n , n = possible_vars_n , p=possible_missense_n/possible_vars_n )[[4]][1])
miss_misssynon_ratio_control_CI_high <- c(miss_misssynon_ratio_control_CI_high, binom.test(x = possible_missense_n , n = possible_vars_n , p=possible_missense_n/possible_vars_n )[[4]][2])
##### SUMMARIZE RESULTS
gene_stat_data <- data.frame(stat_gene = stat_gene, miss_misssynon_ratio_test = miss_misssynon_ratio_test, miss_misssynon_ratio_control = miss_misssynon_ratio_control,
pval_missynon_ratio = pval_missynon_ratio, miss_per_length = miss_per_length, pval_miss_per_length= pval_miss_per_length,
ratio_ci_low = miss_misssynon_ratio_test_CI_low,  ratio_ci_high = miss_misssynon_ratio_test_CI_high,
ratio_control_ci_low = miss_misssynon_ratio_control_CI_low,  ratio_control_ci_high = miss_misssynon_ratio_control_CI_high)
gene_stat_data <- gene_stat_data[order(gene_stat_data$stat_gene), ]
gene_stat_data_paired_bar <- data.frame(stat_gene = c(stat_gene, stat_gene),
ratio = c(miss_misssynon_ratio_control, miss_misssynon_ratio_test),
ci_low = c(miss_misssynon_ratio_control_CI_low, miss_misssynon_ratio_test_CI_low),
ci_high = c(miss_misssynon_ratio_control_CI_high, miss_misssynon_ratio_test_CI_high),
test_cont = c(rep("control", time = length(stat_gene)), rep("test", time = length(stat_gene))))
gene_stat_data_paired_bar <- gene_stat_data_paired_bar[order(gene_stat_data_paired_bar$stat_gene), ]
gene_stat_data_per_len <- data.frame(stat_gene = c(stat_gene, "All"), miss_per_length = c(miss_per_length, dim(rare_miss_vars)[1]/sum(data$length)))
gene_stat_data_per_len <- gene_stat_data_per_len[order(gene_stat_data_per_len$stat_gene), ]
library(ggplot2)
library("wesanderson")
miss_sysnon_paired_bar <-  ggplot(gene_stat_data_paired_bar, aes(x = stat_gene, y = ratio, fill = test_cont))+
geom_col(position = "dodge") +
scale_fill_manual(values = c(wes_palette("Moonrise3")[4], wes_palette("Moonrise3")[3]))+
geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width=.2, position = position_dodge(1))+
coord_cartesian(ylim = c(0, 1.15))+
theme_classic()+
geom_signif(y_position = c(1.1), xmin = c(1.8), xmax = c(2.2), annotation = c("***"), tip_length = 0.01)+ #E1E4
geom_signif(y_position = c(0.925), xmin = c(2.8), xmax = c(3.2), annotation = c("***"), tip_length = 0.01)+ #E2
geom_signif(y_position = c(0.85), xmin = c(3.8), xmax = c(4.2), annotation = c("***"), tip_length = 0.01)+ #E5
  geom_signif(y_position = c(0.65), xmin = c(5.8), xmax = c(6.2), annotation = c("***"), tip_length = 0.01)#E7

pdf(file = "all_vars_conservation.pdf", width = 5, height = 3.5)
miss_sysnon_paired_bar
dev.off()
#miss_sysnon_paired_bar





     