library(maftools)
### see https://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html#4_overview_of_the_package
require(TCGAbiolinks)
#require(TCGAbiolinksGUI)
require(DT)
require(dplyr)
require(SummarizedExperiment)
library(here)
library(stringr)
library(seqinr)
library(deconstructSigs)




reference_hpv16 <- read.fasta(file = "reference_hpv16_snpeff.fasta")

##### truncate a sample name to make it in normal TCGA format
truncate <- function(x){
  return(substr(x, 1, 12))
}
##### truncate a sample name to make it in normal TCGA format
truncate_new <- function(x){
  return(str_sub(x, start = -12L, end = -3L))
}
##### 
#vcf_data <- read.table("all_tumor_hpv_snps_all_genes.vcf", header = TRUE, comment.char = "#")
#vcf_data$tumor_id <- truncate_new(as.character(vcf_data$Sample1))

vcf_data <- read.table("all_tumors_hpv_snps_all_genes.vcf", header = FALSE, comment.char = "#")
names(vcf_data) <- c("chromosome", "POS", "blank_A", "REF", "ALT", "blank_B", "filter", "effect_raw", "scheme", "stats", "study_id")
vcf_data$study_id <- str_sub(vcf_data$study_id, start = 1L, end = 10L)


##### VCF normal #######
vcf_normal <- vcf_data[, names(vcf_data) %in% c("POS", "REF")]
vcf_normal <- unique(vcf_normal)
vcf_normal <- vcf_normal[order(vcf_normal$POS), ]

##### VCF Variants #######
vcf_variants <- vcf_data[, names(vcf_data) %in% c("ALT", "POS", "study_id")]

################# get normal reference genome for blank DF

############z
positions_altered <- as.character(1:7906)
all_altered_tumors <- names(table(vcf_data$study_id))
### BLANK DF
blank_df_initial <- data.frame(POS = positions_altered, REF = as.character(reference_hpv16$NC_001526.4))
row.names(blank_df_initial) <- blank_df_initial$POS
for(i in all_altered_tumors){
  blank_df_initial[[i]] <- reference_hpv16$NC_001526.4
}
#blank_df_initial <- blank_df_initial[, !(names(blank_df_initial) %in% c("REF"))]
##########
##########
genotype_df <- blank_df_initial

for(i in 1:length(vcf_variants$tumor_id)){
  genotype_df[as.character(vcf_variants$POS[i]), vcf_variants$tumor_id[i]] <- vcf_variants$ALT[i]
}

test <- genotype_df[, !(names(genotype_df) %in% c("POS"))]

genotopye_rows <- data.frame(t(test))
unique_genotypes<- unique(genotopye_rows)

sequences <- c()
for(i in names(test)){
  sequences <- c(sequences, str_c(test[[i]], collapse = ""))
}
sequences <- as.list(sequences)
write.fasta(sequences, names(test), file.out = "HPV16_genotypes_UNCseq.fasta", open = "w", nbchar = 60, as.string = FALSE)

#### SECTION TO DO WITH ONLY CODING VARIANTS
var_effect <- c()
viral_gene <- c()
effect_coding <- c()
for(i in 1:length(vcf_data$chromosome)){
  var_effect <- c(var_effect, str_split(vcf_data$effect_raw[[i]][1], "[|]")[[1]][2])
  viral_gene <- c(viral_gene, str_split(vcf_data$effect_raw[[i]][1], "[|]")[[1]][4])
  effect_coding <- c(effect_coding, str_split(vcf_data$effect_raw[[i]][1], "[|]")[[1]][11])
}

vcf_data_with_effect <- vcf_data
vcf_data_with_effect$effect <- var_effect
vcf_data_with_effect$effect_coding <- effect_coding
vcf_data_with_effect$viral_gene <- viral_gene

vcf_data_with_effect$unique_mut_id <- paste(vcf_data_with_effect$viral_gene, vcf_data_with_effect$effect_coding, sep = "_")

library(stringr)
#stat <- c()
#for(i in vcf_data_with_effect$study_id){
#  holder <- str_split(i, " ")[[1]][1]
#  stat <- c(stat, holder)
#}
#vcf_data_with_effect$stat <- stat

depth <- c()
vaf <- c()

for(ii in as.character(vcf_data_with_effect$stats)){
  depth <- c(depth, str_split(ii, "[:]")[[1]][4])
  holder <- str_split(ii, "[:]")[[1]][7]
  vaf<- c(vaf, as.numeric(str_split(holder, "[%]")[[1]][1]))
}

vcf_data_with_effect$deth <- depth
vcf_data_with_effect$vaf <- vaf


hpv16_ref_seq_vect <- toupper(reference_hpv16$NC_001526.4)

hpv16_ref_seq_vect <- as.vector(hpv16_ref_seq_vect)

hpv16_ref_seq_vect_times3 <- c(hpv16_ref_seq_vect, hpv16_ref_seq_vect, hpv16_ref_seq_vect)

preceed_POS <- c()
follows_POS <- c()
for(i in vcf_data_with_effect$POS){
  
  preceed_POS_position <-  (i - 1) + 7906
  preceed_POS <- c(preceed_POS, hpv16_ref_seq_vect_times3[preceed_POS_position])
  
  follows_POS_position <-  (i + 1) + 7906
  follows_POS <- c(follows_POS, hpv16_ref_seq_vect_times3[follows_POS_position])
  
}

vcf_data_with_effect$preceed_POS <- preceed_POS
vcf_data_with_effect$follows_POS <- follows_POS


sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref, 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt")

trinuc96_ref <- names(sigs.input)

bases <- c("A", "T", "C", "G")
complements <- c("T", "A", "G", "C")
names(complements) <- bases



trinuc_alt <-c()
trinuc_alt_rc <- c()
for(i in 1:length(vcf_data_with_effect$POS)){
  trinuc_alt <-c(trinuc_alt, paste(vcf_data_with_effect$preceed_POS[i], "[", vcf_data_with_effect$REF[i], ">", vcf_data_with_effect$ALT[i], "]", vcf_data_with_effect$follows_POS[i], sep = ""))
  trinuc_alt_rc <-c(trinuc_alt_rc, paste(complements[vcf_data_with_effect$follows_POS[i]], "[", complements[vcf_data_with_effect$REF[i]], ">", complements[vcf_data_with_effect$ALT[i]], "]", complements[vcf_data_with_effect$preceed_POS[i]], sep = ""))
}

trinuc96_alts <- c()

for(i in 1:length(vcf_data_with_effect$POS)){
  if(trinuc_alt[i] %in% trinuc96_ref){
    trinuc96_alts <- c(trinuc96_alts, trinuc_alt[i])
  }else{
    trinuc96_alts <- c(trinuc96_alts, trinuc_alt_rc[i])
  }
}

vcf_data_with_effect$trinuc96_alts <- trinuc96_alts

vcf_data_with_effect_all <- vcf_data_with_effect
save(vcf_data_with_effect_all, file = "vcf_data_with_effect_all.RData")

###
### FILTER SUBCLONAL
###
###

vcf_data_with_effect_subclonal <- subset(vcf_data_with_effect, vcf_data_with_effect$vaf <= 95)
vcf_data_with_effect_subclonal <- subset(vcf_data_with_effect_subclonal, vcf_data_with_effect_subclonal$vaf > 5)


vcf_data_with_effect <- subset(vcf_data_with_effect, vcf_data_with_effect$vaf >= 90)

L2_269_muts_effects <- c("L2_p.Ser269Pro", "L2_p.Ser269Thr", "L2_p.Ser269Tyr", "L2_p.Ser269Ala")
L2_269_vars <- subset(vcf_data_with_effect, vcf_data_with_effect$unique_mut_id %in% L2_269_muts_effects)
save(L2_269_vars, file = "L2_269_vars.RData")

###
### FILTER SYNONOMOS
###
###

vcf_data_with_effect <- subset(vcf_data_with_effect, !(vcf_data_with_effect$effect %in% c("synonymous_variant")) )
#### filter frequent vars
vcf_data_with_effect$var_id_prot_nuc <- paste(vcf_data_with_effect$unique_mut_id, vcf_data_with_effect$trinuc96_alts, sep = "_")
table_var_ids <- table(vcf_data_with_effect$var_id_prot_nuc)
vcf_data_with_effect$var_freq <- table_var_ids[vcf_data_with_effect$var_id_prot_nuc]

vcf_data_with_effect_rare <- subset(vcf_data_with_effect, vcf_data_with_effect$var_freq <= 20)
vcf_data_with_effect_common <- subset(vcf_data_with_effect, vcf_data_with_effect$var_freq > 20)

vcf_data_with_effect_coding_rare <- vcf_data_with_effect_rare 
save(vcf_data_with_effect_coding_rare, file = "vcf_data_with_effect_coding_rare.RData")

vcf_data_with_effect_coding_common <- vcf_data_with_effect_common
save(vcf_data_with_effect_coding_common, file = "vcf_data_with_effect_coding_common.RData")


my_tric_data <- table(c(vcf_data_with_effect$trinuc96_alts, names(sigs.input)), c(vcf_data_with_effect$study_id, rep("TEST", times = 96)))
my_tric_data <- t(my_tric_data)
my_tric_data <-my_tric_data[!(row.names(my_tric_data) %in% c("TEST")), ]
mytest <- rbind(sigs.input, colSums(my_tric_data))

all_tumors_signatures = whichSignatures(tumor.ref = mytest, 
                           signatures.ref = signatures.cosmic, 
                           sample.id = 3, 
                           contexts.needed = TRUE,
                           tri.counts.method = 'default')

plotSignatures(all_tumors_signatures)
makePie(all_tumors_signatures) 





### RARE
my_tric_data <- table(c(vcf_data_with_effect_rare$trinuc96_alts, names(sigs.input)), c(vcf_data_with_effect_rare$study_id, rep("TEST", times = 96)))
my_tric_data <- t(my_tric_data)
my_tric_data <-my_tric_data[!(row.names(my_tric_data) %in% c("TEST")), ]
mytest <- rbind(sigs.input, colSums(my_tric_data))

all_tumors_signatures_rare = whichSignatures(tumor.ref = mytest, 
                                        signatures.ref = signatures.cosmic, 
                                        sample.id = 3, 
                                        contexts.needed = TRUE,
                                        tri.counts.method = 'default')

plotSignatures(all_tumors_signatures_rare)
makePie(all_tumors_signatures_rare) 
pdf(file = "all_tumors_signatures_rare_SIG_PIE.pdf", height = 4.5, width = 4.5)
print(makePie(all_tumors_signatures_rare))
dev.off()

### COMMON
my_tric_data <- table(c(vcf_data_with_effect_common$trinuc96_alts, names(sigs.input)), c(vcf_data_with_effect_common$study_id, rep("TEST", times = 96)))
my_tric_data <- t(my_tric_data)
my_tric_data <-my_tric_data[!(row.names(my_tric_data) %in% c("TEST")), ]
mytest <- rbind(sigs.input, colSums(my_tric_data))

all_tumors_signatures_common = whichSignatures(tumor.ref = mytest, 
                                             signatures.ref = signatures.cosmic, 
                                             sample.id = 3, 
                                             contexts.needed = TRUE,
                                             tri.counts.method = 'default')

plotSignatures(all_tumors_signatures_common)
makePie(all_tumors_signatures_common) 



###
### BRIN IN N CODIN POLY DATA
###

load("n_clonal_missense.RData")
load("far_clade_parsimony.RData")
med_val <- median(n_clonal_missense_vars$n_clonal_missense)
high_n_coding <- subset(n_clonal_missense_vars, n_clonal_missense_vars$n_clonal_missense > 9)
low_n_coding <- subset(n_clonal_missense_vars, n_clonal_missense_vars$n_clonal_missense <= 9)




high_n_coding_tumor_vars <- subset(vcf_data_with_effect_rare, (vcf_data_with_effect_rare$study_id %in% high_n_coding$study_id))


my_tric_data <- table(c(vcf_data_with_effect_rare$trinuc96_alts, names(sigs.input)), c(vcf_data_with_effect_rare$study_id, rep("TEST", times = 96)))
my_tric_data <- t(my_tric_data)
my_tric_data <-my_tric_data[!(row.names(my_tric_data) %in% c("TEST")), ]
my_tric_data <-my_tric_data[row.names(my_tric_data) %in% high_n_coding$study_id, ]
mytest <- rbind(sigs.input, colSums(my_tric_data))

high_n_coding_tumors_signatures = whichSignatures(tumor.ref = mytest, 
                                        signatures.ref = signatures.cosmic, 
                                        sample.id = 3, 
                                        contexts.needed = TRUE,
                                        tri.counts.method = 'default')

plotSignatures(high_n_coding_tumors_signatures)
makePie(high_n_coding_tumors_signatures)
pdf(file = "high_n_coding_tumors_signatures_uncommon_SIG_PIE.pdf", height = 4.5, width = 4.5)
print(makePie(high_n_coding_tumors_signatures))
dev.off()



low_n_coding_tumor_vars <- subset(vcf_data_with_effect_rare, (vcf_data_with_effect_rare$study_id %in% low_n_coding$study_id))


my_tric_data <- table(c(vcf_data_with_effect_rare$trinuc96_alts, names(sigs.input)), c(vcf_data_with_effect_rare$study_id, rep("TEST", times = 96)))
my_tric_data <- t(my_tric_data)
my_tric_data <-my_tric_data[!(row.names(my_tric_data) %in% c("TEST")), ]
my_tric_data <-my_tric_data[row.names(my_tric_data) %in% low_n_coding$study_id, ]
mytest <- rbind(sigs.input, colSums(my_tric_data))

low_n_coding_tumors_signatures = whichSignatures(tumor.ref = mytest, 
                                                  signatures.ref = signatures.cosmic, 
                                                  sample.id = 3, 
                                                  contexts.needed = TRUE,
                                                  tri.counts.method = 'default')

plotSignatures(low_n_coding_tumors_signatures)
makePie(low_n_coding_tumors_signatures)
pdf(file = "low_n_coding_tumors_signatures_uncommon_SIG_PIE.pdf", height = 4.5, width = 4.5)
print(makePie(low_n_coding_tumors_signatures))
dev.off()




apobec_data <- vcf_data_with_effect_all 
#apobec_data <- subset(apobec_data, apobec_data$vaf > 5)
apobec_data$subclonal <- apobec_data$vaf <= 95
apobec_data$composite <- paste(apobec_data$REF, apobec_data$ALT, sep = "")
apobec_data$apobec <- str_sub(apobec_data$trinuc96_alts, start = 1L, end = 7L) %in% c("T[C>G]A", "T[C>G]C", "T[C>G]T")  ### simplistic Sig 13 APOBEC3B
apobec_data$apobec_2 <- str_sub(apobec_data$trinuc96_alts, start = 1L, end = 7L) %in% c("T[C>T]A", "T[C>T]C", "T[C>T]T")  ### simplistic Sig 2 APOBEC3B

apobec_status  <- c()
trics <- str_sub(apobec_data$trinuc96_alts, start = 1L, end = 7L)
for(i in trics){
  if(i %in%  c("T[C>G]A", "T[C>G]C", "T[C>G]T")){
    apobec_status  <- c(apobec_status, "Sig13")
  }else if(i %in% c("T[C>T]A", "T[C>T]C", "T[C>T]T")){
    apobec_status  <- c(apobec_status, "Sig2")
  }else{
    apobec_status  <- c(apobec_status, "Other")
  }
}

apobec_data$apobec_status <- apobec_status


apobec_data$n_coding_status <- apobec_data$study_id %in% high_n_coding$sample_id
apobec_data$subclonal_apobec <- apobec_data$subclonal & apobec_data$apobec

high_var_tumors_vars <- subset(apobec_data, apobec_data$n_coding_status == TRUE)
low_var_tumors_vars <- subset(apobec_data, apobec_data$n_coding_status == FALSE)

clonal_vars <- subset(apobec_data, apobec_data$subclonal == FALSE)
subclonal_vars <- subset(apobec_data, apobec_data$subclonal == TRUE)
subclonal_apo_vars <- subset(apobec_data, apobec_data$subclonal_apobec == TRUE)


n_subclonal_apo_vars <- data.frame(table(subclonal_apo_vars$study_id))
names(n_subclonal_apo_vars) <- c("study_id", "n_subclonal_apo")
n_subclonal_apo_vars$n_coding_status <- n_subclonal_apo_vars$study_id %in% high_n_coding$sample_id

n_subclonal_apo_vars_high_coding <- subset(n_subclonal_apo_vars, n_subclonal_apo_vars$n_coding_status == TRUE)
n_subclonal_apo_vars_low_coding <- subset(n_subclonal_apo_vars, n_subclonal_apo_vars$n_coding_status == FALSE)

subclonal_apobec_data <- subset(apobec_data, apobec_data$subclonal == TRUE)

n_sub_apo <- data.frame(table(subclonal_apobec_data$study_id))

n_sub_apo$high_n_poly <- n_sub_apo$Var1 %in% high_n_coding_tumor_vars$study_id

n_sub_apo_high_n_coding <- subset(n_sub_apo, n_sub_apo$high_n_poly == TRUE)
n_sub_apo_low_n_coding <- subset(n_sub_apo, n_sub_apo$high_n_poly == FALSE)

print(t.test(n_sub_apo_high_n_coding$Freq, n_sub_apo_low_n_coding$Freq))


prop.test(c(table(subclonal_vars$apobec)[2], table(clonal_vars$apobec)[2]), c(sum(table(subclonal_vars$apobec)), sum(table(clonal_vars$apobec))))

prop.test(c(table(subclonal_vars$apobec_2  )[2], table(clonal_vars$apobec_2)[2]), c(sum(table(subclonal_vars$apobec_2)), sum(table(clonal_vars$apobec_2))))

apobec_data_lt95 <- apobec_data
save(apobec_data_lt95, file = "all_vars_lt95_apobec_annotated.RData")

apobec_data <- subset(apobec_data, apobec_data$vaf > 5)
library(ggplot2)
library(wesanderson)

#apobec_data <- subset(apobec_data, apobec_data$effect != "synonymous_variant")
test <-ggplot(apobec_data, aes(x=subclonal, fill=apobec))
mybar <- test+ geom_bar(position = "fill") + theme_classic() + scale_fill_manual(values = c(wes_palette("Moonrise2"),wes_palette("Moonrise3")[1:3], wes_palette("Moonrise3")[5]))+
  labs(y="APOBEC Mutations (Proportion)", x = "Sub-clonal")

pdf(file = "bar_abopec_subclonal.pdf", width = 3.5, height = 4)
mybar
dev.off()



#apobec_data <- subset(apobec_data, apobec_data$effect != "synonymous_variant")
test <-ggplot(apobec_data, aes(x=subclonal, fill=apobec_2))
mybar_2 <- test+ geom_bar(position = "fill") + theme_classic() + scale_fill_manual(values = c(wes_palette("Moonrise2")[1],wes_palette("Moonrise2")[3]))+
  labs(y="APOBEC2 Mutations (Proportion)", x = "Sub-clonal")

pdf(file = "bar_abopec_2_subclonal.pdf", width = 3.5, height = 4)
mybar_2
dev.off()




#apobec_data <- subset(apobec_data, apobec_data$effect != "synonymous_variant")
test <-ggplot(apobec_data, aes(x=subclonal, fill=apobec_status))
mybar_2 <- test+ geom_bar(position = "fill") + theme_classic() + scale_fill_manual(values = c(wes_palette("Moonrise2")[1],wes_palette("Moonrise2")[2], wes_palette("Moonrise2")[3]))+
  labs(y="APOBEC Mutations (Proportion)", x = "Sub-clonal")

pdf(file = "bar_abopec_status_subclonal.pdf", width = 3.5, height = 4)
mybar_2
dev.off()


















