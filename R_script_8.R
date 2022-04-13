library(stringr)
library(here)
library(dplyr)
library(pheatmap)
###bring in sample data
load("merged_genomics_genotpye_prediction.RData")
sample_data <- merged_genomics_genotype_prediction
#load("clinical_out_LT30.RData")
#merged_genomics_genotype_prediction <- subset(merged_genomics_genotype_prediction, merged_genomics_genotype_prediction$study_id %in% clinical_out_LT30$study_id)
merged_genomics_genotype_prediction$file_string <- str_sub(merged_genomics_genotype_prediction$sample.y, start = 1L, end = -7L)


### bring in genomic data
load("raw_counts.RData")
load("abundance.RData")
#raw_counts <- abundance
names(raw_counts) <- str_sub(names(raw_counts), start = 1L, end = -7L)
load("hpv16_raw_reads.RData")

hs_hpv_counts <- rbind(raw_counts, hpv16_raw_reads)
hs_hpv_counts <- hs_hpv_counts[ , names(hs_hpv_counts) %in% merged_genomics_genotype_prediction$file_string]
hs_hpv_counts <- hs_hpv_counts[ ,merged_genomics_genotype_prediction$file_string]

name_dictionary <- merged_genomics_genotype_prediction$study_id_dot
names(name_dictionary) <- merged_genomics_genotype_prediction$file_string

new_names <- c() 
for(i in names(hs_hpv_counts)){
  new_names <- c(new_names, name_dictionary[[i]])
}
names(hs_hpv_counts) <- new_names

total_mapped_dictionary <- colSums(hs_hpv_counts)
names(total_mapped_dictionary) <- new_names


hpv16_tx_reads <- hpv16_raw_reads
hpv16_tx_reads <- hpv16_tx_reads[ , names(hpv16_tx_reads) %in% merged_genomics_genotype_prediction$file_string]
hpv16_tx_reads <- hpv16_tx_reads[ , merged_genomics_genotype_prediction$file_string]
new_names <- c() 
for(i in names(hpv16_tx_reads)){
  new_names <- c(new_names, name_dictionary[[i]])
}
names(hpv16_tx_reads) <- new_names
hpv16_tx_reads <- hpv16_tx_reads * 1000000
hpv16_tx_reads <- hpv16_tx_reads +1


hpv16_tx_cpm <- hpv16_tx_reads
for(i in names(hpv16_tx_cpm)){
  #print(total_mapped_dictionary[i])
  hpv16_tx_cpm[[i]] <- hpv16_tx_cpm[[i]] / total_mapped_dictionary[i]
}

hpv16_tx_Lcpm_plot <- log(hpv16_tx_cpm, base = 2)


vifi_cluster_data <- read.table("all_tumor_clusters_TUMOR_ID.txt", header = FALSE, sep = "")
#vifi_cluster_data <- subset(vifi_cluster_data, vifi_cluster_data$V5 > 1)
vifi_cluster_data <- subset(vifi_cluster_data, vifi_cluster_data$V4 > 25)
sum_clusters <- data.frame(table(vifi_cluster_data$V7))
sum_clusters <- subset(sum_clusters, sum_clusters$Var1 %in% merged_genomics_genotype_prediction$file_string)
sum_clusters$sample_id <- name_dictionary[as.character(sum_clusters$Var1)]





#sum_clusters <- subset(sum_clusters, sum_clusters$sample_id %in% hpv16_pos_tumors$sample_id)
#sum_clusters <- subset(sum_clusters, sum_clusters$Freq > 1)






t_scaled <- t(hpv16_tx_cpm)
t_scaled <- data.frame(t_scaled)
#E1E2 <- (t_scaled$E2 + t_scaled$E5 )/2
E1E2 <- (t_scaled$E2 + t_scaled$E5 )/2
E6E7 <- (t_scaled$E6 + t_scaled$E7)/2
ratio <- log(E6E7/E1E2, base = 2)
all_hpv <- log(rowSums(t_scaled), base = 2)
E5 <- t_scaled$E5
ratio_data <- data.frame(sample_id = row.names(t_scaled), ratio = ratio, E5 = E5, all_hpv = all_hpv)
ratio_data <- ratio_data[order(ratio_data$ratio), ]
ratio_data <- subset(ratio_data, ratio_data$sample_id %in% names(hpv16_tx_Lcpm_plot))


hpv16_tx_Lcpm_plot <- hpv16_tx_Lcpm_plot[, ratio_data$sample_id]

integrated_annotaion_df <- data.frame(integrated = as.numeric(names(hpv16_tx_Lcpm_plot) %in% sum_clusters$sample_id), E6E7toE2E5ratio = ratio_data$ratio, all_hpv = all_hpv)
far_clade <- subset(merged_genomics_genotype_prediction, merged_genomics_genotype_prediction$prediction == TRUE)
row.names(integrated_annotaion_df) <- names(hpv16_tx_Lcpm_plot)
integrated_annotaion_df$integrated_ratio <- as.factor(integrated_annotaion_df$E6E7toE2E5ratio > -0.304)
integrated_annotaion_df$far_clade <- as.factor(row.names(integrated_annotaion_df) %in% far_clade$study_id_dot)

integrated_annotaion_df <- integrated_annotaion_df[, c("integrated",  "E6E7toE2E5ratio","integrated_ratio", "far_clade")]

pheatmap(as.matrix(hpv16_tx_Lcpm_plot), annotation_col = integrated_annotaion_df, cluster_cols = FALSE, cluster_rows = FALSE, clustering_distance_cols = "euclidean", clustering_method = "ward.D2")

chisq.test(integrated_annotaion_df$far_clade, integrated_annotaion_df$integrated)
chisq.test(integrated_annotaion_df$far_clade, integrated_annotaion_df$integrated_ratio)


integrated_annotaion_df$study_id_dot <- row.names(integrated_annotaion_df)
integrated_annotaion_df <- merge(integrated_annotaion_df, merged_genomics_genotype_prediction, by = "study_id_dot")

annot_far <- subset(integrated_annotaion_df, integrated_annotaion_df$integrated_ratio == TRUE)
annot_near <- subset(integrated_annotaion_df, integrated_annotaion_df$integrated_ratio == FALSE)
wilcox.test(annot_far$E6E7toE2E5ratio, annot_near$E6E7toE2E5ratio)
wilcox.test(annot_far$NFkB_TPS, annot_near$NFkB_TPS)

save(integrated_annotaion_df, file = "integrated_annotaion_df.RData")

clade_data <- integrated_annotaion_df[, c("study_id_dot","far_clade", "integrated")]

library(reshape2)


box_data <- hpv16_tx_Lcpm_plot
box_data <- data.frame(t(box_data))
box_data$study_id_dot <- row.names(box_data)
box_data <- merge(box_data, clade_data,  by = "study_id_dot")


melted_box <- melt(box_data, id.vars = c("far_clade", "study_id_dot", "integrated"))

names(melted_box) <- c("far_clade","study_id_dot","integrated", "gene", "expression")
melted_box$gene <- as.factor(melted_box$gene)

melted_box$integrated <- as.logical(melted_box$integrated)

library(ggplot2)
library(wesanderson)

my_box <- ggplot(melted_box, aes(x = gene, y = expression, fill = integrated))+
  geom_boxplot()

far_box_data <- subset(box_data, box_data$far_clade == TRUE)
near_box_data <- subset(box_data, box_data$far_clade == FALSE)

wilcox.test(far_box_data$E5, near_box_data$E5)

