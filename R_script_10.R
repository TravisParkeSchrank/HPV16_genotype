library(here)
library(stringr)
library(ggplot2)

hla_data <- read.table("out_hla_data.csv", header = TRUE,sep = ",")
hla_data  <- subset(hla_data, hla_data$sample_TN %in% c("T"))
hla_data$hla_count <- apply(hla_data, 1, function(x){length(unique(x))}) - 2
hpv_pos <- read.csv("104_hpv16_pos_cases.csv")
hla_data$hpv_status <- hla_data$sample_ID %in% hpv_pos$sample
hla_data <- subset(hla_data, hla_data$hpv_status == TRUE)


#load("netMHCpan_reults.RData")



#pept_len <- c()
#for(i in netMHCpan_results$Peptide){
#  pept_len <- c(pept_len, str_length(i))
#}

#netMHCpan_results$pept_len <- pept_len

#netMHCpan_results$start_position <- netMHCpan_results$Position
#netMHCpan_results$end_position  <- netMHCpan_results$Position + (netMHCpan_results$pept_len - 1)


#netMHCpan_results_len <- netMHCpan_results
#save(netMHCpan_results_len, file = "netMHCpan_reults_len.RData")


load("netMHCpan_reults_len.RData")
netMHC_pan_cases <- subset(netMHCpan_results_len, netMHCpan_results_len$sample_ID %in% hpv_pos$sample)
load("far_clade_parsimony.RData")

netMHC_pan_cases$far_clade <- netMHC_pan_cases$sample_ID %in% far_clade_parsimony

nM_cutoff = 325
results_HA <- subset(netMHC_pan_cases, netMHC_pan_cases$nM_aff < nM_cutoff)

hold_results <- subset(netMHC_pan_cases, netMHC_pan_cases$gene %in% c("L2"))
hold_results_HA <- subset(hold_results, hold_results$nM_aff < nM_cutoff)
L2_plot <- ggplot(hold_results_HA, aes(x=Position, fill = far_clade, color = far_clade))+ geom_histogram(position="identity", alpha= 0.5)

hold_results <- subset(netMHC_pan_cases, netMHC_pan_cases$gene %in% c("L1"))
hold_results_HA <- subset(hold_results, hold_results$nM_aff < nM_cutoff)
L1_plot <- ggplot(hold_results_HA, aes(x=Position, fill = far_clade, color = far_clade))+ geom_histogram(position="identity", alpha= 0.5)

hold_results <- subset(netMHC_pan_cases, netMHC_pan_cases$gene %in% c("E1"))
hold_results_HA <- subset(hold_results, hold_results$nM_aff < nM_cutoff)
E1_plot <- ggplot(hold_results_HA, aes(x=Position, fill = far_clade, color = far_clade))+ geom_histogram(position="identity", alpha= 0.5)

hold_results <- subset(netMHC_pan_cases, netMHC_pan_cases$gene %in% c("E2"))
hold_results_HA <- subset(hold_results, hold_results$nM_aff < nM_cutoff)
E2_plot <- ggplot(hold_results_HA, aes(x=Position, fill = far_clade, color = far_clade))+ geom_histogram(position="identity", alpha= 0.5)

hold_results <- subset(netMHC_pan_cases, netMHC_pan_cases$gene %in% c("E5"))
hold_results_HA <- subset(hold_results, hold_results$nM_aff < nM_cutoff)
E5_plot <- ggplot(hold_results_HA, aes(x=Position, fill = far_clade, color = far_clade))+ geom_histogram(position="identity", alpha= 0.5)

hold_results <- subset(netMHC_pan_cases, netMHC_pan_cases$gene %in% c("E1E4"))
hold_results_HA <- subset(hold_results, hold_results$nM_aff < nM_cutoff)
E1E4_plot <- ggplot(hold_results_HA, aes(x=Position, fill = far_clade, color = far_clade))+ geom_histogram(position="identity", alpha= 0.5)

hold_results <- subset(netMHC_pan_cases, netMHC_pan_cases$gene %in% c("E6"))
hold_results_HA <- subset(hold_results, hold_results$nM_aff < nM_cutoff)
E6_plot <- ggplot(hold_results_HA, aes(x=Position, fill = far_clade, color = far_clade))+ geom_histogram(position="identity", alpha= 0.5)

hold_results <- subset(netMHC_pan_cases, netMHC_pan_cases$gene %in% c("E7"))
hold_results_HA <- subset(hold_results, hold_results$nM_aff < nM_cutoff)
E7_plot <- ggplot(hold_results_HA, aes(x=Position, fill = far_clade, color = far_clade))+ geom_histogram(position="identity", alpha= 0.5)

n_HA_pept_data <- data.frame(table(results_HA$sample_ID))
names(n_HA_pept_data) <- c("study_id", "n_HA_pept")
save(n_HA_pept_data, file = "n_HA_pept_data.RData")

