library(maftools)
unc_maf = read.maf(maf = "grouped_maf.maf", removeDuplicatedVariants = FALSE)
unc_df <- unc_maf@data
unc_gene_summary_df <- getGeneSummary(unc_maf)
#plot_me <- unc_gene_summary_df$Hugo_Symbol[1:50]
#plot_me <- c(plot_me, "TP53")

hpv_maf <- subsetMaf(maf = unc_maf, tsb =  clin_data$study_id , mafObj = TRUE)
hpv_gene_summary_df <- getGeneSummary(hpv_maf)
plot_me <- unc_gene_summary_df$Hugo_Symbol[1:25]
plot_me <- c(plot_me, "TP53")
plot_me_HLA <- c("HLA-A", "HLA-B", "HLA-C")


load("merged_hpv16_p16_BOT_TON.RData")
clin_data <- merged_hpv16_p16_BOT_TON
load("far_clade_parsimony.RData")

clin_data$far_clade <- clin_data$study_id %in% far_clade_parsimony 



near_df <- subset(clin_data, clin_data$far_clade == FALSE)
far_df <- subset(clin_data, clin_data$far_clade == TRUE)

near_maf <- subsetMaf(maf = unc_maf, tsb =  near_df$study_id , mafObj = TRUE)
near_gene_summary_df <- getGeneSummary(near_maf)

far_maf <- subsetMaf(maf = unc_maf, tsb = far_df$study_id , mafObj = TRUE)
far_gene_summary_df <- getGeneSummary(far_maf)

#plot_me <- c(near_gene_summary_df$Hugo_Symbol[1:10], far_gene_summary_df$Hugo_Symbol[1:10], "TP53")

library(RColorBrewer)
library(RColorBrewer)
col = RColorBrewer::brewer.pal(n = 8, name = 'Set1')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')
bg <- rgb(0,0,0, 0.0)

#quartz()
#pdf(file = "coOncoplot_het.pdf", height = 12, width = 8)
coOncoplot(near_maf, far_maf, m1Name = "A1", m2Name = "Divergent", genes = plot_me_HLA, colors = col, bgCol = "lightgrey", borderCol = bg, showSampleNames = TRUE)
#dev.off()

print(mafCompare(near_maf, far_maf, m1Name = "near_A1", m2Name = "divergent_A1", useCNV = FALSE, minMut = 1))
comparison_df <- mafCompare(near_maf, far_maf, m1Name = "near_A1", m2Name = "divergent_A1", useCNV = FALSE, minMut = 1)
write.table(comparison_df$results, file = "gene_comparison_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)


