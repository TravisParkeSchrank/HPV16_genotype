library(ggbio)
library(GenomicRanges)
library(stringr)



in_data <- read.table("idx_stats_summary.csv", header = TRUE, sep = ",")

total_align_list <- c()
prop_test_p <- c()
log10_prop_other_hpv <- c()
log10_prop_hpv16 <- c()
pvals <- c()
log10_prop_hpv16_unmapped <- c()



for(i in 1:length(in_data$sample)){
	total_alignmnets = sum(in_data$hpv16_aligns[i], in_data$human_aligns[i], in_data$other_hpv_aligns[i])
	total_align_list <- c(total_align_list, total_alignmnets)
	log10_prop_other_hpv <- c(log10_prop_other_hpv, log((in_data$other_hpvs_aligns[i]+1)/total_alignmnets, base = 10))
	log10_prop_hpv16 <- c(log10_prop_hpv16, log((in_data$hpv16_aligns[i]+1)/total_alignmnets, base = 10))
	hold_p <- prop.test(c(in_data$hpv16_aligns[i], in_data$other_hpvs_aligns[i]), c(total_alignmnets, total_alignmnets))[3][[1]]
	#print(prop.test(c(in_data$hpv16_aligns[i], in_data$other_hpvs_aligns[i]), c(total_alignmnets, total_alignmnets)))
	pvals <- c(pvals, hold_p)
	
	log10_prop_hpv16_unmapped <- c(log10_prop_hpv16_unmapped, log(in_data$hpv16_not_aligned[i]/in_data$hpv16_aligns[i], base = 10))
	
	
}

in_data$total_alignmnets <- total_align_list
in_data$log10_prop_other_hpv <- log10_prop_other_hpv
in_data$log10_prop_hpv16 <- log10_prop_hpv16
in_data$pvals <- pvals
in_data$log10_prop_hpv16_unmapped <- log10_prop_hpv16_unmapped


plot_data <- data.frame(

log_ratio = c(log10_prop_hpv16, log10_prop_other_hpv),
group = c(rep("hpv16", times = length(log10_prop_hpv16)), rep("other_hpv", times = length(log10_prop_other_hpv)) )

)

quartz()
plot(as.factor(plot_data$group), plot_data$log_ratio)



#load("/Volumes/SDSSD_D/ViFi_integration_HPV16_OPSCC/out_tumor_hpv_integr_data_lisle.RData")


#new_samples   <- c()

#for (i in 1:length(out_tumor_hpv_integr_data$sample)){
#new_samples   <- c(new_samples, str_split(out_tumor_hpv_integr_data$sample[i], "_")[[1]][1])
#	
#}

#out_tumor_hpv_integr_data$sample <- new_samples

#merged <- merge(in_data, out_tumor_hpv_integr_data, by = "sample")

library(ggplot2)
library(wesanderson)




box_data <- plot_data 


bp <-ggplot(box_data, aes(x=group, y=log_ratio, fill=group)) +
   geom_boxplot(outlier.shape=NA)+
   geom_jitter(width = 0.15, pch = 16)+
   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
   scale_fill_manual(values=c(wes_palette(n=4, name="Royal1")[4],wes_palette(n=4, name="Royal1")[1])) 

pdf(file="AAA_HPV16_Positivity_BOX_TEST.pdf", width = 4, height = 3)
print(bp)
dev.off()



##### PLOT THE RELATIVE RATIO OF HPV16 READS to NONSPECIFIC HPV READS

hpv_to_other_logR <- log(((in_data$hpv16_aligns + 1)/(in_data$other_hpvs_aligns+1)), base = 10)


in_data$hpv_to_other_logR <- hpv_to_other_logR

hist_viral_ratio <- ggplot(in_data, aes(x=hpv_to_other_logR ))+
geom_histogram(binwidth = 0.15, fill = "white", color = "black") #fill = wes_palette(n=4, name="Royal1")[4])

pdf(file="AAA_HPV16_Positivity_HIST_TEST.pdf", width = 2.5, height = 2.5)
print(hist_viral_ratio)
dev.off()

low_hpv16_other_ratio <- subset(in_data, in_data$hpv_to_other_logR <3)

##### PLOT THE RELATIVE HPV COPY NUMBER

 

hist_viral_cn <- ggplot(in_data, aes(x=log10_prop_hpv16))+
geom_histogram(binwidth = 0.15, fill = "white", color = "black") #fill = wes_palette(n=4, name="Royal1")[4])

pdf(file="AAA_HPV16_COPY_NUMBER_HIST_TEST.pdf", width = 2.5, height = 2.5)
print(hist_viral_cn)
dev.off()





viral_cn_stats <- in_data

save(viral_cn_stats, file = "viral_cn_stats.RData")


hpv16_pos_cases <- subset(viral_cn_stats, viral_cn_stats$log10_prop_hpv16 > -4)

write.csv(hpv16_pos_cases, file = "104_hpv16_pos_cases.csv")





