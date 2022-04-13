
load("far_clade_parsimony.RData")

HLA_data <- read.csv("out_hla_data.csv")
data <- read.csv("104_hpv16_pos_cases.csv", header = TRUE)
merged <- data
merged$far_clade <- merged$study_id %in% far_clade_parsimony
merged <- merge(merged, HLA_data, by = "study_id")
merged <- subset(merged, merged$sample_TN %in% c("T"))


binned_hla_data <- data.frame(study_id = rep(merged$study_id, times = 2), hla_A = c(merged$A1, merged$A2), far_clade = rep(merged$far_clade, times = 2))
binned_hla_data$hla_B <- c(merged$B1, merged$B2)
binned_hla_data$hla_C <- c(merged$C1, merged$C2)

my_bar_pal <- c(wes_palette(n=5, name="IsleofDogs2")[1], wes_palette(n=5, name="IsleofDogs2")[3])
test_int <-ggplot(merge_int_subclonal, aes(x=far_clade, fill=highest_conf_int))
mybar_int <- test_int+ geom_bar(position = "fill") + theme_classic() + scale_fill_manual(values = my_bar_pal)+
  labs(y="Integrated", x = "Far Clade")
print(chisq.test(table(merge_int_subclonal$far_clade, merge_int_subclonal$highest_conf_int)))

pdf(file = "phylo_integration_bar.pdf", height = 2, width = 3)
mybar_int
dev.off()

##
##
## Subclonal 
##

data <- read.csv("clinical_data_update_7_21_21_ed_last_contact.csv", header = TRUE)
merged <- data
merged <- subset(merged, merged$OP_subsite %in% c("BOT", "tonsil"))
merged$far_clade <- merged$study_id %in% far_clade_parsimony

load("vcf_data_with_effect_all.RData")
subclonal_vars <- vcf_data_with_effect_all
subclonal_vars <- subset(subclonal_vars, subclonal_vars$vaf < 85)
subclonal <- subset(subclonal_vars, subclonal_vars$vaf > 15)

n_sub <- table(subclonal$study_id)
n_sub <- data.frame(n_sub)
names(n_sub) <- c("study_id", "n_subclonal")

merged_subclonal <- merge(merged, n_sub, by = "study_id", all.x = TRUE)
merged_subclonal$n_subclonal[is.na(merged_subclonal$n_subclonal)] <- 0
merged_subclonal$log_subclonal <- log(merged_subclonal$n_subclonal + 1, base = 10)

load("out_sublineage_data.RData")

A1_data <- subset(out_sublineage_data, out_sublineage_data$sublineage_JC69 %in% c("A1"))

merged_subclonal$A1_status <- merged_subclonal$study_id %in% A1_data$study_id
merged_subclonal$has_subclone <- merged_subclonal$n_subclonal > 0

#merged_subclonal <- subset(merged_subclonal, merged_subclonal$smoking_exposure <= 50) 

plot(as.factor(merged_subclonal$far_clade), merged_subclonal$log_subclonal)
points(as.factor(merged_subclonal$far_clade), merged_subclonal$log_subclonal)
print(wilcox.test(n_subclonal ~ far_clade, data = merged_subclonal))


library(stringr)
library(ggsignif)
library(ggplot2)
library("wesanderson")
library("ggpubr")


library(survival);
library(stringr)
library(survminer)
library(wesanderson)
surv_pal <- c(wes_palette(n=4, name="Darjeeling2")[3],wes_palette(n=4, name="Darjeeling2")[2])
surv_pal_sublin <- c(wes_palette(n=4, name="Royal1")[1],wes_palette(n=5, name="IsleofDogs2")[3])

my_box_pal <- c(wes_palette("Moonrise3")[5], wes_palette("Moonrise3")[1])

subclonal_panel <-ggplot(merged_subclonal, aes(x=far_clade, y = n_subclonal)) + 
  geom_boxplot(aes(y=n_subclonal), outlier.alpha = 0.0) +
  geom_jitter(height = 0.2, width = 0.15, aes(color=far_clade), cex = 3, alpha = 0.5) +
  #geom_dotplot(binaxis = "y", stackdir = "center",  aes(y=log_subclonal, fill = far_clade), position=position_dodge(0.05), binwidth = 0.02255)+
  theme_classic() +
  scale_color_manual(values = surv_pal) + 
  #coord_cartesian(ylim = c(  min(merged$plt_smo[!is.na(merged$plt_smo)])*1 , max(merged$plt_smo[!is.na(merged$plt_smo)])*1.1)) +
  #geom_signif(y_position = c(max(merged$plt_smo[!is.na(merged$plt_smo)])*1.05), xmin = c(1), xmax = c(2), annotation = c("NS"), tip_length = 0.03)+
  labs(y="N. Subclonal Poly.", x = "Clade")
print(wilcox.test(n_subclonal ~ far_clade, merged_subclonal))

my_bar_pal <- c(wes_palette(n=5, name="IsleofDogs2")[1], wes_palette(n=5, name="IsleofDogs2")[3])
test <-ggplot(merged_subclonal, aes(x=far_clade, fill=has_subclone))
mybar <- test+ geom_bar(position = "fill") + theme_classic() + scale_fill_manual(values = my_bar_pal)+
  labs(y="Has Subclonal", x = "Far Clade")
print(chisq.test(table(merged_subclonal$far_clade, merged_subclonal$has_subclone)))

pdf(file = "phylo_subclone_bar.pdf", height = 2, width = 3)
mybar
dev.off()
