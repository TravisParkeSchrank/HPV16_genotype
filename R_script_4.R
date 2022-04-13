
library(Gviz)
library(GenomicRanges)
library(Biostrings)
library("reshape2")
library("phylotools")
library("phytools")
library("phangorn")
library("geiger")
library("ggplot2")
library("ggtree")
library("treeio")
library("dplyr")
library("wesanderson")
library(RColorBrewer)
library(phangorn)
library(ape)
library(here)
library(treeman)
library(caper)

my_data <- read.phyDat(here("sample_HPV16_genotypes.fasta"), format = "fasta")

dm  <- dist.ml(my_data)
treeUPGMA  <- upgma(dm)
tree  <- NJ(dm)
plot(tree, main="UPGMA")
treePars  <- optim.parsimony(tree, my_data)

include_tumors <- row.names(as.matrix(dm))


library(Gviz)
library(GenomicRanges)
library(Biostrings)
library("reshape2")
library("phylotools")
library("phytools")
library("phangorn")
library("geiger")
library("ggplot2")
library("ggtree")
library("treeio")
library("dplyr")
library("wesanderson")
library("RColorBrewer")
library("ggnewscale")

HPV_copy_data <- read.csv("idx_stats_summary.csv")

HPV_copy_data$label <- HPV_copy_data$sample
HPV_copy_data$logRatioCN <- log(HPV_copy_data$hpv16_aligns/HPV_copy_data$human_aligns, base = 2)
HPV_copy_data$logRatioCN <- HPV_copy_data$logRatioCN + 12
row.names(HPV_copy_data) <- HPV_copy_data$sample

HPV_copy_data <- HPV_copy_data[include_tumors, ]




p <- ggtree(treePars)
p$data$myvar <- 1:212


holder <- p$data
ordered_samples <- holder$label[1:107]
HPV_copy_data <- HPV_copy_data[match(ordered_samples, HPV_copy_data$label), ]
p$data$viral_CN <- c(HPV_copy_data$logRatioCN, rep(2, times = 212-107) )

##
##
## SUBLINEAGE DATA
##
##
load("out_sublineage_data.RData")

sublin_dictionary <- out_sublineage_data$sublineage_JC69
names(sublin_dictionary) <- out_sublineage_data$study_id
assign_sublin <- sublin_dictionary[p$data$label[1:107]]
p$data$sublineage <- c(as.character(assign_sublin), rep("", times = 212-107))

short_dar <- wes_palette("Darjeeling2")[1:4]
sublineage_names <- c("A1","A2","A3","A4","B1","B2","B3","B4","C1","C2","C3","C4","D1","D2","D3","D4")
#pal <- wes_palette("Royal1", 16, type = "continuous")
pal <- c("grey",
         brewer.pal(n=9, "Blues")[6:8],
         rev(brewer.pal(n=9, "Purples")[5:8]),
         rev(brewer.pal(n=9, "Greens")[6:9]),
         brewer.pal(n=5, "Reds")[2:5]
)

pal <- alpha(pal, alpha = 0.85)
names(pal) <- sublineage_names
row.names(out_sublineage_data) <- out_sublineage_data$study_id
plot_sublineage <- data.frame(out_sublineage_data)[, names(out_sublineage_data) %in% c("sublineage_JC69")]
plot_sublineage <- data.frame(plot_sublineage)
row.names(plot_sublineage) <- out_sublineage_data$study_id


p2 <- gheatmap(p, plot_sublineage , font.size=3, colnames_angle= 90, offset = 0, width  = 0.1, colnames_offset_y = 1, colnames_position = "top") + 
  scale_fill_manual(breaks = names(pal), values=pal) + coord_cartesian(clip = "off")



#### N CODIN MUTATNTS
load("n_clonal_missense.RData")
n_muts_lib <- n_clonal_missense_vars$n_clonal_missense
names(n_muts_lib) <- n_clonal_missense_vars$study_id

n_muts <- rep(0, times = 212)
count = 0
for(i in p$data$label){
  count <- count +1
  if(i %in% names(n_muts_lib)){
    n_muts[count] <- n_muts_lib[[i]]
  }
}

p$data$n_muts <- n_muts
p$data$n_muts[is.na(p$data$n_muts)] <- 0

n_muts_pal <- rev(wes_palette(100, name = "Zissou1", type = "continuous")[1:70])

plot_n_muts <- n_clonal_missense_vars
plot_n_muts <- plot_n_muts[, names(plot_n_muts) %in% c("n_clonal_missense")]
plot_n_muts <- data.frame(plot_n_muts)
row.names(plot_n_muts) <- n_clonal_missense_vars$study_id
#plot_n_muts$plot_n_muts <- as.factor(plot_n_muts$plot_n_muts)

p3 <- p2 + new_scale_fill()
p4 <- gheatmap(p3, plot_n_muts, font.size=3, colnames_angle=90, offset = 2, width  = 0.1, colnames_offset_y = 1, colnames_position = "top")  + 
  scale_fill_steps(low = n_muts_pal[1] , high = brewer.pal(n=9, "Blues")[9], space="Lab", guide = "coloursteps")+ coord_cartesian(clip = "off")



plot_HPV_copy <- HPV_copy_data$logRatioCN
plot_HPV_copy <- data.frame(plot_HPV_copy)
row.names(plot_HPV_copy) <- HPV_copy_data$sample

p5 <- p4 + new_scale_fill()
#p6 <- gheatmap(p5, plot_HPV_copy, font.size=3, colnames_angle=90, offset = 4, width  = 0.1, colnames_offset_y = 1, colnames_position = "top")  + 
#  scale_fill_continuous()+ coord_cartesian(clip = "off")

p6 <- gheatmap(p5, plot_HPV_copy, font.size=3, colnames_angle=90, offset = 4, width  = 0.1, colnames_offset_y = 1, colnames_position = "top")  + 
    scale_fill_steps(low = "white", high = wes_palette("Cavalcanti1")[5], space="Lab", guide = "coloursteps")+ coord_cartesian(clip = "off")

load("missense_vars.RData")
missense_vars$var_id <- paste(missense_vars$gene, missense_vars$aa_change, sep = "_")
freq_data <- data.frame(table(missense_vars$var_id))
freq_data <- freq_data[rev(order(freq_data$Freq)), ]
freq_data <- head(freq_data, n = 15)

zero_mat <- matrix(0, 15, 107)
row.names(zero_mat) <- as.character(freq_data$Var1)
colnames(zero_mat) <- HPV_copy_data$sample

missense_vars_data <- subset(missense_vars, missense_vars$study_id %in% HPV_copy_data$sample)
missense_vars_data <- subset(missense_vars_data, missense_vars_data$var_id %in% freq_data$Var1)

for(i in 1:length(missense_vars_data$var_id)){
  zero_mat[missense_vars_data$var_id[i], missense_vars_data$study_id[i]] <- 1
}

var_matrix <- zero_mat
var_matrix <- var_matrix == 1
plot_vars <- data.frame(var_matrix)
#plot_vars <- t(plot_vars)
plot_vars <- data.frame(plot_vars)
plot_vars <- data.frame(t(plot_vars))

p7 <- p6 + new_scale_fill()
p8 <- gheatmap(p7, plot_vars, font.size=3, colnames_angle=90, offset = 6.5, width  = 0.75, colnames_offset_y = 1, colnames_position = "top")  + 
  scale_fill_manual(breaks = c("TRUE", "FALSE"), values = c(wes_palette("Moonrise2")[4], wes_palette("Moonrise1")[3]))+ coord_cartesian(clip = "off")


pdf(file = "phylo_waterfall.pdf", height = 7, width = 7)
print(p8)
dev.off()


load("out_sublineage_data.RData")
out_sublineage_data_op <- out_sublineage_data
out_sublineage_data_op$origin <- rep("OP", times = length(out_sublineage_data_op$study_id))

load("out_sublineage_cervix_data.RData")
out_sublineage_data_uc <- out_sublineage_data
out_sublineage_data_uc$origin <- rep("UC", times = length(out_sublineage_data_uc$study_id))


sublin_data_plot <- rbind(out_sublineage_data_op, out_sublineage_data_uc)

library(ggplot2)

test <-ggplot(sublin_data_plot, aes(x=origin, fill=sublineage_JC69))
mybar <- test+ geom_bar(position = "fill") + theme_classic() + scale_fill_manual(breaks = names(pal), values=pal) +
  labs(y="Sublineage", x = "Anat. Origin")
pdf(file = "sublineage_vs_UC_OP.pdf", width = 3.5, height = 4)
mybar
dev.off()

