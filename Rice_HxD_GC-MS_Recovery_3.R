#load required packages and functions
#NOTE: Refer to "Rice_HxD_Metabolomics" repository for the functions
library(VennDiagram)
library(grid)
library(RColorBrewer)
library(pcaMethods)
library(gplots)
library(outliers)
library(stats)
library(nortest)
library(plyr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(pheatmap)
source("H:/LMFLawas/PhD/Experiments/Protocols/R script/RemoveFactors_function.R")
source("H:/LMFLawas/PhD/Experiments/Protocols/R script/func_normalize.R")
source('H:/LMFLawas/PhD/Experiments/Protocols/R script/func_find_one_outlier.R')
source('H:/LMFLawas/PhD/Experiments/Protocols/R script/func_replace_outlier.R')
source('H:/LMFLawas/PhD/Experiments/Protocols/R script/func_hist_outlier.R')
source("H:/LMFLawas/PhD/Experiments/Protocols/R script/func_log_transform.R")
source("H:/LMFLawas/PhD/Experiments/Protocols/R script/func_log2_median_transform.R")


#NOTE: early stress means mild stress; late stress means severe stress
#NOTE: for initial data pre-processing (e.g. log and median transformation) and related data, refer to the other GitHub repository


###################################################
#DEVELOPING SEED#
###################################################



##################################################
#PRINCIPAL COMPONENT ANALYSIS#
##################################################


#PCA - controls, late stress, recovery

pca_res_combined_controls_LS_R_developing_seed = pca(data_combined_log2_mean_developing_seed_2, method = "ppca", nPcs = 5, scale = "pareto", center = TRUE)   
pca_res_scores_combined_controls_LS_R_developing_seed = scores(pca_res_combined_controls_LS_R_developing_seed)    #scores
pca_res_loadings_combined_controls_LS_R_developing_seed = loadings(pca_res_combined_controls_LS_R_developing_seed)  #loadings

#scores + additional info as data frame to use in score plot - ggplot
scores_combined_controls_LS_R_developing_seed = as.data.frame(pca_res_scores_combined_controls_LS_R_developing_seed)
scores_combined_controls_LS_R_developing_seed$Timepoint = interaction(data_combined_log2_mean_developing_seed$Treatment, data_combined_log2_mean_developing_seed$Timepoint)
scores_combined_controls_LS_R_developing_seed$Cultivar = data_combined_log2_mean_developing_seed$Cultivar
scores_combined_controls_LS_R_developing_seed$Label = gsub("G", "", rownames(scores_combined_controls_LS_R_developing_seed))


#biplot
pdf("biplot_pca_combined_controls_LS_R_developing_seed.pdf")
biplot(pca_res_combined_controls_LS_R_developing_seed, cex = 0.8, main = "Combined_Controls, Late stress, Recovery \nDeveloping seed")
dev.off()

#screeplot
pdf("screeplot_pca_combined_cotrols_LS_R_developing_seed.pdf")
text(barplot(pca_res_combined_controls_LS_R_developing_seed@R2*100, names.arg = c("PC1", "PC2", "PC3", "PC4", "PC5"), main = "Combined_Controls, Late stress, Recovery \nDeveloping seed", ylab = "Variance (%)", ylim = c(0, 50)), 0, round(pca_res_combined_controls_LS_R_developing_seed@R2*100, 3), pos = 3)
box()
dev.off()


#score plot

#symbols - legends at top
#in assigning colors and shapes to timepoints and cultivars, check the order of factor levels first
#colors similar to FL plot
#with labels below symbols
png("score_plot_pca_combined_controls_LS_R_developing_seed.png", width = 8*300, height = 8*300, res = 300)
ggplot(scores_combined_controls_LS_R_developing_seed, 
       aes(x = PC1, 
           y = PC2, 
           color = Timepoint, 
           fill = Timepoint, 
           shape = Cultivar, 
           label = Label)) +
  geom_point(size = 7) +
  geom_text(vjust = 2, color = "black") +
  scale_color_manual(breaks = c("Control.Late stress", "Control.12h Rewatering", "Control.36h Rewatering", "Control.60h Rewatering",
                                "Heat & drought.Late stress", "Heat & drought.12h Rewatering", "Heat & drought.36h Rewatering", 
                                "Heat & drought.60h Rewatering"),
                     values = c("#fed976", "#bd0026", "#feb24c", "#e31a1c", "#fd8d3c", "#fc4e2a", "#ffeda0", "#800026"),
                     labels = c("Control - Late stress", "Control - 12h rewatering", "Control - 36h rewatering", 
                                "Control - 60h rewatering", "Late stress", "12h rewatering", "36h rewatering", "60h rewatering")) +
  scale_fill_manual(breaks = c("Control.Late stress", "Control.12h Rewatering", "Control.36h Rewatering", "Control.60h Rewatering",
                               "Heat & drought.Late stress", "Heat & drought.12h Rewatering", "Heat & drought.36h Rewatering", 
                               "Heat & drought.60h Rewatering"),
                    values = c("#fed976", "#bd0026", "#feb24c", "#e31a1c", "#fd8d3c", "#fc4e2a", "#ffeda0", "#800026"),
                    labels = c("Control - Late stress", "Control - 12h rewatering", "Control - 36h rewatering", 
                               "Control - 60h rewatering", "Late stress", "12h rewatering", "36h rewatering", "60h rewatering")) +
  scale_shape_manual(breaks = c("N22", "Dular", "Anjali"),
                     values = c(22, 24, 21),
                     labels = c("N22", "Dular", "Anjali")) +
  scale_x_continuous(limits = c(-8, 8), breaks = seq(-8, 8, by = 2)) +
  scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, by = 2)) +
  xlab(paste("PC1 (", round(pca_res_combined_controls_LS_R_developing_seed@R2[1], 4)*100, "%)", sep = "")) +
  ylab(paste("PC2 (", round(pca_res_combined_controls_LS_R_developing_seed@R2[2], 4)*100, "%)", sep = "")) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_blank(),
        legend.justification = c(1, 0),
        legend.box.just = "top",
        legend.position = c(0.995, 0.83),
        legend.background = element_blank(),
        legend.text = element_text(size = 10),
        legend.key = element_blank(),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.spacing = unit(16.5, "lines"),
        legend.key.height = unit(1.3, "lines")) +
  guides(shape = guide_legend(order = 1, override.aes = list(size = 5)),
         fill = guide_legend(override.aes = list(size = 5),
                             ncol = 2),
         color = guide_legend(ncol = 2))
dev.off()


#export as table then upload in flag leaf R script to make a multipanel figure of PCAs of all organs
write.table(scores_combined_controls_LS_R_developing_seed, file = "PCA_scores_combined_developing_seed.txt", sep = "\t", quote = F)


###################################################
#comparison of EGF controls
###################################################


#relative levels of metabolites under control conditions
control_order = c("GNCLS", "GNC12", "GNC36", "GNC60", "GDCLS", "GDC12", "GDC36", "GDC60", "GACLS", "GAC12", "GAC36", "GAC60")
data_combined_log2_mean_developing_seed_controls = data_combined_log2_mean_developing_seed_2[control_order, ]
data_combined_log2_mean_developing_seed_controls = t(data_combined_log2_mean_developing_seed_controls)


#heatmap including all metabolites
heatmap.2(data_combined_log2_mean_developing_seed_controls, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(4, 20), 
          main = "Developing seed - Combined \nRelative levels, Controls", 
          cexCol = 1, 
          cexRow = 0.9, 
          density.info = "none",
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          srtCol = 90)



#compare controls of rewatering time points with the earliest control timepoint (control corresponding to late stress)


#CLS vs. C12

#N22

#data distribution
shapiro_res_combined_N22_CLS_C12_developing_seed = t(data.frame(lapply(data_combined_N22_log2median_developing_seed[which(data_combined_N22_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Control.12h Rewatering")), 3:ncol(data_combined_N22_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_CLS_C12_developing_seed) = "p.value"
shapiro_res_combined_N22_CLS_C12_developing_seed = as.data.frame(shapiro_res_combined_N22_CLS_C12_developing_seed)
shapiro_res_combined_N22_CLS_C12_developing_seed$Distribution = ifelse(shapiro_res_combined_N22_CLS_C12_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_CLS_C12_developing_seed$Distribution == "Normal")     #30 metabolites
sum(shapiro_res_combined_N22_CLS_C12_developing_seed$Distribution == "Non-normal") #37 metabolites

#wilcox test
wilcox_res_N22_CLS_C12_developing_seed_combined = t(data.frame(lapply(data_combined_N22_log2median_developing_seed[ , 3:ncol(data_combined_N22_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_N22_log2median_developing_seed$Timepoint, subset = data_combined_N22_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Control.12h Rewatering"))$p.value)))
colnames(wilcox_res_N22_CLS_C12_developing_seed_combined) = "p.value"
rownames(wilcox_res_N22_CLS_C12_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_N22_CLS_C12_developing_seed_combined = as.data.frame(wilcox_res_N22_CLS_C12_developing_seed_combined)
wilcox_res_N22_CLS_C12_developing_seed_combined$Significance = with(wilcox_res_N22_CLS_C12_developing_seed_combined,
                                                                    ifelse(p.value <= 0.001, "***", ifelse(
                                                                      p.value <= 0.01, "**", ifelse(
                                                                        p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_CLS_C12_developing_seed_combined$Significance == "*") #8 metabolites
sum(wilcox_res_N22_CLS_C12_developing_seed_combined$Significance == "**") #2 metabolites
sum(wilcox_res_N22_CLS_C12_developing_seed_combined$Significance == "***") #0 metabolite
sum(wilcox_res_N22_CLS_C12_developing_seed_combined$Significance != "ns") #10 metabolites
sum(wilcox_res_N22_CLS_C12_developing_seed_combined$Significance == "ns") #57 metabolites


#Dular

#data distribution
shapiro_res_combined_Dular_CLS_C12_developing_seed = t(data.frame(lapply(data_combined_Dular_log2median_developing_seed[which(data_combined_Dular_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Control.12h Rewatering")), 3:ncol(data_combined_Dular_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_CLS_C12_developing_seed) = "p.value"
shapiro_res_combined_Dular_CLS_C12_developing_seed = as.data.frame(shapiro_res_combined_Dular_CLS_C12_developing_seed)
shapiro_res_combined_Dular_CLS_C12_developing_seed$Distribution = ifelse(shapiro_res_combined_Dular_CLS_C12_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_CLS_C12_developing_seed$Distribution == "Normal")     #50 metabolites
sum(shapiro_res_combined_Dular_CLS_C12_developing_seed$Distribution == "Non-normal") #17 metabolites

#wilcox test
wilcox_res_Dular_CLS_C12_developing_seed_combined = t(data.frame(lapply(data_combined_Dular_log2median_developing_seed[ , 3:ncol(data_combined_Dular_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_developing_seed$Timepoint, subset = data_combined_Dular_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Control.12h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_CLS_C12_developing_seed_combined) = "p.value"
rownames(wilcox_res_Dular_CLS_C12_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_Dular_CLS_C12_developing_seed_combined = as.data.frame(wilcox_res_Dular_CLS_C12_developing_seed_combined)
wilcox_res_Dular_CLS_C12_developing_seed_combined$Significance = with(wilcox_res_Dular_CLS_C12_developing_seed_combined,
                                                                      ifelse(p.value <= 0.001, "***", ifelse(
                                                                        p.value <= 0.01, "**", ifelse(
                                                                          p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_CLS_C12_developing_seed_combined$Significance == "*") #14 metabolites
sum(wilcox_res_Dular_CLS_C12_developing_seed_combined$Significance == "**") #8 metabolites
sum(wilcox_res_Dular_CLS_C12_developing_seed_combined$Significance == "***") #2 metabolites
sum(wilcox_res_Dular_CLS_C12_developing_seed_combined$Significance != "ns") #24 metabolites
sum(wilcox_res_Dular_CLS_C12_developing_seed_combined$Significance == "ns") #43 metabolites


#Anjali

#data distribution
shapiro_res_combined_Anjali_CLS_C12_developing_seed = t(data.frame(lapply(data_combined_Anjali_log2median_developing_seed[which(data_combined_Anjali_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Control.12h Rewatering")), 3:ncol(data_combined_Anjali_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_CLS_C12_developing_seed) = "p.value"
shapiro_res_combined_Anjali_CLS_C12_developing_seed = as.data.frame(shapiro_res_combined_Anjali_CLS_C12_developing_seed)
shapiro_res_combined_Anjali_CLS_C12_developing_seed$Distribution = ifelse(shapiro_res_combined_Anjali_CLS_C12_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_CLS_C12_developing_seed$Distribution == "Normal")     #44 metabolites
sum(shapiro_res_combined_Anjali_CLS_C12_developing_seed$Distribution == "Non-normal") #23 metabolites

#wilcox test
wilcox_res_Anjali_CLS_C12_developing_seed_combined = t(data.frame(lapply(data_combined_Anjali_log2median_developing_seed[ , 3:ncol(data_combined_Anjali_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_developing_seed$Timepoint, subset = data_combined_Anjali_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Control.12h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_CLS_C12_developing_seed_combined) = "p.value"
rownames(wilcox_res_Anjali_CLS_C12_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_Anjali_CLS_C12_developing_seed_combined = as.data.frame(wilcox_res_Anjali_CLS_C12_developing_seed_combined)
wilcox_res_Anjali_CLS_C12_developing_seed_combined$Significance = with(wilcox_res_Anjali_CLS_C12_developing_seed_combined,
                                                                       ifelse(p.value <= 0.001, "***", ifelse(
                                                                         p.value <= 0.01, "**", ifelse(
                                                                           p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_CLS_C12_developing_seed_combined$Significance == "*") #5 metabolites
sum(wilcox_res_Anjali_CLS_C12_developing_seed_combined$Significance == "**") #1 metabolite
sum(wilcox_res_Anjali_CLS_C12_developing_seed_combined$Significance == "***") #1 metabolite
sum(wilcox_res_Anjali_CLS_C12_developing_seed_combined$Significance != "ns") #7 metabolites
sum(wilcox_res_Anjali_CLS_C12_developing_seed_combined$Significance == "ns") #60 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_CLS_C12_developing_seed_combined = cbind(wilcox_res_N22_CLS_C12_developing_seed_combined, 
                                                    wilcox_res_Dular_CLS_C12_developing_seed_combined$Significance, 
                                                    wilcox_res_Anjali_CLS_C12_developing_seed_combined$Significance)
wilcox_sig_CLS_C12_developing_seed_combined = wilcox_sig_CLS_C12_developing_seed_combined[ , -1]
colnames(wilcox_sig_CLS_C12_developing_seed_combined)[1] = "wilcox_res_N22_CLS_C12_developing_seed_combined$Significance"
wilcox_sig_CLS_C12_developing_seed_combined = data.frame(lapply(wilcox_sig_CLS_C12_developing_seed_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_CLS_C12_developing_seed_combined) = rownames(wilcox_res_N22_CLS_C12_developing_seed_combined)
wilcox_sig_CLS_C12_developing_seed_combined$Nonsig = rowSums(wilcox_sig_CLS_C12_developing_seed_combined == " ")  #count number of nonsignificant
wilcox_sig_final_CLS_C12_developing_seed_combined = wilcox_sig_CLS_C12_developing_seed_combined[-which(wilcox_sig_CLS_C12_developing_seed_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 34 metabolites

wilcox_sig_log2mean_CLS_C12_developing_seed_combined = data_combined_log2_mean_developing_seed_controls[which(rownames(data_combined_log2_mean_developing_seed_controls) %in% rownames(wilcox_sig_final_CLS_C12_developing_seed_combined)), ]  #relative control levels of metabolites significant in at least one of the comparisons

#overview
heatmap.2(wilcox_sig_log2mean_CLS_C12_developing_seed_combined, 
          Rowv = TRUE,
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(5, 20), 
          main = "Developing seed-Combined-Wilcoxon \nRelative levels, Sig in CLS vs C12",
          cexCol = 1.5,
          cexRow = 0.9, 
          density.info = "none",
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          srtCol = 90, 
          adjCol = c(0.9, 0.5), 
          sepcolor = "black",
          colsep = c(seq(0, 12)),
          rowsep = c(seq(0, 34)),
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))


####################

#CLS vs. C36

#N22

#data distribution
shapiro_res_combined_N22_CLS_C36_developing_seed = t(data.frame(lapply(data_combined_N22_log2median_developing_seed[which(data_combined_N22_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Control.36h Rewatering")), 3:ncol(data_combined_N22_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_CLS_C36_developing_seed) = "p.value"
shapiro_res_combined_N22_CLS_C36_developing_seed = as.data.frame(shapiro_res_combined_N22_CLS_C36_developing_seed)
shapiro_res_combined_N22_CLS_C36_developing_seed$Distribution = ifelse(shapiro_res_combined_N22_CLS_C36_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_CLS_C36_developing_seed$Distribution == "Normal")     #43 metabolites
sum(shapiro_res_combined_N22_CLS_C36_developing_seed$Distribution == "Non-normal") #24 metabolites

#wilcox test
wilcox_res_N22_CLS_C36_developing_seed_combined = t(data.frame(lapply(data_combined_N22_log2median_developing_seed[ , 3:ncol(data_combined_N22_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_N22_log2median_developing_seed$Timepoint, subset = data_combined_N22_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Control.36h Rewatering"))$p.value)))
colnames(wilcox_res_N22_CLS_C36_developing_seed_combined) = "p.value"
rownames(wilcox_res_N22_CLS_C36_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_N22_CLS_C36_developing_seed_combined = as.data.frame(wilcox_res_N22_CLS_C36_developing_seed_combined)
wilcox_res_N22_CLS_C36_developing_seed_combined$Significance = with(wilcox_res_N22_CLS_C36_developing_seed_combined,
                                                                    ifelse(p.value <= 0.001, "***", ifelse(
                                                                      p.value <= 0.01, "**", ifelse(
                                                                        p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_CLS_C36_developing_seed_combined$Significance == "*") #4 metabolites
sum(wilcox_res_N22_CLS_C36_developing_seed_combined$Significance == "**") #17 metabolites
sum(wilcox_res_N22_CLS_C36_developing_seed_combined$Significance == "***") #6 metabolites
sum(wilcox_res_N22_CLS_C36_developing_seed_combined$Significance != "ns") #27 metabolites
sum(wilcox_res_N22_CLS_C36_developing_seed_combined$Significance == "ns") #40 metabolites


#Dular

#data distribution
shapiro_res_combined_Dular_CLS_C36_developing_seed = t(data.frame(lapply(data_combined_Dular_log2median_developing_seed[which(data_combined_Dular_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Control.36h Rewatering")), 3:ncol(data_combined_Dular_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_CLS_C36_developing_seed) = "p.value"
shapiro_res_combined_Dular_CLS_C36_developing_seed = as.data.frame(shapiro_res_combined_Dular_CLS_C36_developing_seed)
shapiro_res_combined_Dular_CLS_C36_developing_seed$Distribution = ifelse(shapiro_res_combined_Dular_CLS_C36_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_CLS_C36_developing_seed$Distribution == "Normal")     #55 metabolites
sum(shapiro_res_combined_Dular_CLS_C36_developing_seed$Distribution == "Non-normal") #12 metabolites

#wilcox test
wilcox_res_Dular_CLS_C36_developing_seed_combined = t(data.frame(lapply(data_combined_Dular_log2median_developing_seed[ , 3:ncol(data_combined_Dular_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_developing_seed$Timepoint, subset = data_combined_Dular_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Control.36h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_CLS_C36_developing_seed_combined) = "p.value"
rownames(wilcox_res_Dular_CLS_C36_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_Dular_CLS_C36_developing_seed_combined = as.data.frame(wilcox_res_Dular_CLS_C36_developing_seed_combined)
wilcox_res_Dular_CLS_C36_developing_seed_combined$Significance = with(wilcox_res_Dular_CLS_C36_developing_seed_combined,
                                                                      ifelse(p.value <= 0.001, "***", ifelse(
                                                                        p.value <= 0.01, "**", ifelse(
                                                                          p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_CLS_C36_developing_seed_combined$Significance == "*") #13 metabolites
sum(wilcox_res_Dular_CLS_C36_developing_seed_combined$Significance == "**") #6 metabolites
sum(wilcox_res_Dular_CLS_C36_developing_seed_combined$Significance == "***") #8 metabolites
sum(wilcox_res_Dular_CLS_C36_developing_seed_combined$Significance != "ns") #27 metabolites
sum(wilcox_res_Dular_CLS_C36_developing_seed_combined$Significance == "ns") #40 metabolites


#Anjali

#data distribution
shapiro_res_combined_Anjali_CLS_C36_developing_seed = t(data.frame(lapply(data_combined_Anjali_log2median_developing_seed[which(data_combined_Anjali_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Control.36h Rewatering")), 3:ncol(data_combined_Anjali_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_CLS_C36_developing_seed) = "p.value"
shapiro_res_combined_Anjali_CLS_C36_developing_seed = as.data.frame(shapiro_res_combined_Anjali_CLS_C36_developing_seed)
shapiro_res_combined_Anjali_CLS_C36_developing_seed$Distribution = ifelse(shapiro_res_combined_Anjali_CLS_C36_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_CLS_C36_developing_seed$Distribution == "Normal")     #47 metabolites
sum(shapiro_res_combined_Anjali_CLS_C36_developing_seed$Distribution == "Non-normal") #20 metabolites

#wilcox test
wilcox_res_Anjali_CLS_C36_developing_seed_combined = t(data.frame(lapply(data_combined_Anjali_log2median_developing_seed[ , 3:ncol(data_combined_Anjali_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_developing_seed$Timepoint, subset = data_combined_Anjali_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Control.36h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_CLS_C36_developing_seed_combined) = "p.value"
rownames(wilcox_res_Anjali_CLS_C36_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_Anjali_CLS_C36_developing_seed_combined = as.data.frame(wilcox_res_Anjali_CLS_C36_developing_seed_combined)
wilcox_res_Anjali_CLS_C36_developing_seed_combined$Significance = with(wilcox_res_Anjali_CLS_C36_developing_seed_combined,
                                                                       ifelse(p.value <= 0.001, "***", ifelse(
                                                                         p.value <= 0.01, "**", ifelse(
                                                                           p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_CLS_C36_developing_seed_combined$Significance == "*") #9 metabolites
sum(wilcox_res_Anjali_CLS_C36_developing_seed_combined$Significance == "**") #2 metabolites
sum(wilcox_res_Anjali_CLS_C36_developing_seed_combined$Significance == "***") #3 metabolites
sum(wilcox_res_Anjali_CLS_C36_developing_seed_combined$Significance != "ns") #14 metabolites
sum(wilcox_res_Anjali_CLS_C36_developing_seed_combined$Significance == "ns") #53 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_CLS_C36_developing_seed_combined = cbind(wilcox_res_N22_CLS_C36_developing_seed_combined, 
                                                    wilcox_res_Dular_CLS_C36_developing_seed_combined$Significance,
                                                    wilcox_res_Anjali_CLS_C36_developing_seed_combined$Significance)
wilcox_sig_CLS_C36_developing_seed_combined = wilcox_sig_CLS_C36_developing_seed_combined[ , -1]
colnames(wilcox_sig_CLS_C36_developing_seed_combined)[1] = "wilcox_res_N22_CLS_C36_developing_seed_combined$Significance"
wilcox_sig_CLS_C36_developing_seed_combined = data.frame(lapply(wilcox_sig_CLS_C36_developing_seed_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_CLS_C36_developing_seed_combined) = rownames(wilcox_res_N22_CLS_C36_developing_seed_combined)
wilcox_sig_CLS_C36_developing_seed_combined$Nonsig = rowSums(wilcox_sig_CLS_C36_developing_seed_combined == " ")  #count number of nonsignificant
wilcox_sig_final_CLS_C36_developing_seed_combined = wilcox_sig_CLS_C36_developing_seed_combined[-which(wilcox_sig_CLS_C36_developing_seed_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 41 metabolites

wilcox_sig_log2mean_CLS_C36_developing_seed_combined = data_combined_log2_mean_developing_seed_controls[which(rownames(data_combined_log2_mean_developing_seed_controls) %in% rownames(wilcox_sig_final_CLS_C36_developing_seed_combined)), ]  #relative control levels of metabolites significant in at least one of the comparisons

#overview
heatmap.2(wilcox_sig_log2mean_CLS_C36_developing_seed_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599),
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(5, 20), 
          main = "Developing seed-Combined-Wilcoxon \nRelative levels, Sig in CLS vs C36", 
          cexCol = 1.5, 
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA,
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          srtCol = 90, 
          adjCol = c(0.9, 0.5), 
          sepcolor = "black",
          colsep = c(seq(0, 36)),
          rowsep = c(seq(0, 41)),
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))


####################

#CLS vs. C60

#N22

#data distribution
shapiro_res_combined_N22_CLS_C60_developing_seed = t(data.frame(lapply(data_combined_N22_log2median_developing_seed[which(data_combined_N22_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Control.60h Rewatering")), 3:ncol(data_combined_N22_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_CLS_C60_developing_seed) = "p.value"
shapiro_res_combined_N22_CLS_C60_developing_seed = as.data.frame(shapiro_res_combined_N22_CLS_C60_developing_seed)
shapiro_res_combined_N22_CLS_C60_developing_seed$Distribution = ifelse(shapiro_res_combined_N22_CLS_C60_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_CLS_C60_developing_seed$Distribution == "Normal")     #45 metabolites
sum(shapiro_res_combined_N22_CLS_C60_developing_seed$Distribution == "Non-normal") #22 metabolites

#wilcox test
wilcox_res_N22_CLS_C60_developing_seed_combined = t(data.frame(lapply(data_combined_N22_log2median_developing_seed[ , 3:ncol(data_combined_N22_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_N22_log2median_developing_seed$Timepoint, subset = data_combined_N22_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Control.60h Rewatering"))$p.value)))
colnames(wilcox_res_N22_CLS_C60_developing_seed_combined) = "p.value"
rownames(wilcox_res_N22_CLS_C60_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_N22_CLS_C60_developing_seed_combined = as.data.frame(wilcox_res_N22_CLS_C60_developing_seed_combined)
wilcox_res_N22_CLS_C60_developing_seed_combined$Significance = with(wilcox_res_N22_CLS_C60_developing_seed_combined,
                                                                    ifelse(p.value <= 0.001, "***", ifelse(
                                                                      p.value <= 0.01, "**", ifelse(
                                                                        p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_CLS_C60_developing_seed_combined$Significance == "*") #5 metabolites
sum(wilcox_res_N22_CLS_C60_developing_seed_combined$Significance == "**") #14 metabolites
sum(wilcox_res_N22_CLS_C60_developing_seed_combined$Significance == "***") #20 metabolites
sum(wilcox_res_N22_CLS_C60_developing_seed_combined$Significance != "ns") #39 metabolites
sum(wilcox_res_N22_CLS_C60_developing_seed_combined$Significance == "ns") #28 metabolites


#Dular

#data distribution
shapiro_res_combined_Dular_CLS_C60_developing_seed = t(data.frame(lapply(data_combined_Dular_log2median_developing_seed[which(data_combined_Dular_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Control.60h Rewatering")), 3:ncol(data_combined_Dular_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_CLS_C60_developing_seed) = "p.value"
shapiro_res_combined_Dular_CLS_C60_developing_seed = as.data.frame(shapiro_res_combined_Dular_CLS_C60_developing_seed)
shapiro_res_combined_Dular_CLS_C60_developing_seed$Distribution = ifelse(shapiro_res_combined_Dular_CLS_C60_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_CLS_C60_developing_seed$Distribution == "Normal")     #56 metabolites
sum(shapiro_res_combined_Dular_CLS_C60_developing_seed$Distribution == "Non-normal") #11 metabolites

#wilcox test
wilcox_res_Dular_CLS_C60_developing_seed_combined = t(data.frame(lapply(data_combined_Dular_log2median_developing_seed[ , 3:ncol(data_combined_Dular_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_developing_seed$Timepoint, subset = data_combined_Dular_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Control.60h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_CLS_C60_developing_seed_combined) = "p.value"
rownames(wilcox_res_Dular_CLS_C60_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_Dular_CLS_C60_developing_seed_combined = as.data.frame(wilcox_res_Dular_CLS_C60_developing_seed_combined)
wilcox_res_Dular_CLS_C60_developing_seed_combined$Significance = with(wilcox_res_Dular_CLS_C60_developing_seed_combined,
                                                                      ifelse(p.value <= 0.001, "***", ifelse(
                                                                        p.value <= 0.01, "**", ifelse(
                                                                          p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_CLS_C60_developing_seed_combined$Significance == "*") #8 metabolites
sum(wilcox_res_Dular_CLS_C60_developing_seed_combined$Significance == "**") #11 metabolites
sum(wilcox_res_Dular_CLS_C60_developing_seed_combined$Significance == "***") #19 metabolites
sum(wilcox_res_Dular_CLS_C60_developing_seed_combined$Significance != "ns") #38 metabolites
sum(wilcox_res_Dular_CLS_C60_developing_seed_combined$Significance == "ns") #29 metabolites


#Anjali

#data distribution
shapiro_res_combined_Anjali_CLS_C60_developing_seed = t(data.frame(lapply(data_combined_Anjali_log2median_developing_seed[which(data_combined_Anjali_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Control.60h Rewatering")), 3:ncol(data_combined_Anjali_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_CLS_C60_developing_seed) = "p.value"
shapiro_res_combined_Anjali_CLS_C60_developing_seed = as.data.frame(shapiro_res_combined_Anjali_CLS_C60_developing_seed)
shapiro_res_combined_Anjali_CLS_C60_developing_seed$Distribution = ifelse(shapiro_res_combined_Anjali_CLS_C60_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_CLS_C60_developing_seed$Distribution == "Normal")     #49 metabolites
sum(shapiro_res_combined_Anjali_CLS_C60_developing_seed$Distribution == "Non-normal") #18 metabolites

#wilcox test
wilcox_res_Anjali_CLS_C60_developing_seed_combined = t(data.frame(lapply(data_combined_Anjali_log2median_developing_seed[ , 3:ncol(data_combined_Anjali_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_developing_seed$Timepoint, subset = data_combined_Anjali_log2median_developing_seed$Timepoint %in% c("Control.Late stress", "Control.60h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_CLS_C60_developing_seed_combined) = "p.value"
rownames(wilcox_res_Anjali_CLS_C60_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_Anjali_CLS_C60_developing_seed_combined = as.data.frame(wilcox_res_Anjali_CLS_C60_developing_seed_combined)
wilcox_res_Anjali_CLS_C60_developing_seed_combined$Significance = with(wilcox_res_Anjali_CLS_C60_developing_seed_combined,
                                                                       ifelse(p.value <= 0.001, "***", ifelse(
                                                                         p.value <= 0.01, "**", ifelse(
                                                                           p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_CLS_C60_developing_seed_combined$Significance == "*") #7 metabolites
sum(wilcox_res_Anjali_CLS_C60_developing_seed_combined$Significance == "**") #10 metabolites
sum(wilcox_res_Anjali_CLS_C60_developing_seed_combined$Significance == "***") #11 metabolites
sum(wilcox_res_Anjali_CLS_C60_developing_seed_combined$Significance != "ns") #28 metabolites
sum(wilcox_res_Anjali_CLS_C60_developing_seed_combined$Significance == "ns") #39 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_CLS_C60_developing_seed_combined = cbind(wilcox_res_N22_CLS_C60_developing_seed_combined, 
                                                    wilcox_res_Dular_CLS_C60_developing_seed_combined$Significance, 
                                                    wilcox_res_Anjali_CLS_C60_developing_seed_combined$Significance)
wilcox_sig_CLS_C60_developing_seed_combined = wilcox_sig_CLS_C60_developing_seed_combined[ , -1]
colnames(wilcox_sig_CLS_C60_developing_seed_combined)[1] = "wilcox_res_N22_CLS_C60_developing_seed_combined$Significance"
wilcox_sig_CLS_C60_developing_seed_combined = data.frame(lapply(wilcox_sig_CLS_C60_developing_seed_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_CLS_C60_developing_seed_combined) = rownames(wilcox_res_N22_CLS_C60_developing_seed_combined)
wilcox_sig_CLS_C60_developing_seed_combined$Nonsig = rowSums(wilcox_sig_CLS_C60_developing_seed_combined == " ")  #count number of nonsignificant
wilcox_sig_final_CLS_C60_developing_seed_combined = wilcox_sig_CLS_C60_developing_seed_combined[-which(wilcox_sig_CLS_C60_developing_seed_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 53 metabolites

wilcox_sig_log2mean_CLS_C60_developing_seed_combined = data_combined_log2_mean_developing_seed_controls[which(rownames(data_combined_log2_mean_developing_seed_controls) %in% rownames(wilcox_sig_final_CLS_C60_developing_seed_combined)), ]  #relative control levels of metabolites significant in at least one of the comparisons

#overview
heatmap.2(wilcox_sig_log2mean_CLS_C60_developing_seed_combined,
          Rowv = TRUE, 
          Colv = FALSE,
          dendrogram = "row",
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(5, 20), 
          main = "Developing seed-Combined-Wilcoxon \nRelative levels, Sig in CLS vs C60", 
          cexCol = 1.5, 
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8),
          lwid = c(1.5, 5), 
          srtCol = 90, 
          adjCol = c(0.9, 0.5), 
          sepcolor = "black", 
          colsep = c(seq(0, 60)),
          rowsep = c(seq(0, 53)),
          sepwidth = c(0.01, 0.01),
          hclustfun = function(x) hclust(x, method = "average"))


####################

#significance of wilcoxon test - all cultivars, all four control timepoints

wilcox_sig_CLS_C12_C36_C60_developing_seed_combined = cbind(wilcox_res_N22_CLS_C12_developing_seed_combined,
                                                            wilcox_res_N22_CLS_C36_developing_seed_combined$Significance,
                                                            wilcox_res_N22_CLS_C60_developing_seed_combined$Significance,
                                                            wilcox_res_Dular_CLS_C12_developing_seed_combined$Significance,
                                                            wilcox_res_Dular_CLS_C36_developing_seed_combined$Significance,
                                                            wilcox_res_Dular_CLS_C60_developing_seed_combined$Significance,
                                                            wilcox_res_Anjali_CLS_C12_developing_seed_combined$Significance,
                                                            wilcox_res_Anjali_CLS_C36_developing_seed_combined$Significance,
                                                            wilcox_res_Anjali_CLS_C60_developing_seed_combined$Significance)
wilcox_sig_CLS_C12_C36_C60_developing_seed_combined = wilcox_sig_CLS_C12_C36_C60_developing_seed_combined[ , -1]
colnames(wilcox_sig_CLS_C12_C36_C60_developing_seed_combined)[1] = "wilcox_res_N22_CLS_C12_developing_seed_combined$Significance"
wilcox_sig_CLS_C12_C36_C60_developing_seed_combined = data.frame(lapply(wilcox_sig_CLS_C12_C36_C60_developing_seed_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_CLS_C12_C36_C60_developing_seed_combined) = rownames(wilcox_res_N22_CLS_C12_developing_seed_combined)
wilcox_sig_CLS_C12_C36_C60_developing_seed_combined$Nonsig = rowSums(wilcox_sig_CLS_C12_C36_C60_developing_seed_combined == " ")  #count number of nonsignificant
wilcox_sig_final_CLS_C12_C36_C60_developing_seed_combined = wilcox_sig_CLS_C12_C36_C60_developing_seed_combined[-which(wilcox_sig_CLS_C12_C36_C60_developing_seed_combined$Nonsig == "9"), ] #metabolites significant in at least one of the comparisons - 57 metabolites

wilcox_sig_log2mean_CLS_C12_C36_C60_developing_seed_combined = data_combined_log2_mean_developing_seed_controls[which(rownames(data_combined_log2_mean_developing_seed_controls) %in% rownames(wilcox_sig_final_CLS_C12_C36_C60_developing_seed_combined)), ]  #relative control levels of metabolites significant in at least one of the comparisons


#significance as cell notes
wilcox_sig_final_CLS_C12_C36_C60_developing_seed_combined_2 = wilcox_sig_final_CLS_C12_C36_C60_developing_seed_combined[ , -10]
wilcox_sig_final_CLS_C12_C36_C60_developing_seed_combined_2$N22_Late_stress = ""
wilcox_sig_final_CLS_C12_C36_C60_developing_seed_combined_2$Dular_Late_stress = ""
wilcox_sig_final_CLS_C12_C36_C60_developing_seed_combined_2$Anjali_Late_stress = ""
wilcox_sig_final_CLS_C12_C36_C60_developing_seed_combined_2 = wilcox_sig_final_CLS_C12_C36_C60_developing_seed_combined_2[ , c(10, 1, 2, 3, 11, 4, 5, 6, 12, 7, 8, 9)]


#visualization
png("wilcox_sig_log2mean_controls_developing_seed_combined.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_log2mean_CLS_C12_C36_C60_developing_seed_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(7.5, 18),
          lmat = rbind(c(0, 4, 5), c(0, 1, 0), c(3, 2, 0)),
          labCol = rep(c("Control - Severe stress", "Control - 12 h RW", "Control - 36 h RW", "Control - 60 h RW"), 3),
          cexCol = 1, 
          cexRow = 1.1, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 0.2, 11), 
          lwid = c(0.8, 5, 1.5),  
          srtCol = 45, 
          cellnote = wilcox_sig_final_CLS_C12_C36_C60_developing_seed_combined_2, 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 12)), 
          rowsep = c(seq(0, 57)), 
          sepwidth = c(0.01, 0.01), 
          ColSideColors = c(rep("#66c2a5", 4), rep("#fc8d62", 4), rep("#8da0cb", 4)),
          hclustfun = function(x) hclust(x, method = "average"))
legend(x = 0.0065, y = 1.1, xpd = TRUE, legend = "N22", bty = "n", cex = 1.2)
legend(x = 0.153, y = 1.1, xpd = TRUE, legend = "Dular", bty = "n", cex = 1.2)
legend(x = 0.305, y = 1.1, xpd = TRUE, legend = "Anjali", bty = "n", cex = 1.2)
legend(x = 0.87, y = 1.04, xpd = TRUE, legend = "Relative content", bty = "n", cex = 0.8)
dev.off()


####################


#analysis of rewatering time points

#late stress vs. recovery - NOTE: consider only the stress-responsive metabolites (i.e. metabolites which have significant difference in LS/CLS comparison)

#late stress vs. recovery

#log2-fold change - late stress vs. recovery
latestress_recovery_order = c("GNLS", "GNR12", "GNR36", "GNR60", "GDLS", "GDR12", "GDR36", "GDR60", "GALS", "GAR12", "GAR36", "GAR60")
data_combined_latestress_recovery_log2_mean_developing_seed = data_combined_log2_mean_developing_seed_2[latestress_recovery_order, ]  #late stress and recovery data in preferred order
data_combined_latestress_recovery_log2_mean_developing_seed = t(data_combined_latestress_recovery_log2_mean_developing_seed)

#log2-fc per cultivar
log2fc_combined_N22_latestress_recovery_developing_seed = data_combined_latestress_recovery_log2_mean_developing_seed[ , 2:4] - data_combined_latestress_recovery_log2_mean_developing_seed[ , 1]
log2fc_combined_Dular_latestress_recovery_developing_seed = data_combined_latestress_recovery_log2_mean_developing_seed[ , 6:8] - data_combined_latestress_recovery_log2_mean_developing_seed[ , 5]
log2fc_combined_Anjali_latestress_recovery_developing_seed = data_combined_latestress_recovery_log2_mean_developing_seed[ , 10:12] - data_combined_latestress_recovery_log2_mean_developing_seed[ , 9]

#combine log2-fc of all cultivars in one data set
log2fc_combined_latestress_recovery_developing_seed = cbind(log2fc_combined_N22_latestress_recovery_developing_seed, 
                                                            log2fc_combined_Dular_latestress_recovery_developing_seed, 
                                                            log2fc_combined_Anjali_latestress_recovery_developing_seed)
colnames(log2fc_combined_latestress_recovery_developing_seed)[1:3] = c(paste(colnames(log2fc_combined_latestress_recovery_developing_seed)[1:3], "- GNLS"))
colnames(log2fc_combined_latestress_recovery_developing_seed)[4:6] = c(paste(colnames(log2fc_combined_latestress_recovery_developing_seed)[4:6], "- GDLS"))
colnames(log2fc_combined_latestress_recovery_developing_seed)[7:9] = c(paste(colnames(log2fc_combined_latestress_recovery_developing_seed)[7:9], "- GALS"))

#overview
heatmap.2(log2fc_combined_latestress_recovery_developing_seed, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none",
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(7, 20), 
          main = "Developing seed - Combined - Log2 FC, Recovery/Late Stress",
          cexCol = 1, 
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          srtCol = 90)


#late stress vs. 12h rewatering

#N22 - 28 metabolites which have significant difference in LS/CLS comparison
data_combined_N22_log2median_stress_responsive_developing_seed = data_combined_N22_log2median_developing_seed
colnames(data_combined_N22_log2median_stress_responsive_developing_seed)[3:ncol(data_combined_N22_log2median_stress_responsive_developing_seed)] = reduced_metabolite_list_final_developing_seed_2013$Name
#only metabolites which have significant difference in LS/C comparison
data_combined_N22_log2median_stress_responsive_developing_seed = data_combined_N22_log2median_stress_responsive_developing_seed[ , c(1, 2, (which(colnames(data_combined_N22_log2median_stress_responsive_developing_seed) %in% rownames(wilcox_sig_log2fc_CLS_LS_developing_seed_combined))))]

#data distribution
shapiro_res_combined_stress_responsive_N22_latestress_recovery12_developing_seed = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_developing_seed[which(data_combined_N22_log2median_stress_responsive_developing_seed$Timepoint %in% c("Late stress", "12h Rewatering")), 3:ncol(data_combined_N22_log2median_stress_responsive_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_stress_responsive_N22_latestress_recovery12_developing_seed) = "p.value"
shapiro_res_combined_stress_responsive_N22_latestress_recovery12_developing_seed = as.data.frame(shapiro_res_combined_stress_responsive_N22_latestress_recovery12_developing_seed)
shapiro_res_combined_stress_responsive_N22_latestress_recovery12_developing_seed$Distribution = ifelse(shapiro_res_combined_stress_responsive_N22_latestress_recovery12_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_stress_responsive_N22_latestress_recovery12_developing_seed$Distribution == "Normal")     #23 metabolites
sum(shapiro_res_combined_stress_responsive_N22_latestress_recovery12_developing_seed$Distribution == "Non-normal") #5 metabolites

#wilcox test
wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_developing_seed[ , 3:ncol(data_combined_N22_log2median_stress_responsive_developing_seed)], function(x) wilcox.test(x ~ data_combined_N22_log2median_stress_responsive_developing_seed$Timepoint, subset = data_combined_N22_log2median_stress_responsive_developing_seed$Timepoint %in% c("Late stress", "12h Rewatering"))$p.value)))
colnames(wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined) = "p.value"
rownames(wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined) = rownames(wilcox_sig_log2fc_CLS_LS_developing_seed_combined)
wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined = as.data.frame(wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined)
wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined$Significance = with(wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined,
                                                                                     ifelse(p.value <= 0.001, "***", ifelse(
                                                                                       p.value <= 0.01, "**", ifelse(
                                                                                         p.value < 0.05, "*", "ns"))))
sum(wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined$Significance == "*") #1 metabolite
sum(wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined$Significance == "**") #2 metabolites
sum(wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined$Significance == "***") #0 metabolite
sum(wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined$Significance != "ns") #3 metabolites
sum(wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined$Significance == "ns") #25 metabolites


#Dular - 28 metabolites which have significant difference in LS/CLS comparison
data_combined_Dular_log2median_stress_responsive_developing_seed = data_combined_Dular_log2median_developing_seed
colnames(data_combined_Dular_log2median_stress_responsive_developing_seed)[3:ncol(data_combined_Dular_log2median_stress_responsive_developing_seed)] = reduced_metabolite_list_final_developing_seed_2013$Name
#only metabolites which have significant difference in LS/C comparison
data_combined_Dular_log2median_stress_responsive_developing_seed = data_combined_Dular_log2median_stress_responsive_developing_seed[ , c(1, 2, (which(colnames(data_combined_Dular_log2median_stress_responsive_developing_seed) %in% rownames(wilcox_sig_log2fc_CLS_LS_developing_seed_combined))))]

#data distribution
shapiro_res_combined_stress_responsive_Dular_latestress_recovery12_developing_seed = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_developing_seed[which(data_combined_Dular_log2median_stress_responsive_developing_seed$Timepoint %in% c("Late stress", "12h Rewatering")), 3:ncol(data_combined_Dular_log2median_stress_responsive_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_stress_responsive_Dular_latestress_recovery12_developing_seed) = "p.value"
shapiro_res_combined_stress_responsive_Dular_latestress_recovery12_developing_seed = as.data.frame(shapiro_res_combined_stress_responsive_Dular_latestress_recovery12_developing_seed)
shapiro_res_combined_stress_responsive_Dular_latestress_recovery12_developing_seed$Distribution = ifelse(shapiro_res_combined_stress_responsive_Dular_latestress_recovery12_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_stress_responsive_Dular_latestress_recovery12_developing_seed$Distribution == "Normal")     #19 metabolites
sum(shapiro_res_combined_stress_responsive_Dular_latestress_recovery12_developing_seed$Distribution == "Non-normal") #9 metabolites

#wilcox test
wilcox_res_stress_responsive_Dular_LS_R12_developing_seed_combined = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_developing_seed[ , 3:ncol(data_combined_Dular_log2median_stress_responsive_developing_seed)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_stress_responsive_developing_seed$Timepoint, subset = data_combined_Dular_log2median_stress_responsive_developing_seed$Timepoint %in% c("Late stress", "12h Rewatering"))$p.value)))
colnames(wilcox_res_stress_responsive_Dular_LS_R12_developing_seed_combined) = "p.value"
rownames(wilcox_res_stress_responsive_Dular_LS_R12_developing_seed_combined) = rownames(wilcox_sig_log2fc_CLS_LS_developing_seed_combined)
wilcox_res_stress_responsive_Dular_LS_R12_developing_seed_combined = as.data.frame(wilcox_res_stress_responsive_Dular_LS_R12_developing_seed_combined)
wilcox_res_stress_responsive_Dular_LS_R12_developing_seed_combined$Significance = with(wilcox_res_stress_responsive_Dular_LS_R12_developing_seed_combined,
                                                                                       ifelse(p.value <= 0.001, "***", ifelse(
                                                                                         p.value <= 0.01, "**", ifelse(
                                                                                           p.value < 0.05, "*", "ns"))))
sum(wilcox_res_stress_responsive_Dular_LS_R12_developing_seed_combined$Significance == "*") #1 metabolite
sum(wilcox_res_stress_responsive_Dular_LS_R12_developing_seed_combined$Significance == "**") #1 metabolite
sum(wilcox_res_stress_responsive_Dular_LS_R12_developing_seed_combined$Significance == "***") #2 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R12_developing_seed_combined$Significance != "ns") #4 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R12_developing_seed_combined$Significance == "ns") #24 metabolites


#Anjali - 28 metabolites which have significant difference in LS/CLS comparison
data_combined_Anjali_log2median_stress_responsive_developing_seed = data_combined_Anjali_log2median_developing_seed
colnames(data_combined_Anjali_log2median_stress_responsive_developing_seed)[3:ncol(data_combined_Anjali_log2median_stress_responsive_developing_seed)] = reduced_metabolite_list_final_developing_seed_2013$Name
#only metabolites which have significant difference in LS/C comparison
data_combined_Anjali_log2median_stress_responsive_developing_seed = data_combined_Anjali_log2median_stress_responsive_developing_seed[ , c(1, 2, (which(colnames(data_combined_Anjali_log2median_stress_responsive_developing_seed) %in% rownames(wilcox_sig_log2fc_CLS_LS_developing_seed_combined))))]

#data distribution
shapiro_res_combined_stress_responsive_Anjali_latestress_recovery12_developing_seed = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_developing_seed[which(data_combined_Anjali_log2median_stress_responsive_developing_seed$Timepoint %in% c("Late stress", "12h Rewatering")), 3:ncol(data_combined_Anjali_log2median_stress_responsive_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery12_developing_seed) = "p.value"
shapiro_res_combined_stress_responsive_Anjali_latestress_recovery12_developing_seed = as.data.frame(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery12_developing_seed)
shapiro_res_combined_stress_responsive_Anjali_latestress_recovery12_developing_seed$Distribution = ifelse(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery12_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery12_developing_seed$Distribution == "Normal")     #19 metabolites
sum(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery12_developing_seed$Distribution == "Non-normal") #9 metabolites

#wilcox test
wilcox_res_stress_responsive_Anjali_LS_R12_developing_seed_combined = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_developing_seed[ , 3:ncol(data_combined_Anjali_log2median_stress_responsive_developing_seed)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_stress_responsive_developing_seed$Timepoint, subset = data_combined_Anjali_log2median_stress_responsive_developing_seed$Timepoint %in% c("Late stress", "12h Rewatering"))$p.value)))
colnames(wilcox_res_stress_responsive_Anjali_LS_R12_developing_seed_combined) = "p.value"
rownames(wilcox_res_stress_responsive_Anjali_LS_R12_developing_seed_combined) = rownames(wilcox_sig_log2fc_CLS_LS_developing_seed_combined)
wilcox_res_stress_responsive_Anjali_LS_R12_developing_seed_combined = as.data.frame(wilcox_res_stress_responsive_Anjali_LS_R12_developing_seed_combined)
wilcox_res_stress_responsive_Anjali_LS_R12_developing_seed_combined$Significance = with(wilcox_res_stress_responsive_Anjali_LS_R12_developing_seed_combined,
                                                                                        ifelse(p.value <= 0.001, "***", ifelse(
                                                                                          p.value <= 0.01, "**", ifelse(
                                                                                            p.value < 0.05, "*", "ns"))))
sum(wilcox_res_stress_responsive_Anjali_LS_R12_developing_seed_combined$Significance == "*") #0 metabolite
sum(wilcox_res_stress_responsive_Anjali_LS_R12_developing_seed_combined$Significance == "**") #3 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R12_developing_seed_combined$Significance == "***") #1 metabolite
sum(wilcox_res_stress_responsive_Anjali_LS_R12_developing_seed_combined$Significance != "ns") #4 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R12_developing_seed_combined$Significance == "ns") #24 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_stress_responsive_LS_R12_developing_seed_combined = cbind(wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined,
                                                                     wilcox_res_stress_responsive_Dular_LS_R12_developing_seed_combined$Significance, wilcox_res_stress_responsive_Anjali_LS_R12_developing_seed_combined$Significance)
wilcox_sig_stress_responsive_LS_R12_developing_seed_combined = wilcox_sig_stress_responsive_LS_R12_developing_seed_combined[ , -1]
colnames(wilcox_sig_stress_responsive_LS_R12_developing_seed_combined)[1] = "wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined$Significance"
wilcox_sig_stress_responsive_LS_R12_developing_seed_combined = data.frame(lapply(wilcox_sig_stress_responsive_LS_R12_developing_seed_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_stress_responsive_LS_R12_developing_seed_combined) = rownames(wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined)
wilcox_sig_stress_responsive_LS_R12_developing_seed_combined$Nonsig = rowSums(wilcox_sig_stress_responsive_LS_R12_developing_seed_combined == " ")  #count number of nonsignificant
wilcox_sig_stress_responsive_final_LS_R12_developing_seed_combined = wilcox_sig_stress_responsive_LS_R12_developing_seed_combined[-which(wilcox_sig_stress_responsive_LS_R12_developing_seed_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 7 metabolites

#log2fc of metabolites significant in at least one of the cultivars
wilcox_sig_stress_responsive_log2fc_LS_R12_developing_seed_combined = log2fc_combined_latestress_recovery_developing_seed[which(rownames(log2fc_combined_latestress_recovery_developing_seed) %in% rownames(wilcox_sig_stress_responsive_final_LS_R12_developing_seed_combined)), c(1, 4, 7)]  

#overview - with HCL - Euclidean distance and average linkage
heatmap.2(wilcox_sig_stress_responsive_log2fc_LS_R12_developing_seed_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)),
          margins = c(5, 20), 
          main = "Developing seed-Combined-Wilcoxon-Stress-responsive \nLog2 FC, EGF, 12h RW/Late stress", 
          cexCol = 1.5, 
          cexRow = 0.9, 
          density.info = "none",
          key.xlab = NA,
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          labCol = c("N22", "Dular", "Anjali"), 
          srtCol = 0, 
          adjCol = c(0.5, 0.5), 
          cellnote = wilcox_sig_stress_responsive_final_LS_R12_developing_seed_combined, 
          notecol = "black",
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 3)), 
          rowsep = c(seq(0, 7)), 
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))


#venn diagram - common among cultivars

#N22
wilcox_sig_stress_responsive_N22_LS_R12_developing_seed_combined = wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined[-which(wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_N22_LS_R12_developing_seed_combined = as.data.frame(wilcox_sig_stress_responsive_log2fc_LS_R12_developing_seed_combined[match(rownames(wilcox_sig_stress_responsive_N22_LS_R12_developing_seed_combined), rownames(wilcox_sig_stress_responsive_log2fc_LS_R12_developing_seed_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_developing_seed_combined) = colnames(wilcox_sig_stress_responsive_log2fc_LS_R12_developing_seed_combined)[1]
wilcox_sig_stress_responsive_log2fc_N22_LS_R12_developing_seed_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_developing_seed_combined$`GNR12 - GNLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_developing_seed_combined$Up.Down == "Up")    #0 metabolite
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_developing_seed_combined$Up.Down == "Down")  #3 metabolites

#Dular
wilcox_sig_stress_responsive_Dular_LS_R12_developing_seed_combined = wilcox_res_stress_responsive_Dular_LS_R12_developing_seed_combined[-which(wilcox_res_stress_responsive_Dular_LS_R12_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_developing_seed_combined = as.data.frame(wilcox_sig_stress_responsive_log2fc_LS_R12_developing_seed_combined[match(rownames(wilcox_sig_stress_responsive_Dular_LS_R12_developing_seed_combined), rownames(wilcox_sig_stress_responsive_log2fc_LS_R12_developing_seed_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_developing_seed_combined) = colnames(wilcox_sig_stress_responsive_log2fc_LS_R12_developing_seed_combined)[2]
wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_developing_seed_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_developing_seed_combined$`GDR12 - GDLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_developing_seed_combined$Up.Down == "Up")    #0 metabolite
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_developing_seed_combined$Up.Down == "Down")  #4 metabolites

#Anjali
wilcox_sig_stress_responsive_Anjali_LS_R12_developing_seed_combined = wilcox_res_stress_responsive_Anjali_LS_R12_developing_seed_combined[-which(wilcox_res_stress_responsive_Anjali_LS_R12_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_developing_seed_combined = as.data.frame(wilcox_sig_stress_responsive_log2fc_LS_R12_developing_seed_combined[match(rownames(wilcox_sig_stress_responsive_Anjali_LS_R12_developing_seed_combined), rownames(wilcox_sig_stress_responsive_log2fc_LS_R12_developing_seed_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_developing_seed_combined) = colnames(wilcox_sig_stress_responsive_log2fc_LS_R12_developing_seed_combined)[3]
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_developing_seed_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_developing_seed_combined$`GAR12 - GALS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_developing_seed_combined$Up.Down == "Up")    #1 metabolites
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_developing_seed_combined$Up.Down == "Down")  #3 metabolites


#venn - increased - 12h rewatering/late stress
vennlist_up_wilcox_stress_responsive_LS_R12_developing_seed_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_developing_seed_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R12_developing_seed_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_developing_seed_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_developing_seed_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_developing_seed_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_developing_seed_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_stress_responsive_LS_R12_developing_seed_combined) = c("N22", "Dular", "Anjali")
venn_up_stress_responsive_wilcox_LS_R12_developing_seed_combined = venn(vennlist_up_wilcox_stress_responsive_LS_R12_developing_seed_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_stress_responsive_R12_latestress_developing_seed_combined = venn.diagram(vennlist_up_wilcox_stress_responsive_LS_R12_developing_seed_combined, category.names = c(paste(names(vennlist_up_wilcox_stress_responsive_LS_R12_developing_seed_combined[1]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R12_developing_seed_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R12_developing_seed_combined[2]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R12_developing_seed_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R12_developing_seed_combined[3]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R12_developing_seed_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "A", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label A
intersect_stress_responsive_wilcox_LS_R12_venn_up_developing_seed_combined = attr(venn_up_stress_responsive_wilcox_LS_R12_developing_seed_combined, "intersections")  #list of unique and common metabolites
venn_up_stress_responsive_wilcox_metabolites_LS_R12_developing_seed_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R12_venn_up_developing_seed_combined, rbind))  #convert from list to matrix
colnames(venn_up_stress_responsive_wilcox_metabolites_LS_R12_developing_seed_combined) = venn_up_stress_responsive_wilcox_metabolites_LS_R12_developing_seed_combined[1, ]
venn_up_stress_responsive_wilcox_metabolites_LS_R12_developing_seed_combined = venn_up_stress_responsive_wilcox_metabolites_LS_R12_developing_seed_combined[-1, ]


#venn - decreased - 12h rewatering/late stress
vennlist_down_wilcox_stress_responsive_LS_R12_developing_seed_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_developing_seed_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R12_developing_seed_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_developing_seed_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_developing_seed_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_developing_seed_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_developing_seed_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_stress_responsive_LS_R12_developing_seed_combined) = c("N22", "Dular", "Anjali")
venn_down_stress_responsive_wilcox_LS_R12_developing_seed_combined = venn(vennlist_down_wilcox_stress_responsive_LS_R12_developing_seed_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_stress_responsive_R12_latestress_developing_seed_combined = venn.diagram(vennlist_down_wilcox_stress_responsive_LS_R12_developing_seed_combined, category.names = c(paste(names(vennlist_down_wilcox_stress_responsive_LS_R12_developing_seed_combined[1]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R12_developing_seed_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R12_developing_seed_combined[2]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R12_developing_seed_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R12_developing_seed_combined[3]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R12_developing_seed_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "B", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label B
intersect_stress_responsive_wilcox_LS_R12_venn_down_developing_seed_combined = attr(venn_down_stress_responsive_wilcox_LS_R12_developing_seed_combined, "intersections")  #list of unique and common metabolites
venn_down_stress_responsive_wilcox_metabolites_LS_R12_developing_seed_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R12_venn_down_developing_seed_combined, rbind))  #convert from list to matrix
colnames(venn_down_stress_responsive_wilcox_metabolites_LS_R12_developing_seed_combined) = venn_down_stress_responsive_wilcox_metabolites_LS_R12_developing_seed_combined[1, ]
venn_down_stress_responsive_wilcox_metabolites_LS_R12_developing_seed_combined = venn_down_stress_responsive_wilcox_metabolites_LS_R12_developing_seed_combined[-1, ]


write.table(venn_up_stress_responsive_wilcox_metabolites_LS_R12_developing_seed_combined, file = "wilcox_stress_responsive_sig_up_metabolites_developing_seed_R12_latestress_combined.txt", sep = "\t", quote = F)
write.table(venn_down_stress_responsive_wilcox_metabolites_LS_R12_developing_seed_combined, file = "wilcox_stress_responsive_sig_down_metabolites_developing_seed_R12_latestress_combined.txt", sep = "\t", quote = F)


####################


#late stress vs. 36h rewatering

#N22

#data distribution
shapiro_res_combined_stress_responsive_N22_latestress_recovery36_developing_seed = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_developing_seed[which(data_combined_N22_log2median_stress_responsive_developing_seed$Timepoint %in% c("Late stress", "36h Rewatering")), 3:ncol(data_combined_N22_log2median_stress_responsive_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_stress_responsive_N22_latestress_recovery36_developing_seed) = "p.value"
shapiro_res_combined_stress_responsive_N22_latestress_recovery36_developing_seed = as.data.frame(shapiro_res_combined_stress_responsive_N22_latestress_recovery36_developing_seed)
shapiro_res_combined_stress_responsive_N22_latestress_recovery36_developing_seed$Distribution = ifelse(shapiro_res_combined_stress_responsive_N22_latestress_recovery36_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_stress_responsive_N22_latestress_recovery36_developing_seed$Distribution == "Normal")     #23 metabolites
sum(shapiro_res_combined_stress_responsive_N22_latestress_recovery36_developing_seed$Distribution == "Non-normal") #5 metabolites

#wilcox test
wilcox_res_stress_responsive_N22_LS_R36_developing_seed_combined = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_developing_seed[ , 3:ncol(data_combined_N22_log2median_stress_responsive_developing_seed)], function(x) wilcox.test(x ~ data_combined_N22_log2median_stress_responsive_developing_seed$Timepoint, subset = data_combined_N22_log2median_stress_responsive_developing_seed$Timepoint %in% c("Late stress", "36h Rewatering"))$p.value)))
colnames(wilcox_res_stress_responsive_N22_LS_R36_developing_seed_combined) = "p.value"
rownames(wilcox_res_stress_responsive_N22_LS_R36_developing_seed_combined) = rownames(wilcox_sig_log2fc_CLS_LS_developing_seed_combined)
wilcox_res_stress_responsive_N22_LS_R36_developing_seed_combined = as.data.frame(wilcox_res_stress_responsive_N22_LS_R36_developing_seed_combined)
wilcox_res_stress_responsive_N22_LS_R36_developing_seed_combined$Significance = with(wilcox_res_stress_responsive_N22_LS_R36_developing_seed_combined,
                                                                                     ifelse(p.value <= 0.001, "***", ifelse(
                                                                                       p.value <= 0.01, "**", ifelse(
                                                                                         p.value < 0.05, "*", "ns"))))
sum(wilcox_res_stress_responsive_N22_LS_R36_developing_seed_combined$Significance == "*") #4 metabolites
sum(wilcox_res_stress_responsive_N22_LS_R36_developing_seed_combined$Significance == "**") #2 metabolites
sum(wilcox_res_stress_responsive_N22_LS_R36_developing_seed_combined$Significance == "***") #7 metabolite
sum(wilcox_res_stress_responsive_N22_LS_R36_developing_seed_combined$Significance != "ns") #13 metabolites
sum(wilcox_res_stress_responsive_N22_LS_R36_developing_seed_combined$Significance == "ns") #15 metabolites


#Dular

#data distribution
shapiro_res_combined_stress_responsive_Dular_latestress_recovery36_developing_seed = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_developing_seed[which(data_combined_Dular_log2median_stress_responsive_developing_seed$Timepoint %in% c("Late stress", "36h Rewatering")), 3:ncol(data_combined_Dular_log2median_stress_responsive_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_stress_responsive_Dular_latestress_recovery36_developing_seed) = "p.value"
shapiro_res_combined_stress_responsive_Dular_latestress_recovery36_developing_seed = as.data.frame(shapiro_res_combined_stress_responsive_Dular_latestress_recovery36_developing_seed)
shapiro_res_combined_stress_responsive_Dular_latestress_recovery36_developing_seed$Distribution = ifelse(shapiro_res_combined_stress_responsive_Dular_latestress_recovery36_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_stress_responsive_Dular_latestress_recovery36_developing_seed$Distribution == "Normal")     #23 metabolites
sum(shapiro_res_combined_stress_responsive_Dular_latestress_recovery36_developing_seed$Distribution == "Non-normal") #5 metabolites

#wilcox test
wilcox_res_stress_responsive_Dular_LS_R36_developing_seed_combined = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_developing_seed[ , 3:ncol(data_combined_Dular_log2median_stress_responsive_developing_seed)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_stress_responsive_developing_seed$Timepoint, subset = data_combined_Dular_log2median_stress_responsive_developing_seed$Timepoint %in% c("Late stress", "36h Rewatering"))$p.value)))
colnames(wilcox_res_stress_responsive_Dular_LS_R36_developing_seed_combined) = "p.value"
rownames(wilcox_res_stress_responsive_Dular_LS_R36_developing_seed_combined) = rownames(wilcox_sig_log2fc_CLS_LS_developing_seed_combined)
wilcox_res_stress_responsive_Dular_LS_R36_developing_seed_combined = as.data.frame(wilcox_res_stress_responsive_Dular_LS_R36_developing_seed_combined)
wilcox_res_stress_responsive_Dular_LS_R36_developing_seed_combined$Significance = with(wilcox_res_stress_responsive_Dular_LS_R36_developing_seed_combined,
                                                                                       ifelse(p.value <= 0.001, "***", ifelse(
                                                                                         p.value <= 0.01, "**", ifelse(
                                                                                           p.value < 0.05, "*", "ns"))))
sum(wilcox_res_stress_responsive_Dular_LS_R36_developing_seed_combined$Significance == "*") #2 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R36_developing_seed_combined$Significance == "**") #7 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R36_developing_seed_combined$Significance == "***") #3 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R36_developing_seed_combined$Significance != "ns") #12 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R36_developing_seed_combined$Significance == "ns") #16 metabolites


#Anjali

#data distribution
shapiro_res_combined_stress_responsive_Anjali_latestress_recovery36_developing_seed = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_developing_seed[which(data_combined_Anjali_log2median_stress_responsive_developing_seed$Timepoint %in% c("Late stress", "36h Rewatering")), 3:ncol(data_combined_Anjali_log2median_stress_responsive_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery36_developing_seed) = "p.value"
shapiro_res_combined_stress_responsive_Anjali_latestress_recovery36_developing_seed = as.data.frame(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery36_developing_seed)
shapiro_res_combined_stress_responsive_Anjali_latestress_recovery36_developing_seed$Distribution = ifelse(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery36_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery36_developing_seed$Distribution == "Normal")     #17 metabolites
sum(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery36_developing_seed$Distribution == "Non-normal") #11 metabolites

#wilcox test
wilcox_res_stress_responsive_Anjali_LS_R36_developing_seed_combined = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_developing_seed[ , 3:ncol(data_combined_Anjali_log2median_stress_responsive_developing_seed)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_stress_responsive_developing_seed$Timepoint, subset = data_combined_Anjali_log2median_stress_responsive_developing_seed$Timepoint %in% c("Late stress", "36h Rewatering"))$p.value)))
colnames(wilcox_res_stress_responsive_Anjali_LS_R36_developing_seed_combined) = "p.value"
rownames(wilcox_res_stress_responsive_Anjali_LS_R36_developing_seed_combined) = rownames(wilcox_sig_log2fc_CLS_LS_developing_seed_combined)
wilcox_res_stress_responsive_Anjali_LS_R36_developing_seed_combined = as.data.frame(wilcox_res_stress_responsive_Anjali_LS_R36_developing_seed_combined)
wilcox_res_stress_responsive_Anjali_LS_R36_developing_seed_combined$Significance = with(wilcox_res_stress_responsive_Anjali_LS_R36_developing_seed_combined,
                                                                                        ifelse(p.value <= 0.001, "***", ifelse(
                                                                                          p.value <= 0.01, "**", ifelse(
                                                                                            p.value < 0.05, "*", "ns"))))
sum(wilcox_res_stress_responsive_Anjali_LS_R36_developing_seed_combined$Significance == "*") #3 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R36_developing_seed_combined$Significance == "**") #2 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R36_developing_seed_combined$Significance == "***") #6 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R36_developing_seed_combined$Significance != "ns") #11 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R36_developing_seed_combined$Significance == "ns") #17 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_stress_responsive_LS_R36_developing_seed_combined = cbind(wilcox_res_stress_responsive_N22_LS_R36_developing_seed_combined,
                                                                     wilcox_res_stress_responsive_Dular_LS_R36_developing_seed_combined$Significance, 
                                                                     wilcox_res_stress_responsive_Anjali_LS_R36_developing_seed_combined$Significance)
wilcox_sig_stress_responsive_LS_R36_developing_seed_combined = wilcox_sig_stress_responsive_LS_R36_developing_seed_combined[ , -1]
colnames(wilcox_sig_stress_responsive_LS_R36_developing_seed_combined)[1] = "wilcox_res_stress_responsive_N22_LS_R36_developing_seed_combined$Significance"
wilcox_sig_stress_responsive_LS_R36_developing_seed_combined = data.frame(lapply(wilcox_sig_stress_responsive_LS_R36_developing_seed_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_stress_responsive_LS_R36_developing_seed_combined) = rownames(wilcox_res_stress_responsive_N22_LS_R36_developing_seed_combined)
wilcox_sig_stress_responsive_LS_R36_developing_seed_combined$Nonsig = rowSums(wilcox_sig_stress_responsive_LS_R36_developing_seed_combined == " ")  #count number of nonsignificant
wilcox_sig_stress_responsive_final_LS_R36_developing_seed_combined = wilcox_sig_stress_responsive_LS_R36_developing_seed_combined[-which(wilcox_sig_stress_responsive_LS_R36_developing_seed_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 17 metabolites

#log2fc of metabolites significant in at least one of the cultivars
wilcox_sig_stress_responsive_log2fc_LS_R36_developing_seed_combined = log2fc_combined_latestress_recovery_developing_seed[which(rownames(log2fc_combined_latestress_recovery_developing_seed) %in% rownames(wilcox_sig_stress_responsive_final_LS_R36_developing_seed_combined)), c(2, 5, 8)]  

#overview - with HCL - Euclidean distance and average linkage
heatmap.2(wilcox_sig_stress_responsive_log2fc_LS_R36_developing_seed_combined, 
          Rowv = TRUE,
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599),
          breaks = c(seq(-3, -0.01, length.out = 300), seq(0.01, 3, length.out = 300)), 
          margins = c(5, 20), 
          main = "Developing seed-Combined-Wilcoxon-Stress-responsive \nLog2 FC, EGF, 36h RW/Late stress", 
          cexCol = 1.5, 
          cexRow = 0.9,
          density.info = "none",
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          labCol = c("N22", "Dular", "Anjali"),
          srtCol = 0, 
          adjCol = c(0.5, 0.5), 
          cellnote = wilcox_sig_stress_responsive_final_LS_R36_developing_seed_combined, 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 3)), 
          rowsep = c(seq(0, 17)),
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))


#venn diagram - common among cultivars

#N22
wilcox_sig_stress_responsive_N22_LS_R36_developing_seed_combined = wilcox_res_stress_responsive_N22_LS_R36_developing_seed_combined[-which(wilcox_res_stress_responsive_N22_LS_R36_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_N22_LS_R36_developing_seed_combined = as.data.frame(wilcox_sig_stress_responsive_log2fc_LS_R36_developing_seed_combined[match(rownames(wilcox_sig_stress_responsive_N22_LS_R36_developing_seed_combined), rownames(wilcox_sig_stress_responsive_log2fc_LS_R36_developing_seed_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_developing_seed_combined) = colnames(wilcox_sig_stress_responsive_log2fc_LS_R36_developing_seed_combined)[1]
wilcox_sig_stress_responsive_log2fc_N22_LS_R36_developing_seed_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_developing_seed_combined$`GNR36 - GNLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_developing_seed_combined$Up.Down == "Up")    #1 metabolite
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_developing_seed_combined$Up.Down == "Down")  #12 metabolites

#Dular
wilcox_sig_stress_responsive_Dular_LS_R36_developing_seed_combined = wilcox_res_stress_responsive_Dular_LS_R36_developing_seed_combined[-which(wilcox_res_stress_responsive_Dular_LS_R36_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_developing_seed_combined = as.data.frame(wilcox_sig_stress_responsive_log2fc_LS_R36_developing_seed_combined[match(rownames(wilcox_sig_stress_responsive_Dular_LS_R36_developing_seed_combined), rownames(wilcox_sig_stress_responsive_log2fc_LS_R36_developing_seed_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_developing_seed_combined) = colnames(wilcox_sig_stress_responsive_log2fc_LS_R36_developing_seed_combined)[2]
wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_developing_seed_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_developing_seed_combined$`GDR36 - GDLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_developing_seed_combined$Up.Down == "Up")    #1 metabolite
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_developing_seed_combined$Up.Down == "Down")  #11 metabolites

#Anjali
wilcox_sig_stress_responsive_Anjali_LS_R36_developing_seed_combined = wilcox_res_stress_responsive_Anjali_LS_R36_developing_seed_combined[-which(wilcox_res_stress_responsive_Anjali_LS_R36_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_developing_seed_combined = as.data.frame(wilcox_sig_stress_responsive_log2fc_LS_R36_developing_seed_combined[match(rownames(wilcox_sig_stress_responsive_Anjali_LS_R36_developing_seed_combined), rownames(wilcox_sig_stress_responsive_log2fc_LS_R36_developing_seed_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_developing_seed_combined) = colnames(wilcox_sig_stress_responsive_log2fc_LS_R36_developing_seed_combined)[3]
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_developing_seed_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_developing_seed_combined$`GAR36 - GALS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_developing_seed_combined$Up.Down == "Up")    #3 metabolites
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_developing_seed_combined$Up.Down == "Down")  #8 metabolites


#venn - increased - 36h rewatering/late stress
vennlist_up_wilcox_stress_responsive_LS_R36_developing_seed_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_developing_seed_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R36_developing_seed_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_developing_seed_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_developing_seed_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_developing_seed_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_developing_seed_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_stress_responsive_LS_R36_developing_seed_combined) = c("N22", "Dular", "Anjali")
venn_up_stress_responsive_wilcox_LS_R36_developing_seed_combined = venn(vennlist_up_wilcox_stress_responsive_LS_R36_developing_seed_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_stress_responsive_R36_latestress_developing_seed_combined = venn.diagram(vennlist_up_wilcox_stress_responsive_LS_R36_developing_seed_combined, category.names = c(paste(names(vennlist_up_wilcox_stress_responsive_LS_R36_developing_seed_combined[1]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R36_developing_seed_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R36_developing_seed_combined[2]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R36_developing_seed_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R36_developing_seed_combined[3]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R36_developing_seed_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "C", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label C
intersect_stress_responsive_wilcox_LS_R36_venn_up_developing_seed_combined = attr(venn_up_stress_responsive_wilcox_LS_R36_developing_seed_combined, "intersections")  #list of unique and common metabolites
venn_up_stress_responsive_wilcox_metabolites_LS_R36_developing_seed_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R36_venn_up_developing_seed_combined, rbind))  #convert from list to matrix
colnames(venn_up_stress_responsive_wilcox_metabolites_LS_R36_developing_seed_combined) = venn_up_stress_responsive_wilcox_metabolites_LS_R36_developing_seed_combined[1, ]
venn_up_stress_responsive_wilcox_metabolites_LS_R36_developing_seed_combined = venn_up_stress_responsive_wilcox_metabolites_LS_R36_developing_seed_combined[-1, ]


#venn - decreased - 36h rewatering/late stress
vennlist_down_wilcox_stress_responsive_LS_R36_developing_seed_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_developing_seed_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R36_developing_seed_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_developing_seed_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_developing_seed_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_developing_seed_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_developing_seed_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_stress_responsive_LS_R36_developing_seed_combined) = c("N22", "Dular", "Anjali")
venn_down_stress_responsive_wilcox_LS_R36_developing_seed_combined = venn(vennlist_down_wilcox_stress_responsive_LS_R36_developing_seed_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_stress_responsive_R36_latestress_developing_seed_combined = venn.diagram(vennlist_down_wilcox_stress_responsive_LS_R36_developing_seed_combined, category.names = c(paste(names(vennlist_down_wilcox_stress_responsive_LS_R36_developing_seed_combined[1]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R36_developing_seed_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R36_developing_seed_combined[2]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R36_developing_seed_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R36_developing_seed_combined[3]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R36_developing_seed_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "D", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label D
intersect_stress_responsive_wilcox_LS_R36_venn_down_developing_seed_combined = attr(venn_down_stress_responsive_wilcox_LS_R36_developing_seed_combined, "intersections")  #list of unique and common metabolites
venn_down_stress_responsive_wilcox_metabolites_LS_R36_developing_seed_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R36_venn_down_developing_seed_combined, rbind))  #convert from list to matrix
colnames(venn_down_stress_responsive_wilcox_metabolites_LS_R36_developing_seed_combined) = venn_down_stress_responsive_wilcox_metabolites_LS_R36_developing_seed_combined[1, ]
venn_down_stress_responsive_wilcox_metabolites_LS_R36_developing_seed_combined = venn_down_stress_responsive_wilcox_metabolites_LS_R36_developing_seed_combined[-1, ]


write.table(venn_up_stress_responsive_wilcox_metabolites_LS_R36_developing_seed_combined, file = "wilcox_stress_responsive_sig_up_metabolites_developing_seed_R36_latestress_combined.txt", sep = "\t", quote = F)
write.table(venn_down_stress_responsive_wilcox_metabolites_LS_R36_developing_seed_combined, file = "wilcox_stress_responsive_sig_down_metabolites_developing_seed_R36_latestress_combined.txt", sep = "\t", quote = F)


####################


#late stress vs. 60h rewatering

#N22

#data distribution
shapiro_res_combined_stress_responsive_N22_latestress_recovery60_developing_seed = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_developing_seed[which(data_combined_N22_log2median_stress_responsive_developing_seed$Timepoint %in% c("Late stress", "60h Rewatering")), 3:ncol(data_combined_N22_log2median_stress_responsive_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_stress_responsive_N22_latestress_recovery60_developing_seed) = "p.value"
shapiro_res_combined_stress_responsive_N22_latestress_recovery60_developing_seed = as.data.frame(shapiro_res_combined_stress_responsive_N22_latestress_recovery60_developing_seed)
shapiro_res_combined_stress_responsive_N22_latestress_recovery60_developing_seed$Distribution = ifelse(shapiro_res_combined_stress_responsive_N22_latestress_recovery60_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_stress_responsive_N22_latestress_recovery60_developing_seed$Distribution == "Normal")     #21 metabolites
sum(shapiro_res_combined_stress_responsive_N22_latestress_recovery60_developing_seed$Distribution == "Non-normal") #7 metabolites

#wilcox test
wilcox_res_stress_responsive_N22_LS_R60_developing_seed_combined = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_developing_seed[ , 3:ncol(data_combined_N22_log2median_stress_responsive_developing_seed)], function(x) wilcox.test(x ~ data_combined_N22_log2median_stress_responsive_developing_seed$Timepoint, subset = data_combined_N22_log2median_stress_responsive_developing_seed$Timepoint %in% c("Late stress", "60h Rewatering"))$p.value)))
colnames(wilcox_res_stress_responsive_N22_LS_R60_developing_seed_combined) = "p.value"
rownames(wilcox_res_stress_responsive_N22_LS_R60_developing_seed_combined) = rownames(wilcox_sig_log2fc_CLS_LS_developing_seed_combined)
wilcox_res_stress_responsive_N22_LS_R60_developing_seed_combined = as.data.frame(wilcox_res_stress_responsive_N22_LS_R60_developing_seed_combined)
wilcox_res_stress_responsive_N22_LS_R60_developing_seed_combined$Significance = with(wilcox_res_stress_responsive_N22_LS_R60_developing_seed_combined,
                                                                                     ifelse(p.value <= 0.001, "***", ifelse(
                                                                                       p.value <= 0.01, "**", ifelse(
                                                                                         p.value < 0.05, "*", "ns"))))
sum(wilcox_res_stress_responsive_N22_LS_R60_developing_seed_combined$Significance == "*") #1 metabolite
sum(wilcox_res_stress_responsive_N22_LS_R60_developing_seed_combined$Significance == "**") #4 metabolites
sum(wilcox_res_stress_responsive_N22_LS_R60_developing_seed_combined$Significance == "***") #8 metabolites
sum(wilcox_res_stress_responsive_N22_LS_R60_developing_seed_combined$Significance != "ns") #13 metabolites
sum(wilcox_res_stress_responsive_N22_LS_R60_developing_seed_combined$Significance == "ns") #15 metabolites


#Dular

#data distribution
shapiro_res_combined_stress_responsive_Dular_latestress_recovery60_developing_seed = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_developing_seed[which(data_combined_Dular_log2median_stress_responsive_developing_seed$Timepoint %in% c("Late stress", "60h Rewatering")), 3:ncol(data_combined_Dular_log2median_stress_responsive_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_stress_responsive_Dular_latestress_recovery60_developing_seed) = "p.value"
shapiro_res_combined_stress_responsive_Dular_latestress_recovery60_developing_seed = as.data.frame(shapiro_res_combined_stress_responsive_Dular_latestress_recovery60_developing_seed)
shapiro_res_combined_stress_responsive_Dular_latestress_recovery60_developing_seed$Distribution = ifelse(shapiro_res_combined_stress_responsive_Dular_latestress_recovery60_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_stress_responsive_Dular_latestress_recovery60_developing_seed$Distribution == "Normal")     #22 metabolites
sum(shapiro_res_combined_stress_responsive_Dular_latestress_recovery60_developing_seed$Distribution == "Non-normal") #6 metabolites

#wilcox test
wilcox_res_stress_responsive_Dular_LS_R60_developing_seed_combined = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_developing_seed[ , 3:ncol(data_combined_Dular_log2median_stress_responsive_developing_seed)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_stress_responsive_developing_seed$Timepoint, subset = data_combined_Dular_log2median_stress_responsive_developing_seed$Timepoint %in% c("Late stress", "60h Rewatering"))$p.value)))
colnames(wilcox_res_stress_responsive_Dular_LS_R60_developing_seed_combined) = "p.value"
rownames(wilcox_res_stress_responsive_Dular_LS_R60_developing_seed_combined) = rownames(wilcox_sig_log2fc_CLS_LS_developing_seed_combined)
wilcox_res_stress_responsive_Dular_LS_R60_developing_seed_combined = as.data.frame(wilcox_res_stress_responsive_Dular_LS_R60_developing_seed_combined)
wilcox_res_stress_responsive_Dular_LS_R60_developing_seed_combined$Significance = with(wilcox_res_stress_responsive_Dular_LS_R60_developing_seed_combined,
                                                                                       ifelse(p.value <= 0.001, "***", ifelse(
                                                                                         p.value <= 0.01, "**", ifelse(
                                                                                           p.value < 0.05, "*", "ns"))))
sum(wilcox_res_stress_responsive_Dular_LS_R60_developing_seed_combined$Significance == "*") #2 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R60_developing_seed_combined$Significance == "**") #6 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R60_developing_seed_combined$Significance == "***") #12 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R60_developing_seed_combined$Significance != "ns") #20 metabolites
sum(wilcox_res_stress_responsive_Dular_LS_R60_developing_seed_combined$Significance == "ns") #8 metabolites


#Anjali

#data distribution
shapiro_res_combined_stress_responsive_Anjali_latestress_recovery60_developing_seed = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_developing_seed[which(data_combined_Anjali_log2median_stress_responsive_developing_seed$Timepoint %in% c("Late stress", "60h Rewatering")), 3:ncol(data_combined_Anjali_log2median_stress_responsive_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery60_developing_seed) = "p.value"
shapiro_res_combined_stress_responsive_Anjali_latestress_recovery60_developing_seed = as.data.frame(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery60_developing_seed)
shapiro_res_combined_stress_responsive_Anjali_latestress_recovery60_developing_seed$Distribution = ifelse(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery60_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery60_developing_seed$Distribution == "Normal")     #17 metabolites
sum(shapiro_res_combined_stress_responsive_Anjali_latestress_recovery60_developing_seed$Distribution == "Non-normal") #11 metabolites

#wilcox test
wilcox_res_stress_responsive_Anjali_LS_R60_developing_seed_combined = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_developing_seed[ , 3:ncol(data_combined_Anjali_log2median_stress_responsive_developing_seed)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_stress_responsive_developing_seed$Timepoint, subset = data_combined_Anjali_log2median_stress_responsive_developing_seed$Timepoint %in% c("Late stress", "60h Rewatering"))$p.value)))
colnames(wilcox_res_stress_responsive_Anjali_LS_R60_developing_seed_combined) = "p.value"
rownames(wilcox_res_stress_responsive_Anjali_LS_R60_developing_seed_combined) = rownames(wilcox_sig_log2fc_CLS_LS_developing_seed_combined)
wilcox_res_stress_responsive_Anjali_LS_R60_developing_seed_combined = as.data.frame(wilcox_res_stress_responsive_Anjali_LS_R60_developing_seed_combined)
wilcox_res_stress_responsive_Anjali_LS_R60_developing_seed_combined$Significance = with(wilcox_res_stress_responsive_Anjali_LS_R60_developing_seed_combined,
                                                                                        ifelse(p.value <= 0.001, "***", ifelse(
                                                                                          p.value <= 0.01, "**", ifelse(
                                                                                            p.value < 0.05, "*", "ns"))))
sum(wilcox_res_stress_responsive_Anjali_LS_R60_developing_seed_combined$Significance == "*") #2 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R60_developing_seed_combined$Significance == "**") #4 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R60_developing_seed_combined$Significance == "***") #9 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R60_developing_seed_combined$Significance != "ns") #15 metabolites
sum(wilcox_res_stress_responsive_Anjali_LS_R60_developing_seed_combined$Significance == "ns") #13 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_stress_responsive_LS_R60_developing_seed_combined = cbind(wilcox_res_stress_responsive_N22_LS_R60_developing_seed_combined,
                                                                     wilcox_res_stress_responsive_Dular_LS_R60_developing_seed_combined$Significance, 
                                                                     wilcox_res_stress_responsive_Anjali_LS_R60_developing_seed_combined$Significance)
wilcox_sig_stress_responsive_LS_R60_developing_seed_combined = wilcox_sig_stress_responsive_LS_R60_developing_seed_combined[ , -1]
colnames(wilcox_sig_stress_responsive_LS_R60_developing_seed_combined)[1] = "wilcox_res_stress_responsive_N22_LS_R60_developing_seed_combined$Significance"
wilcox_sig_stress_responsive_LS_R60_developing_seed_combined = data.frame(lapply(wilcox_sig_stress_responsive_LS_R60_developing_seed_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_stress_responsive_LS_R60_developing_seed_combined) = rownames(wilcox_res_stress_responsive_N22_LS_R60_developing_seed_combined)
wilcox_sig_stress_responsive_LS_R60_developing_seed_combined$Nonsig = rowSums(wilcox_sig_stress_responsive_LS_R60_developing_seed_combined == " ")  #count number of nonsignificant
wilcox_sig_stress_responsive_final_LS_R60_developing_seed_combined = wilcox_sig_stress_responsive_LS_R60_developing_seed_combined[-which(wilcox_sig_stress_responsive_LS_R60_developing_seed_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 22 metabolites

#log2fc of metabolites significant in at least one of the cultivars
wilcox_sig_stress_responsive_log2fc_LS_R60_developing_seed_combined = log2fc_combined_latestress_recovery_developing_seed[which(rownames(log2fc_combined_latestress_recovery_developing_seed) %in% rownames(wilcox_sig_stress_responsive_final_LS_R60_developing_seed_combined)), c(3, 6, 9)]  

#overview - with HCL - Euclidean distance and average linkage
heatmap.2(wilcox_sig_stress_responsive_log2fc_LS_R60_developing_seed_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-3, -0.01, length.out = 300), seq(0.01, 3, length.out = 300)), 
          margins = c(5, 20),  
          main = "Developing seed-Combined-Wilcoxon-Stress-responsive \nLog2 FC, EGF, 60h RW/Late stress", 
          cexCol = 1.5, 
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA,
          lhei = c(1, 8),
          lwid = c(1.5, 5), 
          labCol = c("N22", "Dular", "Anjali"), 
          srtCol = 0, 
          adjCol = c(0.5, 0.5),
          cellnote = wilcox_sig_stress_responsive_final_LS_R60_developing_seed_combined,
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 3)), 
          rowsep = c(seq(0, 22)),
          sepwidth = c(0.01, 0.01),
          hclustfun = function(x) hclust(x, method = "average"))


#venn diagram - common among cultivars

#N22
wilcox_sig_stress_responsive_N22_LS_R60_developing_seed_combined = wilcox_res_stress_responsive_N22_LS_R60_developing_seed_combined[-which(wilcox_res_stress_responsive_N22_LS_R60_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_N22_LS_R60_developing_seed_combined = as.data.frame(wilcox_sig_stress_responsive_log2fc_LS_R60_developing_seed_combined[match(rownames(wilcox_sig_stress_responsive_N22_LS_R60_developing_seed_combined), rownames(wilcox_sig_stress_responsive_log2fc_LS_R60_developing_seed_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_developing_seed_combined) = colnames(wilcox_sig_stress_responsive_log2fc_LS_R60_developing_seed_combined)[1]
wilcox_sig_stress_responsive_log2fc_N22_LS_R60_developing_seed_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_developing_seed_combined$`GNR60 - GNLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_developing_seed_combined$Up.Down == "Up")    #1 metabolite
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_developing_seed_combined$Up.Down == "Down")  #12 metabolites

#Dular
wilcox_sig_stress_responsive_Dular_LS_R60_developing_seed_combined = wilcox_res_stress_responsive_Dular_LS_R60_developing_seed_combined[-which(wilcox_res_stress_responsive_Dular_LS_R60_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_developing_seed_combined = as.data.frame(wilcox_sig_stress_responsive_log2fc_LS_R60_developing_seed_combined[match(rownames(wilcox_sig_stress_responsive_Dular_LS_R60_developing_seed_combined), rownames(wilcox_sig_stress_responsive_log2fc_LS_R60_developing_seed_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_developing_seed_combined) = colnames(wilcox_sig_stress_responsive_log2fc_LS_R60_developing_seed_combined)[2]
wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_developing_seed_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_developing_seed_combined$`GDR60 - GDLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_developing_seed_combined$Up.Down == "Up")    #2 metabolites
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_developing_seed_combined$Up.Down == "Down")  #18 metabolites

#Anjali
wilcox_sig_stress_responsive_Anjali_LS_R60_developing_seed_combined = wilcox_res_stress_responsive_Anjali_LS_R60_developing_seed_combined[-which(wilcox_res_stress_responsive_Anjali_LS_R60_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_developing_seed_combined = as.data.frame(wilcox_sig_stress_responsive_log2fc_LS_R60_developing_seed_combined[match(rownames(wilcox_sig_stress_responsive_Anjali_LS_R60_developing_seed_combined), rownames(wilcox_sig_stress_responsive_log2fc_LS_R60_developing_seed_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_developing_seed_combined) = colnames(wilcox_sig_stress_responsive_log2fc_LS_R60_developing_seed_combined)[3]
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_developing_seed_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_developing_seed_combined$`GAR60 - GALS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_developing_seed_combined$Up.Down == "Up")    #1 metabolite
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_developing_seed_combined$Up.Down == "Down")  #14 metabolites


#venn - increased - 60h rewatering/late stress
vennlist_up_wilcox_stress_responsive_LS_R60_developing_seed_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_developing_seed_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R60_developing_seed_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_developing_seed_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_developing_seed_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_developing_seed_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_developing_seed_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_stress_responsive_LS_R60_developing_seed_combined) = c("N22", "Dular", "Anjali")
venn_up_stress_responsive_wilcox_LS_R60_developing_seed_combined = venn(vennlist_up_wilcox_stress_responsive_LS_R60_developing_seed_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_stress_responsive_R60_latestress_developing_seed_combined = venn.diagram(vennlist_up_wilcox_stress_responsive_LS_R60_developing_seed_combined, category.names = c(paste(names(vennlist_up_wilcox_stress_responsive_LS_R60_developing_seed_combined[1]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R60_developing_seed_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R60_developing_seed_combined[2]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R60_developing_seed_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R60_developing_seed_combined[3]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R60_developing_seed_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "E", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label E
intersect_stress_responsive_wilcox_LS_R60_venn_up_developing_seed_combined = attr(venn_up_stress_responsive_wilcox_LS_R60_developing_seed_combined, "intersections")  #list of unique and common metabolites
venn_up_stress_responsive_wilcox_metabolites_LS_R60_developing_seed_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R60_venn_up_developing_seed_combined, rbind))  #convert from list to matrix
colnames(venn_up_stress_responsive_wilcox_metabolites_LS_R60_developing_seed_combined) = venn_up_stress_responsive_wilcox_metabolites_LS_R60_developing_seed_combined[1, ]
venn_up_stress_responsive_wilcox_metabolites_LS_R60_developing_seed_combined = venn_up_stress_responsive_wilcox_metabolites_LS_R60_developing_seed_combined[-1, ]


#venn - decreased - 60h rewatering/late stress
vennlist_down_wilcox_stress_responsive_LS_R60_developing_seed_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_developing_seed_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R60_developing_seed_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_developing_seed_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_developing_seed_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_developing_seed_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_developing_seed_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_stress_responsive_LS_R60_developing_seed_combined) = c("N22", "Dular", "Anjali")
venn_down_stress_responsive_wilcox_LS_R60_developing_seed_combined = venn(vennlist_down_wilcox_stress_responsive_LS_R60_developing_seed_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_stress_responsive_R60_latestress_developing_seed_combined = venn.diagram(vennlist_down_wilcox_stress_responsive_LS_R60_developing_seed_combined, category.names = c(paste(names(vennlist_down_wilcox_stress_responsive_LS_R60_developing_seed_combined[1]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R60_developing_seed_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R60_developing_seed_combined[2]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R60_developing_seed_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R60_developing_seed_combined[3]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R60_developing_seed_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "F", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label F
intersect_stress_responsive_wilcox_LS_R60_venn_down_developing_seed_combined = attr(venn_down_stress_responsive_wilcox_LS_R60_developing_seed_combined, "intersections")  #list of unique and common metabolites
venn_down_stress_responsive_wilcox_metabolites_LS_R60_developing_seed_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R60_venn_down_developing_seed_combined, rbind))  #convert from list to matrix
colnames(venn_down_stress_responsive_wilcox_metabolites_LS_R60_developing_seed_combined) = venn_down_stress_responsive_wilcox_metabolites_LS_R60_developing_seed_combined[1, ]
venn_down_stress_responsive_wilcox_metabolites_LS_R60_developing_seed_combined = venn_down_stress_responsive_wilcox_metabolites_LS_R60_developing_seed_combined[-1, ]


write.table(venn_up_stress_responsive_wilcox_metabolites_LS_R60_developing_seed_combined, file = "wilcox_stress_responsive_sig_up_metabolites_developing_seed_R60_latestress_combined.txt", sep = "\t", quote = F)
write.table(venn_down_stress_responsive_wilcox_metabolites_LS_R60_developing_seed_combined, file = "wilcox_stress_responsive_sig_down_metabolites_developing_seed_R60_latestress_combined.txt", sep = "\t", quote = F)


####################


#rownames same - check before combining
#significance of wilcoxon test - all cultivars, all three RW timepoints
wilcox_sig_stress_responsive_LS_R12_R36_R60_developing_seed_combined = cbind(wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined,
                                                                             wilcox_res_stress_responsive_N22_LS_R36_developing_seed_combined$Significance,
                                                                             wilcox_res_stress_responsive_N22_LS_R60_developing_seed_combined$Significance,
                                                                             wilcox_res_stress_responsive_Dular_LS_R12_developing_seed_combined$Significance,
                                                                             wilcox_res_stress_responsive_Dular_LS_R36_developing_seed_combined$Significance,
                                                                             wilcox_res_stress_responsive_Dular_LS_R60_developing_seed_combined$Significance,
                                                                             wilcox_res_stress_responsive_Anjali_LS_R12_developing_seed_combined$Significance,
                                                                             wilcox_res_stress_responsive_Anjali_LS_R36_developing_seed_combined$Significance,
                                                                             wilcox_res_stress_responsive_Anjali_LS_R60_developing_seed_combined$Significance)
wilcox_sig_stress_responsive_LS_R12_R36_R60_developing_seed_combined = wilcox_sig_stress_responsive_LS_R12_R36_R60_developing_seed_combined[ , -1]
colnames(wilcox_sig_stress_responsive_LS_R12_R36_R60_developing_seed_combined)[1] = "wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined$Significance"
wilcox_sig_stress_responsive_LS_R12_R36_R60_developing_seed_combined = data.frame(lapply(wilcox_sig_stress_responsive_LS_R12_R36_R60_developing_seed_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_stress_responsive_LS_R12_R36_R60_developing_seed_combined) = rownames(wilcox_res_stress_responsive_N22_LS_R12_developing_seed_combined)
wilcox_sig_stress_responsive_LS_R12_R36_R60_developing_seed_combined$Nonsig = rowSums(wilcox_sig_stress_responsive_LS_R12_R36_R60_developing_seed_combined == " ")  #count number of nonsignificant
wilcox_sig_stress_responsive_final_LS_R12_R36_R60_developing_seed_combined = wilcox_sig_stress_responsive_LS_R12_R36_R60_developing_seed_combined[-which(wilcox_sig_stress_responsive_LS_R12_R36_R60_developing_seed_combined$Nonsig == "9"), ] #metabolites significant in at least one of the comparisons - 24 metabolites

#log2fc values of metabolites with significant difference in any one of the comparisons - 12, 36, 60h recovery vs. late stress in any cultivar
wilcox_sig_stress_responsive_log2fc_LS_R12_R36_R60_developing_seed_combined = log2fc_combined_latestress_recovery_developing_seed[which(rownames(log2fc_combined_latestress_recovery_developing_seed) %in% rownames(wilcox_sig_stress_responsive_final_LS_R12_R36_R60_developing_seed_combined)), ] 


#visualization
png("wilcox_sig_stress_responsive_log2fc_LS_RW_developing_seed_combined.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_stress_responsive_log2fc_LS_R12_R36_R60_developing_seed_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-3, -0.01, length.out = 300), seq(0.01, 3, length.out = 300)), 
          margins = c(20, 20),  
          lmat = rbind(c(0, 4, 5), c(0, 1, 0), c(3, 2, 0)),
          labCol = rep(c("12 h RW/Severe stress", "36 h RW/Severe stress", "60 h RW/Severe stress"), 3),
          cexCol = 1, 
          cexRow = 1.2, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 0.2, 11), 
          lwid = c(0.8, 5, 2.1),
          srtCol = 45, 
          cellnote = wilcox_sig_stress_responsive_final_LS_R12_R36_R60_developing_seed_combined, 
          notecol = "black", 
          notecex = 1.8, 
          sepcolor = "black", 
          colsep = c(seq(0, 9)), 
          rowsep = c(seq(0, 24)), 
          sepwidth = c(0.01, 0.01), 
          ColSideColors = c(rep("#66c2a5", 3), rep("#fc8d62", 3), rep("#8da0cb", 3)),
          hclustfun = function(x) hclust(x, method = "average"))
legend(x = -0.02, y = 1.1, xpd = TRUE, legend = "N22", bty = "n", cex = 1.2)
legend(x = 0.09, y = 1.1, xpd = TRUE, legend = "Dular", bty = "n", cex = 1.2)
legend(x = 0.21, y = 1.1, xpd = TRUE, legend = "Anjali", bty = "n", cex = 1.2)
legend(x = 0.83, y = 1.04, xpd = TRUE, legend = expression(paste("Log" ["2"], "-fold change")), bty = "n", cex = 0.8)
legend(x = -0.15, y = 1.13, xpd = TRUE, legend = c("G"), bty = "n", cex = 1.05, text.font = 2)
dev.off()


####################

#all venn diagrams

#multiplot
png("venn_up_down_developing_seed_stress_responsive_RW_latestress_combined.png", width = 16*300, height = 24*300, res = 300)
grid.arrange(gTree(children = venn_up_stress_responsive_R12_latestress_developing_seed_combined), 
             gTree(children = venn_down_stress_responsive_R12_latestress_developing_seed_combined), 
             gTree(children = venn_up_stress_responsive_R36_latestress_developing_seed_combined),
             gTree(children = venn_down_stress_responsive_R36_latestress_developing_seed_combined),
             gTree(children = venn_up_stress_responsive_R60_latestress_developing_seed_combined),
             gTree(children = venn_down_stress_responsive_R60_latestress_developing_seed_combined),
             ncol = 2, nrow = 3)
dev.off()


####################


#recovery vs. corresponding control - NOTE: consider all metabolites, not just the stress-responsive

#log2-fold change - control vs. rewatering
control_recovery_order = c("GNC12", "GNC36", "GNC60", "GNR12", "GNR36", "GNR60", "GDC12", "GDC36", "GDC60", "GDR12", "GDR36", "GDR60", "GAC12", "GAC36", "GAC60", "GAR12", "GAR36", "GAR60")
data_combined_control_recovery_log2_mean_developing_seed = data_combined_log2_mean_developing_seed_2[control_recovery_order, ]  #control and recovery data in preferred order
data_combined_control_recovery_log2_mean_developing_seed = t(data_combined_control_recovery_log2_mean_developing_seed)

log2fc_combined_N22_control_recovery_developing_seed = data_combined_control_recovery_log2_mean_developing_seed[ , 4:6] - data_combined_control_recovery_log2_mean_developing_seed[ , 1:3]
log2fc_combined_Dular_control_recovery_developing_seed = data_combined_control_recovery_log2_mean_developing_seed[ , 10:12] - data_combined_control_recovery_log2_mean_developing_seed[ , 7:9]
log2fc_combined_Anjali_control_recovery_developing_seed = data_combined_control_recovery_log2_mean_developing_seed[ , 16:18] - data_combined_control_recovery_log2_mean_developing_seed[ , 13:15]

#combine in one data set
log2fc_combined_control_recovery_developing_seed = cbind(log2fc_combined_N22_control_recovery_developing_seed, log2fc_combined_Dular_control_recovery_developing_seed, log2fc_combined_Anjali_control_recovery_developing_seed)   #combine log2-fc of all cultivars in one dataset
#rename columns
colnames(log2fc_combined_control_recovery_developing_seed) = c(paste(colnames(log2fc_combined_control_recovery_developing_seed), "-", sub("R", "C", colnames(log2fc_combined_control_recovery_developing_seed))))

#heatmap including all metabolites
heatmap.2(log2fc_combined_control_recovery_developing_seed, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none",
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(7, 20),  
          main = "Developing seed - Combined - Log2 FC, Recovery/Control", 
          cexCol = 1, 
          cexRow = 0.9, 
          density.info = "none",
          key.xlab = NA,
          key.title = NA, 
          lhei = c(1, 8),
          lwid = c(1.5, 5), srtCol = 90)


####################

#control vs. 12h rewatering

#N22

#data distribution
shapiro_res_combined_N22_control_recovery12_developing_seed = t(data.frame(lapply(data_combined_N22_log2median_developing_seed[which(data_combined_N22_log2median_developing_seed$Timepoint %in% c("Control.12h Rewatering", "12h Rewatering")), 3:ncol(data_combined_N22_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_control_recovery12_developing_seed) = "p.value"
shapiro_res_combined_N22_control_recovery12_developing_seed = as.data.frame(shapiro_res_combined_N22_control_recovery12_developing_seed)
shapiro_res_combined_N22_control_recovery12_developing_seed$Distribution = ifelse(shapiro_res_combined_N22_control_recovery12_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_control_recovery12_developing_seed$Distribution == "Normal")     #42 metabolites
sum(shapiro_res_combined_N22_control_recovery12_developing_seed$Distribution == "Non-normal") #25 metabolites

#wilcox test
wilcox_res_N22_C12_R12_developing_seed_combined = t(data.frame(lapply(data_combined_N22_log2median_developing_seed[ , 3:ncol(data_combined_N22_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_N22_log2median_developing_seed$Timepoint, subset = data_combined_N22_log2median_developing_seed$Timepoint %in% c("Control.12h Rewatering", "12h Rewatering"))$p.value)))
colnames(wilcox_res_N22_C12_R12_developing_seed_combined) = "p.value"
rownames(wilcox_res_N22_C12_R12_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_N22_C12_R12_developing_seed_combined = as.data.frame(wilcox_res_N22_C12_R12_developing_seed_combined)
wilcox_res_N22_C12_R12_developing_seed_combined$Significance = with(wilcox_res_N22_C12_R12_developing_seed_combined,
                                                                    ifelse(p.value <= 0.001, "***", ifelse(
                                                                      p.value <= 0.01, "**", ifelse(
                                                                        p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_C12_R12_developing_seed_combined$Significance == "*") #11 metabolites
sum(wilcox_res_N22_C12_R12_developing_seed_combined$Significance == "**") #5 metabolites
sum(wilcox_res_N22_C12_R12_developing_seed_combined$Significance == "***") #5 metabolites
sum(wilcox_res_N22_C12_R12_developing_seed_combined$Significance != "ns") #21 metabolites
sum(wilcox_res_N22_C12_R12_developing_seed_combined$Significance == "ns") #46 metabolites


#Dular

#data distribution
shapiro_res_combined_Dular_control_recovery12_developing_seed = t(data.frame(lapply(data_combined_Dular_log2median_developing_seed[which(data_combined_Dular_log2median_developing_seed$Timepoint %in% c("Control.12h Rewatering", "12h Rewatering")), 3:ncol(data_combined_Dular_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_control_recovery12_developing_seed) = "p.value"
shapiro_res_combined_Dular_control_recovery12_developing_seed = as.data.frame(shapiro_res_combined_Dular_control_recovery12_developing_seed)
shapiro_res_combined_Dular_control_recovery12_developing_seed$Distribution = ifelse(shapiro_res_combined_Dular_control_recovery12_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_control_recovery12_developing_seed$Distribution == "Normal")     #46 metabolites
sum(shapiro_res_combined_Dular_control_recovery12_developing_seed$Distribution == "Non-normal") #21 metabolites

#wilcox test
wilcox_res_Dular_C12_R12_developing_seed_combined = t(data.frame(lapply(data_combined_Dular_log2median_developing_seed[ , 3:ncol(data_combined_Dular_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_developing_seed$Timepoint, subset = data_combined_Dular_log2median_developing_seed$Timepoint %in% c("Control.12h Rewatering", "12h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_C12_R12_developing_seed_combined) = "p.value"
rownames(wilcox_res_Dular_C12_R12_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_Dular_C12_R12_developing_seed_combined = as.data.frame(wilcox_res_Dular_C12_R12_developing_seed_combined)
wilcox_res_Dular_C12_R12_developing_seed_combined$Significance = with(wilcox_res_Dular_C12_R12_developing_seed_combined,
                                                                      ifelse(p.value <= 0.001, "***", ifelse(
                                                                        p.value <= 0.01, "**", ifelse(
                                                                          p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_C12_R12_developing_seed_combined$Significance == "*") #6 metabolites
sum(wilcox_res_Dular_C12_R12_developing_seed_combined$Significance == "**") #13 metabolites
sum(wilcox_res_Dular_C12_R12_developing_seed_combined$Significance == "***") #11 metabolites
sum(wilcox_res_Dular_C12_R12_developing_seed_combined$Significance != "ns") #30 metabolites
sum(wilcox_res_Dular_C12_R12_developing_seed_combined$Significance == "ns") #37 metabolites


#Anjali

#data distribution
shapiro_res_combined_Anjali_control_recovery12_developing_seed = t(data.frame(lapply(data_combined_Anjali_log2median_developing_seed[which(data_combined_Anjali_log2median_developing_seed$Timepoint %in% c("Control.12h Rewatering", "12h Rewatering")), 3:ncol(data_combined_Anjali_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_control_recovery12_developing_seed) = "p.value"
shapiro_res_combined_Anjali_control_recovery12_developing_seed = as.data.frame(shapiro_res_combined_Anjali_control_recovery12_developing_seed)
shapiro_res_combined_Anjali_control_recovery12_developing_seed$Distribution = ifelse(shapiro_res_combined_Anjali_control_recovery12_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_control_recovery12_developing_seed$Distribution == "Normal")     #47 metabolites
sum(shapiro_res_combined_Anjali_control_recovery12_developing_seed$Distribution == "Non-normal") #20 metabolites

#wilcox test
wilcox_res_Anjali_C12_R12_developing_seed_combined = t(data.frame(lapply(data_combined_Anjali_log2median_developing_seed[ , 3:ncol(data_combined_Anjali_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_developing_seed$Timepoint, subset = data_combined_Anjali_log2median_developing_seed$Timepoint %in% c("Control.12h Rewatering", "12h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_C12_R12_developing_seed_combined) = "p.value"
rownames(wilcox_res_Anjali_C12_R12_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_Anjali_C12_R12_developing_seed_combined = as.data.frame(wilcox_res_Anjali_C12_R12_developing_seed_combined)
wilcox_res_Anjali_C12_R12_developing_seed_combined$Significance = with(wilcox_res_Anjali_C12_R12_developing_seed_combined,
                                                                       ifelse(p.value <= 0.001, "***", ifelse(
                                                                         p.value <= 0.01, "**", ifelse(
                                                                           p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_C12_R12_developing_seed_combined$Significance == "*") #7 metabolites
sum(wilcox_res_Anjali_C12_R12_developing_seed_combined$Significance == "**") #5 metabolites
sum(wilcox_res_Anjali_C12_R12_developing_seed_combined$Significance == "***") #11 metabolites
sum(wilcox_res_Anjali_C12_R12_developing_seed_combined$Significance != "ns") #23 metabolites
sum(wilcox_res_Anjali_C12_R12_developing_seed_combined$Significance == "ns") #44 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_C12_R12_developing_seed_combined = cbind(wilcox_res_N22_C12_R12_developing_seed_combined, 
                                                    wilcox_res_Dular_C12_R12_developing_seed_combined$Significance, 
                                                    wilcox_res_Anjali_C12_R12_developing_seed_combined$Significance)
wilcox_sig_C12_R12_developing_seed_combined = wilcox_sig_C12_R12_developing_seed_combined[ , -1]
colnames(wilcox_sig_C12_R12_developing_seed_combined)[1] = "wilcox_res_N22_C12_R12_developing_seed_combined$Significance"
wilcox_sig_C12_R12_developing_seed_combined = data.frame(lapply(wilcox_sig_C12_R12_developing_seed_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_C12_R12_developing_seed_combined) = rownames(wilcox_res_N22_C12_R12_developing_seed_combined)
wilcox_sig_C12_R12_developing_seed_combined$Nonsig = rowSums(wilcox_sig_C12_R12_developing_seed_combined == " ")  #count number of nonsignificant
wilcox_sig_final_C12_R12_developing_seed_combined = wilcox_sig_C12_R12_developing_seed_combined[-which(wilcox_sig_C12_R12_developing_seed_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 42 metabolites

wilcox_sig_log2fc_C12_R12_developing_seed_combined = log2fc_combined_control_recovery_developing_seed[which(rownames(log2fc_combined_control_recovery_developing_seed) %in% rownames(wilcox_sig_final_C12_R12_developing_seed_combined)), c(1, 4, 7)]  #log2FC values of metabolites significant in at least one of the cultivars

#overview - with HCA - Euclidean distance and average linkage 
heatmap.2(wilcox_sig_log2fc_C12_R12_developing_seed_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(5, 20), 
          main = "Developing seed-Combined-Wilcoxon \nLog2 FC, EGF, 12h RW/Control 12h", 
          cexCol = 1.5, 
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          labCol = c("N22", "Dular", "Anjali"), 
          srtCol = 0, 
          adjCol = c(0.5, 0.5),
          cellnote = wilcox_sig_final_C12_R12_developing_seed_combined, 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black",
          colsep = c(seq(0, 3)), 
          rowsep = c(seq(0, 42)), 
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))


#venn diagram - common among cultivars

#N22
wilcox_sig_N22_C12_R12_developing_seed_combined = wilcox_res_N22_C12_R12_developing_seed_combined[-which(wilcox_res_N22_C12_R12_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_N22_C12_R12_developing_seed_combined = as.data.frame(wilcox_sig_log2fc_C12_R12_developing_seed_combined[match(rownames(wilcox_sig_N22_C12_R12_developing_seed_combined), rownames(wilcox_sig_log2fc_C12_R12_developing_seed_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_N22_C12_R12_developing_seed_combined) = colnames(wilcox_sig_log2fc_C12_R12_developing_seed_combined)[1]
wilcox_sig_log2fc_N22_C12_R12_developing_seed_combined$Up.Down = ifelse(wilcox_sig_log2fc_N22_C12_R12_developing_seed_combined$`GNR12 - GNC12` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_N22_C12_R12_developing_seed_combined$Up.Down == "Up")    #13 metabolites
sum(wilcox_sig_log2fc_N22_C12_R12_developing_seed_combined$Up.Down == "Down")  #8 metabolites

#Dular
wilcox_sig_Dular_C12_R12_developing_seed_combined = wilcox_res_Dular_C12_R12_developing_seed_combined[-which(wilcox_res_Dular_C12_R12_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Dular_C12_R12_developing_seed_combined = as.data.frame(wilcox_sig_log2fc_C12_R12_developing_seed_combined[match(rownames(wilcox_sig_Dular_C12_R12_developing_seed_combined), rownames(wilcox_sig_log2fc_C12_R12_developing_seed_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Dular_C12_R12_developing_seed_combined) = colnames(wilcox_sig_log2fc_C12_R12_developing_seed_combined)[2]
wilcox_sig_log2fc_Dular_C12_R12_developing_seed_combined$Up.Down = ifelse(wilcox_sig_log2fc_Dular_C12_R12_developing_seed_combined$`GDR12 - GDC12` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Dular_C12_R12_developing_seed_combined$Up.Down == "Up")    #24 metabolites
sum(wilcox_sig_log2fc_Dular_C12_R12_developing_seed_combined$Up.Down == "Down")  #6 metabolites

#Anjali
wilcox_sig_Anjali_C12_R12_developing_seed_combined = wilcox_res_Anjali_C12_R12_developing_seed_combined[-which(wilcox_res_Anjali_C12_R12_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Anjali_C12_R12_developing_seed_combined = as.data.frame(wilcox_sig_log2fc_C12_R12_developing_seed_combined[match(rownames(wilcox_sig_Anjali_C12_R12_developing_seed_combined), rownames(wilcox_sig_log2fc_C12_R12_developing_seed_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Anjali_C12_R12_developing_seed_combined) = colnames(wilcox_sig_log2fc_C12_R12_developing_seed_combined)[3]
wilcox_sig_log2fc_Anjali_C12_R12_developing_seed_combined$Up.Down = ifelse(wilcox_sig_log2fc_Anjali_C12_R12_developing_seed_combined$`GAR12 - GAC12` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Anjali_C12_R12_developing_seed_combined$Up.Down == "Up")    #10 metabolites
sum(wilcox_sig_log2fc_Anjali_C12_R12_developing_seed_combined$Up.Down == "Down")  #13 metabolites


#venn - increased - 12h rewatering/control 12h
vennlist_up_wilcox_C12_R12_developing_seed_combined = list(rownames(wilcox_sig_log2fc_N22_C12_R12_developing_seed_combined)[wilcox_sig_log2fc_N22_C12_R12_developing_seed_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Dular_C12_R12_developing_seed_combined)[wilcox_sig_log2fc_Dular_C12_R12_developing_seed_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Anjali_C12_R12_developing_seed_combined)[wilcox_sig_log2fc_Anjali_C12_R12_developing_seed_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_C12_R12_developing_seed_combined) = c("N22", "Dular", "Anjali")
venn_up_wilcox_C12_R12_developing_seed_combined = venn(vennlist_up_wilcox_C12_R12_developing_seed_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_R12_C12_developing_seed_combined = venn.diagram(vennlist_up_wilcox_C12_R12_developing_seed_combined, category.names = c(paste(names(vennlist_up_wilcox_C12_R12_developing_seed_combined[1]), " (", length(vennlist_up_wilcox_C12_R12_developing_seed_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_C12_R12_developing_seed_combined[2]), " (", length(vennlist_up_wilcox_C12_R12_developing_seed_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_C12_R12_developing_seed_combined[3]), " (", length(vennlist_up_wilcox_C12_R12_developing_seed_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "A", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label A
intersect_wilcox_C12_R12_venn_up_developing_seed_combined = attr(venn_up_wilcox_C12_R12_developing_seed_combined, "intersections")  #list of unique and common metabolites
venn_up_wilcox_metabolites_C12_R12_developing_seed_combined = t(ldply(intersect_wilcox_C12_R12_venn_up_developing_seed_combined, rbind))  #convert from list to matrix
colnames(venn_up_wilcox_metabolites_C12_R12_developing_seed_combined) = venn_up_wilcox_metabolites_C12_R12_developing_seed_combined[1, ]
venn_up_wilcox_metabolites_C12_R12_developing_seed_combined = venn_up_wilcox_metabolites_C12_R12_developing_seed_combined[-1, ]


#venn - decreased - 12h rewatering/control 12h
vennlist_down_wilcox_C12_R12_developing_seed_combined = list(rownames(wilcox_sig_log2fc_N22_C12_R12_developing_seed_combined)[wilcox_sig_log2fc_N22_C12_R12_developing_seed_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Dular_C12_R12_developing_seed_combined)[wilcox_sig_log2fc_Dular_C12_R12_developing_seed_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Anjali_C12_R12_developing_seed_combined)[wilcox_sig_log2fc_Anjali_C12_R12_developing_seed_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_C12_R12_developing_seed_combined) = c("N22", "Dular", "Anjali")
venn_down_wilcox_C12_R12_developing_seed_combined = venn(vennlist_down_wilcox_C12_R12_developing_seed_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_R12_C12_developing_seed_combined = venn.diagram(vennlist_down_wilcox_C12_R12_developing_seed_combined, category.names = c(paste(names(vennlist_down_wilcox_C12_R12_developing_seed_combined[1]), " (", length(vennlist_down_wilcox_C12_R12_developing_seed_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_C12_R12_developing_seed_combined[2]), " (", length(vennlist_down_wilcox_C12_R12_developing_seed_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_C12_R12_developing_seed_combined[3]), " (", length(vennlist_down_wilcox_C12_R12_developing_seed_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "B", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label B
intersect_wilcox_C12_R12_venn_down_developing_seed_combined = attr(venn_down_wilcox_C12_R12_developing_seed_combined, "intersections")  #list of unique and common metabolites
venn_down_wilcox_metabolites_C12_R12_developing_seed_combined = t(ldply(intersect_wilcox_C12_R12_venn_down_developing_seed_combined, rbind))  #convert from list to matrix
colnames(venn_down_wilcox_metabolites_C12_R12_developing_seed_combined) = venn_down_wilcox_metabolites_C12_R12_developing_seed_combined[1, ]
venn_down_wilcox_metabolites_C12_R12_developing_seed_combined = venn_down_wilcox_metabolites_C12_R12_developing_seed_combined[-1, ]


write.table(venn_up_wilcox_metabolites_C12_R12_developing_seed_combined, file = "wilcox_sig_up_metabolites_developing_seed_R12_C12_combined.txt", sep = "\t", quote = F)
write.table(venn_down_wilcox_metabolites_C12_R12_developing_seed_combined, file = "wilcox_sig_down_metabolites_developing_seed_R12_C12_combined.txt", sep = "\t", quote = F)


####################


#control 36h vs. 36h rewatering

#N22

#data distribution
shapiro_res_combined_N22_control_recovery36_developing_seed = t(data.frame(lapply(data_combined_N22_log2median_developing_seed[which(data_combined_N22_log2median_developing_seed$Timepoint %in% c("Control.36h Rewatering", "36h Rewatering")), 3:ncol(data_combined_N22_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_control_recovery36_developing_seed) = "p.value"
shapiro_res_combined_N22_control_recovery36_developing_seed = as.data.frame(shapiro_res_combined_N22_control_recovery36_developing_seed)
shapiro_res_combined_N22_control_recovery36_developing_seed$Distribution = ifelse(shapiro_res_combined_N22_control_recovery36_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_control_recovery36_developing_seed$Distribution == "Normal")     #42 metabolites
sum(shapiro_res_combined_N22_control_recovery36_developing_seed$Distribution == "Non-normal") #25 metabolites

#wilcox test
wilcox_res_N22_C36_R36_developing_seed_combined = t(data.frame(lapply(data_combined_N22_log2median_developing_seed[ , 3:ncol(data_combined_N22_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_N22_log2median_developing_seed$Timepoint, subset = data_combined_N22_log2median_developing_seed$Timepoint %in% c("Control.36h Rewatering", "36h Rewatering"))$p.value)))
colnames(wilcox_res_N22_C36_R36_developing_seed_combined) = "p.value"
rownames(wilcox_res_N22_C36_R36_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_N22_C36_R36_developing_seed_combined = as.data.frame(wilcox_res_N22_C36_R36_developing_seed_combined)
wilcox_res_N22_C36_R36_developing_seed_combined$Significance = with(wilcox_res_N22_C36_R36_developing_seed_combined,
                                                                    ifelse(p.value <= 0.001, "***", ifelse(
                                                                      p.value <= 0.01, "**", ifelse(
                                                                        p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_C36_R36_developing_seed_combined$Significance == "*") #5 metabolites
sum(wilcox_res_N22_C36_R36_developing_seed_combined$Significance == "**") #5 metabolites
sum(wilcox_res_N22_C36_R36_developing_seed_combined$Significance == "***") #7 metabolites
sum(wilcox_res_N22_C36_R36_developing_seed_combined$Significance != "ns") #17 metabolites
sum(wilcox_res_N22_C36_R36_developing_seed_combined$Significance == "ns") #50 metabolites


#Dular

#data distribution
shapiro_res_combined_Dular_control_recovery36_developing_seed = t(data.frame(lapply(data_combined_Dular_log2median_developing_seed[which(data_combined_Dular_log2median_developing_seed$Timepoint %in% c("Control.36h Rewatering", "36h Rewatering")), 3:ncol(data_combined_Dular_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_control_recovery36_developing_seed) = "p.value"
shapiro_res_combined_Dular_control_recovery36_developing_seed = as.data.frame(shapiro_res_combined_Dular_control_recovery36_developing_seed)
shapiro_res_combined_Dular_control_recovery36_developing_seed$Distribution = ifelse(shapiro_res_combined_Dular_control_recovery36_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_control_recovery36_developing_seed$Distribution == "Normal")     #55 metabolites
sum(shapiro_res_combined_Dular_control_recovery36_developing_seed$Distribution == "Non-normal") #12 metabolites

#wilcox test
wilcox_res_Dular_C36_R36_developing_seed_combined = t(data.frame(lapply(data_combined_Dular_log2median_developing_seed[ , 3:ncol(data_combined_Dular_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_developing_seed$Timepoint, subset = data_combined_Dular_log2median_developing_seed$Timepoint %in% c("Control.36h Rewatering", "36h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_C36_R36_developing_seed_combined) = "p.value"
rownames(wilcox_res_Dular_C36_R36_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_Dular_C36_R36_developing_seed_combined = as.data.frame(wilcox_res_Dular_C36_R36_developing_seed_combined)
wilcox_res_Dular_C36_R36_developing_seed_combined$Significance = with(wilcox_res_Dular_C36_R36_developing_seed_combined,
                                                                      ifelse(p.value <= 0.001, "***", ifelse(
                                                                        p.value <= 0.01, "**", ifelse(
                                                                          p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_C36_R36_developing_seed_combined$Significance == "*") #9 metabolites
sum(wilcox_res_Dular_C36_R36_developing_seed_combined$Significance == "**") #5 metabolites
sum(wilcox_res_Dular_C36_R36_developing_seed_combined$Significance == "***") #7 metabolites
sum(wilcox_res_Dular_C36_R36_developing_seed_combined$Significance != "ns") #21 metabolites
sum(wilcox_res_Dular_C36_R36_developing_seed_combined$Significance == "ns") #46 metabolites


#Anjali

#data distribution
shapiro_res_combined_Anjali_control_recovery36_developing_seed = t(data.frame(lapply(data_combined_Anjali_log2median_developing_seed[which(data_combined_Anjali_log2median_developing_seed$Timepoint %in% c("Control.36h Rewatering", "36h Rewatering")), 3:ncol(data_combined_Anjali_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_control_recovery36_developing_seed) = "p.value"
shapiro_res_combined_Anjali_control_recovery36_developing_seed = as.data.frame(shapiro_res_combined_Anjali_control_recovery36_developing_seed)
shapiro_res_combined_Anjali_control_recovery36_developing_seed$Distribution = ifelse(shapiro_res_combined_Anjali_control_recovery36_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_control_recovery36_developing_seed$Distribution == "Normal")     #46 metabolites
sum(shapiro_res_combined_Anjali_control_recovery36_developing_seed$Distribution == "Non-normal") #21 metabolites

#wilcox test
wilcox_res_Anjali_C36_R36_developing_seed_combined = t(data.frame(lapply(data_combined_Anjali_log2median_developing_seed[ , 3:ncol(data_combined_Anjali_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_developing_seed$Timepoint, subset = data_combined_Anjali_log2median_developing_seed$Timepoint %in% c("Control.36h Rewatering", "36h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_C36_R36_developing_seed_combined) = "p.value"
rownames(wilcox_res_Anjali_C36_R36_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_Anjali_C36_R36_developing_seed_combined = as.data.frame(wilcox_res_Anjali_C36_R36_developing_seed_combined)
wilcox_res_Anjali_C36_R36_developing_seed_combined$Significance = with(wilcox_res_Anjali_C36_R36_developing_seed_combined,
                                                                       ifelse(p.value <= 0.001, "***", ifelse(
                                                                         p.value <= 0.01, "**", ifelse(
                                                                           p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_C36_R36_developing_seed_combined$Significance == "*") #10 metabolites
sum(wilcox_res_Anjali_C36_R36_developing_seed_combined$Significance == "**") #3 metabolites
sum(wilcox_res_Anjali_C36_R36_developing_seed_combined$Significance == "***") #7 metabolites
sum(wilcox_res_Anjali_C36_R36_developing_seed_combined$Significance != "ns") #20 metabolites
sum(wilcox_res_Anjali_C36_R36_developing_seed_combined$Significance == "ns") #47 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_C36_R36_developing_seed_combined = cbind(wilcox_res_N22_C36_R36_developing_seed_combined, 
                                                    wilcox_res_Dular_C36_R36_developing_seed_combined$Significance, 
                                                    wilcox_res_Anjali_C36_R36_developing_seed_combined$Significance)
wilcox_sig_C36_R36_developing_seed_combined = wilcox_sig_C36_R36_developing_seed_combined[ , -1]
colnames(wilcox_sig_C36_R36_developing_seed_combined)[1] = "wilcox_res_N22_C36_R36_developing_seed_combined$Significance"
wilcox_sig_C36_R36_developing_seed_combined = data.frame(lapply(wilcox_sig_C36_R36_developing_seed_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_C36_R36_developing_seed_combined) = rownames(wilcox_res_N22_C36_R36_developing_seed_combined)
wilcox_sig_C36_R36_developing_seed_combined$Nonsig = rowSums(wilcox_sig_C36_R36_developing_seed_combined == " ")  #count number of nonsignificant
wilcox_sig_final_C36_R36_developing_seed_combined = wilcox_sig_C36_R36_developing_seed_combined[-which(wilcox_sig_C36_R36_developing_seed_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 35 metabolites

wilcox_sig_log2fc_C36_R36_developing_seed_combined = log2fc_combined_control_recovery_developing_seed[which(rownames(log2fc_combined_control_recovery_developing_seed) %in% rownames(wilcox_sig_final_C36_R36_developing_seed_combined)), c(2, 5, 8)]  #log2FC values of metabolites significant in at least one of the cultivars

#overview - with HCA - Euclidean distance and average linkage 
heatmap.2(wilcox_sig_log2fc_C36_R36_developing_seed_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row",
          trace = "none", 
          col = bluered(599),
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(5, 20), 
          main = "Developing seed-Combined-Wilcoxon \nLog2 FC, EGF, 36h RW/Control 36h", 
          cexCol = 1.5, 
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          labCol = c("N22", "Dular", "Anjali"), 
          srtCol = 0, 
          adjCol = c(0.5, 0.5), 
          cellnote = wilcox_sig_final_C36_R36_developing_seed_combined, 
          notecol = "black",
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 3)), 
          rowsep = c(seq(0, 35)), 
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))


#venn diagram - common among cultivars

#N22
wilcox_sig_N22_C36_R36_developing_seed_combined = wilcox_res_N22_C36_R36_developing_seed_combined[-which(wilcox_res_N22_C36_R36_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_N22_C36_R36_developing_seed_combined = as.data.frame(wilcox_sig_log2fc_C36_R36_developing_seed_combined[match(rownames(wilcox_sig_N22_C36_R36_developing_seed_combined), rownames(wilcox_sig_log2fc_C36_R36_developing_seed_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_N22_C36_R36_developing_seed_combined) = colnames(wilcox_sig_log2fc_C36_R36_developing_seed_combined)[1]
wilcox_sig_log2fc_N22_C36_R36_developing_seed_combined$Up.Down = ifelse(wilcox_sig_log2fc_N22_C36_R36_developing_seed_combined$`GNR36 - GNC36` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_N22_C36_R36_developing_seed_combined$Up.Down == "Up")    #5 metabolites
sum(wilcox_sig_log2fc_N22_C36_R36_developing_seed_combined$Up.Down == "Down")  #12 metabolites

#Dular
wilcox_sig_Dular_C36_R36_developing_seed_combined = wilcox_res_Dular_C36_R36_developing_seed_combined[-which(wilcox_res_Dular_C36_R36_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Dular_C36_R36_developing_seed_combined = as.data.frame(wilcox_sig_log2fc_C36_R36_developing_seed_combined[match(rownames(wilcox_sig_Dular_C36_R36_developing_seed_combined), rownames(wilcox_sig_log2fc_C36_R36_developing_seed_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Dular_C36_R36_developing_seed_combined) = colnames(wilcox_sig_log2fc_C36_R36_developing_seed_combined)[2]
wilcox_sig_log2fc_Dular_C36_R36_developing_seed_combined$Up.Down = ifelse(wilcox_sig_log2fc_Dular_C36_R36_developing_seed_combined$`GDR36 - GDC36` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Dular_C36_R36_developing_seed_combined$Up.Down == "Up")    #17 metabolites
sum(wilcox_sig_log2fc_Dular_C36_R36_developing_seed_combined$Up.Down == "Down")  #4 metabolites

#Anjali
wilcox_sig_Anjali_C36_R36_developing_seed_combined = wilcox_res_Anjali_C36_R36_developing_seed_combined[-which(wilcox_res_Anjali_C36_R36_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Anjali_C36_R36_developing_seed_combined = as.data.frame(wilcox_sig_log2fc_C36_R36_developing_seed_combined[match(rownames(wilcox_sig_Anjali_C36_R36_developing_seed_combined), rownames(wilcox_sig_log2fc_C36_R36_developing_seed_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Anjali_C36_R36_developing_seed_combined) = colnames(wilcox_sig_log2fc_C36_R36_developing_seed_combined)[3]
wilcox_sig_log2fc_Anjali_C36_R36_developing_seed_combined$Up.Down = ifelse(wilcox_sig_log2fc_Anjali_C36_R36_developing_seed_combined$`GAR36 - GAC36` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Anjali_C36_R36_developing_seed_combined$Up.Down == "Up")    #6 metabolites
sum(wilcox_sig_log2fc_Anjali_C36_R36_developing_seed_combined$Up.Down == "Down")  #14 metabolites


#venn - increased - 36h rewatering/control 36h
vennlist_up_wilcox_C36_R36_developing_seed_combined = list(rownames(wilcox_sig_log2fc_N22_C36_R36_developing_seed_combined)[wilcox_sig_log2fc_N22_C36_R36_developing_seed_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Dular_C36_R36_developing_seed_combined)[wilcox_sig_log2fc_Dular_C36_R36_developing_seed_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Anjali_C36_R36_developing_seed_combined)[wilcox_sig_log2fc_Anjali_C36_R36_developing_seed_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_C36_R36_developing_seed_combined) = c("N22", "Dular", "Anjali")
venn_up_wilcox_C36_R36_developing_seed_combined = venn(vennlist_up_wilcox_C36_R36_developing_seed_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_R36_C36_developing_seed_combined = venn.diagram(vennlist_up_wilcox_C36_R36_developing_seed_combined, category.names = c(paste(names(vennlist_up_wilcox_C36_R36_developing_seed_combined[1]), " (", length(vennlist_up_wilcox_C36_R36_developing_seed_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_C36_R36_developing_seed_combined[2]), " (", length(vennlist_up_wilcox_C36_R36_developing_seed_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_C36_R36_developing_seed_combined[3]), " (", length(vennlist_up_wilcox_C36_R36_developing_seed_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "C", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label C
intersect_wilcox_C36_R36_venn_up_developing_seed_combined = attr(venn_up_wilcox_C36_R36_developing_seed_combined, "intersections")  #list of unique and common metabolites
venn_up_wilcox_metabolites_C36_R36_developing_seed_combined = t(ldply(intersect_wilcox_C36_R36_venn_up_developing_seed_combined, rbind))  #convert from list to matrix
colnames(venn_up_wilcox_metabolites_C36_R36_developing_seed_combined) = venn_up_wilcox_metabolites_C36_R36_developing_seed_combined[1, ]
venn_up_wilcox_metabolites_C36_R36_developing_seed_combined = venn_up_wilcox_metabolites_C36_R36_developing_seed_combined[-1, ]


#venn - decreased - 36h rewatering/control 36h
vennlist_down_wilcox_C36_R36_developing_seed_combined = list(rownames(wilcox_sig_log2fc_N22_C36_R36_developing_seed_combined)[wilcox_sig_log2fc_N22_C36_R36_developing_seed_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Dular_C36_R36_developing_seed_combined)[wilcox_sig_log2fc_Dular_C36_R36_developing_seed_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Anjali_C36_R36_developing_seed_combined)[wilcox_sig_log2fc_Anjali_C36_R36_developing_seed_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_C36_R36_developing_seed_combined) = c("N22", "Dular", "Anjali")
venn_down_wilcox_C36_R36_developing_seed_combined = venn(vennlist_down_wilcox_C36_R36_developing_seed_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_R36_C36_developing_seed_combined = venn.diagram(vennlist_down_wilcox_C36_R36_developing_seed_combined, category.names = c(paste(names(vennlist_down_wilcox_C36_R36_developing_seed_combined[1]), " (", length(vennlist_down_wilcox_C36_R36_developing_seed_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_C36_R36_developing_seed_combined[2]), " (", length(vennlist_down_wilcox_C36_R36_developing_seed_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_C36_R36_developing_seed_combined[3]), " (", length(vennlist_down_wilcox_C36_R36_developing_seed_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "D", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label D
intersect_wilcox_C36_R36_venn_down_developing_seed_combined = attr(venn_down_wilcox_C36_R36_developing_seed_combined, "intersections")  #list of unique and common metabolites
venn_down_wilcox_metabolites_C36_R36_developing_seed_combined = t(ldply(intersect_wilcox_C36_R36_venn_down_developing_seed_combined, rbind))  #convert from list to matrix
colnames(venn_down_wilcox_metabolites_C36_R36_developing_seed_combined) = venn_down_wilcox_metabolites_C36_R36_developing_seed_combined[1, ]
venn_down_wilcox_metabolites_C36_R36_developing_seed_combined = venn_down_wilcox_metabolites_C36_R36_developing_seed_combined[-1, ]


write.table(venn_up_wilcox_metabolites_C36_R36_developing_seed_combined, file = "wilcox_sig_up_metabolites_developing_seed_R36_C36_combined.txt", sep = "\t", quote = F)
write.table(venn_down_wilcox_metabolites_C36_R36_developing_seed_combined, file = "wilcox_sig_down_metabolites_developing_seed_R36_C36_combined.txt", sep = "\t", quote = F)


####################


#control 60h vs. 60h rewatering

#N22

#data distribution
shapiro_res_combined_N22_control_recovery60_developing_seed = t(data.frame(lapply(data_combined_N22_log2median_developing_seed[which(data_combined_N22_log2median_developing_seed$Timepoint %in% c("Control.60h Rewatering", "60h Rewatering")), 3:ncol(data_combined_N22_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_control_recovery60_developing_seed) = "p.value"
shapiro_res_combined_N22_control_recovery60_developing_seed = as.data.frame(shapiro_res_combined_N22_control_recovery60_developing_seed)
shapiro_res_combined_N22_control_recovery60_developing_seed$Distribution = ifelse(shapiro_res_combined_N22_control_recovery60_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_control_recovery60_developing_seed$Distribution == "Normal")     #37 metabolites
sum(shapiro_res_combined_N22_control_recovery60_developing_seed$Distribution == "Non-normal") #30 metabolites

#wilcox test
wilcox_res_N22_C60_R60_developing_seed_combined = t(data.frame(lapply(data_combined_N22_log2median_developing_seed[ , 3:ncol(data_combined_N22_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_N22_log2median_developing_seed$Timepoint, subset = data_combined_N22_log2median_developing_seed$Timepoint %in% c("Control.60h Rewatering", "60h Rewatering"))$p.value)))
colnames(wilcox_res_N22_C60_R60_developing_seed_combined) = "p.value"
rownames(wilcox_res_N22_C60_R60_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_N22_C60_R60_developing_seed_combined = as.data.frame(wilcox_res_N22_C60_R60_developing_seed_combined)
wilcox_res_N22_C60_R60_developing_seed_combined$Significance = with(wilcox_res_N22_C60_R60_developing_seed_combined,
                                                                    ifelse(p.value <= 0.001, "***", ifelse(
                                                                      p.value <= 0.01, "**", ifelse(
                                                                        p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_C60_R60_developing_seed_combined$Significance == "*") #5 metabolites
sum(wilcox_res_N22_C60_R60_developing_seed_combined$Significance == "**") #5 metabolites
sum(wilcox_res_N22_C60_R60_developing_seed_combined$Significance == "***") #3 metabolites
sum(wilcox_res_N22_C60_R60_developing_seed_combined$Significance != "ns") #13 metabolites
sum(wilcox_res_N22_C60_R60_developing_seed_combined$Significance == "ns") #54 metabolites


#Dular

#data distribution
shapiro_res_combined_Dular_control_recovery60_developing_seed = t(data.frame(lapply(data_combined_Dular_log2median_developing_seed[which(data_combined_Dular_log2median_developing_seed$Timepoint %in% c("Control.60h Rewatering", "60h Rewatering")), 3:ncol(data_combined_Dular_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_control_recovery60_developing_seed) = "p.value"
shapiro_res_combined_Dular_control_recovery60_developing_seed = as.data.frame(shapiro_res_combined_Dular_control_recovery60_developing_seed)
shapiro_res_combined_Dular_control_recovery60_developing_seed$Distribution = ifelse(shapiro_res_combined_Dular_control_recovery60_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_control_recovery60_developing_seed$Distribution == "Normal")     #44 metabolites
sum(shapiro_res_combined_Dular_control_recovery60_developing_seed$Distribution == "Non-normal") #23 metabolites

#wilcox test
wilcox_res_Dular_C60_R60_developing_seed_combined = t(data.frame(lapply(data_combined_Dular_log2median_developing_seed[ , 3:ncol(data_combined_Dular_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_developing_seed$Timepoint, subset = data_combined_Dular_log2median_developing_seed$Timepoint %in% c("Control.60h Rewatering", "60h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_C60_R60_developing_seed_combined) = "p.value"
rownames(wilcox_res_Dular_C60_R60_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_Dular_C60_R60_developing_seed_combined = as.data.frame(wilcox_res_Dular_C60_R60_developing_seed_combined)
wilcox_res_Dular_C60_R60_developing_seed_combined$Significance = with(wilcox_res_Dular_C60_R60_developing_seed_combined,
                                                                      ifelse(p.value <= 0.001, "***", ifelse(
                                                                        p.value <= 0.01, "**", ifelse(
                                                                          p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_C60_R60_developing_seed_combined$Significance == "*") #7 metabolites
sum(wilcox_res_Dular_C60_R60_developing_seed_combined$Significance == "**") #2 metabolites
sum(wilcox_res_Dular_C60_R60_developing_seed_combined$Significance == "***") #3 metabolites
sum(wilcox_res_Dular_C60_R60_developing_seed_combined$Significance != "ns") #12 metabolites
sum(wilcox_res_Dular_C60_R60_developing_seed_combined$Significance == "ns") #55 metabolites


#Anjali

#data distribution
shapiro_res_combined_Anjali_control_recovery60_developing_seed = t(data.frame(lapply(data_combined_Anjali_log2median_developing_seed[which(data_combined_Anjali_log2median_developing_seed$Timepoint %in% c("Control.60h Rewatering", "60h Rewatering")), 3:ncol(data_combined_Anjali_log2median_developing_seed)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_control_recovery60_developing_seed) = "p.value"
shapiro_res_combined_Anjali_control_recovery60_developing_seed = as.data.frame(shapiro_res_combined_Anjali_control_recovery60_developing_seed)
shapiro_res_combined_Anjali_control_recovery60_developing_seed$Distribution = ifelse(shapiro_res_combined_Anjali_control_recovery60_developing_seed$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_control_recovery60_developing_seed$Distribution == "Normal")     #49 metabolites
sum(shapiro_res_combined_Anjali_control_recovery60_developing_seed$Distribution == "Non-normal") #18 metabolites

#wilcox test
wilcox_res_Anjali_C60_R60_developing_seed_combined = t(data.frame(lapply(data_combined_Anjali_log2median_developing_seed[ , 3:ncol(data_combined_Anjali_log2median_developing_seed)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_developing_seed$Timepoint, subset = data_combined_Anjali_log2median_developing_seed$Timepoint %in% c("Control.60h Rewatering", "60h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_C60_R60_developing_seed_combined) = "p.value"
rownames(wilcox_res_Anjali_C60_R60_developing_seed_combined) = reduced_metabolite_list_final_developing_seed_2013$Name
wilcox_res_Anjali_C60_R60_developing_seed_combined = as.data.frame(wilcox_res_Anjali_C60_R60_developing_seed_combined)
wilcox_res_Anjali_C60_R60_developing_seed_combined$Significance = with(wilcox_res_Anjali_C60_R60_developing_seed_combined,
                                                                       ifelse(p.value <= 0.001, "***", ifelse(
                                                                         p.value <= 0.01, "**", ifelse(
                                                                           p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_C60_R60_developing_seed_combined$Significance == "*") #8 metabolites
sum(wilcox_res_Anjali_C60_R60_developing_seed_combined$Significance == "**") #3 metabolites
sum(wilcox_res_Anjali_C60_R60_developing_seed_combined$Significance == "***") #4 metabolites
sum(wilcox_res_Anjali_C60_R60_developing_seed_combined$Significance != "ns") #15 metabolites
sum(wilcox_res_Anjali_C60_R60_developing_seed_combined$Significance == "ns") #52 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_C60_R60_developing_seed_combined = cbind(wilcox_res_N22_C60_R60_developing_seed_combined, 
                                                    wilcox_res_Dular_C60_R60_developing_seed_combined$Significance, 
                                                    wilcox_res_Anjali_C60_R60_developing_seed_combined$Significance)
wilcox_sig_C60_R60_developing_seed_combined = wilcox_sig_C60_R60_developing_seed_combined[ , -1]
colnames(wilcox_sig_C60_R60_developing_seed_combined)[1] = "wilcox_res_N22_C60_R60_developing_seed_combined$Significance"
wilcox_sig_C60_R60_developing_seed_combined = data.frame(lapply(wilcox_sig_C60_R60_developing_seed_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_C60_R60_developing_seed_combined) = rownames(wilcox_res_N22_C60_R60_developing_seed_combined)
wilcox_sig_C60_R60_developing_seed_combined$Nonsig = rowSums(wilcox_sig_C60_R60_developing_seed_combined == " ")  #count number of nonsignificant
wilcox_sig_final_C60_R60_developing_seed_combined = wilcox_sig_C60_R60_developing_seed_combined[-which(wilcox_sig_C60_R60_developing_seed_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 27 metabolites

wilcox_sig_log2fc_C60_R60_developing_seed_combined = log2fc_combined_control_recovery_developing_seed[which(rownames(log2fc_combined_control_recovery_developing_seed) %in% rownames(wilcox_sig_final_C60_R60_developing_seed_combined)), c(3, 6, 9)]  #log2FC values of metabolites significant in at least one of the cultivars

#overview - with HCA - Euclidean distance and average linkage 
heatmap.2(wilcox_sig_log2fc_C60_R60_developing_seed_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599),
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(5, 20), 
          main = "Developing seed-Combined-Wilcoxon \nLog2 FC, EGF, 60h RW/Control 60h", 
          cexCol = 1.5, 
          cexRow = 0.9, 
          density.info = "none",
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5),
          labCol = c("N22", "Dular", "Anjali"),
          srtCol = 0, 
          adjCol = c(0.5, 0.5),
          cellnote = wilcox_sig_final_C60_R60_developing_seed_combined,
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black",
          colsep = c(seq(0, 3)), 
          rowsep = c(seq(0, 27)),
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))


#venn diagram - common among cultivars

#N22
wilcox_sig_N22_C60_R60_developing_seed_combined = wilcox_res_N22_C60_R60_developing_seed_combined[-which(wilcox_res_N22_C60_R60_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_N22_C60_R60_developing_seed_combined = as.data.frame(wilcox_sig_log2fc_C60_R60_developing_seed_combined[match(rownames(wilcox_sig_N22_C60_R60_developing_seed_combined), rownames(wilcox_sig_log2fc_C60_R60_developing_seed_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_N22_C60_R60_developing_seed_combined) = colnames(wilcox_sig_log2fc_C60_R60_developing_seed_combined)[1]
wilcox_sig_log2fc_N22_C60_R60_developing_seed_combined$Up.Down = ifelse(wilcox_sig_log2fc_N22_C60_R60_developing_seed_combined$`GNR60 - GNC60` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_N22_C60_R60_developing_seed_combined$Up.Down == "Up")    #8 metabolites
sum(wilcox_sig_log2fc_N22_C60_R60_developing_seed_combined$Up.Down == "Down")  #5 metabolites

#Dular
wilcox_sig_Dular_C60_R60_developing_seed_combined = wilcox_res_Dular_C60_R60_developing_seed_combined[-which(wilcox_res_Dular_C60_R60_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Dular_C60_R60_developing_seed_combined = as.data.frame(wilcox_sig_log2fc_C60_R60_developing_seed_combined[match(rownames(wilcox_sig_Dular_C60_R60_developing_seed_combined), rownames(wilcox_sig_log2fc_C60_R60_developing_seed_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Dular_C60_R60_developing_seed_combined) = colnames(wilcox_sig_log2fc_C60_R60_developing_seed_combined)[2]
wilcox_sig_log2fc_Dular_C60_R60_developing_seed_combined$Up.Down = ifelse(wilcox_sig_log2fc_Dular_C60_R60_developing_seed_combined$`GDR60 - GDC60` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Dular_C60_R60_developing_seed_combined$Up.Down == "Up")    #6 metabolites
sum(wilcox_sig_log2fc_Dular_C60_R60_developing_seed_combined$Up.Down == "Down")  #6 metabolites

#Anjali
wilcox_sig_Anjali_C60_R60_developing_seed_combined = wilcox_res_Anjali_C60_R60_developing_seed_combined[-which(wilcox_res_Anjali_C60_R60_developing_seed_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Anjali_C60_R60_developing_seed_combined = as.data.frame(wilcox_sig_log2fc_C60_R60_developing_seed_combined[match(rownames(wilcox_sig_Anjali_C60_R60_developing_seed_combined), rownames(wilcox_sig_log2fc_C60_R60_developing_seed_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Anjali_C60_R60_developing_seed_combined) = colnames(wilcox_sig_log2fc_C60_R60_developing_seed_combined)[3]
wilcox_sig_log2fc_Anjali_C60_R60_developing_seed_combined$Up.Down = ifelse(wilcox_sig_log2fc_Anjali_C60_R60_developing_seed_combined$`GAR60 - GAC60` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Anjali_C60_R60_developing_seed_combined$Up.Down == "Up")    #10 metabolites
sum(wilcox_sig_log2fc_Anjali_C60_R60_developing_seed_combined$Up.Down == "Down")  #5 metabolites


#venn - increased - 60h rewatering/control 60h
vennlist_up_wilcox_C60_R60_developing_seed_combined = list(rownames(wilcox_sig_log2fc_N22_C60_R60_developing_seed_combined)[wilcox_sig_log2fc_N22_C60_R60_developing_seed_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Dular_C60_R60_developing_seed_combined)[wilcox_sig_log2fc_Dular_C60_R60_developing_seed_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Anjali_C60_R60_developing_seed_combined)[wilcox_sig_log2fc_Anjali_C60_R60_developing_seed_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_C60_R60_developing_seed_combined) = c("N22", "Dular", "Anjali")
venn_up_wilcox_C60_R60_developing_seed_combined = venn(vennlist_up_wilcox_C60_R60_developing_seed_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_R60_C60_developing_seed_combined = venn.diagram(vennlist_up_wilcox_C60_R60_developing_seed_combined, category.names = c(paste(names(vennlist_up_wilcox_C60_R60_developing_seed_combined[1]), " (", length(vennlist_up_wilcox_C60_R60_developing_seed_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_C60_R60_developing_seed_combined[2]), " (", length(vennlist_up_wilcox_C60_R60_developing_seed_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_C60_R60_developing_seed_combined[3]), " (", length(vennlist_up_wilcox_C60_R60_developing_seed_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "E", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label E
intersect_wilcox_C60_R60_venn_up_developing_seed_combined = attr(venn_up_wilcox_C60_R60_developing_seed_combined, "intersections")  #list of unique and common metabolites
venn_up_wilcox_metabolites_C60_R60_developing_seed_combined = t(ldply(intersect_wilcox_C60_R60_venn_up_developing_seed_combined, rbind))  #convert from list to matrix
colnames(venn_up_wilcox_metabolites_C60_R60_developing_seed_combined) = venn_up_wilcox_metabolites_C60_R60_developing_seed_combined[1, ]
venn_up_wilcox_metabolites_C60_R60_developing_seed_combined = venn_up_wilcox_metabolites_C60_R60_developing_seed_combined[-1, ]


#venn - decreased - 60h rewatering/control 60h
vennlist_down_wilcox_C60_R60_developing_seed_combined = list(rownames(wilcox_sig_log2fc_N22_C60_R60_developing_seed_combined)[wilcox_sig_log2fc_N22_C60_R60_developing_seed_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Dular_C60_R60_developing_seed_combined)[wilcox_sig_log2fc_Dular_C60_R60_developing_seed_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Anjali_C60_R60_developing_seed_combined)[wilcox_sig_log2fc_Anjali_C60_R60_developing_seed_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_C60_R60_developing_seed_combined) = c("N22", "Dular", "Anjali")
venn_down_wilcox_C60_R60_developing_seed_combined = venn(vennlist_down_wilcox_C60_R60_developing_seed_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_R60_C60_developing_seed_combined = venn.diagram(vennlist_down_wilcox_C60_R60_developing_seed_combined, category.names = c(paste(names(vennlist_down_wilcox_C60_R60_developing_seed_combined[1]), " (", length(vennlist_down_wilcox_C60_R60_developing_seed_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_C60_R60_developing_seed_combined[2]), " (", length(vennlist_down_wilcox_C60_R60_developing_seed_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_C60_R60_developing_seed_combined[3]), " (", length(vennlist_down_wilcox_C60_R60_developing_seed_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "F", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label F
intersect_wilcox_C60_R60_venn_down_developing_seed_combined = attr(venn_down_wilcox_C60_R60_developing_seed_combined, "intersections")  #list of unique and common metabolites
venn_down_wilcox_metabolites_C60_R60_developing_seed_combined = t(ldply(intersect_wilcox_C60_R60_venn_down_developing_seed_combined, rbind))  #convert from list to matrix
colnames(venn_down_wilcox_metabolites_C60_R60_developing_seed_combined) = venn_down_wilcox_metabolites_C60_R60_developing_seed_combined[1, ]
venn_down_wilcox_metabolites_C60_R60_developing_seed_combined = venn_down_wilcox_metabolites_C60_R60_developing_seed_combined[-1, ]


write.table(venn_up_wilcox_metabolites_C60_R60_developing_seed_combined, file = "wilcox_sig_up_metabolites_developing_seed_R60_C60_combined.txt", sep = "\t", quote = F)
write.table(venn_down_wilcox_metabolites_C60_R60_developing_seed_combined, file = "wilcox_sig_down_metabolites_developing_seed_R60_C60_combined.txt", sep = "\t", quote = F)


####################


#rownames were checked to be the same for the data files to be combined
#significance of wilcoxon test - all cultivars, all three RW timepoints - including LS/C (refer to previous paper and associated GitHub files)
wilcox_sig_C_LS_R12_R36_R60_developing_seed_combined = cbind(wilcox_res_N22_CLS_LS_developing_seed_combined,
                                                             wilcox_res_N22_C12_R12_developing_seed_combined$Significance,
                                                             wilcox_res_N22_C36_R36_developing_seed_combined$Significance,
                                                             wilcox_res_N22_C60_R60_developing_seed_combined$Significance,
                                                             wilcox_res_Dular_CLS_LS_developing_seed_combined$Significance,
                                                             wilcox_res_Dular_C12_R12_developing_seed_combined$Significance,
                                                             wilcox_res_Dular_C36_R36_developing_seed_combined$Significance,
                                                             wilcox_res_Dular_C60_R60_developing_seed_combined$Significance,
                                                             wilcox_res_Anjali_CLS_LS_developing_seed_combined$Significance,
                                                             wilcox_res_Anjali_C12_R12_developing_seed_combined$Significance,
                                                             wilcox_res_Anjali_C36_R36_developing_seed_combined$Significance,
                                                             wilcox_res_Anjali_C60_R60_developing_seed_combined$Significance)
wilcox_sig_C_LS_R12_R36_R60_developing_seed_combined = wilcox_sig_C_LS_R12_R36_R60_developing_seed_combined[ , -1]
colnames(wilcox_sig_C_LS_R12_R36_R60_developing_seed_combined)[1] = "wilcox_res_N22_CLS_LS_developing_seed_combined$Significance"
wilcox_sig_C_LS_R12_R36_R60_developing_seed_combined = data.frame(lapply(wilcox_sig_C_LS_R12_R36_R60_developing_seed_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_C_LS_R12_R36_R60_developing_seed_combined) = rownames(wilcox_res_N22_CLS_LS_developing_seed_combined)
wilcox_sig_C_LS_R12_R36_R60_developing_seed_combined$Nonsig = rowSums(wilcox_sig_C_LS_R12_R36_R60_developing_seed_combined == " ")  #count number of nonsignificant
wilcox_sig_final_C_LS_R12_R36_R60_developing_seed_combined = wilcox_sig_C_LS_R12_R36_R60_developing_seed_combined[-which(wilcox_sig_C_LS_R12_R36_R60_developing_seed_combined$Nonsig == "12"), ] #metabolites significant in at least one of the comparisons - 51 metabolites

#log2fc values of late stress and recovery relative to respective controls (with respect to stress paper and associated GitHub files)
rownames(log2fc_combined_control_stress_developing_seed) == rownames(log2fc_combined_control_recovery_developing_seed)  #double-check if rownames are the same
log2fc_combined_controls_latestress_recovery_developing_seed = cbind(log2fc_combined_control_stress_developing_seed, 
                                                                     log2fc_combined_control_recovery_developing_seed)   #combine in one data set
log2fc_combined_controls_latestress_recovery_developing_seed = log2fc_combined_controls_latestress_recovery_developing_seed[ , c(1, 4:6, 2, 7:9, 3, 10:12)]   #rearrange in preferred order


#log2fc values of metabolites with significant difference in at least one of the comparisons - late stress, 12, 36, 60h recovery vs. corresponding controls in any cultivar
wilcox_sig_log2fc_C_LS_R12_R36_R60_developing_seed_combined = log2fc_combined_controls_latestress_recovery_developing_seed[which(rownames(log2fc_combined_controls_latestress_recovery_developing_seed) %in% rownames(wilcox_sig_final_C_LS_R12_R36_R60_developing_seed_combined)), ]


#rewatering-responsive metabolites only - indicate as red font in heat map
setdiff(rownames(wilcox_sig_log2fc_C_LS_R12_R36_R60_developing_seed_combined), rownames(wilcox_sig_log2fc_CLS_LS_developing_seed_combined))

#visualization
png("wilcox_sig_log2fc_C_LS_RW_developing_seed_combined.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_log2fc_C_LS_R12_R36_R60_developing_seed_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(7, 18.5),
          lmat = rbind(c(0, 4, 5), c(0, 1, 0), c(3, 2, 0)),
          labCol = rep(c("Severe stress/Control", "12 h RW/Control", "36 h RW/Control", "60 h RW/Control"), 3),
          cexCol = 1, 
          cexRow = 1, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 0.2, 11), 
          lwid = c(1, 5, 1.6),
          srtCol = 45, 
          cellnote = wilcox_sig_final_C_LS_R12_R36_R60_developing_seed_combined[ , c(1:12)], 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 12)), 
          rowsep = c(seq(0, 51)), 
          sepwidth = c(0.01, 0.01), 
          ColSideColors = c(rep("#66c2a5", 4), rep("#fc8d62", 4), rep("#8da0cb", 4)),
          hclustfun = function(x) hclust(x, method = "average"),
          colRow = c(rep("red", 4), "black", "red", rep("black", 2), "red", rep("black", 3), "red", rep("black", 2), "red",
                     "black", rep("red", 3), rep("black", 9), rep("red", 4), rep("black", 2), rep("red", 4), rep("black", 2), "red", 
                     rep("black", 3), "red", rep("black", 2), rep("red", 2), "black"))
legend(x = 0.026, y = 1.1, xpd = TRUE, legend = "N22", bty = "n", cex = 1.2)
legend(x = 0.155, y = 1.1, xpd = TRUE, legend = "Dular", bty = "n", cex = 1.2)
legend(x = 0.295, y = 1.1, xpd = TRUE, legend = "Anjali", bty = "n", cex = 1.2)
legend(x = 0.863, y = 1.04, xpd = TRUE, legend = expression(paste("Log" ["2"], "-fold change")), bty = "n", cex = 0.8)
dev.off()


####################

#venn diagrams

#multiplot
png("venn_up_down_developing_seed_RW_controls_combined.png", width = 16*300, height = 24*300, res = 300)
grid.arrange(gTree(children = venn_up_R12_C12_developing_seed_combined), 
             gTree(children = venn_down_R12_C12_developing_seed_combined), 
             gTree(children = venn_up_R36_C36_developing_seed_combined),
             gTree(children = venn_down_R36_C36_developing_seed_combined),
             gTree(children = venn_up_R60_C60_developing_seed_combined),
             gTree(children = venn_down_R60_C60_developing_seed_combined),
             ncol = 2, nrow = 3)
dev.off()


####################


#correlation analysis

#correlation between change in yield and metabolite changes during 60h RW
#correlation between change in proportion of chalky grains (>50% chalk content) and metabolite changes during 60h RW

#only 60h RW was considered since it is the latest sample collection time point


#yield and chalkiness - mean per cultivar per treatment per year already ran previously

stress_R60_order = c("GNLS", "GNR60", "GDLS", "GDR60", "GALS", "GAR60")


#log2fc - 60h RW / LS

#2013
data_LS_R60_log2_mean_developing_seed_2013 = data_log2_mean_developing_seed_2013_2[stress_R60_order, ]  #data with stress and recovery60 values only
data_LS_R60_log2_mean_developing_seed_2013 = t(data_LS_R60_log2_mean_developing_seed_2013)

log2fc_N22_stress_recovery_developing_seed_2013 = data_LS_R60_log2_mean_developing_seed_2013[ , 2] - data_LS_R60_log2_mean_developing_seed_2013[ , 1]
log2fc_Dular_stress_recovery_developing_seed_2013 = data_LS_R60_log2_mean_developing_seed_2013[ , 4] - data_LS_R60_log2_mean_developing_seed_2013[ , 3]
log2fc_Anjali_stress_recovery_developing_seed_2013 = data_LS_R60_log2_mean_developing_seed_2013[ , 6] - data_LS_R60_log2_mean_developing_seed_2013[ , 5]
log2fc_stress_recovery_developing_seed_2013 = as.data.frame(cbind(log2fc_N22_stress_recovery_developing_seed_2013, log2fc_Dular_stress_recovery_developing_seed_2013, log2fc_Anjali_stress_recovery_developing_seed_2013))
colnames(log2fc_stress_recovery_developing_seed_2013) = c("GNR60 - GNLS", "GDR60 - GDLS", "GAR60 - GALS")


#2014
data_LS_R60_log2_mean_developing_seed_2014 = data_log2_mean_developing_seed_2014_2[stress_R60_order, ]  #data with stress and recovery60 values only
data_LS_R60_log2_mean_developing_seed_2014 = t(data_LS_R60_log2_mean_developing_seed_2014)

log2fc_N22_stress_recovery_developing_seed_2014 = data_LS_R60_log2_mean_developing_seed_2014[ , 2] - data_LS_R60_log2_mean_developing_seed_2014[ , 1]
log2fc_Dular_stress_recovery_developing_seed_2014 = data_LS_R60_log2_mean_developing_seed_2014[ , 4] - data_LS_R60_log2_mean_developing_seed_2014[ , 3]
log2fc_Anjali_stress_recovery_developing_seed_2014 = data_LS_R60_log2_mean_developing_seed_2014[ , 6] - data_LS_R60_log2_mean_developing_seed_2014[ , 5]
log2fc_stress_recovery_developing_seed_2014 = as.data.frame(cbind(log2fc_N22_stress_recovery_developing_seed_2014, log2fc_Dular_stress_recovery_developing_seed_2014, log2fc_Anjali_stress_recovery_developing_seed_2014))
colnames(log2fc_stress_recovery_developing_seed_2014) = c("GNR60 - GNLS", "GDR60 - GDLS", "GAR60 - GALS")


#2015
data_LS_R60_log2_mean_developing_seed_2015 = data_log2_mean_developing_seed_2015_2[stress_R60_order, ]  #data with stress and recovery60 values only
data_LS_R60_log2_mean_developing_seed_2015 = t(data_LS_R60_log2_mean_developing_seed_2015)

log2fc_N22_stress_recovery_developing_seed_2015 = data_LS_R60_log2_mean_developing_seed_2015[ , 2] - data_LS_R60_log2_mean_developing_seed_2015[ , 1]
log2fc_Dular_stress_recovery_developing_seed_2015 = data_LS_R60_log2_mean_developing_seed_2015[ , 4] - data_LS_R60_log2_mean_developing_seed_2015[ , 3]
log2fc_Anjali_stress_recovery_developing_seed_2015 = data_LS_R60_log2_mean_developing_seed_2015[ , 6] - data_LS_R60_log2_mean_developing_seed_2015[ , 5]
log2fc_stress_recovery_developing_seed_2015 = as.data.frame(cbind(log2fc_N22_stress_recovery_developing_seed_2015, log2fc_Dular_stress_recovery_developing_seed_2015, log2fc_Anjali_stress_recovery_developing_seed_2015))
colnames(log2fc_stress_recovery_developing_seed_2015) = c("GNR60 - GNLS", "GDR60 - GDLS", "GAR60 - GALS")


#check first if rownames are the same before combining data
rownames(log2fc_stress_recovery_developing_seed_2013) == rownames(log2fc_stress_recovery_developing_seed_2014)
rownames(log2fc_stress_recovery_developing_seed_2013) == rownames(log2fc_stress_recovery_developing_seed_2015)

log2fc_developing_seed_R60_latestress_20131415 = cbind(log2fc_stress_recovery_developing_seed_2013[ , ncol(log2fc_stress_recovery_developing_seed_2013):1], 
                                                       log2fc_stress_recovery_developing_seed_2014[ , ncol(log2fc_stress_recovery_developing_seed_2014):1], 
                                                       log2fc_stress_recovery_developing_seed_2015[ , ncol(log2fc_stress_recovery_developing_seed_2015):1])    #log2fc of 60h RW / late stress - 3-yr data in one data table - order: Anjali, Dular, N22
colnames(log2fc_developing_seed_R60_latestress_20131415) = c("2013 Anjali", "2013 Dular", "2013 N22", "2014 Anjali", "2014 Dular", "2014 N22", "2015 Anjali", "2015 Dular", "2015 N22")
log2fc_developing_seed_R60_latestress_20131415 = log2fc_developing_seed_R60_latestress_20131415[ , match(colnames(yield_chalk_change_HxD_EGF), colnames(log2fc_developing_seed_R60_latestress_20131415))]     #reorder columns to match column order of yield and chalkiness data


#yield, chalkiness, and log2FC in one data table

#check first if rownames are the same before combining data
colnames(yield_chalk_change_HxD_EGF) == colnames(log2fc_developing_seed_R60_latestress_20131415)

yield_chalk_log2fc_R60_LS_developing_seed_HxD = rbind(yield_chalk_change_HxD_EGF, log2fc_developing_seed_R60_latestress_20131415)
yield_chalk_log2fc_R60_LS_developing_seed_HxD = t(yield_chalk_log2fc_R60_LS_developing_seed_HxD)
yield_chalk_log2fc_R60_LS_developing_seed_HxD = as.data.frame(yield_chalk_log2fc_R60_LS_developing_seed_HxD)

#test for normality of data
shapiro_yield_chalk_log2fc_R60_LS_developing_seed_HxD = func_shapiro_test(yield_chalk_log2fc_R60_LS_developing_seed_HxD)
shapiro_yield_chalk_log2fc_R60_LS_developing_seed_HxD$Distribution = ifelse(shapiro_yield_chalk_log2fc_R60_LS_developing_seed_HxD$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_yield_chalk_log2fc_R60_LS_developing_seed_HxD$Distribution == "Normal")   #61 variables
sum(shapiro_yield_chalk_log2fc_R60_LS_developing_seed_HxD$Distribution == "Non-normal")   #9 variables


#correlation test - spearman

#correlation coefficient

#accounts for change in yield/chalkiness
cor_res_spearman_rho_R60_LS_developing_seed_HxD = cor(yield_chalk_log2fc_R60_LS_developing_seed_HxD, method = "spearman")


#correlation with change in yield and chalkiness

#yield
yield_cor_res_spearman_pval_R60_LS_developing_seed_HxD = t(data.frame(lapply(yield_chalk_log2fc_R60_LS_developing_seed_HxD[4:ncol(yield_chalk_log2fc_R60_LS_developing_seed_HxD)], function(x) cor.test(yield_chalk_log2fc_R60_LS_developing_seed_HxD$`Change in yield`, x, method = "spearman")$p.value)))    #significance of correlation
yield_cor_res_spearman_pval_R60_LS_developing_seed_HxD = as.data.frame(yield_cor_res_spearman_pval_R60_LS_developing_seed_HxD)
colnames(yield_cor_res_spearman_pval_R60_LS_developing_seed_HxD) = "p.value"
rownames(yield_cor_res_spearman_pval_R60_LS_developing_seed_HxD) = reduced_metabolite_list_final_developing_seed_2013$Name
yield_cor_res_spearman_pval_R60_LS_developing_seed_HxD$Significance = with(yield_cor_res_spearman_pval_R60_LS_developing_seed_HxD,
                                                                           ifelse(p.value <= 0.001, "***", ifelse(
                                                                             p.value <= 0.01, "**", ifelse(
                                                                               p.value < 0.05, "*", "ns"))))
length(which(yield_cor_res_spearman_pval_R60_LS_developing_seed_HxD$p.value < 0.05))   #15 metabolites
yield_cor_sig_metabolites_R60_LS_developing_seed_HxD = rownames(yield_cor_res_spearman_pval_R60_LS_developing_seed_HxD)[yield_cor_res_spearman_pval_R60_LS_developing_seed_HxD$Significance != "ns"]   #metabolites with significant correlation with change in yield
yield_sig_cor_rho_spearman_R60_LS_developing_seed_HxD = data.frame(cor_res_spearman_rho_R60_LS_developing_seed_HxD[1, which(colnames(cor_res_spearman_rho_R60_LS_developing_seed_HxD) %in% yield_cor_sig_metabolites_R60_LS_developing_seed_HxD)])     #metabolites with significant correlation with change in yield and rho values
yield_sig_cor_rho_spearman_R60_LS_developing_seed_HxD = round(yield_sig_cor_rho_spearman_R60_LS_developing_seed_HxD, 2)   #round to 2 decimals
yield_sig_cor_rho_spearman_R60_LS_developing_seed_HxD = as.matrix(yield_sig_cor_rho_spearman_R60_LS_developing_seed_HxD)  #convert to matrix
yield_sig_cor_rho_spearman_R60_LS_developing_seed_HxD = cbind(yield_sig_cor_rho_spearman_R60_LS_developing_seed_HxD, yield_sig_cor_rho_spearman_R60_LS_developing_seed_HxD)   #make 2-column data frame to be able to sort rownames alphabetically
yield_sig_cor_rho_spearman_R60_LS_developing_seed_HxD = yield_sig_cor_rho_spearman_R60_LS_developing_seed_HxD[order(rownames(yield_sig_cor_rho_spearman_R60_LS_developing_seed_HxD)), ]   #sort metabolites alphabetially
yield_sig_cor_rho_spearman_R60_LS_developing_seed_HxD = as.data.frame(yield_sig_cor_rho_spearman_R60_LS_developing_seed_HxD[ , -2])
colnames(yield_sig_cor_rho_spearman_R60_LS_developing_seed_HxD) = "rho value"
#export table
write.table(yield_sig_cor_rho_spearman_R60_LS_developing_seed_HxD, file = "cor_yield_sig_R60_LS_developing_seed_HxD.txt", sep = "\t", quote = F)


#chalky grains
chalk_50_75_cor_res_spearman_pval_R60_LS_developing_seed_HxD = t(data.frame(lapply(yield_chalk_log2fc_R60_LS_developing_seed_HxD[4:ncol(yield_chalk_log2fc_R60_LS_developing_seed_HxD)], function(x) cor.test(yield_chalk_log2fc_R60_LS_developing_seed_HxD$`Change in chalky_50-75%`, x, method = "spearman")$p.value)))    #significance of correlation
chalk_50_75_cor_res_spearman_pval_R60_LS_developing_seed_HxD = as.data.frame(chalk_50_75_cor_res_spearman_pval_R60_LS_developing_seed_HxD)
colnames(chalk_50_75_cor_res_spearman_pval_R60_LS_developing_seed_HxD) = "p.value"
rownames(chalk_50_75_cor_res_spearman_pval_R60_LS_developing_seed_HxD) = reduced_metabolite_list_final_developing_seed_2013$Name
chalk_50_75_cor_res_spearman_pval_R60_LS_developing_seed_HxD$Significance = with(chalk_50_75_cor_res_spearman_pval_R60_LS_developing_seed_HxD,
                                                                                 ifelse(p.value <= 0.001, "***", ifelse(
                                                                                   p.value <= 0.01, "**", ifelse(
                                                                                     p.value < 0.05, "*", "ns"))))
length(which(chalk_50_75_cor_res_spearman_pval_R60_LS_developing_seed_HxD$p.value < 0.05))   #1 metabolite
chalk_50_75_cor_sig_metabolites_R60_LS_developing_seed_HxD = rownames(chalk_50_75_cor_res_spearman_pval_R60_LS_developing_seed_HxD)[chalk_50_75_cor_res_spearman_pval_R60_LS_developing_seed_HxD$Significance != "ns"]   #metabolites with significant correlation with change in chalk_50_75
chalk_50_75_sig_cor_rho_spearman_R60_LS_developing_seed_HxD = data.frame(cor_res_spearman_rho_R60_LS_developing_seed_HxD[3, which(colnames(cor_res_spearman_rho_R60_LS_developing_seed_HxD) %in% chalk_50_75_cor_sig_metabolites_R60_LS_developing_seed_HxD)])     #metabolites with significant correlation with change in chalk_50_75 and rho values
rownames(chalk_50_75_sig_cor_rho_spearman_R60_LS_developing_seed_HxD) = chalk_50_75_cor_sig_metabolites_R60_LS_developing_seed_HxD
chalk_50_75_sig_cor_rho_spearman_R60_LS_developing_seed_HxD = round(chalk_50_75_sig_cor_rho_spearman_R60_LS_developing_seed_HxD, 2)   #round to 2 decimals
colnames(chalk_50_75_sig_cor_rho_spearman_R60_LS_developing_seed_HxD) = "rho value"
#export table
write.table(chalk_50_75_sig_cor_rho_spearman_R60_LS_developing_seed_HxD, file = "cor_chalk_sig_R60_LS_developing_seed_HxD.txt", sep = "\t", quote = F)
