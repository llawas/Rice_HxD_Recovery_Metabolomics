#load required packages and functions
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


##################################################
#PRINCIPAL COMPONENT ANALYSIS#
##################################################

#PCA of flag leaf, flowering spikelet, and developing seed - control, late (severe) stress, 12 h rewatering (RW), 36 h RW, 60 h RW

#flag leaf - Flowering stage

#exclude early stress
remove_earlystress = c("FAES", "FDES", "FNES")

data_combined_log2_mean_C_LS_R_flag_leaf_FL = data_combined_log2_mean_FL_flag_leaf[!(data_combined_log2_mean_FL_flag_leaf$Code) %in% remove_earlystress, ]    #data without early stress
data_combined_log2_mean_C_LS_R_flag_leaf_FL_2 = data_combined_log2_mean_FL_flag_leaf_2[!rownames(data_combined_log2_mean_FL_flag_leaf_2) %in% remove_earlystress, ]    #data without early stress

#PCA
pca_res_combined_C_LS_R_flag_leaf_FL = pca(data_combined_log2_mean_C_LS_R_flag_leaf_FL_2, 
                                           method = "ppca", 
                                           nPcs = 5, 
                                           scale = "pareto", 
                                           center = TRUE)   
pca_res_scores_combined_C_LS_R_flag_leaf_FL = scores(pca_res_combined_C_LS_R_flag_leaf_FL)    #scores
pca_res_loadings_combined_C_LS_R_flag_leaf_FL = loadings(pca_res_combined_C_LS_R_flag_leaf_FL)  #loadings

#scores + additional info as data frame to use in score plot
scores_combined_C_LS_R_flag_leaf_FL = as.data.frame(pca_res_scores_combined_C_LS_R_flag_leaf_FL)
scores_combined_C_LS_R_flag_leaf_FL$Timepoint = data_combined_log2_mean_C_LS_R_flag_leaf_FL$Timepoint
scores_combined_C_LS_R_flag_leaf_FL$Cultivar = data_combined_log2_mean_C_LS_R_flag_leaf_FL$Cultivar
scores_combined_C_LS_R_flag_leaf_FL$Label = sub(".", "", rownames(scores_combined_C_LS_R_flag_leaf_FL))   #remove first character from rownames

#biplot
pdf("biplot_pca_combined_C_LS_R_flag_leaf_FL.pdf")
biplot(pca_res_combined_C_LS_R_flag_leaf_FL, cex = 0.8, main = "Combined_Control, Late stress, Recovery \nFlag leaf - FL stage")
dev.off()

#screeplot
pdf("screeplot_pca_combined_C_LS_R_flag_leaf_FL.pdf")
text(barplot(pca_res_combined_C_LS_R_flag_leaf_FL@R2*100, names.arg = c("PC1", "PC2", "PC3", "PC4", "PC5"), main = "Combined_Control, Late stress, Recovery \nFlag leaf - FL stage", ylab = "Variance (%)", ylim = c(0, 40)), 0, round(pca_res_combined_C_LS_R_flag_leaf_FL@R2*100, 3), pos = 3)
box()
dev.off()

#colors for PCA score plot
display.brewer.pal(12, "Paired")    #diplays colors - choose from this
colors_pca = brewer.pal(12, "Paired")     
colors_pca     #shows the code for each color

#score plot
PCA_combined_flag_leaf_FL = ggplot(scores_combined_C_LS_R_flag_leaf_FL, 
                                   aes(x = PC1, 
                                       y = PC2, 
                                       color = Timepoint, 
                                       fill = Timepoint, 
                                       shape = Cultivar, 
                                       label = Label)) +
  geom_point(size = 8) +
  scale_color_manual(breaks = c("Control", "Late stress", "12h Rewatering", "36h Rewatering", "60h Rewatering"),
                     values = c("#6A3D9A", "#B15928", "#FF7F00", "#E31A1C", "#FB9A99"),
                     labels = c("Control", "Severe stress", "12 h rewatering", "36 h rewatering", "60 h rewatering")) +
  scale_fill_manual(breaks = c("Control", "Late stress", "12h Rewatering", "36h Rewatering", "60h Rewatering"),
                    values = c("#6A3D9A", "#B15928", "#FF7F00", "#E31A1C", "#FB9A99"),
                    labels = c("Control", "Severe stress", "12 h rewatering", "36 h rewatering", "60 h rewatering")) +
  scale_shape_manual(breaks = c("N22", "Dular", "Anjali"),
                     values = c(22, 24, 21),
                     labels = c("N22  ", "Dular  ", "Anjali")) +
  scale_x_continuous(limits = c(-8, 6), breaks = seq(-8, 6, by = 2)) +
  scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, by = 2)) +
  xlab(paste("PC1 (", round(pca_res_combined_C_LS_R_flag_leaf_FL@R2[1], 3)*100, "%)", sep = "")) +
  ylab(paste("PC2 (", round(pca_res_combined_C_LS_R_flag_leaf_FL@R2[2], 3)*100, "%)", sep = "")) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_blank(),
        legend.justification = c(1, 0),
        legend.box.just = "top",
        legend.position = c(0.995, 0.8),
        legend.background = element_blank(),
        legend.text = element_text(size = 10),
        legend.key = element_blank(),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.spacing = unit(24.5, "lines"),
        legend.key.height = unit(1.3, "lines")) +
  guides(shape = guide_legend(order = 1, override.aes = list(size = 5)),
         fill = guide_legend(override.aes = list(size = 5)))


#flag leaf - Early grain filling (EGF) stage

#PCA
pca_res_combined_controls_LS_R_flag_leaf_EGF = pca(data_combined_log2_mean_EGF_flag_leaf_2, method = "ppca", nPcs = 5, scale = "pareto", center = TRUE)   
pca_res_scores_combined_controls_LS_R_flag_leaf_EGF = scores(pca_res_combined_controls_LS_R_flag_leaf_EGF)    #scores
pca_res_loadings_combined_controls_LS_R_flag_leaf_EGF = loadings(pca_res_combined_controls_LS_R_flag_leaf_EGF)  #loadings

#scores + additional info as data frame to use in score plot
scores_combined_controls_LS_R_flag_leaf_EGF = as.data.frame(pca_res_scores_combined_controls_LS_R_flag_leaf_EGF)
scores_combined_controls_LS_R_flag_leaf_EGF$Timepoint = interaction(data_combined_log2_mean_EGF_flag_leaf$Treatment, data_combined_log2_mean_EGF_flag_leaf$Timepoint)
scores_combined_controls_LS_R_flag_leaf_EGF$Cultivar = data_combined_log2_mean_EGF_flag_leaf$Cultivar
scores_combined_controls_LS_R_flag_leaf_EGF$Label = gsub("G", "", rownames(scores_combined_controls_LS_R_flag_leaf_EGF))

#biplot
pdf("biplot_pca_combined_controls_LS_R_flag_leaf_EGF.pdf")
biplot(pca_res_combined_controls_LS_R_flag_leaf_EGF, cex = 0.8, main = "Combined_Controls, Late stress, Recovery \nFlag leaf - EGF stage")
dev.off()

#screeplot
pdf("screeplot_pca_combined_cotrols_LS_R_flag_leaf_EGF.pdf")
text(barplot(pca_res_combined_controls_LS_R_flag_leaf_EGF@R2*100, names.arg = c("PC1", "PC2", "PC3", "PC4", "PC5"), main = "Combined_Controls, Late stress, Recovery \nFlag leaf - EGF stage", ylab = "Variance (%)", ylim = c(0, 40)), 0, round(pca_res_combined_controls_LS_R_flag_leaf_EGF@R2*100, 3), pos = 3)
box()
dev.off()

#score plot
PCA_combined_flag_leaf_EGF = ggplot(scores_combined_controls_LS_R_flag_leaf_EGF, 
                                    aes(x = PC1, 
                                        y = PC2, 
                                        color = Timepoint, 
                                        fill = Timepoint, 
                                        shape = Cultivar, 
                                        label = Label)) +
  geom_point(size = 8) +
  scale_color_manual(breaks = c("Control.Late stress", "Control.12h Rewatering", "Control.36h Rewatering", "Control.60h Rewatering",
                                "Heat & drought.Late stress", "Heat & drought.12h Rewatering", "Heat & drought.36h Rewatering", 
                                "Heat & drought.60h Rewatering"),
                     values = c("#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FDBF6F", "#FF7F00", "#FB9A99", "#E31A1C"),
                     labels = c("Control - Severe stress", "Control - 12 h rewatering", "Control - 36 h rewatering", 
                                "Control - 60 h rewatering", "Severe stress", "12 h rewatering", "36 h rewatering", "60 h rewatering")) +
  scale_fill_manual(breaks = c("Control.Late stress", "Control.12h Rewatering", "Control.36h Rewatering", "Control.60h Rewatering",
                               "Heat & drought.Late stress", "Heat & drought.12h Rewatering", "Heat & drought.36h Rewatering", 
                               "Heat & drought.60h Rewatering"),
                    values = c("#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FDBF6F", "#FF7F00", "#FB9A99", "#E31A1C"),
                    labels = c("Control - Severe stress", "Control - 12 h rewatering", "Control - 36 h rewatering", 
                               "Control - 60 h rewatering", "Severe stress", "12 h rewatering", "36 h rewatering", "60 h rewatering")) +
  scale_shape_manual(breaks = c("N22", "Dular", "Anjali"),
                     values = c(22, 24, 21),
                     labels = c("N22", "Dular", "Anjali")) +
  scale_x_continuous(limits = c(-8, 8), breaks = seq(-8, 8, by = 2)) +
  scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, by = 2)) +
  xlab(paste("PC1 (", round(pca_res_combined_controls_LS_R_flag_leaf_EGF@R2[1], 3)*100, "%)", sep = "")) +
  ylab(paste("PC2 (", round(pca_res_combined_controls_LS_R_flag_leaf_EGF@R2[2], 3)*100, "%)", sep = "")) +
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
        legend.spacing = unit(15.5, "lines"),
        legend.key.height = unit(1.3, "lines")) +
  guides(shape = guide_legend(order = 1, override.aes = list(size = 5)),
         fill = guide_legend(override.aes = list(size = 5),
                             ncol = 2),
         color = guide_legend(ncol = 2))


#flowering spikelet

#for txt file, refer to Rice_HxD_GC-MS_Recovery_2 GitHub file
PCA_scores_combined_flowering_spikelet = read.table("PCA_scores_combined_flowering_spikelet.txt", header = TRUE, sep = "\t")


#check order of factor levels first before assigning colors/shapes
PCA_combined_flowering_spikelet = ggplot(PCA_scores_combined_flowering_spikelet, 
                                         aes(x = PC1, 
                                             y = PC2, 
                                             color = Timepoint, 
                                             fill = Timepoint, 
                                             shape = Cultivar, 
                                             label = Label)) +
  geom_point(size = 8) +
  scale_color_manual(breaks = c("Control", "Late stress", "12h Rewatering", "36h Rewatering", "60h Rewatering"),
                     values = c("#6A3D9A", "#B15928", "#FF7F00", "#FB9A99", "#E31A1C"),
                     labels = c("Control", "Severe stress", "12 h rewatering", "36 h rewatering", "60 h rewatering")) +
  scale_fill_manual(breaks = c("Control", "Late stress", "12h Rewatering", "36h Rewatering", "60h Rewatering"),
                    values = c("#6A3D9A", "#B15928", "#FF7F00", "#FB9A99", "#E31A1C"),
                    labels = c("Control", "Severe stress", "12 h rewatering", "36 h rewatering", "60 h rewatering")) +
  scale_shape_manual(breaks = c("N22", "Dular", "Anjali"),
                     values = c(22, 24, 21),
                     labels = c("N22", "Dular", "Anjali")) +
  scale_x_continuous(limits = c(-8, 6), breaks = seq(-8, 6, by = 2)) +
  scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, by = 2)) +
  labs(x = "PC1 (31.4%)", y = "PC2 (20.6%)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_blank(),
        legend.justification = c(1, 0),
        legend.box.just = "top",
        legend.position = c(0.995, 0.8),
        legend.background = element_blank(),
        legend.text = element_text(size = 10),
        legend.key = element_blank(),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.spacing = unit(24.5, "lines"),
        legend.key.height = unit(1.3, "lines")) +
  guides(shape = guide_legend(order = 1, override.aes = list(size = 5)),
         fill = guide_legend(override.aes = list(size = 5)))


#developing seed

#for txt file, refer to Rice_HxD_GC-MS_Recovery_3 GitHub file
PCA_scores_combined_developing_seed = read.table("PCA_scores_combined_developing_seed.txt", header = TRUE, sep = "\t")


#check order of factor levels first before assigning colors/shapes
PCA_combined_developing_seed = ggplot(PCA_scores_combined_developing_seed, 
                                      aes(x = PC1, 
                                          y = PC2, 
                                          color = Timepoint, 
                                          fill = Timepoint, 
                                          shape = Cultivar, 
                                          label = Label)) +
  geom_point(size = 8) +
  scale_color_manual(breaks = c("Control.Late stress", "Control.12h Rewatering", "Control.36h Rewatering", "Control.60h Rewatering",
                                "Heat & drought.Late stress", "Heat & drought.12h Rewatering", "Heat & drought.36h Rewatering", 
                                "Heat & drought.60h Rewatering"),
                     values = c("#CAB2D6", "#FFFF99", "#FDBF6F", "#FB9A99", "#6A3D9A", "#B15928", "#FF7F00", "#E31A1C"),
                     labels = c("Control - Severe stress", "Control - 12 h rewatering", "Control - 36 h rewatering", 
                                "Control - 60 h rewatering", "Severe stress", "12 h rewatering", "36 h rewatering", "60 h rewatering")) +
  scale_fill_manual(breaks = c("Control.Late stress", "Control.12h Rewatering", "Control.36h Rewatering", "Control.60h Rewatering",
                               "Heat & drought.Late stress", "Heat & drought.12h Rewatering", "Heat & drought.36h Rewatering", 
                               "Heat & drought.60h Rewatering"),
                    values = c("#CAB2D6", "#FFFF99", "#FDBF6F", "#FB9A99", "#6A3D9A", "#B15928", "#FF7F00", "#E31A1C"),
                    labels = c("Control - Severe stress", "Control - 12 h rewatering", "Control - 36 h rewatering", 
                               "Control - 60 h rewatering", "Severe stress", "12 h rewatering", "36 h rewatering", "60 h rewatering")) +
  scale_shape_manual(breaks = c("N22", "Dular", "Anjali"),
                     values = c(22, 24, 21),
                     labels = c("N22", "Dular", "Anjali")) +
  scale_x_continuous(limits = c(-8, 8), breaks = seq(-8, 8, by = 2)) +
  scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, by = 2)) +
  labs(x = "PC1 (43.1%)", y = "PC2 (24.5%)") +
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
        legend.spacing = unit(15.5, "lines"),
        legend.key.height = unit(1.3, "lines")) +
  guides(shape = guide_legend(order = 1, override.aes = list(size = 5)),
         fill = guide_legend(override.aes = list(size = 5),
                             ncol = 2),
         color = guide_legend(ncol = 2))



#multiplot - all organs
PCA_combined_leaf_spikelet_seed = cowplot::plot_grid(PCA_combined_flag_leaf_FL, NULL, PCA_combined_flag_leaf_EGF, 
                                                     NULL, NULL, NULL,
                                                     PCA_combined_flowering_spikelet, NULL, PCA_combined_developing_seed,
                                                     labels = c("A", "", "B", "", "", "", "C", "", "D"),
                                                     nrow = 3,
                                                     rel_widths = c(1, 0.03, 1),
                                                     rel_heights = c(1, 0.03, 1),
                                                     label_size = 20)
png("PCA_combined_leaf_spikelet_seed.png", width = 16*300, height = 16*300, res = 300)
print(PCA_combined_leaf_spikelet_seed)
dev.off()



###################################################
#FLAG LEAF#
###################################################


#comparison of EGF controls

#relative levels of metabolites under control conditions
#control_order_EGF already ran/called previously
data_combined_log2_mean_flag_leaf_EGF_controls = data_combined_log2_mean_EGF_flag_leaf_2[control_order_EGF, ]
data_combined_log2_mean_flag_leaf_EGF_controls = t(data_combined_log2_mean_flag_leaf_EGF_controls)


#heatmap including all metabolites
heatmap.2(data_combined_log2_mean_flag_leaf_EGF_controls, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(4, 20),  
          main = "Flag leaf-EGF - Combined \nRelative levels, Controls", 
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
shapiro_res_combined_N22_CLS_C12_flag_leaf_EGF = t(data.frame(lapply(data_combined_N22_log2median_EGF_flag_leaf[which(data_combined_N22_log2median_EGF_flag_leaf$Timepoint %in% c("Control.Late stress", "Control.12h Rewatering")), 3:ncol(data_combined_N22_log2median_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_CLS_C12_flag_leaf_EGF) = "p.value"
shapiro_res_combined_N22_CLS_C12_flag_leaf_EGF = as.data.frame(shapiro_res_combined_N22_CLS_C12_flag_leaf_EGF)
shapiro_res_combined_N22_CLS_C12_flag_leaf_EGF$Distribution = ifelse(shapiro_res_combined_N22_CLS_C12_flag_leaf_EGF$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_CLS_C12_flag_leaf_EGF$Distribution == "Normal")     #60 metabolites
sum(shapiro_res_combined_N22_CLS_C12_flag_leaf_EGF$Distribution == "Non-normal") #21 metabolites

#wilcox test
wilcox_res_N22_CLS_C12_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_N22_log2median_EGF_flag_leaf[ , 3:ncol(data_combined_N22_log2median_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_N22_log2median_EGF_flag_leaf$Timepoint, subset = data_combined_N22_log2median_EGF_flag_leaf$Timepoint %in% c("Control.Late stress", "Control.12h Rewatering"))$p.value)))
colnames(wilcox_res_N22_CLS_C12_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_N22_CLS_C12_flag_leaf_EGF_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_N22_CLS_C12_flag_leaf_EGF_combined = as.data.frame(wilcox_res_N22_CLS_C12_flag_leaf_EGF_combined)
wilcox_res_N22_CLS_C12_flag_leaf_EGF_combined$Significance = with(wilcox_res_N22_CLS_C12_flag_leaf_EGF_combined,
                                                                  ifelse(p.value <= 0.001, "***", ifelse(
                                                                    p.value <= 0.01, "**", ifelse(
                                                                      p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_CLS_C12_flag_leaf_EGF_combined$Significance == "*") #6 metabolites
sum(wilcox_res_N22_CLS_C12_flag_leaf_EGF_combined$Significance == "**") #3 metabolites
sum(wilcox_res_N22_CLS_C12_flag_leaf_EGF_combined$Significance == "***") #1 metabolites
sum(wilcox_res_N22_CLS_C12_flag_leaf_EGF_combined$Significance != "ns") #10 metabolites
sum(wilcox_res_N22_CLS_C12_flag_leaf_EGF_combined$Significance == "ns") #71 metabolites


#Dular

#data distribution
shapiro_res_combined_Dular_CLS_C12_flag_leaf_EGF = t(data.frame(lapply(data_combined_Dular_log2median_EGF_flag_leaf[which(data_combined_Dular_log2median_EGF_flag_leaf$Timepoint %in% c("Control.Late stress", "Control.12h Rewatering")), 3:ncol(data_combined_Dular_log2median_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_CLS_C12_flag_leaf_EGF) = "p.value"
shapiro_res_combined_Dular_CLS_C12_flag_leaf_EGF = as.data.frame(shapiro_res_combined_Dular_CLS_C12_flag_leaf_EGF)
shapiro_res_combined_Dular_CLS_C12_flag_leaf_EGF$Distribution = ifelse(shapiro_res_combined_Dular_CLS_C12_flag_leaf_EGF$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_CLS_C12_flag_leaf_EGF$Distribution == "Normal")     #60 metabolites
sum(shapiro_res_combined_Dular_CLS_C12_flag_leaf_EGF$Distribution == "Non-normal") #21 metabolites

#wilcox test
wilcox_res_Dular_CLS_C12_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_Dular_log2median_EGF_flag_leaf[ , 3:ncol(data_combined_Dular_log2median_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_EGF_flag_leaf$Timepoint, subset = data_combined_Dular_log2median_EGF_flag_leaf$Timepoint %in% c("Control.Late stress", "Control.12h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_CLS_C12_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_Dular_CLS_C12_flag_leaf_EGF_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_Dular_CLS_C12_flag_leaf_EGF_combined = as.data.frame(wilcox_res_Dular_CLS_C12_flag_leaf_EGF_combined)
wilcox_res_Dular_CLS_C12_flag_leaf_EGF_combined$Significance = with(wilcox_res_Dular_CLS_C12_flag_leaf_EGF_combined,
                                                                    ifelse(p.value <= 0.001, "***", ifelse(
                                                                      p.value <= 0.01, "**", ifelse(
                                                                        p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_CLS_C12_flag_leaf_EGF_combined$Significance == "*") #3 metabolites
sum(wilcox_res_Dular_CLS_C12_flag_leaf_EGF_combined$Significance == "**") #0 metabolites
sum(wilcox_res_Dular_CLS_C12_flag_leaf_EGF_combined$Significance == "***") #1 metabolites
sum(wilcox_res_Dular_CLS_C12_flag_leaf_EGF_combined$Significance != "ns") #4 metabolites
sum(wilcox_res_Dular_CLS_C12_flag_leaf_EGF_combined$Significance == "ns") #77 metabolites


#Anjali

#data distribution
shapiro_res_combined_Anjali_CLS_C12_flag_leaf_EGF = t(data.frame(lapply(data_combined_Anjali_log2median_EGF_flag_leaf[which(data_combined_Anjali_log2median_EGF_flag_leaf$Timepoint %in% c("Control.Late stress", "Control.12h Rewatering")), 3:ncol(data_combined_Anjali_log2median_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_CLS_C12_flag_leaf_EGF) = "p.value"
shapiro_res_combined_Anjali_CLS_C12_flag_leaf_EGF = as.data.frame(shapiro_res_combined_Anjali_CLS_C12_flag_leaf_EGF)
shapiro_res_combined_Anjali_CLS_C12_flag_leaf_EGF$Distribution = ifelse(shapiro_res_combined_Anjali_CLS_C12_flag_leaf_EGF$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_CLS_C12_flag_leaf_EGF$Distribution == "Normal")     #55 metabolites
sum(shapiro_res_combined_Anjali_CLS_C12_flag_leaf_EGF$Distribution == "Non-normal") #26 metabolites

#wilcox test
wilcox_res_Anjali_CLS_C12_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_Anjali_log2median_EGF_flag_leaf[ , 3:ncol(data_combined_Anjali_log2median_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_EGF_flag_leaf$Timepoint, subset = data_combined_Anjali_log2median_EGF_flag_leaf$Timepoint %in% c("Control.Late stress", "Control.12h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_CLS_C12_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_Anjali_CLS_C12_flag_leaf_EGF_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_Anjali_CLS_C12_flag_leaf_EGF_combined = as.data.frame(wilcox_res_Anjali_CLS_C12_flag_leaf_EGF_combined)
wilcox_res_Anjali_CLS_C12_flag_leaf_EGF_combined$Significance = with(wilcox_res_Anjali_CLS_C12_flag_leaf_EGF_combined,
                                                                     ifelse(p.value <= 0.001, "***", ifelse(
                                                                       p.value <= 0.01, "**", ifelse(
                                                                         p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_CLS_C12_flag_leaf_EGF_combined$Significance == "*") #4 metabolites
sum(wilcox_res_Anjali_CLS_C12_flag_leaf_EGF_combined$Significance == "**") #1 metabolite
sum(wilcox_res_Anjali_CLS_C12_flag_leaf_EGF_combined$Significance == "***") #0 metabolite
sum(wilcox_res_Anjali_CLS_C12_flag_leaf_EGF_combined$Significance != "ns") #5 metabolites
sum(wilcox_res_Anjali_CLS_C12_flag_leaf_EGF_combined$Significance == "ns") #76 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_CLS_C12_flag_leaf_EGF_combined = cbind(wilcox_res_N22_CLS_C12_flag_leaf_EGF_combined, 
                                                  wilcox_res_Dular_CLS_C12_flag_leaf_EGF_combined$Significance, 
                                                  wilcox_res_Anjali_CLS_C12_flag_leaf_EGF_combined$Significance)
wilcox_sig_CLS_C12_flag_leaf_EGF_combined = wilcox_sig_CLS_C12_flag_leaf_EGF_combined[ , -1]
colnames(wilcox_sig_CLS_C12_flag_leaf_EGF_combined)[1] = "wilcox_res_N22_CLS_C12_flag_leaf_EGF_combined$Significance"
wilcox_sig_CLS_C12_flag_leaf_EGF_combined = data.frame(lapply(wilcox_sig_CLS_C12_flag_leaf_EGF_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_CLS_C12_flag_leaf_EGF_combined) = rownames(wilcox_res_N22_CLS_C12_flag_leaf_EGF_combined)
wilcox_sig_CLS_C12_flag_leaf_EGF_combined$Nonsig = rowSums(wilcox_sig_CLS_C12_flag_leaf_EGF_combined == " ")  #count number of nonsignificant
wilcox_sig_final_CLS_C12_flag_leaf_EGF_combined = wilcox_sig_CLS_C12_flag_leaf_EGF_combined[-which(wilcox_sig_CLS_C12_flag_leaf_EGF_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 14 metabolites

wilcox_sig_log2mean_CLS_C12_flag_leaf_EGF_combined = data_combined_log2_mean_flag_leaf_EGF_controls[which(rownames(data_combined_log2_mean_flag_leaf_EGF_controls) %in% rownames(wilcox_sig_final_CLS_C12_flag_leaf_EGF_combined)), ]  #relative control levels of metabolites significant in at least one of the comparisons

#overview
heatmap.2(wilcox_sig_log2mean_CLS_C12_flag_leaf_EGF_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)),
          margins = c(5, 20),  
          main = "Flag leaf-EGF-Combined-Wilcoxon \nRelative levels, Sig in CLS vs C12", 
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
          rowsep = c(seq(0, 14)), 
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))


####################

#CLS vs. C36

#N22

#data distribution
shapiro_res_combined_N22_CLS_C36_flag_leaf_EGF = t(data.frame(lapply(data_combined_N22_log2median_EGF_flag_leaf[which(data_combined_N22_log2median_EGF_flag_leaf$Timepoint %in% c("Control.Late stress", "Control.36h Rewatering")), 3:ncol(data_combined_N22_log2median_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_CLS_C36_flag_leaf_EGF) = "p.value"
shapiro_res_combined_N22_CLS_C36_flag_leaf_EGF = as.data.frame(shapiro_res_combined_N22_CLS_C36_flag_leaf_EGF)
shapiro_res_combined_N22_CLS_C36_flag_leaf_EGF$Distribution = ifelse(shapiro_res_combined_N22_CLS_C36_flag_leaf_EGF$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_CLS_C36_flag_leaf_EGF$Distribution == "Normal")     #63 metabolites
sum(shapiro_res_combined_N22_CLS_C36_flag_leaf_EGF$Distribution == "Non-normal") #18 metabolites

#wilcox test
wilcox_res_N22_CLS_C36_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_N22_log2median_EGF_flag_leaf[ , 3:ncol(data_combined_N22_log2median_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_N22_log2median_EGF_flag_leaf$Timepoint, subset = data_combined_N22_log2median_EGF_flag_leaf$Timepoint %in% c("Control.Late stress", "Control.36h Rewatering"))$p.value)))
colnames(wilcox_res_N22_CLS_C36_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_N22_CLS_C36_flag_leaf_EGF_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_N22_CLS_C36_flag_leaf_EGF_combined = as.data.frame(wilcox_res_N22_CLS_C36_flag_leaf_EGF_combined)
wilcox_res_N22_CLS_C36_flag_leaf_EGF_combined$Significance = with(wilcox_res_N22_CLS_C36_flag_leaf_EGF_combined,
                                                                  ifelse(p.value <= 0.001, "***", ifelse(
                                                                    p.value <= 0.01, "**", ifelse(
                                                                      p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_CLS_C36_flag_leaf_EGF_combined$Significance == "*") #4 metabolites
sum(wilcox_res_N22_CLS_C36_flag_leaf_EGF_combined$Significance == "**") #5 metabolites
sum(wilcox_res_N22_CLS_C36_flag_leaf_EGF_combined$Significance == "***") #0 metabolite
sum(wilcox_res_N22_CLS_C36_flag_leaf_EGF_combined$Significance != "ns") #9 metabolites
sum(wilcox_res_N22_CLS_C36_flag_leaf_EGF_combined$Significance == "ns") #72 metabolites


#Dular

#data distribution
shapiro_res_combined_Dular_CLS_C36_flag_leaf_EGF = t(data.frame(lapply(data_combined_Dular_log2median_EGF_flag_leaf[which(data_combined_Dular_log2median_EGF_flag_leaf$Timepoint %in% c("Control.Late stress", "Control.36h Rewatering")), 3:ncol(data_combined_Dular_log2median_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_CLS_C36_flag_leaf_EGF) = "p.value"
shapiro_res_combined_Dular_CLS_C36_flag_leaf_EGF = as.data.frame(shapiro_res_combined_Dular_CLS_C36_flag_leaf_EGF)
shapiro_res_combined_Dular_CLS_C36_flag_leaf_EGF$Distribution = ifelse(shapiro_res_combined_Dular_CLS_C36_flag_leaf_EGF$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_CLS_C36_flag_leaf_EGF$Distribution == "Normal")     #60 metabolites
sum(shapiro_res_combined_Dular_CLS_C36_flag_leaf_EGF$Distribution == "Non-normal") #21 metabolites

#wilcox test
wilcox_res_Dular_CLS_C36_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_Dular_log2median_EGF_flag_leaf[ , 3:ncol(data_combined_Dular_log2median_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_EGF_flag_leaf$Timepoint, subset = data_combined_Dular_log2median_EGF_flag_leaf$Timepoint %in% c("Control.Late stress", "Control.36h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_CLS_C36_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_Dular_CLS_C36_flag_leaf_EGF_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_Dular_CLS_C36_flag_leaf_EGF_combined = as.data.frame(wilcox_res_Dular_CLS_C36_flag_leaf_EGF_combined)
wilcox_res_Dular_CLS_C36_flag_leaf_EGF_combined$Significance = with(wilcox_res_Dular_CLS_C36_flag_leaf_EGF_combined,
                                                                    ifelse(p.value <= 0.001, "***", ifelse(
                                                                      p.value <= 0.01, "**", ifelse(
                                                                        p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_CLS_C36_flag_leaf_EGF_combined$Significance == "*") #7 metabolites
sum(wilcox_res_Dular_CLS_C36_flag_leaf_EGF_combined$Significance == "**") #3 metabolites
sum(wilcox_res_Dular_CLS_C36_flag_leaf_EGF_combined$Significance == "***") #2 metabolites
sum(wilcox_res_Dular_CLS_C36_flag_leaf_EGF_combined$Significance != "ns") #12 metabolites
sum(wilcox_res_Dular_CLS_C36_flag_leaf_EGF_combined$Significance == "ns") #69 metabolites


#Anjali

#data distribution
shapiro_res_combined_Anjali_CLS_C36_flag_leaf_EGF = t(data.frame(lapply(data_combined_Anjali_log2median_EGF_flag_leaf[which(data_combined_Anjali_log2median_EGF_flag_leaf$Timepoint %in% c("Control.Late stress", "Control.36h Rewatering")), 3:ncol(data_combined_Anjali_log2median_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_CLS_C36_flag_leaf_EGF) = "p.value"
shapiro_res_combined_Anjali_CLS_C36_flag_leaf_EGF = as.data.frame(shapiro_res_combined_Anjali_CLS_C36_flag_leaf_EGF)
shapiro_res_combined_Anjali_CLS_C36_flag_leaf_EGF$Distribution = ifelse(shapiro_res_combined_Anjali_CLS_C36_flag_leaf_EGF$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_CLS_C36_flag_leaf_EGF$Distribution == "Normal")     #64 metabolites
sum(shapiro_res_combined_Anjali_CLS_C36_flag_leaf_EGF$Distribution == "Non-normal") #17 metabolites

#wilcox test
wilcox_res_Anjali_CLS_C36_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_Anjali_log2median_EGF_flag_leaf[ , 3:ncol(data_combined_Anjali_log2median_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_EGF_flag_leaf$Timepoint, subset = data_combined_Anjali_log2median_EGF_flag_leaf$Timepoint %in% c("Control.Late stress", "Control.36h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_CLS_C36_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_Anjali_CLS_C36_flag_leaf_EGF_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_Anjali_CLS_C36_flag_leaf_EGF_combined = as.data.frame(wilcox_res_Anjali_CLS_C36_flag_leaf_EGF_combined)
wilcox_res_Anjali_CLS_C36_flag_leaf_EGF_combined$Significance = with(wilcox_res_Anjali_CLS_C36_flag_leaf_EGF_combined,
                                                                     ifelse(p.value <= 0.001, "***", ifelse(
                                                                       p.value <= 0.01, "**", ifelse(
                                                                         p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_CLS_C36_flag_leaf_EGF_combined$Significance == "*") #4 metabolites
sum(wilcox_res_Anjali_CLS_C36_flag_leaf_EGF_combined$Significance == "**") #5 metabolite
sum(wilcox_res_Anjali_CLS_C36_flag_leaf_EGF_combined$Significance == "***") #0 metabolite
sum(wilcox_res_Anjali_CLS_C36_flag_leaf_EGF_combined$Significance != "ns") #9 metabolites
sum(wilcox_res_Anjali_CLS_C36_flag_leaf_EGF_combined$Significance == "ns") #72 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_CLS_C36_flag_leaf_EGF_combined = cbind(wilcox_res_N22_CLS_C36_flag_leaf_EGF_combined, 
                                                  wilcox_res_Dular_CLS_C36_flag_leaf_EGF_combined$Significance, 
                                                  wilcox_res_Anjali_CLS_C36_flag_leaf_EGF_combined$Significance)
wilcox_sig_CLS_C36_flag_leaf_EGF_combined = wilcox_sig_CLS_C36_flag_leaf_EGF_combined[ , -1]
colnames(wilcox_sig_CLS_C36_flag_leaf_EGF_combined)[1] = "wilcox_res_N22_CLS_C36_flag_leaf_EGF_combined$Significance"
wilcox_sig_CLS_C36_flag_leaf_EGF_combined = data.frame(lapply(wilcox_sig_CLS_C36_flag_leaf_EGF_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_CLS_C36_flag_leaf_EGF_combined) = rownames(wilcox_res_N22_CLS_C36_flag_leaf_EGF_combined)
wilcox_sig_CLS_C36_flag_leaf_EGF_combined$Nonsig = rowSums(wilcox_sig_CLS_C36_flag_leaf_EGF_combined == " ")  #count number of nonsignificant
wilcox_sig_final_CLS_C36_flag_leaf_EGF_combined = wilcox_sig_CLS_C36_flag_leaf_EGF_combined[-which(wilcox_sig_CLS_C36_flag_leaf_EGF_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 23 metabolites

wilcox_sig_log2mean_CLS_C36_flag_leaf_EGF_combined = data_combined_log2_mean_flag_leaf_EGF_controls[which(rownames(data_combined_log2_mean_flag_leaf_EGF_controls) %in% rownames(wilcox_sig_final_CLS_C36_flag_leaf_EGF_combined)), ]  #relative control levels of metabolites significant in at least one of the comparisons

#overview
heatmap.2(wilcox_sig_log2mean_CLS_C36_flag_leaf_EGF_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)),
          margins = c(5, 20), 
          main = "Flag leaf-EGF-Combined-Wilcoxon \nRelative levels, Sig in CLS vs C36", 
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
          rowsep = c(seq(0, 23)), 
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))


####################

#CLS vs. C60

#N22

#data distribution
shapiro_res_combined_N22_CLS_C60_flag_leaf_EGF = t(data.frame(lapply(data_combined_N22_log2median_EGF_flag_leaf[which(data_combined_N22_log2median_EGF_flag_leaf$Timepoint %in% c("Control.Late stress", "Control.60h Rewatering")), 3:ncol(data_combined_N22_log2median_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_CLS_C60_flag_leaf_EGF) = "p.value"
shapiro_res_combined_N22_CLS_C60_flag_leaf_EGF = as.data.frame(shapiro_res_combined_N22_CLS_C60_flag_leaf_EGF)
shapiro_res_combined_N22_CLS_C60_flag_leaf_EGF$Distribution = ifelse(shapiro_res_combined_N22_CLS_C60_flag_leaf_EGF$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_CLS_C60_flag_leaf_EGF$Distribution == "Normal")     #58 metabolites
sum(shapiro_res_combined_N22_CLS_C60_flag_leaf_EGF$Distribution == "Non-normal") #23 metabolites

#wilcox test
wilcox_res_N22_CLS_C60_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_N22_log2median_EGF_flag_leaf[ , 3:ncol(data_combined_N22_log2median_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_N22_log2median_EGF_flag_leaf$Timepoint, subset = data_combined_N22_log2median_EGF_flag_leaf$Timepoint %in% c("Control.Late stress", "Control.60h Rewatering"))$p.value)))
colnames(wilcox_res_N22_CLS_C60_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_N22_CLS_C60_flag_leaf_EGF_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_N22_CLS_C60_flag_leaf_EGF_combined = as.data.frame(wilcox_res_N22_CLS_C60_flag_leaf_EGF_combined)
wilcox_res_N22_CLS_C60_flag_leaf_EGF_combined$Significance = with(wilcox_res_N22_CLS_C60_flag_leaf_EGF_combined,
                                                                  ifelse(p.value <= 0.001, "***", ifelse(
                                                                    p.value <= 0.01, "**", ifelse(
                                                                      p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_CLS_C60_flag_leaf_EGF_combined$Significance == "*") #1 metabolite
sum(wilcox_res_N22_CLS_C60_flag_leaf_EGF_combined$Significance == "**") #2 metabolites
sum(wilcox_res_N22_CLS_C60_flag_leaf_EGF_combined$Significance == "***") #2 metabolites
sum(wilcox_res_N22_CLS_C60_flag_leaf_EGF_combined$Significance != "ns") #5 metabolites
sum(wilcox_res_N22_CLS_C60_flag_leaf_EGF_combined$Significance == "ns") #76 metabolites


#Dular

#data distribution
shapiro_res_combined_Dular_CLS_C60_flag_leaf_EGF = t(data.frame(lapply(data_combined_Dular_log2median_EGF_flag_leaf[which(data_combined_Dular_log2median_EGF_flag_leaf$Timepoint %in% c("Control.Late stress", "Control.60h Rewatering")), 3:ncol(data_combined_Dular_log2median_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_CLS_C60_flag_leaf_EGF) = "p.value"
shapiro_res_combined_Dular_CLS_C60_flag_leaf_EGF = as.data.frame(shapiro_res_combined_Dular_CLS_C60_flag_leaf_EGF)
shapiro_res_combined_Dular_CLS_C60_flag_leaf_EGF$Distribution = ifelse(shapiro_res_combined_Dular_CLS_C60_flag_leaf_EGF$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_CLS_C60_flag_leaf_EGF$Distribution == "Normal")     #57 metabolites
sum(shapiro_res_combined_Dular_CLS_C60_flag_leaf_EGF$Distribution == "Non-normal") #24 metabolites

#wilcox test
wilcox_res_Dular_CLS_C60_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_Dular_log2median_EGF_flag_leaf[ , 3:ncol(data_combined_Dular_log2median_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_EGF_flag_leaf$Timepoint, subset = data_combined_Dular_log2median_EGF_flag_leaf$Timepoint %in% c("Control.Late stress", "Control.60h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_CLS_C60_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_Dular_CLS_C60_flag_leaf_EGF_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_Dular_CLS_C60_flag_leaf_EGF_combined = as.data.frame(wilcox_res_Dular_CLS_C60_flag_leaf_EGF_combined)
wilcox_res_Dular_CLS_C60_flag_leaf_EGF_combined$Significance = with(wilcox_res_Dular_CLS_C60_flag_leaf_EGF_combined,
                                                                    ifelse(p.value <= 0.001, "***", ifelse(
                                                                      p.value <= 0.01, "**", ifelse(
                                                                        p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_CLS_C60_flag_leaf_EGF_combined$Significance == "*") #8 metabolites
sum(wilcox_res_Dular_CLS_C60_flag_leaf_EGF_combined$Significance == "**") #5 metabolites
sum(wilcox_res_Dular_CLS_C60_flag_leaf_EGF_combined$Significance == "***") #3 metabolites
sum(wilcox_res_Dular_CLS_C60_flag_leaf_EGF_combined$Significance != "ns") #16 metabolites
sum(wilcox_res_Dular_CLS_C60_flag_leaf_EGF_combined$Significance == "ns") #65 metabolites


#Anjali

#data distribution
shapiro_res_combined_Anjali_CLS_C60_flag_leaf_EGF = t(data.frame(lapply(data_combined_Anjali_log2median_EGF_flag_leaf[which(data_combined_Anjali_log2median_EGF_flag_leaf$Timepoint %in% c("Control.Late stress", "Control.60h Rewatering")), 3:ncol(data_combined_Anjali_log2median_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_CLS_C60_flag_leaf_EGF) = "p.value"
shapiro_res_combined_Anjali_CLS_C60_flag_leaf_EGF = as.data.frame(shapiro_res_combined_Anjali_CLS_C60_flag_leaf_EGF)
shapiro_res_combined_Anjali_CLS_C60_flag_leaf_EGF$Distribution = ifelse(shapiro_res_combined_Anjali_CLS_C60_flag_leaf_EGF$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_CLS_C60_flag_leaf_EGF$Distribution == "Normal")     #49 metabolites
sum(shapiro_res_combined_Anjali_CLS_C60_flag_leaf_EGF$Distribution == "Non-normal") #32 metabolites

#wilcox test
wilcox_res_Anjali_CLS_C60_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_Anjali_log2median_EGF_flag_leaf[ , 3:ncol(data_combined_Anjali_log2median_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_EGF_flag_leaf$Timepoint, subset = data_combined_Anjali_log2median_EGF_flag_leaf$Timepoint %in% c("Control.Late stress", "Control.60h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_CLS_C60_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_Anjali_CLS_C60_flag_leaf_EGF_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_Anjali_CLS_C60_flag_leaf_EGF_combined = as.data.frame(wilcox_res_Anjali_CLS_C60_flag_leaf_EGF_combined)
wilcox_res_Anjali_CLS_C60_flag_leaf_EGF_combined$Significance = with(wilcox_res_Anjali_CLS_C60_flag_leaf_EGF_combined,
                                                                     ifelse(p.value <= 0.001, "***", ifelse(
                                                                       p.value <= 0.01, "**", ifelse(
                                                                         p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_CLS_C60_flag_leaf_EGF_combined$Significance == "*") #10 metabolites
sum(wilcox_res_Anjali_CLS_C60_flag_leaf_EGF_combined$Significance == "**") #9 metabolites
sum(wilcox_res_Anjali_CLS_C60_flag_leaf_EGF_combined$Significance == "***") #0 metabolite
sum(wilcox_res_Anjali_CLS_C60_flag_leaf_EGF_combined$Significance != "ns") #19 metabolites
sum(wilcox_res_Anjali_CLS_C60_flag_leaf_EGF_combined$Significance == "ns") #62 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_CLS_C60_flag_leaf_EGF_combined = cbind(wilcox_res_N22_CLS_C60_flag_leaf_EGF_combined, 
                                                  wilcox_res_Dular_CLS_C60_flag_leaf_EGF_combined$Significance, 
                                                  wilcox_res_Anjali_CLS_C60_flag_leaf_EGF_combined$Significance)
wilcox_sig_CLS_C60_flag_leaf_EGF_combined = wilcox_sig_CLS_C60_flag_leaf_EGF_combined[ , -1]
colnames(wilcox_sig_CLS_C60_flag_leaf_EGF_combined)[1] = "wilcox_res_N22_CLS_C60_flag_leaf_EGF_combined$Significance"
wilcox_sig_CLS_C60_flag_leaf_EGF_combined = data.frame(lapply(wilcox_sig_CLS_C60_flag_leaf_EGF_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_CLS_C60_flag_leaf_EGF_combined) = rownames(wilcox_res_N22_CLS_C60_flag_leaf_EGF_combined)
wilcox_sig_CLS_C60_flag_leaf_EGF_combined$Nonsig = rowSums(wilcox_sig_CLS_C60_flag_leaf_EGF_combined == " ")  #count number of nonsignificant
wilcox_sig_final_CLS_C60_flag_leaf_EGF_combined = wilcox_sig_CLS_C60_flag_leaf_EGF_combined[-which(wilcox_sig_CLS_C60_flag_leaf_EGF_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 30 metabolites

wilcox_sig_log2mean_CLS_C60_flag_leaf_EGF_combined = data_combined_log2_mean_flag_leaf_EGF_controls[which(rownames(data_combined_log2_mean_flag_leaf_EGF_controls) %in% rownames(wilcox_sig_final_CLS_C60_flag_leaf_EGF_combined)), ]  #relative control levels of metabolites significant in at least one of the comparisons

#overview
heatmap.2(wilcox_sig_log2mean_CLS_C60_flag_leaf_EGF_combined,
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(5, 20),  
          main = "Flag leaf-EGF-Combined-Wilcoxon \nRelative levels, Sig in CLS vs C60", 
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
          rowsep = c(seq(0, 30)),
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))


####################

#significance of wilcoxon test - all cultivars, all four control timepoints
wilcox_sig_CLS_C12_C36_C60_flag_leaf_EGF_combined = cbind(wilcox_res_N22_CLS_C12_flag_leaf_EGF_combined,
                                                          wilcox_res_N22_CLS_C36_flag_leaf_EGF_combined$Significance,
                                                          wilcox_res_N22_CLS_C60_flag_leaf_EGF_combined$Significance,
                                                          wilcox_res_Dular_CLS_C12_flag_leaf_EGF_combined$Significance,
                                                          wilcox_res_Dular_CLS_C36_flag_leaf_EGF_combined$Significance,
                                                          wilcox_res_Dular_CLS_C60_flag_leaf_EGF_combined$Significance,
                                                          wilcox_res_Anjali_CLS_C12_flag_leaf_EGF_combined$Significance,
                                                          wilcox_res_Anjali_CLS_C36_flag_leaf_EGF_combined$Significance,
                                                          wilcox_res_Anjali_CLS_C60_flag_leaf_EGF_combined$Significance)
wilcox_sig_CLS_C12_C36_C60_flag_leaf_EGF_combined = wilcox_sig_CLS_C12_C36_C60_flag_leaf_EGF_combined[ , -1]
colnames(wilcox_sig_CLS_C12_C36_C60_flag_leaf_EGF_combined)[1] = "wilcox_res_N22_CLS_C12_flag_leaf_EGF_combined$Significance"
wilcox_sig_CLS_C12_C36_C60_flag_leaf_EGF_combined = data.frame(lapply(wilcox_sig_CLS_C12_C36_C60_flag_leaf_EGF_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_CLS_C12_C36_C60_flag_leaf_EGF_combined) = rownames(wilcox_res_N22_CLS_C12_flag_leaf_EGF_combined)
wilcox_sig_CLS_C12_C36_C60_flag_leaf_EGF_combined$Nonsig = rowSums(wilcox_sig_CLS_C12_C36_C60_flag_leaf_EGF_combined == " ")  #count number of nonsignificant
wilcox_sig_final_CLS_C12_C36_C60_flag_leaf_EGF_combined = wilcox_sig_CLS_C12_C36_C60_flag_leaf_EGF_combined[-which(wilcox_sig_CLS_C12_C36_C60_flag_leaf_EGF_combined$Nonsig == "9"), ] #metabolites significant in at least one of the comparisons - 45 metabolites

wilcox_sig_log2mean_CLS_C12_C36_C60_flag_leaf_EGF_combined = data_combined_log2_mean_flag_leaf_EGF_controls[which(rownames(data_combined_log2_mean_flag_leaf_EGF_controls) %in% rownames(wilcox_sig_final_CLS_C12_C36_C60_flag_leaf_EGF_combined)), ]  #relative control levels of metabolites significant in at least one of the comparisons


#significance as cell notes
wilcox_sig_final_CLS_C12_C36_C60_flag_leaf_EGF_combined_2 = wilcox_sig_final_CLS_C12_C36_C60_flag_leaf_EGF_combined[ , -10]
wilcox_sig_final_CLS_C12_C36_C60_flag_leaf_EGF_combined_2$N22_Late_stress = ""
wilcox_sig_final_CLS_C12_C36_C60_flag_leaf_EGF_combined_2$Dular_Late_stress = ""
wilcox_sig_final_CLS_C12_C36_C60_flag_leaf_EGF_combined_2$Anjali_Late_stress = ""
wilcox_sig_final_CLS_C12_C36_C60_flag_leaf_EGF_combined_2 = wilcox_sig_final_CLS_C12_C36_C60_flag_leaf_EGF_combined_2[ , c(10, 1, 2, 3, 11, 4, 5, 6, 12, 7, 8, 9)]


#visualization
png("wilcox_sig_log2mean_controls_flag_leaf_EGF_combined.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_log2mean_CLS_C12_C36_C60_flag_leaf_EGF_combined, 
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
          cexRow = 1.2, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 0.2, 11), 
          lwid = c(0.8, 5, 1.5), 
          srtCol = 45, 
          cellnote = wilcox_sig_final_CLS_C12_C36_C60_flag_leaf_EGF_combined_2, 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 12)), 
          rowsep = c(seq(0, 45)), 
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


#late stress vs. recovery - NOTE: consider only the stress-responsive metabolites (i.e. metabolites which have significant difference in LS/C comparison)

#Flowering stage

#late stress vs. recovery

#log2-fold change - late stress vs. recovery
latestress_recovery_order = c("FNLS", "FNR12", "FNR36", "FNR60", "FDLS", "FDR12", "FDR36", "FDR60", "FALS", "FAR12", "FAR36", "FAR60")
data_combined_latestress_recovery_log2_mean_FL_flag_leaf = data_combined_log2_mean_FL_flag_leaf_2[latestress_recovery_order, ]  #late stress and recovery data in preferred order
data_combined_latestress_recovery_log2_mean_FL_flag_leaf = t(data_combined_latestress_recovery_log2_mean_FL_flag_leaf)

log2fc_combined_N22_latestress_recovery_FL_flag_leaf = data_combined_latestress_recovery_log2_mean_FL_flag_leaf[ , 2:4] - data_combined_latestress_recovery_log2_mean_FL_flag_leaf[ , 1]
log2fc_combined_Dular_latestress_recovery_FL_flag_leaf = data_combined_latestress_recovery_log2_mean_FL_flag_leaf[ , 6:8] - data_combined_latestress_recovery_log2_mean_FL_flag_leaf[ , 5]
log2fc_combined_Anjali_latestress_recovery_FL_flag_leaf = data_combined_latestress_recovery_log2_mean_FL_flag_leaf[ , 10:12] - data_combined_latestress_recovery_log2_mean_FL_flag_leaf[ , 9]

log2fc_combined_latestress_recovery_FL_flag_leaf = cbind(log2fc_combined_N22_latestress_recovery_FL_flag_leaf, 
                                                         log2fc_combined_Dular_latestress_recovery_FL_flag_leaf, 
                                                         log2fc_combined_Anjali_latestress_recovery_FL_flag_leaf)   #combine log2-fc of all cultivars in one dataset
colnames(log2fc_combined_latestress_recovery_FL_flag_leaf)[1:3] = c(paste(colnames(log2fc_combined_latestress_recovery_FL_flag_leaf)[1:3], "- FNLS"))
colnames(log2fc_combined_latestress_recovery_FL_flag_leaf)[4:6] = c(paste(colnames(log2fc_combined_latestress_recovery_FL_flag_leaf)[4:6], "- FDLS"))
colnames(log2fc_combined_latestress_recovery_FL_flag_leaf)[7:9] = c(paste(colnames(log2fc_combined_latestress_recovery_FL_flag_leaf)[7:9], "- FALS"))

#overview
heatmap.2(log2fc_combined_latestress_recovery_FL_flag_leaf, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(7, 20),  
          main = "Flag leaf_FL - Combined - Log2 FC, Recovery/Late Stress", 
          cexCol = 1, 
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          srtCol = 90)


#late stress vs. 12h rewatering

#N22

#data distribution 
data_combined_N22_log2median_stress_responsive_FL_flag_leaf = data_combined_N22_log2median_FL_flag_leaf
colnames(data_combined_N22_log2median_stress_responsive_FL_flag_leaf)[3:83] = reduced_metabolite_list_final_flag_leaf_2013$Name
#only metabolites which have significant difference in LS/C comparison
data_combined_N22_log2median_stress_responsive_FL_flag_leaf = data_combined_N22_log2median_stress_responsive_FL_flag_leaf[ , c(1, 2, (which(colnames(data_combined_N22_log2median_stress_responsive_FL_flag_leaf) %in% rownames(wilcox_sig_log2fc_C_LS_flag_leaf_FL_combined))))]

shapiro_res_combined_N22_latestress_recovery12_FL_flag_leaf_2 = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_FL_flag_leaf[which(data_combined_N22_log2median_stress_responsive_FL_flag_leaf$Timepoint %in% c("Late stress", "12h Rewatering")), 3:ncol(data_combined_N22_log2median_stress_responsive_FL_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_latestress_recovery12_FL_flag_leaf_2) = "p.value"
shapiro_res_combined_N22_latestress_recovery12_FL_flag_leaf_2 = as.data.frame(shapiro_res_combined_N22_latestress_recovery12_FL_flag_leaf_2)
shapiro_res_combined_N22_latestress_recovery12_FL_flag_leaf_2$Distribution = ifelse(shapiro_res_combined_N22_latestress_recovery12_FL_flag_leaf_2$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_latestress_recovery12_FL_flag_leaf_2$Distribution == "Normal")     #28 metabolites
sum(shapiro_res_combined_N22_latestress_recovery12_FL_flag_leaf_2$Distribution == "Non-normal") #27 metabolites

#wilcox test
wilcox_res_N22_LS_R12_flag_leaf_FL_combined_2 = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_FL_flag_leaf[ , 3:ncol(data_combined_N22_log2median_stress_responsive_FL_flag_leaf)], function(x) wilcox.test(x ~ data_combined_N22_log2median_stress_responsive_FL_flag_leaf$Timepoint, subset = data_combined_N22_log2median_stress_responsive_FL_flag_leaf$Timepoint %in% c("Late stress", "12h Rewatering"))$p.value)))
colnames(wilcox_res_N22_LS_R12_flag_leaf_FL_combined_2) = "p.value"
rownames(wilcox_res_N22_LS_R12_flag_leaf_FL_combined_2) = rownames(wilcox_sig_log2fc_C_LS_flag_leaf_FL_combined)
wilcox_res_N22_LS_R12_flag_leaf_FL_combined_2 = as.data.frame(wilcox_res_N22_LS_R12_flag_leaf_FL_combined_2)
wilcox_res_N22_LS_R12_flag_leaf_FL_combined_2$Significance = with(wilcox_res_N22_LS_R12_flag_leaf_FL_combined_2,
                                                                  ifelse(p.value <= 0.001, "***", ifelse(
                                                                    p.value <= 0.01, "**", ifelse(
                                                                      p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_LS_R12_flag_leaf_FL_combined_2$Significance == "*") #1 metabolites
sum(wilcox_res_N22_LS_R12_flag_leaf_FL_combined_2$Significance == "**") #2 metabolites
sum(wilcox_res_N22_LS_R12_flag_leaf_FL_combined_2$Significance == "***") #1 metabolites
sum(wilcox_res_N22_LS_R12_flag_leaf_FL_combined_2$Significance != "ns") #4 metabolites
sum(wilcox_res_N22_LS_R12_flag_leaf_FL_combined_2$Significance == "ns") #51 metabolites 


#Dular

#data distribution 
data_combined_Dular_log2median_stress_responsive_FL_flag_leaf = data_combined_Dular_log2median_FL_flag_leaf
colnames(data_combined_Dular_log2median_stress_responsive_FL_flag_leaf)[3:83] = reduced_metabolite_list_final_flag_leaf_2013$Name
#only metabolites which have significant difference in LS/C comparison
data_combined_Dular_log2median_stress_responsive_FL_flag_leaf = data_combined_Dular_log2median_stress_responsive_FL_flag_leaf[ , c(1, 2, (which(colnames(data_combined_Dular_log2median_stress_responsive_FL_flag_leaf) %in% rownames(wilcox_sig_log2fc_C_LS_flag_leaf_FL_combined))))]

shapiro_res_combined_Dular_latestress_recovery12_FL_flag_leaf_2 = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_FL_flag_leaf[which(data_combined_Dular_log2median_stress_responsive_FL_flag_leaf$Timepoint %in% c("Late stress", "12h Rewatering")), 3:ncol(data_combined_Dular_log2median_stress_responsive_FL_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_latestress_recovery12_FL_flag_leaf_2) = "p.value"
shapiro_res_combined_Dular_latestress_recovery12_FL_flag_leaf_2 = as.data.frame(shapiro_res_combined_Dular_latestress_recovery12_FL_flag_leaf_2)
shapiro_res_combined_Dular_latestress_recovery12_FL_flag_leaf_2$Distribution = ifelse(shapiro_res_combined_Dular_latestress_recovery12_FL_flag_leaf_2$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_latestress_recovery12_FL_flag_leaf_2$Distribution == "Normal")     #39 metabolites
sum(shapiro_res_combined_Dular_latestress_recovery12_FL_flag_leaf_2$Distribution == "Non-normal") #16 metabolites

#wilcox test
wilcox_res_Dular_LS_R12_flag_leaf_FL_combined_2 = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_FL_flag_leaf[ , 3:ncol(data_combined_Dular_log2median_stress_responsive_FL_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_stress_responsive_FL_flag_leaf$Timepoint, subset = data_combined_Dular_log2median_stress_responsive_FL_flag_leaf$Timepoint %in% c("Late stress", "12h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_LS_R12_flag_leaf_FL_combined_2) = "p.value"
rownames(wilcox_res_Dular_LS_R12_flag_leaf_FL_combined_2) = rownames(wilcox_sig_log2fc_C_LS_flag_leaf_FL_combined)
wilcox_res_Dular_LS_R12_flag_leaf_FL_combined_2 = as.data.frame(wilcox_res_Dular_LS_R12_flag_leaf_FL_combined_2)
wilcox_res_Dular_LS_R12_flag_leaf_FL_combined_2$Significance = with(wilcox_res_Dular_LS_R12_flag_leaf_FL_combined_2,
                                                                    ifelse(p.value <= 0.001, "***", ifelse(
                                                                      p.value <= 0.01, "**", ifelse(
                                                                        p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_LS_R12_flag_leaf_FL_combined_2$Significance == "*") #3 metabolites
sum(wilcox_res_Dular_LS_R12_flag_leaf_FL_combined_2$Significance == "**") #0 metabolites
sum(wilcox_res_Dular_LS_R12_flag_leaf_FL_combined_2$Significance == "***") #0 metabolites
sum(wilcox_res_Dular_LS_R12_flag_leaf_FL_combined_2$Significance != "ns") #3 metabolites
sum(wilcox_res_Dular_LS_R12_flag_leaf_FL_combined_2$Significance == "ns") #52 metabolites


#Anjali

#data distribution 
data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf = data_combined_Anjali_log2median_FL_flag_leaf
colnames(data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf)[3:83] = reduced_metabolite_list_final_flag_leaf_2013$Name
#only metabolites which have significant difference in LS/C comparison
data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf = data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf[ , c(1, 2, (which(colnames(data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf) %in% rownames(wilcox_sig_log2fc_C_LS_flag_leaf_FL_combined))))]

shapiro_res_combined_Anjali_latestress_recovery12_FL_flag_leaf_2 = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf[which(data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf$Timepoint %in% c("Late stress", "12h Rewatering")), 3:ncol(data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_latestress_recovery12_FL_flag_leaf_2) = "p.value"
shapiro_res_combined_Anjali_latestress_recovery12_FL_flag_leaf_2 = as.data.frame(shapiro_res_combined_Anjali_latestress_recovery12_FL_flag_leaf_2)
shapiro_res_combined_Anjali_latestress_recovery12_FL_flag_leaf_2$Distribution = ifelse(shapiro_res_combined_Anjali_latestress_recovery12_FL_flag_leaf_2$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_latestress_recovery12_FL_flag_leaf_2$Distribution == "Normal")     #36 metabolites
sum(shapiro_res_combined_Anjali_latestress_recovery12_FL_flag_leaf_2$Distribution == "Non-normal") #19 metabolites

#wilcox test
wilcox_res_Anjali_LS_R12_flag_leaf_FL_combined_2 = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf[ , 3:ncol(data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf$Timepoint, subset = data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf$Timepoint %in% c("Late stress", "12h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_LS_R12_flag_leaf_FL_combined_2) = "p.value"
rownames(wilcox_res_Anjali_LS_R12_flag_leaf_FL_combined_2) = rownames(wilcox_sig_log2fc_C_LS_flag_leaf_FL_combined)
wilcox_res_Anjali_LS_R12_flag_leaf_FL_combined_2 = as.data.frame(wilcox_res_Anjali_LS_R12_flag_leaf_FL_combined_2)
wilcox_res_Anjali_LS_R12_flag_leaf_FL_combined_2$Significance = with(wilcox_res_Anjali_LS_R12_flag_leaf_FL_combined_2,
                                                                     ifelse(p.value <= 0.001, "***", ifelse(
                                                                       p.value <= 0.01, "**", ifelse(
                                                                         p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_LS_R12_flag_leaf_FL_combined_2$Significance == "*") #4 metabolites
sum(wilcox_res_Anjali_LS_R12_flag_leaf_FL_combined_2$Significance == "**") #1 metabolites
sum(wilcox_res_Anjali_LS_R12_flag_leaf_FL_combined_2$Significance == "***") #0 metabolites
sum(wilcox_res_Anjali_LS_R12_flag_leaf_FL_combined_2$Significance != "ns") #5 metabolites
sum(wilcox_res_Anjali_LS_R12_flag_leaf_FL_combined_2$Significance == "ns") #50 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_LS_R12_flag_leaf_FL_combined_2 = cbind(wilcox_res_N22_LS_R12_flag_leaf_FL_combined_2,
                                                  wilcox_res_Dular_LS_R12_flag_leaf_FL_combined_2$Significance, 
                                                  wilcox_res_Anjali_LS_R12_flag_leaf_FL_combined_2$Significance)
wilcox_sig_LS_R12_flag_leaf_FL_combined_2 = wilcox_sig_LS_R12_flag_leaf_FL_combined_2[ , -1]
colnames(wilcox_sig_LS_R12_flag_leaf_FL_combined_2)[1] = "wilcox_res_N22_LS_R12_flag_leaf_FL_combined$Significance"
wilcox_sig_LS_R12_flag_leaf_FL_combined_2 = data.frame(lapply(wilcox_sig_LS_R12_flag_leaf_FL_combined_2, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_LS_R12_flag_leaf_FL_combined_2) = rownames(wilcox_res_N22_LS_R12_flag_leaf_FL_combined_2)
wilcox_sig_LS_R12_flag_leaf_FL_combined_2$Nonsig = rowSums(wilcox_sig_LS_R12_flag_leaf_FL_combined_2 == " ")  #count number of nonsignificant
wilcox_sig_final_LS_R12_flag_leaf_FL_combined_2 = wilcox_sig_LS_R12_flag_leaf_FL_combined_2[-which(wilcox_sig_LS_R12_flag_leaf_FL_combined_2$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars

wilcox_sig_log2fc_LS_R12_flag_leaf_FL_combined_2 = log2fc_combined_latestress_recovery_FL_flag_leaf[which(rownames(log2fc_combined_latestress_recovery_FL_flag_leaf) %in% rownames(wilcox_sig_final_LS_R12_flag_leaf_FL_combined_2)), c(1, 4, 7)]  #log2FC values of metabolites significant in at least one of the cultivars


#venn diagram - common among cultivars

#N22
wilcox_sig_stress_responsive_N22_LS_R12_flag_leaf_FL_combined = wilcox_res_N22_LS_R12_flag_leaf_FL_combined_2[-which(wilcox_res_N22_LS_R12_flag_leaf_FL_combined_2$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_LS_R12_flag_leaf_FL_combined_2[match(rownames(wilcox_sig_stress_responsive_N22_LS_R12_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_LS_R12_flag_leaf_FL_combined_2)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_LS_R12_flag_leaf_FL_combined_2)[1]
wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined$`FNR12 - FNLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined$Up.Down == "Up")    #1 metabolite
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined$Up.Down == "Down")  #3 metabolites

#Dular
wilcox_sig_stress_responsive_Dular_LS_R12_flag_leaf_FL_combined = wilcox_res_Dular_LS_R12_flag_leaf_FL_combined_2[-which(wilcox_res_Dular_LS_R12_flag_leaf_FL_combined_2$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_LS_R12_flag_leaf_FL_combined_2[match(rownames(wilcox_sig_stress_responsive_Dular_LS_R12_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_LS_R12_flag_leaf_FL_combined_2)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_LS_R12_flag_leaf_FL_combined_2)[2]
wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined$`FDR12 - FDLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined$Up.Down == "Up")    #2 metabolites
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined$Up.Down == "Down")  #1 metabolite

#Anjali
wilcox_sig_stress_responsive_Anjali_LS_R12_flag_leaf_FL_combined = wilcox_res_Anjali_LS_R12_flag_leaf_FL_combined_2[-which(wilcox_res_Anjali_LS_R12_flag_leaf_FL_combined_2$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_LS_R12_flag_leaf_FL_combined_2[match(rownames(wilcox_sig_stress_responsive_Anjali_LS_R12_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_LS_R12_flag_leaf_FL_combined_2)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_LS_R12_flag_leaf_FL_combined_2)[3]
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined$`FAR12 - FALS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined$Up.Down == "Up")    #1 metabolite
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined$Up.Down == "Down")  #4 metabolites


#venn - increased - 12h rewatering/late stress
vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined) = c("N22", "Dular", "Anjali")
venn_up_stress_responsive_wilcox_LS_R12_flag_leaf_FL_combined = venn(vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_stress_responsive_R12_latestress_flag_leaf_FL_combined = venn.diagram(vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined, category.names = c(paste(names(vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined[1]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined[2]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined[3]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "A", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label A
intersect_stress_responsive_wilcox_LS_R12_venn_up_flag_leaf_FL_combined = attr(venn_up_stress_responsive_wilcox_LS_R12_flag_leaf_FL_combined, "intersections")  #list of unique and common metabolites
venn_up_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_FL_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R12_venn_up_flag_leaf_FL_combined, rbind))  #convert from list to matrix
colnames(venn_up_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_FL_combined) = venn_up_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_FL_combined[1, ]
venn_up_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_FL_combined = venn_up_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_FL_combined[-1, ]


#venn - decreased - 12h rewatering/late stress
vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined) = c("N22", "Dular", "Anjali")
venn_down_stress_responsive_wilcox_LS_R12_flag_leaf_FL_combined = venn(vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_stress_responsive_R12_latestress_flag_leaf_FL_combined = venn.diagram(vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined, category.names = c(paste(names(vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined[1]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined[2]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined[3]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_FL_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "B", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label B
intersect_stress_responsive_wilcox_LS_R12_venn_down_flag_leaf_FL_combined = attr(venn_down_stress_responsive_wilcox_LS_R12_flag_leaf_FL_combined, "intersections")  #list of unique and common metabolites
venn_down_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_FL_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R12_venn_down_flag_leaf_FL_combined, rbind))  #convert from list to matrix
colnames(venn_down_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_FL_combined) = venn_down_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_FL_combined[1, ]
venn_down_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_FL_combined = venn_down_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_FL_combined[-1, ]


write.table(venn_up_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_FL_combined, file = "wilcox_stress_responsive_sig_up_metabolites_flag_leaf_FL_R12_latestress_combined.txt", sep = "\t", quote = F)
write.table(venn_down_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_FL_combined, file = "wilcox_stress_responsive_sig_down_metabolites_flag_leaf_FL_R12_latestress_combined.txt", sep = "\t", quote = F)


####################

#late stress vs. 36h rewatering

#N22

#data distribution 
shapiro_res_combined_N22_latestress_recovery36_FL_flag_leaf_2 = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_FL_flag_leaf[which(data_combined_N22_log2median_stress_responsive_FL_flag_leaf$Timepoint %in% c("Late stress", "36h Rewatering")), 3:ncol(data_combined_N22_log2median_stress_responsive_FL_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_latestress_recovery36_FL_flag_leaf_2) = "p.value"
shapiro_res_combined_N22_latestress_recovery36_FL_flag_leaf_2 = as.data.frame(shapiro_res_combined_N22_latestress_recovery36_FL_flag_leaf_2)
shapiro_res_combined_N22_latestress_recovery36_FL_flag_leaf_2$Distribution = ifelse(shapiro_res_combined_N22_latestress_recovery36_FL_flag_leaf_2$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_latestress_recovery36_FL_flag_leaf_2$Distribution == "Normal")     #31 metabolites
sum(shapiro_res_combined_N22_latestress_recovery36_FL_flag_leaf_2$Distribution == "Non-normal") #24 metabolites

#wilcox test
wilcox_res_N22_LS_R36_flag_leaf_FL_combined_2 = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_FL_flag_leaf[ , 3:ncol(data_combined_N22_log2median_stress_responsive_FL_flag_leaf)], function(x) wilcox.test(x ~ data_combined_N22_log2median_stress_responsive_FL_flag_leaf$Timepoint, subset = data_combined_N22_log2median_stress_responsive_FL_flag_leaf$Timepoint %in% c("Late stress", "36h Rewatering"))$p.value)))
colnames(wilcox_res_N22_LS_R36_flag_leaf_FL_combined_2) = "p.value"
rownames(wilcox_res_N22_LS_R36_flag_leaf_FL_combined_2) = rownames(wilcox_sig_log2fc_C_LS_flag_leaf_FL_combined)
wilcox_res_N22_LS_R36_flag_leaf_FL_combined_2 = as.data.frame(wilcox_res_N22_LS_R36_flag_leaf_FL_combined_2)
wilcox_res_N22_LS_R36_flag_leaf_FL_combined_2$Significance = with(wilcox_res_N22_LS_R36_flag_leaf_FL_combined_2,
                                                                  ifelse(p.value <= 0.001, "***", ifelse(
                                                                    p.value <= 0.01, "**", ifelse(
                                                                      p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_LS_R36_flag_leaf_FL_combined_2$Significance == "*") #6 metabolites
sum(wilcox_res_N22_LS_R36_flag_leaf_FL_combined_2$Significance == "**") #5 metabolites
sum(wilcox_res_N22_LS_R36_flag_leaf_FL_combined_2$Significance == "***") #4 metabolites
sum(wilcox_res_N22_LS_R36_flag_leaf_FL_combined_2$Significance != "ns") #15 metabolites
sum(wilcox_res_N22_LS_R36_flag_leaf_FL_combined_2$Significance == "ns") #40 metabolites 


#Dular

#data distribution 
shapiro_res_combined_Dular_latestress_recovery36_FL_flag_leaf_2 = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_FL_flag_leaf[which(data_combined_Dular_log2median_stress_responsive_FL_flag_leaf$Timepoint %in% c("Late stress", "36h Rewatering")), 3:ncol(data_combined_Dular_log2median_stress_responsive_FL_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_latestress_recovery36_FL_flag_leaf_2) = "p.value"
shapiro_res_combined_Dular_latestress_recovery36_FL_flag_leaf_2 = as.data.frame(shapiro_res_combined_Dular_latestress_recovery36_FL_flag_leaf_2)
shapiro_res_combined_Dular_latestress_recovery36_FL_flag_leaf_2$Distribution = ifelse(shapiro_res_combined_Dular_latestress_recovery36_FL_flag_leaf_2$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_latestress_recovery36_FL_flag_leaf_2$Distribution == "Normal")     #35 metabolites
sum(shapiro_res_combined_Dular_latestress_recovery36_FL_flag_leaf_2$Distribution == "Non-normal") #20 metabolites

#wilcox test
wilcox_res_Dular_LS_R36_flag_leaf_FL_combined_2 = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_FL_flag_leaf[ , 3:ncol(data_combined_Dular_log2median_stress_responsive_FL_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_stress_responsive_FL_flag_leaf$Timepoint, subset = data_combined_Dular_log2median_stress_responsive_FL_flag_leaf$Timepoint %in% c("Late stress", "36h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_LS_R36_flag_leaf_FL_combined_2) = "p.value"
rownames(wilcox_res_Dular_LS_R36_flag_leaf_FL_combined_2) = rownames(wilcox_sig_log2fc_C_LS_flag_leaf_FL_combined)
wilcox_res_Dular_LS_R36_flag_leaf_FL_combined_2 = as.data.frame(wilcox_res_Dular_LS_R36_flag_leaf_FL_combined_2)
wilcox_res_Dular_LS_R36_flag_leaf_FL_combined_2$Significance = with(wilcox_res_Dular_LS_R36_flag_leaf_FL_combined_2,
                                                                    ifelse(p.value <= 0.001, "***", ifelse(
                                                                      p.value <= 0.01, "**", ifelse(
                                                                        p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_LS_R36_flag_leaf_FL_combined_2$Significance == "*") #5 metabolites
sum(wilcox_res_Dular_LS_R36_flag_leaf_FL_combined_2$Significance == "**") #6 metabolites
sum(wilcox_res_Dular_LS_R36_flag_leaf_FL_combined_2$Significance == "***") #7 metabolites
sum(wilcox_res_Dular_LS_R36_flag_leaf_FL_combined_2$Significance != "ns") #18 metabolites
sum(wilcox_res_Dular_LS_R36_flag_leaf_FL_combined_2$Significance == "ns") #37 metabolites


#Anjali

#data distribution 
shapiro_res_combined_Anjali_latestress_recovery36_FL_flag_leaf_2 = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf[which(data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf$Timepoint %in% c("Late stress", "36h Rewatering")), 3:ncol(data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_latestress_recovery36_FL_flag_leaf_2) = "p.value"
shapiro_res_combined_Anjali_latestress_recovery36_FL_flag_leaf_2 = as.data.frame(shapiro_res_combined_Anjali_latestress_recovery36_FL_flag_leaf_2)
shapiro_res_combined_Anjali_latestress_recovery36_FL_flag_leaf_2$Distribution = ifelse(shapiro_res_combined_Anjali_latestress_recovery36_FL_flag_leaf_2$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_latestress_recovery36_FL_flag_leaf_2$Distribution == "Normal")     #42 metabolites
sum(shapiro_res_combined_Anjali_latestress_recovery36_FL_flag_leaf_2$Distribution == "Non-normal") #13 metabolites

#wilcox test
wilcox_res_Anjali_LS_R36_flag_leaf_FL_combined_2 = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf[ , 3:ncol(data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf$Timepoint, subset = data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf$Timepoint %in% c("Late stress", "36h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_LS_R36_flag_leaf_FL_combined_2) = "p.value"
rownames(wilcox_res_Anjali_LS_R36_flag_leaf_FL_combined_2) = rownames(wilcox_sig_log2fc_C_LS_flag_leaf_FL_combined)
wilcox_res_Anjali_LS_R36_flag_leaf_FL_combined_2 = as.data.frame(wilcox_res_Anjali_LS_R36_flag_leaf_FL_combined_2)
wilcox_res_Anjali_LS_R36_flag_leaf_FL_combined_2$Significance = with(wilcox_res_Anjali_LS_R36_flag_leaf_FL_combined_2,
                                                                     ifelse(p.value <= 0.001, "***", ifelse(
                                                                       p.value <= 0.01, "**", ifelse(
                                                                         p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_LS_R36_flag_leaf_FL_combined_2$Significance == "*") #9 metabolites
sum(wilcox_res_Anjali_LS_R36_flag_leaf_FL_combined_2$Significance == "**") #10 metabolites
sum(wilcox_res_Anjali_LS_R36_flag_leaf_FL_combined_2$Significance == "***") #6 metabolites
sum(wilcox_res_Anjali_LS_R36_flag_leaf_FL_combined_2$Significance != "ns") #25 metabolites
sum(wilcox_res_Anjali_LS_R36_flag_leaf_FL_combined_2$Significance == "ns") #30 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_LS_R36_flag_leaf_FL_combined_2 = cbind(wilcox_res_N22_LS_R36_flag_leaf_FL_combined_2, 
                                                  wilcox_res_Dular_LS_R36_flag_leaf_FL_combined_2$Significance, 
                                                  wilcox_res_Anjali_LS_R36_flag_leaf_FL_combined_2$Significance)
wilcox_sig_LS_R36_flag_leaf_FL_combined_2 = wilcox_sig_LS_R36_flag_leaf_FL_combined_2[ , -1]
colnames(wilcox_sig_LS_R36_flag_leaf_FL_combined_2)[1] = "wilcox_res_N22_LS_R36_flag_leaf_FL_combined$Significance"
wilcox_sig_LS_R36_flag_leaf_FL_combined_2 = data.frame(lapply(wilcox_sig_LS_R36_flag_leaf_FL_combined_2, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_LS_R36_flag_leaf_FL_combined_2) = rownames(wilcox_res_N22_LS_R36_flag_leaf_FL_combined_2)
wilcox_sig_LS_R36_flag_leaf_FL_combined_2$Nonsig = rowSums(wilcox_sig_LS_R36_flag_leaf_FL_combined_2 == " ")  #count number of nonsignificant
wilcox_sig_final_LS_R36_flag_leaf_FL_combined_2 = wilcox_sig_LS_R36_flag_leaf_FL_combined_2[-which(wilcox_sig_LS_R36_flag_leaf_FL_combined_2$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars

wilcox_sig_log2fc_LS_R36_flag_leaf_FL_combined_2 = log2fc_combined_latestress_recovery_FL_flag_leaf[which(rownames(log2fc_combined_latestress_recovery_FL_flag_leaf) %in% rownames(wilcox_sig_final_LS_R36_flag_leaf_FL_combined_2)), c(2, 5, 8)]  #log2FC values of metabolites significant in at least one of the cultivars


#venn diagram - common among cultivars

#N22
wilcox_sig_stress_responsive_N22_LS_R36_flag_leaf_FL_combined = wilcox_res_N22_LS_R36_flag_leaf_FL_combined_2[-which(wilcox_res_N22_LS_R36_flag_leaf_FL_combined_2$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_LS_R36_flag_leaf_FL_combined_2[match(rownames(wilcox_sig_stress_responsive_N22_LS_R36_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_LS_R36_flag_leaf_FL_combined_2)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_LS_R36_flag_leaf_FL_combined_2)[1]
wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined$`FNR36 - FNLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined$Up.Down == "Up")    #2 metabolites
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined$Up.Down == "Down")  #13 metabolites

#Dular
wilcox_sig_stress_responsive_Dular_LS_R36_flag_leaf_FL_combined = wilcox_res_Dular_LS_R36_flag_leaf_FL_combined_2[-which(wilcox_res_Dular_LS_R36_flag_leaf_FL_combined_2$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_LS_R36_flag_leaf_FL_combined_2[match(rownames(wilcox_sig_stress_responsive_Dular_LS_R36_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_LS_R36_flag_leaf_FL_combined_2)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_LS_R36_flag_leaf_FL_combined_2)[2]
wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined$`FDR36 - FDLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined$Up.Down == "Up")    #6 metabolites
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined$Up.Down == "Down")  #12 metabolites

#Anjali
wilcox_sig_stress_responsive_Anjali_LS_R36_flag_leaf_FL_combined = wilcox_res_Anjali_LS_R36_flag_leaf_FL_combined_2[-which(wilcox_res_Anjali_LS_R36_flag_leaf_FL_combined_2$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_LS_R36_flag_leaf_FL_combined_2[match(rownames(wilcox_sig_stress_responsive_Anjali_LS_R36_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_LS_R36_flag_leaf_FL_combined_2)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_LS_R36_flag_leaf_FL_combined_2)[3]
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined$`FAR36 - FALS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined$Up.Down == "Up")    #8 metabolites
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined$Up.Down == "Down")  #17 metabolites


#venn - increased - 36h rewatering/late stress
vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined) = c("N22", "Dular", "Anjali")
venn_up_stress_responsive_wilcox_LS_R36_flag_leaf_FL_combined = venn(vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_stress_responsive_R36_latestress_flag_leaf_FL_combined = venn.diagram(vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined, category.names = c(paste(names(vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined[1]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined[2]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined[3]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "C", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label C
intersect_stress_responsive_wilcox_LS_R36_venn_up_flag_leaf_FL_combined = attr(venn_up_stress_responsive_wilcox_LS_R36_flag_leaf_FL_combined, "intersections")  #list of unique and common metabolites
venn_up_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_FL_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R36_venn_up_flag_leaf_FL_combined, rbind))  #convert from list to matrix
colnames(venn_up_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_FL_combined) = venn_up_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_FL_combined[1, ]
venn_up_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_FL_combined = venn_up_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_FL_combined[-1, ]


#venn - decreased - 36h rewatering/late stress
vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined) = c("N22", "Dular", "Anjali")
venn_down_stress_responsive_wilcox_LS_R36_flag_leaf_FL_combined = venn(vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_stress_responsive_R36_latestress_flag_leaf_FL_combined = venn.diagram(vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined, category.names = c(paste(names(vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined[1]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined[2]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined[3]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_FL_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "D", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label D
intersect_stress_responsive_wilcox_LS_R36_venn_down_flag_leaf_FL_combined = attr(venn_down_stress_responsive_wilcox_LS_R36_flag_leaf_FL_combined, "intersections")  #list of unique and common metabolites
venn_down_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_FL_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R36_venn_down_flag_leaf_FL_combined, rbind))  #convert from list to matrix
colnames(venn_down_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_FL_combined) = venn_down_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_FL_combined[1, ]
venn_down_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_FL_combined = venn_down_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_FL_combined[-1, ]


write.table(venn_up_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_FL_combined, file = "wilcox_stress_responsive_sig_up_metabolites_flag_leaf_FL_R36_latestress_combined.txt", sep = "\t", quote = F)
write.table(venn_down_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_FL_combined, file = "wilcox_stress_responsive_sig_down_metabolites_flag_leaf_FL_R36_latestress_combined.txt", sep = "\t", quote = F)


####################

#late stress vs. 60h rewatering

#N22

#data distribution 
shapiro_res_combined_N22_latestress_recovery60_FL_flag_leaf_2 = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_FL_flag_leaf[which(data_combined_N22_log2median_stress_responsive_FL_flag_leaf$Timepoint %in% c("Late stress", "60h Rewatering")), 3:ncol(data_combined_N22_log2median_stress_responsive_FL_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_latestress_recovery60_FL_flag_leaf_2) = "p.value"
shapiro_res_combined_N22_latestress_recovery60_FL_flag_leaf_2 = as.data.frame(shapiro_res_combined_N22_latestress_recovery60_FL_flag_leaf_2)
shapiro_res_combined_N22_latestress_recovery60_FL_flag_leaf_2$Distribution = ifelse(shapiro_res_combined_N22_latestress_recovery60_FL_flag_leaf_2$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_latestress_recovery60_FL_flag_leaf_2$Distribution == "Normal")     #35 metabolites
sum(shapiro_res_combined_N22_latestress_recovery60_FL_flag_leaf_2$Distribution == "Non-normal") #20 metabolites

#wilcox test
wilcox_res_N22_LS_R60_flag_leaf_FL_combined_2 = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_FL_flag_leaf[ , 3:ncol(data_combined_N22_log2median_stress_responsive_FL_flag_leaf)], function(x) wilcox.test(x ~ data_combined_N22_log2median_stress_responsive_FL_flag_leaf$Timepoint, subset = data_combined_N22_log2median_stress_responsive_FL_flag_leaf$Timepoint %in% c("Late stress", "60h Rewatering"))$p.value)))
colnames(wilcox_res_N22_LS_R60_flag_leaf_FL_combined_2) = "p.value"
rownames(wilcox_res_N22_LS_R60_flag_leaf_FL_combined_2) = rownames(wilcox_sig_log2fc_C_LS_flag_leaf_FL_combined)
wilcox_res_N22_LS_R60_flag_leaf_FL_combined_2 = as.data.frame(wilcox_res_N22_LS_R60_flag_leaf_FL_combined_2)
wilcox_res_N22_LS_R60_flag_leaf_FL_combined_2$Significance = with(wilcox_res_N22_LS_R60_flag_leaf_FL_combined_2,
                                                                  ifelse(p.value <= 0.001, "***", ifelse(
                                                                    p.value <= 0.01, "**", ifelse(
                                                                      p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_LS_R60_flag_leaf_FL_combined_2$Significance == "*") #3 metabolites
sum(wilcox_res_N22_LS_R60_flag_leaf_FL_combined_2$Significance == "**") #6 metabolites
sum(wilcox_res_N22_LS_R60_flag_leaf_FL_combined_2$Significance == "***") #5 metabolites
sum(wilcox_res_N22_LS_R60_flag_leaf_FL_combined_2$Significance != "ns") #14 metabolites
sum(wilcox_res_N22_LS_R60_flag_leaf_FL_combined_2$Significance == "ns") #41 metabolites 


#Dular

#data distribution 
shapiro_res_combined_Dular_latestress_recovery60_FL_flag_leaf_2 = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_FL_flag_leaf[which(data_combined_Dular_log2median_stress_responsive_FL_flag_leaf$Timepoint %in% c("Late stress", "60h Rewatering")), 3:ncol(data_combined_Dular_log2median_stress_responsive_FL_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_latestress_recovery60_FL_flag_leaf_2) = "p.value"
shapiro_res_combined_Dular_latestress_recovery60_FL_flag_leaf_2 = as.data.frame(shapiro_res_combined_Dular_latestress_recovery60_FL_flag_leaf_2)
shapiro_res_combined_Dular_latestress_recovery60_FL_flag_leaf_2$Distribution = ifelse(shapiro_res_combined_Dular_latestress_recovery60_FL_flag_leaf_2$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_latestress_recovery60_FL_flag_leaf_2$Distribution == "Normal")     #37 metabolites
sum(shapiro_res_combined_Dular_latestress_recovery60_FL_flag_leaf_2$Distribution == "Non-normal") #18 metabolites

#wilcox test
wilcox_res_Dular_LS_R60_flag_leaf_FL_combined_2 = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_FL_flag_leaf[ , 3:ncol(data_combined_Dular_log2median_stress_responsive_FL_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_stress_responsive_FL_flag_leaf$Timepoint, subset = data_combined_Dular_log2median_stress_responsive_FL_flag_leaf$Timepoint %in% c("Late stress", "60h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_LS_R60_flag_leaf_FL_combined_2) = "p.value"
rownames(wilcox_res_Dular_LS_R60_flag_leaf_FL_combined_2) = rownames(wilcox_sig_log2fc_C_LS_flag_leaf_FL_combined)
wilcox_res_Dular_LS_R60_flag_leaf_FL_combined_2 = as.data.frame(wilcox_res_Dular_LS_R60_flag_leaf_FL_combined_2)
wilcox_res_Dular_LS_R60_flag_leaf_FL_combined_2$Significance = with(wilcox_res_Dular_LS_R60_flag_leaf_FL_combined_2,
                                                                    ifelse(p.value <= 0.001, "***", ifelse(
                                                                      p.value <= 0.01, "**", ifelse(
                                                                        p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_LS_R60_flag_leaf_FL_combined_2$Significance == "*") #12 metabolites
sum(wilcox_res_Dular_LS_R60_flag_leaf_FL_combined_2$Significance == "**") #5 metabolites
sum(wilcox_res_Dular_LS_R60_flag_leaf_FL_combined_2$Significance == "***") #5 metabolites
sum(wilcox_res_Dular_LS_R60_flag_leaf_FL_combined_2$Significance != "ns") #22 metabolites
sum(wilcox_res_Dular_LS_R60_flag_leaf_FL_combined_2$Significance == "ns") #33 metabolites


#Anjali

#data distribution 
shapiro_res_combined_Anjali_latestress_recovery60_FL_flag_leaf_2 = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf[which(data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf$Timepoint %in% c("Late stress", "60h Rewatering")), 3:ncol(data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_latestress_recovery60_FL_flag_leaf_2) = "p.value"
shapiro_res_combined_Anjali_latestress_recovery60_FL_flag_leaf_2 = as.data.frame(shapiro_res_combined_Anjali_latestress_recovery60_FL_flag_leaf_2)
shapiro_res_combined_Anjali_latestress_recovery60_FL_flag_leaf_2$Distribution = ifelse(shapiro_res_combined_Anjali_latestress_recovery60_FL_flag_leaf_2$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_latestress_recovery60_FL_flag_leaf_2$Distribution == "Normal")     #44 metabolites
sum(shapiro_res_combined_Anjali_latestress_recovery60_FL_flag_leaf_2$Distribution == "Non-normal") #11 metabolites

#wilcox test
wilcox_res_Anjali_LS_R60_flag_leaf_FL_combined_2 = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf[ , 3:ncol(data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf$Timepoint, subset = data_combined_Anjali_log2median_stress_responsive_FL_flag_leaf$Timepoint %in% c("Late stress", "60h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_LS_R60_flag_leaf_FL_combined_2) = "p.value"
rownames(wilcox_res_Anjali_LS_R60_flag_leaf_FL_combined_2) = rownames(wilcox_sig_log2fc_C_LS_flag_leaf_FL_combined)
wilcox_res_Anjali_LS_R60_flag_leaf_FL_combined_2 = as.data.frame(wilcox_res_Anjali_LS_R60_flag_leaf_FL_combined_2)
wilcox_res_Anjali_LS_R60_flag_leaf_FL_combined_2$Significance = with(wilcox_res_Anjali_LS_R60_flag_leaf_FL_combined_2,
                                                                     ifelse(p.value <= 0.001, "***", ifelse(
                                                                       p.value <= 0.01, "**", ifelse(
                                                                         p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_LS_R60_flag_leaf_FL_combined_2$Significance == "*") #8 metabolites
sum(wilcox_res_Anjali_LS_R60_flag_leaf_FL_combined_2$Significance == "**") #6 metabolites
sum(wilcox_res_Anjali_LS_R60_flag_leaf_FL_combined_2$Significance == "***") #4 metabolites
sum(wilcox_res_Anjali_LS_R60_flag_leaf_FL_combined_2$Significance != "ns") #18 metabolites
sum(wilcox_res_Anjali_LS_R60_flag_leaf_FL_combined_2$Significance == "ns") #37 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_LS_R60_flag_leaf_FL_combined_2 = cbind(wilcox_res_N22_LS_R60_flag_leaf_FL_combined_2, 
                                                  wilcox_res_Dular_LS_R60_flag_leaf_FL_combined_2$Significance, 
                                                  wilcox_res_Anjali_LS_R60_flag_leaf_FL_combined_2$Significance)
wilcox_sig_LS_R60_flag_leaf_FL_combined_2 = wilcox_sig_LS_R60_flag_leaf_FL_combined_2[ , -1]
colnames(wilcox_sig_LS_R60_flag_leaf_FL_combined_2)[1] = "wilcox_res_N22_LS_R60_flag_leaf_FL_combined$Significance"
wilcox_sig_LS_R60_flag_leaf_FL_combined_2 = data.frame(lapply(wilcox_sig_LS_R60_flag_leaf_FL_combined_2, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_LS_R60_flag_leaf_FL_combined_2) = rownames(wilcox_res_N22_LS_R60_flag_leaf_FL_combined_2)
wilcox_sig_LS_R60_flag_leaf_FL_combined_2$Nonsig = rowSums(wilcox_sig_LS_R60_flag_leaf_FL_combined_2 == " ")  #count number of nonsignificant
wilcox_sig_final_LS_R60_flag_leaf_FL_combined_2 = wilcox_sig_LS_R60_flag_leaf_FL_combined_2[-which(wilcox_sig_LS_R60_flag_leaf_FL_combined_2$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars

wilcox_sig_log2fc_LS_R60_flag_leaf_FL_combined_2 = log2fc_combined_latestress_recovery_FL_flag_leaf[which(rownames(log2fc_combined_latestress_recovery_FL_flag_leaf) %in% rownames(wilcox_sig_final_LS_R60_flag_leaf_FL_combined_2)), c(3, 6, 9)]  #log2FC values of metabolites significant in at least one of the cultivars


#venn diagram - common among cultivars

#N22
wilcox_sig_stress_responsive_N22_LS_R60_flag_leaf_FL_combined = wilcox_res_N22_LS_R60_flag_leaf_FL_combined_2[-which(wilcox_res_N22_LS_R60_flag_leaf_FL_combined_2$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_LS_R60_flag_leaf_FL_combined_2[match(rownames(wilcox_sig_stress_responsive_N22_LS_R60_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_LS_R60_flag_leaf_FL_combined_2)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_LS_R60_flag_leaf_FL_combined_2)[1]
wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined$`FNR60 - FNLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined$Up.Down == "Up")    #1 metabolite
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined$Up.Down == "Down")  #13 metabolites

#Dular
wilcox_sig_stress_responsive_Dular_LS_R60_flag_leaf_FL_combined = wilcox_res_Dular_LS_R60_flag_leaf_FL_combined_2[-which(wilcox_res_Dular_LS_R60_flag_leaf_FL_combined_2$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_LS_R60_flag_leaf_FL_combined_2[match(rownames(wilcox_sig_stress_responsive_Dular_LS_R60_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_LS_R60_flag_leaf_FL_combined_2)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_LS_R60_flag_leaf_FL_combined_2)[2]
wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined$`FDR60 - FDLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined$Up.Down == "Up")    #5 metabolites
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined$Up.Down == "Down")  #17 metabolites

#Anjali
wilcox_sig_stress_responsive_Anjali_LS_R60_flag_leaf_FL_combined = wilcox_res_Anjali_LS_R60_flag_leaf_FL_combined_2[-which(wilcox_res_Anjali_LS_R60_flag_leaf_FL_combined_2$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_LS_R60_flag_leaf_FL_combined_2[match(rownames(wilcox_sig_stress_responsive_Anjali_LS_R60_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_LS_R60_flag_leaf_FL_combined_2)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_LS_R60_flag_leaf_FL_combined_2)[3]
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined$`FAR60 - FALS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined$Up.Down == "Up")    #6 metabolites
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined$Up.Down == "Down")  #12 metabolites


#venn - increased - 60h rewatering/late stress
vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined) = c("N22", "Dular", "Anjali")
venn_up_stress_responsive_wilcox_LS_R60_flag_leaf_FL_combined = venn(vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_stress_responsive_R60_latestress_flag_leaf_FL_combined = venn.diagram(vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined, category.names = c(paste(names(vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined[1]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined[2]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined[3]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "E", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label E
intersect_stress_responsive_wilcox_LS_R60_venn_up_flag_leaf_FL_combined = attr(venn_up_stress_responsive_wilcox_LS_R60_flag_leaf_FL_combined, "intersections")  #list of unique and common metabolites
venn_up_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_FL_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R60_venn_up_flag_leaf_FL_combined, rbind))  #convert from list to matrix
colnames(venn_up_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_FL_combined) = venn_up_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_FL_combined[1, ]
venn_up_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_FL_combined = venn_up_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_FL_combined[-1, ]


#venn - decreased - 60h rewatering/late stress
vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined) = c("N22", "Dular", "Anjali")
venn_down_stress_responsive_wilcox_LS_R60_flag_leaf_FL_combined = venn(vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_stress_responsive_R60_latestress_flag_leaf_FL_combined = venn.diagram(vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined, category.names = c(paste(names(vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined[1]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined[2]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined[3]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_FL_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "F", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label F
intersect_stress_responsive_wilcox_LS_R60_venn_down_flag_leaf_FL_combined = attr(venn_down_stress_responsive_wilcox_LS_R60_flag_leaf_FL_combined, "intersections")  #list of unique and common metabolites
venn_down_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_FL_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R60_venn_down_flag_leaf_FL_combined, rbind))  #convert from list to matrix
colnames(venn_down_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_FL_combined) = venn_down_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_FL_combined[1, ]
venn_down_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_FL_combined = venn_down_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_FL_combined[-1, ]


write.table(venn_up_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_FL_combined, file = "wilcox_stress_responsive_sig_up_metabolites_flag_leaf_FL_R60_latestress_combined.txt", sep = "\t", quote = F)
write.table(venn_down_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_FL_combined, file = "wilcox_stress_responsive_sig_down_metabolites_flag_leaf_FL_R60_latestress_combined.txt", sep = "\t", quote = F)


####################

#significance of wilcoxon test - all cultivars, all three RW timepoints
wilcox_sig_LS_R12_R36_R60_flag_leaf_FL_combined_2 = cbind(wilcox_res_N22_LS_R12_flag_leaf_FL_combined_2,
                                                          wilcox_res_N22_LS_R36_flag_leaf_FL_combined_2$Significance,
                                                          wilcox_res_N22_LS_R60_flag_leaf_FL_combined_2$Significance,
                                                          wilcox_res_Dular_LS_R12_flag_leaf_FL_combined_2$Significance,
                                                          wilcox_res_Dular_LS_R36_flag_leaf_FL_combined_2$Significance,
                                                          wilcox_res_Dular_LS_R60_flag_leaf_FL_combined_2$Significance,
                                                          wilcox_res_Anjali_LS_R12_flag_leaf_FL_combined_2$Significance,
                                                          wilcox_res_Anjali_LS_R36_flag_leaf_FL_combined_2$Significance,
                                                          wilcox_res_Anjali_LS_R60_flag_leaf_FL_combined_2$Significance)
wilcox_sig_LS_R12_R36_R60_flag_leaf_FL_combined_2 = wilcox_sig_LS_R12_R36_R60_flag_leaf_FL_combined_2[ , -1]
colnames(wilcox_sig_LS_R12_R36_R60_flag_leaf_FL_combined_2)[1] = "wilcox_res_N22_LS_R12_flag_leaf_FL_combined$Significance"
wilcox_sig_LS_R12_R36_R60_flag_leaf_FL_combined_2 = data.frame(lapply(wilcox_sig_LS_R12_R36_R60_flag_leaf_FL_combined_2, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_LS_R12_R36_R60_flag_leaf_FL_combined_2) = rownames(wilcox_res_N22_LS_R12_flag_leaf_FL_combined_2)
wilcox_sig_LS_R12_R36_R60_flag_leaf_FL_combined_2$Nonsig = rowSums(wilcox_sig_LS_R12_R36_R60_flag_leaf_FL_combined_2 == " ")  #count number of nonsignificant
wilcox_sig_final_LS_R12_R36_R60_flag_leaf_FL_combined_3 = wilcox_sig_LS_R12_R36_R60_flag_leaf_FL_combined_2[-which(wilcox_sig_LS_R12_R36_R60_flag_leaf_FL_combined_2$Nonsig == "9"), ] #metabolites significant in at least one of the comparisons - 41 metabolites

#log2fc values of metabolites with significant difference in any one of the comparisons - 12, 36, 60h recovery vs. late stress in any cultivar
wilcox_sig_log2fc_LS_R12_R36_R60_flag_leaf_FL_combined_2 = log2fc_combined_latestress_recovery_FL_flag_leaf[which(rownames(log2fc_combined_latestress_recovery_FL_flag_leaf) %in% rownames(wilcox_sig_final_LS_R12_R36_R60_flag_leaf_FL_combined_3)), ] 


#visualization
png("wilcox_sig_stress_responsive_log2fc_LS_RW_flag_leaf_FL_combined.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_log2fc_LS_R12_R36_R60_flag_leaf_FL_combined_2, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(7.5, 23), 
          lmat = rbind(c(0, 4, 5), c(0, 1, 0), c(3, 2, 0)),
          labCol = rep(c("12 h RW/Severe stress", "36 h RW/Severe stress", "60 h RW/Severe stress"), 3),
          cexCol = 1, 
          cexRow = 1.2, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 0.2, 11), 
          lwid = c(0.8, 5, 1.5),
          srtCol = 45, 
          cellnote = wilcox_sig_final_LS_R12_R36_R60_flag_leaf_FL_combined_3, 
          notecol = "black", 
          notecex = 1.8, 
          sepcolor = "black", 
          colsep = c(seq(0, 9)), 
          rowsep = c(seq(0, 41)), 
          sepwidth = c(0.01, 0.01), 
          ColSideColors = c(rep("#66c2a5", 3), rep("#fc8d62", 3), rep("#8da0cb", 3)),
          hclustfun = function(x) hclust(x, method = "average"))
legend(x = -0.009, y = 1.1, xpd = TRUE, legend = "N22", bty = "n", cex = 1.2)
legend(x = 0.105, y = 1.1, xpd = TRUE, legend = "Dular", bty = "n", cex = 1.2)
legend(x = 0.22, y = 1.1, xpd = TRUE, legend = "Anjali", bty = "n", cex = 1.2)
legend(x = 0.865, y = 1.04, xpd = TRUE, legend = expression(paste("Log" ["2"], "-fold change")), bty = "n", cex = 0.8)
legend(x = -0.15, y = 1.13, xpd = TRUE, legend = c("G"), bty = "n", cex = 1.05, text.font = 2)
dev.off()


####################

#all venn diagrams

#multiplot
png("venn_up_down_flag_leaf_FL_stress_responsive_RW_latestress_combined.png", width = 16*300, height = 24*300, res = 300)
grid.arrange(gTree(children = venn_up_stress_responsive_R12_latestress_flag_leaf_FL_combined), 
             gTree(children = venn_down_stress_responsive_R12_latestress_flag_leaf_FL_combined), 
             gTree(children = venn_up_stress_responsive_R36_latestress_flag_leaf_FL_combined),
             gTree(children = venn_down_stress_responsive_R36_latestress_flag_leaf_FL_combined),
             gTree(children = venn_up_stress_responsive_R60_latestress_flag_leaf_FL_combined),
             gTree(children = venn_down_stress_responsive_R60_latestress_flag_leaf_FL_combined),
             ncol = 2, nrow = 3)
dev.off()


####################


#recovery vs. corresponding control - NOTE: consider all metabolites, not just the stress-responsive

#log2-fold change - control vs. rewatering
control_recovery_order_FL = c("FNC", "FNR12", "FNR36", "FNR60", "FDC", "FDR12", "FDR36", "FDR60", "FAC", "FAR12", "FAR36", "FAR60") #line command already ran previously
data_combined_control_recovery_log2_mean_FL_flag_leaf = data_combined_log2_mean_FL_flag_leaf_2[control_recovery_order_FL, ]  #control and recovery data in preferred order
data_combined_control_recovery_log2_mean_FL_flag_leaf = t(data_combined_control_recovery_log2_mean_FL_flag_leaf)

log2fc_combined_N22_control_recovery_FL_flag_leaf = data_combined_control_recovery_log2_mean_FL_flag_leaf[ , 2:4] - data_combined_control_recovery_log2_mean_FL_flag_leaf[ , 1]
log2fc_combined_Dular_control_recovery_FL_flag_leaf = data_combined_control_recovery_log2_mean_FL_flag_leaf[ , 6:8] - data_combined_control_recovery_log2_mean_FL_flag_leaf[ , 5]
log2fc_combined_Anjali_control_recovery_FL_flag_leaf = data_combined_control_recovery_log2_mean_FL_flag_leaf[ , 10:12] - data_combined_control_recovery_log2_mean_FL_flag_leaf[ , 9]

#combine in one data set
log2fc_combined_control_recovery_FL_flag_leaf = cbind(log2fc_combined_N22_control_recovery_FL_flag_leaf, 
                                                      log2fc_combined_Dular_control_recovery_FL_flag_leaf, 
                                                      log2fc_combined_Anjali_control_recovery_FL_flag_leaf)   #combine log2-fc of all cultivars in one dataset
#rename columns
colnames(log2fc_combined_control_recovery_FL_flag_leaf)[1:3] = c(paste(colnames(log2fc_combined_control_recovery_FL_flag_leaf)[1:3], "- FNC"))
colnames(log2fc_combined_control_recovery_FL_flag_leaf)[4:6] = c(paste(colnames(log2fc_combined_control_recovery_FL_flag_leaf)[4:6], "- FDC"))
colnames(log2fc_combined_control_recovery_FL_flag_leaf)[7:9] = c(paste(colnames(log2fc_combined_control_recovery_FL_flag_leaf)[7:9], "- FAC"))

#heatmap including all metabolites
heatmap.2(log2fc_combined_control_recovery_FL_flag_leaf, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(7, 20),  main = "Flag leaf_FL - Combined - Log2 FC, Recovery/Control", 
          cexCol = 1, 
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          srtCol = 90)


####################

#control vs. 12h rewatering


#N22

#data distribution
shapiro_res_combined_N22_control_recovery12_FL_flag_leaf = t(data.frame(lapply(data_combined_N22_log2median_FL_flag_leaf[which(data_combined_N22_log2median_FL_flag_leaf$Timepoint %in% c("Control", "12h Rewatering")), 3:ncol(data_combined_N22_log2median_FL_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_control_recovery12_FL_flag_leaf) = "p.value"
shapiro_res_combined_N22_control_recovery12_FL_flag_leaf = as.data.frame(shapiro_res_combined_N22_control_recovery12_FL_flag_leaf)
shapiro_res_combined_N22_control_recovery12_FL_flag_leaf$Distribution = ifelse(shapiro_res_combined_N22_control_recovery12_FL_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_control_recovery12_FL_flag_leaf$Distribution == "Normal")     #63 metabolites
sum(shapiro_res_combined_N22_control_recovery12_FL_flag_leaf$Distribution == "Non-normal") #18 metabolites

#wilcox test
wilcox_res_N22_C_R12_flag_leaf_FL_combined = t(data.frame(lapply(data_combined_N22_log2median_FL_flag_leaf[ , 3:ncol(data_combined_N22_log2median_FL_flag_leaf)], function(x) wilcox.test(x ~ data_combined_N22_log2median_FL_flag_leaf$Timepoint, subset = data_combined_N22_log2median_FL_flag_leaf$Timepoint %in% c("Control", "12h Rewatering"))$p.value)))
colnames(wilcox_res_N22_C_R12_flag_leaf_FL_combined) = "p.value"
rownames(wilcox_res_N22_C_R12_flag_leaf_FL_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_N22_C_R12_flag_leaf_FL_combined = as.data.frame(wilcox_res_N22_C_R12_flag_leaf_FL_combined)
wilcox_res_N22_C_R12_flag_leaf_FL_combined$Significance = with(wilcox_res_N22_C_R12_flag_leaf_FL_combined,
                                                               ifelse(p.value <= 0.001, "***", ifelse(
                                                                 p.value <= 0.01, "**", ifelse(
                                                                   p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_C_R12_flag_leaf_FL_combined$Significance == "*") #8 metabolites
sum(wilcox_res_N22_C_R12_flag_leaf_FL_combined$Significance == "**") #9 metabolites
sum(wilcox_res_N22_C_R12_flag_leaf_FL_combined$Significance == "***") #13 metabolites
sum(wilcox_res_N22_C_R12_flag_leaf_FL_combined$Significance != "ns") #30 metabolites
sum(wilcox_res_N22_C_R12_flag_leaf_FL_combined$Significance == "ns") #51 metabolites


#Dular

#data distribution
shapiro_res_combined_Dular_control_recovery12_FL_flag_leaf = t(data.frame(lapply(data_combined_Dular_log2median_FL_flag_leaf[which(data_combined_Dular_log2median_FL_flag_leaf$Timepoint %in% c("Control", "12h Rewatering")), 3:ncol(data_combined_Dular_log2median_FL_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_control_recovery12_FL_flag_leaf) = "p.value"
shapiro_res_combined_Dular_control_recovery12_FL_flag_leaf = as.data.frame(shapiro_res_combined_Dular_control_recovery12_FL_flag_leaf)
shapiro_res_combined_Dular_control_recovery12_FL_flag_leaf$Distribution = ifelse(shapiro_res_combined_Dular_control_recovery12_FL_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_control_recovery12_FL_flag_leaf$Distribution == "Normal")     #65 metabolites
sum(shapiro_res_combined_Dular_control_recovery12_FL_flag_leaf$Distribution == "Non-normal") #16 metabolites

#wilcox test
wilcox_res_Dular_C_R12_flag_leaf_FL_combined = t(data.frame(lapply(data_combined_Dular_log2median_FL_flag_leaf[ , 3:ncol(data_combined_Dular_log2median_FL_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_FL_flag_leaf$Timepoint, subset = data_combined_Dular_log2median_FL_flag_leaf$Timepoint %in% c("Control", "12h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_C_R12_flag_leaf_FL_combined) = "p.value"
rownames(wilcox_res_Dular_C_R12_flag_leaf_FL_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_Dular_C_R12_flag_leaf_FL_combined = as.data.frame(wilcox_res_Dular_C_R12_flag_leaf_FL_combined)
wilcox_res_Dular_C_R12_flag_leaf_FL_combined$Significance = with(wilcox_res_Dular_C_R12_flag_leaf_FL_combined,
                                                                 ifelse(p.value <= 0.001, "***", ifelse(
                                                                   p.value <= 0.01, "**", ifelse(
                                                                     p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_C_R12_flag_leaf_FL_combined$Significance == "*") #11 metabolites
sum(wilcox_res_Dular_C_R12_flag_leaf_FL_combined$Significance == "**") #11 metabolites
sum(wilcox_res_Dular_C_R12_flag_leaf_FL_combined$Significance == "***") #19 metabolites
sum(wilcox_res_Dular_C_R12_flag_leaf_FL_combined$Significance != "ns") #41 metabolites
sum(wilcox_res_Dular_C_R12_flag_leaf_FL_combined$Significance == "ns") #40 metabolites


#Anjali

#data distribution
shapiro_res_combined_Anjali_control_recovery12_FL_flag_leaf = t(data.frame(lapply(data_combined_Anjali_log2median_FL_flag_leaf[which(data_combined_Anjali_log2median_FL_flag_leaf$Timepoint %in% c("Control", "12h Rewatering")), 3:ncol(data_combined_Anjali_log2median_FL_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_control_recovery12_FL_flag_leaf) = "p.value"
shapiro_res_combined_Anjali_control_recovery12_FL_flag_leaf = as.data.frame(shapiro_res_combined_Anjali_control_recovery12_FL_flag_leaf)
shapiro_res_combined_Anjali_control_recovery12_FL_flag_leaf$Distribution = ifelse(shapiro_res_combined_Anjali_control_recovery12_FL_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_control_recovery12_FL_flag_leaf$Distribution == "Normal")     #49 metabolites
sum(shapiro_res_combined_Anjali_control_recovery12_FL_flag_leaf$Distribution == "Non-normal") #32 metabolites

#wilcox test
wilcox_res_Anjali_C_R12_flag_leaf_FL_combined = t(data.frame(lapply(data_combined_Anjali_log2median_FL_flag_leaf[ , 3:ncol(data_combined_Anjali_log2median_FL_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_FL_flag_leaf$Timepoint, subset = data_combined_Anjali_log2median_FL_flag_leaf$Timepoint %in% c("Control", "12h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_C_R12_flag_leaf_FL_combined) = "p.value"
rownames(wilcox_res_Anjali_C_R12_flag_leaf_FL_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_Anjali_C_R12_flag_leaf_FL_combined = as.data.frame(wilcox_res_Anjali_C_R12_flag_leaf_FL_combined)
wilcox_res_Anjali_C_R12_flag_leaf_FL_combined$Significance = with(wilcox_res_Anjali_C_R12_flag_leaf_FL_combined,
                                                                  ifelse(p.value <= 0.001, "***", ifelse(
                                                                    p.value <= 0.01, "**", ifelse(
                                                                      p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_C_R12_flag_leaf_FL_combined$Significance == "*") #11 metabolites
sum(wilcox_res_Anjali_C_R12_flag_leaf_FL_combined$Significance == "**") #16 metabolites
sum(wilcox_res_Anjali_C_R12_flag_leaf_FL_combined$Significance == "***") #15 metabolites
sum(wilcox_res_Anjali_C_R12_flag_leaf_FL_combined$Significance != "ns") #42 metabolites
sum(wilcox_res_Anjali_C_R12_flag_leaf_FL_combined$Significance == "ns") #39 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_C_R12_flag_leaf_FL_combined = cbind(wilcox_res_N22_C_R12_flag_leaf_FL_combined, 
                                               wilcox_res_Dular_C_R12_flag_leaf_FL_combined$Significance, 
                                               wilcox_res_Anjali_C_R12_flag_leaf_FL_combined$Significance)
wilcox_sig_C_R12_flag_leaf_FL_combined = wilcox_sig_C_R12_flag_leaf_FL_combined[ , -1]
colnames(wilcox_sig_C_R12_flag_leaf_FL_combined)[1] = "wilcox_res_N22_C_R12_flag_leaf_FL_combined$Significance"
wilcox_sig_C_R12_flag_leaf_FL_combined = data.frame(lapply(wilcox_sig_C_R12_flag_leaf_FL_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_C_R12_flag_leaf_FL_combined) = rownames(wilcox_res_N22_C_R12_flag_leaf_FL_combined)
wilcox_sig_C_R12_flag_leaf_FL_combined$Nonsig = rowSums(wilcox_sig_C_R12_flag_leaf_FL_combined == " ")  #count number of nonsignificant
wilcox_sig_final_C_R12_flag_leaf_FL_combined = wilcox_sig_C_R12_flag_leaf_FL_combined[-which(wilcox_sig_C_R12_flag_leaf_FL_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 56 metabolites

wilcox_sig_log2fc_C_R12_flag_leaf_FL_combined = log2fc_combined_control_recovery_FL_flag_leaf[which(rownames(log2fc_combined_control_recovery_FL_flag_leaf) %in% rownames(wilcox_sig_final_C_R12_flag_leaf_FL_combined)), c(1, 4, 7)]  #log2FC values of metabolites significant in at least one of the cultivars

#overview - with HCA - Euclidean distance and average linkage 
png("wilcox_sig_log2fc_C_R12_flag_leaf_FL_combined.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_log2fc_C_R12_flag_leaf_FL_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(5, 20),  
          main = "Flag leaf-Combined-Wilcoxon \nLog2 FC, FL, 12h RW/Control", 
          cexCol = 1.5, 
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          labCol = c("N22", "Dular", "Anjali"), 
          srtCol = 0, adjCol = c(0.5, 0.5), 
          cellnote = wilcox_sig_final_C_R12_flag_leaf_FL_combined, 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 3)), 
          rowsep = c(seq(0, 56)), 
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))
dev.off()


####################

#common metabolites with significant difference - LS vs C and R12 vs C 
intersect(rownames(wilcox_sig_log2fc_C_LS_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_C_R12_flag_leaf_FL_combined))
#cultivar-specific - Anjali - 33 metabolites - example only
intersect(rownames(wilcox_res_Anjali_C_LS_flag_leaf_FL_combined)[wilcox_res_Anjali_C_LS_flag_leaf_FL_combined$Significance != "ns"], rownames(wilcox_res_Anjali_C_R12_flag_leaf_FL_combined)[wilcox_res_Anjali_C_R12_flag_leaf_FL_combined$Significance != "ns"])  


#metabolites with significant difference only in LS vs C
setdiff(rownames(wilcox_sig_log2fc_C_LS_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_C_R12_flag_leaf_FL_combined))
#cultivar-specific - Anjali - 11 metabolites - example only
setdiff(rownames(wilcox_res_Anjali_C_LS_flag_leaf_FL_combined)[wilcox_res_Anjali_C_LS_flag_leaf_FL_combined$Significance != "ns"], rownames(wilcox_res_Anjali_C_R12_flag_leaf_FL_combined)[wilcox_res_Anjali_C_R12_flag_leaf_FL_combined$Significance != "ns"])


#metabolites with significant difference only in R12 vs C
setdiff(rownames(wilcox_sig_log2fc_C_R12_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_C_LS_flag_leaf_FL_combined))
#cultivar-specific - Anjali - 9 metabolites - example only
setdiff(rownames(wilcox_res_Anjali_C_R12_flag_leaf_FL_combined)[wilcox_res_Anjali_C_R12_flag_leaf_FL_combined$Significance != "ns"], rownames(wilcox_res_Anjali_C_LS_flag_leaf_FL_combined)[wilcox_res_Anjali_C_LS_flag_leaf_FL_combined$Significance != "ns"])


####################


#venn diagram - common among cultivars

#N22
wilcox_sig_N22_C_R12_flag_leaf_FL_combined = wilcox_res_N22_C_R12_flag_leaf_FL_combined[-which(wilcox_res_N22_C_R12_flag_leaf_FL_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_N22_C_R12_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_C_R12_flag_leaf_FL_combined[match(rownames(wilcox_sig_N22_C_R12_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_C_R12_flag_leaf_FL_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_N22_C_R12_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_C_R12_flag_leaf_FL_combined)[1]
wilcox_sig_log2fc_N22_C_R12_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_log2fc_N22_C_R12_flag_leaf_FL_combined$`FNR12 - FNC` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_N22_C_R12_flag_leaf_FL_combined$Up.Down == "Up")    #16 metabolites
sum(wilcox_sig_log2fc_N22_C_R12_flag_leaf_FL_combined$Up.Down == "Down")  #14 metabolites

#Dular
wilcox_sig_Dular_C_R12_flag_leaf_FL_combined = wilcox_res_Dular_C_R12_flag_leaf_FL_combined[-which(wilcox_res_Dular_C_R12_flag_leaf_FL_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Dular_C_R12_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_C_R12_flag_leaf_FL_combined[match(rownames(wilcox_sig_Dular_C_R12_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_C_R12_flag_leaf_FL_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Dular_C_R12_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_C_R12_flag_leaf_FL_combined)[2]
wilcox_sig_log2fc_Dular_C_R12_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_log2fc_Dular_C_R12_flag_leaf_FL_combined$`FDR12 - FDC` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Dular_C_R12_flag_leaf_FL_combined$Up.Down == "Up")    #25 metabolites
sum(wilcox_sig_log2fc_Dular_C_R12_flag_leaf_FL_combined$Up.Down == "Down")  #16 metabolites

#Anjali
wilcox_sig_Anjali_C_R12_flag_leaf_FL_combined = wilcox_res_Anjali_C_R12_flag_leaf_FL_combined[-which(wilcox_res_Anjali_C_R12_flag_leaf_FL_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Anjali_C_R12_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_C_R12_flag_leaf_FL_combined[match(rownames(wilcox_sig_Anjali_C_R12_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_C_R12_flag_leaf_FL_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Anjali_C_R12_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_C_R12_flag_leaf_FL_combined)[3]
wilcox_sig_log2fc_Anjali_C_R12_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_log2fc_Anjali_C_R12_flag_leaf_FL_combined$`FAR12 - FAC` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Anjali_C_R12_flag_leaf_FL_combined$Up.Down == "Up")    #23 metabolites
sum(wilcox_sig_log2fc_Anjali_C_R12_flag_leaf_FL_combined$Up.Down == "Down")  #19 metabolites


#venn - increased - 12h rewatering/control
vennlist_up_wilcox_C_R12_flag_leaf_FL_combined = list(rownames(wilcox_sig_log2fc_N22_C_R12_flag_leaf_FL_combined)[wilcox_sig_log2fc_N22_C_R12_flag_leaf_FL_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Dular_C_R12_flag_leaf_FL_combined)[wilcox_sig_log2fc_Dular_C_R12_flag_leaf_FL_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Anjali_C_R12_flag_leaf_FL_combined)[wilcox_sig_log2fc_Anjali_C_R12_flag_leaf_FL_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_C_R12_flag_leaf_FL_combined) = c("N22", "Dular", "Anjali")
venn_up_wilcox_C_R12_flag_leaf_FL_combined = venn(vennlist_up_wilcox_C_R12_flag_leaf_FL_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_R12_control_flag_leaf_FL_combined = venn.diagram(vennlist_up_wilcox_C_R12_flag_leaf_FL_combined, category.names = c(paste(names(vennlist_up_wilcox_C_R12_flag_leaf_FL_combined[1]), " (", length(vennlist_up_wilcox_C_R12_flag_leaf_FL_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_C_R12_flag_leaf_FL_combined[2]), " (", length(vennlist_up_wilcox_C_R12_flag_leaf_FL_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_C_R12_flag_leaf_FL_combined[3]), " (", length(vennlist_up_wilcox_C_R12_flag_leaf_FL_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "A", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label A
intersect_wilcox_C_R12_venn_up_flag_leaf_FL_combined = attr(venn_up_wilcox_C_R12_flag_leaf_FL_combined, "intersections")  #list of unique and common metabolites
venn_up_wilcox_metabolites_C_R12_flag_leaf_FL_combined = t(ldply(intersect_wilcox_C_R12_venn_up_flag_leaf_FL_combined, rbind))  #convert from list to matrix
colnames(venn_up_wilcox_metabolites_C_R12_flag_leaf_FL_combined) = venn_up_wilcox_metabolites_C_R12_flag_leaf_FL_combined[1, ]
venn_up_wilcox_metabolites_C_R12_flag_leaf_FL_combined = venn_up_wilcox_metabolites_C_R12_flag_leaf_FL_combined[-1, ]


#venn - decreased - 12h rewatering/control
vennlist_down_wilcox_C_R12_flag_leaf_FL_combined = list(rownames(wilcox_sig_log2fc_N22_C_R12_flag_leaf_FL_combined)[wilcox_sig_log2fc_N22_C_R12_flag_leaf_FL_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Dular_C_R12_flag_leaf_FL_combined)[wilcox_sig_log2fc_Dular_C_R12_flag_leaf_FL_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Anjali_C_R12_flag_leaf_FL_combined)[wilcox_sig_log2fc_Anjali_C_R12_flag_leaf_FL_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_C_R12_flag_leaf_FL_combined) = c("N22", "Dular", "Anjali")
venn_down_wilcox_C_R12_flag_leaf_FL_combined = venn(vennlist_down_wilcox_C_R12_flag_leaf_FL_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_R12_control_flag_leaf_FL_combined = venn.diagram(vennlist_down_wilcox_C_R12_flag_leaf_FL_combined, category.names = c(paste(names(vennlist_down_wilcox_C_R12_flag_leaf_FL_combined[1]), " (", length(vennlist_down_wilcox_C_R12_flag_leaf_FL_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_C_R12_flag_leaf_FL_combined[2]), " (", length(vennlist_down_wilcox_C_R12_flag_leaf_FL_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_C_R12_flag_leaf_FL_combined[3]), " (", length(vennlist_down_wilcox_C_R12_flag_leaf_FL_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "B", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label B
intersect_wilcox_C_R12_venn_down_flag_leaf_FL_combined = attr(venn_down_wilcox_C_R12_flag_leaf_FL_combined, "intersections")  #list of unique and common metabolites
venn_down_wilcox_metabolites_C_R12_flag_leaf_FL_combined = t(ldply(intersect_wilcox_C_R12_venn_down_flag_leaf_FL_combined, rbind))  #convert from list to matrix
colnames(venn_down_wilcox_metabolites_C_R12_flag_leaf_FL_combined) = venn_down_wilcox_metabolites_C_R12_flag_leaf_FL_combined[1, ]
venn_down_wilcox_metabolites_C_R12_flag_leaf_FL_combined = venn_down_wilcox_metabolites_C_R12_flag_leaf_FL_combined[-1, ]


write.table(venn_up_wilcox_metabolites_C_R12_flag_leaf_FL_combined, file = "wilcox_sig_up_metabolites_flag_leaf_FL_R12_control_combined.txt", sep = "\t", quote = F)
write.table(venn_down_wilcox_metabolites_C_R12_flag_leaf_FL_combined, file = "wilcox_sig_down_metabolites_flag_leaf_FL_R12_control_combined.txt", sep = "\t", quote = F)


####################


#Control vs. 36h rewatering

#N22

#data distribution
shapiro_res_combined_N22_control_recovery36_FL_flag_leaf = t(data.frame(lapply(data_combined_N22_log2median_FL_flag_leaf[which(data_combined_N22_log2median_FL_flag_leaf$Timepoint %in% c("Control", "36h Rewatering")), 3:ncol(data_combined_N22_log2median_FL_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_control_recovery36_FL_flag_leaf) = "p.value"
shapiro_res_combined_N22_control_recovery36_FL_flag_leaf = as.data.frame(shapiro_res_combined_N22_control_recovery36_FL_flag_leaf)
shapiro_res_combined_N22_control_recovery36_FL_flag_leaf$Distribution = ifelse(shapiro_res_combined_N22_control_recovery36_FL_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_control_recovery36_FL_flag_leaf$Distribution == "Normal")     #58 metabolites
sum(shapiro_res_combined_N22_control_recovery36_FL_flag_leaf$Distribution == "Non-normal") #23 metabolites

#wilcox test
wilcox_res_N22_C_R36_flag_leaf_FL_combined = t(data.frame(lapply(data_combined_N22_log2median_FL_flag_leaf[ , 3:ncol(data_combined_N22_log2median_FL_flag_leaf)], function(x) wilcox.test(x ~ data_combined_N22_log2median_FL_flag_leaf$Timepoint, subset = data_combined_N22_log2median_FL_flag_leaf$Timepoint %in% c("Control", "36h Rewatering"))$p.value)))
colnames(wilcox_res_N22_C_R36_flag_leaf_FL_combined) = "p.value"
rownames(wilcox_res_N22_C_R36_flag_leaf_FL_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_N22_C_R36_flag_leaf_FL_combined = as.data.frame(wilcox_res_N22_C_R36_flag_leaf_FL_combined)
wilcox_res_N22_C_R36_flag_leaf_FL_combined$Significance = with(wilcox_res_N22_C_R36_flag_leaf_FL_combined,
                                                               ifelse(p.value <= 0.001, "***", ifelse(
                                                                 p.value <= 0.01, "**", ifelse(
                                                                   p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_C_R36_flag_leaf_FL_combined$Significance == "*") #12 metabolites
sum(wilcox_res_N22_C_R36_flag_leaf_FL_combined$Significance == "**") #3 metabolites
sum(wilcox_res_N22_C_R36_flag_leaf_FL_combined$Significance == "***") #10 metabolites
sum(wilcox_res_N22_C_R36_flag_leaf_FL_combined$Significance != "ns") #25 metabolites
sum(wilcox_res_N22_C_R36_flag_leaf_FL_combined$Significance == "ns") #56 metabolites


#Dular

#data distribution
shapiro_res_combined_Dular_control_recovery36_FL_flag_leaf = t(data.frame(lapply(data_combined_Dular_log2median_FL_flag_leaf[which(data_combined_Dular_log2median_FL_flag_leaf$Timepoint %in% c("Control", "36h Rewatering")), 3:ncol(data_combined_Dular_log2median_FL_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_control_recovery36_FL_flag_leaf) = "p.value"
shapiro_res_combined_Dular_control_recovery36_FL_flag_leaf = as.data.frame(shapiro_res_combined_Dular_control_recovery36_FL_flag_leaf)
shapiro_res_combined_Dular_control_recovery36_FL_flag_leaf$Distribution = ifelse(shapiro_res_combined_Dular_control_recovery36_FL_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_control_recovery36_FL_flag_leaf$Distribution == "Normal")     #58 metabolites
sum(shapiro_res_combined_Dular_control_recovery36_FL_flag_leaf$Distribution == "Non-normal") #23 metabolites

#wilcox test
wilcox_res_Dular_C_R36_flag_leaf_FL_combined = t(data.frame(lapply(data_combined_Dular_log2median_FL_flag_leaf[ , 3:ncol(data_combined_Dular_log2median_FL_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_FL_flag_leaf$Timepoint, subset = data_combined_Dular_log2median_FL_flag_leaf$Timepoint %in% c("Control", "36h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_C_R36_flag_leaf_FL_combined) = "p.value"
rownames(wilcox_res_Dular_C_R36_flag_leaf_FL_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_Dular_C_R36_flag_leaf_FL_combined = as.data.frame(wilcox_res_Dular_C_R36_flag_leaf_FL_combined)
wilcox_res_Dular_C_R36_flag_leaf_FL_combined$Significance = with(wilcox_res_Dular_C_R36_flag_leaf_FL_combined,
                                                                 ifelse(p.value <= 0.001, "***", ifelse(
                                                                   p.value <= 0.01, "**", ifelse(
                                                                     p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_C_R36_flag_leaf_FL_combined$Significance == "*") #11 metabolites
sum(wilcox_res_Dular_C_R36_flag_leaf_FL_combined$Significance == "**") #8 metabolites
sum(wilcox_res_Dular_C_R36_flag_leaf_FL_combined$Significance == "***") #9 metabolites
sum(wilcox_res_Dular_C_R36_flag_leaf_FL_combined$Significance != "ns") #28 metabolites
sum(wilcox_res_Dular_C_R36_flag_leaf_FL_combined$Significance == "ns") #53 metabolites


#Anjali

#data distribution
shapiro_res_combined_Anjali_control_recovery36_FL_flag_leaf = t(data.frame(lapply(data_combined_Anjali_log2median_FL_flag_leaf[which(data_combined_Anjali_log2median_FL_flag_leaf$Timepoint %in% c("Control", "36h Rewatering")), 3:ncol(data_combined_Anjali_log2median_FL_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_control_recovery36_FL_flag_leaf) = "p.value"
shapiro_res_combined_Anjali_control_recovery36_FL_flag_leaf = as.data.frame(shapiro_res_combined_Anjali_control_recovery36_FL_flag_leaf)
shapiro_res_combined_Anjali_control_recovery36_FL_flag_leaf$Distribution = ifelse(shapiro_res_combined_Anjali_control_recovery36_FL_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_control_recovery36_FL_flag_leaf$Distribution == "Normal")     #56 metabolites
sum(shapiro_res_combined_Anjali_control_recovery36_FL_flag_leaf$Distribution == "Non-normal") #25 metabolites

#wilcox test
wilcox_res_Anjali_C_R36_flag_leaf_FL_combined = t(data.frame(lapply(data_combined_Anjali_log2median_FL_flag_leaf[ , 3:ncol(data_combined_Anjali_log2median_FL_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_FL_flag_leaf$Timepoint, subset = data_combined_Anjali_log2median_FL_flag_leaf$Timepoint %in% c("Control", "36h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_C_R36_flag_leaf_FL_combined) = "p.value"
rownames(wilcox_res_Anjali_C_R36_flag_leaf_FL_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_Anjali_C_R36_flag_leaf_FL_combined = as.data.frame(wilcox_res_Anjali_C_R36_flag_leaf_FL_combined)
wilcox_res_Anjali_C_R36_flag_leaf_FL_combined$Significance = with(wilcox_res_Anjali_C_R36_flag_leaf_FL_combined,
                                                                  ifelse(p.value <= 0.001, "***", ifelse(
                                                                    p.value <= 0.01, "**", ifelse(
                                                                      p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_C_R36_flag_leaf_FL_combined$Significance == "*") #8 metabolites
sum(wilcox_res_Anjali_C_R36_flag_leaf_FL_combined$Significance == "**") #9 metabolites
sum(wilcox_res_Anjali_C_R36_flag_leaf_FL_combined$Significance == "***") #12 metabolites
sum(wilcox_res_Anjali_C_R36_flag_leaf_FL_combined$Significance != "ns") #29 metabolites
sum(wilcox_res_Anjali_C_R36_flag_leaf_FL_combined$Significance == "ns") #52 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_C_R36_flag_leaf_FL_combined = cbind(wilcox_res_N22_C_R36_flag_leaf_FL_combined, 
                                               wilcox_res_Dular_C_R36_flag_leaf_FL_combined$Significance, 
                                               wilcox_res_Anjali_C_R36_flag_leaf_FL_combined$Significance)
wilcox_sig_C_R36_flag_leaf_FL_combined = wilcox_sig_C_R36_flag_leaf_FL_combined[ , -1]
colnames(wilcox_sig_C_R36_flag_leaf_FL_combined)[1] = "wilcox_res_N22_C_R36_flag_leaf_FL_combined$Significance"
wilcox_sig_C_R36_flag_leaf_FL_combined = data.frame(lapply(wilcox_sig_C_R36_flag_leaf_FL_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_C_R36_flag_leaf_FL_combined) = rownames(wilcox_res_N22_C_R36_flag_leaf_FL_combined)
wilcox_sig_C_R36_flag_leaf_FL_combined$Nonsig = rowSums(wilcox_sig_C_R36_flag_leaf_FL_combined == " ")  #count number of nonsignificant
wilcox_sig_final_C_R36_flag_leaf_FL_combined = wilcox_sig_C_R36_flag_leaf_FL_combined[-which(wilcox_sig_C_R36_flag_leaf_FL_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 45 metabolites

wilcox_sig_log2fc_C_R36_flag_leaf_FL_combined = log2fc_combined_control_recovery_FL_flag_leaf[which(rownames(log2fc_combined_control_recovery_FL_flag_leaf) %in% rownames(wilcox_sig_final_C_R36_flag_leaf_FL_combined)), c(2, 5, 8)]  #log2FC values of metabolites significant in at least one of the cultivars

#overview - with HCA - Euclidean distance and average linkage 
png("wilcox_sig_log2fc_C_R36_flag_leaf_FL_combined.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_log2fc_C_R36_flag_leaf_FL_combined,
          Rowv = TRUE,
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(5, 20), 
          main = "Flag leaf-Combined-Wilcoxon \nLog2 FC, FL, 36h RW/Control", 
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
          cellnote = wilcox_sig_final_C_R36_flag_leaf_FL_combined, 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 3)), 
          rowsep = c(seq(0, 45)), 
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))
dev.off()


####################

#common metabolites with significant difference - R12 vs C and R36 vs C 
intersect(rownames(wilcox_sig_log2fc_C_R12_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_C_R36_flag_leaf_FL_combined))

#significant only in R12 vs C
setdiff(rownames(wilcox_sig_log2fc_C_R12_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_C_R36_flag_leaf_FL_combined))

#significant only in R36 vs C
setdiff(rownames(wilcox_sig_log2fc_C_R36_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_C_R12_flag_leaf_FL_combined))

####################


#venn diagram - common among cultivars

#N22
wilcox_sig_N22_C_R36_flag_leaf_FL_combined = wilcox_res_N22_C_R36_flag_leaf_FL_combined[-which(wilcox_res_N22_C_R36_flag_leaf_FL_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_N22_C_R36_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_C_R36_flag_leaf_FL_combined[match(rownames(wilcox_sig_N22_C_R36_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_C_R36_flag_leaf_FL_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_N22_C_R36_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_C_R36_flag_leaf_FL_combined)[1]
wilcox_sig_log2fc_N22_C_R36_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_log2fc_N22_C_R36_flag_leaf_FL_combined$`FNR36 - FNC` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_N22_C_R36_flag_leaf_FL_combined$Up.Down == "Up")    #10 metabolites
sum(wilcox_sig_log2fc_N22_C_R36_flag_leaf_FL_combined$Up.Down == "Down")  #15 metabolites

#Dular
wilcox_sig_Dular_C_R36_flag_leaf_FL_combined = wilcox_res_Dular_C_R36_flag_leaf_FL_combined[-which(wilcox_res_Dular_C_R36_flag_leaf_FL_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Dular_C_R36_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_C_R36_flag_leaf_FL_combined[match(rownames(wilcox_sig_Dular_C_R36_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_C_R36_flag_leaf_FL_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Dular_C_R36_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_C_R36_flag_leaf_FL_combined)[2]
wilcox_sig_log2fc_Dular_C_R36_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_log2fc_Dular_C_R36_flag_leaf_FL_combined$`FDR36 - FDC` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Dular_C_R36_flag_leaf_FL_combined$Up.Down == "Up")    #15 metabolites
sum(wilcox_sig_log2fc_Dular_C_R36_flag_leaf_FL_combined$Up.Down == "Down")  #13 metabolites

#Anjali
wilcox_sig_Anjali_C_R36_flag_leaf_FL_combined = wilcox_res_Anjali_C_R36_flag_leaf_FL_combined[-which(wilcox_res_Anjali_C_R36_flag_leaf_FL_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Anjali_C_R36_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_C_R36_flag_leaf_FL_combined[match(rownames(wilcox_sig_Anjali_C_R36_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_C_R36_flag_leaf_FL_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Anjali_C_R36_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_C_R36_flag_leaf_FL_combined)[3]
wilcox_sig_log2fc_Anjali_C_R36_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_log2fc_Anjali_C_R36_flag_leaf_FL_combined$`FAR36 - FAC` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Anjali_C_R36_flag_leaf_FL_combined$Up.Down == "Up")    #16 metabolites
sum(wilcox_sig_log2fc_Anjali_C_R36_flag_leaf_FL_combined$Up.Down == "Down")  #13 metabolites


#venn - increased - 36h rewatering/control
vennlist_up_wilcox_C_R36_flag_leaf_FL_combined = list(rownames(wilcox_sig_log2fc_N22_C_R36_flag_leaf_FL_combined)[wilcox_sig_log2fc_N22_C_R36_flag_leaf_FL_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Dular_C_R36_flag_leaf_FL_combined)[wilcox_sig_log2fc_Dular_C_R36_flag_leaf_FL_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Anjali_C_R36_flag_leaf_FL_combined)[wilcox_sig_log2fc_Anjali_C_R36_flag_leaf_FL_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_C_R36_flag_leaf_FL_combined) = c("N22", "Dular", "Anjali")
venn_up_wilcox_C_R36_flag_leaf_FL_combined = venn(vennlist_up_wilcox_C_R36_flag_leaf_FL_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_R36_control_flag_leaf_FL_combined = venn.diagram(vennlist_up_wilcox_C_R36_flag_leaf_FL_combined, category.names = c(paste(names(vennlist_up_wilcox_C_R36_flag_leaf_FL_combined[1]), " (", length(vennlist_up_wilcox_C_R36_flag_leaf_FL_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_C_R36_flag_leaf_FL_combined[2]), " (", length(vennlist_up_wilcox_C_R36_flag_leaf_FL_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_C_R36_flag_leaf_FL_combined[3]), " (", length(vennlist_up_wilcox_C_R36_flag_leaf_FL_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "C", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label C
intersect_wilcox_C_R36_venn_up_flag_leaf_FL_combined = attr(venn_up_wilcox_C_R36_flag_leaf_FL_combined, "intersections")  #list of unique and common metabolites
venn_up_wilcox_metabolites_C_R36_flag_leaf_FL_combined = t(ldply(intersect_wilcox_C_R36_venn_up_flag_leaf_FL_combined, rbind))  #convert from list to matrix
colnames(venn_up_wilcox_metabolites_C_R36_flag_leaf_FL_combined) = venn_up_wilcox_metabolites_C_R36_flag_leaf_FL_combined[1, ]
venn_up_wilcox_metabolites_C_R36_flag_leaf_FL_combined = venn_up_wilcox_metabolites_C_R36_flag_leaf_FL_combined[-1, ]


#venn - decreased - 36h rewatering/control
vennlist_down_wilcox_C_R36_flag_leaf_FL_combined = list(rownames(wilcox_sig_log2fc_N22_C_R36_flag_leaf_FL_combined)[wilcox_sig_log2fc_N22_C_R36_flag_leaf_FL_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Dular_C_R36_flag_leaf_FL_combined)[wilcox_sig_log2fc_Dular_C_R36_flag_leaf_FL_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Anjali_C_R36_flag_leaf_FL_combined)[wilcox_sig_log2fc_Anjali_C_R36_flag_leaf_FL_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_C_R36_flag_leaf_FL_combined) = c("N22", "Dular", "Anjali")
venn_down_wilcox_C_R36_flag_leaf_FL_combined = venn(vennlist_down_wilcox_C_R36_flag_leaf_FL_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_R36_control_flag_leaf_FL_combined = venn.diagram(vennlist_down_wilcox_C_R36_flag_leaf_FL_combined, category.names = c(paste(names(vennlist_down_wilcox_C_R36_flag_leaf_FL_combined[1]), " (", length(vennlist_down_wilcox_C_R36_flag_leaf_FL_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_C_R36_flag_leaf_FL_combined[2]), " (", length(vennlist_down_wilcox_C_R36_flag_leaf_FL_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_C_R36_flag_leaf_FL_combined[3]), " (", length(vennlist_down_wilcox_C_R36_flag_leaf_FL_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "D", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label D
intersect_wilcox_C_R36_venn_down_flag_leaf_FL_combined = attr(venn_down_wilcox_C_R36_flag_leaf_FL_combined, "intersections")  #list of unique and common metabolites
venn_down_wilcox_metabolites_C_R36_flag_leaf_FL_combined = t(ldply(intersect_wilcox_C_R36_venn_down_flag_leaf_FL_combined, rbind))  #convert from list to matrix
colnames(venn_down_wilcox_metabolites_C_R36_flag_leaf_FL_combined) = venn_down_wilcox_metabolites_C_R36_flag_leaf_FL_combined[1, ]
venn_down_wilcox_metabolites_C_R36_flag_leaf_FL_combined = venn_down_wilcox_metabolites_C_R36_flag_leaf_FL_combined[-1, ]


write.table(venn_up_wilcox_metabolites_C_R36_flag_leaf_FL_combined, file = "wilcox_sig_up_metabolites_flag_leaf_FL_R36_control_combined.txt", sep = "\t", quote = F)
write.table(venn_down_wilcox_metabolites_C_R36_flag_leaf_FL_combined, file = "wilcox_sig_down_metabolites_flag_leaf_FL_R36_control_combined.txt", sep = "\t", quote = F)


####################


#Control vs. 60h rewatering

#N22

#data distribution
shapiro_res_combined_N22_control_recovery60_FL_flag_leaf = t(data.frame(lapply(data_combined_N22_log2median_FL_flag_leaf[which(data_combined_N22_log2median_FL_flag_leaf$Timepoint %in% c("Control", "60h Rewatering")), 3:ncol(data_combined_N22_log2median_FL_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_control_recovery60_FL_flag_leaf) = "p.value"
shapiro_res_combined_N22_control_recovery60_FL_flag_leaf = as.data.frame(shapiro_res_combined_N22_control_recovery60_FL_flag_leaf)
shapiro_res_combined_N22_control_recovery60_FL_flag_leaf$Distribution = ifelse(shapiro_res_combined_N22_control_recovery60_FL_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_control_recovery60_FL_flag_leaf$Distribution == "Normal")     #64 metabolites
sum(shapiro_res_combined_N22_control_recovery60_FL_flag_leaf$Distribution == "Non-normal") #17 metabolites

#wilcox test
wilcox_res_N22_C_R60_flag_leaf_FL_combined = t(data.frame(lapply(data_combined_N22_log2median_FL_flag_leaf[ , 3:ncol(data_combined_N22_log2median_FL_flag_leaf)], function(x) wilcox.test(x ~ data_combined_N22_log2median_FL_flag_leaf$Timepoint, subset = data_combined_N22_log2median_FL_flag_leaf$Timepoint %in% c("Control", "60h Rewatering"))$p.value)))
colnames(wilcox_res_N22_C_R60_flag_leaf_FL_combined) = "p.value"
rownames(wilcox_res_N22_C_R60_flag_leaf_FL_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_N22_C_R60_flag_leaf_FL_combined = as.data.frame(wilcox_res_N22_C_R60_flag_leaf_FL_combined)
wilcox_res_N22_C_R60_flag_leaf_FL_combined$Significance = with(wilcox_res_N22_C_R60_flag_leaf_FL_combined,
                                                               ifelse(p.value <= 0.001, "***", ifelse(
                                                                 p.value <= 0.01, "**", ifelse(
                                                                   p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_C_R60_flag_leaf_FL_combined$Significance == "*") #15 metabolites
sum(wilcox_res_N22_C_R60_flag_leaf_FL_combined$Significance == "**") #4 metabolites
sum(wilcox_res_N22_C_R60_flag_leaf_FL_combined$Significance == "***") #12 metabolites
sum(wilcox_res_N22_C_R60_flag_leaf_FL_combined$Significance != "ns") #31 metabolites
sum(wilcox_res_N22_C_R60_flag_leaf_FL_combined$Significance == "ns") #50 metabolites


#Dular

#data distribution
shapiro_res_combined_Dular_control_recovery60_FL_flag_leaf = t(data.frame(lapply(data_combined_Dular_log2median_FL_flag_leaf[which(data_combined_Dular_log2median_FL_flag_leaf$Timepoint %in% c("Control", "60h Rewatering")), 3:ncol(data_combined_Dular_log2median_FL_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_control_recovery60_FL_flag_leaf) = "p.value"
shapiro_res_combined_Dular_control_recovery60_FL_flag_leaf = as.data.frame(shapiro_res_combined_Dular_control_recovery60_FL_flag_leaf)
shapiro_res_combined_Dular_control_recovery60_FL_flag_leaf$Distribution = ifelse(shapiro_res_combined_Dular_control_recovery60_FL_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_control_recovery60_FL_flag_leaf$Distribution == "Normal")     #58 metabolites
sum(shapiro_res_combined_Dular_control_recovery60_FL_flag_leaf$Distribution == "Non-normal") #23 metabolites

#wilcox test
wilcox_res_Dular_C_R60_flag_leaf_FL_combined = t(data.frame(lapply(data_combined_Dular_log2median_FL_flag_leaf[ , 3:ncol(data_combined_Dular_log2median_FL_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_FL_flag_leaf$Timepoint, subset = data_combined_Dular_log2median_FL_flag_leaf$Timepoint %in% c("Control", "60h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_C_R60_flag_leaf_FL_combined) = "p.value"
rownames(wilcox_res_Dular_C_R60_flag_leaf_FL_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_Dular_C_R60_flag_leaf_FL_combined = as.data.frame(wilcox_res_Dular_C_R60_flag_leaf_FL_combined)
wilcox_res_Dular_C_R60_flag_leaf_FL_combined$Significance = with(wilcox_res_Dular_C_R60_flag_leaf_FL_combined,
                                                                 ifelse(p.value <= 0.001, "***", ifelse(
                                                                   p.value <= 0.01, "**", ifelse(
                                                                     p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_C_R60_flag_leaf_FL_combined$Significance == "*") #13 metabolites
sum(wilcox_res_Dular_C_R60_flag_leaf_FL_combined$Significance == "**") #8 metabolites
sum(wilcox_res_Dular_C_R60_flag_leaf_FL_combined$Significance == "***") #4 metabolites
sum(wilcox_res_Dular_C_R60_flag_leaf_FL_combined$Significance != "ns") #25 metabolites
sum(wilcox_res_Dular_C_R60_flag_leaf_FL_combined$Significance == "ns") #56 metabolites


#Anjali

#data distribution
shapiro_res_combined_Anjali_control_recovery60_FL_flag_leaf = t(data.frame(lapply(data_combined_Anjali_log2median_FL_flag_leaf[which(data_combined_Anjali_log2median_FL_flag_leaf$Timepoint %in% c("Control", "60h Rewatering")), 3:ncol(data_combined_Anjali_log2median_FL_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_control_recovery60_FL_flag_leaf) = "p.value"
shapiro_res_combined_Anjali_control_recovery60_FL_flag_leaf = as.data.frame(shapiro_res_combined_Anjali_control_recovery60_FL_flag_leaf)
shapiro_res_combined_Anjali_control_recovery60_FL_flag_leaf$Distribution = ifelse(shapiro_res_combined_Anjali_control_recovery60_FL_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_control_recovery60_FL_flag_leaf$Distribution == "Normal")     #64 metabolites
sum(shapiro_res_combined_Anjali_control_recovery60_FL_flag_leaf$Distribution == "Non-normal") #17 metabolites

#wilcox test
wilcox_res_Anjali_C_R60_flag_leaf_FL_combined = t(data.frame(lapply(data_combined_Anjali_log2median_FL_flag_leaf[ , 3:ncol(data_combined_Anjali_log2median_FL_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_FL_flag_leaf$Timepoint, subset = data_combined_Anjali_log2median_FL_flag_leaf$Timepoint %in% c("Control", "60h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_C_R60_flag_leaf_FL_combined) = "p.value"
rownames(wilcox_res_Anjali_C_R60_flag_leaf_FL_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_Anjali_C_R60_flag_leaf_FL_combined = as.data.frame(wilcox_res_Anjali_C_R60_flag_leaf_FL_combined)
wilcox_res_Anjali_C_R60_flag_leaf_FL_combined$Significance = with(wilcox_res_Anjali_C_R60_flag_leaf_FL_combined,
                                                                  ifelse(p.value <= 0.001, "***", ifelse(
                                                                    p.value <= 0.01, "**", ifelse(
                                                                      p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_C_R60_flag_leaf_FL_combined$Significance == "*") #11 metabolites
sum(wilcox_res_Anjali_C_R60_flag_leaf_FL_combined$Significance == "**") #7 metabolites
sum(wilcox_res_Anjali_C_R60_flag_leaf_FL_combined$Significance == "***") #11 metabolites
sum(wilcox_res_Anjali_C_R60_flag_leaf_FL_combined$Significance != "ns") #29 metabolites
sum(wilcox_res_Anjali_C_R60_flag_leaf_FL_combined$Significance == "ns") #52 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_C_R60_flag_leaf_FL_combined = cbind(wilcox_res_N22_C_R60_flag_leaf_FL_combined, 
                                               wilcox_res_Dular_C_R60_flag_leaf_FL_combined$Significance, 
                                               wilcox_res_Anjali_C_R60_flag_leaf_FL_combined$Significance)
wilcox_sig_C_R60_flag_leaf_FL_combined = wilcox_sig_C_R60_flag_leaf_FL_combined[ , -1]
colnames(wilcox_sig_C_R60_flag_leaf_FL_combined)[1] = "wilcox_res_N22_C_R60_flag_leaf_FL_combined$Significance"
wilcox_sig_C_R60_flag_leaf_FL_combined = data.frame(lapply(wilcox_sig_C_R60_flag_leaf_FL_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_C_R60_flag_leaf_FL_combined) = rownames(wilcox_res_N22_C_R60_flag_leaf_FL_combined)
wilcox_sig_C_R60_flag_leaf_FL_combined$Nonsig = rowSums(wilcox_sig_C_R60_flag_leaf_FL_combined == " ")  #count number of nonsignificant
wilcox_sig_final_C_R60_flag_leaf_FL_combined = wilcox_sig_C_R60_flag_leaf_FL_combined[-which(wilcox_sig_C_R60_flag_leaf_FL_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 50 metabolites

wilcox_sig_log2fc_C_R60_flag_leaf_FL_combined = log2fc_combined_control_recovery_FL_flag_leaf[which(rownames(log2fc_combined_control_recovery_FL_flag_leaf) %in% rownames(wilcox_sig_final_C_R60_flag_leaf_FL_combined)), c(3, 6, 9)]  #log2FC values of metabolites significant in at least one of the cultivars

#overview - with HCA - Euclidean distance and average linkage 
png("wilcox_sig_log2fc_C_R60_flag_leaf_FL_combined.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_log2fc_C_R60_flag_leaf_FL_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)),
          margins = c(5, 20), 
          main = "Flag leaf-Combined-Wilcoxon \nLog2 FC, FL, 60h RW/Control", 
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
          cellnote = wilcox_sig_final_C_R60_flag_leaf_FL_combined, 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 3)), 
          rowsep = c(seq(0, 50)), 
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))
dev.off()


####################

#common metabolites with significant difference - R36 vs C and R60 vs C 
intersect(rownames(wilcox_sig_log2fc_C_R36_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_C_R60_flag_leaf_FL_combined))

#significant only in R36 vs C
setdiff(rownames(wilcox_sig_log2fc_C_R36_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_C_R60_flag_leaf_FL_combined))

#significant only in R60 vs C
setdiff(rownames(wilcox_sig_log2fc_C_R60_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_C_R36_flag_leaf_FL_combined))

####################


#venn diagram - common among cultivars

#N22
wilcox_sig_N22_C_R60_flag_leaf_FL_combined = wilcox_res_N22_C_R60_flag_leaf_FL_combined[-which(wilcox_res_N22_C_R60_flag_leaf_FL_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_N22_C_R60_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_C_R60_flag_leaf_FL_combined[match(rownames(wilcox_sig_N22_C_R60_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_C_R60_flag_leaf_FL_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_N22_C_R60_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_C_R60_flag_leaf_FL_combined)[1]
wilcox_sig_log2fc_N22_C_R60_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_log2fc_N22_C_R60_flag_leaf_FL_combined$`FNR60 - FNC` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_N22_C_R60_flag_leaf_FL_combined$Up.Down == "Up")    #13 metabolites
sum(wilcox_sig_log2fc_N22_C_R60_flag_leaf_FL_combined$Up.Down == "Down")  #18 metabolites

#Dular
wilcox_sig_Dular_C_R60_flag_leaf_FL_combined = wilcox_res_Dular_C_R60_flag_leaf_FL_combined[-which(wilcox_res_Dular_C_R60_flag_leaf_FL_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Dular_C_R60_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_C_R60_flag_leaf_FL_combined[match(rownames(wilcox_sig_Dular_C_R60_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_C_R60_flag_leaf_FL_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Dular_C_R60_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_C_R60_flag_leaf_FL_combined)[2]
wilcox_sig_log2fc_Dular_C_R60_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_log2fc_Dular_C_R60_flag_leaf_FL_combined$`FDR60 - FDC` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Dular_C_R60_flag_leaf_FL_combined$Up.Down == "Up")    #10 metabolites
sum(wilcox_sig_log2fc_Dular_C_R60_flag_leaf_FL_combined$Up.Down == "Down")  #15 metabolites

#Anjali
wilcox_sig_Anjali_C_R60_flag_leaf_FL_combined = wilcox_res_Anjali_C_R60_flag_leaf_FL_combined[-which(wilcox_res_Anjali_C_R60_flag_leaf_FL_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Anjali_C_R60_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_C_R60_flag_leaf_FL_combined[match(rownames(wilcox_sig_Anjali_C_R60_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_C_R60_flag_leaf_FL_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Anjali_C_R60_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_C_R60_flag_leaf_FL_combined)[3]
wilcox_sig_log2fc_Anjali_C_R60_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_log2fc_Anjali_C_R60_flag_leaf_FL_combined$`FAR60 - FAC` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Anjali_C_R60_flag_leaf_FL_combined$Up.Down == "Up")    #14 metabolites
sum(wilcox_sig_log2fc_Anjali_C_R60_flag_leaf_FL_combined$Up.Down == "Down")  #15 metabolites


#venn - increased - 60h rewatering/control
vennlist_up_wilcox_C_R60_flag_leaf_FL_combined = list(rownames(wilcox_sig_log2fc_N22_C_R60_flag_leaf_FL_combined)[wilcox_sig_log2fc_N22_C_R60_flag_leaf_FL_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Dular_C_R60_flag_leaf_FL_combined)[wilcox_sig_log2fc_Dular_C_R60_flag_leaf_FL_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Anjali_C_R60_flag_leaf_FL_combined)[wilcox_sig_log2fc_Anjali_C_R60_flag_leaf_FL_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_C_R60_flag_leaf_FL_combined) = c("N22", "Dular", "Anjali")
venn_up_wilcox_C_R60_flag_leaf_FL_combined = venn(vennlist_up_wilcox_C_R60_flag_leaf_FL_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_R60_control_flag_leaf_FL_combined = venn.diagram(vennlist_up_wilcox_C_R60_flag_leaf_FL_combined, category.names = c(paste(names(vennlist_up_wilcox_C_R60_flag_leaf_FL_combined[1]), " (", length(vennlist_up_wilcox_C_R60_flag_leaf_FL_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_C_R60_flag_leaf_FL_combined[2]), " (", length(vennlist_up_wilcox_C_R60_flag_leaf_FL_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_C_R60_flag_leaf_FL_combined[3]), " (", length(vennlist_up_wilcox_C_R60_flag_leaf_FL_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "E", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label E
intersect_wilcox_C_R60_venn_up_flag_leaf_FL_combined = attr(venn_up_wilcox_C_R60_flag_leaf_FL_combined, "intersections")  #list of unique and common metabolites
venn_up_wilcox_metabolites_C_R60_flag_leaf_FL_combined = t(ldply(intersect_wilcox_C_R60_venn_up_flag_leaf_FL_combined, rbind))  #convert from list to matrix
colnames(venn_up_wilcox_metabolites_C_R60_flag_leaf_FL_combined) = venn_up_wilcox_metabolites_C_R60_flag_leaf_FL_combined[1, ]
venn_up_wilcox_metabolites_C_R60_flag_leaf_FL_combined = venn_up_wilcox_metabolites_C_R60_flag_leaf_FL_combined[-1, ]


#venn - decreased - 60h rewatering/control
vennlist_down_wilcox_C_R60_flag_leaf_FL_combined = list(rownames(wilcox_sig_log2fc_N22_C_R60_flag_leaf_FL_combined)[wilcox_sig_log2fc_N22_C_R60_flag_leaf_FL_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Dular_C_R60_flag_leaf_FL_combined)[wilcox_sig_log2fc_Dular_C_R60_flag_leaf_FL_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Anjali_C_R60_flag_leaf_FL_combined)[wilcox_sig_log2fc_Anjali_C_R60_flag_leaf_FL_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_C_R60_flag_leaf_FL_combined) = c("N22", "Dular", "Anjali")
venn_down_wilcox_C_R60_flag_leaf_FL_combined = venn(vennlist_down_wilcox_C_R60_flag_leaf_FL_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_R60_control_flag_leaf_FL_combined = venn.diagram(vennlist_down_wilcox_C_R60_flag_leaf_FL_combined, category.names = c(paste(names(vennlist_down_wilcox_C_R60_flag_leaf_FL_combined[1]), " (", length(vennlist_down_wilcox_C_R60_flag_leaf_FL_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_C_R60_flag_leaf_FL_combined[2]), " (", length(vennlist_down_wilcox_C_R60_flag_leaf_FL_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_C_R60_flag_leaf_FL_combined[3]), " (", length(vennlist_down_wilcox_C_R60_flag_leaf_FL_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "F", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label F
intersect_wilcox_C_R60_venn_down_flag_leaf_FL_combined = attr(venn_down_wilcox_C_R60_flag_leaf_FL_combined, "intersections")  #list of unique and common metabolites
venn_down_wilcox_metabolites_C_R60_flag_leaf_FL_combined = t(ldply(intersect_wilcox_C_R60_venn_down_flag_leaf_FL_combined, rbind))  #convert from list to matrix
colnames(venn_down_wilcox_metabolites_C_R60_flag_leaf_FL_combined) = venn_down_wilcox_metabolites_C_R60_flag_leaf_FL_combined[1, ]
venn_down_wilcox_metabolites_C_R60_flag_leaf_FL_combined = venn_down_wilcox_metabolites_C_R60_flag_leaf_FL_combined[-1, ]


write.table(venn_up_wilcox_metabolites_C_R60_flag_leaf_FL_combined, file = "wilcox_sig_up_metabolites_flag_leaf_FL_R60_control_combined.txt", sep = "\t", quote = F)
write.table(venn_down_wilcox_metabolites_C_R60_flag_leaf_FL_combined, file = "wilcox_sig_down_metabolites_flag_leaf_FL_R60_control_combined.txt", sep = "\t", quote = F)


####################


#rownames were checked to be the same for the data files to be combined
#significance of wilcoxon test - all cultivars, all three RW timepoints and LS/C (refer to previous paper and associated GitHub files)
wilcox_sig_C_LS_R12_R36_R60_flag_leaf_FL_combined = cbind(wilcox_res_N22_C_LS_flag_leaf_FL_combined,
                                                          wilcox_res_N22_C_R12_flag_leaf_FL_combined$Significance,
                                                          wilcox_res_N22_C_R36_flag_leaf_FL_combined$Significance,
                                                          wilcox_res_N22_C_R60_flag_leaf_FL_combined$Significance,
                                                          wilcox_res_Dular_C_LS_flag_leaf_FL_combined$Significance,
                                                          wilcox_res_Dular_C_R12_flag_leaf_FL_combined$Significance,
                                                          wilcox_res_Dular_C_R36_flag_leaf_FL_combined$Significance,
                                                          wilcox_res_Dular_C_R60_flag_leaf_FL_combined$Significance,
                                                          wilcox_res_Anjali_C_LS_flag_leaf_FL_combined$Significance,
                                                          wilcox_res_Anjali_C_R12_flag_leaf_FL_combined$Significance,
                                                          wilcox_res_Anjali_C_R36_flag_leaf_FL_combined$Significance,
                                                          wilcox_res_Anjali_C_R60_flag_leaf_FL_combined$Significance)
wilcox_sig_C_LS_R12_R36_R60_flag_leaf_FL_combined = wilcox_sig_C_LS_R12_R36_R60_flag_leaf_FL_combined[ , -1]
colnames(wilcox_sig_C_LS_R12_R36_R60_flag_leaf_FL_combined)[1] = "wilcox_res_N22_C_LS_flag_leaf_FL_combined$Significance"
wilcox_sig_C_LS_R12_R36_R60_flag_leaf_FL_combined = data.frame(lapply(wilcox_sig_C_LS_R12_R36_R60_flag_leaf_FL_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_C_LS_R12_R36_R60_flag_leaf_FL_combined) = rownames(wilcox_res_N22_C_LS_flag_leaf_FL_combined)
wilcox_sig_C_LS_R12_R36_R60_flag_leaf_FL_combined$Nonsig = rowSums(wilcox_sig_C_LS_R12_R36_R60_flag_leaf_FL_combined == " ")  #count number of nonsignificant
wilcox_sig_final_C_LS_R12_R36_R60_flag_leaf_FL_combined = wilcox_sig_C_LS_R12_R36_R60_flag_leaf_FL_combined[-which(wilcox_sig_C_LS_R12_R36_R60_flag_leaf_FL_combined$Nonsig == "12"), ] #metabolites significant in at least one of the comparisons - 70 metabolites


#log2fc values of late stress and recovery relative to control (with reference to stress paper and associated GitHub files)
log2fc_combined_control_latestress_recovery_FL_flag_leaf = cbind(log2fc_combined_control_latestress_FL_flag_leaf, 
                                                                 log2fc_combined_control_recovery_FL_flag_leaf)
log2fc_combined_control_latestress_recovery_FL_flag_leaf = log2fc_combined_control_latestress_recovery_FL_flag_leaf[ , c(1, 4:6, 2, 7:9, 3, 10:12)]   #rearrange in preferred order


#log2fc values of metabolites with significant difference in at least one of the comparisons - late stress, 12, 36, 60h recovery vs. control in any cultivar
wilcox_sig_log2fc_C_LS_R12_R36_R60_flag_leaf_FL_combined = log2fc_combined_control_latestress_recovery_FL_flag_leaf[which(rownames(log2fc_combined_control_latestress_recovery_FL_flag_leaf) %in% rownames(wilcox_sig_final_C_LS_R12_R36_R60_flag_leaf_FL_combined)), ]


#rewatering-responsive metabolites only - indicate as red font in heat map
setdiff(rownames(wilcox_sig_log2fc_C_LS_R12_R36_R60_flag_leaf_FL_combined), rownames(wilcox_sig_log2fc_C_LS_flag_leaf_FL_combined))

#visualization
png("wilcox_sig_log2fc_C_LS_RW_flag_leaf_FL_combined.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_log2fc_C_LS_R12_R36_R60_flag_leaf_FL_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-3, -0.01, length.out = 300), seq(0.01, 3, length.out = 300)), 
          margins = c(7, 17),
          lmat = rbind(c(0, 4, 5), c(0, 1, 0), c(3, 2, 0)),
          labCol = rep(c("Severe stress/Control", "12 h RW/Control", "36 h RW/Control", "60 h RW/Control"), 3),
          cexCol = 1, 
          cexRow = 1, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 0.2, 11), 
          lwid = c(0.8, 5, 2.1),
          srtCol = 45, 
          cellnote = wilcox_sig_final_C_LS_R12_R36_R60_flag_leaf_FL_combined, 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 12)), 
          rowsep = c(seq(0, 70)), 
          sepwidth = c(0.01, 0.01), 
          ColSideColors = c(rep("#66c2a5", 4), rep("#fc8d62", 4), rep("#8da0cb", 4)),
          hclustfun = function(x) hclust(x, method = "average"),
          colRow = c(rep("black", 7), "red", rep("black", 2), "red", rep("black", 9), "red", rep("black", 2), rep("red", 2), 
                     rep("black", 2), "red", rep("black", 2), rep("red", 2), rep("black", 4), "red", rep("black", 2), "red", 
                     rep("black", 5), "red", rep("black", 2), "red", rep("black", 3), "red", rep("black", 5), "red", rep("black", 3),
                     "red", rep("black", 7)))
legend(x = -0.01, y = 1.1, xpd = TRUE, legend = "N22", bty = "n", cex = 1.2)
legend(x = 0.12, y = 1.1, xpd = TRUE, legend = "Dular", bty = "n", cex = 1.2)
legend(x = 0.26, y = 1.1, xpd = TRUE, legend = "Anjali", bty = "n", cex = 1.2)
legend(x = 0.83, y = 1.04, xpd = TRUE, legend = expression(paste("Log" ["2"], "-fold change")), bty = "n", cex = 0.8)
dev.off()


####################

#venn diagrams

#multiplot
png("venn_up_down_flag_leaf_FL_RW_control_combined.png", width = 16*300, height = 24*300, res = 300)
grid.arrange(gTree(children = venn_up_R12_control_flag_leaf_FL_combined), 
             gTree(children = venn_down_R12_control_flag_leaf_FL_combined), 
             gTree(children = venn_up_R36_control_flag_leaf_FL_combined),
             gTree(children = venn_down_R36_control_flag_leaf_FL_combined),
             gTree(children = venn_up_R60_control_flag_leaf_FL_combined),
             gTree(children = venn_down_R60_control_flag_leaf_FL_combined),
             ncol = 2, nrow = 3)
dev.off()


####################


#Early grain filling stage

#late stress vs. recovery

#log2-fold change - late stress vs. recovery
latestress_recovery_order_EGF = c("GNLS", "GNR12", "GNR36", "GNR60", "GDLS", "GDR12", "GDR36", "GDR60", "GALS", "GAR12", "GAR36", "GAR60")
data_combined_latestress_recovery_log2_mean_EGF_flag_leaf = data_combined_log2_mean_EGF_flag_leaf_2[latestress_recovery_order_EGF, ]  #late stress and recovery data in preferred order
data_combined_latestress_recovery_log2_mean_EGF_flag_leaf = t(data_combined_latestress_recovery_log2_mean_EGF_flag_leaf)

#log2-fc per cultivar
log2fc_combined_N22_latestress_recovery_EGF_flag_leaf = data_combined_latestress_recovery_log2_mean_EGF_flag_leaf[ , 2:4] - data_combined_latestress_recovery_log2_mean_EGF_flag_leaf[ , 1]
log2fc_combined_Dular_latestress_recovery_EGF_flag_leaf = data_combined_latestress_recovery_log2_mean_EGF_flag_leaf[ , 6:8] - data_combined_latestress_recovery_log2_mean_EGF_flag_leaf[ , 5]
log2fc_combined_Anjali_latestress_recovery_EGF_flag_leaf = data_combined_latestress_recovery_log2_mean_EGF_flag_leaf[ , 10:12] - data_combined_latestress_recovery_log2_mean_EGF_flag_leaf[ , 9]

#combine log2-fc of all cultivars in one data set
log2fc_combined_latestress_recovery_EGF_flag_leaf = cbind(log2fc_combined_N22_latestress_recovery_EGF_flag_leaf, log2fc_combined_Dular_latestress_recovery_EGF_flag_leaf, log2fc_combined_Anjali_latestress_recovery_EGF_flag_leaf)
colnames(log2fc_combined_latestress_recovery_EGF_flag_leaf)[1:3] = c(paste(colnames(log2fc_combined_latestress_recovery_EGF_flag_leaf)[1:3], "- GNLS"))
colnames(log2fc_combined_latestress_recovery_EGF_flag_leaf)[4:6] = c(paste(colnames(log2fc_combined_latestress_recovery_EGF_flag_leaf)[4:6], "- GDLS"))
colnames(log2fc_combined_latestress_recovery_EGF_flag_leaf)[7:9] = c(paste(colnames(log2fc_combined_latestress_recovery_EGF_flag_leaf)[7:9], "- GALS"))

#heatmap including all metabolites
heatmap.2(log2fc_combined_latestress_recovery_EGF_flag_leaf, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none",
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(7, 20), 
          main = "Flag leaf_EGF - Combined - Log2 FC, Recovery/Late Stress", 
          cexCol = 1, 
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          srtCol = 90)


#consider only the stress-responsive metabolites (i.e. metabolites which have significant difference in LS/CLS comparison)

#late stress vs. 12h rewatering


#N22 - 51 metabolites which have significant difference in LS/CLS comparison
data_combined_N22_log2median_stress_responsive_EGF_flag_leaf = data_combined_N22_log2median_EGF_flag_leaf
colnames(data_combined_N22_log2median_stress_responsive_EGF_flag_leaf)[3:ncol(data_combined_N22_log2median_stress_responsive_EGF_flag_leaf)] = reduced_metabolite_list_final_flag_leaf_2013$Name
#only metabolites which have significant difference in LS/C comparison
data_combined_N22_log2median_stress_responsive_EGF_flag_leaf = data_combined_N22_log2median_stress_responsive_EGF_flag_leaf[ , c(1, 2, (which(colnames(data_combined_N22_log2median_stress_responsive_EGF_flag_leaf) %in% rownames(wilcox_sig_log2fc_CLS_LS_flag_leaf_EGF_combined))))]

#data distribution
shapiro_res_combined_N22_latestress_recovery12_EGF_flag_leaf = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_EGF_flag_leaf[which(data_combined_N22_log2median_stress_responsive_EGF_flag_leaf$Timepoint %in% c("Late stress", "12h Rewatering")), 3:ncol(data_combined_N22_log2median_stress_responsive_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_latestress_recovery12_EGF_flag_leaf) = "p.value"
shapiro_res_combined_N22_latestress_recovery12_EGF_flag_leaf = as.data.frame(shapiro_res_combined_N22_latestress_recovery12_EGF_flag_leaf)
shapiro_res_combined_N22_latestress_recovery12_EGF_flag_leaf$Distribution = ifelse(shapiro_res_combined_N22_latestress_recovery12_EGF_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_latestress_recovery12_EGF_flag_leaf$Distribution == "Normal")     #36 metabolites
sum(shapiro_res_combined_N22_latestress_recovery12_EGF_flag_leaf$Distribution == "Non-normal") #15 metabolites

#wilcox test
wilcox_res_N22_LS_R12_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_EGF_flag_leaf[ , 3:ncol(data_combined_N22_log2median_stress_responsive_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_N22_log2median_stress_responsive_EGF_flag_leaf$Timepoint, subset = data_combined_N22_log2median_stress_responsive_EGF_flag_leaf$Timepoint %in% c("Late stress", "12h Rewatering"))$p.value)))
colnames(wilcox_res_N22_LS_R12_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_N22_LS_R12_flag_leaf_EGF_combined) = rownames(wilcox_sig_log2fc_CLS_LS_flag_leaf_EGF_combined)
wilcox_res_N22_LS_R12_flag_leaf_EGF_combined = as.data.frame(wilcox_res_N22_LS_R12_flag_leaf_EGF_combined)
wilcox_res_N22_LS_R12_flag_leaf_EGF_combined$Significance = with(wilcox_res_N22_LS_R12_flag_leaf_EGF_combined,
                                                                 ifelse(p.value <= 0.001, "***", ifelse(
                                                                   p.value <= 0.01, "**", ifelse(
                                                                     p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_LS_R12_flag_leaf_EGF_combined$Significance == "*") #3 metabolites
sum(wilcox_res_N22_LS_R12_flag_leaf_EGF_combined$Significance == "**") #2 metabolites
sum(wilcox_res_N22_LS_R12_flag_leaf_EGF_combined$Significance == "***") #0 metabolites
sum(wilcox_res_N22_LS_R12_flag_leaf_EGF_combined$Significance != "ns") #5 metabolites
sum(wilcox_res_N22_LS_R12_flag_leaf_EGF_combined$Significance == "ns") #46 metabolites 


#Dular - 51 metabolites which have significant difference in LS/CLS comparison
data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf = data_combined_Dular_log2median_EGF_flag_leaf
colnames(data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf)[3:ncol(data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf)] = reduced_metabolite_list_final_flag_leaf_2013$Name
#only metabolites which have significant difference in LS/C comparison
data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf = data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf[ , c(1, 2, (which(colnames(data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf) %in% rownames(wilcox_sig_log2fc_CLS_LS_flag_leaf_EGF_combined))))]

#data distribution
shapiro_res_combined_Dular_latestress_recovery12_EGF_flag_leaf = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf[which(data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf$Timepoint %in% c("Late stress", "12h Rewatering")), 3:ncol(data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_latestress_recovery12_EGF_flag_leaf) = "p.value"
shapiro_res_combined_Dular_latestress_recovery12_EGF_flag_leaf = as.data.frame(shapiro_res_combined_Dular_latestress_recovery12_EGF_flag_leaf)
shapiro_res_combined_Dular_latestress_recovery12_EGF_flag_leaf$Distribution = ifelse(shapiro_res_combined_Dular_latestress_recovery12_EGF_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_latestress_recovery12_EGF_flag_leaf$Distribution == "Normal")     #33 metabolites
sum(shapiro_res_combined_Dular_latestress_recovery12_EGF_flag_leaf$Distribution == "Non-normal") #18 metabolites

#wilcox test
wilcox_res_Dular_LS_R12_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf[ , 3:ncol(data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf$Timepoint, subset = data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf$Timepoint %in% c("Late stress", "12h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_LS_R12_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_Dular_LS_R12_flag_leaf_EGF_combined) = rownames(wilcox_sig_log2fc_CLS_LS_flag_leaf_EGF_combined)
wilcox_res_Dular_LS_R12_flag_leaf_EGF_combined = as.data.frame(wilcox_res_Dular_LS_R12_flag_leaf_EGF_combined)
wilcox_res_Dular_LS_R12_flag_leaf_EGF_combined$Significance = with(wilcox_res_Dular_LS_R12_flag_leaf_EGF_combined,
                                                                   ifelse(p.value <= 0.001, "***", ifelse(
                                                                     p.value <= 0.01, "**", ifelse(
                                                                       p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_LS_R12_flag_leaf_EGF_combined$Significance == "*") #4 metabolites
sum(wilcox_res_Dular_LS_R12_flag_leaf_EGF_combined$Significance == "**") #1 metabolites
sum(wilcox_res_Dular_LS_R12_flag_leaf_EGF_combined$Significance == "***") #0 metabolite
sum(wilcox_res_Dular_LS_R12_flag_leaf_EGF_combined$Significance != "ns") #5 metabolites
sum(wilcox_res_Dular_LS_R12_flag_leaf_EGF_combined$Significance == "ns") #46 metabolites 


#Anjali
#Anjali - 51 metabolites which have significant difference in LS/CLS comparison
data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf = data_combined_Anjali_log2median_EGF_flag_leaf
colnames(data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf)[3:ncol(data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf)] = reduced_metabolite_list_final_flag_leaf_2013$Name
#only metabolites which have significant difference in LS/C comparison
data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf = data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf[ , c(1, 2, (which(colnames(data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf) %in% rownames(wilcox_sig_log2fc_CLS_LS_flag_leaf_EGF_combined))))]

#data distribution
shapiro_res_combined_Anjali_latestress_recovery12_EGF_flag_leaf = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf[which(data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf$Timepoint %in% c("Late stress", "12h Rewatering")), 3:ncol(data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_latestress_recovery12_EGF_flag_leaf) = "p.value"
shapiro_res_combined_Anjali_latestress_recovery12_EGF_flag_leaf = as.data.frame(shapiro_res_combined_Anjali_latestress_recovery12_EGF_flag_leaf)
shapiro_res_combined_Anjali_latestress_recovery12_EGF_flag_leaf$Distribution = ifelse(shapiro_res_combined_Anjali_latestress_recovery12_EGF_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_latestress_recovery12_EGF_flag_leaf$Distribution == "Normal")     #32 metabolites
sum(shapiro_res_combined_Anjali_latestress_recovery12_EGF_flag_leaf$Distribution == "Non-normal") #19 metabolites

#wilcox test
wilcox_res_Anjali_LS_R12_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf[ , 3:ncol(data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf$Timepoint, subset = data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf$Timepoint %in% c("Late stress", "12h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_LS_R12_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_Anjali_LS_R12_flag_leaf_EGF_combined) = rownames(wilcox_sig_log2fc_CLS_LS_flag_leaf_EGF_combined)
wilcox_res_Anjali_LS_R12_flag_leaf_EGF_combined = as.data.frame(wilcox_res_Anjali_LS_R12_flag_leaf_EGF_combined)
wilcox_res_Anjali_LS_R12_flag_leaf_EGF_combined$Significance = with(wilcox_res_Anjali_LS_R12_flag_leaf_EGF_combined,
                                                                    ifelse(p.value <= 0.001, "***", ifelse(
                                                                      p.value <= 0.01, "**", ifelse(
                                                                        p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_LS_R12_flag_leaf_EGF_combined$Significance == "*") #7 metabolites
sum(wilcox_res_Anjali_LS_R12_flag_leaf_EGF_combined$Significance == "**") #2 metabolites
sum(wilcox_res_Anjali_LS_R12_flag_leaf_EGF_combined$Significance == "***") #1 metabolite
sum(wilcox_res_Anjali_LS_R12_flag_leaf_EGF_combined$Significance != "ns") #10 metabolites
sum(wilcox_res_Anjali_LS_R12_flag_leaf_EGF_combined$Significance == "ns") #41 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_LS_R12_flag_leaf_EGF_combined = cbind(wilcox_res_N22_LS_R12_flag_leaf_EGF_combined, 
                                                 wilcox_res_Dular_LS_R12_flag_leaf_EGF_combined$Significance, 
                                                 wilcox_res_Anjali_LS_R12_flag_leaf_EGF_combined$Significance)
wilcox_sig_LS_R12_flag_leaf_EGF_combined = wilcox_sig_LS_R12_flag_leaf_EGF_combined[ , -1]
colnames(wilcox_sig_LS_R12_flag_leaf_EGF_combined)[1] = "wilcox_res_N22_LS_R12_flag_leaf_EGF_combined$Significance"
wilcox_sig_LS_R12_flag_leaf_EGF_combined = data.frame(lapply(wilcox_sig_LS_R12_flag_leaf_EGF_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_LS_R12_flag_leaf_EGF_combined) = rownames(wilcox_res_N22_LS_R12_flag_leaf_EGF_combined)
wilcox_sig_LS_R12_flag_leaf_EGF_combined$Nonsig = rowSums(wilcox_sig_LS_R12_flag_leaf_EGF_combined == " ")  #count number of nonsignificant
wilcox_sig_final_LS_R12_flag_leaf_EGF_combined = wilcox_sig_LS_R12_flag_leaf_EGF_combined[-which(wilcox_sig_LS_R12_flag_leaf_EGF_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 17 metabolites

#log2fc of metabolites significant in at least one of the cultivars
wilcox_sig_log2fc_LS_R12_flag_leaf_EGF_combined = log2fc_combined_latestress_recovery_EGF_flag_leaf[which(rownames(log2fc_combined_latestress_recovery_EGF_flag_leaf) %in% rownames(wilcox_sig_final_LS_R12_flag_leaf_EGF_combined)), c(1, 4, 7)]  

#overview - with HCA - Euclidean distance and average linkage 
png("wilcox_sig_stress_responsive_log2fc_LS_R12_flag_leaf_EGF_combined.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_log2fc_LS_R12_flag_leaf_EGF_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(5, 20),  
          main = "Flag leaf-Combined-Wilcoxon-Stress-responsive \nLog2 FC, EGF, 12h RW/Late stress", 
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
          cellnote = wilcox_sig_final_LS_R12_flag_leaf_EGF_combined, 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 3)), 
          rowsep = c(seq(0, 22)), 
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))
dev.off()


#venn diagram - common among cultivars

#N22
wilcox_sig_stress_responsive_N22_LS_R12_flag_leaf_EGF_combined = wilcox_res_N22_LS_R12_flag_leaf_EGF_combined[-which(wilcox_res_N22_LS_R12_flag_leaf_EGF_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_LS_R12_flag_leaf_EGF_combined[match(rownames(wilcox_sig_stress_responsive_N22_LS_R12_flag_leaf_EGF_combined), rownames(wilcox_sig_log2fc_LS_R12_flag_leaf_EGF_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_LS_R12_flag_leaf_EGF_combined)[1]
wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined$`GNR12 - GNLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined$Up.Down == "Up")    #0 metabolite
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined$Up.Down == "Down")  #5 metabolites

#Dular
wilcox_sig_stress_responsive_Dular_LS_R12_flag_leaf_EGF_combined = wilcox_res_Dular_LS_R12_flag_leaf_EGF_combined[-which(wilcox_res_Dular_LS_R12_flag_leaf_EGF_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_LS_R12_flag_leaf_EGF_combined[match(rownames(wilcox_sig_stress_responsive_Dular_LS_R12_flag_leaf_EGF_combined), rownames(wilcox_sig_log2fc_LS_R12_flag_leaf_EGF_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_LS_R12_flag_leaf_EGF_combined)[2]
wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined$`GDR12 - GDLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined$Up.Down == "Up")    #3 metabolites
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined$Up.Down == "Down")  #2 metabolites

#Anjali
wilcox_sig_stress_responsive_Anjali_LS_R12_flag_leaf_EGF_combined = wilcox_res_Anjali_LS_R12_flag_leaf_EGF_combined[-which(wilcox_res_Anjali_LS_R12_flag_leaf_EGF_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_LS_R12_flag_leaf_EGF_combined[match(rownames(wilcox_sig_stress_responsive_Anjali_LS_R12_flag_leaf_EGF_combined), rownames(wilcox_sig_log2fc_LS_R12_flag_leaf_EGF_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_LS_R12_flag_leaf_EGF_combined)[3]
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined$`GAR12 - GALS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined$Up.Down == "Up")    #0 metabolite
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined$Up.Down == "Down")  #10 metabolites


#venn - increased - 12h rewatering/late stress
vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined) = c("N22", "Dular", "Anjali")
venn_up_stress_responsive_wilcox_LS_R12_flag_leaf_EGF_combined = venn(vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_stress_responsive_R12_latestress_flag_leaf_EGF_combined = venn.diagram(vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined, category.names = c(paste(names(vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined[1]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined[2]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined[3]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "A", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label A
intersect_stress_responsive_wilcox_LS_R12_venn_up_flag_leaf_EGF_combined = attr(venn_up_stress_responsive_wilcox_LS_R12_flag_leaf_EGF_combined, "intersections")  #list of unique and common metabolites
venn_up_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_EGF_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R12_venn_up_flag_leaf_EGF_combined, rbind))  #convert from list to matrix
colnames(venn_up_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_EGF_combined) = venn_up_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_EGF_combined[1, ]
venn_up_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_EGF_combined = venn_up_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_EGF_combined[-1, ]


#venn - decreased - 12h rewatering/late stress
vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R12_flag_leaf_FL_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R12_flag_leaf_FL_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R12_flag_leaf_FL_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined) = c("N22", "Dular", "Anjali")
venn_down_stress_responsive_wilcox_LS_R12_flag_leaf_EGF_combined = venn(vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_stress_responsive_R12_latestress_flag_leaf_EGF_combined = venn.diagram(vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined, category.names = c(paste(names(vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined[1]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined[2]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined[3]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R12_flag_leaf_EGF_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "B", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label B
intersect_stress_responsive_wilcox_LS_R12_venn_down_flag_leaf_EGF_combined = attr(venn_down_stress_responsive_wilcox_LS_R12_flag_leaf_EGF_combined, "intersections")  #list of unique and common metabolites
venn_down_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_EGF_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R12_venn_down_flag_leaf_EGF_combined, rbind))  #convert from list to matrix
colnames(venn_down_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_EGF_combined) = venn_down_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_EGF_combined[1, ]
venn_down_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_EGF_combined = venn_down_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_EGF_combined[-1, ]


write.table(venn_up_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_EGF_combined, file = "wilcox_stress_responsive_sig_up_metabolites_flag_leaf_EGF_R12_latestress_combined.txt", sep = "\t", quote = F)
write.table(venn_down_stress_responsive_wilcox_metabolites_LS_R12_flag_leaf_EGF_combined, file = "wilcox_stress_responsive_sig_down_metabolites_flag_leaf_EGF_R12_latestress_combined.txt", sep = "\t", quote = F)


####################

#late stress vs. 36h rewatering


#N22

#data distribution
shapiro_res_combined_N22_latestress_recovery36_EGF_flag_leaf = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_EGF_flag_leaf[which(data_combined_N22_log2median_stress_responsive_EGF_flag_leaf$Timepoint %in% c("Late stress", "36h Rewatering")), 3:ncol(data_combined_N22_log2median_stress_responsive_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_latestress_recovery36_EGF_flag_leaf) = "p.value"
shapiro_res_combined_N22_latestress_recovery36_EGF_flag_leaf = as.data.frame(shapiro_res_combined_N22_latestress_recovery36_EGF_flag_leaf)
shapiro_res_combined_N22_latestress_recovery36_EGF_flag_leaf$Distribution = ifelse(shapiro_res_combined_N22_latestress_recovery36_EGF_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_latestress_recovery36_EGF_flag_leaf$Distribution == "Normal")     #39 metabolites
sum(shapiro_res_combined_N22_latestress_recovery36_EGF_flag_leaf$Distribution == "Non-normal") #12 metabolites

#wilcox test
wilcox_res_N22_LS_R36_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_EGF_flag_leaf[ , 3:ncol(data_combined_N22_log2median_stress_responsive_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_N22_log2median_stress_responsive_EGF_flag_leaf$Timepoint, subset = data_combined_N22_log2median_stress_responsive_EGF_flag_leaf$Timepoint %in% c("Late stress", "36h Rewatering"))$p.value)))
colnames(wilcox_res_N22_LS_R36_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_N22_LS_R36_flag_leaf_EGF_combined) = rownames(wilcox_sig_log2fc_CLS_LS_flag_leaf_EGF_combined)
wilcox_res_N22_LS_R36_flag_leaf_EGF_combined = as.data.frame(wilcox_res_N22_LS_R36_flag_leaf_EGF_combined)
wilcox_res_N22_LS_R36_flag_leaf_EGF_combined$Significance = with(wilcox_res_N22_LS_R36_flag_leaf_EGF_combined,
                                                                 ifelse(p.value <= 0.001, "***", ifelse(
                                                                   p.value <= 0.01, "**", ifelse(
                                                                     p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_LS_R36_flag_leaf_EGF_combined$Significance == "*") #9 metabolites
sum(wilcox_res_N22_LS_R36_flag_leaf_EGF_combined$Significance == "**") #5 metabolites
sum(wilcox_res_N22_LS_R36_flag_leaf_EGF_combined$Significance == "***") #9 metabolites
sum(wilcox_res_N22_LS_R36_flag_leaf_EGF_combined$Significance != "ns") #23 metabolites
sum(wilcox_res_N22_LS_R36_flag_leaf_EGF_combined$Significance == "ns") #28 metabolites 


#Dular

#data distribution
shapiro_res_combined_Dular_latestress_recovery36_EGF_flag_leaf = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf[which(data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf$Timepoint %in% c("Late stress", "36h Rewatering")), 3:ncol(data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_latestress_recovery36_EGF_flag_leaf) = "p.value"
shapiro_res_combined_Dular_latestress_recovery36_EGF_flag_leaf = as.data.frame(shapiro_res_combined_Dular_latestress_recovery36_EGF_flag_leaf)
shapiro_res_combined_Dular_latestress_recovery36_EGF_flag_leaf$Distribution = ifelse(shapiro_res_combined_Dular_latestress_recovery36_EGF_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_latestress_recovery36_EGF_flag_leaf$Distribution == "Normal")     #31 metabolites
sum(shapiro_res_combined_Dular_latestress_recovery36_EGF_flag_leaf$Distribution == "Non-normal") #20 metabolites

#wilcox test
wilcox_res_Dular_LS_R36_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf[ , 3:ncol(data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf$Timepoint, subset = data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf$Timepoint %in% c("Late stress", "36h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_LS_R36_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_Dular_LS_R36_flag_leaf_EGF_combined) = rownames(wilcox_sig_log2fc_CLS_LS_flag_leaf_EGF_combined)
wilcox_res_Dular_LS_R36_flag_leaf_EGF_combined = as.data.frame(wilcox_res_Dular_LS_R36_flag_leaf_EGF_combined)
wilcox_res_Dular_LS_R36_flag_leaf_EGF_combined$Significance = with(wilcox_res_Dular_LS_R36_flag_leaf_EGF_combined,
                                                                   ifelse(p.value <= 0.001, "***", ifelse(
                                                                     p.value <= 0.01, "**", ifelse(
                                                                       p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_LS_R36_flag_leaf_EGF_combined$Significance == "*") #7 metabolites
sum(wilcox_res_Dular_LS_R36_flag_leaf_EGF_combined$Significance == "**") #4 metabolites
sum(wilcox_res_Dular_LS_R36_flag_leaf_EGF_combined$Significance == "***") #7 metabolite
sum(wilcox_res_Dular_LS_R36_flag_leaf_EGF_combined$Significance != "ns") #18 metabolites
sum(wilcox_res_Dular_LS_R36_flag_leaf_EGF_combined$Significance == "ns") #33 metabolites 


#Anjali

#data distribution
shapiro_res_combined_Anjali_latestress_recovery36_EGF_flag_leaf = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf[which(data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf$Timepoint %in% c("Late stress", "36h Rewatering")), 3:ncol(data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_latestress_recovery36_EGF_flag_leaf) = "p.value"
shapiro_res_combined_Anjali_latestress_recovery36_EGF_flag_leaf = as.data.frame(shapiro_res_combined_Anjali_latestress_recovery36_EGF_flag_leaf)
shapiro_res_combined_Anjali_latestress_recovery36_EGF_flag_leaf$Distribution = ifelse(shapiro_res_combined_Anjali_latestress_recovery36_EGF_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_latestress_recovery36_EGF_flag_leaf$Distribution == "Normal")     #41 metabolites
sum(shapiro_res_combined_Anjali_latestress_recovery36_EGF_flag_leaf$Distribution == "Non-normal") #10 metabolites

#wilcox test
wilcox_res_Anjali_LS_R36_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf[ , 3:ncol(data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf$Timepoint, subset = data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf$Timepoint %in% c("Late stress", "36h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_LS_R36_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_Anjali_LS_R36_flag_leaf_EGF_combined) = rownames(wilcox_sig_log2fc_CLS_LS_flag_leaf_EGF_combined)
wilcox_res_Anjali_LS_R36_flag_leaf_EGF_combined = as.data.frame(wilcox_res_Anjali_LS_R36_flag_leaf_EGF_combined)
wilcox_res_Anjali_LS_R36_flag_leaf_EGF_combined$Significance = with(wilcox_res_Anjali_LS_R36_flag_leaf_EGF_combined,
                                                                    ifelse(p.value <= 0.001, "***", ifelse(
                                                                      p.value <= 0.01, "**", ifelse(
                                                                        p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_LS_R36_flag_leaf_EGF_combined$Significance == "*") #11 metabolites
sum(wilcox_res_Anjali_LS_R36_flag_leaf_EGF_combined$Significance == "**") #5 metabolites
sum(wilcox_res_Anjali_LS_R36_flag_leaf_EGF_combined$Significance == "***") #4 metabolite
sum(wilcox_res_Anjali_LS_R36_flag_leaf_EGF_combined$Significance != "ns") #20 metabolites
sum(wilcox_res_Anjali_LS_R36_flag_leaf_EGF_combined$Significance == "ns") #31 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_LS_R36_flag_leaf_EGF_combined = cbind(wilcox_res_N22_LS_R36_flag_leaf_EGF_combined, 
                                                 wilcox_res_Dular_LS_R36_flag_leaf_EGF_combined$Significance,
                                                 wilcox_res_Anjali_LS_R36_flag_leaf_EGF_combined$Significance)
wilcox_sig_LS_R36_flag_leaf_EGF_combined = wilcox_sig_LS_R36_flag_leaf_EGF_combined[ , -1]
colnames(wilcox_sig_LS_R36_flag_leaf_EGF_combined)[1] = "wilcox_res_N22_LS_R36_flag_leaf_EGF_combined$Significance"
wilcox_sig_LS_R36_flag_leaf_EGF_combined = data.frame(lapply(wilcox_sig_LS_R36_flag_leaf_EGF_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_LS_R36_flag_leaf_EGF_combined) = rownames(wilcox_res_N22_LS_R36_flag_leaf_EGF_combined)
wilcox_sig_LS_R36_flag_leaf_EGF_combined$Nonsig = rowSums(wilcox_sig_LS_R36_flag_leaf_EGF_combined == " ")  #count number of nonsignificant
wilcox_sig_final_LS_R36_flag_leaf_EGF_combined = wilcox_sig_LS_R36_flag_leaf_EGF_combined[-which(wilcox_sig_LS_R36_flag_leaf_EGF_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 32 metabolites

#log2fc of metabolites significant in at least one of the cultivars
wilcox_sig_log2fc_LS_R36_flag_leaf_EGF_combined = log2fc_combined_latestress_recovery_EGF_flag_leaf[which(rownames(log2fc_combined_latestress_recovery_EGF_flag_leaf) %in% rownames(wilcox_sig_final_LS_R36_flag_leaf_EGF_combined)), c(2, 5, 8)]  

#overview - with HCA - Euclidean distance and average linkage
png("wilcox_sig_stress_responsive_log2fc_LS_R36_flag_leaf_EGF_combined.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_log2fc_LS_R36_flag_leaf_EGF_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(5, 20), 
          main = "Flag leaf-Combined-Wilcoxon-Stress-responsive \nLog2 FC, EGF, 36h RW/Late stress", 
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
          cellnote = wilcox_sig_final_LS_R36_flag_leaf_EGF_combined, 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black",
          colsep = c(seq(0, 3)), 
          rowsep = c(seq(0, 32)), 
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))
dev.off()


#venn diagram - common among cultivars

#N22
wilcox_sig_stress_responsive_N22_LS_R36_flag_leaf_EGF_combined = wilcox_res_N22_LS_R36_flag_leaf_EGF_combined[-which(wilcox_res_N22_LS_R36_flag_leaf_EGF_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_LS_R36_flag_leaf_EGF_combined[match(rownames(wilcox_sig_stress_responsive_N22_LS_R36_flag_leaf_EGF_combined), rownames(wilcox_sig_log2fc_LS_R36_flag_leaf_EGF_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_LS_R36_flag_leaf_EGF_combined)[1]
wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined$`GNR36 - GNLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined$Up.Down == "Up")    #4 metabolites
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined$Up.Down == "Down")  #19 metabolites

#Dular
wilcox_sig_stress_responsive_Dular_LS_R36_flag_leaf_EGF_combined = wilcox_res_Dular_LS_R36_flag_leaf_EGF_combined[-which(wilcox_res_Dular_LS_R36_flag_leaf_EGF_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_LS_R36_flag_leaf_EGF_combined[match(rownames(wilcox_sig_stress_responsive_Dular_LS_R36_flag_leaf_EGF_combined), rownames(wilcox_sig_log2fc_LS_R36_flag_leaf_EGF_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_LS_R36_flag_leaf_EGF_combined)[2]
wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined$`GDR36 - GDLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined$Up.Down == "Up")    #2 metabolites
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined$Up.Down == "Down")  #16 metabolites

#Anjali
wilcox_sig_stress_responsive_Anjali_LS_R36_flag_leaf_EGF_combined = wilcox_res_Anjali_LS_R36_flag_leaf_EGF_combined[-which(wilcox_res_Anjali_LS_R36_flag_leaf_EGF_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_LS_R36_flag_leaf_EGF_combined[match(rownames(wilcox_sig_stress_responsive_Anjali_LS_R36_flag_leaf_EGF_combined), rownames(wilcox_sig_log2fc_LS_R36_flag_leaf_EGF_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_LS_R36_flag_leaf_EGF_combined)[3]
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined$`GAR36 - GALS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined$Up.Down == "Up")    #5 metabolites
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined$Up.Down == "Down")  #15 metabolites


#venn - increased - 36h rewatering/late stress
vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined) = c("N22", "Dular", "Anjali")
venn_up_stress_responsive_wilcox_LS_R36_flag_leaf_EGF_combined = venn(vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_stress_responsive_R36_latestress_flag_leaf_EGF_combined = venn.diagram(vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined, category.names = c(paste(names(vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined[1]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined[2]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined[3]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "C", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label C
intersect_stress_responsive_wilcox_LS_R36_venn_up_flag_leaf_EGF_combined = attr(venn_up_stress_responsive_wilcox_LS_R36_flag_leaf_EGF_combined, "intersections")  #list of unique and common metabolites
venn_up_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_EGF_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R36_venn_up_flag_leaf_EGF_combined, rbind))  #convert from list to matrix
colnames(venn_up_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_EGF_combined) = venn_up_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_EGF_combined[1, ]
venn_up_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_EGF_combined = venn_up_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_EGF_combined[-1, ]


#venn - decreased - 36h rewatering/late stress
vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R36_flag_leaf_FL_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R36_flag_leaf_FL_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R36_flag_leaf_FL_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined) = c("N22", "Dular", "Anjali")
venn_down_stress_responsive_wilcox_LS_R36_flag_leaf_EGF_combined = venn(vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_stress_responsive_R36_latestress_flag_leaf_EGF_combined = venn.diagram(vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined, category.names = c(paste(names(vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined[1]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined[2]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined[3]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R36_flag_leaf_EGF_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "D", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label D
intersect_stress_responsive_wilcox_LS_R36_venn_down_flag_leaf_EGF_combined = attr(venn_down_stress_responsive_wilcox_LS_R36_flag_leaf_EGF_combined, "intersections")  #list of unique and common metabolites
venn_down_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_EGF_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R36_venn_down_flag_leaf_EGF_combined, rbind))  #convert from list to matrix
colnames(venn_down_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_EGF_combined) = venn_down_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_EGF_combined[1, ]
venn_down_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_EGF_combined = venn_down_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_EGF_combined[-1, ]


write.table(venn_up_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_EGF_combined, file = "wilcox_stress_responsive_sig_up_metabolites_flag_leaf_EGF_R36_latestress_combined.txt", sep = "\t", quote = F)
write.table(venn_down_stress_responsive_wilcox_metabolites_LS_R36_flag_leaf_EGF_combined, file = "wilcox_stress_responsive_sig_down_metabolites_flag_leaf_EGF_R36_latestress_combined.txt", sep = "\t", quote = F)


####################


#late stress vs. 60h rewatering

#N22

#data distribution
shapiro_res_combined_N22_latestress_recovery60_EGF_flag_leaf = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_EGF_flag_leaf[which(data_combined_N22_log2median_stress_responsive_EGF_flag_leaf$Timepoint %in% c("Late stress", "60h Rewatering")), 3:ncol(data_combined_N22_log2median_stress_responsive_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_latestress_recovery60_EGF_flag_leaf) = "p.value"
shapiro_res_combined_N22_latestress_recovery60_EGF_flag_leaf = as.data.frame(shapiro_res_combined_N22_latestress_recovery60_EGF_flag_leaf)
shapiro_res_combined_N22_latestress_recovery60_EGF_flag_leaf$Distribution = ifelse(shapiro_res_combined_N22_latestress_recovery60_EGF_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_latestress_recovery60_EGF_flag_leaf$Distribution == "Normal")     #36 metabolites
sum(shapiro_res_combined_N22_latestress_recovery60_EGF_flag_leaf$Distribution == "Non-normal") #15 metabolites

#wilcox test
wilcox_res_N22_LS_R60_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_N22_log2median_stress_responsive_EGF_flag_leaf[ , 3:ncol(data_combined_N22_log2median_stress_responsive_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_N22_log2median_stress_responsive_EGF_flag_leaf$Timepoint, subset = data_combined_N22_log2median_stress_responsive_EGF_flag_leaf$Timepoint %in% c("Late stress", "60h Rewatering"))$p.value)))
colnames(wilcox_res_N22_LS_R60_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_N22_LS_R60_flag_leaf_EGF_combined) = rownames(wilcox_sig_log2fc_CLS_LS_flag_leaf_EGF_combined)
wilcox_res_N22_LS_R60_flag_leaf_EGF_combined = as.data.frame(wilcox_res_N22_LS_R60_flag_leaf_EGF_combined)
wilcox_res_N22_LS_R60_flag_leaf_EGF_combined$Significance = with(wilcox_res_N22_LS_R60_flag_leaf_EGF_combined,
                                                                 ifelse(p.value <= 0.001, "***", ifelse(
                                                                   p.value <= 0.01, "**", ifelse(
                                                                     p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_LS_R60_flag_leaf_EGF_combined$Significance == "*") #5 metabolites
sum(wilcox_res_N22_LS_R60_flag_leaf_EGF_combined$Significance == "**") #7 metabolites
sum(wilcox_res_N22_LS_R60_flag_leaf_EGF_combined$Significance == "***") #10 metabolites
sum(wilcox_res_N22_LS_R60_flag_leaf_EGF_combined$Significance != "ns") #22 metabolites
sum(wilcox_res_N22_LS_R60_flag_leaf_EGF_combined$Significance == "ns") #29 metabolites 


#Dular

#data distribution
shapiro_res_combined_Dular_latestress_recovery60_EGF_flag_leaf = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf[which(data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf$Timepoint %in% c("Late stress", "60h Rewatering")), 3:ncol(data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_latestress_recovery60_EGF_flag_leaf) = "p.value"
shapiro_res_combined_Dular_latestress_recovery60_EGF_flag_leaf = as.data.frame(shapiro_res_combined_Dular_latestress_recovery60_EGF_flag_leaf)
shapiro_res_combined_Dular_latestress_recovery60_EGF_flag_leaf$Distribution = ifelse(shapiro_res_combined_Dular_latestress_recovery60_EGF_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_latestress_recovery60_EGF_flag_leaf$Distribution == "Normal")     #38 metabolites
sum(shapiro_res_combined_Dular_latestress_recovery60_EGF_flag_leaf$Distribution == "Non-normal") #13 metabolites

#wilcox test
wilcox_res_Dular_LS_R60_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf[ , 3:ncol(data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf$Timepoint, subset = data_combined_Dular_log2median_stress_responsive_EGF_flag_leaf$Timepoint %in% c("Late stress", "60h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_LS_R60_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_Dular_LS_R60_flag_leaf_EGF_combined) = rownames(wilcox_sig_log2fc_CLS_LS_flag_leaf_EGF_combined)
wilcox_res_Dular_LS_R60_flag_leaf_EGF_combined = as.data.frame(wilcox_res_Dular_LS_R60_flag_leaf_EGF_combined)
wilcox_res_Dular_LS_R60_flag_leaf_EGF_combined$Significance = with(wilcox_res_Dular_LS_R60_flag_leaf_EGF_combined,
                                                                   ifelse(p.value <= 0.001, "***", ifelse(
                                                                     p.value <= 0.01, "**", ifelse(
                                                                       p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_LS_R60_flag_leaf_EGF_combined$Significance == "*") #11 metabolites
sum(wilcox_res_Dular_LS_R60_flag_leaf_EGF_combined$Significance == "**") #6 metabolites
sum(wilcox_res_Dular_LS_R60_flag_leaf_EGF_combined$Significance == "***") #4 metabolite
sum(wilcox_res_Dular_LS_R60_flag_leaf_EGF_combined$Significance != "ns") #21 metabolites
sum(wilcox_res_Dular_LS_R60_flag_leaf_EGF_combined$Significance == "ns") #30 metabolites


#Anjali

#data distribution
shapiro_res_combined_Anjali_latestress_recovery60_EGF_flag_leaf = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf[which(data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf$Timepoint %in% c("Late stress", "60h Rewatering")), 3:ncol(data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_latestress_recovery60_EGF_flag_leaf) = "p.value"
shapiro_res_combined_Anjali_latestress_recovery60_EGF_flag_leaf = as.data.frame(shapiro_res_combined_Anjali_latestress_recovery60_EGF_flag_leaf)
shapiro_res_combined_Anjali_latestress_recovery60_EGF_flag_leaf$Distribution = ifelse(shapiro_res_combined_Anjali_latestress_recovery60_EGF_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_latestress_recovery60_EGF_flag_leaf$Distribution == "Normal")     #38 metabolites
sum(shapiro_res_combined_Anjali_latestress_recovery60_EGF_flag_leaf$Distribution == "Non-normal") #13 metabolites

#wilcox test
wilcox_res_Anjali_LS_R60_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf[ , 3:ncol(data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf$Timepoint, subset = data_combined_Anjali_log2median_stress_responsive_EGF_flag_leaf$Timepoint %in% c("Late stress", "60h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_LS_R60_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_Anjali_LS_R60_flag_leaf_EGF_combined) = rownames(wilcox_sig_log2fc_CLS_LS_flag_leaf_EGF_combined)
wilcox_res_Anjali_LS_R60_flag_leaf_EGF_combined = as.data.frame(wilcox_res_Anjali_LS_R60_flag_leaf_EGF_combined)
wilcox_res_Anjali_LS_R60_flag_leaf_EGF_combined$Significance = with(wilcox_res_Anjali_LS_R60_flag_leaf_EGF_combined,
                                                                    ifelse(p.value <= 0.001, "***", ifelse(
                                                                      p.value <= 0.01, "**", ifelse(
                                                                        p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_LS_R60_flag_leaf_EGF_combined$Significance == "*") #6 metabolites
sum(wilcox_res_Anjali_LS_R60_flag_leaf_EGF_combined$Significance == "**") #2 metabolites
sum(wilcox_res_Anjali_LS_R60_flag_leaf_EGF_combined$Significance == "***") #2 metabolite
sum(wilcox_res_Anjali_LS_R60_flag_leaf_EGF_combined$Significance != "ns") #10 metabolites
sum(wilcox_res_Anjali_LS_R60_flag_leaf_EGF_combined$Significance == "ns") #41 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_LS_R60_flag_leaf_EGF_combined = cbind(wilcox_res_N22_LS_R60_flag_leaf_EGF_combined, 
                                                 wilcox_res_Dular_LS_R60_flag_leaf_EGF_combined$Significance, 
                                                 wilcox_res_Anjali_LS_R60_flag_leaf_EGF_combined$Significance)
wilcox_sig_LS_R60_flag_leaf_EGF_combined = wilcox_sig_LS_R60_flag_leaf_EGF_combined[ , -1]
colnames(wilcox_sig_LS_R60_flag_leaf_EGF_combined)[1] = "wilcox_res_N22_LS_R60_flag_leaf_EGF_combined$Significance"
wilcox_sig_LS_R60_flag_leaf_EGF_combined = data.frame(lapply(wilcox_sig_LS_R60_flag_leaf_EGF_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_LS_R60_flag_leaf_EGF_combined) = rownames(wilcox_res_N22_LS_R60_flag_leaf_EGF_combined)
wilcox_sig_LS_R60_flag_leaf_EGF_combined$Nonsig = rowSums(wilcox_sig_LS_R60_flag_leaf_EGF_combined == " ")  #count number of nonsignificant
wilcox_sig_final_LS_R60_flag_leaf_EGF_combined = wilcox_sig_LS_R60_flag_leaf_EGF_combined[-which(wilcox_sig_LS_R60_flag_leaf_EGF_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 31 metabolites

#log2fc of metabolites significant in at least one of the cultivars
wilcox_sig_log2fc_LS_R60_flag_leaf_EGF_combined = log2fc_combined_latestress_recovery_EGF_flag_leaf[which(rownames(log2fc_combined_latestress_recovery_EGF_flag_leaf) %in% rownames(wilcox_sig_final_LS_R60_flag_leaf_EGF_combined)), c(3, 6, 9)]  

#overview - with HCA - Euclidean distance and average linkage 
png("wilcox_sig_stress_responsive_log2fc_LS_R60_flag_leaf_EGF_combined.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_log2fc_LS_R60_flag_leaf_EGF_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none",
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(5, 20),  
          main = "Flag leaf-Combined-Wilcoxon-Stress-responsive \nLog2 FC, EGF, 60h RW/Late stress", 
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
          cellnote = wilcox_sig_final_LS_R60_flag_leaf_EGF_combined,
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 3)), 
          rowsep = c(seq(0, 31)), 
          sepwidth = c(0.01, 0.01),
          hclustfun = function(x) hclust(x, method = "average"))
dev.off()


#venn diagram - common among cultivars

#N22
wilcox_sig_stress_responsive_N22_LS_R60_flag_leaf_EGF_combined = wilcox_res_N22_LS_R60_flag_leaf_EGF_combined[-which(wilcox_res_N22_LS_R60_flag_leaf_EGF_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_LS_R60_flag_leaf_EGF_combined[match(rownames(wilcox_sig_stress_responsive_N22_LS_R60_flag_leaf_EGF_combined), rownames(wilcox_sig_log2fc_LS_R60_flag_leaf_EGF_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_LS_R60_flag_leaf_EGF_combined)[1]
wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined$`GNR60 - GNLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined$Up.Down == "Up")    #4 metabolites
sum(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined$Up.Down == "Down")  #18 metabolites

#Dular
wilcox_sig_stress_responsive_Dular_LS_R60_flag_leaf_EGF_combined = wilcox_res_Dular_LS_R60_flag_leaf_EGF_combined[-which(wilcox_res_Dular_LS_R60_flag_leaf_EGF_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_LS_R60_flag_leaf_EGF_combined[match(rownames(wilcox_sig_stress_responsive_Dular_LS_R60_flag_leaf_EGF_combined), rownames(wilcox_sig_log2fc_LS_R60_flag_leaf_EGF_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_LS_R60_flag_leaf_EGF_combined)[2]
wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined$`GDR60 - GDLS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined$Up.Down == "Up")    #4 metabolites
sum(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined$Up.Down == "Down")  #17 metabolites

#Anjali
wilcox_sig_stress_responsive_Anjali_LS_R60_flag_leaf_EGF_combined = wilcox_res_Anjali_LS_R60_flag_leaf_EGF_combined[-which(wilcox_res_Anjali_LS_R60_flag_leaf_EGF_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined = as.data.frame(wilcox_sig_log2fc_LS_R60_flag_leaf_EGF_combined[match(rownames(wilcox_sig_stress_responsive_Anjali_LS_R60_flag_leaf_EGF_combined), rownames(wilcox_sig_log2fc_LS_R60_flag_leaf_EGF_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined) = colnames(wilcox_sig_log2fc_LS_R60_flag_leaf_EGF_combined)[3]
wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined$Up.Down = ifelse(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined$`GAR60 - GALS` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined$Up.Down == "Up")    #3 metabolites
sum(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined$Up.Down == "Down")  #7 metabolites


#venn - increased - 60h rewatering/late stress
vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined$Up.Down == "Up"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined) = c("N22", "Dular", "Anjali")
venn_up_stress_responsive_wilcox_LS_R60_flag_leaf_EGF_combined = venn(vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_stress_responsive_R60_latestress_flag_leaf_EGF_combined = venn.diagram(vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined, category.names = c(paste(names(vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined[1]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined[2]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined[3]), " (", length(vennlist_up_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "E", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label E
intersect_stress_responsive_wilcox_LS_R60_venn_up_flag_leaf_EGF_combined = attr(venn_up_stress_responsive_wilcox_LS_R60_flag_leaf_EGF_combined, "intersections")  #list of unique and common metabolites
venn_up_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_EGF_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R60_venn_up_flag_leaf_EGF_combined, rbind))  #convert from list to matrix
colnames(venn_up_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_EGF_combined) = venn_up_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_EGF_combined[1, ]
venn_up_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_EGF_combined = venn_up_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_EGF_combined[-1, ]


#venn - decreased - 60h rewatering/late stress
vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined = list(rownames(wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_N22_LS_R60_flag_leaf_FL_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Dular_LS_R60_flag_leaf_FL_combined$Up.Down == "Down"], rownames(wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined)[wilcox_sig_stress_responsive_log2fc_Anjali_LS_R60_flag_leaf_FL_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined) = c("N22", "Dular", "Anjali")
venn_down_stress_responsive_wilcox_LS_R60_flag_leaf_EGF_combined = venn(vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_stress_responsive_R60_latestress_flag_leaf_EGF_combined = venn.diagram(vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined, category.names = c(paste(names(vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined[1]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined[2]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined[3]), " (", length(vennlist_down_wilcox_stress_responsive_LS_R60_flag_leaf_EGF_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "F", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label F
intersect_stress_responsive_wilcox_LS_R60_venn_down_flag_leaf_EGF_combined = attr(venn_down_stress_responsive_wilcox_LS_R60_flag_leaf_EGF_combined, "intersections")  #list of unique and common metabolites
venn_down_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_EGF_combined = t(ldply(intersect_stress_responsive_wilcox_LS_R60_venn_down_flag_leaf_EGF_combined, rbind))  #convert from list to matrix
colnames(venn_down_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_EGF_combined) = venn_down_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_EGF_combined[1, ]
venn_down_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_EGF_combined = venn_down_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_EGF_combined[-1, ]


write.table(venn_up_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_EGF_combined, file = "wilcox_stress_responsive_sig_up_metabolites_flag_leaf_EGF_R60_latestress_combined.txt", sep = "\t", quote = F)
write.table(venn_down_stress_responsive_wilcox_metabolites_LS_R60_flag_leaf_EGF_combined, file = "wilcox_stress_responsive_sig_down_metabolites_flag_leaf_EGF_R60_latestress_combined.txt", sep = "\t", quote = F)


####################


#significance of wilcoxon test - all cultivars, all three RW timepoints
wilcox_sig_LS_R12_R36_R60_flag_leaf_EGF_combined = cbind(wilcox_res_N22_LS_R12_flag_leaf_EGF_combined,
                                                         wilcox_res_N22_LS_R36_flag_leaf_EGF_combined$Significance,
                                                         wilcox_res_N22_LS_R60_flag_leaf_EGF_combined$Significance,
                                                         wilcox_res_Dular_LS_R12_flag_leaf_EGF_combined$Significance,
                                                         wilcox_res_Dular_LS_R36_flag_leaf_EGF_combined$Significance,
                                                         wilcox_res_Dular_LS_R60_flag_leaf_EGF_combined$Significance,
                                                         wilcox_res_Anjali_LS_R12_flag_leaf_EGF_combined$Significance,
                                                         wilcox_res_Anjali_LS_R36_flag_leaf_EGF_combined$Significance,
                                                         wilcox_res_Anjali_LS_R60_flag_leaf_EGF_combined$Significance)
wilcox_sig_LS_R12_R36_R60_flag_leaf_EGF_combined = wilcox_sig_LS_R12_R36_R60_flag_leaf_EGF_combined[ , -1]
colnames(wilcox_sig_LS_R12_R36_R60_flag_leaf_EGF_combined)[1] = "wilcox_res_N22_LS_R12_flag_leaf_EGF_combined$Significance"
wilcox_sig_LS_R12_R36_R60_flag_leaf_EGF_combined = data.frame(lapply(wilcox_sig_LS_R12_R36_R60_flag_leaf_EGF_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_LS_R12_R36_R60_flag_leaf_EGF_combined) = rownames(wilcox_res_N22_LS_R12_flag_leaf_EGF_combined)
wilcox_sig_LS_R12_R36_R60_flag_leaf_EGF_combined$Nonsig = rowSums(wilcox_sig_LS_R12_R36_R60_flag_leaf_EGF_combined == " ")  #count number of nonsignificant
wilcox_sig_final_LS_R12_R36_R60_flag_leaf_EGF_combined = wilcox_sig_LS_R12_R36_R60_flag_leaf_EGF_combined[-which(wilcox_sig_LS_R12_R36_R60_flag_leaf_EGF_combined$Nonsig == "9"), ] #metabolites significant in at least one of the comparisons - 39 metabolites

#log2fc values of metabolites with significant difference in any one of the comparisons - 12, 36, 60h recovery vs. late stress in any cultivar
wilcox_sig_log2fc_LS_R12_R36_R60_flag_leaf_EGF_combined = log2fc_combined_latestress_recovery_EGF_flag_leaf[which(rownames(log2fc_combined_latestress_recovery_EGF_flag_leaf) %in% rownames(wilcox_sig_final_LS_R12_R36_R60_flag_leaf_EGF_combined)), ] 


#visualization
png("wilcox_sig_stress_responsive_log2fc_LS_RW_flag_leaf_EGF_combined.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_log2fc_LS_R12_R36_R60_flag_leaf_EGF_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-3, -0.01, length.out = 300), seq(0.01, 3, length.out = 300)), 
          margins = c(7.5, 20),
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
          cellnote = wilcox_sig_final_LS_R12_R36_R60_flag_leaf_EGF_combined, 
          notecol = "black", 
          notecex = 1.8, 
          sepcolor = "black", 
          colsep = c(seq(0, 9)), 
          rowsep = c(seq(0, 39)), 
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
png("venn_up_down_flag_leaf_EGF_stress_responsive_RW_latestress_combined.png", width = 16*300, height = 24*300, res = 300)
grid.arrange(gTree(children = venn_up_stress_responsive_R12_latestress_flag_leaf_EGF_combined), 
             gTree(children = venn_down_stress_responsive_R12_latestress_flag_leaf_EGF_combined), 
             gTree(children = venn_up_stress_responsive_R36_latestress_flag_leaf_EGF_combined),
             gTree(children = venn_down_stress_responsive_R36_latestress_flag_leaf_EGF_combined),
             gTree(children = venn_up_stress_responsive_R60_latestress_flag_leaf_EGF_combined),
             gTree(children = venn_down_stress_responsive_R60_latestress_flag_leaf_EGF_combined),
             ncol = 2, nrow = 3)
dev.off()


####################


#recovery vs. corresponding control - NOTE: consider all metabolites, not just the stress-responsive

#log2-fold change - control vs. rewatering
control_recovery_order_EGF = c("GNC12", "GNC36", "GNC60", "GNR12", "GNR36", "GNR60", "GDC12", "GDC36", "GDC60", "GDR12", "GDR36", "GDR60", "GAC12", "GAC36", "GAC60", "GAR12", "GAR36", "GAR60")
data_combined_control_recovery_log2_mean_EGF_flag_leaf = data_combined_log2_mean_EGF_flag_leaf_2[control_recovery_order_EGF, ]  #control and recovery data in preferred order
data_combined_control_recovery_log2_mean_EGF_flag_leaf = t(data_combined_control_recovery_log2_mean_EGF_flag_leaf)

log2fc_combined_N22_control_recovery_EGF_flag_leaf = data_combined_control_recovery_log2_mean_EGF_flag_leaf[ , 4:6] - data_combined_control_recovery_log2_mean_EGF_flag_leaf[ , 1:3]
log2fc_combined_Dular_control_recovery_EGF_flag_leaf = data_combined_control_recovery_log2_mean_EGF_flag_leaf[ , 10:12] - data_combined_control_recovery_log2_mean_EGF_flag_leaf[ , 7:9]
log2fc_combined_Anjali_control_recovery_EGF_flag_leaf = data_combined_control_recovery_log2_mean_EGF_flag_leaf[ , 16:18] - data_combined_control_recovery_log2_mean_EGF_flag_leaf[ , 13:15]

#combine in one data set
log2fc_combined_control_recovery_EGF_flag_leaf = cbind(log2fc_combined_N22_control_recovery_EGF_flag_leaf, 
                                                       log2fc_combined_Dular_control_recovery_EGF_flag_leaf, 
                                                       log2fc_combined_Anjali_control_recovery_EGF_flag_leaf)   #combine log2-fc of all cultivars in one dataset
#rename columns
colnames(log2fc_combined_control_recovery_EGF_flag_leaf) = c(paste(colnames(log2fc_combined_control_recovery_EGF_flag_leaf), "-", sub("R", "C", colnames(log2fc_combined_control_recovery_EGF_flag_leaf))))

#heatmap including all metabolites
heatmap.2(log2fc_combined_control_recovery_EGF_flag_leaf, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          col = bluered(599),
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(7, 20),
          main = "Flag leaf_EGF - Combined - Log2 FC, Recovery/Control", 
          cexCol = 1, 
          cexRow = 0.9, 
          density.info = "none",
          key.xlab = NA,
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          srtCol = 90)


####################


#control 12h vs. 12h rewatering

#N22

#data distribution
shapiro_res_combined_N22_control_recovery12_EGF_flag_leaf = t(data.frame(lapply(data_combined_N22_log2median_EGF_flag_leaf[which(data_combined_N22_log2median_EGF_flag_leaf$Timepoint %in% c("Control.12h Rewatering", "12h Rewatering")), 3:ncol(data_combined_N22_log2median_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_control_recovery12_EGF_flag_leaf) = "p.value"
shapiro_res_combined_N22_control_recovery12_EGF_flag_leaf = as.data.frame(shapiro_res_combined_N22_control_recovery12_EGF_flag_leaf)
shapiro_res_combined_N22_control_recovery12_EGF_flag_leaf$Distribution = ifelse(shapiro_res_combined_N22_control_recovery12_EGF_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_control_recovery12_EGF_flag_leaf$Distribution == "Normal")     #58 metabolites
sum(shapiro_res_combined_N22_control_recovery12_EGF_flag_leaf$Distribution == "Non-normal") #23 metabolites

#wilcox test
wilcox_res_N22_C12_R12_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_N22_log2median_EGF_flag_leaf[ , 3:ncol(data_combined_N22_log2median_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_N22_log2median_EGF_flag_leaf$Timepoint, subset = data_combined_N22_log2median_EGF_flag_leaf$Timepoint %in% c("Control.12h Rewatering", "12h Rewatering"))$p.value)))
colnames(wilcox_res_N22_C12_R12_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_N22_C12_R12_flag_leaf_EGF_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_N22_C12_R12_flag_leaf_EGF_combined = as.data.frame(wilcox_res_N22_C12_R12_flag_leaf_EGF_combined)
wilcox_res_N22_C12_R12_flag_leaf_EGF_combined$Significance = with(wilcox_res_N22_C12_R12_flag_leaf_EGF_combined,
                                                                  ifelse(p.value <= 0.001, "***", ifelse(
                                                                    p.value <= 0.01, "**", ifelse(
                                                                      p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_C12_R12_flag_leaf_EGF_combined$Significance == "*") #7 metabolites
sum(wilcox_res_N22_C12_R12_flag_leaf_EGF_combined$Significance == "**") #5 metabolites
sum(wilcox_res_N22_C12_R12_flag_leaf_EGF_combined$Significance == "***") #4 metabolites
sum(wilcox_res_N22_C12_R12_flag_leaf_EGF_combined$Significance != "ns") #16 metabolites
sum(wilcox_res_N22_C12_R12_flag_leaf_EGF_combined$Significance == "ns") #65 metabolites


#Dular

#data distribution
shapiro_res_combined_Dular_control_recovery12_EGF_flag_leaf = t(data.frame(lapply(data_combined_Dular_log2median_EGF_flag_leaf[which(data_combined_Dular_log2median_EGF_flag_leaf$Timepoint %in% c("Control.12h Rewatering", "12h Rewatering")), 3:ncol(data_combined_Dular_log2median_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_control_recovery12_EGF_flag_leaf) = "p.value"
shapiro_res_combined_Dular_control_recovery12_EGF_flag_leaf = as.data.frame(shapiro_res_combined_Dular_control_recovery12_EGF_flag_leaf)
shapiro_res_combined_Dular_control_recovery12_EGF_flag_leaf$Distribution = ifelse(shapiro_res_combined_Dular_control_recovery12_EGF_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_control_recovery12_EGF_flag_leaf$Distribution == "Normal")     #58 metabolites
sum(shapiro_res_combined_Dular_control_recovery12_EGF_flag_leaf$Distribution == "Non-normal") #23 metabolites

#wilcox test
wilcox_res_Dular_C12_R12_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_Dular_log2median_EGF_flag_leaf[ , 3:ncol(data_combined_Dular_log2median_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_EGF_flag_leaf$Timepoint, subset = data_combined_Dular_log2median_EGF_flag_leaf$Timepoint %in% c("Control.12h Rewatering", "12h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_C12_R12_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_Dular_C12_R12_flag_leaf_EGF_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_Dular_C12_R12_flag_leaf_EGF_combined = as.data.frame(wilcox_res_Dular_C12_R12_flag_leaf_EGF_combined)
wilcox_res_Dular_C12_R12_flag_leaf_EGF_combined$Significance = with(wilcox_res_Dular_C12_R12_flag_leaf_EGF_combined,
                                                                    ifelse(p.value <= 0.001, "***", ifelse(
                                                                      p.value <= 0.01, "**", ifelse(
                                                                        p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_C12_R12_flag_leaf_EGF_combined$Significance == "*") #7 metabolites
sum(wilcox_res_Dular_C12_R12_flag_leaf_EGF_combined$Significance == "**") #4 metabolites
sum(wilcox_res_Dular_C12_R12_flag_leaf_EGF_combined$Significance == "***") #6 metabolites
sum(wilcox_res_Dular_C12_R12_flag_leaf_EGF_combined$Significance != "ns") #17 metabolites
sum(wilcox_res_Dular_C12_R12_flag_leaf_EGF_combined$Significance == "ns") #64 metabolites


#Anjali

#data distribution
shapiro_res_combined_Anjali_control_recovery12_EGF_flag_leaf = t(data.frame(lapply(data_combined_Anjali_log2median_EGF_flag_leaf[which(data_combined_Anjali_log2median_EGF_flag_leaf$Timepoint %in% c("Control.12h Rewatering", "12h Rewatering")), 3:ncol(data_combined_Anjali_log2median_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_control_recovery12_EGF_flag_leaf) = "p.value"
shapiro_res_combined_Anjali_control_recovery12_EGF_flag_leaf = as.data.frame(shapiro_res_combined_Anjali_control_recovery12_EGF_flag_leaf)
shapiro_res_combined_Anjali_control_recovery12_EGF_flag_leaf$Distribution = ifelse(shapiro_res_combined_Anjali_control_recovery12_EGF_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_control_recovery12_EGF_flag_leaf$Distribution == "Normal")     #50 metabolites
sum(shapiro_res_combined_Anjali_control_recovery12_EGF_flag_leaf$Distribution == "Non-normal") #31 metabolites

#wilcox test
wilcox_res_Anjali_C12_R12_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_Anjali_log2median_EGF_flag_leaf[ , 3:ncol(data_combined_Anjali_log2median_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_EGF_flag_leaf$Timepoint, subset = data_combined_Anjali_log2median_EGF_flag_leaf$Timepoint %in% c("Control.12h Rewatering", "12h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_C12_R12_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_Anjali_C12_R12_flag_leaf_EGF_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_Anjali_C12_R12_flag_leaf_EGF_combined = as.data.frame(wilcox_res_Anjali_C12_R12_flag_leaf_EGF_combined)
wilcox_res_Anjali_C12_R12_flag_leaf_EGF_combined$Significance = with(wilcox_res_Anjali_C12_R12_flag_leaf_EGF_combined,
                                                                     ifelse(p.value <= 0.001, "***", ifelse(
                                                                       p.value <= 0.01, "**", ifelse(
                                                                         p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_C12_R12_flag_leaf_EGF_combined$Significance == "*") #7 metabolites
sum(wilcox_res_Anjali_C12_R12_flag_leaf_EGF_combined$Significance == "**") #2 metabolites
sum(wilcox_res_Anjali_C12_R12_flag_leaf_EGF_combined$Significance == "***") #7 metabolites
sum(wilcox_res_Anjali_C12_R12_flag_leaf_EGF_combined$Significance != "ns") #16 metabolites
sum(wilcox_res_Anjali_C12_R12_flag_leaf_EGF_combined$Significance == "ns") #65 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_C12_R12_flag_leaf_EGF_combined = cbind(wilcox_res_N22_C12_R12_flag_leaf_EGF_combined, 
                                                  wilcox_res_Dular_C12_R12_flag_leaf_EGF_combined$Significance, 
                                                  wilcox_res_Anjali_C12_R12_flag_leaf_EGF_combined$Significance)
wilcox_sig_C12_R12_flag_leaf_EGF_combined = wilcox_sig_C12_R12_flag_leaf_EGF_combined[ , -1]
colnames(wilcox_sig_C12_R12_flag_leaf_EGF_combined)[1] = "wilcox_res_N22_C12_R12_flag_leaf_EGF_combined$Significance"
wilcox_sig_C12_R12_flag_leaf_EGF_combined = data.frame(lapply(wilcox_sig_C12_R12_flag_leaf_EGF_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_C12_R12_flag_leaf_EGF_combined) = rownames(wilcox_res_N22_C12_R12_flag_leaf_EGF_combined)
wilcox_sig_C12_R12_flag_leaf_EGF_combined$Nonsig = rowSums(wilcox_sig_C12_R12_flag_leaf_EGF_combined == " ")  #count number of nonsignificant
wilcox_sig_final_C12_R12_flag_leaf_EGF_combined = wilcox_sig_C12_R12_flag_leaf_EGF_combined[-which(wilcox_sig_C12_R12_flag_leaf_EGF_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 32 metabolites

wilcox_sig_log2fc_C12_R12_flag_leaf_EGF_combined = log2fc_combined_control_recovery_EGF_flag_leaf[which(rownames(log2fc_combined_control_recovery_EGF_flag_leaf) %in% rownames(wilcox_sig_final_C12_R12_flag_leaf_EGF_combined)), c(1, 4, 7)]  #log2FC values of metabolites significant in at least one of the cultivars

#overview - with HCA - Euclidean distance and average linkage 
png("wilcox_sig_log2fc_C12_R12_flag_leaf_EGF_combined.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_log2fc_C12_R12_flag_leaf_EGF_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(5, 20), 
          main = "Flag leaf-Combined-Wilcoxon \nLog2 FC, EGF, 12h RW/Control 12h", 
          cexCol = 1.5, 
          cexRow = 0.9, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 8), 
          lwid = c(1.5, 5), 
          labCol = c("N22", "Dular", "Anjali"), 
          srtCol = 0, adjCol = c(0.5, 0.5), 
          cellnote = wilcox_sig_final_C12_R12_flag_leaf_EGF_combined, 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 3)), 
          rowsep = c(seq(0, 32)),
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))
dev.off()


#venn diagram - common among cultivars

#N22
wilcox_sig_N22_C12_R12_flag_leaf_EGF_combined = wilcox_res_N22_C12_R12_flag_leaf_EGF_combined[-which(wilcox_res_N22_C12_R12_flag_leaf_EGF_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_N22_C12_R12_flag_leaf_EGF_combined = as.data.frame(wilcox_sig_log2fc_C12_R12_flag_leaf_EGF_combined[match(rownames(wilcox_sig_N22_C12_R12_flag_leaf_EGF_combined), rownames(wilcox_sig_log2fc_C12_R12_flag_leaf_EGF_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_N22_C12_R12_flag_leaf_EGF_combined) = colnames(wilcox_sig_log2fc_C12_R12_flag_leaf_EGF_combined)[1]
wilcox_sig_log2fc_N22_C12_R12_flag_leaf_EGF_combined$Up.Down = ifelse(wilcox_sig_log2fc_N22_C12_R12_flag_leaf_EGF_combined$`GNR12 - GNC12` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_N22_C12_R12_flag_leaf_EGF_combined$Up.Down == "Up")    #6 metabolites
sum(wilcox_sig_log2fc_N22_C12_R12_flag_leaf_EGF_combined$Up.Down == "Down")  #10 metabolites

#Dular
wilcox_sig_Dular_C12_R12_flag_leaf_EGF_combined = wilcox_res_Dular_C12_R12_flag_leaf_EGF_combined[-which(wilcox_res_Dular_C12_R12_flag_leaf_EGF_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Dular_C12_R12_flag_leaf_EGF_combined = as.data.frame(wilcox_sig_log2fc_C12_R12_flag_leaf_EGF_combined[match(rownames(wilcox_sig_Dular_C12_R12_flag_leaf_EGF_combined), rownames(wilcox_sig_log2fc_C12_R12_flag_leaf_EGF_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Dular_C12_R12_flag_leaf_EGF_combined) = colnames(wilcox_sig_log2fc_C12_R12_flag_leaf_EGF_combined)[2]
wilcox_sig_log2fc_Dular_C12_R12_flag_leaf_EGF_combined$Up.Down = ifelse(wilcox_sig_log2fc_Dular_C12_R12_flag_leaf_EGF_combined$`GDR12 - GDC12` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Dular_C12_R12_flag_leaf_EGF_combined$Up.Down == "Up")    #10 metabolites
sum(wilcox_sig_log2fc_Dular_C12_R12_flag_leaf_EGF_combined$Up.Down == "Down")  #7 metabolites

#Anjali
wilcox_sig_Anjali_C12_R12_flag_leaf_EGF_combined = wilcox_res_Anjali_C12_R12_flag_leaf_EGF_combined[-which(wilcox_res_Anjali_C12_R12_flag_leaf_EGF_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Anjali_C12_R12_flag_leaf_EGF_combined = as.data.frame(wilcox_sig_log2fc_C12_R12_flag_leaf_EGF_combined[match(rownames(wilcox_sig_Anjali_C12_R12_flag_leaf_EGF_combined), rownames(wilcox_sig_log2fc_C12_R12_flag_leaf_EGF_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Anjali_C12_R12_flag_leaf_EGF_combined) = colnames(wilcox_sig_log2fc_C12_R12_flag_leaf_EGF_combined)[3]
wilcox_sig_log2fc_Anjali_C12_R12_flag_leaf_EGF_combined$Up.Down = ifelse(wilcox_sig_log2fc_Anjali_C12_R12_flag_leaf_EGF_combined$`GAR12 - GAC12` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Anjali_C12_R12_flag_leaf_EGF_combined$Up.Down == "Up")    #5 metabolites
sum(wilcox_sig_log2fc_Anjali_C12_R12_flag_leaf_EGF_combined$Up.Down == "Down")  #11 metabolites


#venn - increased - 12h rewatering/control 12h
vennlist_up_wilcox_C12_R12_flag_leaf_EGF_combined = list(rownames(wilcox_sig_log2fc_N22_C12_R12_flag_leaf_EGF_combined)[wilcox_sig_log2fc_N22_C12_R12_flag_leaf_EGF_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Dular_C12_R12_flag_leaf_EGF_combined)[wilcox_sig_log2fc_Dular_C12_R12_flag_leaf_EGF_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Anjali_C12_R12_flag_leaf_EGF_combined)[wilcox_sig_log2fc_Anjali_C12_R12_flag_leaf_EGF_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_C12_R12_flag_leaf_EGF_combined) = c("N22", "Dular", "Anjali")
venn_up_wilcox_C12_R12_flag_leaf_EGF_combined = venn(vennlist_up_wilcox_C12_R12_flag_leaf_EGF_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_R12_C12_flag_leaf_EGF_combined = venn.diagram(vennlist_up_wilcox_C12_R12_flag_leaf_EGF_combined, category.names = c(paste(names(vennlist_up_wilcox_C12_R12_flag_leaf_EGF_combined[1]), " (", length(vennlist_up_wilcox_C12_R12_flag_leaf_EGF_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_C12_R12_flag_leaf_EGF_combined[2]), " (", length(vennlist_up_wilcox_C12_R12_flag_leaf_EGF_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_C12_R12_flag_leaf_EGF_combined[3]), " (", length(vennlist_up_wilcox_C12_R12_flag_leaf_EGF_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "A", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label A
intersect_wilcox_C12_R12_venn_up_flag_leaf_EGF_combined = attr(venn_up_wilcox_C12_R12_flag_leaf_EGF_combined, "intersections")  #list of unique and common metabolites
venn_up_wilcox_metabolites_C12_R12_flag_leaf_EGF_combined = t(ldply(intersect_wilcox_C12_R12_venn_up_flag_leaf_EGF_combined, rbind))  #convert from list to matrix
colnames(venn_up_wilcox_metabolites_C12_R12_flag_leaf_EGF_combined) = venn_up_wilcox_metabolites_C12_R12_flag_leaf_EGF_combined[1, ]
venn_up_wilcox_metabolites_C12_R12_flag_leaf_EGF_combined = venn_up_wilcox_metabolites_C12_R12_flag_leaf_EGF_combined[-1, ]


#venn - decreased - 12h rewatering/control 12h
vennlist_down_wilcox_C12_R12_flag_leaf_EGF_combined = list(rownames(wilcox_sig_log2fc_N22_C12_R12_flag_leaf_EGF_combined)[wilcox_sig_log2fc_N22_C12_R12_flag_leaf_EGF_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Dular_C12_R12_flag_leaf_EGF_combined)[wilcox_sig_log2fc_Dular_C12_R12_flag_leaf_EGF_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Anjali_C12_R12_flag_leaf_EGF_combined)[wilcox_sig_log2fc_Anjali_C12_R12_flag_leaf_EGF_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_C12_R12_flag_leaf_EGF_combined) = c("N22", "Dular", "Anjali")
venn_down_wilcox_C12_R12_flag_leaf_EGF_combined = venn(vennlist_down_wilcox_C12_R12_flag_leaf_EGF_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_R12_C12_flag_leaf_EGF_combined = venn.diagram(vennlist_down_wilcox_C12_R12_flag_leaf_EGF_combined, category.names = c(paste(names(vennlist_down_wilcox_C12_R12_flag_leaf_EGF_combined[1]), " (", length(vennlist_down_wilcox_C12_R12_flag_leaf_EGF_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_C12_R12_flag_leaf_EGF_combined[2]), " (", length(vennlist_down_wilcox_C12_R12_flag_leaf_EGF_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_C12_R12_flag_leaf_EGF_combined[3]), " (", length(vennlist_down_wilcox_C12_R12_flag_leaf_EGF_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "B", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label B
intersect_wilcox_C12_R12_venn_down_flag_leaf_EGF_combined = attr(venn_down_wilcox_C12_R12_flag_leaf_EGF_combined, "intersections")  #list of unique and common metabolites
venn_down_wilcox_metabolites_C12_R12_flag_leaf_EGF_combined = t(ldply(intersect_wilcox_C12_R12_venn_down_flag_leaf_EGF_combined, rbind))  #convert from list to matrix
colnames(venn_down_wilcox_metabolites_C12_R12_flag_leaf_EGF_combined) = venn_down_wilcox_metabolites_C12_R12_flag_leaf_EGF_combined[1, ]
venn_down_wilcox_metabolites_C12_R12_flag_leaf_EGF_combined = venn_down_wilcox_metabolites_C12_R12_flag_leaf_EGF_combined[-1, ]


write.table(venn_up_wilcox_metabolites_C12_R12_flag_leaf_EGF_combined, file = "wilcox_sig_up_metabolites_flag_leaf_EGF_R12_C12_combined.txt", sep = "\t", quote = F)
write.table(venn_down_wilcox_metabolites_C12_R12_flag_leaf_EGF_combined, file = "wilcox_sig_down_metabolites_flag_leaf_EGF_R12_C12_combined.txt", sep = "\t", quote = F)


####################

#Control 36h vs. 36h rewatering

#N22

#data distribution
shapiro_res_combined_N22_control_recovery36_EGF_flag_leaf = t(data.frame(lapply(data_combined_N22_log2median_EGF_flag_leaf[which(data_combined_N22_log2median_EGF_flag_leaf$Timepoint %in% c("Control.36h Rewatering", "36h Rewatering")), 3:ncol(data_combined_N22_log2median_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_control_recovery36_EGF_flag_leaf) = "p.value"
shapiro_res_combined_N22_control_recovery36_EGF_flag_leaf = as.data.frame(shapiro_res_combined_N22_control_recovery36_EGF_flag_leaf)
shapiro_res_combined_N22_control_recovery36_EGF_flag_leaf$Distribution = ifelse(shapiro_res_combined_N22_control_recovery36_EGF_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_control_recovery36_EGF_flag_leaf$Distribution == "Normal")     #62 metabolites
sum(shapiro_res_combined_N22_control_recovery36_EGF_flag_leaf$Distribution == "Non-normal") #19 metabolites

#wilcox test
wilcox_res_N22_C36_R36_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_N22_log2median_EGF_flag_leaf[ , 3:ncol(data_combined_N22_log2median_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_N22_log2median_EGF_flag_leaf$Timepoint, subset = data_combined_N22_log2median_EGF_flag_leaf$Timepoint %in% c("Control.36h Rewatering", "36h Rewatering"))$p.value)))
colnames(wilcox_res_N22_C36_R36_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_N22_C36_R36_flag_leaf_EGF_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_N22_C36_R36_flag_leaf_EGF_combined = as.data.frame(wilcox_res_N22_C36_R36_flag_leaf_EGF_combined)
wilcox_res_N22_C36_R36_flag_leaf_EGF_combined$Significance = with(wilcox_res_N22_C36_R36_flag_leaf_EGF_combined,
                                                                  ifelse(p.value <= 0.001, "***", ifelse(
                                                                    p.value <= 0.01, "**", ifelse(
                                                                      p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_C36_R36_flag_leaf_EGF_combined$Significance == "*") #14 metabolites
sum(wilcox_res_N22_C36_R36_flag_leaf_EGF_combined$Significance == "**") #3 metabolites
sum(wilcox_res_N22_C36_R36_flag_leaf_EGF_combined$Significance == "***") #1 metabolites
sum(wilcox_res_N22_C36_R36_flag_leaf_EGF_combined$Significance != "ns") #18 metabolites
sum(wilcox_res_N22_C36_R36_flag_leaf_EGF_combined$Significance == "ns") #63 metabolites


#Dular

#data distribution
shapiro_res_combined_Dular_control_recovery36_EGF_flag_leaf = t(data.frame(lapply(data_combined_Dular_log2median_EGF_flag_leaf[which(data_combined_Dular_log2median_EGF_flag_leaf$Timepoint %in% c("Control.36h Rewatering", "36h Rewatering")), 3:ncol(data_combined_Dular_log2median_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_control_recovery36_EGF_flag_leaf) = "p.value"
shapiro_res_combined_Dular_control_recovery36_EGF_flag_leaf = as.data.frame(shapiro_res_combined_Dular_control_recovery36_EGF_flag_leaf)
shapiro_res_combined_Dular_control_recovery36_EGF_flag_leaf$Distribution = ifelse(shapiro_res_combined_Dular_control_recovery36_EGF_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_control_recovery36_EGF_flag_leaf$Distribution == "Normal")     #62 metabolites
sum(shapiro_res_combined_Dular_control_recovery36_EGF_flag_leaf$Distribution == "Non-normal") #19 metabolites

#wilcox test
wilcox_res_Dular_C36_R36_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_Dular_log2median_EGF_flag_leaf[ , 3:ncol(data_combined_Dular_log2median_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_EGF_flag_leaf$Timepoint, subset = data_combined_Dular_log2median_EGF_flag_leaf$Timepoint %in% c("Control.36h Rewatering", "36h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_C36_R36_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_Dular_C36_R36_flag_leaf_EGF_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_Dular_C36_R36_flag_leaf_EGF_combined = as.data.frame(wilcox_res_Dular_C36_R36_flag_leaf_EGF_combined)
wilcox_res_Dular_C36_R36_flag_leaf_EGF_combined$Significance = with(wilcox_res_Dular_C36_R36_flag_leaf_EGF_combined,
                                                                    ifelse(p.value <= 0.001, "***", ifelse(
                                                                      p.value <= 0.01, "**", ifelse(
                                                                        p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_C36_R36_flag_leaf_EGF_combined$Significance == "*") #9 metabolites
sum(wilcox_res_Dular_C36_R36_flag_leaf_EGF_combined$Significance == "**") #3 metabolites
sum(wilcox_res_Dular_C36_R36_flag_leaf_EGF_combined$Significance == "***") #2 metabolites
sum(wilcox_res_Dular_C36_R36_flag_leaf_EGF_combined$Significance != "ns") #14 metabolites
sum(wilcox_res_Dular_C36_R36_flag_leaf_EGF_combined$Significance == "ns") #67 metabolites


#Anjali

#data distribution
shapiro_res_combined_Anjali_control_recovery36_EGF_flag_leaf = t(data.frame(lapply(data_combined_Anjali_log2median_EGF_flag_leaf[which(data_combined_Anjali_log2median_EGF_flag_leaf$Timepoint %in% c("Control.36h Rewatering", "36h Rewatering")), 3:ncol(data_combined_Anjali_log2median_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_control_recovery36_EGF_flag_leaf) = "p.value"
shapiro_res_combined_Anjali_control_recovery36_EGF_flag_leaf = as.data.frame(shapiro_res_combined_Anjali_control_recovery36_EGF_flag_leaf)
shapiro_res_combined_Anjali_control_recovery36_EGF_flag_leaf$Distribution = ifelse(shapiro_res_combined_Anjali_control_recovery36_EGF_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_control_recovery36_EGF_flag_leaf$Distribution == "Normal")     #66 metabolites
sum(shapiro_res_combined_Anjali_control_recovery36_EGF_flag_leaf$Distribution == "Non-normal") #15 metabolites

#wilcox test
wilcox_res_Anjali_C36_R36_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_Anjali_log2median_EGF_flag_leaf[ , 3:ncol(data_combined_Anjali_log2median_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_EGF_flag_leaf$Timepoint, subset = data_combined_Anjali_log2median_EGF_flag_leaf$Timepoint %in% c("Control.36h Rewatering", "36h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_C36_R36_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_Anjali_C36_R36_flag_leaf_EGF_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_Anjali_C36_R36_flag_leaf_EGF_combined = as.data.frame(wilcox_res_Anjali_C36_R36_flag_leaf_EGF_combined)
wilcox_res_Anjali_C36_R36_flag_leaf_EGF_combined$Significance = with(wilcox_res_Anjali_C36_R36_flag_leaf_EGF_combined,
                                                                     ifelse(p.value <= 0.001, "***", ifelse(
                                                                       p.value <= 0.01, "**", ifelse(
                                                                         p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_C36_R36_flag_leaf_EGF_combined$Significance == "*") #3 metabolites
sum(wilcox_res_Anjali_C36_R36_flag_leaf_EGF_combined$Significance == "**") #5 metabolites
sum(wilcox_res_Anjali_C36_R36_flag_leaf_EGF_combined$Significance == "***") #9 metabolites
sum(wilcox_res_Anjali_C36_R36_flag_leaf_EGF_combined$Significance != "ns") #17 metabolites
sum(wilcox_res_Anjali_C36_R36_flag_leaf_EGF_combined$Significance == "ns") #64 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_C36_R36_flag_leaf_EGF_combined = cbind(wilcox_res_N22_C36_R36_flag_leaf_EGF_combined, 
                                                  wilcox_res_Dular_C36_R36_flag_leaf_EGF_combined$Significance,
                                                  wilcox_res_Anjali_C36_R36_flag_leaf_EGF_combined$Significance)
wilcox_sig_C36_R36_flag_leaf_EGF_combined = wilcox_sig_C36_R36_flag_leaf_EGF_combined[ , -1]
colnames(wilcox_sig_C36_R36_flag_leaf_EGF_combined)[1] = "wilcox_res_N22_C36_R36_flag_leaf_EGF_combined$Significance"
wilcox_sig_C36_R36_flag_leaf_EGF_combined = data.frame(lapply(wilcox_sig_C36_R36_flag_leaf_EGF_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_C36_R36_flag_leaf_EGF_combined) = rownames(wilcox_res_N22_C36_R36_flag_leaf_EGF_combined)
wilcox_sig_C36_R36_flag_leaf_EGF_combined$Nonsig = rowSums(wilcox_sig_C36_R36_flag_leaf_EGF_combined == " ")  #count number of nonsignificant
wilcox_sig_final_C36_R36_flag_leaf_EGF_combined = wilcox_sig_C36_R36_flag_leaf_EGF_combined[-which(wilcox_sig_C36_R36_flag_leaf_EGF_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 31 metabolites

wilcox_sig_log2fc_C36_R36_flag_leaf_EGF_combined = log2fc_combined_control_recovery_EGF_flag_leaf[which(rownames(log2fc_combined_control_recovery_EGF_flag_leaf) %in% rownames(wilcox_sig_final_C36_R36_flag_leaf_EGF_combined)), c(2, 5, 8)]  #log2FC values of metabolites significant in at least one of the cultivars

#overview - with HCA - Euclidean distance and average linkage 
png("wilcox_sig_log2fc_C36_R36_flag_leaf_EGF_combined.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_log2fc_C36_R36_flag_leaf_EGF_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row",
          trace = "none", 
          col = bluered(599),
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)),
          margins = c(5, 20), 
          main = "Flag leaf-Combined-Wilcoxon \nLog2 FC, EGF, 36h RW/Control 36h", 
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
          cellnote = wilcox_sig_final_C36_R36_flag_leaf_EGF_combined, 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 3)), 
          rowsep = c(seq(0, 31)), 
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))
dev.off()


#venn diagram - common among cultivars

#N22
wilcox_sig_N22_C36_R36_flag_leaf_EGF_combined = wilcox_res_N22_C36_R36_flag_leaf_EGF_combined[-which(wilcox_res_N22_C36_R36_flag_leaf_EGF_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_N22_C36_R36_flag_leaf_EGF_combined = as.data.frame(wilcox_sig_log2fc_C36_R36_flag_leaf_EGF_combined[match(rownames(wilcox_sig_N22_C36_R36_flag_leaf_EGF_combined), rownames(wilcox_sig_log2fc_C36_R36_flag_leaf_EGF_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_N22_C36_R36_flag_leaf_EGF_combined) = colnames(wilcox_sig_log2fc_C36_R36_flag_leaf_EGF_combined)[1]
wilcox_sig_log2fc_N22_C36_R36_flag_leaf_EGF_combined$Up.Down = ifelse(wilcox_sig_log2fc_N22_C36_R36_flag_leaf_EGF_combined$`GNR36 - GNC36` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_N22_C36_R36_flag_leaf_EGF_combined$Up.Down == "Up")    #7 metabolites
sum(wilcox_sig_log2fc_N22_C36_R36_flag_leaf_EGF_combined$Up.Down == "Down")  #11 metabolites

#Dular
wilcox_sig_Dular_C36_R36_flag_leaf_EGF_combined = wilcox_res_Dular_C36_R36_flag_leaf_EGF_combined[-which(wilcox_res_Dular_C36_R36_flag_leaf_EGF_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Dular_C36_R36_flag_leaf_EGF_combined = as.data.frame(wilcox_sig_log2fc_C36_R36_flag_leaf_EGF_combined[match(rownames(wilcox_sig_Dular_C36_R36_flag_leaf_EGF_combined), rownames(wilcox_sig_log2fc_C36_R36_flag_leaf_EGF_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Dular_C36_R36_flag_leaf_EGF_combined) = colnames(wilcox_sig_log2fc_C36_R36_flag_leaf_EGF_combined)[2]
wilcox_sig_log2fc_Dular_C36_R36_flag_leaf_EGF_combined$Up.Down = ifelse(wilcox_sig_log2fc_Dular_C36_R36_flag_leaf_EGF_combined$`GDR36 - GDC36` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Dular_C36_R36_flag_leaf_EGF_combined$Up.Down == "Up")    #3 metabolites
sum(wilcox_sig_log2fc_Dular_C36_R36_flag_leaf_EGF_combined$Up.Down == "Down")  #11 metabolites

#Anjali
wilcox_sig_Anjali_C36_R36_flag_leaf_EGF_combined = wilcox_res_Anjali_C36_R36_flag_leaf_EGF_combined[-which(wilcox_res_Anjali_C36_R36_flag_leaf_EGF_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Anjali_C36_R36_flag_leaf_EGF_combined = as.data.frame(wilcox_sig_log2fc_C36_R36_flag_leaf_EGF_combined[match(rownames(wilcox_sig_Anjali_C36_R36_flag_leaf_EGF_combined), rownames(wilcox_sig_log2fc_C36_R36_flag_leaf_EGF_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Anjali_C36_R36_flag_leaf_EGF_combined) = colnames(wilcox_sig_log2fc_C36_R36_flag_leaf_EGF_combined)[3]
wilcox_sig_log2fc_Anjali_C36_R36_flag_leaf_EGF_combined$Up.Down = ifelse(wilcox_sig_log2fc_Anjali_C36_R36_flag_leaf_EGF_combined$`GAR36 - GAC36` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Anjali_C36_R36_flag_leaf_EGF_combined$Up.Down == "Up")    #10 metabolites
sum(wilcox_sig_log2fc_Anjali_C36_R36_flag_leaf_EGF_combined$Up.Down == "Down")  #7 metabolites


#venn - increased - 36h rewatering/control 36h
vennlist_up_wilcox_C36_R36_flag_leaf_EGF_combined = list(rownames(wilcox_sig_log2fc_N22_C36_R36_flag_leaf_EGF_combined)[wilcox_sig_log2fc_N22_C36_R36_flag_leaf_EGF_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Dular_C36_R36_flag_leaf_EGF_combined)[wilcox_sig_log2fc_Dular_C36_R36_flag_leaf_EGF_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Anjali_C36_R36_flag_leaf_EGF_combined)[wilcox_sig_log2fc_Anjali_C36_R36_flag_leaf_EGF_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_C36_R36_flag_leaf_EGF_combined) = c("N22", "Dular", "Anjali")
venn_up_wilcox_C36_R36_flag_leaf_EGF_combined = venn(vennlist_up_wilcox_C36_R36_flag_leaf_EGF_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_R36_C36_flag_leaf_EGF_combined = venn.diagram(vennlist_up_wilcox_C36_R36_flag_leaf_EGF_combined, category.names = c(paste(names(vennlist_up_wilcox_C36_R36_flag_leaf_EGF_combined[1]), " (", length(vennlist_up_wilcox_C36_R36_flag_leaf_EGF_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_C36_R36_flag_leaf_EGF_combined[2]), " (", length(vennlist_up_wilcox_C36_R36_flag_leaf_EGF_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_C36_R36_flag_leaf_EGF_combined[3]), " (", length(vennlist_up_wilcox_C36_R36_flag_leaf_EGF_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "C", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label C
intersect_wilcox_C36_R36_venn_up_flag_leaf_EGF_combined = attr(venn_up_wilcox_C36_R36_flag_leaf_EGF_combined, "intersections")  #list of unique and common metabolites
venn_up_wilcox_metabolites_C36_R36_flag_leaf_EGF_combined = t(ldply(intersect_wilcox_C36_R36_venn_up_flag_leaf_EGF_combined, rbind))  #convert from list to matrix
colnames(venn_up_wilcox_metabolites_C36_R36_flag_leaf_EGF_combined) = venn_up_wilcox_metabolites_C36_R36_flag_leaf_EGF_combined[1, ]
venn_up_wilcox_metabolites_C36_R36_flag_leaf_EGF_combined = venn_up_wilcox_metabolites_C36_R36_flag_leaf_EGF_combined[-1, ]


#venn - decreased - 36h rewatering/control 36h
vennlist_down_wilcox_C36_R36_flag_leaf_EGF_combined = list(rownames(wilcox_sig_log2fc_N22_C36_R36_flag_leaf_EGF_combined)[wilcox_sig_log2fc_N22_C36_R36_flag_leaf_EGF_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Dular_C36_R36_flag_leaf_EGF_combined)[wilcox_sig_log2fc_Dular_C36_R36_flag_leaf_EGF_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Anjali_C36_R36_flag_leaf_EGF_combined)[wilcox_sig_log2fc_Anjali_C36_R36_flag_leaf_EGF_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_C36_R36_flag_leaf_EGF_combined) = c("N22", "Dular", "Anjali")
venn_down_wilcox_C36_R36_flag_leaf_EGF_combined = venn(vennlist_down_wilcox_C36_R36_flag_leaf_EGF_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_R36_C36_flag_leaf_EGF_combined = venn.diagram(vennlist_down_wilcox_C36_R36_flag_leaf_EGF_combined, category.names = c(paste(names(vennlist_down_wilcox_C36_R36_flag_leaf_EGF_combined[1]), " (", length(vennlist_down_wilcox_C36_R36_flag_leaf_EGF_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_C36_R36_flag_leaf_EGF_combined[2]), " (", length(vennlist_down_wilcox_C36_R36_flag_leaf_EGF_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_C36_R36_flag_leaf_EGF_combined[3]), " (", length(vennlist_down_wilcox_C36_R36_flag_leaf_EGF_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "D", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label D
intersect_wilcox_C36_R36_venn_down_flag_leaf_EGF_combined = attr(venn_down_wilcox_C36_R36_flag_leaf_EGF_combined, "intersections")  #list of unique and common metabolites
venn_down_wilcox_metabolites_C36_R36_flag_leaf_EGF_combined = t(ldply(intersect_wilcox_C36_R36_venn_down_flag_leaf_EGF_combined, rbind))  #convert from list to matrix
colnames(venn_down_wilcox_metabolites_C36_R36_flag_leaf_EGF_combined) = venn_down_wilcox_metabolites_C36_R36_flag_leaf_EGF_combined[1, ]
venn_down_wilcox_metabolites_C36_R36_flag_leaf_EGF_combined = venn_down_wilcox_metabolites_C36_R36_flag_leaf_EGF_combined[-1, ]


write.table(venn_up_wilcox_metabolites_C36_R36_flag_leaf_EGF_combined, file = "wilcox_sig_up_metabolites_flag_leaf_EGF_R36_C36_combined.txt", sep = "\t", quote = F)
write.table(venn_down_wilcox_metabolites_C36_R36_flag_leaf_EGF_combined, file = "wilcox_sig_down_metabolites_flag_leaf_EGF_R36_C36_combined.txt", sep = "\t", quote = F)


####################


#Control 60h vs. 60h rewatering

#N22

#data distribution
shapiro_res_combined_N22_control_recovery60_EGF_flag_leaf = t(data.frame(lapply(data_combined_N22_log2median_EGF_flag_leaf[which(data_combined_N22_log2median_EGF_flag_leaf$Timepoint %in% c("Control.60h Rewatering", "60h Rewatering")), 3:ncol(data_combined_N22_log2median_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_N22_control_recovery60_EGF_flag_leaf) = "p.value"
shapiro_res_combined_N22_control_recovery60_EGF_flag_leaf = as.data.frame(shapiro_res_combined_N22_control_recovery60_EGF_flag_leaf)
shapiro_res_combined_N22_control_recovery60_EGF_flag_leaf$Distribution = ifelse(shapiro_res_combined_N22_control_recovery60_EGF_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_N22_control_recovery60_EGF_flag_leaf$Distribution == "Normal")     #62 metabolites
sum(shapiro_res_combined_N22_control_recovery60_EGF_flag_leaf$Distribution == "Non-normal") #19 metabolites

#wilcox test
wilcox_res_N22_C60_R60_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_N22_log2median_EGF_flag_leaf[ , 3:ncol(data_combined_N22_log2median_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_N22_log2median_EGF_flag_leaf$Timepoint, subset = data_combined_N22_log2median_EGF_flag_leaf$Timepoint %in% c("Control.60h Rewatering", "60h Rewatering"))$p.value)))
colnames(wilcox_res_N22_C60_R60_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_N22_C60_R60_flag_leaf_EGF_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_N22_C60_R60_flag_leaf_EGF_combined = as.data.frame(wilcox_res_N22_C60_R60_flag_leaf_EGF_combined)
wilcox_res_N22_C60_R60_flag_leaf_EGF_combined$Significance = with(wilcox_res_N22_C60_R60_flag_leaf_EGF_combined,
                                                                  ifelse(p.value <= 0.001, "***", ifelse(
                                                                    p.value <= 0.01, "**", ifelse(
                                                                      p.value < 0.05, "*", "ns"))))
sum(wilcox_res_N22_C60_R60_flag_leaf_EGF_combined$Significance == "*") #13 metabolites
sum(wilcox_res_N22_C60_R60_flag_leaf_EGF_combined$Significance == "**") #6 metabolites
sum(wilcox_res_N22_C60_R60_flag_leaf_EGF_combined$Significance == "***") #3 metabolites
sum(wilcox_res_N22_C60_R60_flag_leaf_EGF_combined$Significance != "ns") #22 metabolites
sum(wilcox_res_N22_C60_R60_flag_leaf_EGF_combined$Significance == "ns") #59 metabolites


#Dular

#data distribution
shapiro_res_combined_Dular_control_recovery60_EGF_flag_leaf = t(data.frame(lapply(data_combined_Dular_log2median_EGF_flag_leaf[which(data_combined_Dular_log2median_EGF_flag_leaf$Timepoint %in% c("Control.60h Rewatering", "60h Rewatering")), 3:ncol(data_combined_Dular_log2median_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Dular_control_recovery60_EGF_flag_leaf) = "p.value"
shapiro_res_combined_Dular_control_recovery60_EGF_flag_leaf = as.data.frame(shapiro_res_combined_Dular_control_recovery60_EGF_flag_leaf)
shapiro_res_combined_Dular_control_recovery60_EGF_flag_leaf$Distribution = ifelse(shapiro_res_combined_Dular_control_recovery60_EGF_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Dular_control_recovery60_EGF_flag_leaf$Distribution == "Normal")     #51 metabolites
sum(shapiro_res_combined_Dular_control_recovery60_EGF_flag_leaf$Distribution == "Non-normal") #30 metabolites

#wilcox test
wilcox_res_Dular_C60_R60_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_Dular_log2median_EGF_flag_leaf[ , 3:ncol(data_combined_Dular_log2median_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Dular_log2median_EGF_flag_leaf$Timepoint, subset = data_combined_Dular_log2median_EGF_flag_leaf$Timepoint %in% c("Control.60h Rewatering", "60h Rewatering"))$p.value)))
colnames(wilcox_res_Dular_C60_R60_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_Dular_C60_R60_flag_leaf_EGF_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_Dular_C60_R60_flag_leaf_EGF_combined = as.data.frame(wilcox_res_Dular_C60_R60_flag_leaf_EGF_combined)
wilcox_res_Dular_C60_R60_flag_leaf_EGF_combined$Significance = with(wilcox_res_Dular_C60_R60_flag_leaf_EGF_combined,
                                                                    ifelse(p.value <= 0.001, "***", ifelse(
                                                                      p.value <= 0.01, "**", ifelse(
                                                                        p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Dular_C60_R60_flag_leaf_EGF_combined$Significance == "*") #10 metabolites
sum(wilcox_res_Dular_C60_R60_flag_leaf_EGF_combined$Significance == "**") #11 metabolites
sum(wilcox_res_Dular_C60_R60_flag_leaf_EGF_combined$Significance == "***") #6 metabolites
sum(wilcox_res_Dular_C60_R60_flag_leaf_EGF_combined$Significance != "ns") #27 metabolites
sum(wilcox_res_Dular_C60_R60_flag_leaf_EGF_combined$Significance == "ns") #54 metabolites


#Anjali

#data distribution
shapiro_res_combined_Anjali_control_recovery60_EGF_flag_leaf = t(data.frame(lapply(data_combined_Anjali_log2median_EGF_flag_leaf[which(data_combined_Anjali_log2median_EGF_flag_leaf$Timepoint %in% c("Control.60h Rewatering", "60h Rewatering")), 3:ncol(data_combined_Anjali_log2median_EGF_flag_leaf)], function(x) shapiro.test(x)$p.value)))
colnames(shapiro_res_combined_Anjali_control_recovery60_EGF_flag_leaf) = "p.value"
shapiro_res_combined_Anjali_control_recovery60_EGF_flag_leaf = as.data.frame(shapiro_res_combined_Anjali_control_recovery60_EGF_flag_leaf)
shapiro_res_combined_Anjali_control_recovery60_EGF_flag_leaf$Distribution = ifelse(shapiro_res_combined_Anjali_control_recovery60_EGF_flag_leaf$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_res_combined_Anjali_control_recovery60_EGF_flag_leaf$Distribution == "Normal")     #63 metabolites
sum(shapiro_res_combined_Anjali_control_recovery60_EGF_flag_leaf$Distribution == "Non-normal") #18 metabolites

#wilcox test
wilcox_res_Anjali_C60_R60_flag_leaf_EGF_combined = t(data.frame(lapply(data_combined_Anjali_log2median_EGF_flag_leaf[ , 3:ncol(data_combined_Anjali_log2median_EGF_flag_leaf)], function(x) wilcox.test(x ~ data_combined_Anjali_log2median_EGF_flag_leaf$Timepoint, subset = data_combined_Anjali_log2median_EGF_flag_leaf$Timepoint %in% c("Control.60h Rewatering", "60h Rewatering"))$p.value)))
colnames(wilcox_res_Anjali_C60_R60_flag_leaf_EGF_combined) = "p.value"
rownames(wilcox_res_Anjali_C60_R60_flag_leaf_EGF_combined) = reduced_metabolite_list_final_flag_leaf_2013$Name
wilcox_res_Anjali_C60_R60_flag_leaf_EGF_combined = as.data.frame(wilcox_res_Anjali_C60_R60_flag_leaf_EGF_combined)
wilcox_res_Anjali_C60_R60_flag_leaf_EGF_combined$Significance = with(wilcox_res_Anjali_C60_R60_flag_leaf_EGF_combined,
                                                                     ifelse(p.value <= 0.001, "***", ifelse(
                                                                       p.value <= 0.01, "**", ifelse(
                                                                         p.value < 0.05, "*", "ns"))))
sum(wilcox_res_Anjali_C60_R60_flag_leaf_EGF_combined$Significance == "*") #13 metabolites
sum(wilcox_res_Anjali_C60_R60_flag_leaf_EGF_combined$Significance == "**") #8 metabolites
sum(wilcox_res_Anjali_C60_R60_flag_leaf_EGF_combined$Significance == "***") #14 metabolites
sum(wilcox_res_Anjali_C60_R60_flag_leaf_EGF_combined$Significance != "ns") #35 metabolites
sum(wilcox_res_Anjali_C60_R60_flag_leaf_EGF_combined$Significance == "ns") #46 metabolites


#metabolites with Wilcoxon test significant results

wilcox_sig_C60_R60_flag_leaf_EGF_combined = cbind(wilcox_res_N22_C60_R60_flag_leaf_EGF_combined, 
                                                  wilcox_res_Dular_C60_R60_flag_leaf_EGF_combined$Significance, 
                                                  wilcox_res_Anjali_C60_R60_flag_leaf_EGF_combined$Significance)
wilcox_sig_C60_R60_flag_leaf_EGF_combined = wilcox_sig_C60_R60_flag_leaf_EGF_combined[ , -1]
colnames(wilcox_sig_C60_R60_flag_leaf_EGF_combined)[1] = "wilcox_res_N22_C60_R60_flag_leaf_EGF_combined$Significance"
wilcox_sig_C60_R60_flag_leaf_EGF_combined = data.frame(lapply(wilcox_sig_C60_R60_flag_leaf_EGF_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_C60_R60_flag_leaf_EGF_combined) = rownames(wilcox_res_N22_C60_R60_flag_leaf_EGF_combined)
wilcox_sig_C60_R60_flag_leaf_EGF_combined$Nonsig = rowSums(wilcox_sig_C60_R60_flag_leaf_EGF_combined == " ")  #count number of nonsignificant
wilcox_sig_final_C60_R60_flag_leaf_EGF_combined = wilcox_sig_C60_R60_flag_leaf_EGF_combined[-which(wilcox_sig_C60_R60_flag_leaf_EGF_combined$Nonsig == "3"), ] #metabolites significant in at least one of the cultivars - 47 metabolites

wilcox_sig_log2fc_C60_R60_flag_leaf_EGF_combined = log2fc_combined_control_recovery_EGF_flag_leaf[which(rownames(log2fc_combined_control_recovery_EGF_flag_leaf) %in% rownames(wilcox_sig_final_C60_R60_flag_leaf_EGF_combined)), c(3, 6, 9)]  #log2FC values of metabolites significant in at least one of the cultivars

#overview - with HCA - Euclidean distance and average linkage 
png("wilcox_sig_log2fc_C60_R60_flag_leaf_EGF_combined.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_log2fc_C60_R60_flag_leaf_EGF_combined, 
          Rowv = TRUE, 
          Colv = FALSE,
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-2, -0.01, length.out = 300), seq(0.01, 2, length.out = 300)), 
          margins = c(5, 20), 
          main = "Flag leaf-Combined-Wilcoxon \nLog2 FC, EGF, 60h RW/Control 60h", 
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
          cellnote = wilcox_sig_final_C60_R60_flag_leaf_EGF_combined, 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 3)), 
          rowsep = c(seq(0, 47)), 
          sepwidth = c(0.01, 0.01), 
          hclustfun = function(x) hclust(x, method = "average"))
dev.off()


#venn diagram - common among cultivars

#N22
wilcox_sig_N22_C60_R60_flag_leaf_EGF_combined = wilcox_res_N22_C60_R60_flag_leaf_EGF_combined[-which(wilcox_res_N22_C60_R60_flag_leaf_EGF_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_N22_C60_R60_flag_leaf_EGF_combined = as.data.frame(wilcox_sig_log2fc_C60_R60_flag_leaf_EGF_combined[match(rownames(wilcox_sig_N22_C60_R60_flag_leaf_EGF_combined), rownames(wilcox_sig_log2fc_C60_R60_flag_leaf_EGF_combined)), 1])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_N22_C60_R60_flag_leaf_EGF_combined) = colnames(wilcox_sig_log2fc_C60_R60_flag_leaf_EGF_combined)[1]
wilcox_sig_log2fc_N22_C60_R60_flag_leaf_EGF_combined$Up.Down = ifelse(wilcox_sig_log2fc_N22_C60_R60_flag_leaf_EGF_combined$`GNR60 - GNC60` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_N22_C60_R60_flag_leaf_EGF_combined$Up.Down == "Up")    #13 metabolites
sum(wilcox_sig_log2fc_N22_C60_R60_flag_leaf_EGF_combined$Up.Down == "Down")  #9 metabolites

#Dular
wilcox_sig_Dular_C60_R60_flag_leaf_EGF_combined = wilcox_res_Dular_C60_R60_flag_leaf_EGF_combined[-which(wilcox_res_Dular_C60_R60_flag_leaf_EGF_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Dular_C60_R60_flag_leaf_EGF_combined = as.data.frame(wilcox_sig_log2fc_C60_R60_flag_leaf_EGF_combined[match(rownames(wilcox_sig_Dular_C60_R60_flag_leaf_EGF_combined), rownames(wilcox_sig_log2fc_C60_R60_flag_leaf_EGF_combined)), 2])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Dular_C60_R60_flag_leaf_EGF_combined) = colnames(wilcox_sig_log2fc_C60_R60_flag_leaf_EGF_combined)[2]
wilcox_sig_log2fc_Dular_C60_R60_flag_leaf_EGF_combined$Up.Down = ifelse(wilcox_sig_log2fc_Dular_C60_R60_flag_leaf_EGF_combined$`GDR60 - GDC60` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Dular_C60_R60_flag_leaf_EGF_combined$Up.Down == "Up")    #14 metabolites
sum(wilcox_sig_log2fc_Dular_C60_R60_flag_leaf_EGF_combined$Up.Down == "Down")  #13 metabolites

#Anjali
wilcox_sig_Anjali_C60_R60_flag_leaf_EGF_combined = wilcox_res_Anjali_C60_R60_flag_leaf_EGF_combined[-which(wilcox_res_Anjali_C60_R60_flag_leaf_EGF_combined$Significance == "ns"), ]    #significantly changed metabolites
wilcox_sig_log2fc_Anjali_C60_R60_flag_leaf_EGF_combined = as.data.frame(wilcox_sig_log2fc_C60_R60_flag_leaf_EGF_combined[match(rownames(wilcox_sig_Anjali_C60_R60_flag_leaf_EGF_combined), rownames(wilcox_sig_log2fc_C60_R60_flag_leaf_EGF_combined)), 3])   #log2fc values of significantly changed metabolites
colnames(wilcox_sig_log2fc_Anjali_C60_R60_flag_leaf_EGF_combined) = colnames(wilcox_sig_log2fc_C60_R60_flag_leaf_EGF_combined)[3]
wilcox_sig_log2fc_Anjali_C60_R60_flag_leaf_EGF_combined$Up.Down = ifelse(wilcox_sig_log2fc_Anjali_C60_R60_flag_leaf_EGF_combined$`GAR60 - GAC60` > 0, yes = "Up", no = "Down") #categorize if FC is up/down
sum(wilcox_sig_log2fc_Anjali_C60_R60_flag_leaf_EGF_combined$Up.Down == "Up")    #26 metabolites
sum(wilcox_sig_log2fc_Anjali_C60_R60_flag_leaf_EGF_combined$Up.Down == "Down")  #9 metabolites


#venn - increased - 60h rewatering/control 60h
vennlist_up_wilcox_C60_R60_flag_leaf_EGF_combined = list(rownames(wilcox_sig_log2fc_N22_C60_R60_flag_leaf_EGF_combined)[wilcox_sig_log2fc_N22_C60_R60_flag_leaf_EGF_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Dular_C60_R60_flag_leaf_EGF_combined)[wilcox_sig_log2fc_Dular_C60_R60_flag_leaf_EGF_combined$Up.Down == "Up"], rownames(wilcox_sig_log2fc_Anjali_C60_R60_flag_leaf_EGF_combined)[wilcox_sig_log2fc_Anjali_C60_R60_flag_leaf_EGF_combined$Up.Down == "Up"])
names(vennlist_up_wilcox_C60_R60_flag_leaf_EGF_combined) = c("N22", "Dular", "Anjali")
venn_up_wilcox_C60_R60_flag_leaf_EGF_combined = venn(vennlist_up_wilcox_C60_R60_flag_leaf_EGF_combined)    #plots venn diagram and lists results
grid.newpage()
venn_up_R60_C60_flag_leaf_EGF_combined = venn.diagram(vennlist_up_wilcox_C60_R60_flag_leaf_EGF_combined, category.names = c(paste(names(vennlist_up_wilcox_C60_R60_flag_leaf_EGF_combined[1]), " (", length(vennlist_up_wilcox_C60_R60_flag_leaf_EGF_combined$N22), ")", sep = ""), paste(names(vennlist_up_wilcox_C60_R60_flag_leaf_EGF_combined[2]), " (", length(vennlist_up_wilcox_C60_R60_flag_leaf_EGF_combined$Dular), ")", sep = ""), paste(names(vennlist_up_wilcox_C60_R60_flag_leaf_EGF_combined[3]), " (", length(vennlist_up_wilcox_C60_R60_flag_leaf_EGF_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "E", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label E
intersect_wilcox_C60_R60_venn_up_flag_leaf_EGF_combined = attr(venn_up_wilcox_C60_R60_flag_leaf_EGF_combined, "intersections")  #list of unique and common metabolites
venn_up_wilcox_metabolites_C60_R60_flag_leaf_EGF_combined = t(ldply(intersect_wilcox_C60_R60_venn_up_flag_leaf_EGF_combined, rbind))  #convert from list to matrix
colnames(venn_up_wilcox_metabolites_C60_R60_flag_leaf_EGF_combined) = venn_up_wilcox_metabolites_C60_R60_flag_leaf_EGF_combined[1, ]
venn_up_wilcox_metabolites_C60_R60_flag_leaf_EGF_combined = venn_up_wilcox_metabolites_C60_R60_flag_leaf_EGF_combined[-1, ]


#venn - decreased - 60h rewatering/control 60h
vennlist_down_wilcox_C60_R60_flag_leaf_EGF_combined = list(rownames(wilcox_sig_log2fc_N22_C60_R60_flag_leaf_EGF_combined)[wilcox_sig_log2fc_N22_C60_R60_flag_leaf_EGF_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Dular_C60_R60_flag_leaf_EGF_combined)[wilcox_sig_log2fc_Dular_C60_R60_flag_leaf_EGF_combined$Up.Down == "Down"], rownames(wilcox_sig_log2fc_Anjali_C60_R60_flag_leaf_EGF_combined)[wilcox_sig_log2fc_Anjali_C60_R60_flag_leaf_EGF_combined$Up.Down == "Down"])
names(vennlist_down_wilcox_C60_R60_flag_leaf_EGF_combined) = c("N22", "Dular", "Anjali")
venn_down_wilcox_C60_R60_flag_leaf_EGF_combined = venn(vennlist_down_wilcox_C60_R60_flag_leaf_EGF_combined)    #plots venn diagram and lists results
grid.newpage()
venn_down_R60_C60_flag_leaf_EGF_combined = venn.diagram(vennlist_down_wilcox_C60_R60_flag_leaf_EGF_combined, category.names = c(paste(names(vennlist_down_wilcox_C60_R60_flag_leaf_EGF_combined[1]), " (", length(vennlist_down_wilcox_C60_R60_flag_leaf_EGF_combined$N22), ")", sep = ""), paste(names(vennlist_down_wilcox_C60_R60_flag_leaf_EGF_combined[2]), " (", length(vennlist_down_wilcox_C60_R60_flag_leaf_EGF_combined$Dular), ")", sep = ""), paste(names(vennlist_down_wilcox_C60_R60_flag_leaf_EGF_combined[3]), " (", length(vennlist_down_wilcox_C60_R60_flag_leaf_EGF_combined$Anjali), ")", sep = "")), filename = NULL, fill = c ("orchid", "yellow", "light blue"), cex = 4.5, fontfamily = "Arial", cat.cex = 2.5, cat.default.pos = c("text"), cat.pos = c(12, -14, 0), cat.dist = c(0.1, 0.1, -0.07), cat.fontfamily = "Arial", cat.fontface = "bold", euler.d = FALSE, scaled = FALSE, main = "F", main.pos = c(0.05, 0.95), main.fontface = "bold", main.fontfamily = "Arial", main.cex = 3)   #with total number of metabolites per cultivar
#fig label F
intersect_wilcox_C60_R60_venn_down_flag_leaf_EGF_combined = attr(venn_down_wilcox_C60_R60_flag_leaf_EGF_combined, "intersections")  #list of unique and common metabolites
venn_down_wilcox_metabolites_C60_R60_flag_leaf_EGF_combined = t(ldply(intersect_wilcox_C60_R60_venn_down_flag_leaf_EGF_combined, rbind))  #convert from list to matrix
colnames(venn_down_wilcox_metabolites_C60_R60_flag_leaf_EGF_combined) = venn_down_wilcox_metabolites_C60_R60_flag_leaf_EGF_combined[1, ]
venn_down_wilcox_metabolites_C60_R60_flag_leaf_EGF_combined = venn_down_wilcox_metabolites_C60_R60_flag_leaf_EGF_combined[-1, ]#


write.table(venn_up_wilcox_metabolites_C60_R60_flag_leaf_EGF_combined, file = "wilcox_sig_up_metabolites_flag_leaf_EGF_R60_C60_combined.txt", sep = "\t", quote = F)
write.table(venn_down_wilcox_metabolites_C60_R60_flag_leaf_EGF_combined, file = "wilcox_sig_down_metabolites_flag_leaf_EGF_R60_C60_combined.txt", sep = "\t", quote = F)


####################


#rownames were checked to be the same for the data files to be combined
#significance of wilcoxon test - all cultivars, all three RW timepoints - including LS/CLS (refer to previous paper and associated GitHub files)
wilcox_sig_C_LS_R12_R36_R60_flag_leaf_EGF_combined = cbind(wilcox_res_N22_CLS_LS_flag_leaf_EGF_combined,
                                                           wilcox_res_N22_C12_R12_flag_leaf_EGF_combined$Significance,
                                                           wilcox_res_N22_C36_R36_flag_leaf_EGF_combined$Significance,
                                                           wilcox_res_N22_C60_R60_flag_leaf_EGF_combined$Significance,
                                                           wilcox_res_Dular_CLS_LS_flag_leaf_EGF_combined$Significance,
                                                           wilcox_res_Dular_C12_R12_flag_leaf_EGF_combined$Significance,
                                                           wilcox_res_Dular_C36_R36_flag_leaf_EGF_combined$Significance,
                                                           wilcox_res_Dular_C60_R60_flag_leaf_EGF_combined$Significance,
                                                           wilcox_res_Anjali_CLS_LS_flag_leaf_EGF_combined$Significance,
                                                           wilcox_res_Anjali_C12_R12_flag_leaf_EGF_combined$Significance,
                                                           wilcox_res_Anjali_C36_R36_flag_leaf_EGF_combined$Significance,
                                                           wilcox_res_Anjali_C60_R60_flag_leaf_EGF_combined$Significance)
wilcox_sig_C_LS_R12_R36_R60_flag_leaf_EGF_combined = wilcox_sig_C_LS_R12_R36_R60_flag_leaf_EGF_combined[ , -1]
colnames(wilcox_sig_C_LS_R12_R36_R60_flag_leaf_EGF_combined)[1] = "wilcox_res_N22_CLS_LS_flag_leaf_EGF_combined$Significance"
wilcox_sig_C_LS_R12_R36_R60_flag_leaf_EGF_combined = data.frame(lapply(wilcox_sig_C_LS_R12_R36_R60_flag_leaf_EGF_combined, function(x) {gsub("ns", " ", x)}))
rownames(wilcox_sig_C_LS_R12_R36_R60_flag_leaf_EGF_combined) = rownames(wilcox_res_N22_CLS_LS_flag_leaf_EGF_combined)
wilcox_sig_C_LS_R12_R36_R60_flag_leaf_EGF_combined$Nonsig = rowSums(wilcox_sig_C_LS_R12_R36_R60_flag_leaf_EGF_combined == " ")  #count number of nonsignificant
wilcox_sig_final_C_LS_R12_R36_R60_flag_leaf_EGF_combined = wilcox_sig_C_LS_R12_R36_R60_flag_leaf_EGF_combined[-which(wilcox_sig_C_LS_R12_R36_R60_flag_leaf_EGF_combined$Nonsig == "12"), ] #metabolites significant in at least one of the comparisons - 66 metabolites

#log2fc values of late stress and recovery relative to respective controls (with respect to stress paper and associated GitHub files)
rownames(log2fc_combined_control_latestress_EGF_flag_leaf) == rownames(log2fc_combined_control_recovery_EGF_flag_leaf)  #double-check if rownames are the same
log2fc_combined_controls_latestress_recovery_EGF_flag_leaf = cbind(log2fc_combined_control_latestress_EGF_flag_leaf, 
                                                                   log2fc_combined_control_recovery_EGF_flag_leaf)   #combine in one data set
log2fc_combined_controls_latestress_recovery_EGF_flag_leaf = log2fc_combined_controls_latestress_recovery_EGF_flag_leaf[ , c(1, 4:6, 2, 7:9, 3, 10:12)]   #rearrange in preferred order


#log2fc values of metabolites with significant difference in at least one of the comparisons - late stress, 12, 36, 60h recovery vs. corresponding controls in any cultivar
wilcox_sig_log2fc_C_LS_R12_R36_R60_flag_leaf_EGF_combined = log2fc_combined_controls_latestress_recovery_EGF_flag_leaf[which(rownames(log2fc_combined_controls_latestress_recovery_EGF_flag_leaf) %in% rownames(wilcox_sig_final_C_LS_R12_R36_R60_flag_leaf_EGF_combined)), ]


#rewatering-responsive metabolites only - indicate as red font in heat map
setdiff(rownames(wilcox_sig_log2fc_C_LS_R12_R36_R60_flag_leaf_EGF_combined), rownames(wilcox_sig_log2fc_CLS_LS_flag_leaf_EGF_combined))

#visualization
png("wilcox_sig_log2fc_C_LS_RW_flag_leaf_EGF_combined.png", width = 8*300, height = 8*300, res = 300)
heatmap.2(wilcox_sig_log2fc_C_LS_R12_R36_R60_flag_leaf_EGF_combined, 
          Rowv = TRUE, 
          Colv = FALSE, 
          dendrogram = "row", 
          trace = "none", 
          col = bluered(599), 
          breaks = c(seq(-3, -0.01, length.out = 300), seq(0.01, 3, length.out = 300)), 
          margins = c(7, 15),
          lmat = rbind(c(0, 4, 5), c(0, 1, 0), c(3, 2, 0)),
          labCol = rep(c("Severe stress/Control", "12 h RW/Control", "36 h RW/Control", "60 h RW/Control"), 3),
          cexCol = 1, 
          cexRow = 1, 
          density.info = "none", 
          key.xlab = NA, 
          key.title = NA, 
          lhei = c(1, 0.2, 11), 
          lwid = c(1.1, 5, 2.2),
          srtCol = 45, 
          cellnote = wilcox_sig_final_C_LS_R12_R36_R60_flag_leaf_EGF_combined[ , c(1:12)], 
          notecol = "black", 
          notecex = 1.5, 
          sepcolor = "black", 
          colsep = c(seq(0, 12)), 
          rowsep = c(seq(0, 66)), 
          sepwidth = c(0.01, 0.01), 
          ColSideColors = c(rep("#66c2a5", 4), rep("#fc8d62", 4), rep("#8da0cb", 4)),
          hclustfun = function(x) hclust(x, method = "average"),
          colRow = c(rep("black", 8), rep("red", 2), rep("black", 2), "red", rep("black", 8), "red", "black", "red", rep("black", 5),
                     "red", "black", rep("red", 2), rep("black", 3), "red", rep("black", 6), "red", rep("black", 2), "red", 
                     rep("black", 2), "red", rep("black", 3), rep("red", 2), rep("black", 8), "red", rep("black", 2)))
legend(x = 0.025, y = 1.1, xpd = TRUE, legend = "N22", bty = "n", cex = 1.2)
legend(x = 0.16, y = 1.1, xpd = TRUE, legend = "Dular", bty = "n", cex = 1.2)
legend(x = 0.3, y = 1.1, xpd = TRUE, legend = "Anjali", bty = "n", cex = 1.2)
legend(x = 0.83, y = 1.04, xpd = TRUE, legend = expression(paste("Log" ["2"], "-fold change")), bty = "n", cex = 0.8)
dev.off()


####################

#venn diagrams

#multiplot
png("venn_up_down_flag_leaf_EGF_RW_controls_combined.png", width = 16*300, height = 24*300, res = 300)
grid.arrange(gTree(children = venn_up_R12_C12_flag_leaf_EGF_combined), 
             gTree(children = venn_down_R12_C12_flag_leaf_EGF_combined), 
             gTree(children = venn_up_R36_C36_flag_leaf_EGF_combined),
             gTree(children = venn_down_R36_C36_flag_leaf_EGF_combined),
             gTree(children = venn_up_R60_C60_flag_leaf_EGF_combined),
             gTree(children = venn_down_R60_C60_flag_leaf_EGF_combined),
             ncol = 2, nrow = 3)
dev.off()


####################


#correlation analysis

#correlation between change in yield and metabolite changes during 60h RW
#correlation between change in proportion of chalky grains (>50% chalk content) and metabolite changes during 60h RW

#only 60h RW was considered since it is the latest sample collection time point


#FL stage

#yield and chalkiness - mean per cultivar per treatment per year already ran previously (refer to related GitHub files)

#log2fc - 60h RW / LS

#check first if rownames are the same before combining data
rownames(log2fc_stress_recovery_FL_flag_leaf_2013) == rownames(log2fc_stress_recovery_FL_flag_leaf_2014)
rownames(log2fc_stress_recovery_FL_flag_leaf_2013) == rownames(log2fc_stress_recovery_FL_flag_leaf_2015)

log2fc_flag_leaf_FL_R60_latestress_20131415 = cbind(log2fc_stress_recovery_FL_flag_leaf_2013[ , c(9, 6, 3)], 
                                                    log2fc_stress_recovery_FL_flag_leaf_2014[ , c(9, 6, 3)], 
                                                    log2fc_stress_recovery_FL_flag_leaf_2015[ , c(9, 6, 3)])    #log2fc of 60h RW / late stress - 3-yr data in one data table - order: Anjali, Dular, N22
colnames(log2fc_flag_leaf_FL_R60_latestress_20131415) = c("2013 Anjali", "2013 Dular", "2013 N22", "2014 Anjali", "2014 Dular", "2014 N22", "2015 Anjali", "2015 Dular", "2015 N22")
log2fc_flag_leaf_FL_R60_latestress_20131415 = log2fc_flag_leaf_FL_R60_latestress_20131415[ , match(colnames(yield_chalk_change_HxD_FL), colnames(log2fc_flag_leaf_FL_R60_latestress_20131415))]     #reorder columns to match column order of yield and chalkiness data


#yield, chalkiness, and log2FC in one data table

#check first if rownames are the same before combining data
colnames(yield_chalk_change_HxD_FL) == colnames(log2fc_flag_leaf_FL_R60_latestress_20131415)

yield_chalk_log2fc_R60_LS_flag_leaf_FL_HxD = rbind(yield_chalk_change_HxD_FL, log2fc_flag_leaf_FL_R60_latestress_20131415)
yield_chalk_log2fc_R60_LS_flag_leaf_FL_HxD = t(yield_chalk_log2fc_R60_LS_flag_leaf_FL_HxD)
yield_chalk_log2fc_R60_LS_flag_leaf_FL_HxD = as.data.frame(yield_chalk_log2fc_R60_LS_flag_leaf_FL_HxD)

#test for normality of data
shapiro_yield_chalk_log2fc_R60_LS_flag_leaf_FL_HxD = func_shapiro_test(yield_chalk_log2fc_R60_LS_flag_leaf_FL_HxD)
shapiro_yield_chalk_log2fc_R60_LS_flag_leaf_FL_HxD$Distribution = ifelse(shapiro_yield_chalk_log2fc_R60_LS_flag_leaf_FL_HxD$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_yield_chalk_log2fc_R60_LS_flag_leaf_FL_HxD$Distribution == "Normal")   #78 variables
sum(shapiro_yield_chalk_log2fc_R60_LS_flag_leaf_FL_HxD$Distribution == "Non-normal")   #6 variables


#correlation test - spearman

#correlation coefficient

#accounts for change in yield/chalkiness
cor_res_spearman_rho_R60_LS_flag_leaf_FL_HxD = cor(yield_chalk_log2fc_R60_LS_flag_leaf_FL_HxD, method = "spearman")


#correlation with change in yield and chalkiness

#yield
yield_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD = t(data.frame(lapply(yield_chalk_log2fc_R60_LS_flag_leaf_FL_HxD[4:84], function(x) cor.test(yield_chalk_log2fc_R60_LS_flag_leaf_FL_HxD$`Change in yield`, x, method = "spearman")$p.value)))    #significance of correlation
yield_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD = as.data.frame(yield_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD)
colnames(yield_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD) = "p.value"
rownames(yield_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD) = reduced_metabolite_list_final_flag_leaf_2013$Name
yield_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD$Significance = with(yield_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD,
                                                                        ifelse(p.value <= 0.001, "***", ifelse(
                                                                          p.value <= 0.01, "**", ifelse(
                                                                            p.value < 0.05, "*", "ns"))))
length(which(yield_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD$p.value < 0.05))   #3 metabolites
yield_cor_sig_metabolites_R60_LS_flag_leaf_FL_HxD = rownames(yield_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD)[yield_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD$Significance != "ns"]   #metabolites with significant correlation with change in yield
yield_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD = data.frame(cor_res_spearman_rho_R60_LS_flag_leaf_FL_HxD[1, which(colnames(cor_res_spearman_rho_R60_LS_flag_leaf_FL_HxD) %in% yield_cor_sig_metabolites_R60_LS_flag_leaf_FL_HxD)])     #metabolites with significant correlation with change in yield and rho values
yield_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD = round(yield_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD, 2)   #round to 2 decimals
yield_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD = as.matrix(yield_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD)  #convert to matrix
yield_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD = cbind(yield_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD, yield_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD)   #make 2-column data frame to be able to sort rownames alphabetically
yield_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD = yield_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD[order(rownames(yield_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD)), ]   #sort metabolites alphabetially
yield_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD = as.data.frame(yield_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD[ , -2])
colnames(yield_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD) = "rho value"
#export table
write.table(yield_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD, file = "cor_yield_sig_R60_LS_flag_leaf_FL_HxD.txt", sep = "\t", quote = F)


#chalky grains
chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD = t(data.frame(lapply(yield_chalk_log2fc_R60_LS_flag_leaf_FL_HxD[4:84], function(x) cor.test(yield_chalk_log2fc_R60_LS_flag_leaf_FL_HxD$`Change in chalky_50-75%`, x, method = "spearman")$p.value)))    #significance of correlation
chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD = as.data.frame(chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD)
colnames(chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD) = "p.value"
rownames(chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD) = reduced_metabolite_list_final_flag_leaf_2013$Name
chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD$Significance = with(chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD,
                                                                              ifelse(p.value <= 0.001, "***", ifelse(
                                                                                p.value <= 0.01, "**", ifelse(
                                                                                  p.value < 0.05, "*", "ns"))))
length(which(chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD$p.value < 0.05))   #3 metabolites
chalk_50_75_cor_sig_metabolites_R60_LS_flag_leaf_FL_HxD = rownames(chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD)[chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_FL_HxD$Significance != "ns"]   #metabolites with significant correlation with change in chalk_50_75
chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD = data.frame(cor_res_spearman_rho_R60_LS_flag_leaf_FL_HxD[3, which(colnames(cor_res_spearman_rho_R60_LS_flag_leaf_FL_HxD) %in% chalk_50_75_cor_sig_metabolites_R60_LS_flag_leaf_FL_HxD)])     #metabolites with significant correlation with change in chalk_50_75 and rho values
chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD = round(chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD, 2)   #round to 2 decimals
chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD = as.matrix(chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD)  #convert to matrix
chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD = cbind(chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD, chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD)   #make 2-column data frame to be able to sort rownames alphabetically
chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD = chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD[order(rownames(chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD)), ]   #sort metabolites alphabetially
chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD = as.data.frame(chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD[ , -2])
colnames(chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD) = "rho value"
#export table
write.table(chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_FL_HxD, file = "cor_chalk_sig_R60_LS_flag_leaf_FL_HxD.txt", sep = "\t", quote = F)


####################

#correlation analysis

#Early grain filling stage

#yield and chalkiness - mean per cultivar per treatment per year already ran previously

#R60 vs. late stress
stress_R60_EGF_order = c("GNLS", "GNR60", "GDLS", "GDR60", "GALS", "GAR60")

#2013
data_LS_R60_log2_mean_flag_leaf_EGF_2013 = data_log2_mean_EGF_flag_leaf_2013_2[stress_R60_EGF_order, ]  #data with stress and recovery60 values only
data_LS_R60_log2_mean_flag_leaf_EGF_2013 = t(data_LS_R60_log2_mean_flag_leaf_EGF_2013)

#log2-fold change - 60h RW / LS
log2fc_N22_stress_recovery_flag_leaf_EGF_2013 = data_LS_R60_log2_mean_flag_leaf_EGF_2013[ , 2] - data_LS_R60_log2_mean_flag_leaf_EGF_2013[ , 1]
log2fc_Dular_stress_recovery_flag_leaf_EGF_2013 = data_LS_R60_log2_mean_flag_leaf_EGF_2013[ , 4] - data_LS_R60_log2_mean_flag_leaf_EGF_2013[ , 3]
log2fc_Anjali_stress_recovery_flag_leaf_EGF_2013 = data_LS_R60_log2_mean_flag_leaf_EGF_2013[ , 6] - data_LS_R60_log2_mean_flag_leaf_EGF_2013[ , 5]
log2fc_stress_recovery_flag_leaf_EGF_2013 = as.data.frame(cbind(log2fc_N22_stress_recovery_flag_leaf_EGF_2013, 
                                                                log2fc_Dular_stress_recovery_flag_leaf_EGF_2013, 
                                                                log2fc_Anjali_stress_recovery_flag_leaf_EGF_2013))
colnames(log2fc_stress_recovery_flag_leaf_EGF_2013) = c("GNR60 - GNLS", "GDR60 - GDLS", "GAR60 - GALS")


#2014
data_LS_R60_log2_mean_flag_leaf_EGF_2014 = data_log2_mean_EGF_flag_leaf_2014_2[stress_R60_EGF_order, ]  #data with stress and recovery60 values only
data_LS_R60_log2_mean_flag_leaf_EGF_2014 = t(data_LS_R60_log2_mean_flag_leaf_EGF_2014)

#log2-fold change - 60h RW / LS
log2fc_N22_stress_recovery_flag_leaf_EGF_2014 = data_LS_R60_log2_mean_flag_leaf_EGF_2014[ , 2] - data_LS_R60_log2_mean_flag_leaf_EGF_2014[ , 1]
log2fc_Dular_stress_recovery_flag_leaf_EGF_2014 = data_LS_R60_log2_mean_flag_leaf_EGF_2014[ , 4] - data_LS_R60_log2_mean_flag_leaf_EGF_2014[ , 3]
log2fc_Anjali_stress_recovery_flag_leaf_EGF_2014 = data_LS_R60_log2_mean_flag_leaf_EGF_2014[ , 6] - data_LS_R60_log2_mean_flag_leaf_EGF_2014[ , 5]
log2fc_stress_recovery_flag_leaf_EGF_2014 = as.data.frame(cbind(log2fc_N22_stress_recovery_flag_leaf_EGF_2014, 
                                                                log2fc_Dular_stress_recovery_flag_leaf_EGF_2014, 
                                                                log2fc_Anjali_stress_recovery_flag_leaf_EGF_2014))
colnames(log2fc_stress_recovery_flag_leaf_EGF_2014) = c("GNR60 - GNLS", "GDR60 - GDLS", "GAR60 - GALS")


#2015
data_LS_R60_log2_mean_flag_leaf_EGF_2015 = data_log2_mean_EGF_flag_leaf_2015_2[stress_R60_EGF_order, ]  #data with stress and recovery60 values only
data_LS_R60_log2_mean_flag_leaf_EGF_2015 = t(data_LS_R60_log2_mean_flag_leaf_EGF_2015)

#log2-fold change - 60h RW / LS
log2fc_N22_stress_recovery_flag_leaf_EGF_2015 = data_LS_R60_log2_mean_flag_leaf_EGF_2015[ , 2] - data_LS_R60_log2_mean_flag_leaf_EGF_2015[ , 1]
log2fc_Dular_stress_recovery_flag_leaf_EGF_2015 = data_LS_R60_log2_mean_flag_leaf_EGF_2015[ , 4] - data_LS_R60_log2_mean_flag_leaf_EGF_2015[ , 3]
log2fc_Anjali_stress_recovery_flag_leaf_EGF_2015 = data_LS_R60_log2_mean_flag_leaf_EGF_2015[ , 6] - data_LS_R60_log2_mean_flag_leaf_EGF_2015[ , 5]
log2fc_stress_recovery_flag_leaf_EGF_2015 = as.data.frame(cbind(log2fc_N22_stress_recovery_flag_leaf_EGF_2015, 
                                                                log2fc_Dular_stress_recovery_flag_leaf_EGF_2015, 
                                                                log2fc_Anjali_stress_recovery_flag_leaf_EGF_2015))
colnames(log2fc_stress_recovery_flag_leaf_EGF_2015) = c("GNR60 - GNLS", "GDR60 - GDLS", "GAR60 - GALS")


#log2-fold change - 60h RW / LS

#check first if rownames are the same before combining data
rownames(log2fc_stress_recovery_flag_leaf_EGF_2013) == rownames(log2fc_stress_recovery_flag_leaf_EGF_2014)
rownames(log2fc_stress_recovery_flag_leaf_EGF_2013) == rownames(log2fc_stress_recovery_flag_leaf_EGF_2015)

log2fc_flag_leaf_EGF_R60_latestress_20131415 = cbind(log2fc_stress_recovery_flag_leaf_EGF_2013[ , ncol(log2fc_stress_recovery_flag_leaf_EGF_2013):1], log2fc_stress_recovery_flag_leaf_EGF_2014[ , ncol(log2fc_stress_recovery_flag_leaf_EGF_2014):1], log2fc_stress_recovery_flag_leaf_EGF_2015[ , ncol(log2fc_stress_recovery_flag_leaf_EGF_2015):1])    #log2fc of 60h RW / late stress - 3-yr data in one data table - order: Anjali, Dular, N22
colnames(log2fc_flag_leaf_EGF_R60_latestress_20131415) = c("2013 Anjali", "2013 Dular", "2013 N22", "2014 Anjali", "2014 Dular", "2014 N22", "2015 Anjali", "2015 Dular", "2015 N22")
log2fc_flag_leaf_EGF_R60_latestress_20131415 = log2fc_flag_leaf_EGF_R60_latestress_20131415[ , match(colnames(yield_chalk_change_HxD_EGF), colnames(log2fc_flag_leaf_EGF_R60_latestress_20131415))]     #reorder columns to match column order of yield and chalkiness data


#yield, chalkiness, and log2FC in one data table

#check first if rownames are the same before combining data
colnames(yield_chalk_change_HxD_EGF) == colnames(log2fc_flag_leaf_EGF_R60_latestress_20131415)

yield_chalk_log2fc_R60_LS_flag_leaf_EGF_HxD = rbind(yield_chalk_change_HxD_EGF, log2fc_flag_leaf_EGF_R60_latestress_20131415)
yield_chalk_log2fc_R60_LS_flag_leaf_EGF_HxD = t(yield_chalk_log2fc_R60_LS_flag_leaf_EGF_HxD)
yield_chalk_log2fc_R60_LS_flag_leaf_EGF_HxD = as.data.frame(yield_chalk_log2fc_R60_LS_flag_leaf_EGF_HxD)

#test for normality of data
shapiro_yield_chalk_log2fc_R60_LS_flag_leaf_EGF_HxD = func_shapiro_test(yield_chalk_log2fc_R60_LS_flag_leaf_EGF_HxD)
shapiro_yield_chalk_log2fc_R60_LS_flag_leaf_EGF_HxD$Distribution = ifelse(shapiro_yield_chalk_log2fc_R60_LS_flag_leaf_EGF_HxD$p.value < 0.05, yes = "Non-normal", no = "Normal")
sum(shapiro_yield_chalk_log2fc_R60_LS_flag_leaf_EGF_HxD$Distribution == "Normal")   #80 variables
sum(shapiro_yield_chalk_log2fc_R60_LS_flag_leaf_EGF_HxD$Distribution == "Non-normal")   #4 variables


#correlation test - spearman

#correlation coefficient

#accounts for change in yield/chalkiness
cor_res_spearman_rho_R60_LS_flag_leaf_EGF_HxD = cor(yield_chalk_log2fc_R60_LS_flag_leaf_EGF_HxD, method = "spearman")


#correlation with change in yield and chalkiness

#yield
yield_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD = t(data.frame(lapply(yield_chalk_log2fc_R60_LS_flag_leaf_EGF_HxD[4:84], function(x) cor.test(yield_chalk_log2fc_R60_LS_flag_leaf_EGF_HxD$`Change in yield`, x, method = "spearman")$p.value)))    #significance of correlation
yield_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD = as.data.frame(yield_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD)
colnames(yield_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD) = "p.value"
rownames(yield_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD) = reduced_metabolite_list_final_flag_leaf_2013$Name
yield_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD$Significance = with(yield_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD,
                                                                         ifelse(p.value <= 0.001, "***", ifelse(
                                                                           p.value <= 0.01, "**", ifelse(
                                                                             p.value < 0.05, "*", "ns"))))
length(which(yield_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD$p.value < 0.05))   #7 metabolites
yield_cor_sig_metabolites_R60_LS_flag_leaf_EGF_HxD = rownames(yield_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD)[yield_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD$Significance != "ns"]   #metabolites with significant correlation with change in yield
yield_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD = data.frame(cor_res_spearman_rho_R60_LS_flag_leaf_EGF_HxD[1, which(colnames(cor_res_spearman_rho_R60_LS_flag_leaf_EGF_HxD) %in% yield_cor_sig_metabolites_R60_LS_flag_leaf_EGF_HxD)])     #metabolites with significant correlation with change in yield and rho values
yield_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD = round(yield_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD, 2)   #round to 2 decimals
yield_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD = as.matrix(yield_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD)  #convert to matrix
yield_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD = cbind(yield_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD, yield_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD)   #make 2-column data frame to be able to sort rownames alphabetically
yield_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD = yield_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD[order(rownames(yield_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD)), ]   #sort metabolites alphabetially
yield_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD = as.data.frame(yield_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD[ , -2])
colnames(yield_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD) = "rho value"
#export table
write.table(yield_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD, file = "cor_yield_sig_R60_LS_flag_leaf_EGF_HxD.txt", sep = "\t", quote = F)


#chalky grains
chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD = t(data.frame(lapply(yield_chalk_log2fc_R60_LS_flag_leaf_EGF_HxD[4:84], function(x) cor.test(yield_chalk_log2fc_R60_LS_flag_leaf_EGF_HxD$`Change in chalky_50-75%`, x, method = "spearman")$p.value)))    #significance of correlation
chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD = as.data.frame(chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD)
colnames(chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD) = "p.value"
rownames(chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD) = reduced_metabolite_list_final_flag_leaf_2013$Name
chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD$Significance = with(chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD,
                                                                               ifelse(p.value <= 0.001, "***", ifelse(
                                                                                 p.value <= 0.01, "**", ifelse(
                                                                                   p.value < 0.05, "*", "ns"))))
length(which(chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD$p.value < 0.05))   #2 metabolites
chalk_50_75_cor_sig_metabolites_R60_LS_flag_leaf_EGF_HxD = rownames(chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD)[chalk_50_75_cor_res_spearman_pval_R60_LS_flag_leaf_EGF_HxD$Significance != "ns"]   #metabolites with significant correlation with change in chalk_50_75
chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD = data.frame(cor_res_spearman_rho_R60_LS_flag_leaf_EGF_HxD[3, which(colnames(cor_res_spearman_rho_R60_LS_flag_leaf_EGF_HxD) %in% chalk_50_75_cor_sig_metabolites_R60_LS_flag_leaf_EGF_HxD)])     #metabolites with significant correlation with change in chalk_50_75 and rho values
chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD = round(chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD, 2)   #round to 2 decimals
chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD = as.matrix(chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD)  #convert to matrix
chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD = cbind(chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD, chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD)   #make 2-column data frame to be able to sort rownames alphabetically
chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD = chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD[order(rownames(chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD)), ]   #sort metabolites alphabetially
chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD = as.data.frame(chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD[ , -2])
colnames(chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD) = "rho value"
#export table
write.table(chalk_50_75_sig_cor_rho_spearman_R60_LS_flag_leaf_EGF_HxD, file = "cor_chalk_sig_R60_LS_flag_leaf_EGF_HxD.txt", sep = "\t", quote = F)