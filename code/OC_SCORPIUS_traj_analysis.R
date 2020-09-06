rm(list=ls())
library(SCORPIUS)
library(Seurat)

library(dplyr)
library(future)
library(future.apply)
library(harmony)
library(monocle)
plan("multiprocess", workers = 4) 
options(future.globals.maxSize = 40000 * 1024^2) 

oct<-readRDS("OC_mature.rds")
srt <- NormalizeData(oct)
srt<-FindVariableFeatures(srt,selection.method = "vst", nfeatures = 2000)
vargenes<-srt@assays$RNA@var.features
expression <- t(as.matrix(srt@assays$RNA@data))
expression<-expression[,vargenes]
dim(expression)
group_name <- srt@meta.data$cell_group
group_name
####Reduce dimensionality of the dataset###
?reduce_dimensionality
set.seed(123)
space <- reduce_dimensionality(expression, dist="pearson", ndim = 3,num_landmarks = 500)
draw_trajectory_plot(space, progression_group = group_name, contour = TRUE)

traj <- infer_trajectory(space,k = 3)
traj$time
###draw the trajectory plot 
draw_trajectory_plot(
  space, 
  progression_group = group_name,
  path = traj$path,
  contour = TRUE
)

###Finding candidate marker genes###
?gene_importances
#We search for genes whose expression is seems to be a function of the trajectory timeline that was inferred, as such genes might be good candidate marker genes for dendritic cell maturation.
gimp <- gene_importances(expression, traj$time, num_permutations = 10, num_threads = 8)
gene_sel <- gimp[gimp$pvalue<0.05,]
expr_sel <- expression[,gene_sel$gene]
reverse_traj <- reverse_trajectory(traj)
draw_trajectory_heatmap(expr_sel, reverse_traj$time, group_name,show_labels_row = F)

###grouped into modules
?extract_modules
modules <- extract_modules(scale_quantile(expr_sel), reverse_traj$time, verbose = FALSE)
draw_trajectory_heatmap(expr_sel, reverse_traj$time, group_name, modules,show_labels_row = T)

write.csv(modules,file="OC_analysis_SCORPIUS/modules.csv")

###draw the top 50 genes###
gene_sel <- gimp[1:50,]
expr_sel <- expression[,gene_sel$gene]
reverse_traj <- reverse_trajectory(traj)
draw_trajectory_heatmap(expr_sel, reverse_traj$time, group_name,show_labels_row = T)

###grouped into modules
?extract_modules
modules <- extract_modules(scale_quantile(expr_sel), reverse_traj$time, verbose = FALSE)
draw_trajectory_heatmap(expr_sel, reverse_traj$time, group_name, modules,show_labels_row = T)

save.image(file="./OC_analysis_SCORPIUS/OC_SCORPIUS.RData")

