rm(list = ls())
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
plan("multiprocess", workers = 6) 
options(future.globals.maxSize = 40000 * 1024^2) 
library(devtools)
library(harmony)
library(tidyverse)
library(dyno)
library(SingleCellExperiment)
library(destiny)
#dynwrap::test_docker_installation(detailed = TRUE)
#data("fibroblast_reprogramming_treutlein")
getwd()
#setwd("")
sc20<-readRDS("~/Documents/hu_scRNA/repeated/BC20_OS.rds")
data.countbc20<-as.matrix(sc20@assays$RNA@counts)
sim <- SingleCellExperiment(assays = List(counts = data.countbc20))
#sim
#sim <- sim[tmp, ]
#dim(sim)
#log_mat <- log1p(GetAssayData(sc22, "counts"))
#so <- SetAssayData(sc22, "data", new.data = log_mat)
#sce <- as.SingleCellExperiment(so)
#dim(sim)
#geneFilter <- apply(assays(sim)$counts,1,function(x){
# sum(x >= 1) >= round(dim(sim)[2]/4)
})
#table(geneFilter)
sim <- sim[sc20@assays$RNA@var.features, ]
dim(sim)


pc20<-log1p(assays(sim)$counts)
pca <- prcomp(t(pc20), scale. = FALSE)
rd1 <- pca$x[,1:15]
plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)


library(destiny, quietly = TRUE)
dm <- DiffusionMap(t(pc20))
rd2 <- cbind(DC1 = dm$DC1, DC2 = dm$DC2)

plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)

plot(dm,1:2,col=sc20$harmony_cell_subgroup)


####performe the reduction with the diffusemap####


plot_dm_3D <- function(dm=dm, dc=c(1:3), condition=condition, colours=colours){
  cond <- factor(condition)
  col <- factor(condition)
  levels(col) <- colours
  col <- as.vector(col)
  DCs <- paste("DC",dc, sep="")
  
  data <- data.frame(
    dm@eigenvectors[,DCs[1]], 
    dm@eigenvectors[,DCs[2]], 
    dm@eigenvectors[,DCs[3]]
  )
  colnames(data) <- DCs
  
  plot3d(
    data,
    col=col,
    size=6.5,
    box = FALSE
  )
  
}

library(ggsci)
cors<-pal_d3()(10)
cors
#BiocManager::install("rgl")
#install.packages("rgl")

library(rgl)
plot_dm_3D(
  dm=dm, 
  dc=c(1:3),
  condition=sc20$harmony_cell_subgroup, 
  colour=cors
)


rgl.postscript("bc20_2.svg", fmt = "svg", drawText = TRUE )
plot_eigenVal <- function(dm=dm){
  linepad <- .5
  plot(
    eigenvalues(dm), 
    ylim = 0:1, 
    pch = 20, 
    xlab ='Diffusion component (DC)', 
    ylab ='Eigenvalue'
  )
}

plot_eigenVal(
  dm=dm
)

#table(cl1)
library(slingshot)

get_lineage <- function(dm=dm, 
                        dim=c(0:20), 
                        condition=condition, 
                        start=start, end=end, 
                        shrink.method="cosine"){
  
  data <- data.frame(
    dm@eigenvectors[,dim]
  )
  crv <- slingshot(
    data, 
    condition, 
    start.clus = start, 
    end.clus=end,
    maxit=100000,
    shrink.method=shrink.method,shrink=1
    # shrink.method="tricube"
  )
  
  return(crv)
}
?slingshot
#plot(dm,1:2,col=cl1)

sc20_lineage <- get_lineage(
  dm=dm, 
  dim=c(1:2), 
  condition=sc20$harmony_cell_subgroup,
  #shrink.method="tricube",
  start="Chondrocyte",
  end=c("Osteoblast_4"),
  shrink.method="cosine"
)




source("/Users/peizhanchen/Downloads/mouse-xx-xy-C1-e00dc2c8fe9e1bb28c34f6776979c882f0f8e81c/step5-Slingshot/functions.R")
sc20_pseudotime <- get_pseudotime(sc20_lineage, wthres=0.9)


summary(sc20_pseudotime)
dev.off()
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(10,'Spectral'))(100)
plotcol <- colors[cut(sc20_pseudotime, breaks=10)]

reducedDims(sim) <- SimpleList(PCA = rd1, DiffMap = rd2)
plot(reducedDims(sim)$DiffMap,  pch=16, asp = 0.5,col = sc20$harmony_cell_subgroup)

plot(reducedDims(sim)$DiffMap,  pch=16, asp = 0.5,col = plotcol)

plot(dm,1:2,col=sc20$harmony_cell_subgroup)
lines(SlingshotDataSet(sc20_lineage), lwd=1, col=c('white',"red","black"))






####identify the genes change along with the pseudotime
pseudotime_sc20<-sc20_pseudotime
pseudotime_lin <- pseudotime_sc20[,"curve1"]
max_pseudotime <- max(pseudotime_lin, na.rm = TRUE)
pseudotime_lin1_percent <- (pseudotime_lin*100)/max_pseudotime

pseudotime_lin <- pseudotime_sc20[,"curve2"]
max_pseudotime <- max(pseudotime_lin, na.rm = TRUE)
pseudotime_lin2_percent <- (pseudotime_lin*100)/max_pseudotime

pseudotime_sc20[,"curve1"] <- pseudotime_lin1_percent
pseudotime_sc20[,"curve2"] <- pseudotime_lin2_percent

pseudotime_sc20<-as.data.frame(pseudotime_sc20)
# select the ptime values 
ptime_1 <- pseudotime_sc20$curve1

# get cells in that lineage
lineage_cells_1 <- rownames(pseudotime_sc20)[!is.na(ptime_1)]
lineage_cells_1[1:5]
# remove values for cells not in the lineage
ptime_1 <- ptime_1[!is.na(ptime_1)]

# just test variable genes to save some time
genes_to_test <- VariableFeatures(sc20)[1:3000]

# get log normalized data to test
cnts <- pc20[genes_to_test, lineage_cells_1]

library(gam)
# fit a GAM with a loess term for pseudotime
gam.pval <- apply(cnts, 1, function(z){
  d <- data.frame(z = z, 
                  ptime = ptime_1)
  tmp <- suppressWarnings(gam(z ~ lo(ptime_1), data=d))
  p <- summary(tmp)[4][[1]][1, 5]
  p
})

# adjust pvalues 
res <- tibble(
  id = names(gam.pval),
  pvals = gam.pval,
  qval = p.adjust(gam.pval, method = "fdr")) %>% 
  arrange(qval)



library(ComplexHeatmap)
#BiocManager::install("ComplexHeatmap")

# get log normalized counts 
to_plot <- as.matrix(pc20[res$id[1:100], lineage_cells_1])

# arrange cells by pseudotime
ptime_order <- colnames(to_plot)[order(ptime_1)]
to_plot<-to_plot[,ptime_order]

cell_group<-as.vector(sc20@meta.data[ptime_order,]$harmony_cell_subgroup)
cell_group
time<-ptime_1[order(ptime_1)]

annotations<-cbind(cell_group,time)
annotations<-as.data.frame(annotations)
rownames(annotations)<-ptime_order
annotations$time<-as.numeric(annotations$time)
# add useful annotations
col_fun = colorRamp2(c(0, 50, 100), c("white", "grey", "black"))

ha <- HeatmapAnnotation(df = annotations,col = list(cell_group=c("Chondrocyte"="#1F77B4FF","Osteoblast_pro_1"="#FF7F0EFF","Osteoblast_pro_2"="#2CA02CFF",
                                                                 "Osteoblast_1"="#D62728FF","Osteoblast_2"="#9467BDFF","Osteoblast_3"="#8C564BFF","Osteoblast_4"="#E377C2FF"),
                                                    time=col_fun))
#ha <- HeatmapAnnotation(df = annotations)

Heatmap(to_plot,
        column_order = ptime_order,
        show_column_names = FALSE,
        show_row_names = T,
        row_names_gp = gpar(fontsize = 6),
        top_annotation = ha)

?Heatmap

######cell fate 2######
# select the ptime values 
ptime_2 <- pseudotime_sc20$curve2

# get cells in that lineage
lineage_cells_2 <- rownames(pseudotime_sc20)[!is.na(ptime_2)]
length(lineage_cells_2)
# remove values for cells not in the lineage
ptime_2 <- ptime_2[!is.na(ptime_2)]

# just test variable genes to save some time
genes_to_test <- VariableFeatures(sc20)[1:3000]

# get log normalized data to test
cnts <- pc20[genes_to_test, lineage_cells_2]

library(gam)
# fit a GAM with a loess term for pseudotime
gam.pval <- apply(cnts, 1, function(z){
  d <- data.frame(z = z, 
                  ptime = ptime_2)
  tmp <- suppressWarnings(gam(z ~ lo(ptime_2), data=d))
  p <- summary(tmp)[4][[1]][1, 5]
  p
})

# adjust pvalues 
res <- tibble(
  id = names(gam.pval),
  pvals = gam.pval,
  qval = p.adjust(gam.pval, method = "fdr")) %>% 
  arrange(qval)



library(ComplexHeatmap)
#BiocManager::install("ComplexHeatmap")

# get log normalized counts 
to_plot <- as.matrix(pc20[res$id[1:100], lineage_cells_2])

# arrange cells by pseudotime
ptime_order <- colnames(to_plot)[order(ptime_2)]
to_plot<-to_plot[,ptime_order]

cell_group<-as.vector(sc20@meta.data[ptime_order,]$harmony_cell_subgroup)
cell_group
time<-ptime_2[order(ptime_2)]

annotations<-cbind(cell_group,time)
annotations<-as.data.frame(annotations)
rownames(annotations)<-ptime_order
annotations$time<-as.numeric(annotations$time)
# add useful annotations
cors
table(sc20$harmony_cell_subgroup)
library(circlize)
col_fun = colorRamp2(c(0, 50, 100), c("white", "grey", "black"))

ha <- HeatmapAnnotation(df = annotations,col = list(cell_group=c("Chondrocyte"="#1F77B4FF","Osteoblast_pro_1"="#FF7F0EFF","Osteoblast_pro_2"="#2CA02CFF",
                                                                 "Osteoblast_1"="#D62728FF","Osteoblast_2"="#9467BDFF","Osteoblast_3"="#8C564BFF","Osteoblast_4"="#E377C2FF"),
                                                    time=col_fun))
ha
Heatmap(to_plot,
        column_order = ptime_order,
        show_column_names = FALSE,
        show_row_names = T,
        row_names_gp = gpar(fontsize = 6),
        top_annotation = ha)



getwd()
setwd("~/Documents/hu_scRNA/repeated/")
save.image(file="BC20_trajectory.RData")
load(file="BC20_trajectory.RData")
