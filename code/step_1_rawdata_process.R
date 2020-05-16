rm(list = ls())
library(Seurat)##version 3.1.1
library(dplyr)
library(future)
library(future.apply)
library(DoubletFinder)
plan("multiprocess", workers = 4) ###compute cores
options(future.globals.maxSize = 40000 * 1024^2)

getwd()
list.files(path = "1_raw_data",pattern = "^BC")
####BC2####
BC2<-Read10X(data.dir = "1_raw_data/BC2/")
bc2_object<- CreateSeuratObject(counts = BC2, project = "BC2", min.cells = 0, min.features = 200)
bc2_object[["percent.mt"]] <- PercentageFeatureSet(bc2_object, pattern = "^MT-")
hist(bc2_object[["percent.mt"]]$percent.mt)
VlnPlot(bc2_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(bc2_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bc2_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
CombinePlots(plots = list(plot1,plot2))
##select the cells with 300 genes at least and 4000 at most, the percent of mitochondrion genes is less than 10%
BC2_val<- subset(bc2_object, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 10)
colname<-paste("BC2_",colnames(BC2_val),sep="")
BC2_val<-RenameCells(object = BC2_val,colname)
rm(list = c("BC2","bc2_object"))

###Detect the doublets###
BC2_val <- NormalizeData(object = BC2_val)
BC2_val <- ScaleData(object = BC2_val)
BC2_val <- FindVariableFeatures(object = BC2_val)
BC2_val <- RunPCA(BC2_val, verbose = FALSE)
BC2_val <- RunTSNE(BC2_val,dims = 1:10, do.fast = T) 
BC2_val <- FindNeighbors(BC2_val,dims = 1:10)
BC2_val <- FindClusters(BC2_val, resolution = 0.6)

DimPlot(BC2_val,reduction = "tsne")
sweep.res.list_BC2_val <- paramSweep_v3(BC2_val, PCs = 1:10)
sweep.stats_BC2_val <- summarizeSweep(sweep.res.list_BC2_val, GT = FALSE)
bcmvn_BC2_val<-find.pK(sweep.stats_BC2_val)
pK_value <- as.numeric(as.character(bcmvn_BC2_val$pK[bcmvn_BC2_val$BCmetric == max(bcmvn_BC2_val$BCmetric)]))
annotations <- BC2_val@meta.data$RNA_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)  ####get the estimated number
nExp_poi <- round(0.104*length(BC2_val@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BC2_val <- doubletFinder_v3(BC2_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BC2_val <- doubletFinder(BC2_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)  
BC2_val@meta.data$Doublet <-BC2_val@meta.data$DF.classifications_0.25_0.005_545

BC2_val@meta.data[1:2,]
site<-rep("Insitu",nrow(BC2_val@meta.data))
BC2_val<-AddMetaData(object = BC2_val,metadata = site,col.name = "Site")

save(BC2_val,file = paste("2_output/BC2_db_remove",".RData",sep = ""))


####BC3####
BC3<-Read10X(data.dir = "1_raw_data/BC3")
bc3_object<- CreateSeuratObject(counts = BC3, project = "BC3", min.cells = 0, min.features = 200)
bc3_object[["percent.mt"]] <- PercentageFeatureSet(bc3_object, pattern = "^MT-")
hist(bc3_object[["percent.mt"]]$percent.mt)
VlnPlot(bc3_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(bc3_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bc3_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))
##select the cells with 200 genes at least and 3000 at most, the percent of mitochondrion genes is less than 10%
BC3_val<-subset(bc3_object, subset = nFeature_RNA > 300 & nFeature_RNA < 3000 & percent.mt < 10)
colname<-paste("BC3_",colnames(BC3_val),sep="")
BC3_val<-RenameCells(object = BC3_val,colname)
rm(list = c("BC3","bc3_object"))

###Detect the doublets###
BC3_val <- NormalizeData(object = BC3_val)
BC3_val <- ScaleData(object = BC3_val)
BC3_val <- FindVariableFeatures(object = BC3_val)
BC3_val <- RunPCA(BC3_val, verbose = FALSE)
BC3_val <- RunTSNE(BC3_val,dims = 1:10, do.fast = T) 
BC3_val <- FindNeighbors(BC3_val,dims = 1:10)
BC3_val <- FindClusters(BC3_val, resolution = 0.8)

DimPlot(BC3_val,reduction = "tsne")
sweep.res.list_BC3_val <- paramSweep_v3(BC3_val, PCs = 1:10)
sweep.stats_BC3_val <- summarizeSweep(sweep.res.list_BC3_val, GT = FALSE)
bcmvn_BC3_val<-find.pK(sweep.stats_BC3_val)
pK_value <- as.numeric(as.character(bcmvn_BC3_val$pK[bcmvn_BC3_val$BCmetric == max(bcmvn_BC3_val$BCmetric)]))
annotations <- BC3_val@meta.data$RNA_snn_res.0.8
homotypic.prop <- modelHomotypic(annotations)  ####get the estimated number
nExp_poi <- round(0.104*length(BC3_val@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BC3_val <- doubletFinder_v3(BC3_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BC3_val <- doubletFinder(BC3_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)  
BC3_val@meta.data$Doublet <-BC3_val@meta.data$DF.classifications_0.25_0.005_743

BC3_val@meta.data[1:2,]
site<-rep("Insitu",nrow(BC3_val@meta.data))
BC3_val<-AddMetaData(object = BC3_val,metadata = site,col.name = "Site")

save(BC3_val,file = paste("2_output/BC3_db_remove",".RData",sep = ""))

rm(list=ls())
####BC5####
BC5<-Read10X(data.dir = "1_raw_data/BC5")
bc5_object<- CreateSeuratObject(counts = BC5, project = "BC5", min.cells = 0, min.features = 200)
bc5_object[["percent.mt"]] <- PercentageFeatureSet(bc5_object, pattern = "^MT-")
hist(bc5_object[["percent.mt"]]$percent.mt)
VlnPlot(bc5_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(bc5_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bc5_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))
table((bc5_object[["percent.mt"]]$percent.mt)>10)
BC5_val<- subset(bc5_object, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 10)
colname<-paste("BC5_",colnames(BC5_val),sep="")
BC5_val<-RenameCells(object = BC5_val,colname)
rm(list = c("BC5","bc5_object"))

###Detect the doublets###
BC5_val <- NormalizeData(object = BC5_val)
BC5_val <- ScaleData(object = BC5_val)
BC5_val <- FindVariableFeatures(object = BC5_val)
BC5_val <- RunPCA(BC5_val, verbose = FALSE)
BC5_val <- RunTSNE(BC5_val,dims = 1:10, do.fast = T) 
BC5_val <- FindNeighbors(BC5_val,dims = 1:10)
BC5_val <- FindClusters(BC5_val, resolution = 0.6)

DimPlot(BC5_val,reduction = "tsne")
sweep.res.list_BC5_val <- paramSweep_v3(BC5_val, PCs = 1:10)
sweep.stats_BC5_val <- summarizeSweep(sweep.res.list_BC5_val, GT = FALSE)
bcmvn_BC5_val<-find.pK(sweep.stats_BC5_val)
pK_value <- as.numeric(as.character(bcmvn_BC5_val$pK[bcmvn_BC5_val$BCmetric == max(bcmvn_BC5_val$BCmetric)]))
annotations <- BC5_val@meta.data$RNA_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)  ####get the estimated number
nExp_poi <- round(0.0946*length(BC5_val@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BC5_val <- doubletFinder_v3(BC5_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BC5_val <- doubletFinder(BC5_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)  
BC5_val@meta.data$Doublet <-BC5_val@meta.data$DF.classifications_0.25_0.19_1481

BC5_val@meta.data[1:2,]
site<-rep("Insitu",nrow(BC5_val@meta.data))
BC5_val<-AddMetaData(object = BC5_val,metadata = site,col.name = "Site")

save(BC5_val,file = paste("2_output/BC5_db_remove",".RData",sep = ""))


rm(list=ls())

#####BC6###
BC6<-Read10X(data.dir = "1_raw_data/BC6")
bc6_object<- CreateSeuratObject(counts = BC6, project = "BC6", min.cells = 0, min.features = 200)
bc6_object[["percent.mt"]] <- PercentageFeatureSet(bc6_object, pattern = "^MT-")
hist(bc6_object[["percent.mt"]]$percent.mt)
VlnPlot(bc6_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(bc6_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bc6_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))
table((bc6_object[["percent.mt"]]$percent.mt)>10)
BC6_val<- subset(bc6_object, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 10)
colname<-paste("BC6_",colnames(BC6_val),sep="")
BC6_val<-RenameCells(object = BC6_val,colname)
rm(list = c("BC6","bc6_object"))

###Detect the doublets###
BC6_val <- NormalizeData(object = BC6_val)
BC6_val <- ScaleData(object = BC6_val)
BC6_val <- FindVariableFeatures(object = BC6_val)
BC6_val <- RunPCA(BC6_val, verbose = FALSE)
BC6_val <- RunTSNE(BC6_val,dims = 1:10, do.fast = T) 
BC6_val <- FindNeighbors(BC6_val,dims = 1:10)
BC6_val <- FindClusters(BC6_val, resolution = 0.6)

DimPlot(BC6_val,reduction = "tsne")
sweep.res.list_BC6_val <- paramSweep_v3(BC6_val, PCs = 1:10)
sweep.stats_BC6_val <- summarizeSweep(sweep.res.list_BC6_val, GT = FALSE)
bcmvn_BC6_val<-find.pK(sweep.stats_BC6_val)
pK_value <- as.numeric(as.character(bcmvn_BC6_val$pK[bcmvn_BC6_val$BCmetric == max(bcmvn_BC6_val$BCmetric)]))
annotations <- BC6_val@meta.data$RNA_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)  ####get the estimated number
nExp_poi <- round(0.0792*length(BC6_val@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BC6_val <- doubletFinder_v3(BC6_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BC6_val <- doubletFinder(BC6_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)  
BC6_val@meta.data$Doublet <-BC6_val@meta.data$DF.classifications_0.25_0.22_1447

BC6_val@meta.data[1:2,]
site<-rep("Insitu",nrow(BC6_val@meta.data))
BC6_val<-AddMetaData(object = BC6_val,metadata = site,col.name = "Site")

save(BC6_val,file = paste("2_output/BC6_db_remove",".RData",sep = ""))

rm(list=ls())


####BC10####
BC10<-Read10X(data.dir = "1_raw_data/BC10")
bc10_object<- CreateSeuratObject(counts = BC10, project = "BC10", min.cells = 0, min.features = 200)
bc10_object[["percent.mt"]] <- PercentageFeatureSet(bc10_object, pattern = "^MT-")
hist(bc10_object[["percent.mt"]]$percent.mt)
VlnPlot(bc10_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(bc10_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bc10_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))
BC10_val<- subset(bc10_object, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 10)
colname<-paste("BC10_",colnames(BC10_val),sep="")
BC10_val<-RenameCells(object = BC10_val,colname)
rm(list = c("BC10","bc10_object"))

###Detect the doublets###
BC10_val <- NormalizeData(object = BC10_val)
BC10_val <- ScaleData(object = BC10_val)
BC10_val <- FindVariableFeatures(object = BC10_val)
BC10_val <- RunPCA(BC10_val, verbose = FALSE)
BC10_val <- RunTSNE(BC10_val,dims = 1:10, do.fast = T) 
BC10_val <- FindNeighbors(BC10_val,dims = 1:10)
BC10_val <- FindClusters(BC10_val, resolution = 0.6)

DimPlot(BC10_val,reduction = "tsne")
sweep.res.list_BC10_val <- paramSweep_v3(BC10_val, PCs = 1:10)
sweep.stats_BC10_val <- summarizeSweep(sweep.res.list_BC10_val, GT = FALSE)
bcmvn_BC10_val<-find.pK(sweep.stats_BC10_val)
pK_value <- as.numeric(as.character(bcmvn_BC10_val$pK[bcmvn_BC10_val$BCmetric == max(bcmvn_BC10_val$BCmetric)]))
annotations <- BC10_val@meta.data$RNA_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)  ####get the estimated number
nExp_poi <- round(0.088*length(BC10_val@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BC10_val <- doubletFinder_v3(BC10_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BC10_val <- doubletFinder(BC10_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)  
BC10_val@meta.data$Doublet <-BC10_val@meta.data$DF.classifications_0.25_0.13_1257

BC10_val@meta.data[1:2,]
site<-rep("lungmeta",nrow(BC10_val@meta.data))
BC10_val<-AddMetaData(object = BC10_val,metadata = site,col.name = "Site")

save(BC10_val,file = paste("2_output/BC10_db_remove",".RData",sep = ""))

rm(list=ls())


#####BC11#####
BC11<-Read10X(data.dir = "1_raw_data/BC11/")
bc11_object<- CreateSeuratObject(counts = BC11, project = "BC11", min.cells = 0, min.features = 200)
bc11_object[["percent.mt"]] <- PercentageFeatureSet(bc11_object, pattern = "^MT-")
hist(bc11_object[["percent.mt"]]$percent.mt)
VlnPlot(bc11_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0,ncol = 3)
plot1 <- FeatureScatter(bc11_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bc11_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(bc11_object, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1,plot2,plot3))
BC11_val<- subset(bc11_object, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 10)
colname<-paste("BC11_",colnames(BC11_val),sep="")
BC11_val<-RenameCells(object = BC11_val,colname)
rm(list = c("BC11","bc11_object"))

###Detect the doublets###
BC11_val <- NormalizeData(object = BC11_val)
BC11_val <- ScaleData(object = BC11_val)
BC11_val <- FindVariableFeatures(object = BC11_val)
BC11_val <- RunPCA(BC11_val, verbose = FALSE)
BC11_val <- RunTSNE(BC11_val,dims = 1:10, do.fast = T) 
BC11_val <- FindNeighbors(BC11_val,dims = 1:10)
BC11_val <- FindClusters(BC11_val, resolution = 0.6)

DimPlot(BC11_val,reduction = "tsne")
sweep.res.list_BC11_val <- paramSweep_v3(BC11_val, PCs = 1:10)
sweep.stats_BC11_val <- summarizeSweep(sweep.res.list_BC11_val, GT = FALSE)
bcmvn_BC11_val<-find.pK(sweep.stats_BC11_val)
pK_value <- as.numeric(as.character(bcmvn_BC11_val$pK[bcmvn_BC11_val$BCmetric == max(bcmvn_BC11_val$BCmetric)]))
annotations <- BC11_val@meta.data$RNA_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)  ####get the estimated number
nExp_poi <- round(0.128*length(BC11_val@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BC11_val <- doubletFinder_v3(BC11_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BC11_val <- doubletFinder(BC11_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)  
BC11_val@meta.data$Doublet <-BC11_val@meta.data$DF.classifications_0.25_0.08_1078

BC11_val@meta.data[1:2,]
site<-rep("recurrence",nrow(BC11_val@meta.data))
BC11_val<-AddMetaData(object = BC11_val,metadata = site,col.name = "Site")

save(BC11_val,file = paste("2_output/BC11_db_remove",".RData",sep = ""))

rm(list=ls())



####BC16####
BC16<-Read10X(data.dir = "1_raw_data/BC16/")
bc16_object<- CreateSeuratObject(counts = BC16, project = "BC16", min.cells = 0, min.features = 200)
bc16_object[["percent.mt"]] <- PercentageFeatureSet(bc16_object, pattern = "^MT-")
hist(bc16_object[["percent.mt"]]$percent.mt)
VlnPlot(bc16_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(bc16_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bc16_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
CombinePlots(plots = list(plot1,plot2))
BC16_val<- subset(bc16_object, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & percent.mt < 10)
colname<-paste("BC16_",colnames(BC16_val),sep="")
BC16_val<-RenameCells(object = BC16_val,colname)
rm(list = c("BC16","bc16_object"))

###Detect the doublets###
BC16_val <- NormalizeData(object = BC16_val)
BC16_val <- ScaleData(object = BC16_val)
BC16_val <- FindVariableFeatures(object = BC16_val)
BC16_val <- RunPCA(BC16_val, verbose = FALSE)
BC16_val <- RunTSNE(BC16_val,dims = 1:10, do.fast = T) 
BC16_val <- FindNeighbors(BC16_val,dims = 1:10)
BC16_val <- FindClusters(BC16_val, resolution = 0.6)

DimPlot(BC16_val,reduction = "tsne")
sweep.res.list_BC16_val <- paramSweep_v3(BC16_val, PCs = 1:10)
sweep.stats_BC16_val <- summarizeSweep(sweep.res.list_BC16_val, GT = FALSE)
bcmvn_BC16_val<-find.pK(sweep.stats_BC16_val)
pK_value <- as.numeric(as.character(bcmvn_BC16_val$pK[bcmvn_BC16_val$BCmetric == max(bcmvn_BC16_val$BCmetric)]))
annotations <- BC16_val@meta.data$RNA_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)  ####get the estimated number
nExp_poi <- round(0.092*length(BC16_val@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BC16_val <- doubletFinder_v3(BC16_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BC16_val <- doubletFinder(BC16_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)  
BC16_val@meta.data$Doublet <-BC16_val@meta.data$DF.classifications_0.25_0.13_768

BC16_val@meta.data[1:2,]
site<-rep("Insitu",nrow(BC16_val@meta.data))
BC16_val<-AddMetaData(object = BC16_val,metadata = site,col.name = "Site")

save(BC16_val,file = paste("2_output/BC16_db_remove",".RData",sep = ""))

rm(list=ls())





####BC17####
BC17<-Read10X(data.dir = "1_raw_data/BC17/")
BC17_object<- CreateSeuratObject(counts= BC17, project = "BC17", min.cells = 0, min.features = 200)
BC17_object[["percent.mt"]] <- PercentageFeatureSet(BC17_object, pattern = "^MT-")
hist(BC17_object[["percent.mt"]]$percent.mt)
VlnPlot(BC17_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(BC17_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BC17_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))
BC17_val<- subset(BC17_object, subset = nFeature_RNA > 300 & nFeature_RNA < 7500 & percent.mt < 10)
colname<-paste("BC17_",colnames(BC17_val),sep="")
BC17_val<-RenameCells(object = BC17_val,colname)
rm(list = c("BC17","BC17_object"))

###Detect the doublets###
BC17_val <- NormalizeData(object = BC17_val)
BC17_val <- ScaleData(object = BC17_val)
BC17_val <- FindVariableFeatures(object = BC17_val)
BC17_val <- RunPCA(BC17_val, verbose = FALSE)
BC17_val <- RunTSNE(BC17_val,dims = 1:10, do.fast = T) 
BC17_val <- FindNeighbors(BC17_val,dims = 1:10)
BC17_val <- FindClusters(BC17_val, resolution = 0.6)

DimPlot(BC17_val,reduction = "tsne")
sweep.res.list_BC17_val <- paramSweep_v3(BC17_val, PCs = 1:10)
sweep.stats_BC17_val <- summarizeSweep(sweep.res.list_BC17_val, GT = FALSE)
bcmvn_BC17_val<-find.pK(sweep.stats_BC17_val)
pK_value <- as.numeric(as.character(bcmvn_BC17_val$pK[bcmvn_BC17_val$BCmetric == max(bcmvn_BC17_val$BCmetric)]))
annotations <- BC17_val@meta.data$RNA_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)  ####get the estimated number
nExp_poi <- round(0.137*length(BC17_val@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BC17_val <- doubletFinder_v3(BC17_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BC17_val <- doubletFinder(BC17_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)  
BC17_val@meta.data$Doublet <-BC17_val@meta.data$DF.classifications_0.25_0.04_447

BC17_val@meta.data[1:2,]
site<-rep("lungmeta",nrow(BC17_val@meta.data))
BC17_val<-AddMetaData(object = BC17_val,metadata = site,col.name = "Site")

save(BC17_val,file = paste("2_output/BC17_db_remove",".RData",sep = ""))


rm(list=ls())



####BC20####

BC20<-Read10X(data.dir = "1_raw_data/BC20/")
BC20_object<- CreateSeuratObject(counts= BC20, project = "BC20", min.cells = 0, min.features = 200)
BC20_object[["percent.mt"]] <- PercentageFeatureSet(BC20_object, pattern = "^MT-")
hist(BC20_object[["percent.mt"]]$percent.mt)
VlnPlot(BC20_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(BC20_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BC20_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))
BC20_val<- subset(BC20_object, subset = nFeature_RNA > 300 & nFeature_RNA <10000 & percent.mt < 10)
colname<-paste("BC20_",colnames(BC20_val),sep="")
BC20_val<-RenameCells(object = BC20_val,colname)
rm(list = c("BC20","BC20_object"))

###Detect the doublets###
BC20_val <- NormalizeData(object = BC20_val)
BC20_val <- ScaleData(object = BC20_val)
BC20_val <- FindVariableFeatures(object = BC20_val)
BC20_val <- RunPCA(BC20_val, verbose = FALSE)
BC20_val <- RunTSNE(BC20_val,dims = 1:10, do.fast = T) 
BC20_val <- FindNeighbors(BC20_val,dims = 1:10)
BC20_val <- FindClusters(BC20_val, resolution = 0.6)

DimPlot(BC20_val,reduction = "tsne")
sweep.res.list_BC20_val <- paramSweep_v3(BC20_val, PCs = 1:10)
sweep.stats_BC20_val <- summarizeSweep(sweep.res.list_BC20_val, GT = FALSE)
bcmvn_BC20_val<-find.pK(sweep.stats_BC20_val)
pK_value <- as.numeric(as.character(bcmvn_BC20_val$pK[bcmvn_BC20_val$BCmetric == max(bcmvn_BC20_val$BCmetric)]))
annotations <- BC20_val@meta.data$RNA_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)  ####get the estimated number
nExp_poi <- round(0.098*length(BC20_val@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BC20_val <- doubletFinder_v3(BC20_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BC20_val <- doubletFinder(BC20_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)  
BC20_val@meta.data$Doublet <-BC20_val@meta.data$DF.classifications_0.25_0.19_851

BC20_val@meta.data[1:2,]
site<-rep("recurrence",nrow(BC20_val@meta.data))
BC20_val<-AddMetaData(object = BC20_val,metadata = site,col.name = "Site")

save(BC20_val,file = paste("2_output/BC20_db_remove",".RData",sep = ""))


rm(list=ls())


####BC21####

BC21<-Read10X(data.dir = "1_raw_data/BC21/")
BC21_object<- CreateSeuratObject(counts= BC21, project = "BC21", min.cells = 0, min.features = 200)
BC21_object[["percent.mt"]] <- PercentageFeatureSet(BC21_object, pattern = "^MT-")
hist(BC21_object[["percent.mt"]]$percent.mt)
VlnPlot(BC21_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(BC21_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BC21_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))
BC21_val<- subset(BC21_object, subset = nFeature_RNA > 300 & nFeature_RNA < 7500 & percent.mt < 10)
colname<-paste("BC21_",colnames(BC21_val),sep="")
BC21_val<-RenameCells(object = BC21_val,colname)
rm(list = c("BC21","BC21_object"))

###Detect the doublets###
BC21_val <- NormalizeData(object = BC21_val)
BC21_val <- ScaleData(object = BC21_val)
BC21_val <- FindVariableFeatures(object = BC21_val)
BC21_val <- RunPCA(BC21_val, verbose = FALSE)
BC21_val <- RunTSNE(BC21_val,dims = 1:10, do.fast = T) 
BC21_val <- FindNeighbors(BC21_val,dims = 1:10)
BC21_val <- FindClusters(BC21_val, resolution = 0.6)

DimPlot(BC21_val,reduction = "tsne")
sweep.res.list_BC21_val <- paramSweep_v3(BC21_val, PCs = 1:10)
sweep.stats_BC21_val <- summarizeSweep(sweep.res.list_BC21_val, GT = FALSE)
bcmvn_BC21_val<-find.pK(sweep.stats_BC21_val)
pK_value <- as.numeric(as.character(bcmvn_BC21_val$pK[bcmvn_BC21_val$BCmetric == max(bcmvn_BC21_val$BCmetric)]))
annotations <- BC21_val@meta.data$RNA_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)  ####get the estimated number
nExp_poi <- round(0.117*length(BC21_val@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BC21_val <- doubletFinder_v3(BC21_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BC21_val <- doubletFinder(BC21_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)  
BC21_val@meta.data$Doublet <-BC21_val@meta.data$DF.classifications_0.25_0.17_453

BC21_val@meta.data[1:2,]
site<-rep("Insitu",nrow(BC21_val@meta.data))
BC21_val<-AddMetaData(object = BC21_val,metadata = site,col.name = "Site")

save(BC21_val,file = paste("2_output/BC21_db_remove",".RData",sep = ""))


rm(list=ls())






####BC22####

BC22<-Read10X(data.dir = "1_raw_data/BC22/")
BC22_object<- CreateSeuratObject(counts= BC22, project = "BC22", min.cells = 0, min.features = 200)
BC22_object[["percent.mt"]] <- PercentageFeatureSet(BC22_object, pattern = "^MT-")
hist(BC22_object[["percent.mt"]]$percent.mt)
VlnPlot(BC22_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(BC22_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BC22_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))
table((BC22_object[["percent.mt"]]$percent.mt)>10)
BC22_val<- subset(BC22_object, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & percent.mt < 10)
colname<-paste("BC22_",colnames(BC22_val),sep="")
BC22_val<-RenameCells(object = BC22_val,colname)
rm(list = c("BC22","BC22_object"))



###Detect the doublets###
BC22_val <- NormalizeData(object = BC22_val)
BC22_val <- ScaleData(object = BC22_val)
BC22_val <- FindVariableFeatures(object = BC22_val)
BC22_val <- RunPCA(BC22_val, verbose = FALSE)
BC22_val <- RunTSNE(BC22_val,dims = 1:10, do.fast = T) 
BC22_val <- FindNeighbors(BC22_val,dims = 1:10)
BC22_val <- FindClusters(BC22_val, resolution = 0.6)

DimPlot(BC22_val,reduction = "tsne")
sweep.res.list_BC22_val <- paramSweep_v3(BC22_val, PCs = 1:10)
sweep.stats_BC22_val <- summarizeSweep(sweep.res.list_BC22_val, GT = FALSE)
bcmvn_BC22_val<-find.pK(sweep.stats_BC22_val)
pK_value <- as.numeric(as.character(bcmvn_BC22_val$pK[bcmvn_BC22_val$BCmetric == max(bcmvn_BC22_val$BCmetric)]))
annotations <- BC22_val@meta.data$RNA_snn_res.0.6
homotypic.prop <- modelHomotypic(annotations)  ####get the estimated number
nExp_poi <- round(homotypic.prop*length(BC22_val@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BC22_val <- doubletFinder_v3(BC22_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BC22_val <- doubletFinder(BC22_val, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)  
BC22_val@meta.data$Doublet <-BC22_val@meta.data$DF.classifications_0.25_0.18_547

BC22_val@meta.data[1:2,]
site<-rep("Insitu",nrow(BC22_val@meta.data))
BC22_val<-AddMetaData(object = BC22_val,metadata = site,col.name = "Site")

save(BC22_val,file = paste("2_output/BC22_db_remove",".RData",sep = ""))

#rm(list=ls())



