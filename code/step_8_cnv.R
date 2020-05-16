rm(list = ls())
gc()
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
plan("multiprocess", workers = 8) 
options(future.globals.maxSize = 40000 * 1024^2) 
####test the harmony
###
library(devtools)
library(harmony)

tumor_chon<-read.csv("./chon_clean_cells_0425.csv",header = T,row.names = 1)
tumor_osteo<-read.csv("./Obsteoblast_clean_cells_0425.csv",header = T,row.names = 1)
tnk<-read.csv("./tandnk_0425.csv",header = T,row.names = 1)
myeloid<-read.csv("./myleoid_cell_0425.csv",header = T,row.names = 1)
tnkHarmony_group<-rep("tnk",dim(tnk)[1])
tnk$harmony_cell_subgroup<-tnkHarmony_group

colnames(tumor_chon)
chon_cnv<-tumor_chon[,c("orig.ident","Harmony_group","harmony_cell_subgroup")]
osteo_cnv<-tumor_osteo[,c("orig.ident","Harmony_group","harmony_cell_subgroup")]
tnk_cnv<-tnk[,c("orig.ident","Harmony_group","harmony_cell_subgroup")]
myeloid_cnv<-myeloid[,c("orig.ident","Harmony_group","harmony_cell_subgroup")]

cells_for_cnv<-rbind(chon_cnv,osteo_cnv,tnk_cnv,myeloid_cnv)
dim(cells_for_cnv)
cells_for_cnv[1:5,]
cells_for_cnv$cellname<-rownames(cells_for_cnv)
table(cells_for_cnv$harmony_cell_subgroup)
library(stringr)
celltype<-str_split(cells_for_cnv$harmony_cell_subgroup,pattern = "[_]",simplify = T)[,1]
table(celltype)
cells_for_cnv$cell_type<-celltype

cells_for_cnv[1:5,1:5]
write.csv(cells_for_cnv,file="cells_for_cnv.csv")


#sce_cnv<-sce_inter[,rownames(cells_for_cnv)]
#saveRDS(sce_cnv,file="sce_cnv.rds")
setwd("~/Documents/hu_scRNA/NC_Revised/osteosarcoma_project/4_CNV/")
sce_cnv<-readRDS(file="sce_cnv.rds")
table(rownames(sce_cnv@meta.data)==rownames(cells_for_cnv))
getwd()

#rm(BC20,BC202,blas,BC22,chon_clean,DC_cell,DC_Mat,em_lung_H,exp,exprMat,osb_clean,osb_row,osb_row_2,sce_inter)
gc()
library(infercnv)
#setwd("./4_CNV/")
getwd()
sce_cnv$cell_type<-cells_for_cnv$cell_type
sce_cnv@meta.data[1:5,]
table(sce_cnv@meta.data$nCount_RNA>2000,sce_cnv$orig.ident)

###select the cells with high quality####with UMI >2000
tmp<-sce_cnv@meta.data$nCount_RNA>2000
sce_cnv<-sce_cnv[,tmp]
sce_cnv
library(Seurat)

sce_split<-SplitObject(sce_cnv,split.by = "orig.ident")

subject<-names(sce_split)
subject
gene.order<-read.table("./gene_position_updated.txt")
gene.order<-gene.order[!duplicated(gene.order$V1),]
rownames(gene.order)<-gene.order$V1
gene.order<-gene.order[,-1]

for (i in 1:length(subject)){
  dir.create(path=subject[i])
  sce1<-sce_split[[i]]
  matrix<-as.matrix(sce1@assays$RNA@counts)
  ann<-as.data.frame(ifelse(sce1@meta.data[,"cell_type"]=="chon","maglignant_Chon",
                            ifelse(sce1@meta.data[,"cell_type"]=="Osteoblast","maglignant_osteoblast",
                                   ifelse(sce1@meta.data[,"cell_type"]=="myl","myeloid",
                                          ifelse(sce1@meta.data[,"cell_type"]=="tnk","tnk","")))),
                     row.names=rownames(sce1@meta.data))
  
  
  ann<-as.matrix(ann)
  colnames(ann)<-"cell_type"
  inferobj<-CreateInfercnvObject(raw_counts_matrix=matrix,
                                 annotations_file=ann,
                                 gene_order_file=gene.order,
                                 ref_group_names=c("myeloid","tnk"))
  outdir_1<-paste0("./",subject[i])
  
  inferobj <- infercnv::run(inferobj,
                            cutoff=0.1,
                            out_dir=outdir_1, 
                            cluster_by_groups=F, 
                            plot_steps = F,
                            no_prelim_plot = T,
                            denoise=TRUE,num_ref_groups = 2,
                            HMM=T,
                            num_threads=8,analysis_mode='subclusters',
                            png_res = 300)#, ,hclust_method='ward.D2'
}