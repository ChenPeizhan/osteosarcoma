####load the data####
#rm(list = ls())
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
plan("multiprocess", workers = 6) 
options(future.globals.maxSize = 40000 * 1024^2) 
####test the harmony
###
library(devtools)
library(harmony)
sce.big<-readRDS(file="2_output/merge_data_db_remove_harmony_group.rds")
table(sce.big@meta.data$Harmony_group)
tmp<-sce.big@meta.data$Harmony_group=="6.Myeloid"
table(tmp)
myl<-sce.big[,tmp]
rm(sce.big)

myl<-NormalizeData(myl,verbose = FALSE) 
myl<-FindVariableFeatures(myl,selection.method = "vst", nfeatures = 3000)
myl<-ScaleData(myl,vars.to.regress = "percent.mt",verbose = FALSE)
myl<-RunPCA(myl,verbose = T)
p1 <- DimPlot(object = myl, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = myl, features = "PC_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p1,p2))
myl<-RunHarmony(myl,"orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(myl, 'harmony')
dim(harmony_embeddings)
p3 <- DimPlot(object = myl, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p4 <- VlnPlot(object = myl, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p3,p4))


myl <- myl %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.6) %>% 
  identity()

myl<- FindClusters(myl,resolution = 1)
#rm(fibro)
DimPlot(myl, reduction = "tsne", pt.size = .1, split.by = 'Site')
DimPlot(myl, reduction = "umap", label = TRUE, pt.size = .1)

marker<-FindMarkers(myl,ident.1 = 8,only.pos = T)
FeaturePlot(myl,features = c("CD74","CD14","CD3D","MS4A1"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl,features = c("COL1A1","CD74","CD14","DCN","CTSK"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl,features = c("JCHAIN","CPA3","MS4A1"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl,features = c("CD3D","CD3E","NKG7","LUM"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl,features = c("COL1A1","LUM","IBSP","CTSK"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl,features = c("CD163","MRC1","TOP2A","HBB","GZMB"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl,features = c("IFIT1","IFNG","CCL2","CXCL1"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl,features = c("FCGR3A","S100A8","CD79A","MACRO"),reduction = "tsne",max.cutoff =3,label =T)
######remove the cells#####
tmp<-ifelse((Idents(myl)=="15"|Idents(myl)=="13"|Idents(myl)=="14"|Idents(myl)=="18"|Idents(myl)=="8"|Idents(myl)=="21"|Idents(myl)=="20"),TRUE,FALSE)
table(tmp)
myl_to_remove<-myl[,tmp]
table(Idents(myl_to_remove))
myl_to_remove<-RenameIdents(myl_to_remove,"15"="Osteoblast","13"="Osteoblast","14"="T/NK","18"="plasma","8"="doublets","21"="doublets","20"="lung_cells")
table(Idents(myl_to_remove))
myl_to_remove@meta.data[["harmony_cell_subgroup"]]<-paste("myl_to_remove",Idents(myl_to_remove),sep="_")
myl_to_remove@meta.data[1:5,]
write.csv(myl_to_remove@meta.data,file ="3_Figures_Tables/myl_to_remove0420.csv")



tmp<-ifelse((Idents(myl)=="15"|Idents(myl)=="13"|Idents(myl)=="14"|Idents(myl)=="18"|Idents(myl)=="8"|Idents(myl)=="21"|Idents(myl)=="20"),FALSE,TRUE)
table(tmp)
myl_2<-myl[,tmp]

myl_2<-NormalizeData(myl_2,verbose = FALSE) 
myl_2<-FindVariableFeatures(myl_2,selection.method = "vst", nfeatures = 3000)
myl_2<-ScaleData(myl_2,vars.to.regress = "percent.mt",verbose = FALSE)
myl_2<-RunPCA(myl_2,verbose = T)
p1 <- DimPlot(object = myl_2, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = myl_2, features = "PC_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p1,p2))
myl_2<-RunHarmony(myl_2,"orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(myl_2, 'harmony')
dim(harmony_embeddings)
p3 <- DimPlot(object = myl_2, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p4 <- VlnPlot(object = myl_2, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p3,p4))


myl_2 <- myl_2 %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.6) %>% 
  identity()

#myl_2<-FindClusters(myl_2,resolution = 1)
DimPlot(myl_2, reduction = "tsne", pt.size = .1, split.by = 'Site',label = T)
DimPlot(myl_2, reduction = "tsne", pt.size = .1, split.by = 'orig.ident',label = T)
DimPlot(myl_2, reduction = "umap", label = TRUE, pt.size = .1)
DimPlot(myl_2, reduction = "tsne", pt.size = .1,cols = cors)



FeaturePlot(myl_2,features = c("CD74","CD14","CD3D","MS4A1"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("CCL3","CXCL2","STAB1","S100A9"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("JCHAIN","CPA3","MS4A1","LTF"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("CD3D","CD3E","NKG7","LUM"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("CD74","CD14","CD3D","MS4A1"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("CXCL9","IL7R"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("FCER1A", "CD14","CD163","MRC1"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("MRC1","CD14","LYZ","CXCL8"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("CCL2","S100A8","S100A9"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("CSF1R","CD79A","IFIT1","IRF7"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("CLEC9A","PTCRA"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("DCLK1","COL1A1","IBSP"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("CD163","MRC1","SPP1","FCN1","FABP4"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("FCER1A","SELL","CCL2","CCL3"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("A2M","TGFB1","TREM2","LY6Z"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("CD14","FCGR3A","CDH1"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("FGF23","IBSP","SPP1","VEGFA"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("IL10","CCL19","ITGAX"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("S100A12","CSF3R","BST1"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("CD32","CD44","CD55"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("CD302","CD14","FCGR3A","IL2RA"),reduction = "umap",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("IL6","ATF3","MAF","FN1","TGFBI"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("CCL18","CCL20","BIRC3"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("TNF","CD68","CXCL2","PTPRC"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("CD14","FCGR3A","PTPRC"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(myl_2,features = c("CD14","CD1C","THBD","CLEC9A"),reduction = "tsne",max.cutoff =3,label =F)
FeaturePlot(myl_2,features = c("G0S2","CCR7","MRC1","CD163"),reduction = "tsne",max.cutoff =3,label =F)

markers<-FindAllMarkers(myl_2,only.pos = T,logfc.threshold = 0.1,min.pct = 0.1)

#markers_25<-FindAllMarkers(myl_2,only.pos = T,min.pct = 0.25)


myl_2<-RenameIdents(myl_2,"0"="CD14+_monocytes","7"="FCGR3A+_monocytes","6"="M1_like_TAM","1"="M2_like_TAM","2"="M2_like_TAM","5"="M2_like_TAM","9"="M2_like_TAM","10"="M2_like_TAM",
                    "8"="IFNg_activated_macrophage","13"="FABP4+_macrophage","3"="Neutrophil","4"="FCER1A+_DC",
                   "11"="CD141+/CLEC9A+_DC",
                    "12"="CCR7+_DC")

table(Idents(myl_2),myl_2$Site)

library(ggsci)
cors<-pal_d3()(10)
DimPlot(myl_2, reduction = "tsne", pt.size = .1, split.by = 'Site',cols = cors)
DimPlot(myl_2, reduction = "tsne", pt.size = .1,label = F,cols = cors)
DimPlot(myl_2, reduction = "tsne", pt.size = .1, split.by = 'orig.ident',ncol = 3,cols = cors)
DimPlot(myl_2, reduction = "umap", label = F, pt.size = .1,cols = cors)

VlnPlot(myl_2,features = c("CD74", "CD14","CD163","MRC1"),pt.size = 0,ncol = 2)
VlnPlot(myl_2,features = c("FCGR3A", "CCR7","CLEC9A","IFIT1"),pt.size = 0,ncol = 2)
VlnPlot(myl_2,features = c("CCL3", "FABP4","S100A8","G0S2"),pt.size = 0,ncol = 2)
table(Idents(myl_2))
markers<-FindAllMarkers(myl_2,only.pos = T,min.pct = 0.25)
markers_2<-FindAllMarkers(myl_2,only.pos = T,min.pct = 0.25,test.use = "roc")


top20<-markers %>% group_by(cluster) %>% top_n(n=20,wt=avg_logFC)
top10<-markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_logFC)
DoHeatmap(myl_2, features = top10$gene,group.colors = cors,label = F)
myl_2


markers_2<-markers[markers$p_val_adj<0.05,]
tab1<-table(Idents(myl_2),myl_2$tumor_type)
tab2<-prop.table(tab1,2)*100
tab2<-data.frame(tab2)
tab2<-tab2[20:1,]
library(ggplot2)
ggplot(tab2,aes(x=Var2,y=Freq,fill=Var1))+ 
  geom_bar(stat='identity',position='stack',alpha=.5)+ 
  labs(title='Bar with Stack',x='',y='')+ 
  theme(legend.justification = 'right', 
        legend.position = 'right', 
        legend.key.height = unit(0.1,'cm'),
        panel.background = element_blank(),
        )+ 
  scale_fill_manual(values=cors)


tab1<-table(Idents(myl_2),myl_2$Site)
tab2<-prop.table(tab1,2)*100
tab2<-data.frame(tab2)


ggplot(tab2,aes(x=Var2,y=Freq,fill=Var1))+ 
  geom_bar(stat='identity',position='stack',alpha=.5)+ 
  labs(title='Bar with Stack',x='',y='')+ 
  theme(legend.justification = 'right', 
        legend.position = 'right', 
        legend.key.height = unit(0.1,'cm'),
        panel.background = element_blank(),
  )+ 
  scale_fill_manual(values=cors)


######perform the heatmap of the differntialtion expression genes####
markers.to.plot <- c("CD74","CD14","FCGR3A","CD163","MRC1","IRF1","IDO1","MAF","FABP4","IFIT1","S100A8","CD1C","FCER1A","CLEC9A","CCR7")
DotPlot(myl_2, features = rev(markers.to.plot), cols = c("blue","#E64B35FF"), dot.scale = 6,scale = T,)

####GSEA analyssis###
#library(splatter)
library(msigdbr)
library(singleseqgset)
#install.packages("heatmap3")
library(heatmap3)
#install_github("arc85/singleseqgset")

h.human <- msigdbr(species="Homo sapiens",category="H")

h.names <- unique(h.human$gs_name)

h.sets <- vector("list",length=length(h.names))
names(h.sets) <- h.names

for (i in names(h.sets)) {
  h.sets[[i]] <- pull(h.human[h.human$gs_name==i,"gene_symbol"])
}


table(Idents(myl_2))
myl_2$real_id<-Idents(myl_2)
logfc.data <- logFC(cluster.ids=myl_2@meta.data$real_id,expr.mat=myl_2@assays$RNA@data)
names(logfc.data)
gse.res <- wmw_gsea(expr.mat=myl_2@assays$RNA@data,cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],gene.sets=h.sets)
res.stats <- gse.res[["GSEA_statistics"]]
res.stats
res.pvals <- gse.res[["GSEA_p_values"]]
res.pvals <- apply(res.pvals,2,p.adjust,method="fdr") #Correct for multiple comparisons

#res.stats[order(res.stats[,1],decreasing=TRUE)[1:10],] #Top gene sets enriched by z scores
ac<-data.frame(cluster=factor(colnames(res.stats)))
ac
rownames(ac)=colnames(res.stats)

#paste(rownames(ac),cors,sep="=")

ann_colors = list(cluster = c("CD14+_monocytes"="#1F77B4FF","FCGR3A+_monocytes"="#FF7F0EFF",
                              "M1_like_TAM"="#2CA02CFF","M2_like_TAM"="#D62728FF",
                              "IFNg_activated_macrophage"="#9467BDFF","FABP4+_macrophage"="#8C564BFF",
                              "Neutrophil"="#E377C2FF","FCER1A+_DC"="#7F7F7FFF",
                              "CD141+/CLEC9A+_DC"="#BCBD22FF","CCR7+_DC"="#17BECFFF"))

dev.off()
pheatmap::pheatmap(res.stats,fontsize_row = 8,annotation_col = ac,annotation_legend = F,annotation_colors = ann_colors,scale="row",cluster_rows = F,
                   cluster_cols = F)


######GSVA analysis of the #######
signature<-read.csv(file="~/Documents/public_data/Macrophage.csv",header=T,as.is = T)
gmt_in_use<-list()
gmt_in_use[["M1_actite"]]<-signature$M1_UP
gmt_in_use[["M2_actite"]]<-signature$M2_UP
gmt_in_use[["IFNg"]]<-signature$GSE1925_CTRL_VS_IFNG_PRIMED_MACROPHAGE_24H_IFNG_STIM_DN



table(Idents(myl_2))
macro_tmp<-Idents(myl_2)=="M2_like_TAM"|Idents(myl_2)=="M1_like_TAM"|Idents(myl_2)=="FABP4+_macrophage"|Idents(myl_2)=="IFNg_activated_macrophage"
table(macro_tmp)

macro<-myl_2[,macro_tmp]
exprMat<-GetAssayData(macro,slot="data")
exprMat<-as.matrix(exprMat)
res_es <- gsva(exprMat, gmt_in_use, min.sz=1, max.sz=500, verbose=FALSE, parallel.sz=8)
write.csv(res_es,file="three_subgroup_macrophae.csv")
table(Idents(macro))
cluster_infor<-sort(Idents(macro))
cluster_infor[1:5]
res_es2<-res_es[,names(cluster_infor)]
res_es2[,1:5]
res_es2<-as.data.frame(res_es2)
ac<-as.data.frame(cluster_infor)
table(rownames(ac)==colnames(res_es2))
levels(cluster_infor)
cors
table(Idents(myl_2))
top_anno <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#2CA02CFF","#D62728FF","#9467BDFF","#8C564BFF")), 
                       labels = levels(cluster_infor),
                       labels_gp = gpar(cex = 0, col = "white"))) 
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

res_es2<-as.matrix(res_es2)
?scale
res_es3<-scale(res_es2)
cluster_infor<-sort(Idents(macro))
Heatmap(res_es3,
        cluster_rows = FALSE,
        cluster_columns = T,
        show_column_names = FALSE,
        show_row_names = T,
        top_annotation = top_anno,
        column_split = cluster_infor)
pheatmap::pheatmap(res_es2,annotation_col = ac,annotation_legend = T,cluster_cols = T,annotation_names_col = T,cluster_rows = F,scale = "column")


#####下面对DC的数据进行分析###

DC_tmp<-Idents(myl_2)=="FCER1A+_DC"|Idents(myl_2)=="CD141+/CLEC9A+_DC"|Idents(myl_2)=="CCR7+_DC"
table(DC_tmp)

DC_cell<-myl_2[,DC_tmp]
DC_Mat<-GetAssayData(DC_cell,slot="data")
DC_Mat<-as.matrix(DC_Mat)

h.human <- msigdbr(species="Homo sapiens",category="C7")
h.names <- unique(h.human$gs_name)

h.sets <- vector("list",length=length(h.names))
names(h.sets) <- h.names
for (i in names(h.sets)) {
  h.sets[[i]] <- pull(h.human[h.human$gs_name==i,"gene_symbol"])
}


table(Idents(DC_cell))
DC_cell$real_id<-Idents(DC_cell)
logfc.data <- logFC(cluster.ids=DC_cell@meta.data$real_id,expr.mat=DC_cell@assays$RNA@data)
names(logfc.data)
gse.res <- wmw_gsea(expr.mat=DC_cell@assays$RNA@data,cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],gene.sets=h.sets_2)
res.stats <- gse.res[["GSEA_statistics"]]

matt2<-res.stats[order(res.stats$`CCR7+_DC`,decreasing = T),]
dim(matt2)
matt3<-matt2[1:10,]
matt4<-res.stats[order(res.stats$`CD141+/CLEC9A+_DC`,decreasing = T),]
dim(matt4)
matt4<-matt4[1:10,]
matt5<-res.stats[order(res.stats$`FCER1A+_DC`,decreasing = T),]
dim(matt5)
matt5<-matt5[1:10,]
names<-unique(c(rownames(matt3),rownames(matt4),rownames(matt5)))

mmmm<-res.stats[names,]

Heatmap(as.matrix(mmmm),
        cluster_rows = FALSE,
        cluster_columns = F,
        show_column_names = T,
        show_row_names = T,
        top_annotation = top_anno,row_names_max_width = unit(1, "cm"))



pheatmap::pheatmap(mmmm,fontsize_row = 8,annotation_legend = F,annotation_colors = ann_colors,cluster_rows = F,scale="row",
                   cluster_cols = F)





table(Idents(myl_2))
myl_2@meta.data[["harmony_cell_subgroup"]]<-paste("myl_to_keep",Idents(myl_2),sep="_")
myl_2@meta.data[1:5,]
myl_2
write.csv(myl_2@meta.data,file ="3_Figures_Tables/myl_2.csv")


#####myeloid cells distrubution###
table(sce.big$orig.ident)
table(myl_2$orig.ident,Idents(myl_2))

myl_pro<-data.frame(table(sce.big$orig.ident),table(myl_2$orig.ident))[,c(1,2,4)]
colnames(myl_pro)<-c("sampleID","total_cell","myeloid_cell")
myl_pro<-mutate(myl_pro,total_ratio=myl_pro$myeloid_cell/myl_pro$total_cell)
myl_pro$sampleID<-factor(myl_pro$sampleID,levels=c("BC2","BC3","BC5","BC6","BC21","BC16","BC22","BC10","BC17","BC11","BC20"))


specific<-cbind(table(sce.big$orig.ident),table(myl_2$orig.ident,Idents(myl_2)))

specific<-data.frame(specific)
specific$sampleID<-rownames(specific)
specific<-mutate(specific,neutrophil_ratio=specific$Neutrophil/specific$V1)
specific$sampleID<-factor(specific$sampleID,levels=c("BC2","BC3","BC5","BC6","BC21","BC16","BC22","BC10","BC17","BC11","BC20"))


library(ggplot2)
ggplot(data=myl_pro,aes(sampleID,total_ratio))+
  geom_bar(stat = "identity",color="white",fill="navy")+
  theme_light()


ggplot(data=specific,aes(sampleID,neutrophil_ratio))+
  geom_bar(stat = "identity",color="white",fill="navy")+
  theme_light()




save.image(file="./myeloid_0422.RData")
load(file="./myeloid_0422.RData")
DimPlot()
write.csv(myl_2@meta.data,file="myleoid_cell_0425.csv")

table(Idents(myl_2),myl_2$orig.ident)
prop.table(table(Idents(myl_2),myl_2$orig.ident),2)

