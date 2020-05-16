####load the data####
rm(list = ls())
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
plan("multiprocess", workers = 6) 
options(future.globals.maxSize = 40000 * 1024^2) 
library(devtools)
library(harmony)

DimPlot(sce_inter,reduction = "tsne")
table(Idents(sce_inter))
sce_inter@meta.data[1:5,]
sce_inter@meta.data$Harmony_group<-Idents(sce_inter)

sce.big<-readRDS(file="2_output/merge_data_db_remove_harmony_group.rds")
table(sce_inter@meta.data$Harmony_group)
tmp<-ifelse((sce_inter@meta.data$Harmony_group=="1.Osteoblast"|sce_inter@meta.data$Harmony_group=="2.Osteoblast_proli"|sce_inter@meta.data$Harmony_group=="3.Chondrocyte"),TRUE,FALSE)
#FeaturePlot(sce_inter,features = c("SOX9","ACAN","COL2A1"),reduction = "tsne",label = T)

table(tmp)
osb<-sce_inter[,tmp]


osb<-NormalizeData(osb,verbose = FALSE) 
osb<-FindVariableFeatures(osb,selection.method = "vst", nfeatures = 3000)
osb<-ScaleData(osb,vars.to.regress = "percent.mt",verbose = FALSE)
osb<-RunPCA(osb,verbose = T)
p1 <- DimPlot(object = osb, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = osb, features = "PC_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p1,p2))
osb<-RunHarmony(osb,"orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(osb, 'harmony')
dim(harmony_embeddings)
p3 <- DimPlot(object = osb, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p4 <- VlnPlot(object = osb, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p3,p4))


osb<- osb %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.6) %>% 
  identity()
#osb<-FindClusters(osb,resolution = 0.6)
DimPlot(osb, reduction = "tsne", pt.size = .1, split.by = 'Site')
DimPlot(osb, reduction = "umap", label = TRUE, pt.size = .1)
DimPlot(osb, reduction = "umap", label = TRUE, pt.size = .1,split.by = 'Site')
DimPlot(osb, reduction = "tsne", label = TRUE, pt.size = .1)

####chondrocytes
FeaturePlot(osb,features = c("SOX9","TOP2A","ACAN","COL2A1"),reduction = "tsne",max.cutoff =3,label =F)
####MSC
FeaturePlot(osb,features = c("CXCL12","LEPR","RUNX1"),reduction = "tsne",max.cutoff =3,label =F)
####OSteobalst
FeaturePlot(osb,features = c("RUNX2","COL1A1","TNC","SPP1"),reduction = "tsne",max.cutoff =3,label =F)
FeaturePlot(osb,features = c("EPCAM","PECAM1","HBB","JCHAIN"),reduction = "tsne",max.cutoff =3,label =T)
markers<-FindAllMarkers(osb,only.pos = T)
top10<-markers %>% group_by(cluster) %>% top_n(n=10,wt = avg_logFC)


####remove the doublets and red blood cells according to the markers
tmp<-ifelse((Idents(osb)==5|Idents(osb)==13|Idents(osb)==11|Idents(osb)==10|Idents(osb)==16|Idents(osb)==17),FALSE,TRUE)
table(tmp)
####analysis the clean osb####
osb_clean<-osb[,tmp]
osb_clean<-NormalizeData(osb_clean,verbose = FALSE) 
osb_clean<-FindVariableFeatures(osb_clean,selection.method = "vst", nfeatures = 3000)
osb_clean<-ScaleData(osb_clean,vars.to.regress = "percent.mt",verbose = FALSE)
osb_clean<-RunPCA(osb_clean,verbose = T,npcs = 30)
p1 <- DimPlot(object = osb_clean, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = osb_clean, features = "PC_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p1,p2))
osb_clean<-RunHarmony(osb_clean,"orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(osb_clean, 'harmony')
dim(harmony_embeddings)
p3 <- DimPlot(object = osb_clean, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p4 <- VlnPlot(object = osb_clean, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p3,p4))


osb_clean<- osb_clean %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) 

osb_clean<-FindClusters(osb_clean,resolution = 0.8)
Idents(osb_clean)<-osb_clean$RNA_snn_res.0.8
DimPlot(osb_clean, reduction = "tsne", pt.size = .1,label = T)
DimPlot(osb_clean, reduction = "tsne", label = TRUE, pt.size = .1,split.by = "Site")
DimPlot(osb_clean, reduction = "umap", label = TRUE, pt.size = .1,split.by = 'Site')


####chondrocytes#
FeaturePlot(osb_clean,features = c("SOX9","TOP2A","ACAN","COL2A1"),reduction = "tsne",max.cutoff =3,label =T)
####OSteobalst####
FeaturePlot(osb_clean,features = c("RUNX2","COL1A1","TNC","SPP1"),reduction = "umap",max.cutoff =3,label =F)
FeaturePlot(osb_clean,features = c("CXCL12","MMP13","SP7"),reduction = "umap",max.cutoff =3,label =F)
FeaturePlot(osb_clean,features = c("COL2A1","SOX9","ACAN"),reduction = "umap",max.cutoff =3,label =F)
FeaturePlot(osb_clean,features = c("LUM","COL3A1","SERPINF1","RUNX2"),reduction = "tsne",max.cutoff =3,label =F)


osb_clean<-RenameIdents(osb_clean,"10"="Chondrocyte","2"="Chondrocyte","15"="Osteoblast_pro_1","4"="Osteoblast_pro_1",
                        "3"="Osteoblast_pro_2","14"="Osteoblast_pro_2",
                        "1"="Osteoblast_1","8"="Osteoblast_1","19"="Osteoblast_1","22"="Osteoblast_1","17"="Osteoblast_1",
                         "6"="Osteoblast_2","9"="Osteoblast_2",
                        "7"="Osteoblast_3","18"="Osteoblast_4","5"="Osteoblast_4","13"="Osteoblast_4",
                        "20"="Osteoblast_4","0"="Osteoblast_4","12"="Osteoblast_4","21"="Osteoblast_4","16"="Osteoblast_4","11"="Osteoblast_4")


table(Idents(osb_clean))
table(osb_clean$harmony_cell_subgroup)
library(ggsci)
cors<-pal_d3()(10)
cors

DimPlot(osb_clean, reduction = "tsne", pt.size = .1, split.by = 'Site',cols = cors)
DimPlot(osb_clean, reduction = "tsne", pt.size = .1,label = F,cols = cors)
DimPlot(osb_clean, reduction = "umap", pt.size = .1,label = F,cols = cors)
DimPlot(osb_clean, reduction = "tsne", pt.size = .1, split.by = 'orig.ident',ncol = 3,cols = cors)

library(RColorBrewer)
FeaturePlot(osb_clean,features = c("COL1A1","RUNX2","TNC","COL2A1","ACAN","SOX9"),reduction = "tsne",max.cutoff =8,label =F,cols =  brewer.pal(9,"Reds"))
#FeaturePlot(osb_clean,features = c(),reduction = "tsne",max.cutoff =8,label =F,cols =  brewer.pal(9,"RdPu"))

display.brewer.all()
brewer.pal(9,"RdBu")[9:1]
markers_total_pos<-FindAllMarkers(osb_clean,only.pos = T,min.pct = 0.1,logfc.threshold = 0.1)
markers_total_total<-FindAllMarkers(osb_clean,only.pos = F,min.pct = 0.1,logfc.threshold = 0.1)




#####compare the genes between the groups####
####compare the osteoblast cells from different site of the cells####
table(Idents(osb_clean))
osb_clean[["harmony_cell_subgroup"]]<-Idents(osb_clean)
#rm(osb)
osb_clean@meta.data[1:5,]
tmp<-Idents(osb_clean)=="Chondrocyte"
tmp<-!tmp
table(tmp)
blas<-osb_clean[,tmp]
table(Idents(blas))
blas@meta.data[1:5,]
Idents(blas)<-blas@meta.data$Site
library(DESeq2)
lungmeta_findMarkers<-FindMarkers(blas,ident.1 = "lungmeta",ident.2 = "Insitu",verbose = T,only.pos = F,min.pct = 0.1,logfc.threshold = 0.1)

blas.cells <- log1p(AverageExpression(blas, verbose = FALSE)$RNA)

####lung vs. in situ###
library(ggplot2)
library(cowplot)
?top_n
top10<-rownames(lungmeta_findMarkers[order(lungmeta_findMarkers$avg_logFC,decreasing = T),])[1:10]
botom10<-rownames(lungmeta_findMarkers[order(lungmeta_findMarkers$avg_logFC,decreasing = F),])[1:10]

?AverageExpression
genes.to.label = c(top10,botom10)
library(ggplot2)
dev.new()
p1<-ggplot(blas.cells, aes(Insitu, lungmeta)) + geom_point(size=0.5) + 
  ggtitle("Osteoblast (lung metastasis vs in situ) tumor")+theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
LabelPoints(plot = p1, points = genes.to.label, repel = TRUE,colour="red")

####enrichment in the lung###
lungmeta_findMarkers$gene<-rownames(lungmeta_findMarkers)
lungmeta_marker<-lungmeta_findMarkers %>% filter(p_val_adj<0.05) %>% filter(avg_logFC>0)

library(clusterProfiler)
library(msigdbr)

m_t2g_C5 <- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "BP") %>% 
  dplyr::select(gs_name, gene_symbol)

m_t2g_H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

?msigdbr
m_t2g_C7 <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, gene_symbol)
m_t2g_C2 <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = c("CP:KEGG")) %>% 
  dplyr::select(gs_name, gene_symbol) ###signating pathways
###group1##
em_lung <- enricher(lungmeta_marker$gene, TERM2GENE=m_t2g_C5)
em_lung_H <- enricher(lungmeta_marker$gene, TERM2GENE=m_t2g_H)
###simplify(em_group1)
head(em_lung)
head(em_lung_H)
em_lung<-simplify(em_lung)
barplot(em_lung ,showCategory = 15)
barplot(em_lung_H,showCategory = 15)
write.csv(em_lung,file="3_Figures_Tables/Figure_20417/lung_enrichment.csv")
write.csv(em_lung_H,file="3_Figures_Tables/Figure_20417//lung_enrichment_H.csv")

#####recurrence###
table(Idents(blas))
recurrence_findMarkers<-FindMarkers(blas,ident.1 = "recurrence",ident.2 = "Insitu",verbose = T,only.pos = F,min.pct = 0.1,logfc.threshold = 0.1)
top10<-rownames(recurrence_findMarkers[order(recurrence_findMarkers$avg_logFC,decreasing = T),])[1:10]
botom10<-rownames(recurrence_findMarkers[order(recurrence_findMarkers$avg_logFC,decreasing = F),])[1:10]

?AverageExpression
genes.to.label = c(top10,botom10)

dev.off()
p1<-ggplot(blas.cells, aes(Insitu, recurrence)) + geom_point(size=0.5) + ggtitle("Osteoblast (recurrence vs in situ)")+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
LabelPoints(plot = p1, points = genes.to.label, repel = TRUE,colour="red")

###GO analysis of differentially expressed genes##

recurrence_findMarkers$gene<-rownames(recurrence_findMarkers)
recurrence_marker<-recurrence_findMarkers %>% filter(p_val_adj<0.05) %>% filter(avg_logFC>0)

em_recurrence <- enricher(recurrence_marker$gene, TERM2GENE=m_t2g_C5)
em_recurrence_H <- enricher(recurrence_marker$gene, TERM2GENE=m_t2g_H)

head(em_recurrence)
head(em_recurrence_H)
barplot(em_recurrence,showCategory = 15)
barplot(em_recurrence_H,showCategory = 50)
write.csv(em_recurrence,file="3_Figures_Tables/Figure_20417//recurrence_enrichment.csv")
write.csv(em_recurrence_H,file="3_Figures_Tables/Figure_20417//recurrence_enrichment_H.csv")


cellID<-bitr(recurrence_marker$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
em_rec<- enrichGO(cellID$ENTREZID, 'org.Hs.eg.db',ont="BP")
head(em_rec)
em_rec<-simplify(em_rec)
em_rec
barplot(em_rec,showCategory = 20)
write.csv(em_rec,file="3_Figures_Tables/Figure_20417//recurrence_enrichmentBP.csv")

cellID_2<-bitr(lungmeta_marker$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
em_lung<- enrichGO(cellID_2$ENTREZID, 'org.Hs.eg.db',ont="BP")
head(em_lung)
em_lung<-simplify(em_lung)
barplot(em_lung,showCategory = 20)
write.csv(em_lung,file="3_Figures_Tables/Figure_20417//lung_enrichmentBP.csv")



####select the chrocytes####
table(Idents(osb_clean))
tmp_chro<-ifelse((Idents(osb_clean)=="Chondrocyte"),TRUE,FALSE)
table(tmp_chro)

chon<-osb_clean[,tmp_chro]


chon<-NormalizeData(chon,verbose = FALSE) 
chon<-FindVariableFeatures(chon,selection.method = "vst", nfeatures = 3000)
chon<-ScaleData(chon,vars.to.regress = "percent.mt",verbose = FALSE)
chon<-RunPCA(chon,verbose = T)
p1 <- DimPlot(object = chon, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = chon, features = "PC_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p1,p2))
chon<-RunHarmony(chon,"orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(chon, 'harmony')
dim(harmony_embeddings)
p3 <- DimPlot(object = chon, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p4 <- VlnPlot(object = chon, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p3,p4))


chon<- chon %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) 
 
chon<-FindClusters(chon,resolution = 1.2)
DimPlot(chon, reduction = "tsne", pt.size = .1,label = T)
DimPlot(chon, reduction = "umap", pt.size = .1, split.by = 'Site')
DimPlot(chon, reduction = "tsne", pt.size = .1, split.by = 'orig.ident')
#tab<-table(chon@meta.data$orig.ident)
#write.csv(tab,file="chon_from_sample.csv")
DimPlot(chon, reduction = "umap", label = TRUE, pt.size = .1,split.by = "Site")
DimPlot(chon, reduction = "umap", label = TRUE, pt.size = .1)
FeaturePlot(chon,features = c("ACAN"),reduction = "tsne",max.cutoff =3,label =T,ncol = 3)
FeaturePlot(chon,features = c("SOX9","TOP2A","ACAN","COL2A1","COL1A1"),reduction = "tsne",max.cutoff =3,label =T,ncol = 3)
FeaturePlot(chon,features = c("PTH1R","IHH","MEF2C","RUNX2","SOX9","ACAN","COL10A1","COL2A1"),reduction = "tsne",max.cutoff =3,label =F)
FeaturePlot(chon,features = c("COL2A1","CD74","CDH1"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(chon,features = c("RUNX2","SP7","ALPL","SPP1","SOX9","ACAN","GREM1","PTH1R","COL10A1"),reduction = "tsne",max.cutoff =3,label =T)


tmp_chon<-ifelse(Idents(chon)=="15"|Idents(chon)=="17",FALSE,TRUE)###remove the contamination cells
table(tmp_chon)
chon_clean<-chon[,tmp_chon]

chon_clean<-NormalizeData(chon_clean,verbose = FALSE) 
chon_clean<-FindVariableFeatures(chon_clean,selection.method = "vst", nfeatures = 3000)
rownames(chon_clean)[1:5]
chon_clean<-ScaleData(chon_clean,vars.to.regress = "percent.mt",verbose = FALSE,features = rownames(chon_clean))
chon_clean<-RunPCA(chon_clean,verbose = T)
p1 <- DimPlot(object = chon_clean, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = chon_clean, features = "PC_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p1,p2))
chon_clean<-RunHarmony(chon_clean,"orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(chon_clean, 'harmony')
dim(harmony_embeddings)
p3 <- DimPlot(object = chon_clean, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p4 <- VlnPlot(object = chon_clean, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p3,p4))


chon_clean<- chon_clean %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.6) %>% 
  identity()

chon_clean<-FindClusters(chon_clean,resolution = 0.6)

DimPlot(chon_clean, reduction = "tsne", pt.size = .1, split.by = 'Site',label = T)
DimPlot(chon_clean, reduction = "tsne", pt.size = .1, split.by = 'orig.ident',label = T)
DimPlot(chon_clean, reduction = "tsne", pt.size = .1,label = T,cols = cors)
DimPlot(chon_clean, reduction = "umap", pt.size = .1,label = T)
DimPlot(chon_clean, reduction = "umap", pt.size = .1,label = T,split.by = "Site")
tab<-table(chon_clean@meta.data$orig.ident)
write.csv(tab,file="chon_from_sample.csv")
table(osb_clean$orig.ident)
Idents(chon_clean)<-chon_clean$RNA_snn_res.0.6
FeaturePlot(chon_clean,features = c("SOX9","TOP2A","ACAN","COL2A1"),reduction = "umap",max.cutoff =3,label =T)
FeaturePlot(chon_clean,features = c("PTH1R","IHH","MEF2C","RUNX2","SOX9","ACAN","COL10A1","COL2A1","GREM1","ALPL","SP7","SPP1"),reduction = "umap",max.cutoff =3,label =T)
FeaturePlot(chon_clean,features = c("IHH","COL9A1","TOP2A","UBE2C"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(chon_clean,features = c("COL10A1","WIF1","IHH","TGFBI"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(chon_clean,features = c("PCNA","TOP2A","TYMS","MKI67"),reduction = "tsne",max.cutoff =3,label =F)
FeaturePlot(chon_clean,features = c("CPE","CPM","COL9A1","CXCL12","THY1"),reduction = "tsne",max.cutoff =3,label =F)
FeaturePlot(chon_clean,features = c("SYT8","COL9A1","RARRES1","TIMP1","IFITM2"),reduction = "tsne",max.cutoff =3,label =F)
FeaturePlot(chon_clean,features = c("SPP1","RUNX2","SP7","SOX9","COL2A1","COL1A1"),reduction = "tsne",max.cutoff =3,label =F)
FeaturePlot(chon_clean,features = c("OSR2","CPM","PRG4","SCARA3","COL1A1"),reduction = "tsne",max.cutoff =3,label =F)
FeaturePlot(chon_clean,features = c("SP7","ALPL","OCN","RUNX2"),reduction = "tsne",max.cutoff =3,label =F)
FeaturePlot(chon_clean,features = c("MEF2C","PTH1R", "IHH","COL9A1","ALPL"),reduction = "tsne",max.cutoff =3,label =F)
chon_markers<-FindAllMarkers(chon_clean,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)


library(singleseqgset)
h.human <- msigdbr(species="Homo sapiens",category="H")
h.names <- unique(h.human$gs_name)
h.sets <- vector("list",length=length(h.names))
names(h.sets) <- h.names
for (i in names(h.sets)) {
  h.sets[[i]] <- pull(h.human[h.human$gs_name==i,"gene_symbol"])
}

table(Idents(chon_clean))
chon_clean$real_id<-Idents(chon_clean)
logfc.data <- logFC(cluster.ids=chon_clean@meta.data$real_id,expr.mat=chon_clean@assays$RNA@data)
names(logfc.data)
gse.res <- wmw_gsea(expr.mat=chon_clean@assays$RNA@data,cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],gene.sets=h.sets)
res.stats <- gse.res[["GSEA_statistics"]]
res.stats

ac<-data.frame(cluster=factor(colnames(res.stats)))
ac
rownames(ac)=colnames(res.stats)
table(Idents(chon_clean))
ann_colors = list(cluster = c("Chondro_Pro"="#1F77B4FF","Chondro_pre_1"="#FF7F0EFF", "Chondro_pre_2"="#2CA02CFF",
                              "Chondro_hyper"="#D62728FF"))

dev.off()
pheatmap::pheatmap(res.stats,fontsize_row = 8,annotation_col = ac,annotation_legend = F,annotation_colors = ann_colors,cluster_rows = F,
                   cluster_cols = F,scale = "row")


Idents(chon_clean)<-chon_clean$RNA_snn_res.0.6
chon_clean<-RenameIdents(chon_clean,"7"="Chondro_Pro","4"="Chondro_hyper_1","5"="Chondro_hyper_2","8"="Chondro_hyper_2","10"="Chondro_hyper_2","9"="Chondro_hyper_2","0"="Chondro_hyper_2","6"="Chondro_hyper_2",
                         "2"="Chondro_trans_diff","3"="Chondro_trans_diff","1"="Chondro_trans_diff","11"="Chondro_trans_diff")
table(Idents(chon_clean))


chon_clean@meta.data[["harmony_cell_subgroup"]]<-paste("chon_clean",Idents(chon_clean),sep="_")
chon_clean@meta.data[1:5,]
write.csv(chon_clean@meta.data,file ="3_Figures_Tables/Figure_20417/chon_clean_0419.csv")
table(Idents(chon_clean))








library(ggsci)
cors<-pal_d3()(10)
cors

DimPlot(chon_clean, reduction = "tsne", pt.size = .1, split.by = 'Site',cols = cors)
DimPlot(chon_clean, reduction = "tsne", pt.size = .1,label = F,cols = cors)
DimPlot(chon_clean, reduction = "umap", pt.size = .1,label = F,cols = cors)
DimPlot(chon_clean, reduction = "tsne", pt.size = .1, split.by = 'orig.ident')
DimPlot(chon_clean, reduction = "umap", label = TRUE, pt.size = .1,cols = cors,split.by = "Site")
tabe<-table(chon_clean$orig.ident,Idents(chon_clean))
prop.table(tabe,1)
write.csv(prop.table(tabe,1),file="./3_Figures_Tables/Figure_20417//distinct_chondro_subgroup_in_samples.csv")


####select the osteoblast####
table(Idents(osb_clean))

tmp_osb<-ifelse((Idents(osb_clean)=="Osteoblast_pro_1"|Idents(osb_clean)=="Osteoblast_pro_2"|Idents(osb_clean)=="Osteoblast_1"|Idents(osb_clean)=="Osteoblast_2"|
                   Idents(osb_clean)=="Osteoblast_3"|Idents(osb_clean)=="Osteoblast_4"),TRUE,FALSE)
table(tmp_osb)

osb_row<-osb_clean[,tmp_osb]
DimPlot(osb_row, reduction = "tsne", pt.size = .1, split.by = 'Site',label = T)
DimPlot(osb_row_2, reduction = "tsne", pt.size = .1,label = T)
DimPlot(osb_row_2, reduction = "umap", pt.size = .1,label = T)
DimPlot(osb_row, reduction = "tsne", pt.size = .1, split.by = 'orig.ident')

FeaturePlot(osb_row,features = c("COL1A1","COL2A1","SOX9","ACAN"),reduction = "tsne",label = T)
FeaturePlot(osb_row,features = c("ACAN","IBSP","TCN"),reduction = "tsne",label = T)
FeaturePlot(osb_row,features = c("CD74","CD14","CD3D","MS4A1"),reduction = "tsne",label = T)
FeaturePlot(osb_row,features = c("CTSK","EPCAM","MYH11","CXCL12"),reduction = "tsne",label = T)
FeaturePlot(osb_row,features = c("CD44","CTSC","ID3","RUNX2","COL1A1"),reduction = "tsne",label = T)
FeaturePlot(osb_row,features = c("COL1A1","SPP1","ACAN"),reduction = "umap",label = T)
FeaturePlot(osb_row,features = c("CD14","CD74",'FCER1A'),reduction = "umap",label = T)
FeaturePlot(osb_row,features = c("CDH11","ACTA2","RUNX2","DCN","VIM","LUM"),reduction = "umap",label = T,max.cutoff = 3)
FeaturePlot(osb_row,features = c("LY6A","THY1","CD44","PDGFRA","CD34"),reduction = "umap",label = T)
FeaturePlot(osb_row,features = c("SPP1","NT5E","CSPG4"),reduction = "umap",label = T)

osb_row[["cell_type_seurat"]]<-Idents(osb_row)

table(osb_row$Harmony_group,osb_row$harmony_cell_subgroup)

write.csv(osb_row@meta.data,file="osteoblast_clean_0420.csv")


####functional analysis####
markers_osb_pos<-FindAllMarkers(osb_row,logfc.threshold = 0.1,min.pct = 0.1,only.pos = T)
markers_osb_total<-FindAllMarkers(osb_row,logfc.threshold = 0.1,min.pct = 0.1,only.pos = F)
markers_osb_roc<-FindAllMarkers(osb_row,test.use = "roc",only.pos = T)

FeaturePlot(osb_row,features = c("PCNA","NEAT1","IFITM5"),reduction = "tsne",label = T)
FeaturePlot(osb_row,features = c("RUNX2","SP7","CTSK"),reduction = "tsne",label = T,max.cutoff = 2)
FeaturePlot(osb_row,features = c("SPP1","PTH1R","ALPL","COL1A1"),reduction = "tsne",label = T)
FeaturePlot(osb_row,features = c("FTL","FN1","CD74"),reduction = "tsne",label = T,max.cutoff = 5)
FeaturePlot(osb_row,features = c("PTH1R","ALPL","COL1A1","SPP1"),reduction = "tsne",label = T)

####GESA plot of the H_markers
library(singleseqgset)
h.human <- msigdbr(species="Homo sapiens",category="H")

h.names <- unique(h.human$gs_name)

h.sets <- vector("list",length=length(h.names))
names(h.sets) <- h.names

for (i in names(h.sets)) {
  h.sets[[i]] <- pull(h.human[h.human$gs_name==i,"gene_symbol"])
}


table(Idents(osb_row))
osb_row$real_id<-Idents(osb_row)
logfc.data <- logFC(cluster.ids=osb_row@meta.data$real_id,expr.mat=osb_row@assays$RNA@data)
names(logfc.data)
gse.res <- wmw_gsea(expr.mat=osb_row@assays$RNA@data,cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],gene.sets=h.sets)
res.stats <- gse.res[["GSEA_statistics"]]
res.stats

ac<-data.frame(cluster=factor(colnames(res.stats)))
ac
rownames(ac)=colnames(res.stats)

ann_colors = list(cluster = c("Osteoblast_pro_1"="#1F77B4FF","Osteoblast_pro_2"="#FF7F0EFF", "Osteoblast_1"="#2CA02CFF",
                              "Osteoblast_2"="#D62728FF","Osteoblast_3"="#9467BDFF","Osteoblast_4"="#8C564BFF"))

dev.off()
pheatmap::pheatmap(res.stats,fontsize_row = 8,annotation_col = ac,annotation_legend = F,annotation_colors = ann_colors,scale = "row",cluster_rows = F,
                   cluster_cols = F)






#####cell proportion among the samples
tab1<-table(osb_row$orig.ident,Idents(osb_row))
write.csv(tab1,file = "./3_Figures_Tables/Figure_20417/osteoblast_cells_among_samples.csv")
tab2<-prop.table(tab1,1)
tab2
write.csv(tab2,file = "./3_Figures_Tables/Figure_20417/osteoblast_cells_proportion_among_samples.csv")

markers_osb_pos<-FindAllMarkers(osb_row,logfc.threshold = 0.25,min.pct = 0.25,only.pos = T)
top15<-markers_osb_pos %>% group_by(cluster) %>% top_n(n=15,wt=avg_logFC)
table(Idents(osb_row))
aver_osb_row<-AverageExpression(osb_row)
top15$gene

cc<-aver_osb_row$RNA[unique(top15$gene),]

pheatmap::pheatmap(cc,scale = "row",cluster_rows = F,cluster_cols = F,fontsize_row = 6)


osb_row@meta.data[1:5,]
write.csv(osb_row@meta.data,file="cell_subgroup_osteoblast_cells0420.csv")



####analysis the cells between the chondrocytes and the osteoblast####
table(Idents(osb_clean))
osb_clean@meta.data[1:5,]

type<-ifelse(Idents(osb_clean)=="Chondrocyte","Chondrocyte","Osteoblast")
table(type)
Idents(osb_clean)<-type
table(Idents(osb_clean))


markers<-FindMarkers(osb_clean,ident.1 = "Chondrocyte",ident.2 = "Osteoblast",logfc.threshold = 0.25,min.pct = 0.25)
total_2_cells <- log1p(AverageExpression(osb_clean, verbose = FALSE)$RNA)

####lung vs. in situ###
library(ggplot2)
library(cowplot)
?top_n
top10<-rownames(markers[order(markers$avg_logFC,decreasing = T),])[1:10]
botom10<-rownames(markers[order(markers$avg_logFC,decreasing = F),])[1:10]

?AverageExpression
genes.to.label = c(top10,botom10)

dev.new()
p1<-ggplot(total_2_cells, aes(Osteoblast, Chondrocyte)) + geom_point(size=0.5) + ggtitle("Chondrocyte vs. Osteoblast")+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
LabelPoints(plot = p1, points = genes.to.label, repel = TRUE,colour="red")

markers$gene<-rownames(markers)
markers_chondro<-markers %>% filter(p_val_adj<0.05) %>% filter(avg_logFC>0)
markers_osteoblast<-markers %>% filter(p_val_adj<0.05) %>% filter(avg_logFC<0)
library(clusterProfiler)
library(msigdbr)

m_t2g_C5 <- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "BP") %>% 
  dplyr::select(gs_name, gene_symbol)

m_t2g_H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

?msigdbr
m_t2g_C7 <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, gene_symbol)
m_t2g_C2 <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = c("CP:KEGG")) %>% 
  dplyr::select(gs_name, gene_symbol) ###signating pathways
###group1##
markers_chondro_enrich<- enricher(markers_chondro$gene, TERM2GENE=m_t2g_C5)
markers_osteoblast_enrich <- enricher(markers_osteoblast$gene, TERM2GENE=m_t2g_C5)
###simplify(em_group1)
barplot(markers_chondro_enrich,showCategory = 30)
barplot(markers_osteoblast_enrich,showCategory = 30)
write.csv(markers_chondro_enrich,file="3_Figures_Tables/Figure_20417/chondrocyte_enrichment.csv")
write.csv(markers_osteoblast_enrich,file="3_Figures_Tables/Figure_20417/osteoblast_enrichment_H.csv")


####we get the BC20 and BC22 cells for the trajectory analysis####

table(osb_clean$orig.ident)
table(Idents(osb_clean))
BC20_tmp<-osb_clean$orig.ident=="BC20"
BC22_tmp<-osb_clean$orig.ident=="BC22"

BC20<-osb_clean[,BC20_tmp]
BC22<-osb_clean[,BC22_tmp]
table(Idents(BC20))

table(BC20$harmony_cell_subgroup)
table(BC22$harmony_cell_subgroup)
length(chon_clean@assays$RNA@var.features)

#bc20<-BC20@meta.data
#dim(bc20)

table(chon_clean$harmony_cell_subgroup)


BC20$cell_group_row<-Idents(BC20)
write.csv(BC20@meta.data,file="BC20_meta_infor_0507.csv")
BC22$cell_group_row<-Idents(BC22)
write.csv(BC22@meta.data,file="BC22_meta_infor_0507.csv")

table(Idents(BC22))

bc20<-as.matrix(BC20@meta.data)
class(bc20)
chon_infor<-chon_clean@meta.data
chon_infor_bc20<-chon_infor[chon_infor$orig.ident=="BC20",]
table(rownames(chon_infor_bc20)%in%rownames(bc20))
bc20[rownames(chon_infor_bc20),14]<-chon_infor_bc20$harmony_cell_subgroup

dim(bc20)
class(bc20)
BC20
table(rownames(BC20@meta.data)==rownames(bc20))
head(BC20@meta.data)
BC20$detailed_cell<-bc20[,14]
head(BC20@meta.data)

bc22<-as.matrix(BC22@meta.data)
class(bc22)
chon_infor<-chon_clean@meta.data
chon_infor_bc22<-chon_infor[chon_infor$orig.ident=="BC22",]
table(rownames(chon_infor_bc22)%in%rownames(bc22))
bc22[rownames(chon_infor_bc22),14]<-chon_infor_bc22$harmony_cell_subgroup

table(rownames(BC22@meta.data)==rownames(bc22))
head(BC22@meta.data)
BC22$detailed_cell<-bc22[,14]
head(BC22@meta.data)
table(BC20$detailed_cell,BC20$harmony_cell_subgroup)



####perform the monocle analysis of the condrocyte and osteoblast cell in BC20 and BC22####
data.countbc20<-as.matrix(BC20@assays$RNA@counts)
data.countbc22<-as.matrix(BC22@assays$RNA@counts)
dim(data.countbc20)[2]
dim(data.countbc22)[2]
# filter out genes with zero expression across 25% of cells
keep <- rowSums(data.countbc20> 0) >= round(dim(data.countbc20)[2]/4)
sum(keep)
data.countbc20 <- data.countbc20[keep,] 
#dim(data.countbc20)
# filter out genes with very low average nonzero expression across all cells
gene_ann <- data.frame(
  gene_short_name = row.names(data.countbc20), 
  row.names = row.names(data.countbc20))

sample_ann <- BC20@meta.data
dim(sample_ann)
library(monocle)
fd <- new("AnnotatedDataFrame",data=gene_ann)
pd<-new("AnnotatedDataFrame",data=sample_ann)

sc_cds_2 <- newCellDataSet(data.countbc20,  phenoData = pd,featureData =fd,expressionFamily = negbinomial.size(),lowerDetectionLimit=0.1)
sc_cds_2 <- estimateSizeFactors(sc_cds_2)
sc_cds_2 <- estimateDispersions(sc_cds_2)
sc_cds_2 <- detectGenes(sc_cds_2, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(sc_cds_2), num_cells_expressed >= 5))

# vVerifying normalization
#library(reshape2)
#L <- log(exprs(sc_cds_2[expressed_genes,]))
#melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))



disp_table <- dispersionTable(sc_cds_2) # 挑有差异的
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1) # 挑表达量不太低的
sc_cds_2 <- setOrderingFilter(sc_cds_2, unsup_clustering_genes$gene_id)  # 准备聚类基因名单
plot_ordering_genes(sc_cds_2) 

plot_pc_variance_explained(sc_cds_2, return_all = F)
sc_cds_2 <- reduceDimension(sc_cds_2, max_components = 2, num_dim = 5,reduction_method = 'tSNE', 
                            residualModelFormulaStr = "~num_genes_expressed",norm_method = 'log',
                            verbose = T)

sc_cds_2 <- clusterCells(sc_cds_2, num_clusters = 5)
p1<-plot_cell_clusters(sc_cds_2, 1, 2, color = "Cluster")
p1
pData(sc_cds_2)[1:5,]
p2<-plot_cell_clusters(sc_cds_2, 1, 2, color = "detailed_cell")
CombinePlots(list(p1,p2))



diff_test_res <- differentialGeneTest(sc_cds_2[expressed_genes,],
                                      fullModelFormulaStr = "~as.factor(cell_group_row)",cores = 3,
                                      verbose = T)

table(Idents(BC20))

ordering_genes <- row.names(subset(diff_test_res, qval < 1e-2))

####根据variance来选择排序####
library(dplyr)
disp_table <- dispersionTable(sc_cds_2)
tmp<-disp_table$mean_expression >= 0.1&disp_table$dispersion_empirical >= 1*disp_table$dispersion_fit
ordering_genes_3 <- as.character(disp_table[tmp,]$gene_id)
####end####


length(ordering_genes)
sc_cds2 <- setOrderingFilter(sc_cds_2, ordering_genes)
plot_ordering_genes(sc_cds2)
sc_cds2<- reduceDimension(sc_cds2, max_components = 2, reduction_method  = "DDRTree")
sc_cds2 <- orderCells(sc_cds2)

plot_cell_trajectory(sc_cds2, color_by = "Pseudotime",show_branch_points = F)
plot_cell_trajectory(sc_cds2, color_by = "State",show_branch_points = F)
plot_cell_trajectory(sc_cds2, color_by = "cell_group_row",show_branch_points = F)+facet_wrap(~cell_group_row)
plot_cell_trajectory(sc_cds2, color_by = "detailed_cell",show_branch_points = F)+facet_wrap(~detailed_cell)
plot_cell_trajectory(sc_cds2, color_by = "harmony_cell_subgroup",show_branch_points = F)
plot_cell_trajectory(sc_cds2, markers = c("ACAN","COL1A1","COL2A1","SOX9","RUNX2"),use_color_gradient = T,show_branch_points = F)
plot_cell_trajectory(sc_cds2, markers = c("SOX9","ACAN","RUNX2"),use_color_gradient = T,show_branch_points = F)


BEAM_res <- BEAM(sc_cds2[expressed_genes, ], branch_point = 1, cores = 6)


BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
ccc<-plot_genes_branched_heatmap(sc_cds2[row.names(subset(BEAM_res,
                                                          qval < 1e-10)),],
                                 branch_point = 1,
                                 num_clusters = 4,
                                 cores = 4,
                                 use_gene_short_name = T,return_heatmap = T,
                                 show_rownames = T)

clusters<-ccc$annotation_row
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)

write.csv(clustering, file="./5.OSteoblast/gene_BC20_cluster.csv")

###Functional analysis###
library(msigdbr)
library(clusterProfiler)
msigdbr_show_species()  
m_df <- msigdbr(species = "Homo sapiens") 

table(m_df$gs_cat)
head(m_df,2)%>%as.data.frame()

m_t2g_C5 <- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "BP") %>% 
  dplyr::select(gs_name, gene_symbol)
?msigdbr
#m_t2g_C7 <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
# dplyr::select(gs_name, gene_symbol)
m_t2g_C2 <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = c("CP:KEGG")) %>% 
  dplyr::select(gs_name, gene_symbol) ###signating pathways
###group1##
gene1<-rownames(clustering)[clustering[,1]==1]
em_group1 <- enricher(gene1, TERM2GENE=m_t2g_C5)
em_group_pathway_1 <- enricher(gene1, TERM2GENE=m_t2g_C2)
###simplify(em_group1)
head(em_group1)
head(em_group_pathway_1)
barplot(em_group1,showCategory = 20)
barplot(em_group_pathway_1)
write.csv(em_group1,file="5.OSteoblast/heatmap_cluster1_enrichment_BC20.csv")
write.csv(em_group_pathway_1,file="5.OSteoblast/heatmap_cluster1_enrichment_pathway_BC20.csv")

###group1##
gene2<-rownames(clustering)[clustering[,1]==2]
em_group2 <- enricher(gene2, TERM2GENE=m_t2g_C5)
em_group_pathway_2 <- enricher(gene2, TERM2GENE=m_t2g_C2)
###simplify(em_group1)
head(em_group2)
head(em_group_pathway_2)
barplot(em_group2,showCategory = 20)
barplot(em_group_pathway_2)
write.csv(em_group2,file="5.OSteoblast/heatmap_cluster2_enrichment_BC20.csv")
write.csv(em_group_pathway_2,file="5.OSteoblast/heatmap_cluster2_enrichment_pathway_BC20.csv")
###group3##
gene3<-rownames(clustering)[clustering[,1]==3]
em_group3 <- enricher(gene3, TERM2GENE=m_t2g_C5)
em_group_pathway_3 <- enricher(gene3, TERM2GENE=m_t2g_C2)
###simplify(em_group1)
head(em_group3)
head(em_group_pathway_3)
barplot(em_group3,showCategory = 50)
barplot(em_group_pathway_1)
write.csv(em_group3,file="5.OSteoblast/heatmap_cluster3_enrichment_BC20.csv")
write.csv(em_group_pathway_3,file="5.OSteoblast/heatmap_cluster3_enrichment_pathway_BC20.csv")
###group4##
gene4<-rownames(clustering)[clustering[,1]==4]
em_group4 <- enricher(gene4, TERM2GENE=m_t2g_C5)
em_group_pathway_4 <- enricher(gene4, TERM2GENE=m_t2g_C2)
###simplify(em_group1)
head(em_group4)
head(em_group_pathway_4)
barplot(em_group4,showCategory = 50)
barplot(em_group_pathway_4)
write.csv(em_group4,file="5.OSteoblast/heatmap_cluster4_enrichment_BC20.csv")
write.csv(em_group_pathway_4,file="5.OSteoblast/heatmap_cluster4_enrichment_pathway_BC20.csv")










####分析BC22####



####perform the monocle analysis of the condrocyte and osteoblast cell in BC22####
data.countbc22<-as.matrix(BC22@assays$RNA@counts)
dim(data.countbc22)[2]
# filter out genes with zero expression across 25% of cells
keep <- rowSums(data.countbc22> 0) >= round(dim(data.countbc22)[2]/4)
sum(keep)
data.countbc22 <- data.countbc22[keep,] 
#dim(data.countbc20)
# filter out genes with very low average nonzero expression across all cells
gene_ann <- data.frame(
  gene_short_name = row.names(data.countbc22), 
  row.names = row.names(data.countbc22))

sample_ann <- BC22@meta.data
dim(sample_ann)
library(monocle)
fd <- new("AnnotatedDataFrame",data=gene_ann)
pd<-new("AnnotatedDataFrame",data=sample_ann)

sc_cds_2 <- newCellDataSet(data.countbc22,  phenoData = pd,featureData =fd,expressionFamily = negbinomial.size(),lowerDetectionLimit=0.1)
sc_cds_2 <- estimateSizeFactors(sc_cds_2)
sc_cds_2 <- estimateDispersions(sc_cds_2)
sc_cds_2 <- detectGenes(sc_cds_2, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(sc_cds_2), num_cells_expressed >= 5))



disp_table <- dispersionTable(sc_cds_2) 
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1) 
sc_cds_2 <- setOrderingFilter(sc_cds_2, unsup_clustering_genes$gene_id)  
plot_ordering_genes(sc_cds_2) 

plot_pc_variance_explained(sc_cds_2, return_all = F)
sc_cds_2 <- reduceDimension(sc_cds_2, max_components = 2, num_dim = 5,reduction_method = 'tSNE', 
                            residualModelFormulaStr = "~num_genes_expressed",norm_method = 'log',
                            verbose = T)

sc_cds_2 <- clusterCells(sc_cds_2, num_clusters = 5)
p1<-plot_cell_clusters(sc_cds_2, 1, 2, color = "Cluster")
p1
pData(sc_cds_2)[1:5,]
p2<-plot_cell_clusters(sc_cds_2, 1, 2, color = "cell_group_row")
CombinePlots(list(p1,p2))



diff_test_res <- differentialGeneTest(sc_cds_2[expressed_genes,],
                                      fullModelFormulaStr = "~as.factor(cell_group_row)",cores = 3,
                                      verbose = T)
ordering_genes <- row.names(subset(diff_test_res, qval < 1e-2))
sc_cds2 <- setOrderingFilter(sc_cds_2, ordering_genes)
plot_ordering_genes(sc_cds2)
sc_cds2<- reduceDimension(sc_cds2, max_components = 2, reduction_method  = "DDRTree")
sc_cds2 <- orderCells(sc_cds2,root_state = 3)

plot_cell_trajectory(sc_cds2, color_by = "Pseudotime",show_branch_points = T)
plot_cell_trajectory(sc_cds2, color_by = "State",show_branch_points = F)
plot_cell_trajectory(sc_cds2, color_by = "detailed_cell",show_branch_points = F)+facet_wrap(~detailed_cell)
plot_cell_trajectory(sc_cds2, color_by = "harmony_cell_subgroup",show_branch_points = F)

plot_cell_trajectory(sc_cds2, markers = c("ACAN","COL1A1","COL2A1","SOX9","RUNX2"),use_color_gradient = T,show_branch_points = F)
plot_cell_trajectory(sc_cds2, markers = c(),use_color_gradient = T,show_branch_points = F)


BEAM_res <- BEAM(sc_cds2[expressed_genes, ], branch_point = 1, cores = 4)



BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
ccc<-plot_genes_branched_heatmap(sc_cds2[row.names(subset(BEAM_res,
                                                          qval < 1e-10)),],
                                 branch_point = 1,
                                 num_clusters = 4,
                                 cores = 4,
                                 use_gene_short_name = T,return_heatmap = T,
                                 show_rownames = T)

clusters<-ccc$annotation_row
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)

write.csv(clustering, file="./5.OSteoblast/gene_BC20_cluster.csv")



###Functional analysis###
library(msigdbr)
library(clusterProfiler)
msigdbr_show_species()  
m_df <- msigdbr(species = "Homo sapiens") 

table(m_df$gs_cat)
head(m_df,2)%>%as.data.frame()

m_t2g_C5 <- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "BP") %>% 
  dplyr::select(gs_name, gene_symbol)
?msigdbr
#m_t2g_C7 <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
# dplyr::select(gs_name, gene_symbol)
m_t2g_C2 <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = c("CP:KEGG")) %>% 
  dplyr::select(gs_name, gene_symbol) ###signating pathways
###group1##
gene1<-rownames(clustering)[clustering[,1]==1]
em_group1 <- enricher(gene1, TERM2GENE=m_t2g_C5)
em_group_pathway_1 <- enricher(gene1, TERM2GENE=m_t2g_C2)
###simplify(em_group1)
head(em_group1)
head(em_group_pathway_1)
barplot(em_group1,showCategory = 20)
barplot(em_group_pathway_1)
write.csv(em_group1,file="5.OSteoblast/heatmap_cluster1_enrichment_BC22.csv")
write.csv(em_group_pathway_1,file="5.OSteoblast/heatmap_cluster1_enrichment_pathway_BC22.csv")

###group1##
gene2<-rownames(clustering)[clustering[,1]==2]
em_group2 <- enricher(gene2, TERM2GENE=m_t2g_C5)
em_group_pathway_2 <- enricher(gene2, TERM2GENE=m_t2g_C2)
###simplify(em_group1)
head(em_group2)
head(em_group_pathway_2)
barplot(em_group2,showCategory = 20)
barplot(em_group_pathway_2)
write.csv(em_group2,file="5.OSteoblast/heatmap_cluster2_enrichment_BC22.csv")
write.csv(em_group_pathway_2,file="5.OSteoblast/heatmap_cluster2_enrichment_pathway_BC22.csv")
###group3##
gene3<-rownames(clustering)[clustering[,1]==3]
em_group3 <- enricher(gene3, TERM2GENE=m_t2g_C5)
em_group_pathway_3 <- enricher(gene3, TERM2GENE=m_t2g_C2)
###simplify(em_group1)
head(em_group3)
head(em_group_pathway_3)
barplot(em_group3,showCategory = 50)
barplot(em_group_pathway_1)
write.csv(em_group3,file="5.OSteoblast/heatmap_cluster3_enrichment_BC22.csv")
write.csv(em_group_pathway_3,file="5.OSteoblast/heatmap_cluster3_enrichment_pathway_BC22.csv")
###group4##
gene4<-rownames(clustering)[clustering[,1]==4]
em_group4 <- enricher(gene4, TERM2GENE=m_t2g_C5)
em_group_pathway_4 <- enricher(gene4, TERM2GENE=m_t2g_C2)
###simplify(em_group1)
head(em_group4)
head(em_group_pathway_4)
barplot(em_group4,showCategory = 50)
barplot(em_group_pathway_4)
write.csv(em_group4,file="5.OSteoblast/heatmap_cluster4_enrichment_BC22.csv")
write.csv(em_group_pathway_4,file="5.OSteoblast/heatmap_cluster4_enrichment_pathway_BC22.csv")


save.image(file = "OS_analysis_0420.RData")
load(file = "OS_analysis_0420.RData")
BiocManager::install("infercnv")



tab1<-table(Idents(osb_clean),osb_clean$orig.ident)
tab1<-tab1[,c(5,9,10,11,1,3,7,2,4,6,8)]
tab1
tab2<-data.frame(tab1)
tab2
number<-ggplot(tab2,aes(x=Var2,y=Freq,fill=Var1))+ 
  geom_bar(stat='identity',position='stack',alpha=.5)+ 
  labs(title='',x='sample ID',y='Number')+ 
  theme(legend.justification = 'right', 
        legend.position = 'right', 
        legend.key.height = unit(0.1,'cm'),
        panel.background = element_blank(),
        axis.line=element_line(size=0.5,colour="black")
  )+ 
  scale_fill_manual(values=cors)
number
library(ggThemeAssist)
ggThemeAssistGadget(number)


tab1
tab2<-prop.table(tab1,2)*100
tab2<-data.frame(tab2)

ggplot(tab2,aes(x=Var2,y=Freq,fill=Var1))+ 
  geom_bar(stat='identity',position='stack',alpha=.5)+ 
  labs(title='',x='sampleID',y='Cell Proportion (%)')+ 
  theme(legend.justification = 'right', 
        legend.position = 'right', 
        legend.key.height = unit(0.1,'cm'),
        panel.background = element_blank(),
        axis.line=element_line(size=0.5,colour="black")
  )+ 
  scale_fill_manual(values=cors)


write.csv(osb_clean@meta.data,file="osteosarcoma_0515.csv")




