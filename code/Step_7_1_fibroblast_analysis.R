rm(list = ls())
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
plan("multiprocess", workers = 4) 
options(future.globals.maxSize = 40000 * 1024^2) 
####test the harmony
###
library(devtools)
library(harmony)
sce.big<-readRDS(file="2_output/merge_data_db_remove_harmony_group.rds")
table(sce.big@meta.data$Harmony_group)
tmp<-sce.big@meta.data$Harmony_group=="7.Fibroblast"
table(tmp)
fibro<-sce.big[,tmp]
rm(sce.big)

fibro<-NormalizeData(fibro,verbose = FALSE) 
fibro<-FindVariableFeatures(fibro,selection.method = "vst", nfeatures = 3000)
fibro<-ScaleData(fibro,vars.to.regress = "percent.mt",verbose = FALSE,features = rownames(fibro))
fibro<-RunPCA(fibro,verbose = T)
p1 <- DimPlot(object = fibro, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = fibro, features = "PC_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p1,p2))
fibro<-RunHarmony(fibro,"orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(fibro, 'harmony')
dim(harmony_embeddings)
p3 <- DimPlot(object = fibro, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p4 <- VlnPlot(object = fibro, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p3,p4))


fibro <- fibro %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) 

fibro<- FindClusters(fibro,resolution = 1.0)
DimPlot(fibro, reduction = "tsne", pt.size = .1, 
        split.by = 'Site',label = T)
DimPlot(fibro, reduction = "tsne", label = TRUE, pt.size = .1)
table(Idents(fibro))



FeaturePlot(fibro,features = c("CD74","CD14","CD3D","MS4A1","VEGFA","JCHAIN","CPA3","HBB","COL2A1","SOX9","COL9A1","IBSP"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(fibro,features = c("CD3D","VWF","COL13A1","COL1A1","COL2A1"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(fibro,features = c("DCN","THY1","ACTA2","FN1","SOX9","CXCL12","MYH11"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(fibro,features = c("SOX9","CD44","ACTA2","THY1","ANGPT1","SCX","IBSP","SPP1"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(fibro,features = c("CSPG4","NT5E","RUNX2","PTH1R","CLIP","DES"),reduction = "tsne",max.cutoff =3,label =T)

#markers<-FindAllMarkers(fibro,only.pos = T)
#table(hi)
#cells<-names(hi)[hi==4]
#DimPlot(fibro, reduction = "tsne", label = TRUE, pt.size = .1,cells.highlight = cells)


####remove the microglical cells and the doublets##
tmp<-ifelse((Idents(fibro)=="9"|Idents(fibro)=="12"|Idents(fibro)=="13"|Idents(fibro)=="14"),FALSE,TRUE)
table(tmp)

####get the pure cells##
fibro<-fibro[,tmp]
fibro<-NormalizeData(fibro,verbose = FALSE) 
fibro<-FindVariableFeatures(fibro,selection.method = "vst", nfeatures = 2000)
fibro<-ScaleData(fibro,vars.to.regress = "percent.mt",verbose = FALSE)
fibro<-RunPCA(fibro,verbose = T)
p1 <- DimPlot(object = fibro, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = fibro, features = "PC_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p1,p2))
fibro<-RunHarmony(fibro,"orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(fibro, 'harmony')
dim(harmony_embeddings)
p3 <- DimPlot(object = fibro, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p4 <- VlnPlot(object = fibro, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p3,p4))
DimPlot(object = fibro, reduction = "harmony", pt.size = .1, split.by = "Site")

fibro <- fibro %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) 

fibro<-FindClusters(fibro,resolution = 0.8)



#fibro<- FindClusters(fibro,resolution = 1.2)
DimPlot(fibro, reduction = "tsne", pt.size = .1, split.by = 'Site')
DimPlot(fibro, reduction = "tsne", pt.size = .1,label = T)
DimPlot(fibro, reduction = "tsne", pt.size = .1, split.by = 'orig.ident')
DimPlot(fibro, reduction = "umap", label = TRUE, pt.size = .1,split.by = 'Site')



FeaturePlot(fibro,features = c("DCN","THY1","CD24","SFTPA1"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(fibro,features = c("DES","CXCL12","COL14A1","MME","THY1"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(fibro,features = c("DCN","THY1","COL1A1","COL1A2","SFRP2","CLEC11A"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(fibro,features = c("COL13A1","COL14A1","DES","FABP5","MYLK","ACTA2"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(fibro,features = c("TAGLN","ACTA2","MSLN"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(fibro,features = c("GSN","PI16","CYGB",'FABP5',"APOD"),reduction = "tsne",max.cutoff =3,label =T)


marker10<-FindMarkers(fibro,ident.1 = 10,only.pos = T)
Idents(fibro)<-fibro$RNA_snn_res.0.8
FeaturePlot(fibro,features = c("CD74","CXCL12","ACTA2","DES"),reduction = "tsne",max.cutoff =6,label =F)

fibro<-RenameIdents(fibro,"2"="COL14A1+ fibroblast","4"="Smooth muscle cell","6"="Smooth muscle cell","9"="Smooth muscle cell",
                    "0"="Myofibroblast","1"="Myofibroblast","5"="Myofibroblast","3"="Myofibroblast","7"="Myofibroblast","8"="Myofibroblast","10"="Myofibroblast")
###1 with DES from myofibroblast, 
DefaultAssay(fibro)<-"RNA"
markers<-FindAllMarkers(fibro,only.pos = T, min.pct = 0.25,logfc.threshold = 0.25)
markers_total<-FindAllMarkers(fibro,only.pos = F,logfc.threshold = 0.25,min.pct = 0.25)

library(ggsci)
cors<-pal_d3()(10)
cors

DimPlot(fibro, reduction = "tsne", pt.size = .1, split.by = 'Site',cols = cors)
DimPlot(fibro, reduction = "tsne", pt.size = .1,label = F,cols = cors)
DimPlot(fibro, reduction = "tsne", pt.size = .1, split.by = 'orig.ident',ncol = 3,cols=cors)
DimPlot(fibro, reduction = "umap", label = TRUE, pt.size = .1,cols = cors)


###pan-CAF
VlnPlot(fibro,features = c("COL1A1","LUM","COL14A1","CXCL12"),cols = cors,pt.size = 0,ncol = 2)
VlnPlot(fibro,features = c("DES","ACTA2","SPP1","MYL9"),cols = cors,pt.size = 0,ncol = 2)


####Distrubition###


tab2<-table(Idents(fibro),fibro$Site)
tab2<-prop.table(tab2,2)*100
tab2
tab2<-as.data.frame(tab2)

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


markers_2<-markers[markers$p_val_adj<0.05,]
#grep(markers_2$gene,pattern = "[-]")
#marker_gene<-gsub("[.]","-",markers_2$gene)
#markers_2$marker_gene<-marker_gene

top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
top20$gene

DoHeatmap(fibro, features = top20$gene,group.colors = cors)



#####singleseqgset#####
library(singleseqgset)
h.human <- msigdbr(species="Homo sapiens",category="H")
h.names <- unique(h.human$gs_name)
h.sets <- vector("list",length=length(h.names))
names(h.sets) <- h.names
for (i in names(h.sets)) {
  h.sets[[i]] <- pull(h.human[h.human$gs_name==i,"gene_symbol"])
}

table(Idents(fibro))
fibro$real_id<-Idents(fibro)
logfc.data <- logFC(cluster.ids=fibro@meta.data$real_id,expr.mat=fibro@assays$RNA@data)
names(logfc.data)
gse.res <- wmw_gsea(expr.mat=fibro@assays$RNA@data,cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],gene.sets=h.sets)
res.stats <- gse.res[["GSEA_statistics"]]
res.stats

ac<-data.frame(cluster=factor(colnames(res.stats)))
ac
rownames(ac)=colnames(res.stats)
table(Idents(fibro))
cors
ann_colors = list(cluster = c("COL14A1+ fibroblast"="#1F77B4FF","Smooth muscle cell"="#FF7F0EFF", "Myofibroblast"="#2CA02CFF"))

dev.off()
pheatmap::pheatmap(res.stats,fontsize_row = 8,annotation_col = ac,annotation_legend = T,annotation_colors = ann_colors,cluster_rows = F,annotation_names_col = F,
                   cluster_cols = F,scale = "row")


save.image(file="fibrobalst_0424.RData")
load(file="fibrobalst_0424.RData")

table(Idents(fibro))

tab1<-table(Idents(fibro),fibro$orig.ident)
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
#library(ggThemeAssist)
#ggThemeAssistGadget(number)


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

write.csv(fibro@meta.data,file="fibroblast.csv")
