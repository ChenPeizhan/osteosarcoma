rm(list = ls())
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(harmony)
plan("multiprocess", workers = 4) 
options(future.globals.maxSize = 40000 * 1024^2) 

####Subgroup aaasyas###
#load("./harmony_integrated_0413.RData")
sce.big<-readRDS(file="2_output/merge_data_db_remove_harmony_group.rds")
table(sce.big@meta.data$Harmony_group)
tmp<-sce.big@meta.data$Harmony_group=="5.T/NK"
table(tmp)
tnk<-sce.big[,tmp]
rm(sce.big)
tnk<-NormalizeData(tnk,verbose = T) 
tnk<-FindVariableFeatures(tnk,selection.method = "vst", nfeatures = 3000)
tnk<-ScaleData(tnk,verbose = FALSE)
tnk<-RunPCA(tnk,verbose = T)
p1 <- DimPlot(object = tnk, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = tnk, features = "PC_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p1,p2))
tnk<-RunHarmony(tnk,"orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(tnk, 'harmony')
dim(harmony_embeddings)
p3 <- DimPlot(object = tnk, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p4 <- VlnPlot(object = tnk, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p3,p4))


tnk <- tnk %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  RunTSNE(reduction = "harmony", dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50)

tnk<- FindClusters(tnk,resolution = 0.8)
DimPlot(tnk, reduction = "tsne", pt.size = .1,label = T)
markers<-FindAllMarkers(tnk,only.pos = T)

tnk<-RenameIdents(tnk,"0"="CD4-/CD8-","1"="CD4-/CD8-","5"="CD4-/CD8-","8"="CD4-/CD8-","4"="T-reg","11"="Plasma",
                  "2"="CD8+T","6"="NK","7"="NK","3"="CD4+T","10"="CD4+T","9"="Proliferating")

table(Idents(tnk))
FeaturePlot(tnk,features = c("TOP2A"),reduction = "tsne",max.cutoff =3)
FeaturePlot(tnk,features = c("CD4","CD8A"),reduction = "tsne",max.cutoff =3)
FeaturePlot(tnk,features = c("TIGIT","PDCD1","CTLA4","FOXP3"),reduction = "tsne",max.cutoff =3)
FeaturePlot(tnk,features = c("TIGIT","PDCD1","CTLA4","FOXP3","HBB"),reduction = "tsne",max.cutoff =3)
FeaturePlot(tnk,features = c("NKG7","GZMB","GZMA","GNLY","CD3D","CD2","JCHAIN","CPA3","MS4A1"),reduction = "tsne",max.cutoff =3)
FeaturePlot(tnk,features = c("FCGR3A","CD3D","NKG7","NCAM1","CD8A","CD4"),reduction = "tsne",max.cutoff =3)


library(ggsci)
cors<-pal_d3()(8)
cors

###Figure 6A
DimPlot(tnk,reduction = "tsne",label = T,cols = cors,label.size = 0)
DimPlot(tnk,reduction = "umap",label = T,cols = cors,label.size = 0)
###Figure S6A
DimPlot(tnk,reduction = "tsne",split.by = "orig.ident",ncol = 3,cols = cors,label.size = 0)
DimPlot(tnk,reduction = "tsne",split.by = "Site",ncol = 3,cols = cors,label.size = 0)


###Figure 6B
write.csv(table(Idents(tnk),tnk@meta.data$Site),file="./3_Figures_Tables/Figure 6 and Supplementary Figure 6/tnk_cell_distribution.csv")


markers.to.plot <- c("CD2","IL7R","CD8A","CD4","NKG7","GZMB","GNLY","FOXP3","CTLA4","TIGIT","TNFRSF4","JCHAIN")
DotPlot(tnk, features = rev(markers.to.plot), cols = c("lightgrey","blue"), dot.scale = 8,col.max = 1)

###Figure 6D

DefaultAssay(tnk)<-"RNA"
VlnPlot(tnk,features=c("CD4","CD3D","CD8A","IL7R"),cols = cors,pt.size = 0,adjust = 1,ncol=2,y.max = 6)
VlnPlot(tnk,features=c("NKG7","GNLY","TOP2A","JCHAIN"),cols = cors,pt.size = 0,adjust = 1,ncol=2,y.max = 6)
VlnPlot(tnk,features=c("FOXP3","TIGIT","JCHAIN"),cols = cors,pt.size = 0,adjust = 1,ncol=2)


FeaturePlot(tnk,features = c("FOXP3","TIGIT","CTLA4","TNFRSF4"),reduction = "tsne",max.cutoff =3)
FeaturePlot(tnk,features = c("LUM","CD8A","CD4","FOXP3","CD25"),reduction = "tsne",max.cutoff =3)
FeaturePlot(tnk,features = c("CD3D","FCGR3A","CD56",),reduction = "tsne",max.cutoff =3)
FeaturePlot(tnk,features = c("DCL1"),reduction = "tsne",max.cutoff =3)
####Heatmap 

markers<-FindAllMarkers(tnk,only.pos = T)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(tnk, features = top10$gene,group.colors = cors)

tnk@meta.data[["harmony_cell_subgroup"]]<-paste("T/NK",Idents(tnk),sep="_")
tnk@meta.data[1:5,]
write.csv(tnk@meta.data,file ="3_Figures_Tables/TNK_cell_subgroup.csv")
table(Idents(tnk))


#####GSVA analysis######
library(clusterProfiler)
library(msigdbr)
library(GSVA)
gmt<-read.gmt( "~/Documents/public_data/T_cell.gmt.txt") ###gene list was derived from the Nat Commun. 2017 May 5;8:15081. 
t_tmp<-Idents(tnk)=="CD4-/CD8-"|Idents(tnk)=="T-reg"|Idents(tnk)=="CD8+T"|Idents(tnk)=="CD4+T"|Idents(tnk)=="Proliferating"
table(t_tmp)
t<-tnk[,t_tmp]
gmt$gene[gmt$ont=="Cytotoxic"]
table(gmt$ont)

class(gmt)
gmt_in_use<-list()
gmt_in_use[["Cytotoxic"]]<-gmt$gene[gmt$ont=="Cytotoxic"]
gmt_in_use[["Exhausted"]]<-gmt$gene[gmt$ont=="Exhausted"]
gmt_in_use[["Regulatory"]]<-gmt$gene[gmt$ont=="Regulatory"]
gmt_in_use[["Naive"]]<-gmt$gene[gmt$ont=="Naive"]
gmt_in_use[["Costimulatory"]]<-gmt$gene[gmt$ont=="Costimulatory"]
gmt_in_use[["G1_S"]]<-gmt$gene[gmt$ont=="G1_S"]
gmt_in_use[["G2_M"]]<-gmt$gene[gmt$ont=="G2_M"]



library(Seurat)
exprMat<-GetAssayData(t,slot="data")
exprMat<-GetAssayData(t,slot="counts")
exprMat<-as.matrix(exprMat)
?gsva
res_es <- gsva(exprMat, gmt_in_use, min.sz=1, max.sz=500, verbose=FALSE, parallel.sz=8)
class(res_es)

write.csv(res_es,file="T_cell_matrix.csv")
res_es<-read.csv("./T_cell_matrix.csv",header=T,row.names = 1)
dim(res_es)
Idents(t)

res_es_2<-data.frame(t(res_es))
res_es_2$group<-factor(Idents(t))

ac<-data.frame(cluster=factor(Idents(t)))
table(rownames(res_es_2)==rownames(ac))
library(pheatmap)
t_reg_matrix<-res_es[,ac$cluster=="T-reg"]
p1<-pheatmap(t_reg_matrix,annotation_legend = F,cluster_cols = T,annotation_names_col = F,cluster_rows = F,scale = "column")
t_reg_matrix<-t_reg_matrix[,p1$tree_col$order]

t_cd8_matrix<-res_es[,ac$cluster=="CD8+T"]
p2<-pheatmap(t_cd8_matrix,annotation_legend = F,cluster_cols = T,annotation_names_col = F,cluster_rows = F,scale = "column")
t_cd8_matrix<-t_cd8_matrix[,p2$tree_col$order]
t_cd4_matrix<-res_es[,ac$cluster=="CD4+T"]
p3<-pheatmap(t_cd4_matrix,annotation_legend = F,cluster_cols = T,annotation_names_col = F,cluster_rows = F,scale = "column")
t_cd4_matrix<-t_cd4_matrix[,p3$tree_col$order]

t_pro_matrix<-res_es[,ac$cluster=="Proliferating"]
p4<-pheatmap(t_pro_matrix,annotation_legend = F,cluster_cols = T,annotation_names_col = F,cluster_rows = F,scale = "column")
t_pro_matrix<-t_pro_matrix[,p4$tree_col$order]

t_DN_matrix<-res_es[,ac$cluster=="CD4-/CD8-"]
p5<-pheatmap(t_DN_matrix,annotation_legend = F,cluster_cols = T,annotation_names_col = F,cluster_rows = F,scale = "column")
t_DN_matrix<-t_DN_matrix[,p5$tree_col$order]

total<-cbind(t_cd8_matrix,t_cd4_matrix,t_reg_matrix,t_DN_matrix,t_pro_matrix)
dim(total)

ac<-c(rep("CD8+T",dim(t_cd8_matrix)[2]),rep("CD4+T",dim(t_cd4_matrix)[2]),rep("T-reg",dim(t_reg_matrix)[2]),
           rep("DN_T",dim(t_DN_matrix)[2]), rep("T_pro",dim(t_pro_matrix)[2]))
ac<-as.data.frame(ac)
rownames(ac)<-colnames(total)


pheatmap::pheatmap(total,annotation_col = ac,annotation_legend = T,cluster_cols = F,annotation_names_col = T,cluster_rows = F,scale = "column")

exprMat[1:5,1:5]

exp<-exprMat[,colnames(total)]
#pheatmap::pheatmap(exp,annotation_col = ac,annotation_legend = T,cluster_cols = F,annotation_names_col = T,cluster_rows = F,scale = "column")


gene<-as.vector(gmt$gene[1:30])
gene

exp2<-as.data.frame(exp)[gene,]
pheatmap::pheatmap(exp2,annotation_col = ac,annotation_legend = T,cluster_cols = F,annotation_names_col = T,cluster_rows = F,color = colorRampPalette(c("white", "firebrick3"))(5))


#####proportion of the cell distribution######
table(sce.big$orig.ident)
table(tnk$orig.ident)
table(tnk$orig.ident,Idents(tnk))

tnk_pro<-data.frame(table(sce.big$orig.ident),table(tnk$orig.ident))[,c(1,2,4)]
colnames(tnk_pro)<-c("sampleID","total_cell","tnk_cell")
tnk_pro<-mutate(tnk_pro,total_ratio=tnk_pro$tnk_cell/tnk_pro$total_cell)
tnk_pro$sampleID<-factor(tnk_pro$sampleID,levels=c("BC2","BC3","BC5","BC6","BC16","BC21","BC22","BC10","BC17","BC11","BC20"))


specific<-cbind(table(sce.big$orig.ident),table(tnk$orig.ident,Idents(tnk)))
specific<-data.frame(specific)
specific$sampleID<-rownames(specific)
specific<-mutate(specific,double_negative_ratio=specific$CD4..CD8./specific$V1)
specific$sampleID<-factor(specific$sampleID,levels=c("BC2","BC3","BC5","BC6","BC16","BC21","BC22","BC10","BC17","BC11","BC20"))


library(ggplot2)
ggplot(data=tnk_pro,aes(sampleID,total_ratio))+
  geom_bar(stat = "identity",color="white",fill="navy")+
  theme_light()


ggplot(data=specific,aes(sampleID,double_negative_ratio))+
  geom_bar(stat = "identity",color="white",fill="navy")+
  theme_light()



tab1<-table(Idents(tnk),tnk$orig.ident)
tab2<-prop.table(tab1,2)*100
tab2<-data.frame(tab2)
tab2$Var2<-factor(tab2$Var2,levels=c("BC2","BC3","BC5","BC6","BC16","BC21","BC22","BC10","BC17","BC11","BC20"))
#tab2<-tab2[20:1,]

library(ggplot2)
ggplot(tab2,aes(x=Var2,y=Freq,fill=Var1))+ 
  geom_bar(stat='identity',position='stack',alpha=.5)+ 
  labs(title='Cells_proportion_by_patient',x='',y='')+ 
  theme(legend.justification = 'right', 
        legend.position = 'right', 
        legend.key.height = unit(0.1,'cm'),
        panel.background = element_blank(),
  )+ 
  scale_fill_manual(values=cors)

save.image(file="t&nk0421.RData")
load(file="t&nk0421.RData")
write.csv(tnk@meta.data,file="tandnk_0425.csv")




tab1<-table(Idents(tnk),tnk$orig.ident)

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


