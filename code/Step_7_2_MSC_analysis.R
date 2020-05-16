
rm(list = ls())
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
plan("multiprocess", workers = 4) 
options(future.globals.maxSize = 40000 * 1024^2) 
####test the harmony
###
#load("./harmony_integrated_0413.RData")
#DimPlot(sce_inter,reduction = "tsne")
#DimPlot(sce_inter,reduction = "umap")


library(devtools)
library(harmony)
sce.big<-readRDS(file="2_output/merge_data_db_remove_harmony_group.rds")
table(sce.big@meta.data$Harmony_group)
tmp<-sce.big@meta.data$Harmony_group=="9.MSC"
table(tmp)
msc<-sce.big[,tmp]
rm(sce.big)

msc<-NormalizeData(msc,verbose = FALSE) 
msc<-FindVariableFeatures(msc,selection.method = "vst", nfeatures = 2000)
msc<-ScaleData(msc,vars.to.regress = "percent.mt",verbose = FALSE)
msc<-RunPCA(msc,verbose = T)
p1 <- DimPlot(object = msc, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = msc, features = "PC_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p1,p2))
msc<-RunHarmony(msc,"orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(msc, 'harmony')
dim(harmony_embeddings)
p3 <- DimPlot(object = msc, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p4 <- VlnPlot(object = msc, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p3,p4))


msc <- msc %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  RunTSNE(reduction = "harmony", dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50)


msc<- FindClusters(msc,resolution = 0.8)
DimPlot(msc, reduction = "tsne", pt.size = .1, 
        split.by = 'Site',label = T)
DimPlot(msc, reduction = "tsne", label = TRUE, pt.size = .1)
FeaturePlot(msc,features = c("CD74","CD14","CD3D","MS4A1","VEGFA","JCHAIN","CPA3","HBB","COL2A1","SOX9","COL9A1","IBSP"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(msc,features = c("DCN","THY1","ACTA2","FN1","SOX9","CXCL12","MYH11"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(msc,features = c("THY1","ENG","VCAM1","NT5E","RGS5","WISP2","EGR1"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(msc,features = c("CLEC11A","GREM1","STMN1","SPP1"),reduction = "tsne",max.cutoff =5,label =T)
FeaturePlot(msc,features = c("PTX3"),reduction = "tsne",max.cutoff =5,label =T)


msc<-RenameIdents(msc,"2"="NT5E+MSC","0"="WISP2+MSC","1"="WISP2+MSC","6"="WISP2+MSC","3"="CLEC11A+MSC",
                    "4"="CLEC11A+MSC","5"="CLEC11A+MSC")

library(ggsci)
cors<-pal_d3()(10)
cors

markers<-FindAllMarkers(msc,only.pos = T,min.pct = 0.1,logfc.threshold = 0.1)
markers_total<-FindAllMarkers(msc,only.pos = F,min.pct = 0.1,logfc.threshold = 0.1)

msc@meta.data[1:5,]
DimPlot(msc, reduction = "tsne", pt.size = .1, split.by = 'Site',cols = cors)
DimPlot(msc, reduction = "tsne", pt.size = .1,label = F,cols = cors,split.by = "orig.ident",ncol = 3)
DimPlot(msc, reduction = "tsne", pt.size = .1, split.by = 'Site',cols = cors)
DimPlot(msc, reduction = "tsne", pt.size = .1,label = F,cols = cors)

msc2<-msc
#msc2<-ScaleData(msc2,vars.to.regress = "percent.mt",verbose = FALSE,features = top20$gene)
top20<-markers %>% group_by(cluster) %>% top_n(n=20,wt=avg_logFC)
DoHeatmap(msc2, features = top20$gene,group.colors = cors)
top50<-markers %>% group_by(cluster) %>% top_n(n=50,wt=avg_logFC)
DoHeatmap(msc,top20$gene)
tab2<-table(Idents(msc),msc$Site)
tab2<-prop.table(tab2,2)*100
tab2
tab2<-data.frame(tab2)
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

VlnPlot(msc,features = c("NT5E","WISP2","CLEC11A","B2M"),cols = cors,pt.size = 0,ncol = 2)
VlnPlot(msc,features = c("VEGFA","TGFBI","SFRP2","CXCL14"),cols = cors,pt.size = 0,ncol = 2)



tab1<-table(Idents(msc),msc$orig.ident)
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

save.image(file="MSC0515.RData")



