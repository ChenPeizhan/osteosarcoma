####merge the scRNA-data####
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
plan("multiprocess", workers = 6) ###set the compute core
options(future.globals.maxSize = 60000 * 1024^2)

load("2_output/BC10_db_remove.RData") 
load("2_output/BC11_db_remove.RData")  
load("2_output/BC16_db_remove.RData")
load("2_output/BC17_db_remove.RData") 
load("2_output/BC2_db_remove.RData")
load("2_output/BC20_db_remove.RData")
load("2_output/BC21_db_remove.RData")                                                                            
load("2_output/BC22_db_remove.RData")                                                                            
load("2_output/BC3_db_remove.RData") 
load("2_output/BC5_db_remove.RData")                                                                           
load("2_output/BC6_db_remove.RData")

####get the cells to use####
BC10_pass<-BC10_val@meta.data[BC10_val$Doublet=="Singlet",c("Site","Doublet")]
BC11_pass<-BC11_val@meta.data[BC11_val$Doublet=="Singlet",c("Site","Doublet")]
BC16_pass<-BC16_val@meta.data[BC16_val$Doublet=="Singlet",c("Site","Doublet")]
BC17_pass<-BC17_val@meta.data[BC17_val$Doublet=="Singlet",c("Site","Doublet")]
BC2_pass<-BC2_val@meta.data[BC2_val$Doublet=="Singlet",c("Site","Doublet")]
BC20_pass<-BC20_val@meta.data[BC20_val$Doublet=="Singlet",c("Site","Doublet")]
BC21_pass<-BC21_val@meta.data[BC21_val$Doublet=="Singlet",c("Site","Doublet")]
BC22_pass<-BC22_val@meta.data[BC22_val$Doublet=="Singlet",c("Site","Doublet")]
BC3_pass<-BC3_val@meta.data[BC3_val$Doublet=="Singlet",c("Site","Doublet")]
BC5_pass<-BC5_val@meta.data[BC5_val$Doublet=="Singlet",c("Site","Doublet")]
BC6_pass<-BC6_val@meta.data[BC6_val$Doublet=="Singlet",c("Site","Doublet")]


###merge the data###
sce.merge<- merge(BC2_val,y=c(BC3_val,BC5_val,BC6_val,BC16_val,BC21_val,BC22_val,BC11_val,BC20_val,BC10_val,BC17_val),project = "scTOTAL")

cells_to_use<-rbind(BC10_pass,BC11_pass,BC16_pass,BC17_pass,BC2_pass,BC20_pass,
                    BC21_pass,BC22_pass,BC3_pass,BC5_pass,BC6_pass)
rm(list=c("BC2_val","BC3_val","BC5_val","BC6_val","BC16_val","BC21_val","BC22_val","BC11_val","BC20_val","BC10_val","BC17_val"))
dim(cells_to_use)
rm(BC10_pass,BC11_pass,BC16_pass,BC17_pass,BC2_pass,BC20_pass,
                    BC21_pass,BC22_pass,BC3_pass,BC5_pass,BC6_pass)
tmp<-rownames(sce.merge@meta.data)%in%rownames(cells_to_use)
table(tmp)

####save the final data for further analysis####
scedata<-sce.merge[,tmp]
saveRDS(scedata,file="2_output/merge_data_db_remove.rds") ###the final data for next step analysis

###load the data
sce_data<-readRDS(file="2_output/merge_data_db_remove.rds")

####intergration with the harmony method####
library(harmony)
gc()
sce_inter<-NormalizeData(sce_data,verbose = T) 
sce_inter<-FindVariableFeatures(sce_inter,selection.method = "vst", nfeatures = 3000)
sce_inter<-ScaleData(sce_inter,verbose = FALSE)
sce_inter<-RunPCA(sce_inter,verbose = T,npcs = 50)
ElbowPlot(sce_inter,ndims = 100)
p1 <- DimPlot(object = sce_inter, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = sce_inter, features = "PC_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p1,p2))
sce_inter<-RunHarmony(sce_inter,"orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(sce_inter, 'harmony')
dim(harmony_embeddings)
p3 <- DimPlot(object = sce_inter, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p4 <- VlnPlot(object = sce_inter, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p3,p4))


sce_inter <- sce_inter %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  RunTSNE(reduction = "harmony", dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50)

sce_inter<-FindClusters(sce_inter,resolution = 0.5)
table(Idents(sce_inter))
Idents(sce_inter)<-sce_inter$seurat_clusters
sce_inter@meta.data[1:5,]
DimPlot(sce_inter,reduction = "tsne",label = T)
DimPlot(sce_inter,reduction = "umap",label = T)
table(sce_inter@meta.data$RNA_snn_res.0.5)

####define the subgroups#####
sce_inter<-RenameIdents(sce_inter,
                        "0"="1.Osteoblast","3"="1.Osteoblast","6"="1.Osteoblast","15"="1.Osteoblast","17"="1.Osteoblast","21"="1.Osteoblast","22"="1.Osteoblast","25"="1.Osteoblast",
                        "8"="2.Osteoblast_proli","12"="2.Osteoblast_proli",
                        "2"="3.Chondrocyte","20"="3.Chondrocyte",
                        "5"="4.Osteoclast","13"="4.Osteoclast",
                        "7"="5.T/NK",
                        "1"="6.Myeloid","9"="6.Myeloid","19"="6.Myeloid","24"="6.Myeloid","18"="6.Myeloid",
                        "4"="7.Fibroblast","16"="7.Fibroblast","26"="7.Fibroblast",
                        "11"="8.Pericytes",
                        "14"="9.MSC",
                        "23"="10.Myoblast",
                        "10"="11.Endothelial",
                        "27"="6.Myeloid")

markers<-FindAllMarkers(sce_inter,only.pos = T,min.pct = 0.25,test.use = "roc")
markers_2<-FindAllMarkers(sce_inter,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
                          
                        
library(ggsci)
##define the color
cors<-c("#E64B35FF", "#7E6148FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#BB0021FF","#4DBBD5FF", "#B09C85FF","#808180FF", "#1B1919FF")
##Figure 1A##
DimPlot(sce_inter,reduction="tsne",cols=cors)
DimPlot(sce_inter,reduction="umap",cols=cors)
###Figure S1A
DimPlot(sce_inter,reduction="tsne",cols = cors,split.by = "orig.ident",ncol = 3)
DimPlot(sce_inter,reduction="umap",cols = cors,split.by = "orig.ident",ncol = 3)
###Figure 1B

DimPlot(sce_inter,reduction="tsne",cols = cors,split.by = "Site")
DimPlot(sce_inter,reduction="umap",cols = cors,label = F)
DimPlot(sce_inter,reduction="umap",cols = cors,split.by = "orig.ident",ncol = 3)

###Figure S1C##
cell_in_group<-table(Idents(sce_inter),sce_inter@meta.data$orig.ident)
cell_in_group
tab2<-prop.table(cell_in_group,2)
tab2
write.csv(cell_in_group,file="3_Figures_Tables/Figure 1 and Supplementary Figure 1/Harmony_Supplementary_1_cell_in_group.csv")
write.csv(tab2,file="3_Figures_Tables/Figure 1 and Supplementary Figure 1/percent_of_cells.csv")

###Figure 1B####
###Osteoblast and Fibroblast##
VlnPlot(sce_inter,features=c("COL1A1","LUM","CDH11","ACTA2"),cols = cors,pt.size = 0,adjust = 1,ncol = 2)

VlnPlot(sce_inter,features=c("DCN","MME","CXCL14","CXCL12"),cols = cors,pt.size = 0,adjust = 1,ncol = 2)
###T&NK cells
VlnPlot(sce_inter,features=c("CD3D","CD2","RGS5","IGFBP7"),cols = cors,pt.size = 0,adjust = 1,ncol=2)
###Osteoclast cells and endothelial cells###
VlnPlot(sce_inter,features=c("CTSK","MMP9","PECAM1","VWF"),cols = cors,pt.size = 0,adjust = 1,ncol=2)
###myeloid cells
VlnPlot(sce_inter,features=c("FCGR3A","CD14","CD74","HLA-DRA"),cols = cors,pt.size = 0,adjust = 1,ncol=2)

###supplementary Figure
VlnPlot(sce_inter,features=c("PECAM1","DCN","THY1","MYL1"),cols = cors,pt.size = 0,adjust = 1,ncol = 2)
VlnPlot(sce_inter,features=c("IBSP","TOP2A","PCNA","MKI67"),cols = cors,pt.size = 0,adjust = 1,ncol = 2)
VlnPlot(sce_inter,features=c("PTH1R","RUNX2","DCN","LUM"),cols = cors,pt.size = 0,adjust = 1,ncol = 2)
VlnPlot(sce_inter,features=c("CTSK","ACAN","SOX9","COL2A1"),cols = cors,pt.size = 0,adjust = 1,ncol = 2)
VlnPlot(sce_inter,features=c("IL7R","NKG7","GNLY","CD14"),cols = cors,pt.size = 0,adjust = 1,ncol = 2)
VlnPlot(sce_inter,features=c("ACTA2","RGS5","CXCL12","MME"),cols = cors,pt.size = 0,adjust = 1,ncol = 2)
VlnPlot(sce_inter,features=c("SFRP2","MYLPF","MYL1"),cols = cors,pt.size = 0,adjust = 1,ncol = 2)

###Figure 1C##
markers.to.plot <- c("COL1A1","SPP1","LUM","COL2A1","SOX9","ACAN","CTSK","MMP9","IL7R","CD3D","NKG7","CD74","CD14","FCGR3A","RGS5","ACTA2","CXCL12","SFRP2","MME","MYLPF","VWF")
DotPlot(sce_inter, features = rev(markers.to.plot), cols = c("blue","#E64B35FF"), dot.scale = 6,scale = T)
?DotPlot
###Figure 1D###
tab<-table(Idents(sce_inter),sce_inter$orig.ident)
tab2<-prop.table(tab,2)
tab2<-as.data.frame(tab2)
library(ggplot2)
ggplot(tab2,aes(x=Var2,y=Freq,fill=Var1))+ 
  geom_bar(stat='identity',position='stack',alpha=.5)+ 
  labs(title='Bar with Stack',x='',y='')+ 
  theme(legend.justification = 'right', 
        legend.position = 'top', 
        legend.key.height = unit(0.1,'cm'),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank())+ 
  scale_fill_manual(values=cors)


table(Idents(sce_inter))



#####correlation between the cell types######
pheat=list()
for (i in unique(Idents(sce_inter))){
  matrix<-as.matrix(sce_inter@assays$RNA@data[,Idents(sce_inter)==i])
  matrix=rowMeans(matrix)
  pheat[[i]]=matrix
}

pheat_2<-do.call(cbind,pheat)


####correlation of the correlated cells based on the HVGs####
vargenes<-sce_inter@assays$RNA@var.features
pheat_3<-pheat_2[vargenes,]
pheat_3<-pheat_3[,c(6,9,11,2,8,4,3,7,10,1,5)]
ac<-data.frame(cluster=factor(colnames(pheat_3)))
ac
rownames(ac)=colnames(pheat_3)
?pheatmap
ccc<-cors[1:11]
ann_colors = list(cluster = c("1.Osteoblast"="#E64B35FF", "2.Osteoblast_proli"="#7E6148FF", "3.Chondrocyte"="#00A087FF","4.Osteoclast"= "#3C5488FF","5.T/NK"="#F39B7FFF","6.Myeloid"= "#8491B4FF",
 "7.Fibroblast"="#91D1C2FF", "8.Pericytes"="#BB0021FF","9.MSC"= "#4DBBD5FF","10.Myoblast"= "#B09C85FF","11.Endothelial"= "#808180FF"))

pheatmap::pheatmap(corr,fontsize_row = 10,annotation_col = ac,annotation_legend = T,annotation_colors = ann_colors,cluster_rows = T,
                   cluster_cols = T, color = colorRampPalette(c("white", "red"))(10),display_numbers = T)


####correlation of the cells between the samples###
table(sce_inter$orig.ident)
pheat_2=list()
for (i in unique(sce_inter$orig.ident)){
  sample<-sce_inter[,sce_inter$orig.ident==i]
  table(Idents(sample))
  pheat=list()
  for (j in unique(Idents(sce_inter))){
    matrix<-as.matrix(sce_inter@assays$RNA@data[,Idents(sce_inter)==j])
    matrix=rowMeans(matrix)
    pheat[[j]]=matrix}
  pheat_2[[i]]<-do.call(cbind,pheat)
}


name<-colnames(pheat_2[[1]])
name
for (i in c(7)){
  matrix<-sapply(pheat_2, function(v) return(v[,i]))
  pheat_4<-matrix[vargenes,]
  library(corrplot)
  corr<-cor(pheat_4)
  dim(pheat_4)
  pdf(file=paste("t_nk",".pdf",sep=""))
  corrplot(corr,method="color",order="hclust", number.font = 5,diag = TRUE,hclust.method="average",addCoef.col = "blue", col = colorRampPalette(c("blue","white", "red"))(50))
  dev.off()
}





####correlation between the cells of the 

top10 <- ost.markers_total %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DefaultAssay(sce_inter)<-"RNA"
DoHeatmap(sce_inter_2, features = roc,group.colors = cors)


library(dplyr)

man<-sce_inter@meta.data %>% group_by(orig.ident) %>% summarise(nfea.mean=mean(nFeature_RNA))
man<-sce_inter@meta.data %>% group_by(orig.ident) %>% summarise(nfea.mean=median(nFeature_RNA))###from 564 to 4371#
man_2<-sce_inter@meta.data %>% group_by(orig.ident) %>% summarise(nfea.mean=median(nCount_RNA))##1041 to 16249###

library(ggplot2)
dt<-sce_inter@meta.data[,c(1,2,3,4)]
ggplot(dt, aes(orig.ident,nFeature_RNA))+
  geom_violin(aes(fill=orig.ident),trim = F)+
  #geom_boxplot(width=0.2)+
  theme_light()


gg<-ggplot(dt, aes(orig.ident,nCount_RNA))+
  geom_violin(aes(fill=orig.ident),trim = F)+
  #geom_boxplot(width=0.2)+
  theme_light()
library(ggThemeAssist)
ggThemeAssistGadget(gg)

#######save the cell distrubution data##
sce.big<-readRDS(file="2_output/merge_data_db_remove.rds")
table(rownames(sce.big@meta.data)==rownames(sce_inter@meta.data))
table(Idents(sce_inter))
sce_inter@meta.data$Harmony_group=Idents(sce_inter)
sce_inter@meta.data[1:5,]

sce.big@meta.data<-sce_inter@meta.data[,c("orig.ident", "nCount_RNA", "nFeature_RNA",
                                          "percent.mt", "Site","Harmony_group")]
sce.big@meta.data[1:5,]
saveRDS(sce.big,file="2_output/merge_data_db_remove_harmony_group_0427.rds")

rm(osb,chon,chon_clean)
rm(sce.big)

save.image(file="./harmony_integrated_0427.RData")
load("./harmony_integrated_0427.RData")

table(sce_inter$Site)


table(Idents(sce_inter))
tab1<-table(Idents(sce_inter),sce_inter$orig.ident)
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


#rm(list=ls())
