#rm(list = ls())
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(harmony)
library(monocle)
plan("multiprocess", workers = 4) 
options(future.globals.maxSize = 40000 * 1024^2) 

####Subgroup aaasyas###
sce.big<-readRDS(file="2_output/merge_data_db_remove_harmony_group.rds")
table(sce.big@meta.data$Harmony_group)
tmp<-sce.big@meta.data$Harmony_group=="4.Osteoclast"
table(tmp)
oct<-sce.big[,tmp]
rm(sce.big)
oct<-NormalizeData(oct,verbose = T) 
oct<-FindVariableFeatures(oct,selection.method = "vst", nfeatures = 3000)
oct<-ScaleData(oct,verbose = FALSE)
oct<-RunPCA(oct,verbose = T)
p1 <- DimPlot(object = oct, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = oct, features = "PC_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p1,p2))
oct<-RunHarmony(oct,"orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(oct, 'harmony')
dim(harmony_embeddings)
p3 <- DimPlot(object = oct, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p4 <- VlnPlot(object = oct, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p3,p4))


oct <- oct %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

#tnk<- FindClusters(tnk,resolution = 0.8)
DimPlot(oct, reduction = "tsne", pt.size = .1,label = T,split.by = "Site")
DimPlot(oct, reduction = "tsne", pt.size = .1,label = T)
markers<-FindAllMarkers(oct,only.pos = T)
top10<-markers%>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
####we remove those T, mast cells and Doublets ###317 cells##
tmp<-ifelse((Idents(oct)=="8"|Idents(oct)=="9"|Idents(oct)=="10"|Idents(oct)=="11"),FALSE,TRUE)
table(tmp)
FeaturePlot(oct,features = c("CD3D","CPA3","IFIT1","TOP2A","CD74","CD14"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(oct,features = c("ACP5","CTSK"),reduction = "tsne",max.cutoff =3,label =T)




###re--run the process##
oct<-oct[,tmp]
rm(sce.big)
oct<-NormalizeData(oct,verbose = T) 
oct<-FindVariableFeatures(oct,selection.method = "vst", nfeatures = 3000)
oct<-ScaleData(oct,verbose = FALSE)
oct<-RunPCA(oct,verbose = T)
p1 <- DimPlot(object = oct, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = oct, features = "PC_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p1,p2))
oct<-RunHarmony(oct,"orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(oct, 'harmony')
dim(harmony_embeddings)
p3 <- DimPlot(object = oct, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p4 <- VlnPlot(object = oct, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
CombinePlots(plots=list(p3,p4))


oct <- oct %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

library(RColorBrewer)
FeaturePlot(oct,features = c("CD74","CTSK","TOP2A","CD14","ACP5","MMP9"),reduction = "tsne",max.cutoff =3,label =T)
FeaturePlot(oct,features = c("CTSK","CD14","ACP5","CD74"),reduction = "umap",max.cutoff =8,label =T)
FeaturePlot(oct,features = c("CTSK","CD14","ACP5","CD74","TOP2A","MMP9"),reduction = "tsne",max.cutoff =8,label =F,cols =  brewer.pal(11,"PuOr"))
DimPlot(oct, label = T,reduction = "umap")

Idents(oct)<-oct@meta.data$RNA_snn_res.0.5
oct<-RenameIdents(oct,"2"="OC_progenitor","4"="OC_progenitor","6"="OC_progenitor","9"="OC_progenitor",
                  "3"="OC_immature","1"="OC_immature","7"="OC_immature","8"="OC_immature","10"="OC_immature",
                 "0"="OC_mature", "5"="OC_mature")


markers<-FindAllMarkers(oct,only.pos = T)

library(ggsci)
cors<-pal_d3()(8)
cors
DimPlot(oct, label = F,reduction = "tsne",cols = cors,split.by = "Site")



#####trajectory analysis####
library(monocle)
matrix<-as.matrix(oct@assays$RNA@counts)
dim(matrix)
matrix[1:50,1:10]
gene_ann <- data.frame(
  gene_short_name = row.names(matrix), 
  row.names = row.names(matrix)
)
oct[["cell_group"]]<-Idents(oct)
sample_ann <- oct@meta.data

fd <- new("AnnotatedDataFrame",data=gene_ann)
pd<-new("AnnotatedDataFrame",data=sample_ann)
#?newCellDataSet
sc_cds_2 <- newCellDataSet(matrix,  phenoData = pd,featureData =fd,expressionFamily = negbinomial.size(),lowerDetectionLimit=0.1)
sc_cds_2 <- estimateSizeFactors(sc_cds_2)
sc_cds_2 <- estimateDispersions(sc_cds_2)
####QC#####

sc_cds_2 <- detectGenes(sc_cds_2, min_expr = 3)
expressed_genes <- row.names(subset(fData(sc_cds_2), num_cells_expressed >= 5))

#####cluster analysis######

disp_table <- dispersionTable(sc_cds_2) # 
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1) #
sc_cds_2 <- setOrderingFilter(sc_cds_2, unsup_clustering_genes$gene_id) 
plot_ordering_genes(sc_cds_2) 

plot_pc_variance_explained(sc_cds_2, return_all = F)

sc_cds_2 <- reduceDimension(sc_cds_2, max_components = 2, num_dim = 6,reduction_method = 'tSNE', 
                            residualModelFormulaStr = "~num_genes_expressed+orig.ident",norm_method = 'log',
                            verbose = T)
sc_cds_2 <- clusterCells(sc_cds_2, num_clusters = 5)
p1<-plot_cell_clusters(sc_cds_2, 1, 2, color = "Cluster")
p2<-plot_cell_clusters(sc_cds_2, 1, 2, color = "cell_group")
CombinePlots(list(p1,p2))
plot_cell_clusters(sc_cds_2, 1, 2, color = "cell_group",markers = c("CD74","CD14","TOP2A","CTSK","ACP5","MMP9"))
p1<-DimPlot(oct,reduction = "tsne",label = T)
p2<-plot_cell_clusters(sc_cds_2, 1, 2, color = "cell_group")
CombinePlots(list(p1,p2))
#save.image(file="oct_0318.RData")
#load(file="oct_0318.RData")
diff_test_res <- differentialGeneTest(sc_cds_2[expressed_genes,],
                                      fullModelFormulaStr = "~as.factor(cell_group)",cores = 3,
                                      verbose = T) #+num_genes_expressed+orig.ident
ordering_genes <- row.names(subset(diff_test_res, qval < 1e-2)) #1e-200


#ordering_genes<-markers %>% group_by(cluster) %>% top_n(n=200,wt=avg_logFC)

sc_cds2 <- setOrderingFilter(sc_cds_2, ordering_genes)
plot_ordering_genes(sc_cds2)
sc_cds2<- reduceDimension(sc_cds2, max_components = 2, reduction_method  = "DDRTree")
sc_cds2 <- orderCells(sc_cds2)
sc_cds2 <- orderCells(sc_cds2,root_state = 8)
pData(sc_cds_2)[1:5,]
plot_cell_trajectory(sc_cds2, color_by = "cell_group",show_branch_points = F)+facet_wrap(~cell_group)
plot_cell_trajectory(sc_cds2, color_by = "State",show_branch_points = T)+facet_wrap(~Site)
plot_cell_trajectory(sc_cds2, markers = c("ACP5"),use_color_gradient = T,show_branch_points = F)
plot_cell_trajectory(sc_cds2, markers = c("CTSK"),use_color_gradient = T,show_branch_points = F)
plot_cell_trajectory(sc_cds2, markers = c("CD74"),use_color_gradient = T,show_branch_points = F)
plot_cell_trajectory(sc_cds2, markers = c("TOP2A"),use_color_gradient = T,show_branch_points = F)
plot_cell_trajectory(sc_cds2, markers = c("CD14"),use_color_gradient = T,show_branch_points = F)
plot_cell_trajectory(sc_cds2, color_by = "Pseudotime",show_branch_points = F)


####genes along with the pseudotime###

diff_test_res_pso <- differentialGeneTest(sc_cds2[expressed_genes,],cores = 3,
                                          fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res_pso, qval < 1e-10))
dev.off()
pseudoplot<-plot_pseudotime_heatmap(sc_cds2[sig_gene_names,],
                        cores = 4,num_clusters = 4,
                        show_rownames = F,norm_method = "log",return_heatmap = T)

TF<-read.csv("~/Documents/public_data/Homo_sapiens_TF.txt",sep="\t",header=T)
sig_gene_names_2 <- row.names(subset(diff_test_res_pso, qval < 1e-20))
Tftmp<-sig_gene_names_2%in%TF$Symbol
table(Tftmp)
sig_TFs<-sig_gene_names_2[Tftmp]
dev.off()

pseudoplot_2<-plot_pseudotime_heatmap(sc_cds2[sig_TFs,],
                                    cores = 4,num_clusters = 3,
                                    show_rownames = T,norm_method = "log",return_heatmap = T)


####get the gene names in each cluster##
clusters <- cutree(pseudoplot$tree_row, k = 4)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)


library(msigdbr)
library(clusterProfiler)
msigdbr_show_species()  
m_df <- msigdbr(species = "Homo sapiens") 


m_t2g_C5 <- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "BP") %>% 
  dplyr::select(gs_name, gene_symbol)
?msigdbr
m_t2g_C7 <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, gene_symbol)
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
write.csv(em_group1,file="6.Figure_3/heatmap_cluster1_enrichment.csv")


gene2<-rownames(clustering)[clustering[,1]==2]
em_group2<- enricher(gene2, TERM2GENE=m_t2g_C5)
head(em_group2)
barplot(em_group2,showCategory = 20)
write.csv(em_group2,file="6.Figure_3/heatmap_cluster2_enrichment.csv")

gene3<-rownames(clustering)[clustering[,1]==3]
em_group3<- enricher(gene3, TERM2GENE=m_t2g_C5)
head(em_group3)
barplot(em_group3,showCategory = 20)
write.csv(em_group3,file="6.Figure_3/heatmap_cluster3_enrichment.csv")


gene4<-rownames(clustering)[clustering[,1]==4]
em_group4<- enricher(gene4, TERM2GENE=m_t2g_C5)
head(em_group4)
barplot(em_group4,showCategory = 20)
write.csv(em_group4,file="6.Figure_3/heatmap_cluster4_enrichment.csv")


save.image(file = "./OCT_analysis.RData")
load(file="./OCT_analysis.RData")

#oct@meta.data[["harmony_cell_subgroup"]]<-paste("oct",Idents(oct),sep="_")
#oct@meta.data[1:5,]
#write.csv(oct@meta.data,file ="3_Figures_Tables/osteoclast_oct_subgroup.csv")

table(Idents(oct))

tab1<-table(Idents(oct),oct$orig.ident)
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

write.csv(oct@meta.data,file="OC_sample_0515.csv")
 