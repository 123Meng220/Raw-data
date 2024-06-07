
if (T) {
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("figures")
  dir.create("origin_datas/GEO",recursive = T)
}
library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)

library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(harmony)
options(stringsAsFactors = F)
# # # 
# GSE184880_name=list.files('origin_datas/GEO/GSE184880_RAW/')
# #########
# gsm.list=unique(substr(GSE184880_name,1,10))
# for (f in gsm.list) {
#   dir.create(paste0('origin_datas/GEO/GSE184880_RAW/',f))
# }
# 
# ########
# for (d in gsm.list) {
#   #d=gsm.list[1]
#   files<-list.files(paste0('origin_datas/GEO/GSE184880_RAW/',d))
#   for (f in files){
#     newname<-stringr::str_split_fixed(f,'[.]',2)[,2]
#     file.rename(paste0('origin_datas/GEO/GSE184880_RAW/',d,'/',f),
#                 paste0('origin_datas/GEO/GSE184880_RAW/',d,'/',newname))
#   }
# }

#01.######
dir.create('results/01.cell_annotation')
##############
GSE184880.sample=read.delim('origin_datas/GEO/GSE184880_sample.txt')
GSE184880.sample$Sample_title=gsub('[1234567]','',GSE184880.sample$Sample_title)
head(GSE184880.sample)
GSE184880.sample$tumor_stage=gsub('[ABC123]','',GSE184880.sample$tumor_stage)

dir_name=list.files('origin_datas/GEO/GSE184880_RAW/')
datalist=list()
for (i in 1:length(dir_name)){
  dir.10x = paste0("origin_datas/GEO/GSE184880_RAW/",dir_name[i])
  list.files(dir.10x)
  my.data <- Read10X(data.dir = dir.10x)
  datalist[[i]] = CreateSeuratObject(counts = my.data, project = dir_name[i], min.cells = 3, min.features = 200)
  Samples1=dir_name[i]
  tissue=GSE184880.sample$Sample_title[i]
  stage=GSE184880.sample$tumor_stage[i]
  datalist[[i]] = AddMetaData(datalist[[i]] , Samples1,col.name = "Samples")
  datalist[[i]] = AddMetaData(datalist[[i]] , tissue,col.name = "tissue")
  datalist[[i]] = AddMetaData(datalist[[i]] , stage,col.name = "stage")
}
names(datalist)=dir_name
rm(my.data)

####
for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")# 
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")# 
  datalist[[i]] <- sce
  rm(sce)
}
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
rm(datalist)
#
raw_meta=sce@meta.data
raw_count <- table(raw_meta$Samples)
raw_count
sum(raw_count)#  64659
pearplot_befor<-VlnPlot(sce,group.by ='Samples',
                        features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                        pt.size = 0,
                        ncol = 3)
pearplot_befor
ggsave('results/01.cell_annotation/pearplot_befor.pdf',pearplot_befor,height = 5,width = 10)

#
sce=subset(sce, subset=nFeature_RNA>200 & nFeature_RNA<8000 & percent.mt<15)
clean_meta=sce@meta.data
clean_count <- table(clean_meta$Samples)
clean_count
sum(clean_count)#45448
pearplot_after <- VlnPlot(sce,group.by ='Samples',
                          features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                          pt.size = 0,
                          ncol = 3)
pearplot_after
ggsave('results/01.cell_annotation/pearplot_after.pdf',pearplot_after,height = 6,width = 15)

##1.1 #########
library(harmony)
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(sce, features = rownames(sce))

#
sce <- RunPCA(sce, features = VariableFeatures(sce))
colnames(sce@meta.data)

##
sce = RunHarmony(sce, group.by.vars="Samples", lambda=0.5, max.iter.harmony=30)


pca.plot=ElbowPlot(sce,ndims = 50) #
pdf('results/01.cell_annotation/PCA_plot.pdf',height = 5,width = 5)
pca.plot+theme(text = element_text(family = 'Times',size = 12))
dev.off()
sce <- RunUMAP(sce, dims=1:20, reduction="harmony")

after_batch=DimPlot(sce,group.by='Samples',
                    reduction="umap",
                    label = "T",
                    pt.size = 0.2,
                    label.size = 0)+
  theme(text = element_text(family = 'Times',size = 12))
# after_batch=LabelClusters(after_batch,id = 'Samples',family='Times')
ggsave('results/01.cell_annotation/after_batch.pdf',after_batch,height = 5,width = 6)




##1.2 ####
library(clustree)
sce <- FindNeighbors(sce, dims = 1:20, reduction="harmony")
sce <- FindClusters(object = sce,resolution = .1)
DefaultAssay(sce) <- "RNA"
colnames(sce@meta.data)
length(table(sce@meta.data$seurat_clusters))

sce=subset(sce, subset=seurat_clusters%in% c(0:10))
table(sce@meta.data$seurat_clusters)
p=DimPlot(sce,group.by='seurat_clusters',reduction="umap",label = F,pt.size = 0.2)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(text = element_text(family = 'Times',size = 12),panel.grid = element_blank())
p=LabelClusters(p,id = 'seurat_clusters',family='Times')
p

ggsave('results/01.cell_annotation/cluster_umap.pdf',p,height = 5,width = 5.5)


##1.3 #####
#
Logfc = 0.25
#
Minpct = 0.25
DefaultAssay(sce) <- "RNA"
colnames(sce@meta.data)
Idents(sce)<-'seurat_clusters'

sce.markers <- FindAllMarkers(object = sce,logfc.threshold = Logfc, min.pct = Minpct,only.pos = T)
head(sce.markers)
sce.markers["pct.diff"]=sce.markers$pct.1-sce.markers$pct.2
sce.markers <- sce.markers[sce.markers$p_val_adj<0.05,]
table(sce.markers$cluster)
length(unique(sce.markers$gene))#762
head(sce.markers)
table(sce.markers$cluster)
write.csv(sce.markers,'results/01.cell_annotation/diff_marker_gene.csv')
### 
sce.markers=read.csv('results/01.cell_annotation/diff_marker_gene.csv')
Top5 <- sce.markers %>% group_by(cluster) %>% slice_max(n =5, order_by = avg_logFC)
length(Top5$gene)
length(unique(Top5$gene))
###
diff.marker.dotplot=DotPlot(object = sce, features = unique(Top5$gene),
                            cols=c("snow", "blue"),scale = T)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('')+ylab('')+coord_flip()
diff.marker.dotplot
ggsave('results/01.cell_annotation/cluster_diffmarker.pdf',diff.marker.dotplot,height = 10,width = 10)

VlnPlot(sce, features = c("S100A9",'LYZ'),pt.size = 0,group.by  = 'seurat_clusters')

marker <- data.frame(cluster = 0:10,cell = 0:10)
marker[marker$cluster %in% c(0),2] <- 'Fibroblasts'
marker[marker$cluster %in% c(1,2),2] <- 'T/NK cells'
marker[marker$cluster %in% c(3),2] <- 'Epithelial cells'
marker[marker$cluster %in% c(4),2] <- 'Myeloid cells'
marker[marker$cluster %in% c(5),2] <- 'Ovarian stromal cell'
marker[marker$cluster %in% c(6),2] <- 'Endothelial cell'
marker[marker$cluster %in% c(7),2] <- 'Plasma cell'
marker[marker$cluster %in% c(8),2] <- 'proliferating T cells'
marker[marker$cluster %in% c(9),2] <- 'SMC/myofibroblasts'
marker[marker$cluster %in% c(10),2] <- 'B cells'
marker
sce@meta.data$cell_type <- sapply(sce@meta.data$seurat_clusters,function(x){marker[x,2]})
my.cols=brewer.pal(11,"Set3")[-9]
cell_type_umap=DimPlot(sce,group.by='cell_type',reduction="umap",label = F,pt.size = 0.2,cols =my.cols)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),text = element_text(family = 'Times',size = 12))
cell_type_umap=LabelClusters(cell_type_umap,id = 'cell_type',family='Times')
cell_type_umap
ggsave('results/01.cell_annotation/cell_type_umap.pdf',cell_type_umap,height = 5,width = 5.5)
table(sce$cell_type)


# 0:T cell[NKG7,GZMK,CD3D]
# 1:Hepatocyte[APOA2,ALB]
# 2,3:Macrophage[C1QC,S100A9,LYZ]
# 4:Endothelial cell[FCN3,CLDN5]
# 5:Cycling cells[MKI67,TOP2A]
# 6:Fibroblast[ACTA2,TAGLN,DCN]
# 7:Epithelial cell[KRT8,KRT18]
# 8:B cell[CD79A,MS4A1]
# 9:Plasma cell[IGHA1,IGLC3]

marker_gene=c( 'GNLY','NKG7','CD3D','CD8A',
               'COL1A1','DCN',
               'KRT8','EPCAM',
               'STAR', 'LYZ','C1QA','CD14',
               'VWF','CLDN5','PECAM1',
               'ACTA2','TAGLN','MYH11',
               'MKI67','TOP2A' ,
               'IGLC2','JCHAIN',
               'CD79A','MS4A1'
)
Idents(sce)='cell_type'
marker.dot=DotPlot(object = sce, features = marker_gene,
                   col.min = .5,cols=c("white", "#1B9E77"),scale = T)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size=12),text = element_text(family = 'Times',size = 12)) +
  xlab('')+ylab('')+coord_flip()
marker.dot
ggsave('results/01.cell_annotation/marker_dotplot.pdf',marker.dot,height = 6,width = 6.5)
table(sce$cell_type)

saveRDS(sce,file = 'results/01.cell_annotation/sce.rds')

VlnPlot(sce, features = c('CD2','NKG7','CD3D','CD8A'),pt.size = 0,group.by  = 'seurat_clusters')

##1.4 #####
#####
cell_freq2=melt(prop.table(table(sce$cell_type,sce$tissue),margin=2))
cell_freq2
colnames(cell_freq2)<-c('cell_type','Tissue','Freq')
cell_prop_fig2=ggplot(cell_freq2,aes(x=reorder(cell_type,Freq),y=Freq,fill=Tissue))+
  scale_fill_manual(values = pal_npg()(10))+
  geom_bar(stat="identity", position="dodge")+
  xlab('')+ylab('Proportion')+theme_classic()+
  theme(text = element_text(family = 'Times',size=12),
        axis.text.x = element_text(angle = 30,hjust = 1),
        legend.position = 'top',legend.text =element_text(family = 'Times'),
        legend.title = element_text(family = 'Times'))
cell_prop_fig2
            


cell_prop_fig1=ggplot(cell_freq2[cell_freq2$Tissue=='Cancer_HGSOC',],aes(x="", y=Freq, fill=cell_type)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  scale_fill_manual(values = my.cols)+ 
  geom_text(aes(label = paste0(round(Freq/sum(Freq)*100, 1), "%")),
            position = position_stack(vjust = 0.5))+
  theme_classic()+ggtitle('Cancer_HGSOC')

ggplot(cell_freq2[cell_freq2$Tissue=='Normal ovarian tissue',],aes(x="", y=Freq, fill=cell_type)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +scale_fill_igv()+
  geom_text(aes(label = paste0(round(Freq/sum(Freq)*100, 1), "%")),
            position = position_stack(vjust = 0.5))+
  theme_classic()+ggtitle('Normal ovarian tissue')



# #####
cell_freq3=melt(prop.table(table(sce$cell_type,sce$stage),margin=2))
cell_freq3
colnames(cell_freq3)<-c('cell_type','stage','Freq')
cell_prop_fig3=ggplot(cell_freq3,aes(x=reorder(cell_type,Freq),y=Freq,fill=stage))+
  scale_fill_manual(values = c("#FF8888FF","#D4D4FFFF","#8787FFFF","#4500ACFF"))+
  geom_bar(stat="identity", position="dodge")+
  xlab('')+ylab('Proportion')+theme_classic()+
  theme(text = element_text(family = 'Times',size=12),
        axis.text.x = element_text(angle = 30,hjust = 1),
        legend.position = 'top',legend.text =element_text(family = 'Times'),
        legend.title = element_text(family = 'Times'))
cell_prop_fig3


pdf('results/01.cell_annotation/cell_prop.pdf',height = 5,width = 12,onefile = F)
mg_merge_plot(cell_prop_fig2,cell_prop_fig1,widths = c(1,1.2))
dev.off()


pdf('results/01.cell_annotation/Fig1.pdf',height = 12,width = 12,onefile = F)
mg_merge_plot(mg_merge_plot(cell_type_umap,marker.dot,labels = c('A','B')),
              mg_merge_plot(cell_prop_fig2,cell_prop_fig1,widths = c(1,1.2),labels = c('C','D')),
              nrow=2,heights = c(1.2,1))
dev.off()



#02.############## 
dir.create('results/02.T_NK_cell')
# Tcell=subset(sce,subset=cell_type %in% 'T/NK cells' & tissue=='Cancer_HGSOC')
# Tcell <- NormalizeData(Tcell)
# Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
# Tcell <- ScaleData(Tcell, features = rownames(Tcell))
# Tcell <- RunPCA(Tcell, features = VariableFeatures(Tcell))
# colnames(Tcell@meta.data)
# ElbowPlot(Tcell,ndims = 50)
# dev.off()
# ###
# Tcell <- RunUMAP(Tcell, dims=1:20, reduction="harmony")
# Tcell <- FindNeighbors(Tcell, dims = 1:20, reduction="harmony")
# saveRDS(Tcell,file = 'results/02.T_NK_cell/Tcell.rds')

Tcell=readRDS('results/02.T_NK_cell/Tcell.rds')
Tcell <- FindClusters(object = Tcell,resolution = .4)
DefaultAssay(Tcell) <- "RNA"
colnames(Tcell@meta.data)
table(Tcell@meta.data$seurat_clusters)
Tcell=subset(Tcell,seurat_clusters%in%c(0,1,2,3,5))
DimPlot(Tcell,group.by='seurat_clusters',reduction="umap",label = T,pt.size = 0.4)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),text = element_text(family = 'Times',size = 12))


Idents(Tcell)<-'seurat_clusters'
Tcell.sce.markers <- FindAllMarkers(object = Tcell,logfc.threshold = 0.25, min.pct = 0.25,only.pos = T)
head(Tcell.sce.markers)
Tcell.sce.markers["pct.diff"]=Tcell.sce.markers$pct.1-Tcell.sce.markers$pct.2
Tcell.sce.markers <- Tcell.sce.markers[Tcell.sce.markers$p_val_adj<0.05,]
table(Tcell.sce.markers$cluster)
length(unique(Tcell.sce.markers$gene))#762
head(Tcell.sce.markers)
table(Tcell.sce.markers$cluster)
write.csv(Tcell.sce.markers,'results/02.T_NK_cell/Tcell_diff_gene.csv')

##2.1 ########
marker <- data.frame(cluster = 0:5,cell = 0:5)
marker[marker$cluster %in% c(0,1,5),2] <- 'cytotoxic NK/T cells'
marker[marker$cluster %in% c(2),2] <- 'naive T cells'
marker[marker$cluster %in% c(3),2] <- 'Tregs'
marker[marker$cluster %in% c(4),2] <- ''

marker
Tcell@meta.data$subcluster <- sapply(Tcell@meta.data$seurat_clusters,function(x){marker[x,2]})

cell_type_umap2=DimPlot(Tcell,group.by='subcluster',reduction="umap",label = F,pt.size = 1,cols = c("#E7298A","#66A61E","#E6AB02"))+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),text = element_text(family = 'Times',size = 12))
cell_type_umap2=LabelClusters(cell_type_umap2,id = 'subcluster',family='Times')
cell_type_umap2
ggsave('results/02.T_NK_cell/T_cell_type_umap.pdf',cell_type_umap2,height = 5,width = 7)



# cytotoxic  :0,1,5,
p=VlnPlot(Tcell, features = c('CD8A','CD8B','GZMK','NKG7'),pt.size = 0,group.by  = 'subcluster',ncol = 4)
ggsave('results/02.T_NK_cell/vlnplot_1.pdf',p,height = 3,width = 12)
#naive 2
p=VlnPlot(Tcell, features = c('CCR7','TCF7'),pt.size = 0,group.by  = 'subcluster')
ggsave('results/02.T_NK_cell/vlnplot_2.pdf',p,height = 3,width = 5)
#3:Regulatory T(Treg) cell	FOXP3,IL2RA
p=VlnPlot(Tcell, features = c('FOXP3','IL2RA','RTKN2'),pt.size = 0,group.by  = 'subcluster')
ggsave('results/02.T_NK_cell/vlnplot_3.pdf',p,height = 3,width = 7)

p=VlnPlot(Tcell, features = c('CD8A','GZMK','NKG7','CCR7','TCF7','FOXP3','IL2RA','RTKN2'),
          pt.size = 0,group.by  = 'subcluster',ncol = 8)
ggsave('results/02.T_NK_cell/vlnplot_merge.pdf',p,height = 4,width = 15)


##2.2 ######
bar = Tcell@meta.data %>% group_by(stage, subcluster) %>% count()
Type_label = c('Stage I',  'Stage II', 'Stage III')
bar$stage = factor(bar$stage, levels=Type_label)
bar = bar %>% group_by(stage) %>% mutate(percent=100*n/sum(n))

p=ggplot(data=bar, aes(x=subcluster, y=percent, fill=stage,label = sprintf("%.2f", percent)))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c("#FDB462","#B3DE69","#FCCDE5"))+theme_classic()+
  ggtitle("Percent(%)")+xlab('')+
  geom_text(position = position_dodge(width = 0.9), vjust = -0.5, size = 4)+
  theme(axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"),
        legend.text=element_text(family = 'Times', size=12), 
        legend.title=element_blank(), text = element_text(family = 'Times'))
ggsave('results/02.T_NK_cell/cell_barplot.pdf',p,height = 5,width = 8)


# ####
# Idents(Tcell)<-'subcluster'
# Tcell.sce.markers <- FindAllMarkers(object = Tcell,logfc.threshold = 0.5, min.pct = 0.25,only.pos = T)
# head(Tcell.sce.markers)
# Tcell.sce.markers["pct.diff"]=Tcell.sce.markers$pct.1-Tcell.sce.markers$pct.2
# Tcell.sce.markers <- Tcell.sce.markers[Tcell.sce.markers$p_val_adj<0.05,]
# table(Tcell.sce.markers$cluster)
# length(unique(Tcell.sce.markers$gene))#762
# head(Tcell.sce.markers)
# table(Tcell.sce.markers$cluster)
# write.csv(Tcell.sce.markers,'results/02.T_NK_cell/T_subcluster_diff_gene.csv')



##2.3 ########
countexp<-Tcell@assays$RNA@counts
countexp<-data.frame(as.matrix(countexp))
library(AUCell)
library(GSEABase)
cells_rankings <- AUCell_buildRankings(as.matrix(countexp)) #rank
####
file_paths=list.files('origin_datas/immune_gmt/')
read_and_process <- function(file_path) {
  df <- read.gmt(paste0('origin_datas/immune_gmt/',file_path))
  return(df)
}
list_of_dfs <- lapply(file_paths, read_and_process)
list_of_immune=sapply(list_of_dfs, function(x){subset(x,select='gene',drop=TRUE)})
names(list_of_immune)=str_split_fixed(file_paths,'[.]',2)[,1]
geneSets <- list_of_immune
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) #calc


signature_exp <- data.frame(getAUC(cells_AUC))
Tcell@assays$HALLMARK$score<-signature_exp
Tcell$group=paste0(Tcell$subcluster,'_',Tcell$stage)

df = as.data.frame(t(Tcell@assays[["HALLMARK"]][["score"]]))
rownames(df) <- gsub(".", "-", rownames(df), fixed = TRUE)
df = df[names(Tcell$group),]
df[1:5,1:5]


df$group <- Tcell$group
avg_df =aggregate(df[,1:ncol(df)-1],list(df$group),mean)
rownames(avg_df) = avg_df$Group.1
avg_df=avg_df[,-1]
avg_df <- as.data.frame(t(avg_df))
rownames(avg_df)
head(avg_df)


pdf('results/02.T_NK_cell/immune_pathway_heatmap.pdf',height = 7,width = 8,onefile = F)
pheatmap(avg_df[,1:3], 
         show_colnames = T,scale='row',name = 'AUCell score',
         color=colorRampPalette(c('#1A5592','white',"orange"))(100),
         cluster_cols = F, cluster_rows = T,
         fontsize_row = 14,fontsize_col = 14)
dev.off()



# ####
# 
table(Tcell$subcluster)
cytotoxic_T=subset(Tcell,subcluster=='cytotoxic NK/T cells')
DotPlot(object = cytotoxic_T, features =c('GZMA','GZMB','GZMM','GZMH','GZMK','PRF1'),
        cols=c("skyblue", "orange"),scale = T,group.by ='stage')+
  RotatedAxis()+ xlab('')+ylab('')+#coord_flip()+
  theme(plot.title = element_text(hjust = 0.5))
ggsave('results/02.T_NK_cell/cytotoxic_gene_dotplot.pdf',height = 4.5,width = 6)




#03.###########
dir.create('results/03.Epithelial')
# sce=readRDS('results/01.cell_annotation/sce.rds')
# table(sce$cell_type,sce$tissue)
# Epithelial=subset(sce,subset=cell_type %in% 'Epithelial cells' & tissue %in% 'Cancer_HGSOC')
# Epithelial <- NormalizeData(Epithelial)
# Epithelial <- FindVariableFeatures(Epithelial, selection.method = "vst", nfeatures = 2000)
# Epithelial <- ScaleData(Epithelial, features = rownames(Epithelial))
# Epithelial <- RunPCA(Epithelial, features = VariableFeatures(Epithelial))
# colnames(Epithelial@meta.data)
# ElbowPlot(Epithelial,ndims = 50)
# dev.off()
# ###降维聚类
# Epithelial <- RunUMAP(Epithelial, dims=1:20, reduction="harmony")
# Epithelial <- FindNeighbors(Epithelial, dims = 1:20, reduction="harmony")
# saveRDS(Epithelial,file = 'results/03.Epithelial/Epithelial.rds')


Epithelial=readRDS('results/03.Epithelial/Epithelial.rds')
Epithelial <- FindClusters(object = Epithelial,resolution = .2)
DefaultAssay(Epithelial) <- "RNA"
colnames(Epithelial@meta.data)
length(table(Epithelial@meta.data$seurat_clusters))

Epithelial=subset(Epithelial,seurat_clusters%in%c(0,2,3))
table(Epithelial$seurat_clusters)
p=DimPlot(Epithelial,group.by='seurat_clusters',reduction="umap",label = F,pt.size = 0.2)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(text = element_text(family = 'Times',size = 12),panel.grid = element_blank())
p=LabelClusters(p,id = 'seurat_clusters',family='Times')
p


Idents(Epithelial)<-'seurat_clusters'
Epithelial.sce.markers <- FindAllMarkers(object = Epithelial,logfc.threshold = 0.25, min.pct = 0.25,only.pos = T)
head(Epithelial.sce.markers)
Epithelial.sce.markers["pct.diff"]=Epithelial.sce.markers$pct.1-Epithelial.sce.markers$pct.2
Epithelial.sce.markers <- Epithelial.sce.markers[Epithelial.sce.markers$p_val_adj<0.05,]
table(Epithelial.sce.markers$cluster)
length(unique(Epithelial.sce.markers$gene))#762
head(Epithelial.sce.markers)
table(Epithelial.sce.markers$cluster)
write.csv(Epithelial.sce.markers,'results/03.Epithelial/Epithelial_diff_gene.csv')

Top5 <- Epithelial.sce.markers %>% group_by(cluster) %>% slice_max(n =5, order_by = avg_logFC)
length(Top5$gene)
length(unique(Top5$gene))
###
# pdf('results/03.Epithelial/Epithelial_diff_gene.pdf',height = 7,width = 5,onefile = F)
DotPlot(object = Epithelial, features = unique(Top5$gene), cols=c("snow", "blue"),scale = T )+
  RotatedAxis()+ xlab('')+ylab('')+coord_flip()+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


##3.1 ####
marker <- data.frame(cluster = c(0:3),cell =c(0:3))
marker[marker$cluster %in% c(0),2] <- 'C1'
marker[marker$cluster %in% c(1),2] <- ' '
marker[marker$cluster %in% c(2),2] <- 'C2'
marker[marker$cluster %in% c(3),2] <- 'C3'
marker
Epithelial@meta.data$subcluster <- sapply(Epithelial@meta.data$seurat_clusters,function(x){marker[x,2]})


cell_type_umap3=DimPlot(Epithelial,group.by='subcluster',reduction="umap",label = F,pt.size = 0.4)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),text = element_text(family = 'Times',size = 12))
cell_type_umap3=LabelClusters(cell_type_umap3,id = 'subcluster',family='Times')
cell_type_umap3
ggsave('results/03.Epithelial/subcluster_umap.pdf',cell_type_umap3,height = 5,width = 6)


##3.2 ####
bar = Epithelial@meta.data %>% group_by(stage, subcluster) %>% count()
Type_label = c('Stage I',  'Stage II', 'Stage III')
bar$stage = factor(bar$stage, levels=Type_label)
bar = bar %>% group_by(stage) %>% mutate(percent=100*n/sum(n))

ggplot(data=bar, aes(x=subcluster, y=percent, fill=stage,label = sprintf("%.2f", percent)))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c("#FDB462","#B3DE69","#FCCDE5"))+theme_classic()+
  ggtitle("Percent(%)")+
  geom_text(position = position_dodge(width = 0.9), vjust = -0.5, size = 4)+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank())
ggsave('results/03.Epithelial/subcluster_barplot.pdf',height = 5,width = 8)

##3.3 ####
Idents(Epithelial)<-'subcluster'
Epithelial.sce.markers <- FindAllMarkers(object = Epithelial,logfc.threshold = 0.5, min.pct = 0.25,only.pos = T)
head(Epithelial.sce.markers)
Epithelial.sce.markers["pct.diff"]=Epithelial.sce.markers$pct.1-Epithelial.sce.markers$pct.2
Epithelial.sce.markers <- Epithelial.sce.markers[Epithelial.sce.markers$p_val_adj<0.05,]
table(Epithelial.sce.markers$cluster)
length(unique(Epithelial.sce.markers$gene))#762
head(Epithelial.sce.markers)
table(Epithelial.sce.markers$cluster)
write.csv(Epithelial.sce.markers,'results/03.Epithelial/Epithelial_subcluster_degs.csv')



Epithelial.marker.genesets.list=split(x=Epithelial.sce.markers,f=Epithelial.sce.markers$cluster)
Epithelial.marker.genesets.list=sapply(Epithelial.marker.genesets.list, function(x){subset(x,select='gene',drop=TRUE)})
library(topGO)
library(clusterProfiler)
library(Rgraphviz)
library(org.Hs.eg.db)
entrez_id = mapIds(x = org.Hs.eg.db, keys = Epithelial.marker.genesets.list$Epithelial_C1,
                   keytype = "SYMBOL",column = "ENTREZID")
entrez_id = na.omit(entrez_id)
c1.erich.go.BP = enrichGO(gene = entrez_id,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05)
c1.erich.go.BP.res=c1.erich.go.BP@result
write.xlsx(c1.erich.go.BP.res,'results/03.Epithelial/Epithelial_c1_enrichment.xlsx',overwrite = T)

c1.erich.go.BP.res=c1.erich.go.BP.res[c1.erich.go.BP.res$p.adjust<0.05,c('Description','GeneRatio','Count','p.adjust')]
c1.erich.go.BP.res=c1.erich.go.BP.res %>% slice_min(n =10, order_by = p.adjust) 
head(c1.erich.go.BP.res)
p1=ggplot(data=c1.erich.go.BP.res,aes(x=Count,y=reorder(Description,Count), color = -log10(p.adjust))) +
  geom_point(aes(size=Count),show.legend = T) +
  scale_color_gradient(low = "blue", high = "orange")+
  geom_segment(aes(yend=Description,xend=0),size=.5) +
  labs(x='Count',y='')+theme_bw()+
  theme(text = element_text(family = 'Times',size = 14))+
  scale_y_discrete(labels=function(y)str_wrap(y,width = 35))
p1


entrez_id = mapIds(x = org.Hs.eg.db, keys = Epithelial.marker.genesets.list$Epithelial_C2,
                   keytype = "SYMBOL",column = "ENTREZID")
entrez_id = na.omit(entrez_id)
c2.erich.go.BP = enrichGO(gene = entrez_id,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05)
c2.erich.go.BP.res=c2.erich.go.BP@result
write.xlsx(c2.erich.go.BP.res,'results/03.Epithelial/Epithelial_c2_enrichment.xlsx',overwrite = T)

c2.erich.go.BP.res=c2.erich.go.BP.res[c2.erich.go.BP.res$p.adjust<0.05,c('Description','GeneRatio','Count','p.adjust')]
c2.erich.go.BP.res=c2.erich.go.BP.res %>% slice_min(n =10, order_by = p.adjust) 
head(c2.erich.go.BP.res)
p2=ggplot(data=c2.erich.go.BP.res,aes(x=Count,y=reorder(Description,Count), color = -log10(p.adjust))) +
  geom_point(aes(size=Count),show.legend = T) +
  scale_color_gradient(low = "blue", high = "orange")+
  geom_segment(aes(yend=Description,xend=0),size=.5) +
  labs(x='Count',y='')+theme_bw()+
  theme(text = element_text(family = 'Times',size = 14))+
  scale_y_discrete(labels=function(y)str_wrap(y,width = 35))
p2


entrez_id = mapIds(x = org.Hs.eg.db, keys = Epithelial.marker.genesets.list$Epithelial_C3,
                   keytype = "SYMBOL",column = "ENTREZID")
entrez_id = na.omit(entrez_id)
c3.erich.go.BP = enrichGO(gene = entrez_id,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05)
c3.erich.go.BP.res=c3.erich.go.BP@result
write.xlsx(c3.erich.go.BP.res,'results/03.Epithelial/Epithelial_c3_enrichment.xlsx',overwrite = T)

c3.erich.go.BP.res=c3.erich.go.BP.res[c3.erich.go.BP.res$p.adjust<0.05,c('Description','GeneRatio','Count','p.adjust')]
c3.erich.go.BP.res=c3.erich.go.BP.res %>% slice_min(n =10, order_by = p.adjust) 
head(c3.erich.go.BP.res)
p3=ggplot(data=c3.erich.go.BP.res,aes(x=Count,y=reorder(Description,Count), color = -log10(p.adjust))) +
  geom_point(aes(size=Count),show.legend = T) +
  scale_color_gradient(low = "blue", high = "orange")+
  geom_segment(aes(yend=Description,xend=0),size=.5) +
  labs(x='Count',y='')+theme_bw()+
  theme(text = element_text(family = 'Times',size = 14))+
  scale_y_discrete(labels=function(y)str_wrap(y,width = 35))
p3

p=mg_merge_plot(p1,p2,p3,nrow = 3,labels = LETTERS[3:5])
ggsave('results/03.Epithelial/subcluster_GOBP_enrichment.pdf',p,height = 16,width = 8)


##3.4 ###########
countexp<-Epithelial@assays$RNA@counts
countexp<-data.frame(as.matrix(countexp))
gmtFile<-"h.all.v2023.1.Hs.symbols.gmt"
library(AUCell)
library(GSEABase)
cells_rankings <- AUCell_buildRankings(as.matrix(countexp)) #rank
geneSets <- getGmt(gmtFile) #signature read
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) #calc
signature_exp <- data.frame(getAUC(cells_AUC))
Epithelial@assays$HALLMARK$score<-signature_exp

df = as.data.frame(t(Epithelial@assays[["HALLMARK"]][["score"]]))
rownames(df) <- gsub(".", "-", rownames(df), fixed = TRUE)
df = df[names(Epithelial$subcluster),]
df[1:5,1:5]
colnames(df)=gsub('HALLMARK_','',colnames(df))
# colnames(df)=gsub('_',' ',colnames(df))
Epithelial.meta=Epithelial@meta.data
Epithelial.meta=Epithelial.meta[,c('stage','subcluster')]
head(Epithelial.meta)

Epithelial.meta$group=paste0(Epithelial.meta$subcluster,' ',Epithelial.meta$stage)
df$group <- Epithelial.meta$group
avg_df =aggregate(df[,1:ncol(df)-1],list(df$group),mean)
rownames(avg_df) = avg_df$Group.1
avg_df=avg_df[,-1]
avg_df <- as.data.frame(t(avg_df))
rownames(avg_df)
head(avg_df)
cell.anno=data.frame(cell=rep(c('Epithelial_C1','Epithelial_C2','Epithelial_C3'),c(3,3,3)))
rownames(cell.anno)=colnames(avg_df)

 pdf('results/03.Epithelial/hallmark_heatmap_ALL.pdf',height = 15,width = 10,onefile = F)
pheatmap(avg_df[,rownames(cell.anno)], show_colnames = T,scale='row',
         color=colorRampPalette(c('#1A5592','white',"orange"))(100),
         cluster_cols = F, cluster_rows = T,
         annotation_col  = cell.anno,
         gaps_col = c(3,6),#gaps_row = c(3,9,12,19,26,31,37),
         fontsize_row = 10,fontsize_col = 12)
dev.off()




####
df_c1=df[rownames(Epithelial.meta[Epithelial.meta$subcluster=='Epithelial_C1',]),]
df_c1$stage <- Epithelial.meta$stage[Epithelial.meta$subcluster=='Epithelial_C1']
avg_df =aggregate(df_c1[,1:ncol(df_c1)-1],list(df_c1$stage),mean)
rownames(avg_df) = avg_df$Group.1
avg_df=avg_df[,-1]
avg_df <- as.data.frame(t(avg_df))
rownames(avg_df)
head(avg_df)

pdf('results/03.Epithelial/hallmark_heatmap_c1.pdf',height = 15,width = 6,onefile = F)
pheatmap(avg_df, show_colnames = T,scale='row',
         main = 'Epithelial_C1',name = 'AUCell',
         color=colorRampPalette(c('#1A5592','snow',"#FB8072"))(100),
         cluster_cols = F, cluster_rows = T,
         fontsize_row = 10,fontsize_col = 12)
dev.off()

####
df_c2=df[rownames(Epithelial.meta[Epithelial.meta$subcluster=='Epithelial_C2',]),]
df_c2$stage <- Epithelial.meta$stage[Epithelial.meta$subcluster=='Epithelial_C2']
avg_df =aggregate(df_c2[,1:ncol(df_c2)-1],list(df_c2$stage),mean)
rownames(avg_df) = avg_df$Group.1
avg_df=avg_df[,-1]
avg_df <- as.data.frame(t(avg_df))
rownames(avg_df)
head(avg_df)
pdf('results/03.Epithelial/hallmark_heatmap_c2.pdf',height = 15,width = 6,onefile = F)
pheatmap(avg_df, show_colnames = T,scale='row',
         main = 'Epithelial_C2',name = 'AUCell',
         color=colorRampPalette(c('#1A5592','snow',"#FB8072"))(100),
         cluster_cols = F, cluster_rows = T,
         fontsize_row = 10,fontsize_col = 12)
dev.off()

####
df_c3=df[rownames(Epithelial.meta[Epithelial.meta$subcluster=='Epithelial_C3',]),]
df_c3$stage <- Epithelial.meta$stage[Epithelial.meta$subcluster=='Epithelial_C3']
avg_df =aggregate(df_c3[,1:ncol(df_c3)-1],list(df_c3$stage),mean)
rownames(avg_df) = avg_df$Group.1
avg_df=avg_df[,-1]
avg_df <- as.data.frame(t(avg_df))
rownames(avg_df)
head(avg_df)

pdf('results/03.Epithelial/hallmark_heatmap_c3.pdf',height = 15,width = 6,onefile = F)
pheatmap(avg_df, show_colnames = T,scale='row',
         main = 'Epithelial_C3',name = 'AUCell',
         color=colorRampPalette(c('#1A5592','snow',"#FB8072"))(100),
         cluster_cols = F, cluster_rows = T,
         fontsize_row = 10,fontsize_col = 12)
dev.off()




#04.##############
dir.create('results/04.Myeloid')
# sce=readRDS('results/01.cell_annotation/sce.rds')
# table(sce$cell_type,sce$tissue)
# Myeloid=subset(sce,subset=cell_type %in% 'Myeloid cells' & tissue %in% 'Cancer_HGSOC')
# Myeloid <- NormalizeData(Myeloid)
# Myeloid <- FindVariableFeatures(Myeloid, selection.method = "vst", nfeatures = 2000)
# Myeloid <- ScaleData(Myeloid, features = rownames(Myeloid))
# Myeloid <- RunPCA(Myeloid, features = VariableFeatures(Myeloid))
# colnames(Myeloid@meta.data)
# ElbowPlot(Myeloid,ndims = 50)
# dev.off()
# ###
# Myeloid <- RunUMAP(Myeloid, dims=1:20, reduction="harmony")
# Myeloid <- FindNeighbors(Myeloid, dims = 1:20, reduction="harmony")
# saveRDS(Myeloid,file = 'results/04.Myeloid/Myeloid.rds')


Myeloid=readRDS('results/04.Myeloid/Myeloid.rds')
Myeloid <- FindClusters(object = Myeloid,resolution = .4)
DefaultAssay(Myeloid) <- "RNA"
colnames(Myeloid@meta.data)
table(Myeloid@meta.data$seurat_clusters)
length(table(Myeloid@meta.data$seurat_clusters))
Myeloid=subset(Myeloid,seurat_clusters%in%c(0,1,4))
table(Myeloid$seurat_clusters)
p=DimPlot(Myeloid,group.by='seurat_clusters',reduction="umap",label = F,pt.size = 0.2)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(text = element_text(family = 'Times',size = 12),panel.grid = element_blank())
p=LabelClusters(p,id = 'seurat_clusters',family='Times')
p



DotPlot(object = Myeloid, features = c('FCN1','S100A8',
                                       'CD68','CD163',
                                       'CD86'),
        cols=c("white", "#1B9E77"),scale = T,col.min = .5)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size=12),
        text = element_text(family = 'Times',size = 12)) +
  xlab('')+ylab('')#+coord_flip()



##4.1 ####


marker <- data.frame(cluster = c(0:5),cell =c(0:5))
marker[marker$cluster %in% c(0),2] <- 'M1 Macrophage'
marker[marker$cluster %in% c(1),2] <- 'Monocyte'
marker[marker$cluster %in% c(2,3,5),2] <- ' '
marker[marker$cluster %in% c(4),2] <- 'M0 Macrophage'
marker
Myeloid@meta.data$subcluster <- sapply(Myeloid@meta.data$seurat_clusters,function(x){marker[x,2]})
table(Myeloid@meta.data$subcluster)

Myeloid$seurat_clusters=paste0('Myeloid_C',Myeloid$seurat_clusters)
cell_type_umap3=DimPlot(Myeloid,group.by='subcluster',reduction="umap",label = F,pt.size = 0.4,cols = c("#8DD3C7","#D95F02" ,"#E7298A"))+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),text = element_text(family = 'Times',size = 12))
cell_type_umap3=LabelClusters(cell_type_umap3,id = 'subcluster',family='Times')
cell_type_umap3
ggsave('results/04.Myeloid/subcluster_umap.pdf',cell_type_umap3,height = 4,width = 5)

DotPlot(object = Myeloid, features = c('CD68','CD163','CD86','FCN1','S100A8'),
        cols=c("white", "#1B9E77"),scale = T,group.by = 'subcluster',col.min = 0)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size=12),
        text = element_text(family = 'Times',size = 12)) +
  xlab('')+ylab('')+coord_flip()
ggsave('results/04.Myeloid/subcluster_marker.pdf',height = 4,width = 5)

##4.2 ####
bar = Myeloid@meta.data %>% group_by(stage, subcluster) %>% count()
bar=bar[bar$subcluster!=' ',]
Type_label = c('Stage I',  'Stage II', 'Stage III')
bar$stage = factor(bar$stage, levels=Type_label)
bar = bar %>% group_by(stage) %>% mutate(percent=100*n/sum(n))

ggplot(data=bar, aes(x=subcluster, y=percent, fill=stage,label = sprintf("%.2f", percent)))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c("#FDB462","#B3DE69","#FCCDE5"))+theme_classic()+
  ggtitle("Percent(%)")+
  geom_text(position = position_dodge(width = 0.9), vjust = -0.5, size = 4)+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank())
ggsave('results/04.Myeloid/subcluster_barplot.pdf',height = 4,width = 7)

##4.3############
Idents(Myeloid)<-'subcluster'
Myeloid.sce.markers <- FindAllMarkers(object = Myeloid,logfc.threshold = 0.25, min.pct = 0.25,only.pos = T)
head(Myeloid.sce.markers)
Myeloid.sce.markers["pct.diff"]=Myeloid.sce.markers$pct.1-Myeloid.sce.markers$pct.2
Myeloid.sce.markers <- Myeloid.sce.markers[Myeloid.sce.markers$p_val_adj<0.05,]
table(Myeloid.sce.markers$cluster)
length(unique(Myeloid.sce.markers$gene))#762
head(Myeloid.sce.markers)
table(Myeloid.sce.markers$cluster)
write.csv(Myeloid.sce.markers,'results/04.Myeloid/Myeloid_subcluster_degs.csv')


Myeloid.marker.genesets.list=split(x=Myeloid.sce.markers,f=Myeloid.sce.markers$cluster)
Myeloid.marker.genesets.list=sapply(Myeloid.marker.genesets.list, function(x){subset(x,select='gene',drop=TRUE)})
library(topGO)
library(clusterProfiler)
library(Rgraphviz)
library(org.Hs.eg.db)
entrez_id = mapIds(x = org.Hs.eg.db, keys = Myeloid.marker.genesets.list$Monocyte,
                   keytype = "SYMBOL",column = "ENTREZID")
entrez_id = na.omit(entrez_id)
Monocyte.erich.go.BP = enrichGO(gene = entrez_id,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05)
Monocyte.erich.go.BP.res=Monocyte.erich.go.BP@result
write.xlsx(Monocyte.erich.go.BP.res,'results/04.Myeloid/Myeloid_Monocyte_enrichment.xlsx',overwrite = T)###logfc=0.5

Monocyte.erich.go.BP.res=Monocyte.erich.go.BP.res[Monocyte.erich.go.BP.res$p.adjust<0.05,c('Description','GeneRatio','Count','p.adjust')]
Monocyte.erich.go.BP.res=Monocyte.erich.go.BP.res %>% slice_min(n =10, order_by = p.adjust) 
head(Monocyte.erich.go.BP.res)
p1=ggplot(data=Monocyte.erich.go.BP.res,aes(x=Count,y=reorder(Description,Count), color = -log10(p.adjust))) +
  geom_point(aes(size=Count),show.legend = T) +
  scale_color_gradient(low = "#8DD3C7", high = "#E6AB02")+
  geom_segment(aes(yend=Description,xend=0),size=.5) +
  labs(x='Count',y='')+theme_bw()+ggtitle('Monocyte')+
  theme(text = element_text(family = 'Times',size = 14))+
  scale_y_discrete(labels=function(y)str_wrap(y,width = 35))
p1


entrez_id = mapIds(x = org.Hs.eg.db, keys = Myeloid.marker.genesets.list$`M0 Macrophage`,
                   keytype = "SYMBOL",column = "ENTREZID")
entrez_id = na.omit(entrez_id)
M0_Mac.erich.go.BP = enrichGO(gene = entrez_id,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05)
M0_Mac.erich.go.BP.res=M0_Mac.erich.go.BP@result
write.xlsx(M0_Mac.erich.go.BP.res,'results/04.Myeloid/Myeloid_M0_Mac_enrichment.xlsx',overwrite = T)###logfc=0.25

M0_Mac.erich.go.BP.res=M0_Mac.erich.go.BP.res[M0_Mac.erich.go.BP.res$p.adjust<0.05,c('Description','GeneRatio','Count','p.adjust')]
M0_Mac.erich.go.BP.res=M0_Mac.erich.go.BP.res %>% slice_min(n =10, order_by = p.adjust) 
head(M0_Mac.erich.go.BP.res)
p2=ggplot(data=M0_Mac.erich.go.BP.res,aes(x=Count,y=reorder(Description,Count), color = -log10(p.adjust))) +
  geom_point(aes(size=Count),show.legend = T) +
  scale_color_gradient(low = "#8DD3C7", high = "#E6AB02")+
  geom_segment(aes(yend=Description,xend=0),size=.5) +
  labs(x='Count',y='')+theme_bw()+ggtitle('M0 Macrophage')+
  theme(text = element_text(family = 'Times',size = 14))+
  scale_y_discrete(labels=function(y)str_wrap(y,width = 35))
p2



entrez_id = mapIds(x = org.Hs.eg.db, keys = Myeloid.marker.genesets.list$`M1 Macrophage`,
                   keytype = "SYMBOL",column = "ENTREZID")
entrez_id = na.omit(entrez_id)
M1_Mac.erich.go.BP = enrichGO(gene = entrez_id,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05)
M1_Mac.erich.go.BP.res=M1_Mac.erich.go.BP@result
write.xlsx(M1_Mac.erich.go.BP.res,'results/04.Myeloid/Myeloid_M1_Mac_enrichment.xlsx',overwrite = T)###logfc=0.25

M1_Mac.erich.go.BP.res=M1_Mac.erich.go.BP.res[M1_Mac.erich.go.BP.res$p.adjust<0.05,c('Description','GeneRatio','Count','p.adjust')]
M1_Mac.erich.go.BP.res=M1_Mac.erich.go.BP.res %>% slice_min(n =10, order_by = p.adjust) 
head(M1_Mac.erich.go.BP.res)
p3=ggplot(data=M1_Mac.erich.go.BP.res,aes(x=Count,y=reorder(Description,Count), color = -log10(p.adjust))) +
  geom_point(aes(size=Count),show.legend = T) +
  scale_color_gradient(low = "#8DD3C7", high = "#E6AB02")+
  geom_segment(aes(yend=Description,xend=0),size=.5) +
  labs(x='Count',y='')+theme_bw()+ggtitle('M1 Macrophage')+
  theme(text = element_text(family = 'Times',size = 14))+
  scale_y_discrete(labels=function(y)str_wrap(y,width = 35))
p3

pdf('results/04.Myeloid/marker_enrichment.pdf',height = 6,width = 18,onefile = F)
mg_merge_plot(p1,p2,p3,ncol = 3,labels = LETTERS[4:6])
dev.off()
save(Myeloid.marker.genesets.list,file = 'results/04.Myeloid/Myeloid.marker.genesets.list.RData')




##4.4 ####
#####
Immunomodulator_and_chemokines=read.delim
Immunomodulator_and_chemokines=Immunomodulator_and_chemokines[1:150,]
head(Immunomodulator_and_chemokines)
table(Immunomodulator_and_chemokines$category)
MHC.genes=Immunomodulator_and_chemokines$g[Immunomodulator_and_chemokines$category=='MHC']

# 
# Idents(Myeloid)=rev(unique(Myeloid$stage))#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10687432/
# Idents(Myeloid)=rev(unique(Myeloid$stage))
####
Myeloid_m1=subset(Myeloid,subcluster=='M1 Macrophage')
p=DotPlot(object = Myeloid_m1, features =c('IL1B','CCL2',"CCL7","CCL8",'CXCL9','CXCL10','CXCL11',MHC.genes), 
        cols=c("skyblue", "orange"),scale = T,group.by = 'stage')+
  RotatedAxis()+ xlab('')+ylab('')+coord_flip()+
  theme(plot.title = element_text(hjust = 0.5))
ggsave('results/04.Myeloid/chemokine_MHC_plot.pdf',p,height = 6,width = 5)
p


##
countexp<-Myeloid_m1@assays$RNA@counts
countexp<-data.frame(as.matrix(countexp))
library(AUCell)
library(GSEABase)
cells_rankings <- AUCell_buildRankings(as.matrix(countexp)) #rank
####
file_paths=list.files('origin_datas/migration_gmt/')
read_and_process <- function(file_path) {
  df <- read.gmt(paste0('origin_datas/migration_gmt/',file_path))
  return(df)
}
list_of_dfs <- lapply(file_paths, read_and_process)
list_of_immune=sapply(list_of_dfs, function(x){subset(x,select='gene',drop=TRUE)})
names(list_of_immune)=str_split_fixed(file_paths,'[.]',2)[,1]
geneSets <- list_of_immune
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) #calc


signature_exp <- data.frame(getAUC(cells_AUC))
Myeloid_m1@assays$HALLMARK$score<-signature_exp

df = as.data.frame(t(Myeloid_m1@assays[["HALLMARK"]][["score"]]))
rownames(df) <- gsub(".", "-", rownames(df), fixed = TRUE)
df = df[names(Myeloid_m1$stage),]
df[1:5,1:5]


df$group <- Myeloid_m1$stage
avg_df =aggregate(df[,1:ncol(df)-1],list(df$group),mean)
rownames(avg_df) = avg_df$Group.1
avg_df=avg_df[,-1]
avg_df <- as.data.frame(t(avg_df))
rownames(avg_df)
head(avg_df)


pdf('results/04.Myeloid/migration_pathway_heatmap.pdf',height = 6,width = 7,onefile = F)
pheatmap(avg_df,
         show_colnames = T,scale='row',name = 'AUCell score',
         color=colorRampPalette(c('#1A5592','white',"orange"))(100),
         cluster_cols = F, cluster_rows = T,
         fontsize_row = 10,fontsize_col = 14)
dev.off()


#05.########
dir.create('results/05.cellchat')
library(CellChat)

colnames(Tcell@meta.data)
table(Tcell$subcluster)
table(Epithelial$subcluster)
Epithelial$subcluster=paste0('Epithelial_',Epithelial$subcluster)
table(Myeloid$subcluster)

Tcell.df=Tcell@meta.data[,'subcluster',drop=F]
Epithelial.df=Epithelial@meta.data[,'subcluster',drop=F]
Myeloid.df=Myeloid@meta.data[,'subcluster',drop=F]
my.data.df=rbind(Tcell.df,Epithelial.df,Myeloid.df)
head(my.data.df)
dim(my.data.df)

sce=readRDS('results/01.cell_annotation/sce.rds')
table(sce$cell_type)
my.data=sce[,rownames(my.data.df)]
rm(sce)
my.data=AddMetaData(object = my.data,metadata = my.data.df,col.name = 'subcluster')
table(my.data$subcluster)
library(CellChat)

table(my.data$subcluster)
my.data1=subset(my.data,stage=='Stage I')
my.data2=subset(my.data,stage=='Stage II')
my.data3=subset(my.data,stage=='Stage III')

cellchat.stage1 <- createCellChat(object = my.data1, meta = my.data1@meta.data, group.by = "subcluster")
rm(my.data1)
cellchat.stage1@DB <- subsetDB(CellChatDB = CellChatDB.human, search = "Secreted Signaling")
cellchat.stage1 <- subsetData(cellchat.stage1)
cellchat.stage1@data.signaling
#
cellchat.stage1 <- identifyOverExpressedGenes(cellchat.stage1)
cellchat.stage1 <- identifyOverExpressedInteractions(cellchat.stage1)
cellchat.stage1 <- projectData(cellchat.stage1, PPI.human)
# # #
# #
# # #
cellchat.stage1 <- computeCommunProb(cellchat.stage1)
cellchat.stage1 <- computeCommunProbPathway(cellchat.stage1)
cellchat.stage1 <- aggregateNet(cellchat.stage1)
cellchat.stage1 <- netAnalysis_computeCentrality(cellchat.stage1)
# # 
# # 
# library(NMF)
# selectK(cellchat.stage1, pattern = "outgoing")
cellchat.stage1 <- identifyCommunicationPatterns(cellchat.stage1, pattern = "outgoing", k = 5)


cellchat.stage2 <- createCellChat(object = my.data2, meta = my.data2@meta.data, group.by = "subcluster")
rm(my.data2)
cellchat.stage2@DB <- subsetDB(CellChatDB = CellChatDB.human, search = "Secreted Signaling")
cellchat.stage2 <- subsetData(cellchat.stage2)
cellchat.stage2@data.signaling
#
cellchat.stage2 <- identifyOverExpressedGenes(cellchat.stage2)
cellchat.stage2 <- identifyOverExpressedInteractions(cellchat.stage2)
cellchat.stage2 <- projectData(cellchat.stage2, PPI.human)
# # #
# #
# # #
cellchat.stage2 <- computeCommunProb(cellchat.stage2)
cellchat.stage2 <- computeCommunProbPathway(cellchat.stage2)
cellchat.stage2 <- aggregateNet(cellchat.stage2)
cellchat.stage2 <- netAnalysis_computeCentrality(cellchat.stage2)
# # 
# # 
# library(NMF)
# selectK(cellchat.stage2, pattern = "outgoing")
cellchat.stage2 <- identifyCommunicationPatterns(cellchat.stage2, pattern = "outgoing", k = 5)



cellchat.stage3 <- createCellChat(object = my.data3, meta = my.data3@meta.data, group.by = "subcluster")
rm(my.data3)
cellchat.stage3@DB <- subsetDB(CellChatDB = CellChatDB.human, search = "Secreted Signaling")
cellchat.stage3 <- subsetData(cellchat.stage3)
cellchat.stage3@data.signaling
#
cellchat.stage3 <- identifyOverExpressedGenes(cellchat.stage3)
cellchat.stage3 <- identifyOverExpressedInteractions(cellchat.stage3)
cellchat.stage3 <- projectData(cellchat.stage3, PPI.human)
# # #
# #
# # #
cellchat.stage3 <- computeCommunProb(cellchat.stage3)
cellchat.stage3 <- computeCommunProbPathway(cellchat.stage3)
cellchat.stage3 <- aggregateNet(cellchat.stage3)
cellchat.stage3 <- netAnalysis_computeCentrality(cellchat.stage3)
# # 
# # 
# library(NMF)
# selectK(cellchat.stage3, pattern = "outgoing")
cellchat.stage3 <- identifyCommunicationPatterns(cellchat.stage3, pattern = "outgoing", k = 5)



########
object.list <- list(stage1 = cellchat.stage1, stage2 = cellchat.stage2,stage3 = cellchat.stage3)
cellchat_merge <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(cellchat,file = 'results/05.cellchat/cellchat_merge.rds')


cellchat_merge=readRDS('results/05.cellchat/cellchat_merge.rds')
gg1 <- compareInteractions(cellchat_merge, show.legend = F, group = c(1:3))
gg2 <- compareInteractions(cellchat_merge, show.legend = F, group = c(1:3), measure = "weight")
gg1 + gg2





# 
# 
# 
gg1 <- rankNet(cellchat_merge, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(1,3),color.use = c("#7570B3","#FB8072"))
gg2 <- rankNet(cellchat_merge, mode = "comparison", stacked = F, do.stat = TRUE,comparison = c(1,3),color.use = c("#7570B3","#FB8072"))
pdf('results/05.cellchat/cellchat_barplot.pdf',height = 6,width = 9,onefile = F)
gg1 + gg2
dev.off()



# [1] "Tregs"                "cytotoxic NK/T cells" "naive T cells"        "Epithelial_C1"        "Epithelial_C2"        "Epithelial_C3"       
# [7] "Monocyte"             "M2 Macrophage"        "M1 Macrophage"        "M0 Macrophage"   


unique(cellchat_merge@idents[["stage1"]])
pdf('results/05.cellchat/cellchat_dotplot_1.pdf',height = 8,width = 4.5,onefile = F)
netVisual_bubble(cellchat_merge,sources.use = 5, 
             targets.use = c(1,8,9),
             comparison = c(1,3), angle.x = 45)
dev.off()

pdf('results/05.cellchat/cellchat_dotplot_2.pdf',height = 8,width = 4.5,onefile = F)
netVisual_bubble(cellchat_merge,sources.use = c(6) , 
             targets.use = c(1:2),
             comparison = c(1,3), angle.x = 45)
dev.off()

pdf('results/05.cellchat/cellchat_dotplot_3.pdf',height = 8,width = 4.5,onefile = F)
netVisual_bubble(cellchat_merge,sources.use = c(7) , 
             targets.use = c(1:2),
             comparison = c(1,3), angle.x = 45)
dev.off()


pdf('results/05.cellchat/cellchat_dotplot_4.pdf',height = 8,width = 6,onefile = F)
netVisual_bubble(cellchat_merge,sources.use = c(5:7) , 
                 targets.use = 'cytotoxic NK/T cells',
                 comparison = c(1,3), angle.x = 45)
dev.off()


pdf('results/05.cellchat/cellchat_dotplot_5.pdf',height = 8,width = 6,onefile = F)
netVisual_bubble(cellchat_merge, targets.use= c(5:7) , 
                 sources.use = 'cytotoxic NK/T cells',
                 comparison = c(1,3), angle.x = 45)
dev.off()
