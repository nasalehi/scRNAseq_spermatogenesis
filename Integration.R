library(dplyr)
library(Seurat)
library(pheatmap)
library(ggplot2)
library(cowplot)
library(stats4)
library(splines)
library(VGAM)
library(irlba)
library(DDRTree)
library(reshape2)
library(tidyr)
library(Biobase)
library(monocle)
library(monocle3)
install.packages("readxl")
packageVersion('monocle')
#packageVersion('monocle3')
setwd("/home/salehi/Projects/scRNAseq/Data")
setwd("/home/nsalehi/scRNAseq/Data")
setwd("C:/Najmeh/Education/MyProject/scRNAseq/Data/")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

abc
######################################Initial Data############################################################################
spermatogonia.data <- Read10X(data.dir = "GSE109037/AdultHumanSpermatogonia-reps1-3/")
Spermatocytes.data <- Read10X(data.dir = "GSE109037/AdultHuman-Spermatocytes_Reps1-2")
Spermatids.data <- Read10X(data.dir = "GSE109037/AdultHuman-Spermatids_Reps1-2")
Spermatogenesis.data <- Read10X(data.dir = "GSE109037/AdultHumanSpermatogenesis-reps1-3")
C1.data <- read.delim("GSE106487/C1.txt")
dim(spermatogonia.data)
#33694 11104
dim(Spermatocytes.data)
#33694  4884
dim(Spermatids.data)
#33694  7434
dim(Spermatogenesis.data)
#33694  7434
dim(C1.data)
#24153  3059
spermatogonia <- CreateSeuratObject(counts = spermatogonia.data, project = "spermatogonia", min.cells = 3, min.features = 200)
Spermatocytes <- CreateSeuratObject(counts = Spermatocytes.data, project = "Spermatocytes", min.cells = 3, min.features = 200)
Spermatids <- CreateSeuratObject(counts = Spermatids.data, project = "Spermatids", min.cells = 3, min.features = 200)
Spermatogenesis <- CreateSeuratObject(counts = Spermatogenesis.data , project = "Spermatogenesis", min.cells = 3, min.features = 200)
C1<- CreateSeuratObject(counts = C1.data , project = "C1", min.cells = 3, min.features = 200)
dim(spermatogonia)
#28320 11104
dim(Spermatocytes)
#27405  4884
dim(Spermatids)
#27898  7434
dim(Spermatogenesis)
#29119  7132
dim(C1)
#22911  3046
spermatogonia[["percent.mt"]] <- PercentageFeatureSet(spermatogonia, pattern = "^MT-")
Spermatocytes[["percent.mt"]] <- PercentageFeatureSet(Spermatocytes, pattern = "^MT-")
Spermatids[["percent.mt"]] <- PercentageFeatureSet(Spermatids, pattern = "^MT-")
Spermatogenesis[["percent.mt"]] <- PercentageFeatureSet(Spermatogenesis, pattern = "^MT-")
C1[["percent.mt"]] <- PercentageFeatureSet(C1, pattern = "^MT-")
head(spermatogonia@meta.data, 5)
head(Spermatocytes@meta.data, 5)
head(Spermatids@meta.data, 5)
head(Spermatogenesis@meta.data, 5)
head(C1@meta.data, 5)
g1 <- merge(x = spermatogonia, y = list(Spermatocytes, Spermatids,Spermatogenesis))
######################################QC######################################
pdf("./Results2/quality1_1_vlnplot.pdf", width=20)
lapply(c(g1,Spermatogenesis,spermatogonia,Spermatocytes,Spermatids),VlnPlot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

outliers <- function(x) {
  qnt <- quantile(x, probs=c(.25, .75))
  H <- 1.5 * IQR(x)
  c((qnt[1] - H),(qnt[2] + H))
  
}
#sort(spermatogonia@meta.data$nFeature_RNA)[1:50]
#sort(spermatogonia@meta.data$nFeature_RNA,decreasing = TRUE)[1:30]
#sort(spermatogonia@meta.data$percent.mt,decreasing = TRUE)[1:45]
#outliers(spermatogonia@meta.data$nFeature_RNA)
#outliers(spermatogonia@meta.data$percent.mt)
spermatogonia <- subset(spermatogonia, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 20)
dim(spermatogonia)
#28320 10998

#sort(Spermatocytes@meta.data$nFeature_RNA)[1:30]
#sort(Spermatocytes@meta.data$nFeature_RNA,decreasing = TRUE)[1:30]
#sort(Spermatocytes@meta.data$percent.mt,decreasing = TRUE)[1:50]
Spermatocytes <- subset(Spermatocytes, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 10)
dim(Spermatocytes)
#27405  4797


#sort(Spermatids@meta.data$nFeature_RNA)[1:30]
#sort(Spermatids@meta.data$nFeature_RNA,decreasing = TRUE)[1:20]
#sort(Spermatids@meta.data$percent.mt,decreasing = TRUE)[1:50]
#outliers(Spermatids@meta.data$nFeature_RNA)
#outliers(Spermatids@meta.data$percent.mt)
Spermatids <- subset(Spermatids, subset = nFeature_RNA > 1200 & nFeature_RNA < 9000 & percent.mt < 15)
dim(Spermatids)
#27898  7382

#sort(Spermatogenesis@meta.data$nFeature_RNA)[1:30]
#sort(Spermatogenesis@meta.data$nFeature_RNA,decreasing = TRUE)[1:30]
#sort(Spermatogenesis@meta.data$percent.mt,decreasing = TRUE)[1:50]
Spermatogenesis <- subset(Spermatogenesis, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 25)
dim(Spermatogenesis)
#29119  7037

g1 <- merge(x = spermatogonia, y = list(Spermatocytes, Spermatids,Spermatogenesis))
pdf("./Results2/quality1_2_vlnplot.pdf", width=20)
lapply(c(g1,Spermatogenesis,spermatogonia,Spermatocytes,Spermatids),VlnPlot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

Idents(C1)<-rep("C1",dim(C1)[2])
pdf("./Results2/quality2_1_vlnplot.pdf", width=20)
VlnPlot(C1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

C1 <- subset(C1, subset = nFeature_RNA > 800 & nFeature_RNA < 12000)
dim(C1)

pdf("./Results2/quality2_2_vlnplot.pdf", width=20)
VlnPlot(C1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
######################################Normalization######################################
all.list <- c(Spermatogenesis,spermatogonia,Spermatocytes, Spermatids,C1)
names(all.list) <- c("Spermatogenesis", "spermatogonia","Spermatocytes", "Spermatids", "C1")

for (i in names(all.list)) {
  all.list[[i]] <- NormalizeData(all.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
}
######################################Variable############################################################################

for (i in names(all.list)) {
  all.list[[i]] <- FindVariableFeatures(all.list[[i]], selection.method = "vst", nfeatures = 2000)
}
######################################intersec############################################################################
length(Reduce(intersect, list(Spermatogenesis@assays$RNA@data@Dimnames[[1]],spermatogonia@assays$RNA@data@Dimnames[[1]],Spermatocytes@assays$RNA@data@Dimnames[[1]],Spermatids@assays$RNA@data@Dimnames[[1]],C1@assays$RNA@data@Dimnames[[1]])))
# 18659
length(Reduce(union, list(Spermatogenesis@assays$RNA@data@Dimnames[[1]],spermatogonia@assays$RNA@data@Dimnames[[1]],Spermatocytes@assays$RNA@data@Dimnames[[1]],Spermatids@assays$RNA@data@Dimnames[[1]],C1@assays$RNA@data@Dimnames[[1]])))
#33551
######################################removing############################################################################
save.image(file = "scRNAseq1.RData")
remove(spermatogonia)
remove(Spermatocytes)
remove(Spermatids)
remove(Spermatogenesis)
remove(C1)
remove(spermatogonia.data)
remove(Spermatocytes.data)
remove(Spermatids.data)
remove(Spermatogenesis.data )
remove(C1.data)
remove(g1)
######################################INTEGRATION############################################################################
all.anchors <- FindIntegrationAnchors(object.list=all.list,dims = 1:36)
all.integrated <- IntegrateData(anchorset = all.anchors, dims = 1:36)
dim(all.integrated )
######################################removing############################################################################
save.image(file = "scRNAseq2_36.RData")
load("scRNAseq2.RData")
remove(all.list)
remove(all.anchors)
######################################Dimention Reduction#########################################################
DefaultAssay(all.integrated) <- "integrated"
all.integrated <- ScaleData(all.integrated, verbose = FALSE)
all.integrated <- RunPCA(all.integrated, npcs = 50, verbose = FALSE)
all.integrated <- RunUMAP(all.integrated, reduction = "pca", dims = 1:36)
dim(all.integrated)
cellTypes <- Idents(all.integrated)
##Clusters
all.integrated <- FindNeighbors(all.integrated, dims = 1:35)
all.integrated <- FindClusters(all.integrated, resolution = 0.2)
class<- Idents(all.integrated)
table(Idents(all.integrated))


all.integrated@reductions[["umap"]]@cell.embeddings <- all.integrated@reductions[["umap"]]@cell.embeddings*-1
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
new.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,16)
Idents(all.integrated) <- plyr::mapvalues(x = Idents(all.integrated), from = current.cluster.ids, to = new.cluster.ids)

#save.image(file = "scRNAseq3_36.RData")
load("scRNAseq3_36.RData")
class2<- Idents(all.integrated)


current.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,16)
new.cluster.ids <- c("Undiff. SPG1", "Undiff. SPG2", "Round SPT2", "Round SPT1", "Pachytene SPC", "Zygotene SPC",
                      "Early round SPT", "Diplotene SPC", "Leptotene SPC", "Diff.ing SPG1", "Elongating SPT", 
                     "Diff.ed SPG","Diff.ing SPG2", "Somatic cells", "Somatic cells", "Somatic cells")
Idents(all.integrated) <- plyr::mapvalues(x = Idents(all.integrated), from = current.cluster.ids, to = new.cluster.ids)
cellTypes3<-Idents(all.integrated)
DefaultAssay(all.integrated) <- "RNA"

my_levels <-  c("Somatic cells","Undiff. SPG1", "Undiff. SPG2","Diff.ing SPG1","Diff.ing SPG2", "Diff.ed SPG","Leptotene SPC", "Pachytene SPC", "Zygotene SPC", "Diplotene SPC", "Early round SPT", "Round SPT1", "Round SPT2", "Elongating SPT")
# Re-level object@ident
all.integrated@active.ident <- factor(x = all.integrated@active.ident, levels = my_levels)


pdf("./Results2/integrated15.pdf", width=10)
Idents(all.integrated) <- cellTypes
current.cluster.ids <- c("C1", "Spermatids", "Spermatocytes", "Spermatogenesis", "spermatogonia")
new.cluster.ids <- c("Spermatogenesis2", "Spermatid", "Spermatocyte", "Spermatogenesis1", "Spermatogonia")
Idents(all.integrated) <- plyr::mapvalues(x = Idents(all.integrated), from = current.cluster.ids, to = new.cluster.ids)

DimPlot(all.integrated, reduction = "umap", group.by = 'ident', repel = TRUE, 
        cols= c( "blue3","lightskyblue","thistle","darkseagreen","lightpink1"),
        order = c("Spermatogonia", "Spermatocyte", "Spermatid", "Spermatogenesis1","Spermatogenesis2" ))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
Idents(all.integrated) <- class
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
new.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,16)
Idents(all.integrated) <- plyr::mapvalues(x = Idents(all.integrated), from = current.cluster.ids, to = new.cluster.ids)
DimPlot(all.integrated, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Results2/integrated12.pdf", width=10)
DimPlot(all.integrated, reduction = "umap", group.by = 'ident', repel = TRUE, 
        cols= c( "gray","gray","gray","gray","lightskyblue"),
        order = c("Spermatogenesis1", "Spermatocyte", "Spermatid","Spermatogonia","Spermatogenesis2"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))

DimPlot(all.integrated, reduction = "umap", group.by = 'ident', repel = TRUE, 
        cols= c( "gray","gray","gray","gray","blue3"),
        order = c( "Spermatogenesis2","Spermatocyte", "Spermatid", "Spermatogenesis1","Spermatogonia" ))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))

DimPlot(all.integrated, reduction = "umap", group.by = 'ident', repel = TRUE, 
        cols= c( "gray","gray","gray","gray","lightpink1"),
        order = c("Spermatogonia","Spermatocyte", "Spermatid","Spermatogenesis2", "Spermatogenesis1"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))

DimPlot(all.integrated, reduction = "umap", group.by = 'ident', repel = TRUE, 
        cols= c( "gray","gray","gray","gray","darkseagreen"),
        order = c("Spermatocyte","Spermatogonia", "Spermatid","Spermatogenesis2", "Spermatogenesis1"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))

DimPlot(all.integrated, reduction = "umap", group.by = 'ident', repel = TRUE, 
        cols= c( "gray","gray","gray","gray","thistle"),
        order = c("Spermatid","Spermatocyte","Spermatogonia", "Spermatogenesis2", "Spermatogenesis1"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))

dev.off()


pdf("./Results2/UMAP_dim.pdf", width=8)
DimPlot(all.integrated[,1:10998], reduction = "umap",group.by = 'ident')
DimPlot(all.integrated[,10999:15795], reduction = "umap",group.by = 'ident')
DimPlot(all.integrated[,15795:23177], reduction = "umap",group.by = 'ident')
DimPlot(all.integrated[,23178:30214], reduction = "umap",group.by = 'ident')
DimPlot(all.integrated[,30215:33185], reduction = "umap",group.by = 'ident')
dev.off()

prop.table(table(Idents(all.integrated)))
clusdata <-table(Idents(all.integrated), cellTypes)
write.csv(clusdata, file= "./Results2/clusdata_integrated.csv")

pdf("./Results2/PCA_all.pdf", width=10)
VizDimLoadings(all.integrated, dims = 1:2, reduction = "pca")
Idents(all.integrated) <- cellTypes
DimPlot(all.integrated, reduction = "pca",group.by = 'ident')
Idents(all.integrated) <- class
DimPlot(all.integrated, reduction = "pca",group.by = 'ident')
DimHeatmap(all.integrated, dims = 1:3, cells = 500, balanced = TRUE)
ElbowPlot(all.integrated, ndims = 50)
dev.off()

all.integrated <-RunTSNE(all.integrated, dims = 1:35)
pdf("./Results2/tSNE_all.pdf", width=10)
Idents(all.integrated) <- cellTypes
DimPlot(all.integrated, reduction = "tsne")
Idents(all.integrated) <- class
DimPlot(all.integrated, reduction = "tsne")
dev.off()

######################################MARKERS############################################################################
load("scRNAseq3.RData")
DefaultAssay(all.integrated) <- "RNA"
pdf("./Results2/markers_allfeature2.pdf", width=10)
FeaturePlot(all.integrated, features = c("DDX4","GPR125"))
FeaturePlot(all.integrated, features = c("GFRA1","RET","STRA8","NANOS2","ZBTB16","SALL4","POU3F1","UTF1","NANOS3","FGFR3","BMPR1B",""))
FeaturePlot(all.integrated, features = c("KIT","GPR125","MKI67"))
FeaturePlot(all.integrated, features = c("STRA8","GPR125","MKI67"))
FeaturePlot(all.integrated, features = c("ALDH1A1","CYP26B1","ALDH1A2"))
FeaturePlot(all.integrated, features = c("ALDH1A2","ALDH1A3","SPC7"))
dev.off()

pdf("./Results2/STAG3_2.pdf", width=10)
FeaturePlot(all.integrated, features = c("STAG3","REC8","RAD21L1", "RAD21"))
DotPlot(all.integrated, features = c("STAG3","REC8","RAD21L1", "RAD21"),col.min = 0) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

DotPlot(all.integrated, features = c("STAG3","ESRRA", "ESR1"),col.min = 0) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))

pdf("./Results2/Marker_germcells3.pdf", width=15)
FeaturePlot(all.integrated, features = c("DDX4"),cols=c("lightgrey","brown"))
dev.off()
pdf("./Results2/newMarkers.pdf", width=10)
FeaturePlot(all.integrated, features = c("ZNF428", "UTF1", "C19orf84","ID4", "CCNI", "RPS19","TKTL1", "ASB9", "DMRT1", "HIST1H4C","PTMA", "TUBA1B", "HIST1H4C", "PTMA", "GAGE2A", "TEX101","SYCP1","MEIOB", "CETN3","SPATA8", "RNFT1", "PPP3R2", "CLGN","PRKCDBP", "GLIPR1L1","PRKCDBP","LDHC", "LINC00643","GOLGA6L2","CDRT15", "GOLGA6L2, "FAM24A", "EQTN"
                                         
                                         
                                         
                                         
),ncol = 3)
######################################Somatic cell
pdf("./Results2/Markers_somatic3.pdf", width=15)
FeaturePlot(all.integrated, features = c("CYP26B1","ALDH1A1", "MYH11", "INSL3", "CD68","CD163"),cols=c("lightgrey","darkblue")) + RotatedAxis()
VlnPlot(all.integrated, c("CYP26B1","ALDH1A1", "MYH11", "INSL3", "CD68","CD163"))
dev.off()
pdf("./Results2/Markers_somatic_4.pdf", width=15)
DotPlot(all.integrated, features = c("CYP26B1","ALDH1A1", "MYH11", "INSL3", "CD68","CD163"),cols=c("lightgrey","darkblue")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
DotPlot(all.integrated, features = c("GDNF","BMP4","FGF9","GATA4","GATA1","AMH","SOX9","SRY","NR5A1","CYP26B1", "ALDH1A1","DMRT1","MYH11","DLK1","CD68","CD163"),cols=c("lightgrey","darkblue")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
DotPlot(all.integrated, features = c("PHF7","ICAM1","PLIN2","BAX", "CST3","FATE1","KRT18","PTGDS"),cols=c("lightgrey","darkblue")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Results2/Markers_somatic_5.pdf", width=15)
DotPlot(all.integrated, features = c("CYP26B1","INSL3", "MYH11", "ALDH1A1", "CD68","CD163"),cols=c("lightgrey","darkblue")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Results2/Markers_sertoli.pdf", width=15)
FeaturePlot(all.integrated, features = c("AMH","SOX9","MYH11","DLK1","CD163","CD68","ALDH1A1","CYP26B1","ALDH1A2"))
VlnPlot(all.integrated, c("AMH","SOX9","MYH11","DLK1","CD163","CD68","ALDH1A1","CYP26B1","ALDH1A2"))
dev.off()
#######################################Spermatogonia
pdf("./Results2/Markers_spermatogonia_2.pdf", width=15)
FeaturePlot(all.integrated, features = c("POU2F2","CD9", "ENO2", "EPCAM","FGFR3","FMR1","GFRA1"))
VlnPlot(all.integrated, c("POU2F2","CD9", "ENO2", "EPCAM","FGFR3","FMR1","GFRA1"))
FeaturePlot(all.integrated, features = c("ID4","ZBTB16","SPOCD1","UCHL1","CDH1","UTF1","ZKSCAN2","PIWIL4"))
VlnPlot(all.integrated, c("ID4","UCHL1","ZBTB16","SPOCD1","CDH1","UTF1","ZKSCAN2","PIWIL4"))
FeaturePlot(all.integrated, features = c("SALL4","PHF13","ITGA6", "HMGA1", "EXOSC10","LMO4", "MAGEA4" ))
VlnPlot(all.integrated, c("SALL4","PHF13","ITGA6", "HMGA1", "EXOSC10","LMO4", "MAGEA4" ))
FeaturePlot(all.integrated, features = c("ELAVL2","PAX7","DMRT1","SSX2","SSX4", "KIT","STRA8", "SAGE1", "TRAPPC6A", "SSX1", "SYCP3"))
VlnPlot(all.integrated, c("ELAVL2","PAX7","DMRT1","SSX2","SSX4", "KIT","STRA8", "SAGE1", "TRAPPC6A", "SSX1", "SYCP3"))
dev.off()

pdf("./Results2/Markers_spermatogonia_8.pdf", width=12)
FeaturePlot(all.integrated, features = c("NANOS2","PIWIL4","GFRA1","SALL4","FGFR3","UTF1","LMO4","HMGA1","MAGEA4" ))
VlnPlot(all.integrated, c("NANOS2","PIWIL4","GFRA1","SALL4","FGFR3","UTF1","LMO4","HMGA1","MAGEA4"))
dev.off()

pdf("./Results2/Markers_spermatogonia_9.pdf", width=15)
DotPlot(all.integrated, features = c("PIWIL4","NANOS2","GFRA1","SALL4","FGFR3","UTF1","LMO4","HMGA1","MAGEA4" )) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Results2/Markers_spermatogonia_10.pdf", width=15)
DotPlot(all.integrated, features = c("PIWIL4","NANOS2","GFRA1","SALL4","MAGEA4","HMGA1" )) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Results2/Markers_spermatogonia_11.pdf", width=15)
DotPlot(all.integrated, features = c("ID4", "KIT", "POU2F2","OCT2" )) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
############################################spermatocyte
pdf("./Results2/Markers_spermatocyte.pdf", width=15)
FeaturePlot(all.integrated, features = c("SPO11","DMC1","RAD51AP2","SYCP3","PIWIL1","PGK2","ACR","OVOL2","NME8"))
VlnPlot(all.integrated, c("SPO11","DMC1","RAD51AP2","SYCP3","PIWIL1","PGK2","ACR","OVOL2","NME8"))
dev.off()
pdf("./Results2/Markers_spermatocyte_5.pdf", width=15)
DotPlot(all.integrated, features = c("DMC1", "SPO11", "RAD51AP2","PIWIL1","SYCP3", "OVOL2","NME8" ),cols=c("lightgrey","red")) + RotatedAxis()+
theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Results2/Markers_spermatocyte_6.pdf", width=15)
FeaturePlot(all.integrated, features = c("DMC1", "RAD51AP2","PIWIL1","SYCP3", "OVOL2","NME8"),cols=c("lightgrey","red"))
VlnPlot(all.integrated, c("DMC1", "RAD51AP2","PIWIL1","SYCP3", "OVOL2","NME8"))
dev.off()
pdf("./Results2/Markers_spermatocyte_7.pdf", width=15)
DotPlot(all.integrated, features = c("DMC1", "RAD51AP2","PIWIL1","SYCP3", "OVOL2" ),cols=c("lightgrey","red")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

##############################################spermatid
pdf("./Results2/Markers_spermatid_4.pdf", width=15)
FeaturePlot(all.integrated, features = c("TEX29","TXNDC2","SUN5","SPEM1","ACR","PGK2"),cols=c("lightgrey","darkgreen"))
VlnPlot(all.integrated, c("TEX29","TXNDC2","SUN5","SPEM1","ACR","PGK2"))
dev.off()
pdf("./Results2/Markers_spermatid_5.pdf", width=15)
DotPlot(all.integrated, features = c("SYCP3","OVOL2", "NME8","TEX29","TXNDC2","SUN5","SPEM1","ACR","PGK2"),cols=c("lightgrey","darkgreen")) + RotatedAxis()+
theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
pdf("./Results2/Markers_spermatid_7.pdf", width=15)
DotPlot(all.integrated, features = c("TEX29","SUN5","SPEM1","ACR","PGK2"),cols=c("lightgrey","darkgreen")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Results2/ACE2.pdf", width=10)
FeaturePlot(all.integrated, features = c("ACE2"))
VlnPlot(all.integrated, c("ACE2"))
dev.off()
theme(panel.border = element_rect(colour = "black", fill=NA, size=2))


pdf("./Results2/umap_remarks4.pdf", width=10)
DimPlot(object = all.integrated, reduction = "umap", label = TRUE, repel = TRUE)+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
######################################DEG##############################################################
load("scRNAseq3_36.RData")
##On integrated data
Idents(all.integrated) <- class 
all.markers <- FindAllMarkers(all.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(all.markers, file= "./Results2/DEG_all_integrated.csv")
#
all.markers2 <- FindAllMarkers(all.integrated, logfc.threshold = 0.25, test.use = "roc")
write.csv(all.markers2, file= "./Results2/DEG_all_integrated_ROC.csv")

#########################on initial data
DefaultAssay(all.integrated) <- "RNA"
Idents(all.integrated) <- class 
all.markers3 <- FindAllMarkers(all.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(all.markers3, file= "./Results2/DEG_all_integrated2.csv")
all.markers4 <- FindAllMarkers(all.integrated, logfc.threshold = 0.25, test.use = "roc")
write.csv(all.markers4, file= "./Results2/DEG_all_integrated_ROC2.csv")

pdf("./Results2/DEG_all_integrated_heatmap.pdf", width=10)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(all.markers, features = top10$gene) + NoLegend()
dev.off()

Idents(all.integrated) <- cellTypes2
all.markers10 <- FindAllMarkers(all.integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(all.markers10, file= "./Results2/DEG_all_integrated10.csv")
all.markers6 <- FindAllMarkers(all.integrated, logfc.threshold = 0.25, test.use = "roc")
write.csv(all.markers6, file= "./Results2/DEG_all_integrated_ROC3.csv")

#####CellTypes
spermatogonia1.markers <- FindMarkers(all.integrated, ident.1 = c("0","1"), min.pct = 0.25)
write.csv(spermatogonia1.markers, file= "./Results2/DEG_spermatogonia1.csv")

spermatogonia2.markers <- FindMarkers(all.integrated, ident.1 = c("9","12"), min.pct = 0.25)
write.csv(spermatogonia2.markers, file= "./Results2/DEG_spermatogonia23.csv")

Spermatocytes.markers <- FindMarkers(all.integrated, ident.1 = c("4","5","7","8"), min.pct = 0.25)
write.csv(Spermatocytes.markers, file= "./Results2/DEG_Spermatocytes.csv")

Spermatids.markers <- FindMarkers(all.integrated, ident.1 = c("2","3","6","10"), min.pct = 0.25)
write.csv(Spermatids.markers, file= "./Results2/DEG_Spermatidss.csv")
Spermatids.markers <- FindMarkers(all.integrated, ident.1 = c("2","3"), min.pct = 0.25)
write.csv(Spermatids.markers, file= "./Results2/DEG_2Spermatidss.csv")

cluster.markers <- FindMarkers(all.integrated, ident.1 = "leptotene",  min.pct = 0.25)
write.csv(cluster.markers, file= "./Results2/DEG_leptotene.csv")

spermatogonia12.markers <- FindMarkers(all.integrated, ident.1 = c("0","1"), logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
write.csv(spermatogonia12.markers, file= "./Results2/DEG_spermatogonia12.csv")

spermatogonia22.markers <- FindMarkers(all.integrated, ident.1 = c("9","11","12"), logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
write.csv(spermatogonia22.markers, file= "./Results2/DEG_spermatogonia22.csv")

Spermatocytes2.markers <- FindMarkers(all.integrated, ident.1 = c("4","5","7","8"), logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
write.csv(Spermatocytes2.markers , file= "./Results2/DEG_Spermatocytes2.csv")

Spermatids2.markers <- FindMarkers(all.integrated, ident.1 = c("2","3","6","10"), logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
write.csv(Spermatids2.markers , file= "./Results2/DEG_Spermatids2.csv")
###each type vs previous type
cluster1.markers <- FindMarkers(all.integrated, ident.1 = 1, ident.2 = 0, min.pct = 0.25)
write.csv(cluster1.markers, file= "./Results2/cluster1.markers.csv")

cluster12.markers <- FindMarkers(all.integrated, ident.1 = 12, ident.2 = 1, min.pct = 0.25)
write.csv(cluster12.markers, file= "./Results2/cluster12.markers.csv")

cluster9.markers <- FindMarkers(all.integrated, ident.1 = 9, ident.2 = 1, min.pct = 0.25)
write.csv(cluster9.markers, file= "./Results2/cluster9.markers.csv")

cluster11.markers <- FindMarkers(all.integrated, ident.1 = 11, ident.2 = 9, min.pct = 0.25)
write.csv(cluster11.markers, file= "./Results2/cluster11-9.markers.csv")
cluster11.markers <- FindMarkers(all.integrated, ident.1 = 11, ident.2 = 12, min.pct = 0.25)
write.csv(cluster11.markers, file= "./Results2/cluster11-12.markers.csv")

cluster8.markers <- FindMarkers(all.integrated, ident.1 = 8, ident.2 = 11, min.pct = 0.25)
write.csv(cluster8.markers, file= "./Results2/cluster8_11.markers.csv")

cluster8.markers <- FindMarkers(all.integrated, ident.1 = 8, ident.2 = 11, min.pct = 0.25)
write.csv(cluster8.markers, file= "./Results2/cluster8_11.markers.csv")
cluster8.markers <- FindMarkers(all.integrated, ident.1 = 8, ident.2 = 12, min.pct = 0.25)
write.csv(cluster8.markers, file= "./Results2/cluster8_12.markers.csv")

cluster5.markers <- FindMarkers(all.integrated, ident.1 = 5, ident.2 = 8, min.pct = 0.25)
write.csv(cluster5.markers, file= "./Results2/cluster5.markers.csv")

cluster4.markers <- FindMarkers(all.integrated, ident.1 = 4, ident.2 = 5, min.pct = 0.25)
write.csv(cluster4.markers, file= "./Results2/cluster4.markers.csv")

cluster7.markers <- FindMarkers(all.integrated, ident.1 = 7, ident.2 = 4, min.pct = 0.25)
write.csv(cluster7.markers, file= "./Results2/cluster7.markers.csv")

cluster6.markers <- FindMarkers(all.integrated, ident.1 = 6, ident.2 = 7, min.pct = 0.25)
write.csv(cluster6.markers, file= "./Results2/cluster6.markers.csv")

cluster3.markers <- FindMarkers(all.integrated, ident.1 = 3, ident.2 = 6, min.pct = 0.25)
write.csv(cluster3.markers, file= "./Results2/cluster3.markers.csv")

cluster2.markers <- FindMarkers(all.integrated, ident.1 = 2, ident.2 = 3, min.pct = 0.25)
write.csv(cluster2.markers, file= "./Results2/cluster2.markers.csv")

cluster10.markers <- FindMarkers(all.integrated, ident.1 = 10, ident.2 = 2, min.pct = 0.25)
write.csv(cluster10.markers, file= "./Results2/cluster10.markers.csv")

#

data_1 <- readxl::read_excel("Results2/DEG/summary6.xlsx")
pdf("./Results2/DEG/DEG_all3.pdf", width=10)
positions <- c("UnDiff. SPG1", "UnDiff. SPG2","Diff.ing SPG1","Diff.ing SPG2", "Diff.ed SPG","Leptotene SPC", "Zygotene SPC", "Pachytene SPC",  "Diplotene SPC", "Early round SPT", "Round SPT1", "Round SPT2", "Elongating SPT")
data_1 %>%
  ggplot(aes(x = Cell_Types, y = DEG_No, fill = DEG))+ scale_x_discrete(limits = positions) +
  geom_bar(stat = "identity")+ 
  theme(text = element_text(size=20))+
  coord_flip()
dev.off()
pdf("./Results2/enrichment5.pdf", width=13)
ggplot(df_sp, aes(x = sample, y = BP, fill = value)) +
  geom_tile(aes(fill = value)) + labs(fill = "Benjamini Score") + 
  scale_y_discrete(limits = BioP) + scale_x_discrete(limits = positions)+
  scale_fill_gradientn(colours = c("blue","white","red"))+ theme(text = element_text(size=20), axis.text.x = element_text(angle = 90))
dev.off()

data_2 <- readxl::read_excel("Results2/DEG/summary7.xlsx")
BioP <- c("translation","cell cycle","chromosome organization", 
          "DNA metabolic process", "cellular macromolecular complex subunit organization", 
           "spermatogenesis","meiosis",
         "spermatid development","cell wall macromolecule catabolic process", 
         "sperm motility", "nucleus organization", "spermatid nucleus differentiation"
)
positions <- c("UnDiff. SPG1", "UnDiff. SPG2","Diff.ing SPG1","Diff.ing SPG2", "Diff.ed SPG","Leptotene SPC", "Zygotene SPC", "Pachytene SPC", "Diplotene SPC", "Early round SPT", "Round SPT1", "Round SPT2", "Elongating SPT")

df_sp = data_2 %>% gather(sample, value, -BP)


top10 <- all.markers7%>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf("./Results2/DEg_Heatmap.pdf", width=10)
DoHeatmap(all.integrated, features = top10$gene) + NoLegend()
dev.off()
######################################Enrichment##############################################################
library("RDAVIDWebService")
david<-DAVIDWebService$new(email="nsalehi@ut.ac.ir")
data(demoList1)
result<-addList(david, demoList1,idType="AFFYMETRIX_3PRIME_IVT_ID",listName="demoList1", listType="Gene")

data_2 <- read.csv("Results2/DEG/DEG_all_integrated9.csv")
USPG1 <- subset(data_2 , cluster== "Undiff. Spermatogonia1" & avg_logFC > 0.5)
USPG1$gene
length(USPG1$gene)
qw
######################################Monocle2######################################
data <- as(as.matrix(all.integrated@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = all.integrated@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              #lowerDetectionLimit = 0.5,
                              expressionFamily = uninormal())
#View data
head(pData(monocle_cds))
head(fData(monocle_cds))
#Run ordering algorithm
var_genes <- all.integrated[["RNA"]]@var.features
ordering_genes <- var_genes
monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
print(dim(exprs(monocle_cds)))
remove(all.integrated)
## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
monocle_cds <- reduceDimension(monocle_cds,norm_method="none", 
                               reduction_method="DDRTree",
                               max_components=2,
                               scaling=TRUE,
                               verbose=TRUE,
                               pseudo_expr=0)
###
# First decide what you want to color your cells by
print(head(pData(monocle_cds)))

## order cells change colors and theta to match your plot
monocle_cds <- orderCells(monocle_cds)

GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$seurat_clusters)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

monocle_cds <- orderCells(monocle_cds, root_state = GM_state(monocle_cds))

pdf("./Results2/trajectory.pdf", width=10)
plot_cell_trajectory(monocle_cds, 
                     color_by = "seurat_clusters",
                     theta = -15,
                     show_branch_points = FALSE,
                     show_tree = TRUE,
                     cell_size = 4) + theme(legend.position = "right")
plot_cell_trajectory(monocle_cds, 
                     color_by = "seurat_clusters") + theme(legend.position = "right")
plot_cell_trajectory(monocle_cds, 
                     color_by = "State") + theme(legend.position = "State")
plot_cell_trajectory(monocle_cds, 
                     color_by = "Pseudotime") + theme(legend.position = "State")
dev.off()

######################################Monocle3 alpha######################################
monocle_cds <- updateCDS(monocle_cds)
monocle_cds <- partitionCells(monocle_cds)
monocle_cds <- learnGraph(monocle_cds,  RGE_method = 'SimplePPT')
pdf("./Results2/trajectory2_2.pdf", width=10)
plot_cell_trajectory(monocle_cds,
                     color_by = "seurat_clusters",  cell_size = 0.01)+ theme(legend.position = "right")
dev.off()

get_correct_root_state <- function(monocle_cds, cell_phenotype, root_type){
  cell_ids <- which(pData(monocle_cds)[, cell_phenotype] == root_type)
  
  closest_vertex <-
    monocle_cds@auxOrderingData[[monocle_cds@rge_method]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(monocle_cds), ])
  root_pr_nodes <-
    V(monocle_cds@minSpanningTree)$name[as.numeric(names
                                                   (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
MPP_node_ids = get_correct_root_state(monocle_cds,
                                      cell_phenotype = "seurat_clusters", "0")
monocle_cds <- orderCells(monocle_cds, root_pr_nodes = MPP_node_ids)

pdf("./Results2/trajectory3.pdf", width=10)
plot_cell_trajectory(monocle_cds)
dev.off()
######################################Monocle3######################################
setwd("C:/Najmeh/Education/MyProject/scRNAseq/Data/")
load("scRNAseq3_36.RData")
library(monocle3)
library(htmlwidgets)
nPC <-36
cluster.res <- 0.2
Dim = "2D"
input.dir <- ("./Results2/Monocle")
output.dir <- sprintf("./Results2/Monocle")

DefaultAssay(all.integrated) <- "integrated"
Idents(all.integrated) <- cellTypes2
if (Dim = "3D"){
  print ("Running UMAP 3D")
  all.integrated <- RunUMAP(object = all.integrated, reduction = "pca", dims = 1:nPC, n.components = 3)
  print("Clustering 3D")
  all.integrated <- FindNeighbors(object=all.integrated, dims=1:nPC)
  all.integrated <- FindClusters(object=all.integrated, resolution=cluster.res)
  all.integrated[[sprintf("ClusterNames_%.1f_%dPC", cluster.res, nPC)]] <- Idents(object = all.integrated)
}

# part one, gene annotations
gene_annotation <- as.data.frame(rownames(all.integrated@reductions[["pca"]]@feature.loadings), row.names = rownames(all.integrated@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"
head(gene_annotation)

# part two, cell information
cell_metadata <- as.data.frame(all.integrated@assays[["RNA"]]@counts@Dimnames[[2]], row.names = all.integrated@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"
head(cell_metadata)

# part three, counts sparse matrix
New_matrix <- all.integrated@assays[["RNA"]]@data
New_matrix <- New_matrix[rownames(all.integrated@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix
New_matrix[1:10,1:10]

### Construct the basic cds object
cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)

### Construct and assign the made up partition
recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

### Assign the cluster info
list_cluster <- all.integrated@meta.data[[sprintf("ClusterNames_%s_%sPC", cluster.res, nPC)]]
list_cluster <-all.integrated@meta.data$seurat_clusters
list_cluster <-Idents(all.integrated)
names(list_cluster) <- all.integrated@assays[["RNA"]]@data@Dimnames[[2]]
cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

### Could be a space-holder, but essentially fills out louvain parameters
cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
###
#p1<-all.integrated@reductions[["umap"]]@cell.embeddings
#p2 <- p1 [ -2 <p1[,"UMAP_1"] & p1[,"UMAP_1"] < 0,]
#p3 <- p2 [0 <p2[,"UMAP_2"] & p2[,"UMAP_2"] < 2,]
#p4 <- row.names(p3)
#p5 <-p1[!rownames(p1) %in% p4, ]
#p1[rownames(p1) %in% p4,1]=-9
#p1[rownames(p1) %in% p4,1]=4

#p12 <- p1 [ -2 <p1[,"UMAP_1"] & p1[,"UMAP_1"] < 0,]
#p13 <- p12 [-6 <p12[,"UMAP_2"] & p12[,"UMAP_2"] < -4,]
#p14 <- row.names(p13)
#p1[rownames(p1) %in% p14,1]= p1[rownames(p1) %in% p14,1]*0.9
#p1[rownames(p1) %in% p14,2]= p1[rownames(p1) %in% p14,2]*0.9
### Assign UMAP coordinate
#cds_from_seurat@reducedDims@listData[["UMAP"]] <-p1
cds_from_seurat@reducedDims@listData[["UMAP"]] <-all.integrated@reductions[["umap"]]@cell.embeddings

#cds_from_seurat@reducedDims@listData[["UMAP"]] <- cds_from_seurat@reducedDims@listData[["UMAP"]]*-1 

### Assign feature loading for downstream module analysis
cds_from_seurat@preprocess_aux$gene_loadings <- all.integrated@reductions[["pca"]]@feature.loadings


remove(all.integrated)
### Learn graph, this step usually takes a significant period of time for larger samples
print("Learning graph, which can take a while depends on the sample")
cds_from_seurat <- learn_graph(cds_from_seurat)
cds_from_seurat <- learn_graph(cds_from_seurat, verbose = FALSE, 
                               learn_graph_control = list(geodesic_distance_ratio=0.11))


### Plot cluster info with trajectory
print("Plotting clusters")
Dim = "2D"
if (Dim == "2D") {
  pdf(sprintf("%s/clusters.with.trajectory12.%s.pdf", output.dir, Dim), width = 10, height = 10)
  clus <- plot_cells(cds_from_seurat, 
                     color_cells_by = 'cluster',
                     label_groups_by_cluster=TRUE,
                     label_leaves=FALSE,
                     group_label_size=5,
                     label_branch_points=TRUE) + 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) +
    theme(legend.position = "right")
  print(clus)
  dev.off()
}
save.image(file = "monocle3_4.RData")
################################################Monocle3
setwd("C:/Najmeh/Education/MyProject/scRNAseq/Data/")
load("monocle3_4.RData")



root_cells = choose_cells(cds_from_seurat ,reduction_method ="UMAP", return_list = TRUE)
cds_from_seurat <- order_cells(cds_from_seurat, root_cells = root_cells)

cds_from_seurat_sub <- choose_graph_segments(cds_from_seurat)

pdf("./Results2/Monocle/trajectory13_gd.11.pdf", width = 10, height = 10)
plot_cells(cds_from_seurat,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  theme(legend.position = "right")

dev.off()
###
data_1 <- readxl::read_excel("Results2/WGCNA/all_BC.xlsx")
features=c(data_1[1:12,2])
colData(cds_from_seurat)$cell.type <-list_cluster 
my_levels <-  c("Somatic cells", "Undiff. Spermatogonia", "Diff.ing Spermatogonia", "Diff.ed Spermatogonia", "leptotene",
                "Zygotene", "Pachytene","Diplotene","Early round spermatids","Round spermatids","Elongating spermatids")
AFD_lineage_cds <- cds_from_seurat[rowData(cds_from_seurat)$gene_short_name %in% features$Gene,
                       colData(cds_from_seurat)$cell.type %in% my_levels]

pdf("./Results2/Monocle/Stime7.pdf",width = 7)
plot_genes_in_pseudotime(AFD_lineage_cds,ncol = 2,cell_size = 0.5,
                         color_cells_by = "pseudotime",panel_order=features$Gene,
                         min_expr=0.5)
dev.off()
###
data_2 <- readxl::read_excel("Results2/DEG/TopDEGs.xlsx")
features=c(data_2[1:14,2])
colData(cds_from_seurat)$cell.type <-list_cluster 
my_levels <-  c("Somatic cells", "Undiff. Spermatogonia", "Diff.ing Spermatogonia", "Diff.ed Spermatogonia", "leptotene",
                "Zygotene", "Pachytene","Diplotene","Early round spermatids","Round spermatids","Elongating spermatids")
AFD_lineage_cds <- cds_from_seurat[rowData(cds_from_seurat)$gene_short_name %in% features$Gene,
                                   colData(cds_from_seurat)$cell.type %in% my_levels]

pdf("./Results2/Monocle/Stime12.pdf",width = 10)
plot_genes_in_pseudotime(AFD_lineage_cds,ncol = 4,cell_size = 0.5,
                         color_cells_by = "pseudotime",panel_order=features$Gene,
                         min_expr=0.5)
dev.off()
###
data_3 <- readxl::read_excel("Results2/WGCNA/all_BC.xlsx")
features=c(data_3[1:10,2])
colData(cds_from_seurat)$cell.type <-list_cluster 
my_levels <-  c("Somatic cells", "Undiff. Spermatogonia", "Diff.ing Spermatogonia", "Diff.ed Spermatogonia", "leptotene",
                "Zygotene", "Pachytene","Diplotene","Early round spermatids","Round spermatids","Elongating spermatids")
AFD_lineage_cds <- cds_from_seurat[rowData(cds_from_seurat)$gene_short_name %in% features$Gene,
                                   colData(cds_from_seurat)$cell.type %in% my_levels]

pdf("./Results2/Monocle/Stime16.pdf",width = 17, height = 6)
plot_genes_in_pseudotime(AFD_lineage_cds,ncol = 5,cell_size = 0.5,
                         color_cells_by = "pseudotime",panel_order=features$Gene,
                         min_expr=0.5)
dev.off()
###
data_4 <- readxl::read_excel("Results2/WGCNA/blue-brown_BC.xlsx")
features=c(data_4[1:10,2])
colData(cds_from_seurat)$cell.type <-list_cluster 
my_levels <-  c("Somatic cells", "Undiff. Spermatogonia", "Diff.ing Spermatogonia", "Diff.ed Spermatogonia", "leptotene",
                "Zygotene", "Pachytene","Diplotene","Early round spermatids","Round spermatids","Elongating spermatids")
AFD_lineage_cds <- cds_from_seurat[rowData(cds_from_seurat)$gene_short_name %in% features$Gene,
                                   colData(cds_from_seurat)$cell.type %in% my_levels]

pdf("./Results2/Monocle/Stime17.pdf",width = 17, height = 6)
plot_genes_in_pseudotime(AFD_lineage_cds,ncol = 5,cell_size = 0.5,
                         color_cells_by = "pseudotime",panel_order=features$Gene,
                         min_expr=0.5)
dev.off()
###
data_5 <- readxl::read_excel("Results2/WGCNA/blue-turquoise-yellow_BC.xlsx")
features=c(data_5[1:4,2])
colData(cds_from_seurat)$cell.type <-list_cluster 
my_levels <-  c("Somatic cells", "Undiff. Spermatogonia", "Diff.ing Spermatogonia", "Diff.ed Spermatogonia", "leptotene",
                "Zygotene", "Pachytene","Diplotene","Early round spermatids","Round spermatids","Elongating spermatids")
AFD_lineage_cds <- cds_from_seurat[rowData(cds_from_seurat)$gene_short_name %in% features$Gene,
                                   colData(cds_from_seurat)$cell.type %in% my_levels]

pdf("./Results2/Monocle/Stime18.pdf",width = 17, height = 6)
plot_genes_in_pseudotime(AFD_lineage_cds,ncol = 5,cell_size = 0.5,
                         color_cells_by = "pseudotime",panel_order=features$Gene,
                         min_expr=0.5)
dev.off()
features=c("STAG3", "REC8", "RAD21L1","RAD21")
colData(cds_from_seurat)$cell.type <-list_cluster 
my_levels <-  c("Somatic cells", "Undiff. Spermatogonia", "Diff.ing Spermatogonia", "Diff.ed Spermatogonia", "leptotene",
                "Zygotene", "Pachytene","Diplotene","Early round spermatids","Round spermatids","Elongating spermatids")
AFD_lineage_cds <- cds_from_seurat[rowData(cds_from_seurat)$gene_short_name %in% features,
                                   colData(cds_from_seurat)$cell.type %in% my_levels]

pdf("./Results2/MonSTAG3.pdf",width = 7)
plot_genes_in_pseudotime(AFD_lineage_cds,ncol = 2,cell_size = 0.5,
                         color_cells_by = "pseudotime",panel_order=features,
                         min_expr=0.5)
dev.off()
avd
#############################################WGCNA_all######################################
setwd("C:/Najmeh/Education/MyProject/scRNAseq/Data/")
library(WGCNA)
load("scRNAseq3_36.RData")  
DefaultAssay(all.integrated) <- "RNA"
DefaultAssay(all.integrated) <- "integrated"

library(doParallel)
registerDoParallel(cores=4)
library(WGCNA)
options(stringsAsFactors = F)

#nfeatures = 2000
datExpr <- t(as.matrix(GetAssayData(all.integrated)))
dim(datExpr)
head(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

#Save for all data in IPM server
#save.image(file = "scRNAseq5.RData")

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
pdf("./Results2/WGCNA/Scale-free_7.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
####
#power = 6
net <- blockwiseModules(datExpr, power = 6,
                        corType = "pearson", 
                        networkType = "signed", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = F, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "TOM_all",
                        verbose = 3)
#Big data
bwnet = blockwiseModules(datExpr, maxBlockSize = 2000,
                         power = 8, TOMType = "signed", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = F,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "TOM_all",
                         verbose = 3)


table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
#save.image(file = "scRNAseq6.RData")
# Convert labels to colors for plotting
mergedColors = net$colors
# Plot the dendrogram and the module colors underneath
pdf("./Results2/WGCNA/dendrogram_pearson_7.pdf")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

df <- data.frame(grey=colnames(datExpr)[net$colors == "grey"]) #notice I created a data.frame vs. the vector you called
blue=colnames(datExpr)[net$colors == "blue"]
brown=colnames(datExpr)[net$colors == "brown"]
green=colnames(datExpr)[net$colors == "green"]
turquoise=colnames(datExpr)[net$colors == "turquoise"]
red=colnames(datExpr)[net$colors == "red"]
yellow=colnames(datExpr)[net$colors == "yellow"]

#METHOD 1
df$blue <- c(blue, rep(NA, nrow(df)-length(blue)))#keep as integer
df$brown <- c(brown, rep(NA, nrow(df)-length(brown)))
df$green <- c(green, rep(NA, nrow(df)-length(green)))
df$turquoise <- c(turquoise, rep(NA, nrow(df)-length(turquoise)))
df$red <- c(red, rep(NA, nrow(df)-length(red)))
df$yellow <- c(yellow, rep(NA, nrow(df)-length(yellow)))
str(df) 

write.csv(df, "./Results2/WGCNA/coE_Pearson_7.csv")

grey.eigengene <- unlist(net$MEs[paste0("ME", "grey")])
names(grey.eigengene) <- rownames(datExpr)
all.integrated$grey.eigengene <- grey.eigengene

blue.eigengene <- unlist(net$MEs[paste0("ME", "blue")])
names(blue.eigengene) <- rownames(datExpr)
all.integrated$blue.eigengene <- blue.eigengene

brown.eigengene <- unlist(net$MEs[paste0("ME", "brown")])
names(brown.eigengene) <- rownames(datExpr)
all.integrated$brown.eigengene <- brown.eigengene

green.eigengene <- unlist(net$MEs[paste0("ME", "green")])
names(green.eigengene) <- rownames(datExpr)
all.integrated$green.eigengene <- green.eigengene

turquoise.eigengene <- unlist(net$MEs[paste0("ME", "turquoise")])
names(turquoise.eigengene) <- rownames(datExpr)
all.integrated$turquoise.eigengene <- turquoise.eigengene

red.eigengene <- unlist(net$MEs[paste0("ME", "red")])
names(red.eigengene) <- rownames(datExpr)
all.integrated$red.eigengene <- red.eigengene

yellow.eigengene <- unlist(net$MEs[paste0("ME", "yellow")])
names(yellow.eigengene) <- rownames(datExpr)
all.integrated$yellow.eigengene <- yellow.eigengene

net$colors


pdf("./Results2/WGCNA/FeaturePlot_7.pdf")
FeaturePlot(all.integrated, features = "blue.eigengene" ,cols = c("grey", "blue"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
FeaturePlot(all.integrated, features ="brown.eigengene" ,cols = c("grey", "brown"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
FeaturePlot(all.integrated, features ="green.eigengene",cols = c("grey","green"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
FeaturePlot(all.integrated, features = "turquoise.eigengene" ,cols = c("grey", "turquoise"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
FeaturePlot(all.integrated, features = "red.eigengene" ,cols = c("grey","red"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
FeaturePlot(all.integrated, features = "yellow.eigengene" ,cols = c("grey","yellow"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
FeaturePlot(all.integrated, features = "grey.eigengene" ,cols = c("grey","Orange"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

modNames = substring(names(MEs), 3)

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
plotTOM = dissTOM^7;
diag(plotTOM) = NA;
sizeGrWindow(9,9)

pdf("./Results2/WGCNA/heatmap_6.pdf")
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()

# Recalculate module eigengenes
#MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
MET =MEs
sizeGrWindow(5,7.5);
par(cex = 0.9)
pdf("./Results2/WGCNA/Eigengene_6.pdf")
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
dev.off()

sizeGrWindow(6,6);
par(cex = 1.0)
pdf("./Results2/WGCNA/Eigengene_dendrogram_6.pdf")
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()

par(cex = 1.0)
pdf("./Results2/WGCNA/Eigengene_heatmap_6.pdf")
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

j1 <- max.col(p1, "first")
WGCNA <- names(p1)[j1]
Idents(all.integrated) <- WGCNA

pdf("./Results2/WGCNA/FeaturePlot_6.pdf")
DimPlot(all.integrated, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE, 
        cols= c("blue","green" ,"grey" , "brown", "turquoise", "red" , "yellow"))+ NoLegend()
dev.off()

f <- function(module){
  eigengene <- unlist(net$MEs[paste0("ME", module)])
  means <- tapply(eigengene, Idents(all.integrated), mean, na.rm = T)
  return(means)
}
modules <- c("blue", "brown", "green", "turquoise", "red","yellow", "")

Idents(all.integrated)<- class2
pdf("./Results2/WGCNA/clusters_WGCNA_7.pdf")
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = "Seurat Cluster",
        ylab = "WGCNA Module Eigengene")
axis(1, at = 1:18, labels = 1:18, cex.axis = 0.8)
matpoints(plotdat, col = modules, pch = 21)
dev.off()

save.image(file = "scRNAseq5.RData")
load("scRNAseq5.RData")  
DefaultAssay(all.integrated) <- "integrated"

TOM = TOMsimilarityFromExpr(datExpr, power = 6);
modules <- c( "blue", "turquoise", "yellow")
probes = colnames(datExpr)
inModule = is.finite(match(moduleLabels , modules));
modProbes = probes[inModule];
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
moduleColors=moduleLabels
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("./Results2/WGCNA/CytoscapeInput-edges-2", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("./Results2/WGCNA/CytoscapeInput-nodes-2", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

sdds

######################################tmp_WGCNA######################################
setwd("C:/Najmeh/Education/MyProject/scRNAseq/Data/")
load("scRNAseq4.RData") 
library(WGCNA)
options(stringsAsFactors = F)
#spermatogonia1
spermatogonia1  <- FindVariableFeatures(spermatogonia1 , selection.method = "vst", nfeatures = 2000)
datExpr <- t(as.matrix(GetAssayData(spermatogonia1)))[,VariableFeatures(spermatogonia1 )] 
dim(datExpr)
head(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

#########Save for all data in IPM server
#save.image(file = "scRNAseq5.RData")

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
pdf("./Results2/WGCNA/Scale-free_spermatogonia1.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.96,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#power = 6
net <- blockwiseModules(datExpr, power = 6,
                        corType = "pearson", # use robust correlation
                        networkType = "signed", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = F, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "TOM_spermatogonia1",
                        verbose = 3)
#Big data
bwnet = blockwiseModules(datExpr, maxBlockSize = 2000,
                         power = 8, TOMType = "signed", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = F,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "TOM_all",
                         verbose = 3)


table(net$colors)
#save.image(file = "scRNAseq6.RData")
# Convert labels to colors for plotting
mergedColors = net$colors
# Plot the dendrogram and the module colors underneath
pdf("./Results2/WGCNA/dendrogram_pearson_spermatogonia1.pdf")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

df <- data.frame(grey=colnames(datExpr)[net$colors == "grey"]) #notice I created a data.frame vs. the vector you called
blue=colnames(datExpr)[net$colors == "blue"]
brown=colnames(datExpr)[net$colors == "brown"]
#green=colnames(datExpr)[net$colors == "green"]
turquoise=colnames(datExpr)[net$colors == "turquoise"]
#red=colnames(datExpr)[net$colors == "red"]
yellow=colnames(datExpr)[net$colors == "yellow"]

#METHOD 1
df$blue <- c(blue, rep(NA, nrow(df)-length(blue)))#keep as integer
df$brown <- c(brown, rep(NA, nrow(df)-length(brown)))
#df$green <- c(green, rep(NA, nrow(df)-length(green)))
df$turquoise <- c(turquoise, rep(NA, nrow(df)-length(turquoise)))
#df$red <- c(red, rep(NA, nrow(df)-length(red)))
df$yellow <- c(yellow, rep(NA, nrow(df)-length(yellow)))
str(df) 

write.csv(df, "./Results2/WGCNA/coE_Pearson_spermatogonia1.csv")

grey.eigengene <- unlist(net$MEs[paste0("ME", "grey")])
names(grey.eigengene) <- rownames(datExpr)
all.integrated$grey.eigengene <- grey.eigengene

blue.eigengene <- unlist(net$MEs[paste0("ME", "blue")])
names(blue.eigengene) <- rownames(datExpr)
all.integrated$blue.eigengene <- blue.eigengene

brown.eigengene <- unlist(net$MEs[paste0("ME", "brown")])
names(brown.eigengene) <- rownames(datExpr)
all.integrated$brown.eigengene <- brown.eigengene

#green.eigengene <- unlist(net$MEs[paste0("ME", "green")])
#names(green.eigengene) <- rownames(datExpr)
#all.integrated$green.eigengene <- green.eigengene

turquoise.eigengene <- unlist(net$MEs[paste0("ME", "turquoise")])
names(turquoise.eigengene) <- rownames(datExpr)
all.integrated$turquoise.eigengene <- turquoise.eigengene

#red.eigengene <- unlist(net$MEs[paste0("ME", "red")])
#names(red.eigengene) <- rownames(datExpr)
#all.integrated$red.eigengene <- red.eigengene

yellow.eigengene <- unlist(net$MEs[paste0("ME", "yellow")])
names(yellow.eigengene) <- rownames(datExpr)
all.integrated$yellow.eigengene <- yellow.eigengene

net$colors


pdf("./Results2/WGCNA/FeaturePlot_spermatogonia1.pdf")
FeaturePlot(all.integrated, features = "blue.eigengene" ,cols = c("grey", "blue"))
FeaturePlot(all.integrated, features ="brown.eigengene" ,cols = c("grey", "brown"))
FeaturePlot(all.integrated, features ="green.eigengene",cols = c("grey","green"))
FeaturePlot(all.integrated, features = "turquoise.eigengene" ,cols = c("grey", "turquoise"))
FeaturePlot(all.integrated, features = "red.eigengene" ,cols = c("grey","red"))
FeaturePlot(all.integrated, features = "yellow.eigengene" ,cols = c("grey","yellow"))
FeaturePlot(all.integrated, features = "grey.eigengene" ,cols = c("grey","Orange"))
dev.off()

p1<-net$MEs
j1 <- max.col(p1, "first")
WGCNA <- names(p1)[j1]
Idents(all.integrated) <- WGCNA

pdf("./Results2/WGCNA/FeaturePlot_5_2.pdf")
DimPlot(all.integrated, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE, 
        cols= c("blue","green" ,"grey" , "brown", "turquoise", "red" , "yellow"))+ NoLegend()
dev.off()

f <- function(module){
  eigengene <- unlist(net$MEs[paste0("ME", module)])
  means <- tapply(eigengene, Idents(all.integrated), mean, na.rm = T)
  return(means)
}
modules <- c("blue", "brown", "green", "turquoise", "red","yellow")

pdf("./Results2/WGCNA/clusters_WGCNA_5.pdf")
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = "Seurat Cluster",
        ylab = "WGCNA Module Eigengene")
axis(1, at = 1:19, labels = 0:18, cex.axis = 0.8)
matpoints(plotdat, col = modules, pch = 21)
dev.off()

spermatogoniaExpr2 <-t(as.matrix(GetAssayData(spermatogonia2)))
dim(spermatogoniaExpr2)
SpermatocytesExpr <-t(as.matrix(GetAssayData(Spermatocytes)))
dim(SpermatocytesExpr)
SpermatidsExpr <-t(as.matrix(GetAssayData(Spermatids)))
dim(SpermatidsExpr)
#b
spermatogonia <- FindVariableFeatures(spermatogonia, selection.method = "vst", nfeatures = 2000)
spermatogoniaExpr <-t(as.matrix(GetAssayData(spermatogonia)))[,VariableFeatures(spermatogonia)] 
dim(spermatogoniaExpr)

spermatogonia1 <- FindVariableFeatures(spermatogonia1, selection.method = "vst", nfeatures = 2000)
spermatogoniaExpr1 <-t(as.matrix(GetAssayData(spermatogonia1)))[,VariableFeatures(spermatogonia1)] 
dim(spermatogoniaExpr1)

spermatogonia2 <- FindVariableFeatures(spermatogonia2, selection.method = "vst", nfeatures = 2000)
spermatogoniaExpr2 <-t(as.matrix(GetAssayData(spermatogonia2)))[,VariableFeatures(spermatogonia2)] 
dim(spermatogoniaExpr2)


Spermatocytes <- FindVariableFeatures(Spermatocytes, selection.method = "vst", nfeatures = 2000)
SpermatocytesExpr <-t(as.matrix(GetAssayData(Spermatocytes)))[,VariableFeatures(Spermatocytes)] 
dim(SpermatocytesExpr)

Spermatids <- FindVariableFeatures(Spermatids, selection.method = "vst", nfeatures = 2000)
SpermatidsExpr <-t(as.matrix(GetAssayData(Spermatids)))[,VariableFeatures(Spermatids)] 
dim(SpermatidsExpr)

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(spermatogoniaExpr, powerVector = powers, verbose = 5)
sft1 = pickSoftThreshold(spermatogoniaExpr1, powerVector = powers, verbose = 5)
sft2 = pickSoftThreshold(spermatogoniaExpr2, powerVector = powers, verbose = 5)
sft3 = pickSoftThreshold(SpermatocytesExpr, powerVector = powers, verbose = 5)
sft4 = pickSoftThreshold(SpermatidsExpr, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

pdf("./Results2/Scale-free_variable2.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

plot(sft1$fitIndices[,1], -sign(sft1$fitIndices[,3])*sft1$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft1$fitIndices[,1], -sign(sft1$fitIndices[,3])*sft1$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft1$fitIndices[,1], sft1$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft1$fitIndices[,1], sft1$fitIndices[,5], labels=powers, cex=cex1,col="red")

plot(sft2$fitIndices[,1], -sign(sft2$fitIndices[,3])*sft2$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft2$fitIndices[,1], -sign(sft2$fitIndices[,3])*sft2$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft2$fitIndices[,1], sft2$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft2$fitIndices[,1], sft2$fitIndices[,5], labels=powers, cex=cex1,col="red")

plot(sft3$fitIndices[,1], -sign(sft3$fitIndices[,3])*sft3$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft3$fitIndices[,1], -sign(sft3$fitIndices[,3])*sft3$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft3$fitIndices[,1], sft3$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft3$fitIndices[,1], sft3$fitIndices[,5], labels=powers, cex=cex1,col="red")

plot(sft4$fitIndices[,1], -sign(sft4$fitIndices[,3])*sft4$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft4$fitIndices[,1], -sign(sft4$fitIndices[,3])*sft4$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft4$fitIndices[,1], sft4$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft4$fitIndices[,1], sft4$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


net1 <- blockwiseModules(spermatogoniaExpr1, power = 3,
                        corType = "pearson", # use robust correlation
                        networkType = "signed", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "TOM_spermatogonia",
                        verbose = 3)
table(net1$colors)
mergedColors = net1$colors

net2 <- blockwiseModules(spermatogoniaExpr2, power = 5,
                        corType = "pearson", # use robust correlation
                        networkType = "signed", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "TOM_spermatogonia2",
                        verbose = 3)
table(net2$colors)
mergedColors2 = net2$colors

net3 <- blockwiseModules(SpermatocytesExpr, power = 2,
                        corType = "pearson", # use robust correlation
                        networkType = "signed", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "TOM_Spermatocytes",
                        verbose = 3)
table(net3$colors)
mergedColors3 = net3$colors

net4 <- blockwiseModules(SpermatidsExpr, power = 2,
                        corType = "pearson", # use robust correlation
                        networkType = "signed", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "TOM_Spermatids",
                        verbose = 3)
table(net4$colors)
mergedColors4 = net4$colors
pdf("./Results2/dendrogram_pearson2.pdf")
plotDendroAndColors(net1$dendrograms[[1]], mergedColors[net1$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
plotDendroAndColors(net2$dendrograms[[1]], mergedColors2[net2$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
plotDendroAndColors(net3$dendrograms[[1]], mergedColors3[net3$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
plotDendroAndColors(net4$dendrograms[[1]], mergedColors4[net4$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#f <- function(module){
#  eigengene <- unlist(net$MEs[paste0("ME", module)])
#  means <- tapply(eigengene, Idents(all.integrated), mean, na.rm = T)
#  return(means)
#}
#modules <- c("blue", "brown", "green", "turquoise", "yellow")
#plotdat <- sapply(modules, f)
#matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = "Seurat Cluster",
#        ylab = "WGCNA Module Eigengene")
#axis(1, at = 1:19, labels = 0:18, cex.axis = 0.8)
#matpoints(plotdat, col = modules, pch = 21)
#################################################disaese#################################################
library("readxl")
library(data.table)
setwd("C:/Najmeh/Education/MyProject/scRNAseq/Data/")
load("scRNAseq3.RData")
DefaultAssay(all.integrated) <- "RNA"
Idents(all.integrated) <- class 
###DEGall
DEGall <- fread("./Results2/DEG/DEG_all_integrated8.csv",header=T, data.table = F)
colnames(DEGall) <- c("Symbol","p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster","gene")
dim(DEGall)
DEGall<- subset(DEGall,avg_logFC > 0.5 & p_val_adj < 0.05)
dim(DEGall)

#WGCNAall
WGCNAall <- fread("./Results2/WGCNA/coE_Pearson_5.csv",header=T, data.table = T)

#################oligospermia
#OMIM
oligospermia1 <- read_xlsx("./Results2/Disease/Oligospermia or oligozoospermia or oligoasthenospermia/OMIM-Oligospermia or oligozoospermia.xlsx", skip=4, col_names=T)#, col_types=c("text"), range=cell_cols("F")#)
lapply(c(unique(oligospermia1[,3]),unique(oligospermia1[,6])),length)
OligoSym <-c(unique(oligospermia1[,6]))
OligoSym <- na.omit(OligoSym$'Approved Symbol')
length(OligoSym)


#GAD
oligospermia_GAD <- read_xlsx("./Results2/Disease/Oligospermia or oligozoospermia or oligoasthenospermia/GAD.xlsx", col_names=T)  
OligoSym_GAD <-c(unique(oligospermia_GAD[,"GENE"]))
OligoSym_GAD <-na.omit(OligoSym_GAD$"GENE")
length(OligoSym_GAD )

#ClinVar
oligospermia_CV <- read_xlsx("./Results2/Disease/Oligospermia or oligozoospermia or oligoasthenospermia/ClinVar.xlsx", col_names=T)  
oligospermia_CV <- c(unique(oligospermia_CV[,"Gene(s)"]))
oligospermia_CV <-na.omit(oligospermia_CV$"Gene(s)")
length(oligospermia_CV)

#total
p1<- unique(c(OligoSym_GAD,OligoSym))
length(p1)
length(intersect(OligoSym_GAD,OligoSym))

#Disease-DEG
p2 <-intersect(DEGall[,8],p1)
p3 <-subset( DEGall, gene %in% p2)
length(row.names(p3))
length(unique(p4[,"gene"]))
write.csv(p3,"./Results2/Disease/Oligospermia or oligozoospermia or oligoasthenospermia/Oligo_DEG.csv")

#Disease-WGCNA
intersect(WGCNAall$grey,p1)
#8
intersect(WGCNAall$blue,p1)
#7
intersect(WGCNAall$brown ,p1)
#2
intersect(WGCNAall$green  ,p1)
#1
intersect(WGCNAall$turquoise  ,p1)
#8
intersect(WGCNAall$red   ,p1)
#1
intersect(WGCNAall$yellow   ,p1)
#2


#################Azoospermia
#OMIM
Azoospermia <- read_xlsx("./Results2/Disease/Azoospermia/OMIM-Gene-Map-Retrieval.xlsx", skip=4, col_names=T)#, col_types=c("text"), range=cell_cols("F")#)
lapply(c(unique(Azoospermia[,3]),unique(Azoospermia[,6])),length)
AzooSym <-c(unique(Azoospermia[,6]))
AzooSym <- na.omit(AzooSym$'Approved Symbol')
length(AzooSym)


#GAD
Azoospermia_GAD <- read_xlsx("./Results2/Disease/Azoospermia/GAD.xlsx", col_names=T)  
AzooSym_GAD <-c(unique(Azoospermia_GAD[,"GENE"]))
AzooSym_GAD <-na.omit(AzooSym_GAD$"GENE")
length(AzooSym_GAD )

#total
p1<- unique(c(AzooSym_GAD ,AzooSym))
length(p1)
length(intersect(AzooSym_GAD ,AzooSym))

#Disease-DEG
p2 <-intersect(DEGall[,8],p1)
p3 <-subset( DEGall, gene %in% p2)
length(row.names(p3))
length(unique(p3[,"gene"]))
write.csv(p3,"./Results2/Disease/Azoospermia/Azoo_DEG.csv")


Asthenozoospermia 
#################Asthenozoospermia 
#OMIM
Asthenozoospermia <- read_xlsx("./Results2/Disease/Asthenozoospermia or asthenospermia or oligoasthenospermia/OMIM-Gene-Map-Retrieval.xlsx", skip=4, col_names=T)#, col_types=c("text"), range=cell_cols("F")#)
lapply(c(unique(Asthenozoospermia[,3]),unique(Asthenozoospermia[,6])),length)
AsthSym <-c(unique(Asthenozoospermia[,6]))
AsthSym <- na.omit(AsthSym$'Approved Symbol')
length(AsthSym)


#GAD
Asthenozoospermia_GAD <- read_xlsx("./Results2/Disease/Asthenozoospermia or asthenospermia or oligoasthenospermia/GAD.xlsx", col_names=T)  
AsthSym_GAD <-c(unique(Asthenozoospermia_GAD[,"GENE"]))
AsthSym_GAD <-na.omit(AsthSym_GAD$"GENE")
length(AsthSym_GAD )

#total
p1<- unique(c(AsthSym_GAD ,AsthSym))
length(p1)
length(intersect(AsthSym_GAD ,AsthSym))

#Disease-DEG
p2 <-intersect(DEGall[,8],p1)
p3 <-subset( DEGall, gene %in% p2)
length(row.names(p3))
length(unique(p3[,"gene"]))
write.csv(p3,"./Results2/Disease/Asthenozoospermia or asthenospermia or oligoasthenospermia/Asth_DEG.csv")

#################male infertility 
#OMIM
MInfertility <- read_xlsx("./Results2/Disease/male infetility/OMIM-Gene-Map-Retrieval1.xlsx", skip=4, col_names=T)#, col_types=c("text"), range=cell_cols("F")#)
lapply(c(unique(MInfertility[,3]),unique(MInfertility[,6])),length)
MISym <-c(unique(MInfertility[,6]))
MISym <- na.omit(MISym$'Approved Symbol')
length(MISym)


#GAD
MInfertility_GAD <- read_xlsx("./Results2/Disease/male infetility/GAD1.xlsx", col_names=T)  
MI_GAD <-c(unique(MInfertility_GAD[,"GENE"]))
MI_GAD <-na.omit(MI_GAD$"GENE")
length(MI_GAD )

#total
p1<- unique(c(MI_GAD ,MISym))
length(p1)
length(intersect(MI_GAD ,MISym))

#Disease-DEG
p2 <-intersect(DEGall[,8],p1)
p3 <-subset( DEGall, gene %in% p2)
length(row.names(p3))
length(unique(p3[,"gene"]))
write.csv(p3,"./Results2/Disease/male infetility/MI_DEG2.csv")

######################################Co-expression Network######################################
data_1 <- readxl::read_excel("Results2/WGCNA/blue-brown_BC.xlsx")
features=c(data_1[1:10,2])

RidgePlot(all.integrated, features = features$Gene, ncol = 2)
Idents(all.integrated)<-cellTypes3

my_levels <-  c("Somatic cells","Undiff. SPG1", "Undiff. SPG2","Diff.ing SPG1","Diff.ing SPG2", "Diff.ed SPG","Leptotene SPC", "Pachytene SPC", "Zygotene SPC", "Diplotene SPC", "Early round SPT", "Round SPT1", "Round SPT2", "Elongating SPT")
# Re-level object@ident
all.integrated@active.ident <- factor(x = all.integrated@active.ident, levels = my_levels)

pdf("./Results2/WGCNA/TopBet6.pdf", width=15)
FeaturePlot(all.integrated, features =features$Gene)
VlnPlot(all.integrated, features$Gene)
dev.off()
pdf("./Results2/WGCNA/TopBet14.pdf", width=15)
DotPlot(all.integrated, features = features$Gene) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
