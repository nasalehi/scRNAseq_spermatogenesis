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

setwd("....")


######################################Dimention Reduction& Clustering##########################################
load("integratedData.RData")
###DR
DefaultAssay(all.integrated) <- "integrated"
all.integrated <- ScaleData(all.integrated, verbose = FALSE)
all.integrated <- RunPCA(all.integrated, npcs = 50, verbose = FALSE)
all.integrated <- RunUMAP(all.integrated, reduction = "pca", dims = 1:35)
dim(all.integrated)
dataIDs <- Idents(all.integrated)
###Clustering
all.integrated <- FindNeighbors(all.integrated, dims = 1:35)
all.integrated <- FindClusters(all.integrated, resolution = 0.2)
table(Idents(all.integrated))

current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
new.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,16)
Idents(all.integrated) <- plyr::mapvalues(x = Idents(all.integrated), from = current.cluster.ids, to = new.cluster.ids)
classIDs<- Idents(all.integrated)

pdf("./integrated_UMAP.pdf", width=10)
Idents(all.integrated) <- dataIDs
DimPlot(all.integrated, reduction = "umap", group.by = 'ident', repel = TRUE, 
        cols= c( "blue3","lightskyblue","thistle","darkseagreen","lightpink1"),
        order = c("Spermatogonia", "Spermatocyte", "Spermatid", "Spermatogenesis1","Spermatogenesis2" ))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))

Idents(all.integrated) <- classIDs
DimPlot(all.integrated, reduction = "umap", group.by = 'ident',label = TRUE, repel = TRUE,label.size = 6)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()


save.image(file = "integrated_UMAP.RData")
#####################################################MARKERS###################################################
load("integrated_UMAP.RData")
DefaultAssay(all.integrated) <- "RNA"

pdf("./Markers_somatic.pdf", width=15)
FeaturePlot(all.integrated, features = c("CYP26B1","INSL3", "MYH11", "ALDH1A1", "CD68","CD163"),cols=c("lightgrey","darkblue")) + RotatedAxis()
VlnPlot(all.integrated, c("CYP26B1","INSL3", "MYH11", "ALDH1A1", "CD68","CD163"))
DotPlot(all.integrated, features =  c("CYP26B1","INSL3", "MYH11", "ALDH1A1", "CD68","CD163"),cols=c("lightgrey","darkblue")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Markers_spermatogonia.pdf", width=12)
FeaturePlot(all.integrated, features = c("NANOS2","PIWIL4","GFRA1","SALL4","MAGEA4","HMGA1" ),cols=c("lightgrey","blue"))
VlnPlot(all.integrated, c("NANOS2","PIWIL4","GFRA1","SALL4","MAGEA4","HMGA1" ))
DotPlot(all.integrated, features = c("NANOS2","PIWIL4","GFRA1","SALL4","MAGEA4","HMGA1" ),cols=c("lightgrey","blue")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Results2/Markers_spermatocyte.pdf", width=15)
FeaturePlot(all.integrated, features = c("DMC1", "RAD51AP2","PIWIL1","SYCP3", "OVOL2"),cols=c("lightgrey","red"))
VlnPlot(all.integrated, c("DMC1", "RAD51AP2","PIWIL1","SYCP3", "OVOL2"))
DotPlot(all.integrated, features = c("DMC1", "RAD51AP2","PIWIL1","SYCP3", "OVOL2"),cols=c("lightgrey","red")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Results2/Markers_spermatid.pdf", width=15)
FeaturePlot(all.integrated, features = c("TEX29","SUN5","SPEM1","SYCP3","OVOL2","ACR","PGK2"),cols=c("lightgrey","darkgreen"))
VlnPlot(all.integrated, c("TEX29","SUN5","SPEM1","SYCP3","OVOL2","ACR","PGK2"))
DotPlot(all.integrated, features = c("TEX29","SUN5","SPEM1","SYCP3","OVOL2","ACR","PGK2"),cols=c("lightgrey","darkgreen")) + RotatedAxis()+
theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

pdf("./Results2/Marker_germcells.pdf", width=15)
FeaturePlot(all.integrated, features = c("DDX4"),cols=c("lightgrey","brown"))
dev.off()
                                  
  
current.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,16)
new.cluster.ids <- c("Undiff. SPG1", "Undiff. SPG2", "Round SPT2", "Round SPT1", "Pachytene SPC", "Zygotene SPC",
                      "Early round SPT", "Diplotene SPC", "Leptotene SPC", "Diff.ing SPG1", "Elongating SPT", 
                     "Diff.ed SPG","Diff.ing SPG2", "Somatic cells", "Somatic cells", "Somatic cells")
Idents(all.integrated) <- plyr::mapvalues(x = Idents(all.integrated), from = current.cluster.ids, to = new.cluster.ids)
cellIDs<-Idents(all.integrated)
my_levels <-  c("Somatic cells","Undiff. SPG1", "Undiff. SPG2","Diff.ing SPG1","Diff.ing SPG2", "Diff.ed SPG","Leptotene SPC", "Pachytene SPC", "Zygotene SPC", "Diplotene SPC", "Early round SPT", "Round SPT1", "Round SPT2", "Elongating SPT")
# Re-level object@ident
all.integrated@active.ident <- factor(x = all.integrated@active.ident, levels = my_levels)                                       


save.image(file = "integrated_UMAP2.RData")
#####################################################DEGs###################################################
load("integrated_UMAP2.RData")
DefaultAssay(all.integrated) <- "RNA"

all.markers <- FindAllMarkers(all.integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(all.markers, file= "./DEGs_integrated.csv")


