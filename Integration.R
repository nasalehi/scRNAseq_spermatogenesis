---
### integration using Seurat3.2 standard integration flow
### outputs .rds file containing the integrated seurat object.
---

library(dplyr)
library(Seurat)
setwd("....")


######################################Initial Data######################################
###Loading data
spermatogonia.data <- Read10X(data.dir = "GSE109037/AdultHumanSpermatogonia-reps1-3/")
Spermatocytes.data <- Read10X(data.dir = "GSE109037/AdultHuman-Spermatocytes_Reps1-2")
Spermatids.data <- Read10X(data.dir = "GSE109037/AdultHuman-Spermatids_Reps1-2")
Spermatogenesis1.data <- Read10X(data.dir = "GSE109037/AdultHumanSpermatogenesis-reps1-3")
Spermatogenesis2.data <- read.delim("GSE106487/Spermatogenesis2.txt")
dim(spermatogonia.data)
#33694 11104
dim(Spermatocytes.data)
#33694  4884
dim(Spermatids.data)
#33694  7434
dim(Spermatogenesis1.data)
#33694  7434
dim(Spermatogenesis2.data)
#24153  3059
spermatogonia <- CreateSeuratObject(counts = spermatogonia.data, project = "spermatogonia", min.cells = 3, min.features = 200)
Spermatocytes <- CreateSeuratObject(counts = Spermatocytes.data, project = "Spermatocytes", min.cells = 3, min.features = 200)
Spermatids <- CreateSeuratObject(counts = Spermatids.data, project = "Spermatids", min.cells = 3, min.features = 200)
Spermatogenesis1 <- CreateSeuratObject(counts = Spermatogenesis1.data , project = "Spermatogenesis1", min.cells = 3, min.features = 200)
Spermatogenesis2 <- CreateSeuratObject(counts = Spermatogenesis2.data , project = "Spermatogenesis2", min.cells = 3, min.features = 200)
dim(spermatogonia)
#28320 11104
dim(Spermatocytes)
#27405  4884
dim(Spermatids)
#27898  7434
dim(Spermatogenesis1)
#29119  7132
dim(Spermatogenesis2)
#22911  3046
spermatogonia[["percent.mt"]] <- PercentageFeatureSet(spermatogonia, pattern = "^MT-")
Spermatocytes[["percent.mt"]] <- PercentageFeatureSet(Spermatocytes, pattern = "^MT-")
Spermatids[["percent.mt"]] <- PercentageFeatureSet(Spermatids, pattern = "^MT-")
Spermatogenesis1[["percent.mt"]] <- PercentageFeatureSet(Spermatogenesis1, pattern = "^MT-")
Spermatogenesis2[["percent.mt"]] <- PercentageFeatureSet(Spermatogenesis2, pattern = "^MT-")
head(spermatogonia@meta.data, 5)
head(Spermatocytes@meta.data, 5)
head(Spermatids@meta.data, 5)
head(Spermatogenesis1@meta.data, 5)
head(Spermatogenesis2@meta.data, 5)
######################################QC######################################
pdf("./quality1_vlnplot.pdf", width=20)
lapply(c(Spermatogenesis2,Spermatogenesis1,spermatogonia,Spermatocytes,Spermatids),VlnPlot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

spermatogonia <- subset(spermatogonia, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 20)
dim(spermatogonia)
#28320 10998

Spermatocytes <- subset(Spermatocytes, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 10)
dim(Spermatocytes)
#27405  4797

Spermatids <- subset(Spermatids, subset = nFeature_RNA > 1200 & nFeature_RNA < 9000 & percent.mt < 15)
dim(Spermatids)
#27898  7382

Spermatogenesis1 <- subset(Spermatogenesis1, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 25)
dim(Spermatogenesis)
#29119  7037

Spermatogenesis2 <- subset(Spermatogenesis2, subset = nFeature_RNA > 800 & nFeature_RNA < 12000)
dim(Spermatogenesis2)

pdf("./quality2_vlnplot.pdf", width=20)
lapply(c(Spermatogenesis2,Spermatogenesis1,spermatogonia,Spermatocytes,Spermatids),VlnPlot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

######################################Normalization######################################
all.list <- c(Spermatogenesis1,spermatogonia,Spermatocytes, Spermatids,Spermatogenesis2)
names(all.list) <- c("Spermatogenesis1", "spermatogonia","Spermatocytes", "Spermatids", "Spermatogenesis2")

for (i in names(all.list)) {
  all.list[[i]] <- NormalizeData(all.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
}

######################################Variable############################################

for (i in names(all.list)) {
  all.list[[i]] <- FindVariableFeatures(all.list[[i]], selection.method = "vst", nfeatures = 2000)
}


######################################INTEGRATION###########################################
n=35
all.anchors <- FindIntegrationAnchors(object.list=all.list,dims = 1:n)
all.integrated <- IntegrateData(anchorset = all.anchors, dims = 1:n)
dim(all.integrated )
######################################Saving###############################################
save.image(file = "integratedData.RData")
