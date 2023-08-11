### This is a script for the bone marrow dataset ###

### Road libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(SoupX)
library(stringr)
library(harmony)
library(SingleR)
library(celldex)
library(ggplot2)
library(DoMultiBarHeatmap)


### Set up the Seurat Object and pre-process the data
## 1) Naive
## Load the dataset and remove ambient RNA
# Load the naive dataset
naive_toc = Read10X(data.dir = "./filtered_feature_bc_matrix")
naive_tod = Read10X(data.dir = "./raw_feature_bc_matrix")
# Profile the Soup
naive_sc = SoupChannel(naive_tod, naive_toc, calcSoupProfile = FALSE)
naive_sc = estimateSoup(naive_sc)
# Add extra meta data to the SoupChannel object
naive_metadata <-read.csv(file = "./clusters.csv")
naive_sc = setClusters(naive_sc, setNames(naive_metadata$Cluster, rownames(naive_metadata)))
# Estimate rho
naive_sc = autoEstCont(naive_sc)
# Clean the data
naive_out = adjustCounts(naive_sc)

## Initialize the Seurat object and pre-process the data (quality control) 
# Set up the Seurat object
naive <- CreateSeuratObject(counts = naive_out, project = "0 dpi", min.cells = 5)
naive$group <- "Naive"
# Remove long non-coding RNA contamination
counts.naive <- GetAssayData(naive, assay = "RNA")
counts.naive <- counts.naive[-(which(rownames(counts.naive) %in% c("Gm42418", "AY036118"))),]
naive <- subset(naive, features = rownames(counts.naive))
# Calculate mitochondrial QC metrics
naive[["percent.MT"]] <- PercentageFeatureSet(naive, pattern = "^mt-")
# Visualize QC metrics as a violin plot before filtering low-quality cells
VlnPlot(naive, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
plot1 <- FeatureScatter(naive, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot2 <- FeatureScatter(naive, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# Filter low-quality cells by adjusting these quality metrics: nFeature_RNA and percent.MT
naive <- subset(naive, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.MT < 10)
# Visualize QC metrics as a violin plot after filtering low-quality cells
VlnPlot(naive, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
plot1 <- FeatureScatter(naive, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot2 <- FeatureScatter(naive, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 

## 2) 14dpi
## Load the dataset and remove ambient RNA
# Load the 14dpi dataset
dpi14_toc = Read10X(data.dir = "./filtered_feature_bc_matrix")
dpi14_tod = Read10X(data.dir = "./raw_feature_bc_matrix")
# Profile the Soup
dpi14_sc = SoupChannel(dpi14_tod, dpi14_toc, calcSoupProfile = FALSE)
dpi14_sc = estimateSoup(dpi14_sc)
# Add extra meta data to the SoupChannel object
dpi14_metadata <-read.csv(file = "./clusters.csv")
dpi14_sc = setClusters(dpi14_sc, setNames(dpi14_metadata$Cluster, rownames(dpi14_metadata)))
# Estimate rho
dpi14_sc = autoEstCont(dpi14_sc)
# Clean the data
dpi14_out = adjustCounts(dpi14_sc)

## Initialize the Seurat object and pre-process the data (quality control) 
# Set up the Seurat object
dpi14 <- CreateSeuratObject(counts = dpi14_out, project = "14 dpi", min.cells = 5)
dpi14$group <- "14dpi"
# Remove long non-coding RNA contamination
counts.dpi14 <- GetAssayData(dpi14, assay = "RNA")
counts.dpi14 <- counts.dpi14[-(which(rownames(counts.dpi14) %in% c("Gm42418", "AY036118"))),]
dpi14 <- subset(dpi14, features = rownames(counts.dpi14))
# Calculate mitochondrial QC metrics
dpi14[["percent.MT"]] <- PercentageFeatureSet(dpi14, pattern = "^mt-")
# Visualize QC metrics as a violin plot before filtering low-quality cells
VlnPlot(dpi14, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
plot1 <- FeatureScatter(dpi14, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot2 <- FeatureScatter(dpi14, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 
# Filter low-quality cells by adjusting these quality metrics: nFeature_RNA and percent.MT
dpi14 <- subset(dpi14, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.MT < 10)
# Visualize QC metrics as a violin plot after filtering low-quality cells
VlnPlot(dpi14, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
plot1 <- FeatureScatter(dpi14, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot2 <- FeatureScatter(dpi14, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 

## Merge naive and 14dpi datasets
# Merge two Seurat objects 
evdata <- merge(x= naive, y = dpi14, add.cell.ids = c("Naive", "14dpi"), project = "T.b.b.BM")
# Normalize the counts
evdata <- NormalizeData(evdata) 


### Cell-cycle scoring and regression
# Acquire cell cycle markers from the published paper; Tirosh et al, 2015 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# S phase genes
s.genes = str_to_title(s.genes)
# G2/M phage genes 
g2m.genes = str_to_title(g2m.genes)
# Score cells for cell cycle
evdata <- CellCycleScoring(evdata, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# Regress out the difference between the G2M and S phase scores
evdata$CC.Difference <- evdata$S.Score - evdata$G2M.Score 
# Scale the counts
evdata <- ScaleData(evdata, vars.to.regress = "CC.Difference", features = rownames(evdata))


### Cell clustering
# Identify the most variable genes and perform linear dimensional reduction
evdata <- FindVariableFeatures(evdata) %>% RunPCA(verbose = FALSE) 
# Correct the batch effects
evdata <- evdata %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
options(repr.plot.height = 5, repr.plot.width = 12)
# Display Elbow plot
ElbowPlot(evdata, n = 50)
# Determine the dimensionality of the dataset based on Elbow plot
# Run non-linear dimensional reduction
# Cluster cells
evdata <- evdata %>% 
  RunUMAP(reduction = "harmony", dims = 1:26) %>%  
  FindNeighbors(reduction = "harmony", dims = 1:26) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
options(repr.plot.height = 4, repr.plot.width = 10)
# Visualize individual clusters in a dimensional reduction plot
DimPlot(evdata, reduction = "umap", group.by = "orig.ident", pt.size = .1)
# Visualize individual clusters split by sample group in a dimensional reduction plot
DimPlot(evdata, reduction = "umap",group.by = "orig.ident", split.by = "orig.ident", pt.size = .1)


### Select 1st neutrophil clusters
## 1) Manual annotation of cell type using neutrophil markers identified in the literature
# 1-1) Pre-neutrophil markers
pre_np <- c("Runx1", "Gfi1", "Cebpe",  "Stmn1", "H2afx", "Ube2c", "Tubb5", "Ptma", "Ube2s")         
FeaturePlot(evdata, features = pre_np, label = T)
# 1-2) Immature neutrophil markers-1      
immature_np_1 <- c("Gfi1", "Cebpe", "Orm1", "Lcn2", "Ngp", "Ltf", "Cybb", "Ly6g", "Cd177")         
FeaturePlot(evdata, features = immature_np_1, label = T)
# 1-3) Immature neutrophil markers-2          
immature_np_2 <- c("Thbs1", "Mmp8", "Mmp25", "Mmp9", "Acvrl1")         
FeaturePlot(evdata, features = immature_np_2, label = T)
# 1-4) Mature neutrophil markers      
mature_np <- c("Cebpd", "Spi1", "Il1b", "Ccl6", "Junb", "Jund", "Marcks", "Wfdc17", "Irf7")         
FeaturePlot(evdata, features = mature_np, label = T)
# 1-5) Monocyte markers 
mo <- c("Ccr2", "Cx3cr1")
FeaturePlot(evdata, features = mo, label = T)

## Find all markers of cluster 16 to identify cell type
cluster16.markers <- FindMarkers(evdata, ident.1 = 16, min.pct = 0.25)
write.csv(cluster16.markers, file = "./cluster16.markers.csv")

## 2) Automated cell type annotation
# Load built-in references to annotate cells
evdata.diet = DietSeurat(evdata)
pred.evdata = as.SingleCellExperiment(evdata.diet)
Immgen.red = ImmGenData()
# Perform prediction
pred.evdata.r <- SingleR(test = pred.evdata, ref = Immgen.red, assay.type.test=1,
                         labels = Immgen.red$label.main)
evdata[["SingleR.labels"]] <- pred.evdata.r$labels
# Visualize individual clusters in a dimensional reduction plot
DimPlot(evdata, reduction = "umap", group.by = "SingleR.labels", label = T, label.size=3, repel=T)

## Subset the Seurat object
neutrophil <- subset(evdata, idents = c(0, 1, 2, 3, 4, 16))
# <Description of how clusters were annotated as neutrophils:
# SingleR annotated cluster 0, 1, 2, 3 and 4 as neutrophils.
# However, the top markers belonging to cluster 16, linked to cluster 4, presented by the FindMarkers() function were neutrophil markers.
# So, cluster 16 was also annotated as neutrophils.>

## Cell clustering
# Normalize the counts
# Identify the most variable genes
# Scale the counts
# Perform linear dimensional reduction
all.genes <- rownames(neutrophil)
neutrophil <- NormalizeData(neutrophil) %>% FindVariableFeatures() %>% ScaleData(features=all.genes) %>% RunPCA(verbose = FALSE)
# Display Elbow Plot
ElbowPlot(neutrophil, n=50)
# Display dimensional reduction heatmap
DimHeatmap(neutrophil, dims = 1:10, cells = 500, balanced = TRUE)
# Determine the dimensionality of the dataset based on Elbow plot
# Run non-linear dimensional reduction
# Cluster cells
neutrophil <- neutrophil %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
neutrophil <- neutrophil %>% 
  RunUMAP(reduction = "harmony", dims = 1:7) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:7) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
# Visualize neutrophil clusters in a dimensional reduction plot
DimPlot(neutrophil, label=T)
# Visualize neutrophil clusters split by sample group in a dimensional reduction plot
DimPlot(neutrophil, split.by = "orig.ident", label = T)

## Identify expression levels of the same neutrophil markers which were used for assigning cell type identity in neutrophil clusters
# 1-1) Pre-neutrophil markers
FeaturePlot(neutrophil, features = pre_np, label = T)
# 1-2) Immature neutrophil markers-1          
FeaturePlot(neutrophil, features = immature_np_1, label = T)
# 1-3) Immature neutrophil markers-2 
FeaturePlot(neutrophil, features = immature_np_2, label = T)
# 1-4) Mature neutrophil markers       
FeaturePlot(neutrophil, features = mature_np, label = T)


### Select 2nd neutrophil clusters
# Remove GMP clusters from neutrophil clusters
factor(Idents(neutrophil))
levels(Idents(neutrophil))
neutrophil_2 <- neutrophil
# Identify expression levels of the granulocyte-monocyte progenitor (GMP) markers in neutrophil clusters  
gmp <- c("Mpo", "Elane", "Ctsg", "Prtn3", "Prss57", "Ctsc") 
FeaturePlot(neutrophil, features = gmp, label = T)
# Subset the Seurat object
neutrophil_2 <- subset(neutrophil_2, idents = c(0, 1, 2, 3, 4, 5, 6, 7))
# Assign new cluster IDs
new.cluster.ids <- c("BM_N4", "BM_N3", "BM_N3", "BM_N2", "BM_N4", "BM_N2", "BM_N1", "BM_N1")
names(new.cluster.ids) <- levels(neutrophil_2)
neutrophil_2 <- RenameIdents(neutrophil_2, new.cluster.ids)
levels(neutrophil_2) <- c("BM_N1", "BM_N2", "BM_N3", "BM_N4")
# Add new identity as SeuratObject level
neutrophil_2$Identity <- Idents(neutrophil_2)
head(neutrophil_2)
# Number of neutrophils by sample group
table(Idents(neutrophil_2), neutrophil_2$orig.ident)

## Visualize final neutrophil clusters split by cluster identity in a dimensional reduction plot
DimPlot(neutrophil_2, group.by = "Identity", split.by = "orig.ident", label = F, pt.size=1)

## Visualize final neutrophil clusters split by sample group in a dimensional reduction plot
# Split the Seurat object by sample group
neutrophil_2_naive <- subset(neutrophil_2, subset = orig.ident == "0 dpi")
neutrophil_2_14dpi <- subset(neutrophil_2, subset = orig.ident == "14 dpi")
# Assign new cluster IDs
# 1) Naive
new.cluster.ids <- c("Naive", "Naive", "Naive", "Naive")
names(new.cluster.ids) <- levels(neutrophil_2_naive)
neutrophil_2_naive <- RenameIdents(neutrophil_2_naive, new.cluster.ids)
neutrophil_2_naive$Identity <- Idents(neutrophil_2_naive)
# 2) 14 dpi
new.cluster.ids <- c("14 dpi", "14 dpi", "14 dpi", "14 dpi")
names(new.cluster.ids) <- levels(neutrophil_2_14dpi)
neutrophil_2_14dpi <- RenameIdents(neutrophil_2_14dpi, new.cluster.ids)
neutrophil_2_14dpi$Identity <- Idents(neutrophil_2_14dpi)
# Dimensional reduction plot 
DimPlot(neutrophil_2_naive, pt.size=1, cols=c('Naive'="#008080"))
DimPlot(neutrophil_2_14dpi,  pt.size=1, cols=c('14 dpi'="#FF8000"))


### Dot plot showing the expression level of neutrophil markers in final neutrophil clusters
final_np <- c("Cebpe", "Cd177", "S100a8", "Ly6g", "Ly6c2", "Itgam",
              "Mmp8", "Mmp9", "Csf1r", "Cxcr4", "Ccl6", "Cxcr2", "Il1b")
DotPlot(neutrophil_2, features = final_np, 
        cols = c("blue", "red"), dot.scale = 8) + theme(axis.text = element_text(size = 15, face = "italic"),
                                                        legend.title = element_text(size = 15),
                                                        legend.text = element_text(size = 10),
                                                        axis.title.x = element_text(size = 15),
                                                        axis.title.y.left = element_text(size = 15))


### Perform differential expression genes (DEGs) test

# Find top 10 genes per cluster
neutrophil_2_markers <- FindAllMarkers(neutrophil_2, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- neutrophil_2_markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)

# Visualize the top 10 genes per cluster in a heatmap
DoMultiBarHeatmap(object = neutrophil_2, features = as.character(top10$gene),
                  group.by="Identity", additional.group.by="group", size=4) + 
  scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", 
                       high = rev(c('#b2182b','#ef8a26','#fddbc7')),
                       midpoint = 0, guide = "colourbar",
                       aesthetics = "fill") + 
  theme(axis.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 15),legend.text = element_text(size = 10))


### Number of neutrophil subpopulations by sample group
Idents(neutrophil_2) = "Identity"
neutrophil_2$celltype.orig <- paste(Idents(neutrophil_2), neutrophil_2$orig.ident, sep = "_")
Idents(neutrophil_2) <- "celltype.orig"
DimPlot(neutrophil_2, reduction = "umap", group.by = "celltype.orig")
table(neutrophil_2@meta.data$celltype.orig)