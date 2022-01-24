# Seurat_mt_batch_cc.R
# Last updated on NOV.23 2020
# used for scReQTL (SRR10018149, SRR10018150,SRR10018151 dataset)

# load package
suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(tidyverse, quietly = TRUE))
suppressMessages(library(Seurat, quietly = TRUE))
suppressMessages(library(SingleR, quietly = TRUE))



if (packageVersion("Seurat") < "3.0.0") {
  stop(paste0("You have Seurat version", packageVersion("Seurat"), "installed. Please make sure you have Seurat version > 3.0.0"))
}

# Load all data-----------------------------------------------------------------------------------------
#For these samples: 3k<nfeatures<8k and mt_percent<7
# load filtered Seurat data
SRR10018149 <- readRDS('SRR10018149_Seurat_clustered_singleR.rds')
SRR10018150 <- readRDS('SRR10018150_Seurat_clustered_singleR.rds')
SRR10018151 <- readRDS('SRR10018151_Seurat_clustered_singleR.rds')
SRR10018152 <- readRDS('SRR10018152_Seurat_clustered_singleR.rds')

# cell cycle calculation
s.genes <- fread('stage_s.txt')
s.genes <- s.genes$Gene
g2m.genes <- fread('stage_G2M.txt')
g2m.genes <- g2m.genes$Gene
#-------------------------------------------------------------------------------------------------------

# Data pre-processing-----------------------------------------------------------------------------------
# remove mitochondrial genes
SRR10018149_outMT <- SRR10018149@assays[['RNA']]@counts
SRR10018149_outMT <- SRR10018149_outMT[-grep(pattern = '^MT', row.names(SRR10018149_outMT)),  ]
SRR10018149 <- CreateSeuratObject(counts = SRR10018149_outMT,project = "SRR10018149")
rm(SRR10018149_outMT)
SRR10018150_outMT <- SRR10018150@assays[['RNA']]@counts
SRR10018150_outMT <- SRR10018150_outMT[-grep(pattern = '^MT', row.names(SRR10018150_outMT)),  ]
SRR10018150 <- CreateSeuratObject(counts = SRR10018150_outMT,project = "SRR10018150")
rm(SRR10018150_outMT)
SRR10018151_outMT <- SRR10018151@assays[['RNA']]@counts
SRR10018151_outMT <- SRR10018151_outMT[-grep(pattern = '^MT', row.names(SRR10018151_outMT)),  ]
SRR10018151 <- CreateSeuratObject(counts = SRR10018151_outMT,project = "SRR10018151")
rm(SRR10018151_outMT)
SRR10018152_outMT <- SRR10018152@assays[['RNA']]@counts
SRR10018152_outMT <- SRR10018152_outMT[-grep(pattern = '^MT', row.names(SRR10018152_outMT)),  ]
SRR10018152 <- CreateSeuratObject(counts = SRR10018152_outMT,project = "SRR10018152")
rm(SRR10018152_outMT)
# normalization and scale the dataset (only use nUMI to regress which is default)
# Set 3000 variable features(default: variable.features.n = 3000)
## our data has been regressed (nUMI) based on umi_tools
SRR10018149 <- SCTransform(SRR10018149)
SRR10018150 <- SCTransform(SRR10018150)
SRR10018151 <- SCTransform(SRR10018151)
SRR10018152 <- SCTransform(SRR10018152)
#-------------------------------------------------------------------------------------------------------

#integrate three samples together
list <- list(SRR10018152, SRR10018151, SRR10018150, SRR10018149)
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 2000)
options(future.globals.maxSize = 9768*1024^2)
list <- PrepSCTIntegration(object.list = list, anchor.features = features,
                           verbose = FALSE)
anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT",
                                  anchor.features = features, verbose = FALSE)
samples_comb <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)

png("samples_comb_feature_distribution_vlnplot.png", width = 850, height = 400)
VlnPlot(object = samples_comb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# saveRDS(samples_comb, 'samples_comb_outMT_perPCA.rds')

# cluster
DefaultAssay(samples_comb) <- "integrated"
samples_comb <- RunPCA(samples_comb, verbose = FALSE)
samples_comb <- RunUMAP(samples_comb, verbose = FALSE, dim = 1:30)
samples_comb<- FindNeighbors(samples_comb, dims = 1:30) #,  k.param = 10
samples_comb<- FindClusters(samples_comb, resolution = 0.2)

# plot cluster information
png("samples_comb_clusters_Seurat_umap.png", width = 450, height = 400)
DimPlot(samples_comb, group.by = "orig.ident",reduction = "umap", label = TRUE)
dev.off()

# annotation cell type
rna_re <- BlueprintEncodeData()
cluster <- samples_comb@active.ident
b <- GetAssayData(samples_comb)
result_cluster <- SingleR(test = b, ref = rna_re, labels = rna_re$label.fine, method="cluster", clusters = cluster)
samples_comb[["SingleR.cluster.labels"]] <-
  result_cluster$labels[match(samples_comb[[]]["seurat_clusters"]$seurat_clusters, rownames(result_cluster))]

png("samples_comb_clusters_SingleR_umap.png", width = 450, height = 400)
DimPlot(samples_comb, group.by =  "SingleR.cluster.labels", reduction = "umap", label = TRUE)
dev.off()

png("samples_comb_clusters_SingleR_sample_split_umap.png", width = 850, height = 400)
DimPlot(samples_comb, group.by= "SingleR.cluster.labels", split.by = "orig.ident",reduction = "umap", label = TRUE) + NoLegend()
dev.off()

#split the file
list_new <- SplitObject(samples_comb, split.by = "orig.ident")
SRR10018149 <- list_new$SRR10018149
SRR10018150 <- list_new$SRR10018150
SRR10018151 <- list_new$SRR10018151
SRR10018152 <- list_new$SRR10018152

# Identifying clusters individually after batch effect removal
# SRR10018149-----------------------------------------------------------------------------------------------
# Assign cell cycle scores and visualize 
DefaultAssay(SRR10018149) <- "integrated"
SRR10018149 <- CellCycleScoring(SRR10018149, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) #  (use head(SRR10018149[[]]) to check)
SRR10018149 <- RunPCA(SRR10018149)
SRR10018149 <- RunUMAP(SRR10018149, dims = 1:30)
png("SRR10018149_cell_cycle_effect_umap.png", width = 450, height = 400)
DimPlot(SRR10018149, group.by = "Phase", reduction = "umap",label = TRUE)
dev.off()

# Regress out the cell cycle effect
SRR10018149 <- ScaleData(SRR10018149, vars.to.regress = c("S.Score", "G2M.Score"),
                 use.umi = FALSE, do.scale = FALSE, do.center =  FALSE )
SRR10018149 <- RunPCA(SRR10018149)
SRR10018149<- RunUMAP(SRR10018149, dims = 1:30)#, spread = 1, min.dist = 0.0001, n.epochs = 200, n.neighbors = 10)

png("SRR10018149_cell_cycle_effect_removed_umap.png", width = 450, height = 400)
DimPlot(SRR10018149, group.by = "Phase", reduction = "umap",label = TRUE)
dev.off()

# Cluster the cells
SRR10018149<- FindNeighbors(SRR10018149, dims = 1:30)#, k.param = 10)
SRR10018149<- FindClusters(SRR10018149, resolution = 0.2)

png("SRR10018149_elbow_plot.png", width = 450, height = 400)
ElbowPlot(SRR10018149)
dev.off()

# annotation
b <- GetAssayData(SRR10018149)
res_re <- BlueprintEncodeData()
cluster <- SRR10018149@active.ident
result_cluster <- SingleR(test = b, ref = res_re, labels = res_re$label.fine, method="cluster", clusters = cluster)

png("SRR10018149_SingleR_cell_label_scores.png", width = 450, height = 400)
plotScoreHeatmap(result_cluster)
dev.off()

SRR10018149[["SingleR.cluster.labels"]] <-
  result_cluster$labels[match(SRR10018149[[]]["seurat_clusters"]$seurat_clusters, rownames(result_cluster))]

png("SRR10018149_clusters_SingleR_umap.png", width = 450, height = 400)
DimPlot(SRR10018149, group.by =  "SingleR.cluster.labels", reduction = "umap", label = TRUE)
dev.off()

png("SRR10018149_single_cell_feature_exp_heatmap.png", width = 450, height = 400)
DoHeatmap(subset(SRR10018149, downsample = 100), features = features, size = 3, group.by = "SingleR.cluster.labels")
dev.off()

# write cluster information
ide <- as.data.frame(SRR10018149@active.ident)
fwrite(ide, "SRR10018149_cluster_sample.txt", row.names = T, sep = '\t', quote = F)
cluster <- data.frame(SRR10018149[["SingleR.cluster.labels"]])
fwrite(cluster, "SRR10018149_cluster_sample_named.txt", row.names = T, sep = ',', quote = F)

# pca used for SRR10018149
SRR10018149_pca <- t(as.data.frame(SRR10018149@reductions[["pca"]]@cell.embeddings))
fwrite(SRR10018149_pca, "SRR10018149_pca_matrix.txt", row.names = T, sep = '\t', quote = F)
#-------------------------------------------------------------------------------------------------------

#SRR10018150----------------------------------------------------------------------------------------------------
DefaultAssay(SRR10018150) <- "integrated"
SRR10018150 <- CellCycleScoring(SRR10018150, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) #  (use head(SRR10018150[[]]) to check)
SRR10018150 <- RunPCA(SRR10018150)
SRR10018150 <- RunUMAP(SRR10018150, dims = 1:30)

png("SRR10018150_cell_cycle_effect_umap.png", width = 450, height = 400)
DimPlot(object = SRR10018150, group.by = "Phase", reduction = "umap",label = TRUE)
dev.off()

# Sale SRR10018150 sample by the score of S and G2M
SRR10018150 <- ScaleData(SRR10018150, vars.to.regress = c("S.Score", "G2M.Score"),
                use.umi = FALSE, do.scale = FALSE, do.center =  FALSE )
SRR10018150 <- RunPCA(SRR10018150)
SRR10018150<- RunUMAP(SRR10018150, dims = 1:30)#, spread = 1, min.dist = 0.0001, n.epochs = 200, n.neighbors = 10)

png("SRR10018150_cell_cycle_effect_removed_umap.png", width = 450, height = 400)
DimPlot(object = SRR10018150, group.by = "Phase", reduction = "umap",label = TRUE)
dev.off()

SRR10018150<- FindNeighbors(SRR10018150, dims = 1:30)#, k.param = 10)
SRR10018150<- FindClusters(SRR10018150, resolution = 0.2)

png("SRR10018150_elbow_plot.png", width = 450, height = 400)
ElbowPlot(SRR10018150)
dev.off()

# write the cluster information
ide <- as.data.frame(SRR10018150@active.ident)
fwrite(ide, "SRR10018150_cluster_sample.txt", row.names = T, sep = '\t', quote = F)
cluster <- data.frame(SRR10018150[["SingleR.cluster.labels"]])
fwrite(cluster, "SRR10018150_cluster_sample_named.txt", row.names = T, sep = ',', quote = F)

# annotation
b <- GetAssayData(SRR10018150)
cluster <- SRR10018150@active.ident
result_cluster <- SingleR(test = b, ref = res_re, labels = res_re$label.fine, method="cluster", clusters = cluster)
SRR10018150[["SingleR.cluster.labels"]] <-
  result_cluster$labels[match(SRR10018150[[]]["seurat_clusters"]$seurat_clusters, rownames(result_cluster))]

png("SRR10018150_SingleR_cell_label_scores.png", width = 450, height = 400)
plotScoreHeatmap(result_cluster)
dev.off()

png("SRR10018150_SingleR_cell_label_scores.png", width = 450, height = 400)
DimPlot(SRR10018150, group.by =  "SingleR.cluster.labels", reduction = "umap", label = TRUE)
dev.off()

png("SRR10018150_SRR10018149_single_cell_feature_exp_heatmap.png", width = 450, height = 400)
DoHeatmap(subset(SRR10018150, downsample = 100), features = features, size = 3, group.by = "SingleR.cluster.labels")
dev.off()

# pca used for SRR10018150
SRR10018150_pca <- t(as.data.frame(SRR10018150@reductions[["pca"]]@cell.embeddings))
fwrite(SRR10018150_pca, "SRR10018150_pca_matrix.txt", row.names = T, sep = '\t', quote = F)
#-------------------------------------------------------------------------------------------------------

#SRR10018151-----------------------------------------------------------------------------------------------------
DefaultAssay(SRR10018151) <- "integrated"
SRR10018151 <- CellCycleScoring(SRR10018151, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) # #  (use head(SRR10018151[[]]) to check)
SRR10018151 <- RunPCA(SRR10018151)
SRR10018151 <- RunUMAP(SRR10018151, dims = 1:30)

png("SRR10018151_cell_cycle_effect_umap.png", width = 450, height = 400)
DimPlot(object = SRR10018151, group.by = "Phase", reduction = "umap",label = TRUE)
dev.off()

# Scale SRR10018151 sample by the score of S and G2M
SRR10018151 <- ScaleData(SRR10018151, vars.to.regress = c("S.Score", "G2M.Score"),
                use.umi = FALSE, do.scale = FALSE, do.center =  FALSE )
SRR10018151 <- RunPCA(SRR10018151)
SRR10018151<- RunUMAP(SRR10018151, dims = 1:30)#, spread = 1, min.dist = 0.0001, n.epochs = 200, n.neighbors = 10)

png("SRR10018151_cell_cycle_effect_removed_umap.png", width = 450, height = 400)
DimPlot(object = SRR10018151, group.by = "Phase", reduction = "umap",label = TRUE)
dev.off()

SRR10018151<- FindNeighbors(SRR10018151, dims = 1:30)#, k.param = 10)
SRR10018151<- FindClusters(SRR10018151, resolution = 0.2)

png("SRR10018151_elbow_plot.png", width = 450, height = 400)
ElbowPlot(SRR10018151)
dev.off()

# write cluster information
ide <- as.data.frame(SRR10018151@active.ident)
fwrite(ide, "SRR10018151_cluster_sample.txt", row.names = T, sep = '\t', quote = F)
cluster <- data.frame(SRR10018151[["SingleR.cluster.labels"]])
fwrite(cluster, "SRR10018151_cluster_sample_named.txt", row.names = T, sep = ',', quote = F)

# annotation
b <- GetAssayData(SRR10018151)
cluster <- SRR10018151@active.ident
result_cluster <- SingleR(test = b, ref = res_re, labels = res_re$label.fine, method="cluster", clusters = cluster)
SRR10018151[["SingleR.cluster.labels"]] <-
  result_cluster$labels[match(SRR10018151[[]]["seurat_clusters"]$seurat_clusters, rownames(result_cluster))]

png("SRR10018151_SingleR_cell_label_scores.png", width = 450, height = 400)
plotScoreHeatmap(result_cluster)
dev.off()

png("SRR10018151_clusters_SingleR_umap.png", width = 450, height = 400)
DimPlot(SRR10018151, group.by =  "SingleR.cluster.labels", reduction = "umap", label = TRUE)
dev.off()

png("SRR10018151_single_cell_feature_exp_heatmap.png", width = 450, height = 400)
DoHeatmap(subset(SRR10018151, downsample = 100), features = features, size = 3, group.by = "SingleR.cluster.labels")
dev.off()

# pca used for SRR10018151
SRR10018151_pca <- t(as.data.frame(SRR10018151@reductions[["pca"]]@cell.embeddings))
fwrite(SRR10018151_pca, "SRR10018151_pca_matrix.txt", row.names = T, sep = '\t', quote = F)
#-------------------------------------------------------------------------------------------------------

#SRR10018152----------------------------------------------------------------------------------------------------
DefaultAssay(SRR10018152) <- "integrated"
SRR10018152 <- CellCycleScoring(SRR10018152, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) #  (use head(SRR10018152[[]]) to check)
SRR10018152 <- RunPCA(SRR10018152)
SRR10018152 <- RunUMAP(SRR10018152, dims = 1:30)

png("SRR10018152_cell_cycle_effect_umap.png", width = 450, height = 400)
DimPlot(object = SRR10018152, group.by = "Phase", reduction = "umap",label = TRUE)
dev.off()

# Sale SRR10018152 sample by the score of S and G2M
SRR10018152 <- ScaleData(SRR10018152, vars.to.regress = c("S.Score", "G2M.Score"),
                     use.umi = FALSE, do.scale = FALSE, do.center =  FALSE )
SRR10018152 <- RunPCA(SRR10018152)
SRR10018152<- RunUMAP(SRR10018152, dims = 1:30)#, spread = 1, min.dist = 0.0001, n.epochs = 200, n.neighbors = 10)

png("SRR10018152_cell_cycle_effect_removed_umap.png", width = 450, height = 400)
DimPlot(object = SRR10018152, group.by = "Phase", reduction = "umap",label = TRUE)
dev.off()

SRR10018152<- FindNeighbors(SRR10018152, dims = 1:30)#, k.param = 10)
SRR10018152<- FindClusters(SRR10018152, resolution = 0.2)

png("SRR10018152_elbow_plot.png", width = 450, height = 400)
ElbowPlot(SRR10018152)
dev.off()

# write the cluster information
ide <- as.data.frame(SRR10018152@active.ident)
fwrite(ide, "SRR10018152_cluster_sample.txt", row.names = T, sep = '\t', quote = F)
cluster <- data.frame(SRR10018152[["SingleR.cluster.labels"]])
fwrite(cluster, "SRR10018152_cluster_sample_named.txt", row.names = T, sep = ',', quote = F)

# annotation
b <- GetAssayData(SRR10018152)
cluster <- SRR10018152@active.ident
result_cluster <- SingleR(test = b, ref = res_re, labels = res_re$label.fine, method="cluster", clusters = cluster)
SRR10018152[["SingleR.cluster.labels"]] <-
  result_cluster$labels[match(SRR10018152[[]]["seurat_clusters"]$seurat_clusters, rownames(result_cluster))]

png("SRR10018152_SingleR_cell_label_scores.png", width = 450, height = 400)
plotScoreHeatmap(result_cluster)
dev.off()

png("SRR10018152_SingleR_cell_label_scores.png", width = 450, height = 400)
DimPlot(SRR10018152, group.by =  "SingleR.cluster.labels", reduction = "umap", label = TRUE)
dev.off()

png("SRR10018152_SRR10018149_single_cell_feature_exp_heatmap.png", width = 450, height = 400)
DoHeatmap(subset(SRR10018152, downsample = 100), features = features, size = 3, group.by = "SingleR.cluster.labels")
dev.off()

# pca used for SRR10018152
SRR10018152_pca <- t(as.data.frame(SRR10018152@reductions[["pca"]]@cell.embeddings))
fwrite(SRR10018152_pca, "SRR10018152_pca_matrix.txt", row.names = T, sep = '\t', quote = F)
#-------------------------------------------------------------------------------------------------------

# save output
saveRDS(SRR10018149, "SRR10018149_ccregress.rds")
saveRDS(SRR10018150, "SRR10018150_ccregress.rds")
saveRDS(SRR10018151, "SRR10018151_ccregress.rds")
saveRDS(SRR10018152, "SRR10018152_ccregress.rds")
