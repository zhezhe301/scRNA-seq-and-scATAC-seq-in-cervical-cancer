# snode1
# module load R/3.5.2
# R
## source /homes/zhezhang/miniconda2/etc/profile.d/conda.sh
## conda activate /liulab/xhuang/miniconda3/envs/MAESTRO0904
## R
## library(MAESTRO,lib.loc='/liulab/xhuang/miniconda3/envs/MAESTRO0904/lib/R/library')
library(Seurat,lib.loc='/homes/zhezhang/software/R/Seurat3')
library(Matrix)
library(ggplot2)
library(MAESTRO,lib.loc='/liulab/xhuang/miniconda3/envs/MAESTRO0904/lib/R/library')

# scRNA pipeline for individual sample
scRNA_qc <- function(pbmc.data = pbmc.data, project=project ){
  pbmc <- CreateSeuratObject(counts = pbmc.data, project = project, min.cells = 3, min.features = 200)
  pbmc
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  # pdf(paste0(project,'_vlnplot.pdf'),width=8,height=4.5)
  p1 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0(project,'_vlnplot.pdf'),p1,width=8,height=4.5)
  # pdf(paste0(project,'feature.pdf'),width=10,height=4.5)
  plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p2 <- CombinePlots(plots = list(plot1, plot2))
  ggsave(paste0(project,'feature.pdf'),p2,width=12,height=4.5)
}
setwd('/liulab/zhezhang/project/CESC/scRNA/scRNA_Analyses/CellRanger_out/intergation')
library(cowplot)
p1ca.rna <- Read10X('/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/TCGA-CESC/cesc/scRNA/scRNA_Analyses/CellRanger_out/190608A_1_ai/outs/filtered_feature_bc_matrix')
p1ln.rna <- Read10X('/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/TCGA-CESC/cesc/scRNA/scRNA_Analyses/CellRanger_out/190608A_2_linba/outs/filtered_feature_bc_matrix')
p2ca.rna <- Read10X('/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/TCGA-CESC/cesc/scRNA/scRNA_Analyses/CellRanger_out/190608B_A/outs/filtered_feature_bc_matrix')
p2ln.rna <- Read10X('/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/TCGA-CESC/cesc/scRNA/scRNA_Analyses/CellRanger_out/190608B_L/outs/filtered_feature_bc_matrix')
p3ca.rna <- Read10X('/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/TCGA-CESC/cesc/scRNA/scRNA_Analyses/CellRanger_out/190608C_A/outs/filtered_feature_bc_matrix')
p3ln.rna <- Read10X('/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/TCGA-CESC/cesc/scRNA/scRNA_Analyses/CellRanger_out/190608C_L/outs/filtered_feature_bc_matrix')
scRNA_qc(p1ca.rna, project = 'P1_CA_RNA')
scRNA_qc(p1ln.rna, project = 'P2_LN_RNA')
scRNA_qc(p1ca.rna, project = 'P2_CA_RNA')
scRNA_qc(p1ln.rna, project = 'P2_LN_RNA')
scRNA_qc(p1ca.rna, project = 'P3_CA_RNA')
scRNA_qc(p1ln.rna, project = 'P3_LN_RNA')
# to observe groups split by patients
colnames(p1ca.rna) <- paste0("p1_",colnames(p1ca.rna))
colnames(p1ln.rna) <- paste0("p1_",colnames(p1ln.rna))
colnames(p2ca.rna) <- paste0("p2_",colnames(p2ca.rna))
colnames(p2ln.rna) <- paste0("p2_",colnames(p2ln.rna))
colnames(p3ca.rna) <- paste0("p3_",colnames(p3ca.rna))
colnames(p3ln.rna) <- paste0("p3_",colnames(p3ln.rna))
gene <- cbind(p1ca.rna,p1ln.rna,p2ca.rna,p2ln.rna,p3ca.rna,p3ln.rna)
# 'gene' matrix: 33694 * 46405
gene <- gene[-grep('^MT-', rownames(gene)),]
# 'gene' matrix: 33681 * 46405
scRNA_qc(gene, project = 'CESC')
pbmc <- CreateSeuratObject(counts = gene, project = "CESC")
pbmc@meta.data$ca_ln <- as.factor(c(rep('ca',6000),rep('ln',10098),rep('ca',5371),rep('ln',9349),rep('ca',4657),rep('ln',10930)))
#pbmc <- CreateSeuratObject(counts = gene, project = "CESC", min.cells = 3, min.features = 200)
# this 'pbmc' matrix: 24618 * 44604
# pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# cancer <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10)
pbmc<- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 11000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 10000)
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
p3 <- CombinePlots(plots = list(plot1, plot2))
ggsave(paste0('integ_withoutcca','feature.pdf'),p3,width=12,height=4.5)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
p4 <- DimPlot(pbmc, reduction = "pca")
ggsave(paste0('integ_withoutcca','_pca.pdf'),p4,width=6,height=4.5)
p5 <- DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
ggsave(paste0('integ_withoutcca','_heaatmap.pdf'),p5,width=6,height=4.5)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.1)
pbmc <- RunUMAP(pbmc, dims = 1:30)
p6 <- DimPlot(pbmc, reduction = "umap")
ggsave(paste0('integ_withoutcca','_umap.pdf'),p6,width=6,height=4.5)
p7 <- DimPlot(pbmc, reduction = "umap",split.by = "orig.ident")
ggsave(paste0('integ_withoutcca','_umap_bypatient.pdf'),p7,width=12,height=4.5)
p8 <- DimPlot(pbmc, reduction = "umap",split.by = "ca_ln")
ggsave(paste0('integ_withoutcca','_umap_byloc.pdf'),p8,width=9,height=4.5)
p9 <- DimPlot(pbmc, reduction = "umap", group.by = "orig.ident")
ggsave(paste0('integ_withoutcca','_umap_mix.pdf'),p9,width=6,height=4.5)
p10 <- DimPlot(pbmc, reduction = "umap", label = TRUE)
ggsave(paste0('integ_withoutcca','_umap_noannotation.pdf'),p10,width=6,height=4.5)





cluster.markers1 <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.1)
cluster.markers2 <- FindAllMarkers(object = pbmc, min.pct = 0.1)

saveRDS(pbmc,'seuratobj.rds')
saveRDS(cluster.markers1,'cluster.markers1.rds')
saveRDS(cluster.markers2,'cluster.markers2.rds')

AssignCellTypeSeurat <- function(cluster.markers, marker.use = markers.CIBERSORT){
  print(unique(cluster.markers$cluster))
  print(as.integer(levels(cluster.markers$cluster)))
  print(as.character(levels(cluster.markers$cluster)))
  CT = sapply(as.character(levels(cluster.markers$cluster)), function(x){
    idx = cluster.markers$cluster==x
    tmpFC = cluster.markers$avg_logFC[idx]
    names(tmpFC) = toupper(cluster.markers$gene[idx])
    
    score_cluster <- NULL
    for (i in 1:length(summary(marker.use[,1]))) {
      g.n <- which(marker.use[,1]==names(summary(marker.use[,1]))[i])
      y <- marker.use[g.n, 2]
      score = sum(tmpFC[y], na.rm = TRUE) / log2(length(y))
      names(score) <- names(summary(marker.use[,1]))[i]
      score_cluster<- c(score_cluster, score)
    }
    # if(max(score_cluster, na.rm = TRUE)>0.1)
    CT = names(score_cluster)[which.max(score_cluster)]
    # else
      # CT = "Others"
  })
  return(CT)
}

data(human.immune.CIBERSORT)
human.immune.CIBERSORT <- rbind(human.immune.CIBERSORT,
                                cbind("Fibroblasts", c("FAP", "PDPN", "MMP2", "PDGFRA", "THY1", "MMP11","PDGFRL", "TGFB3", "COL1A2", "DCN", "COL3A1", "COL6A1")),
                                cbind("Myofibroblasts", c("ACTA2", "MCAM", "MYLK", "MYL9", "IL6", "PDGFA")),
                                cbind("Endothelial", c("PECAM1", "VWF", "ENG")),
                                cbind("Pericyte",c("ACTA1", "CD248","ACTA2","MMP9")),
                                cbind("SCC",c("CDKN2A","PIK3CA","MCM3","MCM5","CDC6","GMNN","NUF2","CCND1","TOP2A","BIRC5","KRT19","TP53")) )

current.cluster.ids = as.integer(levels(cluster.markers$cluster))
new.cluster.ids = AssignCellTypeSeurat(cluster.markers, marker.use = human.immune.CIBERSORT)


####################################
# scATAC 
setwd('/liulab/zhezhang/project/CESC/scATAC/scATAC_Analyses/CellRanger_out/intergation')
p1ca.atac.peak <- Read10X_h5('/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/TCGA-CESC/cesc/scRNA/scATAC_Analyses/CellRanger_out/1-ai/outs/filtered_peak_bc_matrix.h5')
p1ca.atac.tf <- Read10X_h5('/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/TCGA-CESC/cesc/scRNA/scATAC_Analyses/CellRanger_out/1-ai/outs/filtered_tf_bc_matrix.h5')
p1ln.atac.peak <- Read10X_h5('/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/TCGA-CESC/cesc/scRNA/scATAC_Analyses/CellRanger_out/2-linba/outs/filtered_peak_bc_matrix.h5')
p1ln.atac.tf <- Read10X_h5('/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/TCGA-CESC/cesc/scRNA/scATAC_Analyses/CellRanger_out/2-linba/outs/filtered_tf_bc_matrix.h5')
p2ca.atac.peak <- Read10X_h5('/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/TCGA-CESC/cesc/scRNA/scATAC_Analyses/CellRanger_out/190608B-A/outs/filtered_peak_bc_matrix.h5')
p2ca.atac.tf <- Read10X_h5('/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/TCGA-CESC/cesc/scRNA/scATAC_Analyses/CellRanger_out/190608B-A/outs/filtered_tf_bc_matrix.h5')
p2ln.atac.peak <- Read10X_h5('/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/TCGA-CESC/cesc/scRNA/scATAC_Analyses/CellRanger_out/190608B-L/outs/filtered_peak_bc_matrix.h5')
p2ln.atac.tf <- Read10X_h5('/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/TCGA-CESC/cesc/scRNA/scATAC_Analyses/CellRanger_out/190608B-L/outs/filtered_tf_bc_matrix.h5')
p3ca.atac.peak <- Read10X_h5('/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/TCGA-CESC/cesc/scRNA/scATAC_Analyses/CellRanger_out/190608C-A/outs/filtered_peak_bc_matrix.h5')
p3ca.atac.tf <- Read10X_h5('/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/TCGA-CESC/cesc/scRNA/scATAC_Analyses/CellRanger_out/190608C-A/outs/filtered_tf_bc_matrix.h5')
p3ln.atac.peak <- Read10X_h5('/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/TCGA-CESC/cesc/scRNA/scATAC_Analyses/CellRanger_out/190608C-L/outs/filtered_peak_bc_matrix.h5')
p3ln.atac.tf <- Read10X_h5('/tcga/tcga_ATAC_Seq/tcga_ATAC_Seq_processed_data/TCGA-CESC/cesc/scRNA/scATAC_Analyses/CellRanger_out/190608C-L/outs/filtered_tf_bc_matrix.h5')




























