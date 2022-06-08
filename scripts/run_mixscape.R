#!/home/baker02/miniconda3/envs/scRNA_env/bin/Rscript

# execute the analysis pipeline for the ER data
# each step will write intermediate output files to disk

#! How many cores per task?
#SBATCH --cpus-per-task=1
#! How much memory do you need?
#SBATCH --mem=12G
#! How much wallclock time will be required?
#SBATCH --time=10-00:00:00
#! Specify your email address here otherwise you won't recieve emails!
#SBATCH --mail-user=alexander.baker@cruk.cam.ac.uk
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
#SBATCH --no-requeue
#! General partition
#SBATCH -p general

# Data Cleaning packages
library(biomaRt)
library(tidyverse)
library(data.table)
library(hrbrthemes)
library(EnsDb.Hsapiens.v79)

# Load packages.
library(Seurat)
library(scales)
library(ggplot2)
library(argparse)
library(scRNAseq)
library(pheatmap)
library(patchwork)
library(DropletUtils)
library(BiocSingular)
library(SingleCellExperiment)

# Utility Function
extract_degs = function(results, matrix_genes){
    degs = list()
    
    for (gene in names(results)){
        gene_df = as.data.frame(results[[gene]])
        gene_df = gene_df[rownames(gene_df) %in% matrix_genes, ]
        gene_df = gene_df[!is.na(gene_df$padj), ]
        gene_df = gene_df[gene_df$padj <= 0.05, ]
        degs[[gene]] = rownames(gene_df)
    }
    
}

# Setup custom theme for plotting.
ctrl_label = 'Control'
data_dir = '/scratchb/fmlab/baker02/data'
sce_dir = paste0(data_dir, '/SCE/')
results_dir = paste0(data_dir, '/DEGs/results/')
custom_theme = theme(plot.title = element_text(size = 16, hjust = 0.5), legend.key.size = unit(0.7, "cm"), legend.text = element_text(size = 14))

############################### Loading the Franeigh et al 2021 dataset ###############################
parser = ArgumentParser(description='Conducting Differential Expression Analysis on Sub Groups of Franiegh Data')
parser$add_argument('library', help='Experimental Condition of Franiegh Data', type='character')

arguments =  parser$parse_args()
library_used = arguments$library

if (library_used == 'Yuza') {
    sce = readRDS(paste0(sce_dir, 'Yuza.Rds'))
    res = readRDS(paste0(results_dir, 'Yuza.Rds'))
} else if (library_used == 'Edwards') {
    sce = readRDS(paste0(sce_dir, 'Edwards.Rds'))
    res = readRDS(paste0(results_dir, 'Edwards.Rds'))
} else if (library_used == 'Brunello'){
    sce = readRDS(paste0(sce_dir, 'Brunello.Rds'))
    res = readRDS(paste0(results_dir, 'Brunello.Rds'))
} else {
    sce = readRDS(paste0(sce_dir, 'Merged.Rds'))
    res = readRDS(paste0(results_dir, 'Merged.Rds'))
}

dedup_cells = !duplicated(colnames(sce))
sce = sce[, dedup_cells]

# prepping count data and create Seurat Object
sce = CreateSeuratObject(counts=assay(sce, 'counts'), meta.data=as.data.frame(colData(sce)[colnames(sce), ]), paste0(library_used, '_library'), min.cells=3, min.features=200)
DefaultAssay(object = sce) = 'RNA'

print ("Normalizing Data...")
# Identifying the Phase of Cell Cycle a given cell is in 
sce = NormalizeData(object = sce) %>% FindVariableFeatures() %>% ScaleData()
sce = CellCycleScoring(sce, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=TRUE)

# Run Principle Component Analysis (PCA) and Uniform Manifold Approximation and Projection (UMAP) to reduce the dimensionality of the data.
sce = RunPCA(object = sce)
sce = RunUMAP(object = sce, dims = 1:40)

print ("Calculate Perturbation Signatures...")

# Calculating Perturbation Signatures
sce = CalcPerturbSig(object = sce, assay = "RNA", slot = "data", gd.class ="gene", nt.cell.class = ctrl_label,
      reduction = "pca", ndims = 40, num.neighbors = 20, split.by = "random_group", new.assay.name = "PRTB")

# Prepare PRTB assay for dimensionality reduction: Normalize data, find variable features and center data.
DefaultAssay(object = sce) = 'PRTB'
VariableFeatures(object = sce) = VariableFeatures(object = sce[["RNA"]])
sce = ScaleData(object = sce, do.scale = F, do.center = T)

# Run PCA and UMAP to reduce the dimensionality of the data.
sce = RunPCA(object = sce, reduction.key = 'prtbpca', reduction.name = 'prtbpca')
sce = RunUMAP(object = sce, dims = 1:40, reduction = 'prtbpca', reduction.key = 'prtbumap', reduction.name = 'prtbumap')

print ("Run MixScape...")

# Extract Differentially Expressed Genes and Run MixScape on scRNA-Seq CRISPR Screen Data
degs = extract_degs(res, rownames(sce))
sce = RunMixscape(object = sce, assay = "PRTB", slot = "scale.data", labels = "gene", nt.class.name = ctrl_label,
      de.gene.list = degs, provided.de.list = TRUE, min.de.genes = 5, iter.num = 10, de.assay = "RNA", verbose = F)

if (library_used == 'Yuza') {
    saveRDS(sce, paste0(sce_dir, 'Yuza_mixscape.Rds'))
} else if (library_used == 'Edwards') {
    saveRDS(sce, paste0(sce_dir, 'Edwards_mixscape.Rds'))
} else if (library_used == 'Brunello'){
     saveRDS(sce, paste0(sce_dir, 'Brunello_mixscape.Rds'))
} else {
    saveRDS(sce, paste0(sce_dir, 'Merged_mixscape.Rds'))
}