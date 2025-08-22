# Initial
rm(list = ls())

# Set directory
getwd()
setwd('/Users/User/Desktop/OC/')

# Library
library(rjson)
library(tidyverse)
library(TCGAbiolinks)
library(maftools)
library(writexl)
library(magrittr)
library(readxl)
library(SummarizedExperiment)
library(DESeq2)

# Find projects in GDC
proj <- getGDCprojects()
# View(proj)

c_metadata <- fromJSON(file = 'GDC_OC_2025_06_01/gc_cart_metadata.json')
c_sample_sheet <- read_tsv(file = 'GDC_OC_2025_06_01/gc_cart_sample_sheet.tsv')

c_metadata_df <- data.frame(barcode = sapply(c_metadata, function(x) x$associated_entities[[1]]$entity_submitter_id),
                            file_id = sapply(c_metadata, function(x) x$file_id),
                            file_name = sapply(c_metadata, function(x) x$file_name),
                            stringsAsFactors = FALSE)

c_sample_sheet_df <- c_sample_sheet[,c('File ID', 'File Name', 'Tissue Type', 'Tumor Descriptor', 'Project ID')] %>% 
  set_colnames(c('file_id', 'file_name', 'sample_type', 'tumor_descriptor', 'project'))

merge_df <- merge(c_metadata_df, c_sample_sheet_df, by = c('file_id', 'file_name'))

# project <- c_sample_sheet$`Project ID` %>% table() %>% names()
# barcode <- c_metadata_df$barcode
query <- GDCquery(project = project[2],
                  data.category = 'Transcriptome Profiling',
                  data.type = 'Gene Expression Quantification',
                  experimental.strategy = 'RNA-Seq',
                  workflow.type = 'STAR - Counts',
                  barcode = c_metadata_df$barcode)

# GDCdownload(query = query, directory = 'GDC/', files.per.chunk = 10)

oc <- GDCprepare(query = query, directory = 'GDC/', summarizedExperiment = T)

oc_meta <- colData(oc) %>% as.data.frame()
oc_meta_sel <- oc_meta[, c("barcode", "tissue_type", "tissue_or_organ_of_origin")]
oc_meta_sel <- oc_meta_sel[oc_meta_sel$tissue_or_organ_of_origin %in% c('Upper Gum', 'Tongue, NOS',
                                                                        'Overlapping lesion of lip, oral cavity and pharynx',
                                                                        'Mouth, NOS', 'Lymph nodes of head, face and neck',
                                                                        'Hard palate', 'Gum, NOS', 'Floor of mouth, NOS',
                                                                        'Base of tongue, NOS'), ]

oc_meta_sel <- oc_meta_sel[, c("barcode", "tissue_type")]
rownames(oc_meta_sel) <- NULL
oc_meta_sel <- column_to_rownames(oc_meta_sel, var = "barcode")
oc_meta_sel$tissue_type <- as.factor(oc_meta_sel$tissue_type)

genes <- rowData(oc) %>% as.data.frame() %>% select(c('gene_id', 'gene_name', 'gene_type')) %>% set_rownames(NULL)
genes <- genes[genes$gene_type %in% c('protein_coding'), ] %>% select(c('gene_id', 'gene_name'))

oc_count <- assay(oc)
oc_count_sel <- oc_count[, rownames(oc_meta_sel)]
oc_count_sel <- oc_count_sel[rownames(oc_count_sel) %in% genes$gene_id,]
rownames(oc_count_sel) <- genes$gene_name[match(rownames(oc_count_sel), genes$gene_id)]
oc_count_sel_mean <- apply(oc_count_sel, 2, tapply, rownames(oc_count_sel), mean) %>% round()

all(colnames(oc_count_sel) == rownames(oc_meta_sel))

ddsobj <- DESeqDataSetFromMatrix(countData = oc_count_sel_mean,
                                 colData = oc_meta_sel,
                                 design = ~ tissue_type)

filt <- min(matrixStats::count(oc_meta_sel$tissue_type == 'Tumor'),
            matrixStats::count(oc_meta_sel$tissue_type == 'Normal'))
keep <- rowSums(counts(ddsobj) >= 10) >= filt
ddsobj.filt <- ddsobj[keep,]

# DESeq2
dds <- DESeq(object = ddsobj.filt, fitType = 'local')

# Results and filter DEGs
res <- results(object = dds, contrast = c('tissue_type', 'Tumor', 'Normal'))
res_gene <- res %>% data.frame()
deg_gene <- res_gene %>% filter(!is.na(x = padj), padj <= 0.05, abs(log2FoldChange) >= 1)