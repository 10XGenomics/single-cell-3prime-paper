#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#

rm(list=ls()) # clear workspace
# ----------------------------
# load relevant libraries
# ----------------------------
library(Matrix)
library(ggplot2)
library(Rtsne)
library(svd)
library(dplyr)
library(data.table)
library(pheatmap)
# -------------------------------------
# specify paths and load functions
# -------------------------------------
DATA_DIR <- "PATH_TO_DATA_DIRECTORY"        # SPECIFY HERE
PROG_DIR <- "PATH_TO_PROGRAM_DIRECTORY"     # SPECIFY HERE
FIG_DIR <-  "PATH_TO_FIGURE_DIRECTORY"      # SPECIFY HERE
source(file.path(PROG_DIR,'util.R')) 
source(file.path(PROG_DIR,'select_pure_pbmc.R'))
# -----------------------
# load purified data
# -----------------------
pure_pbmcs <- readRDS(file.path(DATA_DIR,'all_pure_pbmc_data.rds'))
all_data <- pure_pbmcs$all_data
all_json <- pure_pbmcs$all_json
all_metrics <- pure_pbmcs$all_metrics
all_mol_info <- pure_pbmcs$all_mol_info
genes <- all_data[[1]]$hg19$gene_symbols
# -------------------------------------------------------------------------
# downsample mapped reads/cell
# so that all samples have the same # of confidently mapped reads/cell
# -------------------------------------------------------------------------
set.seed(1)
rpc <- all_metrics %>% 
  mutate(conf_mapped_rpc=raw_rpc*conf_mapped_frac*good_bc_frac*good_umi_frac) %>%
  select(sample_id, description, conf_mapped_rpc)
tgt_rpc <- floor(min(rpc$conf_mapped_rpc)) # 13995
subsampled_purified_mats <- lapply(1:length(all_data), function(i) { # subsample the matrix to match tgt_rpc
  cat(sprintf("%d...\n", i))
  .downsample_gene_bc_mtx(all_json[[i]], all_data[[i]], all_mol_info[[i]], tgt_rpc, 'conf_mapped_reads')[[1]]
} )
# -------------------------------------------
# run PCA and tSNE on down-sampled data
# the steps take a long time to complete
# -------------------------------------------
# using 10 PCs provides enough information to filter out undesired populations from purified populations
set.seed(1)
all_pure_pca<-lapply(1:length(subsampled_purified_mats),function(i) pure_pca_i<-.do_propack(subsampled_purified_mats[[i]],10))
all_pure_tsne<-lapply(1:length(all_pure_pca),function(i) pure_tsne_i<-Rtsne(all_pure_pca[[i]]$pca,pca=F))
# ------------------------------------
# curate the purified populations
# ------------------------------------
pure_id<-c("CD34+","CD56+ NK","CD4+/CD45RA+/CD25- Naive T", "CD4+/CD25 T Reg","CD8+/CD45RA+ Naive Cytotoxic",
           "CD4+/CD45RO+ Memory","CD8+ Cytotoxic T","CD19+ B","CD4+ T Helper2","CD14+ Monocyte","Dendritic")
sub_idx <-list(data.frame(sample=1, use=(get_pure_pop_idx(genes,pure_id[1],all_pure_pca[[1]],all_pure_tsne[[1]],FIG_DIR))), 
               data.frame(sample=2, use=(get_pure_pop_idx(genes,pure_id[2],all_pure_pca[[2]],all_pure_tsne[[2]],FIG_DIR))),
               data.frame(sample=3, use=(get_pure_pop_idx(genes,pure_id[3],all_pure_pca[[3]],all_pure_tsne[[3]],FIG_DIR))),
               data.frame(sample=4, use=(get_pure_pop_idx(genes,pure_id[4],all_pure_pca[[4]],all_pure_tsne[[4]],FIG_DIR))),
               data.frame(sample=5, use=(get_pure_pop_idx(genes,pure_id[5],all_pure_pca[[5]],all_pure_tsne[[5]],FIG_DIR))),
               data.frame(sample=6, use=(get_pure_pop_idx(genes,pure_id[6],all_pure_pca[[6]],all_pure_tsne[[6]],FIG_DIR))),
               data.frame(sample=7, use=(get_pure_pop_idx(genes,pure_id[7],all_pure_pca[[7]],all_pure_tsne[[7]],FIG_DIR))),
               data.frame(sample=8, use=(get_pure_pop_idx(genes,pure_id[8],all_pure_pca[[8]],all_pure_tsne[[8]],FIG_DIR))),
               data.frame(sample=9, use=(get_pure_pop_idx(genes,pure_id[9],all_pure_pca[[9]],all_pure_tsne[[9]],FIG_DIR))),
               data.frame(sample=10,use=(get_pure_pop_idx(genes,pure_id[10],all_pure_pca[[10]],all_pure_tsne[[10]],FIG_DIR))),
               data.frame(sample=10,use=(get_pure_pop_idx(genes,pure_id[11],all_pure_pca[[10]],all_pure_tsne[[10]],FIG_DIR))))
pure_select_11<-lapply(1:length(sub_idx),function(i) {subsampled_purified_mats[[sub_idx[[i]]$sample[1]]][sub_idx[[i]]$use,]})
# -------------------------------------------------
# build correlation heatmap
# this produces Supp. Fig. 7 in the manuscript
# -------------------------------------------------
pure_use_genes<-which(colSums(do.call(rBind,lapply(pure_select_11,function(x) x)))>1)
pure_select_use_genes<-lapply(1:length(pure_select_11),function(i) pure_select_11[[i]][,pure_use_genes])
pure_avg<-do.call(rbind,lapply(pure_select_use_genes,.train_multinomial))
pure_select_11_cor<-cor(t(pure_avg),method='spearman')
rownames(pure_select_11_cor)<-pure_id
pheatmap(pure_select_11_cor,color=colorRampPalette(c("gray100", 'gray30'))(20),cluster_rows=T,cluster_cols=T)
# ------------------------------------------
# save mean expression for the 11 types
# ------------------------------------------
pure_11 <- list(pure_id=pure_id, pure_avg=pure_avg, pure_use_genes=pure_use_genes)
saveRDS(pure_11,file=file.path(DATA_DIR, "all_pure_select_11types.rds"))
