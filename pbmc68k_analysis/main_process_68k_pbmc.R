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
library(plyr)
library(data.table)
library(pheatmap)
# -------------------------------------
# specify paths and load functions
# -------------------------------------
DATA_DIR <- "PATH_TO_DATA_DIRECTORY"        # SPECIFY HERE
PROG_DIR <- "PATH_TO_PROGRAM_DIRECTORY"     # SPECIFY HERE
RES_DIR  <- "PATH_TO_RESULT_DIRECTORY"      # SPECIFY HERE
source(file.path(PROG_DIR,'util.R'))
# ------------------------------------------------------------
# load 68k PBMC data, 11 purified PBMC data and meta-data
# ------------------------------------------------------------
pbmc_68k <- readRDS(file.path(DATA_DIR,'pbmc68k_data.rds'))
pure_11 <- readRDS(file.path(DATA_DIR,'all_pure_select_11types.rds'))
all_data <- pbmc_68k$all_data
purified_ref_11 <- load_purified_pbmc_types(pure_11,pbmc_68k$ens_genes)
# --------------------------------------------------------------------------------------
# normalize by RNA content (umi counts) and select the top 1000 most variable genes
# --------------------------------------------------------------------------------------
m<-all_data[[1]]$hg19$mat
l<-.normalize_by_umi(m)
m_n<-l$m
df<-.get_variable_gene(m_n)
disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[1000]
df$used<-df$dispersion_norm >= disp_cut_off
# --------------------------------------------------
# plot dispersion vs. mean for the genes
# this produces Supp. Fig. 5c in the manuscript
# --------------------------------------------------
ggplot(df,aes(mean,dispersion,col=used))+geom_point(size=0.5)+scale_x_log10()+scale_y_log10()+
  scale_color_manual(values=c("grey","black"))+theme_classic()
# --------------------------------------------
# use top 1000 variable genes for PCA
# --------------------------------------------
set.seed(0)
m_n_1000<-m_n[,head(order(-df$dispersion_norm),1000)]
pca_n_1000<-.do_propack(m_n_1000,50)
# --------------------------------------------
# generate 2-D tSNE embedding
# this step may take a long time
# --------------------------------------------
tsne_n_1000<-Rtsne(pca_n_1000$pca,pca=F)
tdf_n_1000<-data.frame(tsne_n_1000$Y)
# ---------------------------------------------------------------------------------------------------------------------------
# assign IDs by comparing the transcriptome profile of each cell to the reference profile from purified PBMC populations
# this produces Fig. 3j in the manuscript
# ---------------------------------------------------------------------------------------------------------------------------
m_filt<-m_n_1000
use_genes_n<-order(-df$dispersion_norm)
use_genes_n_id<-all_data[[1]]$hg19$gene_symbols[l$use_genes][order(-df$dispersion_norm)]
use_genes_n_ens<-all_data[[1]]$hg19$genes[l$use_genes][order(-df$dispersion_norm)]
z_1000_11<-.compare_by_cor(m_filt,use_genes_n_ens[1:1000],purified_ref_11)
# reassign IDs, as there're some overlaps in the purified pbmc populations
test<-.reassign_pbmc_11(z_1000_11)
cls_id<-factor(colnames(z_1000_11)[test])
tdf_n_1000$cls_id<-cls_id
# adjust ordering of cells for plotting aesthetics
tdf_mod <- tdf_n_1000[tdf_n_1000$cls_id!='CD4+/CD45RA+/CD25- Naive T',]
tdf_mod <- rbind(tdf_mod,tdf_n_1000[tdf_n_1000$cls_id=='CD4+/CD45RA+/CD25- Naive T',])
tdf_mod_2 <- tdf_mod[tdf_mod$cls_id!='CD56+ NK',]
tdf_mod_2 <- rbind(tdf_mod_2,tdf_mod[tdf_mod$cls_id=='CD56+ NK',])
ggplot(tdf_mod_2,aes(X1,X2,col=cls_id))+geom_point(size=0,alpha=1)+theme_classic()+.set_pbmc_color_11()
# --------------------------------------------------
# use k-means clustering to specify populations
# this produces Fig. 3b in the manuscript
# --------------------------------------------------
set.seed(0)
k_n_1000<-kmeans(pca_n_1000$pca,10,iter.max=150,algorithm="MacQueen")
tdf_n_1000$k<-k_n_1000$cluster
ggplot(tdf_n_1000,aes(X1,X2,col=as.factor(k)))+geom_point(size=0,alpha=0.6)+theme_classic()+
   scale_color_manual(values=c("#FB9A99","#FF7F00","yellow","orchid","grey",
                               "red","dodgerblue2","tan4","green4","#99c9fb"))
# -------------------------------------------------
# find cluster-specific genes and plot markers
# this produces Fig. 3c in the manuscript
# -------------------------------------------------
# plot heatmap with gene clusters and write the genes to file
.plot_heatmap_and_write_genes(n_clust=10, tsne_df=tdf_n_1000, mat=m_n_1000, genes=use_genes_n_id[1:1000],
                              dir=RES_DIR, topic="k10_marker_genes")
markers <- c("HLA-DRA","CD3D","CD8A","NKG7")
.plot_gene_expression(markers=markers, genes=all_data[[1]]$hg19$gene_symbols, mat=all_data[[1]]$hg19$mat,
                      tsne_df=tdf_n_1000, style='type1', title='pbmc68k', dir=RES_DIR)

# ----------------------------------------------
# identify sub-clusters in cluster #9
# this produces Fig. 3g-i in the manuscript
# ----------------------------------------------
c9_m_n_1000<-m_n_1000[tdf_n_1000$k==9,]
pca_c9<-.do_propack(c9_m_n_1000,10)
set.seed(0)
tsne_c9<-Rtsne(pca_c9$pca,pca=F)
km_c9<-kmeans(pca_c9$pca,3,iter.max=150,algorithm="MacQueen")
tdf_c9<-data.frame(tsne_c9$Y,k=km_c9$cluster)
.plot_heatmap_and_write_genes(n_clust=3, tsne_df=tdf_c9, mat=c9_m_n_1000, genes=use_genes_n_id[1:1000],
                              dir=RES_DIR, topic="myeloid_genes_3grps")
markers <- c("FCER1A","FCGR3A","S100A8","CD1C")
c9_m <- all_data[[1]]$hg19$mat[tdf_n_1000$k==9,]
.plot_gene_expression(markers=markers, genes=all_data[[1]]$hg19$gene_symbols, mat=c9_m,
                      tsne_df=tdf_c9, style='type2', title='cluster9', dir=RES_DIR)



