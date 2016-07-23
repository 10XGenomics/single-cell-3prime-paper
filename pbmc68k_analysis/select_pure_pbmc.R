#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#

# ---------------------------
# curate pure cell types
# ---------------------------
get_pure_pop_idx <- function(all_genes_symbols,id,pure_type_pca,pure_type_tsne,DIR=NULL) {
  known_ids <- c("CD19+ B","CD14+ Monocyte","Dendritic","CD56+ NK","CD34+","CD4+/CD25 T Reg","CD4+/CD45RA+/CD25- Naive T",
                "CD4+/CD45RO+ Memory","CD4+ T Helper2","CD8+/CD45RA+ Naive Cytotoxic","CD8+ Cytotoxic T")
  if (id %in% known_ids)  {
    k <- 1; sel <- c(1) # default is to select all populations
  } else {
    stop('Unknown cell type:',id,'\n')
  }
  if (id=="CD19+ B") { check_genes <- c("CD74","CD27")} 
  if (id=="CD14+ Monocyte") { check_genes <- c("FTL");      k <- 5; sel <- c(2,4) }  # select based on FTL
  if (id=="Dendritic") { check_genes <- c("CLEC9A","CD1C"); k <- 5; sel <- c(5) }      # select based on CLEC9A
  if (id=="CD34+")  { check_genes <- "CD34";                k <- 8; sel <- c(1,2,3,7) }# select based on CD34
  if (id=="CD56+ NK") {check_genes <- c("CD3D","NKG7") }
  if (id=="CD4+/CD25 T Reg") { check_genes <- "NKG7" }
  if (id=="CD4+/CD45RA+/CD25- Naive T") { check_genes <- "NKG7" }
  if (id=="CD4+/CD45RO+ Memory") { check_genes <- "NKG7" }
  if (id=="CD4+ T Helper2") { check_genes <- "NKG7" }
  if (id=="CD8+/CD45RA+ Naive Cytotoxic") { check_genes <- "NKG7"} 
  if (id=="CD8+ Cytotoxic T") { check_genes <- c("CD8A","CD3D")} 

  .plot_gene_markers(pure_type_tsne,pure_type_pca,check_genes,all_genes_symbols,id,DIR)
  return(.pick_pop_clusters(pure_type_tsne,pure_type_pca,k,sel,all_genes_symbols,id,DIR))
}
# ----------------------------------------------------
# figure plotting for purified populations
# this produces Supp. Fig. 7b-k of the manuscript
# ----------------------------------------------------
.plot_k_cluster <- function(pure_tsne,cluster,id,dir,marker,no_col=FALSE) {
  tsne_df<-data.frame(pure_tsne$Y)
  if (no_col) { # no color no color bar
    p <- ggplot(tsne_df,aes(X1,X2)) 
  } else {
    p <- ggplot(tsne_df,aes(X1,X2,col=as.factor(cluster))) + 
      scale_colour_brewer(name='val',palette="Set2")
  }
  if (marker=='Sel') {# compute fraction of selected cells for refernce
    fraction <- sum(cluster)/length(cluster) * 100
    tit <- paste(id,' (',fraction,'%)',sep="")
  } else {
    tit <- id
  }
  p <- p + xlab("") + ylab("") + geom_point(size=0.1)+theme_bw()+ggtitle(tit)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text.x = element_blank(), axis.text.y = element_blank()) 
  if (!is.null(dir)) {
    if (no_col) {
      .save_figure(p,dir,id,marker,width=3.2,height=3.5)
    } else {
      .save_figure(p,dir,id,marker,width=4.0,height=3.5)
    }
  }
}
.plot_gene_markers <- function(pure_tsne,pure_pca,markers,all_genes_symbols,id,dir) {
  if (is.null(dir)) {
    return()
  }
  tsne_df<-data.frame(pure_tsne$Y)
  m_n<-pure_pca$m 
  for (marker in markers) {
    marker_val<-m_n[,match(marker,all_genes_symbols[pure_pca$use_genes])]
    midpnt <- 0
    if (marker=='FTL') {midpnt <- -2} 
    if (marker=='CLEC9A') {midpnt <- 1} 
    if (marker=='NKG7') {midpnt <- 6} 
    if (marker=='CD34') {midpnt <- 0.5} 
    if (marker=='NKG7') {midpnt <- -1} 
    if (marker=='CD3D') {midpnt <- -1} 
    if (marker=='CD27') {midpnt <- 1} 
    if (marker=='CD74') {midpnt <- -1} 
    p <- ggplot(tsne_df,aes(X1,X2,col=marker_val))+geom_point(size=0.3)+
      theme_bw()+ggtitle(paste(id,'(',marker,')'))+xlab("")+ylab("")+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.ticks = element_blank(), 
            axis.text.x = element_blank(), axis.text.y = element_blank())+
      scale_colour_gradient2(low="skyblue4",high="red",mid="white",midpoint=midpnt,name='val',
                             guide=guide_colorbar(barwidth = 1,barheight = 5,limits=c(0,max(marker_val)))) 
    .save_figure(p,dir,id,marker,width=4.0,height=3.5)
  }
}
.pick_pop_clusters <- function(pure_tsne,pure_pca,k,select_cluster_ids,all_genes_symbols,id,dir) {
  set.seed(0) # use k-means clustering to help select regions of cells
  k_means<-kmeans(pure_pca$pca,k)
  selection <- k_means$cluster %in% select_cluster_ids
  no_colour <- (k == 1)
  if (!is.null(dir)) {
    .plot_k_cluster(pure_tsne,k_means$cluster,id,dir,'Kmeans',no_col=no_colour)
    .plot_k_cluster(pure_tsne,selection,id,dir,'Sel')
  }
  return(selection)
}
