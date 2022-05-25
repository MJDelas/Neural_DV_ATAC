Clustering of accessible elements
================

# Clustering element dynamics by kmeans after filtering elements

``` r
rm(list=ls())

library(RColorBrewer)
library(tidyverse)
library(factoextra)
library(ggrepel)
library(ComplexHeatmap)
```

## Load data

Vsd with all the samples generated in
`cats-atac_2_deseq_PCA_heatmaps.md`
[here](NeuralDV_Rproject/cats-atac_2_deseq_PCA_heatmaps.md)

And the list of genes that we have filtered

``` r
# read the vsd data from DESeq
vsd_data <- read.csv("/Users/delasj/Documents/BriscoeLab/project_DV_ATAC_reproduce_analysis/outputs_cats-atac_2/consensus_peaks.mLb.vsd.csv") %>%
  remove_rownames() %>%
  column_to_rownames("X")

diffexpr_time_gate <- read.csv("/Users/delasj/Documents/BriscoeLab/project_DV_ATAC_reproduce_analysis/outputs_cats-atac_4/DiffAccessible_elements.csv", stringsAsFactors = FALSE)


vsd_sub <- vsd_data %>%
  rownames_to_column("Peakid") %>% filter(Peakid %in% diffexpr_time_gate$x) %>%
  column_to_rownames("Peakid")
```

## Colors and shapes

``` r
sorted.DayGate <- c("D3_NMP","D4_1","D4_2","D4_M","D4_3",
                    "D5_1","D5_2","D5_M","D5_3",
                    "D6_1","D6_2","D6_M","D6_3")

sorted.samples <- c("D3_0_NMP","D4_0_1","D4_10_1","D5_0_1","D5_10_1","D6_0_1","D6_10_1",
                   "D4_10_2","D4_100_2","D5_10_2","D5_100_2","D6_10_2","D6_100_2",
                   "D4_100_M","D4_500_M","D5_100_M","D5_500_M","D6_100_M","D6_500_M",
                   "D4_100_3","D4_500_3","D5_100_3","D5_500_3","D6_100_3","D6_500_3")

colorJD <- c("#878787","#6da4ba","#f0be56","#ec936f","#5bb357",
             "#477d92","#e5a114","#e3602b","#009640",
             "#2e525e","#9f7113","#ab4117","#044a23")

shapes4_manual = c(18,15,16,17) # these are block
shapes5_manual = c(25,21,22,23,24) # these are filled
shapes4_fill_manual = c(23,21,22,24)

# Annotated heatmap with selected colors
hm_colors = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
```

## Shape data

Calculate average

``` r
vsd_sub_ave <- vsd_sub %>%
  rownames_to_column("order") %>%
  gather(sample, vsd, starts_with("D")) %>%
  mutate(condition = gsub('_R[1-3]', '', sample)) %>%
  group_by(order, condition) %>%
  summarise(ave_vsd = mean(vsd)) %>%
  spread(condition, ave_vsd) %>%
  ungroup() %>%
  column_to_rownames("order")
```

    ## `summarise()` has grouped output by 'order'. You can override using the
    ## `.groups` argument.

## kmeans cluster and re-cluster

The number of clusters at both steps has been tested in multiple
iterations to find a stable grouping. Using seed for reproducibility.

The function will save all the results to the folder.

``` r
outdir=("/Users/delasj/Documents/BriscoeLab/project_DV_ATAC_reproduce_analysis/outputs_cats-atac_5_kmeans/")

seed_options <- 10

atac_clustering <- function(x){
  
  #set seed to you can reproduce
  set.seed(seed_options)
  
  ## kmeans
  vsd_kclusterd_2 <- kmeans(vsd_sub_ave,30,iter.max = 100, nstart = 50)
  clustering_qc <- data.frame("Tot.withinss"=vsd_kclusterd_2$tot.withinss,
                              "ifault"=vsd_kclusterd_2$ifault)
  
  write.table(clustering_qc, file=paste0(outdir,"Stats40_iter_",x,".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  vsd_kclusterd <- vsd_kclusterd_2
  
  # get cluster centers to plot and get distances
  clusters_k30_plot <- vsd_kclusterd$centers %>% as.data.frame() %>%
    rownames_to_column("cluster") %>%
    gather(Sample,centers, starts_with("D")) %>%
    separate(Sample, into = c("Day","SAG","Gate"), remove = FALSE) %>%
    mutate(DayGate=paste(Day,Gate,sep = "_")) %>%
    mutate(DayGate = factor(DayGate, levels = sorted.DayGate),
           Sample = factor(Sample, levels = sorted.samples)) 


# plot center per cluster
  plot_gg <- ggplot(clusters_k30_plot, aes(x=Sample,y=centers, group=cluster))+
    geom_line() +
    geom_point(aes(color=DayGate)) +
    scale_color_manual(values = colorJD) +
    facet_wrap(~cluster) +
    theme_classic(base_size=8) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
    ggsave(paste0(outdir,"Plots_30Cluster_expression_Iter_",x,".pdf"), plot=plot_gg,
           width=15, height=5, units="in", useDingbats=FALSE)

  
  #plot_centers
# get cluster centers
  clusters_k30 <- vsd_kclusterd$centers
  write.table(clusters_k30, file=paste0(outdir,"Cluster30Centers_iter_",x,".txt"), sep = "\t", quote = FALSE, row.names = FALSE)

  
  clusters_intervals_k30 <- as.data.frame(vsd_kclusterd$cluster) %>%
    rownames_to_column()
  colnames(clusters_intervals_k30) <- c("Interval","Cluster")


  clusters_intervals_k30 <- clusters_intervals_k30[order(clusters_intervals_k30[,2] ),]

  write.csv(clusters_intervals_k30, file=paste0(outdir,"Genes_Cluster30Centers_iter_",x,".csv"), quote = FALSE, row.names = FALSE)
  
  ## Computing correlation based distances
  dist_cluster <- get_dist(clusters_k30, method = "pearson")
  dist_cluster_hm <- as.matrix(dist_cluster)
  
  #calculate hclust
  hcl=hclust(dist_cluster)
  
  # print hierchical cluster and cutree for each iteration
  pdf(paste0(outdir,"Plot_hclust_Iter_",x,".pdf")) 

  k=9 # how many ReClusters
  plot(hcl, hang= -1) 
  reclusters <- rect.hclust(hcl, k, border = "red")
  
  dev.off() 


# visualize the re-clustering 
  # visualize the clusters and mega clusters
  ReClustering_dirty <- do.call(rbind.data.frame, reclusters)
  
  ReClustering <- ReClustering_dirty %>% rownames_to_column("ReCluster") %>% 
    gather("col_del", "cluster", 2:(ncol(ReClustering_dirty)+1)) %>%
    dplyr::select("cluster","ReCluster") %>%
    unique() %>%
    arrange(cluster) %>%
    remove_rownames() %>%
    column_to_rownames("cluster")
  
  #color clusters
  
  Nclusters <- c(1:k) %>% as.character() # make vector
  Ncolors <- colorRampPalette(brewer.pal(12, "Set3"))(k) # get colors

  
  names(Ncolors) <- Nclusters # named vector
  
  hm_colors = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
  
  ann_color_cluster <- list(
    ReCluster = c(Ncolors))
  
  # make Complex heatmap annotation 
  heatmap_ann_row <- rowAnnotation(df=ReClustering, col=ann_color_cluster)
  heatmap_ann <- HeatmapAnnotation(df=ReClustering, col=ann_color_cluster,
                                   show_legend = FALSE)
  
  
  
  # print hierchical cluster and cutree for each iteration
  pdf(paste0(outdir,"Plot",k,"_HeatmapReCluster_Iter_",x,".pdf"), width = 10, height = 9) 
  
  
  # now with colors of the reclustering
  draw(Heatmap(cor(t(clusters_k30)), name="Cluster cor", col=hm_colors,
          left_annotation = heatmap_ann_row, top_annotation = heatmap_ann,
          clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
          clustering_method_columns = 'complete',
          clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
          clustering_method_rows = 'complete',
          column_dend_height = unit(2, "cm"), row_dend_width = unit(2, "cm")))
  
  dev.off()
  

 # export cluster and recluter
  reclusteres_cluster_order <- vsd_kclusterd$cluster %>% as.data.frame() %>%
    rename("cluster"=".") %>%
    rownames_to_column("order") %>%
    merge(ReClustering %>% rownames_to_column("cluster")) %>%
    arrange(order)
  write.csv(reclusteres_cluster_order, file = paste0(outdir,"Output",k,"_ReClusters_Iter_",x,".csv"), 
            row.names = FALSE, quote = FALSE)
  


  ## get the number of elements per recluster for violin plot representation 
  
  n.hcluster <- reclusteres_cluster_order %>%
    group_by(ReCluster) %>%
    summarise(n = n()) %>%
    mutate(label = sprintf('n = %s', n))
  
  
  
  ## plot zcores per ReCluster

  # scale the vsd data
  heat <- t(scale(t(vsd_sub_ave)))
  
  ## Violin   
  reclustered_hclust_zcoresconditions <- heat %>% as.data.frame() %>%
    rownames_to_column("order") %>%
    merge(reclusteres_cluster_order, by="order" ) %>%
    gather(sample, zscore, starts_with("D")) %>%
    separate(sample, into = c("Day","SAG","Gate"), remove = FALSE) %>%
    mutate(DayGate=paste(Day,Gate,sep = "_")) %>%
    mutate(DayGate = factor(DayGate, levels = sorted.DayGate),
           sample = factor(sample, levels = sorted.samples)) 
  

  plot_violins <- ggplot(reclustered_hclust_zcoresconditions, aes(x = sample, y = zscore, fill=DayGate)) + 
    geom_violin(size=0.2) +
    scale_fill_manual(values = colorJD) +
    facet_wrap(~ReCluster) + xlab('') + ylab('Z-Score') +
    geom_text_repel(data = n.hcluster, aes(x = 0.75, y = 4, label = label), inherit.aes = F) +
    geom_hline(yintercept = 0, col = 'red', linetype = 'dashed') +
    theme_bw(base_size=10)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = 'none')

    ggsave(paste0(outdir,"Plots_",k,"ReCluster_expression_Iter_",x,".pdf"), plot=plot_violins,
           width=10, height=5, units="in", useDingbats=FALSE)
}
```

Run the function with the chosen seed.

``` r
atac_clustering(seed_options)
```

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1083000)

``` r
sessionInfo()
```

    ## R version 3.6.3 (2020-02-29)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Catalina 10.15.7
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
    ## 
    ## attached base packages:
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] ComplexHeatmap_2.2.0 ggrepel_0.9.1        factoextra_1.0.7    
    ##  [4] forcats_0.5.1        stringr_1.4.0        dplyr_1.0.8         
    ##  [7] purrr_0.3.4          readr_2.1.2          tidyr_1.2.0         
    ## [10] tibble_3.1.6         ggplot2_3.3.5        tidyverse_1.3.1     
    ## [13] RColorBrewer_1.1-3  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.8.3        circlize_0.4.14     lubridate_1.8.0    
    ##  [4] png_0.1-7           assertthat_0.2.1    digest_0.6.29      
    ##  [7] utf8_1.2.2          R6_2.5.1            cellranger_1.1.0   
    ## [10] backports_1.4.1     reprex_2.0.1        evaluate_0.15      
    ## [13] httr_1.4.2          pillar_1.7.0        GlobalOptions_0.1.2
    ## [16] rlang_1.0.2         readxl_1.4.0        rstudioapi_0.13    
    ## [19] GetoptLong_1.0.5    rmarkdown_2.13      labeling_0.4.2     
    ## [22] textshaping_0.3.6   munsell_0.5.0       broom_0.7.12       
    ## [25] compiler_3.6.3      modelr_0.1.8        xfun_0.30          
    ## [28] systemfonts_1.0.4   pkgconfig_2.0.3     shape_1.4.6        
    ## [31] htmltools_0.5.2     tidyselect_1.1.2    fansi_1.0.3        
    ## [34] crayon_1.5.1        tzdb_0.3.0          dbplyr_2.1.1       
    ## [37] withr_2.5.0         jsonlite_1.8.0      gtable_0.3.0       
    ## [40] lifecycle_1.0.1     DBI_1.1.2           magrittr_2.0.3     
    ## [43] scales_1.1.1        cli_3.2.0           stringi_1.7.6      
    ## [46] farver_2.1.0        fs_1.5.2            xml2_1.3.3         
    ## [49] ragg_1.2.2          ellipsis_0.3.2      generics_0.1.2     
    ## [52] vctrs_0.4.0         rjson_0.2.20        tools_3.6.3        
    ## [55] glue_1.6.2          hms_1.1.1           parallel_3.6.3     
    ## [58] fastmap_1.1.0       yaml_2.3.5          clue_0.3-60        
    ## [61] colorspace_2.0-3    cluster_2.1.2       rvest_1.0.2        
    ## [64] knitr_1.38          haven_2.4.3
