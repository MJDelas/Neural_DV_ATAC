Filter intervals below threshold
================

``` r
rm(list=ls())

library(tidyverse)
```

The default atacseq pipeline threshold includes very small peaks. The
option –macs\_fdr 0.00001 did not filter as expected so subsequently
filter peaks as follows.

### Set dirs

``` r
workingdir="/Users/delasj/Documents/BriscoeLab/project_DV_ATAC_reproduce_analysis/"
subworkinput="inputs_cats-atac_1/"
outdir="outputs_cats-atac_1/"

ifelse(!dir.exists(file.path(workingdir,outdir)), dir.create(file.path(workingdir,outdir)), "Directory exists")
```

    ## [1] "Directory exists"

## Load data

The boolean file from nf-core/atacseq
`cats_atac/results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.boolean.annotatePeaks.txt`

``` r
peaks_boolean <- read.delim(file=paste0(workingdir,subworkinput,"consensus_peaks.mLb.clN.boolean.annotatePeaks.txt"),header=TRUE, stringsAsFactors = FALSE)
```

Clean empty columns and column names

``` r
peaks_boolean_clean <- peaks_boolean %>%
  select(chr,start,end,interval_id,num_peaks,num_samples, matches("qval$")) %>%
  gather(variable, qval, matches("^D[0-9]"))

peaks_boolean_clean$qval <- as.numeric(peaks_boolean_clean$qval)
```

    ## Warning: NAs introduced by coercion

## Filter peaks that have qval \> 0.00001

Use qval \>5 as it is -log10

Peak is completely eliminated if it does not pass filter in any sample

``` r
peaks_filter_gather <- peaks_boolean_clean %>%
  filter(qval > 5)

peaks_filter_wide <- peaks_filter_gather %>%
  select(interval_id,variable,qval) %>%
  spread(variable,qval)

write.table(peaks_filter_wide, file =paste0(workingdir,outdir,"consensus_peakfdr_filtered.csv"), quote = FALSE, row.names = FALSE)
```

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] forcats_0.5.1   stringr_1.4.0   dplyr_1.0.8     purrr_0.3.4    
    ## [5] readr_2.1.2     tidyr_1.2.0     tibble_3.1.6    ggplot2_3.3.5  
    ## [9] tidyverse_1.3.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.2 xfun_0.30        haven_2.4.3      colorspace_2.0-3
    ##  [5] vctrs_0.4.0      generics_0.1.2   htmltools_0.5.2  yaml_2.3.5      
    ##  [9] utf8_1.2.2       rlang_1.0.2      pillar_1.7.0     glue_1.6.2      
    ## [13] withr_2.5.0      DBI_1.1.2        dbplyr_2.1.1     modelr_0.1.8    
    ## [17] readxl_1.4.0     lifecycle_1.0.1  munsell_0.5.0    gtable_0.3.0    
    ## [21] cellranger_1.1.0 rvest_1.0.2      evaluate_0.15    knitr_1.38      
    ## [25] tzdb_0.3.0       fastmap_1.1.0    fansi_1.0.3      broom_0.7.12    
    ## [29] backports_1.4.1  scales_1.1.1     jsonlite_1.8.0   fs_1.5.2        
    ## [33] hms_1.1.1        digest_0.6.29    stringi_1.7.6    grid_3.6.3      
    ## [37] cli_3.2.0        tools_3.6.3      magrittr_2.0.3   crayon_1.5.1    
    ## [41] pkgconfig_2.0.3  ellipsis_0.3.2   xml2_1.3.3       reprex_2.0.1    
    ## [45] lubridate_1.8.0  assertthat_0.2.1 rmarkdown_2.13   httr_1.4.2      
    ## [49] rstudioapi_0.13  R6_2.5.1         compiler_3.6.3
