# Neural DV ATAC
Scripts for manuscript 

Custom code to reproduce all figures. 

All figures generated using R can be visualized as markdown (.md files). 


## Sequencing processing
- CaTS-ATAC sequencing processing (using [nf-core](https://nf-co.re/atacseq)): 
    - [pipeline script](sh/run_cats-atac.sh)
    - [fdr interval filtering](NeuralDV_Rproject/cats-atac_1_filter_fdr.md)
- RNA-seq processing (using [nf-core](https://nf-co.re/rnaseq)): 
    - [pipeline script](sh/run_rnaseq.sh) 
    - [removal of PCR duplicates](R/R_geneCounts.R)
- Published ChIP-seq processing (using [nf-core](https://nf-co.re/chipseq)) : 
    - [summary and scripts](docs/chip-seq_processing.md)
- Bulk ATAC-seq processing (using [nf-core](https://nf-co.re/atacseq)):
    - [pipeline script](sh/run_atac-ifoxa2.sh)

## CaTS-ATAC analysis
- Comparison of accessibility across conditions (using [Deseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
    - [Differential accessibility between SAG concentrations](NeuralDV_Rproject/cats-atac_3_cross_condition_diffacc.md)
- Global accessibility analysis (using [Deseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)):
    - [PCA and correlation heatmaps](NeuralDV_Rproject/cats-atac_2_deseq_PCA_heatmaps.md)
    - [Differentially accessible elements](NeuralDV_Rproject/cats-atac_4_DiffAcc_filter.md)
    - [kmeans clustering](NeuralDV_Rproject/cats-atac_5_kmeans.md)
- Heatmap visualization of ATAC-seq and ChIP-seq (using [deeptools](https://deeptools.readthedocs.io/en/develop/)):
    - [CaTS-ATAC and ChIP](sh/run_deeptools_atac_chip.sh)
    - [FOXA2 ChIP](sh/run_deeptools_atac_Foxa2chip.sh)
- Metaplot of ChIPseq over ATAC peaks, and ATAC clusters over ChIPseq peaks

## RNA-seq analysis
- Normalized counts export and key genes

## Footprinting analysis
- Footprinting score calculations (using [TOBIAS](https://github.com/loosolab/TOBIAS)) 
    - [Scripts](docs/tobias_allscripts.md)
- Most variable footprints across p0-M
- Most variable footprints across all: p3-specific

## Inducible Foxa2 ATACseq
- Differential accessibility in inducible Foxa2 (using [Deseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
    - [Analysis and plots](NeuralDV_Rproject/ifoxa2-atac_1_maplot.md)
- Intersect 
- Plotting intersects 