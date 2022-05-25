# Neural DV ATAC
Scripts for manuscript 

Custom code to reproduce all figures. 

All figures generated using R can be visualized as markdown (.md files). 


Make links here: 
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
- Comparison of accessibility across conditions
    - [Differential accessibility between SAG concentrations](NeuralDV_Rproject/cats-atac_3_cross_condition_diffacc.md)
- Global accessibility analysis: PCA and differentially accessible elements selection
    - [PCA and correlation heatmaps](NeuralDV_Rproject/cats-atac_2_deseq_PCA_heatmaps.md)
- Heatmap visualization of ATAC-seq and ChIP-seq (using [deeptools](https://deeptools.readthedocs.io/en/develop/)):
    - script
- Metaplot of ChIPseq over ATAC peaks, and ATAC clusters over ChIPseq peaks
- Footprinting score calculations (using [TOBIAS](https://github.com/loosolab/TOBIAS)) 
- Most variable footprints across p0-M
- Most variable footprints across all: p3-specific
- Differential accessibility in inducible Foxa2
    - do next 
- Intersect 
- Plotting intersects 