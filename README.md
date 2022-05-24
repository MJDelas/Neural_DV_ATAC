# Neural DV ATAC
Scripts for manuscript 

Custom code to reproduce all figures. 

All figures generated using R can be visualized as markdown (.md files). 


Make links here: 
- [Test script](1_test.md)
- CaTS-ATAC sequencing processing (using [nf-core](https://nf-co.re/atacseq)): 
    - [pipeline script](sh/run_cats-atac.sh)
    - [fdr interval filtering](Neural_DV_ATAC/blob/main/NeuralDV_Rproject/cats-atac_1_filter_fdr.md)
- RNA-seq processing (using [nf-core](https://nf-co.re/rnaseq)): 
    - [pipeline script](sh/run_rnaseq.sh) 
    - [removal of PCR duplicates](R/R_geneCounts.R)
- Published ChIP-seq processing (using [nf-core](https://nf-co.re/chipseq)) : [summary and scripts](docs/chip-seq_processing.md)
- Bulk ATAC-seq processing (using nf-core)
- Comparison of accessibility in cell type across SAG concentrations
- Global accessibility analysis: PCA and differentially accessible elements selection
- Heatmap visualization of ATAC-seq and ChIP-seq
- Footprinting score calculations (using TOBIAS) 
- Most variable footprints across p0-M
- Most variable footprints across all: p3-specific
- Differential accessibility in inducible Foxa2
- Intersect 
- Plotting intersects 