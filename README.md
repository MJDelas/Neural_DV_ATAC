# Neural DV ATAC
Scripts for manuscript 

Custom code to reproduce all figures. 

All figures generated using R can be visualized as markdown (.md files). 


Make links here: 
- [Test script](1_test.md)
- CaTS-ATAC sequencing processing (using [nf-core](https://nf-co.re/atacseq)): [script](sh/run_cats-atac.sh)
- RNA-seq processing (using [nf-core](https://nf-co.re/rnaseq)): [script](sh/run_rnaseq.sh) followed by [removal of PCR duplicates](R/R_geneCounts.R)
- Published ChIP-seq processing (using nf-core)
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