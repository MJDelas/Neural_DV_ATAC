
## Processing of published ChIP-seq summary

### Neural embryoids with SAG: NKX6.1, NKX2.2, OLIG2

Source: Nishi et al (2015) 

PMID: 26293298, DOI: 10.1242/dev.124636

GEO: GSE65462

`nextflow run nf-core/chipseq --input design.csv --single_end --genome mm10 --email joaquina.delas@crick.ac.uk -profile crick -r 1.1.0 -resume`

### Neural embryoids with SAG: SOX2, FOXA2

Source: Peterson et al.(2012) 

PMID: 23249739, DOI: 10.1101/gad.207142.112

GEO: GSE42594

```
nextflow run nf-core/chipseq \
        --input design.csv \
        --single_end \
        --genome mm10 \
        --email joaquina.delas@crick.ac.uk \
        -profile crick \
        -r 1.1.0 \
	--skip_diff_analysis \
        -resume
``` 

Available trackhubs: instructions in the [briscoelab wiki](https://briscoelab.github.io/wiki/BriscoeLab_trackhubs.html)


### Endoderm differentiation: FOXA2

Source: Cernilogar et al (2019)

PMID: PMID31350899, DOI: 10.1093/nar/gkz627

GEO: GSE116258

Due to lack of input, these sampels were processed with the ATAC-seq pipeline, as follows:

```
nextflow run nf-core/atacseq \
        --input design.csv \
        --single_end \
        --mito_name false \
        --genome mm10 \
        --email joaquina.delas@crick.ac.uk \
        -profile crick \
        -r 1.1.0 \
        -resume
```

