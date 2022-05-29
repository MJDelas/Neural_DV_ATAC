## Data download from GEO

Use [nf-core/fetchngs](https://nf-co.re/fetchngs) to download all three datasets from GEO.

### CaTS-ATAC

```
nextflow run nf-core/fetchngs \
    --input  input_files/geo_cats.txt \
    --outdir cats-atac \
    --email joaquina.delas@crick.ac.uk \
    -profile crick \
    -r 1.6 \
    -resume
```

`geo_cats.txt`
```
GSE204690
```

### Fixed RNA-seq

```
nextflow run nf-core/fetchngs \
    --input  geo_fixedRNA.txt \
    --outdir fixed-rna \
    --email joaquina.delas@crick.ac.uk \
    -profile crick \
    -r 1.6 \
    -resume
```

`geo_fixedRNA.txt`
```
GSE204920
```

### Inducible Foxa2 bulk ATAC-seq

```
nextflow run nf-core/fetchngs \
    --input  geo_iFoxa2_ATAC.txt \
    --outdir ifoxa2-atac \
    --email joaquina.delas@crick.ac.uk \
    -profile crick \
    -r 1.6 \
    -resume
```

`geo_iFoxa2_ATAC.txt`
```
GSE204664
```