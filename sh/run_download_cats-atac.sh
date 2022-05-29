export NXF_WORK="/camp/svc/scratch/briscoej/joaquina/nf-core_fetch"

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/21.10.3
ml Singularity/2.6.0-foss-2016b


## UPDATE PIPELINE
nextflow pull nf-core/fetchngs

## RUN PIPELINE
nextflow run nf-core/fetchngs \
    --input  geo_cats.csv \
    --outdir cats-atac \
    --email joaquina.delas@crick.ac.uk \
    -profile crick \
    -r 1.6 \
    -resume