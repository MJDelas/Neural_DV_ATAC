export NXF_WORK="/camp/svc/scratch/briscoej/joaquina/nf-core_atacseq"

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/20.01.0
ml Singularity/2.6.0-foss-2016b


## UPDATE PIPELINE
nextflow pull nf-core/atacseq

## RUN PIPELINE
nextflow run nf-core/atacseq \
    --input  input_files/design_cats-atac.csv \
    --genome mm10 \
    --gtf '/camp/lab/briscoej/working/Joaquina/hpc_camp/delasj_reference_files//mm10_gtf_noheader/mm10.refGene.gtf' \
    --macs_fdr 0.00001 \
    --skip_diff_analysis \
    --min_reps_consensus 2 \
    --email joaquina.delas@crick.ac.uk \
    -profile crick \
    -r 1.2.0 \
    -resume