#!/bin/bash 
export NXF_SINGULARITY_CACHEDIR=/camp/lab/briscoej/working/Joaquina/hpc_camp/sing/atac
export NXF_WORK=/camp/svc/scratch/briscoej/joaquina/nf-core_atac_iFoxa2

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/21.04.0
ml Singularity/2.6.0-foss-2016b
ml CAMP_proxy

## UPDATE PIPLINE
nextflow pull nf-core/atacseq

## RUN PIPELINE
nextflow run nf-core/atacseq \
	--input input_files/design_atac_ifoxa2.csv \
        --genome mm10 \
        --macs_fdr 0.00001 \
    	--skip_diff_analysis \
    	--min_reps_consensus 2 \
        --email joaquina.delas@crick.ac.uk \
        -profile crick \
        -r 1.2.0 \
        -resume