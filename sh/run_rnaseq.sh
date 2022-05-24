export NXF_WORK="/camp/svc/scratch/briscoej/joaquina/nf-core_rnaseq"

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/20.01.0
ml Singularity/2.6.0-foss-2016b


## UPDATE PIPLINE
nextflow pull nf-core/rnaseq

## RUN PIPELINE
nextflow run nf-core/rnaseq \
    -profile crick \
    --reads '/camp/lab/briscoej/working/Joaquina/hpc_camp/RNA_Full/fastqgz_clean/D*read{1,2}.fastq.gz' \
    --star_index '/camp/lab/briscoej/working/Joaquina/hpc_camp/delasj_reference_files/star_ercc_mm10_genome' \
    --gtf '/camp/lab/briscoej/working/Joaquina/hpc_camp/delasj_reference_files/mm10_gtf_noheader/mm10.refGene_wERCC.gtf' \
    --fc_group_features_type gene_id \
    --fc_extra_attributes gene_id \
    --email joaquina.delas@crick.ac.uk \
    -r 1.4.2 \
    -resume