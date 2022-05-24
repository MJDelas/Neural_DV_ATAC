#!/usr/bin/env Rscript
##### R script to count features on all samples and print results 


##########
# Load packages

library(Rsubread)
library(dplyr)
library(tidyr)
library(tibble)


###########
## load inputs

targetDir="/camp/lab/briscoej/working/Joaquina/hpc_camp/RNA_DV_mm10ercc/results/markDuplicates"
#ann=read.delim("/camp/lab/briscoej/working/Joaquina/hpc_camp/delasj_reference_files/mm10.refGene_wHeader.gtf",
#               header = T, stringsAsFactors=FALSE)             


#########
# run features counts


counts_bams <- lapply(list.files(path = targetDir, pattern="*_read1Aligned.sortedByCoord.out.markDups.bam$", full.names=TRUE),function(x){
  data_all <- featureCounts(x,
                        annot.ext="/camp/lab/briscoej/working/Joaquina/hpc_camp//delasj_reference_files/mm10_gtf_header/mm10.refGene_wHeader_wERCC.gtf", # 
                        isGTFAnnotationFile=TRUE,
                        GTF.featureType=,
                        GTF.attrType="gene_id",
                        isPairedEnd=TRUE, # paired end
                        strandSpecific=0, # default for GTF
                        ignoreDup=TRUE # this is my choice to remove PCR duplicates
                  )
  data_all
  })

# get the feature counts
counts_bams_features <- lapply(c(1:length(counts_bams)), function(z){
  merge(counts_bams[[z]]$annotation,counts_bams[[z]]$counts, by.x="GeneID", by.y="row.names")
})

counts_features_clean <- Reduce(function(...)
  merge(...,by=c("GeneID","Chr","Start","End","Strand","Length")),counts_bams_features)

write.table(counts_features_clean, file = "Counts_DupR.featureCounts_custom.txt",
            row.names = F, quote = F)

# get the summaries
# Initial list:
counts_bams_summaries <- list()
for(i in 1:length(counts_bams)){
  counts_bams_summaries[[length(counts_bams_summaries)+1]] <- counts_bams[[i]]$stat
}

counts_summaries_clean <- Reduce(function(...)
    merge(...,by="Status"),counts_bams_summaries)

write.table(counts_summaries_clean, file = "Counts_DupR.featureCounts_custom.txt.summary",
            row.names = F, quote = F)

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

RLogFile <- "R_sessionInfo.log"
if (file.exists(RLogFile) == FALSE) {
    sink(RLogFile)
    a <- sessionInfo()
    print(a)
    sink()
}

