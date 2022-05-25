#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --job-name=deepmatrix_endochip
#SBATCH --mem=16G
#SBATCH -n 4
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joaquina.delas@crick.ac.uk
#SBATCH --output=deepmatrix_endochip.o
#SBATCH --error=deepmatrix_endochip.e

computeMatrix reference-point -S cats_atacATAC_FullFilter_qval/results/bwa/mergedReplicate/bigwig/D3_0_NMP.mRp.clN.bigWig \
				cats_atacATAC_FullFilter_qval/results/bwa/mergedReplicate/bigwig/D6_0_1.mRp.clN.bigWig \
				cats_atacATAC_FullFilter_qval/results/bwa/mergedReplicate/bigwig/D6_10_2.mRp.clN.bigWig \
				cats_atacATAC_FullFilter_qval/results/bwa/mergedReplicate/bigwig/D6_100_M.mRp.clN.bigWig \
				cats_atacATAC_FullFilter_qval/results/bwa/mergedReplicate/bigwig/D6_500_3.mRp.clN.bigWig \
                cats_atacpublic_data/Peterson_more/results/bwa/mergedLibrary/bigwig/Foxa2-IP_R1.mLb.clN.bigWig \
                cats_atacpublic_data/Cernilogar2019_PMID31350899/results/bwa/mergedReplicate/bigwig/Foxa2-d3F.mRp.clN.bigWig \
                cats_atacpublic_data/Cernilogar2019_PMID31350899/results/bwa/mergedReplicate/bigwig/Foxa2-d5FS.mRp.clN.bigWig \
			-R cats_atacATAC_deeptools_filteredclusters/1_deeptools_ATAC_ChIP/kmeans_clusters.bed \
            --referencePoint center \
            -a 1500 \
            -b 1500 \
            --missingDataAsZero \
            --sortRegions no \
            --skipZeros \
			-o endochip_on_atac_clusters.mat.gz

plotHeatmap -m endochip_on_atac_clusters.mat.gz -o Heatmap_endoChIP_ATAC_colors.pdf \
            --refPointLabel ATAC_centre \
            --colorMap Greys Blues Oranges Reds Greens Greys Greys Purples Greys RdPu RdPu \
            --zMax 2.0 2.0 2.0 2.0 2.0 0.60 0.15 0.15
