#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --job-name=deepmatrix_chip
#SBATCH --mem=16G
#SBATCH -n 4
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joaquina.delas@crick.ac.uk
#SBATCH --output=deepmatrix_chip.o
#SBATCH --error=deepmatrix_chip.e

computeMatrix reference-point -S cats_atac/results/bwa/mergedReplicate/bigwig/D3_0_NMP.mRp.clN.bigWig \
                cats_atac/results/bwa/mergedReplicate/bigwig/D4_0_1.mRp.clN.bigWig \
				cats_atac/results/bwa/mergedReplicate/bigwig/D4_10_2.mRp.clN.bigWig \
				cats_atac/results/bwa/mergedReplicate/bigwig/D4_100_M.mRp.clN.bigWig \
				cats_atac/results/bwa/mergedReplicate/bigwig/D4_500_3.mRp.clN.bigWig \
				cats_atac/results/bwa/mergedReplicate/bigwig/D5_0_1.mRp.clN.bigWig \
				cats_atac/results/bwa/mergedReplicate/bigwig/D5_10_2.mRp.clN.bigWig \
				cats_atac/results/bwa/mergedReplicate/bigwig/D5_100_M.mRp.clN.bigWig \
				cats_atac/results/bwa/mergedReplicate/bigwig/D5_500_3.mRp.clN.bigWig \
                cats_atac/results/bwa/mergedReplicate/bigwig/D6_10_1.mRp.clN.bigWig \
				cats_atac/results/bwa/mergedReplicate/bigwig/D6_10_2.mRp.clN.bigWig \
				cats_atac/results/bwa/mergedReplicate/bigwig/D6_100_M.mRp.clN.bigWig \
				cats_atac/results/bwa/mergedReplicate/bigwig/D6_500_3.mRp.clN.bigWig \
                public_data/Nishi2015_PMID26293298/results/bwa/mergedLibrary/bigwig/Nkx61_IP_R1.mLb.clN.bigWig \
                public_data/Nishi2015_PMID26293298/results/bwa/mergedLibrary/bigwig/Olig2_IP_R1.mLb.clN.bigWig \
                public_data/Nishi2015_PMID26293298/results/bwa/mergedLibrary/bigwig/Nkx22_IP_R2.mLb.clN.bigWig \
                sox2_bcat/results/bwa/mergedLibrary/bigwig/WT_SOX2_IP_R1.mLb.clN.bigWig \
                public_data/Peterson_more/results/bwa/mergedLibrary/bigwig/SOX2_SAG_72h_R1.mLb.clN.bigWig \
            -R cats_atac/beds_diffaccess_clusters/Intervals_Cluster_1.bed \
                cats_atac/beds_diffaccess_clusters/Intervals_Cluster_3.bed \
                cats_atac/beds_diffaccess_clusters/Intervals_Cluster_2.bed \
                cats_atac/beds_diffaccess_clusters/Intervals_Cluster_9.bed \
                cats_atac/beds_diffaccess_clusters/Intervals_Cluster_6.bed \
                cats_atac/beds_diffaccess_clusters/Intervals_Cluster_7.bed \
                cats_atac/beds_diffaccess_clusters/Intervals_Cluster_8.bed \
                cats_atac/beds_diffaccess_clusters/Intervals_Cluster_4.bed \
                cats_atac/beds_diffaccess_clusters/Intervals_Cluster_5.bed \
			--outFileSortedRegions kmeans_clusters.bed \
            --referencePoint center \
            -a 1500 \
            -b 1500 \
            --missingDataAsZero \
            --sortRegions no \
            --skipZeros \
			-o chip_on_atac_clusters.mat.gz

plotHeatmap -m chip_on_atac_clusters.mat.gz -o Heatmap_ChIP_ATAC_colors.pdf \
            --refPointLabel ATAC_centre \
            --colorMap Greys Blues Oranges Reds Greens Blues Oranges Reds Greens Blues Oranges Reds Greens Oranges Reds Greens Greys Greys \
            --zMax 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 1.5 1.5 1.5 0.9 1.0