#!/bin/bash

clusterpaths=(outputs_cats-atac_5_kmeans//Intervals_Cluster_1.bed \
                outputs_cats-atac_5_kmeans//Intervals_Cluster_7.bed \
                outputs_cats-atac_5_kmeans//Intervals_Cluster_2.bed \
                outputs_cats-atac_5_kmeans//Intervals_Cluster_3.bed \
                outputs_cats-atac_5_kmeans//Intervals_Cluster_4.bed \
                outputs_cats-atac_5_kmeans//Intervals_Cluster_5.bed \
                outputs_cats-atac_5_kmeans//Intervals_Cluster_6.bed \
                outputs_cats-atac_5_kmeans//Intervals_Cluster_8.bed \
                outputs_cats-atac_5_kmeans//Intervals_Cluster_9.bed)

foxa2paths=(outputs_ifoxa2-atac_1/Foxa2_up.bed \
        outputs_ifoxa2-atac_1//Foxa2_down.bed )


for cluster in "${clusterpaths[@]}"; do
    cleanCluster=${cluster%%.bed}
    cleanCluster=${cleanCluster##*/Intervals_}
    for foxa2atac in "${foxa2paths[@]}"; do
        cleanfoxa2atac=${foxa2atac%%_peaks.broadPeak}
        cleanfoxa2atac=${cleanfoxa2atac%%.bed}
        cleanfoxa2atac=${cleanfoxa2atac##*/}

        echo "#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --job-name=run_intersects_${cleanfoxa2atac}_${cleanCluster}
#SBATCH --mem=8G
#SBATCH -n 4
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joaquina.delas@crick.ac.uk
#SBATCH --output=run_intersects_${cleanfoxa2atac}_${cleanCluster}.o
#SBATCH --error=run_intersects_${cleanfoxa2atac}_${cleanCluster}.e" > run_intersects_${cleanfoxa2atac}_${cleanCluster}.sh


        echo "bedtools intersect -wa -a $cluster -b $foxa2atac > ${cleanfoxa2atac}__within__${cleanCluster}.bed" >> run_intersects_${cleanfoxa2atac}_${cleanCluster}.sh

        echo "sbatch run_intersects_${cleanfoxa2atac}_${cleanCluster}.sh" >> go.sh
    done
done

echo "Done: run 'go.sh' to submit jobs"
