# Running TOBIAS footprinting software

Footprint score were calculated using [TOBIAS](https://github.com/loosolab/TOBIAS), installed following the documentation on their wiki. 

The scripts below will reproduced the analysis performed in this paper. 
To use this tool in further work, we are also developing a nextflow pipeline: [nf-tobias](https://github.com/luslab/briscoe-nf-tobias)

## 1. TOBIAS ATACorrect

This script will generate a script per sample and `go.sh` to run all in parallel. 

```#!/bin/bash
GENOME='cats_atac/results/reference_genome/genome.fa'
BED='cats_atac/results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.bed'
BLACKLIST='~/.nextflow/assets/nf-core/atacseq/assets/blacklists/mm10-blacklist.bed'

mkdir ATACorrect_mergedReps

echo "#!/bin/bash" > go_correct.sh

for i in cats_atac/results/bwa/mergedReplicate/*.bam; do

cleanName=${i%%.mRp.clN.sorted.bam}
cleanName=${cleanName##*/}

echo $cleanName


echo "#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --job-name=ACorrect_${cleanName}
#SBATCH --mem=16G
#SBATCH -n 32
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joaquina.delas@crick.ac.uk
#SBATCH --output=ACorr_${cleanName}.o
#SBATCH --error=ACorr_${cleanName}.e" > RunATACorrect_$cleanName.sh;

echo "TOBIAS ATACorrect --bam $i --genome $GENOME --peaks $BED --blacklist $BLACKLIST --outdir ATACorrect_mergedReps --cores 32 " >> RunATACorrect_$cleanName.sh;

echo "sbatch RunATACorrect_$cleanName.sh" >> go_correct.sh;

done

echo "Done: run 'go_correct.sh' to submit jobs"
```

## 2. TOBIAS Footprint

This will generate a script per sample to run ATACFootprint on the corrected files from the previous step

```
#!/bin/bash
BED='cats_atac/results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.bed'

mkdir ATACFootprint_mergedReps

echo "#!/bin/bash" > go_footprint.sh

for i in ATACorrect_mergedReps/*_corrected.bw; do

cleanName=${i%%.mRp.clN.sorted_corrected.bw}
cleanName=${cleanName##*/}

echo $cleanName

echo "#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --job-name=AFootprint_${cleanName}
#SBATCH --mem=16G
#SBATCH -n 32
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joaquina.delas@crick.ac.uk
#SBATCH --output=ACorr_${cleanName}.o
#SBATCH --error=ACorr_${cleanName}.e" > RunFootprint_$cleanName.sh;
echo "TOBIAS FootprintScores --signal $i \
	--regions $BED  \
	--output ATACFootprint_mergedReps/${cleanName}_footprints.bw --cores 32" >> RunFootprint_$cleanName.sh;

echo "sbatch RunFootprint_$cleanName.sh" >> go_footprint.sh;

done

echo "Done: run 'go_footprint.sh' to submit jobs"
```

## 3. TOBIAS BINDetect

Script to run BINDetect on all samples. It requires very high memory. 

We used motifs from the following four databases concatenated into `motifs_archetypes.meme`
- `HOCOMOCOv11_core_HUMAN_mono_meme_format.meme`
- `HOCOMOCOv11_core_MOUSE_mono_meme_format.meme`
- `JASPAR2018_CORE_vertebrates_non-redundant_pfms.meme`
- `jolma2013.meme`

```
#!/bin/bash
#SBATCH --partition=hmem
#SBATCH --job-name=bindectect
#SBATCH --mem=1504G
#SBATCH -n 32
#SBATCH --time=72:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joaquina.delas@crick.ac.uk
#SBATCH --output=bindectect.o
#SBATCH --error=bindectect.e


GENOME='cats_atacl/results/reference_genome/genome.fa'
BED='cats_atac/results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.bed'

DATABASE='motifs_archetypes.meme'

mkdir BINDetect_25conditions_arch

QT_QPA_PLATFORM=offscreen TOBIAS BINDetect --motifs $DATABASE \
 	--signal ATACFootprint_mergedReps/D6_100_3_footprints.bw \
	ATACFootprint_mergedReps/D6_500_3_footprints.bw  \
	ATACFootprint_mergedReps/D6_0_1_footprints.bw  \
	ATACFootprint_mergedReps/D6_500_M_footprints.bw  \
	ATACFootprint_mergedReps/D4_100_3_footprints.bw  \
	ATACFootprint_mergedReps/D5_500_M_footprints.bw  \
	ATACFootprint_mergedReps/D6_100_M_footprints.bw  \
	ATACFootprint_mergedReps/D5_500_3_footprints.bw  \
	ATACFootprint_mergedReps/D4_500_3_footprints.bw  \
	ATACFootprint_mergedReps/D4_100_M_footprints.bw  \
	ATACFootprint_mergedReps/D4_10_1_footprints.bw  \
	ATACFootprint_mergedReps/D5_10_2_footprints.bw  \
	ATACFootprint_mergedReps/D4_500_M_footprints.bw  \
	ATACFootprint_mergedReps/D3_0_NMP_footprints.bw  \
	ATACFootprint_mergedReps/D6_100_2_footprints.bw  \
	ATACFootprint_mergedReps/D5_10_1_footprints.bw  \
	ATACFootprint_mergedReps/D5_100_3_footprints.bw  \
	ATACFootprint_mergedReps/D4_100_2_footprints.bw  \
	ATACFootprint_mergedReps/D6_10_2_footprints.bw  \
	ATACFootprint_mergedReps/D4_0_1_footprints.bw  \
	ATACFootprint_mergedReps/D6_10_1_footprints.bw  \
	ATACFootprint_mergedReps/D4_10_2_footprints.bw  \
	ATACFootprint_mergedReps/D5_100_M_footprints.bw  \
	ATACFootprint_mergedReps/D5_0_1_footprints.bw  \
	ATACFootprint_mergedReps/D5_100_2_footprints.bw \
	--genome $GENOME \
	--peaks  $BED \
	--outdir BINDetect_25conditions_arch \
	--cond_names D6_100_3 D6_500_3 D6_0_1 D6_500_M D4_100_3 D5_500_M D6_100_M D5_500_3 D4_500_3 D4_100_M D4_10_1 D5_10_2 D4_500_M D3_0_NMP D6_100_2 D5_10_1 D5_100_3 D4_100_2 D6_10_2 D4_0_1 D6_10_1 D4_10_2 D5_100_M D5_0_1 D5_100_2 \
	--cores 32
```