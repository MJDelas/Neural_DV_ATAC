# Running TOBIAS footprinting software

Footprint score were calculated using [TOBIAS](https://github.com/loosolab/TOBIAS), installed following the documentation on their wiki. 

The scripts used were as follows

## TOBIAS ATACorrect

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