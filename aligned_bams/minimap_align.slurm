#!/bin/bash
#SBATCH --job-name=petrina_align_array
#SBATCH --partition=main
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --account=tp_612_1653
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daphnew@usc.edu
#SBATCH --array=11-18
#SBATCH --output=logs/align_%a.out
#SBATCH --error=logs/align_%a.err

module load minimap2
module load samtools

FASTQ_BASE="/project2/tp_612_1653/for_Tirzah/Petrina_Kan_20260714/no_sample_id/20260714_1720_P2I-00347-A_PBM43482_cde9b208/fastq_pass"
REF="/project2/tp_612_1653/resources/ont/hg38/hg38.fa"
OUTDIR="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/aligned_bams"

mkdir -p ${OUTDIR}
mkdir -p logs

BC=$(printf "barcode%02d" ${SLURM_ARRAY_TASK_ID})

echo "Aligning ${BC}..."

minimap2 -ax map-ont -t 16 ${REF} ${FASTQ_BASE}/${BC}/*.fastq.gz \
    | samtools sort -@ 16 -o ${OUTDIR}/${BC}.sorted.bam -
samtools index ${OUTDIR}/${BC}.sorted.bam


echo "### Alignment DONE for ${BC}"
