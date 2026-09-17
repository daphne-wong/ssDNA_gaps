#!/bin/bash
#SBATCH --job-name=split_bam_into_chunks
#SBATCH --partition=main
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --account=tp_612_1653
#SBATCH --mail-user=daphnew@usc.edu
#SBATCH --mail-type=ALL

module load samtools

#samplenames="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/samplenames.txt"
#sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${samplenames})
#sample_bam="${bam_dir}/${sample}.sorted.merged.bam"

bam_dir="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/aligned_bams"
out_dir="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/chunk_bams"
mkdir -p ${out_dir}

# Group chromosomes into 4 roughly balanced chunks by size
declare -a CHUNKS=(
  "chr1 chr2 chr3"
  "chr4 chr5 chr6 chr7"
  "chr8 chr9 chr10 chr11 chr12 chr13"
  "chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
)



for i in "${!CHUNKS[@]}"; do
    samtools view -b \
        ${bam_dir}/OV3_Cisp1_masked.bam \
        ${CHUNKS[$i]} \
    | samtools sort \
        -@ ${SLURM_CPUS_PER_TASK} \
        -o ${out_dir}/OV3_Cisp1_masked_chunk${i}.sorted.bam -

    samtools index \
        -@ ${SLURM_CPUS_PER_TASK} \
        ${out_dir}/OV3_Cisp1_masked_chunk${i}.sorted.bam
done
