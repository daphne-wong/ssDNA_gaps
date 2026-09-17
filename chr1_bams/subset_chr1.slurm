#!/bin/bash
#SBATCH --job-name=subset_chr1
#SBATCH --partition=main
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --account=tp_612_1653
#SBATCH --mail-user=daphnew@usc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --array=1-4

module load samtools

samples=("UWB_Cisp1" "UWB_Cisp0" "OV3_Cisp1_masked" "OV3_Cisp0")
sample=${samples[$((SLURM_ARRAY_TASK_ID - 1))]}

out_dir="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/chr1_bams"
merged_dir="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/merged_bams"

in_bam="${merged_dir}/${sample}_detect_q30_full.bam"
out_bam="${out_dir}/${sample}_detect_chr1.bam"

samtools view -h -@ ${SLURM_CPUS_PER_TASK} -b ${in_bam} -o ${out_bam} chr1
samtools index -@ "${SLURM_CPUS_PER_TASK}" "${out_bam}"
samtools quickcheck -v "${out_bam}"
