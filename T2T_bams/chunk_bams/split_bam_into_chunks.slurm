#!/bin/bash
#SBATCH --job-name=split_bam_into_chunks
#SBATCH --partition=main
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --account=tp_612_1653
#SBATCH --array=1-N
#SBATCH --output=logs/split_%A_%a.out
#SBATCH --mail-user=daphnew@usc.edu
#SBATCH --mail-type=ALL

set -euo pipefail
module load samtools

samplenames="/project2/tp_612_1653/for_Tirzah/ssDNA_gaps/samplenames.txt"
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${samplenames})
[[ -z "${sample}" ]] && { echo "No sample for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }

bam_dir="/project2/tp_612_1653/for_Tirzah/ssDNA_gaps/T2T_bams/aligned_T2T"
out_dir="/project2/tp_612_1653/for_Tirzah/ssDNA_gaps/T2T_bams/chunk_bams"
mkdir -p ${out_dir}
in_bam=${bam_dir}/${sample}.sorted.merged.bam

# RefSeq T2T-CHM13v2.0 names; chunks = chr1-3 | chr4-7 | chr8-13 | chr14-22,X,Y
declare -a CHUNKS=(
"NC_060925.1 NC_060926.1 NC_060927.1"
"NC_060928.1 NC_060929.1 NC_060930.1 NC_060931.1"
"NC_060932.1 NC_060933.1 NC_060934.1 NC_060935.1 NC_060936.1 NC_060937.1"
"NC_060938.1 NC_060939.1 NC_060940.1 NC_060941.1 NC_060942.1 NC_060943.1 NC_060944.1 NC_060945.1 NC_060946.1 NC_060947.1 NC_060948.1"
)

for i in "${!CHUNKS[@]}"; do
    out_bam=${out_dir}/${sample}_chunk${i}.sorted.bam
    samtools view -@ ${SLURM_CPUS_PER_TASK} -b -o ${out_bam} ${in_bam} ${CHUNKS[$i]}
    samtools index -@ ${SLURM_CPUS_PER_TASK} ${out_bam}
    n=$(samtools view -c ${out_bam})
    echo "${sample} chunk${i}: ${n} records"
    [[ ${n} -eq 0 ]] && echo "WARNING: chunk${i} is empty - check chromosome names"
done
