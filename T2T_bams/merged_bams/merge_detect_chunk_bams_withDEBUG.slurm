#!/bin/bash
#SBATCH --job-name=merge_detect_chunks_DB
#SBATCH --partition=main
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --account=tp_612_1653
#SBATCH --mail-user=daphnew@usc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --array=1-4
#SBATCH --output=logs/merge_withDB_%A_%a.out
#SBATCH --error=logs/merge_withDB_%A_%a.err

module load samtools

samples=("UWB_Cisp1" "UWB_Cisp0" "OV3_Cisp1" "OV3_Cisp0")
sample=${samples[$((SLURM_ARRAY_TASK_ID - 1))]}

detect_dir="/project2/tp_612_1653/for_Tirzah/ssDNA_gaps/T2T_bams/detect_chunks"
merged_dir="/project2/tp_612_1653/for_Tirzah/ssDNA_gaps/T2T_bams/merged_bams"
mkdir -p "${merged_dir}"

unsorted_out="${merged_dir}/${sample}_detect_q20_unsorted.bam"
final_out="${merged_dir}/${sample}_detect_q20_full.bam"

echo "### Starting merge for ${sample}"
echo $(date)

samtools merge -f -@ ${SLURM_CPUS_PER_TASK} "${unsorted_out}" ${detect_dir}/${sample}_chunk*_detect_q20.sorted.bam
if [ $? -ne 0 ]; then
    echo "samtools merge FAILED for ${sample}"
    echo $(date)
    exit 1
fi

samtools sort -@ ${SLURM_CPUS_PER_TASK} -o "${final_out}" "${unsorted_out}"
if [ $? -ne 0 ]; then
    echo "samtools sort FAILED for ${sample}"
    echo $(date)
    exit 1
fi

samtools index "${final_out}"
if [ $? -ne 0 ]; then
    echo "samtools index FAILED for ${sample}"
    echo $(date)
    exit 1
fi

rm -f "${unsorted_out}"

echo "### Merge DONE for ${sample}"
echo $(date)
