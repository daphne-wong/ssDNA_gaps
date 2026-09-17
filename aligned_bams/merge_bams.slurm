#!/bin/bash
#SBATCH --job-name=merge_bams
#SBATCH --partition=main
#SBATCH --output=logs/merge_bams_%A_%a.out
#SBATCH --error=logs/merge_bams_%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB
#SBATCH --account=tp_612_1653
#SBATCH --mail-user=daphnew@usc.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --array=1-4

module load samtools

samplelist="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/aligned_bams/samplesheet.txt"
sample=$(tail -n +2 ${samplelist} | cut -d, -f1 | sort -u | sed -n "${SLURM_ARRAY_TASK_ID}p")
#sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${samplelist})
barcodes=$(awk -F, -v s="${sample}" '$1==s {print $2}' ${samplelist})

out_dir="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/aligned_bams"

mkdir -p logs
mkdir -p ${out_dir}

# BAM directories
merged_bam="${out_dir}/${sample}.sorted.merged.bam"


echo "[$(date)] Starting merge for sample: ${sample}"
echo "Barcodes: ${barcodes}"
echo "Output:   ${merged_bam}"

am_inputs=()
missing=0
for bc in ${barcodes}; do
    bam="${out_dir}/demux_bams/${bc}.sorted.bam"
    if [[ -f "${bam}" ]]; then
        bam_inputs+=("${bam}")
    else
        echo "[$(date)] Missing BAM: ${bam}"
        missing=1
    fi
done

if [[ ${missing} -eq 0 && ${#bam_inputs[@]} -gt 0 ]]; then
    samtools merge -@ "${SLURM_CPUS_PER_TASK}" "${merged_bam}" "${bam_inputs[@]}"
    samtools index "${merged_bam}"
    echo "[$(date)] Merge complete for ${sample}"

    echo "[$(date)] Running flagstat on merged BAM..."
    samtools flagstat -@ "${SLURM_CPUS_PER_TASK}" "${merged_bam}" | tee "${out_dir}/${sample}.flagstat.txt"
else
    echo "[$(date)] Skipping merge for ${sample} due to missing BAM(s)"
fi


