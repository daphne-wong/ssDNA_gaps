#!/bin/bash
#SBATCH --job-name=modkit_summary
#SBATCH --partition=main
#SBATCH --account=tp_612_1653
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --array=1-4
#SBATCH --output=logs/modkit_summary_%A_%a.out
#SBATCH --error=logs/modkit_summary_%A_%a.err

set -euo pipefail


#### User settings ####

samplelist="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/samplenames.txt"
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${samplelist}")

bam="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/merged_bams/${sample}_detect_q30_full.bam"
bam_name=$(basename ${bam})

out_dir="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/modkit_summaries"
mkdir -p "${out_dir}"



# Modkit
modkit_bin="/project2/tp_612_1653/resources/ont/modkit/bin/modkit"
probability_threshold=0.5



#### DEBUG ####

if [[ -z "${bam}" ]]; then
    echo "ERROR: No BAM found for array task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

if [[ ! -f "${bam}" ]]; then
    echo "ERROR: BAM does not exist: ${bam}" >&2
    exit 1
fi

if [[ ! -x "${modkit_bin}" ]]; then
    echo "ERROR: modkit executable not found: ${modkit_bin}" >&2
    exit 1
fi



#### Run modkit ####

raw_summary="${out_dir}/${sample}_modkit_summary.txt"
tsv_summary="${out_dir}/${sample}_analogue_counts.tsv"


echo "Input BAM: ${bam}"
echo "Probability threshold: ${probability_threshold}"
echo "Output TSV: ${tsv_summary}"

"${modkit_bin}" modbam summary "${bam}" \
    --filter-threshold "${probability_threshold}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    > "${raw_summary}"


#### Extract EdU and BrdU counts ####

awk -v sample="${sample}" \
    -v bam="${bam_name}" \
    -v threshold="${probability_threshold}" \
    'BEGIN {
        OFS = "\t"
        print "sample", "bam", "analogue", "code", \
              "probability_threshold", "called_bases", \
              "called_fraction", "all_calls", "all_fraction"
    }
    $2 == "b" {
        print sample, bam, "BrdU", $2, threshold, $3, $4, $5, $6
    }
    $2 == "e" {
        print sample, bam, "EdU", $2, threshold, $3, $4, $5, $6
    }' "${raw_summary}" > "${tsv_summary}"

echo "Finished ${sample}"
