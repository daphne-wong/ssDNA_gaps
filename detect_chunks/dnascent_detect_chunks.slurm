#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --account=tp_612_1653
#SBATCH --job-name=dnascent_detect_bam_chunk_masked
#SBATCH --output=logs/dnascent_detect_chunk_masked_%A_%a.out
#SBATCH --error=logs/dnascent_detect_chunk_masked_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=2-00:00:00
#SBATCH --array=1-4
#SBATCH --mail-user=daphnew@usc.edu
#SBATCH --mail-type=ALL

module load apptainer
module load openjdk
module load cuda
module load samtools

export APPTAINER_BINDPATH=/project2
export APPTAINER_CACHEDIR=/project2/tp_612_1653/resources/apptainer_cache
export NXF_SINGULARITY_CACHEDIR=/project2/tp_612_1653/for_Tirzah/dnascent_petrina/local_singularity_cache

dnascent="/project2/tp_612_1653/software/dnascent-4.2.1/DNAscent.sif"

#bam_dir="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/aligned_bams"
bam_dir="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/chunk_bams"
out_dir="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/detect_chunks"

mkdir -p ${out_dir}

#samplenames="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/samplenames.txt"
samplenames="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/detect_chunks/chunklist_masked.txt"
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${samplenames})
#sample_bam="${bam_dir}/${sample}.sorted.merged.bam"
sample_bam="${bam_dir}/${sample}.sorted.bam"

dnascent_index="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/index.dnascent"
ref_genome="/project2/tp_612_1653/resources/ont/hg38/hg38.fa"


echo "Running DNAscent detect on ${sample}..."


# Run DNAscent via Apptainer/Singularity
apptainer run --nv ${dnascent} detect -b ${sample_bam} -r ${ref_genome} -i ${dnascent_index} -o ${out_dir}/${sample}_detect_q30.bam -q 30 -l 5000 -t 8 --GPU ${CUDA_VISIBLE_DEVICES}


# Run DNAscent via local
# dnascent="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/dw_dnascent/bin/DNAscent"
# ${dnascent} detect -b ${sample_bam} -r ${ref_genome} -i ${dnascent_index} -o ${out_dir}/${sample}_detect_q30.bam -q 30 -l 5000 -t 8 --GPU ${CUDA_VISIBLE_DEVICES}

# Stop if DNAscent failed
if [ $? -ne 0 ]; then
    echo "DNAscent detect failed for ${sample}"
    exit 1
fi

echo "Sorting output BAM..."

samtools sort \
    -@ ${SLURM_CPUS_PER_TASK} \
    -o ${out_dir}/${sample}_detect_q30.sorted.bam \
    ${out_dir}/${sample}_detect_q30.bam

samtools index \
    -@ ${SLURM_CPUS_PER_TASK} \
    ${out_dir}/${sample}_detect_q30.sorted.bam
