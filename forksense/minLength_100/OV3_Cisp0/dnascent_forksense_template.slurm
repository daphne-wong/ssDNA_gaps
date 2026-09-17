#!/bin/bash
#SBATCH --partition=main
#SBATCH --account=tp_612_1653
#SBATCH --job-name=dnascent_forksense
#SBATCH --output=%A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --array=1-4


# THIS MUST BE RUN IN SEPARATE DIRECTORIES
# BED OUTPUTS DEFAULT TO DIR WHERE COMMAND WAS SUBMITTED, SO RUNNING AN ARRAY FROM ONE DIR WILL CAUSE THEM TO OVERWRITE EACH OTHER

module load apptainer
module load openjdk
module load cuda
module load bzip2

export APPTAINER_BINDPATH=/project2
export APPTAINER_CACHEDIR=/project2/tp_612_1653/resources/apptainer_cache
export NXF_SINGULARITY_CACHEDIR=/project2/tp_612_1653/for_Tirzah/dnascent_petrina/local_singularity_cache

dnascent="/project2/tp_612_1653/software/dnascent-4.2.1/DNAscent.sif"

bam_dir="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/merged_bams"

samplenames="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/samplenames.txt"
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${samplenames})

sample_bam="${bam_dir}/${sample}_detect_q30_full.bam"

out_dir="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/forksense/minLength_100/${sample}"
mkdir -p ${out_dir}

#dnascent="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/dw_dnascent/bin/DNAscent"

#${dnascent} forkSense \

apptainer run --nv ${dnascent} forkSense \
    -d ${sample_bam} \
    -o ${out_dir}/${sample}.forkSense \
    --order EdU,BrdU \
    -t ${SLURM_CPUS_PER_TASK} \
    --markAnalogues \
    --markOrigins \
    --markTerminations \
    --markForks \
    --makeSignatures

