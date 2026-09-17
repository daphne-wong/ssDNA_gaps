#!/bin/bash
#SBATCH --partition=main
#SBATCH --account=tp_612_1653
#SBATCH --job-name=dnascent_meIODy
#SBATCH --output=logs/meIODy_%A_%a.out
#SBATCH --error=logs/meIODy_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --mail-user=daphnew@usc.edu
#SBATCH --mail-type=ALL

#module load apptainer
#module load openjdk
module load bzip2

#export APPTAINER_BINDPATH=/project2
#export APPTAINER_CACHEDIR=/project2/tp_612_1653/resources/apptainer_cache
#export NXF_SINGULARITY_CACHEDIR=/project2/tp_612_1653/for_Tirzah/dnascent_petrina/local_singularity_cache

#dnascent="/project2/tp_612_1653/software/dnascent-4.2.1/DNAscent.sif"

samplenames="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/samplenames.txt"
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${samplenames})

dnascent="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/dw_dnascent/bin/DNAscent"

forksense_outputs="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/forksense/minLength_100"
detect_bam="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/merged_bams/${sample}_detect_q30_full.bam"

out_file="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/meIODy/${sample}_iod.dnascent"


${dnascent} meIODy \
    -l ${forksense_outputs}/${sample}/leftForks_DNAscent_forkSense.bed \
    -r ${forksense_outputs}/${sample}/rightForks_DNAscent_forkSense.bed \
    --origin ${forksense_outputs}/${sample}/origins_DNAscent_forkSense.bed \
    --termination ${forksense_outputs}/${sample}/terminations_DNAscent_forkSense.bed \
    -d ${detect_bam} \
    --tPulse1 5. \
    --tPulse2 10. \
    -o ${out_file} \
    -t ${SLURM_CPUS_PER_TASK}
