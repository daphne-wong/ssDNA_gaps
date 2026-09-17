#!/bin/bash
#SBATCH --partition=main
#SBATCH --account=tp_612_1653
#SBATCH --job-name=gaps_py
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --array=1-4


module load legacy/CentOS7
module load gcc/11.3.0
module load openblas/0.3.20
module load python/3.9.12

source /project2/tp_612_1653/for_Tirzah/dnascent_petrina/venvs/dnascent_flank/bin/activate
python -c "import pysam; print(pysam.__file__); print(pysam.__version__)"
python -c 'import pysam; print(pysam.__file__); print(pysam.__version__); print(hasattr(pysam.AlignedSegment, "modified_bases"))'

bam_dir="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/merged_bams"

samplenames="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/samplenames.txt"
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${samplenames})
sample_bam="${bam_dir}/${sample}_detect_q30_full.bam"

ref_genome="/project2/tp_612_1653/resources/ont/hg38/hg38.fa"

py_script="/project2/tp_612_1653/for_Tirzah/dnascent_petrina/gaps_python/dnascent_gaps.py"

python ${py_script} -i ${sample_bam} -o /project2/tp_612_1653/for_Tirzah/dnascent_petrina/gaps_python/${sample}
