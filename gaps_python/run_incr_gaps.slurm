#!/bin/bash
#SBATCH --partition=main
#SBATCH --account=tp_612_1653
#SBATCH --job-name=gaps_incr_py
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --array=1-4

# Load modules
module load legacy/CentOS7
module load gcc/11.3.0
module load openblas/0.3.20
module load python/3.9.12

# Activate python venv
source /project2/tp_612_1653/for_Tirzah/ssDNA_gaps/venvs/dnascent_flank/bin/activate

python -c "import pysam; print(pysam.__file__); print(pysam.__version__)"
python -c 'import pysam; print(pysam.__file__); print(pysam.__version__); print(hasattr(pysam.AlignedSegment, "modified_bases"))'


# Set directory paths
bam_dir="/project2/tp_612_1653/for_Tirzah/ssDNA_gaps/T2T_bams/merged_bams"
out_dir="/project2/tp_612_1653/for_Tirzah/ssDNA_gaps/gaps_python/20261001_T2T_gaps"

mkdir -p ${out_dir}


# Set sample variables
samplenames="/project2/tp_612_1653/for_Tirzah/ssDNA_gaps/samplenames.txt"
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${samplenames})
sample_bam="${bam_dir}/${sample}_detect_q20_full.bam"


# Start gap calling
echo "Starting gaps_incremental analysis for ${sample}..."
echo "Using sample bam file: ${sample_bam}"

ref_genome="/project2/tp_612_1653/resources/ont/T2T-CHM13v2/T2T-CHM13v2.fa"
py_script="/project2/tp_612_1653/for_Tirzah/ssDNA_gaps/gaps_python/gaps_incremental_annotated.py"


python ${py_script} -i ${sample_bam} -o /project2/tp_612_1653/for_Tirzah/ssDNA_gaps/gaps_python/${sample} --event-types EdU-BrdU-EdU
