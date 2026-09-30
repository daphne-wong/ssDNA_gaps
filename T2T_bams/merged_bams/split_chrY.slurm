#!/bin/bash
#SBATCH --job-name=split_chrY
#SBATCH --account=tp_612_1653
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/split_chrY_%A_%a.out
#SBATCH --error=logs/split_chrY_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=daphnew@usc.edu
#SBATCH --array=1-4

module load samtools

samples=("UWB_Cisp1" "UWB_Cisp0" "OV3_Cisp1" "OV3_Cisp0")
sample=${samples[$((SLURM_ARRAY_TASK_ID - 1))]}

merged_dir="/project2/tp_612_1653/for_Tirzah/ssDNA_gaps/T2T_bams/merged_bams"
chrY_dir="${merged_dir}/chrY"
mkdir -p "${chrY_dir}"

in_bam="${merged_dir}/${sample}_detect_q20_full.bam"
chrY_bam="${chrY_dir}/${sample}_chrY.bam"
tmp_bam="${merged_dir}/${sample}_detect_q20_full.noChrY.tmp.bam"
chrY_name="NC_060948.1"
t=${SLURM_CPUS_PER_TASK}

fail() { echo "$1 FAILED for ${sample}"; date; exit 1; }

echo "### Starting chrY split for ${sample}"
echo "$(date)"

# Make sure input is indexed and actually contains chrY
[ -f "${in_bam}.bai" ] || samtools index "${in_bam}" || fail "samtools index (input)"
samtools idxstats "${in_bam}" | cut -f1 | grep -qx "${chrY_name}" || fail "finding ${chrY_name} in header"

# 1) chrY reads -> own BAM
samtools view -b -@ ${t} -o "${chrY_bam}" "${in_bam}" "${chrY_name}" || fail "samtools view (chrY)"
samtools index "${chrY_bam}" || fail "samtools index (chrY)"

# 2) All other contigs -> temp BAM
keep=$(samtools idxstats "${in_bam}" | cut -f1 | grep -vx -e "${chrY_name}" -e '*')
samtools view -b -@ ${t} -o "${tmp_bam}" "${in_bam}" ${keep} || fail "samtools view (non-chrY)"

# 3) Check read counts add up before overwriting
n_total=$(samtools view -c -@ ${t} "${in_bam}")
n_y=$(samtools view -c -@ ${t} "${chrY_bam}")
n_rest=$(samtools view -c -@ ${t} "${tmp_bam}")
echo "total: ${n_total}  chrY: ${n_y}  non-chrY: ${n_rest}"
[ $((n_y + n_rest)) -eq ${n_total} ] || fail "read count check (original BAM left untouched)"

# 4) Replace merged BAM with non-chrY version and re-index
mv "${tmp_bam}" "${in_bam}" || fail "mv"
samtools index "${in_bam}" || fail "samtools index (non-chrY)"

echo "### chrY split DONE for ${sample}"
echo "$(date)"
