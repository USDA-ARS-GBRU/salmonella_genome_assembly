#!/bin/bash

#SBATCH -N 1 -n 48
#SBATCH --partition=atlas
#SBATCH --mail-type=END,FAIL
#SBATCH --account=gbru_fy21_salm_compgen
#SBATCH --job-name=run_salm_genome
#SBATCH --output=run_salm_genome_%j.out
#SBATCH --time=7-0

# Load all needed modules
module load openjdk
module load apptainer
module load nextflow

# Directory with Nextflow files
cd /project/gbru_fy21_salm_compgen/annette/salmonella_genome_assembly

date;hostname;pwd
start_time="$(date -u +%s)"

echo "starting nextflow job at $start_time"
report="out_${SLURM_JOB_ID}.report.html"
trace="out_${SLURM_JOB_ID}.trace.txt"

nextflow run main.nf -with-apptainer -resume -with-report $report -with-trace $trace

end_time="$(date -u +%s)"
elapsed="$(($end_time-$start_time))"
echo "Total of $elapsed seconds elapsed for nextflow pipeline"
