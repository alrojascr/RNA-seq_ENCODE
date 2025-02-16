#!/bin/bash
#SBATCH --job-name=FastQC                              # Job name
#SBATCH --output=fastqc_%j.log                         # Standard output log file
#SBATCH --error=fastqc_%j.err                          # Error log file
#SBATCH --nodelist=pujnodo0
#SBATCH --ntasks=1                                     # Number of tasks
#SBATCH --cpus-per-task=4                              # CPUs per task
#SBATCH --mem=12G                                      # Memory allocation
#SBATCH --partition=debug                              # Slurm partition name (debug partition)

# Directories
FASTQ_DIR="/opt/data/HPC01A/alexis_rojasc/OSA_RNA-seq/Control_quality/FastQC_Raw"
OUTPUT_DIR="/opt/data/HPC01A/alexis_rojasc/OSA_RNA-seq/Control_quality/FastQC_Raw/FastQC_Output"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Record the start time of the job
job_start=$(date +"%Y-%m-%d %H:%M:%S")
echo "Job started at: $job_start"

# Process each .fastq file in the directory
for fastq_file in "$FASTQ_DIR"/*.fastq; do
    # Extract the base name of the file (without directory and extension)
    base_name=$(basename "$fastq_file" .fastq)

    # Record the start time for the current file
    file_start=$(date +"%Y-%m-%d %H:%M:%S")
    echo "Processing $fastq_file started at: $file_start"

    # Run FastQC on the current file
    fastqc -o "$OUTPUT_DIR" -t $SLURM_CPUS_PER_TASK "$fastq_file"

    # Record the end time for the current file
    file_end=$(date +"%Y-%m-%d %H:%M:%S")
    echo "Completed: $base_name at $file_end"
done

# Record the end time of the job
job_end=$(date +"%Y-%m-%d %H:%M:%S")
echo "Job completed at: $job_end"
