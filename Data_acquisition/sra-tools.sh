#!/bin/bash
#SBATCH --job-name=sra_download                        # Job name
#SBATCH --output=sra_%j.log                            # Standard output log file, %j is replaced by job ID
#SBATCH --error=sra_%j.err                             # Error log file
#SBATCH --ntasks=1                                     # Number of tasks
#SBATCH --cpus-per-task=4                              # CPUs per task
#SBATCH --mem=8G                                       # Memory allocation
#SBATCH --partition=debug                              # Slurm partition name (debug partition)


# Set the base path for logs and output
LOG_PATH="/opt/data/HPC01A/alexis_rojasc/OSA_RNA-seq/Dataset"

# Start time of the script
START_TIME=$(date +%s)

echo "Job started at: $(date)"

# Create a directory for the downloaded data
OUTPUT_DIR="${LOG_PATH}/SRP449460"
mkdir -p ${OUTPUT_DIR}

# Navigate to the output directory
cd ${OUTPUT_DIR}

# Prefetch: Download the SRA dataset (SRP319371) from NCBI SRA
echo "Starting data download from SRA using prefetch..."
prefetch SRP449460

# Check if prefetch was successful
if [ $? -eq 0 ]; then
    echo "Data download completed successfully with prefetch."
else
    echo "Error during the prefetch download."
    exit 1
fi

# Fastq-dump: Convert downloaded .sra files to FASTQ format
echo "Converting SRA files to FASTQ format using fastq-dump..."
fastq-dump --split-3 ${LOG_PATH}/SRP449460/*/*.sra

# Check if fastq-dump was successful
if [ $? -eq 0 ]; then
    echo "FASTQ conversion completed successfully."
else
    echo "Error during FASTQ conversion."
    exit 1
fi

# End time of the script
END_TIME=$(date +%s)

# Calculate execution time
EXECUTION_TIME=$((END_TIME - START_TIME))

# Output execution time
echo "Job completed at: $(date)"
echo "Total execution time: $((EXECUTION_TIME / 60)) minutes and $((EXECUTION_TIME % 60)) seconds."

# End of script
echo "Process completed."

