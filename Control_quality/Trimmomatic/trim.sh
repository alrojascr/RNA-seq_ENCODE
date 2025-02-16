#!/bin/bash
#SBATCH --job-name=Trimmomatic                      
#SBATCH --output=trim_%j.log                        
#SBATCH --error=trim_%j.err                           
#SBATCH --nodelist=pujnodo0
#SBATCH --ntasks=1                                    
#SBATCH --cpus-per-task=10                            
#SBATCH --mem=18G                             
#SBATCH --partition=debug                              

# Directories
INPUT_DIR="/opt/data/HPC01A/alexis_rojasc/OSA_RNA-seq/Control_quality/FastQC_Raw"
OUTPUT_DIR="/opt/data/HPC01A/alexis_rojasc/OSA_RNA-seq/Control_quality/Trimmomatic"
ADAPTERS="/opt/data/HPC01A/alexis_rojasc/OSA_RNA-seq/Control_quality/Trimmomatic/TruSeq3.fa"

# Create output directories for trimmed and unpaired reads
TRIMMED_DIR="${OUTPUT_DIR}/Trimmed"
UNPAIRED_DIR="${OUTPUT_DIR}/Unpaired"

mkdir -p "${TRIMMED_DIR}"
mkdir -p "${UNPAIRED_DIR}"

# Start timing the execution
start_time=$(date +%s)

# Loop through all fastq files in the input directory
for file1 in "${INPUT_DIR}"/*_1.fastq; do
  # Extract the sample name from the filename (assuming _1 and _2 are part of the filenames)
  sample_name=$(basename "$file1" _1.fastq)

  # Define the corresponding reverse read file
  file2="${INPUT_DIR}/${sample_name}_2.fastq"

  # Define the output file paths for paired-end data
  output_forward="${TRIMMED_DIR}/${sample_name}_1_trimmed.fastq"
  output_reverse="${TRIMMED_DIR}/${sample_name}_2_trimmed.fastq"
  output_unpaired_forward="${UNPAIRED_DIR}/${sample_name}_1_unpaired.fastq"
  output_unpaired_reverse="${UNPAIRED_DIR}/${sample_name}_2_unpaired.fastq"

  # Run Trimmomatic with paired-end parameters
  trimmomatic PE -phred33 \
    "${file1}" "${file2}" \
    "${output_forward}" "${output_unpaired_forward}" \
    "${output_reverse}" "${output_unpaired_reverse}" \
    ILLUMINACLIP:${ADAPTERS}:2:30:10 \
    LEADING:5 \
    TRAILING:5 \
    SLIDINGWINDOW:4:30 \
    MINLEN:30

done

# End timing the execution
end_time=$(date +%s)

# Calculate the elapsed time
elapsed_time=$((end_time - start_time))

# Print the elapsed time
echo "Execution time: ${elapsed_time} seconds"

