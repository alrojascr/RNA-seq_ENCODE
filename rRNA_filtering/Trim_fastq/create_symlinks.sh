# Source directory containing the original .fastq files
source_dir="/opt/data/HPC01A/alexis_rojasc/OSA_RNA-seq/Control_quality/Trimmomatic/Trimmed"

# Destination directory for the symbolic links
target_dir="/opt/data/HPC01A/alexis_rojasc/OSA_RNA-seq/RNA_ribosomal/Trim_fastq"

# List all .fastq files in the source directory

files=("$source_dir"/*.fastq)

# Create symbolic links with original names
for file in "${files[@]}"; do
    # Extract the filename (without the directory)
    filename=$(basename "$file")
    modified_filename=${filename/_trimmed.fastq/.fastq}
    
    # Create symbolic link
    ln -s "$file" "${target_dir}/${modified_filename}"
    echo "Symbolic link created: $file -> ${target_dir}/${modified_filename}"
done
