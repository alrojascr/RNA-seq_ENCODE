# RNA-Seq Pipeline Using ENCODE

This pipeline provides a comprehensive framework for analyzing RNA sequencing (RNA-seq) data with tools and standards defined by ENCODE. It covers every step, from data acquisition to quantification and differential expression analysis, incorporating tools designed to ensure high-quality results. 

## Table of Contents
- [Data Acquisition](#data-acquisition)
- [Control Quality](#control-quality)
- [Ribosomal RNA Filtering](#ribosomal-rna-filtering)
- [Mapping](#mapping)
- [Mapping Quality](#mapping-quality)
- [Gene Biotype](#gene-biotype)
- [Quantification](#quantification)
- [MultiQC](#multiqc)
- [How to Run?](#how-to-run)
- [Acknowledgements](#acknowledgements)

## Mapping
**Description:**  
The RNA-seq reads are mapped to a reference genome using a highly efficient aligner. This step produces the alignment of sequencing reads to the reference genome, essential for the next steps in quantification and differential expression analysis.

The reference genomes used in this step are those deposited in the [GENCODE](https://www.gencodegenes.org/) page, which provides high-quality, functional annotations for human and mouse genomes. Please ensure you download the relevant reference genome and its annotations from the GENCODE website.

**Tools:**

- **STAR**  
  - **Description:** STAR is a fast and accurate RNA-seq read aligner. It is capable of aligning reads to large genomes and transcriptomes.  
  - **GitHub Repository:** [STAR GitHub](https://github.com/alexdobin/STAR)  
  - **Command Example:**  
    ```bash
    STAR --genomeDir genome_index --readFilesIn input.fastq --runThreadN 4 --outFileNamePrefix aligned_
    ```

- **RSEM**  
  - **Description:** RSEM is used to perform accurate quantification of transcript abundance from RNA-seq data.  
  - **GitHub Repository:** [RSEM GitHub](https://github.com/deweylab/RSEM)  
  - **Command Example:**  
    ```bash
    rsem-calculate-expression --paired-end aligned_reads.bam reference transcripts
    ```

## Mapping Quality
**Description:**  
This section assesses the quality of the alignments and checks for any potential issues with the mapped data. This includes evaluating read coverage, mapping efficiency, and detecting any issues with duplicates or misaligned reads.

**Tools:**

- **Samtools**  
  - **Description:** A suite of tools for manipulating alignments in the SAM, BAM, and CRAM formats.  
  - **GitHub Repository:** [Samtools GitHub](https://github.com/samtools/samtools)  
  - **Command Example:**  
    ```bash
    samtools flagstat aligned_reads.bam
    ```

- **RSeQC**  
  - **Description:** A quality control toolset for RNA-seq data that provides a wide range of mapping quality checks.  
  - **GitHub Repository:** [RSeQC GitHub](https://github.com/hasherm/seqc)  
  - **Command Example:**  
    ```bash
    python read_distribution.py -i aligned_reads.bam -r reference.gtf
    ```

- **Preseq**  
  - **Description:** A tool for evaluating sequencing depth, assessing library complexity, and predicting the impact of sequencing depth.  
  - **GitHub Repository:** [Preseq GitHub](https://github.com/smithlabcode/preseq)  
  - **Command Example:**  
    ```bash
    preseq c_curve -B input.bam > preseq_output.txt
    ```

- **Picard**  
  - **Description:** A set of command-line tools for working with high-throughput sequencing data, including various mapping quality metrics.  
  - **GitHub Repository:** [Picard GitHub](https://github.com/broadinstitute/picard)  
  - **Command Example:**  
    ```bash
    java -jar picard.jar CollectAlignmentSummaryMetrics I=aligned_reads.bam O=alignment_metrics.txt
    ```

- **Qualimap**  
  - **Description:** A tool for assessing the quality of alignment in RNA-seq data, with a focus on visualization of mapping quality and coverage.  
  - **GitHub Repository:** [Qualimap GitHub](https://github.com/ualib-ros/Qualimap)  
  - **Command Example:**  
    ```bash
    qualimap bamqc -bam aligned_reads.bam -outdir ./qualimap_output
    ```

## Gene Biotype
**Description:**  
This step predicts the types of genes present in the RNA-seq data, helping to categorize genes into various biotypes, such as protein-coding, non-coding, or others. This can provide further insight into gene expression profiles.

**Tools:**

- **FeatureCounts**  
  - **Description:** A tool used to assign reads to genomic features (such as genes, exons, etc.) based on a reference annotation. It also allows for the identification and classification of gene biotypes.  
  - **GitHub Repository:** [FeatureCounts GitHub](https://github.com/subreadteam/subread)  
  - **Command Example:**  
    ```bash
    featureCounts -a annotation.gtf -o gene_biotype_counts.txt aligned_reads.bam
    ```

## Quantification
**Description:**  
In this step, the expression levels of transcripts are quantified. The aligned reads are counted to estimate gene expression levels, which can be used in differential expression analysis.

**Tools:**

- **RSEM**  
  - **Description:** RSEM is used again in this section to calculate gene and transcript expression levels based on the aligned reads.  
  - **GitHub Repository:** [RSEM GitHub](https://github.com/deweylab/RSEM)  
  - **Command Example:**  
    ```bash
    rsem-calculate-expression --paired-end aligned_reads.bam reference transcripts
    ```

## MultiQC
**Description:**  
After the RNA-seq analysis is completed, MultiQC provides a comprehensive summary report. It aggregates quality control metrics from all steps of the pipeline (e.g., FastQC, STAR, RSEM) into a single, easy-to-understand report, providing an overview of the analysis and highlighting potential issues.

**Tools:**

- **MultiQC**  
  - **Description:** A tool for aggregating quality control metrics from various tools into a single HTML report for easy visualization and interpretation.  
  - **GitHub Repository:** [MultiQC GitHub](https://github.com/ewels/MultiQC)  
  - **Command Example:**  
    ```bash
    multiqc ./
    ```

## How to Run?

All scripts for each section of the RNA-seq pipeline are designed to be executed on a high-performance computing cluster using the **SLURM** system. To run any of the scripts, use the following command:

```bash
sbatch <script_name>.sh
```

## Acknowledgements

This pipeline has benefited from the resources and support of the **High Performance Computing (HPC) Cluster at Pontificia Universidad Javeriana**. Special thanks to **Juan Guillermo Torres Hurtado** for his continuous assistance and services that made this computational work possible.

For any questions or issues related to the pipeline, please feel free to open an issue or contact the contributors.
