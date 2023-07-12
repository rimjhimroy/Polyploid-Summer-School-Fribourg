# Polyploid Genotyping Protocol
Created by: Rimjhim Roy Choudhury, Armel Salmon, and Christian Parisod  
Date: 2023-07-18

## Table of Contents
- [Introduction](#introduction)
- [Overview of the Protocol](#overview-of-the-protocol)
- [Requirements](#requirements)
- [Data access and Resource allocation](#data-access-and-resource-allocation)
- [How to submit jobs on the IBU cluster](#how-to-submit-jobs-on-the-ibu-cluster)
- [Step 1: Setting up the Conda Environment](#step-1-setting-up-the-conda-environment)
- [Step 2: Preparing the Reference Genome](#step-2-preparing-the-reference-genome)
- [Step 3: Preparing the Truth Set](#step-3-preparing-the-truth-set)
- [Step 4: Genotype Calling with GATK](#step-4-genotype-calling-with-gatk)
- [Step 5: Genotype Calling with FreeBayes](#step-5-genotype-calling-with-freebayes)
- [Step 6: Genotype Calling with Updog (using unfiltered GATK called sites as input)](#step-6-genotype-calling-with-updog-using-unfiltered-gatk-called-sites-as-input)
- [Step 7: Comparing to the Truth Set using vcfeval](#step-7-comparing-to-the-truth-set-using-vcfeval)
- [Conclusion](#conclusion)
- [References](#references)

## Introduction
Genotype calling is the process of determining the genetic variants present in a sample. The most common type of genetic variants are single nucleotide polymorphisms (SNPs), which are single base pair differences between the sample and a reference genome. Genotype calling is a crucial step in many genomic analyses, including genome-wide association studies (GWAS), phylogenetic analyses, and population genetics.

Polyploid genotyping poses several challenges due to the inherent complexity of polyploid genomes.
- **Allelic complexity:** Polyploid genomes contain multiple copies of each chromosome, resulting in increased allelic complexity. Each locus can have multiple alleles, making it difficult to accurately determine the genotype for each allele. Resolving and distinguishing between these multiple alleles is a major challenge in polyploid genotyping.
- **Allele dosage ambiguity:** In polyploids, determining the dosage of each allele is challenging. Differentiating between homozygous and heterozygous genotypes becomes complex when multiple copies of each allele are present. Determining the true dosage for each allele at a given locus is a critical task in polyploid genotyping.
- **Read alignment complexities:** Aligning short sequencing reads to a polyploid reference genome can be challenging due to the presence of multiple similar sequences. The ambiguity in read alignments, especially when dealing with highly similar sequences, can lead to incorrect genotype calls. Furthermore, sequencing read cannot distinguish identical allele copies in the absence of physical linkage with other heterozygous alleles. 
- **Lower per-allele coverage:** Determining allele-specific copy number becomes progressively challenging as the ploidy increases and the per-allele coverage decreases. The reduced coverage per allele also introduces additional uncertainty in distinguishing true genetic variations from sequencing errors.  

## Overview of the Protocol
In this protocol, we will perform genotype calling on polyploid genomes using GATK, FreeBayes, and Updog. 

For this session of practicals, we have generated polyploid datasets by using two reference sequences: `At_Chr4_140-150.fas` and `At_Chr4_35-45.fas` (also known as the knob region).

Next, we created randomly generated 4x sequences with 0.3% polymorphic sites based on the reference sequences. These sequences were stored in the files `At_Chr4_140-150_4x.fas` and `At_Chr4_35-45_4x.fas`. We recorded the SNP positions within these 4x sequences in the files `At_Chr4_140-150_4x_snp_pos.txt` and `At_Chr4_35-45_4x_snp_pos.txt`.
<!-- Armel, can you give more details on how you generated the 4x sequences? -->

Additionally, to simulate realistic sequencing data, we generated insilico Illumina NovaSeq paired reads at different coverage levels. The coverage levels used were 150x, 100x, 50x, 20x, and 5x. These reads represent the experimental data that will be used in further analyses.
<!-- Armel, can you give more details on how you simulated the reads? -->

To align the generated reads to the reference sequences, we employed the Bowtie2 software with specific parameters (`--sensitive -N 1 --no-contain --no-overlap --no-discordant`). The resulting mapped reads were stored in BAM files, which serve as the final files required for this session.

This summer school protocol introduces three widely used tools: GATK (Genome Analysis Toolkit), FreeBayes, and UpDog, which provide powerful solutions for variant calling and genotyping. 
 
- GATK (McKenna et al., 2010), developed by the Broad Institute, offers a comprehensive suite of algorithms and best-practice pipelines for accurate variant discovery  (DePristo et al., 2011; Van der Auwera & O'Connor, 2020). It employs a sophisticated local de novo assembly approach known as the HaplotypeCaller (Poplin et al., 2017), reconstructing the most likely haplotypes in the target regions and calling variants based on these assembled haplotypes.
   - Key Features:
     - Comprehensive variant calling and genotyping.
     - Variant quality score recalibration to filter out false positives.
     - Utilizes the Genome Analysis Toolkit's best-practice pipelines for data analysis.
 
- FreeBayes (Garrison & Marth, 2012) uses Bayesian inference to assess the probability of a variant at each genomic position. It is a haplotype-based variant caller that detects variants by analyzing the actual sequences of aligned reads, rather than relying solely on their alignment. FreeBayes can handle various types of variants and is especially effective in detecting heterozygous and low-frequency variants.
   - Key Features:
     - Bayesian framework for accurate variant calling.
     - Handles various types of genetic variants.
     - Suitable for detecting low-frequency and heterozygous variants.
 
- UpDog (Gerard et. al., 2021) is a novel tool that is designed for genotyping polyploids. It addresses various challenges commonly encountered in NGS data analysis, such as allele bias, overdispersion, and sequencing errors. The name "updog" stands for "Using Parental Data for Offspring Genotyping" and reflects the initial focus of the method on full-sib populations. However, Updog has evolved to be applicable to a wider range of populations, making it a versatile tool for polyploid genotyping from NGS data.
   - Key Features:
     - Specific focus on variant calling in polyploid genomes.
     - Probabilistic model for accurate variant identification.
     - Robust handling of different ploidy levels.
<div style="page-break-after: always;"></div>

## Requirements
- Conda environment created for the summer school
- Polyploid reference genome fasta files
- Simulated reads at different coverage levels: 5x, 20x, 50x, 100x, and 150x and their corresponding BAM files mapped to the reference genome
- Truth set of known variants

## Data access and Resource allocation
Please note the following information regarding data access and resource allocation on the cluster:

1. Data Access:
   - The data for the Summer school practical sessions are located in the directory: /data/courses/pploidy/
   - Students are in the Unix group `pploidystudents` have read-only permissions for these data.

2. User Quotas:
   - Each student and user has a 4TB quota to write into the directory /data/users/loginname.
   - Each student and user has a 20GB quota to write into the directory /home/loginname.
   - Use the `lsquota` command to view the individual quota usage.

3. Cluster Resources:
   - The Slurm partition for the cluster is named `pploidy`.
   - The `pploidy` partition has 128 CPUs reserved.
   - The `pploidy` partition has 512 GB of RAM.


## How to submit jobs on the IBU cluster
You can create a job script to submit jobs on the cluster. The job script contains the commands that you would like to run on the cluster. You can submit the job script using the `sbatch` command. The job will be added to the queue and will be executed when the resources are available. You can check the status of the job using the `squeue` or `sacct` commands. You can cancel the job using the `scancel` command.

1. Example on how to create a job submission script (e.g. `submit.sh`) with the following content:
   ```
   #!/bin/bash
   #SBATCH --job-name=polyploid_genotyping
   #SBATCH --output=polyploid_genotyping_%j.out
   #SBATCH --error=polyploid_genotyping_%j.err
   #SBATCH --time=00:30:00
   #SBATCH --mem=10G
   #SBATCH --cpus-per-task=1
   #SBATCH --partition=pploidy

   # Load the miniconda module
   module load Conda/miniconda/latest

   # Activate the Conda environment created for the summer school
   conda activate -p /mnt/Users/choudhuryr/Polyploid_summer_school_2023/envs/polyploid_env

   # Run the command
   gatk MarkDuplicates -I input.bam -O deduplicated.bam -M metrics.txt
   ```
2. Submit the job using the `sbatch` command:
   ```
   sbatch submit.sh
   ```
3. Check the status of the job using the `squeue` command:
   ```
   squeue -u <your-username>
   ```
4. Check the status of the job using the `sacct` command:

   ```
   sacct -j <job-id>
   ```
5. Cancel the job using the `scancel` command:
   ```
   scancel <job-id>
   ```
Alternatively, you can run the command directly uing the --wrap option:
```
sbatch -J polyploid_genotyping -o polyploid_genotyping_%j.out -e polyploid_genotyping_%j.err --time=00:30:00 --mem=10G -c 1 --wrap="gatk MarkDuplicates -I input.bam -O deduplicated.bam -M metrics.txt"
```
<div style="page-break-after: always;"></div>

## Step 1: Setting up the Conda Environment
1. Load the miniconda module
   ```
   module load Conda/miniconda/latestS
   ```

2. Activate the Conda environment created for the summer school
   ```
   conda activate -p /mnt/Users/choudhuryr/Polyploid_summer_school_2023/envs/polyploid_env
   ```
  
## Step 2: Preparing the Reference Genome
1. Index the reference genome using samtools:
   ```
   samtools faidx reference.fasta
   ```
2. Create a sequence dictionary for the reference genome using Picard:
   ```
   gatk CreateSequenceDictionary -R reference.fasta -O reference.dict
   ```
  
## Step 3: Preparing the Truth Set
1. Convert the truth set of known variants to VCF format using the script provided in the scripts directory:
   ```
   Rscript --vanilla scripts/convert_truthset.R --snp_pos_file truthset.txt --ref reference.fasta --out truthset.vcf
   ```
  
## Step 4: Genotype Calling with GATK

1. Mark duplicates in the aligned BAM files using GATK:
   ```
   gatk MarkDuplicates -I input.bam -O deduplicated.bam -M metrics.txt
   ```

2. Add read groups to the BAM files using GATK:
   ```
   gatk AddOrReplaceReadGroups -I deduplicated.bam -O dedup_rg.bam -LB samplename -PL illumina -PU samplename -SM samplename
   ```

3. Perform variant calling using GATK HaplotypeCaller:
   ```
   gatk HaplotypeCaller -R reference.fasta -I dedup_rg.bam -O raw_variants.vcf -ploidy 4 --native-pair-hmm-threads 1 -stand-call-conf 10
   ```

4. Apply filtering on the raw variant calls using GATK VariantFiltration, adjusting the parameters as needed:
   ```
   gatk VariantFiltration -R reference.fasta -V raw_variants.vcf -O filtered_variants.vcf \
   -filter 'QD < 2.0' --filter-name 'QD2' \
   -filter 'QUAL < 30' --filter-name 'QUAL30' \
   -filter 'GQ < 5' --filter-name 'GQ5' \
   -filter 'FS > 60.0' --filter-name 'FS60' \
   -filter 'SOR > 3.0' --filter-name 'SOR3' \
   -filter 'MQ < 40.0' --filter-name 'MQ40' \
   -filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' \
   -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8'
   ```
  
## Step 5: Genotype Calling with FreeBayes

1. Perform variant calling using FreeBayes:
   ```
   freebayes -f reference.fasta -p 4 --genotype-qualities dedup_rg.bam > raw_variants.vcf
   ```

2. Apply filtering on the raw variant calls using BCFtools, adjusting the parameters as needed:
   ```
   bcftools filter -i 'QUAL > 1 && GQ > 1 && SAF > 0 & SAR > 0 && RPR > 1 && RPL > 1' -s FAIL -o filtered_variants.vcf raw_variants.vcf
   ```
  
## Step 6: Genotype Calling with Updog (using unfiltered GATK called sites as input)
1. Use the variant calls generated by GATK HaplotypeCaller as input for Updog.

2. Use the script in the scripts folder to genotype and filter the variants using Updog:
   ```
   Rscript --vanilla scripts/updog.R -i gatk_variants.vcf -s sample_coverge -o -p 4 -n 10 -o updog_genotyped.vcf
   ```
  
## Step 7: Comparing to the Truth Set using vcfeval

1. Obtain a truth set of known variants in VCF format.

2. Compare the variant calls from each method to the truth set using rtg vcfeval:
   ```
   rtg vcfeval -b truth_set.vcf.gz -c method_calls.vcf.gz -o output_dir -t reference_SDF --sample-ploidy 4  --ref-overlap
   ```
3. Compare the variant calls from each method to the truth set using rtg vcfeval at the base call level, without considering the genotypes
   ```
   rtg vcfeval -b truth_set.vcf.gz -c method_calls.vcf.gz -o output_dir -t reference_SDF --sample-ploidy 4  --ref-overlap --squash-ploidy --output-mode "annotate" --XXcom.rtg.vcf.eval.flag-alternates=true
   ```
   - **TP**: True Positives - The number of variant calls that are correctly identified as true positives compared to the truth set.
   - **FP**: False Positives - The number of variant calls that are incorrectly identified as positive (i.e., present) when they are not present in the truth set.
   - **FN**: False Negatives - The number of variant calls that are present in the truth set but are not identified as positive (i.e., missed) by the variant calling methods.
   - **Precision**: The proportion of true positive calls out of all the positive calls made by the variant calling methods (TP / (TP + FP)).
   - **Sensitivity**: The proportion of true positive calls out of all the variants present in the truth set, representing the ability of the method to identify true positive variants. (TP / (TP + FN)).
   - **F-measure**: The harmonic mean of precision and recall, providing a single metric that balances the trade-off between precision and recall. (2 * ((precision * recall) / (precision + recall)))
   - **Specificity**: The ability of the method to correctly identify true negative calls, calculated as 1 - false positive rate. (TN / (TN + FP))
   - **Accuracy**: The overall correctness of the variant calls, calculated as (TP + TN) / (TP + FP + TN + FN).

4. Analyze and interpret the comparison results to evaluate the performance of each method. Consider the precision, sensitivity, F-measure, specificity, and accuracy values to assess the accuracy and reliability of the variant calling methods.

5. Generate the ROC curve using rtg rocplot:
   ```
   rtg rocplot -o output_dir/rocplot.pdf output_dir/tp-baseline.vcf.gz output_dir/tp-call.vcf.gz
   ```
   - **ROC curve**: A graphical plot of the true positive rate (sensitivity) against the false positive rate (1 - specificity) at various threshold settings. The ROC curve provides a visual representation of the trade-off between sensitivity and specificity for each method. The closer the curve is to the top left corner, the better the performance of the method.


## Conclusion
This protocol provides a step-by-step guide for genotyping polyploid genomes using GATK, FreeBayes, and Updog. By following these steps, you will be able to compare the variant calls to a truth set and assess the performance of each method. Adjust the parameters and filtering criteria based on your specific needs and data characteristics. 
 
We acknowledge the significant contribution of Cooke et al. (2022) in their paper titled "Benchmarking small-variant genotyping in polyploids." This study presents a comprehensive evaluation and comparison of various genotyping methods tailored specifically for polyploid organisms. By systematically analyzing the strengths and limitations of each method, Cooke et al. provide valuable insights into their performance in accurately determining genotypes within polyploid populations. 
 
Enjoy your journey into polyploid genomics!

## References
Cooke, D.P., Wedge, D.C., Lunter, G. Benchmarking small-variant genotyping in polyploids (2022). *Genome research*, 32(2):403-408. 

McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., ... & DePristo, M. A. (2010). The genome analysis toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome research*, 20(9), 1297-1303.

DePristo, M. A., Banks, E., Poplin, R., Garimella, K. V., Maguire, J. R., Hartl, C., ... & Daly, M. J. (2011). A framework for variation discovery and genotyping using next-generation DNA sequencing data. *Nature genetics*, 43(5), 491-498.

Van der Auwera GA, Carneiro MO, Hartl C, Poplin R, Del Angel G, Levy-Moonshine A, Jordan T, Shakir K, Roazen D, Thibault J, Banks E, Garimella KV, Altshuler D, Gabriel S, DePristo MA. From FastQ data to high confidence variant calls: the Genome Analysis Toolkit best practices pipeline. *Current protocols in bioinformatics*. 2013;43(1110):11.10.1-11.10.33.

Poplin R, Ruano-Rubio V, DePristo MA, Fennell TJ, Carneiro MO, Auwera GA, Kling DE, Gauthier LD, Levy-Moonshine A, Roazen D, Shakir K, Thibault J, Chandran S, Whelan CW, Lek M, Gabriel S, Daly MJ, Neale BM, MacArthur DG, Banks E. Scaling accurate genetic variant discovery to tens of thousands of samples. *bioRxiv*. doi: https://doi.org/10.1101/201178

Garrison, E., & Marth, G. (2012). Haplotype-based variant detection from short-read sequencing. *arXiv preprint arXiv:1207.3907*. doi: 
https://doi.org/10.48550/arXiv.1207.3907

Gerard, D., Ferrão, L. F., Garcia, A. A., & Stephens, M. (2018). Genotyping polyploids from messy sequencing data. *Genetics*, 210(3), 789-807.

