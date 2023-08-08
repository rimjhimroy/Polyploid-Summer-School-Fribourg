# Polyploid SNP calling Protocol

Created by: Rimjhim Roy Choudhury, Armel Salmon, and Christian Parisod  
Date: 2023-07-18

## Table of Contents

- [Introduction](#introduction)
- [Overview of the Protocol](#overview-of-the-protocol)
- [Requirements](#requirements)
- [Data access and Resource allocation](#data-access-and-resource-allocation)
- [How to submit jobs on the IBU cluster](#how-to-submit-jobs-on-the-ibu-cluster)
- [Step 1: Preparing the Reference Genome](#step-1-preparing-the-reference-genome-already-done-for-you)
- [Step 2: Preparing the bam files](#step-2-preparing-the-bam-files)
- [Step3: Genotype Calling with GATK](#step-3-genotype-calling-with-gatk)
- [Step 4: Genotype Calling with FreeBayes](#step-4-genotype-calling-with-freebayes)
- [Step 5: Genotype Calling with Octopus](#step-5-genotype-calling-with-octopus)
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

In this protocol, we will perform genotype calling on polyploid genomes using GATK, FreeBayes, Octopus and Updog.

For this session of practicals, we have generated polyploid datasets by using two reference sequences: `At_Chr4_140-150.fas` and `At_Chr4_35-45.fas` (also known as the knob region).

Next, we created randomly generated 4x sequences with 0.3% polymorphic sites based on the reference sequences using a custom python script set up to randomly select 0.3% of site to mutate from reference sequences (substitution model: C<->G; A<->T) with 1/4 of Aaaa; 1/2 of AAaa; and 1/4 AAAa. These sequences were stored in the files `At_Chr4_140-150_4x.fas` and `At_Chr4_35-45_4x.fas`. We recorded the SNP positions within these 4x sequences in the files `At_Chr4_140-150_4x_snp_pos.txt` and `At_Chr4_35-45_4x_snp_pos.txt`. VCF files were generated from aligned sequences thanks to the *snp-sites* software (<https://github.com/sanger-pathogens/snp-sites>).

Additionally, to simulate realistic sequencing data, we generated *in silico* Illumina NovaSeq paired reads at different coverage levels. The coverage levels used were 150x, 100x, 50x, 20x, and 5x. These reads represent the experimental data that will be used in further analyses.

Simulated reads were generated using the *InSilicoSeq* software (<https://insilicoseq.readthedocs.io/en/latest/>) that can be installed with conda (`conda install -c bioconda insilicoseq` or `conda create -p /your_path/insilicoseq insilicoseq`); to activate conda, then run (`conda activate /your_path/insilicoseq`).
For example, to generate and simulate a 100x coverage of a 1Mb long sequence with 150bp long paired-reads: (`iss generate --n_reads 670000 --draft At_Chr4_35-45_4x.fas --model novaseq --output At_Chr4_35-45_4x_reads_100x`). This software generates reads and simulates sequencing errors using an "Error Model" that can be set by the user. We used the prebuilt error model (for Novaseq sequencing) with the command line above.

To align the generated reads to the reference sequences, we employed the *Bowtie2* software (<https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>) with specific parameters (`--sensitive -N 1 --no-contain --no-overlap --no-discordant`). The resulting mapped reads were stored in BAM files, which serve as the final files required for this session. To visualize alignments (BAM files), you can use the *Tablet* software
(<https://ics.hutton.ac.uk/tablet/download-tablet/>)

This summer school protocol introduces three widely used tools: *GATK* (Genome Analysis Toolkit), *FreeBayes*, *Octopus*, and *UpDog*, which provide powerful solutions for variant calling and genotyping.

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

- Octopus (Cooke et al., 2021) is a variant caller specifically designed to handle various ploidies. It utilizes a polymorphic Bayesian genotyping model within a unified haplotype-aware framework. It demonstrates the ability to handle polyploid genomes, including those with higher ploidies such as tetraploid and hexaploid, and showcases superior performance compared to other popular tools in terms of genotyping accuracy.
  - Key Features:
    - Demonstrates high accuracy in calling germline variants in individuals
    - Can detect SNPs, indels, and small complex replacements such as microinversions.
    - Determine phased genotypes of arbitrary ploidy, even those containing somatic mutations

- UpDog (Gerard et. al., 2021) is a novel tool that is designed for genotyping polyploids. It addresses various challenges commonly encountered in NGS data analysis, such as allele bias, overdispersion, and sequencing errors. The name "updog" stands for "Using Parental Data for Offspring Genotyping" and reflects the initial focus of the method on full-sib populations. However, Updog has evolved to be applicable to a wider range of populations, making it a versatile tool for polyploid genotyping from NGS data.
  - Key Features:
    - Specific focus on variant calling in polyploid genomes.
    - Probabilistic model for accurate variant identification.
    - Robust handling of different ploidy levels.

<div style="page-break-after: always;"></div>

## Requirements

- Download the zip file: <https://drive.switch.ch/index.php/s/1dj6yPey42rCmyB>, and unzip in your renku session
- Polyploid reference genome fasta files `At_Chr4_140-150.fas` and `At_Chr4_35-45.fas` in ref/ folder
- Sorted Bam files of mapped simulated reads at different coverage levels: 5x, 20x, 50x, 100x, and 150x, mapped to the `At_Chr4_140-150` and `At_Chr4_35-45` reference fasta sequences in mapping/bowtie_accurate/ folder
- Truth set vcf of known variants `At_Chr4_35-45_4x_truth_set.vcf.gz` and `At_Chr4_35-45_4x_truth_set.vcf.gz` in the truth_set/ folder
  
## Step 1: Preparing the Reference Genome (already done for you)

1. Index the reference genome using samtools:

   ```bash
   samtools faidx ref/reference.fasta
   ```

2. Create a sequence dictionary for the reference genome using Picard:

   ```bash
   gatk CreateSequenceDictionary -R ref/reference.fasta -O ref/reference.dict
   ```

## Step 2: Preparing the bam files

1. Mark duplicates in the aligned BAM files using GATK. Replace the input.bam and deduplicated.bam file names with the appropriate ones based on the reference fasta and the coverage level:

   ```bash
   gatk MarkDuplicates -I mapping/bowtie_accurate/input.bam -O mapping/bowtie_accurate/deduplicated.bam -M mapping/bowtie_accurate/metrics.txt
   ```

2. Add read groups to the BAM files using GATK. Do not forget to replace the file names with the appropriate ones based on the reference fasta and the coverage level:

   ```bash
   gatk AddOrReplaceReadGroups -I mapping/bowtie_accurate/deduplicated.bam -O mapping/bowtie_accurate/dedup_rg.bam --CREATE_INDEX true -LB samplename -PL illumina -PU samplename -SM samplename
   ```

## Step3: Genotype Calling with GATK

1. Perform variant calling using GATK HaplotypeCaller:

   ```bash
   mkdir -p vcf/gatk
   gatk HaplotypeCaller -R ref/reference.fasta -I mapping/bowtie_accurate/dedup_rg.bam -O vcf/gatk/gatk_raw_variants.vcf -ploidy 4 --native-pair-hmm-threads 1 -stand-call-conf 10
   ```

2. Apply filtering on the raw variant calls using GATK VariantFiltration, adjusting the parameters as needed:

   ```bash
   gatk VariantFiltration -R ref/reference.fasta -V vcf/gatk/gatk_raw_variants.vcf -O vcf/gatk/gatk_filtered_variants.vcf \
   -filter 'QD < 2.0' --filter-name 'QD2' \
   -filter 'QUAL < 30' --filter-name 'QUAL30' \
   -filter 'GQ < 5' --filter-name 'GQ5' \
   -filter 'FS > 60.0' --filter-name 'FS60' \
   -filter 'SOR > 3.0' --filter-name 'SOR3' \
   -filter 'MQ < 40.0' --filter-name 'MQ40' \
   -filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' \
   -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8'
   ```
  
## Step 4: Genotype Calling with FreeBayes

1. Perform variant calling using FreeBayes:

   ```bash
   mkdir -p vcf/freebayes
   freebayes -f ref/reference.fasta -p 4 --genotype-qualities mapping/bowtie_accurate/dedup_rg.bam > vcf/freebayes/freebayes_raw_variants.vcf
   ```

2. Apply filtering on the raw variant calls using BCFtools, adjusting the parameters as needed:

   ```bash
   bcftools filter -i 'QUAL > 1 && GQ > 1 && SAF > 0 && SAR > 0 && RPR > 1 && RPL > 1' -s FAIL -o vcf/freebayes/freebayes_filtered_variants.vcf vcf/freebayes/freebayes_raw_variants.vcf
   ```

## Step 5: Genotype Calling with Octopus

1. Perform variant calling using Octopus:

   ```bash
   export PATH=/home/rstudio/octopus/bin/:$PATH
   mkdir -p vcf/octopus
   octopus --reference ref/reference.fasta -I mapping/bowtie_accurate/sorted_dedup_rg.bam -o vcf/octopus/octopus_variants.vcf -P 4 --threads 8 --max-genotypes 20000
   ```
  
  **Note:** The highest coverage file will run for 15-20 mins.

2. Call longer haplotypes using Octopus:

   ```bash
   octopus --reference ref/reference.fasta -I mapping/bowtie_accurate/sorted_dedup_rg.bam -o vcf/octopus/octopus_variants_longhaps.vcf -P 4 --threads 8 --disable-denovo-variant-discovery --source-candidates vcf/octopus/octopus_variants.vcf --use-filtered-source-candidates  -x 400 --lagging-level OPTIMISTIC --backtrack-level AGGRESSIVE 
   ```
  
  **Note:** The highest coverage file will run for 1 hour 30 mins.

## Step 6: Genotype Calling with Updog (using unfiltered GATK called sites as input)

1. Use the variant calls generated by GATK HaplotypeCaller as input for Updog.

2. remove sites with "\*" in ALT field from the vcf file as R package `VariantAnnotation` cannot handle it properly.
What are the \* in the vcf file? <https://gatk.broadinstitute.org/hc/en-us/articles/360035531912-Spanning-or-overlapping-deletions-allele->

```bash
grep -v "*" output/gatk/gatk_raw_variants.vcf > output/gatk/gatk_raw_variants_fix.vcf
```

3. Use the script in the scripts folder to genotype and filter the variants using Updog, replace `sample_coverge` to 5x, 20x, 50x, 100x, and 150x, depending on the coverage of the GATK vcf file. Change `sample_coverge` with "5x", "20x", "50x", "100x", and "150x" depending on the input GATK vcf coverage. The script will output a filtered vcf file in the vcf/ folder:

   ```bash
   mkdir -p vcf/updog
   conda activate snp_calling
   cd vcf/updog
   Rscript --vanilla scripts/updog.R -i ../gatk/gatk_raw_variants_fix.vcf-s sample_coverge -p 4 -n 8
   ```
  
## Step 7: Comparing to the Truth Set using vcfeval

1. Create a reference SDF file using the reference genome fasta file:

   ```
   rtg format -o ref/reference_SDF ref/reference.fasta 
   ```

2. gzip all the called vcf files

      ```
      bgzip vcf/gatk/gatk_filtered_variants.vcf
      tabix -p vcf vcf/gatk/gatk_filtered_variants.vcf.gz
      bgzip vcf/freebayes/freebayes_filtered_variants.vcf
      tabix -p vcf vcf/freebayes/freebayes_filtered_variants.vcf.gz
      bgzip vcf/octopus/octopus_variants_longhaps.vcf
      tabix -p vcf vcf/octopus/octopus_variants_longhaps.vcf.gz
      bgzip vcf/updog/updog_filtered_variants.vcf
      tabix -p vcf vcf/updog/updog_filtered_variants.vcf.gz
      ```

3. Compare the variant calls from each method to the truth set using rtg vcfeval, change the `method_calls.vcf.gz` with the final vcf file of each method, truth_set.vcf.gz with the truth set vcf file for the reference fasta, and output_dir with an appropriate name based on the reference, the method used, and the coverage:

   ```
   rtg vcfeval -b truth_set.vcf.gz -c method_calls.vcf.gz -o output_dir -t reference_SDF --sample-ploidy 4  --ref-overlap
   ```

4. Compare the variant calls from each method to the truth set using rtg vcfeval at the base call level, without considering the genotypes

   ```
   rtg vcfeval -b truth_set.vcf.gz -c method_calls.vcf.gz -o output_dir -t reference_SDF --sample-ploidy 4  --ref-overlap --squash-ploidy --output-mode "annotate" 
   ```

   - **TP**: True Positives - The number of variant calls that are correctly identified as true positives compared to the truth set.
   - **FP**: False Positives - The number of variant calls that are incorrectly identified as positive (i.e., present) when they are not present in the truth set.
   - **FN**: False Negatives - The number of variant calls that are present in the truth set but are not identified as positive (i.e., missed) by the variant calling methods.
   - **Precision**: The proportion of true positive calls out of all the positive calls made by the variant calling methods (TP / (TP + FP)).
   - **Sensitivity**: The proportion of true positive calls out of all the variants present in the truth set, representing the ability of the method to identify true positive variants. (TP / (TP + FN)).
   - **F-measure**: The harmonic mean of precision and recall, providing a single metric that balances the trade-off between precision and recall. (2 \*((precision\* recall) / (precision + recall)))
   - **Specificity**: The ability of the method to correctly identify true negative calls, calculated as 1 - false positive rate. (TN / (TN + FP))
   - **Accuracy**: The overall correctness of the variant calls, calculated as (TP + TN) / (TP + FP + TN + FN).

5. Analyze and interpret the comparison results to evaluate the performance of each method. Consider the precision, sensitivity, F-measure, specificity, and accuracy values to assess the accuracy and reliability of the variant calling methods.

6. Generate the ROC curve using rtg rocplot:

   ```
   rtg rocplot --png output_dir/roc.png output_dir/snp_roc.tsv.gz
   ```

   - **ROC curve**: A graphical plot of the true positive rate (sensitivity) against the false positive rate (1 - specificity) at various threshold settings. The ROC curve provides a visual representation of the trade-off between sensitivity and specificity for each method. The closer the curve is to the top left corner, the better the performance of the method.
   **Note**: Take a look at the different options of `rtg vcfeval` specially `--all-records` if you are using a filtered vcf file and `--vcf-score-field` which you can change and subsequently generate the roc plots for different thresholds.

[Need Preprocessed output vcf files? Download the zip file: <https://drive.switch.ch/index.php/s/ZJvt7XPrUBieFrT>, and unzip in your renku session]

## Conclusion

This protocol provides a step-by-step guide for genotyping polyploid genomes using GATK, FreeBayes, and Updog. By following these steps, you will be able to compare the variant calls to a truth set and assess the performance of each method. Adjust the parameters and filtering criteria based on your specific needs and data characteristics.

We acknowledge the significant contribution of Cooke et al. (2022) in their paper titled "Benchmarking small-variant genotyping in polyploids." This study presents a comprehensive evaluation and comparison of various genotyping methods tailored specifically for polyploid organisms. By systematically analyzing the strengths and limitations of each method, Cooke et al. provide valuable insights into their performance in accurately determining genotypes within polyploid populations.

Enjoy your journey into polyploid genomics!

## References

Cooke, D.P., Wedge, D.C., Lunter, G. Benchmarking small-variant genotyping in polyploids (2022). *Genome research*, 32(2):403-408.
  
Cooke, D. P., Wedge, D. C., & Lunter, G. (2021). A unified haplotype-based method for accurate and comprehensive variant calling. Nature biotechnology, 39(7):885–892.
  
McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., ... & DePristo, M. A. (2010). The genome analysis toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome research*, 20(9), 1297-1303.
  
DePristo, M. A., Banks, E., Poplin, R., Garimella, K. V., Maguire, J. R., Hartl, C., ... & Daly, M. J. (2011). A framework for variation discovery and genotyping using next-generation DNA sequencing data. *Nature genetics*, 43(5), 491-498.
  
Van der Auwera GA, Carneiro MO, Hartl C, Poplin R, Del Angel G, Levy-Moonshine A, Jordan T, Shakir K, Roazen D, Thibault J, Banks E, Garimella KV, Altshuler D, Gabriel S, DePristo MA. From FastQ data to high confidence variant calls: the Genome Analysis Toolkit best practices pipeline. *Current protocols in bioinformatics*. 2013;43(1110):11.10.1-11.10.33.
  
Poplin R, Ruano-Rubio V, DePristo MA, Fennell TJ, Carneiro MO, Auwera GA, Kling DE, Gauthier LD, Levy-Moonshine A, Roazen D, Shakir K, Thibault J, Chandran S, Whelan CW, Lek M, Gabriel S, Daly MJ, Neale BM, MacArthur DG, Banks E. Scaling accurate genetic variant discovery to tens of thousands of samples. *bioRxiv*. doi: <https://doi.org/10.1101/201178>
  
Garrison, E., & Marth, G. (2012). Haplotype-based variant detection from short-read sequencing. *arXiv preprint arXiv:1207.3907*. doi:
<https://doi.org/10.48550/arXiv.1207.3907>
  
Gerard, D., Ferrão, L. F., Garcia, A. A., & Stephens, M. (2018). Genotyping polyploids from messy sequencing data. *Genetics*, 210(3), 789-807.
