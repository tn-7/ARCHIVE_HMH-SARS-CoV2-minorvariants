# Data files and analysis scripts for "Within-host diversity of SARS-CoV-2 in the context of hospital-associated genomic surveillance"

This repository contains plain text files summarizing deep sequencing data used to identify minor variants (also known as low-frequency variants or intrahost single-nucleotide variants) in samples from Houston Methodist Hospital. Additional files contain information about the characteristics of the samples or the patients they were collected from. Analyses of the distributions and characteristics of minor variants and their associations with patient characteristics were carried out in R, described and annotated here in R markdown files.

To generate minor variant data, upstream of this analysis, the variant calling pipeline timo (https://github.com/GhedinLab/timo) was run on raw fastq files from HMH. The output of timo is a .csv file for each sample, with one line corresponding to each nucleotide position in the genome, with information about the number of reads supporting consensus and minor alleles at each site. The timo files were filtered for all sites with at least 100x read coverage with a minor variant present at any frequency, as a minimum detection limit, and were concatenated into minor_sites_100x_all_20220507.csv to contain minor variant information from all samples in a single file. This file serves as the starting point for further analyses in which minor variants are filtered according to various criteria and levels of stringency, described in the paper and in the R markdown files.

The scripts must be run sequentially, 01-04, to generate the figures and statistics in the paper. Some figures were subsequently annotated in Adobe Illustrator for additional clarity. 

01. Sample quality control and overall tallying of minor variants. 01.5 summary statistics regarding sequencing depth 
02. Association of minor variant richness with patient characteristics
03. Characterization of most prevalent minor variants
04. Re-running analyses in (02) on alternate datasets with different thresholds
