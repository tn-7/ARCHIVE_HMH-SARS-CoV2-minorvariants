# Within-host diversity of SARS-CoV-2 in the context of hospital-associated genomic surveillance

Alexandra Mushegian et al

This repository contains the R script and data files used to generate the figures and analysis presented in the manuscript.

## Requirements

Analysis requires installation of R, Rstudio, and the following packages: tidyverse, cowplot, gggenes, nlme, randomForest, pROC, glmnet 

Analyses were run with the following R session info:

```
R version 4.1.2 (2021-11-01)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 11.6.3

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.7          pillar_1.6.4        compiler_4.1.2      plyr_1.8.6          tools_4.1.2         digest_0.6.29      
 [7] evaluate_0.14       lifecycle_1.0.1     tibble_3.1.6        gtable_0.3.0        pkgconfig_2.0.3     rlang_0.4.12       
[13] DBI_1.1.1           cli_3.1.0           rstudioapi_0.13     yaml_2.2.1          xfun_0.28           fastmap_1.1.0      
[19] dplyr_1.0.7         knitr_1.36          pROC_1.18.0         generics_0.1.1      vctrs_0.3.8         grid_4.1.2         
[25] cowplot_1.1.1       tidyselect_1.1.1    glue_1.5.1          R6_2.5.1            fansi_0.5.0         rmarkdown_2.11     
[31] ggplot2_3.3.5       purrr_0.3.4         magrittr_2.0.1      scales_1.1.1        ellipsis_0.3.2      htmltools_0.5.2    
[37] randomForest_4.6-14 assertthat_0.2.1    colorspace_2.0-2    utf8_1.2.2          munsell_0.5.0       crayon_1.4.2
```

## Running the analysis

Open houston_analysis.Rmd in Rstudio, set the working directory to the directory containing the data files, and run the script. Figures will be generated within the R Notebook. 

The first step of the analysis joins the split files in `./variant_lists` into three main data files: `variant_sites_all.csv `, `consensus_and_variant_sites.csv`, and `replicated_samples.csv`

## Files

`variant_sites_all.csv `, `consensus_and_variant_sites.csv`, and `replicated_samples.csv` were created from outputs of the timo pipeline, run on sample fastq files from Bioproject PRJNA767388. Timo creates .csv files with sequencing depth and minor allele frequencies at every nucleotide position in the SARS-CoV-2 genome. 

From the output files of all samples, only locations with at least 100x sequencing coverage and consensus or minor changes were collected into a single file (`variant_sites_all.csv ` for sites with minor variants and  `consensus_and_variant_sites.csv` for sites with either consensus changes from the reference or minor variants). Further filtering of variant sites according to depth or frequency or other criteria occurs in the R script.

`replicated_samples.csv` contains the concatenated entire timo files of samples that were sequenced twice as technical replicates. 

Additional input files in this folder contain information on sites of particular interest in the genome, sample collection or patient information, SARS-CoV-2 reference genome annotations, or other quality control parameters, as explained in the script and in the manuscript's Materials and Methods.
