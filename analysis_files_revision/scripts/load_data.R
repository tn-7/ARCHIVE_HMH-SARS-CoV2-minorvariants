# Tung Nguyen
# January 11, 2022

# This file will load all large datasets 
# needed for the startup of each Rmd Notebook after 01_load_new_samples.Rmd.


#### DERIVED FROM 01.RMD and used for later analysis
mcov_info = fread("sample_date_and_run.csv", data.table=F) %>% 
  mutate(COLLECTION_DT=as.Date(COLLECTION_DT, "%m/%d/%y"), 
         MCoVNumber = str_remove(mcov_id, "-")) %>% arrange(COLLECTION_DT) %>% 
  filter(!duplicated(PatientID))


genes <- fread("ntpos_gene_update.csv", data.table = F)
gene_names <- genes %>% pull(gene_id) %>% unique()
genes$gene_id <- factor(genes$gene_id, levels = gene_names)


#all samples sequenced in each run (including those from outside of this study period) -- sometimes relevant for run QC purposes
runs_all_samples <- fread("run_samples.csv", data.table = F) %>% 
  mutate(MCoVNumber=str_remove(`Sample ID`,"-"), run=Run) %>% 
  select(run, MCoVNumber)


d1 = fread('coverage_levels_20220507.csv', data.table = F) %>% 
  filter(duplicated(samplename)) %>% pull(samplename)
d2 = fread('coverage_levels_20220507.csv', data.table = F) %>% 
  filter(samplename %in% d1) %>% 
  filter(!duplicated(samplename)) %>% pull(filename)
coverage_levels<-fread('coverage_levels_20220507.csv', data.table = F) %>% 
  filter(!filename %in% d2) %>% filter(samplename != "MCoV30904") # gets rid of first sample of duplicated sample


nc <- fread("houston_nextclade.tsv", sep='\t', data.table = F) %>% 
  mutate(MCoVNumber=regmatches(seqName, regexpr("[M,R,S,O]CoV.[0-9]+", seqName)) %>% 
           str_remove("-") %>% str_remove("_")) %>% filter(!duplicated(MCoVNumber))
nextclade_bad_samples<-nc %>% filter(qc.overallStatus %in% c("bad")) %>% pull(MCoVNumber)


mcov_samples_filtered <- fread("full_lineage_report_20220507.tsv", data.table=F) %>% 
  mutate(MCoVNumber=mcov_reformat(taxon)) %>% 
  filter(startsWith(MCoVNumber, "MCoV")) %>% 
  filter(MCoVNumber!="MCoV30904") %>% filter(!taxon %in% dup65.66) %>% 
  filter(MCoVNumber %in% mcov_info$MCoVNumber) %>% 
  left_join(mcov_info) %>% select(taxon, MCoVNumber, lineage, scorpio_call, qc_status, 
                                  COLLECTION_DT, INSTRUMENT, INSTRUMENT_RESULT, 
                                  run=run_group) %>% 
  left_join(coverage_levels, by=c("MCoVNumber"="samplename")) %>% 
  filter(!run %in% runs_to_drop) %>% 
  filter(qc_status=="pass") %>% 
  filter(!MCoVNumber %in% nextclade_bad_samples) %>% 
  filter(scorpio_call!="Omicron (BA.1-like)") %>% 
 ##### #main coverage criterion for fair comparisons: X depth over Y percent of the genome
  filter(fraction_100x_coverage>=0.98) %>% droplevels()

minor_variant_sites_threshold <- fread('minor_variants_filtered_100x0.01_50.csv', data.table=F)

n_var <- minor_variant_sites_threshold %>% group_by(MCoVNumber) %>% tally() %>% 
  arrange(desc(n)) #this tally doesn't include any samples with 0 variants, so need to join to original list
samples_n_var<-mcov_samples_filtered %>% left_join(n_var) %>% arrange(COLLECTION_DT) %>% 
  mutate(n_var=tidyr::replace_na(n,0)) 

###### Important for 03.RMD

primers = fread("nCoV-2019.artic_v3.primer.txt", sep="\t", header=FALSE, data.table=F) %>% 
  select(start = V2, end = V3)
primer_positions_v3<-as.numeric()
for (i in 1:nrow(primers)){
  primer_positions_v3<-c(primer_positions_v3, primers[i,]$start:primers[i,]$end)
}

primers = fread("nCoV-2019.artic_v4.primer.bed", sep="\t", 
                header=FALSE, data.table = F) %>% select(start = V2, end = V3)
primer_positions_v4<-as.numeric()
for (i in 1:nrow(primers)){
  primer_positions_v4<-c(primer_positions_v4, primers[i,]$start:primers[i,]$end)
}

primers = fread("V4.1.bed", sep="\t", header=FALSE, data.table=F) %>% 
  select(start = V2, end = V3)
primer_positions_v4.1 <- as.numeric()
for (i in 1:nrow(primers)){
  primer_positions_v4.1 <- c(primer_positions_v4.1, primers[i,]$start:primers[i,]$end)
}

problem_sites = fread("problematic_sites_sarsCov2_v8-20211027.vcf", skip = 88)

primer_positions_all <- c(primer_positions_v3, primer_positions_v4, primer_positions_v4.1,
                          problem_sites$POS) %>% unique()

all_replicates_table_filt <- fread("replicated_samples.csv", data.table=F) %>% 
  filter(MCoVNumber %in% samples_n_var$MCoVNumber)

minors_table<-all_replicates_table_filt %>% 
  filter(major.original %in% c("A","C","T","G")) %>% 
  filter(minor.original %in% c("A","C","T","G")) %>%
  filter(totalcount.original>=100 & minorfreq.original>=0.01 & totalcount.original*minorfreq.original>=50 & 
           tolower(binocheck.original)!="false") %>% filter((!ntpos %in% primer_positions_all)) %>% 
  filter(!ntpos %in% 1:265) %>% filter(!ntpos>29674)

  #was the same minor variant found in the second replicate? if not, set minor frequency to 0 in second rep
minors_table<-minors_table %>% 
  mutate(minorfreq.reseq=if_else(minor.original==minor.reseq, minorfreq.reseq, 0)) %>% 
  mutate(detected_minor_in_repl=if_else(minor.original==minor.reseq,"yes","no"))


##### CLEANING PATIENT DATA FOR 02.RMD
patient_data <- fread("sample_and_patient_data.csv",data.table=F) %>% 
  mutate(MCoVNumber=str_remove(mcov_id, "-")) %>% mutate(collection_date=as.Date(COLLECTION_DT, "%m/%d/%y")) %>% 
  mutate(collection_month=format(as.Date(collection_date), "%Y-%m")) %>%
  mutate(CT=ifelse(INSTRUMENT_RESULT<50, INSTRUMENT_RESULT, NA_integer_)) %>% 
  mutate(vaccine_status=if_else(Vaccine_Status=="No vaccine",0,1)) %>%
  mutate(age18under=if_else(Age_Group=="00-17",1,0)) %>% 
  mutate(age18to54=if_else(Age_Group=="18-54",1,0)) %>% 
  mutate(age55plus=if_else(Age_Group=="55-64"|Age_Group=="65+",1,0)) %>%
  select(MCoVNumber, collection_date, collection_month, run=run_group, CT, 
         ordering_clinic=ORDERING_CLINIC_TYPE,pui=PUI, age18under, age18to54, 
         age55plus, sex=SEX, ethnicity=Ethnicity, obesity=Obesity_YN, 
         chronic_lung_disease=Chronic_Lung_Disease_YN, 
         chronic_liver_disease=Chronic_Liver_Disease_YN, 
         hcw=IS_SURVEILLANCE, chronic_heart_disease=Chronic_Heart_Disease_YN,
         chronic_kidney_disease=Chronic_Kidney_Disease_YN, 
         hypertension=Hypertension_YN, diabetes=Diabetes_YN, 
         cancer=Cancer_YN, hiv=HIV_YN, transplant_patient=Transplant_Patient, 
         vaccine_status, admitted_hospital=Admitted_YN, highest_level=HIGHEST_LEVEL_OF_CARE, 
         max_respiratory_support=MaxRespiratorySupport, mAb=mAb_YN, plasma=Plasma_YN) %>%
  mutate(surveillance = if_else(hcw == "Yes Surveillance",1,0))

factor_columns <- c("collection_month","run","ordering_clinic", "pui", "age18under", 
                    "age18to54", "age55plus","sex","ethnicity","obesity","surveillance", 
                    "chronic_lung_disease","chronic_liver_disease","chronic_heart_disease",
                    "chronic_kidney_disease","hypertension","diabetes","cancer","hiv",
                    "transplant_patient","vaccine_status","admitted_hospital",
                    "highest_level","max_respiratory_support","mAb","plasma") 

patient_data[factor_columns] <- lapply(patient_data[factor_columns], factor)

sample_type = fread("sample_type_PUI.csv", data.table = F) %>% 
  mutate(MCoVNumber=str_remove(mcov_id, "-"))
sample_duration =  fread("timestamp_sample_RNA_extraction_processing.csv",
                        data.table = F) %>% 
  mutate(MCoVNumber=str_remove(mcov_id, "-"))

patient_data[factor_columns]<-lapply(patient_data[factor_columns], factor)

patient_var_tmp = fread("processing/var_aa_ct.txt", data.table = F) %>% 
  left_join(patient_data) %>% mutate(high_counts = n_var > 14) %>% 
 select(-pui) %>% left_join(sample_type) %>% 
  left_join(sample_duration %>% select(-COLLECTION_DT)) %>%
  mutate(vocAlpha=if_else(str_starts(scorpio_call, "Alpha"),1,0), 
         vocDelta=if_else(str_starts(scorpio_call, "Delta"),1,0)) %>% 
  mutate(vocAlpha=as.factor(vocAlpha), vocDelta=as.factor(vocDelta))

#### ACTUAL PATIENT ANALYSIS FOR 02 RMD.
patient_var = patient_var_tmp %>% filter(INSTRUMENT_RESULT < 26)# & 
                                           #minor %in% c("A","C", "T", "G") & 
                                           #major %in% c("A","C","T","G"))
# the reason you get minor or major that isn't defined hence attempting the filter above
# with nucleotides is because 853 samples have no minor variants n_var == 0. So you want
# to keep that to represent the zeros.

patient_var_30 = patient_var %>% filter(n_var < 30)


patient_counts = patient_var %>% select(MCoVNumber,lineage,Duration,COLLECTION_DT:high_counts, 
                                        vocAlpha, vocDelta,PUI) %>% unique
patient_counts_30 = patient_counts %>% filter(n_var < 30)


# Spectra analysis post filter from RMD 03
hmap<-table(patient_var_30$major[patient_var_30$major!=""], 
            patient_var_30$minor[patient_var_30$minor!=""]) %>% data.frame() %>%
  ggplot(aes(x=Var1, y=Var2, fill=Freq/sum(Freq))) + geom_tile(colour = "black") + # grid color
  scale_fill_gradient(low = "white",
                      high = "darkblue") +
  theme_minimal() +   labs(fill = "Fraction",
                           x = "Consensus allele",
                           y = "Minor allele", title="Post-filter minor variants")

