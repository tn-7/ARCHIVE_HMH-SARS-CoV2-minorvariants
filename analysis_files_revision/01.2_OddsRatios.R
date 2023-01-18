# Tung Nguyen
# 10 August 2022
# EDA of the HMH Cohort
install.packages(c('data.table', 'tidyverse'))
setwd('~/git/hmh/')
library(data.table); library(tidyverse)
#### PART 1 - what kind of PANGO variants are there in this cohort? ####
# What is the time distribution?
# Anyone who is an outlier?
#### Patient meta-data ####
patients = fread('TN-Houston_ScienceAdvances_MS/analysis_files_revision/sample_and_patient_data.csv',
                 data.table=F)
patients_tmp = patients %>% filter(INSTRUMENT_RESULT > 0)
patients_tmp = patients %>% filter(INSTRUMENT_RESULT > 0 &
                                     INSTRUMENT_RESULT < 26)
plot(density(patients_tmp$INSTRUMENT_RESULT))
abline(v=26, col = 'red')
# Count the number of different
plot(density(patients_tmp$COLLECTION_DT))
# Count the respirator support
plot(density(patients_tmp$))
barplot(colSums(patients_tmp %>%
                  select(Obesity_YN:Admitted_YN,-Vaccine_Status)), las=2,
        cex.names = 0.4)