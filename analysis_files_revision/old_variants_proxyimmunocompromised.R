# Tung Nguyen
# Plot the distribution of old variants over time
# August 15, 2022

# Load files from the PANGO
#install.packages("ggpubr")
setwd('~/github/hmh/old_variants')
library(tidyverse); library(data.table); library(ggpubr); library(lubridate)
##### LOAD FUNCTIONS ######
# Calculating outliers by Tukey's 1.5 * IQR method
iqr_outlier = function(date_vector) {
  lower = quantile(date_vector, 0.25, type = 1)
  upper = quantile(date_vector, 0.75, type = 1)
  iqr = upper - lower
  outlier_upper = lower + 1.5*iqr
  outlier_lower = upper - 1.5*iqr
  interval = outlier_lower %--% outlier_upper
  return(interval)
}



##### LOAD DATA ######
lineages = fread('220815-full_lineage_report.tsv', data.table = F)







# Combine with date file with the dates

patients = fread('../TN-Houston_ScienceAdvances_MS/analysis_files_revision/sample_and_patient_data.csv',
                 data.table=F) %>% mutate(COLLECTION_DT = as.Date(COLLECTION_DT,
  format = '%m/%d/%y'))

omicron = fread('~/github/hmh/omicron_dates/omicron_metadata.tsv', data.table = F) %>%
  rename(id = mcov_id) %>% mutate(COLLECTION_DT = as.Date(COLLECTION_DT))

# Cross-reference Inclusion table (just the first column)
# rename the mcovid to just a number
lineages = lineages %>% mutate(id = taxon %>% parse_number %>% abs)
patients = patients %>% mutate(id = mcov_id %>% parse_number %>% abs) #%>% 
  #filter(INSTRUMENT_RESULT < 26)
# note that the omicron data set is basically also a 'patients' data set
omicronplus_lineages = lineages %>% select(id, lineage, scorpio_call) %>% 
  
  left_join(
    patients %>% 
    select(colnames(omicron)) %>% 
      cbind(source = "prior") %>%
      rbind(cbind(omicron, source="omicron"))
    ) %>% 
  
  na.omit(cols = COLLECTION_DT)# stuff we have sequences for

tmp = omicronplus_lineages

table(omicronplus_lineages$source)
# note check for errors not in the other
#tmp = patients %>% left_join(lineages) %>% select(id, lineage, scorpio_call, COLLECTION_DT) %>%
#  na.omit()

# Plot out the date histograms
tmp %>% select(lineage) %>% unique() %>% nrow() # 258 different lineages
sum(tmp$scorpio_call != "") # 6742 samples have the Greek name spec 'Scorpio Call'
tmp$scorpio_call %>% unique %>% length # 17 Greek letter assignments
omicronplus_counts = tmp$lineage %>% table %>% as.data.frame() %>% 
  rename(Lineage = 1) %>% arrange(-Freq) %>% 
  mutate(Lineage = fct_relevel(Lineage, as.character(Lineage)))

ggbarplot(data = omicronplus_counts %>% filter(Freq > 200), x = "Lineage", y = "Freq") +
  rotate_x_text(angle = 90) + scale_y_continuous(breaks=seq(0,max(omicronplus_counts$Freq),1000)) +
  grids(linetype = "dashed")

# PLOT STACK PLOTS OF DENSITY PLOTS ON TOP OF EACH OTHER 12 of them and the null distribution
# Ridge plots





# Let's just do Alpha and Delta
# alpha = tmp %>% filter(grepl("Alpha",scorpio_call)) %>% select(!scorpio_call) %>%
#   cbind(VOC="Alpha")
# delta = tmp %>% filter(grepl("Delta",scorpio_call)) %>% select(!scorpio_call) %>% 
#   cbind(VOC="Delta")
# 
# alpha_delta = rbind(alpha, delta) %>% mutate(Date = as.Date(COLLECTION_DT, format = "%m/%d/%y")) %>% 
#   select(!COLLECTION_DT)

# ggplot(alpha_delta, aes(x = Date, y = (..count..), fill = VOC)) + 
#   geom_histogram() +
#   scale_x_date(
#                breaks = seq(min(alpha_delta$Date)-5, max(alpha_delta$Date)+5, 30),
#                limits = c(as.Date("2021-01-01"), as.Date("2021-12-31"))) +
#   theme(axis.text.x = element_text(angle = 90))

### RIDGE PLOTS with Dates ####
min_cut = 10
prioritized = tmp %>% filter(lineage %in% (omicronplus_counts %>% filter(Freq > min_cut) %>%
                                             select(Lineage) %>% .[,1]))
table(prioritized$source)


ggplot(data = prioritized, 
       aes(x = COLLECTION_DT, y = lineage, fill=source, color = source, point_color = source)) +
  
  geom_density_ridges(
    jittered_points = T, alpha = 0.01,
    trim = TRUE,
    scale = 1, 
    size = 0.2, 
    rel_min_height = 0.01) +   theme_ridges(center = TRUE)

######### PULL IN THE ENTIRE HMH COHORT DATA FOR BASELINE ########
date_range = omicronplus_lineages$COLLECTION_DT %>% as.character %>% range 
date_interval = date_range[1] %--% date_range[2] 

HMH_all = fread('~/github/hmh/omicron_dates/HMH_omicron_id_date2.tsv',
                data.table = F) %>% rename(id = 1, COLLECTION_DT = 2,
                                           lineage = 3) %>% 
  mutate(COLLECTION_DT = as.Date(COLLECTION_DT)) %>% 
  filter(COLLECTION_DT %within% date_interval)


ggplot(data = HMH_all, 
       aes(x = COLLECTION_DT, y = lineage)) +
  geom_density_ridges(
    jittered_points = T, alpha = 0.01,
    trim = TRUE,
    scale = 1, 
    size = 0.2, 
    rel_min_height = 0.01) +   theme_ridges(center = TRUE)

hmh_all_counts = table(HMH_all$lineage) %>%  as.data.frame() %>% rename(Lineage = 1) %>% 
  arrange(-Freq) %>% 
  mutate(Lineage = fct_relevel(Lineage, as.character(Lineage)))

  
ggbarplot(data = hmh_all_counts %>% filter(Freq > 200), x = "Lineage", y = "Freq") +
  rotate_x_text(angle = 90) + scale_y_continuous(breaks=seq(0,max(omicronplus_counts$Freq),1000)) +
  grids(linetype = "dashed")

# GET ALL US SAMPLES
# between our date ranges
#US_all = fread()


USA_all = fread('~/github/hmh/omicron_dates/test_gisaid_USA.tsv',
                data.table = F) %>% rename(id = 1, COLLECTION_DT = 2,
                                           lineage = 3) %>% 
  mutate(COLLECTION_DT = as.Date(COLLECTION_DT)) %>% 
  filter(COLLECTION_DT %within% date_interval)

ggplot(data = USA_all[1:10000,], 
       aes(x = COLLECTION_DT, y = lineage)) +
  geom_density_ridges(
    jittered_points = T, alpha = 0.01,
    trim = TRUE,
    scale = 1, 
    size = 0.2, 
    rel_min_height = 0.01) + theme_ridges(center = TRUE)

USA_all_counts = table(USA_all$lineage) %>%  as.data.frame() %>% rename(Lineage = 1) %>% 
  arrange(-Freq) %>% 
  mutate(Lineage = fct_relevel(Lineage, as.character(Lineage)))


ggbarplot(data = USA_all_counts %>% .[1:40,], x = "Lineage", y = "Freq") +
  rotate_x_text(angle = 90) + scale_y_continuous(breaks=seq(0,max(omicronplus_counts$Freq),1000)) +
  grids(linetype = "dashed")

###### COMBINE HMH + USA DATA ######
top = 30
myVOI = USA_all_counts %>% select(Lineage) %>% .[,1] %>% 
  .[1:top] %>% 
  as.character()
HMH_USA = HMH_all %>% mutate(source = "GISAID_HMH") %>% 
  rbind(USA_all %>% mutate(source = "GISAID_USA"))

HMH_USA_filtered = HMH_USA %>% filter(lineage %in% myVOI)

# Look at only the top hits in HMH data-set
ggplot(data = HMH_USA_filtered,
       aes(x = COLLECTION_DT, y = lineage, fill = source, alpha = 0.01)) +
  geom_density_ridges() +
  theme_ridges(center = TRUE)

###### COMBINE OMICRONPLUS + USA DATA ########

nonoutlier_lineage_intervals = HMH_USA_filtered %>% group_by(lineage,source) %>% 
  summarize(test = iqr_outlier(COLLECTION_DT)) %>% group_by(lineage) %>% 
  summarize(union_interval = union(test[1],test[2]))

# !!! defining those with infections of older-than-expected strains, as possible 
# 'immuno-compromised' candidates.

jurassic_candidates = omicronplus_lineages %>% left_join(nonoutlier_lineage_intervals) %>% 
  filter(COLLECTION_DT > int_end(union_interval)) %>% 
  mutate(time_discrepancy = COLLECTION_DT-as.Date(int_end(union_interval))) %>%
           arrange(-time_discrepancy)

jurassic_candidates

earlybird_candidates = omicronplus_lineages %>% left_join(nonoutlier_lineage_intervals) %>% 
  filter(COLLECTION_DT < int_start(union_interval)) %>% 
  mutate(time_discrepancy = COLLECTION_DT-as.Date(int_start(union_interval))) %>%
  arrange(time_discrepancy)

earlybird_candidates

# plot distribution of Jurassic & early birds
ggplot(data = jurassic_candidates, 
       aes(x = time_discrepancy)) +
  geom_histogram(bins = 100)

ggplot(data = earlybird_candidates, 
       aes(x = time_discrepancy)) +
  geom_histogram(bins = 100)

# RIDGE PLOT OF THE OUTLIERS WHERE Y IS STRAIN AND X IS THE DISTRIBUTION
outlier_candidates = rbind(jurassic_candidates, earlybird_candidates)

ggplot(data = outlier_candidates,
       aes(x = time_discrepancy, y = lineage, fill = source, alpha = 0.01)) +
  geom_density_ridges() +
  theme_ridges(center = TRUE)

ggplot(data = outlier_candidates %>% filter(abs(time_discrepancy) > 15),
       aes(x = time_discrepancy, y = lineage, fill = source, alpha = 0.01, 
           height = ..count..)) +
  geom_density_ridges(stat = "binline") +
  theme_ridges(center = TRUE)


##
outliers = rbind(jurassic_candidates, earlybird_candidates) %>%
  filter(abs(time_discrepancy) > 30) %>%
  select(lineage) %>% unique() %>% .[,1]

# select the strains, plot them in violin plot and label the outliers
ggplot(HMH_USA %>% filter(lineage %in% outliers), 
       aes(COLLECTION_DT, lineage)) +
  geom_boxplot() #%>% 
  #geom_text_repel(aes(label = lineage), 
#                  na.rm = TRUE, show.legend = F)

# check quality metrics
# count minor variants in each of these
# cross check with filtered


######### SCRATCH ##########
# try it for HMH all
# jurassic_candidates = HMH_all %>% left_join(nonoutlier_lineage_intervals) %>% 
#   filter(COLLECTION_DT > int_end(union_interval)) %>% 
#   mutate(time_discrepancy = COLLECTION_DT-as.Date(int_end(union_interval))) %>%
#   arrange(-time_discrepancy)
