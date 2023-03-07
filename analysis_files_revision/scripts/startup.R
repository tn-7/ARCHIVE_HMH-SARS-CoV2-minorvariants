# Tung Nguyen
# January 11, 2022

# This file will load all the packages, common functions & small variables 
# needed for the startup of each Rmd Notebook.

clear_env = function(){
  rm(list = ls(all.names = TRUE))
  pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
}

clear_env()

libraries = c("tidyverse", "cowplot", "glmnet", 
              "ggforce", "gggenes", "viridis",
             "ggpubr", "data.table", "ggrepel", "Biostrings", "ggpubr", 
             "pheatmap", "ggsci", "ggExtra", "lubridate", "scales",
             "magrittr", "arrow","tableone", "broom", "glmmTMB", "foreach",
             "doParallel", "ggforestplot", "readr", "ComplexHeatmap",
             "tidyHeatmap", "GetoptLong");

if (!require("pacman")) install.packages("pacman")
pacman::p_load(char=libraries)



# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Biostrings")


replace_na = tidyr::replace_na
plot_grid = cowplot::plot_grid

mcov_reformat<-function(badly_formatted_name) {
  MCoVNumber<-regmatches(badly_formatted_name, 
                         regexpr("[M,R,S,O]CoV.[0-9]+", 
                                 badly_formatted_name)) %>% 
  str_remove("-") %>% str_remove("_")
  return(MCoVNumber)
}

# HMH directed us to discard these samples
dup65.66<-c("MCoV-55544_S733", "MCoV-55545_S734", "MCoV-55546_S735", 
            "MCoV-55547_S736", "MCoV-55548_S737", "MCoV-55549_S738", 
          	"MCoV-55550_S739", "MCoV-55551_S740", "MCoV-55552_S741", 
          	"MCoV-55553_S742", "MCoV-55554_S743", "MCoV-55555_S744", 
          	"MCoV-55556_S745", "MCoV-55557_S746", "MCoV-55558_S747", 
          	"MCoV-55559_S748", "MCoV-55560_S749", "MCoV-55561_S750", 
          	"MCoV-55562_S751", "MCoV-55563_S752", "MCoV-55564_S753", 
          	"MCoV-55565_S754", "MCoV-55566_S755", "MCoV-55567_S756")

runs_to_drop = readRDS("processing/runs_to_drop.rds")

select = dplyr::select