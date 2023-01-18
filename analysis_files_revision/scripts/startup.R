# Tung Nguyen
# January 11, 2022

# This file will load all the packages, common functions & small variables 
# needed for the startup of each Rmd Notebook.


libraries = c("tidyverse", "cowplot", "pROC", "epitools", "glmnet", "nlme", 
              "ggforce", "gggenes", "viridis", "GLMMadaptive",
             "ggpubr", "data.table", "ggrepel", "Biostrings", "ggpubr", 
             "pheatmap", "vegan", "sjPlot", "sjlabelled", "sjmisc",
             "broom", "ggsci", "ggExtra");
invisible(suppressPackageStartupMessages(lapply(libraries, require, character.only = TRUE)));
options(dplyr.summarise.inform = FALSE)

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

runs_to_drop = c("Run_20", "Run_89", "Run_21", "Run_58", "Run_76", "Run_71", 
                 "Run_75", "Run_86", "Run_74", "Run_70", "Run_78", "Run_62", 
                 "Run_13", "Run_85", "Run_65", "Run_34", "Run_87", "Run_11", 
                 "Run_83", "Run_16", "Run_77", "Run_81")

select = dplyr::select