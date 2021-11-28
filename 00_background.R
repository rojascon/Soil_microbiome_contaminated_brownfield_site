################################################################################
#
#           Soil Microbiome Composition and Tolerance in Brownfield Soils
#                      
#     Palacios Mejia et.al. 2021. Soil microbial community composition and 
#             tolerance to contaminants in an urban brownfield site
#
#                         Created By: Connie A. Rojas
#                                 On: 19 Mar 2021
#                         Modified By: Connie A. Rojas
#                         Last updated: 28 Nov 2021
#
################################################################################

##CODE FOR: configuring R workspace and printing R version and package versions
#for reader

################################################################################
#             1.  Configure the workspace for subsequent R project scripts                 
################################################################################

#set conditions for R session
rm(list=ls());
options(scipen=999);
options(stringsAsFactors = FALSE) ;

#load necessary packages
library(pacman);
pacman::p_load("car","MASS","dplyr","tidyr","reshape2","vegan","ggplot2",
               "lme4","lmtest","multcomp","grid","phyloseq", "pheatmap",
               "stringr","gridExtra","sjPlot","stringi",
               "ggvegan","Hmisc","fossil");


################################################################################
#             2. Communicate the R version and package versions to reader    
#             sessionInfo() command will provide this information 
################################################################################
print("This code was developed with R version 3.6.2");

print("The packages used and their versions were: pheatmap_1.0.12 | 
phyloseq_1.30.0| multcomp_1.4-15| lmtest_0.9-38| lme4_1.1-26| ggplot2_3.3.3| 
vegan_2.5-7| reshape2_1.4.4| tidyr_1.1.2| dplyr_1.0.3| MASS_7.3-53| car_3.0-10| 
pacman_0.5.1| stringr_1.4.0| gridExtra_2.3| stringi_1.5.3|
      ggvegan_0.1-0 |Hmisc_4.4-2");

