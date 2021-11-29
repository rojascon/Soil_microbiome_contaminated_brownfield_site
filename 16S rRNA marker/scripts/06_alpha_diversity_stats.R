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

##CODE FOR: 
# A) running model to assess the effects of pollutants on
#             soil microbiome alpha diversity 
# B) make plots of soil microbiome alpha diversity ~ pollutant concentrations

source(file="00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load alpha-diversity data
################################################################################
load("16S rRNA marker/data/05_alpha_diversity.Rdata");

#remove samples without pollutant concentrations
tempam=alpha_meta[complete.cases(alpha_meta$Arsenic),];


################################################################################
#             2. Linear model #1 (all samples, exclude neighbor lot)
#           soil microbiome diversity ~ sample depth (sample core random term)
################################################################################
#exclude neighbor lot samples
am=alpha_meta[alpha_meta$sample_type!="neighbor_lot",];

#define model
m1=lmer(log(Chao2)~ sample_depth + (1|core), 
        data=am, na.action=na.fail);

#evaluate model fit via residuals
plot(fitted(m1), resid(m1));

#summary of result
summary(m1)$coeff; 

#assess statistical significance of predictors
Anova(m1);


################################################################################
#             3. Linear model #2 (only surface samples, include neighbor lot)
#             soil microbiome diversity ~ contamination category (Y/N)
################################################################################
#exclude non-neighbor lot samples
am2=alpha_meta[alpha_meta$sample_depth==0.15,];

#define model
m2=lmer(log(Chao2)~ contam_type + (1|core), 
        data=am2, na.action=na.fail);

#evaluate model fit via residuals
plot(fitted(m2), resid(m2));

#summary of result
summary(m2)$coeff; 

#assess statistical significance of predictors
Anova(m2);


################################################################################
#             4. Linear model #3 (only samples with pollutant concentrations)
#                 soil microbiome diversity ~ various pollutants
################################################################################

#define model
m3=lmer(log(Chao2)~log(Arsenic)+log(Cobalt)+log(Chromium)+log(Lead)+
        log(Benzo_a_pyrene) + sample_depth+ (1|core), 
      data=tempam,na.action = "na.fail");

#evaluate model fit via residuals
plot(fitted(m3), resid(m3));

#summary of result
summary(m3)$coeff; 

#assess statistical significance of predictors
Anova(m3);


################################################################################
#         5. Make plots of estimated soil microbiome richness 
#                     ~ pollutant concentrations
################################################################################

#for 16S rRNA data, none of the hydrocarbons or heavy metals significantly predicted
# soil microbiome richness, so there is nothing to plot here


################################################################################
#             6. save plots of estimated soil microbiome alpha diversity 
################################################################################

#no plots to save
