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
load("18S rRNA marker/data/05_alpha_diversity.Rdata");

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
        data=am, na.action=na.fail)

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
#exclude neighbor lot samples
am2=alpha_meta[alpha_meta$sample_depth==0.15,];

#define model
m2=lmer(log(Chao2)~ contam_type + (1|core), 
        data=am2, na.action=na.fail)

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

# only for pollutants that were statistically significant according to model #3
#for 18S rRNA data = Lead

pred1<-ggeffects::ggpredict(m3, terms = c("Lead"));

p1=ggplot(pred1, aes(x=x, y=predicted, ymin=0)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  theme_classic()+
  labs(y="Predicted Chao 2 Richness",
       x="Lead Concentration (mg/kg) ",
       title="18S rRNA marker")+
  theme(axis.title = element_text(size = 12, face="bold"), 
        axis.text = element_text(size = 11),
        plot.title=element_text(size = 12, face="bold")); plot(p1);


################################################################################
#             6. save plots of estimated soil microbiome alpha diversity
################################################################################

ggsave(filename="06_Chao2_plot.pdf",
       device="pdf",path="18S rRNA marker/figures",
       plot=p1,
       width=4.5,
       height=4,
       units="in",
       dpi=500);

