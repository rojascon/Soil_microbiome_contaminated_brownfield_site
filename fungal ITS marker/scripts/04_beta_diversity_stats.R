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

##CODE FOR: running beta-diversity analyses -- PERMANOVAs and CCA analysis
# to identify salient drivers of soil microbiome structure

source(file="00_background.R"); #load necessary packages and specifications


################################################################################
#             1.  Load tidyd metadata and ASV table
################################################################################
load("fungal ITS marker/data/01_brownfield_fITS_ASVtable_filtered.Rdata");
load("fungal ITS marker/data/01_brownfield_fITS_metadata_formatted.Rdata");

#remove the single 6.10m sample because we cannot have an N=1 for that depth
asvdf=asvdf[,!colnames(asvdf) %in% "FUN_ITS_1_55_S55_L001"]
meta=meta[meta$sample_tube!=55,];


################################################################################
#             2. calculate Jaccard distance matrices for statistical tests           
################################################################################

#transpose ASV table so ASVs are columns
asvdf=asvdf[,order(colnames(asvdf))]
asvf=as.matrix(t(asvdf));

###JACCARD distance - all samples [excluding neighbor lot]
tempdf=meta[meta$sample_type!="neighbor_lot",];
asvff=asvf[rownames(asvf) %in% tempdf$SampleID,]
jac=(asvff>0)*1;
jac.dist=vegdist(jac, method="jaccard");

###JACCARD distance - surface samples only [including neighbor lot]
tempdf2=meta[meta$sample_depth==0.15,];
asvf2=asvf[rownames(asvf) %in% tempdf2$SampleID,]
jac2=(asvf2>0)*1;
jac.dist2=vegdist(jac2, method="jaccard");

####JACCARD distance - samples with known pollutant values only
tempdf3=meta[complete.cases(meta$Arsenic),]
asvf3=asvf[rownames(asvf) %in% tempdf3$SampleID,]
jac3=(asvf3>0)*1;
jac.dist3=vegdist(jac3, method="jaccard");


################################################################################
#             3. Run PERMANOVA tests   
################################################################################

#test 1: ALL SAMPLES
# soil microbiome structure ~ sample depth
adonis(jac.dist~sample_depth, strata=tempdf$core, 
       data=tempdf,permutations=999);

#test 2: SURFACE SAMPLES ONLY
#soil microbiome structure ~ contamination category (Y/N)
adonis(jac.dist2~contam_type, strata=tempdf2$core,
       data=tempdf2,permutations=999);

#test 3: SAMPLES WITH KNOWN POLLUTANT VALUES ONLY
#soil microbiome structure~ sample depth + pollutant concentrations
adonis(jac.dist3~Cobalt+Chromium+Lead+Benzo_a_pyrene+Arsenic+sample_depth, 
       strata=tempdf3$core, data=tempdf3,permutations=999,by="margin");


################################################################################
#             4. Run Constrained Correspondence Analysis (CCA)
################################################################################

#shorten sample names in jaccard distance matrix
rownames(jac3)=meta$sample_tube[match(rownames(jac3),meta$SampleID)];

#shorten meta data to only pollutants of interest
cca.meta=tempdf3[,c(2,3,6,9,10,12,25)]; #2 is sample tube label, 3 is depth
rownames(cca.meta)=cca.meta$sample_tube; cca.meta$sample_tube=NULL;

brown.cca <- cca(jac3~ Cobalt+Chromium+Lead+
                   Benzo_a_pyrene+Arsenic+sample_depth, 
                 data=cca.meta, na.action="na.exclude");


################################################################################
#             5. Visualize output of CCA (rudimentary plot)
################################################################################

#see output of CCA analysis
print(summary(brown.cca));   #proportion of variance explained by CCA axes
print(anova.cca(brown.cca, by="terms"));   #which  pollutant predictors are significant?

#plot initial rudimentary plot of CCA output
plot(brown.cca, 
     xlim=c(-3,3), ylim=c(-4,4),
     xlab="CCA1 (23.88%)",
     ylab="CCA2 (20.85%)",
     display=c("ob","cn","wa"));

#a slightly more improved CCA plot; write down the CCA coordinates
# for pollutants that were statistically significant based on CCA output
tp=autoplot(brown.cca)+
  lims(x = c(-5, 4))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=2),
        legend.position="right",
        legend.text=element_text(size=13),
        legend.title=element_text(size=13, face="bold"),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(size=13),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=13, face="bold")); plot(tp)

# my coordinates
#Cobalt= 0, -4
#Arsenic= -1,-4
#depth= 2.5, -4.5


################################################################################
#             6. Visualize output of CCA (manual ggplot plot)
################################################################################

#extract CCA data
ggcca=fortify(brown.cca);
ggcca=ggcca[ggcca$Score=="sites",]; colnames(ggcca)[2]="sample_tube"
cm=merge(ggcca, meta, by="sample_tube")

#set color palette for dots representing sample depth
mycols=c("#fdb462","#bc80bd","#01665e");

#build plot with ggplot2
cca1=ggplot(cm, aes(CCA1, CCA2))+
    geom_point(mapping=aes(fill=as.character(sample_depth)),
               size = 2.5,
               shape=21)+
    scale_fill_manual(values=mycols)+
    labs(y="CCA2 (20.85%)",
         x="CCA1 (23.88%)",
         fill="Depth (m)",
         title="fungal ITS marker")+
    lims(x = c(-8, 3),y=c(-4,4))+
    ##add the pollutant names (specify x and y coordinates)
    annotate(geom="text", x=0.6, y=-4, label="Cobalt", size=3.5)+
    annotate(geom="text", x=-1, y=-4, label="Arsenic",size=3.5)+
    annotate(geom="text", x=-3.3, y=1.6, label="Depth",size=3.5)+
    #add the biplot arrows (specify xend and yend coordinates)
    geom_segment(aes(x = 0, y = 0, xend =0.6, yend = -3.8),
                 arrow = arrow(length = unit(0.3, "cm")), size=0.2)+
    geom_segment(aes(x = 0, y = 0, xend =-1, yend = -3.8),
                 arrow = arrow(length = unit(0.3, "cm")), size=0.2)+
    geom_segment(aes(x = 0, y = 0, xend =-3, yend = 1.3),
                 arrow = arrow(length = unit(0.3, "cm")), size=0.2)+
  #continue tidying plot
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=2),
          legend.position="right",
          legend.text=element_text(size=13),
          legend.title=element_text(size=13, face="bold"),
          plot.title=element_text(size=14, face="bold"),
          axis.text.x=element_text(size=13),
          axis.title.x=element_text(size=13, face="bold"), 
          axis.text.y=element_text(size=13),
          axis.title.y=element_text(size=13, face="bold")); plot(cca1)

##save plot
ggsave(filename="04_CCA_plot.pdf",
       device="pdf",path="fungal ITS marker/figures",
       plot=cca1,
       width=6,
       height=4.5,
       units="in",
       dpi=500);

