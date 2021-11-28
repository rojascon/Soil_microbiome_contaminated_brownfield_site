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

##CODE FOR: generating rarefaction curve of ASV richness after samples were
#subsampled to 370 sequences each

#NOTE: there might be a bit of a lag to run code from beginning to end
#because of large amount of data displayed

source(file="00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Generate rarefaction curve data using mothur software                 
################################################################################
#In script 05_get_alphadiversity we included code for how to subsample our
#samples to 370 sequences each using the mothur software

#this way all samples have the same # of sequences and we can make comparisons
#Now, we will make a rarefaction plot to double check that our 45,000 sequence
#cutoff was appropriate

#Step1: open data file "fungal ITS marker/
#                 data/05_output_mothur.txt" in Excel and save as .txt
#for some reason, mothur insists on this step ^^^

#Step2: open mothur and run command:
#       rarefaction.single(shared=05_output_mothur.txt, calc=sobs, freq=1)
#Code will take ~20 mins to run

#Step3: rename output file as "08_rarefaction_mothur.txt" and place in the same
#data directory


################################################################################
#             2. Clean and arrange data for ggplot2                 
################################################################################
#read in rarefaction data output by mothur
mrare=read.table("fungal ITS marker/data/08_rarefaction_mothur.txt", 
                 sep="\t", header=TRUE);

#load sample metadata
load("fungal ITS marker/data/01_brownfield_fITS_metadata_formatted.Rdata");

#clean rarefaction data 
raref=mrare[, grepl("X0.03.*",names(mrare))]; 
colnames(raref)=gsub("X0.03.", "", colnames(raref)); 

raref=cbind(mrare$numsampled,raref);
colnames(raref)[1]="numsampled";

#melt data for ggplot2 (timepoint|sample|numOTUs)
rf<- reshape2::melt(raref, id.vars="numsampled",value.name = "ASV_Richness");
colnames(rf)=c("TimeSampled", "SampleID", "ASV_Richness");
rare_meta=merge(rf, meta, by="SampleID");


################################################################################
#             3. Make sample rarefaction curves using ggplot2                 
################################################################################

#set color-palette (1 color for each sample depth)
depth_col=c("#91cf60","#d53e4f","#1f78b4"); 

#set tick marks for x-axis 
reads=seq(0, nrow(mrare), by=50);

#plot! 
rfc=ggplot(data=rare_meta) + 
  geom_line(mapping=aes(x=TimeSampled, 
                        y=ASV_Richness, 
                        col=as.character(sample_depth), 
                        group=SampleID))+ 
  scale_x_continuous(breaks=reads)+
  lims(y=c(0,60))+
  scale_colour_manual(values=depth_col)+
  labs(x = "Number of Reads",
       y = "ASV Richness",
       colour= "Depth (m)")+
  ggtitle("Rarefaction Curves -Fungal ITS profiles")+
  theme_bw() +
  guides(color=guide_legend(override.aes = list(size=3)))+
  theme(panel.background= element_rect(colour = "black", size=2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="right", 
        legend.text = element_text(size=14),
        legend.title=element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=14, face="bold"),
        plot.title =  element_text(size=12, face="bold"));

#plot(rfc);

################################################################################
#             4. save plot              
################################################################################

ggsave(filename="08_rarefaction_curves.pdf",
       device="pdf",path="fungal ITS marker/figures",
       plot=rfc,
       width=7.7,
       height=6,
       units="in",
       dpi=500);

