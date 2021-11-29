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
# making heatmaps of PICRUSt2 functional gene pathway abundances for
#         contaminated vs. uncontaminated surface soil microbiomes

source(file="00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Load gene abundance table from PICRUSt2 and sample meta-data
################################################################################

load("16S rRNA marker/data/01_brownfield_16S_metadata_formatted.Rdata");

#download PICRUSt2 file from Google Drive (too big for GitHub) and deposit in
# 16S rRNA marker/data
# https://drive.google.com/file/d/1coLoZVnMS3TSDTYUcV_HW_LLDZ3GERl4/view?usp=sharing

pie=read.delim("16S rRNA marker/data/00_PICRUSt2_output.txt", 
               sep="\t",stringsAsFactors=FALSE,quote="", header=T);

#remove first column from picrust file
pie$pathway=NULL;
pie2=pie;

#convert picrust abundances to relative abundances
pie2[,-1] <- lapply(pie2[,-1], function(x) (x/sum(x))*100);
#colSums(pie2[-1])


################################################################################
#                     2. Plot abundant gene pathways
#         in contaminated vs uncontaminated soils (surface samples only)
################################################################################

#keep pathways >0.65% relative abundance across samples
pie2$AVG=rowMeans(pie2[,-1]);
pie2=pie2[pie2$AVG>0.65,];
pie2$AVG=NULL;

#subset df to only surface samples
snames=meta$SampleID[meta$sample_depth==0.15];
snames=c(snames, "description");
pie3=pie2[,colnames(pie2) %in% snames]

#melt data frame for ggplot
pbar<-reshape2::melt(pie3, id.vars="description",value.name = "abun");
colnames(pbar)[1:2]=c("pathway","SampleID")
pbar=merge(pbar, meta, by="SampleID");

#make a heatmap of pathway abundances
heat1=ggplot(data=pbar, aes(x=SampleID, y=pathway, fill= abun)) + 
  geom_tile() +
  facet_grid(~contam_type,scales="free_x")+ 
  scale_fill_gradient(low="white", high="#88419d") +
  labs(y="PICRUSt2 pathway",
       x="",
       fill="Relative Abundance (%)")+
  theme_bw()+
  theme(legend.position="top",
        legend.text=element_text(size=12),
        legend.title=element_text(size=13, face="bold"),
        axis.text.x=element_blank(),
        axis.title.x=element_text(size=13, face="bold"), 
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=13, face="bold"))+
  theme(strip.text.x = element_text(size = 13)); 
plot(heat1);


################################################################################
#                     3. save plot
################################################################################
ggsave(filename="08_heatmap_picrust.pdf",
       device="pdf",path="16S rRNA marker/figures",
       plot=heat1,
       width=11,
       height=9,
       units="in",
       dpi=500);

