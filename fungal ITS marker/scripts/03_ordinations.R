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

##CODE FOR: generating PCoA ordinations of the soil microbiome 
# based on Jaccard distances

source(file="00_background.R"); #load necessary packages and specifications


################################################################################
#             1.  Load tidyd metadata and ASV table
################################################################################

load("fungal ITS marker/data/01_brownfield_fITS_ASVtable_filtered.Rdata");
load("fungal ITS marker/data/01_brownfield_fITS_metadata_formatted.Rdata");

#remove singleton and doubleton sequences
asvdf$seq=rowSums(asvdf); asvdf=asvdf[asvdf$seq>2,];
asvdf$seq=NULL;


################################################################################
#             2. calculate Jaccard distance matrix           
################################################################################

#transpose ASV table so ASVs are columns
asvf=as.matrix(t(asvdf));

###JACCARD distance [turn ASV table to presence/absence first]
jac=(asvf>0)*1;
print(rowSums(jac));
jac.dist=vegdist(jac, method="jaccard");


################################################################################
#             3. plot PCoA ordination color-coded by sample depth
#                       include neighbor lot samples
################################################################################

#calculate coordinates for PCoA
pcoa_dec=cmdscale(jac.dist, eig=TRUE);
pcoa=as.data.frame(pcoa_dec$points);
colnames(pcoa)=c("Axis1","Axis2");
pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "SampleID");
pcoa_met=merge(pcoa,meta,by="SampleID"); 
pcoa_met$sample_depth=as.character(pcoa_met$sample_depth);

#add a "neighbor lot" entry to the sample depth category
pcoa_met$sample_depth2=pcoa_met$sample_depth;
pcoa_met$sample_depth2[pcoa_met$sample_type=="neighbor_lot"]="neighbor_lot";

#factors
pcoa_met$sample_depth2=factor(pcoa_met$sample_depth2, 
                              levels=c("0.15", "1.52","4.6","6.1",
                                       "neighbor_lot"))

#calculate % explained by PC1 and PC2
pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; 
ax1=format(pcoa_per[1], digits=2, nsmall=2);
ax2=format(pcoa_per[2], digits=2, nsmall=2);

#set color-palette for plot
depth_col=c("#91cf60","#d53e4f","#1f78b4","#f781bf","#878787"); 

#plot the PCoA-color coded by sample depth
pcoa1=ggplot(pcoa_met, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=sample_depth2),
             size = 4,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="Depth (m)",
       title="fungal ITS marker (all samples)")+
  scale_fill_manual(values=depth_col)+
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
        axis.title.y=element_text(size=13, face="bold"))

plot(pcoa1);


################################################################################
#             4. plot PCoA ordination contaminated  vs not contaminated
#                          SURFACE SAMPLES ONLY
################################################################################

#retain only surface sample data
pcoa_met2=pcoa_met[pcoa_met$sample_depth=="0.15",];


#add a "neighbor_lot" entry to the contam_type category
pcoa_met2$contam_type=as.character(pcoa_met2$contam_type);
pcoa_met2$contam_type[pcoa_met2$sample_type=="neighbor_lot"]="neighbor_lot";
pcoa_met2$contam_type=factor(pcoa_met2$contam_type, 
                             levels=c("contaminated","not_contaminated",
                                      "neighbor_lot"));

#plot the PCoA contaminated vs not contaminated [surface samples only]
pcoa2=ggplot(pcoa_met2, aes(Axis1, Axis2))+
  geom_point(mapping=aes(fill=contam_type),
             size = 4,
             shape=21)+
  labs(y=paste("PC2 (",ax2,"%)",sep=""),
       x=paste("PC1 (",ax1,"%)",sep=""),
       fill="",
       title="fungal ITS marker (surface samples)")+
  scale_fill_manual(values=c("#fdc086","#beaed4","#525252"))+
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
        axis.title.y=element_text(size=13, face="bold"));

plot(pcoa2);


################################################################################
#             5. Save PCoA ordinations                  
################################################################################

mypcas=arrangeGrob(pcoa1, pcoa2, nrow=1); 
ggsave(filename="03_jaccard_pcoas.pdf",
       device="pdf",path="fungal ITS marker/figures",
       plot=mypcas,
       width=10.8,
       height=3.5,
       units="in",
       dpi=500);






