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

##CODE FOR: generating stacked barplots of soil microbiome composition at the
##eukaryotic phylum, and order taxonomic levels for each soil depth

source(file="00_background.R"); #load necessary packages and specifications


################################################################################
#             1.  Load metadata, ASV table, and taxonomy table    
#               after tidying them up (see script 01_prep_data.R)
################################################################################

load("18S rRNA marker/data/01_brownfield_18S_ASVtable_filtered.Rdata");
load("18S rRNA marker/data/01_brownfield_18S_taxonomy_reduced.Rdata");
load("18S rRNA marker/data/01_brownfield_18S_metadata_formatted.Rdata");

#attach taxonomy to the ASV table 
asv_tax=merge(asvdf,mytax,by="row.names"); 
colnames(asv_tax)[1]="ASV";

#make a vector of your sample names for later
samples=meta$SampleID;


################################################################################
#             2. Create Phylum level composition barplots by sample depth
################################################################################
# ###################### exclude neighbor lot samples #####################

#remove neighbor lot samples
samples2=meta$SampleID[meta$sample_type!="neighbor_lot"];

#select eukaryotic taxonomic rank
phylum=asv_tax[,colnames(asv_tax) %in% c(samples2, "Phylum")];
colnames(phylum)[ncol(phylum)]="taxa";

#calculate taxa relative abundances (e.g. proportions)
phylum=aggregate(.~taxa, phylum, sum);  
phylum[,-1] <- lapply(phylum[,-1], function(x) (x/sum(x))*100);
#print(colSums(phylum[-1]));

#keep phyla >0.1% relative abundance across samples
phylum$AVG=rowMeans(phylum[,-1]);
phylum=phylum[phylum$AVG>0.1,];
phylum$AVG=NULL;

#denote the rest of phyla as "Other"
newrow=c(NA, 100-colSums(phylum[2:ncol(phylum)])); 
phylum=rbind(phylum, newrow); 
phylum$taxa=as.character(phylum$taxa);
phylum[nrow(phylum),1]="Other";

#melt data frame for ggplot
pbar<-reshape2::melt(phylum, id.vars="taxa",value.name = "abun");
colnames(pbar)[2]="SampleID";
pbar=merge(pbar, meta, by="SampleID");

#set color-palette
phy_col=c("#1B9E77","#D95F02", "#7570B3", "#E7298A" ,"#66A61E" ,"#E6AB02", "#A6761D" , 
          "#253494", "#FC8D62", "#969696","#8DA0CB", "#ae017e", "#A6D854","#666666" ,
          "#E5C494");

#create plot; 
barphy=ggplot(data=pbar, 
              aes(x=SampleID,y=abun, fill=taxa))+
  geom_bar(stat = "identity")+
  facet_grid(~sample_depth, scales="free_x")+ 
  theme_bw()+ 
  labs(x = "",
       y = "Relative Abundance (%)",
       fill="Eukaryotic Phylum")+
  scale_fill_manual(values=phy_col)+
  theme(legend.position="right", 
        legend.text = element_text(size=12),
        legend.title = element_text(size=14, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x = element_text(size=14, face="bold"),
        strip.text = element_text(size =10, face="bold"));

plot(barphy);

##save image 
ggsave(filename="02_barplot_phylum.pdf",
       device="pdf",path="18S rRNA marker/figures",
       plot=barphy,
       width=8,
       height=5,
       units="in",
       dpi=500);


################################################################################
#             3. Create Order level composition barplots by sample depth               
################################################################################
# ###################### exclude neighbor lot samples #####################

#select Eukaryotic taxonomic rank 
order=asv_tax[,which(names(asv_tax) 
                     %in% c(samples2, "Order"))];
colnames(order)[ncol(order)]="taxa";

#calculate taxa relative abundances 
order=aggregate(.~taxa, order, sum);  
order[,-1] <- lapply(order[,-1], function(x) (x/sum(x))*100);
#print(colSums(order[-1]));

#keep orders >0.4% relative abundance across samples
order$AVG=rowMeans(order[,-1]);
order=order[order$AVG>0.424,]; #0.46 mothur
order$AVG=NULL;

#denote the rest of phyla as "Other"
newrow=c(NA, 100-colSums(order[2:ncol(order)])); 
order=rbind(order, newrow); 
order$taxa=as.character(order$taxa);
order[nrow(order),1]="Other";

#melt data frame for ggplot
obar<-reshape2::melt(order, id.vars="taxa",value.name = "abun");
colnames(obar)[2]="SampleID";
obar=merge(obar, meta, by="SampleID");

#set color-palette
ord_col=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
          "#117777","#77CCCC",
          "#737373","#88CCAA","#777711", "#fc8d59", "#fed976",
          "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788",
          "black","grey");

#create plot; 
barord=ggplot(data=obar, 
              aes(x=SampleID,y=abun, fill=taxa))+
  geom_bar(stat = "identity")+
  facet_grid(~sample_depth, scales="free_x")+ 
  theme_bw()+ 
  labs(x = "",
       y = "Relative Abundance (%)",
       fill="Eukaryotic Order")+
  scale_fill_manual(values=ord_col)+
  theme(legend.position="right", 
        legend.text = element_text(size=12),
        legend.title = element_text(size=14, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x = element_text(size=14, face="bold"),
        strip.text = element_text(size =10, face="bold"));

plot(barord);

##save image 
ggsave(filename="02_barplot_order.pdf",
       device="pdf",path="18S rRNA marker/figures",
       plot=barord,
       width=10.5,
       height=6,
       units="in",
       dpi=500);


