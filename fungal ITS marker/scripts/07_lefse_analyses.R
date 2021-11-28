#################################################################################
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
#A) generating table of microbial genera abundances for LEfSe 
#       LEfSe will identify the genera that are differentially enriched 
#       between contaminated and uncomtaminated surface samples
#B) running LEfSe on the Galaxy platform
#C) plotting the output from LEfSe in the form of diverging plots

source(file="00_background.R"); #load packages & specifications


################################################################################
#             1. Load filtered ASV abundance table, ASV taxonomy table and 
#                         filtered metadata table                 
################################################################################
load("fungal ITS marker/data/01_brownfield_fITS_ASVtable_filtered.Rdata");
load("fungal ITS marker/data/01_brownfield_fITS_metadata_formatted.Rdata");
load("fungal ITS marker/data/01_brownfield_fITS_taxonomy_reduced.Rdata");


################################################################################
#             2. Calculate relative abundances of each Genus   
#               and only retain those >0.1% average relative abundance 
################################################################################

#attach taxonomy to the ASV table 
asv_tax=merge(asvdf,mytax[,c(1,6)],by="row.names"); 

#delete the first row, and all columns that are not samples or "Genus"
asv_tax$Row.names=NULL; asv_tax$Kingdom=NULL;

#aggregate genus-level abundances (e.g. sum abundances of taxa classified
#as the same genus)
tempdf=aggregate(.~Genus, asv_tax, sum);  
rownames(tempdf)=tempdf$Genus; tempdf$Genus=NULL

#calculate relative abundances
tempdf=as.matrix(tempdf);
tempdf=prop.table(tempdf,2);
tempdf=data.frame(tempdf);

#only keep genera with > 0.1% mean abundances across samples
tempdf$rM=rowMeans(tempdf); tempdf=tempdf[tempdf$rM>0.001,];


################################################################################
#             3. Make data frame for LEfSe                 
################################################################################
#surface samples only, include neighbor lot
mys=meta$SampleID[meta$sample_depth==0.15]
df1=tempdf[ , (names(tempdf) %in% mys)]
colnames(df1)=meta$contam_type[
  match(colnames(df1),meta$SampleID)];

#save df for LEfSe
write.table(df1, "fungal ITS marker/data/07_contam_surface_lefse.txt", 
            col.names=NA,row.names=TRUE,sep="\t");


##open file in Excel, and add "contamination" to first cell


################################################################################
#             4. Run LEfSe on the Galaxy platform
#               https://huttenhower.sph.harvard.edu/galaxy/
################################################################################

#Step1: Upload data frames from above by going to Get Data > 
# upload file from your Computer

#Step2: Go to LEfSe > Format Data for LEfSe 
#Select whether the vectors (features and meta-data information) are listed in 
#rows or columns: select ROWS
#Select which row to use as class: use #1 contamination or #1 depth
#click EXECUTE

#Step3: Go to LDA Effect Size (LEfSe)
#set the Alpha value for the factorial Kruskal-Wallis test among classes = 0.05
#click EXECUTE

#Step4: compile all tables output by LEFSE into one file = "07_lefse_output.csv"
#colnames are Genus	| Mean	| group	| LDA	| p_value
#save only rows with LDA > 3
#Place that file in the data folder


################################################################################
#                            6. Plot LEfSe data
#             taxa enriched in contaminated vs uncontaminated soils
################################################################################
#read in LefSe output
lefse=read.csv("fungal ITS marker/data/07_lefse_output.csv", header=T);

#convert all non-contaminated effect sizes to -values for diverging plots
lefse$LDA[lefse$group=="not_contaminated"]=
  (lefse$LDA[lefse$group=="not_contaminated"])*-1

##set up color palette and LDA tickmark labels for plot
contam_col=c("navajowhite2","#bf812d");
LDAticks=c(-4.5, -3, 0, 3);
lefse$Genus=forcats::fct_rev(factor(lefse$Genus));
lefse$group=factor(lefse$group, levels=c("not_contaminated","contaminated"))

#plot diverging plot
lefsep<-ggplot(lefse, aes(x=Genus, 
                           y=LDA)) +
  geom_point(stat='identity', 
             aes(colour=group),
             size = 2)+
  scale_y_continuous(breaks=LDAticks)+
  scale_colour_manual(values=contam_col,
                      labels=c("Not Contaminated","Contaminated"))+
  coord_flip() +
  theme_bw() + 
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "black", size=0.6)+
  labs(y="LDA Effect Size",
       x="", 
       colour="Enriched in") +
  theme(legend.position="bottom",
        legend.text=element_text(size=9),
        legend.title=element_text(size=9, face="bold"),
        axis.ticks = element_blank(),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=8, face="bold"),
        axis.text.x=element_text(size=10),
        strip.text = element_text(size =8, face="bold"))+
  guides(colour = guide_legend(nrow=2));

plot(lefsep);


################################################################################
#             7. Save LEfSe plot
################################################################################
ggsave(filename="07_lefse_plot.pdf",
       device="pdf",path="fungal ITS marker/figures",
       plot=lefsep,
       width=3.6,
       height=5,
       units="in",
       dpi=500);



