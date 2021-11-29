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
#A) generating table of microbial phyla, order, and genera abundances for LEfSe 
#       LEfSe will identify the taxa that are differentially enriched 
#       between the different sample depths
#B) running LEfSe on the Galaxy platform

source(file="00_background.R"); #load packages & specifications


################################################################################
#             1. Load filtered ASV abundance table, ASV taxonomy table and 
#                         filtered metadata table                 
################################################################################
load("16S rRNA marker/data/01_brownfield_16S_ASVtable_filtered.Rdata");
load("16S rRNA marker/data/01_brownfield_16S_metadata_formatted.Rdata");
load("16S rRNA marker/data/01_brownfield_16S_taxonomy_reduced.Rdata");


################################################################################
#           2. Calculate relative abundances of each phyla
#         and only retain those >0.1% average relative abundance 
################################################################################

#attach taxonomy to the ASV table 
asv_tax=merge(asvdf,mytax[,c(1,2)],by="row.names"); 

#delete the first row, and unecessary Kingdom columnn
asv_tax$Row.names=NULL; asv_tax$Kingdom=NULL;

#aggregate genus-level abundances (e.g. sum abundances of taxa classified
#as the same genus)
tempdf=aggregate(.~Phylum, asv_tax, sum);  
rownames(tempdf)=tempdf$Phylum; tempdf$Phylum=NULL

#calculate relative abundances
tempdf=as.matrix(tempdf);
tempdf=prop.table(tempdf,2);
tempdf=data.frame(tempdf);

#only keep genera with > 0.1% mean abundances across samples
tempdf$rM=rowMeans(tempdf); tempdfA=tempdf[tempdf$rM>0.001,];


################################################################################
#           2. Calculate relative abundances of each order
#         and only retain those >0.1% average relative abundance 
################################################################################

#attach taxonomy to the ASV table 
asv_tax=merge(asvdf,mytax[,c(1,4)],by="row.names"); 

#delete the first row, and unecessary Kingdom columnn
asv_tax$Row.names=NULL; asv_tax$Kingdom=NULL;

#aggregate genus-level abundances (e.g. sum abundances of taxa classified
#as the same genus)
tempdf=aggregate(.~Order, asv_tax, sum);  
rownames(tempdf)=tempdf$Order; tempdf$Order=NULL

#calculate relative abundances
tempdf=as.matrix(tempdf);
tempdf=prop.table(tempdf,2);
tempdf=data.frame(tempdf);

#only keep orders with > 0.1% mean abundances across samples
tempdf$rM=rowMeans(tempdf); tempdfB=tempdf[tempdf$rM>0.001,];


################################################################################
#           2. Calculate relative abundances of each genera
#         and only retain those >0.1% average relative abundance 
################################################################################

#attach taxonomy to the ASV table 
asv_tax=merge(asvdf,mytax[,c(1,6)],by="row.names"); 

#delete the first row, and unecessary Kingdom columnn
asv_tax$Row.names=NULL; asv_tax$Kingdom=NULL;

#aggregate genus-level abundances (e.g. sum abundances of taxa classified
#as the same genus)
tempdf=aggregate(.~Genus, asv_tax, sum);  
rownames(tempdf)=tempdf$Genus; tempdf$Genus=NULL

#calculate relative abundances
tempdf=as.matrix(tempdf);
tempdf=prop.table(tempdf,2);
tempdf=data.frame(tempdf);

#only keep orders with > 0.1% mean abundances across samples
tempdf$rM=rowMeans(tempdf); tempdfC=tempdf[tempdf$rM>0.001,];


################################################################################
#             3. Make data frames for LEfSe                 
################################################################################
#all samples, do not include neighbor lot
breaks=c(0.15, 1.6, 5, 7);
meta$depth2=cut(meta$sample_depth,
                breaks=breaks,
                include.lowest=TRUE,
                labels = c("0.15-1.5","3.05-4.60","6"));

mys2=meta$SampleID[meta$sample_type!="neighbor_lot"]

#add depth info to each sample
dfA=tempdfA[ , (names(tempdfA) %in% mys2)];
colnames(dfA)=meta$depth2[
  match(colnames(dfA),meta$SampleID)];

dfB=tempdfB[ , (names(tempdfB) %in% mys2)];
colnames(dfB)=meta$depth2[
  match(colnames(dfB),meta$SampleID)];

dfC=tempdfC[ , (names(tempdfC) %in% mys2)];
colnames(dfC)=meta$depth2[
  match(colnames(dfC),meta$SampleID)];

#save these dfs for LEfSe
write.table(dfA, "16S rRNA marker/data/07_depth_phylum_lefse.txt", 
            col.names=NA,row.names=TRUE,sep="\t");

write.table(dfB, "16S rRNA marker/data/07_depth_order_lefse.txt", 
            col.names=NA,row.names=TRUE,sep="\t");

write.table(dfC, "16S rRNA marker/data/07_depth_genus_lefse.txt", 
            col.names=NA,row.names=TRUE,sep="\t");


##open files in Excel, and in the very first cell, add the factor
#that you are evaluating: "depth"


################################################################################
#             4. Run LEfSe on the Galaxy platform
#               https://huttenhower.sph.harvard.edu/galaxy/
################################################################################

#Step1: Upload data frames from above by going to Get Data > 
# upload file from your Computer

#Step2: Go to LEfSe > Format Data for LEfSe 
#Select whether the vectors (features and meta-data information) are listed in 
#rows or columns: select ROWS
#Select which row to use as class: use #1 depth
#click EXECUTE

#Step3: Go to LDA Effect Size (LEfSe)
#set the Alpha value for the factorial Kruskal-Wallis test among classes = 0.05
#click EXECUTE

#Step4: compile all tables output by LEFSE into one table for
# supplementary material


