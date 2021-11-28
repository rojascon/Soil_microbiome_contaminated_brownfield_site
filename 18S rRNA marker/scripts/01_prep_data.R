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

##CODE FOR: formatting sample metadata, ASV abundance table, and 
#ASV taxonomy table for analysis

source(file="00_background.R"); #load necessary packages and specifications


################################################################################
#             1.  Load metadata, ASV table, and taxonomy table              
################################################################################

asv=read.table("18S rRNA marker/data/00_brownfield_18S_ASVtable.txt", 
               header=T, sep="\t", row.names=1);
tax=read.table("18S rRNA marker/data/00_brownfield_18S_taxonomy.txt", 
               header=T, sep="\t");
meta=read.table("18S rRNA marker/data/00_brownfield_18S_metadata.txt", 
                header=T, sep="\t");


################################################################################
#             2.  Tidy ASV table and ASV taxonomy table              
################################################################################

#remove species column because its mostly NAs
mytax=tax[,c(1:7)];
rownames(mytax)=mytax$X18S_seq_number; mytax$X18S_seq_number=NULL;
save(mytax, file="18S rRNA marker/data/01_brownfield_18S_taxonomy_reduced.Rdata");

#reorder ASV table column names 
asvdf=asv[, order(colnames(asv))];
save(asvdf, file="18S rRNA marker/data/01_brownfield_18S_ASVtable_filtered.Rdata");


################################################################################
#             3.  Format metadata factors             
################################################################################

##this metadata file was modified from the original, raw metadata file for this
#project (see 00_raw_metadata.csv)
#we removed unecessary columns and columns of heavy metals that did not 
#show variation across samples

meta$sample_type=factor(meta$sample_type, levels=c("sample","neighbor_lot"));
meta$contam_type=factor(meta$contam_type);
meta=meta[order(meta$SampleID),];
save(meta, file="18S rRNA marker/data/01_brownfield_18S_metadata_formatted.Rdata");

