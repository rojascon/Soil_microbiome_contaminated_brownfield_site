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
#             1.  Load metadata, ASV table, and ASV taxonomy table              
################################################################################

asv=read.table("fungal ITS marker/data/00_brownfield_fITS_ASVtable.txt", 
               header=T, sep="\t", row.names=1);
tax=read.table("fungal ITS marker/data/00_brownfield_fITS_taxonomy.txt", 
               header=T, sep="\t");
meta=read.table("fungal ITS marker/data/00_brownfield_fITS_metadata.txt", 
                header=T, sep="\t");


################################################################################
#             2.  Tidy ASV table and taxonomy table              
################################################################################

#make ASV name the rownames
mytax=tax;
rownames(mytax)=mytax$ITS_seq_number; mytax$ITS_seq_number=NULL;
save(mytax, file=
       "fungal ITS marker/data/01_brownfield_fITS_taxonomy_reduced.Rdata");

#remove ASVs that are unclassified from 16S ASV table
#these ASVs were already removed from the taxonomy table
asvdf=asv[rownames(asv) %in% rownames(mytax),];
asvdf=asvdf[, order(colnames(asvdf))];

#remove samples that did not yield many sequences (<200 sequences)
colSums(asvdf);
bad=c(12,13,14,16,17,18,19,2,21,22,24,26,27,29,32,34,36,
      37,38,39,4,41,42,43,44,45,46,49,50,52,53,54,7,8,9);

meta=meta[!meta$sample_tube %in% bad,];
asvdf=asvdf[, colnames(asvdf) %in% meta$SampleID];
save(asvdf, file="fungal ITS marker/data/01_brownfield_fITS_ASVtable_filtered.Rdata");


################################################################################
#             3.  Format metadata factors             
################################################################################

##this metadata file was modified from the original, raw metadata file for this
#project (see 00_raw_metadata.csv)
#we removed unecessary columns and columns of heavy metals that did not 
#show variation across samples
#above, we also removed the samples that did not yield many sequences

meta$sample_type=factor(meta$sample_type, levels=c("sample","neighbor_lot"));
meta$contam_type=factor(meta$contam_type);
meta=meta[order(meta$SampleID),];
save(meta, file="fungal ITS marker/data/01_brownfield_fITS_metadata_formatted.Rdata");

