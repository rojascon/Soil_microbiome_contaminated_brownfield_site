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
#A) subsampling samples to 45,000 sequences each using mothur software 
#(to reduce biases due to differences in sampling depth)

#B) calculating 2 metrics of soil microbiome alpha-diversity  
#     Observed Richness, Chao2 Richness

source(file="00_background.R"); #load necessary packages and specifications


################################################################################
#             1. Prepara data for subsampling to 45,000 sequences    
#               this number was chosen because it was the 3rd lowest number of 
#                     sequences found in a biological sample 
################################################################################

#load ASV table and ASV taxonomy
load("16S rRNA marker/data/01_brownfield_16S_ASVtable_filtered.Rdata");
load("16S rRNA marker/data/01_brownfield_16S_taxonomy_reduced.Rdata");

#view # of sequences each sample has
A=colSums(asvdf);
sort(A);

#45,000 appears like a good cutoff (only losing 4 samples)

#make a df of ASV sequences and ASV labels
asvdf=as.data.frame(t(asvdf));
names=data.frame("seqs" = colnames(asvdf), 
                 "names" = sprintf("SP%s",
                                   seq(1:ncol(asvdf))));

#format ASV table for mothur
#step1: replace colnames of ASV table with ASV labels from above
#step2:add the label, Group, and numOTUs columns to ASV table
colnames(asvdf)=names$names[match(names(asvdf),names$seqs)];
tempdf=data.frame("label" = rep(0.03,nrow(asvdf)), 
                  "Group" = row.names(asvdf), 
                  "numOtus" = rep(ncol(asvdf),nrow(asvdf)));
mothur=cbind(tempdf, asvdf);
write.table(mothur, file="16S rRNA marker/data/05_input_mothur.txt", 
            row.names=F, sep="\t");


################################################################################
#           2. Run mothur to rarefy samples to the same sequencing depth 
###############################################################################

#step1: for some reason, must open the above .txt file in Excel, and 
#           save as Tab delimited file (.txt). Overwrite the original. 

#step2: download and open mothur software (Schloss et.al 2009)
# https://github.com/mothur/mothur/releases/tag/v1.46.1

#step3: run the command: 
#           sub.sample(shared=05_input_mothur.txt, size=45000, persample=true)
#step4: rename output file as "05_output_mothur.txt") and place in the data directory

#PASTE Output from mothur here
# X16S_1_2_S2_L001 contains 31296. Eliminating.
# X16S_1_32_S32_L001 contains 28242. Eliminating.
# X16S_1_34_S34_L001 contains 37496. Eliminating.
# X16S_1_52_S52_L001 contains 42529. Eliminating.


################################################################################
#             3. Clean file output from mothur
################################################################################

#read in output from mothur
mothur2=read.table("16S rRNA marker/data/05_output_mothur.txt", sep="\t",header=T);
rownames(mothur2)=mothur2$Group
mothur2=mothur2[,-c(1:3)];
asvdf=as.data.frame(t(mothur2));


################################################################################
#             4. Calculate soil microbiome alpha-diversity 
################################################################################
#make a phyloseq object and calculate Observed Richness
ps<- phyloseq(otu_table(mothur2, taxa_are_rows=FALSE));

alpha=estimate_richness(ps,split = TRUE, 
                        measures = c("Observed"));

alpha$SampleID=rownames(alpha); rownames(alpha)=NULL;

#calculate Chao 2 richness using the fossil package
chao.values<-c();

for(i in 1:ncol(asvdf)){
  #save sampleID
  ss=as.character(colnames(asvdf)[i]);
  #subset to one sample
  df=asvdf[,i];
  #calculate Chao2 richness for that sample
  mval=chao2(df, taxa.row = TRUE);
  chao.values<- rbind(chao.values,mval);
  rownames(chao.values)[i]=ss;
};

#combine data from all samples into one data frame
chao.values=as.data.frame(chao.values);
chao.values$SampleID=rownames(chao.values);
colnames(chao.values)[1]="Chao2";
alpha=merge(alpha, chao.values,by="SampleID");

#append sample metadata
load("16S rRNA marker/data/01_brownfield_16S_metadata_formatted.Rdata");
alpha_meta=merge(alpha,meta, by="SampleID")


################################################################################
#             5. save data frame of alpha diversity metrics 
################################################################################
save(alpha_meta, file="16S rRNA marker/data/05_alpha_diversity.Rdata");


