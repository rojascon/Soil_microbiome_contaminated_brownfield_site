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

##CODE FOR: identifying collinear microbial community predictors, specifically
## pollutant values that covary with one another
# only pollutants that do not covary with one another will be included in
# statistical analyses

source(file="00_background.R"); #load necessary packages and specifications


################################################################################
#             1.  Load tidyd metadata 
################################################################################

load("16S rRNA marker/data/01_brownfield_16S_metadata_formatted.Rdata");


################################################################################
#             2. Identify heavy metals, hydrocarbons, PAHs that covary 
#                        (e.g. are collinear)
################################################################################

#make a matrix of only the heavy metals, hydrocarbon, PAHs 
#[and remove rows with NAs]
poll=meta[,6:(ncol(meta)-1)];
poll=poll[complete.cases(poll),];

#compute correlation matrix 
#coeff= correlation coefficient aka the strength of the correlation, 
#pvals=pvalues for that correlation
res2 <- rcorr(as.matrix(poll));
coeff=data.frame(res2[[1]]);
pvals=data.frame(res2[[3]]);
View(coeff); View(pvals);

##remove weak and nonsignificant correlations and View resulting matrix
coeff2=coeff; 
coeff2[pvals>0.01]=NA;  #remove nonsignificant correlations alpha>0.01
coeff2[coeff2<0.5]=NA;  #remove weak correlations r < 0.05
View(coeff2);

##based on matrix above, there are several pollutants that covary and will be
#excluded from all statistical analysises


################################################################################
#             3. Remove the collinear heavy metals, hydrocarbons, PAHs 
#                         and rerun correlation matrix
################################################################################

##only keeping pollutants that do not covary
poll2=poll[,c(1,3,4,5,7,11,19,20)]; 
res2 <- rcorr(as.matrix(poll2));
coeff2b=data.frame(res2[[1]]);
pvals2b=data.frame(res2[[3]]);

#remove nonsigificant correlations (alpha>0.01)
coeff2b[pvals2b>0.01]=NA; 
View(coeff2b);   

#the reduced data frame from above is good- most pollutants are not correlated
## with one another
#TPH_do is correlated with lead, so we will remove
#cadmium is correlaed with cobalt, so we will remove
#benzo_a_anthracene is correlated with arsenic, so we will remove
coeff2b=coeff2b[c(1,3,4,5,8), c(1,3,4,5,8)]
print("our final microbial community predictors are:");
print(colnames(coeff2b));

