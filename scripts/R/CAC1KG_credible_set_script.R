#### The purpose of this script is to generate credible sets for GWAS association signals
#### Last updated 4/17/2023 by Chani Hodonsky
#### Adapted from Doris Wong's script which was adapted from Paul DeVries's original script.

#### NOTE: ALL variants must have a p-value (i.e., remove NA values prior to running this script) 
#### NOTE: This script is designed to be run by locus, so individual files containing all SNPs per locus are required; as is, only one 95% credible set will be generated across each input file.

library(tidyverse)

snp_files<-list.files("/project/cphg-millerlab/chani/scripts/CAC_1KG", pattern="_CS_input_filtered")

credible_set <- function(assoc){
  
  
  ## calculate Bayes Factors based on association P-values
  assoc$bf  <- exp((qnorm(assoc$"P"/2))^2/2)
  ## calculate Posterior Prob in the region based on Bayes Factors
  assoc$Posterior_Prob  <- assoc$bf/sum(assoc$bf)
  
  ## build 95% credible set based on Posterior Prob
  assoc  <- assoc[order(assoc$Posterior_Prob, decreasing=T),]
  assoc$cpp <- 0
  for(i in 1:dim(assoc)[1]){
    assoc$cpp[i]  <- sum(assoc[1:i,]$Posterior_Prob)
  }
  
  ## identify the first line in which the cumulative Posterior Prob is over 95%
  index <- min(which(assoc$cpp>=0.95))
  cred_set <- assoc[1:index,]
  cred_set$bf  <- NULL
  cred_set$cpp <- NULL
  cred_set$Posterior_Prob <- round(cred_set$Posterior_Prob, 3)
  
  return(cred_set)
  print(head(cred_set))  
}

cac_snps<-lapply(snp_files, function(x){ 
read_delim(x, " ") %>% filter(!is.na(P))
})

names(cac_snps)<-snp_files

cac_snps_cred_set<-lapply(cac_snps, function(x){ 
credible_set(x)
})

cac_snps_cred_set_all<-bind_rows(cac_snps_cred_set, .id = "column_label")

write_delim(cac_snps_cred_set_all, "/project/cphg-millerlab/chani/CAC_1KG/CAC_EUR_AFR_cred_set_all_loci_25kb_vs_100kb_Apr11.txt", '\t', append=T)

