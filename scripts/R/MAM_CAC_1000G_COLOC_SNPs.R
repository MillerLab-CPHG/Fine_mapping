###
## Finemapping SNPs 
library(tidyverse) 

library(coloc)

library(data.table)

library(dplyr)

library(tidyr)

library(plyr)

setwd('/scratch/dw2vr/CAC_1000G/CAC_1000G_coloc/CAC_1000G_MAM_coloc')

MAM<-fread('/scratch/dw2vr/STARNET_dat/STARNET_MAM.bed' )

MAM_SNPS<-MAM$V4

cac_gwas<-fread('/scratch/dw2vr/CAC_1000G/EA_All_FUMA_CAC1000G_Filtered_01072020_eur_af_2')

cac_snps<-cac_gwas$V8

cac_gwas<-cac_gwas%>%
  filter(V8 %in% MAM_SNPS)

MAM<-MAM%>%
  filter(V4 %in% cac_snps )

MAM_filtered <- split( MAM , f = MAM$V8 )



result_list<-llply(MAM_filtered, function(dat){
  dat%>%
    distinct(V4, .keep_all = TRUE)%>%
    drop_na()
}
)

result_list<-result_list[sapply(result_list, function(x) dim(x)[1]) > 5]


##### finemapping snps 

dat_all<-data.frame()

for (i in 1:(length(result_list))){
  
  print(result_list[[i]]$V8[1])
  
  snps<-result_list[[i]]$V4
  
  cac_gwas_genes<-cac_gwas%>%
    filter(V4 %in% snps)%>%
    distinct(V4, .keep_all=TRUE)%>%
    drop_na()
  
  dataset1=
    list(snp=cac_gwas_genes$V4,N=cac_gwas_genes$V13,type="quant", 
         sdY=sd(cac_gwas_genes$V9), MAF=cac_gwas_genes$V19, beta=cac_gwas_genes$V9, 
         varbeta=((cac_gwas_genes$V10^2)))
  
  dat1<-(finemap.abf(dataset1,
                     p1 = 1e-04))%>%
    arrange(desc(SNP.PP))%>%
    drop_na()%>%
    mutate(gene=result_list[[i]]$V8[1])
  
  dat_all<-rbind(dat_all, dat1)

}


write_delim(dat_all, 'MAM_CAC_1000G_PP_SNPs.txt', '\t')



### Genomewide Mapping 

setwd('/scratch/dw2vr/CAC_1000G/CAC_1000G_coloc/CAC_1000G_MAM_coloc')

dat_all_snps<-data.frame()

for (i in 1:(length(result_list))){
  
  print(result_list[[i]]$V8[1])
  
  snps<-as.vector(result_list[[i]]$V4)
  
  cac_gwas_genes<-cac_gwas%>%
    filter(V8 %in% snps)%>%
    distinct(V8, .keep_all = TRUE)%>%
    drop_na(V8)
  
  dat_snps<-coloc.signals(
    dataset1=list(snp=cac_gwas_genes$V8,N=cac_gwas_genes$V9,type="quant", 
                  MAF=cac_gwas_genes$V12, beta=cac_gwas_genes$V5, varbeta=((cac_gwas_genes$V6^2))),
    dataset2=list(snp=result_list[[i]]$V4, type="quant", MAF=result_list[[i]]$V7,N=600, 
                  beta=result_list[[i]]$V9,varbeta=((result_list[[i]]$V10)^2)),
    method = c("single"),
    mode = c("iterative"),
    p1 = 1e-04,
    p2 = 1e-04,
    p12=1e-5,
    maxhits = 3,
    pthr = 1e-06)
  
  dat_snps<-data.frame(dat_snps$summary)%>%
    mutate(Gene=result_list[[i]]$V8[1])
  
  dat_all_snps<-rbind(dat_all_snps, dat_snps)
  
}

write_delim(dat_all_snps, 'MAM_CAC_1000G_ALL_single.txt', '\t' )

mam_dat_sig_snps<-dat_all_snps%>%
  filter( PP.H4.abf > 0.50)

genes<-data.frame(mam_dat_sig_snps$Gene)

write_delim(genes, 'mam_coloc_sig_genes.txt', '\t')

MAM_filt_genes<-MAM%>%
  filter(V8 %in% genes$mam_dat_sig_snps.Gene)

genes_vect<-mam_dat_sig_snps$Gene

genes_vect<-sort(genes_vect)

MAM_filt_split_genes<-split(MAM_filt_genes, f=MAM_filt_genes$V8)

MAM_sig_snps<-llply(MAM_filt_split_genes, function(dat){
  dat$V4})

for(i in names(MAM_sig_snps)){
  write_delim(data.frame(MAM_sig_snps[[i]]), paste0(i,"_snps.txt"), '\t')
}

MAM_filt_genes_pos<-MAM_filt_genes[,c(1,2,3,8,4)]

test<-inner_join(genes, MAM_filt_genes_pos, by=c('mam_dat_sig_snps.Gene'='V8'))

test<-test%>%
  distinct(mam_dat_sig_snps.Gene, .keep_all = TRUE)

View(test)

test_split<-split(test, f=test$V1)

test_split<-llply(test_split, function(dat){
  dat<-dat$mam_dat_sig_snps.Gene})

for(i in names(test_split)){
  write_delim(data.frame(test_split[[i]]), paste0("chr",i,"_snps.txt"), '\t')
}





snp_lists <- list.files(path="/scratch/dw2vr/CAC_1000G/CAC_1000G_coloc/CAC_1000G_MAM_coloc/LD_COR", pattern="snplist$")

ld_cor<-list.files(path="/scratch/dw2vr/CAC_1000G/CAC_1000G_coloc/CAC_1000G_MAM_coloc/LD_COR", pattern="LD_COR.ld$")


for (i in 1:length(snp_lists)){
  assign(paste0(genes_vect[i], "_snplist"),  read_delim(paste0("LD_COR/", snp_lists[i]), '\t',col_names=FALSE)) #read in files using the fread function from the data.table package
}

for (i in 1:length(ld_cor)){
  assign(paste0(genes_vect[i], "_ld_cor"),  read_delim(paste0("LD_COR/", ld_cor[i]), '\t', col_names=FALSE)) #read in files using the fread function from the data.table package
}

sig_snps<-ls(pattern="snplist$")

sig_snps_all<-data.frame()

for (t in sig_snps){
  i=get(t)
  sig_snps_all<-rbind(sig_snps_all, i)
}

filter(MAM, V4 %in% sig_snps_all$X1)


##formatting ld cor matrices

ld_cor_list<-ls(pattern='_ld_cor')

for (t in 1:length(genes_vect)){
  i<-get(paste0(genes_vect[t], "_snplist"))$X1
  dat<-get(paste0(genes_vect[t], "_ld_cor"))
  colnames(dat)<-i
  rownames(dat)<-i
  dat<-as.matrix(dat)
  assign(paste0(genes_vect[t], "_LD_COR.mat"), dat)
}

dat_all<-data.frame()

for (t in 1:length(genes_vect)){
  
  print(genes_vect[t])
  
  i<-get(paste0(genes_vect[t], "_snplist"))$X1
  
  cac_gwas_genes<-cac_gwas%>%
    filter(V8 %in% i)%>%
    distinct(V8, .keep_all = TRUE)%>%
    drop_na(V8)
  
  MAM_genes<-MAM%>%
    filter( V4 %in% i)%>%
    distinct(V4, .keep_all = TRUE)%>%
    drop_na(V4)
  
  dat2<-coloc.signals(
    dataset1=list(snp=cac_gwas_genes$V8,N=cac_gwas_genes$V9,type="quant", 
                  MAF=cac_gwas_genes$V12, beta=cac_gwas_genes$V5, varbeta=((cac_gwas_genes$V6^2))),
    dataset2=list(snp=MAM_genes$V4, type="quant", MAF=MAM_genes$V7, N=600, 
                  beta=MAM_genes$V9,varbeta=(MAM_genes$V10)^2),
    LD=get(paste0(genes_vect[t], "_LD_COR.mat")),
    method = c("cond"),
    mode = c("iterative"),
    p1 = 1e-04,
    p2 = 1e-04,
    p12=1e-5,
    maxhits = 3,
    r2thr=0.01,
    pthr = 1e-05) 
  
  dat2<-data.frame(dat2$summary)%>%
    mutate(Gene=genes_vect[t])
  
  dat_all<-rbind(dat_all, dat2)
  
  dat_all<-dat_all%>%
    distinct()
  
} 

View(dat_all)

write_delim(dat_all, "MAM_1000G_CAC_LD_COR_cond.txt", '\t')

View(dat_all)




