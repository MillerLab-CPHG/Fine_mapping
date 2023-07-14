### The purpose of this script is to run coloc, which is described here: https://github.com/chr1swallace/coloc  
library(coloc)
library(data.table)
library(dplyr)
library(tidyverse)
setwd('/path/to/coloc/folder')

qtl <- fread("/path/to/qtl/results/eQTL_results.txt")

gwas <- fread("/path/to/gwas/results/GWAS_sumstats.txt.gz")
setnames(gwas, c("snp", "p_value", "ref", "alt", "REF_AF", "beta", "standard_error"))
gwas <- gwas[snp %in% uva_qtl$snp & REF_AF>0 & REF_AF<1,]

        gene_list <- qtl[!duplicated(gene_name), c("gene_name", "ensg")]
        genes <- qtl[, c("gene_name","snp","beta","N","ALT_AF","pos","standard_error")]

dat_all<-data.table()

for (t in 1:length(gene_list$gene_name))        {

        qtl_dat <- genes[gene_name==gene_list[t, 1]]
        qtl_dat <- qtl_dat[!duplicated(snp)]

        gwas_dat <- gwas[snp %in% qtl_dat$snp,]
        gwas_dat <- gwas_dat[!duplicated(snp)]

        if(!(is.null(gwas_dat))) {
        print(gene_list[t, gene_name])
        }
        if(nrow(gwas_dat)==0)   {
        print("Next gene"); next
        }

        qtl_dat <- qtl_dat[snp %in% gwas_dat$snp,]
        setorder(gwas_dat, snp); setorder(qtl_dat, snp)
        gwas_dat[, position:=qtl_dat$pos]
        dat2<-coloc.signals(
### s refers to case proportions for case-control GWAS studies
        dataset1=list(snp=gwas_dat$snp, type="cc", MAF=gwas_dat$REF_AF, beta=gwas_dat$beta, varbeta=gwas_dat$standard_error^2, position=gwas_dat$position, s=0.138),
### N refers to total number of GWAS participants for quantitative traits, only use one or the other
#       dataset1=list(snp=gwas_dat$snp, type="quant", MAF=gwas_dat$REF_AF, N=500000, beta=gwas_dat$beta, varbeta=gwas_dat$standard_error^2, position=gwas_dat$position),
        dataset2=list(snp=qtl_dat$snp, type="quant", MAF=qtl_dat$ALT_AF, N=qtl_dat$N, beta=qtl_dat$beta, varbeta=qtl_dat$standard_error^2, position=qtl_dat$pos),
        method = c("single"), mode = c("iterative"),
        p1 = 1e-04, p2 = 1e-04, p12= 1e-05,
        maxhits = 3)

  dat2<-as.data.table(dat2$summary)
  dat2[, gene_name:=gene_list[t, 1]]
  dat_all<-rbind(dat_all, dat2)

}

write.table(dat_all, "Study_name_coloc_date.txt", sep="\t", quote=F, row.names=F)

