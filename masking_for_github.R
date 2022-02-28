#!/usr/bin/env Rscript
rm(list=ls())
library(tidyverse)


####Processing

#coordinates to know who are the admix americans and the ref americans
all=read.table("Samples_information.txt")

america=c(grep(0, all[,2]))
america_ref=c(grep(3, all[,2]))

orden=read.table("order.txt")
fam=read.table("dataset.fam")

fam=fam[order(match(fam$V2,orden$V1)),]


#subset only the american individuals that are the ones that we need
#for the downstream analysis

fam=subset(fam, fam$V2 %in% c(unique(all[which(all$V2==0),1]),unique(all[which(all$V2==3),1])) )


####Functions 
rowMax <- function(x)apply(x,1, max)

mi_fun <- function(w) {
  rowSums(viterbi1[,c(w-1,w)], na.rm=TRUE)
}




##################SNP possibilities to transfrom 0,1 code to A,C,G,T

acgt=c("A","G","C","T")

combinaciones=as.data.frame(expand.grid(acgt,acgt)) 

#repeated calls 
x=c(combinaciones[,1]==combinaciones[,2])

#remove them
combinaciones=combinaciones[which(x==0, arr.ind = TRUE),]

#we have 12 possible combinations
#and now we need to find which one of those are 

combinaciones$Var3=paste(combinaciones$Var1,combinaciones$Var2,sep = "")



for (y in 1:22){
  
  
  posterior1=read.table(paste("chrom",y,"_rfmix.2.ForwardBackward.txt",sep = ""),sep = "")
  viterbi1=read.table(paste("chrom",y,"_rfmix.2.Viterbi.txt",sep = ""),sep = "")
  
  
  #Process RFMIX output. I need to get rid of the low posterior prob. for each call. (<0.9)
  
  #creamos grupos en los que va cada asignacion de SNPs
  groups <- rep(1:ncol(viterbi1), each=3)
  
  #This function just pick the maximum
  
  #we apply the function before per group of rows. now we have a table with the prob
  #of the assigned ancestry in the viterbi call
  good_pp=sapply(unique(groups), function(i)rowMax(posterior1[,which(groups==i)]))
  good_pp=as.data.frame(good_pp)
  
  good_pp[good_pp<0.9]=0
  
  
  viterbi1[which(good_pp==0, arr.ind = TRUE)] = 0
  
  
  write.table(viterbi1, file=paste("fragment_call_chr",y,".txt", sep = ""), row.names = F, col.names = F, quote = F)
  write.table(viterbi1,"fragment_call_all_chr.txt", col.names=FALSE ,row.names=FALSE, quote=F, sep = "\t", append = T)
  
  
  viterbi1=viterbi1[,1:america[length(america)]]
  
  
  alleles=read.table(paste("chrom",y,"_rfmix.allelesRephased2_sep.txt",sep=""),sep="")
  
  snp_coding=read.table(paste("../phased_chr/chrom1_snp_coding",sep=""), header = T)
  snp_coding$both=paste(snp_coding$REF,snp_coding$ALT,sep = "")
  
  
  for (ref_alt in 1:nrow(combinaciones)) {
    
    coordenadas=which(snp_coding$both ==combinaciones[ref_alt,3], arr.ind = TRUE) 
    alleles[coordenadas,]=ifelse(alleles[coordenadas,]==0,as.character(combinaciones[ref_alt,1]),as.character(combinaciones[ref_alt,2]))
    
  }
  
  alleles_ref=alleles[,america_ref]
  alleles_america=alleles[,america]
  
  
  #now the masking, every position that is not a 3 on the viterbi file should be
  #removed 
  
  alleles_america[which(viterbi1==0, arr.ind = TRUE)] = 0
  alleles_america[which(viterbi1==1, arr.ind = TRUE)] = 0
  alleles_america[which(viterbi1==2, arr.ind = TRUE)] = 0
  
  
  #Now it is MASKED!!! wiii
  
  # we now need to create 2 types of masking 
  #1)diploid masking (only keeping calls that are coming from the american ancestry
  #in both of the sides of the chromosome)
  #2)psuedo-haploid: duplication the haployd and creating the double of individuals 
  
  
  
  #######Pseudo-haploid: we need to give the same treatment to all the samples from the
  #americas. even the non-admixed. 
  
  
  pseudo=cbind(alleles_america, alleles_ref)
  
  #here we duplicate
  pseudo1=pseudo[, rep(seq_len(ncol(pseudo)), each = 2)]
  
  
  
  write.table(pseudo1, file=paste("masking/pseudo_haploid_chr",y,".txt", sep = ""), row.names = F, col.names = F, quote = F)
  write.table(pseudo1,"masking/pseudo_haploid_all_chr.txt", col.names=FALSE ,row.names=FALSE, quote=F, sep = "\t", append = T)
  
  #let's create the ped file
  pseudo_ped=pseudo %>%  t() %>% as.data.frame()
  
  pseudo_ped=pseudo_ped[, rep(seq_len(ncol(pseudo_ped)), each = 2)]
  
  
  fam_pseudo=fam[rep(seq_len(nrow(fam)), each = 2), ]
  
  fam_pseudo[seq(1,nrow(fam_pseudo),2),2]= paste0(fam_pseudo[seq(1,nrow(fam_pseudo),2),2],"_1")
  fam_pseudo[seq(2,nrow(fam_pseudo),2),2]= paste0(fam_pseudo[seq(2,nrow(fam_pseudo),2),2],"_2")
  
  
  pseudo_ped= cbind(fam_pseudo,pseudo_ped)
  
  write.table(pseudo_ped, file=paste("masking/pseudo_haploid_chr",y,".ped", sep = ""), row.names = F, col.names = F, quote = F)
  
  map=read.table( paste("dataset_chr",y,".map",sep=""))
  map_rfmix=read.table(paste("phased_chr/referenceFromMyDataCentimorgan_Chr",y,".map",sep=""))
  
  map=subset(map, map$V4 %in% map_rfmix$V4)
  
  write.table(map, file=paste("masking/pseudo_haploid_chr",y,".map", sep = ""), row.names = F, col.names = F, quote = F)
  
  
  ##############################
  #now we can do the diploid masking
  
  haplotype2=seq(2, ncol(viterbi1),2)
  
  # in ther viterbi file we have these 4 options: 0,1,2,3
  #when we sum per individual the number of each side in one SNP, only when both
  #calls are coming from the american ancestry will be 6 
  
  
  
  
  diploid=as.data.frame(sapply(haplotype2, mi_fun)) 
  
  diploid=diploid[, rep(seq_len(ncol(diploid)), each = 2)]
  
  
  
  #now we remove everything that is not 6 
  alleles_america1=alleles_america
  alleles_america1[which(diploid==0, arr.ind = TRUE)] = 0
  alleles_america1[which(diploid==1, arr.ind = TRUE)] = 0
  alleles_america1[which(diploid==2, arr.ind = TRUE)] = 0
  alleles_america1[which(diploid==3, arr.ind = TRUE)] = 0
  alleles_america1[which(diploid==4, arr.ind = TRUE)] = 0
  alleles_america1[which(diploid==5, arr.ind = TRUE)] = 0
  
  
  # merge with the ref
  
  
  diploid=cbind(alleles_america1, alleles_ref)
  
  write.table(diploid, file=paste("masking/diploid_chr",y,".txt", sep = ""), row.names = F, col.names = F, quote = F)
  
  
  diploid_ped=diploid %>%  t() %>% as.data.frame()
  
  
  diploid_ped_h1=diploid_ped[ seq(1,nrow(diploid_ped),2),]
  diploid_ped_h2=diploid_ped[ seq(2,nrow(diploid_ped),2),]
  
  
  
  neworder <- order(c(2*(seq_along(diploid_ped_h1) - 1) + 1,
                      2*seq_along(diploid_ped_h2)))
  diploid_ped2 = cbind(diploid_ped_h1, diploid_ped_h2)[,neworder]
  
  
  
  diploid_ped2=cbind(fam,diploid_ped2)
  
  write.table(diploid_ped2, file=paste("masking/diploid_chr",y,".ped", sep = ""), row.names = F, col.names = F, quote = F)
  write.table(map, file=paste("masking/diploid_chr",y,".map", sep = ""), row.names = F, col.names = F, quote = F)
}
#####
