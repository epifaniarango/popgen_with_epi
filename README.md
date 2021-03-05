**Identity by descent**

Before calculating the IBD fragments, we need to phase. 

## 1. Phasing with BEAGLE
```
plink --bfile  yourfile  --allow-no-sex --recode vcf-iid --alleleACGT --out  dataset
```

Compress and index the file

```
bgzip dataset.vcf
tabix -p vcf dataset.vcf.gz
```
Split the vcf file by chromosomes:

```
for chr in `bcftools view -h dataset.vcf.gz | perl -ne '
 if (/^##contig=<ID=([^,]+)/) { if (length($1)<=2) {
   print "$1\n"
 } }'`; do
  bcftools view -Oz -r $chr dataset.vcf.gz  > splitted_$chr.vcf.gz &
done
```

Move the files to a different folder
```
mkdir BEAGLE_files
mv *.vcf.gz BEAGLE_files
cd BEAGLE_files
```


BEAGLES requires a genetic map (R code)
```

setwd("~/BEAGLE_files/")

mapBeagle<-read.table("dataset.map", sep="")
mapBeagle[,2]<-"." # little tweak to add a dot between values. 
mapBeagle<-mapBeagle[,1:4]
for (i in 1:24){ temp<-mapBeagle[which(mapBeagle[,1]==i),] 
write.table(temp, paste(c("referenceFromMyData_Chr",i,".map"), collapse=""), sep="\t", col.names = F, row.names = F) }

```
Now we have the genetic map in morgans but BEAGLE wants them in cM

```
for chromosome in {1..22}; do
awk '{$3 = $3 * 100; print}' referenceFromMyData_Chr${chromosome}.map > referenceFromMyDataCentimorgan_Chr${chromosome}.map
done
```
Now I activate my conda environment where I installed BEAGLE (I called it BEAGLE). Be aware of the number of cores that you want to use! (nthreads)
```
conda activate BEAGLE

for run in {1..3}; do
  for chromosome in {1..22}; do
	  seed=$RANDOM
      beagle gt=splitted_${chromosome}.vcf.gz  map=referenceFromMyDataCentimorgan_Chr${chromosome}.map window=20 seed=$seed out=RUN${run}_BeaglePhased${chromosome} nthreads=16
done
done
```


## 2. IBD fragments with Refined-IBD

You can download the programme like this (check first if there is a new update)

```
wget "http://faculty.washington.edu/browning/refined-ibd/refined-ibd.12Jul18.a0b.jar"
```
Let's do it!

```
for run in {1..3}; do for chromosome in {1..22}; do java -Xmx100g -jar /path_to_the_folder_where_the_jarfile_is/refined-ibd.17Jan20.102.jar gt=/BEAGLE_files/RUN${run}_Beagle5Phased${chromosome}.vcf.gz map=/BEAGLE_files/referenceFromMyDataCentimorgan_Chr${chromosome}.map window=20 trim=0.3 out=RefinedIBD_RUN${run}.$chromosome nthreads=8; done; done

```

```
gunzip *.ibd.gz

wget https://faculty.washington.edu/browning/refined-ibd/merge-ibd-segments.17Jan20.102.jar

```
With this programme, we can merge the runs. It is a nice option if you don't plan to do any downstream analysis because when you merge we lose information about the phase of each SNP. 


These are the specifications of the program:

usage: cat [in] | java -jar merge-ibd-segments.12Jul18.a0b.jar [vcf] [map] [gap] [discord] > [out]

where
  [in]      = IBD output file from Refined IBD analysis
  [vcf]     = Phased input VCF file from Refined IBD analysis
  [map]     = PLINK format genetic map file with centimorgan (cM) distances
  [gap]     = max length of gap between IBD segments (cM) ####we could try with 1 
  [discord] = max number of genotypes in IBD gap that are inconsistent with IBD
  [out]     = IBD output file with IBD segments after merging


```
for chromosome in {1..22}; do
cat RefinedIBD_RUN1.${chromosome}.ibd RefinedIBD_RUN2.${chromosome}.ibd RefinedIBD_RUN3.${chromosome}.ibd > Chr${chromosome}.123.ibd
cat Chr${chromosome}.123.ibd | java -jar /path_to_the_folder_where_the_jarfile_is/merge-ibd-segments.17Jan20.102.jar /BEAGLE_files/RUN2_Beagle5Phased${chromosome}.vcf.gz  /BEAGLE_files/referenceFromMyDataCentimorgan_Chr${chromosome}.map 0.5 1 > Chr${chromosome}_gap0.5.IBD.Merged
done
```


```
cat *gap0.5.IBD.Merged > all.refinedIBD_gap0.5.Merged
```

## 3. Plotting and analysis
For doing the same approach as in Ioannidis2020, on this script we account for fragments over 7 cM and sometimes 2 pairs of individuals share more than 1 fragment, we need to remove those connections. 
Rscrpipt for processing and plotting. In this case I am only going to show the processing that I did for only fragments over 7cM:
```
setwd("~/your_folder/")
ibd<-read.table("all.refinedIBD0.5.Merged", as.is=T)
colnames(ibd)<-c("firstID","firstHapIndex","secondID","secondHapIndex", "chromosome", "start","end","LOD","length")

ibd=ibd[ibd$length>7, ]

ibd=ibd[!duplicated(ibd[c(1,2)],),]

#infoID.csv is a file with the geographical information of your samples
infoID<-read.csv("infoID.csv",header=T, as.is=T , comment.char = "", fill=T)  #info file each line one individual

minimuminfo<-infoID[,c(2,1)]   # select the columns which have the sample ID name and the corresponding population
colnames(minimuminfo)[1]<-"firstID"
ibdmatch<-merge(ibd,minimuminfo,all.x=TRUE)   # associate the population source for the first sample ID of the couple
colnames(ibdmatch)[10]<-"source1"
colnames(minimuminfo)[1]<-"secondID"
ibdmatch2<-merge(ibdmatch,minimuminfo,all.x=TRUE) # associate the population source for the second sample ID of the couple
colnames(ibdmatch2)[11]<-"source2"

write.table(ibdmatch2, "refinedIBD_merged_withinfo.txt", row.names = F, quote=F)
ibd<-ibdmatch2
ibd=read.table("refinedIBD_merged_withinfo.txt", header = T)

### table sharing per population

#a text file with one column and one row per population in your desired order
popp=read.table("pops.txt")


poporder<-unique(infoID$PopName)

poporder=poporder[match(popp$V1,poporder)]

pops<-table(infoID$PopName)

perpop<-matrix(NA,length(poporder),11)
colnames(perpop)<-c("population","samplesize","numberSharingTot","numbersharingWithin","numberSharingOut","FreqSharingTot","FreqsharingWithin","FreqSharingOut","Mean_lengthsharingWithin","totallenghtsharing","howmanypops")
perpop[,1]<-poporder
perpop[,2]<-pops[poporder]

for (i in 1:nrow(perpop)){
  popp<-poporder[i]
  within<-which(ibd$source1%in%popp & ibd$source2%in%popp)
  tempWithin<-ibd[within,]
  tempTOT<-ibd[union(which(ibd$source1%in%popp),which(ibd$source2%in%popp)),]
  tempOut<-tempTOT[-which(tempTOT$source1 == tempTOT$source2),]
  perpop[i,3]<-nrow(tempTOT)
  perpop[i,4]<-nrow(tempWithin)
  perpop[i,5]<-nrow(tempOut)
  perpop[i,6]<-nrow(tempTOT)/as.numeric(perpop[i,2])
  perpop[i,7]<-nrow(tempWithin)/as.numeric(perpop[i,2])
  perpop[i,8]<-nrow(tempOut)/as.numeric(perpop[i,2])
  perpop[i,9]<-mean(tempTOT$length)
  popvarie<-c(tempTOT$source1,tempTOT$source2)
  perpop[i,10]<-sum(tempTOT$length)
  perpop[i,11]<-length(unique(popvarie))
}

library(reshape2)
perpop2 <- melt(perpop[,c(1,3,4,5)], id.vars='population')
write.table(perpop,"popInfoIBDsharing.txt",sep="\t", row.names = F, quote=F)
perpop<-read.table("popInfoIBDsharing.txt",sep="\t",header=T, as.is=T,comment.char = "", fill=T, quote="")



#------------------------------------------------------------------
### SECTION 1 visualize exchange between populations as number of events
#------------------------------------------------------------------

# create a file to plot a symmetric matrix of exchange
bestmirror<-ibd
bestmirror$source1<-ibd$source2
bestmirror$source2<-ibd$source1
bestdouble<-rbind(ibd,bestmirror)

#------------------------------------------------------------------
# Matrices of exchange between populations
#------------------------------------------------------------------

# matrix with the total number of shared blocks
matrixIBD<-matrix(NA, length(poporder),length(poporder), dimnames=list(poporder, poporder))

for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
    temp<-ibd[union(which(ibd$source1==pop.i),which(ibd$source2==pop.i)),]
    if (pop.i==pop.k){
      matrixIBD[i,k]<- length(which(temp$source1==temp$source2))
    } else {
      tempp<-rbind(temp[which(temp$source1==pop.k),],temp[which(temp$source2==pop.k),])
      matrixIBD[i,k]<-nrow(tempp)
    }
  }
}
write.table(matrixIBD,"matrix_refinedIBD_merge_sharing.txt", sep="\t")



# make a matrix with the average length of blocks
matrixIBDAverageLength<-matrix(NA, length(poporder),length(poporder), dimnames=list(poporder, poporder))

for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
    temp<-ibd[union(which(ibd$source1==pop.i),which(ibd$source2==pop.i)),]
    if (pop.i==pop.k){
      temp2<-temp[(which(temp$source1==temp$source2)),]
      matrixIBDAverageLength[i,k]<- mean(temp2$length)
    } else {
      tempp<-rbind(temp[which(temp$source1==pop.k),],temp[which(temp$source2==pop.k),])
      matrixIBDAverageLength[i,k]<-mean(tempp$length)
    }
  }
}
write.table(matrixIBDAverageLength,"matrix_IBD_averageLength.txt", sep="\t")


# make a matrix with the TOTAL length of blocks
matrixIBDTotLength<-matrix(NA, length(poporder),length(poporder), dimnames=list(poporder, poporder))

for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
    temp<-ibd[union(which(ibd$source1==pop.i),which(ibd$source2==pop.i)),]
    if (pop.i==pop.k){
      temp2<-temp[(which(temp$source1==temp$source2)),]
      matrixIBDTotLength[i,k]<- sum(temp2$length)
    } else {
      tempp<-rbind(temp[which(temp$source1==pop.k),],temp[which(temp$source2==pop.k),])
      matrixIBDTotLength[i,k]<-sum(tempp$length)
    }
  }
}
write.table(matrixIBDTotLength,"matrix_IBD_totalLength.txt", sep="\t")

pops<-table(infoID$PopName)
pops<-pops[which(pops>0)]
pops<-pops[rownames(matrixIBD)]

#adjust for population size
matrixIBDadjustpopsize<-matrixIBD
for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
    matrixIBDadjustpopsize[i,k]<- matrixIBD[i,k]/(pops[i]*pops[k])
  }
}
write.table(matrixIBDadjustpopsize,"matrix_IBDsharingAdjustPopSize.txt", sep="\t")


matrixIBDadjustpopsizelength<-matrixIBDTotLength
for (i in 1:length(poporder)){
  for (k in 1:length(poporder)){
    pop.i=poporder[i]
    pop.k=poporder[k]
    matrixIBDadjustpopsizelength[i,k]<- matrixIBDTotLength[i,k]/(pops[i]*pops[k])
  }
}
write.table(matrixIBDadjustpopsizelength,"matrix_IBDsharingAdjustPopSizelength.txt", sep="\t")

# ---------------
# Melt everything for ggplot

library(reshape2)
library(ggplot2)
meltIBD<-melt(matrixIBD)
colnames(meltIBD)<-c("source1", "source2", "n_sharing")

meltIBDaverage<-melt(matrixIBDAverageLength)
meltIBD$averageLength<-meltIBDaverage$value
meltIBDlength<-melt(matrixIBDTotLength)
meltIBD$totalLength<-meltIBDlength$value
meltIBDadjuxt<-melt(matrixIBDadjustpopsize)
meltIBD$sharingadjust<-meltIBDadjuxt$value
meltIBDlengthadjuxt<-melt(matrixIBDadjustpopsizelength)
meltIBD$lengthadjust<-meltIBDlengthadjuxt$value


meltIBD2<-meltIBD[-(which(meltIBD$source1==meltIBD$source2)),] #exclude same pop sharing
meltIBD2$prob=meltIBD2$n_sharing/(333/2)



# you can filter for more significant pairs, like pairs that share more than once, or more than the median
meltIBD22<-meltIBD2[which(meltIBD2$n_sharing!=0),]
meltIBD22$prob=meltIBD22$n_sharing/(333/2)
meltIBD22<-meltIBD2[which(meltIBD2$n_sharing>1),] # more than once
#meltIBD22<-meltIBD22[which(meltIBD22$sharingadjust>median(meltIBD22$sharingadjust)),] #more than the median

gg<-ggplot(meltIBD22,aes(x=source1, y=source2, fill=sharingadjust, size=lengthadjust))+
  geom_point(shape=21)+
  scale_fill_gradient( low="cyan3",high="darkorchid")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  scale_x_discrete(limits=poporder)+
  scale_y_discrete(limits=poporder)
gg


ggsave("matrixAdjustSize_length_sharing_IBDmerged_.pdf", useDingbats=FALSE, width = 12, height = 10) # Figure 4A

gg<-ggplot(meltIBD22,aes(x=source1, y=source2, fill=sharingadjust, size=n_sharing))+
  geom_point(shape=21)+
  theme_bw() +
  scale_fill_gradient(name="Population size adjustment", low="cyan3",high="darkorchid")+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  scale_x_discrete(limits=poporder)+
  scale_y_discrete(limits=poporder)+
  scale_size(name = "Pairs of sharing")
gg


ggsave("matrixAdjustSize_Nsharing_IBDmerged_.pdf", useDingbats=FALSE, width = 12, height = 10) 

# ------------------------------------------------------------------------
# geographic distances
# ------------------------------------------------------------------------


infoID<-infoID[which(infoID$PopName%in%poporder),]
pops<-table(infoID$PopName)
pops<-pops[poporder]

library(fields)
lista<-(cbind(as.numeric(infoID$lon), as.numeric(infoID$lat)))  #with longitude and latitude coordinates
MatrixGeo<-rdist.earth (lista, miles=FALSE)  #matrix of distance in km between locations
rownames(MatrixGeo)<-infoID$PopName
colnames(MatrixGeo)<-infoID$PopName


geomelt<-melt(MatrixGeo)



# ------------------------------------------------------------------------
### SECTION 2



meltIBD3<-meltIBD2[which(meltIBD2$source1%in%geomelt$Var1),]
meltIBD3<-meltIBD3[which(meltIBD3$source2%in%geomelt$Var1),]
meltIBD3<-meltIBD3[which(meltIBD3$source2!=meltIBD3$source1),]

geomelt<-geomelt[which(geomelt$Var1!=geomelt$Var2),]
geomelt=geomelt[!duplicated(geomelt[c(1,2)]),]

meltIBD3$geodist<-as.character(geomelt$value)
meltIBD3<-meltIBD3[which(meltIBD3$n_sharing!=0),]

meltIBD3$lengthGeo<-meltIBD3$lengthadjust*as.numeric(meltIBD3$geodist)
meltIBD3$sharingGeo<-meltIBD3$sharingadjust*as.numeric(meltIBD3$geodist)

#meltIBD33<-meltIBD3[which(meltIBD3$n_sharing>1),]
#meltIBD33<-meltIBD33[which(meltIBD33$sharingadjust>median(meltIBD33$sharingadjust)),]

library(maps)
library('geosphere')

col.1 <- adjustcolor("mediumorchid", alpha=0.4)
col.2 <- adjustcolor("mediumorchid4", alpha=0.4)
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE)
edge.col <- edge.pal(100)


#Now I want to do the same but without the connections with Spain and also not ploting North and MesoAmerica

spain=c("Spain")
`%out%` <- function(a,b) ! a %in% b
meltIBD3=subset(meltIBD3, meltIBD3$source1 %out% spain & meltIBD3$source2 %out% spain)




pdf("mapSharingNetworkRefinedIBD2_zoom.pdf")



map(database = "world", regions = ".",ylim=c(-60,20), xlim=c(-90,-30), col="grey90", fill=TRUE,  lwd=0.1)
#points(x=infoID$lon, y=infoID$lat, pch=19,  cex=0.5, col="mediumturquoise")

for(i in 1:nrow(meltIBD3))  {
  node1 <- infoID[infoID$PopName == as.character(meltIBD3[i,]$source1),]
  node2 <- infoID[infoID$PopName == as.character(meltIBD3[i,]$source2),]
  
  arc <- gcIntermediate(as.numeric(c(node1[1,]$lon, node1[1,]$lat)), 
                        as.numeric(c(node2[1,]$lon, node2[1,]$lat)), 
                        n=1, addStartEnd=TRUE )
  edge.ind <- round(round(15*meltIBD3[i,]$sharingadjust / max(meltIBD3$sharingadjust)))
  
  lines(arc, col=edge.col[edge.ind], lwd=edge.ind)
}

text(x=as.numeric(infoID$lon), y=as.numeric(infoID$lat), labels=infoID$PopName,  cex=0.1, col="black")

dev.off()






```
