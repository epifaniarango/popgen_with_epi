# Local ancestry analysis for masking

For local ancestry analysis, we need first to do the haplotyple estimation or phasing. I will start describing the analysis from a PLINK file where you should have the reference populations and the target populations. You should install these softwares. (I do everything thorugh conda):

- PLINK
- BCFTOOLS
- BEAGLE
- RFMIX

The scripts are not all mine. I picked some codes from https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes and did a few changes.

In my case, I am studing South American populations, sometimes it is not so easy to select a reference panel. If you are in a similar situation I have this short R script based on ADMIXTURE to select the reference panel in these kind of situations.

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

for chromosome in {1..22}; do
	seed=$RANDOM
    beagle gt=splitted_${chromosome}.vcf.gz  map=referenceFromMyDataCentimorgan_Chr${chromosome}.map window=20 seed=$seed out=BeaglePhased${chromosome} nthreads=16
done
```

## 2. Local ancestry inference with RFMIX

The order of the individuals is important when we use RFMIX. I put the admixed individuals first and then the reference panel, this file I called it (samples_order.txt). This file is just a plain text with one column, the name of the individuals on each row on the desired order. We also need to put all the chromose back together.

```
for chr in {1..22}; do bcftools view -S samples_order.txt BeaglePhased$chr.vcf.gz > BeaglePhased$chr.vcf ; done 

for chr in {2..22}; do tabix -p vcf BeaglePhased$chr.vcf.gz ; done
bcftools concat BeaglePhased{1..22}.vcf.gz -Oz -o dataset_phased.vcf.gz
```
