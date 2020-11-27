# Local ancestry analysis for masking

For local ancestry analysis, we need first to do the haplotyple estimation or phasing. I will start describing the analysis from a PLINK file where you should have the reference populations and the target populations. You should install these softwares. (I do everything thorugh conda):

- PLINK
- BCFTOOLS
- BEAGLE
- RFMIX
- VCFTOOLS

The scripts are not all mine. I picked some codes from https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes and did a few changes.



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
First we need to put all the chromose back together and prepared 2 vcf files: one with the target individuals and one with the reference panel. For that you will need a plain text file with one column, each row with the name of an individual, one for the target and one for the reference individuals (ref_filter.txt and admix_ind_america.txt). I selected the reference panel based on a previous run of admixture with a simple Rscript. 

```
bcftools concat BeaglePhased{1..22}.vcf.gz -Oz -o dataset_phased.vcf.gz
bgzip dataset_phased.vcf
tabix -p vcf dataset_phased.vcf.gz

vcftools --vcf dataset_phased.vcf --keep ref_filter.txt --recode --out dataset_reference
vcftools --vcf dataset_phased.vcf --keep admix_ind_america.txt --recode --out dataset_admix

```
The reference VCF file needs to be sorted. In my case I have individuals from Spain, Yoruba and America in the reference. I put first the American, then the Spanish and then the Yoruba individuals on a plain text file in the order that I wanted (

The map file that we created before is not valid for RFMIX, we just need a very simple R code:
```
#Prepare the map file
map<-read.table("dataset.map", sep="")
map=subset(map, map$V1 <23) #only autosomes
map=map[,c(1,4,3)]
map[,3]=map[,3]*100
write.table(map, "ref_map_rfmix_cM.map")
```





