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
A class file is also need. It is just basically a file where you tell RFMIX to which group that individual belongs (target, ref1, ref2, ...)
```
bcftools query -l ref_combination1.recode.vcf | sort > samples_REF.txt
```
Rcode:

```
#I had already data.frames with the individuals belonging to each group (spain, yoruba, america_ref)
sample_file=read.table("samples_REF.txt") 
sample_file[,2]=NA
sample_file[,2]=as.character(sample_file[,2])

for (rowID in 1:nrow(sample_file)){
  if( sample_file[rowID,1]%in% spain$V2) {
    sample_file[rowID,2]="Spain"} 
    else {if(sample_file[rowID,1]%in% yoruba$V2) {
    sample_file[rowID,2]="Yoruba"
    }else{if(sample_file[rowID,1]%in% america_ref$V2) {
      sample_file[rowID,2]="America"} }}}


write.table(sample_file, "sample_ref_file.txt", row.names = F, col.names = F, quote = F)
```


The map file that we created before is not valid for RFMIX, we just need a very simple R code:
```
#Prepare the map file
map<-read.table("dataset.map", sep="")
map=subset(map, map$V1 <23) #only autosomes
map=map[,c(1,4,3)]
map[,3]=map[,3]*100
write.table(map, "ref_map_rfmix_cM.map")
```

### Now we can start!!!
I run it on our server, it took around 24h using 8 cores with more or less 400 individuals. 

```
for chr in {1..22}; do rfmix -f dataset_admix.recode.vcf.gz -r dataset_reference.recode.vcf.gz -m sample_ref_file.txt -g ref_map_rfmix_cM.map -o result_combination1_chr$chr -c 0.2 -r 0.2 -w 0.2 -e 1 --reanalyze-reference  -n 5 -G 11 --chromosome=$chr -n 8 ; done
```




