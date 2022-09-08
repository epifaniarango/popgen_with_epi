# Local ancestry analysis for masking

This script is designed to mask African and European ancestry from American samples. It could be used for other populations with small adjutments. Our primary focus in this study is to understand the precolonial history of the Southern Cone and Mapuche populations. For this purpose, the African and European ancestries that came in colonial times represent a confounding effect. A possible solution to this issue is a process called masking. The masking process represents a methodological challenge; in the present note, we detail our steps and the different quality checks we performed. Before masking, we performed a PCA test including the admixed American individuals and European and African samples (Fig. 1). 

![PCA_continents_noAncient](https://user-images.githubusercontent.com/60963543/189084917-5a7cd5d4-d73f-4c17-8eb7-06513acfc91e.png)

***Figure 1:*** PCA analysis.

The position along PC1 and PC2 is proportional to the precentage of European ancestry calculated with ADMIXTURE (1)(Fig. 2).

![pc1_admix](https://user-images.githubusercontent.com/60963543/189091642-a870ad28-53b5-4143-848e-f99e866994af.png)
***Figure 2:*** Correlation between European ancestry at K=8 and position at PC1 and PC2

By performing a local ancestry analysis (RFMIX in our case, for further details, refer to the method section), we estimated the origins of chromosomal segments in admixed individuals. The local ancestry software performs supervised analysis as we need to provide a set of reference populations with known ancestry to start the analysis. For local ancestry analysis, we need first to do the haplotyple estimation or phasing. I will start describing the analysis from a PLINK file where you should have the reference populations and the target populations. You should install the following softwares (I do everything thorugh conda):

- PLINK
- BCFTOOLS
- BEAGLE
- RFMIX v1.5.4
- VCFTOOLS

(!) I picked some codes from https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes and did a few changes.

## 1. Phasing with BEAGLE
```
plink --bfile  yourfile  --allow-no-sex --recode vcf-iid --alleleACGT --out  dataset

#We will use this later
for chr in {1..22}; do \
plink --bfile yourfile --allow-no-sex --chr $chr --recode --out dataset_chr${chr}; \
done
```
RFMIX v1.5.4 does not accept invariant sites. 
```
bcftools view -v snps -o dataset_noinvariant.vcf dataset.vcf
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
mkdir phased_chr
mv *.vcf.gz  phased_chr
cd  phased_chr
```


BEAGLES requires a genetic map (R code)
```

setwd("~/phased_chr/")

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
    beagle gt=splitted_${chromosome}.vcf.gz  map=referenceFromMyDataCentimorgan_Chr${chromosome}.map window=20 seed=$seed out=chrom${chromosome}_phased nthreads=16
done

for chr in {1..22}; do tabix -p vcf chrom${chr}_phased.vcf.gz ; done

```
## 2. Local ancestry inference with RFMIX v1.5.4 

As we wanted to differentiate the American ancestry fragments from the African and European ones, we had to provide the software with three proxies for each population.The two reference panels for African and European admixture included the Yoruba and Spanish individuals. For the American panel, we selected “unadmixed” American individuals that were established with the following criteria: 

***1)*** not having African or European ancestry (according to K=8 Admixture and the run with the highest likelihood, 0.999 American ancestry)
***2)*** being non-significant for the statistics f4(Unadmixed American Population, Target Individual, Han, San) for Unadmixed American Populations being Karitiana, Mixe, and Xavante, populations, usually chosen to represent unadmixed Native American populations (2, 3).

(!) The design of the reference panel is not provided on this script.

To run this script you need to have 2 files already prepared on your own with a simple script:

#### - Order: a txt file with one column and the name of each individual per row in the desigened order. 1) Admixed samples need to go first 2) European individuals 3) African individuals 4) American samples without admixture as the reference panel of the Americas (My file is called order.txt, so you can get an idea of how it should look)

#### - Sample information: the same as the previous file but adding a second column where the admixed Americans will be coded as 0, Europeans as 1, Africans as 2 and the reference Americans as 3. Don't change the code, or the order, my script will not work.(sample_information.txt)

### 2.1 Prepare input files
The plink samples order follow usually an alphabethical order. We need to change this and follow the order that we already pre-defined. After this we need to change vcf format to the RFMIX format (alleles).
```
for chr in {1..22}; do bcftools view -S order.txt /phased_chr/chrom${chr}_phased.vcf.gz >  phased_chr/chrom${chr}_phased_order.vcf ; done
for chr in {1..22}; do cat chrom${chr}_phased_order.vcf  | grep -v '#' | cut -f10- | tr -d  '\t|' > chrom${chr}.alleles ; done

```
A class file is also need. It is just basically a file where you tell RFMIX to which group that individual belongs (target, ref1, ref2, ...). We will create that file based on the second file I asked you to create. RFMIX asks you to create a code of each HAPLOTYPE, not invidiual. Rcode:

```
all=read.table("sample_information.txt")
all=all[rep(seq_len(nrow(all)), each = 2), ]
write.table(t(all$V2),"sample_file.txt",col.names = F, row.names = F, quote = F, sep = "\t")
```

We will need this at some point during the downstream analysis. I am not keeeping track on which folder I am (at least not always, be carefull).
```
for i in `seq 22`; do

 grep -v "^##" chrom${i}_phased_order.vcf | cut -f1-3 > snps_${i}
done
```

### 2.2 Running the program
I ran it on our server, it took around 24h per chr with 10 cores. If you change some of the options, it would take less. Everything that you need to know to get the programm in your computer is here (https://github.com/indraniel/rfmix). You need to run it from the folder where the programme is. 

```
cd /your_path/RFMix_v1.5.4 
for chr in {1..22}; do python2 RunRFMix.py PopPhased -G 11 --num-threads 10 -n 5 --forward-backward --use-reference-panels-in-EM -e 2 -w 0.2 /your_path/phased_chr/chrom${chr}.alleles /your_path/sample_file.txt /your_path/phased_chr/referenceFromMyDataCentimorgan_Chr${chr}.map -o /your_path/chrom${chr}_rfmix ; done
```
When the analysis is done, run these simple commands. We will need them on the downstream analysis:

```
for chr in {1..22}; do cat chrom${chr}_rfmix.allelesRephased2.txt |  sed 's/./& /g' > chrom${chr}_rfmix.allelesRephased2_sep.txt ; done

for chr in {1..22}; do grep -v "^##" /phased_chr/chrom${chr}_phased_order.vcf | cut -f4-5 > /phased_chr/chrom${chr}_snp_coding ; done

mkdir splitted_plink

for chr in {1..22}; do vcftools --vcf  /phased_chr/chrom${chr}_phased_order.vcf --plink --out splitted_plink/splitted_${chr}; done
```

### 2.3 Processing the output. MASKING

RFMIX creates several output files, one of them contains the information of the posterior probability of the ancestry assignation for each SNP. The following script selects calls with a higher confidence than 0.9. In literature, I observed two masking methods, so I decided to evaluate which protocol was better (pseudo(3,4) or diploid masking(5)). In the diploid masking, for a SNP call to be kept, both sides of the chromosome should have been assigned with the same ancestry (in our case, Native American ancestry) (Fig. 3).While in the pseudo haploid methods, all calls assigned to the ancestry of interest are kept, as many software do not allow half-calls (a call on one chromosome and a missing call on the second one), the individuals will be pseudohaploid (Fig. 3.B).


![panel5 (dragged)](https://user-images.githubusercontent.com/60963543/189116108-c383b15e-19c4-4eba-8113-7a3a7d7c4eda.png)
***Figure 3:*** Visual representation of masking methods. A) diploid and B) pseudohaploid masking of one individual.

The code will not work if the names of the files are different. I created the script under R v4.0.3, be aware that things might not work with other version. You can use a conda environment to select the desired version. The script is above. I create here 2 types of masking: diploid and pseudohaploid.
```
mkdir masking

Rscript masking_for_github.R

```

When it is done (around 1 hour). The text list.txt you can copy it above,  (you can do the same with the diploid files called like this : diploid_chr",y,".ped)

```
for chr in {1..22}; do plink --file pseudo_haploid_chr${chr}  --allow-no-sex --make-bed --out pseudo_haploid_chr${chr} ; done
plink --allow-no-sex --bfile pseudo_haploid_chr1 --merge-list list.txt --make-bed --out pseudo_haploid
```

The plink file pseudo_haploid contains both the masked individuals and the individuals used as a the American reference panel, in these case.

Thank you Jonas for your help with this script =)

## 3. Comparing performance of both methods
After following both masking protocols, we compared the performance by calculating the percentage of missing data per individual (Fig. 4), and the same f4-statistics mentioned above to compare the performance (Fig. 5).
![Rplot09](https://user-images.githubusercontent.com/60963543/189121417-7c79fd67-3b62-49d6-ac6e-2b213447621a.png)
***Figure 4:*** Percentage of missing data per individual after applying A) diploid and B) pseudohaploid masking.

![f4_comparison](https://user-images.githubusercontent.com/60963543/189122958-8a3df657-c833-498a-8c74-9af7414364f0.png)
***Figure 5:*** Performance of masking methods A) diploid and B) pseudohaploid masking, measured through the statistic f4(Unadmixed American Population, Target Individual, Han, San) being Unadmixed American Populations being Karitiana, Mixe, and Xavante. 

The performance of both methods is the same as both are based on the same local ancestry run (Fig. 5), but pseudo haploid masking keeps a higher number of SNPs (Fig. 4). We decided to only used the pseudohaploid masking for all the downstream analysis. We also confirmed the performance wih different PCA visualizations (Fig. 6).

![image](https://user-images.githubusercontent.com/60963543/189123549-78030dc9-7d30-463d-be08-2fa139043e9b.png)
***Figure 6:*** PCA with masked individuals with “only” Native American ancestry. A) masked American Samples with African and European references. B) masked American samples with European references. C) masked American Samples (Dataset 3.2).![image](https://user-images.githubusercontent.com/60963543/189123567-c35dae9e-cf78-4fe4-996e-8d33cb42f0d8.png)



