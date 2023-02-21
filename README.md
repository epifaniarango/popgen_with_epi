# Local ancestry analysis for masking

This script is designed to mask African and European ancestry from American samples and could be used for other populations with minor adjustments. Our primary focus in this study is to understand the precolonial history of the Southern Cone and Mapuche populations. For this purpose, the African and European ancestries that came in colonial times represent a confounding effect. A possible solution to this issue is a process called masking. The masking process represents a methodological challenge; in the present note, we detail our steps and the different quality checks we performed. Before masking, we performed a PCA test including the admixed American individuals and European and African samples (Fig. 1). On the PC1, we observed a gradient of the American individuals towards the African and mainly towards European individuals (Fig. 1). 


![PCA_continents_noAncient](https://user-images.githubusercontent.com/60963543/189084917-5a7cd5d4-d73f-4c17-8eb7-06513acfc91e.png)

***Figure 1:*** PCA analysis.

We confirmed it by directly comparing the percentage of European ancestry calculated with ADMIXTURE (1) and the position along PC1 and PC2 (Fig. 2). 

![pc1_admix](https://user-images.githubusercontent.com/60963543/189091642-a870ad28-53b5-4143-848e-f99e866994af.png)
***Figure 2:*** Correlation between European ancestry at K=8 and position at PC1 and PC2

By performing a local ancestry analysis (RFMIX in our case, for further details, refer to the method section), we estimated the origins of chromosomal segments in admixed individuals. The local ancestry software performs supervised analysis as we need to provide a set of reference populations with known ancestry to start the analysis. For local ancestry analysis, we need first to do the haplotype estimation or phasing. I will begin describing the analysis from a PLINK file where you should have the reference and target populations. You should install the following software (I do everything through Conda):

- PLINK (2)
- BCFTOOLS (3)
- BEAGLE (4)
- RFMIX v1.5.4 (5)
- VCFTOOLS (6)

(!) I picked some codes from https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes and made a few changes.

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

As we wanted to differentiate the American ancestry fragments from the African and European ones, we had to provide the software with three proxies for each population. The two reference panels for African and European admixture included the Yoruba and Spanish individuals. For the American panel, we selected "unadmixed" American individuals that were established with the following criteria: 

***1)*** not having African or European ancestry (according to K=8 Admixture and the run with the highest likelihood, 0.999 American ancestry)
***2)*** being non-significant for the statistics f4(Unadmixed American Population, Target Individual, Han, San) for Unadmixed American Populations being Karitiana, Mixe, and Xavante, populations, usually chosen to represent unadmixed Native American populations (7, 8).

(!) The design of the reference panel is not provided on this script.

To run this script you need to have 2 files already prepared on your own with a simple script:

#### - Order: a txt file with one column and the name of each individual per row in the designed order. 1) Admixed samples need to go first 2) European individuals 3) African individuals 4) American samples without admixture as the reference panel of the Americas (My file is called order.txt, so you can get an idea of how it should look)

#### - Sample information: the same as the previous file but adding a second column where the admixed Americans will be coded as 0, Europeans as 1, Africans as 2 and the reference Americans as 3. Please don't change the code or the order; my script will not work. (sample_information.txt)

### 2.1 Prepare input files
The plink sample order is usually alphabetical. We need to change this and follow the order that we already pre-defined. After this, we need to change vcf format to the RFMIX format (alleles).
```
for chr in {1..22}; do bcftools view -S order.txt /phased_chr/chrom${chr}_phased.vcf.gz >  phased_chr/chrom${chr}_phased_order.vcf ; done
for chr in {1..22}; do cat chrom${chr}_phased_order.vcf  | grep -v '#' | cut -f10- | tr -d  '\t|' > chrom${chr}.alleles ; done

```
A class file is also needed, where you tell RFMIX to which group that individual belongs (target, ref1, ref2, ...). We will create that file based on the second file I asked you to create. RFMIX asks you to create a code of each HAPLOTYPE, not individual. Rcode:

```
all=read.table("sample_information.txt")
all=all[rep(seq_len(nrow(all)), each = 2), ]
write.table(t(all$V2),"sample_file.txt",col.names = F, row.names = F, quote = F, sep = "\t")
```

We will need this at some point during the downstream analysis. I am not keeping track of which folder I am in (at least not consistently; be careful).
```
for i in `seq 22`; do

 grep -v "^##" chrom${i}_phased_order.vcf | cut -f1-3 > snps_${i}
done
```

### 2.2 Running the program
I ran it on our server, which took around 24h per chromosome with 10 cores. If you change some of the options, it will take less. Everything you need to know to get the program on your computer is here (https://github.com/indraniel/rfmix). You need to run it from the folder where the programme is. 


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
#### 2.2.2 Output comparison with ADMIXTURE
The performance of RFMIX in assigning the Native American component is comparable to that of ADMIXTURE (Fig. 3). Despite this, RFMIX tends to underrepresent the Native American component in North American populations, as demonstrated in Fig. 3. This occurs due to the closer genetic relationship between North American populations and European populations, compared to South American ones.

![comparison](https://user-images.githubusercontent.com/60963543/217549786-5eea35d1-0be2-4de6-9843-71845924874e.jpeg)

***Figure 3:*** Comparison of American ancestry asignation by ADMIXTURE and RFMIX


### 2.3 Processing the output. MASKING

RFMIX creates several output files; one of them contains information on the posterior probability of the ancestry assignation for each SNP. The following script selects calls with higher confidence than 0.9. In the literature, I observed two masking methods, so I decided to evaluate which protocol was better (pseudo(8,9) or diploid masking(10)). In the diploid masking, for a SNP call to be kept, both sides of the chromosome should have been assigned with the same ancestry (in our case, Native American ancestry) (Fig. 4.A). While in the pseudo haploid methods, all calls assigned to the ancestry of interest are kept, as many software do not allow half-calls (a call on one chromosome and a missing call on the second one), the individuals will be pseudohaploidized (Fig. 4.B).


![panel5 (dragged)](https://user-images.githubusercontent.com/60963543/189116108-c383b15e-19c4-4eba-8113-7a3a7d7c4eda.png)
***Figure 4:*** Visual representation of masking methods. A) diploid and B) pseudohaploid masking of one individual.

The code will not work if the names of the files are different. I created the script under R v4.0.3, be aware that things might not work with other versions. You can use a conda environment to select the desired version. The script is above. I create here 2 types of masking: diploid and pseudohaploid.

```
mkdir masking

Rscript masking_for_github.R

```

When it is done (around 1 hour). The text list.txt you can copy it above,  (you can do the same with the diploid files called like this : diploid_chr",y,".ped)

```
for chr in {1..22}; do plink --file pseudo_haploid_chr${chr}  --allow-no-sex --make-bed --out pseudo_haploid_chr${chr} ; done
plink --allow-no-sex --bfile pseudo_haploid_chr1 --merge-list list.txt --make-bed --out pseudo_haploid
```

The plink file pseudo_haploid contains both the masked individuals and the individuals used as the American reference panel in these case.

Thank you Jonas, for your help with this script =)

## 3. Comparing performance of both methods
After following both masking protocols, we compared the performance by calculating the percentage of missing data per individual (Fig.5) and the same f4-statistics mentioned above to compare the performance of the masking (Fig. 6).
![Rplot09](https://user-images.githubusercontent.com/60963543/189121417-7c79fd67-3b62-49d6-ac6e-2b213447621a.png)
***Figure 5:*** Percentage of missing data per individual after applying A) diploid and B) pseudohaploid masking.

![f4_comparison](https://user-images.githubusercontent.com/60963543/189122958-8a3df657-c833-498a-8c74-9af7414364f0.png)
***Figure 6:*** Performance of masking methods A) diploid and B) pseudohaploid masking, measured through the statistic f4(Unadmixed American Population, Target Individual, Han, San) being Unadmixed American Populations being Karitiana, Mixe, and Xavante. 

The performance of both methods is the same as both are based on the same local ancestry run (Fig. 6), but pseudo haploid masking keeps a higher number of SNPs (Fig. 5). We decided to only used the pseudohaploid masking for all the downstream analysis. We also confirmed the performance wih different PCA visualizations (Fig. 7).

![pca_plots](https://user-images.githubusercontent.com/60963543/189308775-164c7121-e466-4cd6-9ccc-6bc3f341c6e9.png)
***Figure 7:*** PCA with masked individuals with “only” Native American ancestry. A) masked American Samples with African and European references. B) masked American samples with European references. C) masked American Samples (Dataset 3.2).

The masking performs well, and we could use the masked samples for several downstream analyses as sensible as D-statistics. I hope it was helpful; email me if you need some help (epifaniarango@gmail.com).

***References***
1. Alexander, D.H., Novembre, J., and Lange, K. (2009). Fast model-based estimation of ancestry in unrelated individuals. Genome Res 19, 1655–1664. 10.1101/gr.094052.109.

2. Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M.A.R., Bender, D., Maller, J., Sklar, P., de Bakker, P.I.W., Daly, M.J., et al. (2007). PLINK: a tool set for whole-genome association and population-based linkage analyses. Am J Hum Genet 81, 559–575. 10.1086/519795.

3. Li H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics (2011) 27(21) 2987-93.
4. Browning, S.R., and Browning, B.L. (2007). Rapid and Accurate Haplotype Phasing and Missing-Data Inference for Whole-Genome Association Studies By Use of Localized Haplotype Clustering. The American Journal of Human Genetics 81, 1084–1097. 10.1086/521987.

5. Maples, B.K., Gravel, S., Kenny, E.E., and Bustamante, C.D. (2013). RFMix: A Discriminative Modeling Approach for Rapid and Robust Local-Ancestry Inference. Am J Hum Genet 93, 278–288. 10.1016/j.ajhg.2013.06.020.
6. Danecek, P., Auton, A., Abecasis, G., Albers, C.A., Banks, E., DePristo, M., Handsaker, R., Lunter, G., Marth, G., Sherry, S.T., et al. (2011). The Variant Call Format and VCFtools. Bioinformatics.
7. Reich, D., Patterson, N., Campbell, D., Tandon, A., Mazieres, S., Ray, N., Parra, M.V., Rojas, W., Duque, C., Mesa, N., et al. (2012). Reconstructing Native American population history. Nature 488, 370–374. 10.1038/nature11258.
8. Capodiferro, M.R., Aram, B., Raveane, A., Rambaldi Migliore, N., Colombo, G., Ongaro, L., Rivera, J., Mendizábal, T., Hernández-Mora, I., Tribaldos, M., et al. (2021). Archaeogenomic distinctiveness of the Isthmo-Colombian area. Cell. 10.1016/j.cell.2021.02.040.
9. Ioannidis, A.G., Blanco-Portillo, J., Sandoval, K., Hagelberg, E., Miquel-Poblete, J.F., Moreno-Mayar, J.V., Rodríguez-Rodríguez, J.E., Quinto-Cortés, C.D., Auckland, K., Parks, T., et al. (2020). Native American gene flow into Polynesia predating Easter Island settlement. Nature 583, 572–577. 10.1038/s41586-020-2487-2.
10. Luisi, P., García, A., Berros, J.M., Motti, J.M.B., Demarchi, D.A., Alfaro, E., Aquilano, E., Argüelles, C., Avena, S., Bailliet, G., et al. (2020). Fine-scale genomic analyses of admixed individuals reveal unrecognized genetic ancestry components in Argentina. PLOS ONE 15, e0233808. 10.1371/journal.pone.0233808.






