
# Dating Admixture

I used ALDER 1.03 to date the admixture date with Spanish colonizers. I installed the programme thorugh conda. 


## 1. Prepare the files

I need to convert my plink files to the EIGENSTRAT format. 

This will depend on your environment installation 
```
conda activate admixtools 
``` 
This is the par_file (par.convert_dataset1). I used to transform the plink files. You need to add your own path and the name of your files (in my case dataset1). 
I always had troubles with the familynames option, so I decided to do it manually with R later.

```
genotypename:    /your_path/dataset1.ped
snpname:         /your_path/dataset1.map
indivname:       /your_path/dataset1.ped
outputformat:    EIGENSTRAT
genotypeoutname: /your_path/dataset1.geno
snpoutname:      /your_path/dataset1.snp
indivoutname:    /your_path/dataset1.ind
familynames:     NO

```
Then you just run the following command
```
convertf -p par.convert_dataset1
```

I add the population name with R and merge all Spanish individuals in a single population:
```
#R
library(tidyverse)
ind=read.table("~/your_path/dataset1.ind")
fam=read.table("~/your_path/dataset1.fam")

ind=cbind(ind[,1:2],fam[,1])
ind$V3=as.character(ind$V3)
vec3=which(grepl("Spanish_",ind$V3))
ind$V3=replace(ind$V3, vec3, "Spanish")

write.table(ind, "~/your_path/dataset1.ind", row.names = F,col.names = F,quote = F)
```

## 2. Create par.files and run ALDER
Here again we have to prepare a par_file (par.ALDER).The best proxy for the American population depends on the target population. For this reason, I tried all the combination on my population list

```

genotypename: /your_path/dataset1.geno
snpname:      /your_path/dataset1.snp
indivname:    /your_path/dataset1.ind
admixpop: Target
refpops:      Spanish;Parent2
chrom:        1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22
checkmap:     NO
raw_outname:  /your_path/LDout_Target_Parent2
``` 

```
#R
setwd("~/your_path/")

ind=read.table("dataset1.ind" )

pops=as.data.frame(unique(ind$V3))

pops=subset(pops, pops$`unique(ind$V3)` != "Spanish")

write.table(pops,"all_pops.txt" , row.names = F, col.names = F, quote = F)


popsvec=as.vector(pops$`unique(ind$V3)`)

for(pop in popsvec){
  
  ko=subset(pops, pops$`unique(ind$V3)` != pop)
  write.table(ko,paste("pops_",pop,".txt", sep = "") , row.names = F, col.names = F, quote = F)
}

```


Now we can create the par.files and run ALDER at the same time. It will take a while, so I put it on the server. 

```
conda activate admixtools #this depends on how you installed the programme and the name of your environments
mypops=$( cat all_pops.txt )

for target in  $mypops; do
	sed s/Target/${target}/g par.MALDER >  par.MALDER_${target}
	parent2=$(cat pops_${target}.txt )
	for parent2 in $parent2; do
		sed s/Parent2/${parent2}/g par.MALDER_${target} >  par.MALDER_${target}_${parent2}
		alder -p par.MALDER_${target}_${parent2} | tee par.MALDER_${target}_${parent2}.logfile
		grep 'success' par.MALDER_${target}_${parent2}.logfile > ${target}_${parent2}_success_out.txt
done
done

```
Even if the run is not succesfull it will output a success_out.txt file empty. To remove all the empty files of the folders:


```
find -type f -empty -delete
```


## 3. Select the best parents 
In my case, it was not difficult as I was mainly intereset in a few populations, so I just picked a partent that with work for all of them. You might try to find a better method if you have a lot of populations.


```
cat *KichwaOrellana_success_out.txt > KichwaOrellanaoutputSuccessMalder.cat.txt
cat *Cree_success_out.txt > CreeoutputSuccessMalder.cat.txt
cat *Chilote_success_out.txt > ChiloteoutputSuccessMalder.cat.txt
cat *Huancas_success_out.txt > HuancasoutputSuccessMalder.cat.txt

mdkir final_combinations
mv *outputSuccessMalder.cat.txt final_combinations/
```

Now you got to R to get the final plots: 

```
library(tidyverse)
library(data.table)
setwd("yout_path/final_combinations/")

my.files <- list.files()
all <- map(my.files, fread)
all2= do.call("rbind",all)
 

#I filter by eye, based on the warnings, pvalue and Z-score

all3=all2[c(5,4,7,8,12,17,20,21,22,23,24,25,26),]

all4 <- all3 %>% separate(V11, c("Generations", "SD"),sep ="/- ") 

all4$Generations = substr(all4$Generations,1,nchar(all4$Generations)-2)

all4$Generations= all4$Generations %>% as.numeric()

all4$SD= all4$SD %>% as.numeric()

ThisYear<-2021  # to have the admixture times in calendar years
generationyeras<-30  # chose your generation time

all4$AgeCalendarYear<-ThisYear-(all4$Generations*generationyeras)
all4$Min<-ThisYear-((all4$Generations + all4$SD)*generationyeras)
all4$Max<-ThisYear-((all4$Generations - all4$SD)*generationyeras)

#This is the population order for my plots. This is additional in case you don't want to have a alphabetical order. It is just a text file with one column and #on each row a population.

order=read.table("~/your_path/uniq_pops.txt")
order=order %>% 
              filter(V1 %in% all4$V4 )  
              
vec1=order(match(all4$V4,order$V1))
all4=all4[ vec1,]


gg =ggplot(all4, aes(x=AgeCalendarYear,y=fct_inorder(V4)))+
  geom_point()+
  geom_errorbarh(aes(xmin=Min,xmax=Max), alpha=0.2)+
  labs(title = "Admixture Times from Spain", xlab("Year"), ylab("Population"))+
  theme_bw() +
  ylab("Population")+
  xlab("Year")

gg
ggsave( "Admixture_time.pdf")

```

Hope it was usefull!
