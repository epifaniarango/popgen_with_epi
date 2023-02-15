# Dstat plots

We have included the most relevant plots to clarify our D-statistic analysis, computed using ADMIXTOOLS (1) to examine fine-scale population dynamics between ancient and modern samples. For further details on the data used, please refer to the main manuscript. To gain a brief understanding ofthe logic behind D-statistics, I recommend checking out this tutorial (https://www.bodkan.net/admixr/articles/tutorial.html) about the R Package admixr (2).


## F3-statistics
To quantify the genetic affinity in terms of shared genetic drift, we used the outgroup f3-statistic, with the Mbuti as the outgroup, i.e. f3(Mbuti; Pop1, Pop2), using qp3Pop (1). The supplementary table (f3 tab) contains all possible combinations of this statistics with Dataset 3.3. Here, we preent 
the combinations that include one of our newly sequenced Southern Chile populations (Figure 1).

![Rplot13](https://user-images.githubusercontent.com/60963543/209326981-e80e9967-87a8-49c5-be95-18057005d112.png)
***Figure 1:*** Outgroup f3-statistics of the form f3(Mbuti, X, Y) being X: A) Lafkenche B) Pehuenche C) Huilliche-Chiloe and Y all the available populations or ancient samples from Dataset 3.3. Populations on the Y axes are ordered based on an decreasing value of the f3 results. 

## F4-statistics
In our initial analysis, we investigated allele sharing between Mapuche groups with the configuration f4(Mbuti, X; Y, Z), where X is any ancient and modern population of the dataset, and Y and Z are combinations of Lafkenche, Pehuenche, and Huilliche-Chiloé. As expected, our newly sequenced samples display a high level of allele sharing with each other. Among the configurations analyzed, the only sample that approached significance was the one containing the two individuals from CC_Conchali_700BP (Figure 2).

![figure2](https://user-images.githubusercontent.com/60963543/209647696-dad8a61a-8ac2-44f3-81f3-9ad5bac6ccfe.jpg)
***Figure 2:***  f4-statistics to check the genetic continuity between modern Mapuche groups against significant allele sharing with other ancient and modern populations. A) f4(Mbuti, X; Lafkenche, Pehuenche) B) f4(Mbuti, X; Pehuenche, Huilliche-Chiloe). C) f4(Mbuti, X; Lafkenche, Huilliche-Chiloe). The statistics were calculated using all the available SNPs (blue) and only transversions (red). On the X axes, we observe the Zscore, and the vertical lines represent |Z|= 3.3 and 3. 

We further investigated the relationship between Mapuche and Conchalí with an f4 in the following format: f4(Mbuti, Lafkenche/Pehuenche/Huilliche-Chiloe; Conchalí, Ancient sample / Modern sample) and  f4(Mbuti, Ancient sample / Modern samples; Conchalí, Lafkenche/Pehuenche/HuillicheChiloe;) (Figure 3). Our findings suggest that  modern Mapuche populations display a greater affinity to Conchalí than any other contemporary group and nearly all the available ancient samples (Figure 3.A-B). Additionally, We observed some noteworthy relationships between Conchalí_700BP and with Brazil_LapaDoSanto_9,600BP (Figure 3.B) as well as a few Far South individuals (Aonikenk_400BP and Selknam_400BP).

![figure3](https://user-images.githubusercontent.com/60963543/209677886-ac54967c-1079-4409-abf2-7320f2600785.jpeg)
***Figure 3:*** f4-statistics to check the preferential affinity of modern Mapuche groups with Conchalí_700BP. The statistics were calculated using all the available SNPs (blue) and only transversions (red). Ancient samples in positions A) Z, and B) X, and modern populations in positions C)Z and D) X. We observe the Zscore and the vertical lines represent |Z|= 3.3 and 3, and on the Y axes, samples are ordered based on an increasing Zscore.

The ancient samples geographically closest to our Mapuche and Hulliche populations are CC_Conchalí_700BP, CC_LosRieles_5100BP, and CC_LosRieles_12000BP. We designed several f4-statistics of the form: f4(Mbuti, Pehuenche/Lafkenche/Huilliche-Chiloe; CC ancient, FS ancient) (Figure 4). Out results consistently show a close relationship between CC_Conchali_700BP and contemporary Mapuche populations across all the tests conducted. We also observed some affinity with CC_LosRieles_5100BP, but it is not significant and only observed when certain Southern Cone samples are in the configuration. Consistent with the scenario of a population replacement after the Early Holocene (3) in Central Chile, CC_LosRieles_12000BP does not display affinity towards later Southern Cone populations but instead with other Early Holocene samples from more distant regions(Figure 4.C).

![figure4_git](https://user-images.githubusercontent.com/60963543/209687540-186b9604-6f28-452f-ac24-b9ed0c8dbace.jpeg)
***Figure 4:*** f4-statistics to check affinity of ancient Central Chile samples with modern Southern Chile (SC).  The statistics were calculated using all the available SNPs (blue) and only transversions (red). A) Z= Conchali_700BP. B) Z= LosRieles_5100BP. C) Z= LosRieles_1200BP. On the X axes, we observe the Zscore and the vertical lines represent |Z|= 3.3 and 3, and on the Y axes, samples are ordered based on an increasing Zscore.

Nakatsuka et al. (4) reported evidence of gene flow between the Late Holocene ancestry related to CC_Conchali_700BP and Late Holocene Far South samples such Aonikenk, Haus, Yámana, and Selknam, using the statistic f4(Mbuti, CC_Conchali_700BP; Middle Holocene FS, Late Holocene FS). We verified that this migration also involved the genetic ancestors of contemporary Mapuche populations with the same configuration (Figure 5).

![figure5](https://user-images.githubusercontent.com/60963543/209690015-a3f888fa-2d46-4d4f-8637-fb14faa1cdef.jpeg)
***Figure 5:***  f4-statistics confirm gene flow between contemporary Mapuche and Late Holocene Central Chile (Conchalí) populations with Late Holocene Far South. The statistics were calculated using all the available SNPs (blue) and only transversions (red). A) FS_LaArcillosa2_5800BP, B) FS_PuntaSantaAna_7300BP and C) FS_Ayayema_4700BP. On the X axes, we observe the Zscore and the vertical lines represent |Z|= 3.3 and 3, and on the Y axes, samples are ordered based on a decreasing Zscore. Only statistics with more than 30,000 SNPs are represented. 

We also applied the same approach to investigate gene-flow with other areas (Argentinian Pampas and/or Central Andes). However, in both cases, we did not find significant eveidence of gene flow (Figure 6 and 7). The relative abundance of ancient samples in the Andes gives more strength to the lack of signal.

![figure6](https://user-images.githubusercontent.com/60963543/209691911-9066976a-5622-4a09-8b4c-92854dd6b8bf.jpeg)
***Figure 6:*** f4-statistics deniying a strong signal of gene flow between contemporary Mapuche and Late Holocene Central Chile (Conchalí) populations with ancient Argentina. The statistics were calculated using all the available SNPs (blue) and only transversions (red). A) Argentina_ArroyoSeco2_7700BP and B)Argentina_LagunaChica_6800BP. On the X axes, we observe the Zscore and the vertical lines represent |Z|= 3.3 and 3, and on the Y axes, samples are ordered based on a decreasing Zscore. Only statistics with more than 30,000 SNPs are represented.

![figure7](https://user-images.githubusercontent.com/60963543/209785023-3018ad3d-fab5-4849-9d17-7ed2d286088d.jpeg)
***Figure 7:*** f4-statistics deniying a strong signal of gene flow between contemporary Mapuche and Late Holocene Central Chile (Conchalí) populations with the Andes. The statistics were calculated using all the available SNPs (blue) and only transversions (red). A) Peru_Cuncaicha_9000BP, B)Peru_Lauricocha_8600BP, C)Peru_Cuncaicha_4200BP, D)Peru_LaGalgada_4100BP, and E)Peru_Lauricocha_3500BP. On the X axes, we observe the Zscore and the vertical lines represent |Z|= 3.3 and 3, and on the Y axes, samples are ordered based on a decreasing Zscore. Only statistics with more than 30,000 SNPs are represented.


To asses the relationship of the Southern Cone with external populations, we used the configuration f4(Mbuti, X; Central/Southern Chile, Argentina/Far South), where X is any ancient or modern sample not from the Southern Cone.  We found consistently non significant values, suggestiong a single origin of the Southern Cone populations with minor contribution from the rest of South America.  There is just a few combinations with a significant results, mainly combinations that involve middle Holocene Argentina (ArroyoSeco2_7700BP and LagunaChica_6800BP) (Table 1). To further explore the raw statistics go to the  Supplementary Table.


***Table1:*** Proportion of significant and non-significant f4-statistics of the form f4(Mbuti, X; Central/Southern Chile, Argentina/Far South)

| F4-statistics   | Number        | Proportion    |
| --------------- |:-------------:|:-------------:| 
| Significant     | 297           |  1.73         | 
| Non-Significant | 16860         |  98.27        | 
| Total           | 17157         |  100          |


The best fitting qpGraph topolgy that includes contemporary Mapuche populations (Fig. 6C, main text) requires two admxiture edges from a branch related to LapaDoSanto into Conchalí and LosRieles5100BP. By removing this admixture edge, the tree topolgy is not supported(Fig.S7, Supplementary Material). However, we found no consistent significant results with the f4 configuration f4(Mbuti, LapaDoSanto: Conchalí/LosRieles5100BP, Southern Cone), just a tendency which raises uncertainties about this additional 'ancient substrate' connection (Fig. 8).

![figure8](https://user-images.githubusercontent.com/60963543/219059601-9b90fa93-90db-4e4b-9684-bdd1e6677bec.jpeg)
***Figure 8:*** f4(Mbuti, LapaDoSanto; Conchalí/LosRieles5100, Southern Cone)  The statistics were calculated using all the available SNPs (blue) and only transversions (red). On the X axes, we observe the Zscore and the vertical lines represent |Z|= 3.3 and 3, and on the Y axes, samples are ordered based on a decreasing Zscore. Only statistics with more than 30,000 SNPs are represented.

***References***
1. Patterson, N., Moorjani, P., Luo, Y., Mallick, S., Rohland, N., Zhan, Y., Genschoreck, T., Webster, T., and Reich, D. (2012). Ancient Admixture in Human History. Genetics 192, 1065–1093. 10.1534/genetics.112.145037.
2. Petr, M., Vernot, B., and Kelso, J. (2019). admixr—R package for reproducible analyses using ADMIXTOOLS. Bioinformatics 35, 3194–3195. 10.1093/bioinformatics/btz030.
3. Posth, C., Nakatsuka, N., Lazaridis, I., Skoglund, P., Mallick, S., Lamnidis, T.C., Rohland, N., Nägele, K., Adamski, N., Bertolini, E., et al. (2018). Reconstructing the Deep Population History of Central and South America. Cell 175, 1185-1197.e22. 10.1016/j.cell.2018.10.027.
4. Nakatsuka, N., Luisi, P., Motti, J.M.B., Salemme, M., Santiago, F., D’Angelo del Campo, M.D., Vecchi, R.J., Espinosa-Parrilla, Y., Prieto, A., Adamski, N., et al. (2020). Ancient genomes in South Patagonia reveal population movements associated with technological shifts and geography. Nat. Commun. 11, 3868. 10.1038/s41467-020-17656-w.


