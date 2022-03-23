# Reduced gene templates for supervised analysis of scale-limited CRISPR-Cas9 fitness screens
![alt text](https://github.com/AleVin1995/Reduced_Templates/blob/main/images/Graphical_abstract.png)

## Description
Many analytical tasks performed on data from genome-wide CRISPR-Cas9 screens are either optionally or necessarily performed in a supervised manner. In this scenario, the screening outcomes observed for large sets of positive/negative controls, i.e. genes that are prior known to be essential/nonessential for cell survival, are adopted as benchmark or template classifiers. The analytical tasks accomplished by supervised methods range from quality control [1-5], to fold-change scaling for inter-screen comparisons and interpretability [6-11], to calling statistical significant essential genes [1,9,11-14].<br/>
However, the size of positive/negative control genes has a significant impact on the scale of the experiment to perform and it proves quite prohibitive for scale-limiteds CRISPR-Cas9 screens (e.g, CRISPR-Cas9 screens using focused libraries, primary cultures, organoids or patient-derived xenografts), where the number of reference genes would become comparable or even larger than that of the genes under investigation.<br/>
Minimal Template Estimator (MinTEs) [15] is a computational framework for assembling gene templates of reduced size allowing supervised analyses of data from scale-limited CRISPR-Cas9 screens, while having a limited, user-defined impact on the overall library size. MinTEs is trained on two large genome-wide pooled CRISPR-Cas9 datasets from Project Score [16] and Project Achilles [10], respectively.<br/>
This repository contains the snakemake pipeline (coming soon..) to derive reduced gene templates from [15].

## Setup
Clone code repository and install dependencies:
```
git clone https://github.com/AleVin1995/Reduced_Templates.git

cd Reduced_Templates

conda env create -f envs/MinTEs.yml
```

## References
[1] Behan, F.M., Iorio, F., Picco, G., Gonçalves, E., Beaver, C.M., Migliardi, G., Santos, R., Rao, Y., Sassi, F., Pinnelli, M., et al. (2019). Prioritization of cancer therapeutic targets using CRISPR–Cas9 screens. Nature 568, 511–516.

[2] Gonçalves, E., Thomas, M., Behan, F.M., Picco, G., Pacini, C., Allen, F., Vinceti, A., Sharma, M., Jackson, D.A., Price, S., et al. (2021). Minimal genome-wide human CRISPR-Cas9 library. Genome Biol. 22, 40.

[3] Hart, T., Brown, K.R., Sircoulomb, F., Rottapel, R., and Moffat, J. (2014). Measuring error rates in genomic perturbation screens: gold standards for human functional genomics. Mol. Syst. Biol. 10, 733.

[4] Koike-Yusa, H., Li, Y., Tan, E.-. P., Velasco-Herrera, M.D.C., and Yusa, K. (2014). Genome-wide recessive genetic screening in mammalian cells with a lentiviral CRISPR-guide RNA library. Nat. Biotechnol. 32.

[5] Tzelepis, K., Koike-Yusa, H., De Braekeleer, E., Li, Y., Metzakopian, E., Dovey, O.M., Mupo, A., Grinkevich, V., Li, M., Mazan, M., et al. (2016). A CRISPR Dropout Screen Identifies Genetic Vulnerabilities and Therapeutic Targets in Acute Myeloid Leukemia. Cell Rep. 17, 1193–1205.

[6] Aguirre, A.J., Meyers, R.M., Weir, B.A., Vazquez, F., Zhang, C.-Z., Ben-David, U., Cook, A., Ha, G., Harrington, W.F., Doshi, M.B., et al. (2016). Genomic Copy Number Dictates a Gene-Independent Cell Response to CRISPR/Cas9 Targeting. Cancer Discov. 6, 914–929.

[7] Munoz, D.M., Cassiani, P.J., Li, L., Billy, E., Korn, J.M., Jones, M.D., Golji, J., Ruddy, D.A., Yu, K., McAllister, G., et al. (2016). CRISPR Screens Provide a Comprehensive Assessment of Cancer Vulnerabilities but Generate False-Positive Hits for Highly Amplified Genomic Regions. Cancer Discov. 6, 900–913.

[8] Tsherniak, A., Vazquez, F., Montgomery, P.G., Weir, B.A., Kryukov, G., Cowley, G.S., Gill, S., Harrington, W.F., Pantel, S., Krill-Burger, J.M., et al. (2017). Defining a Cancer Dependency Map. Cell 170, 564–576.e16.

[9] Dempster, J.M., Pacini, C., Pantel, S., Behan, F.M., Green, T., Krill-Burger, J., Beaver, C.M., Younger, S.T., Zhivich, V., Najgebauer, H., et al. (2019a). Agreement between two large pan-cancer CRISPR-Cas9 gene dependency data sets. Nat. Commun. 10, 5817.

[10] Meyers, R.M., Bryan, J.G., McFarland, J.M., Weir, B.A., Sizemore, A.E., Xu, H., Dharia, N.V., Montgomery, P.G., Cowley, G.S., Pantel, S., et al. (2017). Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells. Nat. Genet. 49, 1779–1784.

[11] Pacini, C., Dempster, J.M., Boyle, I., Gonçalves, E., Najgebauer, H., Karakoc, E., van der Meer, D., Barthorpe, A., Lightfoot, H., Jaaks, P., et al. (2021). Integrated cross-study datasets of genetic dependencies in cancer. Nat. Commun. 12, 1661.

[12] Hart, T., and Moffat, J. (2016). BAGEL: a computational framework for identifying essential genes from pooled library screens. BMC Bioinformatics 17, 164.

[13] Kim, E., and Hart, T. (2021). Improved analysis of CRISPR fitness screens and reduced off-target effects with the BAGEL2 gene essentiality classifier. Genome Med. 13, 2.

[14] Hart, T., Chandrashekhar, M., Aregger, M., Steinhart, Z., Brown, K.R., and MacLeod, G. (2015). High-resolution CRISPR screens reveal fitness genes and genotype-specific cancer liabilities. Cell 163.

[15] Vinceti, Alessandro, Umberto Perron, Lucia Trastulla, and Francesco Iorio. 2022. “Reduced Gene Templates for Supervised Analysis of Scale-Limited CRISPR-Cas9 Fitness Screens.” bioRxiv. https://doi.org/10.1101/2022.02.28.482271.

[16] Dwane, L., Behan, F.M., Gonçalves, E., Lightfoot, H., Yang, W., van der Meer, D., Shepherd, R., Pignatelli, M., Iorio, F., and Garnett, M.J. (2021). Project Score database: a resource for investigating cancer cell dependencies and prioritizing therapeutic targets. Nucleic Acids Res. 49, D1365–D1372.
