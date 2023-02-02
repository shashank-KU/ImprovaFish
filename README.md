
- [Feed-microbiome-host interactions in Atlantic salmon over life stages](#feed-microbiome-host-interactions-in-atlantic-salmon-over-life-stages)
  - [Overview of sampling](#overview-of-sampling)
  - [Overview of the data generated in this study](#overview-of-the-data-generated-in-this-study)
  - [Getting Started](#getting-started)
    - [Step 1. Package dependencies](#step-1-package-dependencies)
    - [Step 2. Metagenomics](#step-2-metagenomics)
    - [Step 3. Metatranscriptomics](#step-3-metatranscriptomics)
    - [Step 4. Metaproteomics](#step-4-metaproteomics)
    - [Step 5. Metabolomics](#step-5-metabolomics)
    - [Step 6. Meta-metatranscriptomics](#step-6-meta-metatranscriptomics)
    - [Step 7.Omics-Integration](#step-7omics-integration)
- [Bugs](#bugs)

# Feed-microbiome-host interactions in Atlantic salmon over life stages
To meet future food demands, more efficient and sustainable animal production systems are needed. Given the crucial importance of the gut microbiota to animal (host) health and nutrition, selective enhancement of beneficial microbes via prebiotics may be a powerful approach for promoting farmed fish welfare and robustness. In this study, we employed three versions of a beta-mannan prebiotic that were fed to Atlantic salmon and explored the combined responses of both gut microbiota and the host from freshwater to seawater life stages. We have used weighted gene co-expression network analysis (WGCNA) of host RNA-seq and microbial 16S rRNA amplicon sequencing data to identify biological interactions between the gut ecosystem and the host. We observed several microbial and host modules through WGCNA that were significantly correlated with life stage, but not with diet. Microbial diversity was highest in the early life of Atlantic salmon and decreased over time. In particular, Lactobacillus and Paraburkholderia were the dominating genera and showed the highest correlation with host modules. Our findings demonstrate that salmon-microbiota interactions are mainly influenced by life stage, while further research is required to determine whether supplementation of selected prebiotics to diet can be used to modulate the salmon gut microbiota for improving host health and production sustainability.

## Overview of sampling
<img width="960" alt="ImprovAFish_slidespptx" src="https://user-images.githubusercontent.com/30895959/214268569-7962d10e-bc81-4b13-a322-e02df8e406e0.png">


## Overview of the data generated in this study
<img width="960" alt="ImprovAFish" src="https://user-images.githubusercontent.com/30895959/213148498-c9ec83fc-ee0d-4e58-9a79-520b1748db95.png">

## Getting Started
### Step 1. Package dependencies
Installing R packages can be done through various sources such as GitHub, the Comprehensive R Archive Network (CRAN), or by following the official website of the package.

To install a package from GitHub, use the devtools package and the `install_github()` function. For example, to install the "ggplot2" package from GitHub, run the following code:
```
library(devtools)
install_github("ggplot2/ggplot2")
```
To install a package from CRAN, use the `install.packages()` function. For example, to install the "dplyr" package from CRAN, run the following code:
```
install.packages("dplyr")

```
Finally, to install a package from its official website, download the package source code, and use the `install.packages()` function with the local file path as the argument. For example, to install the "reshape2" package from its official website, first download the source code, then run the following code:

```
install.packages("path/to/reshape2_package.tar.gz", repos = NULL, type = "source")

```
It is recommended to regularly update the installed packages to ensure compatibility and to benefit from new features and bug fixes.

```Package dependencies
library("ranacapa")
library("phyloseq")
library("ggplot2")
library("stringr")
library("plyr")
library("reshape2")
library("reshape")
library("dplyr")
library("tidyr")
library("doBy")
library("plyr")
library("microbiome")
library("ggpubr")
library("vegan")
library("tidyverse")
library("magrittr")
library("cowplot")
library("dendextend")
library("WGCNA")
library("metagenomeSeq")
library("decontam")
library("RColorBrewer")
library("ampvis2")
```

### Step 2. Metagenomics
Primers were removed from the raw paired-end FASTQ files generated via MiSeq using “cutadapt”. Further, reads were analyzed by QIIME2 (qiime2-2021.8) pipeline through dada2 to infer the ASVs present and their relative abundances across the samples. For bed dust samples, using read quality scores for the dataset, forward and reverse reads were truncated at 280 bp and 260 bp, followed by trimming the 5′ end till 25 bp for both forward and reverse reads, respectively; other quality parameters used dada2 default values for both 16S rRNA gene sequencing. For 16S rRNA gene sequencing, taxonomy was assigned using a pre-trained Naïve Bayes classifier (Silva database, release 138, 99% ASV) were used.
To ensure that our analyses were not confounded by spurious results, we first analyzed the alpha diversity of negative control samples (including PCR negative, extraction control) that produced sequencing reads (Fig. X). The DNA extraction and other negative controls had significantly lower observed richness than all samples (Kruskal-Wallis test, p = 2.1e-05) for bacterial data. Furthermore, profiles were significantly different for bacterial by different diet group (PERMANOVA for Bray-Curtis, p = 1e-04, R2 = 0.0759). Sequencing contaminants (93 of 7,481 bacterial ASVs) were identified based on the prevalence of ASVs in the negative control and removed using the decontam package (default parameters). We then removed the PCR and sequencing controls before downstream analysis. Data analysis was conducted in R (R Core Team, 2017). Initial preprocessing of the ASV table was conducted using the phyloseq package (v1.38.0). Further filtering was done by removing ASVs without phylum-level classification from 16S rRNA gene sequencing data. Sequencing contaminants were identified and removed using the decontam package. To avoid the bias due to sampling depth, the ASVs table was relative normalized for 16S rRNA gene, and finally we end up with 6,644 ASVs.
All downstream analyses were performed on this normalized ASVs table unless mentioned. We used two alpha diversity indices, i.e., observed richness and Shannon diversity index. Furthermore, beta diversity was calculated using weighted and unweighted UniFrac metric and visualized by principal coordinates analysis (PCoA). Alpha and beta diversity was calculated using phyloseq v1.38.0 and visualized with ggplot2 v3.3.5 in R v4.1.1. Comparison of community richness and diversity was assessed by the Kruskal-Wallis test between all the groups, and comparison between the two groups was done by Wilcoxon test with Benjamini-Hochberg FDR multiple test correction. Significance testing between the groups for beta diversity was assessed using permutational multivariate analysis of variance (PERMANOVA) using the “vegan” package.

For all the commands used for downstream analysis for [metagenomics](Metagenomics/Step2_RNA_Analysis.pdf)



### Step 3. Metatranscriptomics
For all the commands used for downstream analysis for [trascriptomics](MetaTranscriptomics/Step1_Host_Analysis.R)

### Step 4. Metaproteomics
For all the commands used for downstream analysis for [metaproteomics](MetaProteomics/Step5_Proteomics.R)

### Step 5. Metabolomics
For all the commands used for downstream analysis for metabolomics. For [Lipidomics](Metabolomics/Step3_Lipidomics_WGCNA.R), [HILIC](Metabolomics/Step3_Metabolites_1_WGCNA.R), [RP](Metabolomics/Step3_Metabolites_2_WGCNA.R)

### Step 6. Meta-metatranscriptomics
For all the commands used for downstream analysis for [metatranscriptomics](MetaTranscriptomics/Step1_Host_Analysis.R)

### Step 7.Omics-Integration
For all the commands used for downstream analysis for [omics-integration](OmicsIntegration/Step6_Integration.R)


![Omics_Integration-1 2](https://user-images.githubusercontent.com/30895959/214545865-b85babd0-1f6a-4ba5-8fd9-8c039bdf5869.png)


# Bugs
To inform us of any bugs, please open a new issue or send us an email to shashank.gupta@nmbu.no
