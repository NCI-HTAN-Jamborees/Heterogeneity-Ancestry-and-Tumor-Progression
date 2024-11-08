# Heterogeneity-Ancestry-and-Tumor-Progression
### Team Lead: Isra Elhussin

### Team Members: Jiaying Lai, Ying Wu, Mohamed Abdalla, Aditi Hazra

### Background and Rationale:
Black patients have an almost 20% higher incidence of colorectal cancer compared to Non-Hispanic White (NHW) patients [ACS Cancer Facts and Figures 2024](https://www.cancer.org/research/cancer-facts-statistics/all-cancer-facts-figures/2024-cancer-facts-figures.html). Compelling evidence suggests there is heterogeneity in tumor microenvironment between Black or African American and European American cancer patients. Why immuno-inflammation pathways and immune cells are upregulated in Black patients remains unknown. Furthermore, the single cell drivers of progression from normal to colorectal adenoma (polyp) to invasive colorectal cancer is understudied in Black patients (Image: SciPro/Getty Images; Reference: [Harvard Health. July 20, 2023](https://www.health.harvard.edu/diseases-and-conditions/they-found-colon-polyps-now-what). However, population descriptors, including the social construct of race, may not adequately capture the complex patterns of continuous human genetic variation. To address this gap, we examined the single cell RNA sequencing (scRNAseq) data by genetic ancestry. [Genetic ancestry](https://nap.nationalacademies.org/read/26902/chapter/1#ii) captures an individual’s family tree by which they inherit DNA from specific ancestors. 
![image](https://github.com/user-attachments/assets/d7e7f0b4-d6bc-4a73-bf41-dfa5d9052b9d)

### Objective: 
We explored the continuum from normal tissue, hyperplasia and premalignant adenomas (colorectal polyps). The primary objective was to elucidate distinct immune signature disparities between normal and premalignant tissues within each racial group (Black versus White) and to ascertain their potential role in tumorigenesis.

### Project Workflow:

Figure 1.  Workflow
![image](https://github.com/user-attachments/assets/23678e30-6083-4564-856e-37fc3daca475)


### Data:
The Vanderbilt University Medical Center (VUMC) scRNA-seq data and spatial transcriptomics data used in this project were accessed using the HTAN Data Portal and Google BigQuery v6. The original data is described by [Chen et al. 2021](https://www.cell.com/cell/fulltext/S0092-8674(21)01381-7?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421013817%3Fshowall%3Dtrue).  We accessed scRNA-seq Discovery (DIS) set and the Validation (VAL) set using [CellxGene](https://cellxgene.cziscience.com/collections/a48f5033-3438-4550-8574-cdff3263fdfd). 


### Data Collection and Processing:
We accessed and extracted data from the HTAN portal and subsequently matched sample IDs with ancestry and clinical annotations. For this project, we included individuals with high ancestry (defined as >70%).  The CellxGene datasets were then filtered based on self-reported race/ancestry, focusing on individuals identified as Black, who were represented solely by normal and premalignant sample types. Additionally, we incorporated demographic variables, including sex (female, male) and age-stratified by the median (≤60 years, >60 years), as well as polyp classification (NL: normal; SSL: sessile serrated adenoma; TA: tubular adenoma). We selected the median age since the age distribution was deviated towards older individuals. 
The single cell data included:

ASC: apoptosis-associated speck-like protein containing a caspase recruitment domain 

ABS: absorptive cells 

TAC: transit amplifying cells 



### Statistical Analysis:
We conducted a comparative analysis of single-cell RNA sequencing data from normal tissue, and premalignant serrated polyp tissue, and tubular adenoma tissue samples from Black and White patients. Wilcoxon rank-sum test was used to compare continuous distributions between two groups.


### Results:
We generated reduced-dimensionality (UMAP) visualization plots and bar plots of single epithelial cells in transcriptome space for normal and premalignant data by ancestry [Results folder].



### Installation:
R Programming (R 4.3.2) and the following packages were used for this project: 

tidyverse

ggplot2

Seurat

cowplot

patchwork

gridExtra

[Seurat v5](https://satijalab.org/seurat/) was used to identify and interpret cellular heterogeneity in normal and premalignant colorectal tissue (Seurat Object_5.0.1).


### Limitations
Diversity and representativeness in the datasets is limited. 


### Future Directions:
We propose to evaluate invasive colorectal cancer scRNAseq data by race and tumor type and compare to the current findings.  We will complete analyses in European individuals. We propose to add covariates analyses including sex and age group for African American and European individuals for normal and premalignant tissue. 


### Acknowledgements
2024 HTAN Data Jamboree

Ino de Brujin

Rowan Beck
