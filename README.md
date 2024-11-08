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
The Vanderbilt University Medical Center (VUMC) scRNA-seq data and spatial transcriptomics data used in this project were accessed using the HTAN Data Portal and Google BigQuery v6. The original data is described by [Chen et al. 2021](https://www.cell.com/cell/fulltext/S0092-8674(21)01381-7?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421013817%3Fshowall%3Dtrue).  We accessed scRNA-seq Discovery (DIS) set, the Validation (VAL) set, non-epithelial set using [CellxGene](https://cellxgene.cziscience.com/collections/a48f5033-3438-4550-8574-cdff3263fdfd). 


### Data Collection and Processing:
We accessed and extracted data from the HTAN portal and subsequently matched sample IDs with ancestry and clinical annotations. For this project, we included individuals with high ancestry (defined as >70%).  The CellxGene datasets were then filtered based on self-reported race/ancestry, focusing on individuals identified as Black, who were represented solely by normal and premalignant sample types. Additionally, we incorporated demographic variables, including sex (female, male) and age-stratified by the median (≤60 years, >60 years), as well as polyp classification (NL: normal; SSL: sessile serrated adenoma; TA: tubular adenoma). We selected the median age since the age distribution was deviated towards older individuals. We examined the epithelial cells which comprises of different mature cell types, which are separated into absorptive cells (enterocytes) and secretory (goblet, enteroendocrine and tuft cell) lineages. 
The VUMC CellxGene single epithelial cell dataset included:

ASC: Adenosquamous carcinoma cells have neoplastic potential, it arises from glandular cells. The coexistance of squamous compents is associated with poorer prognosis. Chen et al: ASCs resembled colonic stem and progenitor cells,

ABS: absorptive cells function in the uptake of metabolites and absorbtion of nutrients

TAC: transit amplifying cells are a unique type of progenitor cells crticial for human tissue regeneration and cell turnover.



### Statistical Analysis:
We conducted a comparative analysis of single-cell RNA sequencing data from normal tissue, and premalignant serrated polyp tissue, and tubular adenoma tissue samples from Black and White patients. Wilcoxon rank-sum test was used to compare continuous distributions between two groups.

### Software requirements:  
1. The R used in the analysis is **R4.3.2** or **R4.4.2**    
2. Package used in the analysis is as following.  
  Seurat  
  tidyverse  
  ggplot2  
  Seurat  
  cowplot  
  patchwork  
  gridExtra

[Seurat v5](https://satijalab.org/seurat/) was used to identify and interpret cellular heterogeneity in normal and premalignant colorectal tissue (Seurat Object_5.0.1).

3. Session info: [*Session info1*](https://github.com/NCI-HTAN-Jamborees/Heterogeneity-Ancestry-and-Tumor-Progression/blob/main/scripts/session_info_val_dis_set.txt) and [*Session info2*](https://github.com/NCI-HTAN-Jamborees/Heterogeneity-Ancestry-and-Tumor-Progression/blob/main/scripts/session_info_non-epithelial_set_analysis.txt) 

### Results:
We generated reduced-dimensionality (UMAP) visualization plots and bar plots of single epithelial cells in transcriptome space for normal and premalignant data by ancestry [Results folder]. The following results are for the VAL (Validation) set. Overall, Black individuals had the highest prevalence of tubular adenomas, which is associated with increased risk of colorectal cancer. 

Figure 1: UMAP of epithelial data by ancestry and tissue type
![image](https://github.com/user-attachments/assets/9022897f-04b0-425e-a989-f34d54e4da71)

Figure 2: UMAP of epithelial data by cell-type
![image](https://github.com/user-attachments/assets/1669e944-df0b-4ae0-8c6d-254ad72f5d9c)

Figure 3: UMAP of epithelial data by polyp type

![image](https://github.com/user-attachments/assets/44310fee-18bb-4a82-b87c-d8837443ec54)

Figure 4: Dotplot epithelial data displaying top 10 upregulated genes in European premalignant compared to European normal tissue.

![image](https://github.com/user-attachments/assets/823dc930-871c-4b25-83f6-7e4eebb74e94)

### Non Epithelial 

Figure 5: UMAP of non-epithelials (immune) cells and bar plot of the cell-type distribution.

![image](https://github.com/user-attachments/assets/1fb39bbf-3861-4c25-87ca-cc1a57dfe710)


Figure 6: Feature plot of T cell markers

![image](https://github.com/user-attachments/assets/0eb89ae2-0914-4229-b5a9-38a1d2620a19)


### Limitations
Diversity and representativeness in the datasets is limited. 


### Future Directions:
We propose to evaluate invasive colorectal cancer scRNAseq data by race and tumor type and compare to the current findings.  We will complete analyses in European individuals. We propose to add covariates analyses including sex and age group for African American and European individuals for normal and premalignant tissue. 


### Acknowledgements
Thank you to the 2024 HTAN Data Jamboree Team and Technical Support Team.

