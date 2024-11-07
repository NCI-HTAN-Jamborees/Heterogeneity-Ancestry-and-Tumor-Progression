# Heterogeneity-Ancestry-and-Tumor-Progression
### Team Lead: Isra Elhussin

### Team Members: Jiaying Lai, Ying Wu, Mohamed Abdalla, Aditi Hazra

### Background and Rationale:
Black patients have an almost 20% higher incidence of colorectal cancer compared to Non-Hispanic White (NHW) patients [ACS Cancer Facts and Figures 2024]. Compelling evidence suggests there are differences in tumor microenvironment between Black or African American and European American cancer patients. However, population descriptors, including the social construct of race, may not adequately capture the complex patterns of continuous human genetic variation. Why immuno-inflammation pathways and immune cells are upregulated in Black patients remains unknown. Furthermore, the progression from colorectal adenoma (polyp) to invasive colorectal cancer is understudied in Black patients (Image: SciPro/Getty Images; Reference: They found colon polyps: Now what? Harvard Health. July 20, 2023: https://www.health.harvard.edu/diseases-and-conditions/they-found-colon-polyps-now-what). We propose to evaluate the spatial transcriptome by genetic ancestry, an individual’s family tree by which they inherit DNA from specific ancestors [Reference National Academies of Science, Engineering, and Medicine (NASEM)]. 
![image](https://github.com/user-attachments/assets/d7e7f0b4-d6bc-4a73-bf41-dfa5d9052b9d)

### Objective: 
We propose to conduct a comparative analysis of single-cell RNA sequencing (scRNAseq) data from normal tissue samples from Black and White patients. We will explore the continuum of hyperplasia and premalignant adenomas (colorectal polyps). The primary objective is to elucidate distinct immune signature disparities between normal and premalignant tissues within each racial group (Black versus White) and to ascertain their potential role in tumorigenesis.

### Project Workflow:

Figure 1.  Workflow
![image](https://github.com/user-attachments/assets/23678e30-6083-4564-856e-37fc3daca475)


### Data:
The Vanderbilt University Medical Center (VUMC) scRNA-seq data and spatial transcriptomics data used in this project were accessed using the HTAN Data Portal and Google BigQuery v6.

Ancestry Data (TBA):
We define high ancestry as 70%.

### Data Collection and Processing:

We accessed and extracted data from the HTAN portal and subsequently matched sample IDs with ancestry and clinical annotations. The dataset was then filtered based on self-reported race/ancestry, focusing on individuals identified as Black, who were represented solely by normal and premalignant sample types. Additionally, we incorporated demographic variables, including sex (female, male) and age-stratified by the median (≤60 years, >60 years), as well as polyp classification (NL, SSL, TL).
