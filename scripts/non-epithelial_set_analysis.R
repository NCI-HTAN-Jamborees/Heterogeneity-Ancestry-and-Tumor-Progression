library(tidyverse)
library(ggplot2)
library(Seurat)
library(cowplot)
library(patchwork)
library(gridExtra)

## setup working dir and results dir
working_dir = getwd()
setwd(working_dir)
# getwd()
output_path = paste0(working_dir,"/results/non-epithelial")
dir.create(output_path, showWarnings = F)
input_path = paste0(working_dir,"/CellxGene")
cat("input_path: ", input_path, "\n" )
cat("out_path: ", output_path,"\n" )

integration_plot = output_path
umap_path = paste0(output_path,"/umap")
feature_path = paste0(output_path,"/feature_plot")
barplot_path = paste0(output_path,"/celltype_plot")
marker_path = paste0(output_path,"/marker_plot")
dir.create(umap_path,showWarnings=F)
dir.create(barplot_path,showWarnings=F)
dir.create(feature_path,showWarnings=F)
dir.create(marker_path,showWarnings=F)
plot_path = output_path

getwd()
#This file is downloaded from https://cellxgene.cziscience.com/collections/a48f5033-3438-4550-8574-cdff3263fdfd
object <- readRDS(paste0(input_path,"/non-epithelial_data.rds"))
colors <- c('lightgrey', 'navy')

cell.size = 0.1
resolution = 0.2
reduction_name = "umap"

# Add cell carcode to meta data
m = object@meta.data %>%
  rownames_to_column(var ="cell_barcode")
# rownames(m) = rownames(object@meta.data)
# colnames(m)
# rownames(m)[1000:1005]
# m$cell_barcode[1000:1005]

# table(m$self_reported_ethnicity)
# table(m$disease)
# m$`HTAN Specimen ID`
# getwd()

# Load metadata file with biospecimens info
biospecimens_df= read.table("metadata/biospecimens.tsv",header=T,sep="\t")
#colnames(biospacemens_df)
dim(biospecimens_df)
#biospecimens_df_need = biospecimens_df[,c("HTAN.Biospecimen.ID","HTAN.Parent.ID", "Collection.Days.from.Index", "Storage.Method", "Processing.Days.from.Index", "Preservation.Method", "Tumor.Tissue.Type")]
#colnames(biospecimens_df_need) = c("HTAN_Biospecimen_ID","HTAN_Participant_ID", "Collection_Days_from_Index", "Storage_Method", "Processing_Days_from_Index", "Preservation_Method", "Tumor_Tissue_Type")
biospecimens_df_need = biospecimens_df[,c("HTAN.Biospecimen.ID","Tumor.Tissue.Type")]
colnames(biospecimens_df_need) = c("HTAN_Biospecimen_ID","Tumor_Tissue_Type")
biospecimens_df_need
cat("There are ",length(unique(biospecimens_df_need$HTAN_Biospecimen_ID)), " biospecimens from  Vanderbilt Colon Molecular Atlas Project")
#head(biospecimens_df_need)
as.data.frame(table(biospecimens_df_need$Tumor_Tissue_Type))
#case_df_need = case_df[,c("HTAN.Participant.ID","Ethnicity","Gender","Race")] #,"Age.at.Diagnosis..years.","Primary.Diagnosis","Year.of.Diagnosis")]
#colnames(case_df_need) =c("HTAN_Participant_ID","Ethnicity","Gender","Race")

m_addbiospecimens = m %>%
  left_join(biospecimens_df_need, by= c("HTAN Specimen ID"  = "HTAN_Biospecimen_ID"))
# rownames(m_addbiospecimens) = rownames(m)

summary_tumor_race = as.data.frame(table(m_addbiospecimens$Tumor_Tissue_Type, m_addbiospecimens$self_reported_ethnicity)) %>%
  pivot_wider(names_from ="Var1", values_from = "Freq") %>%
  rename(
    Race = Var2
  )
write.table(summary_tumor_race,"results/non-epithelial/summary_tumor_race_before_filtering.tab",sep="\t",row.names = F, quote=F)
m_needed = m_addbiospecimens %>%
  filter(self_reported_ethnicity == "European" | self_reported_ethnicity == "African American") %>%
  filter(Tumor_Tissue_Type=="Normal" | Tumor_Tissue_Type=="Premalignant")
m_needed$self_reported_ethnicity = as.character(m_needed$self_reported_ethnicity)
m_needed$self_reported_ethnicity = factor(m_needed$self_reported_ethnicity)
colnames(m_addbiospecimens)

cells_needed = m_needed$cell_barcode
m_need_for_metadata = m_needed[,c("cell_barcode","Tumor_Tissue_Type")]
rownames(m_need_for_metadata) = m_need_for_metadata$cell_barcode
so=object
so = subset(so,cells=cells_needed)
so
so = AddMetaData(so,m_need_for_metadata)
summary_tumor_race_after = as.data.frame(table(m_needed$Tumor_Tissue_Type, m_needed$self_reported_ethnicity)) %>%
  pivot_wider(names_from ="Var1", values_from = "Freq") %>%
  rename(
    Race = Var2
  )
summary_tumor_race_after
write.table(summary_tumor_race_after,"summary_tumor_race_before_filtering.tab",sep="\t",row.names = F, quote=F)




