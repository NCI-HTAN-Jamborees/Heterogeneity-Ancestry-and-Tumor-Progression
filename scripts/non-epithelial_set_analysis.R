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

# Load metadata file with biospecimens info
biospecimens_df= read.table("metadata/biospecimens.tsv",header=T,sep="\t")
dim(biospecimens_df)
biospecimens_df_need = biospecimens_df[,c("HTAN.Biospecimen.ID","Tumor.Tissue.Type")]
colnames(biospecimens_df_need) = c("HTAN_Biospecimen_ID","Tumor_Tissue_Type")
biospecimens_df_need
cat("There are ",length(unique(biospecimens_df_need$HTAN_Biospecimen_ID)), " biospecimens from  Vanderbilt Colon Molecular Atlas Project")
as.data.frame(table(biospecimens_df_need$Tumor_Tissue_Type))

m_addbiospecimens = m %>%
  left_join(biospecimens_df_need, by= c("HTAN Specimen ID"  = "HTAN_Biospecimen_ID"))

summary_tumor_race = as.data.frame(table(m_addbiospecimens$Tumor_Tissue_Type, m_addbiospecimens$self_reported_ethnicity)) %>%
  pivot_wider(names_from ="Var1", values_from = "Freq") %>%
  rename(
    Race = Var2
  )
write.table(summary_tumor_race,"results/non-epithelial/summary_tumor_race_before_filtering.tab",sep="\t",row.names = F, quote=F)

m_needed = m_addbiospecimens %>%
  filter(self_reported_ethnicity == "European" | self_reported_ethnicity == "African American") %>%
  filter(Tumor_Tissue_Type=="Normal" | Tumor_Tissue_Type=="Premalignant") %>%
  filter(Tumor_Type=="NL" | Tumor_Type=="SSL" | Tumor_Type=="TA")
m_needed$self_reported_ethnicity = as.character(m_needed$self_reported_ethnicity)
m_needed$self_reported_ethnicity = factor(m_needed$self_reported_ethnicity)

# colnames(m_addbiospecimens)
m_needed = m_needed %>%
  mutate(age = gsub("-year-old stage","",development_stage)) %>%
  mutate(age_group = ifelse(age >60,">60","<=60"))
table(m_needed$age)

table(m_needed$Tumor_Type, m_needed$self_reported_ethnicity)

cells_needed = m_needed$cell_barcode

m_need_for_metadata = m_needed[,c("cell_barcode","Tumor_Tissue_Type","age_group")]
rownames(m_need_for_metadata) = m_need_for_metadata$cell_barcode

so=object
so = subset(so,cells=cells_needed)
so = AddMetaData(so,m_need_for_metadata)
summary_tumor_race_after = as.data.frame(table(m_needed$Tumor_Tissue_Type, m_needed$self_reported_ethnicity)) %>%
  pivot_wider(names_from ="Var1", values_from = "Freq") %>%
  rename(
    Race = Var2
  )
summary_tumor_race_after
write.table(summary_tumor_race_after,"results/non-epithelial/summary_tumor_race_before_filtering.tab",sep="\t",row.names = F, quote=F)

saveRDS(so,"CellxGene/Filtered_no_Epithelial.rds")

cell.size = 0.1
resolution = 0.2
reduction_name = "umap"
#p1 = DimPlot(so,reduction="umap.integrated.rpca", pt.size =0.1, group.by = c("condition"))
p1 = DimPlot(so,reduction=reduction_name, pt.size =0.1, group.by = c("self_reported_ethnicity"))
p2 = DimPlot(so, reduction=reduction_name,pt.size =0.1, group.by = "Tumor_Tissue_Type")
#p2 = DimPlot(so, reduction=reduction_name,group.by = "rna_clusters")
p3 = DimPlot(so, reduction=reduction_name,pt.size =0.1, split.by="self_reported_ethnicity", group.by = "Cell_Type")
p4 = DimPlot(so, reduction=reduction_name,pt.size =0.1, group.by = "sex")
p5 = DimPlot(so, reduction=reduction_name,pt.size =0.1, group.by = "Tumor_Type")
p6 = DimPlot(so, reduction=reduction_name,pt.size =0.1, group.by = "age_group")
p7 = DimPlot(so, reduction=reduction_name,pt.size =0.1, split.by = "Tumor_Type", group.by = "self_reported_ethnicity")
#p1+p2+p3
#p4+p5+p6


p9 = DimPlot(so, reduction=reduction_name,pt.size =0.1, split.by = "self_reported_ethnicity", group.by = "Cell_Type")

m_so = so@meta.data

long_data_cluster <- m_so[,c("self_reported_ethnicity","Cell_Type")] %>% 
  group_by(self_reported_ethnicity,Cell_Type) %>%
  summarise(count=n()) %>%
  mutate(percentage = count / sum(count) * 100)

long_data_cluster
bar_cluster =ggplot(long_data_cluster, aes(x = self_reported_ethnicity, y = percentage, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(
    data=filter(long_data_cluster,percentage>1),
    aes(label = count), position = position_fill(vjust = 0.5),
    color = "white") +
  labs(title = "Percentage Bar Plot", x = "SampleID", y = "Percentage") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 20),  # Increase legend title size
        legend.text = element_text(size = 15)) +
  guides(fill = guide_legend(override.aes = list(size = 8)))
bar_cluster
ggsave(plot = bar_cluster,height =4,width =5, filename = paste0(umap_path, "/barplot_cluster.pdf"))


p10 = DimPlot(so, reduction=reduction_name,pt.size =0.1, split.by = "Tumor_Type", group.by = "Cell_Type")

so$Tumor_Type_Race <- paste(so$Tumor_Type, so$self_reported_ethnicity, sep = "_")

p11 = DimPlot(so, reduction=reduction_name,pt.size =0.1, split.by = "Tumor_Type_Race", ncol = 2, group.by = "Cell_Type")

ggsave(paste0 (umap_path,"/UMAP_splitRaceTumorType_groupCellType_so.pdf"), p11, width = 12, height = 5)

##################

p1
ggsave(paste0 (umap_path,"/UMAP_splitPolyb_groupRace_so.pdf"), p7, width = 18, height = 5)

combined_plot <- p1 + p2 + p3 + p4 +p5+p6
# Arrange plots into a 2x3 grid
combined_plot <- combined_plot + plot_layout(ncol = 3)
plot(combined_plot) 
# Save the combined plot as a PDF file
ggsave(paste0 (umap_path,"/1.UMAP_so.pdf"), combined_plot, width = 18, height = 10)


table(m_so$development_stage)
summary_polyb_race_after = as.data.frame(table(m_so$Polyp_Type, m_so$self_reported_ethnicity)) %>%
  pivot_wider(names_from ="Var1", values_from = "Freq") %>%
  rename(
    Race = Var2
  )
summary_polyb_race_after


## subset based on Polyp_Type
m_final = m_so %>%
  filter(Polyp_Type != "TV" & Polyp_Type != "UNC")
cells_needed_final = m_final$cell_barcode
so_final = subset(so,cells = cells_needed_final)



#p1 = DimPlot(so,reduction="umap.integrated.rpca", pt.size =0.1, group.by = c("condition"))
p1 = DimPlot(so_final,reduction=reduction_name, pt.size =0.1, group.by = c("self_reported_ethnicity"))
p2 = DimPlot(so_final, reduction=reduction_name,pt.size =0.1, group.by = "Tumor_Tissue_Type")
#p2 = DimPlot(so, reduction=reduction_name,group.by = "rna_clusters")
p3 = DimPlot(so_final, reduction=reduction_name,pt.size =0.1, group.by = "Cell_Type")
p4 = DimPlot(so_final, reduction=reduction_name,pt.size =0.1, group.by = "sex")
p5 = DimPlot(so_final, reduction=reduction_name,pt.size =0.1, group.by = "Polyp_Type")
p6 = DimPlot(so_final, reduction=reduction_name,pt.size =0.1, group.by = "age_group")
p7 = DimPlot(so, reduction=reduction_name,pt.size =0.1, split.by = "Polyp_Type", group.by = "self_reported_ethnicity")
p8 = DimPlot(so_final, reduction=reduction_name,pt.size =0.1, split.by = "Polyp_Type", group.by = "self_reported_ethnicity")

p1+p2+p3
p4+p5+p6

p8
ggsave(paste0 (umap_path,"/1.UMAP_splitPolyb_groupRace_so_NLSSLTA.pdf"), p8, width = 12, height = 5)

combined_plot <- p1 + p2 + p3 + p4 +p5+p6
# Arrange plots into a 2x3 grid
combined_plot <- combined_plot + plot_layout(ncol = 3)
plot(combined_plot) 
# Save the combined plot as a PDF file
ggsave(paste0 (umap_path,"/1.UMAP_so_NLSSLTA.pdf"), combined_plot, width = 18, height = 10)





sessionInfo()

