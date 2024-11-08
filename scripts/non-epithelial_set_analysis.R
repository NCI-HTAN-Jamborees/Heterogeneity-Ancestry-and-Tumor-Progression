library(tidyverse)
library(ggplot2)
library(Seurat)
library(cowplot)
library(patchwork)
library(gridExtra)

## setup working dir and results dir
working_dir = getwd()
setwd(working_dir)
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

############################## Part 1: filtering################################

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

# add biospecimen information
m_addbiospecimens = m %>%
  left_join(biospecimens_df_need, by= c("HTAN Specimen ID"  = "HTAN_Biospecimen_ID"))

summary_tumor_race = as.data.frame(table(m_addbiospecimens$Tumor_Tissue_Type, m_addbiospecimens$self_reported_ethnicity)) %>%
  pivot_wider(names_from ="Var1", values_from = "Freq") %>%
  rename(
    Race = Var2
  )
write.table(summary_tumor_race,"results/non-epithelial/summary_tumor_race_before_filtering.tab",sep="\t",row.names = F, quote=F)

# filter cells
m_needed = m_addbiospecimens %>%
  filter(self_reported_ethnicity == "European" | self_reported_ethnicity == "African American") %>%
  filter(Tumor_Tissue_Type=="Normal" | Tumor_Tissue_Type=="Premalignant") %>%
  filter(Tumor_Type=="NL" | Tumor_Type=="SSL" | Tumor_Type=="TA")
m_needed$self_reported_ethnicity = as.character(m_needed$self_reported_ethnicity)
m_needed$self_reported_ethnicity = factor(m_needed$self_reported_ethnicity)

# add age group
m_needed = m_needed %>%
  mutate(age = gsub("-year-old stage","",development_stage)) %>%
  mutate(age_group = ifelse(age >60,">60","<=60"))
table(m_needed$age)

table(m_needed$Tumor_Type, m_needed$self_reported_ethnicity)

# obtain cell_barcodes after filtering
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

# save fitlered seurat object
saveRDS(so,"CellxGene/Filtered_no_Epithelial.rds")

############################## Part 2: analysis ################################
so <- readRDS("CellxGene/Filtered_no_Epithelial.rds")

cell.size = 0.1
resolution = 0.2
reduction_name = "umap"
p1 = DimPlot(so,reduction=reduction_name, pt.size =0.1, group.by = c("self_reported_ethnicity"))
p2 = DimPlot(so, reduction=reduction_name,pt.size =0.1, group.by = "Tumor_Tissue_Type")
p3 = DimPlot(so, reduction=reduction_name,pt.size =0.1, split.by="self_reported_ethnicity", group.by = "Cell_Type")
p4 = DimPlot(so, reduction=reduction_name,pt.size =0.1, group.by = "sex")
p5 = DimPlot(so, reduction=reduction_name,pt.size =0.1, group.by = "Tumor_Type")
p6 = DimPlot(so, reduction=reduction_name,pt.size =0.1, group.by = "age_group")
p7 = DimPlot(so, reduction=reduction_name,pt.size =0.1, split.by = "Tumor_Type", group.by = "self_reported_ethnicity")
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
ggsave(plot = bar_cluster,height=8,width=6, filename = paste0(umap_path, "/barplot_cluster.pdf"))

p10 = DimPlot(so, reduction=reduction_name,pt.size =0.1, split.by = "Tumor_Type", group.by = "Cell_Type")

so$Tumor_Type_Race <- paste(so$Tumor_Type, so$self_reported_ethnicity, sep = "_")

p11 = DimPlot(so, reduction=reduction_name,pt.size =0.1, split.by = "Tumor_Type_Race", ncol = 2, group.by = "Cell_Type")

ggsave(paste0 (umap_path,"/UMAP_splitRaceTumorType_groupCellType_so.pdf"), p11, width = 12, height = 10)

long_data_cluster <- m_so[,c("self_reported_ethnicity","Tumor_Type","Cell_Type")] %>%
  group_by(self_reported_ethnicity,Tumor_Type,Cell_Type) %>%
  summarise(count=n()) %>%
  mutate(percentage = count / sum(count) * 100)
long_data_cluster
bar_cluster =ggplot(long_data_cluster, aes(x = self_reported_ethnicity, y = percentage, fill =Cell_Type)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(
    data=filter(long_data_cluster,percentage>1),
    aes(label = count), position = position_fill(vjust = 0.5),
    color = "white") +
  labs(title = "Percentage Bar Plot", x = "self_reported_ethnicity", y = "Percentage") +
  facet_wrap(~ Tumor_Type, ncol=3, scales = "free_x")  +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 20),  # Increase legend title size
        legend.text = element_text(size = 15)) +
  guides(fill = guide_legend(override.aes = list(size = 8)))
bar_cluster
ggsave(plot = bar_cluster,height =8,width =6, filename = paste0(umap_path, "/barplot_cluster_race_tumor_celltype.pdf"))

# Plot gene expression using FeaturePlot

# B-cell: CD20 (ENSG00000156738), CD21 (ENSG00000117322), CD19 (ENSG00000177455)
# Fibroblasts: VIM (ENSG00000026025), FAP (ENSG00000078098), FSP1 (ENSG00000042286), SMA (ENSG00000172062)
# T-cell: CD8 (ENSG00000153563), CD4 (ENSG00000010610), CD3 (ENSG00000198851), CD28 (ENSG00000178562)
# MYE: CD68 (ENSG00000129226), CD80 (ENSG00000121594), CD86 (ENSG00000114013), CD206 (ENSG00000260314), CD163 (ENSG00000177575)
# PLA: CD38 (ENSG00000004468), MUM1 (ENSG00000160953)
# MAS: CD117 (ENSG00000157404)
# END: CD31 (ENSG00000261371), CD34 (ENSG00000174059), VEGF (ENSG00000112715)

genes <- c("ENSG00000156738", "ENSG00000117322", "ENSG00000177455")

p12 = FeaturePlot(
  object = so,
  features = genes,
  split.by = "self_reported_ethnicity",
  reduction = "umap" # or "tsne", depending on your reduction
)
ggsave(paste0 (umap_path,"/B_cell_marker.pdf"), p12, width = 8, height = 8)

# T-cell: CD8 (ENSG00000153563), CD4 (ENSG00000010610), CD3 (ENSG00000198851), CD28 (ENSG00000178562)
genes <- c("ENSG00000153563", "ENSG00000010610", "ENSG00000198851", "ENSG00000178562")
p13 = FeaturePlot(
  object = so,
  features = genes,
  split.by = "self_reported_ethnicity",
  reduction = "umap" # or "tsne", depending on your reduction
)
ggsave(paste0 (umap_path,"/T_cell_marker.pdf"), p13, width = 8, height = 10)


table(so$Tumor_Tissue_Type)
so$Tumor_Tissue_Race <- paste(so$Tumor_Tissue_Type, so$self_reported_ethnicity, sep = "_")

table(so$Tumor_Tissue_Race)

Idents(so) <- so$Tumor_Tissue_Race
AA_markers <- FindMarkers(so, ident.2 = "Normal_African American", ident.1 = "Premalignant_African American")
# EU_markers <- FindMarkers(so, ident.2 = "Normal_European", ident.1 = "Premalignant_European")

AA_markers <- AA_markers %>%
  rownames_to_column(var = "gene")

# Step 1: Convert Ensembl IDs to gene symbols
AA_gene_conversion <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = AA_markers$gene,
  keytype = "ENSEMBL",
  columns = "SYMBOL"
)

# Rename columns for clarity
colnames(AA_gene_conversion) <- c("ensembl_gene_id", "gene_name")

# Merge the converted gene names back with the markers data
AA_markers <- AA_markers %>%
  left_join(AA_gene_conversion, by = c("gene" = "ensembl_gene_id"))

# Step 2: Filter based on |avg_log2FC| > 2 and p_val_adj < 0.01
AA_filtered_markers <- AA_markers %>%
  filter(!is.na(gene_name), abs(avg_log2FC) > 2, p_val_adj < 0.01) %>%
  arrange(desc(avg_log2FC)) %>%
  distinct(gene_name, .keep_all = TRUE)

write.csv(AA_filtered_markers, "results/non-epithelial/AA_filtered_markers.csv")

###############################################################

sessionInfo()

