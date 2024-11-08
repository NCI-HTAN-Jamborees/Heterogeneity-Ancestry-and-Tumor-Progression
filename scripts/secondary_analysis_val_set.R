library(tidyverse)
library(ggplot2)
library(Seurat)
library(cowplot)
library(patchwork)
library(gridExtra)

## setup working dir and results dir
working_dir = "/data/wuy24/HTAN"
setwd(working_dir)
getwd()
output_path = paste0(working_dir,"/results/Val_set_of_human_colorectal_tumor_Epithelial")
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


## load object
getwd()
object <- readRDS(paste0(input_path,"/Val_set_of_human_colorectal_tumor_Epithelial.rds"))
colors <- c('lightgrey', 'navy')


## add Tumor.Tissue.Type from biospecimens metadata to seurat object metadata
m = object@meta.data %>%
  rownames_to_column(var ="cell_barcode")
rownames(m) = rownames(object@meta.data)
colnames(m)
rownames(m)[1000:1005]
m$cell_barcode[1000:1005]

table(m$self_reported_ethnicity)
table(m$disease)
m$`HTAN Specimen ID`
getwd()
biospecimens_df= read.table("Heterogeneity-Ancestry-and-Tumor-Progression/metadata/biospecimens.tsv",header=T,sep="\t")
#colnames(biospacemens_df)
dim(biospecimens_df)
#biospecimens_df_need = biospecimens_df[,c("HTAN.Biospecimen.ID","HTAN.Parent.ID", "Collection.Days.from.Index", "Storage.Method", "Processing.Days.from.Index", "Preservation.Method", "Tumor.Tissue.Type")]
#colnames(biospecimens_df_need) = c("HTAN_Biospecimen_ID","HTAN_Participant_ID", "Collection_Days_from_Index", "Storage_Method", "Processing_Days_from_Index", "Preservation_Method", "Tumor_Tissue_Type")
biospecimens_df_need = biospecimens_df[,c("HTAN.Biospecimen.ID","Tumor.Tissue.Type")]
colnames(biospecimens_df_need) = c("HTAN_Biospecimen_ID","Tumor_Tissue_Type")
head(biospecimens_df_need)
cat("There are ",length(unique(biospecimens_df_need$HTAN_Biospecimen_ID)), " biospecimens from  Vanderbilt Colon Molecular Atlas Project")
#head(biospecimens_df_need)
as.data.frame(table(biospecimens_df_need$Tumor_Tissue_Type))
#case_df_need = case_df[,c("HTAN.Participant.ID","Ethnicity","Gender","Race")] #,"Age.at.Diagnosis..years.","Primary.Diagnosis","Year.of.Diagnosis")]
#colnames(case_df_need) =c("HTAN_Participant_ID","Ethnicity","Gender","Race")

m_addbiospecimens = m %>%
  left_join(biospecimens_df_need, by= c("HTAN Specimen ID"  = "HTAN_Biospecimen_ID"))
rownames(m_addbiospecimens) = rownames(m)

summary_tumor_race = as.data.frame(table(m_addbiospecimens$Tumor_Tissue_Type, m_addbiospecimens$self_reported_ethnicity)) %>%
  pivot_wider(names_from ="Var1", values_from = "Freq") %>%
  rename(
    Race = Var2
  )
write.table(summary_tumor_race,paste0(output_path ,"/summary_tumor_race_before_filtering.tab"),sep="\t",row.names = F, quote=F)

m_needed = m_addbiospecimens %>%
  filter(self_reported_ethnicity != "unknown") %>%
  filter(Tumor_Tissue_Type=="Normal" | Tumor_Tissue_Type=="Premalignant")
m_needed$self_reported_ethnicity = as.character(m_needed$self_reported_ethnicity)
m_needed$self_reported_ethnicity = factor(m_needed$self_reported_ethnicity)
colnames(m_addbiospecimens)

m_needed = m_needed %>% 
  mutate(age = gsub("-year-old stage","",development_stage)) %>%
  mutate(age_group = ifelse(age >60,">60","<=60"))
m_needed$age = as.numeric(m_needed$age )
table(m_needed$age)
summary(m_needed$age)

cells_needed = m_needed$cell_barcode
m_need_for_metadata = m_needed[,c("cell_barcode","Tumor_Tissue_Type","age_group")]
rownames(m_need_for_metadata) = m_need_for_metadata$cell_barcode
so=object
so = subset(so,cells=cells_needed)
dim(so@meta.data)
so = AddMetaData(so,m_need_for_metadata)
summary_tumor_race_after = as.data.frame(table(m_needed$Tumor_Tissue_Type, m_needed$self_reported_ethnicity)) %>%
  pivot_wider(names_from ="Var1", values_from = "Freq") %>%
  rename(
    Race = Var2
  )
summary_tumor_race_after
write.table(summary_tumor_race_after,paste0(output_path,"/summary_tumor_race_after_filtering.tab"),sep="\t",row.names = F, quote=F)

#saveRDS(so,"CellxGene/Filtered_Val_set_of_human_colorectal_tumor_Epithelial.rds")

## now generate UMAP for filtered object so in which all the cells from unknown self_reported_ethnicity was removed
cell.size = 0.1
resolution = 0.2
reduction_name = "umap"
#p1 = DimPlot(so,reduction="umap.integrated.rpca", pt.size =0.1, group.by = c("condition"))
p1 = DimPlot(so,reduction=reduction_name, pt.size =0.1, group.by = c("self_reported_ethnicity"))
p2 = DimPlot(so, reduction=reduction_name,pt.size =0.1, group.by = "Tumor_Tissue_Type")
#p2 = DimPlot(so, reduction=reduction_name,group.by = "rna_clusters")
p3 = DimPlot(so, reduction=reduction_name,pt.size =0.1, group.by = "Cell_Type")
p4 = DimPlot(so, reduction=reduction_name,pt.size =0.1, group.by = "sex")
p5 = DimPlot(so, reduction=reduction_name,pt.size =0.1, group.by = "Polyp_Type")
p6 = DimPlot(so, reduction=reduction_name,pt.size =0.1, group.by = "age_group")
p7 = DimPlot(so, reduction=reduction_name,pt.size =0.1, split.by = "Polyp_Type", group.by = "self_reported_ethnicity")
#p1+p2+p3
#p4+p5+p6

p7
ggsave(paste0 (umap_path,"/1.UMAP_splitPolyb_groupRace_so.pdf"), p7, width = 18, height = 5)

combined_plot <- p1 + p2 + p3 + p4 +p5+p6
# Arrange plots into a 2x3 grid
combined_plot <- combined_plot + plot_layout(ncol = 3)
plot(combined_plot) 
# Save the combined plot as a PDF file
ggsave(paste0 (umap_path,"/1.UMAP_so.pdf"), combined_plot, width = 18, height = 10)


m_so = so@meta.data
table(m_so$development_stage)
summary_polyb_race_after = as.data.frame(table(m_so$Polyp_Type, m_so$self_reported_ethnicity)) %>%
  pivot_wider(names_from ="Var1", values_from = "Freq") %>%
  rename(
    Race = Var2
  )
summary_polyb_race_after
write.table(summary_polyb_race_after, paste0(output_path,"/summary_polyb_race_after_filterUnknownCell.tab"),
            sep="\t",row.names=F,quote=F)

## subset based on Polyp_Type and remove TV and UNC group which do not have cells in African American
m_final = m_so %>%
  filter(Polyp_Type != "TV" & Polyp_Type != "UNC") %>%
  mutate(Tumor_PolypType = paste0(Tumor_Tissue_Type, "_",Polyp_Type)) %>%
  mutate(Tumor_CellType= paste0(Tumor_Tissue_Type, "_",Cell_Type)) %>%
  mutate(Race_TumorType = paste0(self_reported_ethnicity, "_", Tumor_Tissue_Type)) %>%
  mutate(Race_PolypType = paste0(self_reported_ethnicity, "_", Polyp_Type)) %>%
  mutate(Race_CellType = paste0(self_reported_ethnicity, "_", Cell_Type))
cells_needed_final = m_final$cell_barcode
so_final = subset(so,cells = cells_needed_final)

m_final_add = m_final[,c("cell_barcode","Tumor_PolypType","Tumor_CellType","Race_TumorType","Race_PolypType","Race_CellType")]
dim(m_final_add)
rownames(m_final_add) = m_final_add$cell_barcode
m_final_add$Race_TumorType = factor(m_final_add$Race_TumorType, levels = c("European_Normal" ,"African American_Normal","European_Premalignant", "African American_Premalignant"))
m_final_add$Race_PolypType = factor(m_final_add$Race_TumorType, levels = c("European_NL" ,"African American_NL","European_SSL", "African American_SSL","European_TA","African American_TA"))

so_final = AddMetaData(so_final,m_final_add)


p1 = DimPlot(so_final,reduction=reduction_name, pt.size =0.1, group.by = c("self_reported_ethnicity"))
p2 = DimPlot(so_final, reduction=reduction_name,pt.size =0.1, group.by = "Tumor_Tissue_Type")
p3 = DimPlot(so_final, reduction=reduction_name,pt.size =0.1, group.by = "Cell_Type")
p4 = DimPlot(so_final, reduction=reduction_name,pt.size =0.1, group.by = "sex")
p5 = DimPlot(so_final, reduction=reduction_name,pt.size =0.1, group.by = "Polyp_Type")
p6 = DimPlot(so_final, reduction=reduction_name,pt.size =0.1, group.by = "age_group")
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


## Generate bar plot and Dimplot split by two factor
p9 = DimPlot(so_final, reduction=reduction_name,pt.size =0.1, split.by = "self_reported_ethnicity", group.by = "Cell_Type")
p9 
ggsave(paste0 (umap_path,"/1.UMAP_so_NLSSLTA_splitRace_groupCelltype.pdf"), p9, width = 10, height = 5)

## barplot for cell_type by Race
long_data_cluster <- m_final[,c("self_reported_ethnicity","Cell_Type")] %>% 
  group_by(self_reported_ethnicity,Cell_Type) %>%
  summarise(count=n()) %>%
  mutate(percentage = count / sum(count) * 100)

long_data_cluster
bar_cluster_celltype =ggplot(long_data_cluster, aes(x = self_reported_ethnicity, y = percentage, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(
    data=filter(long_data_cluster,percentage>1),
    aes(label = count), position = position_fill(vjust = 0.5),
    color = "white") +
  labs(title = "Percentage Bar Plot", x = "self_reported_ethnicity", y = "Percentage") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 20),  # Increase legend title size
        legend.text = element_text(size = 15)) +
  guides(fill = guide_legend(override.aes = list(size = 8)))
bar_cluster_celltype
ggsave(plot = bar_cluster_celltype,height =4,width =5, filename = paste0(umap_path, "/4.barplot_cell_type.pdf"))

## barplot for cell_type/polyp_type by RaceTumor
p10 = DimPlot(so_final, reduction=reduction_name,pt.size =0.1, split.by = "self_reported_ethnicity", group.by = "Cell_Type")
p10
ggsave(paste0 (umap_path,"/1.UMAP_so_NLSSLTA_splitRace_groupCelltype.pdf"), p10, width = 6, height = 5)
p11 = DimPlot(so_final, reduction=reduction_name,pt.size =0.1, split.by = "Race_TumorType"  , group.by = "Cell_Type")
p11
ggsave(paste0 (umap_path,"/1.UMAP_so_NLSSLTA_splitRaceTumor_groupCelltype.pdf"), p11, width = 13, height = 5)
p12 = DimPlot(so_final, reduction=reduction_name,pt.size =0.1, split.by = "Race_TumorType"  , group.by = "Polyp_Type")
p12
ggsave(paste0 (umap_path,"/1.UMAP_so_NLSSLTA_splitRaceTumor_groupPolyptype.pdf"), p11, width = 13, height = 5)
long_data_cluster <- m_final[,c("self_reported_ethnicity","Tumor_Tissue_Type","Polyp_Type")] %>% 
  group_by(self_reported_ethnicity,Tumor_Tissue_Type,Polyp_Type) %>%
  summarise(count=n()) %>%
  mutate(percentage = count / sum(count) * 100)

long_data_cluster

bar_cluster_RaceTumor_Polyptype =ggplot(long_data_cluster, aes(x = self_reported_ethnicity, y = percentage, fill =Polyp_Type)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(
    data=filter(long_data_cluster,percentage>1),
    aes(label = count), position = position_fill(vjust = 0.5),
    color = "white") +
  labs(title = "Percentage Bar Plot", x = "self_reported_ethnicity", y = "Percentage") +
  facet_wrap(~ Tumor_Tissue_Type, ncol=2, scales = "free_x")  +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 20),  # Increase legend title size
        legend.text = element_text(size = 15)) +
  guides(fill = guide_legend(override.aes = list(size = 8)))
bar_cluster_RaceTumor_Polyptype
ggsave(filename = paste0(umap_path, "/4.barplot_cluster_race_tumor_polyp.pdf"),plot =bar_cluster_RaceTumor_Polyptype,height =4,width =5, )


long_data_cluster <- m_final[,c("self_reported_ethnicity","Tumor_Tissue_Type","Cell_Type")] %>% 
  group_by(self_reported_ethnicity,Tumor_Tissue_Type,Cell_Type) %>%
  summarise(count=n()) %>%
  mutate(percentage = count / sum(count) * 100)

long_data_cluster

bar_cluster_raceTumor_celltype =ggplot(long_data_cluster, aes(x = self_reported_ethnicity, y = percentage, fill =Cell_Type)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(
    data=filter(long_data_cluster,percentage>1),
    aes(label = count), position = position_fill(vjust = 0.5),
    color = "white") +
  labs(title = "Percentage Bar Plot", x = "self_reported_ethnicity", y = "Percentage") +
  facet_wrap(~ Tumor_Tissue_Type, ncol=2, scales = "free_x")  +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 20),  # Increase legend title size
        legend.text = element_text(size = 15)) +
  guides(fill = guide_legend(override.aes = list(size = 8)))
bar_cluster_raceTumor_celltype
ggsave(filename = paste0(umap_path, "/4.barplot_cluster_race_tumor_celltype.pdf"),plot = bar_cluster_raceTumor_celltype,height =4,width =5, )


## example for feature plot for umap split by race_normal
features_EP = c("KI67","CK8","CK14","CK18")
features_EP_ensembl = c("ENSG00000148773", "ENSG00000170421","ENSG00000186847","ENSG00000111057")
names(features_EP_ensembl) = features_EP 
features_EP_ensembl
#features_ABS = c()
#features_ASC = c()
#features_TAC = c()
p21 = FeaturePlot(so_final,features=features_EP_ensembl, pt.size=0.1, ncol=1, cols=colors, reduction= reduction_name)
p21
so_final
library(patchwork)
plot_list = list()
for (i in c(1:length(features_EP_ensembl))) {
  gene_name <- names(features_EP_ensembl)[i]
  gene_id = features_EP_ensembl[i]
  plot <- FeaturePlot(so_final, features = gene_id, pt.size = 0.2, cols = colors, reduction = reduction_name) + 
    ggtitle(paste0(gene_name,"\n",gene_id)) + 
    theme(plot.title = element_text(size = 10)) 
  plot_list[[i]] <- plot
}
combined_plot2 <- wrap_plots(plot_list, ncol = 1)
combined_plot2

plot_list_racetumor = list()
for (i in c(1:length(features_EP_ensembl))) {
  gene_name <- names(features_EP_ensembl)[i]
  gene_id = features_EP_ensembl[i]
  plot <- FeaturePlot(so_final, features = gene_id, pt.size = 0.2, cols = colors,
                      split.by = "Race_TumorType",reduction = reduction_name) 
  text_plot <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = gene_name,
             size =5, angle = 270, hjust = 0.5, vjust =1.2) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 16))
  
  # Combine the feature plot and text plot
  combined_plot <-  plot / text_plot  + plot_layout(ncol = 2, widths = c(20, 2))  # Adjust layout as needed
  combined_plot
  plot_list_racetumor[[i]] <- combined_plot
}
featureplot_splitRaceTumor <- wrap_plots(plot_list_racetumor, ncol = 1)
featureplot_splitRaceTumor
ggsave(paste0(feature_path,"/feature_EP_featureplot_splitRaceTumor.pdf"),featureplot_splitRaceTumor,height=14,width =14)

## example for vlnplot for markers

vlnplot = VlnPlot(so_final, features = features_EP_ensembl, pt.size = 0, ncol = 2, y.max=5,group.by = "Race_TumorType") 
ggsave(paste0(feature_path,"/feature_EP_vlnplot_RaceTumor.pdf"),vlnplot,height=10,width=8)

so_final

table(m$donor_id,m$self_reported_ethnicity)
m_sub = m[,c("donor_id","self_reported_ethnicity","sex")] %>%
  distinct()
table(m_sub$self_reported_ethnicity)
table(m_sub$sex,m_sub$self_reported_ethnicity)

# DE analysis
# 1.compare American_african group Premalignant vs Normal Tumor_tissue_type
# 2.compare European group Premalignant vs Normal Tumor_tissue_type
findmarker_dir = paste0(output_path,"/markers")
dir.create(marker_dir,showWarnings = F)
colnames(so_final@meta.data)
table(m_final$Race_TumorType)
so_working = so_final
Idents(so_final) = 'Race_TumorType'

markers_AA_prem_norm <- so_final %>% FindMarkers(assay='RNA', ident.1 ="African American_Premalignant",
                                       ident.2 = "African American_Normal",
                                       logfc.threshold = -Inf, , min.pct=0.25, min.diff.pct = -Inf)
dim(markers_AA_prem_norm)
markers_AA_prem_norm$gene = rownames(markers_AA_prem_norm)
head(markers_AA_prem_norm)
columns(org.Hs.eg.db)
df_entrezid =  data.frame(
  gene_entrezid = mapIds(
    org.Hs.eg.db,
    keys = unique(markers_AA_prem_norm$gene),
    keytype = "ENSEMBL",
    column = "ENTREZID",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  tibble::rownames_to_column("Ensembl")

head(df_entrezid)

AA_marker_mapped_df <- data.frame(
  gene_symbol = mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Hs.eg.db,
    keys = unique(markers_AA_prem_norm$gene),
    # Replace with the type of gene identifiers in your data
    keytype = "ENSEMBL",
    # Replace with the type of gene identifiers you would like to map to
    column = "SYMBOL",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  # If an Ensembl gene identifier doesn't map to a gene symbol, drop that
  # from the data frame
  dplyr::filter(!is.na(gene_symbol)) %>%
  # Make an `Ensembl` column to store the rownames
  tibble::rownames_to_column("Ensembl") %>%
  dplyr::full_join(df_entrezid, by = c("Ensembl" = "Ensembl")) %>%
  # Now let's join the rest of the expression data
  dplyr::full_join(markers_AA_prem_norm, by = c("Ensembl" = "gene")) %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol), Ensembl, gene_symbol))
ensembl2symbol = unique(AA_marker_mapped_df[,c(1,2)])
ensembl_to_symbol <- setNames(ensembl2symbol$gene_symbol, ensembl2symbol$Ensembl)
AA_marker_mapped_df = AA_marker_mapped_df %>% 
  arrange (p_val_adj) 

AA_de = AA_marker_mapped_df %>%
  filter(p_val <0.01) %>%
  filter(avg_log2FC>2 | avg_log2FC<(-2))
write.csv(AA_marker_mapped_df,paste0(findmarker_dir,"/AA_marker_premalignant_vs_normal_all.csv"),row.names=F,quote=F)
write.csv(AA_de,paste0(findmarker_dir,"/AA_marker_premalignant_vs_normal_p0.01_log2fc2.csv"),row.names=F,quote=F)
dim(AA_de)
AA_de
dim(EU_de)
## 
markers_EU_prem_norm <- so_final %>% FindMarkers(assay='RNA', ident.1 ="European_Premalignant",
                                                ident.2 = "European_Normal",
                                                logfc.threshold = -Inf, , min.pct=0.25, min.diff.pct = -Inf)
dim(markers_EU_prem_norm)
markers_EU_prem_norm$gene = rownames(markers_EU_prem_norm)
head(markers_EU_prem_norm)
columns(org.Hs.eg.db)
df_entrezid =  data.frame(
  gene_entrezid = mapIds(
    org.Hs.eg.db,
    keys = unique(markers_EU_prem_norm$gene),
    keytype = "ENSEMBL",
    column = "ENTREZID",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  tibble::rownames_to_column("Ensembl")

head(df_entrezid)

EU_marker_mapped_df <- data.frame(
  gene_symbol = mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Hs.eg.db,
    keys = unique(markers_EU_prem_norm$gene),
    # Replace with the type of gene identifiers in your data
    keytype = "ENSEMBL",
    # Replace with the type of gene identifiers you would like to map to
    column = "SYMBOL",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  # If an Ensembl gene identifier doesn't map to a gene symbol, drop that
  # from the data frame
  dplyr::filter(!is.na(gene_symbol)) %>%
  # Make an `Ensembl` column to store the rownames
  tibble::rownames_to_column("Ensembl") %>%
  dplyr::full_join(df_entrezid, by = c("Ensembl" = "Ensembl")) %>%
  # Now let's join the rest of the expression data
  dplyr::full_join(markers_EU_prem_norm, by = c("Ensembl" = "gene")) %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol), Ensembl, gene_symbol))
ensembl2symbol = unique(EU_marker_mapped_df[,c(1,2)])
ensembl_to_symbol <- setNames(ensembl2symbol$gene_symbol, ensembl2symbol$Ensembl)

EU_marker_mapped_df = EU_marker_mapped_df %>% 
  arrange (p_val_adj) 
dim(EU_marker_mapped_df)

EU_de = EU_marker_mapped_df %>%
  filter(p_val <0.01) %>%
  filter(avg_log2FC>2 | avg_log2FC<(-2))
dim(EU_de)
write.csv(EU_marker_mapped_df,paste0(findmarker_dir,"/EU_marker_premalignant_vs_normal_all.csv"),row.names=F,quote=F)
write.csv(EU_de,paste0(findmarker_dir,"/EU_marker_premalignant_vs_normal_p0.01_log2fc2.csv"),row.names=F, quote=F)

top10_AA = AA_de %>% filter(avg_log2FC>0) %>%
   top_n(n =10, wt = avg_log2FC) %>%
  arrange(p_val_adj)
pd10_AA <- DotPlot(
  object = so_final,
  group.by = "Race_TumorType",
  features = unique(top10_AA$Ensembl)
) +
  scale_x_discrete(labels = function(x) ensembl_to_symbol[x]) +  # Map Ensembl to gene_symbol
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
pd10_AA
ggsave(
  plot = pd10_AA,
  height = 4,
  width = 10,
  filename = paste0(marker_path, "/3.AA_top10Marker_dotplot_name.pdf"))

sum(top10_EU %in% top10_AA)
sum( AA_de %in%  EU_de)
top10_EU = EU_de %>% filter(avg_log2FC>0) %>%
  top_n(n =10, wt = avg_log2FC) %>%
  arrange(p_val_adj)
top10_EU
pd10 <- DotPlot(
  object = so_final,
  group.by = "Race_TumorType",
  features = unique(top10_EU$Ensembl)
) +
  scale_x_discrete(labels = function(x) ensembl_to_symbol[x]) +  # Map Ensembl to gene_symbol
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
pd10

ggsave(
  plot = pd10,
  height = 4,
  width = 10,
  filename = paste0(marker_path, "/3.EU_top10Marker_dotplot_name.pdf"))


#saveRDS(so_final,paste0(output_path,"object_final.rds"))

sessionInfo()

