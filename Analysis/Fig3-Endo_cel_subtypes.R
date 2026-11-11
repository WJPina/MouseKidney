library(Seurat)
library(sceasy)
library(reticulate)
library(cowplot)

## ---------- Subset endothelial cells ----------
endo <- subset(sce, subset = celltype == "Endo")

endo <- NormalizeData(endo) %>% FindVariableFeatures() %>% ScaleData() %>% 
  RunPCA() %>% RunUMAP(dim = 1:50)
ElbowPlot(endo,ndims = 50)
endo <- FindNeighbors(endo, dims = 1:20) %>% FindClusters(dims = 1:20) %>% 
  RunUMAP(dims = 1:20)

DimPlot(endo, group.by = "orig.ident", reduction = "umap") + 
  DimPlot(endo, reduction = "umap", label = T, repel = T)

## ---------- Doublet inspection ----------
KidneyMarkers <- list(
  PT = c("Lrp2", "Slc27a2", "Pvalb", "Pgam2", "Uqcrb"),
  DCT = c("Slc12a3", "Wnk1"),
  DLH = c("Bst1", "Epha7", "Tmem207"),
  ALH = c("Umod", "Slc12a1", "Car2", "Slc25a4"),
  `CD-IC` = c("Atp6v0d2", "Atp6v1g3", "Fxyd4", "Apela"),
  `CD-PC` = c("Aqp2", "Hsd11b2", "Magi2", "Cdkn1c"),
  Podo = c("Mafb", "Nphs2", "Adgrl4"),
  Endo = c("Plvap", "Emcn", "Tek", "Rgs5"),
  Stromal = c("Col1a1", "Acta2", "Pdgfrb"),
  UE = c("Psca", "Upk3a", "Ly6d"),
  RBC = c("Hbb-bs", "Hba-a1", "Hba-a2")
)

DotPlot(endo, features = unique(unname(unlist(KidneyMarkers)))) + 
  theme(axis.text.x = element_text(angle = 45, hjust =1))

## ---------- Doublet removal ----------
endo_filter <- endo[, !endo$seurat_clusters %in% c("7", 14)]
DefaultAssay(endo) <- "RNA"
endo_filter[["RNA"]] <- split(endo_filter[["RNA"]], f = endo_filter$orig.ident)
endo_filter <- NormalizeData(endo_filter) %>% FindVariableFeatures() %>% 
  ScaleData() %>% RunPCA()
endo_filter <- IntegrateLayers(endo_filter, method = CCAIntegration, new.reduction = "CCA")
endo_filter <- FindNeighbors(endo_filter, reduction = "CCA") %>% 
  FindClusters(resolution = 0.8) 
ElbowPlot(endo_filter, ndims = 50)
endo_filter <- RunUMAP(endo_filter, dims = 1:20, 
                       reduction = "CCA", reduction.name = "umap.cca")
DimPlot(endo_filter, reduction = "umap.cca", label = T)

DotPlot(endo_filter, features = unique(unname(unlist(KidneyMarkers)))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ---------- Marker-based annotation ----------
deg <- FindAllMarkers(endo_filter)
deg$pct.diff <- deg$pct.1 - deg$pct.2
deg_filter <- deg %>% 
  filter(p_val_adj < 0.05 & pct.1 > 0.25 & pct.diff > 0.1) %>% 
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n = 10)
deg_list <- split(deg_filter$gene, deg_filter$cluster)

EC_Type_Markers <- list(
  Pan_EC= c("Pecam1", "Cdh5", "Kdr"),
  Arterial = c("Hey2", "Gata6", "Gper1", "Adh7"),
  Arterial_Transitional = c("Plat", "Lpl", "Gata5", "Serpina3g"),
  Capillary = c("Plvap", "Esm1", "Col15a1", "Vwa1"),
  Venular = c("Aqp1", "Lef1", "Slc14a1", "Hpgd")
)

FacetDotPlot(endo_filter, marker_list = deg_list, group_by = "seurat_clusters")
FacetDotPlot(endo_filter, marker_list = EC_Type_Markers, group_by = "seurat_clusters")

endo_filter$cellsubtype <- recode(
  endo_filter$clusters,
  "0" = "Endo",    
  "1" = "Capillary",  
  "2" = "Venule",     
  "3" = "Arteriole"   
)
plot_scanpy_umap(endo_filter, group_by = "cellsubtype")
# saveRDS(endo_filter, "02_processed/scRNA/inhouse_endo_filter.rds")

## ---------- Add annotation back to sce ----------
sce$celltype_fine <- as.character(sce$celltype_fine)
sce$celltype_fine[sce$celltype_fine == "Endo"] <- 'Unknow'
sce$celltype_fine[colnames(endo_filter)] <- as.character(endo_filter$cellsubtype)


## ---------- Plot ----------
EC_Type_Markers <- list(
  Endo= c("Pecam1", "Cdh5", "Kdr"),
  Arteriole = c("Hey2", "Gata6", "Gper1", "Adh7"),
 Venule = c("Aqp1", "Lef1", "Slc14a1", "Hpgd"),
  Capillary = c("Plvap", "Esm1", "Col15a1", "Vwa1")
)
endo_filter$cellsubtype <- factor(endo_filter$cellsubtype,
                                  levels = c("Endo", "Arteriole", "Venule",
                                             "Capillary"))

p1 <- FacetDotPlot(endo_filter, EC_Type_Markers, group_by = "cellsubtype",
             colors =  c("#E0F3DB", "#CCEBC5", "#A8DDB5", "#7BCCC4", "#43A2CA"))
gs(p1, pic_supp, "SFig4 DotPlot", width = 6.5, height = 3.2)


endocolors = c("Endo" = "#2d6269", "Arteriole" = "#f7c387", "Venule" = "#8acbea",
             "Capillary" = "#d28daf")
p1 <- plot_scanpy_umap(endo_filter, group_by = "cellsubtype", label = F, 
                       reduction = "umap.cca", custom_colors = endocolors)
gs(p1, pic_fig3, "Fig3B UMAP_Endo", width = 5.5, height = 6)

p <- DimPlot(endo_filter, group.by = "orig.ident", alpha = 1, raster = F, reduction = "umap.cca") +
  scale_color_manual(values = sample_colors) +
  theme(
    axis.text = element_blank(),  
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14)
  )+
  labs(x="", y="", title = "")
p
gs(p, pic_supp, "SFig5F Endo_batch", width = 5, height = 5)
