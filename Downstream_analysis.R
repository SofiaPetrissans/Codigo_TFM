library(dplyr)
library(Seurat)
packageVersion("Seurat")
library(patchwork)
library(readr)
library(RColorBrewer)
library(ggplot2)
library(harmony)
library(readxl)
library(pheatmap)
library(cowplot)
library(clusterProfiler)
library(org.Mm.eg.db)  # Base de datos para genes murinos
library(enrichplot)

## -- C A R G A R - D A T O S --------------------------------------------------
# Después de aplicar todos los pasos de seurat:

data.harmony <- readRDS("Harmony.rds")

path.guardar <- paste("graficos")

data.harmony$Harmony_Log_res.0.7 <- factor(
  data.harmony$Harmony_Log_res.0.7,
  levels = sort(as.numeric(levels(as.factor(data.harmony$Harmony_Log_res.0.7))))
)
data.harmony$Diet <- factor(data.harmony$Diet, levels = c("NCD", "4d_HFD", "2w_HFD"))

# Colores monis
getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))
col <-  getPalette(30)
col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.0.7)))

# -- DimPlot -------------------------------------------------------------------

DimPlot(data.harmony, group.by = "Harmony_Log_res.0.3", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"DimPlot.Harmony.0.3.png",sep="/"),width=10,height=10)

DimPlot(data.harmony, group.by = "Harmony_Log_res.0.5", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"DimPlot.Harmony.0.5.png",sep="/"),width=10,height=10)

DimPlot(data.harmony, group.by = "Harmony_Log_res.0.7", pt.size = 1,label = T,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"DimPlot.Harmony.0.7.png",sep="/"), plot=DimPlot.Harmony.0.7,width=10,height=10)


# -- Feature Plot --------------------------------------------------------------
FeaturePlot(data.harmony, 
            reduction = "umap",
            features = c("Bmx", "Hey1", "Gja5", "Sema3g",
                         "Vwf", "Vcam1", "Selp", "Prss23",
                         "Apln", "Kdr", "Rgcc", "Igfbp7",
                         "Mki67", "Pcna" 
                         ), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = FALSE) & NoAxes()
ggsave(filename=paste(path.guardar,"FeaturePlot_Res0.7.png",sep="/"),width=14,height=10)


# -- Layers --------------------------------------------------------------------


data.harmony <- SetIdent(data.harmony, value="Harmony_Log_res.0.7")
# Capa 1: Nombres generales
layer1 <- c("Capilar", "Vena", "Vena", "Capilar", "Arteria", 
            "Capilar", "Capilar", "Capilar", "Arteria", "Vena", "Vena")

# Capa 2: ingles
layer2 <-c("Capillary", "Vein", "Vein", "Capillary", "Artery", 
           "Capillary", "Capillary", "Capillary", "Artery", "Vein", "Vein")

# Capa 3: Nombres generales + subCap
layer3 <- c("Capilar Tipo I", "Vena", "Vena", "Capilar Tipo II", "Arteria", 
            "Capilar Tipo III", "Capilar Tipo IV", "Capilar Tipo V", "Arteria", "Vena", "Vena")

layer4 <- c("Capillary Type I", "Vein", "Vein", "Capillary Type II", "Artery",
            "Capillary Type III", "Capillary Type IV", "Capillary Type V", "Artery", "Vein", "Vein")
# Asignar los niveles de clusters al nuevo identificador
names(layer1) <- levels(data.harmony)  # Clusters originales en el Seurat
names(layer2) <- levels(data.harmony)
names(layer3) <- levels(data.harmony)
names(layer4) <- levels(data.harmony)

# Crear nuevas columnas en el metadata con los nombres personalizados
data.harmony@meta.data$Cluster_Layer1 <- factor(layer1[Idents(data.harmony)])
data.harmony@meta.data$Cluster_Layer2 <- factor(layer2[Idents(data.harmony)])
data.harmony@meta.data$Cluster_Layer3 <- factor(layer3[Idents(data.harmony)])
data.harmony@meta.data$Cluster_Layer4 <- factor(layer4[Idents(data.harmony)])

annotation.colors.L1 <- c(
  "Arteria" = "#97D489", 
  "Vena" = "#73A4CA", 
  "Capilar" = "#FCA08D"
)

annotation.colors.L2 <- c(
  "Artery" = "#97D489", 
  "Vein" = "#73A4CA", 
  "Capillary" = "#FCA08D"
)

annotation.colors.L3 <- c(
  "Arteria" = "#97D489", 
  "Vena" = "#73A4CA", 
  "Capilar Tipo I" = "#FFDAB9",
  "Capilar Tipo II" = "#F9B478",
  "Capilar Tipo III" = "#DD4636",
  "Capilar Tipo IV" = "#FF8C69",
  "Capilar Tipo V" = "#CD7054"
)

annotation.colors.L4 <- c(
  "Artery" = "#97D489",
  "Vena" = "#73A4CA",
  "Capillary Type I" = "#FFDAB9",
  "Capillary Type II" = "#F9B478",
  "Capillary Type III" = "#DD4636",
  "Capillary Type IV" = "#FF8C69",
  "Capillary Type V" = "#CD7054"
)

# -- DimPlots Layers -----------------------------------------------------------
## -- Layer 1 ------------------------------------------------------------------
DimPlot(data.harmony, reduction = "umap", group.by = "Cluster_Layer1", 
        label = FALSE, pt.size = 0.3, raster=FALSE)+
  scale_color_manual(values = annotation.colors.L1) +
  NoAxes() + 
  ggtitle("Cell Type") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 7)) +  
  guides(colour = guide_legend(override.aes = list(size = 2))) 

ggsave(filename = paste(path.guardar,"General_Res_07_UMAP.png",sep="/"),width = 6,height = 6)

DimPlot(data.harmony, reduction = "umap", group.by = "Cluster_Layer1",
        split.by = "Phenotype",
        label = FALSE, pt.size = 0.3, raster=FALSE)+
  scale_color_manual(values = annotation.colors.L1) +
  NoAxes() + 
  ggtitle("Cell Type") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuousfr(expand = expansion(mult = 0.15)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 7)) +  
  guides(colour = guide_legend(override.aes = list(size = 2))) 
ggsave(filename = paste(path.guardar,"General_Res_07_Split_Pheno.png",sep="/"),width = 8,height = 6)

DimPlot(data.harmony, reduction = "umap", group.by = "Cluster_Layer1",
        split.by = "Diet",
        label = FALSE, pt.size = 0.3, raster=FALSE)+
  scale_color_manual(values = annotation.colors.L1) +
  NoAxes() + 
  ggtitle("Cell Type") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 7)) +  
  guides(colour = guide_legend(override.aes = list(size = 2))) 
ggsave(filename = paste(path.guardar,"General_Res_07_Split_Diet.png",sep="/"),width = 12,height = 6)


## -- Layer 3 ------------------------------------------------------------------

## Dim`Plot clusters`
DimPlot(data.harmony, reduction = "umap", group.by = "Cluster_Layer3", 
        label = FALSE, pt.size = 0.3, raster=FALSE)+
  scale_color_manual(values = annotation.colors.L3) +
  NoAxes() + 
  ggtitle("Cell Type") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 7)) +  
  guides(colour = guide_legend(override.aes = list(size = 2))) 

ggsave(filename = paste(path.guardar,"SubCap_Res_07_UMAP.png",sep="/"),width = 6,height = 6)

## DimPlot dividido por fenotipos
DimPlot(data.harmony, reduction = "umap", group.by = "Cluster_Layer3",
        split.by = "Phenotype",
        label = FALSE, pt.size = 0.3, raster=FALSE)+
  scale_color_manual(values = annotation.colors.L3) +
  NoAxes() + 
  ggtitle("Cell Type") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 7)) +  
  guides(colour = guide_legend(override.aes = list(size = 2))) 
ggsave(filename = paste(path.guardar,"SubCap_Res_07_Split_Pheno.png",sep="/"),width = 8,height = 6)

## DimPlot dividio por dieta
DimPlot(data.harmony, reduction = "umap", group.by = "Cluster_Layer3",
        split.by = "Diet",
        label = FALSE, pt.size = 0.3, raster=FALSE)+
  scale_color_manual(values = annotation.colors.L3) +
  NoAxes() + 
  ggtitle("Cell Type") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 7)) +  
  guides(colour = guide_legend(override.aes = list(size = 2))) 
ggsave(filename = paste(path.guardar,"SubCap_Res_07_Split_Diet.png",sep="/"),width = 12,height = 6)

## DimPlot dividido por fenotipo y dieta
DimPlot(data.harmony, reduction = "umap", group.by = "Cluster_Layer3", 
        label = FALSE, pt.size = 0.3, raster = FALSE) +
  scale_color_manual(values = annotation.colors.L3) +
  NoAxes() +
  ggtitle("Tipo Celular por Dieta y Genotipo") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 7),
        strip.text = element_text(size = 9, face = "bold")) +  
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  facet_grid(rows = vars(data.harmony$Phenotype), cols = vars(data.harmony$Diet))
ggsave(filename = paste(path.guardar,"cell_type_pheno_diet.png",sep="/"),width = 12,height = 8)

# -- HeatMap -------------------------------------------------------------------

# Lista de genes específicos
marker_genes <- c("Bmx", "Hey1", "Gja5", "Sema3g",
                  "Vwf", "Vcam1", "Selp", "Prss23",
                  "Apln", "Kdr", "Rgcc", "Igfbp7",
                  "Mki67", "Pcna")
data.harmony <- SetIdent(data.harmony, value="Cluster_Layer1")
# Calcula la expresión promedio para cada clúster
avg.exp <- AverageExpression(data.harmony, return.seurat = FALSE)
# Extrae la matriz de expresión promedio
expression.matrix <- avg.exp$RNA 
expression.matrix <- expression.matrix[marker_genes, ]
# Escala los valores (z-score por fila)
expression.matrix.scaled <- t(scale(t(expression.matrix)))

# Define los grupos de genes y crea una anotación
gene.groups <- data.frame(
  Tipos_Celulares = c(rep("Arteria", 4), rep("Vena", 4), rep("Capilar", 6))
)
rownames(gene.groups) <- marker_genes  # Asigna los nombres de los genes a las filas de la anotación

# Extraer los nombres de los clústeres únicos
cluster_labels <- unique(data.harmony@meta.data$Cluster_Layer1)
annotation_cols <- data.frame(
  Cluster = cluster_labels
)
rownames(annotation_cols) <- colnames(expression.matrix)

annotation.colors <- list(
  Tipos_Celulares = c(
    "Arteria" = "#97D489", 
    "Vena" = "#73A4CA", 
    "Capilar" = "#FCA08D"
  )
)

# Genera el heatmap con la anotación
HP_marker_genes <- pheatmap(expression.matrix.scaled,
                            cluster_rows = FALSE,       
                            cluster_cols = FALSE,
                            group_by = "Cluster_Layer1",
                            color = colorRampPalette(c("#6495ED", "white", "#A52A2A"))(100),
                            main = "Heatmap marcadores específicos del tipo celular",
                            annotation_row = gene.groups,      # Anotación de las filas
                            annotation_colors = annotation.colors)
print(HP_marker_genes)
ggsave(filename=paste(path.guardar,"HeatMap_L1_MarkerGenes.png",sep="/"),plot=HP_marker_genes,width=8,height=5)


# -- Enrichment ----------------------------------------------------------------
## -- S U B C A P I L L A R Y --------------------------------------------------

data.harmony <- SetIdent(data.harmony, value="Harmony_Log_res.0.7")
# Filtrar para obtener los clústeres positivos para capilares
capillary_clusters <- WhichCells(data.harmony, ident = c("0", "3", "5", "6", "7"))  # Usa los nombres de tus clústeres
data_capillary <- subset(data.harmony, cells = capillary_clusters)

data_capillary.markers <- FindAllMarkers(data_capillary, only.pos = TRUE)
data_capillary.markers %>%
  group_by(cluster)
# Guardar en tabla de excel
write.csv(data_capillary.markers, "Markers_Harmony_Capillary_Res_07.csv")

cluster_data <- read_delim("Markers_Harmony_Capillary_Res_07.csv", 
                           delim = ",", escape_double = FALSE, trim_ws = TRUE)
gene_list <- split(cluster_data$gene, cluster_data$cluster)

# Extraer los identificadores de los genes relevantes
genes_markers <- data_capillary.markers %>%
  filter(p_val_adj < 0.05) %>%
  dplyr::select(gene)

# Realizar la conversión de símbolos de genes a ENTREZID
genes_markers_converted <- bitr(genes_markers$gene, 
                                fromType = "SYMBOL", 
                                toType = "ENTREZID", 
                                OrgDb = org.Mm.eg.db)

# Ahora solo seleccionamos los genes con un ENTREZID válido
genes_markers_valid <- genes_markers %>%
  filter(gene %in% genes_markers_converted$SYMBOL)

# Agregar la columna 'entrez_id' con los valores correspondientes
genes_markers_valid$entrez_id <- genes_markers_converted$ENTREZID[match(genes_markers_valid$gene, genes_markers_converted$SYMBOL)]

# Realizar el análisis de enriquecimiento de GO para los genes de interés
GO_results <- lapply(gene_list, function(genes){
  enrichGO(
    gene = genes,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
})

annotation.colors.L6 <- c(
  "Capilar Tipo II" = "#F9B478",
  "Capilar Tipo III" = "#DD4636",
  "Capilar Tipo I" = "#ffc3b6",
  "Capilar Tipo IV" = "#FF8C69",
  "Capilar Tipo V" = "#CD7054"
)

go_capil_typeI <- barplot(GO_results[[1]], showCategory = 10) + theme(
  axis.text.y = element_text(size = 14))
go_capil_typeII <- barplot(GO_results[[2]], showCategory = 10) + theme(
  axis.text.y = element_text(size = 14)) # Muestra solo los 10 términos más enriquecidos
go_capil_typeIII <- barplot(GO_results[[3]], showCategory = 10) + theme(
  axis.text.y = element_text(size = 14)) # Muestra solo los 10 términos más enriquecidos
go_capil_typeIV <- barplot(GO_results[[4]], showCategory = 10) + theme(
  axis.text.y = element_text(size = 14)) # Muestra solo los 10 términos más enriquecidos
go_capil_typeV <- barplot(GO_results[[5]], showCategory = 10) + theme(
  axis.text.y = element_text(size = 14))

library(patchwork)

# Combina los gráficos en una única imagen
combined_plot <- (go_capil_typeI | go_capil_typeII) / (go_capil_typeIII | go_capil_typeIV) / (go_capil_typeV)

# Muestra la imagen combinada
print(combined_plot)
ggsave(filename=paste(path.guardar,"GOTermsCap_0_5_6_7.png",sep="/"),plot=combined_plot,width=22,height=15)


# Crear un data frame combinando los resultados de cada cluster
GO_results_df <- do.call(rbind, lapply(names(GO_results), function(cluster) {
  result <- GO_results[[cluster]]
  if (!is.null(result)) { # verificar si hay resultados para el cluster
    data <- as.data.frame(result) # Convertir a data.frame
    data$Cluster <- cluster # agregar columnas para identificar el cluster
    return(data)
    
  }
}))
# Guardar los resultados del enriquecimiento en un archivo CSV
write.csv(GO_results_df, "GO_Enrichment_Capillary_Res.0.7.csv")

# -- GO TERMS ELEGIDOS ---------------------------------------------------------

go_terms <- read_delim("go_terms.csv", delim = ";", 
                       escape_double = FALSE, trim_ws = TRUE)

annotation.colors.L6 <- c(
  "Capilar Tipo I" = "#FFDAB9",
  "Capilar Tipo II" = "#F9B478",
  "Capilar Tipo III" = "#DD4636",
  "Capilar Tipo IV" = "#FF8C69",
  "Capilar Tipo V" = "#CD7054"
)

go_terms$NegativeLogFoldChange <- -log(go_terms$p.adjust)
go_terms$y<-0
ggplot(data=go_terms, aes(x=ID, y=NegativeLogFoldChange, fill=ClusterName)) +
  geom_bar(stat="identity",width = 0.6, position=position_dodge(), alpha = 0.7)+
  scale_fill_manual(values = annotation.colors.L6)+
  geom_text(aes(label = Description, y=0),
            color="black",
            size=4,
            fontface="bold",
            hjust=ifelse(go_terms$y < go_terms$NegativeLogFoldChange / 0, -0.02, 1.03))+
  theme(axis.line = element_line(colour = "black", linewidth = 0.5, linetype = "solid"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        axis.text.x = element_text(color="black",size = 10),
        axis.text.y = element_text(color="black",size = 10),
        strip.text = element_text(color = "black",face="bold",size=14),
        panel.background = element_rect(fill = "white", color = NA),
        strip.background = element_rect(fill = "#D5D5D5", color = NA))+
  labs(title="Up Regulated GO Terms",x="", y = expression(log[2](p-value)))+
  theme(plot.title = element_text(face = "bold",hjust = 0.5,size=20))+
  scale_y_continuous(position = "right")+
  coord_flip()+
  facet_wrap(~ ClusterName,scale="free")
ggsave(paste(path.guardar,"BarPlot_GoTerms.png",sep="/"),width=22,height=15)

# -- CONTROL -------------------------------------------------------------------

data_control <- subset(data.harmony, subset = Phenotype == "Control")

## -- Umaps --------------------------------------------------------------------
DimPlot(data_control, reduction = "umap", group.by = "Cluster_Layer3",
        label = FALSE, pt.size = 0.3, raster=FALSE)+
  scale_color_manual(values = annotation.colors.L3) +
  NoAxes() + 
  ggtitle("Control") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 7)) +  
  guides(colour = guide_legend(override.aes = list(size = 2))) 
ggsave(filename = paste(path.guardar,"SubCap_Res_07_Ctrl.png",sep="/"),width = 6,height = 6)

DimPlot(data_control, reduction = "umap", group.by = "Cluster_Layer3",
        split.by = "Diet",
        label = FALSE, pt.size = 0.3, raster=FALSE)+
  scale_color_manual(values = annotation.colors.L3) +
  NoAxes() + 
  ggtitle("Control") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 7)) +  
  guides(colour = guide_legend(override.aes = list(size = 2))) 
ggsave(filename = paste(path.guardar,"SubCap_Res_07_Split_Diet_Ctrl.png",sep="/"),width = 12,height = 6)

## -- Proporciones -------------------------------------------------------------

metadata <- data_control@meta.data
# Calcula porcentajes por condición y tipo de célula
cell_percentage <- metadata %>%
  group_by(Diet, Cluster_Layer1) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(Diet) %>%
  mutate(percentage = (cell_count / sum(cell_count)) * 100)

# Gráfico de barras apilado
ggplot(cell_percentage, aes(x = Diet, y = percentage, fill = Cluster_Layer1)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_fill_manual(values = annotation.colors.L1) +
  labs(title = "Porcentaje de células por dieta", x = "Diet", y = "% células por cluster", fill = "Tipo Celular") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = paste(path.guardar,"BP_L1_TypeCell_Diet_Ctrl.png",sep="/"),width = 5,height = 7)

# Gráfico de barras apilado NCD vs. 4dHFD

cell_percentage_filtered <- cell_percentage %>%
  filter(Diet %in% c("NCD", "4d_HFD"))

ggplot(cell_percentage_filtered, aes(x = Diet, y = percentage, fill = Cluster_Layer1)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_fill_manual(values = annotation.colors.L1) +
  labs(title = "NCD vs. 4d HFD", x = "Diet", y = "% células por cluster", fill = "Tipo Celular") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = paste(path.guardar,"BP_L1_TypeCell_Ctrl_NCDvs4d.png",sep="/"),width = 5,height = 7)

# Gráfico de barras apilado NCD vs. 2wHFD

cell_percentage_filtered <- cell_percentage %>%
  filter(Diet %in% c("NCD", "2w_HFD"))

ggplot(cell_percentage_filtered, aes(x = Diet, y = percentage, fill = Cluster_Layer1)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_fill_manual(values = annotation.colors.L1) +
  labs(title = "NCD vs. 2w HFD", x = "Diet", y = "% células por cluster", fill = "Tipo Celular") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = paste(path.guardar,"BP_L1_TypeCell_Ctrl_NCDvs2w.png",sep="/"),width = 5,height = 7)

# Gráfico de lineas

cell_counts <- metadata %>%
  group_by(Diet, Cluster_Layer1) %>%
  summarise(cell_count = n(), .groups = "drop")

ggplot(cell_counts, aes(x = Diet, y = cell_count, group = Cluster_Layer1, color = Cluster_Layer1)) +
  geom_line(linewidth = 0.3) +  # Añade líneas para conectar los puntos
  geom_point(size = 2) +  # Añade los puntos en cada dieta
  labs(
    title = "Cambios en tipos celulares a lo largo de las dietas",
    x = "Dietas",
    y = "Cantidad de células",
    color = "Tipo Celular"
  ) + 
  scale_color_manual(values = annotation.colors.L1) +  # Usar la paleta personalizada
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

ggsave(filename = paste(path.guardar,"LinePlot_L1_diet_ctrl.png",sep="/"),width = 7,height = 4)


# -- Mutante -------------------------------------------------------------------
data_mut <- subset(data.harmony, subset = Phenotype == "Mutant")

## -- Umaps --------------------------------------------------------------------

DimPlot(data_mut, reduction = "umap", group.by = "Cluster_Layer3",
        label = FALSE, pt.size = 0.3, raster=FALSE)+
  scale_color_manual(values = annotation.colors.L3) +
  NoAxes() + 
  ggtitle("Mutante") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 7)) +  
  guides(colour = guide_legend(override.aes = list(size = 2))) 
ggsave(filename = paste(path.guardar,"SubCap_Res_07_Mut.png",sep="/"),width = 6,height = 6)

DimPlot(data_mut, reduction = "umap", group.by = "Cluster_Layer3",
        split.by = "Diet",
        label = FALSE, pt.size = 0.3, raster=FALSE)+
  scale_color_manual(values = annotation.colors.L3) +
  NoAxes() + 
  ggtitle("Mutante") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 7)) +  
  guides(colour = guide_legend(override.aes = list(size = 2))) 
ggsave(filename = paste(path.guardar,"SubCap_Res_07_Split_Diet_Mut.png",sep="/"),width = 12,height = 6)

## -- Proporciones -------------------------------------------------------------
metadata <- data_mut@meta.data
# Calcula porcentajes por condición y tipo de célula
cell_percentage <- metadata %>%
  group_by(Diet, Cluster_Layer3) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(Diet) %>%
  mutate(percentage = (cell_count / sum(cell_count)) * 100)

# Gráfico de barras apilado
barplot_percentage <- ggplot(cell_percentage, aes(x = Diet, y = percentage, fill = Cluster_Layer3)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_fill_manual(values = annotation.colors.L3) +
  labs(title = "Porcentaje de células por dieta", x = "Diet", y = "% células por cluster", fill = "Tipo Celular") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = paste(path.guardar,"BP_L3_TypeCell_Diet_Mut.png",sep="/"),width = 5,height = 7)

# Gráfico de barras apilado NCD vs. 4dHFD

cell_percentage_filtered <- cell_percentage %>%
  filter(Diet %in% c("NCD", "4d_HFD"))

barplot_percentage <- ggplot(cell_percentage_filtered, aes(x = Diet, y = percentage, fill = Cluster_Layer3)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_fill_manual(values = annotation.colors.L3) +
  labs(title = "NCD vs. 4d HFD", x = "Diet", y = "% células por cluster", fill = "Tipo Celular") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = paste(path.guardar,"BP_L3_TypeCell_Mut_NCDvs4d.png",sep="/"),width = 5,height = 7)

# Gráfico de barras apilado NCD vs. 2wHFD

cell_percentage_filtered <- cell_percentage %>%
  filter(Diet %in% c("NCD", "2w_HFD"))

barplot_percentage <- ggplot(cell_percentage_filtered, aes(x = Diet, y = percentage, fill = Cluster_Layer3)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_fill_manual(values = annotation.colors.L3) +
  labs(title = "NCD vs. 2w HFD", x = "Diet", y = "% células por cluster", fill = "Tipo Celular") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = paste(path.guardar,"BP_L3_TypeCell_Mut_NCDvs2w.png",sep="/"),width = 5,height = 7)

# Gráfico de lineas

cell_counts <- metadata %>%
  group_by(Diet, Cluster_Layer1) %>%
  summarise(cell_count = n(), .groups = "drop")

ggplot(cell_counts, aes(x = Diet, y = cell_count, group = Cluster_Layer1, color = Cluster_Layer1)) +
  geom_line(linewidth = 0.3) +  # Añade líneas para conectar los puntos
  geom_point(size = 2) +  # Añade los puntos en cada dieta
  labs(
    title = "Cambios en tipos celulares a lo largo de las dietas",
    x = "Dietas",
    y = "Cantidad de células",
    color = "Tipo Celular"
  ) + 
  scale_color_manual(values = annotation.colors.L1) +  # Usar la paleta personalizada
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

ggsave(filename = paste(path.guardar,"LinePlot_L1_diet_Mut.png",sep="/"),width = 7,height = 4)


# -- scProportion --------------------------------------------------------------

## -- Control ------------------------------------------------------------------

library("scProportionTest")
# Crear el objeto de análisis para el fenotipo "Control"
prop_test <- sc_utils(data_control)
# Configurar el análisis de permutación en la condición "Control"
prop_test <- permutation_test(
  prop_test,
  cluster_identity = "Cluster_Layer3",  # Tipo celular
  sample_1 = "NCD",                    # Condición 1 (dieta control)
  sample_2 = "4d_HFD",                 # Condición 2 (dieta experimental)
  sample_identity = "Diet"             # Identificación de las dietas
)

permutation_plot(prop_test)
ggsave(filename = paste(path.guardar,"PermutationPlot_NCDvs4d_Ctrl.png",sep="/"),width = 10,height = 7)

# Convertir `results` a data frame (si no lo es ya)
results_df <- as.data.frame(prop_test@results)

ggplot(results_df, aes(x = permutation.clusters, y = permutation.boot_mean_log2FD, 
                                             fill = permutation.FDR < 0.05)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(
    ymin = permutation.boot_CI_2.5,
    ymax = permutation.boot_CI_97.5
  ), width = 0.3, color = "black") +
  scale_fill_manual(values = c("TRUE" = "#FA8072", "FALSE" = "gray"))+
  labs(
    title = "Log2 Fold Change entre NCD y 4d_HFD",
    x = "Tipo celular",
    y = "Log2 Fold Change",
    fill = "FDR < 0.05"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = paste(path.guardar,"BP_scProp_NCDvs4d_Ctrl.png",sep="/"),width = 5,height = 7)

# NCD vs. 2wHFD

prop_test <- sc_utils(data_control)
# Configurar el análisis de permutación en la condición "Control"
prop_test <- permutation_test(
  prop_test,
  cluster_identity = "Cluster_Layer3",  
  sample_1 = "NCD",                    
  sample_2 = "2w_HFD",                
  sample_identity = "Diet"        
)

permutation_plot(prop_test)
ggsave(filename = paste(path.guardar,"PermutationPlot_NCDvs2w_Ctrl.png",sep="/"),width = 10,height = 7)

# Convertir `results` a data frame (si no lo es ya)
results_df <- as.data.frame(prop_test@results)

ggplot(results_df, aes(x = permutation.clusters, y = permutation.boot_mean_log2FD, 
                       fill = permutation.FDR < 0.05)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(
    ymin = permutation.boot_CI_2.5,
    ymax = permutation.boot_CI_97.5
  ), width = 0.3, color = "black") +
  scale_fill_manual(values = c("TRUE" = "#FA8072", "FALSE" = "gray"))+
  labs(
    title = "Log2 Fold Change entre NCD y 2w_HFD",
    x = "Tipo celular",
    y = "Log2 Fold Change",
    fill = "FDR < 0.05"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = paste(path.guardar,"BP_scProp_NCDvs2w_Ctrl.png",sep="/"),width = 5,height = 7)


## -- Mutante ------------------------------------------------------------------

library("scProportionTest")
# Crear el objeto de análisis para el fenotipo "Control"
prop_test <- sc_utils(data_mut)
# Configurar el análisis de permutación en la condición "Control"
prop_test <- permutation_test(
  prop_test,
  cluster_identity = "Cluster_Layer3",  # Tipo celular
  sample_1 = "NCD",                    # Condición 1 (dieta control)
  sample_2 = "4d_HFD",                 # Condición 2 (dieta experimental)
  sample_identity = "Diet"             # Identificación de las dietas
)

permutation_plot(prop_test)
ggsave(filename = paste(path.guardar,"PermutationPlot_NCDvs4d_Mut.png",sep="/"),width = 10,height = 7)

# Convertir `results` a data frame (si no lo es ya)
results_df <- as.data.frame(prop_test@results)

ggplot(results_df, aes(x = permutation.clusters, y = permutation.boot_mean_log2FD, 
                       fill = permutation.FDR < 0.05)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(
    ymin = permutation.boot_CI_2.5,
    ymax = permutation.boot_CI_97.5
  ), width = 0.3, color = "black") +
  scale_fill_manual(values = c("TRUE" = "#FA8072", "FALSE" = "gray"))+
  labs(
    title = "Log2 Fold Change entre NCD y 4d_HFD",
    x = "Tipo celular",
    y = "Log2 Fold Change",
    fill = "FDR < 0.05"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = paste(path.guardar,"BP_scProp_NCDvs4d_Mut.png",sep="/"),width = 5,height = 7)

# NCD vs. 2wHFD

prop_test <- sc_utils(data_mut)
# Configurar el análisis de permutación en la condición "Control"
prop_test <- permutation_test(
  prop_test,
  cluster_identity = "Cluster_Layer3",  
  sample_1 = "NCD",                    
  sample_2 = "2w_HFD",                
  sample_identity = "Diet"        
)

permutation_plot(prop_test)
ggsave(filename = paste(path.guardar,"PermutationPlot_NCDvs2w_Mut.png",sep="/"),width = 10,height = 7)

# Convertir `results` a data frame (si no lo es ya)
results_df <- as.data.frame(prop_test@results)

ggplot(results_df, aes(x = permutation.clusters, y = permutation.boot_mean_log2FD, 
                       fill = permutation.FDR < 0.05)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(
    ymin = permutation.boot_CI_2.5,
    ymax = permutation.boot_CI_97.5
  ), width = 0.3, color = "black") +
  scale_fill_manual(values = c("TRUE" = "#FA8072", "FALSE" = "gray"))+
  labs(
    title = "Log2 Fold Change entre NCD y 2w_HFD",
    x = "Tipo celular",
    y = "Log2 Fold Change",
    fill = "FDR < 0.05"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = paste(path.guardar,"BP_scProp_NCDvs2w_Mut.png",sep="/"),width = 5,height = 7)


# -- R O G U E -----------------------------------------------------------------
## -- Control ------------------------------------------------------------------

#if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
#devtools::install_github("PaulingLiu/ROGUE")
library(ROGUE)

expr_matrix <- as.matrix(data_control@assays$RNA@data)
expr <- matr.filter(expr_matrix, min.cells = 10, min.genes = 10)
ent.res <- SE_fun(expr)
head(ent.res)

rogue.value <- CalculateRogue(ent.res, platform = "UMI")


rogue.res <- rogue(expr, labels = data_control@meta.data$Cluster_Layer1, 
                   samples = data_control@meta.data$Diet, platform = "UMI", span = 0.6)


# Convertir a formato largo
rogue.long <- rogue.res %>%
  as.data.frame() %>%
  mutate(Sample = rownames(rogue.res)) %>%
  pivot_longer(cols = -Sample, names_to = "Cluster", values_to = "ROGUE")

# Generar boxplot
ggplot(rogue.long, aes(x = Cluster, y = ROGUE, fill = Cluster)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = annotation.colors.L1)

ggsave(filename=paste(path.guardar,"Rogue_L1_CTRL.png",sep="/"), width=5,height=5)



# Control
expr_matrix_control <- as.matrix(data_control@assays$RNA@data)
expr_control <- matr.filter(expr_matrix_control, min.cells = 10, min.genes = 10)
ent.res_control <- SE_fun(expr_control)
rogue.res_control <- rogue(expr_control, labels = data_control@meta.data$Cluster_Layer3, 
                           samples = data_control@meta.data$Diet, platform = "UMI", span = 0.6)

rogue.long.control <- rogue.res_control %>%
  as.data.frame() %>%
  mutate(Sample = rownames(rogue.res_control)) %>%
  pivot_longer(cols = -Sample, names_to = "Cluster", values_to = "ROGUE") %>%
  mutate(Phenotype = "Control")  # Añadir columna para el fenotipo

ggplot(rogue.long.control, aes(x = Cluster, y = ROGUE, fill = Cluster)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = annotation.colors.L3)

ggsave(filename=paste(path.guardar,"Rogue_L2_Ctrl.png",sep="/"), width=5,height=5)


# Mutant
expr_matrix <- as.matrix(data_mut@assays$RNA@data)
expr <- matr.filter(expr_matrix, min.cells = 10, min.genes = 10)
ent.res <- SE_fun(expr)
head(ent.res)

rogue.value <- CalculateRogue(ent.res, platform = "UMI")


rogue.res <- rogue(expr, labels = data_mut@meta.data$Cluster_Layer1, 
                   samples = data_mut@meta.data$Diet, platform = "UMI", span = 0.6)

library(tidyr)
# Convertir a formato largo
rogue.long <- rogue.res %>%
  as.data.frame() %>%
  mutate(Sample = rownames(rogue.res)) %>%
  pivot_longer(cols = -Sample, names_to = "Cluster", values_to = "ROGUE")

# Generar boxplot
ggplot(rogue.long, aes(x = Cluster, y = ROGUE, fill = Cluster)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = annotation.colors.L1)

ggsave(filename=paste(path.guardar,"Rogue_L1_MUT.png",sep="/"), width=5,height=5)


expr_matrix_mut <- as.matrix(data_mut@assays$RNA@data)
expr_mut <- matr.filter(expr_matrix_mut, min.cells = 10, min.genes = 10)
ent.res_mut <- SE_fun(expr_mut)
rogue.res_mut <- rogue(expr_mut, labels = data_mut@meta.data$Cluster_Layer3, 
                       samples = data_mut@meta.data$Diet, platform = "UMI", span = 0.6)

rogue.long.mut <- rogue.res_mut %>%
  as.data.frame() %>%
  mutate(Sample = rownames(rogue.res_mut)) %>%
  pivot_longer(cols = -Sample, names_to = "Cluster", values_to = "ROGUE") %>%
  mutate(Phenotype = "Mutant")  # Añadir columna para el fenotipo
ggplot(rogue.long.mut, aes(x = Cluster, y = ROGUE, fill = Cluster)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = annotation.colors.L3)

ggsave(filename=paste(path.guardar,"Rogue_L2_Mut.png",sep="/"), width=5,height=5)


# Combinar ambos
rogue.long.combined <- bind_rows(rogue.long.control, rogue.long.mut)

boxplot <- ggplot(rogue.long.combined, aes(x = Cluster, y = ROGUE, fill = Phenotype)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Separar los boxplots por fenotipo
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("Control" = "#B9DDF1", "Mutant" = "#98D688")) +  # Colores para los fenotipos
  labs(title = "Valores ROGUE por Cluster y Fenotipo",
       x = "Cluster",
       y = "ROGUE",
       fill = "Phenotype")

ggsave(filename=paste(path.guardar,"Rogue_Merge.png",sep="/"),plot=boxplot,width=6,height=5)


# -- Pie Charts ----------------------------------------------------------------

data_ctrl_diet <- subset(data_mut, subset = Diet == "NCD")
metadata <- data_ctrl_diet@meta.data

# Usar metadata para calcular el conteo y porcentajes
cell_counts <- metadata %>%
  group_by(Cluster_Layer3) %>%
  summarise(cell_count = n(), .groups = "drop")

# Capilares Tipo I ##############################################################
# Separar las células X del resto
data <- data.frame(
  Categorias = c("Capilar Tipo I", "Otras células"),
  Count = c(
    cell_counts %>% filter(Cluster_Layer3 == "Capilar Tipo I") %>% pull(cell_count),
    cell_counts %>% filter(Cluster_Layer3 != "Capilar Tipo I") %>% summarise(sum(cell_count)) %>% pull()
  )
)

# Calcular porcentajes
data$Percentage <- data$Count / sum(data$Count) * 100
data$Label <- paste0(round(data$Percentage, 1), "%")

# Crear el gráfico de pastel
capI <- ggplot(data, aes(x = "", y = Percentage, fill = Categorias)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() + # Remover ejes para que luzca como un pie chart
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), color = "black", size = 10) +
  scale_fill_manual(values = c("Capilar Tipo I" = "#CD0000", "Otras células" = "#C2C2C2")) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Capilares Tipo II ##############################################################
# Separar las células X del resto
data <- data.frame(
  Categorias = c("Capilar Tipo II", "Otras células"),
  Count = c(
    cell_counts %>% filter(Cluster_Layer3 == "Capilar Tipo II") %>% pull(cell_count),
    cell_counts %>% filter(Cluster_Layer3 != "Capilar Tipo II") %>% summarise(sum(cell_count)) %>% pull()
  )
)
# Calcular porcentajes
data$Percentage <- data$Count / sum(data$Count) * 100
data$Label <- paste0(round(data$Percentage, 1), "%")

# Crear el gráfico de pastel
capII <- ggplot(data, aes(x = "", y = Percentage, fill = Categorias)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() + # Remover ejes para que luzca como un pie chart
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), color = "black", size = 10) +
  scale_fill_manual(values = c("Capilar Tipo II" = "#CD0000", "Otras células" = "#C2C2C2")) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Capilares Tipo III ##############################################################
# Separar las células X del resto
data <- data.frame(
  Categorias = c("Capilar Tipo III", "Otras células"),
  Count = c(
    cell_counts %>% filter(Cluster_Layer3 == "Capilar Tipo III") %>% pull(cell_count),
    cell_counts %>% filter(Cluster_Layer3 != "Capilar Tipo III") %>% summarise(sum(cell_count)) %>% pull()
  )
)

# Calcular porcentajes
data$Percentage <- data$Count / sum(data$Count) * 100
data$Label <- paste0(round(data$Percentage, 1), "%")

# Crear el gráfico de pastel
capIII <- ggplot(data, aes(x = "", y = Percentage, fill = Categorias)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() + # Remover ejes para que luzca como un pie chart
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), color = "black", size = 10) +
  scale_fill_manual(values = c("Capilar Tipo III" = "#CD0000", "Otras células" = "#C2C2C2")) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Capilares Tipo IV ##############################################################
# Separar las células X del resto
data <- data.frame(
  Categorias = c("Capilar Tipo IV", "Otras células"),
  Count = c(
    cell_counts %>% filter(Cluster_Layer3 == "Capilar Tipo IV") %>% pull(cell_count),
    cell_counts %>% filter(Cluster_Layer3 != "Capilar Tipo IV") %>% summarise(sum(cell_count)) %>% pull()
  )
)

# Calcular porcentajes
data$Percentage <- data$Count / sum(data$Count) * 100
data$Label <- paste0(round(data$Percentage, 1), "%")

capIV <- ggplot(data, aes(x = "", y = Percentage, fill = Categorias)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() + # Remover ejes para que luzca como un pie chart
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), color = "black", size = 10) +
  scale_fill_manual(values = c("Capilar Tipo IV" = "#CD0000", "Otras células" = "#C2C2C2")) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Capilares Tipo V ##############################################################
# Separar las células X del resto
data <- data.frame(
  Categorias = c("Capilar Tipo V", "Otras células"),
  Count = c(
    cell_counts %>% filter(Cluster_Layer3 == "Capilar Tipo V") %>% pull(cell_count),
    cell_counts %>% filter(Cluster_Layer3 != "Capilar Tipo V") %>% summarise(sum(cell_count)) %>% pull()
  )
)

# Calcular porcentajes
data$Percentage <- data$Count / sum(data$Count) * 100
data$Label <- paste0(round(data$Percentage, 1), "%")

capV <- ggplot(data, aes(x = "", y = Percentage, fill = Categorias)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() + # Remover ejes para que luzca como un pie chart
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), color = "black", size = 10) +
  scale_fill_manual(values = c("Capilar Tipo V" = "#CD0000", "Otras células" = "#C2C2C2")) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

#library(gridExtra)

# Combina los gráficos uno debajo del otro
combined_plot <- grid.arrange(capI, capII, capIII, capIV, capV, ncol = 1)

ggsave(filename=paste(path.guardar,"PieCharts_NCD_MUT.png",sep="/"),plot=combined_plot,width=5,height=20)

# -- Feature plots Capilares ---------------------------------------------------
## Pruebas
FeaturePlot(data.harmony, 
            reduction = "umap",
            features = c("Abcg1", "Lrp1",
                         "Tsc1", "Ezr" 
            ), 
            order = TRUE,
            min.cutoff = 'q9', 
            label = FALSE) & NoAxes()

## Cap Tipo I
plot1 <- FeaturePlot(data.harmony, reduction = "umap", features = "Abcg1", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot2 <- FeaturePlot(data.harmony, reduction = "umap", features = "Tsc1", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot3 <- FeaturePlot(data.harmony, reduction = "umap", features = "Ezr", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()

combined_plot <- plot1 | plot2 | plot3
ggsave(filename = paste(path.guardar, "FeaturePlot_Cap_TI.png", sep = "/"), plot = combined_plot, width = 14, height = 5)

## Cap Tipo II
plot1 <- FeaturePlot(data.harmony, reduction = "umap", features = "Mki67", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot2 <- FeaturePlot(data.harmony, reduction = "umap", features = "Top2a", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot3 <- FeaturePlot(data.harmony, reduction = "umap", features = "Racgap1", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()

combined_plot <- plot1 | plot2 | plot3
ggsave(filename = paste(path.guardar, "FeaturePlot_Cap_TII.png", sep = "/"), plot = combined_plot, width = 14, height = 5)

## Cap Tipo III
plot1 <- FeaturePlot(data.harmony, reduction = "umap", features = "Pgf", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot2 <- FeaturePlot(data.harmony, reduction = "umap", features = "Cxcr4", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot3 <- FeaturePlot(data.harmony, reduction = "umap", features = "Angpt2", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()

combined_plot <- plot1 | plot2 | plot3
ggsave(filename = paste(path.guardar, "FeaturePlot_Cap_TIII.png", sep = "/"), plot = combined_plot, width = 14, height = 5)

## Cap Tipo IV
plot1 <- FeaturePlot(data.harmony, reduction = "umap", features = "Plvap", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot2 <- FeaturePlot(data.harmony, reduction = "umap", features = "Cxcl10", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()

combined_plot <- plot1 | plot2 
ggsave(filename = paste(path.guardar, "FeaturePlot_Cap_TIV.png", sep = "/"), plot = combined_plot, width = 14, height = 5)

## Cap Tipo V
plot1 <- FeaturePlot(data.harmony, reduction = "umap", features = "Irs2", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot2 <- FeaturePlot(data.harmony, reduction = "umap", features = "Ano1", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot3 <- FeaturePlot(data.harmony, reduction = "umap", features = "Cpe", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()

combined_plot <- plot1 | plot2 | plot3

ggsave(filename = paste(path.guardar, "FeaturePlot_Cap_TV.png", sep = "/"), plot = combined_plot, width = 14, height = 5)

# -- Feature Plots Capilares control vs. mut -----------------------------------

## Cap Tipo I
plot1 <- FeaturePlot(data_control, reduction = "umap", features = "Abcg1", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot2 <- FeaturePlot(data_control, reduction = "umap", features = "Tsc1", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot3 <- FeaturePlot(data_control, reduction = "umap", features = "Ezr", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()

combined_plot <- plot1 | plot2 | plot3
ggsave(filename = paste(path.guardar, "FeaturePlot_Cap_TI_CTRL.png", sep = "/"), plot = combined_plot, width = 14, height = 5)


plot1 <- FeaturePlot(data_mut, reduction = "umap", features = "Abcg1", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot2 <- FeaturePlot(data_mut, reduction = "umap", features = "Tsc1", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot3 <- FeaturePlot(data_mut, reduction = "umap", features = "Ezr", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()

combined_plot <- plot1 | plot2 | plot3
ggsave(filename = paste(path.guardar, "FeaturePlot_Cap_TI_MUT.png", sep = "/"), plot = combined_plot, width = 14, height = 5)

## Cap Tipo II
plot1 <- FeaturePlot(data_control, reduction = "umap", features = "Mki67", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot2 <- FeaturePlot(data_control, reduction = "umap", features = "Top2a", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot3 <- FeaturePlot(data_control, reduction = "umap", features = "Racgap1", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()

combined_plot <- plot1 | plot2 | plot3
ggsave(filename = paste(path.guardar, "FeaturePlot_Cap_TII_CTRL.png", sep = "/"), plot = combined_plot, width = 14, height = 5)


plot1 <- FeaturePlot(data_mut, reduction = "umap", features = "Mki67", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot2 <- FeaturePlot(data_mut, reduction = "umap", features = "Top2a", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot3 <- FeaturePlot(data_mut, reduction = "umap", features = "Racgap1", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()

combined_plot <- plot1 | plot2 | plot3
ggsave(filename = paste(path.guardar, "FeaturePlot_Cap_TII_MUT.png", sep = "/"), plot = combined_plot, width = 14, height = 5)

## Cap Tipo III
plot1 <- FeaturePlot(data_control, reduction = "umap", features = "Pgf", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot2 <- FeaturePlot(data_control, reduction = "umap", features = "Cxcr4", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot3 <- FeaturePlot(data_control, reduction = "umap", features = "Angpt2", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()

combined_plot <- plot1 | plot2 | plot3
ggsave(filename = paste(path.guardar, "FeaturePlot_Cap_TIII_CTRL.png", sep = "/"), plot = combined_plot, width = 14, height = 5)


plot1 <- FeaturePlot(data_mut, reduction = "umap", features = "Pgf", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot2 <- FeaturePlot(data_mut, reduction = "umap", features = "Cxcr4", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot3 <- FeaturePlot(data_mut, reduction = "umap", features = "Angpt2", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()

combined_plot <- plot1 | plot2 | plot3
ggsave(filename = paste(path.guardar, "FeaturePlot_Cap_TIII_MUT.png", sep = "/"), plot = combined_plot, width = 14, height = 5)

## Cap Tipo IV
plot1 <- FeaturePlot(data_control, reduction = "umap", features = "Plvap", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot2 <- FeaturePlot(data_control, reduction = "umap", features = "Cxcl10", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()

combined_plot <- plot1 | plot2 
ggsave(filename = paste(path.guardar, "FeaturePlot_Cap_TIV_CTRL.png", sep = "/"), plot = combined_plot, width = 14, height = 5)


plot1 <- FeaturePlot(data_mut, reduction = "umap", features = "Plvap", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot2 <- FeaturePlot(data_mut, reduction = "umap", features = "Cxcl10", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()

combined_plot <- plot1 | plot2 
ggsave(filename = paste(path.guardar, "FeaturePlot_Cap_TIV_MUT.png", sep = "/"), plot = combined_plot, width = 14, height = 5)

# Cap Tipo V
plot1 <- FeaturePlot(data_control, reduction = "umap", features = "Irs2", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot2 <- FeaturePlot(data_control, reduction = "umap", features = "Ano1", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot3 <- FeaturePlot(data_control, reduction = "umap", features = "Cpe", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()

combined_plot <- plot1 | plot2 | plot3
ggsave(filename = paste(path.guardar, "FeaturePlot_Cap_TV_CTRL.png", sep = "/"), plot = combined_plot, width = 14, height = 5)


plot1 <- FeaturePlot(data_mut, reduction = "umap", features = "Irs2", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot2 <- FeaturePlot(data_mut, reduction = "umap", features = "Ano1", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()
plot3 <- FeaturePlot(data_mut, reduction = "umap", features = "Cpe", order = TRUE, min.cutoff = 'q10', label = FALSE) & NoAxes()

combined_plot <- plot1 | plot2 | plot3
ggsave(filename = paste(path.guardar, "FeaturePlot_Cap_TV_MUT.png", sep = "/"), plot = combined_plot, width = 14, height = 5)


# -- Heat MAP Capilares -------------------------------------------------------
data.harmony <- SetIdent(data.harmony, value="Harmony_Log_res.0.7")
# Filtrar para obtener los clústeres positivos para capilares
capillary_clusters <- WhichCells(data.harmony, ident = c("0", "3", "5", "6", "7"))  # Usa los nombres de tus clústeres
data_capillary <- subset(data.harmony, cells = capillary_clusters)
# Lista de genes específicos
marker_genes <- c("Abcg1", "Tsc1", "Ezr",
                  "Mki67", "Top2a", "Racgap1",
                  "Pgf", "Cxcr4", "Angpt2",
                  "Plvap", "Cxcl10",
                  "Irs2", "Ano1", "Cpe"
)
data_capillary <- SetIdent(data_capillary, value="Cluster_Layer3")
# Calcula la expresión promedio para cada clúster
avg.exp <- AverageExpression(data_capillary, return.seurat = FALSE)
# Extrae la matriz de expresión promedio
expression.matrix <- avg.exp$RNA 
expression.matrix <- expression.matrix[marker_genes, ]
# Escala los valores (z-score por fila)
expression.matrix.scaled <- t(scale(t(expression.matrix)))

# Define los grupos de genes y crea una anotación
gene.groups <- data.frame(
  Tipos_Celulares = c(rep("Capilar Tipo I", 3), rep("Capilar Tipo II", 3), rep("Capilar Tipo III", 3),
                      rep("Capilar Tipo IV", 2), rep("Capilar Tipo V", 3))
)
rownames(gene.groups) <- marker_genes  # Asigna los nombres de los genes a las filas de la anotación

# Extraer los nombres de los clústeres únicos
cluster_labels <- unique(data_capillary@meta.data$Cluster_Layer3)
annotation_cols <- data.frame(
  Cluster = cluster_labels
)
rownames(annotation_cols) <- colnames(expression.matrix)

annotation.colors <- list(
  Tipos_Celulares = c(
    "Capilar Tipo I" = "#FFDAB9",
    "Capilar Tipo II" = "#F9B478",
    "Capilar Tipo III" = "#DD4636",
    "Capilar Tipo IV" = "#FF8C69",
    "Capilar Tipo V" = "#CD7054"
  )
)
# Genera el heatmap con la anotación
pheatmap(expression.matrix.scaled,
         cluster_rows = FALSE,       
         cluster_cols = FALSE,
         group_by = "Cluster_Layer3",
         color = colorRampPalette(c("#6495ED", "white", "#A52A2A"))(100),
         main = "Heatmap marcadores específicos del tipo celular",
         annotation_row = gene.groups,
         annotation_colors = annotation.colors
)
ggsave(filename=paste(path.guardar,"HeatMap_Capillary_MarkerGenes.png",sep="/"),width=8,height=5)


# -- UMAP sin integrar por dieta y fenotipo
data.harmony <- RunPCA(data.harmony, features = VariableFeatures(object = data.harmony))

# Función para facilitar la obtención del PCA:

EstimatePCA.Dims<-function(data){
  # Estimate the PCA Dimmensions.
  pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  dim.final <- min(co1, co2)
  return(dim.final)
}

dim.max <- EstimatePCA.Dims(data.harmony)


# RunUMAP
data.harmony <- RunUMAP(data.harmony, dims = 1:dim.max) 

# Extraer las coordenadas UMAP
umap_coordinates <- Embeddings(data.harmony, "umap")

# Añadir estas coordenadas a la metadata de Seurat
data.harmony$UMAP_1 <- umap_coordinates[, 1]
data.harmony$UMAP_2 <- umap_coordinates[, 2]

library(ggplot2)

# Crear un dataframe con las coordenadas UMAP y las variables adicionales (dieta y fenotipo)
umap_df <- data.frame(
  UMAP_1 = data.harmony$UMAP_1,
  UMAP_2 = data.harmony$UMAP_2,
  diet = data.harmony$Diet,  # Asegúrate de que "diet" esté en la metadata
  pheno = data.harmony$Phenotype,  # Asegúrate de que "pheno" esté en la metadata
  cell_type = data.harmony$Cluster_Layer3  # Asegúrate de que "cell_type" esté en la metadata
)

# Crear el gráfico con ggplot
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ diet + pheno) +  # Facet por dieta y fenotipo
  theme_minimal() +
  labs(title = "UMAP de células endoteliales por dieta y fenotipo") +
  theme(legend.position = "right")

ggsave(filename=paste(path.guardar,"UMAP_sin_integrar.png",sep="/"),width=14,height=10)


## -- C E L L - C I C L I N G - S C O R I N G ----------------------------------

## CONTROL

data_control <- SetIdent(data_control, value="Harmony_Log_res.0.7")
# Filtrar para obtener los clústeres positivos para capilares
capillary_clusters <- WhichCells(data_control, ident = c("0", "3", "5", "6", "7"))  # Usa los nombres de tus clústeres
data_capillary <- subset(data_control, cells = capillary_clusters)
# Asignar puntuaciones de ciclo celular
# Seurat incluye listas de genes predeterminadas para las fases S y G2/M
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

data_capillary <- CellCycleScoring(
  object = data_capillary,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE
)


# Visualizar la distribución de las fases del ciclo celular
phase_distribution <- data_capillary@meta.data %>%
  group_by(Phase, Cluster_Layer3) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count))

ggplot(phase_distribution, aes(x = Phase, y = Proportion, fill = Cluster_Layer3)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_fill_manual(values = annotation.colors.L6) +
  labs(title = "Proporción de células por Fase - Control", x = "Fase", fill = "Tipo Celular") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste(path.guardar,"CellPhase_Capillary_CNTRL.png",sep="/"),width = 5,height = 7)


## MUTANTE

data_mut <- SetIdent(data_mut, value="Harmony_Log_res.0.7")
# Filtrar para obtener los clústeres positivos para capilares
capillary_clusters <- WhichCells(data_mut, ident = c("0", "3", "5", "6", "7"))  # Usa los nombres de tus clústeres
data_capillary <- subset(data_mut, cells = capillary_clusters)
# Asignar puntuaciones de ciclo celular
# Seurat incluye listas de genes predeterminadas para las fases S y G2/M
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

data_capillary <- CellCycleScoring(
  object = data_capillary,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE
)

# Visualizar la distribución de las fases del ciclo celular
phase_distribution <- data_capillary@meta.data %>%
  group_by(Phase, Cluster_Layer3) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count))

ggplot(phase_distribution, aes(x = Phase, y = Proportion, fill = Cluster_Layer3)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_fill_manual(values = annotation.colors.L6) +
  labs(title = "Proporción de células por Fase - Mutante", x = "Fase", fill = "Tipo Celular") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste(path.guardar,"CellPhase_Capillary_MUT.png",sep="/"),width = 5,height = 7)

