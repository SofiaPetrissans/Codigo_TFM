# Dividir el objeto Seurat por dieta
seurat_list <- SplitObject(data_control, split.by = "Diet")

# Los resultados serán tres objetos: uno para cada dieta
control_obj <- seurat_list[["NCD"]]
hfd_4d_obj <- seurat_list[["4d_HFD"]]
hfd_2w_obj <- seurat_list[["2w_HFD"]]

# Verifica los niveles de las identidades
levels(Idents(control_obj)) # Debe incluir "Capilar Tipo II", "Capilar Tipo I", etc.

# DEG para Control
deg_control <- FindMarkers(control_obj, ident.1 = "3", ident.2 = NULL)

# DEG para 4 días de HFD
deg_hfd_4d <- FindMarkers(hfd_4d_obj, ident.1 = "3", ident.2 = NULL)

# DEG para 2 semanas de HFD
deg_hfd_2w <- FindMarkers(hfd_2w_obj, ident.1 = "3", ident.2 = NULL)
# Guardar las listas de DEG
write.csv(deg_control, "deg_control.csv", row.names = TRUE)
write.csv(deg_hfd_4d, "deg_hfd_4d.csv", row.names = TRUE)
write.csv(deg_hfd_2w, "deg_hfd_2w.csv", row.names = TRUE)


# Contar genes regulados al alza y a la baja en Control
n_up_control <- sum(deg_control$avg_log2FC > 0)
n_down_control <- sum(deg_control$avg_log2FC < 0)

# Contar genes regulados al alza y a la baja en 4dHFD
n_up_hfd_4d <- sum(deg_hfd_4d$avg_log2FC > 0)
n_down_hfd_4d <- sum(deg_hfd_4d$avg_log2FC < 0)

# Contar genes regulados al alza y a la baja en 2wHFD
n_up_hfd_2w <- sum(deg_hfd_2w$avg_log2FC > 0)
n_down_hfd_2w <- sum(deg_hfd_2w$avg_log2FC < 0)

# Mostrar los resultados
cat("Control: ", n_up_control, "genes upregulated y", n_down_control, "genes downregulated\n")
cat("4dHFD: ", n_up_hfd_4d, "genes upregulated y", n_down_hfd_4d, "genes downregulated\n")
cat("2wHFD: ", n_up_hfd_2w, "genes upregulated y", n_down_hfd_2w, "genes downregulated\n")


# Encontrar genes compartidos entre las condiciones
shared_genes <- Reduce(intersect, list(rownames(deg_control), rownames(deg_hfd_4d), rownames(deg_hfd_2w)))
write.csv(shared_genes, "shared_genes.csv", row.names = TRUE)


# Instala el paquete si no lo tienes
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}
library(VennDiagram)

# Extrae los nombres de los genes
genes_control <- rownames(deg_control)
genes_hfd_4d <- rownames(deg_hfd_4d)
genes_hfd_2w <- rownames(deg_hfd_2w)

# Crear el diagrama de Venn
venn.plot <- venn.diagram(
  x = list(
    Control = genes_control,
    `4dHFD` = genes_hfd_4d,
    `2wHFD` = genes_hfd_2w
  ),
  filename = NULL, # Para que lo devuelva como objeto
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  main = "Diagrama de Venn: Capilar Tipo II"
)

# Guardar el diagrama como imagen
png("venn_diagram_DEG.png", width = 800, height = 650)
grid.draw(venn.plot)
dev.off()


# Seleccionar los top 5 upregulated y downregulated de control
top5_up_control <- head(deg_control[order(-deg_control$avg_log2FC), ], 5)
top5_down_control <- head(deg_control[order(deg_control$avg_log2FC), ], 5)

# Unir los resultados en un solo dataframe
top_genes_control <- rbind(
  data.frame(Gene = rownames(top5_up_control), avg_log2FC = top5_up_control$avg_log2FC, Direction = "Up"),
  data.frame(Gene = rownames(top5_down_control), avg_log2FC = top5_down_control$avg_log2FC, Direction = "Down")
)

# Repite el proceso para 4dHFD
top5_up_hfd_4d <- head(deg_hfd_4d[order(-deg_hfd_4d$avg_log2FC), ], 5)
top5_down_hfd_4d <- head(deg_hfd_4d[order(deg_hfd_4d$avg_log2FC), ], 5)

top_genes_hfd_4d <- rbind(
  data.frame(Gene = rownames(top5_up_hfd_4d), avg_log2FC = top5_up_hfd_4d$avg_log2FC, Direction = "Up"),
  data.frame(Gene = rownames(top5_down_hfd_4d), avg_log2FC = top5_down_hfd_4d$avg_log2FC, Direction = "Down")
)

# Repite el proceso para 2wHFD
top5_up_hfd_2w <- head(deg_hfd_2w[order(-deg_hfd_2w$avg_log2FC), ], 5)
top5_down_hfd_2w <- head(deg_hfd_2w[order(deg_hfd_2w$avg_log2FC), ], 5)

top_genes_hfd_2w <- rbind(
  data.frame(Gene = rownames(top5_up_hfd_2w), avg_log2FC = top5_up_hfd_2w$avg_log2FC, Direction = "Up"),
  data.frame(Gene = rownames(top5_down_hfd_2w), avg_log2FC = top5_down_hfd_2w$avg_log2FC, Direction = "Down")
)

# Crear barplot para Control
ggplot(top_genes_control, aes(x = reorder(Gene, avg_log2FC), y = avg_log2FC, fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 5 Up y Down NCD: Capilar Tipo II", x = "Genes", y = "Log2FC") +
  theme_minimal()
ggsave(filename=paste(path.guardar,"top_genes_NCD.png",sep="/"),width=10,height=6)

# Crear barplot para 4dHFD
ggplot(top_genes_hfd_4d, aes(x = reorder(Gene, avg_log2FC), y = avg_log2FC, fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 5 Up y Down 4dHFD: Capilar Tipo II", x = "Genes", y = "Log2FC") +
  theme_minimal()
ggsave(filename=paste(path.guardar,"top_genes_4dHFD.png",sep="/"),width=10,height=6)

# Crear barplot para 2wHFD
ggplot(top_genes_hfd_2w, aes(x = reorder(Gene, avg_log2FC), y = avg_log2FC, fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 5 Up y Down 2wHFD: Capilar Tipo II", x = "Genes", y = "Log2FC") +
  theme_minimal()
ggsave(filename=paste(path.guardar,"top_genes_2wHFD.png",sep="/"),width=10,height=6)


# Combinar los DEG de cada condición con los genes compartidos
deg_control_shared <- deg_control[rownames(deg_control) %in% shared_genes, ]
deg_hfd_4d_shared <- deg_hfd_4d[rownames(deg_hfd_4d) %in% shared_genes, ]
deg_hfd_2w_shared <- deg_hfd_2w[rownames(deg_hfd_2w) %in% shared_genes, ]

# Crear un único dataframe que combine los log2FC de todas las condiciones
shared_genes_combined <- data.frame(
  Gene = shared_genes,
  log2FC_Control = deg_control_shared$avg_log2FC,
  log2FC_HFD_4d = deg_hfd_4d_shared$avg_log2FC,
  log2FC_HFD_2w = deg_hfd_2w_shared$avg_log2FC
)

# Calcular el promedio de log2FC entre las tres condiciones
shared_genes_combined$Mean_log2FC <- rowMeans(shared_genes_combined[, 2:4], na.rm = TRUE)

# Ordenar por el promedio de log2FC
shared_genes_combined <- shared_genes_combined %>%
  arrange(desc(Mean_log2FC))

# Seleccionar los top 5 upregulados y downregulados
top5_up <- shared_genes_combined %>% slice_max(Mean_log2FC, n = 5)
top5_down <- shared_genes_combined %>% slice_min(Mean_log2FC, n = 5)

# Combinar ambos conjuntos para el barplot
top_genes <- bind_rows(top5_up, top5_down)

ggplot(top_genes, aes(x = reorder(Gene, Mean_log2FC), y = Mean_log2FC, fill = Mean_log2FC > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 5 Up y Down NCD: Capilar Tipo II", x = "Genes", y = "Log2FC") +
  theme_minimal()
ggsave(filename=paste(path.guardar,"comun_ctrl.png",sep="/"),width=10,height=6)

# Genes únicos en cada condición
genes_unicos_control <- setdiff(genes_control, c(genes_hfd_4d, genes_hfd_2w))
genes_unicos_hfd_4d <- setdiff(genes_hfd_4d, c(genes_control, genes_hfd_2w))
genes_unicos_hfd_2w <- setdiff(genes_hfd_2w, c(genes_control, genes_hfd_4d))

# Mostrar el número de genes comunes y únicos
cat("Genes comunes entre las 3 condiciones:", length(genes_comunes), "\n")
cat("Genes únicos en Control:", length(genes_unicos_control), "\n")
cat("Genes únicos en 4dHFD:", length(genes_unicos_hfd_4d), "\n")
cat("Genes únicos en 2wHFD:", length(genes_unicos_hfd_2w), "\n")

# Guardar listas en archivos
write.csv(genes_unicos_control, "genes_unicos_control.csv", row.names = FALSE)
write.csv(genes_unicos_hfd_4d, "genes_unicos_hfd_4d.csv", row.names = FALSE)
write.csv(genes_unicos_hfd_2w, "genes_unicos_hfd_2w.csv", row.names = FALSE)

# Filtrar genes únicos en cada condición
deg_unicos_control <- deg_control[rownames(deg_control) %in% genes_unicos_control, ]
deg_unicos_hfd_4d <- deg_hfd_4d[rownames(deg_hfd_4d) %in% genes_unicos_hfd_4d, ]
deg_unicos_hfd_2w <- deg_hfd_2w[rownames(deg_hfd_2w) %in% genes_unicos_hfd_2w, ]

# Seleccionar los top 5 up y down para los genes únicos
top5_up_unicos_control <- head(deg_unicos_control[order(-deg_unicos_control$avg_log2FC), ], 5)
top5_down_unicos_control <- head(deg_unicos_control[order(deg_unicos_control$avg_log2FC), ], 5)

top5_up_unicos_hfd_4d <- head(deg_unicos_hfd_4d[order(-deg_unicos_hfd_4d$avg_log2FC), ], 5)
top5_down_unicos_hfd_4d <- head(deg_unicos_hfd_4d[order(deg_unicos_hfd_4d$avg_log2FC), ], 5)

top5_up_unicos_hfd_2w <- head(deg_unicos_hfd_2w[order(-deg_unicos_hfd_2w$avg_log2FC), ], 5)
top5_down_unicos_hfd_2w <- head(deg_unicos_hfd_2w[order(deg_unicos_hfd_2w$avg_log2FC), ], 5)

# Crear un dataframe para los genes comunes
data_comunes <- rbind(
  data.frame(Gene = rownames(top5_up_comunes), avg_log2FC = top5_up_comunes$avg_log2FC, Direction = "Up", Category = "Common"),
  data.frame(Gene = rownames(top5_down_comunes), avg_log2FC = top5_down_comunes$avg_log2FC, Direction = "Down", Category = "Common")
)

# Crear dataframes para genes únicos
data_unicos_control <- rbind(
  data.frame(Gene = rownames(top5_up_unicos_control), avg_log2FC = top5_up_unicos_control$avg_log2FC, Direction = "Up", Category = "Unique to Control"),
  data.frame(Gene = rownames(top5_down_unicos_control), avg_log2FC = top5_down_unicos_control$avg_log2FC, Direction = "Down", Category = "Unique to Control")
)

data_unicos_hfd_4d <- rbind(
  data.frame(Gene = rownames(top5_up_unicos_hfd_4d), avg_log2FC = top5_up_unicos_hfd_4d$avg_log2FC, Direction = "Up", Category = "Unique to 4dHFD"),
  data.frame(Gene = rownames(top5_down_unicos_hfd_4d), avg_log2FC = top5_down_unicos_hfd_4d$avg_log2FC, Direction = "Down", Category = "Unique to 4dHFD")
)

data_unicos_hfd_2w <- rbind(
  data.frame(Gene = rownames(top5_up_unicos_hfd_2w), avg_log2FC = top5_up_unicos_hfd_2w$avg_log2FC, Direction = "Up", Category = "Unique to 2wHFD"),
  data.frame(Gene = rownames(top5_down_unicos_hfd_2w), avg_log2FC = top5_down_unicos_hfd_2w$avg_log2FC, Direction = "Down", Category = "Unique to 2wHFD")
)

# Combinar todos los datos
data_all <- rbind(data_comunes, data_unicos_control, data_unicos_hfd_4d, data_unicos_hfd_2w)

ggplot(data_comunes, aes(x = reorder(Gene, avg_log2FC), y = avg_log2FC, fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 5 Up y Down Comunes: Capilar Tipo II", x = "Genes", y = "Log2FC") +
  theme_minimal()
ggsave(filename=paste(path.guardar,"top_genes_comunes.png",sep="/"),width=10,height=6)

# Barplot
ggplot(data_all, aes(x = reorder(Gene, avg_log2FC), y = avg_log2FC, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Category, scales = "free", ncol = 1) + # Separar por categoría
  coord_flip() +
  labs(title = "Top 5 Up y Down: Genes Comunes y Únicos", x = "Genes", y = "Log2FC") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) # Ajusta el tamaño del texto
ggsave(filename=paste(path.guardar,"top_genes_common_unique.png",sep="/"),width=10,height=30)


## -- GO TERMS ----------------------------------------------------------------

enrich_NCD <- enrichGO(
  gene = genes_control,
  OrgDb = org.Mm.eg.db, 
  keyType = "SYMBOL",
  ont = "BP", # Biological Processes
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

go_NCD <- barplot(enrich_NCD, showCategory = 10, title = "GO Enrichment - NCD")
write.csv(as.data.frame(enrich_NCD), "enrich_NCD.csv")

enrich_4dHFD <- enrichGO(
  gene = genes_hfd_4d,
  OrgDb = org.Mm.eg.db, 
  keyType = "SYMBOL",
  ont = "BP", # Biological Processes
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

go_4dHFD <- barplot(enrich_4dHFD, showCategory = 10, title = "GO Enrichment - 4dHFD")
write.csv(as.data.frame(enrich_4dHFD), "enrich_4dHFD.csv")

# Enriquecimiento para genes únicos en 2wHFD
enrich_2wHFD <- enrichGO(
  gene = genes_hfd_2w,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

go_2wHFD <- barplot(enrich_2wHFD, showCategory = 10, title = "GO Enrichment - 2wHFD")
write.csv(as.data.frame(enrich_2wHFD), "enrich_2wHFD.csv")

enrich_common <- enrichGO(
  gene = genes_comunes,
  OrgDb = org.Mm.eg.db, 
  keyType = "SYMBOL",
  ont = "BP", # Biological Processes
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

go_comon <- barplot(enrich_common, showCategory = 10, title = "GO Enrichment - comunes")
write.csv(as.data.frame(enrich_common), "enrich_common.csv")


combined_plot <- (go_NCD| go_4dHFD) / (go_2wHFD | go_comon)

# Muestra la imagen combinada
print(combined_plot)
ggsave(filename=paste(path.guardar,"GOTerms_CapII_GO.png",sep="/"),plot=combined_plot,width=22,height=15)


# -- Mutante -------------------------------------------------------------------

# Dividir el objeto Seurat por dieta
seurat_list <- SplitObject(data_mut, split.by = "Diet")

# Los resultados serán tres objetos: uno para cada dieta
control_obj <- seurat_list[["NCD"]]
hfd_4d_obj <- seurat_list[["4d_HFD"]]
hfd_2w_obj <- seurat_list[["2w_HFD"]]

# DEG para Control
deg_control <- FindMarkers(control_obj, ident.1 = "3", ident.2 = NULL)

# DEG para 4 días de HFD
deg_hfd_4d <- FindMarkers(hfd_4d_obj, ident.1 = "3", ident.2 = NULL)

# DEG para 2 semanas de HFD
deg_hfd_2w <- FindMarkers(hfd_2w_obj, ident.1 = "3", ident.2 = NULL)
# Guardar las listas de DEG
write.csv(deg_control, "deg_control.csv", row.names = TRUE)
write.csv(deg_hfd_4d, "deg_hfd_4d.csv", row.names = TRUE)
write.csv(deg_hfd_2w, "deg_hfd_2w.csv", row.names = TRUE)


# Contar genes regulados al alza y a la baja en Control
n_up_control <- sum(deg_control$avg_log2FC > 0)
n_down_control <- sum(deg_control$avg_log2FC < 0)

# Contar genes regulados al alza y a la baja en 4dHFD
n_up_hfd_4d <- sum(deg_hfd_4d$avg_log2FC > 0)
n_down_hfd_4d <- sum(deg_hfd_4d$avg_log2FC < 0)

# Contar genes regulados al alza y a la baja en 2wHFD
n_up_hfd_2w <- sum(deg_hfd_2w$avg_log2FC > 0)
n_down_hfd_2w <- sum(deg_hfd_2w$avg_log2FC < 0)

# Mostrar los resultados
cat("Control: ", n_up_control, "genes upregulated y", n_down_control, "genes downregulated\n")
cat("4dHFD: ", n_up_hfd_4d, "genes upregulated y", n_down_hfd_4d, "genes downregulated\n")
cat("2wHFD: ", n_up_hfd_2w, "genes upregulated y", n_down_hfd_2w, "genes downregulated\n")


# Encontrar genes compartidos entre las condiciones
shared_genes <- Reduce(intersect, list(rownames(deg_control), rownames(deg_hfd_4d), rownames(deg_hfd_2w)))
write.csv(shared_genes, "shared_genes.csv", row.names = TRUE)


# Instala el paquete si no lo tienes
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}
library(VennDiagram)

# Extrae los nombres de los genes
genes_control <- rownames(deg_control)
genes_hfd_4d <- rownames(deg_hfd_4d)
genes_hfd_2w <- rownames(deg_hfd_2w)

# Crear el diagrama de Venn
venn.plot <- venn.diagram(
  x = list(
    Control = genes_control,
    `4dHFD` = genes_hfd_4d,
    `2wHFD` = genes_hfd_2w
  ),
  filename = NULL, # Para que lo devuelva como objeto
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  main = "Diagrama de Venn: Capilar Tipo II"
)

# Guardar el diagrama como imagen
png("venn_diagram_DEG.png", width = 800, height = 650)
grid.draw(venn.plot)
dev.off()


# Seleccionar los top 5 upregulated y downregulated de control
top5_up_control <- head(deg_control[order(-deg_control$avg_log2FC), ], 5)
top5_down_control <- head(deg_control[order(deg_control$avg_log2FC), ], 5)

# Unir los resultados en un solo dataframe
top_genes_control <- rbind(
  data.frame(Gene = rownames(top5_up_control), avg_log2FC = top5_up_control$avg_log2FC, Direction = "Up"),
  data.frame(Gene = rownames(top5_down_control), avg_log2FC = top5_down_control$avg_log2FC, Direction = "Down")
)

# Repite el proceso para 4dHFD
top5_up_hfd_4d <- head(deg_hfd_4d[order(-deg_hfd_4d$avg_log2FC), ], 5)
top5_down_hfd_4d <- head(deg_hfd_4d[order(deg_hfd_4d$avg_log2FC), ], 5)

top_genes_hfd_4d <- rbind(
  data.frame(Gene = rownames(top5_up_hfd_4d), avg_log2FC = top5_up_hfd_4d$avg_log2FC, Direction = "Up"),
  data.frame(Gene = rownames(top5_down_hfd_4d), avg_log2FC = top5_down_hfd_4d$avg_log2FC, Direction = "Down")
)

# Repite el proceso para 2wHFD
top5_up_hfd_2w <- head(deg_hfd_2w[order(-deg_hfd_2w$avg_log2FC), ], 5)
top5_down_hfd_2w <- head(deg_hfd_2w[order(deg_hfd_2w$avg_log2FC), ], 5)

top_genes_hfd_2w <- rbind(
  data.frame(Gene = rownames(top5_up_hfd_2w), avg_log2FC = top5_up_hfd_2w$avg_log2FC, Direction = "Up"),
  data.frame(Gene = rownames(top5_down_hfd_2w), avg_log2FC = top5_down_hfd_2w$avg_log2FC, Direction = "Down")
)

# Crear barplot para Control
ggplot(top_genes_control, aes(x = reorder(Gene, avg_log2FC), y = avg_log2FC, fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 5 Up y Down NCD: Capilar Tipo II", x = "Genes", y = "Log2FC") +
  theme_minimal()
ggsave(filename=paste(path.guardar,"top_genes_NCD.png",sep="/"),width=10,height=6)

# Crear barplot para 4dHFD
ggplot(top_genes_hfd_4d, aes(x = reorder(Gene, avg_log2FC), y = avg_log2FC, fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 5 Up y Down 4dHFD: Capilar Tipo II", x = "Genes", y = "Log2FC") +
  theme_minimal()
ggsave(filename=paste(path.guardar,"top_genes_4dHFD.png",sep="/"),width=10,height=6)

# Crear barplot para 2wHFD
ggplot(top_genes_hfd_2w, aes(x = reorder(Gene, avg_log2FC), y = avg_log2FC, fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 5 Up y Down 2wHFD: Capilar Tipo II", x = "Genes", y = "Log2FC") +
  theme_minimal()
ggsave(filename=paste(path.guardar,"top_genes_2wHFD.png",sep="/"),width=10,height=6)


# Combinar los DEG de cada condición con los genes compartidos
deg_control_shared <- deg_control[rownames(deg_control) %in% shared_genes, ]
deg_hfd_4d_shared <- deg_hfd_4d[rownames(deg_hfd_4d) %in% shared_genes, ]
deg_hfd_2w_shared <- deg_hfd_2w[rownames(deg_hfd_2w) %in% shared_genes, ]

# Crear un único dataframe que combine los log2FC de todas las condiciones
shared_genes_combined <- data.frame(
  Gene = shared_genes,
  log2FC_Control = deg_control_shared$avg_log2FC,
  log2FC_HFD_4d = deg_hfd_4d_shared$avg_log2FC,
  log2FC_HFD_2w = deg_hfd_2w_shared$avg_log2FC
)

# Calcular el promedio de log2FC entre las tres condiciones
shared_genes_combined$Mean_log2FC <- rowMeans(shared_genes_combined[, 2:4], na.rm = TRUE)

# Ordenar por el promedio de log2FC
shared_genes_combined <- shared_genes_combined %>%
  arrange(desc(Mean_log2FC))

# Seleccionar los top 5 upregulados y downregulados
top5_up <- shared_genes_combined %>% slice_max(Mean_log2FC, n = 5)
top5_down <- shared_genes_combined %>% slice_min(Mean_log2FC, n = 5)

# Combinar ambos conjuntos para el barplot
top_genes <- bind_rows(top5_up, top5_down)

ggplot(top_genes, aes(x = reorder(Gene, Mean_log2FC), y = Mean_log2FC, fill = Mean_log2FC > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 5 Up y Down NCD: Capilar Tipo II", x = "Genes", y = "Log2FC") +
  theme_minimal()
ggsave(filename=paste(path.guardar,"comun_mut.png",sep="/"),width=10,height=6)

# Genes únicos en cada condición
genes_unicos_control <- setdiff(genes_control, c(genes_hfd_4d, genes_hfd_2w))
genes_unicos_hfd_4d <- setdiff(genes_hfd_4d, c(genes_control, genes_hfd_2w))
genes_unicos_hfd_2w <- setdiff(genes_hfd_2w, c(genes_control, genes_hfd_4d))

# Mostrar el número de genes comunes y únicos
cat("Genes comunes entre las 3 condiciones:", length(genes_comunes), "\n")
cat("Genes únicos en Control:", length(genes_unicos_control), "\n")
cat("Genes únicos en 4dHFD:", length(genes_unicos_hfd_4d), "\n")
cat("Genes únicos en 2wHFD:", length(genes_unicos_hfd_2w), "\n")

# Guardar listas en archivos
write.csv(genes_unicos_control, "genes_unicos_control.csv", row.names = FALSE)
write.csv(genes_unicos_hfd_4d, "genes_unicos_hfd_4d.csv", row.names = FALSE)
write.csv(genes_unicos_hfd_2w, "genes_unicos_hfd_2w.csv", row.names = FALSE)

# Filtrar genes únicos en cada condición
deg_unicos_control <- deg_control[rownames(deg_control) %in% genes_unicos_control, ]
deg_unicos_hfd_4d <- deg_hfd_4d[rownames(deg_hfd_4d) %in% genes_unicos_hfd_4d, ]
deg_unicos_hfd_2w <- deg_hfd_2w[rownames(deg_hfd_2w) %in% genes_unicos_hfd_2w, ]

# Seleccionar los top 5 up y down para los genes únicos
top5_up_unicos_control <- head(deg_unicos_control[order(-deg_unicos_control$avg_log2FC), ], 5)
top5_down_unicos_control <- head(deg_unicos_control[order(deg_unicos_control$avg_log2FC), ], 5)

top5_up_unicos_hfd_4d <- head(deg_unicos_hfd_4d[order(-deg_unicos_hfd_4d$avg_log2FC), ], 5)
top5_down_unicos_hfd_4d <- head(deg_unicos_hfd_4d[order(deg_unicos_hfd_4d$avg_log2FC), ], 5)

top5_up_unicos_hfd_2w <- head(deg_unicos_hfd_2w[order(-deg_unicos_hfd_2w$avg_log2FC), ], 5)
top5_down_unicos_hfd_2w <- head(deg_unicos_hfd_2w[order(deg_unicos_hfd_2w$avg_log2FC), ], 5)

# Crear un dataframe para los genes comunes
data_comunes <- rbind(
  data.frame(Gene = rownames(top5_up_comunes), avg_log2FC = top5_up_comunes$avg_log2FC, Direction = "Up", Category = "Common"),
  data.frame(Gene = rownames(top5_down_comunes), avg_log2FC = top5_down_comunes$avg_log2FC, Direction = "Down", Category = "Common")
)

# Crear dataframes para genes únicos
data_unicos_control <- rbind(
  data.frame(Gene = rownames(top5_up_unicos_control), avg_log2FC = top5_up_unicos_control$avg_log2FC, Direction = "Up", Category = "Unique to Control"),
  data.frame(Gene = rownames(top5_down_unicos_control), avg_log2FC = top5_down_unicos_control$avg_log2FC, Direction = "Down", Category = "Unique to Control")
)

data_unicos_hfd_4d <- rbind(
  data.frame(Gene = rownames(top5_up_unicos_hfd_4d), avg_log2FC = top5_up_unicos_hfd_4d$avg_log2FC, Direction = "Up", Category = "Unique to 4dHFD"),
  data.frame(Gene = rownames(top5_down_unicos_hfd_4d), avg_log2FC = top5_down_unicos_hfd_4d$avg_log2FC, Direction = "Down", Category = "Unique to 4dHFD")
)

data_unicos_hfd_2w <- rbind(
  data.frame(Gene = rownames(top5_up_unicos_hfd_2w), avg_log2FC = top5_up_unicos_hfd_2w$avg_log2FC, Direction = "Up", Category = "Unique to 2wHFD"),
  data.frame(Gene = rownames(top5_down_unicos_hfd_2w), avg_log2FC = top5_down_unicos_hfd_2w$avg_log2FC, Direction = "Down", Category = "Unique to 2wHFD")
)

# Combinar todos los datos
data_all <- rbind(data_comunes, data_unicos_control, data_unicos_hfd_4d, data_unicos_hfd_2w)

ggplot(data_comunes, aes(x = reorder(Gene, avg_log2FC), y = avg_log2FC, fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 5 Up y Down Comunes: Capilar Tipo II", x = "Genes", y = "Log2FC") +
  theme_minimal()
ggsave(filename=paste(path.guardar,"top_genes_comunes.png",sep="/"),width=10,height=6)

# Barplot
ggplot(data_all, aes(x = reorder(Gene, avg_log2FC), y = avg_log2FC, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Category, scales = "free", ncol = 1) + # Separar por categoría
  coord_flip() +
  labs(title = "Top 5 Up y Down: Genes Comunes y Únicos", x = "Genes", y = "Log2FC") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) # Ajusta el tamaño del texto
ggsave(filename=paste(path.guardar,"top_genes_common_unique.png",sep="/"),width=10,height=30)


## -- GO TERMS ----------------------------------------------------------------

enrich_NCD <- enrichGO(
  gene = genes_control,
  OrgDb = org.Mm.eg.db, 
  keyType = "SYMBOL",
  ont = "BP", # Biological Processes
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

go_NCD <- barplot(enrich_NCD, showCategory = 10, title = "GO Enrichment - NCD")
write.csv(as.data.frame(enrich_NCD), "enrich_NCD.csv")

enrich_4dHFD <- enrichGO(
  gene = genes_hfd_4d,
  OrgDb = org.Mm.eg.db, 
  keyType = "SYMBOL",
  ont = "BP", # Biological Processes
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

go_4dHFD <- barplot(enrich_4dHFD, showCategory = 10, title = "GO Enrichment - 4dHFD")
write.csv(as.data.frame(enrich_4dHFD), "enrich_4dHFD.csv")

# Enriquecimiento para genes únicos en 2wHFD
enrich_2wHFD <- enrichGO(
  gene = genes_hfd_2w,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

go_2wHFD <- barplot(enrich_2wHFD, showCategory = 10, title = "GO Enrichment - 2wHFD")
write.csv(as.data.frame(enrich_2wHFD), "enrich_2wHFD.csv")

enrich_common <- enrichGO(
  gene = genes_comunes,
  OrgDb = org.Mm.eg.db, 
  keyType = "SYMBOL",
  ont = "BP", # Biological Processes
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

go_comon <- barplot(enrich_common, showCategory = 10, title = "GO Enrichment - comunes")
write.csv(as.data.frame(enrich_common), "enrich_common.csv")


combined_plot <- (go_NCD| go_4dHFD) / (go_2wHFD | go_comon)

# Muestra la imagen combinada
print(combined_plot)
ggsave(filename=paste(path.guardar,"GOTerms_CapII_GO.png",sep="/"),plot=combined_plot,width=22,height=15)

