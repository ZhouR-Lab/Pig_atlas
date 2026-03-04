##Figure1C
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(ggh4x)
#-----------------------------
#-----------------------------
out_dir <- "figures"               
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stage_col <- "stage"               
reduction_use <- "umap"        
pt_size <- 1   

custom_colors <- c(
  "E38"   = "#430155",
  "E80"   = "#31688E",
  "PN0"   = "#20918D",
  "PN28"  = "#35B879",
  "PN180" = "#FEE724"
)
obj[[stage_col]][, 1] <- factor(obj[[stage_col]][, 1], levels = names(custom_colors))
#=========================================================
# 1) UMAP colored by stage
#=========================================================
p_umap <- DimPlot(
  obj,
  reduction = reduction_use,
  group.by = stage_col,
  pt.size = pt_size
)

axis_trunc <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)

p_umap <- p_umap +
  guides(x = axis_trunc, y = axis_trunc) +
  scale_color_manual(values = custom_colors, drop = FALSE) +
  theme(
    axis.title.x = element_text(size = 18, hjust = 0, face = "bold"),
    axis.title.y = element_text(size = 18, hjust = 0, angle = 90, face = "bold"),
    axis.line = element_line(
      arrow = arrow(length = unit(0.3, "cm"), ends = "last", type = "closed"),
      linewidth = 1.5
    ),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_blank(),
    legend.text = element_text(size = 18)
  )

print(p_umap)

ggsave(
  filename = file.path(out_dir, "Fig1C_stage_umap.pdf"),
  plot = p_umap,
  width = 8, height = 7, dpi = 300
)

#=========================================================
# 2) Cell number per stage (barplot)
#=========================================================
meta_data <- obj@meta.data %>%
  mutate(Cluster = .data[[stage_col]]) %>%
  group_by(Cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Cluster = factor(Cluster, levels = names(custom_colors)))

p_bar <- ggplot(meta_data, aes(x = Cluster, y = n, fill = Cluster)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = custom_colors, drop = FALSE) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(
    breaks = c(20000, 40000, 60000),
    labels = scales::number_format(),
    expand = c(0, 0),
    limits = c(0, 80000)
  ) +
  coord_flip() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.5),
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black", hjust = 1),
    axis.title.x = element_blank(),    
    axis.title.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(y = "Cell Number")

print(p_bar)

ggsave(
  filename = file.path(out_dir, "Fig1C_stage_number.pdf"),
  plot = p_bar,
  width = 6, height = 7, dpi = 300
)


##Figure1D
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(ggh4x)
#-----------------------------
#-----------------------------
out_dir <- "figures"               
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tissue_col <- "tissue"               
reduction_use <- "umap"        
pt_size <- 1   

custom_colors <- c(
  "Heart"="#EB8C86",
  "Kidney"="#EED314",
  "Muscle"="#3498CB",
  "Liver"="#339848",
  "Lung"="#895b8a"
)
obj[[tissue_col]][, 1] <- factor(obj[[tissue_col]][, 1], levels = names(custom_colors))
#=========================================================
# 1) UMAP colored by tissue
#=========================================================
p_umap <- DimPlot(
  obj,
  reduction = reduction_use,
  group.by = tissue_col,
  pt.size = pt_size
)

axis_trunc <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)

p_umap <- p_umap +
  guides(x = axis_trunc, y = axis_trunc) +
  scale_color_manual(values = custom_colors, drop = FALSE) +
  theme(
    axis.title.x = element_text(size = 18, hjust = 0, face = "bold"),
    axis.title.y = element_text(size = 18, hjust = 0, angle = 90, face = "bold"),
    axis.line = element_line(
      arrow = arrow(length = unit(0.3, "cm"), ends = "last", type = "closed"),
      linewidth = 1.5
    ),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_blank(),
    legend.text = element_text(size = 18)
  )

print(p_umap)

ggsave(
  filename = file.path(out_dir, "Fig1D_tissue_umap.pdf"),
  plot = p_umap,
  width = 8, height = 7, dpi = 300
)

#=========================================================
# 2) Cell number per tissue (barplot)
#=========================================================
meta_data <- obj@meta.data %>%
  mutate(Cluster = .data[[tissue_col]]) %>%
  group_by(Cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Cluster = factor(Cluster, levels = names(custom_colors)))

p_bar <- ggplot(meta_data, aes(x = Cluster, y = n, fill = Cluster)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = custom_colors, drop = FALSE) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(
    breaks = c(20000, 50000, 80000),
    labels = scales::number_format(),
    expand = c(0, 0),
    limits = c(0, 100000)
  ) +
  coord_flip() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.5),
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black", hjust = 1),
    axis.title.x = element_blank(),    
    axis.title.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(y = "Cell Number")

print(p_bar)

ggsave(
  filename = file.path(out_dir, "Fig1D_tissue_number.pdf"),
  plot = p_bar,
  width = 6, height = 7, dpi = 300
)

##Figure1G
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

#=========================================================
# Inputs
#   - Use a generic object name for GitHub
#=========================================================
# obj <- readRDS("path/to/your_seurat_object.rds")  # example (keep your own loading code)
stopifnot(exists("obj"))

out_dir <- "figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Column names in obj@meta.data (edit if yours are different)
stage_col   <- "stage"
sample_col  <- "orig.ident"
lineage_col <- "cell_lineage"

meta <- obj@meta.data

#=========================================================
# Stage comparisons (stage1 / stage2)
#=========================================================
comparisons <- list(
  c("E80",   "E38",   "E80_over_E38"),
  c("PN0",   "E80",   "PN0_over_E80"),
  c("PN28",  "PN0",   "PN28_over_PN0"),
  c("PN180", "PN28",  "PN180_over_PN28")
)

comparison_colors <- c(
  "E80_over_E38"     = "#31688e",
  "PN0_over_E80"     = "#20918d",
  "PN28_over_PN0"    = "#35b879",
  "PN180_over_PN28"  = "#fee724"
)

#=========================================================
# 1) Proportion per sample per lineage
#=========================================================
sample_props <- meta %>%
  group_by(
    stage  = .data[[stage_col]],
    sample = .data[[sample_col]],
    lineage = .data[[lineage_col]]
  ) %>%
  summarise(cell_count = n(), .groups = "drop_last") %>%
  mutate(sample_total = sum(cell_count)) %>%
  ungroup() %>%
  mutate(prop = cell_count / sample_total)

#=========================================================
# 2) Average proportion per stage per lineage
#=========================================================
avg_props <- sample_props %>%
  group_by(stage, lineage) %>%
  summarise(mean_prop = mean(prop), .groups = "drop")

#=========================================================
# 3) Compute log2 ratio across stage comparisons
#=========================================================
ratio_results <- lapply(comparisons, function(comp) {
  stage1 <- comp[1]
  stage2 <- comp[2]
  comp_name <- comp[3]
  
  props_stage1 <- avg_props %>%
    filter(stage == stage1) %>%
    select(lineage, mean_prop1 = mean_prop)
  
  props_stage2 <- avg_props %>%
    filter(stage == stage2) %>%
    select(lineage, mean_prop2 = mean_prop)
  
  full_join(props_stage1, props_stage2, by = "lineage") %>%
    mutate(
      mean_prop1 = replace_na(mean_prop1, 0),
      mean_prop2 = replace_na(mean_prop2, 0),
      ratio      = log2((mean_prop1 ) / (mean_prop2 )),
      size_value = (mean_prop1 + mean_prop2) / 2,
      comparison = comp_name
    )
}) %>% bind_rows()

ratio_results$comparison <- factor(
  ratio_results$comparison,
  levels = names(comparison_colors),
  ordered = TRUE
)

# Order lineages by median ratio (descending)
lineage_order <- ratio_results %>%
  group_by(lineage) %>%
  summarise(median_ratio = median(ratio, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(median_ratio)) %>%
  pull(lineage)

ratio_results$lineage <- factor(ratio_results$lineage, levels = lineage_order)

#=========================================================
# 4) Plot
#=========================================================
p <- ggplot(ratio_results, aes(x = lineage, y = ratio)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_point(aes(color = comparison, size = size_value), alpha = 1) +
  scale_color_manual(
    name = "Comparison Group",
    values = comparison_colors,
    breaks = names(comparison_colors)
  ) +
  scale_size_continuous(
    name = "Proportion",
    range = c(1, 10),
    breaks = c(0.05, 0.10, 0.20, 0.40),                
    labels = scales::percent_format(accuracy = 1),
    guide = guide_legend(override.aes = list(color = "black"), nrow = 2)
  ) +
  labs(x = NULL, y = "log2(avg prop ratio)") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.title  = element_text(color = "black", size = 22),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    axis.line  = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.5, "cm"),
    legend.title = element_text(color = "black", size = 18),
    legend.text  = element_text(size = 15)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 6)),
    size  = guide_legend(order = 2)
  )

print(p)

#=========================================================
# 5) Save (no hard-coded local paths)
#=========================================================
out_file <- file.path(out_dir, "Fig1G_lineage_ratio.pdf")

if (capabilities("cairo")) {
  ggsave(filename = out_file, plot = p, width = 7.5, height = 4.8, dpi = 300, device = cairo_pdf)
} else {
  ggsave(filename = out_file, plot = p, width = 7.5, height = 4.8, dpi = 300)
}

###Figure2A
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
  library(stringr)


# -----------------------------
# Output
# -----------------------------
out_dir <- "figures/heatmap"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Metadata columns (edit if your object uses different names)
# -----------------------------
tissue_col   <- "tissue"
stage_col    <- "stage"
lineage_col  <- "cell_lineage"
celltype_col <- "new_celltype"

# -----------------------------
# 1) Build composite cluster id and filter low-count clusters
# -----------------------------
obj$cluster_id <- paste(
  obj[[tissue_col]][, 1],
  obj[[stage_col]][, 1],
  obj[[lineage_col]][, 1],
  obj[[celltype_col]][, 1],
  sep = "#"
)

cluster_counts <- table(obj$cluster_id)
keep_clusters  <- names(cluster_counts)[cluster_counts >= 100]
obj <- subset(obj, subset = cluster_id %in% keep_clusters)

# -----------------------------
# 2) HVGs and average expression per cluster_id
# -----------------------------
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

avg_expr <- AverageExpression(
  obj,
  assays   = "RNA",
  features = VariableFeatures(obj),
  group.by = "cluster_id",
  verbose  = FALSE
)$RNA

# -----------------------------
# 3) Pearson correlation matrix
# -----------------------------
cor_mat <- cor(avg_expr, method = "spearman")

# -----------------------------
# 4) Parse cluster_id back to metadata for annotations
# -----------------------------
info <- as.data.frame(str_split_fixed(colnames(avg_expr), "#", 4))
colnames(info) <- c("Tissue", "Stage", "Cell_lineage", "Celltype")
rownames(info) <- colnames(avg_expr)

# Optional: set desired ordering (edit levels if needed)
info$Stage <- factor(info$Stage, levels = c("E38", "E80", "PN0", "PN28", "PN180"))
info$Tissue <- factor(info$Tissue, levels = c("Muscle", "Heart", "Kidney", "Lung", "Liver"))

# -----------------------------
# 5) Colors (short + reproducible)
# -----------------------------
tissue_cols <- c(
  "Muscle" = "#3498CB",
  "Lung"   = "#895b8a",
  "Heart"  = "#EB8C86",
  "Kidney" = "#EED314",
  "Liver"  = "#339848"
)

stage_cols <- c(
  "E38"   = "#430155",
  "E80"   = "#31688E",
  "PN0"   = "#20918D",
  "PN28"  = "#35B879",
  "PN180" = "#FEE724"
)

set.seed(1)
lineage_lvls <- sort(unique(as.character(info$Cell_lineage)))
lineage_cols <- setNames(rand_color(length(lineage_lvls), luminosity = "bright"), lineage_lvls)

set.seed(2)
celltype_lvls <- sort(unique(as.character(info$Celltype)))
celltype_cols <- setNames(rand_color(length(celltype_lvls), luminosity = "light"), celltype_lvls)

# -----------------------------
# 6) Annotations
# -----------------------------
col_ha <- columnAnnotation(
  Tissue      = info$Tissue,
  Stage       = info$Stage,
  Cell_lineage= info$Cell_lineage,
  Celltype    = info$Celltype,
  col = list(
    Tissue       = tissue_cols,
    Stage        = stage_cols,
    Cell_lineage = lineage_cols,
    Celltype     = celltype_cols
  ),
  show_legend = TRUE
)

row_ha <- rowAnnotation(
  Tissue = info$Tissue,
  Stage  = info$Stage,
  col = list(
    Tissue = tissue_cols,
    Stage  = stage_cols
  ),
  show_legend = TRUE
)

# -----------------------------
# 7) Heatmap
# -----------------------------
hm_cols <- c(
  "#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0",
  "#F7F7F7",
  "#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"
)

p <- Heatmap(
  cor_mat,
  col = hm_cols,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  top_annotation = col_ha,
  left_annotation = row_ha,
  row_title = NULL,
  use_raster = TRUE
)

pdf(file.path(out_dir, "Figure2A_heatmap.pdf"),
    width = 9, height = 7, useDingbats = FALSE)
draw(p)
dev.off()

##Figure2B
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
})

#=========================================================
# Output directory (GitHub-friendly)
#=========================================================
out_dir <- "figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#=========================================================
# 1) Average expression matrix (cluster-level)
#=========================================================
quant_mx <- AverageExpression(
  obj,
  assays   = "RNA",
  group.by = "cluster_id"
)$RNA

data <- na.omit(quant_mx)

# Keep top genes by mean expression (max 2000)
data_mean <- apply(data, 1, mean)
sel_mean  <- order(data_mean, decreasing = TRUE)[seq_len(min(2000, length(data_mean)))]
data      <- data[sel_mean, , drop = FALSE]

# Select top variable genes by CV (sd/mean)
data_rv <- apply(data, 1, function(x) sd(x) / mean(x))
ntop    <- 1000
sel_rv  <- order(data_rv, decreasing = TRUE)[seq_len(min(ntop, length(data_rv)))]
data    <- data[sel_rv, , drop = FALSE]

#=========================================================
# 2) PCA
#=========================================================
pca <- prcomp(t(data), scale. = TRUE)

pca_summary <- summary(pca)
pct_var <- pca_summary$importance[2, ] * 100
PCA_var <- setNames(round(pct_var, 1), names(pct_var))

#=========================================================
# 3) Plotting data
#    Expecting sample names like: Tissue#Stage#Cell_lineage#Celltype
#=========================================================
PCAplotdata <- as.data.frame(pca$x[, 1:4, drop = FALSE])

info <- as.data.frame(str_split_fixed(rownames(PCAplotdata), pattern = "#", n = 4))
colnames(info) <- c("Tissue", "Stage", "Cell_lineage", "Celltype")

info$Stage <- factor(info$Stage, levels = c("E38", "E80", "PN0", "PN28", "PN180"))
info$Tissue <- factor(info$Tissue, levels = c("Muscle", "Heart", "Kidney", "Lung", "Liver"))

PCAplotdata <- cbind(PCAplotdata, info)

#=========================================================
# 4) PCA scatter plot
#=========================================================
x_axis <- "PC1"
y_axis <- "PC2"

tissue_colors <- c(
  "Muscle" = "#3498CB",
  "Lung"   = "#895b8a",
  "Heart"  = "#EB8C86",
  "Kidney" = "#EED314",
  "Liver"  = "#339848"
)

p <- ggplot(PCAplotdata, aes(x = .data[[x_axis]], y = .data[[y_axis]])) +
  geom_point(aes(shape = Stage, color = Tissue), size = 4, alpha = 1) +
  scale_color_manual(values = tissue_colors, drop = FALSE) +
  geom_vline(xintercept = 0, color = "gray50", linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "gray50", linetype = "dashed", linewidth = 0.4) +
  scale_x_continuous(expand = c(0.2, 0.2)) +
  scale_y_continuous(expand = c(0.2, 0.2)) +
  labs(
    x = paste0(x_axis, " (", PCA_var[[x_axis]], "%)"),
    y = paste0(y_axis, " (", PCA_var[[y_axis]], "%)")
  ) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", fill = "transparent"),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 15, color = "black"),
    axis.title.x = element_text(size = 18, vjust = 0),
    axis.title.y = element_text(size = 16, color = "black", vjust = 0, margin = margin(r = 20)),
    axis.text.x  = element_text(size = 15, color = "black", hjust = 1),
    axis.text.y  = element_text(size = 15, color = "black")
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 5)),
    shape = guide_legend(override.aes = list(size = 5))
  )

print(p)

#=========================================================
# 5) Save figure (no local absolute path)
#=========================================================
ggsave(
  filename = file.path(out_dir, "Fig2B_PCA.pdf"),
  plot     = p,
  device   = cairo_pdf,
  width    = 7.5,
  height   = 6,
  dpi      = 300
)

###Figure2C
VarExp <- function(counts, meta, threshold, inter){
  counts.center <- t(apply(counts, 1, scale, center=TRUE, scale=FALSE))
  cor.counts <- cor(counts.center)
  dim(cor.counts)
  eigen.counts <- eigen(cor.counts)
  eigen.mat <- eigen.counts$vectors
  eigen.val <- eigen.counts$values
  n.eigen <- length(eigen.val)
  eigen.val.sum <- sum(eigen.val)
  percents.pcs <- eigen.val/eigen.val.sum
  meta <- as.data.frame(meta)
  
  all <- 0
  npc.in <- 0
  for(i in 1:n.eigen){
    all <- all + percents.pcs[i]
    npc.in <- npc.in + 1
    if(all > threshold){break}
  }
  if (npc.in < 3) {npc.in <- 3}
  
  pred.list <- colnames(meta)
  meta <- droplevels(meta)
  
  n.preds <- ncol(meta) + 1
  if(inter) {n.preds <- n.preds + choose(ncol(meta),2)}
  
  ran.pred.list <- c()
  for(i in 1:ncol(meta)){
    ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i],")"))
  }
  ##interactions
  if(inter){
    for(i in 1:(ncol(meta)-1)){
      for(j in (i+1):ncol(meta)){
        ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i], ":", pred.list[j], ")"))
        pred.list <- c(pred.list, paste0(pred.list[i], ":", pred.list[j]))
      }
    }
  }
  formula <- paste(ran.pred.list, collapse = " + ")
  formula <- paste("pc", formula, sep=" ~ ")
  ran.var.mat <- NULL
  for(i in 1:npc.in){
    dat <- cbind(eigen.mat[,i],meta)
    colnames(dat) <- c("pc",colnames(meta))
    Rm1ML <- lme4::lmer(formula, dat, REML = TRUE, verbose = FALSE, na.action = na.omit,
                        control = lmerControl(check.nobs.vs.nlev = "ignore",
                                              check.nobs.vs.rankZ = "ignore",
                                              check.nobs.vs.nRE="ignore"))
    var.vec <- unlist(VarCorr(Rm1ML))
    ran.var.mat <- rbind(ran.var.mat, c(var.vec[pred.list], resid = sigma(Rm1ML)^2))
  }
  ran.var.mat.std <- ran.var.mat/rowSums(ran.var.mat)
  wgt.vec <- eigen.val/eigen.val.sum
  prop.var <- colSums(ran.var.mat.std*wgt.vec[1:npc.in])
  std.prop.var <- prop.var/sum(prop.var)
  std.prop.var
}

quant_mx <- AverageExpression(
  obj,
  assays   = "RNA",
  group.by = "cluster_id"
)$RNA

data <- na.omit(quant_mx)
data.mean<-apply(data, 1, function(x){mean(x)})
select.mean <- order(data.mean, decreasing=TRUE)[seq_len(min(2000, length(data.mean)))]
data<-data[select.mean,]
data.rv <- apply(data, 1, function(x){sd(x)/mean(x)})
ntop <-1000
info<-as.data.frame(str_split_fixed(colnames(data),pattern="#",n=3))
colnames(info)<-c("Tissue","Stage","Cell_lineage") 
rownames(info)<-colnames(data)
var<-VarExp(data,info,10,TRUE)
df <- data.frame(eff=names(var), prop=var)
p<-ggplot(df, aes(x=eff, y=prop,fill=eff))+
  ggplot2::geom_bar(stat="identity")+
  ggplot2::geom_text(aes(label=round(prop,3), y=prop+0.04), size=4) +
  labs(x= "Effects", y= "Weighted average proportion variance")
print(p)

ggsave(
  filename = file.path(out_dir, "Fig2C_PVCA.pdf"),
  plot = p,
  device = cairo_pdf,
  width = 7, height = 6, dpi = 300
)

###Figure2D
tissue_cell_lineage_deg <- obj
df <- tissue_cell_lineage_deg %>%  
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) %>%  
  mutate(change = ifelse(avg_log2FC > 0, "up", "down"))  

plot_data <- df %>%
  filter(group %in% c("Erythroid","Immune" ,"Neural"  , "Epithelial","Stromal", "Muscle" ,"Endothelial")) %>%
  mutate(group = factor(group, levels = c("Endothelial", "Epithelial","Erythroid","Immune" ,"Muscle" ,"Neural"  ,"Stromal" ))) %>%
  group_by(group, compare) %>%
  summarise(gene_count = n(), .groups = 'drop') %>%
  group_by(group) %>%
  mutate(percent = gene_count / sum(gene_count) * 100) %>%
  ungroup() %>%
  mutate(compare = factor(compare, 
                          levels = c("E80_vs_E38", "PN0_vs_E80", 
                                     "PN28_vs_PN0", "PN180_vs_PN28")))

p <- ggplot(plot_data,
            aes(x = group,
                y = percent,
                stratum = compare,
                alluvium = compare,
                fill = compare)) +
  geom_flow(alpha = 0.7, width = 0.7) +
  geom_stratum(width = 0.7) +
  scale_fill_manual(
    name = "Comparison Group",
    values = c(
      "E80_vs_E38" = "#31688e", 
      "PN0_vs_E80" = "#20918d", 
      "PN28_vs_PN0" = "#35b879", 
      "PN180_vs_PN28" = "#fee724"
    )
  ) +
  scale_y_continuous(name = "Percentage of DEGs (%)", 
                     expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.0, 0.0)) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.y = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1, color = "black"),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.1, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  guides(fill = guide_legend(
    override.aes = list(shape = 21, size = 1),
    position = "right",
    direction = "vertical"
  ))
print(p)
ggsave(filename = file.path(out_dir, "Fig2D_cell_lineage_DEG.pdf"),
       plot = p, device = cairo_pdf, width = 11,  height = 7,  dpi = 300)

###Figure2F
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Order
compare_order <- c("E80_vs_E38", "PN0_vs_E80", "PN28_vs_PN0", "PN180_vs_PN28")
lineage_order <- c("Endothelial", "Epithelial", "Erythroid",
                   "Immune", "Muscle", "Neural", "Stromal")

# 2. Filter one tissue (edit "Heart" to "Muscle", "Kidney", "Lung", "Liver")
tissue_df <- obj %>%
  filter(Tissue == "Heart")

# 3. Count DEGs per comparison and lineage
deg_counts <- tissue_df %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) %>%
  group_by(Comparison, CellLineage) %>%
  summarise(DEG_count = n(), .groups = "drop")

# 4. Full combination grid
all_combinations <- expand.grid(
  CellLineage = factor(lineage_order, levels = lineage_order),
  Comparison = factor(compare_order, levels = compare_order),
  stringsAsFactors = FALSE
)

# 5. Merge and fill missing with 0
deg_counts_full <- all_combinations %>%
  left_join(deg_counts, by = c("Comparison", "CellLineage")) %>%
  mutate(DEG_count = replace_na(DEG_count, 0))

# 6. Colors
lineage_colors <- c(
  "Endothelial" = "#4673aa",
  "Epithelial"  = "#e89143",
  "Erythroid"   = "#E3CC52",
  "Immune"      = "#D48AAF",
  "Muscle"      = "#d3381c",
  "Neural"      = "#727171",
  "Stromal"     = "#93ca76"
)

deg_counts_full <- deg_counts_full %>%
  mutate(DEG_count_log10 = log10(DEG_count + 1))

p <- ggplot(deg_counts_full, aes(y = CellLineage, x = DEG_count_log10)) +
  geom_point(
    aes(fill = CellLineage),
    size = 4, shape = 21, color = "black", stroke = 0.5, show.legend = FALSE
  ) +
  scale_fill_manual(values = lineage_colors) +
  facet_grid(
    factor(Comparison, levels = compare_order) ~ .,
    drop = FALSE
  ) +
  labs(
    x = expression(log[10]~"(DEG count + 1)"),
    y = ""
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.1)),
    breaks = log10(c(1, 10, 100, 1000)),
    labels = c("0", "10", "100", "1000")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "gray80", size = 0.5),
    panel.grid.minor = element_line(color = "gray90", size = 0.3),
    panel.grid.major.y = element_line(color = "gray80", size = 0.5),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(size = 14, color = "black"),
    axis.title.y = element_blank(),
    strip.text.y = element_text(angle = 0, size = 12, color = "black"),
    strip.background = element_rect(fill = "grey80", color = "black"),
    panel.spacing = unit(0.8, "lines"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )
print(p)
ggsave(
  filename = file.path(out_dir, "Fig2F_DEG.pdf"),
  plot = p, device = "pdf", width = 3.8, height = 6.1, dpi = 300
)

##Figure3
library(reshape2)
library(igraph)
library(ggrepel)
library(hdWGCNA)
pdf(file.path(out_dir, "Figure3A_hubgene_umap.pdf"), width = 15, height = 15)
options(future.globals.maxSize = 30 * 1024^3)
ModuleUMAPPlot(
  obj,
  edge.alpha = 0.5,
  sample_edges = TRUE,
  keep_grey_edges = FALSE,
  edge_prop = 0.075,
  #label_genes = label_genes,
  label_hubs = 5
)
dev.off()

pdf(file.path(out_dir, "Figure3B_Dendrogram.pdf"), width = 9, height = 7)
PlotDendrogram(obj, main = "Dendrogram")
dev.off()

pdf(file.path(out_dir, "Figure3C_ModuleCorrelogram.pdf"), width = 10, height = 8)
ModuleCorrelogram(obj)
dev.off()

# Figure3D: plot with Seurat's DotPlot function
p <- DotPlot(obj, features = mods, group.by = "cell_type")
p <- p +
  RotatedAxis() +
  scale_color_gradient2(high = "red", mid = "grey95", low = "blue")
p

# Figure3E
obj <- RunEnrichr(obj, dbs = dbs, species = "Sus scrofa")
# enrichr dotplot
p <- EnrichrDotPlot(
  obj,
  database = "GO_Biological_Process_2023",
  n_terms = 2,
  break_ties = TRUE,
  p_adj = FALSE
) +
  scale_color_stepsn(colors = rev(viridis::magma(256))) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = "white", colour = NA)
  )

pdf(file.path(out_dir, "Figure3E_GO_BP.pdf"), width = 11, height = 12, useDingbats = FALSE)
print(p)
dev.off()


#Figure4a
library(Seurat)
library(clusterProfiler)
library(org.Ss.eg.db)
library(ComplexHeatmap)
library(tidyverse)
data<-readRDS("CSI.rds")
info<-read.csv("Module_TF_list.csv",row.names = 1)
x=pheatmap::pheatmap(data,clustering_method = "ward.D2",cutree_rows = 8,cutree_cols = 8)
data<-data[x$tree_row$order,x$tree_row$order]
info<-info[rownames(data),]
info[["Module"]]<-paste0("Module",info$Cluster)
# CSI heatmap
row_ha= rowAnnotation(Module=info$Module,
                      col=list(Module=c("Module1"="#0073C2",
                                        "Module2"="#EFC000",
                                        "Module3"="#868686",
                                        "Module4"="#CD534C",
                                        "Module5"="#003C67",
                                        "Module6"="#8F7700",
                                        "Module7"="#3B3B3B",
                                        "Module8"="#A73030")),show_legend=FALSE)
col_ha= columnAnnotation(Module=info$Module,
                         col=list(Module=c("Module1"="#0073C2",
                                           "Module2"="#EFC000",
                                           "Module3"="#868686",
                                           "Module4"="#CD534C",
                                           "Module5"="#003C67",
                                           "Module6"="#8F7700",
                                           "Module7"="#3B3B3B",
                                           "Module8"="#A73030")),
                         show_legend=FALSE)
p<-ComplexHeatmap::Heatmap(data,col=wolfgang_extra,
                           cluster_rows =T,cluster_columns=F,
                           show_row_names = F,show_column_names = F,
                           row_names_gp=grid::gpar(fontsize=6),
                           show_column_dend=F,show_row_dend=F,
                           cluster_column_slices=T,cluster_row_slices=F,
                           row_split = factor(info$Module,levels =unique(info$Module)),
                           column_split = factor(info$Module,levels =unique(info$Module)),
                           row_title=NULL,use_raster=F,column_title_rot = 45,
                           column_title_gp=grid::gpar(fontsize=8),
                           #border = TRUE, 
                           left_annotation = row_ha,top_annotation = col_ha,
                           heatmap_legend_param = list(title="connection specificity index",direction = "horizontal",
                                                       title_position = "topcenter",
                                                       legend_width = unit(6, "cm")))
# Module TF GO enrichment Analysis 
gene_list<-split(info$TF,data$Module)
go1 <- enrichGO(gene=gene_list$Module1, OrgDb = org.Ss.eg.db,ont= "ALL",keyType = "SYMBOL",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.05, readable= TRUE) 
go2 <- enrichGO(gene=gene_list$Module2, OrgDb = org.Ss.eg.db,ont= "ALL",keyType = "SYMBOL",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.05, readable= TRUE) 
go3 <- enrichGO(gene=gene_list$Module3, OrgDb = org.Ss.eg.db,ont= "ALL",keyType = "SYMBOL",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.05, readable= TRUE) 
go4 <- enrichGO(gene=gene_list$Module4, OrgDb = org.Ss.eg.db,ont= "ALL",keyType = "SYMBOL",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.05, readable= TRUE) 
go5 <- enrichGO(gene=gene_list$Module5, OrgDb = org.Ss.eg.db,ont= "ALL",keyType = "SYMBOL",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.05, readable= TRUE) 
go6 <- enrichGO(gene=gene_list$Module6, OrgDb = org.Ss.eg.db,ont= "ALL",keyType = "SYMBOL",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.05, readable= TRUE) 
go7 <- enrichGO(gene=gene_list$Module7, OrgDb = org.Ss.eg.db,ont= "ALL",keyType = "SYMBOL",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.05, readable= TRUE) 
go8 <- enrichGO(gene=gene_list$Module8, OrgDb = org.Ss.eg.db,ont= "ALL",keyType = "SYMBOL",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.05, readable= TRUE) 

#Figure4b
solarExtra<-c("#3361A5","#248AF3","#14B3FF","#88CEEF","#C1D5DC","#EAD397","#FDB31A","#E42A2A","#A31D1D")
p1<-FeaturePlot(proj,features = "Module1",order = T,split.by = "group")&scale_colour_gradientn(colours = solarExtra)&theme(legend.position = "right")
p2<-FeaturePlot(proj,features = "Module2",order = T,split.by = "group")&scale_colour_gradientn(colours = solarExtra)&theme(legend.position = "right")
p3<-FeaturePlot(proj,features = "Module3",order = T,split.by = "group")&scale_colour_gradientn(colours = solarExtra)&theme(legend.position = "right")
p4<-FeaturePlot(proj,features = "Module4",order = T,split.by = "group")&scale_colour_gradientn(colours = solarExtra)&theme(legend.position = "right")
p5<-FeaturePlot(proj,features = "Module5",order = T,split.by = "group")&scale_colour_gradientn(colours = solarExtra)&theme(legend.position = "right")
p6<-FeaturePlot(proj,features = "Module6",order = T,split.by = "group")&scale_colour_gradientn(colours = solarExtra)&theme(legend.position = "right")
p7<-FeaturePlot(proj,features = "Module7",order = T,split.by = "group")&scale_colour_gradientn(colours = solarExtra)&theme(legend.position = "right")
p8<-FeaturePlot(proj,features = "Module8",order = T,split.by = "group")&scale_colour_gradientn(colours = solarExtra)&theme(legend.position = "right")

#Figure4c
modulescore<-readRDS("F://pig_module_score.rds")
cell_lineage_score <-modulescore%>%group_by(cell_lineage)%>%
  summarise(M1 = mean(Module1_score),
            M2 = mean(Module2_score),
            M3 = mean(Module3_score),
            M4 = mean(Module4_score),
            M5 = mean(Module5_score),
            M6 = mean(Module6_score),
            M7 = mean(Module7_score),
            M8 = mean(Module8_score))%>%
  as.data.frame()
rownames(cell_lineage_score)<-cell_lineage_score$cell_lineage
cell_lineage_score$cell_lineage<-NULL
ComplexHeatmap::Heatmap(t(scale(cell_lineage_score)),col=solarExtra,
                        cluster_rows =T,cluster_columns=F,
                        show_row_names = F,show_column_names = F,
                        row_names_gp=grid::gpar(fontsize=6),
                        show_column_dend=F,show_row_dend=F,
                        cluster_column_slices=T,cluster_row_slices=F,
                        row_title=NULL,use_raster=F,column_title_rot = 45,
                        column_title_gp=grid::gpar(fontsize=8))

#Figure4d
modulescore<-readRDS("pig_module_score.rds")
stage <- subset(modulescore, cell_lineage == "Epithelial") %>%group_by(stage) %>%
  summarise(across(all_of(c("M7","M1","M3","M6","M2","M8","M4","M5")), mean)) %>%
  as.data.frame()
rownames(stage) <- stage$stage
stage$stage <- NULL
stage <- stage[, c("M7","M1","M3","M6","M2","M8","M4","M5")]
colnames(stage) <- module_names[colnames(stage)]
stage_long <- stage %>%rownames_to_column("stage") %>%pivot_longer(cols = everything()[-1],names_to = "Module",values_to = "MeanScore")
stage_long_filtered <- stage_long %>%filter(Module %in% c("M4", "M5"))
stage_long_filtered$stage <- factor(stage_long_filtered$stage,levels = c("E38", "E80", "PN0", "PN28", "PN180"))
p <- ggplot(stage_long_filtered, aes(x = stage, y = MeanScore, group = Module, color = Module)) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  scale_color_manual(values = c("M4" = "#1F78B4", "M5" = "#E31A1C")) +  # 蓝红配色
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        panel.grid = element_blank()) +
  labs(x = "Stage", y = "Module Average Score")

###Figure5C
plot_cell_trajectory(
  obj,
  color_by = "State",
  cell_size = 1,
  alpha = 1
)

###Figure5F
visCluster(
  obj,
  plot.type = "both",
  ht.col.list = list(
    col_range = c(-2, 0, 2)
  ),
  pseudotime_col = pseudotime_col,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  annoTerm.data = enrich,
  add.box = FALSE,
  add.line = TRUE,
  line.side = "left",
  go.col = go_col,
  add.bar = FALSE,
  panel.arg = c(2, 0.2, 4, "grey90", NA),
  term.text.limit = c(7, 10),
  word_wrap = TRUE,
  add_new_line = TRUE
)

dev.off()
######Figure5G
exp_matrix <- GetAssayData(obj, assay = "RNA", slot = "data")
solarExtra_colors <- c(
  "#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
  "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"
)
color_palette <- colorRampPalette(solarExtra_colors)(100)
create_heatmap <- function(data, show_rownames, show_colnames, show_legend, title = "") {
  global_min <- min(data, na.rm = TRUE)
  global_max <- max(data, na.rm = TRUE)
  angle_col_val <- if (show_colnames) 45 else 0
  pheatmap::pheatmap(
    as.matrix(data),
    color = color_palette,
    scale = "row",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    fontsize = 10,
    fontsize_row = 10,
    fontsize_col = 10,
    border_color = "white",
    border_width = 0.0000001,
    angle_col = angle_col_val,
    cellwidth = 20,
    cellheight = 20,
    legend = show_legend,
    main = title,
    silent = TRUE
  )
}

###Figure6A
library(tidyverse)
source("R/plotting_config.R")

obj <- readRDS("celltype_cell_lineage_tissue_deg.rds")

geneset1 <- read.csv("European_gene.csv")
geneset2 <- read.csv("Asian_gene.csv")
geneset1$group <- "European"
geneset2$group <- "Asian"
geneset <- rbind(geneset1, geneset2)

head(geneset)

odds_p <- data.frame()
odds <- data.frame()

for (i in seq(unique(geneset$group))) {
  for (j in seq(unique(cell_lineage_specific$group))) {
    data1 <- cell_lineage_specific %>% filter(group == unique(cell_lineage_specific$group)[j])
    data2 <- geneset[which(geneset$group %in% unique(geneset$group)[i]), ]
    a <- length(data1$feature)
    b <- length(data2$gene_name)
    c <- length(intersect(data1$feature, data2$gene_name))
    d <- (27466 - b - c + a)
    m <- matrix(c(c, a - c, b - c, d), nrow = 2, ncol = 2)
    lf <- fisher.test(m)
    odds_p[i, j] <- lf$p.value
    odds[i, j] <- lf$estimate
  }
}

colnames(odds) <- unique(cell_lineage_specific$group)
rownames(odds) <- unique(geneset$group)
colnames(odds_p) <- unique(cell_lineage_specific$group)
rownames(odds_p) <- unique(geneset$group)

odds <- t(odds) %>% as.data.frame()
odds_p <- t(odds_p) %>% as.data.frame()

info <- as.data.frame(str_split_fixed(rownames(odds), "#", 2))
odds$stage <- info$V2
odds$cell_lineage <- info$V1
odds_p$stage <- info$V2
odds_p$cell_lineage <- info$V1

European_odds <- acast(odds, cell_lineage ~ stage, value.var = "European")
Asian_odds <- acast(odds, cell_lineage ~ stage, value.var = "Asian")
European_odds_p <- acast(odds_p, cell_lineage ~ stage, value.var = "European")
Asian_odds_p <- acast(odds_p, cell_lineage ~ stage, value.var = "Asian")

Asian_odds[is.na(Asian_odds)] <- 0
Asian_odds <- Asian_odds[, c(1, 2, 3, 5, 4)]
Asian_odds_p[is.na(Asian_odds_p)] <- 1
Asian_odds_p <- Asian_odds_p[, c(1, 2, 3, 5, 4)]

source("R/plotting_config.R")

h1 <- BORHeatmap(
  Asian_odds, clusterCols = F, clusterRows = T, labelCols = T, labelRows = T,
  dataColors = cmaps_BOR$solarExtra, limits = c(0, 4),
  dataColorMidPoint = 1, showColDendrogram = T, showRowDendrogram = T,
  row_names_gp = gpar(fontsize = axis_font_size),
  column_names_gp = gpar(fontsize = axis_font_size),
  column_title_gp = gpar(fontsize = font_size + 2),
  legendTitle = "Odds_ratio", column_title = "Asian",
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.2f", Asian_odds[i, j]), x, y, gp = gpar(fontsize = 13))
  }
) + theme_BOR

h2 <- BORHeatmap(
  European_odds, clusterCols = F, clusterRows = T, labelCols = T, labelRows = T,
  dataColors = cmaps_BOR$solarExtra, limits = c(0, 4),
  dataColorMidPoint = 1, showColDendrogram = T, showRowDendrogram = T,
  row_names_gp = gpar(fontsize = axis_font_size),
  column_names_gp = gpar(fontsize = axis_font_size),
  column_title_gp = gpar(fontsize = font_size + 2),
  legendTitle = "Odds_ratio", column_title = "European",
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.2f", European_odds[i, j]), x, y, gp = gpar(fontsize = 13))
  }
) + theme_BOR

##3Figure6F
obj <- readRDS("Myonuclei.rds")
selected_groups <- c("E38", "E80", "PN0", "PN28", "PN180")

muscle_subset_selected <- subset(obj, subset = stage %in% selected_groups)

plot_data <- data.frame(
  Expression = muscle_subset_selected[["RNA"]]@data["MYOT", ],
  Period = muscle_subset_selected$stage,
  stringsAsFactors = FALSE
)
plot_data$Period <- factor(plot_data$Period, levels = selected_groups)

min_non_zero <- min(plot_data$Expression[plot_data$Expression > 0])
max_expr <- max(plot_data$Expression)

library(ggplot2)
library(ggridges)

custom_colors <- c(
  "E38"   = "#430155",
  "E80"   = "#31688E",
  "PN0"   = "#20918D",
  "PN28"  = "#35B879",
  "PN180" = "#FEE724"
)

ridge_plot <- ggplot(plot_data, aes(x = Expression, y = Period, fill = Period)) +
  geom_density_ridges(
    scale = 0.9,
    alpha = 0.9,
    rel_min_height = 0.01
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(
    limits = c(min_non_zero, max_expr),
    expand = c(0, 0)
  ) +
  theme(
    panel.background = element_rect(fill = "grey90"),
    panel.grid.major.y = element_line(color = "black", linewidth = 0.2),
    panel.grid.major.x = element_line(color = "white", linewidth = 0.2),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Expression Level",
    y = NULL,
    title = NULL
  ) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, color = "black", hjust = 0.5),
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )
