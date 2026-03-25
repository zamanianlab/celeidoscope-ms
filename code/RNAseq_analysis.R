#############################################
# Transcriptomic analysis of CELeidoscope fluorescent populations
# Description: TPM normalization, clustering, and visualization of fluorescent reporter expression for Figure 4
#############################################

# Libraries
library(tidyverse)
library(cowplot)
library(viridis)
library(ggdendro)
library(Seurat)
library(Matrix)
library(here)
library(stringr)

# Paths
paths <- list(counts_raw = here("data", "Cel_counts.raw.rds"),
  counts_tpm = here("data", "Cel_counts.tpm.rds"),
  metadata = here("data", "sample_metadata.csv"),
  cengen = here("data", "cengen", "032224_L4_all_cells_Seurat5.rds"),
  proteins = here("data", "protein_counts.tsv"),
  results = here("results"))

dir.create(paths$results, showWarnings = FALSE)

# Parameters
min_tpm <- 100
n_clusters <- 6
score_cutoff <- 0.2

sample_order <- c("N2-Cntr", "ZAM47-Neg", "ZAM47-GFP", "ZAM47-YFP", "ZAM47-mCh", "ZAM47-mKO")

# Functions
make_tpm_matrix <- function(tpm_df) {tpm_df %>%
    select(gene_id, sample, expression) %>%
    distinct() %>%
    pivot_wider(names_from = sample, values_from = expression) %>%
    column_to_rownames("gene_id")}

make_zscore_matrix <- function(df, min_tpm = 100) {m <- as.matrix(df)
  # remove zero variance genes
  m <- m[apply(m, 1, var) != 0, ]
  # filter low expression genes
  m <- m[rowSums(m >= min_tpm) >= 1, ]
  # z-score normalization
  t(scale(t(m)))}

# Load data
gene_count <- readRDS(paths$counts_raw)
gene_tpm   <- readRDS(paths$counts_tpm)
samples    <- read_csv(paths$metadata)

#######
# TPM - SAMPLE CLUSTERING AND GENE CLUSTERING
#######
tpm_matrix <- make_tpm_matrix(gene_tpm)
counts_z   <- make_zscore_matrix(tpm_matrix, min_tpm)

# Sample clustering
sample_dist <- dist(t(counts_z))
sample_clust <- hclust(sample_dist, method = "ward.D2")

# Gene clustering
gene_dist <- dist(counts_z)
gene_clust <- hclust(gene_dist, method = "ward.D2")
gene_order <- gene_clust$order

# Heatmap
samples <- samples %>%
  mutate(sample_id = str_remove(sample_id, "-\\d{4}.*$"))

counts_long <- counts_z %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  mutate(gene_id = factor(gene_id, levels = gene_id[gene_order])) %>%
  pivot_longer(-gene_id, names_to = "sample_id", values_to = "expression") %>%
  mutate(sample_id = str_remove(sample_id, "-\\d{4}.*$")) %>%
  left_join(samples, by = "sample_id") %>%
  mutate(sample_id = factor(sample_id, levels = sample_order))

heatmap_plot <- ggplot(counts_long, aes(gene_id, sample_id)) +
  geom_tile(aes(fill = expression)) +
  scale_fill_viridis(option = "magma", limits = c(-2.1, 3.65)) +
  guides(fill=guide_colourbar(title="z-score"), color = "none") +
  labs(x = "C. elegans Genes", y = "") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
heatmap_plot

# Dendrogram
dend_plot <- ggdendrogram(gene_clust, rotate = FALSE, size = 5) + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) 
dend_plot

# Cluster assignment
dend_cluster <- tibble(gene_id = names(cutree(gene_clust, k = n_clusters)),
  dend_cluster = cutree(gene_clust, k = n_clusters))

counts_long <- left_join(counts_long, dend_cluster, by = "gene_id")

# Cluster plot
cluster_plot <- ggplot(counts_long, aes(sample_id, expression)) + 
  geom_line(aes(group = gene_id), color = "grey70", alpha = 0.1) + 
  stat_summary(aes(group = dend_cluster), fun = mean, geom = "line", linewidth = 1) + 
  facet_wrap(~ dend_cluster, ncol = 6) + 
  ylab("z-score") + xlab(" ") +
  theme_classic() +        
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_blank(),
        strip.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45),
        legend.position = "none")
cluster_plot

#######
# Single-Cell UMAP
#######

# RNA data
cengen <- readRDS(paths$cengen)
DefaultAssay(cengen) <- "RNA"

# Expression matrix (log-normalized counts)
expr_mat <- GetAssayData(cengen, layer = "data")

# UMAP embeddings
umap <- as_tibble(cengen@reductions$umap@cell.embeddings, rownames = "cell")

# Cell metadata
meta <- as_tibble(cengen@meta.data, rownames = "cell")

# Gene sets by fp
fp_clusters <- list(mKO2 = 6, mCherry = 1, YFP = 4, GFP = 5)

fp_genes <- map(fp_clusters, ~{counts_long %>%
    filter(dend_cluster == .x) %>%
    pull(gene_id)})

# Module scores
fp_genes_use <- map(fp_genes, ~intersect(.x, rownames(expr_mat)))

# Mean expression per gene set
fp_scores <- map(fp_genes_use, ~{Matrix::colMeans(expr_mat[.x, , drop = FALSE])})
for (name in names(fp_scores)) {cengen[[paste0(name, "_Score")]] <- fp_scores[[name]]}

# Compile in data frame
cengen_df <- meta %>%
  left_join(umap, by = "cell") %>%
  mutate(mKO2_Score = cengen$mKO2_Score, mCherry_Score = cengen$mCherry_Score,
         YFP_Score = cengen$YFP_Score, GFP_Score = cengen$GFP_Score)

cengen_df <- cengen_df %>%
  mutate(max_score = pmax(mKO2_Score, mCherry_Score, YFP_Score, GFP_Score),
    FP = case_when(mKO2_Score == max_score & mKO2_Score > score_cutoff ~ "mKO2",
      mCherry_Score == max_score & mCherry_Score > score_cutoff ~ "mCherry",
      YFP_Score == max_score & YFP_Score > score_cutoff ~ "YFP",
      GFP_Score == max_score & GFP_Score > score_cutoff ~ "GFP", TRUE ~ "Unassigned"))

# UMAP FP plot
UMAP_fp_plot <- ggplot(cengen_df, aes(UMAP_1.y, UMAP_2.y)) +
  geom_point(data = filter(cengen_df, FP == "Unassigned"), aes(color = FP), size = 0.2, alpha = 0.6) +
  geom_point(data = filter(cengen_df, FP != "Unassigned"), aes(color = FP), size = 0.6, alpha = 0.6) +
  scale_color_manual(values = c("mKO2" = "orange", "mCherry" = "red",
               "YFP" = "yellow", "GFP" = "green", "Unassigned" = "gray"),
               breaks = c("mKO2", "mCherry", "YFP", "GFP", "Unassigned")) +
  labs(x = "UMAP_1", y = "UMAP_2", color = "Fluorescent Protein") +
  theme_classic() +
  theme(legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.box.margin = margin(0,0,0,0)) +
  guides(color = guide_legend(nrow = 1,override.aes = list(size = 2)))
UMAP_fp_plot

# Cell type annotations
muscle <- c("Body_wall_muscle", "Body_wall_muscle_anterior", "Vulval_muscle", "Anal_muscle",
            "Coelomocyte", "Gonadal_sheath_cell", "Spermatheca", "hmc",
            "Spermathecal-uterine_junction_or_uterine_toroid", "Distal_tip_cell", "Uterine_cell")

non_neuronal_types <- c(muscle, "Intestine", "Pharyngeal_muscle", "Pharyngeal_gland_cell", "Epidermis", 
                        "Germline", "Sperm", "Seam_cell", "Excretory_cell", "Excretory_gland_cell", "Vulval_cells",
                        "Rectal_gland", "Arcade_cell", "Marginal_cell",  "Unknown_non_neuronal", 
                        "Unknown_non_neuronal_2", "Unknown_non-neuronal_3")

cengen_df$Major_Class <- ifelse(cengen_df$Cell.type %in% non_neuronal_types, cengen_df$Cell.type, "Neuron")

cengen_df$Cell_Type <- case_when(cengen_df$Cell.type == "Intestine" ~ "Intestine",
  cengen_df$Cell.type == "Pharyngeal_muscle" ~ "Pharyngeal_muscle",
  cengen_df$Cell.type %in% muscle ~ "Muscle/Mesoderm",
  !(cengen_df$Cell.type %in% non_neuronal_types) ~ "Neuron", TRUE ~ "Other")

# UMAP cell type plot
UMAP_ct_plot <- ggplot(cengen_df, aes(UMAP_1.y, UMAP_2.y)) +
  geom_point(data = filter(cengen_df, Cell_Type == "Other"), aes(color = Cell_Type), size = 0.2, alpha = 0.6) +
  geom_point(data = filter(cengen_df, Cell_Type != "Other"), aes(color = Cell_Type), size = 0.6, alpha = 0.6) +
  scale_color_manual(values = c("Neuron" = "orange", "Muscle/Mesoderm" = "red",
               "Intestine" = "yellow", "Pharyngeal_muscle" = "green", "Other" = "gray"),
    breaks = c("Neuron", "Muscle/Mesoderm", "Intestine", "Pharyngeal_muscle", "Other")) +
  labs(x = "UMAP_1", y = "UMAP_2", color = "Cell Type") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.x = unit(0.1, "cm"),
        legend.box.margin = margin(0,0,0,0)) +
  guides(color = guide_legend(nrow = 1,override.aes = list(size = 2)))
UMAP_ct_plot


#######
# PROMOTER-ASSOCIATED GENE EXPRESSION
#######

# Genes associated with promoters of interest
promoter_genes <- c("WBGene00003515",  # myo-3
                    "WBGene00003514",  # myo-2
                    "WBGene00009100",  # rgef-1
                    "WBGene00006915")   # vha-6

# Gene expression by sample
tpm_promoter <- gene_tpm %>% 
  select(gene_id, sample, expression) %>% 
  distinct() %>% 
  mutate(sample = stringr::str_remove(sample, "-\\d{4}.*$")) %>% 
  filter(gene_id %in% promoter_genes) %>% 
  mutate(sample = factor(sample, levels = sample_order), 
          gene_id = factor(gene_id, levels = c("WBGene00009100", "WBGene00003515", "WBGene00006915", "WBGene00003514"), 
                            labels = c("rgef-1", "myo-3", "vha-6", "myo-2")))

# Gene expression plot
gene_plot <- ggplot(tpm_promoter, aes(x = sample, y = expression, fill = gene_id)) +
  geom_col(width = 0.5) +
  scale_fill_manual(values = c("myo-2" = "green", "myo-3" = "red", "vha-6" = "yellow", "rgef-1" = "orange")) +
  facet_grid(gene_id ~ ., scales = "free") +
  theme_classic() +
  ylab("TPM") +
  xlab("") +
  theme(panel.border=element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        legend.position = "none", plot.margin = margin(5.5, 20, 5.5, 5.5),
        strip.text.y = element_text(face = "italic")) 
gene_plot

# Protein alignment
proteins <- read_tsv(paths$proteins,col_names = c("sample", "reference", "count"))

# CPM normalization
proteins_cpm <- proteins %>%
  mutate(sample = str_remove(sample, "_proteins")) %>%
  group_by(sample) %>%
  mutate(CPM = count / sum(count) * 1e6) %>%
  ungroup() %>%
  mutate(reference = factor(reference, levels = c("mKO2", "mCherry", "YFP", "GFP")),
         sample = factor(sample, levels = sample_order))

# FP alignment plot
fp_alignment_plot <- ggplot(proteins_cpm, aes(x = sample, y = CPM, fill = reference)) +
  geom_col(width = 0.5) +
  facet_grid(reference ~ ., scales = "free") +   
  scale_fill_manual(values = c("mCherry" = "red", "YFP" = "yellow", "mKO2" = "orange", "GFP" = "green")) +
  theme_classic() +
  ylab("CPM") +
  xlab("") +
  theme(panel.border=element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = margin(5.5, 5.5, 5.5, 20))
fp_alignment_plot

