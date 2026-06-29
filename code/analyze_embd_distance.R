#close inspection time :) 

library(tidyverse)
library(arrow)
library(stringr)

qc_table    <- read_tsv("results/rbcs_qc/rbcs_qc_table.tsv")
esmc_embeddings  <- read_parquet("results/rbcs_embeddings.parquet")
plant_cad_embeddings <- read_tsv("/workdir/hdd29/paralogs/data/plantcad_embeddings_rbcs_qc_table_henry.tsv")
pheno <- read.csv("data/hifi_species_pheno.csv")
metadata <- read.csv("data/hifi_angiosperm_metadata/hifi_representative_newest_per_species.csv")

species_data <- left_join(pheno, metadata, by = c("species" = "Organism.Name"))
colnames(qc_table)[[1]] <- "geneAccession"
qc_table$Accession <- str_split_i(qc_table$geneAccession, "\\_sp_", 1)

colnames(esmc_embeddings)[[1]] <- "geneAccession"
esmc_embeddings$Accession <- str_split_i(esmc_embeddings$geneAccession, "\\_sp_", 1)

# Rename dim_0 … dim_959 to e1 … e960 — shorter and easier to index numerically
esmc_embeddings  <- esmc_embeddings |> rename_with(~ str_replace(., "dim_(\\d+)", "esmc_emb\\1"))

colnames(plant_cad_embeddings)[[2]] <- "geneAccession"
plant_cad_embeddings$Accession <- str_split_i(plant_cad_embeddings$geneAccession, "\\_sp_", 1)

df <- left_join(qc_table, plant_cad_embeddings, by = c("gene_id" = "pid"))
df <- left_join(df, esmc_embeddings, by = c("gene_id" = "gene_id"))


emb_cols <- grep("emb", colnames(df))

pca <- prcomp(df[,emb_cols], scale.=T,center = T )
summary(pca)
pvar <- pca$sdev^2 / sum(pca$sdev^2)
barplot(pvar[1:20])

pca_scores <- pca$x[,1:20]
colnames(pca_scores)
df <- cbind(df, pca_scores)

plot(df$PC1, df$PC2)
head(df$Accession)
head(species_data[,c("Assembly.Accession", "bio8_median")])

df <- df |>
  left_join(
    species_data |> select(Assembly.Accession, bio8_median, species),
    by = c("Accession" = "Assembly.Accession")
  )

ggplot(df, aes(x = PC1, y = PC2, color = bio8_median)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_viridis_c(name = "BIO8 median") +
  theme_classic()

ggplot(df, aes(x = PC1, y = PC2, color = pid_vs_ref.x)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_viridis_c(name = "pid_vs_ref.x") +
  theme_classic()

ggplot(df, aes(x = PC1, y = PC2, color = species)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_viridis_d(name = "species") +
  theme_classic() +
  theme(legend.position = "none")

df <- na.omit(df)
for (taxon in unique(df$species)) {
  
  p <- ggplot(df, aes(x = PC1, y = PC2)) +
    geom_point(color = "grey85", alpha = 0.5, size = 2) +
    geom_point(
      data = df |> filter(species == taxon),
      aes(color = bio8_median.x),
      alpha = 0.9,
      size = 2
    ) +
    scale_color_viridis_c(name = "BIO8 median") +
    labs(title = taxon) +
    theme_classic()
  
  print(p)
}




plot(df$PC1, df$bio8_median.x)
points(df$PC1[grep("Avena", df$species)], df$bio8_median.x[grep("Avena", df$species)], col="seagreen")
points(df$PC1[grep("Saccharum", df$species)], df$bio8_median.x[grep("Saccharum", df$species)], col="gold")
points(df$PC1[grep("sibiricus", df$species)], df$bio8_median.x[grep("sibiricus", df$species)], col="blue")
points(df$PC1[grep("Lolium", df$species)], df$bio8_median.x[grep("Lolium", df$species)], col="skyblue")
points(df$PC1[grep("Bidens", df$species)], df$bio8_median.x[grep("Bidens", df$species)], col="green")
rho <- cor(df$PC1, df$bio8_median.x)
text(20,25, paste0("rho= ", round(rho,2)), col="red")

abline(18,0.9)

copies_below <- df |>
  filter(bio8_median.x < 18 + 0.9 * PC1) |>
  distinct(gene_id) |>
  arrange(gene_id)

copies_below[1:10,]

writeLines(copies_below$gene_id, "results/hifi_rbcs_copies_cold_cluster.txt")

library(tidyverse)

plot_df <- df |>
  mutate(
    highlight_species = case_when(
      str_detect(species, "Avena") ~ "Avena",
      str_detect(species, "Saccharum") ~ "Saccharum",
      str_detect(species, "sibiricus") ~ "Elymus sibiricus",
      str_detect(species, "Lolium") ~ "Lolium multiflorum",
      str_detect(species, "Bidens") ~ "Bidens hawaiensis",
      TRUE ~ "Other"
    )
  )

rho <- cor(plot_df$PC1, plot_df$bio8_median.x, use = "complete.obs")

ggplot(plot_df, aes(x = PC1, y = bio8_median.x)) +
  geom_point(
    data = plot_df |> filter(highlight_species == "Other"),
    color = "black",
    alpha = 0.80,
    size = 1.6
  ) +
  geom_point(
    data = plot_df |> filter(highlight_species != "Other"),
    aes(color = highlight_species),
    alpha = 0.9,
    size = 2.4
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "black",
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0("rho = ", round(rho, 2)),
    hjust = 1.1,
    vjust = 1.5,
    color = "red",
    size = 5
  ) +
  scale_color_manual(
    name = "Highlighted taxa",
    values = c(
      "Avena" = "seagreen",
      "Saccharum" = "gold",
      "Elymus sibiricus" = "blue",
      "Lolium multiflorum" = "skyblue",
      "Bidens hawaiensis" = "green"
    )
  ) +
  labs(
    x = "PC1",
    y = "BIO8 median",
    title = "PC1 vs BIO8 median",
    subtitle = "Highlighted grass and focal taxa"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

#right group is avenas
#lolium multiflorum   elymus sibiricus  -- perenial grasses 
#avenas, elymus sibiricus 
#also bidens hawaiensis 

#right group is avenas
#lolium multiflorum   elymus sibiricus  -- perenial grasses 
#avenas, elymus sibiricus 
#also bidens hawaiensis 

#taxa that split the group - scilia litarderei, 
#smallanthus sonchifolius 


# Drop metadata columns from embeddings before joining (already in qc_table)
embedding_dims_only <- embeddings |> select(gene_id, starts_with("e"))

rbcs <- qc_table |>
  left_join(embedding_dims_only, by = "gene_id")

prot_emb_cols <- grep("esmc_e", colnames(rbcs))

# find rows where all embedding dims are finite
complete_rows <- complete.cases(rbcs[, prot_emb_cols])
cat("rows with NA embeddings:", sum(!complete_rows), "\n")

prot_emb_pca <- prcomp(rbcs[complete_rows, prot_emb_cols], scale. = TRUE, center = TRUE)
col_vars <- apply(rbcs[, prot_emb_cols], 2, var)
cat("zero-variance columns:", sum(col_vars == 0), "\n")
cat("near-zero variance columns:", sum(col_vars < 1e-10), "\n")


# attach scores back, filling NA for the dropped rows so plots still work
pca_scores <- rbcs |>
  mutate(PC1 = NA_real_, PC2 = NA_real_)
pca_scores[complete_rows, c("PC1", "PC2")] <- prot_emb_pca$x[, 1:2]

hist(rbcs$length)
prot_emb_cols <- grep("esmc_e", colnames(rbcs))
library(irlba)
prot_emb_pca <- prcomp_irlba(as.matrix(rbcs[, prot_emb_cols]), n = 10,
                             scale. = TRUE, center = TRUE)
pvar <- prot_emb_pca$sdev^2 / sum(prot_emb_pca$sdev^2)

plot(pvar)
plot(prot_emb_pca$x[,1],prot_emb_pca$x[,2])

library(ggplot2)

# Attach PC scores and metadata to one dataframe for plotting
pca_scores <- as_tibble(prot_emb_pca$x[, 1:2]) |>
  bind_cols(rbcs |> select(species, length, coverage_ref, aln_filter_passed, fold_filter_passed))

# Percent variance explained labels for axes
pc1_label <- sprintf("PC1 (%.1f%%)", pvar[1] * 100)
pc2_label <- sprintf("PC2 (%.1f%%)", pvar[2] * 100)

# Helper so all four plots share the same base aesthetics
base_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  labs(x = pc1_label, y = pc2_label) +
  theme_bw()

# --- by species ---
base_plot + geom_point(aes(colour = species), alpha = 0.7) +
  labs(title = "rbcS ESM-C embeddings — by species", colour = "Species")

# --- by length ---
base_plot + geom_point(aes(colour = length), alpha = 0.7) +
  scale_colour_viridis_c() +
  labs(title = "rbcS ESM-C embeddings — by length", colour = "Length (aa)")

# --- by coverage vs reference ---
base_plot + geom_point(aes(colour = coverage_ref), alpha = 0.7) +
  scale_colour_viridis_c() +
  labs(title = "rbcS ESM-C embeddings — by coverage vs reference", colour = "Coverage (ref)")

# --- by QC filter status (four combinations) ---
pca_scores <- pca_scores |>
  mutate(qc_status = case_when(
    aln_filter_passed  & fold_filter_passed  ~ "pass both",
    aln_filter_passed  & !fold_filter_passed ~ "pass aln only",
    !aln_filter_passed & fold_filter_passed  ~ "pass fold only",
    TRUE                                     ~ "fail both"
  ))

base_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  labs(x = pc1_label, y = pc2_label) +
  theme_bw()

base_plot + geom_point(aes(colour = qc_status), alpha = 0.7) +
  scale_colour_manual(values = c(
    "pass both"      = "#2166ac",
    "pass aln only"  = "#f4a582",
    "pass fold only" = "#92c5de",
    "fail both"      = "#d6604d"
  )) +
  labs(title = "rbcS ESM-C embeddings — by QC filter status", colour = "QC status")

# now it is time to do the plantcad pcas 

library(tidyverse)
library(arrow)

# --- load plantcad embeddings ---
plantcad_raw <- read_tsv("/workdir/hdd29/paralogs/data/plantcad_embeddings_rbcs_qc_table_henry.tsv")

# rename plantcad_emb_N -> pc_e{N} to distinguish from ESM dims
plantcad_dims <- plantcad_raw |>
  rename(gene_id = pid) |>
  rename_with(~ str_replace(., "plantcad_emb_(\\d+)", "pc_e\\1")) |>
  select(gene_id, starts_with("pc_e"))

# join into main table (rbcs already has ESM embeddings + QC metadata)
rbcs_full <- rbcs |> left_join(plantcad_dims, by = "gene_id")

# QC status for colouring (defined earlier, but recompute to be safe)
rbcs_full <- rbcs_full |>
  mutate(qc_status = case_when(
    aln_filter_passed  & fold_filter_passed  ~ "pass both",
    aln_filter_passed  & !fold_filter_passed ~ "pass aln only",
    !aln_filter_passed & fold_filter_passed  ~ "pass fold only",
    TRUE                                     ~ "fail both"
  ))

qc_colours <- c(
  "pass both"      = "#2166ac",
  "pass aln only"  = "#f4a582",
  "pass fold only" = "#92c5de",
  "fail both"      = "#d6604d"
)

# helper: build a pca_scores tibble with metadata attached
make_pca_scores <- function(pca_obj, metadata_df) {
  as_tibble(pca_obj$x[, 1:2]) |>
    bind_cols(metadata_df |> select(species, length, coverage_ref,
                                    aln_filter_passed, fold_filter_passed,
                                    qc_status))
}

pct_var <- function(pca_obj) pca_obj$sdev^2 / sum(pca_obj$sdev^2)

axis_labels <- function(pca_obj) {
  pv <- pct_var(pca_obj)
  list(
    x = sprintf("PC1 (%.1f%%)", pv[1] * 100),
    y = sprintf("PC2 (%.1f%%)", pv[2] * 100)
  )
}

make_plots <- function(pca_scores, ax, title_prefix) {
  base <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
    labs(x = ax$x, y = ax$y) + theme_bw()
  
  p1 <- base + geom_point(aes(colour = species), alpha = 0.7) +
    labs(title = paste(title_prefix, "— by species"), colour = "Species")
  
  p2 <- base + geom_point(aes(colour = length), alpha = 0.7) +
    scale_colour_viridis_c() +
    labs(title = paste(title_prefix, "— by length"), colour = "Length (aa)")
  
  p3 <- base + geom_point(aes(colour = coverage_ref), alpha = 0.7) +
    scale_colour_viridis_c() +
    labs(title = paste(title_prefix, "— by coverage vs ref"), colour = "Coverage (ref)")
  
  p4 <- base + geom_point(aes(colour = qc_status), alpha = 0.7) +
    scale_colour_manual(values = qc_colours) +
    labs(title = paste(title_prefix, "— by QC status"), colour = "QC status")
  
  list(species = p1, length = p2, coverage = p3, qc = p4)
}

# --- PlantCAD PCA ---
plantcad_emb_cols <- grep("^pc_e", colnames(rbcs_full))
plantcad_pca      <- prcomp(rbcs_full[, plantcad_emb_cols], scale. = TRUE, center = TRUE)
plantcad_scores   <- make_pca_scores(plantcad_pca, rbcs_full)
plantcad_ax       <- axis_labels(plantcad_pca)
plantcad_plots    <- make_plots(plantcad_scores, plantcad_ax, "PlantCAD")

# --- Joint PCA (ESM + PlantCAD concatenated) ---
esm_emb_cols  <- grep("^esmc_e", colnames(rbcs_full))
joint_matrix  <- rbcs_full[, c(esm_emb_cols, plantcad_emb_cols)]
joint_pca     <- prcomp(joint_matrix, scale. = TRUE, center = TRUE)
joint_scores  <- make_pca_scores(joint_pca, rbcs_full)
joint_ax      <- axis_labels(joint_pca)
joint_plots   <- make_plots(joint_scores, joint_ax, "Joint ESM+PlantCAD")

# --- print all plots ---
# PlantCAD
pdf("plot.pdf", height=20, width=40)
print(plantcad_plots$species)
dev.off()
print(plantcad_plots$length)
print(plantcad_plots$coverage)
print(plantcad_plots$qc)

# Joint
pdf("plot.pdf", height=20, width=40)
print(joint_plots$species)
dev.off()
print(joint_plots$length)
print(joint_plots$coverage)
print(joint_plots$qc)

library(patchwork)

# rebuild ESM scores in case they weren't saved as a tibble
esm_scores <- make_pca_scores(prot_emb_pca, rbcs_full)
esm_ax     <- axis_labels(prot_emb_pca)

p_esm <- ggplot(esm_scores, aes(PC1, PC2, colour = qc_status)) +
  geom_point(alpha = 0.7) +
  scale_colour_manual(values = qc_colours) +
  labs(title = "ESM-C", x = esm_ax$x, y = esm_ax$y, colour = "QC status") +
  theme_bw()

p_plantcad <- ggplot(plantcad_scores, aes(PC1, PC2, colour = qc_status)) +
  geom_point(alpha = 0.7) +
  scale_colour_manual(values = qc_colours) +
  labs(title = "PlantCAD", x = plantcad_ax$x, y = plantcad_ax$y, colour = "QC status") +
  theme_bw()

p_joint <- ggplot(joint_scores, aes(PC1, PC2, colour = qc_status)) +
  geom_point(alpha = 0.7) +
  scale_colour_manual(values = qc_colours) +
  labs(title = "Joint ESM-C + PlantCAD", x = joint_ax$x, y = joint_ax$y, colour = "QC status") +
  theme_bw()

p_esm + p_plantcad + p_joint + plot_layout(guides = "collect") +
  plot_annotation(title = "rbcS candidate embeddings — QC filter status")

plot(joint_scores)


# ============================================================================
# Within-species divergence in PC space (first 10 PCs), not raw 960-dim ESM.
# Raw embedding distances are noisy; the leading PCs keep the structured
# variation and drop per-dimension noise.
# ============================================================================
library(tidyverse)

n_pcs_for_distance <- 10

# ---- attach gene_id + species + QC flags to the PC scores --------------------
# prot_emb_pca$x rows are in the same order as rbcs, so we can bind directly.
stopifnot("PCA rows don't match rbcs rows" = nrow(prot_emb_pca$x) == nrow(rbcs))

pc_score_table <- as_tibble(prot_emb_pca$x[, 1:n_pcs_for_distance]) |>
  set_names(paste0("pc", 1:n_pcs_for_distance)) |>
  bind_cols(rbcs |> select(gene_id, species,
                           aln_filter_passed, fold_filter_passed))

# only copies passing BOTH filters contribute to within-species divergence
passing_pc_scores <- pc_score_table |>
  filter(aln_filter_passed, fold_filter_passed)

# ---- per species: max pairwise Euclidean distance in 10-PC space -------------
# Euclidean distance between copy i and j over the first 10 PCs:
#   d_ij = sqrt( sum_{k=1..10} (pc_ik - pc_jk)^2 )
pc_columns <- paste0("pc", 1:n_pcs_for_distance)

max_pc_distance_for_species <- function(one_species_df) {
  pc_matrix <- as.matrix(one_species_df[, pc_columns])
  if (nrow(pc_matrix) < 2) return(0)          # single copy = no divergence
  max(as.numeric(dist(pc_matrix, method = "euclidean")))
}

per_species_pc_divergence <-
  split(passing_pc_scores, passing_pc_scores$species) |>
  map_dfr(function(one_species_df) {
    tibble(
      species          = one_species_df$species[1],
      n_passing_copies = nrow(one_species_df),
      max_pc_distance  = max_pc_distance_for_species(one_species_df)
    )
  })

message("species with divergence values: ", nrow(per_species_pc_divergence))
print(per_species_pc_divergence)

# ---- join temperature (same accession -> binomial -> bio8 bridge as before) --
extract_accession <- function(species_string) {
  pieces <- str_split(species_string, "_", simplify = TRUE)
  str_c(pieces[, 1], "_", pieces[, 2])
}

accession_to_species <- read_csv(accession_species_map_path, show_col_types = FALSE)
colnames(accession_to_species)[[1]] <- "accession"
colnames(accession_to_species)[[2]] <- "species"

bio8_table           <- read_csv(bio8_csv_path, show_col_types = FALSE)

pc_model_table <- per_species_pc_divergence |>
  mutate(accession = extract_accession(species)) |>
  left_join(accession_to_species, by = "accession") |>
  rename(binomial = species.y) |>
  left_join(bio8_table, by = c("binomial" = "species")) |>
  filter(!is.na(bio8_median))

message("species in PC-distance model: ", nrow(pc_model_table))
stopifnot("fewer than 10 species" = nrow(pc_model_table) >= 10)

# ---- model: temp ~ max within-species divergence in 10-PC space --------------
# bio8_median     = median temp of wettest quarter (response)
# max_pc_distance = largest within-species pairwise distance over first 10 PCs
cat("\n=== temp ~ max within-species ESM divergence (first 10 PCs) ===\n")
pc_distance_model <- lm(bio8_median ~ max_pc_distance, data = pc_model_table)
print(summary(pc_distance_model))

# ---- plot: observed temp on Y, divergence on X ------------------------------
ggplot(pc_model_table, aes(x = max_pc_distance, y = bio8_median)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "max within-species divergence (first 10 ESM PCs, Euclidean)",
       y = "bio8 median (temp of wettest quarter)",
       title = "Temperature vs within-species rbcS divergence (PC space)") +
  theme_minimal()

plot(pc_model_table$n_passing_copies, pc_model_table$max_pc_distance)


library(tidyverse)
library(patchwork)
library(ggrepel)

# the `species` column in the scores tibbles is the accession string;
# pull the accession (first two underscore tokens) to bridge to bio8 + binomial
extract_accession <- function(species_string) {
  pieces <- str_split(species_string, "_", simplify = TRUE)
  str_c(pieces[, 1], "_", pieces[, 2])
}

# attach bio8 (continuous colour) and binomial (label) onto a scores tibble
attach_bio8 <- function(scores_df) {
  scores_df |>
    mutate(accession = extract_accession(species)) |>
    left_join(accession_to_species, by = "accession") |>
    rename(binomial = species.y) |>
    left_join(bio8_table, by = c("binomial" = "species"))
}

esm_bio8      <- attach_bio8(esm_scores)
plantcad_bio8 <- attach_bio8(plantcad_scores)
joint_bio8    <- attach_bio8(joint_scores)

# one label per species, placed at its centroid in PC space, so repeated
# copies don't stamp the same binomial dozens of times
species_label_positions <- function(scores_bio8_df) {
  scores_bio8_df |>
    group_by(binomial) |>
    summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")
}

make_bio8_panel <- function(scores_bio8_df, ax, title) {
  ggplot(scores_bio8_df, aes(PC1, PC2)) +
    geom_point(aes(colour = bio8_median), alpha = 0.8) +
    #geom_text_repel(
    #  data = species_label_positions(scores_bio8_df),
    #  aes(label = binomial),
    #  size = 2.5, max.overlaps = Inf, segment.size = 0.2
    #) +
    scale_colour_viridis_c() +
    labs(title = title, x = ax$x, y = ax$y, colour = "bio8") +
    theme_bw()
}

p_esm_bio8      <- make_bio8_panel(esm_bio8,      esm_ax,      "ESM-C")
p_plantcad_bio8 <- make_bio8_panel(plantcad_bio8, plantcad_ax, "PlantCAD")
p_joint_bio8    <- make_bio8_panel(joint_bio8,    joint_ax,    "Joint ESM-C + PlantCAD")
p_esm_bio8
p_plantcad_bio8
p_joint_bio8

plot(esm_bio8$PC2, esm_bio8$bio8_median)
plot(p_plantcad_bio8$PC2, p_plantcad_bio8$bio8_median)


p_esm_bio8 + p_plantcad_bio8 + p_joint_bio8 +
  plot_layout(guides = "collect") +
  plot_annotation(title = "rbcS candidate embeddings — coloured by species bio8")

esm_bio8 <- na.omit(esm_bio8)
key <- esm_bio8$aln_filter_passed & esm_bio8$fold_filter_passed

plot(esm_bio8[key,]$PC1,
     esm_bio8[key,]$bio8_median)


cor(esm_bio8[key,]$PC1,
     esm_bio8[key,]$bio8_median)
cor(esm_bio8[key,]$PC2,
    esm_bio8[key,]$bio8_median)
LSD::heatscatter(esm_bio8[key,]$PC2,
    esm_bio8[key,]$bio8_median)

plantcad_bio8 <- na.omit(plantcad_bio8)
LSD::heatscatter(plantcad_bio8[key,]$PC2,
     plantcad_bio8[key,]$bio8_median,
     xlab="plantcad2 embeddings PC2 (16.6%) ", ylab ="bio8",
     main="Plantcad2 embeddings on 765 rbcS in 95 angiosperms")
correl <- round(cor(plantcad_bio8[key,]$PC2,plantcad_bio8[key,]$bio8_median),2)
text(30,25,paste("r=",correl), col="coral")

cor(plantcad_bio8[key,]$PC1,
    plantcad_bio8[key,]$bio8_median)
cor(plantcad_bio8[key,]$PC2,
    plantcad_bio8[key,]$bio8_median)


LSD::heatscatter(joint_bio8[key,]$PC2,
                 joint_bio8[key,]$bio8_median)


joint_bio8 <- na.omit(joint_bio8)
cor(joint_bio8[key,]$PC1,
    joint_bio8[key,]$bio8_median)
cor(joint_bio8[key,]$PC2,
    joint_bio8[key,]$bio8_median)


cor(esm_bio8[key,]$PC2,
    plantcad_bio8[key,]$PC2)
plot(esm_bio8[key,]$PC1,
    plantcad_bio8[key,]$PC1)
