#!/usr/bin/env Rscript
# test_structural_colocalization.R
#
# Do annotated genes hit by structural (photosynthetic) peptides cluster along
# chromosomes more than expected, using ALL annotated genes as the background?
# Clustering is the NUPT co-association signature.
#
# Statistic (per genome): order genes along each chromosome; a "structural gene"
# is one labelled with a structure category. Count pairs of structural genes that
# are CONSECUTIVE in structural order and lie within `window_in_genes` gene-ranks
# of each other -- i.e. how many times two structural genes sit close together.
# Null: randomly relabel which genes are structural, keeping the SAME count per
# chromosome (so the background is the real gene layout). Compare observed to null.
#
# Per-genome outputs: observed statistic, null mean/sd, z-score, fold-enrichment,
# permutation p-value. The z-score is the per-genome colocalization measure that
# gets tested against the phenotype at the end.

library(tidyverse)

results_directory   <- "results"
bio8_csv_path       <- "data/genecad_species_bio8_medians.csv"
window_in_genes     <- 5      # two structural genes within this many gene-ranks = "close"
number_of_permutations <- 1000
minimum_structural_genes <- 5 # genomes with fewer are too low-power to test
set.seed(1)

# categories that are NOT structural (the housekeeping control + non-hits)
non_structural_categories <- c("null", "uncategorised", "unknown_reference")

# ============================================================================
# 1. load positions + labels per genome
# ============================================================================
position_file_paths <- list.files(
  results_directory, pattern = "^gene_positions\\.tsv$",
  recursive = TRUE, full.names = TRUE
)
stopifnot("no gene_positions.tsv files found (run extract_gene_positions.py)" =
            length(position_file_paths) > 0)

load_one_genome <- function(position_file_path) {
  species_folder <- dirname(position_file_path)
  sample_name    <- basename(species_folder)
  category_path  <- file.path(species_folder, "gene_category.tsv")
  if (!file.exists(category_path)) return(NULL)  # not labelled yet
  
  positions  <- read_tsv(position_file_path, show_col_types = FALSE,
                         col_types = cols(chromosome = col_character()))
  categories <- read_tsv(category_path, show_col_types = FALSE) %>%
    select(gene_id, category)
  
  positions %>%
    left_join(categories, by = "gene_id") %>%
    mutate(
      sample_name    = sample_name,
      is_structural  = !is.na(category) & !(category %in% non_structural_categories)
    )
}

all_genomes <- map(position_file_paths, load_one_genome) %>% compact() %>% bind_rows()
stopifnot("no genomes had both positions and labels" = nrow(all_genomes) > 0)

# ============================================================================
# 2. clumping statistic + permutation null, per genome
# ============================================================================
# count "close" consecutive structural pairs given the structural genes' ranks
count_close_consecutive_pairs <- function(structural_ranks, window) {
  if (length(structural_ranks) < 2) return(0L)
  sorted_ranks <- sort(structural_ranks)
  gaps_between_consecutive <- diff(sorted_ranks)
  sum(gaps_between_consecutive <= window)
}

# observed statistic for one genome: rank genes within chromosome by position
compute_observed_statistic <- function(one_genome, window) {
  one_genome %>%
    group_by(chromosome) %>%
    arrange(start, .by_group = TRUE) %>%
    mutate(rank_on_chromosome = row_number()) %>%
    summarise(
      close_pairs = count_close_consecutive_pairs(
        rank_on_chromosome[is_structural], window),
      .groups = "drop"
    ) %>%
    summarise(total = sum(close_pairs)) %>%
    pull(total)
}

# one null draw: per chromosome, randomly pick the same NUMBER of structural genes
# from that chromosome's genes, then recompute the statistic
one_null_statistic <- function(structural_count_by_chromosome,
                               gene_count_by_chromosome, window) {
  total <- 0L
  for (chromosome_name in names(structural_count_by_chromosome)) {
    number_of_genes      <- gene_count_by_chromosome[[chromosome_name]]
    number_of_structural <- structural_count_by_chromosome[[chromosome_name]]
    if (number_of_structural >= 2) {
      random_ranks <- sample.int(number_of_genes, number_of_structural)
      total <- total + count_close_consecutive_pairs(random_ranks, window)
    }
  }
  total
}

analyse_one_genome <- function(one_genome, window, n_perm) {
  total_structural <- sum(one_genome$is_structural)
  
  # per-chromosome counts drive the null (preserve where structural genes are)
  per_chromosome <- one_genome %>%
    group_by(chromosome) %>%
    summarise(n_genes = n(),
              n_structural = sum(is_structural), .groups = "drop")
  gene_count_by_chromosome       <- setNames(per_chromosome$n_genes, per_chromosome$chromosome)
  structural_count_by_chromosome <- setNames(per_chromosome$n_structural, per_chromosome$chromosome)
  
  observed <- compute_observed_statistic(one_genome, window)
  null_draws <- replicate(
    n_perm,
    one_null_statistic(structural_count_by_chromosome, gene_count_by_chromosome, window)
  )
  
  null_mean <- mean(null_draws)
  null_sd   <- sd(null_draws)
  tibble(
    sample_name        = one_genome$sample_name[1],
    n_genes            = nrow(one_genome),
    n_structural       = total_structural,
    observed_close_pairs = observed,
    null_mean          = null_mean,
    null_sd            = null_sd,
    # z-score: how many SDs above the random expectation (the colocalization score)
    colocalization_z   = if (null_sd > 0) (observed - null_mean) / null_sd else NA_real_,
    # fold: observed / expected close pairs
    fold_enrichment    = if (null_mean > 0) observed / null_mean else NA_real_,
    # one-sided permutation p-value (clustering = more close pairs than null)
    p_value            = (1 + sum(null_draws >= observed)) / (1 + n_perm)
  )
}

colocalization_by_genome <- all_genomes %>%
  group_split(sample_name) %>%
  map(~ analyse_one_genome(.x, window_in_genes, number_of_permutations)) %>%
  bind_rows() %>%
  arrange(desc(colocalization_z))

write_tsv(colocalization_by_genome,
          file.path(results_directory, "colocalization_by_genome.tsv"))
print(colocalization_by_genome, n = Inf)

# ============================================================================
# 3. visualise the fair background for one genome (most-clustered by default)
# ============================================================================
# local structural fraction along each chromosome vs the genome-wide fraction
# (the dashed line). Peaks above the line = local structural clustering.
plot_one_genome_enrichment <- function(sample_to_plot, bin_width_nt = 2e6) {
  one_genome <- all_genomes %>% filter(sample_name == sample_to_plot)
  genome_wide_fraction <- mean(one_genome$is_structural)
  
  binned <- one_genome %>%
    mutate(position_bin = floor(start / bin_width_nt) * bin_width_nt) %>%
    group_by(chromosome, position_bin) %>%
    summarise(structural_fraction = mean(is_structural),
              n_genes = n(), .groups = "drop") %>%
    filter(n_genes >= 5)   # ignore sparse bins
  
  ggplot(binned, aes(x = position_bin / 1e6, y = structural_fraction)) +
    geom_col() +
    geom_hline(yintercept = genome_wide_fraction, linetype = "dashed") +
    facet_wrap(~ chromosome, scales = "free_x") +
    labs(x = "position (Mb)", y = "fraction of genes that are structural",
         title = str_c("Structural-gene density vs all-gene background: ",
                       sample_to_plot),
         subtitle = "dashed line = genome-wide structural fraction") +
    theme_minimal()
}

most_clustered_sample <- colocalization_by_genome %>%
  filter(n_structural >= minimum_structural_genes) %>%
  slice_max(colocalization_z, n = 1) %>%
  pull(sample_name)
if (length(most_clustered_sample) == 1) {
  print(plot_one_genome_enrichment(most_clustered_sample))
}

# summary across genomes: which show significant clustering
plot_z_across_genomes <-
  colocalization_by_genome %>%
  filter(n_structural >= minimum_structural_genes) %>%
  ggplot(aes(x = colocalization_z, y = fct_reorder(sample_name, colocalization_z))) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "colocalization z-score (structural clustering vs random)",
       y = NULL, title = "Structural-gene clustering per genome") +
  theme_minimal()
print(plot_z_across_genomes)

# ============================================================================
# 4. associate colocalization with the phenotype (bio8)
# ============================================================================
make_key_from_sample <- function(sample_name) {
  if (str_detect(sample_name, "_")) {
    pieces <- str_split(sample_name, "_", simplify = TRUE)
    genus <- pieces[, 1]; epithet <- pieces[, 2]
  } else {
    genus <- str_sub(sample_name, 1, 1); epithet <- str_sub(sample_name, 2)
  }
  str_c(str_to_lower(str_sub(genus, 1, 1)), "_", str_to_lower(epithet))
}
make_key_from_binomial <- function(binomial) {
  pieces <- str_split(binomial, " ", simplify = TRUE)
  str_c(str_to_lower(str_sub(pieces[, 1], 1, 1)), "_", str_to_lower(pieces[, 2]))
}

bio8_table <- read_csv(bio8_csv_path, show_col_types = FALSE) %>%
  mutate(join_key = make_key_from_binomial(species))

association_table <- colocalization_by_genome %>%
  filter(n_structural >= minimum_structural_genes, !is.na(colocalization_z)) %>%
  mutate(join_key = map_chr(sample_name, make_key_from_sample)) %>%
  inner_join(bio8_table, by = "join_key") %>%
  filter(!is.na(bio8_median))

message("genomes in the colocalization-vs-bio8 test: ", nrow(association_table))

plot_coloc_vs_temperature <-
  ggplot(association_table, aes(x = bio8_median, y = colocalization_z)) +
  geom_point() + geom_smooth(method = "lm", se = TRUE) +
  geom_text(aes(label = species), size = 2.3, vjust = -0.6, check_overlap = TRUE) +
  labs(x = "bio8 median (mean temp of wettest quarter, deg C)",
       y = "structural colocalization z-score",
       title = "Structural-gene clustering vs temperature") +
  theme_minimal()
print(plot_coloc_vs_temperature)

cat("\n--- colocalization z vs bio8 (Spearman) ---\n")
print(cor.test(association_table$bio8_median,
               association_table$colocalization_z, method = "spearman"))
cat("\n--- linear model: colocalization_z ~ bio8_median ---\n")
print(summary(lm(colocalization_z ~ bio8_median, data = association_table)))

# ---------------------------------------------------------------------------
# CAVEATS:
#  * only annotated genes are tested -> unannotated NUPTs are invisible here
#    (that is the price of the fair within-annotation background).
#  * clustering can be real photosynthetic tandem arrays OR NUPT insertions; to
#    isolate NUPTs, re-run with `non_structural_categories` set so only
#    PLASTID-ENCODED reference subunits count as structural.
#  * window_in_genes and bin_width_nt are tunable; check the per-genome plot.
#  * phylogenetic non-independence still applies to the bio8 association.
# ---------------------------------------------------------------------------
message("\nNOTE: phylogenetic non-independence applies to the bio8 association.")

# ============================================================================
# count discrete colocalized BLOCKS per genome and test against bio8
# a block = a maximal run of structural genes where each consecutive pair
# is within window_in_genes gene-ranks of the next
# ============================================================================

count_colocalized_blocks <- function(structural_ranks, window) {
  if (length(structural_ranks) < 2) return(0L)
  sorted_ranks <- sort(structural_ranks)
  gaps <- diff(sorted_ranks)
  # a new block starts wherever the gap exceeds the window
  # number of blocks = number of "breaks" + 1, but only if there are >= 2
  # structural genes close enough to form at least one block
  in_block <- gaps <= window
  if (!any(in_block)) return(0L)
  # count transitions from FALSE to TRUE (new block starts) plus the first block
  transitions <- sum(diff(c(FALSE, in_block)) == 1)
  transitions
}

blocks_per_genome <- all_genomes %>%
  group_by(sample_name) %>%
  group_modify(~ {
    one_genome <- .x
    # rank genes within each chromosome by position
    ranked <- one_genome %>%
      group_by(chromosome) %>%
      arrange(start, .by_group = TRUE) %>%
      mutate(rank_on_chromosome = row_number()) %>%
      ungroup()
    
    # count blocks per chromosome, sum across genome
    total_blocks <- ranked %>%
      group_by(chromosome) %>%
      summarise(
        blocks_on_chromosome = count_colocalized_blocks(
          rank_on_chromosome[is_structural], window_in_genes),
        .groups = "drop"
      ) %>%
      summarise(total_colocalized_blocks = sum(blocks_on_chromosome)) %>%
      pull(total_colocalized_blocks)
    
    tibble(
      n_structural           = sum(one_genome$is_structural),
      total_colocalized_blocks = total_blocks
    )
  }) %>%
  ungroup()

print(blocks_per_genome %>% arrange(desc(total_colocalized_blocks)))

# join bio8
blocks_analysis_table <- blocks_per_genome %>%
  filter(n_structural >= minimum_structural_genes) %>%
  mutate(join_key = map_chr(sample_name, make_key_from_sample)) %>%
  inner_join(bio8_table, by = "join_key") %>%
  filter(!is.na(bio8_median))

message("genomes in blocks vs bio8 test: ", nrow(blocks_analysis_table))
hist(blocks_analysis_table$total_colocalized_blocks)
# block count is a count variable -- use Spearman and also a Poisson GLM
# Poisson GLM: log(n_blocks) ~ bio8 + offset(log(n_structural))
# the offset controls for the fact that more structural genes = more chance
# of forming blocks regardless of biology
plot_blocks_vs_temperature <-
  ggplot(blocks_analysis_table,
         aes(x = bio8_median, y = total_colocalized_blocks)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "poisson"), se = TRUE) +
  geom_text(aes(label = species), size = 2.3,
            vjust = -0.6, check_overlap = TRUE) +
  labs(x = "bio8 median (mean temp of wettest quarter, deg C)",
       y = "number of colocalized structural-gene blocks",
       title = "Colocalized blocks vs temperature") +
  theme_minimal()
print(plot_blocks_vs_temperature)

cat("\n--- colocalized blocks vs bio8 (Spearman) ---\n")
print(cor.test(blocks_analysis_table$bio8_median,
               blocks_analysis_table$total_colocalized_blocks,
               method = "spearman"))

cat("\n--- Poisson GLM: n_blocks ~ bio8 + offset(log(n_structural)) ---\n")
print(summary(
  glm(total_colocalized_blocks ~ bio8_median +
        offset(log(n_structural)),
      family = poisson,
      data = blocks_analysis_table)
))

# is detection rate (n_structural) itself associated with temperature?
cat("\n--- n_structural vs bio8 (Spearman) ---\n")
print(cor.test(blocks_analysis_table$bio8_median,
               blocks_analysis_table$n_structural,
               method = "spearman"))

ggplot(blocks_analysis_table, aes(y = bio8_median, x = n_structural)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  geom_text(aes(label = species), size = 2.3,
            vjust = -0.6, check_overlap = TRUE) +
  labs(y = "bio8 median", y = "number of structural genes detected",
       title = "Are more Photosynthetic Complex genes detected in warmer species?") +
  theme_minimal()

#!/usr/bin/env Rscript
# test_miniprot_nupt_vs_bio8.R
#
# Uses miniprot hits (ALL genomic loci, annotated or not) to test:
#   1. Do structural peptide copy numbers associate with temperature?
#   2. Do structural peptides have more copies than housekeeping (paired)?
#   3. Does copy number vary by specific complex (PDB structure)?
#
# The housekeeping copy count is the within-genome control for genome size /
# assembly quality / annotation completeness effects.

library(tidyverse)

hits_summary_path    <- "results/gene_hits/hits_summary.tsv"
gene_categories_path <- "data/reference_protein_sequences/gene_categories.tsv"
metadata_path        <- "data/hifi_angiosperm_metadata/hifi_representative_newest_per_species.csv"
bio8_csv_path        <- "data/hifi_species_pheno.csv"

housekeeping_category_label <- "housekeeping_candidates_Zm.peptides"

# ============================================================================
# 1. load hits and attach categories
# ============================================================================
hits <- read_tsv(hits_summary_path, show_col_types = FALSE)

# sanitise gene_name to match gene_categories.tsv convention (| -> _)
hits <- hits %>%
  mutate(gene_name_sanitised = str_replace_all(gene_name, "\\|", "_")
         %>% str_replace_all(" ", "_"))

categories <- read_tsv(gene_categories_path, show_col_types = FALSE) %>%
  select(gene_name, category) %>%
  rename(gene_name_sanitised = gene_name)

hits <- hits %>%
  left_join(categories, by = "gene_name_sanitised") %>%
  mutate(
    peptide_class = case_when(
      is.na(category)                      ~ "unlabelled",
      category == housekeeping_category_label ~ "housekeeping",
      TRUE                                 ~ "structural"
    )
  )

message("total hit loci: ", nrow(hits))
message("class breakdown:")
print(count(hits, peptide_class))

# ============================================================================
# 2. per species x peptide: total copy number (= max copy_number per peptide,
#    which is how miniprot counts copies)
# ============================================================================
per_species_per_peptide <- hits %>%
  filter(peptide_class != "unlabelled") %>%
  group_by(assembly_accession, gene_name_sanitised, category, peptide_class) %>%
  summarise(copy_number = max(copy_number), .groups = "drop")

# per species x class: median copy number across peptides in that class
per_species_class_summary <- per_species_per_peptide %>%
  group_by(assembly_accession, peptide_class) %>%
  summarise(
    median_copy_number = median(copy_number),
    mean_copy_number   = mean(copy_number),
    total_copies       = sum(copy_number),
    n_peptides         = n(),
    .groups = "drop"
  )

# per species x category (PDB complex): median copy number
per_species_category_summary <- per_species_per_peptide %>%
  group_by(assembly_accession, category) %>%
  summarise(
    median_copy_number = median(copy_number),
    n_peptides         = n(),
    .groups = "drop"
  )

# ============================================================================
# 3. join accession -> species name -> bio8
# ============================================================================
metadata <- read_csv(metadata_path, show_col_types = FALSE) %>%
  select(assembly_accession = `Assembly Accession`,
         organism_name      = `Organism Name`) %>%
  # drop version suffix (.1) for a more robust key
  mutate(assembly_accession_base = str_remove(assembly_accession, "\\.\\d+$"))

bio8_table <- read_csv(bio8_csv_path, show_col_types = FALSE) %>%
  mutate(
    # match on (genus-initial, epithet) key
    join_key = {
      pieces <- str_split(species, " ", simplify = TRUE)
      str_c(str_to_lower(str_sub(pieces[, 1], 1, 1)), "_",
            str_to_lower(pieces[, 2]))
    }
  )

# build accession -> bio8 lookup via organism name
accession_to_bio8 <- metadata %>%
  mutate(
    join_key = {
      pieces <- str_split(organism_name, " ", simplify = TRUE)
      str_c(str_to_lower(str_sub(pieces[, 1], 1, 1)), "_",
            str_to_lower(pieces[, 2]))
    }
  ) %>%
  inner_join(bio8_table, by = "join_key") %>%
  select(assembly_accession, assembly_accession_base, species, bio8_median)

# report unmatched assemblies
unmatched <- per_species_class_summary %>%
  filter(!assembly_accession %in% accession_to_bio8$assembly_accession) %>%
  pull(assembly_accession) %>% unique()
if (length(unmatched) > 0) {
  message("WARNING: ", length(unmatched),
          " assemblies did not match bio8 (check metadata):")
  print(unmatched)
}

# ============================================================================
# 4. TEST 1: structural vs housekeeping copy number, paired across species
# ============================================================================
# wide form: one row per species, columns = median copy number per class
class_wide <- per_species_class_summary %>%
  select(assembly_accession, peptide_class, median_copy_number) %>%
  pivot_wider(names_from = peptide_class, values_from = median_copy_number) %>%
  filter(!is.na(structural), !is.na(housekeeping))

cat("\n--- TEST 1: structural vs housekeeping copy number (paired Wilcoxon) ---\n")
cat("n species with both classes: ", nrow(class_wide), "\n")
print(wilcox.test(class_wide$structural, class_wide$housekeeping, paired = TRUE))

# visualise the paired difference
class_wide %>%
  mutate(structural_excess = structural - housekeeping,
         species = accession_to_bio8$species[
           match(assembly_accession, accession_to_bio8$assembly_accession)]) %>%
  ggplot(aes(x = fct_reorder(species, structural_excess),
             y = structural_excess)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(x = NULL, y = "structural minus housekeeping median copy number",
       title = "Structural genes have more copies than housekeeping?") +
  theme_minimal(base_size = 8) -> p_paired
print(p_paired)

# ============================================================================
# 5. TEST 2: copy number vs temperature, separately for structural/housekeeping
# ============================================================================
analysis_table <- per_species_class_summary %>%
  inner_join(accession_to_bio8, by = "assembly_accession") %>%
  filter(!is.na(bio8_median))

message("species in temperature association test: ",
        n_distinct(analysis_table$assembly_accession))

for (cls in c("structural", "housekeeping")) {
  class_data <- analysis_table %>% filter(peptide_class == cls)
  cat("\n--- TEST 2:", cls, "copy number vs bio8 (Spearman) ---\n")
  print(cor.test(class_data$bio8_median,
                 class_data$median_copy_number, method = "spearman"))
}

# scatter: structural (colour) and housekeeping (grey) on same plot
analysis_table %>%
  filter(peptide_class %in% c("structural", "housekeeping")) %>%
  ggplot(aes(x = bio8_median, y = median_copy_number,
             colour = peptide_class)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_colour_manual(values = c("structural" = "steelblue",
                                 "housekeeping" = "grey50")) +
  labs(x = "bio8 median (mean temp of wettest quarter, deg C)",
       y = "median peptide copy number",
       colour = NULL,
       title = "Miniprot copy number vs temperature",
       subtitle = "structural (all complexes pooled) vs housekeeping") +
  theme_minimal() -> p_scatter
print(p_scatter)

# ============================================================================
# 6. TEST 3: per-complex (PDB category) copy number vs temperature
# ============================================================================
# each complex tested separately; BH correction across complexes
run_complex_model <- function(cat_name, data) {
  cat_data <- data %>%
    filter(category == cat_name) %>%
    inner_join(accession_to_bio8, by = "assembly_accession") %>%
    filter(!is.na(bio8_median))
  if (nrow(cat_data) < 10) {
    return(tibble(category = cat_name, n_species = nrow(cat_data),
                  estimate = NA_real_, p_value = NA_real_, note = "too few"))
  }
  model_fit  <- lm(bio8_median ~ median_copy_number, data = cat_data)
  model_tidy <- broom::tidy(model_fit) %>% filter(term == "median_copy_number")
  tibble(category  = cat_name,
         n_species = nrow(cat_data),
         estimate  = model_tidy$estimate,
         std_error = model_tidy$std.error,
         p_value   = model_tidy$p.value,
         note      = "ok")
}

structural_categories <- per_species_category_summary %>%
  filter(category != housekeeping_category_label) %>%
  pull(category) %>% unique()

per_complex_results <- map(structural_categories, run_complex_model,
                           data = per_species_category_summary) %>%
  bind_rows() %>%
  filter(note == "ok") %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

cat("\n--- TEST 3: per-complex copy number vs bio8 ---\n")
print(per_complex_results)

per_complex_results %>%
  filter(!is.na(estimate)) %>%
  ggplot(aes(x = estimate,
             y = fct_reorder(category, estimate),
             alpha = p_adjusted < 0.10)) +
  geom_point() +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std_error,
                     xmax = estimate + 1.96 * std_error),
                 height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.35),
                     labels = c("TRUE" = "FDR < 0.10", "FALSE" = "ns"),
                     name = NULL) +
  labs(x = "coefficient: effect of copy number on bio8",
       y = "complex (PDB id)",
       title = "Per-complex miniprot copy number vs temperature") +
  theme_minimal() -> p_complex
print(p_complex)

per_species_category_summary %>%
  filter(category == "1RCX") %>%
  inner_join(accession_to_bio8, by = "assembly_accession") %>%
  filter(!is.na(bio8_median)) %>%
  ggplot(aes(x = bio8_median, y = median_copy_number)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  geom_text(aes(label = species), size = 2.5,
            vjust = -0.6, check_overlap = TRUE) +
  labs(x = "bio8 median (mean temp of wettest quarter, deg C)",
       y = "median RuBisCO SSu (1RCX) peptide copy number",
       title = "RuBisCO SSu copy number vs temperature") +
  theme_minimal()

# ---------------------------------------------------------------------------
# NOTE: species phylogenetically non-independent; PGLS is the proper test.
# Miniprot copy_number counts ALL loci including pseudogene fragments and
# NUPTs -- not just functional genes. This is the feature, not a bug, if
# you want total NUPT burden; but interpret vs the annotated-gene results.
# ---------------------------------------------------------------------------
message("\nNOTE: phylogenetic non-independence applies to all tests above.")