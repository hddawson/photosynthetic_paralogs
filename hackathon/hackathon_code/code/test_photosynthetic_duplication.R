#!/usr/bin/env Rscript
# test_photosynthetic_duplication.R
#
# Hypothesis: photosynthetic gene FAMILIES are more duplicated in warmer-climate
# species, over and above the genome-wide duplication background. Duplication is
# measured by family size from clustering (mechanism-agnostic: WGD/tandem/etc.
# are all lumped, which is the point).
#
# Design:
#   * label genome genes by category (photosynthetic complex vs housekeeping)
#   * a clustering FAMILY is "photosynthetic" if any member is photosynthetic,
#     "housekeeping" if any member is housekeeping  (dedup so big families aren't
#     counted once per labelled member)
#   * per species: summarise photosynthetic vs housekeeping family sizes
#   * housekeeping is the within-genome control; the photosynthetic-minus-
#     housekeeping contrast removes each genome's background duplication level
#   * test that contrast (and the raw photosynthetic size) against bio8

library(tidyverse)

results_directory <- "results"
bio8_csv_path     <- "data/genecad_species_bio8_medians.csv"

# ---- which category labels count as photosynthetic vs the control -----------
# 'null' is the housekeeping control set. Everything else in gene_categories.tsv
# is a photosynthetic complex (PDB id) or rubisco activase. Edit if you want to
# regroup (e.g. drop the chloroplast ribosome 6ERI or PEP polymerase 8XZV).
housekeeping_category_label <- "null"
# (photosynthetic = any labelled category that is not the housekeeping one)

# ============================================================================
# 1. load per-species gene labels + family sizes
#    unit of analysis: reference PEPTIDE x SPECIES
#    value: family_size of the gene that best-hit that peptide
#           (if a peptide hits >1 family in a species, flag it, take largest)
# ============================================================================
gene_category_paths <- list.files(
  results_directory, pattern = "^gene_category\\.tsv$",
  recursive = TRUE, full.names = TRUE
)
stopifnot("no gene_category.tsv files found" = length(gene_category_paths) > 0)
message("found labels for ", length(gene_category_paths), " species")

load_one_species <- function(gene_category_path) {
  species_folder <- dirname(gene_category_path)
  sample_name    <- basename(species_folder)
  
  per_gene_path <- file.path(species_folder, "per_gene_duplication.tsv")
  stopifnot("missing per_gene_duplication.tsv" = file.exists(per_gene_path))
  
  # gene_category.tsv has one row per GENOME gene; reference_hit = which peptide
  # per_gene_duplication.tsv has gene_id -> family_id, family_size
  categories <- read_tsv(gene_category_path, show_col_types = FALSE) %>%
    filter(category != housekeeping_category_label) %>%   # drop housekeeping
    select(gene_id, reference_hit, category)
  
  family_sizes <- read_tsv(per_gene_path, show_col_types = FALSE) %>%
    select(gene_id, family_id, family_size)
  
  # join: each labelled gene gets its family_size
  categories %>%
    left_join(family_sizes, by = "gene_id") %>%
    mutate(sample_name = sample_name)
}

per_gene_labelled <- map(gene_category_paths, load_one_species) %>% bind_rows()

# ============================================================================
# 2. per species x reference_peptide: take largest family size,
#    flag cases where a peptide hit >1 distinct family
# ============================================================================
multi_family_flags <- per_gene_labelled %>%
  group_by(sample_name, reference_hit) %>%
  summarise(
    n_distinct_families  = n_distinct(family_id),
    max_family_size      = max(family_size, na.rm = TRUE),
    has_multiple_families = n_distinct(family_id) > 1,
    .groups = "drop"
  )

# report peptides that hit multiple families in at least one species
multi_hits <- multi_family_flags %>%
  filter(has_multiple_families) %>%
  count(reference_hit, name = "n_species_with_multiple_families") %>%
  arrange(desc(n_species_with_multiple_families))
if (nrow(multi_hits) > 0) {
  message("peptides hitting >1 family in at least one species (largest taken):")
  print(multi_hits, n = Inf)
}

# ============================================================================
# 3. pivot to wide: one column per reference peptide, one row per species
# ============================================================================
peptide_family_size_wide <- multi_family_flags %>%
  select(sample_name, reference_hit, max_family_size) %>%
  pivot_wider(names_from  = reference_hit,
              values_from = max_family_size,
              values_fill = NA)   # NA = that peptide hit nothing in that species

message("predictor matrix: ", nrow(peptide_family_size_wide), " species x ",
        ncol(peptide_family_size_wide) - 1, " peptide predictors")

# ============================================================================
# 4. join climate + genome-wide duplication fraction
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

genome_wide_fraction <- map(
  list.files(results_directory, pattern = "^per_gene_duplication\\.tsv$",
             recursive = TRUE, full.names = TRUE),
  function(p) {
    s <- basename(dirname(p))
    d <- read_tsv(p, show_col_types = FALSE)
    tibble(sample_name = s, fraction_duplicated = mean(d$is_duplicated))
  }
) %>% bind_rows()

analysis_table <- peptide_family_size_wide %>%
  mutate(join_key = map_chr(sample_name, make_key_from_sample)) %>%
  inner_join(bio8_table, by = "join_key") %>%
  left_join(genome_wide_fraction, by = "sample_name") %>%
  filter(!is.na(bio8_median))

message("species in the model: ", nrow(analysis_table))

# ============================================================================
# 5. MLR: bio8_median ~ all peptide family sizes + fraction_duplicated
# ============================================================================
# drop peptide columns that are NA in >50% of species (too sparse to be useful)
peptide_columns <- setdiff(colnames(peptide_family_size_wide), "sample_name")
coverage_of_peptide <- map_dbl(peptide_columns, ~ mean(!is.na(analysis_table[[.x]])))
peptide_columns_kept <- peptide_columns[coverage_of_peptide >= 0.5]
peptide_columns_dropped <- peptide_columns[coverage_of_peptide < 0.5]
if (length(peptide_columns_dropped) > 0) {
  message("dropping ", length(peptide_columns_dropped),
          " peptides present in <50% of species: ",
          paste(peptide_columns_dropped, collapse = ", "))
}
message("peptides kept as predictors: ", length(peptide_columns_kept))

# complete-case analysis: drop species with any NA in kept predictors
model_data <- analysis_table %>%
  select(species, bio8_median, fraction_duplicated,
         all_of(peptide_columns_kept)) %>%
  drop_na()
message("species after dropping NAs: ", nrow(model_data))
stopifnot("fewer than 10 species remain after NA removal" = nrow(model_data) >= 10)

# build the formula programmatically
predictor_terms <- c(peptide_columns_kept, "fraction_duplicated")
model_formula   <- as.formula(
  paste("bio8_median ~", paste(predictor_terms, collapse = " + "))
)
cat("\nfitting: "); print(model_formula)

full_model <- lm(model_formula, data = model_data)
cat("\n--- MLR: bio8 ~ peptide family sizes + fraction_duplicated ---\n")
print(summary(full_model))

# variance inflation factors: values > 10 flag problematic collinearity
if (requireNamespace("car", quietly = TRUE)) {
  cat("\n--- VIF (collinearity check; >10 is a problem) ---\n")
  print(car::vif(full_model))
} else {
  message("install 'car' for VIF: install.packages('car')")
}

# coefficient plot so you can see which peptides drive the association
coef_table <- broom::tidy(full_model) %>%
  filter(term != "(Intercept)") %>%
  mutate(significant = p.value < 0.05)

plot_coefficients <-
  ggplot(coef_table,
         aes(x = estimate,
             y = fct_reorder(term, estimate),
             colour = significant)) +
  geom_point() +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std.error,
                     xmax = estimate + 1.96 * std.error),
                 height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "steelblue")) +
  labs(x = "coefficient (effect on bio8_median)",
       y = NULL, colour = "p < 0.05",
       title = "MLR coefficients: which peptide family sizes predict temperature?") +
  theme_minimal()
print(plot_coefficients)

# ---------------------------------------------------------------------------
# CAVEATS: species phylogenetically non-independent; PGLS is the proper test.
# With ~47 species and 9+ predictors the model is borderline on df.
# Collinear complexes (VIF > 5) should be combined or one dropped.
# ---------------------------------------------------------------------------
message("\nNOTE: phylogenetic non-independence still applies.")

# ============================================================================
# 6. per-peptide univariate tests
#    for each peptide: fit bio8 ~ family_size + fraction_duplicated
#    only in species where that peptide was detected (non-NA)
# ============================================================================
library(broom)

run_one_peptide_model <- function(peptide_name, data) {
  peptide_data <- data %>%
    select(species, bio8_median, fraction_duplicated,
           family_size = all_of(peptide_name)) %>%
    drop_na()
  
  number_of_species <- nrow(peptide_data)
  if (number_of_species < 10) {
    return(tibble(peptide = peptide_name,
                  n_species = number_of_species,
                  estimate = NA_real_, std_error = NA_real_,
                  statistic = NA_real_, p_value = NA_real_,
                  note = "too few species"))
  }
  
  model_fit <- lm(bio8_median ~ family_size + fraction_duplicated,
                  data = peptide_data)
  model_tidy <- tidy(model_fit) %>% filter(term == "family_size")
  
  tibble(
    peptide    = peptide_name,
    n_species  = number_of_species,
    estimate   = model_tidy$estimate,
    std_error  = model_tidy$std.error,
    statistic  = model_tidy$statistic,
    p_value    = model_tidy$p.value,
    note       = "ok"
  )
}

per_peptide_results <- map(peptide_columns_kept, run_one_peptide_model,
                           data = analysis_table) %>%
  bind_rows() %>%
  filter(note == "ok") %>%
  # Benjamini-Hochberg correction across all peptide tests
  mutate(p_adjusted = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

message("peptides tested: ", nrow(per_peptide_results))
message("peptides significant after BH correction (FDR < 0.10): ",
        sum(per_peptide_results$p_adjusted < 0.10, na.rm = TRUE))
print(per_peptide_results, n = Inf)

# join back the category (PDB complex) for context in the plot
category_lookup <- read_tsv(
  list.files(results_directory, pattern = "^gene_category\\.tsv$",
             recursive = TRUE, full.names = TRUE)[1],
  show_col_types = FALSE) %>%
  select(reference_hit, category) %>%
  distinct() %>%
  rename(peptide = reference_hit)

per_peptide_results <- per_peptide_results %>%
  left_join(category_lookup, by = "peptide")

# coefficient plot: one dot per peptide, coloured by complex, ordered by estimate
plot_per_peptide_coefficients <-
  per_peptide_results %>%
  filter(!is.na(estimate)) %>%
  ggplot(aes(x = estimate,
             y = fct_reorder(peptide, estimate),
             colour = category,
             alpha = p_adjusted < 0.10)) +
  geom_point() +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std_error,
                     xmax = estimate + 1.96 * std_error),
                 height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3),
                     labels = c("TRUE" = "FDR < 0.10", "FALSE" = "ns"),
                     name = NULL) +
  scale_colour_manual(
    values = c(
      "housekeeping_candidates_Zm.peptides" = "grey70",
      # auto-assign colours to every other category present
      setNames(
        scales::hue_pal()(length(setdiff(unique(per_peptide_results$category), "null"))),
        setdiff(sort(unique(per_peptide_results$category)), "null")
      )
    )
  ) +
  labs(x = "coefficient: effect of family size on bio8 (adjusted for genome duplication)",
       y = NULL, colour = "complex",
       title = "Per-peptide: family size vs temperature",
       subtitle = "tested only in species where peptide was detected; FDR < 0.10 highlighted") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6))
print(plot_per_peptide_coefficients)

# top 10 significant peptides, one scatterplot each
top_10_peptides <- per_peptide_results %>%
  filter(note == "ok", !is.na(p_adjusted)) %>%
  slice_min(p_adjusted, n = 10) %>%
  pull(peptide)

for (peptide_name in top_10_peptides) {
  plot_data <- analysis_table %>%
    select(species, bio8_median, fraction_duplicated,
           family_size = all_of(peptide_name)) %>%
    drop_na()
  
  peptide_info <- per_peptide_results %>% filter(peptide == peptide_name)
  subtitle_text <- sprintf("category: %s | coef = %.2f | p = %.3f | n = %d",
                           peptide_info$category, peptide_info$estimate,
                           peptide_info$p_value, peptide_info$n_species)
  
  p <- ggplot(plot_data, aes(x = family_size, y = bio8_median)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    geom_text(aes(label = species), size = 2.3,
              vjust = -0.6, check_overlap = TRUE) +
    labs(x = "family size",
         y = "bio8 median (mean temp of wettest quarter)",
         title = peptide_name,
         subtitle = subtitle_text) +
    theme_minimal()
  print(p)
}

# absolute effect size distribution: housekeeping vs photosynthetic peptides
per_peptide_results %>%
  filter(note == "ok", !is.na(estimate)) %>%
  mutate(peptide_class = if_else(category == "null",
                                 "housekeeping", "photosynthetic")) %>%
  ggplot(aes(x = abs(estimate), fill = peptide_class)) +
  geom_histogram(bins = 30, position = "identity", alpha = 0.6) +
  scale_fill_manual(values = c("housekeeping" = "grey70",
                               "photosynthetic" = "steelblue")) +
  labs(x = "absolute effect size (|coefficient|)",
       y = "number of peptides",
       fill = NULL,
       title = "Effect size distribution: housekeeping vs photosynthetic peptides") +
  theme_minimal()
