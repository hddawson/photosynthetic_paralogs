#!/usr/bin/env Rscript

library(ape)
library(arrow)
library(dplyr)
library(stringr)

# ---------------- paths ---------------- #

GRAFTED_TREE_FILE <- "data/grafted_family_tree_raxml.tre"
TIP_METADATA_FILE <- "data/processed_data.parquet"
GENOME_METADATA_FILE <- "data/hifi_angiosperm_metadata/hifi_representative_newest_per_species.csv"

OUT_TREE <- "results/grafted_family_tree_hifi_subset.tre"
OUT_MATCH_TABLE <- "results/grafted_family_tree_hifi_subset_matches.csv"

# ---------------- helpers ---------------- #

clean_binomial <- function(x) {
  x %>%
    str_squish() %>%
    str_split_fixed("\\s+", 3) %>%
    .[, 1:2] %>%
    apply(1, paste, collapse = " ") %>%
    str_to_lower()
}

safe_label <- function(accession, binomial) {
  binomial_clean <- binomial %>%
    str_squish() %>%
    str_replace_all("\\s+", "_")

  paste(accession, binomial_clean, sep = "|")
}

# ---------------- read files ---------------- #

tr <- read.tree(GRAFTED_TREE_FILE)

tip_df <- read_parquet(TIP_METADATA_FILE) %>%
  mutate(
    FileBasename = as.character(FileBasename),
    Organism = as.character(Organism),
    species_binomial = clean_binomial(Organism)
  ) %>%
  select(FileBasename, Organism, species_binomial)

genome_df <- read.csv(GENOME_METADATA_FILE, check.names = FALSE) %>%
  transmute(
    genome_accession = `Assembly Accession`,
    genome_organism = `Organism Name`,
    species_binomial = clean_binomial(`Organism Name`)
  )

# ---------------- match tree tips to species ---------------- #

tree_tip_df <- tibble(FileBasename = tr$tip.label)

tip_matches <- tree_tip_df %>%
  left_join(tip_df, by = "FileBasename")

# ---------------- merge plastid tips to genome species ---------------- #

matched <- tip_matches %>%
  inner_join(genome_df, by = "species_binomial") %>%
  distinct(FileBasename, .keep_all = TRUE) %>%
  mutate(
    new_tip_label = safe_label(genome_accession, species_binomial)
  )

message("Original tree tips: ", length(tr$tip.label))
message("Tree tips with organism metadata: ", sum(!is.na(tip_matches$Organism)))
message("Genome species: ", n_distinct(genome_df$species_binomial))
message("Matched tree tips to genome species: ", nrow(matched))

stopifnot(nrow(matched) > 0)

# ---------------- subset tree ---------------- #

tips_to_drop <- setdiff(tr$tip.label, matched$FileBasename)
tr_subset <- drop.tip(tr, tips_to_drop)

# Relabel tips in tree order
label_map <- matched$new_tip_label
names(label_map) <- matched$FileBasename

tr_subset$tip.label <- label_map[tr_subset$tip.label]

stopifnot(!any(is.na(tr_subset$tip.label)))
stopifnot(length(tr_subset$tip.label) == nrow(matched))

# ---------------- save outputs ---------------- #

write.tree(tr_subset, file = OUT_TREE)
write.csv(matched, OUT_MATCH_TABLE, row.names = FALSE)

message("Saved subsetted tree to: ", OUT_TREE)
message("Saved match table to: ", OUT_MATCH_TABLE)
