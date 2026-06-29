#!/usr/bin/env Rscript

# Incremental local GBIF batch -> cleaned occurrences -> median BIO8
# No GBIF API calls. This reads your existing GBIF occurrence download.
#
# Examples:
#   Rscript species_bio8_from_batch.R \
#     --species "Quercus rubra" \
#     --gbif-file data/0004558-260409193756587.csv \
#     --bio8-tif data/tifs/bio08_Mean_Temperature_Wettest_Quarter.tif \
#     --out data/species_bio8_medians.csv
#
#   Rscript species_bio8_from_batch.R \
#     --species-file data/species.txt \
#     --gbif-file data/0004558-260409193756587.csv \
#     --bio8-tif data/tifs/bio08_Mean_Temperature_Wettest_Quarter.tif \
#     --out data/species_bio8_medians.csv
#
# Required R packages: data.table, CoordinateCleaner, terra

suppressPackageStartupMessages({
  library(data.table)
  library(CoordinateCleaner)
  library(terra)
})

parse_args <- function() {
  x <- commandArgs(trailingOnly = TRUE)
  args <- list(
    cache_dir = "data/species_bio8_cache",
    out = "data/species_bio8_medians.csv",
    force = FALSE,
    scale = 1,
    min_occurrences = 3L,
    exclude_publisher_key = NA_character_
  )
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (key == "--force") {
      args$force <- TRUE
      i <- i + 1L
    } else if (startsWith(key, "--")) {
      nm <- gsub("-", "_", sub("^--", "", key))
      if (i == length(x)) stop("Missing value for ", key)
      args[[nm]] <- x[[i + 1L]]
      i <- i + 2L
    } else {
      stop("Unexpected argument: ", key)
    }
  }
  args$scale <- as.numeric(args$scale)
  args$min_occurrences <- as.integer(args$min_occurrences)
  args
}

read_species <- function(args) {
  spp <- character()
  if (!is.null(args$species)) spp <- c(spp, trimws(unlist(strsplit(args$species, ","))))
  if (!is.null(args$species_file)) spp <- c(spp, trimws(readLines(args$species_file, warn = FALSE)))
  spp <- spp[nzchar(spp) & !startsWith(spp, "#")]
  spp <- unique(spp)
  if (!length(spp)) stop("Provide --species or --species-file")
  spp
}

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
raw_cache_file <- function(sp, cache_dir) file.path(cache_dir, "raw_occurrences", paste0(safe_name(sp), ".rds"))
clean_cache_file <- function(sp, cache_dir) file.path(cache_dir, "clean_occurrences", paste0(safe_name(sp), ".rds"))
bio8_cache_file <- function(sp, cache_dir) file.path(cache_dir, "bio8_occurrences", paste0(safe_name(sp), ".rds"))

filter_gbif_batch <- function(gbif_file, species_needed, cache_dir, exclude_publisher_key = NA_character_) {
  dir.create(file.path(cache_dir, "raw_occurrences"), recursive = TRUE, showWarnings = FALSE)

  species_file <- tempfile("species_", fileext = ".txt")
  fwrite(data.table(species = species_needed), species_file, col.names = FALSE)

  awk_file <- tempfile("filter_gbif_", fileext = ".awk")
  awk_code <- '
BEGIN {
  FS=OFS="\t"
  while ((getline s < species_file) > 0) {
    gsub(/^[ \t]+|[ \t]+$/, "", s)
    if (s != "") wanted[s] = 1
  }
}
NR == 1 {
  for (i=1; i<=NF; i++) {
    if ($i == "species") species_col = i
    if ($i == "scientificName") sciname_col = i
    if ($i == "decimalLatitude") lat_col = i
    if ($i == "decimalLongitude") lon_col = i
    if ($i == "publishingOrgKey") pub_col = i
  }
  if (!species_col || !sciname_col || !lat_col || !lon_col) {
    print "ERROR: required GBIF columns missing" > "/dev/stderr"
    exit 2
  }
  print "query_species", "species", "scientificName", "decimalLatitude", "decimalLongitude"
  next
}
{
  if (exclude_pub != "" && pub_col && $pub_col == exclude_pub) next
  for (sp in wanted) {
    # GBIF species is usually the clean binomial. scientificName may include authorship,
    # so accept exact match or binomial followed by whitespace.
    if ($species_col == sp || $sciname_col == sp || index($sciname_col, sp " ") == 1) {
      print sp, $species_col, $sciname_col, $lat_col, $lon_col
      break
    }
  }
}
'
  writeLines(awk_code, awk_file)

  exclude <- ifelse(is.na(exclude_publisher_key), "", exclude_publisher_key)
  cmd <- sprintf(
    "awk -v species_file=%s -v exclude_pub=%s -f %s %s",
    shQuote(species_file), shQuote(exclude), shQuote(awk_file), shQuote(gbif_file)
  )

  message("Scanning local GBIF batch for ", length(species_needed), " uncached species...")
  dt <- fread(cmd = cmd, sep = "\t", quote = "", showProgress = TRUE)

  for (sp in species_needed) {
    one <- dt[query_species == sp]
    saveRDS(one, raw_cache_file(sp, cache_dir))
    message("  cached raw occurrences: ", sp, " = ", nrow(one))
  }
  invisible(TRUE)
}

clean_occurrences <- function(sp, cache_dir, force = FALSE) {
  cache_file <- clean_cache_file(sp, cache_dir)
  if (file.exists(cache_file) && !force) return(readRDS(cache_file))

  raw_file <- raw_cache_file(sp, cache_dir)
  if (!file.exists(raw_file)) stop("Missing raw cache for ", sp)
  dt <- as.data.table(readRDS(raw_file))

  if (!nrow(dt)) {
    out <- data.table(species = character(), decimalLongitude = numeric(), decimalLatitude = numeric())
    saveRDS(out, cache_file)
    return(out)
  }

  out <- dt[, .(
    species = sp,
    decimalLongitude = as.numeric(decimalLongitude),
    decimalLatitude = as.numeric(decimalLatitude)
  )]

  out <- out[!is.na(decimalLatitude) & !is.na(decimalLongitude)]
  out <- out[decimalLatitude >= -90 & decimalLatitude <= 90 & decimalLongitude >= -180 & decimalLongitude <= 180]
  out <- out[!(decimalLatitude == 0 & decimalLongitude == 0)]
  out <- unique(out, by = c("species", "decimalLongitude", "decimalLatitude"))

  if (nrow(out)) {
    out <- out[cc_val(out, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged", verbose = FALSE)]
  }
  if (nrow(out)) {
    out <- out[cc_sea(out, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged", verbose = FALSE)]
  }
  if (nrow(out) >= 10) {
    out <- tryCatch(
      as.data.table(cc_outl(
        out,
        lon = "decimalLongitude",
        lat = "decimalLatitude",
        species = "species",
        method = "mad",
        mltpl = 10,
        value = "clean",
        verbose = FALSE
      )),
      error = function(e) {
        warning("Outlier filtering failed for ", sp, ": ", conditionMessage(e))
        out
      }
    )
  }

  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, cache_file)
  out
}

extract_bio8 <- function(sp, clean_dt, bio8, cache_dir, scale = 1, force = FALSE) {
  cache_file <- bio8_cache_file(sp, cache_dir)
  if (file.exists(cache_file) && !force) return(readRDS(cache_file))

  if (!nrow(clean_dt)) {
    out <- data.table(species = character(), decimalLongitude = numeric(), decimalLatitude = numeric(), bio8 = numeric())
    saveRDS(out, cache_file)
    return(out)
  }

  pts <- vect(clean_dt, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
  vals <- terra::extract(bio8, pts)[[2]] * scale
  out <- copy(clean_dt)[, bio8 := vals]
  out <- out[!is.na(bio8) & is.finite(bio8)]

  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, cache_file)
  out
}

args <- parse_args()
spp <- read_species(args)
if (is.null(args$gbif_file)) stop("Provide --gbif-file pointing to your local GBIF batch occurrence file")
if (is.null(args$bio8_tif)) stop("Provide --bio8-tif pointing to the BIO8 raster .tif")
if (!file.exists(args$gbif_file)) stop("GBIF file not found: ", args$gbif_file)
if (!file.exists(args$bio8_tif)) stop("BIO8 tif not found: ", args$bio8_tif)

dir.create(args$cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$cache_dir, "raw_occurrences"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$cache_dir, "clean_occurrences"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$cache_dir, "bio8_occurrences"), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args$out), recursive = TRUE, showWarnings = FALSE)

need_raw <- spp[args$force | !file.exists(vapply(spp, raw_cache_file, character(1), cache_dir = args$cache_dir))]
if (length(need_raw)) {
  filter_gbif_batch(args$gbif_file, need_raw, args$cache_dir, args$exclude_publisher_key)
} else {
  message("All requested species already have raw occurrence caches.")
}

message("Loading BIO8 raster: ", args$bio8_tif)
bio8 <- rast(args$bio8_tif)
if (nlyr(bio8) != 1) bio8 <- bio8[[1]]

results <- vector("list", length(spp))
for (i in seq_along(spp)) {
  sp <- spp[[i]]
  message("\n=== ", i, "/", length(spp), ": ", sp, " ===")
  raw <- as.data.table(readRDS(raw_cache_file(sp, args$cache_dir)))
  clean <- clean_occurrences(sp, args$cache_dir, args$force)
  bio <- extract_bio8(sp, clean, bio8, args$cache_dir, args$scale, args$force)

  results[[i]] <- data.table(
    species = sp,
    n_gbif_raw = nrow(raw),
    n_clean_occurrences = nrow(clean),
    n_bio8_occurrences = nrow(bio),
    bio8_median = if (nrow(bio) >= args$min_occurrences) median(bio$bio8, na.rm = TRUE) else NA_real_
  )
}

summary <- rbindlist(results, fill = TRUE)
fwrite(summary, args$out)
saveRDS(summary, sub("\\.csv$", ".rds", args$out))
message("\nSaved: ", args$out)
print(summary)
