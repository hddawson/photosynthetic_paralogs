#!/usr/bin/env Rscript

# Grep-first local GBIF batch -> cleaned occurrences -> median BIO8
# This does NOT download GBIF data and does NOT R-scan the full occurrence table.
# It shells out to grep/ripgrep/parallel-grep first, then R reads only the matched rows.
#
# Examples:
#   Rscript species_bio8_grep.R \
#     --species "Quercus rubra" \
#     --gbif-file data/0004558-260409193756587.csv \
#     --bio8-tif data/tifs/bio08_Mean_Temperature_Wettest_Quarter.tif \
#     --out data/species_bio8_medians.csv
#
#   Rscript species_bio8_grep.R \
#     --species-file data/species.txt \
#     --gbif-file data/0004558-260409193756587.csv \
#     --bio8-tif data/tifs/bio08_Mean_Temperature_Wettest_Quarter.tif \
#     --grep-engine parallel \
#     --threads 40
#
# Required R packages: data.table, CoordinateCleaner, terra
# Optional command-line tools: rg (ripgrep), GNU parallel

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
    min_occurrences = 1L,
    exclude_publisher_key = NA_character_,
    grep_engine = "grep",        # grep, rg, or parallel
    threads = 8L,
    block = "500M",
    sep = "\t"
  )
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (key == "--force") {
      args$force <- TRUE; i <- i + 1L
    } else if (startsWith(key, "--")) {
      nm <- gsub("-", "_", sub("^--", "", key))
      if (i == length(x)) stop("Missing value for ", key)
      args[[nm]] <- x[[i + 1L]]; i <- i + 2L
    } else stop("Unexpected argument: ", key)
  }
  args$scale <- as.numeric(args$scale)
  args$min_occurrences <- as.integer(args$min_occurrences)
  args$threads <- as.integer(args$threads)
  args
}

read_species <- function(args) {
  spp <- character()
  if (!is.null(args$species)) spp <- c(spp, trimws(unlist(strsplit(args$species, ","))))
  if (!is.null(args$species_file)) spp <- c(spp, trimws(readLines(args$species_file, warn = FALSE)))
  spp <- unique(spp[nzchar(spp) & !startsWith(spp, "#")])
  if (!length(spp)) stop("Provide --species or --species-file")
  spp
}

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
raw_cache_file <- function(sp, cache_dir) file.path(cache_dir, "raw_occurrences", paste0(safe_name(sp), ".rds"))
clean_cache_file <- function(sp, cache_dir) file.path(cache_dir, "clean_occurrences", paste0(safe_name(sp), ".rds"))
bio8_cache_file <- function(sp, cache_dir) file.path(cache_dir, "bio8_occurrences", paste0(safe_name(sp), ".rds"))

have_cmd <- function(x) nzchar(Sys.which(x))
q <- function(x) shQuote(x, type = "sh")

write_patterns <- function(species_needed) {
  # Grep broad candidates fast. R exact-filters afterward.
  # Include a tab-prefixed pattern for exact species column matches in TSV-like GBIF files,
  # plus plain binomial to catch scientificName with authorship.
  patterns <- unique(c(paste0("\t", species_needed, "\t"), species_needed))
  f <- tempfile("species_patterns_", fileext = ".txt")
  writeLines(patterns, f, useBytes = TRUE)
  f
}

grep_candidates_to_file <- function(gbif_file, species_needed, cache_dir, args) {
  dir.create(file.path(cache_dir, "grep_candidates"), recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(cache_dir, "grep_candidates", paste0("candidates_", digest_key(species_needed), ".tsv"))
  if (file.exists(out_file) && !args$force) return(out_file)

  pat_file <- write_patterns(species_needed)
  tmp_body <- tempfile("gbif_grep_body_", fileext = ".tsv")
  on.exit(unlink(c(pat_file, tmp_body)), add = TRUE)

  message("Grep-filtering local GBIF batch for ", length(species_needed), " uncached species...")
  message("Engine: ", args$grep_engine, " | file: ", gbif_file)

  # Header only; cheap.
  header_cmd <- sprintf("head -n 1 %s > %s", q(gbif_file), q(out_file))
  status <- system(header_cmd)
  if (status != 0) stop("Failed to read GBIF header with head")

  exclude <- !is.na(args$exclude_publisher_key) && nzchar(args$exclude_publisher_key)

  if (args$grep_engine == "rg") {
    if (!have_cmd("rg")) stop("--grep-engine rg requested, but rg is not on PATH")
    body_cmd <- sprintf("LC_ALL=C rg -F -f %s --no-heading --color never %s",
                        q(pat_file), q(gbif_file))
    # For a file path argument rg prints only matching lines by default; --no-heading is harmless.
  } else if (args$grep_engine == "parallel") {
    if (!have_cmd("parallel")) stop("--grep-engine parallel requested, but GNU parallel is not on PATH")
    body_cmd <- sprintf(
      "LC_ALL=C parallel --pipepart --block %s -j %d -a %s 'grep -F -f %s'",
      q(args$block), args$threads, q(gbif_file), q(pat_file)
    )
  } else if (args$grep_engine == "grep") {
    body_cmd <- sprintf("LC_ALL=C grep -F -f %s %s", q(pat_file), q(gbif_file))
  } else {
    stop("Unknown --grep-engine: ", args$grep_engine, "; use grep, rg, or parallel")
  }

  if (exclude) body_cmd <- sprintf("%s | LC_ALL=C grep -F -v %s", body_cmd, q(args$exclude_publisher_key))

  # grep exits 1 when no matches. That is not an error here.
  cmd <- sprintf("set +e; %s > %s; code=$?; if [ $code -gt 1 ]; then exit $code; else exit 0; fi",
                 body_cmd, q(tmp_body))
  status <- system(cmd)
  if (status != 0) stop("grep candidate filtering failed")

  # Append body after header. If no matches, output is header-only.
  status <- system(sprintf("cat %s >> %s", q(tmp_body), q(out_file)))
  if (status != 0) stop("Failed to append grep candidates")

  message("Candidate file: ", out_file)
  out_file
}

digest_key <- function(x) {
  # Avoid requiring digest package.
  y <- paste(sort(x), collapse = "|")
  sprintf("%08x", sum(utf8ToInt(y) * seq_along(utf8ToInt(y))) %% .Machine$integer.max)
}

cache_raw_from_candidates <- function(candidate_file, species_needed, cache_dir, sep = "\t") {
  dir.create(file.path(cache_dir, "raw_occurrences"), recursive = TRUE, showWarnings = FALSE)

  cols <- names(fread(candidate_file, nrows = 0, sep = sep, quote = "", showProgress = FALSE))
  required <- c("species", "scientificName", "decimalLatitude", "decimalLongitude")
  missing <- setdiff(required, cols)
  if (length(missing)) stop("GBIF candidate file missing columns: ", paste(missing, collapse = ", "))

  keep <- required
  if ("publishingOrgKey" %in% cols) keep <- c(keep, "publishingOrgKey")

  dt <- fread(candidate_file, sep = sep, quote = "", select = keep, showProgress = TRUE)
  if (!nrow(dt)) {
    for (sp in species_needed) saveRDS(data.table(), raw_cache_file(sp, cache_dir))
    return(invisible(TRUE))
  }

  # Exact-ish filter after broad grep. scientificName may contain authorship.
  for (sp in species_needed) {
    one <- dt[
      species == sp |
        scientificName == sp |
        startsWith(scientificName, paste0(sp, " "))
    ]
    saveRDS(one, raw_cache_file(sp, cache_dir))
    message("  cached raw occurrences: ", sp, " = ", nrow(one))
  }
  invisible(TRUE)
}

clean_occurrences <- function(sp, cache_dir, force = FALSE) {
  cache_file <- clean_cache_file(sp, cache_dir)
  if (file.exists(cache_file) && !force) return(readRDS(cache_file))

  dt <- as.data.table(readRDS(raw_cache_file(sp, cache_dir)))
  if (!nrow(dt)) {
    out <- data.table(species = character(), decimalLongitude = numeric(), decimalLatitude = numeric())
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
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

  if (nrow(out)) out <- out[cc_val(out, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged", verbose = FALSE)]
  if (nrow(out)) out <- out[cc_sea(out, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged", verbose = FALSE)]
  if (nrow(out) >= 10) {
    out <- tryCatch(
      as.data.table(cc_outl(out, lon = "decimalLongitude", lat = "decimalLatitude", species = "species",
                            method = "mad", mltpl = 10, value = "clean", verbose = FALSE)),
      error = function(e) { warning("Outlier filtering failed for ", sp, ": ", conditionMessage(e)); out }
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
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
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
dir.create(dirname(args$out), recursive = TRUE, showWarnings = FALSE)

need_raw <- spp[args$force | !file.exists(vapply(spp, raw_cache_file, character(1), cache_dir = args$cache_dir))]
if (length(need_raw)) {
  candidate_file <- grep_candidates_to_file(args$gbif_file, need_raw, args$cache_dir, args)
  cache_raw_from_candidates(candidate_file, need_raw, args$cache_dir, args$sep)
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
