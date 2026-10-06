## =============================================================================
## Add a migration column to the species list
##
## migration: Partners in Flight Avian Conservation Assessment Database (ACAD)
##   Global, version 2024.05.23, "Mig Status":
##     M -> Migratory (migrant or partial migrant)
##     R -> Resident
##   ACAD follows eBird Status & Trends, which counts partial migrants as
##   migrants.
##
## Reads:  data/spp_names_codes_group_aou_taxa_redlist_diet.csv  (from 0b_prepare_diet.R)
##         data/bird_traits/ACAD Global 2024.05.23.xlsx
##           exported by hand from the ACAD scores app
##           (https://pif.birdconservancy.org/avian-conservation-assessment-database-scores/);
##           the data service behind the app needs an access token, so the file
##           cannot be downloaded by script
## Outputs: data/spp_names_codes_group_aou_taxa_redlist_diet_migration.csv
## =============================================================================

library(tidyverse)
library(readxl)
library(here)

here::i_am("0c_prepare_migration.R")

spp_csv       <- "data/spp_names_codes_group_aou_taxa_redlist_diet.csv"
migration_csv <- "data/spp_names_codes_group_aou_taxa_redlist_diet_migration.csv"
acad_xlsx     <- "data/bird_traits/ACAD Global 2024.05.23.xlsx"

if (!file.exists(acad_xlsx)) {
  stop(acad_xlsx, " not found: export the ACAD Global table from ",
       "https://pif.birdconservancy.org/avian-conservation-assessment-database-scores/")
}

## -----------------------------------------------------------------------------
## Species list and name matching
## -----------------------------------------------------------------------------
spp <- read.csv(spp_csv, stringsAsFactors = FALSE) %>%
  mutate(sci_bbs = ifelse(is.na(genus) | is.na(species), NA_character_, paste(genus, species)))

## Lower-case, drop punctuation/spaces, and treat "Grey"/"Gray" alike, so that
## "Black-Hawk" = "Black Hawk" and "Harris' Hawk" = "Harris's Hawk"
norm_name <- function(x) {
  x <- tolower(x)
  x <- gsub("grey", "gray", x)
  x <- gsub("'s\\b", "", x)
  gsub("[^a-z]", "", x)
}

## Return the row in `src_keys` matching each species, trying each candidate
## name column of `cand` in order (first match wins). Species already matched
## in `idx` keep their match, so calls can be chained.
match_names <- function(cand, src_keys, idx = rep(NA_integer_, nrow(cand))) {
  keys <- norm_name(src_keys)
  for (col in names(cand)) {
    m <- match(norm_name(cand[[col]]), keys)
    idx <- ifelse(is.na(idx) & !is.na(cand[[col]]), m, idx)
  }
  idx
}

## -----------------------------------------------------------------------------
## ACAD Mig Status
## -----------------------------------------------------------------------------
## ACAD names follow the AOS Check-list (7th ed., 63rd supplement): matched on
## common name, then scientific name
acad <- read_excel(acad_xlsx, sheet = "Globals") %>%
  select(common = `Common Name`, sci = `Scientific Name`, mig_status = `Mig Status`)

acad_idx <- match_names(spp %>% transmute(Common.Name, bbs_english), acad$common)
acad_idx <- match_names(spp %>% transmute(sci_bbs, redlistScientificName), acad$sci, acad_idx)

miss <- spp$Common.Name[is.na(acad_idx)]
cat("\nACAD: matched", sum(!is.na(acad_idx)), "of", nrow(spp), "species\n")
if (length(miss) > 0) cat(paste(" - no match:", miss), sep = "\n")

spp_migration <- spp %>%
  select(-sci_bbs) %>%
  mutate(migration = recode(acad$mig_status[acad_idx], M = "Migratory", R = "Resident"))

cat("\nmigration:\n")
print(table(spp_migration$migration, useNA = "ifany"))

write.csv(spp_migration, migration_csv, row.names = FALSE)
cat("\nSaved:", migration_csv, "\n")
