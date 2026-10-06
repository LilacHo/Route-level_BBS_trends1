## =============================================================================
## Add an invasive column to the species list
##
## invasive: status in the conterminous 48 United States from the United States
##   Register of Introduced and Invasive Species, US-RIIS ver. 2.0 (Simpson et
##   al. 2022, https://doi.org/10.5066/P9KFFTOD), lower-48 (L48) list,
##   degreeOfEstablishment:
##     "established (category C3)"        -> Established
##     "invasive (category D2)"           -> Invasive
##     "widespread invasive (category E)" -> Widespread invasive
##   Species not on the L48 list -> Not listed. US-RIIS lists only species that
##   are non-native everywhere in the lower 48 and reproducing there, so species
##   native to part of the lower 48 but introduced elsewhere in it (e.g. House
##   Finch in the East) are Not listed.
##
## Reads:  data/spp_names_codes_group_aou_taxa_redlist_diet_migration.csv  (from 0c_prepare_migration.R)
##         data/bird_traits/United States Register of Introduced and Invasive
##           Species (US-RIIS) (ver. 2.0, November 2022)/USRIISv2csvFormat/USRIISv2_MasterList.csv
## Outputs: data/spp_names_codes_group_aou_taxa_redlist_diet_migration_invasive.csv
## =============================================================================

library(tidyverse)
library(here)

here::i_am("0d_prepare_invasive.R")

spp_csv      <- "data/spp_names_codes_group_aou_taxa_redlist_diet_migration.csv"
invasive_csv <- "data/spp_names_codes_group_aou_taxa_redlist_diet_migration_invasive.csv"
riis_csv     <- file.path("data/bird_traits",
                          "United States Register of Introduced and Invasive Species (US-RIIS) (ver. 2.0, November 2022)",
                          "USRIISv2csvFormat", "USRIISv2_MasterList.csv")

if (!file.exists(riis_csv)) {
  stop(riis_csv, " not found: download US-RIIS ver. 2.0 from https://doi.org/10.5066/P9KFFTOD")
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
## US-RIIS lower-48 birds
## -----------------------------------------------------------------------------
## scientificName includes the authority ("Passer domesticus (Linnaeus, 1758)"),
## so only the binomial is kept for matching
riis <- read_csv(riis_csv, show_col_types = FALSE) %>%
  filter(class == "Aves", locality == "L48") %>%
  transmute(sci = word(scientificName, 1, 2), vernacularName, degreeOfEstablishment)
cat("US-RIIS lower-48 bird records:", nrow(riis), "\n")

riis_idx <- match_names(spp %>% transmute(sci_bbs, redlistScientificName), riis$sci)
riis_idx <- match_names(spp %>% transmute(Common.Name, bbs_english), riis$vernacularName, riis_idx)

riis_status <- c(
  "established (category C3)"        = "Established",
  "invasive (category D2)"           = "Invasive",
  "widespread invasive (category E)" = "Widespread invasive"
)
unmapped <- setdiff(riis$degreeOfEstablishment, names(riis_status))
if (length(unmapped) > 0) stop("Unexpected US-RIIS degreeOfEstablishment: ", paste(unmapped, collapse = ", "))

spp_invasive <- spp %>%
  select(-sci_bbs) %>%
  mutate(invasive = coalesce(unname(riis_status[riis$degreeOfEstablishment[riis_idx]]), "Not listed"))

cat("\ninvasive:\n")
print(table(spp_invasive$invasive, useNA = "ifany"))
listed <- spp_invasive %>% filter(invasive != "Not listed") %>% arrange(invasive, Common.Name)
cat(paste0(" - ", listed$Common.Name, ": ", listed$invasive), sep = "\n")

write.csv(spp_invasive, invasive_csv, row.names = FALSE)
cat("\nSaved:", invasive_csv, "\n")
