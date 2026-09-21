## =============================================================================
## Add AOU numeric codes to species list using wildlifeR, then add bbsBayes2
## taxonomy (order, family, genus, species)
## Reads: data/spp_names_codes_group.csv
## Outputs: data/spp_names_codes_group_aou.csv       (skipped if it already exists)
##          data/spp_names_codes_group_aou_taxa.csv
## =============================================================================

library(tidyverse)
library(bbsBayes2)
library(here)

here::i_am("0_prepare_aou.R")

if (!dir.exists("data")) dir.create("data")

aou_csv  <- "data/spp_names_codes_group_aou.csv"
taxa_csv <- "data/spp_names_codes_group_aou_taxa.csv"

## Load bbsBayes2 species list
strat_data <- load_bbs_data()
bbs_species <- strat_data$species %>%
  filter(unid_combined == TRUE)

if (file.exists(aou_csv)) {
  cat("Found", aou_csv, "- skipping AOU matching.\n")
  spp_df <- read.csv(aou_csv, stringsAsFactors = FALSE)
} else {
  library(wildlifeR)

  ## Load species list
  spp_df <- read.csv("data/spp_names_codes_group.csv", stringsAsFactors = FALSE)

  ## Get AOU numeric codes from wildlifeR
  aou_codes <- wildlifeR::AOU_species_codes %>%
    select(spp.num, alpha.code, name)

  ## Join by 4-letter alpha code
  spp_df <- spp_df %>%
    left_join(aou_codes, by = c("Code" = "alpha.code"))

  ## Report any species without an AOU match
  no_aou <- spp_df %>% filter(is.na(spp.num)) %>% distinct(Common.Name, Code)
  if (nrow(no_aou) > 0) {
    cat("Species with no AOU number in wildlifeR:\n")
    cat(paste(" -", no_aou$Common.Name, "(", no_aou$Code, ")"), sep = "\n")
    cat("\n")
  } else {
    cat("All species matched to AOU numbers.\n")
  }

  ## Match to bbsBayes2 species list
  spp_df <- spp_df %>%
    mutate(
      in_bbs_aou  = spp.num %in% bbs_species$aou,
      in_bbs_name = Common.Name %in% bbs_species$english,
      in_bbs      = in_bbs_aou | in_bbs_name,
      bbs_english = case_when(
        in_bbs_aou  ~ bbs_species$english[match(spp.num, bbs_species$aou)],
        in_bbs_name ~ Common.Name,
        TRUE        ~ NA_character_
      )
    )

  ## Manually fill in AOU numbers for species missing from wildlifeR
  manual_aou <- tribble(
    ~Common.Name,            ~spp.num,
    "Dark-eyed Junco",       5677,
    "Northern Flicker",      4123,
    "Northwestern Crow",     4880,
    "Red-faced Cormorant",   1207,
    "Snow Goose",            1690,
    "Yellow-rumped Warbler",  6556
  )

  for (i in seq_len(nrow(manual_aou))) {
    idx <- which(spp_df$Common.Name == manual_aou$Common.Name[i])
    if (length(idx) > 0) {
      spp_df$spp.num[idx] <- manual_aou$spp.num[i]
      # Re-check against bbsBayes2
      if (manual_aou$spp.num[i] %in% bbs_species$aou) {
        spp_df$in_bbs_aou[idx] <- TRUE
        spp_df$in_bbs[idx] <- TRUE
        spp_df$bbs_english[idx] <- bbs_species$english[match(manual_aou$spp.num[i], bbs_species$aou)]
      }
      cat("  Manual AOU set:", manual_aou$Common.Name[i], "->", manual_aou$spp.num[i], "\n")
    }
  }

  cat("\nTotal species:", nrow(spp_df), "\n")
  cat("Matched by AOU number:", sum(spp_df$in_bbs_aou, na.rm = TRUE), "\n")
  cat("Matched by name only:", sum(!spp_df$in_bbs_aou & spp_df$in_bbs_name, na.rm = TRUE), "\n")
  cat("Not found in bbsBayes2:", sum(!spp_df$in_bbs), "\n")

  ## Save
  write.csv(spp_df, aou_csv, row.names = FALSE)
  cat("\nSaved:", aou_csv, "\n")
}

## Add taxonomy (order, family, genus, species) from bbsBayes2
bbs_taxonomy <- bbs_species %>%
  select(bbs_english = english, order, family, genus, species) %>%
  distinct(bbs_english, .keep_all = TRUE)

spp_taxa <- spp_df %>%
  select(-any_of(c("order", "family", "genus", "species"))) %>%
  left_join(bbs_taxonomy, by = "bbs_english")

cat("\nSpecies with bbsBayes2 taxonomy:", sum(!is.na(spp_taxa$genus)), "of", nrow(spp_taxa), "\n")

write.csv(spp_taxa, taxa_csv, row.names = FALSE)
cat("Saved:", taxa_csv, "\n")
