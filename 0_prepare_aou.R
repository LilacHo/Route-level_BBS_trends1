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

## Check for repeated species codes in the input list
dup_check <- read.csv("data/spp_names_codes_group.csv", stringsAsFactors = FALSE)
dup_rows <- dup_check %>%
  filter(Code %in% Code[duplicated(Code)]) %>%
  arrange(Code)
if (nrow(dup_rows) > 0) {
  cat("Repeated codes in data/spp_names_codes_group.csv (",
      n_distinct(dup_rows$Code), " codes, ", nrow(dup_rows), " rows):\n", sep = "")
  print(dup_rows, row.names = FALSE)
  cat("\n")
} else {
  cat("No repeated codes in data/spp_names_codes_group.csv.\n")
}

## Load bbsBayes2 species list
strat_data <- load_bbs_data()
bbs_species <- strat_data$species %>%
  filter(unid_combined == TRUE)

## These species are listed under both "coastal" and "waterbirds", but only the
## waterbirds SDM rasters exist (data/rcp*_coastal/<code>/ folders are empty),
## so the coastal duplicate rows are dropped.
coastal_dup_codes <- c("BLTU", "LIGU", "MEGU", "ROSA", "SEPL", "SUSC", "WILL")
drop_coastal_dups <- function(df) {
  filter(df, !(Code %in% coastal_dup_codes & Group == "coastal"))
}

if (file.exists(aou_csv)) {
  cat("Found", aou_csv, "- skipping AOU matching.\n")
  spp_df <- read.csv(aou_csv, stringsAsFactors = FALSE)
  n_before <- nrow(spp_df)
  spp_df <- drop_coastal_dups(spp_df)
  if (nrow(spp_df) < n_before) {
    cat("Dropped", n_before - nrow(spp_df), "coastal duplicate rows from", aou_csv, "\n")
    tryCatch(
      write.csv(spp_df, aou_csv, row.names = FALSE),
      error = function(e) warning(aou_csv, " could not be rewritten (open in another ",
                                  "program?); it still has the coastal duplicates.")
    )
  }
} else {
  library(wildlifeR)

  ## Load species list
  spp_df <- read.csv("data/spp_names_codes_group.csv", stringsAsFactors = FALSE) %>%
    drop_coastal_dups()

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

## Add IUCN Red List status (redlistCategory, redlistPopulationTrend)
## Match "genus species" to scientificName in the Red List assessments download
redlist_dir <- "data/redlist_species_data_20260921"
redlist_csv <- "data/spp_names_codes_group_aou_taxa_redlist.csv"

assessments <- read.csv(file.path(redlist_dir, "assessments.csv"), stringsAsFactors = FALSE) %>%
  select(scientificName, redlistCategory, redlistPopulationTrend = populationTrend) %>%
  mutate(redlistScientificName = scientificName) %>%
  distinct(scientificName, .keep_all = TRUE)

## Manual lookup for species whose bbsBayes2 genus/species is not the name IUCN
## uses (newer AOS names, splits/lumps, subspecies). genus/species columns are
## left untouched; only the Red List columns are filled from these names.
manual_redlist <- tribble(
  ~Common.Name,                    ~redlist_name,
  "American Three-toed Woodpecker", "Picoides tridactylus",
  "American Tree Sparrow",         "Passerella arborea",
  "Arizona Woodpecker",            "Leuconotopicus arizonae",
  "Baird's Sparrow",               "Passerculus bairdii",
  "Black Oystercatcher",           "Haematopus ater",
  "Black-headed Gull",             "Larus ridibundus",
  "Black-necked Stilt",            "Himantopus himantopus",
  "Bonaparte's Gull",              "Larus philadelphia",
  "Bullock's Oriole",              "Icterus bullockiorum",
  "Cattle Egret",                  "Bubulcus ibis",
  "Common Redpoll",                "Acanthis flammea",
  "Cordilleran Flycatcher",        "Empidonax occidentalis",
  "Evening Grosbeak",              "Hesperiphona vespertina",
  "Franklin's Gull",               "Larus pipixcan",
  "Green Heron",                   "Butorides striata",
  "Hairy Woodpecker",              "Leuconotopicus villosus",
  "Henslow's Sparrow",             "Passerculus henslowii",
  "Hepatic Tanager",               "Piranga hepatica",
  "Herring Gull",                  "Larus smithsonianus",
  "Hoary Redpoll",                 "Acanthis flammea",
  "Laughing Gull",                 "Larus atricilla",
  "Least Bittern",                 "Ixobrychus exilis",
  "Mew Gull",                      "Larus canus",
  "Mountain Plover",               "Charadrius montanus",
  "Northern Goshawk",              "Accipiter gentilis",
  "Pacific-slope Flycatcher",      "Empidonax difficilis", 
  "Pileated Woodpecker",           "Hylatomus pileatus",
  "Red-cockaded Woodpecker",       "Leuconotopicus borealis",
  "Red-faced Cormorant",           "Urile urile",
  "Sandhill Crane",                "Grus canadensis",
  "Snowy Plover",                  "Charadrius nivosus",
  "Spotted Dove",                  "Spilopelia chinensis",
  "White-headed Woodpecker",       "Leuconotopicus albolarvatus",
  "Wilson's Phalarope",            "Steganopus tricolor",
  "Wilson's Plover",               "Charadrius wilsonia",
  "Winter Wren",                   "Troglodytes hiemalis",
  "Woodhouse's Scrub-Jay",         "Aphelocoma californica",
  "Yellow-billed Magpie",          "Pica nutalli"
)

spp_redlist <- spp_taxa %>%
  mutate(scientificName = ifelse(is.na(genus) | is.na(species), NA_character_,
                                 paste(genus, species))) %>%
  left_join(assessments, by = "scientificName") %>%
  left_join(manual_redlist, by = "Common.Name") %>%
  left_join(assessments %>% rename(redlistCategory_man = redlistCategory,
                                   redlistPopulationTrend_man = redlistPopulationTrend,
                                   redlistScientificName_man = redlistScientificName),
            by = c("redlist_name" = "scientificName")) %>%
  mutate(
    redlistCategory        = coalesce(redlistCategory, redlistCategory_man),
    redlistPopulationTrend = coalesce(redlistPopulationTrend, redlistPopulationTrend_man),
    redlistScientificName  = coalesce(redlistScientificName, redlistScientificName_man)
  ) %>%
  select(-scientificName, -redlist_name, -redlistCategory_man,
         -redlistPopulationTrend_man, -redlistScientificName_man) %>%
  relocate(redlistScientificName, .before = redlistCategory)

no_redlist <- spp_redlist %>% filter(is.na(redlistCategory)) %>% distinct(Common.Name, genus, species)
cat("\nSpecies matched to Red List assessments:",
    sum(!is.na(spp_redlist$redlistCategory)), "of", nrow(spp_redlist), "\n")
cat("Species with no Red List match:\n")
cat(paste(" -", no_redlist$Common.Name, "(", no_redlist$genus, no_redlist$species, ")"), sep = "\n")

write.csv(spp_redlist, redlist_csv, row.names = FALSE)
cat("\nSaved:", redlist_csv, "\n")
