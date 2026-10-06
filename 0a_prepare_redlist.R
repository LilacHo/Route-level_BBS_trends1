## =============================================================================
## Add IUCN Red List status (redlistScientificName, redlistCategory,
## redlistPopulationTrend) to the species list
##
## Species are matched on "genus species" (bbsBayes2 taxonomy) to scientificName
## in a Red List assessments download, then through a manual lookup for names
## IUCN spells differently. Species with no match get NA.
##
## Reads:  data/spp_names_codes_group_aou_taxa.csv  (from 0_prepare_aou.R)
##         data/redlist_species_data_20260921/assessments.csv  (IUCN Red List download)
## Outputs: data/spp_names_codes_group_aou_taxa_redlist.csv
## =============================================================================

library(tidyverse)
library(here)

here::i_am("0a_prepare_redlist.R")

taxa_csv    <- "data/spp_names_codes_group_aou_taxa.csv"
redlist_dir <- "data/redlist_species_data_20260921"
redlist_csv <- "data/spp_names_codes_group_aou_taxa_redlist.csv"

spp_taxa <- read.csv(taxa_csv, stringsAsFactors = FALSE)

assessments <- read.csv(file.path(redlist_dir, "assessments.csv"), stringsAsFactors = FALSE) %>%
  select(scientificName, redlistCategory, redlistPopulationTrend = populationTrend) %>%
  mutate(redlistScientificName = scientificName) %>%
  distinct(scientificName, .keep_all = TRUE)

## Manual lookup for species whose bbsBayes2 genus/species is not the name IUCN
## uses (newer AOS names, splits/lumps, subspecies). genus/species columns are
## left untouched; only the Red List columns are filled from these names.
## Black Oystercatcher: H. bachmani is not in the 2026-09-21 download, so it is
## NA until its assessment is added (H. ater is the South American Blackish
## Oystercatcher).
manual_redlist <- tribble(
  ~Common.Name,                    ~redlist_name,
  "American Three-toed Woodpecker", "Picoides tridactylus",
  "American Tree Sparrow",         "Passerella arborea",
  "Arizona Woodpecker",            "Leuconotopicus arizonae",
  "Baird's Sparrow",               "Passerculus bairdii",
  "Black Oystercatcher",           "Haematopus bachmani",
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
cat("Species matched to Red List assessments:",
    sum(!is.na(spp_redlist$redlistCategory)), "of", nrow(spp_redlist), "\n")
cat("Species with no Red List match:\n")
cat(paste(" -", no_redlist$Common.Name, "(", no_redlist$genus, no_redlist$species, ")"), sep = "\n")

write.csv(spp_redlist, redlist_csv, row.names = FALSE)
cat("\nSaved:", redlist_csv, "\n")
