## =============================================================================
## Add a diet column to the species list
##
## diet: BIRDBASE Primary Diet (Sekercioglu et al. 2025) as published, with
##   "Invertebrate" split by AVONET Trophic.Niche (Tobias et al. 2022):
##     Invertivore      -> Terrestrial invertebrate
##     Aquatic predator -> Aquatic invertebrate
##   Invertebrate-eaters with any other AVONET niche (21 species, all AVONET
##   Omnivore) are assigned by hand from other sources (manual_invert below).
##   BIRDBASE Primary Diet: a food type scoring >= 6 of 10 (Invertebrate,
##   Fruit, Nectar, Seed, Vertebrate (land vertebrates), Fish, Scavenger,
##   Plant); otherwise Carnivore (mixed animal food), Herbivore (mixed plant
##   food) or Omnivore (animal and plant food).
##
## Reads:  data/spp_names_codes_group_aou_taxa_redlist.csv
##         data/bird_traits/  (source files, downloaded on first run)
##           BIRDBASE_v2025.1.xlsx                BIRDBASE
##           AVONET_Supplementary_dataset_1.xlsx  AVONET
## Outputs: data/spp_names_codes_group_aou_taxa_redlist_diet.csv
## =============================================================================

library(tidyverse)
library(readxl)
library(here)

here::i_am("0b_prepare_diet.R")

spp_csv    <- "data/spp_names_codes_group_aou_taxa_redlist.csv"
diet_csv   <- "data/spp_names_codes_group_aou_taxa_redlist_diet.csv"
trait_dir  <- "data/bird_traits"

if (!dir.exists(trait_dir)) dir.create(trait_dir, recursive = TRUE)

## -----------------------------------------------------------------------------
## Download source files (skipped for any file that already exists)
## -----------------------------------------------------------------------------
sources <- tribble(
  ~file,                                  ~url,
  "BIRDBASE_v2025.1.xlsx",                "https://ndownloader.figshare.com/files/55634729",
  "AVONET_Supplementary_dataset_1.xlsx",  "https://ndownloader.figshare.com/files/34480856"
)

options(timeout = max(600, getOption("timeout")))
for (i in seq_len(nrow(sources))) {
  dest <- file.path(trait_dir, sources$file[i])
  if (!file.exists(dest)) {
    cat("Downloading", sources$file[i], "\n")
    download.file(sources$url[i], dest, mode = "wb", quiet = TRUE)
  }
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

report_unmatched <- function(source, idx) {
  miss <- spp$Common.Name[is.na(idx)]
  cat("\n", source, ": matched ", sum(!is.na(idx)), " of ", nrow(spp), " species\n", sep = "")
  if (length(miss) > 0) cat(paste(" - no match:", miss), sep = "\n")
}

## -----------------------------------------------------------------------------
## BIRDBASE Primary Diet
## -----------------------------------------------------------------------------
## Two header rows: the second holds the column names
birdbase <- read_excel(file.path(trait_dir, "BIRDBASE_v2025.1.xlsx"), sheet = "Data", skip = 1,
                       col_types = "text", na = c("", "NA"), .name_repair = "minimal") %>%
  select(english  = `English Name (BirdLife > IOC > Clements>AviList)`,
         birdlife = `HBW/BirdLife International (v9.1)`,
         ioc      = `IOC World Bird List (v15.1)`,
         ebird    = `eBird/Clements (V2024b)`,
         primary_diet = `Primary Diet`)

## eBird/Clements 2024 names, then IOC, BirdLife (Red List) and English names,
## then this manual table of eBird names. Northern Goshawk: eBird has split it;
## the American species is "Astur atricapillus". Barn Owl and Dovekie have no
## scientific name in the species list.
manual_birdbase <- tribble(
  ~Common.Name,                    ~birdbase_name,
  "Northern Goshawk",              "Astur atricapillus",
  "Barn Owl",                      "Tyto alba",
  "Dovekie",                       "Alle alle"
)
birdbase_idx <- match_names(spp %>% transmute(sci_bbs), birdbase$ebird)
birdbase_idx <- match_names(spp %>% transmute(sci_bbs), birdbase$ioc, birdbase_idx)
birdbase_idx <- match_names(spp %>% transmute(redlistScientificName), birdbase$birdlife, birdbase_idx)
birdbase_idx <- match_names(spp %>% transmute(Common.Name, bbs_english), birdbase$english, birdbase_idx)
birdbase_idx <- match_names(spp %>% left_join(manual_birdbase, by = "Common.Name") %>% select(birdbase_name),
                            birdbase$ebird, birdbase_idx)
report_unmatched("BIRDBASE", birdbase_idx)

## -----------------------------------------------------------------------------
## AVONET Trophic.Niche (eBird taxonomy sheet, then BirdLife sheet)
## -----------------------------------------------------------------------------
avonet_file <- file.path(trait_dir, "AVONET_Supplementary_dataset_1.xlsx")
avonet <- bind_rows(
  read_excel(avonet_file, sheet = "AVONET2_eBird", na = "NA") %>% select(sci = Species2, Trophic.Niche),
  read_excel(avonet_file, sheet = "AVONET1_BirdLife", na = "NA") %>% select(sci = Species1, Trophic.Niche)
)

## eBird rows come first, so a name found in both sheets takes the eBird row
avonet_idx <- match_names(spp %>% transmute(sci_bbs, birdbase_ebird = birdbase$ebird[birdbase_idx],
                                            redlistScientificName),
                          avonet$sci)
report_unmatched("AVONET", avonet_idx)

## -----------------------------------------------------------------------------
## diet column
## -----------------------------------------------------------------------------
## BIRDBASE invertebrate-eaters that AVONET calls Omnivore, assigned from the
## All About Birds Food label (AAB), DeGraaf et al. 1985 breeding-season
## foraging substrate, EltonTraits 1.0 percent foraging at the water surface
## (water %) and AVONET habitat. Species whose prey changes with season are
## assigned by their breeding-season (BBS survey period) diet.
manual_invert <- tribble(
  ~Common.Name,                 ~invert,         # evidence
  ## all sources point to land: no water foraging, land substrates
  "American Kestrel",           "Terrestrial",   # AAB Small Animals; air hawker
  "Bendire's Thrasher",         "Terrestrial",   # AAB Insects; ground forager
  "Brewer's Sparrow",           "Terrestrial",   # AAB Insects; ground gleaner
  "Brown Thrasher",             "Terrestrial",   # AAB Omnivore; ground forager
  "Cattle Egret",               "Terrestrial",   # AAB Insects; ground gleaner in pastures
  "Gila Woodpecker",            "Terrestrial",   # AAB Omnivore; bark gleaner
  "Gray Kingbird",              "Terrestrial",   # AAB Insects; air sallier
  "Great-tailed Grackle",       "Terrestrial",   # AAB Insects; ground forager
  "Henslow's Sparrow",          "Terrestrial",   # AAB Insects; ground forager
  "Killdeer",                   "Terrestrial",   # AAB Insects; ground gleaner, water 0%
  "Red-breasted Nuthatch",      "Terrestrial",   # AAB Insects; bark gleaner
  "Sulphur-bellied Flycatcher", "Terrestrial",   # AAB Insects; air sallier
  "Western Bluebird",           "Terrestrial",   # AAB Insects; ground gleaner
  "White-breasted Nuthatch",    "Terrestrial",   # AAB Insects; bark gleaner
  "Yellow-billed Magpie",       "Terrestrial",   # AAB Omnivore; ground gleaner
  ## shorebirds, rails and gulls
  "Buff-breasted Sandpiper",    "Terrestrial",   # dry-grassland forager: ground gleaner, water 0%, AVONET Grassland (AAB: Aquatic invertebrates)
  "Pacific Golden-Plover",      "Terrestrial",   # ground gleaner, water 0%, AVONET Grassland (AAB: Aquatic invertebrates)
  "Long-billed Curlew",         "Terrestrial",   # breeds on grasslands eating grasshoppers/beetles; winters on mudflats (AAB: Aquatic invertebrates)
  "Whimbrel",                   "Aquatic",       # shoreline forager, AAB Aquatic invertebrates; tundra insects only while breeding
  "Yellow Rail",                "Aquatic",       # fresh-marsh forager, AAB Aquatic invertebrates, AVONET Wetland
  "Mew Gull",                   "Aquatic"        # shoreline gleaner, water 60%, AVONET Coastal
)

primary_diet <- birdbase$primary_diet[birdbase_idx]
avonet_niche <- avonet$Trophic.Niche[avonet_idx]
manual_split <- manual_invert$invert[match(spp$Common.Name, manual_invert$Common.Name)]

spp_diet <- spp %>%
  select(-sci_bbs) %>%
  mutate(diet = case_when(
    primary_diet != "Invertebrate"       ~ primary_diet,
    avonet_niche %in% "Invertivore"      ~ "Terrestrial invertebrate",
    avonet_niche %in% "Aquatic predator" ~ "Aquatic invertebrate",
    !is.na(manual_split)                 ~ paste(manual_split, "invertebrate"),
    TRUE                                 ~ primary_diet
  ))

cat("\ndiet:\n")
print(table(spp_diet$diet, useNA = "ifany"))

unsplit <- spp_diet %>% filter(diet == "Invertebrate")
if (nrow(unsplit) > 0) {
  cat("\nBIRDBASE invertebrate-eaters with no AVONET or manual split",
      "(left as \"Invertebrate\"):\n")
  cat(paste0(" - ", unsplit$Common.Name, " (AVONET: ",
             avonet_niche[match(unsplit$Common.Name, spp$Common.Name)], ")"), sep = "\n")
}

write.csv(spp_diet, diet_csv, row.names = FALSE)
cat("\nSaved:", diet_csv, "\n")
