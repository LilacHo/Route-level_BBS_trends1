## =============================================================================
## Add AOU numeric codes to species list using wildlifeR
## Reads: data/spp_names_codes_group.csv
## Outputs: data/spp_names_codes_group_aou.csv
## =============================================================================

library(tidyverse)
library(wildlifeR)
library(bbsBayes2)
library(here)

here::i_am("0_prepare_aou.R")

if (!dir.exists("data")) dir.create("data")

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

## Load bbsBayes2 species list and match
strat_data <- load_bbs_data()
bbs_species <- strat_data$species %>%
  filter(unid_combined == TRUE)

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

cat("\nTotal species:", nrow(spp_df), "\n")
cat("Matched by AOU number:", sum(spp_df$in_bbs_aou, na.rm = TRUE), "\n")
cat("Matched by name only:", sum(!spp_df$in_bbs_aou & spp_df$in_bbs_name, na.rm = TRUE), "\n")
cat("Not found in bbsBayes2:", sum(!spp_df$in_bbs), "\n")

## Save
write.csv(spp_df, "data/spp_names_codes_group_aou.csv", row.names = FALSE)
cat("\nSaved: data/spp_names_codes_group_aou.csv\n")
