# 4e_habitat_category_sig_structure_classes.R
#
# Companion to 4d_habitat_category_trends.R. Where 4d compares mean/weighted
# trend BY category, this script classifies each SPECIES by the STRUCTURE of
# its Contraction/Stable/Expansion compact letter display (CLD) -- i.e. which
# of the three categories are/aren't statistically distinguishable from one
# another for that species -- since the letters themselves aren't ordinal
# (see the project record: letter identity is assigned per-species by the
# underlying pairwise test, not by a fixed rank/direction, so "a" doesn't
# mean "highest" or "Contraction" consistently across species).
#
# Each species' (Contraction.sig, Stable.sig, Expansion.sig) triple is
# classified into one of 11 structure classes by counting, for each pair of
# categories, whether their letter SETS overlap (share >=1 letter, i.e. NOT
# significantly different) -- generalized via set overlap rather than exact
# string match so it isn't hardcoded to the "ab"-only compound letters seen
# in the current data:
#
#   0 NA, 3 pairwise overlaps (all mutually non-different)      -> No difference
#   0 NA, 0 pairwise overlaps (all three mutually different)    -> All three differ
#   0 NA, 1 pairwise overlap  (two tied, third differs)         -> C~S!=E / S~E!=C / C~E!=S
#   0 NA, 2 pairwise overlaps (one category "bridges" the other -> Bridging-Contraction /
#                              two, which are themselves distinct)  Bridging-Stable / Bridging-Expansion
#   1 NA, remaining pair overlaps                                -> No difference
#   1 NA, remaining pair doesn't overlap                         -> Partial: <the differing pair>
#   >=2 NA                                                       -> Insufficient data
#
# C~E!=S ("Contraction and Expansion are themselves statistically
# indistinguishable, Stable is the odd one out") is flagged as
# hypothesis-relevant in the opposite direction from the others: it says the
# two categories the project's hypothesis contrasts aren't detectably
# different for that species. Bridging-Stable (a,ab,b-type) is the only
# bridging class that resolves a direct Contraction-vs-Expansion difference;
# the other two bridging classes resolve a different pair and are weaker
# evidence for this project's specific hypothesis. See the project record
# for the full reasoning and worked examples behind this scheme.
#
# Reads:  output/files/all_species_habitat_contraction_stable_expansion_2010_2025_base.xlsx
#         output/files/all_species_habitat_contraction_stable_expansion_2010_2025_anthro.xlsx
#         (see 4d_habitat_category_trends.R's header for where these come from)
#         data/spp_names_codes_group_aou.csv -- species -> bird Group.
#
# Writes: output/files/habitat_category_sig_structure_class_by_species_2010_2025.csv
#           One row per model x scenario x species: its Group, the three sig
#           letters, and its assigned structure_class.
#         output/files/habitat_category_sig_structure_class_totals_2010_2025.csv
#           One row per model x scenario x structure_class: n species.
#         output/files/habitat_category_sig_structure_class_by_group_2010_2025.csv
#           One row per model x scenario x Group x structure_class: n species.

library(here)
library(tidyverse)
library(readxl)

here::i_am("4e_habitat_category_sig_structure_classes.R")

# Settings -------------------------------------------------------------------
in_files <- c(
  base   = here::here("output", "files", "all_species_habitat_contraction_stable_expansion_2010_2025_base.xlsx"),
  anthro = here::here("output", "files", "all_species_habitat_contraction_stable_expansion_2010_2025_anthro.xlsx")
)
group_lookup_csv <- here::here("data", "spp_names_codes_group_aou.csv")

out_species_csv <- here::here("output", "files", "habitat_category_sig_structure_class_by_species_2010_2025.csv")
out_totals_csv  <- here::here("output", "files", "habitat_category_sig_structure_class_totals_2010_2025.csv")
out_group_csv   <- here::here("output", "files", "habitat_category_sig_structure_class_by_group_2010_2025.csv")

missing_files <- in_files[!file.exists(in_files)]
if (length(missing_files) > 0) {
  stop("Missing input workbook(s): ", paste(missing_files, collapse = ", "),
       " -- see this script's header comment for where these come from.")
}

# Read + combine both workbooks (both RCP sheets each) -----------------------
read_one_workbook <- function(path, model_tag) {
  sheets <- excel_sheets(path)
  map_dfr(sheets, function(sh) {
    d <- read_excel(path, sheet = sh)
    d$scenario <- sh
    d$model <- model_tag
    d
  })
}

all_raw <- map2_dfr(in_files, names(in_files), read_one_workbook)
names(all_raw) <- make.names(names(all_raw))  # "Contraction sig" -> "Contraction.sig" etc.

# Bird group lookup, de-duplicated the same way every other script in this
# project resolves species dual-listed under two groups (keep the
# first-listed group only -- same convention as 4d/1c/2c/gamma_lookup.R).
spp_group <- read.csv(group_lookup_csv, stringsAsFactors = FALSE) %>%
  distinct(Common.Name, Code, .keep_all = TRUE) %>%
  select(Species = Common.Name, Group)

sig_long <- all_raw %>%
  select(Species, model, scenario,
         Contraction.sig, Stable.sig, Expansion.sig) %>%
  left_join(spp_group, by = "Species")

#' Classify one species' (Contraction, Stable, Expansion) CLD letters into
#' one of the 11 structure classes described in this script's header.
#'
#' @param c,s,e single-character-or-compound CLD letter strings (e.g. "a",
#'   "ab"), or NA if that category had no valid test for this species.
#' @return a single structure_class string.
classify_structure <- function(c, s, e) {
  letters_of <- function(x) if (is.na(x)) NA else strsplit(x, "")[[1]]
  shares <- function(x, y) {
    lx <- letters_of(x); ly <- letters_of(y)
    if (any(is.na(lx)) || any(is.na(ly))) NA else length(intersect(lx, ly)) > 0
  }

  n_na <- sum(is.na(c), is.na(s), is.na(e))
  if (n_na >= 2) return("Insufficient data")

  if (n_na == 1) {
    if (is.na(c)) return(if (shares(s, e)) "No difference" else "Partial: Stable!=Expansion")
    if (is.na(e)) return(if (shares(c, s)) "No difference" else "Partial: Contraction!=Stable")
    # is.na(s) -- not observed in current data, handled for completeness
    return(if (shares(c, e)) "No difference" else "Partial: Contraction!=Expansion")
  }

  # All three present.
  cs <- shares(c, s); se <- shares(s, e); ce <- shares(c, e)
  n_overlap <- sum(cs, se, ce)

  if (n_overlap == 3) return("No difference")
  if (n_overlap == 0) return("All three differ")
  if (n_overlap == 1) {
    if (cs) return("C~S!=E")
    if (se) return("S~E!=C")
    return("C~E!=S")
  }
  # n_overlap == 2 -- one category bridges the other two, which are
  # themselves distinct. The "hub" is whichever category is in BOTH
  # overlapping pairs.
  if (cs && se) return("Bridging-Stable")        # Stable overlaps both C and E
  if (cs && ce) return("Bridging-Contraction")   # Contraction overlaps both S and E
  return("Bridging-Expansion")                   # Expansion overlaps both C and S
}

sig_long <- sig_long %>%
  rowwise() %>%
  mutate(structure_class = classify_structure(Contraction.sig, Stable.sig, Expansion.sig)) %>%
  ungroup()

n_unmatched <- sig_long %>% filter(is.na(Group)) %>% distinct(Species) %>% nrow()
if (n_unmatched > 0) {
  message(n_unmatched, " species had no Group match in ", group_lookup_csv,
          " -- excluded from the by-group cross-tab (still included in totals).")
}

class_levels <- c("No difference", "All three differ", "C~S!=E", "S~E!=C", "C~E!=S",
                  "Bridging-Expansion", "Bridging-Contraction", "Bridging-Stable",
                  "Partial: Contraction!=Stable", "Partial: Stable!=Expansion",
                  "Partial: Contraction!=Expansion", "Insufficient data")
sig_long <- sig_long %>%
  mutate(structure_class = factor(structure_class, levels = class_levels))

if (!dir.exists(dirname(out_species_csv))) dir.create(dirname(out_species_csv), recursive = TRUE)
write.csv(sig_long, out_species_csv, row.names = FALSE)
cat("Wrote:", out_species_csv, "\n")

class_totals <- sig_long %>%
  count(model, scenario, structure_class, name = "n", .drop = FALSE)
write.csv(class_totals, out_totals_csv, row.names = FALSE)
cat("Wrote:", out_totals_csv, "\n")

class_by_group <- sig_long %>%
  filter(!is.na(Group)) %>%
  count(model, scenario, Group, structure_class, name = "n", .drop = FALSE)
write.csv(class_by_group, out_group_csv, row.names = FALSE)
cat("Wrote:", out_group_csv, "\n")

# Console summary, one model x scenario at a time ----------------------------
for (m in unique(sig_long$model)) {
  for (sc in unique(sig_long$scenario)) {
    cat("\n=== ", m, " | ", sc, ": structure class totals ===\n", sep = "")
    tot <- class_totals %>% filter(model == m, scenario == sc, n > 0) %>%
      arrange(desc(n)) %>% select(structure_class, n)
    print(as.data.frame(tot), row.names = FALSE)

    cat("\n=== ", m, " | ", sc, ": structure class x bird group ===\n", sep = "")
    tab <- class_by_group %>% filter(model == m, scenario == sc) %>%
      select(Group, structure_class, n) %>%
      pivot_wider(names_from = structure_class, values_from = n, values_fill = 0)
    print(as.data.frame(tab), row.names = FALSE)
  }
}

cat("\nDone.\n")
