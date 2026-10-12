# Merge datasets
# @ Miguel Parra Martinez

here::i_am("scripts/combine_datasets.R")
library(here)
source(here("scripts/helpers.R"))

source(here("scripts/studies/Derosa.R"))
source(here("scripts/studies/Gunjur.R"))
source(here("scripts/studies/Healthy_500fg.R"))
source(here("scripts/studies/Kim.R"))
source(here("scripts/studies/Kim_1st.R"))
source(here("scripts/studies/Lee.R"))
source(here("scripts/studies/Routy.R"))
source(here("scripts/studies/Spencer.R"))
source(here("scripts/studies/Drup.R"))

# Merge data -------------------------------------------------------------------

## Species level ## ------------------------------------------------------------

# Merge metaphlan tables
int.species = merge_phyloseq(
  aggregate_taxa(lee.pseq, "Species"),
  aggregate_taxa(derosa.pseq, "Species"),
  aggregate_taxa(routy.pseq, "Species"),
  aggregate_taxa(spencer.pseq, "Species"),
  aggregate_taxa(gunjur.pseq, "Species"),
  aggregate_taxa(kim.pseq, "Species"),
  aggregate_taxa(kim_1st.pseq, "Species"), 
  aggregate_taxa(drup.pseq, "Species")
)

write_rds(int.species, here("output/data/merge/int.species.rds"))

## SGB level ## ------------------------------------------------------------
# Strip the tree from each object
strip_tree <- function(pseq) {
  phyloseq::phyloseq(
    otu_table(pseq),
    tax_table(pseq),
    sample_data(pseq)
  )
}

int.sgb <- merge_phyloseq(
  strip_tree(lee.pseq),
  strip_tree(derosa.pseq),
  strip_tree(routy.pseq),
  strip_tree(spencer.pseq),
  strip_tree(gunjur.pseq),
  strip_tree(kim.pseq),
  strip_tree(kim_1st.pseq),
  strip_tree(drup.pseq)
)

# Add tree again
mpa_tree = ape::read.tree(here(
  "data/resources/mpa_vOct22_CHOCOPhlAnSGB_202212.nwk"
))
mpa_tree$tip.label = paste0("SGB", mpa_tree$tip.label)
int.sgb = merge_phyloseq(int.sgb, mpa_tree)

write_rds(int.sgb, here("output/data/merge/int.sgb.rds"))


# Adjust Batch effects ------------------------------------------------------------

## Species ------------------------------------------------------------

# Batch correct species
int.species.comp = int.species %>%
  microbiome::transform("compositional")

fit_adjust_batch <- MMUPHin::adjust_batch(
  feature_abd = abundances(int.species.comp),
  batch = "Study",
  covariates = c("Response"),
  data = meta(int.species.comp),
  control = list(verbose = T, diagnostic_plot = NULL)
)

int.species.comp.adj = int.species.comp
otu_table(int.species.comp.adj) = otu_table(
  fit_adjust_batch$feature_abd_adj,
  taxa_are_rows = T
)

write_rds(
  int.species.comp.adj,
  here("output/data/merge/int.species.comp.adj.rds")
)

## SGB ------------------------------------------------------------

## SGBS
int.sgb.comp = int.sgb %>%
  microbiome::transform("compositional")

fit_adjust_batch <- MMUPHin::adjust_batch(
  feature_abd = abundances(int.sgb.comp),
  batch = "Study",
  covariates = c("Response"),
  data = meta(int.sgb.comp),
  control = list(verbose = T, diagnostic_plot = NULL)
)

int.sgb.comp.adj = int.sgb.comp
otu_table(int.sgb.comp.adj) = otu_table(
  fit_adjust_batch$feature_abd_adj,
  taxa_are_rows = T
)

write_rds(
  int.sgb.comp.adj,
  here("output/data/merge/int.sgb.comp.adj.rds")
)




# Functional data ------------------------------------------------------------

# Merge metacyc tables
int.metacyc = merge_phyloseq(
  lee.metacyc,
  derosa.metacyc,
  routy.metacyc,
  spencer.metacyc,
  gunjur.metacyc,
  kim.metacyc,
  kim_1st.metacyc,
  drup.metacyc
)
write_rds(int.metacyc, "output/data/merge/int.metacyc.rda")

# Merge metacyc tables with taxonomy
int.metacyc = merge_phyloseq(
  lee.metacyc.tax,
  derosa.metacyc.tax,
  routy.metacyc.tax,
  spencer.metacyc.tax,
  gunjur.metacyc.tax,
  kim.metacyc.tax,
  kim_1st.metacyc.tax,
  drup.metacyc.tax
)
write_rds(int.metacyc, "output/data/merge/int.metacyc.tax.rda")


# Merge KEGG tables
int.kegg = merge_phyloseq(
  lee.kegg,
  derosa.kegg,
  routy.kegg,
  spencer.kegg,
  gunjur.kegg,
  kim.kegg,
  kim_1st.kegg,
  drup.kegg
)
write_rds(int.kegg, "output/data/merge/int.kegg.rda")

# Merge KEGG tables with taxonomy
int.kegg = merge_phyloseq(
  lee.kegg.tax,
  derosa.kegg.tax,
  routy.kegg.tax,
  spencer.kegg.tax,
  gunjur.kegg.tax,
  kim.kegg.tax,
  kim_1st.kegg.tax,
  drup.kegg.tax
)
write_rds(int.kegg, "output/data/merge/int.kegg.tax.rda")

# Merge eggNOG tables
int.eggnog = merge_phyloseq(
  lee.eggnog,
  derosa.eggnog,
  routy.eggnog,
  spencer.eggnog,
  gunjur.eggnog,
  kim.eggnog,
  kim_1st.eggnog,
  drup.eggnog
) %>%
  subset_taxa(!grepl("\\|", Species))
write_rds(int.eggnog, "output/data/merge/int.eggnog.rda")

