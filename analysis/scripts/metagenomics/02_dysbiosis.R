# Script to compute the dysbiosis score.  
# @ Miguel Parra Martinez
# 2025-01

source("scripts/helpers.R")

library (dysbiosisR)

# Microbe ICB -------------------------------------------------------------

# Load data

tse_eggnog <- read_rds(
  "data/metagenomics/int.eggnog.rda"
) %>%
  mia::convertFromPhyloseq() %>%
  transformAssay(method = "relabundance") %>%
  transformAssay(method = "pa") %>% 
  transformAssay(
    method = "clr",
    assay.type = "relabundance",
    pseudocount = T
  ) 

healthy_eggnog_tse <- read_rds(
  "output/data/metagenomics/healthy/tse_healthy500fg_eggnog_families.rds"
)

## Eggnog ------------------------------------------------------------------

# Remove of metadata that are not the same class, and add type info.
colData(tse_eggnog) <- coldf(tse_eggnog) %>%
  select(Sample, timepoint, Response) %>%
  mutate(Type = "Cancer") %>%
  DataFrame()

colData(healthy_eggnog_tse) <- coldf(healthy_eggnog_tse) %>%
  select(run) %>%
  mutate(Type = "Healthy") %>%
  DataFrame()

# Merge tse objects
merge_eggnog <- mergeSEs(tse_eggnog,
                         healthy_eggnog_tse,
                          join = "inner")

# Transform to phyloseq
merge_eggnog_pseq <- merge_eggnog %>% 
  convertToPhyloseq() %>% 
  microbiome::transform("clr") 

# Generate counts
dist_mat_eggnog <- phyloseq::distance(merge_eggnog_pseq, "euclidean")

dysbiosis_eggnog <- dysbiosisR::euclideanDistCentroids(
  merge_eggnog_pseq,
  dist_mat = dist_mat_eggnog,
  use_squared = TRUE,
  group_col = "Type",
  control_label = "Healthy",
  case_label = "Cancer") %>% 
  filter(Type != "Healthy") %>% 
  mutate(dysbiosis_eggnog = scale(CentroidDist_score)[,1])

dysbiosis_eggnog %>% write_csv("output/data/metagenomics/microbeICB/icb-eggnog-dysbiosis.csv")