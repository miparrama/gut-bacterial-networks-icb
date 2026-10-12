# TOPOSCORE calculation on MSI patients treated with ICB.
# Based on https://www.cell.com/cell/fulltext/S0092-8674(24)00538-5 
# @ Miguel Parra Martinez

source("scripts/helpers.R")

# Setup -------------------------------------------------------------------

species_score <- read_xlsx(
  "data/resources/toposcore_species.xlsx",
  skip = 1,
  sheet = 2
) %>%
  mutate(MGS = str_replace_all(MGS, " ", "")) %>%
  mutate(
    MGS = ifelse(
      MGS == "CandidatusCibiobacterqucibialis",
      "CandidatusCibionibacterquicibialis",
      MGS
    )
  ) %>% # fix typo.
  mutate(
    MGS = case_when(
      #Fix names that have changed in new Metaphlan version
      MGS == "ClostridiaceaeunclassifiedSGB4769" ~ "ClostridiaceaebacteriumMarseilleQ4145",
      MGS == "ClostridialesunclassifiedSGB15145" ~ "EubacterialesunclassifiedSGB15145",
      MGS == "CollinsellaSGB14754" ~ "BacteroidesSGB14754",
      MGS == "Enteroclosteraldensis" ~ "Enteroclosteraldenensis",
      MGS == "ClostridiaceaebacteriumOM086BH" ~ "ClostridiaceaebacteriumMarseilleQ4143",
      MGS == "LachnospiraceaebacteriumWCA3601WT6H" ~ "Walteraintestinalis",
      MGS == "Lactobacillusvaginalis" ~ "Limosilactobacillusvaginalis",
      .default = MGS
    )
  )

# Function to compute tse_score
compute_S_score <- function(tse,
                            species_df,
                            assay_name = "counts",
                            agglomerate_rank = "Species") {
  
  #Akkermansia_muciniphila strain to classify gray zone samples. 
  fuso = tse[rowData(tse)$SGB == "SGB9226",] %>% 
    meltSE()
  
  if (dim(fuso)[1] == 0) {
    fuso = tse[rowData(tse)$Species == "Akkermansia_muciniphila",] %>% 
      meltSE()
  }
  
  tse = agglomerateByRank(tse, rank= agglomerate_rank)
  
  # Subset mgs by signature
  sig1_mgs <- species_df %>% filter(SIG == "SIG1") %>% pull(MGS)
  sig2_mgs <- species_df %>% filter(SIG == "SIG2") %>% pull(MGS)
  
  #Extract species:
  sig1_species <- rowData(tse) %>%
    as.data.frame() %>%
    filter (MGS %in% sig1_mgs) %>% 
    pull (agglomerate_rank)
  
  sig2_species <- rowData(tse) %>%
    as.data.frame() %>%
    filter (MGS %in% sig2_mgs) %>% 
    pull (agglomerate_rank)
  
  
  message(paste0(length(sig1_species), "/37 species of SIG1 found"))
  message(paste0(length(sig2_species), "/45 species of SIG2 found"))
  
  # Extract assay data
  abundance_data <- assay(tse, assay_name)
  
  # Subset abundance data for SIG1 and SIG2
  sig1_data <- abundance_data[rownames(abundance_data) %in% sig1_species, , drop = FALSE]
  sig2_data <- abundance_data[rownames(abundance_data) %in% sig2_species, , drop = FALSE]
  
  # Compute proportions of present species
  sig1_present <- colSums(sig1_data > 0) / 37 # SIG1 proportion of present
  sig2_present <- colSums(sig2_data > 0) / 45 # SIG2 proportion of present
  
  # Calculate S-score
  S <- ((sig2_present - sig1_present) + 1) / 2
  
  # Compute the class based on fuso abundance. 
  fuso$S = S
  fuso = fuso %>% 
    mutate( Sig = case_when (
      S > 0.7911 ~ "Sig2", 
      S < 0.5351 ~ "Sig1",
      .default = "gray_zone")) %>% 
    mutate(Sig = case_when (
      (Sig == "gray_zone") & (counts == 0) ~ "Sig1", 
      (Sig == "gray_zone") & (counts > 15) ~ "Sig1",
      (Sig == "gray_zone") & (counts < 15) & (counts > 0 ) ~ "Sig2",
      .default = Sig)) %>% 
    column_to_rownames("SampleID")
  
  sig = fuso$Sig
  names(sig) <- rownames (fuso)
  
  result = list ("s.score" = S, "toposcore" = sig)
  
  return(result)
}



# ------------------------------------------------------------------------#
# ICB-cohort -------------------------------------------------------------
# ------------------------------------------------------------------------#

# tse_int.species <- read_rds("output/data/metagenomics/microbeICB/tse_species.rds")
tse_int.SGB <- read_rds("output/data/metagenomics/microbeICB/tse_int.SGB.rds")

rowData(tse_int.SGB)$MGS <- rowData(tse_int.SGB)$Species %>% str_replace_all("_", "")
rowData(tse_int.SGB)$SIG <- rowtibble(tse_int.SGB) %>%
  left_join(species_score) %>% 
  pull(SIG)
colData(tse_int.SGB)$s_score = compute_S_score(tse_int.SGB,
                                                   species_score,
                                                   agglomerate_rank = "Species")$s.score

write_rds(tse_int.SGB, "output/data/metagenomics/microbeICB/tse_int.SGB.rds")
