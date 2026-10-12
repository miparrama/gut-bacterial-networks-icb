#-----------------------------------------------------------------------------#
# Project: Bacterial communities
# Script:  04_bacterial_communities.R
# Purpose: Unsupervised clustering of bacterial communities in immunotherapy treated patients
# Author:  Miguel Parra-Martinez
# Date:    2025-12-24
#-----------------------------------------------------------------------------#

source("scripts/helpers.R")

# Setup ----------------------------------------------------------------------

# Filter tse

# Load oral bacteria
# Supplementary file from Kim, H. et al. Multi-omic analysis reveals transkingdom gut dysbiosis in metabolic dysfunction-associated steatotic liver disease. Nature Metabolism 2025 7:7 7 (2025). https://doi.org/10.1038/s42255-025-01318-6
oral_bacteria <- read_xlsx("data/resources/s42255-025-01318-6_tables.xlsx", sheet = 15) %>% 
  mutate(SGB = str_extract(speciesname, "SGB\\d+"))

# Load TSE object
tse_int.sgb <- read_rds("output/data/metagenomics/microbeICB/tse_int.SGB.rds")

# Add oral data. 
rowData(tse_int.sgb)$oral <- ifelse(rowData(tse_int.sgb)$SGB %in% oral_bacteria$SGB, "Oral", NA)

# Alpha diversity metrics. 
index_names = c("inverse_simpson_diversity", "shannon_diversity", "absolute_dominance", "simpson_evenness", "observed_richness", "gini_dominance", "pielou")
index = c("inverse_simpson", "shannon", "absolute", "simpson", "observed", "gini_dominance", "pielou_evenness")
tse_int.sgb <- addAlpha(
  tse_int.sgb, assay.type = "counts", index = index, name = index_names)

## Filter samples. 
metadata_pretreatment = coldf(tse_int.sgb) %>% 
  mutate(host_disease_status = as.factor(host_disease_status)) %>% 
  filter(timepoint == "PreTreatment") %>% 
  mutate(batch = paste(Study, cohort, sep = ".") %>% str_replace(".NA", "")) %>% 
  mutate(n = n(), .by = batch) %>% 
  filter(n > 20) %>% 
  mutate(Study = ifelse(str_detect(Study, "Kim"), "Kim", Study))

metadata = coldf(tse_int.sgb) %>%
  mutate(host_disease_status = as.factor(host_disease_status)) %>% 
  mutate(batch = paste(Study, cohort, sep = ".") %>% str_replace(".NA", "")) %>% 
  mutate(n = n(), .by = batch) %>% 
  filter(n > 20) %>% 
  mutate(Study = ifelse(str_detect(Study, "Kim"), "Kim", Study))

metadata %>% 
  tabyl(host_disease_status) %>% 
  arrange(n)

tse_int.sgb.fil.pretreatment <- tse_int.sgb[,rownames(metadata_pretreatment)]
colData(tse_int.sgb.fil.pretreatment) <- DataFrame(metadata_pretreatment)

tse_int.sgb.fil <- tse_int.sgb[,rownames(metadata)]
colData(tse_int.sgb.fil) <- DataFrame(metadata)

# Same but agglomerating by Species
tse_int.species <- read_rds("output/data/metagenomics/microbeICB/tse_int.SGB.rds") %>% 
  agglomerateByRank("Species") # No filtering by prevalence

tse_int.species.fil.pretreatment <- tse_int.species[,rownames(metadata_pretreatment)]
colData(tse_int.species.fil.pretreatment) <- DataFrame(metadata_pretreatment)

tse_int.species.fil <- tse_int.species[,rownames(metadata)]
colData(tse_int.species.fil) <- DataFrame(metadata)

# Metadata plot:
metadata %>%
  mutate(n_study = n(), .by = Study) %>%
  mutate(Study = fct_reorder(Study, n_study)) %>%
  ggplot(aes(Study)) +
  geom_bar(aes(fill = Response)) +
  coord_flip() +
  theme_classic(base_size = 16) +
  theme(
    legend.position = c(0.98, 0.02),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(
      fill = alpha("white", 0.6),
      color = NA
    ),
  )

## Figure 1.A ##
ggsave(
  "output/figures/figures_manuscript/studies_overview.pdf",
  height = 4,
  width = 5
)

# Summary table ------------------------------------------------------------

# data <- read.csv("sample_data.csv", stringsAsFactors = FALSE)

metadata <- read_csv("output/data/tables_manuscript/metadata_cohort.csv")

library(gt)
library(gtsummary)
clinical.char.plot = metadata %>%
  mutate(
    atb = replace_na(atb, "Unknown") %>% 
      str_to_sentence() %>% 
      fct_relevel(atb, c("Yes", "No", "Unknown")),
    host_disease_status = fct_relevel(host_disease_status,
      c(
        "NSCLC",
        "MELA",
        "RCC",
        "CRC",
        "Other-GI",
        "Other-non-GI"
      )
    ), 
    Response = fct_relevel(Response, "Responder"),
    timepoint = fct_relevel(timepoint, "Pre-treatment")
  ) %>%
  rename(
    `Cancer type` = host_disease_status,
    `Objective response` = Response,
    `Sex` = gender,
    `Collection` = timepoint,
    `Antibiotic use` = atb,
    `Age` = age,
    `Species detected` = observed_richness
  ) %>%
  select(
    Study,
    `Objective response`,
    `Cancer type`,
    `Antibiotic use`,
    Sex,
    `Collection`,
    `Age`
  ) %>%
  tbl_summary(by = Study) %>%
  bold_labels() 
# %>%
#   as_gt()

library(huxtable)
library(openxlsx)

gtsummary::as_hux_xlsx(clinical.char.plot, file = "output/data/tables_manuscript/cohort_summary.xlsx")





# ------------------------------------------------------------ #
# Clustering  ICB cohort (fastspar) --------------------------
# ------------------------------------------------------------ #

fastspar_icbm <- read_tsv ("../fastspar/output/correlation/microbe_icb_2_correlation.tsv")
fastspar_icbm_pval <- read_tsv ("../fastspar/output/pvalues/microbe_icb_2_pvalues.tsv")

fastspar_icbm_mt <- fastspar_icbm %>% 
  column_to_rownames(var = "#OTU ID") %>% 
  as.matrix()

fastspar_icbm_pval <- fastspar_icbm_pval %>% 
  column_to_rownames(var = "#OTU ID") %>% 
  as.matrix()

species_to_sgb <- read_csv("output/data/metagenomics/microbeICB/filter_Species_SGB_df.csv")

rownames(fastspar_icbm_mt) <- species_to_sgb$SGB
colnames(fastspar_icbm_mt) <- species_to_sgb$SGB

rownames(fastspar_icbm_pval) <- species_to_sgb$SGB
colnames(fastspar_icbm_pval) <- species_to_sgb$SGB

fastspar_icbm_mt_fil = fastspar_icbm_mt

# Sparsify matrix
fastspar_icbm_mt_fill = fastspar_icbm_mt

fastspar_icbm_mt_fill[p.adjust(fastspar_icbm_pval, "fdr") > 0.05] <- 0
fastspar_icbm_mt_fill[abs(fastspar_icbm_mt_fill) < 0.2] <- 0


## Generate OTU metadata table ------------------------------------------------------------

# Metadata: 
tax_info_Jan21 <- read_csv("data/resources/taxonomy_info_vJan21_CHOCOPhlAnSGB_202103.csv")
tax_info_Jan21 <- tax_info_Jan21 %>%
  mutate(
    species = str_extract(SGBs, "s__[^|]+") %>%
      str_remove("^s__"),
    species_ncbi = str_split(NCBI_tax, "\\|") %>%
      map_chr(~ .x[7]),
    SGB = str_remove_all(SGBs, "_group") %>% str_extract("(SGB|EUK)\\d+$")
  )

toposcore_species <- read_xlsx(
  "data/resources/toposcore_species.xlsx",
  skip = 1,
  sheet = 2
) %>%
  mutate(species = MGS %>% str_replace_all(" ", "_")) %>%
  left_join(tax_info_Jan21 %>% distinct(species, .keep_all = T)) %>%
  mutate(
    Species = case_when(
      MGS ==
        "Candidatus Cibiobacter qucibialis" ~ "Candidatus Cibionibacter quicibialis",
      MGS ==
        "Clostridiacea eunclassified SGB4769" ~ "Clostridiaceae bacterium Marseille Q4145",
      MGS ==
        "Clostridiales unclassified SGB15145" ~ "Eubacterialesunc lassified SGB15145",
      MGS == "Enterocloster aldensis" ~ "Enterocloster aldenensis",
      MGS ==
        "Clostridiaceae bacterium OM08 6BH" ~ "Clostridiaceae bacterium Marseille Q4143",
      MGS == "Lachnospiraceae bacterium WCA3 601 WT 6H" ~ "Waltera intestinalis",
      MGS == "Lactobacillus vaginalis" ~ "Limosilactobacillus vaginalis",
      .default = MGS
    ) %>%
      str_replace_all(" ", "_")
  ) %>% 
  select(-c(SGB))

tax_info_Oct22<- read_csv("data/resources/taxonomy_info_vOct22_CHOCOPhlAnSGB_202403.csv")
tax_info_Oct22 <- tax_info_Oct22 %>%
  mutate(
    Species = str_extract(SGBs, "s__[^|]+") %>%
      str_remove("^s__"),
    Species_ncbi = str_split(NCBI_tax, "\\|") %>%
      map_chr(~ .x[7]),
    SGB = str_remove_all(SGBs, "_group") %>% str_extract("(SGB|EUK)\\d+$")
  )

fit_lm_sgb <- read_rds("output/data/metagenomics/microbeICB/meta_analysis_SGB.rds")
metanalysis <- fit_lm_sgb$meta_fits %>%
  rename(SGB = feature) %>%
  as_tibble() %>%
  select(SGB, pval, qval.fdr, coef) %>%
  mutate(sign_pval = -log10(pval) * sign(coef)) %>%
  mutate(
    sig = case_when(
      pval < 0.05 & coef > 0 ~ "R_associated",
      pval < 0.05 & coef < 0 ~ "NR_associated",
      .default = NA
    )
  )

metadata_icbm_hm <- data.frame(
  SGB = rownames(fastspar_icbm_mt_fill)
) %>%
  left_join(metanalysis, by = "SGB") %>%
  left_join(
    tax_info_Oct22 %>%
      select(-SGBs)
  ) %>%
  left_join(
    toposcore_species %>% select(-NCBI_tax)
  ) %>%
  left_join(rowtibble(tse_int.sgb) %>% select(oral, SGB)) %>%
  mutate(Species_2 = SGB) %>%
  column_to_rownames("Species_2") 


## Clustering -----------------------------------------------------------------
clustering_distance = "euclidean"
clustering_method = "ward.D2"

set.seed(123)
kmeans.out <- clusterRows(fastspar_icbm_mt_fill, KmeansParam(8))
hclust.out <- clusterRows(
  fastspar_icbm_mt_fill,
  HclustParam(cut.params = list(k = 3), method = clustering_method)
)
hclust.out_small <- clusterRows(fastspar_icbm_mt_fill, HclustParam(cut.params = list (k = 5), method = clustering_method))

metadata_icbm_hm$cluster <- hclust.out
metadata_icbm_hm$cluster_small <- hclust.out_small

# Rename clusters 
metadata_icbm_hm <- metadata_icbm_hm %>% 
  mutate(cluster_small = 
  case_when(
    cluster_small == 1 ~ 5,
    cluster_small == 2 ~ 4,
    cluster_small == 3 ~ 2,
    cluster_small == 4 ~ 1,
    cluster_small == 5 ~ 3
  ))

colors_clusters <- c(
  "1" = "#4DAF4A",
  "2" = "#E41A1C",
  "3" = "#377EB8",
  "4" = "#984EA3",
  "5" = "#999999"
)


## Heatmap --------------------------------------------------------------------
set.seed(123)
hm_icb_fastspar <- Heatmap(
  fastspar_icbm_mt_fill,
  name = "FastSpar\ncorrelation",
  column_title = "Microbe-ICB",
  show_column_names = F,
  show_row_names = F,
  row_names_side = "left",
  row_dend_side = "right",
  row_title = "Bacterial species",
  row_title_side = "left",
  col = colorRamp2(c(-0.6, 0, 0.6), c("blue", "white", "red")),
  clustering_distance_columns = clustering_distance,
  clustering_distance_rows = clustering_distance,
  clustering_method_rows = clustering_method,
  clustering_method_columns = clustering_method,
  top_annotation = columnAnnotation(
    cluster = metadata_icbm_hm$cluster_small,
    "-log10(pval)" = metadata_icbm_hm$sign_pval,
    "pval < 0.05" = metadata_icbm_hm$sig,
    col = list(
      "-log10(pval)" = colorRamp2(c(-4, 0, 4), c("red", "white", "blue")),
      "pval < 0.05" = c("R_associated" = "blue", "NR_associated" = "Red"),
      cluster = colors_clusters
    ),
    na_col = "white"
  ),
  left_annotation = rowAnnotation(
    show_legend = FALSE,
    cluster = metadata_icbm_hm$cluster_small,
    "-log10(pval)" = metadata_icbm_hm$sign_pval,
    "pval < 0.05" = metadata_icbm_hm$sig,
    col = list(
      "-log10(pval)" = colorRamp2(c(-4, 0, 4), c("red", "white", "blue")),
      "pval < 0.05" = c("R_associated" = "blue", "NR_associated" = "Red"),
      cluster = colors_clusters
    ),
    na_col = "white"
  ),
  use_raster = T,
  raster_quality = 5
)

hm_icb_fastspar

## Figure S1.B ## 
CairoPDF("output/figures/figures_manuscript/ICB_fastspar_clusters.pdf",
         width = 10, height = 7)
hm_icb_fastspar
dev.off()

write.csv(fastspar_icbm_mt_fill, "output/data/metagenomics/microbeICB/fastspar_icbm_mt_fill.csv", row.names = TRUE)


## Small Version:  ------------------------------------------------------------

metadata_icbm_hm_small <- metadata_icbm_hm %>% 
  filter(cluster_small %in% c(1, 2))

# Save data
write_csv(metadata_icbm_hm_small, "output/data/metagenomics/microbeICB/species_cluster_1_2.csv")

fastspar_icbm_mt_fill_small <- fastspar_icbm_mt_fill[
  rownames(metadata_icbm_hm_small),
  rownames(metadata_icbm_hm_small)
]

set.seed(123)
hm_icb_fastspar_small <- Heatmap(
  fastspar_icbm_mt_fill_small,
  name = "fastspar\ncorrelation",
  column_title = "Microbe-ICB",
  show_column_names = F,
  show_row_names = F,
  row_names_side = "left",
  row_dend_side = "right",
  row_title = "Bacterial species",
  row_title_side = "left",
  col = colorRamp2(c(-0.6, 0, 0.6), c("blue", "white", "red")),
  clustering_distance_columns = clustering_distance,
  clustering_distance_rows = clustering_distance,
  clustering_method_rows = clustering_method,
  clustering_method_columns = clustering_method,
  top_annotation = columnAnnotation(
    cluster = metadata_icbm_hm_small$cluster_small,
    "-log10(pval)" = metadata_icbm_hm_small$sign_pval,
    "pval < 0.05" = metadata_icbm_hm_small$sig,
    col = list(
      "-log10(pval)" = colorRamp2(c(-4, 0, 4), c("red", "white", "blue")),
      "pval < 0.05" = c("R_associated" = "blue", "NR_associated" = "Red"),
      cluster = colors_clusters
    ),
    na_col = "white"
  ),
  left_annotation = rowAnnotation(
    show_legend = FALSE,
    cluster = metadata_icbm_hm_small$cluster_small,
    "-log10(pval)" = metadata_icbm_hm_small$sign_pval,
    "pval < 0.05" = metadata_icbm_hm_small$sig,
    col = list(
      "-log10(pval)" = colorRamp2(c(-4, 0, 4), c("red", "white", "blue")),
      "pval < 0.05" = c("R_associated" = "blue", "NR_associated" = "Red"),
      cluster = colors_clusters
    ),
    na_col = "white"
  ),
  use_raster = T
)

hm_icb_fastspar_small

## Figure 1.B ## 
CairoPDF("output/figures/figures_manuscript/ICB_fastspar_clusters_small.pdf",
         width = 9, height = 7)
hm_icb_fastspar_small
dev.off()

## Save data meta-analysis ------------------------------------------------------------

fit_lm_sgb <- read_rds("output/data/metagenomics/microbeICB/meta_analysis_SGB.rds")
metanalysis <- fit_lm_sgb$meta_fits %>%
  rename(SGB = feature) %>%
  as_tibble() %>%
  # select(SGB, pval, qval.fdr, coef) %>%
  mutate(sign_pval = -log10(pval) * sign(coef)) %>%
  mutate(
    sig = case_when(
      pval < 0.05 & coef > 0 ~ "R_associated",
      pval < 0.05 & coef < 0 ~ "NR_associated",
      .default = NA
    )
  ) %>% 
  left_join(species_to_sgb, by = "SGB") %>% 
  relocate(Species, .after = SGB) %>% 
  left_join(
    metadata_icbm_hm %>% 
      select(SGB, cluster_small) %>% 
      rename(cluster = cluster_small)
  )

openxlsx::write.xlsx(
  c(list("meta_fits" = metanalysis), fit_lm_sgb$maaslin_fits),
  file = "output/data/tables_manuscript/maaslin_fits_per_study.xlsx",
  asTable = TRUE
)



# Meta-analysis barplot ------------------------------------------------------------

## Top SGB significant ------------------------------------------------------------
# plot SGB significant

# Load results meta-analysis. 
fit_lm_sgb <- read_rds("output/data/metagenomics/microbeICB/meta_analysis_SGB.rds")

fit_lm_sgb.res = fit_lm_sgb$meta_fits %>% 
  mutate(fdr = p.adjust(pval, "fdr"), .after = pval) %>% 
  mutate(val = -log10(pval) * sign(coef)) %>% 
  mutate(SGB = feature) %>% 
  as_tibble() %>% 
  left_join(rowtibble(tse_int.sgb.fil.pretreatment)) %>% 
  relocate(Species) %>% 
  arrange (pval)

fit_lm_sgb.res.fil <- fit_lm_sgb.res %>%
  left_join(metadata_icbm_hm %>% select(SGB, cluster_small)) %>%
  filter(pval < 0.055) %>%
  mutate(n_sgb = n(), .by = Species) %>%
  mutate(
    # Species = ifelse(n_sgb == 2, paste(Species, SGB, sep = "_"), Species)
    Species = ifelse(
      str_detect(Species, SGB),
      Species,
      paste(Species, SGB, sep = "_")
    )
  ) %>%
  mutate(
    Cluster = cluster_small,
    sign_pval = -log10(pval) * sign(coef),
    Species = fct_reorder(Species, sign_pval),
    Association = ifelse(coef > 0, "Beneficial", "Detrimental")
  ) %>%
  mutate(Cluster = ifelse(Cluster %in% c(1, 2), Cluster, "Other"))

# Sign p-value. 
# ------------------------------------------------------------

barplot_pval <- fit_lm_sgb.res.fil %>%
  ggplot(aes(x = sign_pval, y = Species)) +
  # geom_col(aes(fill = SIG)) +
  geom_col(aes(fill=Cluster)) +
  geom_vline(xintercept = c(log10(0.05), -log10(0.05)), linetype = "dashed") +
  labs(
    x = "Signed p-value meta-analysis",
    fill = "Cluster"
  ) +
  # scale_fill_manual(values = colors_vec) +
  # scale_fill_manual(
  #   values = c("SIG2" = "#7EB17E", "SIG1" = "black"),
  #   # na.value = "black"
  # ) +
  scale_fill_manual(
    # values =  map_colors(metadata_icbm_hm$cluster_small, palette_name = "Paired")
    values = colors_clusters
  ) + 
  theme_minimal(base_size = 14)
barplot_pval


# Coefficient per study. 
# ------------------------------------------------------------

# Significant species in meta-analysis.
coef_hm <- map_df(names(fit_lm_sgb$maaslin_fits), function(study) {
  res <- fit_lm_sgb$maaslin_fits[[study]]
  res$Study <- study
  res
}) %>%
  mutate(SGB = feature) %>%
  inner_join(fit_lm_sgb.res.fil %>% select(SGB, Species, sign_pval)) %>%
  mutate(
    pval.res.sign = -log10(pval) * sign(coef),
    # pval.res.sign = -log10(pval.res) * sign(coef.res),
    Species = fct_reorder(Species, sign_pval)
  ) %>%
  # filter(pval.res < 0.05) %>%
  ggplot(aes(Study, Species)) +
  geom_tile(aes(fill = coef)) +
  geom_text(
    aes(
      label = ifelse(pval < 0.05, "x", ""),
      color = ifelse(abs(coef) > 0.013, "white", 'black')
    ),
    show.legend = F,
    # color = "black",
    size = 2.8
  ) +
  labs(
    subtitle = "Coefficient",
    fill = "Coefficient\nMMUPHin"
  ) +
  scale_fill_gradient2(
    limits = c(-0.02, 0.02),
    oob = scales::squish,
  ) +
  scale_color_manual(values = c("black" = "black", "white" = "white")) +
  theme_minimal(base_size = 14)
coef_hm


# Number of bacteria with significance in more than 1 study:
number_studies_significance <- map_df(names(fit_lm_sgb$maaslin_fits), function(study) {
  res <- fit_lm_sgb$maaslin_fits[[study]]
  res$Study <- study
  res
}) %>%
  mutate(SGB = feature) %>%
  inner_join(fit_lm_sgb.res.fil %>% select(SGB, Species, sign_pval)) %>%
  mutate(
    pval.res.sign = -log10(pval) * sign(coef),
    # pval.res.sign = -log10(pval.res) * sign(coef.res),
    Species = fct_reorder(Species, sign_pval)
  )

number_studies_significance %>% 
  filter(pval < 0.05) %>% 
  group_by(feature) %>% 
  summarize (n = n()) %>% 
  ungroup() %>% 
  tabyl(n)


# Prevalence difference 
# ------------------------------------------------------------
prevalence_df = tse_int.sgb.fil.pretreatment %>%
  transformAssay(method = "pa") %>% 
  meltSE(assay.type = "pa") %>%
  rename(SGB = FeatureID, Sample = SampleID) %>%
  inner_join(metadata) %>%
  group_by(SGB, Study, Response) %>%
  summarize(prevalence = mean(pa)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Response", values_from = "prevalence") %>% 
  mutate(difference = Responder - `Non-responder`)

prevalence_hm <- prevalence_df %>% 
  inner_join(fit_lm_sgb.res.fil %>% select(SGB, Species, sign_pval)) %>% 
  filter(!is.na(Species)) %>% 
  ggplot(aes(Study, Species)) +
  geom_tile(aes(fill=difference)) +
  scale_fill_gradient2(
    limits = c(-0.3, 0.3),       # fix color scale range
    oob = scales::squish,    # values outside [-5, 5] are mapped to extremes
    # low = "blue", mid = "white", high = "red", midpoint = 0
  ) +
  labs(
    fill = "Prevalence\ndifference",
    subtitle = "Prevalence difference"
  ) +
  theme_minimal(base_size = 14)

# Prevalence per study. 
# ------------------------------------------------------------
n_study = coltibble(tse_int.sgb.fil.pretreatment) %>% 
  count(Study) %>% 
  rename(n_study = n)

prevalence_per_study <- tse_int.sgb.fil.pretreatment %>%
  meltSE(assay_name = "counts", add.col = T) %>%
  rename(SGB = FeatureID) %>%
  inner_join(fit_lm_sgb.res.fil %>% select(SGB, Species, sign_pval)) %>% 
  mutate(SGB = fct_reorder(SGB, sign_pval)) %>%
  mutate(pa = ifelse(counts > 0, 1, 0)) %>%
  group_by(SGB, Study) %>%
  summarize(n = n(), preval = sum(pa) / n()) %>%
  ungroup() %>%
  distinct(SGB, Study, preval)

# Plot prevalence
prevalence_per_study_p <- prevalence_per_study %>%
  ggplot(aes(y = SGB, x = preval, color = Study)) +
  geom_point() +
  scale_color_manual(
    values = colors_studies
  ) +
  labs(
    subtitle = "Prevalence",
    x = "Prevalence"
  ) +
  theme_bw(base_size = 14)


# Mean abundance per study. 
# ------------------------------------------------------------

mean_relab_per_study <- tse_int.sgb.fil %>%
  meltSE(assay_name = "relabundance_adj", add.col = T) %>%
  rename(SGB = FeatureID) %>%
  inner_join(fit_lm_sgb.res.fil %>% select(SGB, Species, sign_pval)) %>% 
  mutate(SGB = fct_reorder(SGB, sign_pval)) %>%
  group_by(SGB, Study) %>%
  summarize(mean_relab = mean(relabundance_adj, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(SGB, Study, mean_relab) %>%
  replace_na(list(mean_relab = 0))

relab_per_study_p <- mean_relab_per_study %>%
  inner_join(fit_lm_sgb.res.fil %>% select(SGB, Species, sign_pval)) %>%
  mutate(Species = fct_reorder(Species, sign_pval)) %>%
  ggplot(aes(y = Species, x = mean_relab, color = Study)) +
  geom_point() +
  scale_color_manual(
    values = colors_studies
  ) +
  scale_x_log10() +
  # scale_y_discrete(position = "right") +
  labs(
    x = "Mean relative abundance",
    subtitle = "Abundance"
  ) +
  theme_bw(base_size = 14)


# Meta-analysis plot combined
# ------------------------------------------------------------
patch <- (barplot_pval +
  theme(
    axis.text.y = element_text(size = 10)
  )) +
  (coef_hm +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )) +
  (prevalence_hm +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )) +
  (prevalence_per_study_p +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank()
    )) +
  (relab_per_study_p +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )) +
  plot_layout(
    guides = "collect",
    nrow = 1,
    widths = c(2 / 10, 2 / 10, 2 / 10, 2 / 10, 2 / 10)
  ) +
  plot_annotation(
    title = "Antagonistic clusters DRUP-MSI (fastspar)",
    subtitle = "Meta-analysis with MMUPhin adjusting for ATB use and cancer type. Filter p<0.1"
  )
patch

# Figure S1.A
ggsave(
  "output/figures/figures_manuscript/consortia_ranking_microbeICB_significant_0.05_2.pdf",
  patch,
  height = 13,
  width = 18
)

## SGBs of clusters of interest  ------------------------------------------------------------

fit_lm_sgb <- read_rds("output/data/metagenomics/microbeICB/meta_analysis_SGB.rds")

fit_lm_sgb.res = fit_lm_sgb$meta_fits %>% 
  mutate(fdr = p.adjust(pval, "fdr")) %>% 
  mutate(val = -log10(pval) * sign(coef)) %>% 
  mutate(SGB = feature) %>% 
  as_tibble() %>% 
  left_join(rowtibble(tse_int.sgb.fil.pretreatment)) %>% 
  relocate(Species) %>% 
  arrange (pval)

fit_lm_sgb.res.fil <- fit_lm_sgb.res %>%
  left_join(metadata_icbm_hm %>% select(SGB, cluster_small)) %>%
  filter(cluster_small %in% c(1, 2)) %>%
  mutate(
    Species = ifelse(
      str_detect(Species, SGB),
      Species,
      paste(Species, SGB, sep = "_")
    )
  ) %>%
  mutate(
    Cluster = cluster_small,
    sign_pval = -log10(pval) * sign(coef),
    Species = fct_reorder(Species, sign_pval)
  ) %>%
  mutate(
    Cluster = ifelse(Cluster %in% c(1, 2), Cluster, "Other") %>% as.factor()
  )

barplot_pval <- fit_lm_sgb.res.fil %>%
  ggplot(aes(x = sign_pval, y = Species)) +
  geom_col(aes(fill=Cluster)) +
  geom_vline(xintercept = c(log10(0.05), -log10(0.05)), linetype = "dashed") +
  labs(
    x = "Signed p-value meta-analysis",
    y = NULL,
    fill = "Cluster"
  ) +
  scale_fill_manual(
    values = colors_clusters
  ) + 
  theme_minimal(base_size = 14)
barplot_pval

# Significant species in meta-analysis.
coef_hm <- map_df(names(fit_lm_sgb$maaslin_fits), function(study) {
  res <- fit_lm_sgb$maaslin_fits[[study]]
  res$Study <- study
  res
}) %>%
  mutate(SGB = feature) %>%
  inner_join(fit_lm_sgb.res.fil %>% select(SGB, Species, sign_pval)) %>%
  mutate(
    pval.res.sign = -log10(pval) * sign(coef),
    # pval.res.sign = -log10(pval.res) * sign(coef.res),
    Species = fct_reorder(Species, sign_pval)
  ) %>%
  # filter(pval.res < 0.05) %>%
  ggplot(aes(Study, Species)) +
  geom_tile(aes(fill = coef)) +
  geom_text(
    aes(
      label = ifelse(pval < 0.05, "x", ""),
      color = ifelse(abs(coef) > 0.013, "white", 'black')
    ),
    show.legend = F,
    # color = "black",
    size = 2.8
  ) +
  labs(
    subtitle = "Coefficient",
    fill = "Coefficient\nMMUPHin"
  ) +
  scale_fill_gradient2(
    limits = c(-0.02, 0.02),
    oob = scales::squish,
  ) +
  scale_color_manual(values = c("black" = "black", "white" = "white")) +
  theme_minimal(base_size = 14)
coef_hm

# prevalence
prevalence_df = tse_int.sgb.fil.pretreatment %>%
  transformAssay(method = "pa") %>% 
  meltSE(assay.type = "pa") %>%
  rename(SGB = FeatureID, Sample = SampleID) %>%
  inner_join(metadata) %>%
  group_by(SGB, Study, Response) %>%
  summarize(prevalence = mean(pa)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Response", values_from = "prevalence") %>% 
  mutate(difference = Responder - `Non-responder`)

prevalence_hm <- prevalence_df %>%
  inner_join(fit_lm_sgb.res.fil %>% select(SGB, Species, sign_pval, pval)) %>%
  filter(!is.na(Species)) %>%
  ggplot(aes(Study, Species)) +
  geom_tile(aes(fill = difference)) +
  scale_fill_gradient2(
    limits = c(-0.3, 0.3), # fix color scale range
    oob = scales::squish, # values outside [-5, 5] are mapped to extremes
    # low = "blue", mid = "white", high = "red", midpoint = 0
  ) +
  labs(
    fill = "Prevalence\ndifference",
    subtitle = "Prevalence"
  ) +
  theme_minimal(base_size = 14)

# Meta-analysis.
patch <- (barplot_pval +
  theme(
    axis.text.y = element_text(size = 10)
  )) +
  (coef_hm +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    guides(
      fill = guide_colorbar(
        label.theme = element_text(angle = 45),
        label.vjust = 0.5
      )
    )) +
  plot_layout(guides = "collect", nrow = 1, widths = c(1/2, 1/2)) +
  plot_annotation(
    title = "Antagonistic clusters DRUP-MSI (fastspar)",
    subtitle = "Meta-analysis with MMUPhin adjusting for ATB use and cancer type."
  ) &
  theme(legend.position = "bottom") 
patch

# Figure 1.E
ggsave(
  "output/figures/figures_manuscript/consortia_ranking_microbeICB_clusters_interest.pdf",
  patch,
  height = 16,
  width = 8
)



# Clusters response association -------------------------------------------------


metadata_hm <-  metadata_icbm_hm 

sgb_cluster_2 <- metadata_hm %>%
  filter(cluster_small %in% c(2)) %>%
  pull(SGB)

sgb_cluster_1 <- metadata_hm %>%
  filter(cluster_small %in% c(1)) %>%
  pull(SGB)

sgb_cluster_3 <- metadata_hm %>%
  filter(cluster_small %in% c(3)) %>%
  pull(SGB)

sgb_cluster_4 <- metadata_hm %>%
  filter(cluster_small %in% c(4)) %>%
  pull(SGB)

sgb_cluster_5 <- metadata_hm %>%
  filter(cluster_small %in% c(5)) %>%
  pull(SGB)


# Compute community summed abundance, cluster prevalence and community scores. 

colData(tse_int.sgb.fil.pretreatment) <- tse_int.sgb.fil.pretreatment %>%
  coldf() %>%
  mutate(
    abun_cluster_2 = get_abundance_set(
      tse_int.sgb.fil.pretreatment,
      sgb_cluster_2
    )$total_abun,
    preval_cluster_2 = get_abundance_set(
      tse_int.sgb.fil.pretreatment,
      sgb_cluster_2
    )$prevalence,
    abun_cluster_1 = get_abundance_set(
      tse_int.sgb.fil.pretreatment,
      sgb_cluster_1
    )$total_abun,
    preval_cluster_1 = get_abundance_set(
      tse_int.sgb.fil.pretreatment,
      sgb_cluster_1
    )$prevalence,
        abun_cluster_3 = get_abundance_set(
      tse_int.sgb.fil.pretreatment,
      sgb_cluster_3
    )$total_abun,
    preval_cluster_3 = get_abundance_set(
      tse_int.sgb.fil.pretreatment,
      sgb_cluster_3
    )$prevalence,
        abun_cluster_4 = get_abundance_set(
      tse_int.sgb.fil.pretreatment,
      sgb_cluster_4
    )$total_abun,
    preval_cluster_4 = get_abundance_set(
      tse_int.sgb.fil.pretreatment,
      sgb_cluster_4
    )$prevalence,
        abun_cluster_5 = get_abundance_set(
      tse_int.sgb.fil.pretreatment,
      sgb_cluster_5
    )$total_abun,
    preval_cluster_5 = get_abundance_set(
      tse_int.sgb.fil.pretreatment,
      sgb_cluster_5
    )$prevalence
  ) %>%
  mutate(
    ratio_abun = log10(abun_cluster_2 + 1e-6) - log10(abun_cluster_1 + 1e-6),
    ratio_prevalence = preval_cluster_2 - preval_cluster_1
  ) %>%
  DataFrame()

colData(tse_int.sgb.fil) <- tse_int.sgb.fil %>%
  coldf() %>%
  mutate(
    abun_cluster_2 = get_abundance_set(
      tse_int.sgb.fil,
      sgb_cluster_2
    )$total_abun,
    preval_cluster_2 = get_abundance_set(
      tse_int.sgb.fil,
      sgb_cluster_2
    )$prevalence,
    abun_cluster_1 = get_abundance_set(
      tse_int.sgb.fil,
      sgb_cluster_1
    )$total_abun,
    preval_cluster_1 = get_abundance_set(
      tse_int.sgb.fil,
      sgb_cluster_1
    )$prevalence,
        abun_cluster_3 = get_abundance_set(
      tse_int.sgb.fil,
      sgb_cluster_3
    )$total_abun,
    preval_cluster_3 = get_abundance_set(
      tse_int.sgb.fil,
      sgb_cluster_3
    )$prevalence,
        abun_cluster_4 = get_abundance_set(
      tse_int.sgb.fil,
      sgb_cluster_4
    )$total_abun,
    preval_cluster_4 = get_abundance_set(
      tse_int.sgb.fil,
      sgb_cluster_4
    )$prevalence,
        abun_cluster_5 = get_abundance_set(
      tse_int.sgb.fil,
      sgb_cluster_5
    )$total_abun,
    preval_cluster_5 = get_abundance_set(
      tse_int.sgb.fil,
      sgb_cluster_5
    )$prevalence
  ) %>%
  mutate(
    ratio_abun = log10(abun_cluster_2 + 1e-6) - log10(abun_cluster_1 + 1e-6),
    ratio_prevalence = preval_cluster_2 - preval_cluster_1
  ) %>%
  DataFrame()


correlation_stats <- tse_int.sgb.fil %>%
  coltibble() %>%
  split(~Study) %>% 
  map_dfr(function(df) {
    stats <- cor.test(df$preval_cluster_2, df$preval_cluster_1)
    tibble(
      r_val = round(stats$estimate, 2),
      p_val = stats$p.value,
      # Create label for plot
      label = paste0("r = ", round(stats$estimate, 2), 
                     ", p = ", format.pval(stats$p.value, digits = 2))
    )
  }, .id = "Study")


# Correlation Cluster 1 & Cluster 2
correlation_clusters <- tse_int.sgb.fil %>%
  coltibble() %>%
  ggplot(aes(preval_cluster_2, preval_cluster_1)) +
  geom_point(size = 1, col = "#378e60ff", show.legend = F) +
  geom_text(data = correlation_stats, 
            aes(x = 0.5, y = Inf, label = label),
            vjust = 1.5, size = 3.5, inherit.aes = FALSE) +
  labs(
    x = "Fraction of taxa from cluster 2 detected",
    y = "Fraction of taxa from cluster 1 detected",
    color = "Community score"
  ) +
  facet_wrap(~Study, ncol = 4) +
  scale_y_continuous(expand = expansion(c(0.05, 0.07))) + 
  theme_bw(base_size = 14) +
  theme(
    strip.background = theme_few()$strip.background,
    strip.text = element_text(size = 14), 
    panel.grid = element_blank()
  )

correlation_clusters

response_enrichment_clusters <- metadata_icbm_hm %>%
  mutate(
    sig = case_when(
      pval < 0.05 & coef > 0 ~ "R associated (p<0.05)",
      pval < 0.1 & pval > 0.05 & coef > 0 ~ "R associated (p<0.1)",
      pval < 0.05 & coef < 0 ~ "NR associated (p<0.05)",
      pval < 0.1 & pval > 0.05 & coef < 0 ~ "NR associated (p<0.1)",
      .default = NA
    )
  ) %>%
  mutate(
    cluster = cluster_small,
    sig = sig %>%
      fct_na_value_to_level(level = "Non significant") %>%
      fct_relevel(c(
        "R associated (p<0.05)",
        "R associated (p<0.1)",
        "Non significant",
        "NR associated (p<0.1)",
        "NR associated (p<0.05)"
      )) %>% 
      fct_rev()
  ) %>%
  ggplot(aes(cluster, fill = sig)) +
  geom_bar(position = "fill") +
  geom_text(
    aes(label = signif(..count..)),
    stat = "count",
    position = position_fill(vjust = 0.5),
    color = "white",
    size = 5.5
  ) +
  labs(y = "Proportion", x = "Bacterial cluster", fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel = element_blank()
  ) +
  scale_fill_manual(
    values = c(
      "Non significant" = "gray50",
      "R associated (p<0.05)" = "#00BFC4",
      "NR associated (p<0.05)" = "#F8766D",
      "R associated (p<0.1)" = "#007478ff",
      "NR associated (p<0.1)" = "#ba5952ff"
    )
  ) 

final_plot <- correlation_clusters / 
  ((response_enrichment_clusters | patchwork::plot_spacer()) + plot_layout( widths = c(4, 1))) +
  plot_layout(heights = c(3, 3))

final_plot

## Figure 1.C & Figure 1.D ## 
ggsave(
  filename = "output/figures/figures_manuscript/correlation_prevalence_proportion_response.pdf",
  final_plot,
  width = 9,
  height = 8.5
)


# Ratio Cluster 1 Cluster 2 response association ----------------------------------

pval_ratio <- tse_int.sgb.fil.pretreatment %>%
  coltibble() %>%
  mutate(abun = ifelse(ratio_abun > 0, "Higher\nCluster 2", "Higher\nCluster 1")) %>%
  tabyl(abun, Response, Study) %>%
  # lapply(fisher.test)
  map_dfr(
    function(data) {
      df = tidy(fisher.test(data))
      return(df)
    },
    .id = "Study"
  ) %>%
  mutate(
    p_label = sprintf("p = %.3f", p.value),
    group1 = "Higher\nCluster 2",
    group2 = "Higher\nCluster 1"
  )

plot_fisher <- tse_int.sgb.fil.pretreatment %>%
  coltibble() %>%
  # mutate(abun = ifelse(ratio_abun > -3, "neg", "pos")) %>%
  mutate(abun = ifelse(ratio_abun > 0, "Higher\nCluster 2", "Higher\nCluster 1")) %>%
  ggplot(aes(abun, fill = Response)) +
  geom_bar(position = "fill") +
  geom_text(
    aes(label = signif(..count..)),
    stat = "count",
    position = position_fill(vjust = 0.5),
    color = "white",
    size = 5
  ) +
  geom_signif(
    data = pval_ratio,
    aes(xmin = group1, xmax = group2, y_position = 1.03, annotations = p_label),
    manual = TRUE,
    tip_length = 0.03,
    inherit.aes = FALSE,
    # map_signif_level = function(p) sprintf("p = %.2g", p) # Not working with manual.
  ) +
  facet_wrap(~Study) +
  labs(
    x = NULL,
    y = "Proportion"
  ) +
    scale_y_continuous(expand = expansion(c(0.05,0.1))) +
  theme_minimal(base_size = 14)

plot_fisher + facet_wrap(~Study, nrow = 1)

## Figure S4.D ##
ggsave(
  "output/figures/figures_manuscript/fisher_test_ratio_long.pdf",
  plot_fisher + facet_wrap(~Study, nrow = 1),
  height = 2.7,
  width = 15
)






# Commmunity score ------------------------------------------------------------

## ROC plot -------------------------------------------------------------------

auc_ratio_prevalence <- roc_from_col(
  tse_int.sgb.fil.pretreatment,
  "ratio_prevalence",
  direction_roc = ">",
  batch_col = "Study",
  base_size = 14,
  make_plot = F
)$auc_table %>% 
  mutate( model = "prevalence")

mean(auc_ratio_prevalence$auc) # median = 0.58
median(auc_ratio_prevalence$auc) #median = 0.561
sd(auc_ratio_prevalence$auc) # sd = 0.04

auc_ratio_abun <- roc_from_col(
  tse_int.sgb.fil.pretreatment,
  "ratio_abun",
  direction_roc = ">",
  batch_col = "Study",
  base_size = 14,
  make_plot = F
)$auc_table %>% 
  mutate( model = "abundance")

mean(auc_ratio_abun$auc) #mean = 0.6
median(auc_ratio_abun$auc) #median = 0.59
sd(auc_ratio_abun$auc) #sd = 0.04

# Save roc plot:
roc_plot <- roc_from_col(
  tse_int.sgb.fil.pretreatment,
  "ratio_abun",
  direction_roc = ">",
  batch_col = "Study",
  base_size = 14,
  make_plot = T,
  average = T
)$plot

roc_plot <- roc_plot +
  labs(
    title = "Community Score (Abundance)",
    subtitle = paste0(
      "All SGBs.",
      sprintf(
        " Mean AUC=%.3f±%.3f",
        mean(auc_ratio_abun$auc),
        sd(auc_ratio_abun$auc)
      )
    ),
    color = "Study"
  ) +
  scale_color_manual(
    values = colors_studies
  )

roc_plot

## Figure 2.A ## 
ggsave(
  "output/figures/figures_manuscript/roc_plot_community_score_abundance.pdf",
  roc_plot,
  height = 7,
  width = 6.5
)

## ggridges --------------------------------------------------------------------

# pvalues community score (abundance)
studies <- unique(colData(tse_int.sgb.fil.pretreatment)$Study)
perm_fit_res <- map_dfr(
  studies %>% purrr::set_names(),
  function(study) {

    df_study <- tse_int.sgb.fil.pretreatment %>%
      coltibble() %>%
      filter(Study == study) %>%
      mutate(bmi = as.numeric(bmi)) %>%
      replace_na(list(bmi = 0, atb = "Unknown", host_disease_status = "Unknown"))

    # What variables to keep in the model. 
    valid_covariates <-
      c("age", "gender", "atb", "batch", "host_disease_status") %>%
      keep(~ {
        # Check if column exists and has more than 1 unique value
        val_count <- df_study[[.x]] %>% na.omit() %>% unique() %>% length()
        val_count > 1
      })

    formula_str <- paste("ratio_abun ~ Response +", paste(valid_covariates, collapse = " + "))
    
    model <- permuco::lmperm(
      as.formula(formula_str),
      data = df_study,
      np = 5000
    )

    sum <- summary(model) %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "term") %>%
      mutate(formula_used = formula_str)
    
    return(sum)
  },
  .id = "Study", .progress = TRUE
)

# host_disease_status, atb, batch
res_studies <- perm_fit_res %>%
  select(-c(`resampled Pr(<t)`, `resampled Pr(>t)`)) %>%
  mutate(
    p.value = `resampled Pr(>|t|)`,
    p = sprintf("p = %.2g", p.value),
    .after = 5
  )

res_studies %>% 
      filter(term == "ResponseResponder")

# Meta-analysis. 
m.gen <- meta::metagen(
  TE = Estimate,              # Treatment Effect (the Estimate from lmperm)
  seTE = `Std. Error`,        # Standard Error from lmperm
  studlab = Study,            # The column with your cohort names
  data = res_studies %>% filter(term == "ResponseResponder"),
  sm = "MD",                  # Mean Difference
  fixed = FALSE,              # Turn off fixed effect
  random = TRUE,              # Use random effects (handles heterogeneity)
  method.tau = "REML"         # Restricted Maximum Likelihood
)

summary(m.gen)

# Abundance
abundance_ggridge <- tse_int.sgb.fil.pretreatment %>%
  coltibble() %>%
  mutate(median_s_score = median(ratio_abun), .by = Study) %>%
  ggplot(aes(x = ratio_abun, y = Study)) +
  ggridges::geom_density_ridges(
    aes(fill = Response),
    quantile_lines = T,
    quantiles = 2,
    quantile_fun = mean,
    vline_linetype = "dashed",
    scale = 1.4, # Reduces height so they don't overlap the next study
    alpha = 0.75
  ) +
  geom_vline(xintercept = 0) +
  geom_text(
    data = res_studies %>%
      filter(term == "ResponseResponder") %>%
      mutate(signif = ifelse(p.value < 0.1, "yes", "no")),
    aes(label = p),
    x = 5,
    show.legend = F,
    nudge_y = 0.4
  ) +
  labs(x = "Community score (abundance)") +
  # xlim(0, 1) +
  theme_bw(base_size = 17) 

abundance_ggridge

ggsave(
  "output/figures/figures_manuscript/community_score_ggridges.pdf",
  abundance_ggridge,
  width = 10,
  height = 6
)

# Forest plot.
col_vec <- colors_studies[m.gen$studlab]
pdf("output/figures/figures_manuscript/community_score_forest_plot.pdf", width = 10, height = 6)  # was 5, now 6
meta::forest(
  m.gen,
  sortvar = TE,
  xlab = "Mean Difference",
  xlim = c(-2, 1),
  leftcols = c("studlab", "w.random"),
  leftlabs = c("Cohort", "Weight"),
  rightcols = c("pval"),
  rightlabs = c("p-value"),
  print.tau2 = TRUE,
  print.I2 = TRUE,
  print.pval.Q = FALSE,
  test.overall.random = TRUE, # adds p-value to the RE summary line
  # col.diamond = "navy",
  # col.square = "gray40",
  col.square = col_vec,
  col.square.lines = col_vec,
  col.study = col_vec,
  comb.random = TRUE,
  addrows.below.overall = 2
)
dev.off()

## Survival ---------------------------------------------------------------------

data_survival <- tse_int.sgb.fil.pretreatment %>%
  coltibble() %>%
  mutate(
    os_time = ifelse(Study == "Derosa", os_time * 30, os_time)
  ) %>%
  mutate(
    os_event = case_when(
      os_event == "yes" ~ "1",
      os_event == "no" ~ "0",
      .default = os_event
    ),
    pfs_event = case_when(
      pfs_event == "yes" ~ "1",
      pfs_event == "no" ~ "0",
      .default = pfs_event
    )
  ) %>%
  mutate(
    ratio_abun = scale(ratio_abun)[, 1],
    .by = Study
  )
  
data_survival_os <- data_survival %>%
  filter(
    !is.na(os_event) & !is.na(os_time)
  ) %>%
  mutate(
    os_event = as.numeric(os_event),
    os_time = as.numeric(os_time)
  )

data_survival_pfs <- data_survival %>%
  filter(
    !is.na(pfs_event) | !is.na(pfs_time)
  ) %>%
  mutate(
    pfs_event = as.numeric(pfs_event),
    pfs_time = as.numeric(pfs_time)
  )

# Survival
community_score_surv <- data_survival_os %>%
  nest_by(Study) %>%
  summarise(
    res = list(
      broom::tidy(
        coxph(
          Surv(os_time, os_event) ~ ratio_abun +
            host_disease_status +
            atb + 
            age +
            gender,
          data = data
        ),
        exponentiate = TRUE,
        conf.int = TRUE
      )
    )
  ) %>%
  unnest() %>%
  ungroup() %>%
  filter(term %in% c("ratio_abun")) %>%
  mutate(
    label = sprintf("HR=%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
    p_label = ifelse(p.value < 0.001, "p<0.001", sprintf("p=%.3f", p.value))
  ) %>%
  mutate(model = "OS")

community_score_surv

# PFS
community_score_surv_pfs <- data_survival_pfs %>%
  nest_by(Study) %>%
  summarise(
    res = list(
      broom::tidy(
        coxph(
          Surv(pfs_time, pfs_event) ~ ratio_abun +
            host_disease_status +
            atb + 
            age +
            gender,
          data = data
        ),
        exponentiate = TRUE,
        conf.int = TRUE
      )
    )
  ) %>%
  unnest() %>%
  ungroup() %>%
  filter(term %in% c("ratio_abun")) %>%
  mutate(
    label = sprintf("HR=%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
    p_label = ifelse(p.value < 0.001, "p<0.001", sprintf("p=%.3f", p.value))
  ) %>%
  mutate(model = "PFS")

community_score_surv_pfs

# Plot both together

community_score_surv_comb <- bind_rows(
  community_score_surv_pfs,
  community_score_surv
) 

forest_community_score <- community_score_surv_comb %>%
  mutate(
    Study = as.factor(Study) %>%
      fct_relevel(c("Derosa", "Spencer", "DRUP", "Gunjur"))
  ) %>%
  ggplot(aes(y = Study, x = estimate)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "gray50",
    linewidth = 0.8
  ) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0,
    linewidth = 1
  ) +
  geom_point(size = 7.5, shape = 18, show.legend = F) +
  geom_point(size = 6, shape = 18, aes(color = Study), show.legend = F) +
  # Add HR and CI labels
  geom_text(
    aes(x = 0.4, label = label),
    hjust = 0,
    size = 4,
    show.legend = FALSE
  ) +
  # Add p-values
  geom_text(
    aes(x = 0.4, label = p_label),
    hjust = 0,
    vjust = 2.5,
    size = 4,
    show.legend = FALSE
  ) +
  scale_x_log10() +
  scale_color_manual(
    values = colors_studies
  ) +
  labs(
    x = "Adjusted Hazard Ratio (95% CI)",
    y = NULL,
    color = "Model",
    title = "Community Score association with survival"
  ) +
  facet_grid(model ~ ., scales = "free", space = "free_y", switch = "y") +
  coord_cartesian(
    xlim = c(0.4, 2.8)
  ) +
  scale_y_discrete(position = "right") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    strip.text.y.left = element_text(angle = 0)
  ) +
  theme_bw(base_size = 14)

forest_community_score

ggsave(
  "output/figures/figures_manuscript/community_score_survival.pdf",
  forest_community_score,
  width = 6.5,
  height = 5
)






# ------------------------------------------------------------------------------#
# Tax phylotree ----------------------------------------------------------------
# ------------------------------------------------------------------------------#
# Figure S2.A

library(ggtree)
library(ggnewscale)
library(ape)

mpa_tree = ape::read.tree("data/resources/mpa_vOct22_CHOCOPhlAnSGB_202212.nwk")
mpa_tree$tip.label = paste0("SGB", mpa_tree$tip.label)

sgb_cluster_2
sgb_cluster_1

rowData(tse_int.sgb)$cluster <- case_when(
  rownames(tse_int.sgb) %in% sgb_cluster_2 ~ "2",
  rownames(tse_int.sgb) %in% sgb_cluster_1 ~ "1",
  .default = NA
)

tse = tse_int.sgb

# library(bugphyzz)
# library(taxize)
# bp <- importBugphyzz()

### Make tree based on taxonomy annotation

tree_3 <- rowtibble(tse) %>% 
  filter(!is.na(cluster)) %>% 
  select(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SGB")) %>% 
  toTree()


tree_3$tip.label <- rowtibble(tse) %>%
  filter(!is.na(cluster)) %>%
  mutate(SGB = fct_relevel(SGB, tree_3$tip.label)) %>% # Order the SGBs in the same way as in the tree
  arrange(SGB) %>%
  mutate(label = paste(Species, SGB, sep = " ")) %>%
  pull(label)

phylum.anno <- rowtibble(tse) %>%
  mutate(label = paste(Species, SGB, sep = " "), label2 = label) %>%
  filter(label %in% tree_3$tip.label) %>% # Filter species in the tree
  column_to_rownames(var = "label2") %>%
  mutate(node = label) %>%
  select(label, SGB, cluster, Phylum) #! Make sure the label column is the first column 
rownames(phylum.anno) <- rownames(phylum.anno)

p <- ggtree(
  tree_3,
  # layout = "radial",
  layout = "rectangular",
  root.position = 45,
  open.angle = 0,
  size = 0.30
) %<+%
  phylum.anno


p1 <-
  # rotate_tree(p, -90) +
  p +
  ylim(-5, NA) +
  geom_tiplab(size = 1.7, offset = 0.6, align = F)

# Load results meta-analysis. 
metadata_icbm_hm

anno <- rowtibble(tse) %>%
  filter(!is.na(cluster)) %>% 
  mutate(label = paste(Species, SGB, sep = " ")) %>% 
  left_join(
    metadata_icbm_hm %>% 
      select(SGB, sig)
  ) %>% 
  column_to_rownames(var = "label") %>%
  select(Phylum, cluster, sig)

colors_phylum <- c(
  "Actinobacteria" = "#7FC97F", 
  "Firmicutes" = "#BEAED4",
  "Proteobacteria" = "#FDC086" 

)

p2 = gheatmap(
  p1,
  anno,
  offset = 0,
  width = 0.5,
  colnames_angle = 0,
  colnames_offset_y = 4.5,
  font.size = 2,
  colnames = T,
  colnames_position = 'top'
) +
  scale_fill_manual(
    values = c(
      "NR_associated" = "red",
      "R_associated" = "blue",
      colors_clusters,
      colors_phylum
    ),
    na.value = "white"
  )

# Add association with respnose (signed pvalue)

df <- rowtibble(tse) %>%
  filter(!is.na(cluster)) %>%
  mutate(label = paste(Species, SGB, sep = " ")) %>%
  left_join(
    metadata_icbm_hm %>%
      select(SGB, sig, sign_pval)
  ) %>%
  column_to_rownames(var = "label") %>%
  select(sign_pval)

p3 <- p2 + new_scale_fill()
p3 <- gheatmap(
  p3,
  df,
  offset = 0.35,
  width = 0.18,
  colnames_angle = 0,
  colnames_offset_y = 4.5,
  # colnames_offset_x = 5,
  font.size = 2,
  colnames = T,
  colnames_position = 'top'
) +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    limits = c(-4, 4),
    name = "-log10(p-value)"
  )

p3 = p3 +
  theme(
    plot.margin = margin(50, 50, 50, 50), # Increase distance so labels are not cut. 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.position = "right"
  )
p3

## Figure S2.A ##
ggsave(
  'output/figures/figures_manuscript/phylotree_drup_clusters_fastspar_v2.pdf',
  p3,
  width = 10,
  height = 10,
  bg = "transparent"
)

# ------------------------------------------------------------------------------
# Healthy Cohorts --------------------------------------------------------------
# ------------------------------------------------------------------------------
# Figure S2.B

# Community score and Cluster 1 & 2 in healthy cohrots. 

# Load data
tse_healthy500fg <- read_rds("output/data/metagenomics/healthy/tse_healthy500fg.rds") 
tse_DAG3 <- read_rds("output/data/metagenomics/DMP/DAG3_tse.rds")

# Compute community score
colData(tse_healthy500fg) <- tse_healthy500fg %>%
  coldf() %>%
  mutate(
    abun_cluster_2 = get_abundance_set(
      tse_healthy500fg,
      sgb_cluster_2
    )$total_abun,
    preval_cluster_2 = get_abundance_set(
      tse_healthy500fg,
      sgb_cluster_2
    )$prevalence,
    abun_cluster_1 = get_abundance_set(
      tse_healthy500fg,
      sgb_cluster_1
    )$total_abun,
    preval_cluster_1 = get_abundance_set(
      tse_healthy500fg,
      sgb_cluster_1
    )$prevalence,
        abun_cluster_3 = get_abundance_set(
      tse_healthy500fg,
      sgb_cluster_3
    )$total_abun,
    preval_cluster_3 = get_abundance_set(
      tse_healthy500fg,
      sgb_cluster_3
    )$prevalence,
        abun_cluster_4 = get_abundance_set(
      tse_healthy500fg,
      sgb_cluster_4
    )$total_abun,
    preval_cluster_4 = get_abundance_set(
      tse_healthy500fg,
      sgb_cluster_4
    )$prevalence,
        abun_cluster_5 = get_abundance_set(
      tse_healthy500fg,
      sgb_cluster_5
    )$total_abun,
    preval_cluster_5 = get_abundance_set(
      tse_healthy500fg,
      sgb_cluster_5
    )$prevalence
  ) %>%
  mutate(
    ratio_abun = log10(abun_cluster_2 + 1e-6) - log10(abun_cluster_1 + 1e-6),
    ratio_prevalence = preval_cluster_2 - preval_cluster_1
  ) %>%
  DataFrame()

colData(tse_DAG3) <- tse_DAG3 %>%
  coldf() %>%
  mutate(
    abun_cluster_2 = get_abundance_set(
      tse_DAG3,
      sgb_cluster_2
    )$total_abun,
    preval_cluster_2 = get_abundance_set(
      tse_DAG3,
      sgb_cluster_2
    )$prevalence,
    abun_cluster_1 = get_abundance_set(
      tse_DAG3,
      sgb_cluster_1
    )$total_abun,
    preval_cluster_1 = get_abundance_set(
      tse_DAG3,
      sgb_cluster_1
    )$prevalence,
        abun_cluster_3 = get_abundance_set(
      tse_DAG3,
      sgb_cluster_3
    )$total_abun,
    preval_cluster_3 = get_abundance_set(
      tse_DAG3,
      sgb_cluster_3
    )$prevalence,
        abun_cluster_4 = get_abundance_set(
      tse_DAG3,
      sgb_cluster_4
    )$total_abun,
    preval_cluster_4 = get_abundance_set(
      tse_DAG3,
      sgb_cluster_4
    )$prevalence,
        abun_cluster_5 = get_abundance_set(
      tse_DAG3,
      sgb_cluster_5
    )$total_abun,
    preval_cluster_5 = get_abundance_set(
      tse_DAG3,
      sgb_cluster_5
    )$prevalence
  ) %>%
  mutate(
    ratio_abun = log10(abun_cluster_2 + 1e-6) - log10(abun_cluster_1 + 1e-6),
    ratio_prevalence = preval_cluster_2 - preval_cluster_1
  ) %>%
  DataFrame()


# Generate plots
correlation_healthy <- tse_DAG3 %>%
  coltibble() %>%
  mutate(Study = "DAG3") %>%
  bind_rows(
    coltibble(tse_healthy500fg) %>% mutate(Study = "500FG"),
  ) %>% 
  ggplot(aes(preval_cluster_2, preval_cluster_1)) +
  geom_jitter(size = 1, col = "#378e60ff", show.legend = F) +
  stat_cor(label.x = 0.2, size = 3) +
  labs(
    x = "Fraction of taxa from cluster 2 detected",
    y = "Fraction of taxa from cluster 1 detected",
    color = "Community\nscore"
  ) +
  scale_color_viridis() +
  facet_wrap(~Study, ncol = 1) +
  theme_bw(base_size = 15) +
  theme(
    # strip.text = theme_few()$strip.text,
    strip.background = theme_few()$strip.background
  )


community_prevalence_distribution_healthy <- tse_DAG3 %>%
  coltibble() %>%
  mutate(Study = "DAG3") %>%
  bind_rows(
    coltibble(tse_healthy500fg) %>% mutate(Study = "500FG"),
    coltibble(tse_int.sgb.fil) %>% select(-bmi)
  ) %>%
  mutate(mean_preval_score = mean(ratio_prevalence),.by = Study) %>% 
  mutate(Study = fct_reorder(Study, -mean_preval_score)) %>%
  ggplot(aes(x = ratio_prevalence, y = Study)) +
  ggridges::geom_density_ridges(
    aes(fill = Study),
    quantile_lines = T,
    quantiles = 2,
    quantile_fun = mean,
    vline_linetype = "dashed",
    scale = 1.4, # Reduces height so they don't overlap the next study
    alpha = 0.65
  ) +
  geom_vline(xintercept = 0) +
  # scale_x_log10() +
  labs(x = "Community score (prevalence)") +
  scale_fill_manual(
    values = c(
      colors_studies,
      "500FG" = "#ABD9E9",
      "DAG3" = "#ABD9E9"
    )
  ) +
  # xlim(0, 1) +
  theme_bw(base_size = 15)

patch <- correlation_healthy + community_prevalence_distribution_healthy + plot_layout(widths = c(1,2))

patch

## Figure S2.B
ggsave(
  filename = "output/figures/figures_manuscript/healthy_cohorts.pdf",
  patch,
  width = 13,
  height = 6
)

# ------------------------------------------------------------------------------
# Oral bacteria ----------------------------------------------------------------
# ------------------------------------------------------------------------------
# Figure S3


### Check cluster with Oral taxa. ---------------------------------------------

metadata_hm %>% 
  tabyl(oral, cluster_small)

# Oral bacteria.
oral_proportion <- metadata_hm %>%
  mutate(
    cluster = paste("Cluster", cluster_small),
    oral = ifelse(oral == "Oral", "Oral bacteria", oral) %>%
      # fct_relevel("Oral bacteria") %>%
      fct_na_value_to_level(level = "Other") %>% 
      fct_rev()
  ) %>%
  ggplot(aes(cluster, fill = oral)) +
  geom_bar(position = "fill") +
  geom_text(
    aes(label = signif(..count..)),
    stat = "count",
    position = position_fill(vjust = 0.5),
    color = "white",
    size = 6
  ) +
  # facet_wrap(~Study) +
  # coord_flip() +
  labs(y = "Proportion", x = NULL, fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel = element_blank(),
    # legend.position = "bottom"
  ) +
  scale_fill_manual(
    values = c("Other" = "gray50", "Oral bacteria" = "#F8766D")
  ) +
  ggpubr::rotate_x_text(angle = 45)

tile_cluster <- metadata_hm %>% 
  distinct(cluster_small) %>% 
  mutate(cluster_small = as.factor(cluster_small)) %>% 
  ggplot(aes(x = cluster_small, y = "tile")) +
  geom_tile(aes(fill = cluster_small), show.legend = F) +
  scale_fill_manual(values = colors_clusters) +
  theme_void()

fit_lm_sgb.res = fit_lm_sgb$meta_fits %>% 
  mutate(fdr = p.adjust(pval, "fdr")) %>% 
  mutate(val = -log10(pval) * sign(coef)) %>% 
  mutate(SGB = feature) %>% 
  as_tibble() %>% 
  left_join(rowtibble(tse_int.sgb.fil.pretreatment)) %>% 
  relocate(Species) %>% 
  arrange (pval)

fit_lm_sgb.res.fil <- fit_lm_sgb.res %>%
  left_join (metadata_icbm_hm %>% select(SGB, cluster_small)) %>% 
  # filter(pval < 0.1) %>% 
  filter(cluster_small == 3) %>% 
  # filter(oral == "Oral") %>% 
  mutate(n_sgb = n(), .by = Species) %>%
  mutate(
    # Species = ifelse(n_sgb == 2, paste(Species, SGB, sep = "_"), Species)
    Species = ifelse(
      str_detect(Species, SGB),
      Species,
      paste(Species, SGB, sep = "_")
    )
  ) %>%
  mutate(
    Cluster = cluster_small,
    sign_pval = -log10(pval) * sign(coef),
    Species = fct_reorder(Species, sign_pval),
    Association = ifelse(coef > 0, "Beneficial", "Detrimental")
  ) %>% 
  mutate(Cluster = ifelse(Cluster %in% c(3, 4), Cluster, "Other"))

barplot_pval <- fit_lm_sgb.res.fil %>%
  ggplot(aes(x = sign_pval, y = Species)) +
  # geom_col(aes(fill = SIG)) +
  geom_col(aes(fill=oral), show.legend = F) +
  geom_vline(xintercept = c(log10(0.05), -log10(0.05)), linetype = "dashed") +
  labs(
    x = "Signed p-value meta-analysis",
    y = NULL
  ) +
  theme_minimal(base_size = 14)
barplot_pval

# Significant species in meta-analysis.
coef_hm <- map_df(names(fit_lm_sgb$maaslin_fits), function(study) {
  res <- fit_lm_sgb$maaslin_fits[[study]]
  res$Study <- study
  res
}) %>%
  mutate(SGB = feature) %>%
  inner_join(fit_lm_sgb.res.fil %>% select(SGB, Species, sign_pval)) %>%
  mutate(
    pval.res.sign = -log10(pval) * sign(coef),
    # pval.res.sign = -log10(pval.res) * sign(coef.res),
    Species = fct_reorder(Species, sign_pval)
  ) %>%
  # filter(pval.res < 0.05) %>%
  ggplot(aes(Study, Species)) +
  geom_tile(aes(fill = coef)) +
  geom_text(
    aes(
      label = ifelse(pval < 0.05, "x", ""),
      color = ifelse(abs(coef) > 0.013, "white", 'black')
    ),
    show.legend = F,
    # color = "black",
    size = 2.8
  ) +
  labs(
    subtitle = "Coefficient",
    fill = "Coefficient\nMMUPHin",
    y = NULL
  ) +
  scale_fill_gradient2(
    limits = c(-0.02, 0.02),
    oob = scales::squish,
  ) +
  scale_color_manual(values = c("black" = "black", "white" = "white")) +
  theme_minimal(base_size = 14)
coef_hm

# prevalence
prevalence_df = tse_int.sgb.fil.pretreatment %>%
  transformAssay(method = "pa") %>% 
  meltSE(assay.type = "pa") %>%
  rename(SGB = FeatureID, Sample = SampleID) %>%
  inner_join(metadata) %>%
  group_by(SGB, Study, Response) %>%
  summarize(prevalence = mean(pa)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Response", values_from = "prevalence") %>% 
  mutate(difference = Responder - `Non-responder`)

prevalence_hm <- prevalence_df %>% 
  inner_join(fit_lm_sgb.res.fil %>% select(SGB, Species, sign_pval)) %>% 
  filter(!is.na(Species)) %>% 
  ggplot(aes(Study, Species)) +
  geom_tile(aes(fill=difference)) +
  scale_fill_gradient2(
    limits = c(-0.3, 0.3),       # fix color scale range
    oob = scales::squish,    # values outside [-5, 5] are mapped to extremes
    # low = "blue", mid = "white", high = "red", midpoint = 0
  ) +
  labs(
    fill = "Prevalence\ndifference",
    subtitle = "Prevalence"
  ) +
  theme_minimal(base_size = 14)

# Meta-analysis.
patch <- (barplot_pval +
  theme(
    axis.text.y = element_text(size = 10)
  )) +
  (coef_hm +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )) +
  (prevalence_hm +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )) +
  plot_layout(guides = "collect", widths = c(2 / 10, 2 / 10, 2 / 10)) &
  theme(
    legend.position = 'bottom',
    legend.text = element_text(angle = 90, vjust = 0.5, hjust = 0.5)
  )

left_col <- wrap_plots(tile_cluster, oral_proportion, plot_spacer(), ncol = 1, heights = c(0.1,2,1.5))

patch2 <- (left_col | patch) +
  plot_layout( widths = c(1,1))

patch2

## Figure S3.A & Figure S3.B ##
ggsave(
  "output/figures/figures_manuscript/cluster_3_oral_bacteria.pdf",
  patch2,
  height = 10,
  width = 15
)

### PPI and oral bacteria --------------------------------------------------------

# Load data DRUP with ppi use: 
tse_reads <- read_rds("output/data/metagenomics/drup/tse_reads_metaphlan.rds")

# Load data Spencser with ppi use: 
spencer.metadata.ppi = read_csv(
  "../studies/spencer_science_2021/SraRunTable.csv"
) %>% 
  select(Run, PPI) %>% 
  rename(Sample = Run) %>% 
  mutate(PPI_use = ifelse(PPI == 1, "Yes", "No"))

# Load data Gunjur with ppi use: 
gunjur.metadata.ppi = read_xlsx(
  "../studies/gunjur_natmed_2024/41591_2024_2823_MOESM3_ESM.xlsx",
  sheet = 4
) %>%
  left_join(
    read_csv(
      "../studies/gunjur_natmed_2024/ena_metadata.csv"
    ) %>%
      mutate(lane_id = str_replace(run_alias, "SC_RUN_", "")) %>%
      relocate(lane_id)
  ) %>% 
  select(run_accession, PPI_use) %>% 
  rename(Sample = run_accession) %>% 
  mutate(PPI_use = ifelse(PPI_use == TRUE, "Yes", "No"))

# Join the PPI metadata. 
drup.metadata.ppi <- tse_reads %>%
  coltibble() %>%
  filter(PPI_use != "Unknown") %>% 
  select(Sample, PPI_use)

metadata.ppi <- bind_rows(
  drup.metadata.ppi,
  spencer.metadata.ppi,
  gunjur.metadata.ppi
)

ppi_use_all_abundance_clusters <- tse_int.sgb.fil %>%
  coltibble() %>%
  inner_join(metadata.ppi) %>%
  left_join(gunjur.metadata.ppi) %>%
  pivot_longer(
    cols = contains("abun_cluster"),
    names_to = "Cluster",
    values_to = "Community_score"
  ) %>%
  mutate(Community_score = Community_score + 1e-5) %>%
  mutate(
    cluster_name = Cluster %>% str_remove("abun_cluster_"),
    Cluster = str_replace(Cluster, "median_c", "C") %>% str_replace("_", " ")
  ) %>%
  ggplot(aes(x = PPI_use, Community_score)) +
  geom_boxplot(
    aes(fill = cluster_name),
    alpha = 0.5,
    show.legend = F,
    outliers = F
  ) +
  geom_quasirandom(size = 0.2) +
  facet_grid(Study ~ Cluster, scales = "free") +
  geom_signif(
    textsize = 3,
    comparisons = list(c("No", "Yes")),
    map_signif_level = function(p) {
      ifelse(
        p > 0.05,
        "ns",
        ifelse(
          p < 0.001,
          "p < 0.001",
          sprintf("p = %.3f", p)
        )
      )
    }
  ) +
  labs(
    x = "PPI use",
    y = "Cluster abundance"
  ) +
  scale_fill_manual(
    values = colors_clusters
  ) +
  scale_y_log10(expand = expansion(c(0.05, 0.15))) +
  theme_bw(base_size = 14) +
  theme(
    # strip.text = element_text(face = "bold", size = 12)
    strip.background = element_rect(fill = "white", color = "black")
  )

ppi_use_all_abundance_clusters

## Figure S3.C ##
ggsave(
  "output/figures/figures_manuscript/cluster_3_oral_bacteria_ppi_use.pdf",
  ppi_use_all_abundance_clusters,
  height = 5,
  width = 7
)

# ------------------------------------------------------------------------------
# Meta-analysis LOO ------------------------------------------------------------
# ------------------------------------------------------------------------------

# How do single significant species predict response to immunotherapy?

# Load meta-analysis results. 
meta_analysis_list <- read_rds("output/data/metagenomics/microbeICB/meta_analysis_LOO_SGB.rds")

Studies <- coltibble(tse_int.sgb.fil)$Study %>% unique()
aucs_individual <- map_dfr(
  Studies %>% set_names(Studies),
  function(batch) {
    df <- meta_analysis_list[[batch]]$res_meta %>%
      filter(pval < 0.05) %>%
      # mutate(auc = ifelse(coef < 0, 1 - auc, auc)) %>% 
      select(feature, auc, pval, coef) 
    return(df)
  },
  .id = "batch"
)

auc_ratio_abun <- roc_from_col(
  tse_int.sgb.fil.pretreatment,
  "ratio_abun",
  direction_roc = ">",
  batch_col = "Study",
  base_size = 14,
  make_plot = F
)$auc_table %>% 
  mutate( model = "abundance")

community_vs_single_species_boxplot <- aucs_individual %>%
  mutate(
    batch = fct_relevel(
      batch,
      auc_ratio_abun %>% arrange(-auc) %>% pull(batch)
    )
  ) %>%
  ggplot(aes(x = batch, y = auc)) +
  geom_quasirandom(size = 1, color = "gray") +
  geom_boxplot(alpha = 0.4, outliers = F, show.legend = F, aes(fill = batch)) +
  geom_point(
    data = auc_ratio_abun,
    color = "red",
    size = 3
  ) +
  labs(
    x = "Study"
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_fill_manual(values = colors_studies) +
  theme_bw() +
  labs(title = "Community vs Single bacteria", y = "ROC AUC") +
  ggpubr::rotate_x_text(angle = 45)

community_vs_single_species_boxplot

## Figure 2.B ##
ggsave(
  filename = "output/figures/figures_manuscript/community_vs_single_species_model.pdf",
  community_vs_single_species_boxplot,
  width = 4.5,
  height = 5.5
)

# Comparison with community score
community_vs_single_species_ranking <- auc_ratio_abun %>%
  mutate(feature = "Community score") %>% 
  bind_rows(aucs_individual) %>% 
  arrange(-auc) %>% 
  mutate(rank = row_number(), .by = batch) %>% 
  mutate(n = n(), mean_rank = median(rank), .by = feature) %>% 
  distinct(feature, .keep_all = T) %>% 
  ggplot(aes(x = mean_rank, y = n) ) +
  geom_label(data = . %>% filter(feature == "Community score"),
                   aes(label = feature),nudge_y = -0.5, nudge_x = 20) +
  geom_point() +
  geom_point(data = . %>% filter(feature == "Community score"), 
             color = "red") +
  labs(y = "# cohorts detected", x = "Median AUC rank") + 
  theme_bw() 

top_species <- bind_rows(aucs_individual) %>%
  arrange(-auc) %>%
  mutate(rank = row_number(), .by = batch) %>%
  mutate(
    n = n(),
    mean_rank = median(rank),
    mean_auc = median(auc),
    .by = feature
  ) %>%
  arrange(mean_rank) %>%
  filter(n == 7) %>%
  mutate(SGB = feature) %>%
  left_join(
    metadata_icbm_hm_small %>%
      select(SGB, cluster_small) %>%
      mutate(cluster = cluster_small %>% as.factor())
  ) %>%
  left_join(
    tse_int.sgb.fil %>%
      rowtibble() %>%
      select(SGB, Species)
  ) %>%
  mutate(
    Species_sgb = ifelse(str_detect(Species, SGB), Species, paste(Species, SGB, sep = "_"))
  ) %>%
  mutate(
    Species_sgb = Species_sgb %>%
      fct_reorder(mean_auc) %>%
      fct_relevel("Community score", after = Inf)
  ) %>%
  mutate(Association = ifelse(coef > 0, "Response", "Non-Response")) %>%
  ggplot(aes(auc, Species_sgb)) +
  geom_boxplot(aes(fill = cluster)) +
  geom_point(aes(color = batch)) +
  geom_boxplot(
    data = auc_ratio_abun %>%
      mutate(feature = "Community score", Species_sgb = "Community score"),
    fill = "grey85"
  ) +
  geom_point(
    data = auc_ratio_abun %>%
      mutate(feature = "Community score", Species_sgb = "Community score"),
    aes(color = batch)
  ) +
  labs(
    color = "Study",
    fill = "Cluster",
    x = "ROC AUC",
    y = NULL
  ) +
  scale_fill_manual(
    values = colors_clusters,
    na.value = "grey85",
    na.translate = FALSE
  ) +
  scale_color_manual(values = colors_studies) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dashed"
  ) +
  theme_few()


patch_median_rank <- 
  (community_vs_single_species_ranking + theme_bw(base_size = 15)) +
  (top_species + theme_few(base_size = 15)) +
  plot_layout(guides = "collect") 

patch_median_rank

## Figure S4.B & Figure S4.C ##
ggsave(
  filename = "output/figures/figures_manuscript/community_vs_single_species_model_top_species.pdf",
  patch_median_rank,
  width = 15,
  height = 5
)


# Plot AUC distribution of significant species:
# Boxplot :
boxplot_significant_species <- aucs_individual %>%
  filter(pval < 0.01) %>%
  group_by(feature) %>%
  summarise(n = n(), mean_auc = mean(auc), median_auc = median(auc)) %>%
  ungroup() %>%
  # filter(n > 1) %>%
  mutate(SGB = feature) %>%
  left_join(
    metadata_icbm_hm_small %>%
      select(SGB, cluster_small) %>%
      mutate(cluster = cluster_small %>% as.factor())
  ) %>%
  ggplot(aes(x = "", y = mean_auc)) +
  geom_boxplot(outliers = F, fill = "#E6F5C9") +
  geom_quasirandom(aes(color = cluster)) +
  geom_hline(
    yintercept = auc_ratio_abun$auc %>% mean(),
    color = "red",
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 0.5, # Places it at the far right
    y = mean(auc_ratio_abun$auc),
    label = "Community score",
    vjust = -1, # Shifts it up
    hjust = 0, # Aligns it to the right edge
    color = "red"
  ) +
  labs(
    x = "Significant species in\nthe LOSO meta-analysis",
    y = "Mean AUC across studies"
  ) +
  coord_cartesian(ylim = c(0.4, 0.63)) +
  scale_color_manual(
    values = colors_clusters,
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

boxplot_significant_species

## Figure 2.C ##
ggsave(
  filename = "output/figures/figures_manuscript/community_vs_single_species_model_boxplot.pdf",
  boxplot_significant_species,
  width = 4,
  height = 6
)


# ------------------------------------------------------------------------------
# Models performance comaparison -----------------------------------------------
# ------------------------------------------------------------------------------

# Compare aucs of all the models. 

aucs_clusters <- map_dfr(
  metadata_hm %>%
    distinct(cluster_small) %>%
    pull(cluster_small, name = cluster_small),
  function(cluster_name) {

    print(cluster_name)
    species_cluster <- metadata_hm %>%
      filter(cluster_small == cluster_name) %>%
      pull(SGB)

    tse = tse_int.sgb.fil.pretreatment
    colData(tse) <- tse %>%
      coldf() %>%
      mutate(
        abun_cluster = get_abundance_set(
          tse,
          species_cluster
        )$total_abun,
        preval_cluster = get_abundance_set(
          tse,
          species_cluster
        )$prevalence
      ) %>%
      DataFrame()

    roc_preval_cluster <- roc_from_col(
      tse,
      "preval_cluster",
      direction_roc = "<",
      batch_col = "Study",
      base_size = 14,
      make_plot = F
    )$auc_table %>%
      rename(auc_preval = auc)

    roc_abun_cluster <- roc_from_col(
      tse,
      "abun_cluster",
      direction_roc = "<",
      batch_col = "Study",
      base_size = 14,
      make_plot = F
    )$auc_table %>%
      rename(auc_abun = auc)

    mean_auc_preval <- median(roc_preval_cluster$auc_preval)
    mean_auc_abun <- median(roc_abun_cluster$auc_abun)

    roc_cluster <- left_join(
      roc_preval_cluster,
      roc_abun_cluster
    ) %>%
      mutate(
        auc_preval = if (mean_auc_preval < 0.5) 1 - auc_preval else auc_preval,
        auc_abun = if (mean_auc_abun < 0.5) 1 - auc_abun else auc_abun
      )
    
    print(roc_cluster)

    return(roc_cluster)
  },.id = "Cluster", .progress = T
)

aucs_clusters



clusters_interest <- aucs_clusters %>%
  # bind_rows(roc_toposcore$auc_table %>% mutate(model = "TopoScore")) %>%
  pivot_longer(cols = contains("auc"), names_to = "model", values_to = "auc") %>% 
  mutate(
    model_type = case_when(
      str_detect(model, "abun") ~ "Abundance",
      str_detect(model, "preval") ~ "Prevalence",
      .default = "S-score"
    )
  ) %>% 
  mutate(
      model = paste0("Cluster ", Cluster, "\n(", model_type, ")")
  )

auc_table_models <- bind_rows(
  auc_ratio_abun,
  auc_ratio_prevalence,
  auc_s_score
) %>%
  mutate(
    model = case_when(
      model == "prevalence" ~ "Comunity score \n (Prevalence)",
      model == "abundance" ~ "Community score \n (Abundance)",
      .default = model
    ),
    model_type = "Community models"
  )

data_barplot = auc_table_models %>%
  bind_rows(clusters_interest) %>%
  mutate(
    model = fct_relevel(
      model,
      "Cluster 1\n(Prevalence)",
      "Cluster 2\n(Prevalence)",
      "Cluster 3\n(Prevalence)",
      "Cluster 4\n(Prevalence)",
      "Cluster 5\n(Prevalence)",
      "Comunity score \n (Prevalence)",
      "Cluster 1\n(Abundance)",
      "Cluster 2\n(Abundance)",
      "Cluster 3\n(Abundance)",
      "Cluster 4\n(Abundance)",
      "Cluster 5\n(Abundance)",
      "Community score \n (Abundance)",
      "s_score"
    ),
    model_type = fct_relevel(
      model_type,
      "Prevalence",
      "Abundance",
      "Community models"
    )
  ) %>%
  filter(model != "s_score")

# Geom summary.
auc_summary <- data_barplot %>%
  group_by(model) %>%
  summarise(
    mean_auc = mean(auc),
    sd_auc = sd(auc),
    n = n(),
    se_auc = sd_auc / sqrt(n)
  ) %>%
  arrange(model) %>%
  mutate(
    model_name = gsub("\\n", "", model),
    model_name = gsub("\\s+", " ", model_name),
    model_name = sprintf("%s\nAUC=%.3f±%.3f", model_name, mean_auc, sd_auc) %>%
      fct_inorder()
  ) %>% 
  left_join(
    data_barplot %>% distinct(model, model_type, Cluster)
  )

performance_barplot <- auc_summary %>%
  ggplot() +
  geom_col(
    aes(x = model_name, y = mean_auc, fill = Cluster),
    # fill = "grey80",
    width = 0.6,
    alpha = 0.6,
    show.legend = F
  ) +
  geom_errorbar(
    data = auc_summary,
    aes(x = model_name, ymin = mean_auc - se_auc, ymax = mean_auc + se_auc),
    width = 0.2,
    color = "grey25",
    size = 0.8
  ) +
  geom_jitter(
    data = data_barplot %>%
      left_join(auc_summary),
    aes(x = model_name, y = auc, fill = batch),
    width = 0.1,
    shape = 21,
    size = 2,
    color = "grey50",
    show.legend = F
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  # ylim(0.4, 0.7) +
  coord_cartesian(ylim = c(0.45, 0.7)) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 14) +
  labs(y = "AUC Score", x = "", color = "Study") +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(
      # size = 13,
      face = "plain",
      color = "black"
    ),
    strip.placement = "outside"
  ) +
  scale_fill_manual(
    values = c(colors_clusters, map_colors(data_barplot$batch, "Set2")),
    na.value = "grey40"
  ) +
  facet_grid(~model_type, scales = "free", space = "free") +
  ggpubr::rotate_x_text(angle = 45)

performance_barplot


##### With ML models ------------------------------------------------------------

results_EN_pa_200f <- read_rds("output/data/metagenomics/microbeICB/ml/EN_pa_200f.rds")
results_rf_pa_200f <- read_rds("output/data/metagenomics/microbeICB/ml/rf_pa_200f.rds")

# Metacyc
results_EN_pa_200f_metacyc <- read_rds("output/data/metagenomics/microbeICB/ml/metacyc_EN_pa_200f.rds")
results_rf_pa_200f_metacyc <- read_rds("output/data/metagenomics/microbeICB/ml/metacyc_rf_pa_200f.rds")

get_table_results_ml <- function(
  results_ml, 
  model = model,
  features = 200,
  feature_type = feature_type,
  normalization = normalization,
  name = name
) {

  # After your loop completes, create a summary
  results_df <- do.call(
    rbind,
    lapply(results_ml, function(x) {
      data.frame(
        Study = x$fold,
        AUC = x$auc_test,
        AUC_train = x$auc_train,
        F1 = x$f1,
        Balanced_Acc = x$balanced_accuracy,
        Sensitivity = x$sensitivity,
        Specificity = x$specificity,
        N_test = x$n_test
      )
    })
  ) %>%
    mutate_if(is.numeric, ~ round(., 3))

  # Add mean and SD
  results_df <- rbind(
    results_df,
    data.frame(
      Study = "Mean ± SD",
      AUC = sprintf("%.3f ± %.3f", mean(results_df$AUC), sd(results_df$AUC)),
      AUC_train = sprintf("%.3f ± %.3f", mean(results_df$AUC_train), sd(results_df$AUC_train)),
      F1 = sprintf(
        "%.3f ± %.3f",
        mean(results_df$F1, na.rm = TRUE),
        sd(results_df$F1, na.rm = TRUE)
      ),
      Balanced_Acc = sprintf(
        "%.3f ± %.3f",
        mean(results_df$Balanced_Acc),
        sd(results_df$Balanced_Acc)
      ),
      Sensitivity = sprintf(
        "%.3f ± %.3f",
        mean(results_df$Sensitivity),
        sd(results_df$Sensitivity)
      ),
      Specificity = sprintf(
        "%.3f ± %.3f",
        mean(results_df$Specificity),
        sd(results_df$Specificity)
      ),
      N_test = sum(as.numeric(results_df$N_test[-nrow(results_df)]))
    )
  )

  results_df = results_df %>%
    mutate(
      model = NA,
      features = 200,
      feature_type = NA,
      normalization = NA,
      name = ifelse(is.na(name), deparse(substitute(results_ml)), name)
    )

  return(results_df)
}

data_barplot_ml = rbind(
  results_EN_pa_200f %>%
    get_table_results_ml(name = "Elastic Net SGBs"),
  results_rf_pa_200f %>%
    get_table_results_ml(name = "Random Forest SGBs"),
  #Metacyc
  results_EN_pa_200f_metacyc %>%
    get_table_results_ml(name = "Elastic Net Metacyc"), 
  results_rf_pa_200f_metacyc %>% 
    get_table_results_ml(name = "Random Forest Metacyc")
) %>%
  filter(Study != "Mean ± SD") %>%
  mutate(AUC = as.numeric(AUC), AUC_train = as.numeric(AUC_train)) %>% 
  mutate(
    model_type = "Machine Learning"
  )

auc_summary_ml = data_barplot_ml %>% 
  group_by(name) %>%
  summarise(
    mean_auc = mean(AUC),
    sd_auc = sd(AUC),
    n = n(),
    se_auc = sd_auc / sqrt(n)
  ) %>%
  ungroup() %>%
  arrange(name) %>%
  mutate(
    model_name = sprintf("%s\nAUC=%.3f±%.3f", name, mean_auc, sd_auc) %>%
      fct_inorder()
  ) %>% 
  mutate(
    model_type = "Machine Learning"
  )

data_barplot_all <- bind_rows(
  data_barplot_ml %>%
    mutate(auc = AUC, model = name) %>%
    select("Study", "auc", "model", "model_type"),
  data_barplot %>%
    rename(Study = batch) %>%
    select("Study", "auc", "model", "model_type")
) %>%
  mutate(
    model_type = fct_relevel(
      model_type,
      "Prevalence",
      "Abundance",
      "Community models",
      "Machine Learning"
    )
  )

auc_summary_all <- bind_rows(
  auc_summary_ml %>%
    mutate(model = name),
  auc_summary
) %>%
  mutate(
    model_type = fct_relevel(
      model_type,
      "Prevalence",
      "Abundance",
      "Community models",
      "Machine Learning"
    )
  ) %>% 
  mutate(Cluster = as.character(Cluster)) %>% 
  mutate(Cluster = ifelse(model_type == "Machine Learning", as.character(model_name), Cluster))

performance_barplot_all <- auc_summary_all %>% 
  ggplot() +
  geom_col(
    aes(x = model_name, y = mean_auc, fill = Cluster),
    # fill = "grey80",
    width = 0.6,
    alpha = 0.6,
    show.legend = F
  ) +
  geom_errorbar(
    data = auc_summary_all,
    aes(x = model_name, ymin = mean_auc - se_auc, ymax = mean_auc + se_auc),
    width = 0.2,
    color = "grey25",
    size = 0.8
  ) +
  geom_jitter(
    data = data_barplot_all %>%
      left_join(auc_summary_all),
    aes(x = model_name, y = auc, fill = Study),
    width = 0.1,
    shape = 21,
    size = 2,
    color = "grey50"
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  coord_cartesian(ylim = c(0.45, 0.7)) +
  theme_bw(base_size = 14) +
  labs(y = "AUC Score", x = "", color = "Study") +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(
      # size = 13,
      face = "plain",
      color = "black"
    ),
    strip.placement = "outside"
  ) +
  scale_fill_manual(
    values = c(
      colors_clusters,
      colors_studies,
      map_colors(auc_summary_ml$model_name, "Paired")
    ),
    na.value = "grey40"
  ) +
  facet_grid(~model_type, scales = "free", space = "free") +
  ggpubr::rotate_x_text(angle = 45)

performance_barplot_all

## Figure S4.A ## 
ggsave(
  filename = "output/figures/figures_manuscript/model_performance_all.pdf",
  performance_barplot_all,
  width = 12,
  height = 6
)




# Functional differences Community score ----------------------------------------



## Alpha diversity differences ------------------------------------------------------------

tse_int_metacyc <- read_rds("data/metagenomics/int.metacyc.rda") %>% 
  convertFromPhyloseq()

metacyc.names = data.table::fread("data/resources/map_metacyc-pwy_name.txt.gz", header = F) %>% 
  rename(ID = V1, name = V2) %>% 
  as_tibble()

rowData (tse_int_metacyc) <- rowtibble(tse_int_metacyc) %>% 
  rename (ID = Species) %>% 
  left_join (metacyc.names) %>% 
  DataFrame()

# Filter to only samples that are also in Metaphlan pre-treatment 
tse_int_metacyc_fil <- tse_int_metacyc %>% 
  filter_se(Sample %in% colData(tse_int.sgb.fil.pretreatment)$Sample)

colData(tse_int_metacyc_fil) <- coldf(tse_int_metacyc_fil) %>%
  inner_join(
    coldf(tse_int.sgb.fil.pretreatment) %>% select(Sample, ratio_abun)
  ) %>%
  mutate(batch = paste(Study, cohort, sep = ".") %>% str_replace(".NA", "")) %>%
  mutate(Study = ifelse(str_detect(Study, "Kim"), "Kim", Study)) %>%
  mutate(
    dysbiosis_level = case_when(
      ratio_abun < quantile(ratio_abun, 0.4, na.rm = T) ~ "Eubiotic",
      ratio_abun > quantile(ratio_abun, 0.6, na.rm = T) ~ "Dysbiotic",
      .default = "other"
    ),
    by = Study
  ) %>%
  mutate(Sample2 = Sample) %>%
  column_to_rownames(var = "Sample2") %>%
  DataFrame()


# Alpha diversity metrics. 
index_names = c("inverse_simpson_diversity", "shannon_diversity", "absolute_dominance", "simpson_evenness", "observed_richness", "gini_dominance", "pielou")
index = c("inverse_simpson", "shannon", "absolute", "simpson", "observed", "gini_dominance", "pielou_evenness")
tse_int_metacyc_fil <- addAlpha(
  tse_int_metacyc_fil, assay.type = "counts", index = index, name = index_names)


# Richness:
richness_metacyc <- coltibble(tse_int_metacyc_fil) %>%
  filter(dysbiosis_level != "other") %>%
  mutate(
    dysbiosis_type = ifelse(
      dysbiosis_level == "Dysbiotic",
      "High",
      "Low"
    )
  ) %>%
  ggplot(aes(y = observed_richness, x = dysbiosis_type)) +
  geom_quasirandom(size = 0.2) +
  geom_boxplot(outliers = F, aes(fill = Study), alpha = 0.3) +
  geom_signif(
    textsize = 3.2,
    comparisons = list(
      c("High", "Low")
    ),
    tip_length = 0.03,
    map_signif_level = function(p) {
      ifelse(
        p < 0.001,
        "p < 0.001",
        sprintf("p = %.3f", p)
      )
    }
  ) +
  scale_y_continuous(expand = expansion(c(0.05, 0.15))) +
  coord_cartesian(ylim = c(150, NA)) +
  facet_wrap(~Study, nrow = 1) +
  theme_few(base_size = 14) +
  labs(
    y = "Observed Richness\n(Metacyc pathways)",
    x = "Community score"
  ) +
  scale_fill_manual(
    values = colors_studies
  ) +
  ggpubr::rotate_x_text(angle = 45) 


# Diversity:
diversity_metacyc <- coltibble(tse_int_metacyc_fil) %>%
  filter(dysbiosis_level != "other") %>%
  mutate(
    dysbiosis_type = ifelse(
      dysbiosis_level == "Dysbiotic",
      "High",
      "Low"
    )
  ) %>%
  ggplot(aes(y = shannon_diversity, x = dysbiosis_type)) +
  geom_quasirandom(size = 0.2) +
  geom_boxplot(outliers = F, aes(fill = Study), alpha = 0.3) +
  geom_signif(
    textsize = 3.2,
    comparisons = list(
      c("High", "Low")
    ),
    tip_length = 0.03,
    map_signif_level = function(p) {
      ifelse(
        p < 0.001,
        "p < 0.001",
        sprintf("p = %.3f", p)
      )
    }
  ) +
  scale_y_continuous(expand = expansion(c(0.05, 0.15))) +
  coord_cartesian(ylim = c(4.7, NA)) +
  facet_wrap(~Study, nrow = 1) +
  theme_few(base_size = 14) +
  labs(
    y = "Shannon diversity\n(Metacyc pathways)",
    x = "Community score"
  ) +
  scale_fill_manual(
    values = colors_studies
  ) +
  ggpubr::rotate_x_text(angle = 45) 

richness_metacyc
diversity_metacyc

patch_metacyc_community <- diversity_metacyc / richness_metacyc +
  plot_layout(guides = "collect")

patch_metacyc_community



## Dysbiosis and all cohorts. ------------------------------------------------------------


# Load dysbiosis scores for all cohorts.

dysbiosis_eggnog <- read_csv("output/data/metagenomics/microbeICB/icb-eggnog-dysbiosis.csv") %>% 
  select(Sample, CentroidDist_score)

dysbiosis_all_data <- tse_int.sgb.fil.pretreatment %>% 
  coltibble() %>% 
    mutate(
    community_score_level = case_when(
      ratio_abun < quantile(ratio_abun, 0.4, na.rm = T) ~ "Low",
      ratio_abun > quantile(ratio_abun, 0.6, na.rm = T) ~ "High",
      .default = "other"
    ),
    by = Study
  )  %>% 
  left_join(dysbiosis_eggnog) %>% 
  mutate(
    dysbiosis_eggnog = scale(CentroidDist_score)[,1],
    .by = "Study"
  )

dysbiosis_all_boxplot <- dysbiosis_all_data %>% 
  filter(community_score_level != "other") %>%
  ggplot(aes(community_score_level, dysbiosis_eggnog)) +
  geom_quasirandom(size = 0.2) +
  geom_boxplot(outliers = F, aes(fill = Study), alpha = 0.3) +
  geom_signif(
    comparisons = list(
      c("High", "Low")
    ),
    tip_length = 0.03,
    map_signif_level = function(p) {
      ifelse(
        p < 0.001,
        "p < 0.001",
        sprintf("p = %.3f", p)
      )
    } 
  ) +
  scale_y_continuous(expand = expansion(c(0.05, 0.15))) +
  # coord_cartesian(ylim = c(150, NA)) +
  facet_wrap(~Study, nrow = 1) +
  theme_few(base_size = 14) +
  labs(
    y = "Functional dysbiosis",
    x = "Community score"
  ) +
  scale_fill_manual(
    values = colors_studies
  ) +
  ggpubr::rotate_x_text(angle = 45)

dysbiosis_all_boxplot

# ggsave(
#   "output/figures/figures_manuscript/functional_dysbiosis_all.pdf",
#   dysbiosis_all_boxplot,
#   height = 4,
#   width = 8.5
# )

# Plot together with richness and diversity of metacyc pathways.

patch_functional <-
  (((diversity_metacyc +
    labs(x = NULL) +
    theme(
      axis.text.x = element_blank(),
      axis.title.y = element_text(size = 12)
    )) /
    (richness_metacyc +
      theme(
        axis.title.y = element_text(size = 12)
      ))) |
    dysbiosis_all_boxplot) +
  plot_layout(guides = "collect")


## Figure 3.A & Figure 3.B ##
ggsave(
  "output/figures/figures_manuscript/functional_differences_all.pdf",
  patch_functional,
  height = 4.5,
  width = 15
)


# Correlation functional dysbiosis & s-score ------------------------------------- 

# Correlation between dysbiosis and community score in all cohorts.

dysbiosis_cor_all <- dysbiosis_all_data %>%
  select(dysbiosis_eggnog, ratio_abun, Study) %>%
  ggplot(aes(ratio_abun, dysbiosis_eggnog)) +
  geom_point(size = 3.5, fill = "#A7CEE2", shape = 21) +
  stat_cor(
    label.y = max(dysbiosis_all_data$dysbiosis_eggnog, na.rm = TRUE) * 1.1,
    label.x = -3
  ) +
  labs(
    x = "Community score",
    y = "Functional dysbiosis"
  ) +
  theme_bw(base_size = 18)

# Correlation between dysbiosis and community score in all cohorts.

s_score_cor_all <- dysbiosis_all_data %>% 
  select(ratio_abun, s_score, Study) %>% 
  ggplot(aes(ratio_abun, s_score)) +
  geom_point(size = 3.5, fill = "#A7CEE2", shape = 21) +
  stat_cor(
    label.y = max(dysbiosis_all_data$s_score) * 1.1,
    label.x = -3
  ) +
  labs(
    x = "Community score",
    y = "S score"
  ) +
  theme_bw(base_size = 18)

patch_cor_all <-  dysbiosis_cor_all + s_score_cor_all

## Figure S5.A & Figure S5.B ##
ggsave(
  filename = "output/figures/figures_manuscript/community_score_cor_dysbiosis_s_score_all.pdf",
  patch_cor_all,
  width = 12,
  height = 6
)


# Overlap with S-score ------------------------------------------------------------

# Overlap in terms of SGBs between toposcore and clusters. 

library(ggvenn)

# Overlap in terms of Species between toposcore and clusters. 
toposcore_species_sgbs <- read_xlsx(
  "data/resources/toposcore_species.xlsx",
  skip = 1,
  sheet = 2
) %>%
  mutate(species = MGS %>% str_replace_all(" ", "_")) %>%
  left_join(tax_info_Jan21 %>% distinct(species, .keep_all = T)) %>%
  mutate(
    Species = case_when(
      MGS ==
        "Candidatus Cibiobacter qucibialis" ~ "Candidatus Cibionibacter quicibialis",
      MGS ==
        "Clostridiacea eunclassified SGB4769" ~ "Clostridiaceae bacterium Marseille Q4145",
      MGS ==
        "Clostridiales unclassified SGB15145" ~ "Eubacterialesunc lassified SGB15145",
      MGS == "Enterocloster aldensis" ~ "Enterocloster aldenensis",
      MGS ==
        "Clostridiaceae bacterium OM08 6BH" ~ "Clostridiaceae bacterium Marseille Q4143",
      MGS == "Lachnospiraceae bacterium WCA3 601 WT 6H" ~ "Waltera intestinalis",
      MGS == "Lactobacillus vaginalis" ~ "Limosilactobacillus vaginalis",
      .default = MGS
    ) %>%
      str_replace_all(" ", "_")
  )

tax_info_Oct22

cluster_1_species <- tax_info_Oct22 %>%
  filter(SGB %in% sgb_cluster_1) %>%
  pull(Species) %>%
  unique()

cluster_2_species <- tax_info_Oct22 %>%
  filter(SGB %in% sgb_cluster_2) %>%
  pull(Species) %>%
  unique()

sig1_species <- toposcore_species_sgbs %>%
  filter(SIG == "SIG1") %>%
  pull(Species)

sig2_species <- toposcore_species_sgbs %>%
  filter(SIG == "SIG2") %>%
  pull(Species)


ven_good_bacteria <- ggvenn(
  list(
    SIG2 = sig2_species,
    Cluster1 = cluster_1_species
  ),
  fill_color = c("#7ea7b1ff", "#4DAF4A"),
  stroke_size = 0.5,
  set_name_size = 4,
  # text_size = 6,
  auto_scale = F
)

ven_bad_bacteria <- ggvenn(
  list(
    SIG1 = sig1_species,
    Cluster2 = cluster_2_species
  ),
  fill_color = c("#9F9596", "#E41A1C"),
  stroke_size = 0.5,
  set_name_size = 4,
  # text_size = 6,
  auto_scale = F
)

patch_venn <- ven_good_bacteria / ven_bad_bacteria

patch_venn

## Figure S5.C ##
ggsave(
  filename = "output/figures/figures_manuscript/venn_toposcore.pdf",
  patch_venn,
  width = 12,
  height = 6
)

#  Meta-analysis (Dysbiosis vs Eubiosis): ---------------------------------------

tse_int_metacyc_fil_2 <- tse_int_metacyc_fil %>% 
  filter_se(dysbiosis_level != "other") %>% 
  transformAssay(method = "relabundance")

## Adjust batch effects 

fit_adjust_batch_metacyc <- MMUPHin::adjust_batch(
  feature_abd = assay(tse_int_metacyc_fil_2, "relabundance"),
  batch = "Study",
  covariates = c("dysbiosis_level"),
  data = coldf(tse_int_metacyc_fil_2),
  control = list(verbose = T, diagnostic_plot = NULL)
)
assay(tse_int_metacyc_fil_2, "relabundance_adj") = fit_adjust_batch_metacyc$feature_abd_adj


## Meta-analysis 
metadata_metacyc = coldf(tse_int_metacyc_fil_2) %>%
  mutate(host_disease_status = as.factor(host_disease_status))
rownames(metadata_metacyc) = metadata_metacyc$Sample

fit_meta_metacyc_dysbiosis <- MMUPHin::lm_meta(
  feature_abd = assay(tse_int_metacyc_fil_2, "relabundance_adj"),
  batch = "Study",
  exposure = "dysbiosis_level",
  covariates = c("host_disease_status", "atb", "batch", "gender", "age"),
  control = list(
    # normalization = "CLR",
    # transform = "NONE",
    output = "output/maaslin2/MMUPHin_lm_meta_KO"
  ),
  data = metadata_metacyc
)

# Save results:
fit_meta_metacyc_dysbiosis$meta_fits %>%
  arrange(pval) %>%
  as_tibble() %>%
  mutate(fdr = p.adjust(pval, "fdr")) %>%
  mutate(ID = feature) %>%
  left_join(rowtibble(tse_int_metacyc_fil_2)) %>%
  # filter(fdr < 0.05) %>% 
  write_rds("output/data/metagenomics/microbeICB/meta_analysis_metacyc_high_low_community_score.csv")


fit_lm_community_score <- read_rds(
  "output/data/metagenomics/microbeICB/meta_analysis_metacyc_high_low_community_score.csv"
)

fit_lm_community_score %>% 
  filter(fdr < 0.05)
# 298/592 # 0.5033784 % of pathways with fdr < 0.05

# Load meta-analys response. 
fit_lm_metacyc <- read_rds("output/data/metagenomics/microbeICB/meta_analysis_metacyc.rds")

barplot_metacyc_microbiomes <- bind_rows(
    fit_lm_community_score %>%
      arrange(pval) %>%
      filter(coef > 0) %>%
      slice_head(n = 20),
    fit_lm_community_score %>%
      arrange(pval) %>%
      filter(coef < 0) %>%
      slice_head(n = 20)
  ) %>% 
  left_join(
    fit_lm_metacyc$meta_fits %>%
      mutate(
        pval_response = pval,
        sign_pval_response = -log10(pval) * sign(coef)
      ) %>%
      select(feature, pval_response, sign_pval_response)
  ) %>%
  mutate(
    `pval<0.05` = case_when(
      pval_response < 0.05 & sign(sign_pval_response) < 0 ~ "NR_associated",
      pval_response < 0.05 & sign(sign_pval_response) > 0 ~ "R_associated",
      .default = "ns"
    )
  ) %>%
  mutate(ID = feature) %>%
  left_join(rowtibble(tse_int_metacyc_fil_2)) %>%
  mutate(
    # fdr_significant = ifelse(fdr < 0.05, "fdr<0.05", "fdr>0.05"),
    signed.pvalue = -log10(qval.fdr) * sign(coef),
    name = make.unique(name),
    # short_name = str_trunc(name, 30, side = "left"),
    # short_name = fct_reorder(short_name, -signed.pvalue),
    name = name %>% str_to_sentence() %>% fct_reorder(signed.pvalue)
  ) %>%
  ggplot(aes(y = name, x = signed.pvalue)) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = 0,
    ymin = -Inf,
    ymax = Inf,
    fill = "#C22A44",
    alpha = 0.05
  ) +
  annotate(
    "rect",
    xmin = 0,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    fill = "#166B2A",
    alpha = 0.05
  ) +
  geom_segment(
    aes(
      x = signed.pvalue,
      xend = 0,
      y = name,
      yend = name
    ),
    size = 1.5
  ) +
  geom_point(
    aes(fill = sign_pval_response),
    size = 5,
    shape = 21
  ) +
  geom_point(
    data = . %>% filter(`pval<0.05` != "ns") %>% 
      mutate(signed.pvalue = ifelse(sign(sign_pval_response) > 0, signed.pvalue + 5.5, signed.pvalue - 5.5)),
    aes(color = `pval<0.05`),
    size = 2,
    shape = 8
  ) +
  geom_vline(
    xintercept = c(-log10(0.05), log10(0.05)),
    linetype = "dashed",
    alpha = 0.5
  ) +
  labs(
    y = NULL,
    x = "Community score association\n(signed log10 p-value)",
    fill = "Significance",
    subtitle = "Metacyc pathways"
  ) +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish,
    name = "Response association\n(signed log10 p-value)"
  ) +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(
    values = c("NR_associated" = "red", "R_associated" = "blue")
  ) +
  theme_few(base_size = 14) +
  theme(legend.position = "bottom", legend.box = "vertical")

barplot_metacyc_microbiomes +
  theme_few(base_size = 18) +
  theme(legend.position = "bottom", legend.box = "vertical")


# Pathways per bacteria clusters ------------------------------------------------

source("scripts/metagenomics/05_functional_per_bacteria.R")

# patch = barplot_metacyc_microbiomes +
#   (pathways_bacteria_contribution_barplot +
#     theme(legend.position = "bottom", legend.box = "vertical"))

# # Figure 3.C & Figure 3.D ## 
# ggsave(
#   "output/figures/figures_manuscript/metacyc_barplot_clusters.pdf",
#   patch,
#   width = 18,
#   height = 10
# )






# ------------------------------------------------------------------------------
# Metabolomics -----------------------------------------------------------------
# ------------------------------------------------------------------------------
# Metabolomics of the DRUP cohort 

se_metabolites_fil <- read_rds ("data/metabolomics/metabolomics.se.fil.rds")

# #Filter for compoud class:
se_metabolites_fil <- se_metabolites_fil[rowData(se_metabolites_fil)$annotation_level %in% c("1", "2a"),]

# Compute abundance log ratio for all the DRUP samples: 

tse_reads_species <- read_rds("output/data/metagenomics/drup/tse_reads_metaphlan.rds")

# Load score. 
colData(tse_reads_species) <- tse_reads_species %>%
  coldf() %>%
  mutate(
    abun_cluster_2 = get_abundance_set(
      tse_reads_species,
      sgb_cluster_2
    )$total_abun,
    preval_cluster_2 = get_abundance_set(
      tse_reads_species,
      sgb_cluster_2
    )$prevalence,
    abun_cluster_1 = get_abundance_set(
      tse_reads_species,
      sgb_cluster_1
    )$total_abun,
    preval_cluster_1 = get_abundance_set(
      tse_reads_species,
      sgb_cluster_1
    )$prevalence
  ) %>%
  mutate(
    # ratio_abun = abun_cluster_1 - abun_cluster_2,
    ratio_abun = log10(abun_cluster_2 + 1e-6) - log10(abun_cluster_1 + 1e-6),
    # ratio_abun = preval_cluster_2 / (preval_cluster_1 + 0.0001),
    ratio_prevalence = preval_cluster_2 - preval_cluster_1
  ) %>%
  DataFrame()

median_ratio_abun <- colData(tse_reads_species)$ratio_abun %>% median()


# Join the datasets
colData(se_metabolites_fil) <- coldf(se_metabolites_fil) %>%
  left_join(coldf(tse_reads_species) %>% select(Sample, ratio_abun)) %>%
  mutate(
    dysbiosis_level = case_when(
      ratio_abun < quantile(ratio_abun, 0.4, na.rm = T) ~ "Eubiotic",
      ratio_abun > quantile(ratio_abun, 0.6, na.rm = T) ~ "Dysbiotic",
      .default = "other"
    )
  ) %>% 
  mutate(Sample2 = Sample) %>%
  column_to_rownames(var = "Sample2") %>%
  DataFrame()

se_metabolites_fil %>% 
  coltibble() %>% 
  tabyl(tumor_class, dysbiosis_level)

# LME to check association between high dysbiotic and low dysbiotic compounds.
core.lme = se_metabolites_fil %>%
  meltSE(assay.type = "log10", add.row = T, add.col = T) %>%
  filter(annotation_level %in% c("1", "2a", "2b")) %>%
  filter(prevalence > 0.5) %>% 
  filter(!is.na(dysbiosis_level)) %>%
  filter(name != "Acadesine(AICAR)") %>%
  nest_by(compound_id) %>%
  mutate(
    model = list(broom.mixed::tidy(lmerTest::lmer(
      log10 ~ dysbiosis_level + ATB_use + tumor_class  + (1 | patientId),
      data = data
    )))
  ) %>%
  select(-data) %>%
  unnest() %>%
  ungroup() %>%
  filter(term == "dysbiosis_levelEubiotic") %>%
  mutate(fdr = p.adjust(p.value), p = sprintf("p = %.2g", p.value)) %>%
  left_join(rowtibble(se_metabolites_fil))

dysbiosis_vp <- core.lme %>%
  arrange(p.value) %>%
  mutate(
    short_name = str_trunc(name, 30, side = "left"),
  ) %>%
  ggplot(aes(estimate, -log10(p.value))) +
  # annotate ("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#166B2A", alpha = 0.6) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = 0,
    ymin = -Inf,
    ymax = Inf,
    fill =  "#E41A1C",
    alpha = 0.03
  ) +
  annotate(
    "rect",
    xmin = 0,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    fill = "#4DAF4A",
    alpha = 0.03
  ) +
  geom_hline(
    yintercept = c(-log10(0.05), -log10(0.01)),
    alpha = 0.7,
    linetype = "dashed"
  ) +
  geom_point(
    data = . %>% filter(fdr > 0.1)
  ) +
  geom_point(
    data = . %>% filter(fdr < 0.1) %>% mutate(fdr_level = "FDR<0.1"),
    aes(colour = fdr_level)
  ) +
  scale_color_manual(values = c("FDR<0.1" = "red")) +
  geom_text_repel(
    data = . %>%
      # filter(fdr < 0.25),
      filter(p.value < 0.01),
    aes(label = short_name),
    max.overlaps = Inf,
    size = 3
  ) +
  # xlim(-3,NA) +
  labs(x = "Coefficient", subtitle = "Stool metabolomics DRUP cohort", color = NULL) +
  theme_few(base_size = 14) +
  theme(
    # legend.position = "bottom",
    legend.position = c(0.98, 0.02),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )

dysbiosis_vp

ggsave(
  "output/figures/figures_manuscript/community_score_metabolomics_vp.pdf",
  dysbiosis_vp,
  width = 5,
  height = 8
)

dysbiosis_boxplot <- se_metabolites_fil %>%
  meltSE(assay.type = "log10", add.row = T, add.col = T) %>%
  filter(annotation_level %in% c("1", "2a", "2b")) %>%
  filter(!is.na(dysbiosis_level)) %>%
  filter(dysbiosis_level != "other") %>%
  # tabyl(class)
  left_join(core.lme) %>%
  filter(fdr < 0.1) %>%
  mutate(
    class = ifelse(
      dysbiosis_level == "Eubiotic",
      "Low\n(bottom 40%)",
      "High\n(top 40%)"
    )
  ) %>%
  mutate(name = fct_reorder(name, -log10(p.value) * sign(statistic))) %>%
  ggplot(aes(x = class, y = log10)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = log10(lod)), linetype = "dashed") +
  geom_quasirandom(size = 1) +
  facet_wrap(~name, scales = 'free', ncol = 3) +
  theme_few() +
  geom_signif(
    data = . %>% distinct(name, .keep_all = T),
    comparisons = list(c("Low\n(bottom 40%)", "High\n(top 40%)")),
    aes(annotations = p),
    manual = T,
    # map_signif_level = function(p) sprintf("p = %.2g", p)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = "Community score", y = "Abundance (log10)")



# Boxplot response significant metabolites

list_plots <- map(
  core.lme %>%
    # filter(p.value < 0.01) %>%
    pull(name) %>%
    set_names(),
  function(metabolite_name){
    data_model = se_metabolites_fil %>%
      meltSE(assay.type = "log10", add.row = T, add.col = T) %>%
      filter(annotation_level %in% c("1", "2a", "2b")) %>%
      filter(!is.na(dysbiosis_level)) %>%
      filter(dysbiosis_level != "other") %>%
      # tabyl(class)
      left_join(core.lme) %>%
      filter(name == metabolite_name) %>%
      filter(time != "T3") %>%
      mutate(
        class = ifelse(
          dysbiosis_level == "Eubiotic",
          "Low\n(bottom 40%)",
          "High\n(top 40%)"
        )
      ) %>%
      mutate(
        name = fct_reorder(name, -log10(p.value) * sign(statistic)),
        time = ifelse(time == "T1", "Pre-treatment", "On-treatment") %>%
          fct_relevel("Pre-treatment")
      )
    
    wilcox_data <- data_model %>%
      group_by(time) %>%
      summarise(p_value = wilcox.test(log10 ~ Response)$p.value) %>%
      ungroup() %>%
      mutate(
        signif = ifelse(p_value < 0.05, "p<0.05", "ns"),
        color = ifelse(p_value < 0.05, "red", "black"),
        p_label = sprintf("p = %.2g", p_value)
      )

    p = data_model %>%
      ggplot(aes(x = Response, y = log10)) +
      geom_boxplot(aes(fill = Response), outliers = F) +
      geom_hline(aes(yintercept = log10(lod)), linetype = "dashed") +
      # geom_quasirandom(aes(group = Response), dodge.width = 0.75, size = 1) +
      geom_quasirandom() +
      facet_wrap(~time, ncol = 2) +
      theme_few() +
      geom_signif(
        data = data_model %>% left_join(wilcox_data),
        comparisons = list(c("R", "NR")),
        aes(annotations = p_label, color = signif),
        manual = T
        # textsize = 3,
        # tip_length = 0.01
      ) +
      scale_color_manual(values = c("p<0.05" = "red", "ns" = "black")) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
      labs(x = "Response", y = "Abundance (log10)", subtitle = metabolite_name)

      return(p)
  }
)

names(list_plots)

plot_wrap <-
  patchwork::wrap_plots(
    c(
      list_plots$`Urocanic acid`,
      list_plots$`Adenine`,
      list_plots$Tryptophan,
      list_plots$`Indole-3-propionic acid`,
      list_plots$`Indole-3-methyl acetate`,
      list_plots$`5-Methoxytryptamine`,
      list_plots$`3-Methyladipic acid`,
      list_plots$`N-Acetylglutamic acid`,
      list_plots$Glycerol
    )
  ) +
  plot_layout(guides = "collect")

plot_wrap

# Boxplot of top metabolites High vs Low community score. 
boxplot_metabolites <- se_metabolites_fil %>%
  meltSE(assay.type = "log10", add.row = T, add.col = T) %>%
  filter(annotation_level %in% c("1", "2a", "2b")) %>%
  filter(!is.na(dysbiosis_level)) %>%
  filter(dysbiosis_level != "other") %>%
  # tabyl(class)
  left_join(core.lme) %>%
  filter(fdr < 0.1) %>%
  mutate(
    class = ifelse(
      dysbiosis_level == "Eubiotic",
      "Low",
      "High"
    )
  ) %>%
  mutate(name = fct_reorder(name, -log10(p.value) * sign(statistic))) %>%
  ggplot(aes(x = class, y = log10)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = log10(lod)), linetype = "dashed") +
  geom_quasirandom(size = 1) +
  facet_wrap(~name, scales = 'free', nrow = 2) +
  theme_few() +
  geom_signif(
    data = . %>% distinct(name, .keep_all = T),
    comparisons = list(c("Low", "High")),
    aes(annotations = p),
    manual = T,
    # map_signif_level = function(p) sprintf("p = %.2g", p)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = "Community score", y = "Abundance (log10)")

patch <- boxplot_metabolites /
  plot_wrap +
  plot_layout(heights = c(1, 2))

patch

ggsave(
  "output/figures/figures_manuscript/metabolites_intesterest_response_boxplot.pdf",
  patch,
  height = 14,
  width = 12
)

# Perform wixcoxon test for responders vs non-responders at the two timepoints. 
core.response.wilcox = se_metabolites_fil %>%
  meltSE(assay.type = "log10", add.row = T, add.col = T) %>%
  filter(annotation_level %in% c("1", "2a", "2b")) %>%
  filter(!is.na(dysbiosis_level)) %>%
  filter(dysbiosis_level != "other") %>%
  filter(time != "T3") %>%
  mutate(
    class = ifelse(
      dysbiosis_level == "Eubiotic",
      "Low\n(bottom 40%)",
      "High\n(top 40%)"
    )
  ) %>%
  mutate(
    time = ifelse(time == "T1", "Pre-treatment", "On-treatment") %>%
      fct_relevel("Pre-treatment")
  ) %>%
  nest_by(name, time) %>%
  summarise(res = wilcox.test(log10 ~ Response, data = data) %>% tidy()) %>%
  unnest() %>% 
  ungroup()

# -----------------------------------------------------------------------------#
### SCFA ... -------------------------------------------------------------------
# -----------------------------------------------------------------------------#

se_scfa <- read_rds("data/metabolomics/scfa.se.fil.rds")

# Compute abundance log ratio for all the DRUP samples:

tse_reads_species <- read_rds("output/data/metagenomics/drup/tse_reads_metaphlan.rds")

# Load score. 
colData(tse_reads_species) <- tse_reads_species %>%
  coldf() %>%
  mutate(
    abun_cluster_2 = get_abundance_set(
      tse_reads_species,
      sgb_cluster_2
    )$total_abun,
    preval_cluster_2 = get_abundance_set(
      tse_reads_species,
      sgb_cluster_2
    )$prevalence,
    abun_cluster_1 = get_abundance_set(
      tse_reads_species,
      sgb_cluster_1
    )$total_abun,
    preval_cluster_1 = get_abundance_set(
      tse_reads_species,
      sgb_cluster_1
    )$prevalence
  ) %>%
  mutate(
    ratio_abun = log10(abun_cluster_2 + 1e-6) - log10(abun_cluster_1 + 1e-6),
    ratio_prevalence = preval_cluster_2 - preval_cluster_1
  ) %>%
  DataFrame()


# Join the datasets
colData(se_scfa) <- coldf(se_scfa) %>%
  left_join(coldf(tse_reads_species) %>% select(Sample, ratio_abun)) %>%
  mutate(
    dysbiosis_level = case_when(
      ratio_abun < quantile(ratio_abun, 0.4, na.rm = T) ~ "Eubiotic",
      ratio_abun > quantile(ratio_abun, 0.4, na.rm = T) ~ "Dysbiotic",
      .default = "other"
    )
  ) %>% 
  mutate(Sample2 = Sample) %>%
  column_to_rownames(var = "Sample2") %>%
  DataFrame()

# LME to check association between high dysbiotic and low dysbiotic compounds.
scfa.lme = se_scfa %>%
  meltSE(assay.type = "log10", add.row = T, add.col = T) %>%
  filter(!is.na(dysbiosis_level)) %>%
  nest_by(compound_id) %>%
  mutate(
    model = list(broom.mixed::tidy(lmerTest::lmer(
      log10 ~ dysbiosis_level + ATB_use + tumor_class  + (1 | patientId),
      data = data
    )))
  ) %>%
  select(-data) %>%
  unnest() %>%
  ungroup() %>%
  filter(term == "dysbiosis_levelEubiotic") %>%
  mutate(fdr = p.adjust(p.value), p = sprintf("p = %.2g", p.value)) %>%
  left_join(rowtibble(se_scfa))

scfa.response.wilcox = se_scfa %>%
  meltSE(assay.type = "log10", add.row = T, add.col = T) %>%
  filter(!is.na(dysbiosis_level)) %>%
  filter(dysbiosis_level != "other") %>%
  filter(time != "T3") %>%
  mutate(
    class = ifelse(
      dysbiosis_level == "Eubiotic",
      "Low\n(bottom 40%)",
      "High\n(top 40%)"
    )
  ) %>%
  mutate(
    time = ifelse(time == "T1", "Pre-treatment", "On-treatment") %>%
      fct_relevel("Pre-treatment"), 
    name = FeatureID
  ) %>%
  nest_by(name, time) %>%
  summarise(res = wilcox.test(log10 ~ Response, data = data) %>% tidy()) %>%
  unnest() %>% 
  ungroup()

list_plots_scfa <- map(
  scfa.lme %>%
    pull(compound_id) %>%
    set_names(),
  function(metabolite_name) {
    p = se_scfa %>%
      meltSE(assay.type = "log10", add.row = T, add.col = T) %>%
      filter(dysbiosis_level != "other") %>%
      left_join(scfa.lme) %>%
      filter(compound_id == metabolite_name) %>%
      filter(time != "T3") %>%
      mutate(
        class = ifelse(
          dysbiosis_level == "Eubiotic",
          "Low\n(bottom 40%)",
          "High\n(top 40%)"
        )
      ) %>%
      mutate(
        time = as.factor(time)
      ) %>%
      ggplot(aes(x = Response, y = log10)) +
      geom_boxplot(aes(fill = Response), outliers = F) +
      geom_quasirandom() +
      facet_wrap(~time, ncol = 2) +
      theme_few() +
      geom_signif(
        comparisons = list(c("R", "NR")),
        map_signif_level = function(p) sprintf("Wilcoxon, p = %.2g", p)
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      labs(x = "Response", y = "Abundance (log10)", title = metabolite_name)

    return(p)
  }
)

names(list_plots_scfa)

wrap_plots(list_plots_scfa, guides = "collect")


#### Save table metabolomics ----------------------------------------------------

# Join scfa & untaregeted metabolites.

all_metabolites_lme <- bind_rows(
  core.lme %>%
    left_join(
      rowtibble(se_metabolites_fil) %>% select(compound_id, name)
    ) %>%
    mutate(compound_id = name) %>%
    select(compound_id, term, estimate, std.error, statistic, df, p.value, fdr),
  scfa.lme %>%
    select(compound_id, term, estimate, std.error, statistic, df, p.value, fdr)
) %>%
  mutate(term = "community_scoreLow")

all_metabolites_response.wilcox <- bind_rows(
  scfa.response.wilcox,
  core.response.wilcox
)


openxlsx::write.xlsx(
  list(
    "Community_score_association" = all_metabolites_lme,
    "Response_association" = all_metabolites_response.wilcox
  ),
  file = "output/data/tables_manuscript/pathways_and_metabolites_da.xlsx",
  asTable = TRUE
)

# Save metadata -----------------------------------------------------------------

tse_int.sgb.fil %>% 
  coltibble() %>% 
  rename(
    community_score_abundance = ratio_abun,
    community_score_prevalence = ratio_prevalence
  ) + 
  write_csv("output/data/tables_manuscript/metadata_cohort.csv")

