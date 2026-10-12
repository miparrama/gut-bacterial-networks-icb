# Meta-analysis of SGBs associated with disease. 
# @ Miguel Parra Martinez

source("scripts/helpers.R")
library(MMUPHin)

# Meta-analysis SGB ------------------------------------------------------------

# Load data 
tse_int.SGB.fil <- read_rds("output/data/metagenomics/microbeICB/tse_int.SGB.rds") %>% 
  subsetByPrevalent(prevalence = 0.05)

# Meta-analaysis
metadata_pretreatment = coldf(tse_int.SGB.fil) %>% 
  mutate(host_disease_status = as.factor(host_disease_status)) %>% 
  filter(timepoint == "PreTreatment") %>% 
  mutate(batch = paste(Study, cohort, sep = ".") %>% str_replace(".NA", "")) %>% 
  mutate(n = n(), .by = batch) %>% 
  filter(n > 15) %>%
  mutate(Study = ifelse(str_detect(Study, "Kim"), "Kim", Study))

metadata_pretreatment %>% 
  tabyl(batch)

metadata_pretreatment %>% 
  tabyl(Study)

# Filter and update metadata
tse_int.sgb.fil.pretreatment <- tse_int.SGB.fil[,rownames(metadata_pretreatment)]
colData(tse_int.sgb.fil.pretreatment) <- DataFrame(metadata_pretreatment)

fit_lm_sgb <- MMUPHin::lm_meta(
  feature_abd = assay(tse_int.sgb.fil.pretreatment, "relabundance_adj"),
  batch = "Study",
  exposure = "Response",
  covariates = c("host_disease_status", "atb", "batch", "gender", "age"),
  control = list(
    output = "output/maaslin2/MMUPHin_lm_meta_SGB"
  ),
  data = metadata_pretreatment
)

fit_lm_sgb <- write_rds("output/data/metagenomics/microbeICB/meta_analysis_SGB.rds")





# Meta-analysis LOO (SGB) ------------------------------------------------------------

# Run meta-analysis with LOO crossvaliation for all the cohorts. 

# Filter samples 
tse_int.sgb.fil.pretreatment <- tse_int.sgb.fil.pretreatment[,rownames(metadata_pretreatment)]

# Add metadata: 
colData(tse_int.sgb.fil.pretreatment) <- DataFrame(metadata_pretreatment)

# Loop through all the Studies.  to performe meta-analysis. 
Studies <- unique(metadata_pretreatment$Study)

meta_analysis_list <- map(Studies %>% setNames(Studies), function(batch_name) {

  set.seed(1234)
  
  tse_int.species.fil_x <- filter_se(tse_int.sgb.fil.pretreatment, batch != batch_name)
  tse_int.species.fil_y <- filter_se(tse_int.sgb.fil.pretreatment, batch == batch_name)
  
  metadata_x = metadata_pretreatment %>% 
    filter(batch != batch_name)
  
  fit_lm_species <- MMUPHin::lm_meta(
    feature_abd = assay(tse_int.species.fil_x, "relabundance_adj"),
    batch = "batch",
    exposure = "Response",
    covariates = c("host_disease_status", "atb"),
    control = list(
      output = paste0("output/maaslin2/",batch_name,"_as_test_set")
    ),
    data = metadata_x
  )
  
  # Meta-analysis. 
  fit_lm_species.res = fit_lm_species$meta_fits %>% 
    mutate(fdr = p.adjust(pval, "fdr")) %>% 
    mutate(val = -log10(pval) * sign(coef)) %>% 
    as_tibble() %>% 
    arrange (pval)
  
  # Get the auc of the Significant species. 
  fit_lm_species.res <- fit_lm_species.res %>%
    mutate(auc = pmap_dbl(list(feature, coef), function(feat, c) {
      dir <- if (c > 0) "<" else ">"
      roc_from_taxa(tse_int.species.fil_y, feat, direction_roc = dir)$auc_table$auc
    }))
  
  species_beneficial <- fit_lm_species.res %>% 
    filter(pval < 0.05) %>% 
    filter(coef > 0) %>% 
    pull(feature)
  
  species_detrimental <- fit_lm_species.res %>% 
    filter(pval < 0.05) %>% 
    filter(coef < 0) %>% 
    pull(feature)
  
  colData(tse_int.species.fil_y) <-  tse_int.species.fil_y %>% 
    coldf() %>% 
    mutate(
      median_beneficial = get_median_abundance_set(
        tse_int.species.fil_y,
        species_beneficial
      )$mean_abun,
      median_detrimental = get_median_abundance_set(
        tse_int.species.fil_y,
        species_detrimental
      )$mean_abun,
      preval_detrimental = get_median_abundance_set(
        tse_int.species.fil_y,
        species_detrimental
      )$prevalence,
      preval_beneficial = get_median_abundance_set(
        tse_int.species.fil_y,
        species_beneficial
      )$prevalence
    ) %>% 
    mutate(preval_score = ((preval_beneficial - preval_detrimental) + 1) / 2) %>% 
    mutate(ratio = median_beneficial/median_detrimental) %>%
    DataFrame()
  
  roc_preval_clusters <- roc_from_col(
    tse_int.species.fil_y,
    "preval_score",
    direction_roc = "<",
    batch_col = "batch",
    base_size = 14, 
    make_plot = TRUE
  )
  
  list_results <- list(
    "res_meta" = fit_lm_species.res,
    "roc" = roc_preval_clusters,
    "tse" = tse_int.species.fil_y
  )
  
  return(list_results)
  
},.progress = T)

# Save results
meta_analysis_list %>% 
  write_rds("output/data/metagenomics/microbeICB/meta_analysis_LOO_SGB.rds")



# Metacyc ------------------------------------------------------------

metacyc_names <- read_tsv("data/resources/map_metacyc-pwy_name.txt.gz", col_names = c("FeatureID", "Description"))

tse_metacyc <- read_rds(
  "data/metagenomics/int.metacyc.rda"
) %>%
  mia::convertFromPhyloseq() %>%
  filter_se(Sample %in% colData(tse_int.sgb.fil.pretreatment)$Sample) %>% # Filter for pre-treatment samples. 
  transformAssay(method = "relabundance") %>%
  transformAssay(
    method = "clr",
    assay.type = "relabundance",
    pseudocount = T
  ) %>% 
    subsetByPrevalent(prevalence = 0.05) # Filter by prevalence ....

colData(tse_metacyc) <- coldf(tse_metacyc) %>%
  mutate(batch = paste(Study, cohort, sep = ".") %>% str_replace(".NA", "")) %>%
  mutate(Study = ifelse(str_detect(Study, "Kim"), "Kim", Study)) %>% 
  DataFrame()

## Adjust batch effects ------------------------------------------------------

fit_adjust_batch_metacyc <- MMUPHin::adjust_batch(
  feature_abd = assay(tse_metacyc, "relabundance"),
  batch = "Study",
  covariates = c("Response"),
  data = coldf(tse_metacyc),
  control = list(verbose = T, diagnostic_plot = NULL)
)

assay(tse_metacyc, "relabundance_adj") = fit_adjust_batch_metacyc$feature_abd_adj

## Meta-analysis ------------------------------------------------------------

metadata_metacyc = coldf(tse_metacyc) %>%
  left_join(
    coldf(tse_int.sgb.fil.pretreatment) %>%
      select(Sample, batch)
  ) %>%
  mutate(host_disease_status = as.factor(host_disease_status))
rownames(metadata_metacyc) = metadata_metacyc$Sample

fit_lm_metacyc <- MMUPHin::lm_meta(
  feature_abd = assay(tse_metacyc, "relabundance_adj"),
  batch = "Study",
  exposure = "Response",
  covariates = c("host_disease_status", "atb", "batch", "gender", "age"),
  control = list(
    # normalization = "CLR",
    # transform = "NONE",
    output = "output/maaslin2/MMUPHin_lm_meta_KO"
  ),
  data = metadata_metacyc
)

# Save results: 
fit_lm_metacyc %>%
  write_rds("output/data/metagenomics/microbeICB/meta_analysis_metacyc.rds")

# Base Meta-analysis -----------------------------------------------------------

# Meta-analysis without adjusting for covariates for the different types of studies. 

## SGB -------------------------------------------------------------------------
tse_int.sgb.fil.pretreatment <- tse_int.sgb.fil.pretreatment %>% 
  transformAssay(method = "relabundance")

fit_lm_sgb_nocovariates <- MMUPHin::lm_meta(
  feature_abd = assay(tse_int.sgb.fil.pretreatment, "relabundance"),
  batch = "Study",
  exposure = "Response",
  control = list(
    output = "output/maaslin2/MMUPHin_lm_meta_SGB_nocovariates"
  ),
  data =  coldf(tse_int.sgb.fil.pretreatment)
)

# Save results,
write_rds(
  fit_lm_sgb_nocovariates,
  file = "output/data/metagenomics/microbeICB/meta-analysis-nocovariates/sgb.mmuphin.rds"
)

# LOSO meta-analysis. 

studies <- tse_int.sgb.fil.pretreatment$Study %>% unique() %>% set_names()

list_meta_analysis <- map(
  studies,
  function(study) {
    print(study)
    tse_no_study = tse_int.sgb.fil.pretreatment %>% filter_se(Study != study)

    fit_lm_metacyc_nocovariates <- MMUPHin::lm_meta(
      feature_abd = assay(tse_no_study, "relabundance"),
      batch = "Study",
      exposure = "Response",
      control = list(
        output = "output/maaslin2/MMUPHin_lm_meta_SGB_nocovariates"
      ),
      data = coldf(tse_no_study)
    )
    return(fit_lm_metacyc_nocovariates)
  }
)

write_rds(
  list_meta_analysis,
  file = "output/data/metagenomics/microbeICB/meta-analysis-nocovariates/sgb.mmuphin.loso.rds"
)


## Metacyc ---------------------------------------------------------------------

metacyc_names <- read_tsv("data/resources/map_metacyc-pwy_name.txt.gz", col_names = c("FeatureID", "Description"))

tse_metacyc <- read_rds(
  "data/metagenomics/int.metacyc.rda"
) %>%
  mia::convertFromPhyloseq() %>%
  filter_se(Sample %in% colData(tse_int.sgb.fil.pretreatment)$Sample) %>% # Filter for pre-treatment samples. 
  transformAssay(method = "relabundance") %>%
  transformAssay(
    method = "clr",
    assay.type = "relabundance",
    pseudocount = T
  ) %>% 
    subsetByPrevalent(prevalence = 0.05) # Filter by prevalence ....

colData(tse_metacyc) <- coldf(tse_metacyc) %>%
  mutate(batch = paste(Study, cohort, sep = ".") %>% str_replace(".NA", "")) %>%
  mutate(Study = ifelse(str_detect(Study, "Kim"), "Kim", Study)) %>% 
  DataFrame()

fit_lm_metacyc_nocovariates <- MMUPHin::lm_meta(
  feature_abd = assay(tse_metacyc, "relabundance"),
  batch = "Study",
  exposure = "Response",
  control = list(
    output = "output/maaslin2/MMUPHin_lm_meta_metacyc_nocovariates"
  ),
  data = coldf(tse_metacyc)
)

# Save results,
write_rds(
  fit_lm_metacyc_nocovariates,
  file = "output/data/metagenomics/microbeICB/meta-analysis-nocovariates/metacyc.mmuphin.rds"
)

# LOSO meta-analysis. 
studies <- tse_int.sgb.fil.pretreatment$Study %>% unique() %>% set_names()

list_meta_analysis_metacyc <- map(
  studies,
  function(study) {
    print(study)
    tse_no_study = tse_metacyc %>% filter_se(Study != study)

    fit_lm_sgb_nocovariates <- MMUPHin::lm_meta(
      feature_abd = assay(tse_no_study, "relabundance"),
      batch = "Study",
      exposure = "Response",
      control = list(
        output = "output/maaslin2/MMUPHin_lm_metacyc_nocovariates"
      ),
      data = coldf(tse_no_study)
    )
    return(fit_lm_sgb_nocovariates)
  }
)

write_rds(
  list_meta_analysis_metacyc,
  file = "output/data/metagenomics/microbeICB/meta-analysis-nocovariates/metacyc.mmuphin.loso.rds"
)