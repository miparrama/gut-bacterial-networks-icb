# Process of the Kim dataset
# @ Miguel Parra Martinez

here::here("scripts/studies/Kim_1st.R")
library(here)
source(here("scripts", "helpers.R"))

# QC --------------------------------------------------------------------

## Kneaddata ------------------------------------------------------------

kneaddata.data = list.files(
  here("data/studies/kim_natcom_2024/results-kim-1st"),
  pattern = "*-kneaddata.txt",
  recursive = T,
  full.names = T
) %>%
  map_dfr(data.table::fread) %>%
  filter(!is.na(Sample)) %>%
  janitor::clean_names() %>%
  filter(sample != "") %>%
  # as_tibble() %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0))) %>%
  mutate(
    raw = raw_pair1 + raw_pair2,
    trimmed = trimmed_pair1 +
      trimmed_pair2 +
      trimmed_orphan1 +
      trimmed_orphan2,
    decontaminated = decontaminated_hg37dec_v0_1_pair1 +
      decontaminated_hg37dec_v0_1_pair2 +
      decontaminated_hg37dec_v0_1_orphan1 +
      decontaminated_hg37dec_v0_1_orphan2,
    final = final_pair1 +
      final_pair2 +
      final_orphan1 +
      final_orphan2
  ) %>%
  mutate(
    trimmed.reads = raw - trimmed,
    human.reads = trimmed - decontaminated
  ) %>%
  mutate(percent_lowqual = trimmed.reads / raw) %>%
  mutate(percent_human = human.reads / raw) %>%
  mutate(recovery_rate = final / raw) %>%
  # select (raw, trimmed, decontaminated, final) %>%
  mutate(sample = str_replace(sample, "-clean", "")) %>%
  mutate(sample = str_replace(sample, "-kd", "")) %>%
  mutate(clean.reads = final, total.reads = raw) %>% 
  mutate(Study = "Kim") %>%
  filter(!is.na(sample))

# Plot distribution of reads
kd.qc.read_distribution <- kneaddata.data %>% 
  ggplot(aes(y = clean.reads, x = "Kim_1st")) +
  geom_boxplot() + 
  labs(y = "Clean reads") +
  expand_limits(y = 0) +
  theme_bw()

kneaddata.report = kneaddata.data  %>% 
  # gather(Metric, counts, -c(sample, Study)) %>%
  pivot_longer(
    cols = -c(sample, Study),
    names_to = "Metric",
    values_to = "counts"
  ) %>%
  filter(
    Metric %in% c("percent_lowqual", "percent_human", "recovery_rate", "Study")
  )

# What is the recovery rate (% retained for analysis)
kd.qc = kneaddata.report %>%
  filter(Metric %in% c("percent_lowqual", "percent_human", "recovery_rate")) %>%
  mutate(
    Metric = fct_recode(
      Metric,
      Human = "percent_human",
      `Low quality` = "percent_lowqual",
      `Remaining reads` = "recovery_rate"
    )
  ) %>%
  ggplot(aes(x = sample, y = counts, fill = Metric)) +
  geom_col() +
  geom_point(aes(x = sample, y = 0), shape = 15) +
  theme_few(base_size = 11) +
  theme(axis.text.x = element_blank(), legend.position = "bottom") +
  scale_fill_d3() +
  # facet_grid(~Study, scales = "free")+
  scale_y_continuous(labels = scales::percent) +
  xlab("Samples") +
  ylab("Percent of sequenced reads")

# What is the sequenced depth of coverage per sample
kd.qc.box = kneaddata.report %>%
  filter(Metric %in% c("percent_lowqual", "percent_human", "recovery_rate")) %>%
  mutate(
    Metric = fct_recode(
      Metric,
      Human = "percent_human",
      `Low quality` = "percent_lowqual",
      `Remaining \n reads` = "recovery_rate"
    )
  ) %>%
  ggplot(aes(x = Metric, y = counts, fill = Metric)) +
  stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.15) +
  geom_boxplot(width = 0.5) +
  theme_few(base_size = 11) +
  guides(fill = FALSE) +
  scale_fill_d3() +
  facet_grid(~Study, scales = "free") +
  ggpubr::rotate_x_text(90) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
  xlab("") +
  ylab("")

kd.qc | kd.qc.box | kd.qc.read_distribution


# Filter samples with low QC. 
kneaddata.data <- kneaddata.data %>% 
  filter(percent_lowqual < 0.75) %>% 
  filter(total.reads > 5e+06)

# Load data ------------------------------------------------------------

kim_1st.metadata = read_csv(
  here("data/studies/kim_natcom_2024/SraRunTable_kim1st.csv")
)

kim_1st.metadata = read_csv(
  here("data/studies/kim_natcom_2024/SraRunTable_kim1st.csv")
) %>%
  mutate(sample_id = Run, Sample = Run, sampleId = Run) %>%
  column_to_rownames("sample_id") %>%
  mutate(
    Response = if_else(
      response == "R",
      "Responder",
      "Non-responder"
    )
  ) %>%
  mutate(Response = fct_relevel(Response, "Non-responder")) %>%
  mutate(host_disease_status = Cancer) %>%
  mutate(Study = "Kim_1st") %>%
  mutate(gender = if_else(Sex == "F", "female", "male")) %>%
  mutate(age_bin = `Age_bin_x_to_x+5`, age = Age) %>%
  mutate(atb = NA) %>%
  select(
    Sample,
    sampleId,
    Study,
    host_disease_status,
    Response,
    gender,
    age,
    age_bin,
    atb
  ) %>%
  mutate(Response = as.character(Response)) %>%
  mutate(timepoint = "PreTreatment") %>% 
  mutate(therapy = 'anti-PD1') %>%
  left_join(
    kneaddata.data %>%
      select(sample, clean.reads, human.reads, total.reads) %>%
      rename(Sample = sample)
  )

rownames(kim_1st.metadata) <- kim_1st.metadata$Sample

# Make metaphlan4 object
kim_1st.files = list.files(
  here("data/studies/kim_natcom_2024/results-kim-1st/metaphlan/read_stats/"),
  full.names = T
)
names(kim_1st.files) = list.files(
  here("data/studies/kim_natcom_2024/results-kim-1st/metaphlan/read_stats/")
) %>%
  basename() %>%
  str_remove("-reads.txt")

kim_1st.pseq = make_phyloseq4(
  table_fn = kim_1st.files,
  sampledata = kim_1st.metadata
)
write_rds(kim_1st.pseq, here("output/data/taxaprofiles/kim-1st-metaphlan.rda"))

# Import Metacyc abundances
kim_1st.metacyc = import_metacyc(
  here("data/studies/kim_natcom_2024/results-kim-1st/humann/pathway_abundance.txt"),
  kim_1st.metadata,
  taxonomy = F
)
write_rds(kim_1st.metacyc, here("output/data/metacyc/kim-1st-metacyc.rda"))

# Import Metacyc abundances with taxonomy 
kim_1st.metacyc.tax = import_metacyc(
  here("data/studies/kim_natcom_2024/results-kim-1st/humann/pathway_abundance.txt"),
  kim_1st.metadata,
  taxonomy = T
)
write_rds(kim_1st.metacyc.tax, here("output/data/regroup/kim-1st-metacyc.tax.rda"))

# Import Metacyc abundances (cpm)
kim_1st.metacyc.cpm = import_metacyc(
  here("data/studies/kim_natcom_2024/results-kim-1st/humann/pathway_abundance-cpm.txt"),
  kim_1st.metadata,
  taxonomy = F
)
write_rds(kim_1st.metacyc.cpm, here("output/data/metacyc/kim-1st-metacyc-cpm.rda"))

# Import KEGG abundances
kim_1st.kegg = import_metacyc(
  here(
    "data/studies/kim_natcom_2024/results-kim-1st/humann/regroup/kegg_families.txt"
  ),
  kim_1st.metadata,
  type = "kegg",
  taxonomy = F
)
write_rds(kim_1st.kegg, here("output/data/regroup/kim-1st-kegg.rda"))

# Import KEGG abundances with taxonomy 
kim_1st.kegg.tax = import_metacyc(
  here("data/studies/kim_natcom_2024/results-kim-1st/humann/regroup/kegg_families.txt"),
  kim_1st.metadata,
  type = "kegg",
  taxonomy = T
)
write_rds(kim_1st.kegg.tax, here("output/data/regroup/kim-1st-kegg.tax.rda"))

# Import eggNOG abundances
kim_1st.eggnog = import_metacyc(
  here(
    "data/studies/kim_natcom_2024/results-kim-1st/humann/regroup/eggnog_families.txt"
  ),
  kim_1st.metadata,
  type = "kegg",
  taxonomy = F
)
write_rds(kim_1st.eggnog, here("output/data/regroup/kim-1st-eggnog.rda"))

# Import KEGG-compound abundances
kim_1st.compound = import_metacyc(
  here(
    "data/studies/kim_natcom_2024/results-kim-1st/humann/regroup/kegg_compounds.txt"
  ),
  kim_1st.metadata
)
write_rds(kim_1st.compound, here("output/data/regroup/kim-1st-compound.rda"))

# Import KEGG-module abundances
kim_1st.module = import_metacyc(
  here(
    "data/studies/kim_natcom_2024/results-kim-1st/humann/regroup/kegg_modules.txt"
  ),
  kim_1st.metadata
)
write_rds(kim_1st.module, here("output/data/regroup/kim-1st-module.rda"))