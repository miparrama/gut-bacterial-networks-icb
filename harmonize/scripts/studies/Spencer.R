# Process of the Spencer dataset
# @ Miguel Parra Martinez

here::here("scripts/studies/Spencer")
library(here)
source(here("scripts", "helpers.R"))

# QC --------------------------------------------------------------------

## Kneaddata ------------------------------------------------------------

kneaddata.data = list.files(
  here("data/studies/spencer_science_2021/results"),
  pattern = "*-kneaddata.txt",
  recursive = T,
  full.names = T
) %>%
  map_dfr(data.table::fread) %>%
  filter(!is.na(Sample)) %>%
  janitor::clean_names() %>%
  filter(sample != "") %>%
  mutate(sample = str_replace(sample, "-kd", "")) %>%
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
  mutate(Study = "Spencer") %>%
  mutate(sample = str_replace(sample, "-clean", "")) %>%
  mutate(sample = str_replace(sample, "-kd", "")) %>%
  mutate(clean.reads = final, total.reads = raw) 

# Plot distribution of reads
kd.qc.read_distribution <- kneaddata.data %>% 
  ggplot(aes(y = clean.reads, x = "Spencer")) +
  geom_boxplot() + 
  labs(y = "Clean reads") +
  expand_limits(y = 0) +
  theme_bw()

# Percentages in log format for plotting. 
kneaddata.report <- kneaddata.data %>%
  filter(!is.na(sample)) %>%
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


# Load data ------------------------------------------------------------

spencer.metadata = read_csv(here(
  "data/studies/spencer_science_2021/SraRunTable.csv"
))

spencer.metadata = spencer.metadata %>%
  janitor::clean_names() %>%
  mutate(sample_id = run) %>%
  column_to_rownames("sample_id") %>%
  mutate(Response = if_else(response == 1, "Responder", "Non-responder")) %>%
  mutate(Response = fct_relevel(Response, "Non-responder")) %>%
  mutate(Study = "Spencer") %>%
  mutate(host_disease_status = "MELA") %>%
  mutate(sex = if_else(sex == "female", "female", "male")) %>%
  rename(Sample = run, gender = sex, sampleId = sample_name) %>%
  mutate(age = as.numeric(age)) %>%
  mutate(age_bin = cut(age, breaks = seq(10, 100, 5))) %>%
  mutate(age_bin = str_remove(age_bin, "\\(")) %>%
  mutate(age_bin = str_remove(age_bin, "\\]")) %>%
  mutate(age_bin = str_replace(age_bin, ",", "-")) %>%
  mutate(atb = if_else(antibiotics == 1, "yes", "no")) %>%
  mutate(atb = replace_na(atb, "no")) %>%
  mutate(pfs_event = pfsevent, pfs_time = pfs_d) %>%
  select(
    Sample,
    sampleId,
    Study,
    host_disease_status,
    Response,
    gender,
    age,
    age_bin,
    atb,
    bmi,
    pfs_event, 
    pfs_time
  ) %>%
  mutate(Response = as.character(Response)) %>%
  mutate(timepoint = "PreTreatment") %>% 
  mutate(therapy = 'anti-PD1 / anti-PD1 + anti-CTL4') %>%
  left_join(
    kneaddata.data %>%
      select(sample, clean.reads, human.reads, total.reads) %>%
      rename(Sample = sample)
  )

rownames(spencer.metadata) <- spencer.metadata$Sample

# Make metaphlan4 object
spencer.files = list.files(
  here("data/studies/spencer_science_2021/results/metaphlan/read_stats/"),
  full.names = T
)
names(spencer.files) = list.files(here(
  "data/studies/spencer_science_2021/results/metaphlan/read_stats/"
)) %>%
  basename() %>%
  str_remove("-reads.txt")
spencer.pseq = make_phyloseq4(
  table_fn = spencer.files,
  sampledata = spencer.metadata
)
write_rds(spencer.pseq, here("output/data/taxaprofiles/spencer-metaphlan.rda"))

# Import Metacyc abundances
spencer.metacyc = import_metacyc(
  here(
    "data/studies/spencer_science_2021/results/humann/pathway_abundance.txt"
  ),
  spencer.metadata,
  taxonomy = F
)
write_rds(spencer.metacyc, here("output/data/metacyc/spencer-metacyc.rda"))

# Import Metacyc abundances with taxonomy 
spencer.metacyc.tax = import_metacyc(
  here("data/studies/spencer_science_2021/results/humann/pathway_abundance.txt"),
  spencer.metadata,
  taxonomy = T
)
write_rds(spencer.metacyc.tax, here("output/data/regroup/spencer-metacyc.tax.rda"))

# Import Metacyc abundances
spencer.metacyc.cpm = import_metacyc(
  here(
    "data/studies/spencer_science_2021/results/humann/pathway_abundance-cpm.txt"
  ),
  spencer.metadata,
  taxonomy = F
)
write_rds(spencer.metacyc.cpm, here("output/data/metacyc/spencer-metacyc-cpm.rda"))

# Import KEGG abundances
spencer.kegg = import_metacyc(
  here(
    "data/studies/spencer_science_2021/results/humann/regroup/kegg_families.txt"
  ),
  spencer.metadata,
  type = "kegg",
  taxonomy = F
)
write_rds(spencer.kegg, here("output/data/regroup/spencer-kegg.rda"))

# Import KEGG abundances with taxonomy 
spencer.kegg.tax = import_metacyc(
  here("data/studies/spencer_science_2021/results/humann/regroup/kegg_families.txt"),
  spencer.metadata,
  type = "kegg",
  taxonomy = T
)
write_rds(spencer.kegg.tax, here("output/data/regroup/spencer-kegg.tax.rda"))

# Import eggNOG abundances
spencer.eggnog = import_metacyc(
  here(
    "data/studies/spencer_science_2021/results/humann/regroup/eggnog_families.txt"
  ),
  spencer.metadata,
  type = "kegg",
  taxonomy = F
)
write_rds(spencer.eggnog, here("output/data/regroup/spencer-eggnog.rda"))

# Import KEGG-compound abundances
spencer.compound = import_metacyc(
  here(
    "data/studies/spencer_science_2021/results/humann/regroup/kegg_compounds.txt"
  ),
  spencer.metadata
)
write_rds(spencer.compound, here("output/data/regroup/spencer-compound.rda"))

# Import KEGG-module abundances
spencer.module = import_metacyc(
  here(
    "data/studies/spencer_science_2021/results/humann/regroup/kegg_modules.txt"
  ),
  spencer.metadata
)
write_rds(spencer.module, here("output/data/regroup/spencer-module.rda"))
