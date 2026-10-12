# Process of the Routy dataset
# @ Miguel Parra Martinez

here::here("scripts/studies/Routy.R")
library(here)
source(here("scripts", "helpers.R"))

# QC --------------------------------------------------------------------

## Kneaddata ------------------------------------------------------------

kneaddata.data = list.files(
  here("data/studies/routy_science_2018"),
  pattern = "*-kneaddata.txt",
  recursive = T,
  full.names = T
) %>%
  map_dfr(data.table::fread) %>%
  filter(!is.na(Sample)) %>%
  janitor::clean_names() %>%
  filter(sample != "") %>%
  mutate(
    trimmed.reads = raw_single - trimmed_single,
    human.reads = trimmed_single - decontaminated_hg37dec_v0_1_single
  ) %>%
  mutate(percent_lowqual = trimmed.reads / raw_single) %>%
  mutate(percent_human = human.reads / raw_single) %>%
  mutate(recovery_rate = final_single / raw_single) %>%
  mutate(sample = str_replace(sample, "-clean", "")) %>%
  mutate(Study = "Routy") %>%
  mutate(clean.reads = final_single, total.reads = raw_single) %>% 
  filter(!is.na(sample))

# Plot distribution of reads
kd.qc.read_distribution <- kneaddata.data %>% 
  ggplot(aes(y = clean.reads, x = "Routy")) +
  geom_boxplot() + 
  labs(y = "Clean reads") +
  expand_limits(y = 0) +
  theme_bw()

kneaddata.report <- kneaddata.data %>%
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

kneaddata.data

# Load data ------------------------------------------------------------

routy.metadata.raw = read_csv(here(
  "data/studies/routy_science_2018/SraRunTable.csv"
)) 

routy.metadata = routy.metadata.raw %>%
  janitor::clean_names() %>%
  mutate(sample_id = run) %>%
  mutate(
    Response = if_else(
      response %in% c("Progression", "Dead"),
      "Non-responder",
      "Responder"
    ),
    BOR = case_when(
      response %in% c("Dead", "Progression") ~ "PD",
      response == "Partial response" ~ "PR",
      response == "Stable" ~ "SD"
    )
  ) %>%
  rename(atb = antibiotics) %>%
  column_to_rownames("sample_id") %>%
  mutate(Study = "Routy") %>%
  rename(Sample = run, sampleId = sample_name_2, gender = host_sex) %>%
  mutate(age_bin = cut(host_age, breaks = seq(10, 100, 5))) %>%
  mutate(age_bin = str_remove(age_bin, "\\(")) %>%
  mutate(age_bin = str_remove(age_bin, "\\]")) %>%
  mutate(age_bin = str_replace(age_bin, ",", "-")) %>%
  rename(age = host_age) %>%
  select(
    Sample,
    sampleId,
    Study,
    host_disease_status,
    Response,
    gender,
    BOR,
    age,
    age_bin,
    timepoint
  ) %>%
  mutate(Response = as.character(Response)) %>%
  mutate(
    timepoint = if_else(timepoint == "T0", "PreTreatment", "PostTreatment")
  ) %>% 
    mutate(therapy = 'anti-PD1') %>%
  left_join(
    kneaddata.data %>%
      select(sample, clean.reads, human.reads, total.reads) %>%
      rename(Sample = sample)
  )

rownames(routy.metadata) <- routy.metadata$Sample

# Make metaphlan4 object
routy.files = list.files(
  here("data/studies/routy_science_2018/results/metaphlan/read_stats/"),
  full.names = T
)
names(routy.files) = list.files(here(
  "data/studies/routy_science_2018/results/metaphlan/read_stats/"
)) %>%
  basename() %>%
  str_remove("-reads.txt")
routy.pseq = make_phyloseq4(table_fn = routy.files, sampledata = routy.metadata)
write_rds(routy.pseq, here("output/data/taxaprofiles/routy-metaphlan.rda"))

# Import Metacyc abundances
routy.metacyc = import_metacyc(
  here("data/studies/routy_science_2018/results/humann/pathway_abundance.txt"),
  routy.metadata,
  taxonomy = F
)
write_rds(routy.metacyc, here("output/data/metacyc/routy-metacyc.rda"))

# Import Metacyc abundances with taxonomy 
routy.metacyc.tax = import_metacyc(
  here("data/studies/routy_science_2018/results/humann/pathway_abundance.txt"),
  routy.metadata,
  taxonomy = T
)
write_rds(routy.metacyc.tax, here("output/data/regroup/routy-metacyc.tax.rda"))

# Import Metacyc abundances (cpm)
routy.metacyc.cpm = import_metacyc(
  here("data/studies/routy_science_2018/results/humann/pathway_abundance-cpm.txt"),
  routy.metadata,
  taxonomy = F
)
write_rds(routy.metacyc.cpm, here("output/data/metacyc/routy-metacyc-cpm.rda"))


# Import KEGG abundances 
routy.kegg = import_metacyc(
  here(
    "data/studies/routy_science_2018/results/humann/regroup/kegg_families.txt"
  ),
  routy.metadata,
  type = "kegg",
  taxonomy = F
)
write_rds(routy.kegg, here("output/data/regroup/routy-kegg.rda"))

# Import KEGG abundances with taxonomy 
routy.kegg.tax = import_metacyc(
  here(
    "data/studies/routy_science_2018/results/humann/regroup/kegg_families.txt"
  ),
  routy.metadata,
  type = "kegg",
  taxonomy = T
)
write_rds(routy.kegg.tax, here("output/data/regroup/routy-kegg.tax.rda"))

# Import eggNOG abundances
routy.eggnog = import_metacyc(
  here(
    "data/studies/routy_science_2018/results/humann/regroup/eggnog_families.txt"
  ),
  routy.metadata,
  type = "kegg",
  taxonomy = F
)
write_rds(routy.eggnog, here("output/data/regroup/routy-eggnog.rda"))

# Import KEGG-compound abundances
routy.compound = import_metacyc(
  here(
    "data/studies/routy_science_2018/results/humann/regroup/kegg_compounds.txt"
  ),
  routy.metadata
)
write_rds(routy.compound, here("output/data/regroup/routy-compound.rda"))

# Import KEGG-module abundances
routy.module = import_metacyc(
  here(
    "data/studies/routy_science_2018/results/humann/regroup/kegg_modules.txt"
  ),
  routy.metadata
)
write_rds(routy.module, here("output/data/regroup/routy-module.rda"))