# Process of the Kim dataset
# @ Miguel Parra Martinez

here::here("scripts/studies/Kim.R")
library(here)
source(here("scripts", "helpers.R"))

# QC --------------------------------------------------------------------

## Kneaddata ------------------------------------------------------------

kneaddata.data = list.files(
  here("data/studies/kim_natcom_2024/results-kim"),
  pattern = "*-kneaddata.txt",
  recursive = T,
  full.names = T
) %>%
  map_dfr(data.table::fread) %>%
  filter(!is.na(Sample)) %>%
  janitor::clean_names() %>%
  filter(sample != "") %>%
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
  mutate(Study = "Kim") %>%
  filter(!is.na(sample)) %>%
  mutate(sample = str_replace(sample, "-kd", "")) %>%
  mutate(clean.reads = final, total.reads = raw) 

# Plot distribution of reads
kd.qc.read_distribution <- kneaddata.data %>% 
  ggplot(aes(y = clean.reads, x = "Kim_2nd")) +
  geom_boxplot() + 
  labs(y = "Clean reads") +
  expand_limits(y = 0) +
  theme_bw()


kneaddata.report <- kneaddata.data %>% 
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


## Humann ------------------------------------------------------------

# Overview of counts
humann.logs = list.files(
  here("data/studies/kim_natcom_2024/results-kim"),
  pattern = "-report.txt",
  recursive = T,
  full.names = T
)
names(humann.logs) = list.files(
  here("data/studies/kim_natcom_2024/results-kim"),
  pattern = "-report.txt",
  recursive = T,
  full.names = F
) %>%
  basename() %>%
  str_remove("-report.txt")

humann.logs.df = humann.logs %>%
  map_dfr(.f = data.table::fread, .id = "Sample") %>%
  janitor::clean_names() %>%
  filter(!is.na(total_reads)) %>%
  rename(Sample = sample) %>%
  as_tibble()

humann.report = humann.logs.df %>%
  na.omit() %>%
  janitor::clean_names() %>%
  rename(Sample = number_samples) %>%
  mutate(
    total_assigned = total_reads /
      (total_nucleotide_aligned + total_translated_aligned)
  ) %>%
  mutate(nucleotide_percent = total_nucleotide_aligned / total_reads) %>%
  mutate(translated = total_translated_aligned / total_reads) %>%
  mutate(unassigned_percent = 1 - translated) %>%
  mutate(translated_percent = 1 - (unassigned_percent + nucleotide_percent))

humann.report %>%
  rowwise() %>%
  mutate(sums = sum(nucleotide_percent, translated_percent)) %>%
  pull(sums) %>%
  sd()

p1 = humann.report %>%
  gather(Metric, percent, -Sample) %>%
  filter(
    Metric %in%
      c("nucleotide_percent", "unassigned_percent", "translated_percent")
  ) %>%
  mutate(
    Metric = fct_recode(
      Metric,
      Nucleotide = "nucleotide_percent",
      Translated = "translated_percent",
      `Unassigned reads` = "unassigned_percent"
    )
  ) %>%
  mutate(percent = as.numeric(percent)) %>%
  ggplot(aes(x = Sample, y = percent, fill = Metric)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_blank(), legend.position = "bottom") +
  # jcolors::scale_fill_jcolors(palette = "default") +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
  xlab("Sample") +
  ylab("Percent of sequenced reads")

# Plot boxplot side-by-side
p2 = humann.report %>%
  gather(Metric, percent, -Sample) %>%
  filter(
    Metric %in%
      c("nucleotide_percent", "unassigned_percent", "translated_percent")
  ) %>%
  mutate(
    Metric = fct_recode(
      Metric,
      Nucleotide = "nucleotide_percent",
      Translated = "translated_percent",
      `Unassigned \n reads` = "unassigned_percent"
    )
  ) %>%
  mutate(percent = as.numeric(percent)) %>%
  ggplot(aes(x = Metric, y = percent, fill = Metric)) +
  stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.15) +
  geom_boxplot(width = 0.5) +
  guides(fill = "none") +
  theme_few() +
  scale_fill_brewer(palette = "Paired") +
  #jcolors::scale_fill_jcolors(palette = "default") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
  xlab("") +
  ylab("")

# Plot together
p = (p1 + p2) +
  plot_layout(widths = c(3, 1), guides = "collect") &
  theme(legend.position = 'bottom')
p


# Load data ------------------------------------------------------------

kim.metadata = read_csv(
  here("data/studies/kim_natcom_2024/SraRunTable.csv")
)

kim.metadata = read_csv(
  here("data/studies/kim_natcom_2024/SraRunTable.csv")
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
  mutate(Study = "Kim_2nd") %>%
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
  mutate(therapy = 'anti-PDL1') %>%
  left_join(
    kneaddata.data %>%
      select(sample, clean.reads, human.reads, total.reads) %>%
      rename(Sample = sample)
  )
rownames(kim.metadata) <- kim.metadata$Sample

# Make metaphlan4 object
kim.files = list.files(
  here("data/studies/kim_natcom_2024/results-kim/metaphlan/read_stats/"),
  full.names = T
)
names(kim.files) = list.files(
  here("data/studies/kim_natcom_2024/results-kim/metaphlan/read_stats/")
) %>%
  basename() %>%
  str_remove("-reads.txt")

kim.pseq = make_phyloseq4(
  table_fn = kim.files,
  sampledata = kim.metadata
)
write_rds(kim.pseq, here("output/data/taxaprofiles/kim-metaphlan.rda"))

# Import Metacyc abundances
kim.metacyc = import_metacyc(
  here("data/studies/kim_natcom_2024/results-kim/humann/pathway_abundance.txt"),
  kim.metadata,
  taxonomy = F
)
write_rds(kim.metacyc, here("output/data/metacyc/kim-metacyc.rda"))

# Import Metacyc abundances with taxonomy 
kim.metacyc.tax = import_metacyc(
  here("data/studies/kim_natcom_2024/results-kim/humann/pathway_abundance.txt"),
  kim.metadata,
  taxonomy = T
)
write_rds(kim.metacyc.tax, here("output/data/regroup/kim-metacyc.tax.rda"))

# Import Metacyc abundances (cpm)
kim.metacyc.cpm = import_metacyc(
  here("data/studies/kim_natcom_2024/results-kim/humann/pathway_abundance-cpm.txt"),
  kim.metadata,
  taxonomy = F
)
write_rds(kim.metacyc.cpm, here("output/data/metacyc/kim-metacyc-cpm.rda"))

# Import KEGG abundances
kim.kegg = import_metacyc(
  here(
    "data/studies/kim_natcom_2024/results-kim/humann/regroup/kegg_families.txt"
  ),
  kim.metadata,
  type = "kegg",
  taxonomy = F
)
write_rds(kim.kegg, here("output/data/regroup/kim-kegg.rda"))

# Import KEGG abundances with taxonomy 
kim.kegg.tax = import_metacyc(
  here("data/studies/kim_natcom_2024/results-kim/humann/regroup/kegg_families.txt"),
  kim.metadata,
  type = "kegg",
  taxonomy = T
)
write_rds(kim.kegg.tax, here("output/data/regroup/kim-kegg.tax.rda"))

# Import eggNOG abundances
kim.eggnog = import_metacyc(
  here(
    "data/studies/kim_natcom_2024/results-kim/humann/regroup/eggnog_families.txt"
  ),
  kim.metadata,
  type = "kegg",
  taxonomy = F
)
write_rds(kim.eggnog, here("output/data/regroup/kim-eggnog.rda"))

# Import KEGG-compound abundances
kim.compound = import_metacyc(
  here(
    "data/studies/kim_natcom_2024/results-kim/humann/regroup/kegg_compounds.txt"
  ),
  kim.metadata
)
write_rds(kim.compound, here("output/data/regroup/kim-compound.rda"))

# Import KEGG-module abundances
kim.module = import_metacyc(
  here(
    "data/studies/kim_natcom_2024/results-kim/humann/regroup/kegg_modules.txt"
  ),
  kim.metadata
)
write_rds(kim.module, here("output/data/regroup/kim-module.rda"))