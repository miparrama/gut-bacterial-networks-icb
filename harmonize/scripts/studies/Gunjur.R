# Process of the gunjur dataset
# @ Miguel Parra Martinez

here::here("scripts/studies/Gunjur.R")
library(here)
source(here("scripts", "helpers.R"))

# QC --------------------------------------------------------------------

## Kneaddata ------------------------------------------------------------

kneaddata.data = list.files(
  here("data/studies/gunjur_natmed_2024/"),
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
  mutate(Study = "Gunjur") %>%
  filter(!is.na(sample)) %>% 
  mutate(clean.reads = final, total.reads = raw) 

# Plot distribution of reads
kd.qc.read_distribution <- kneaddata.data %>% 
  ggplot(aes(y = clean.reads, x = "Gunjur")) +
  geom_boxplot() + 
  labs(y = "Clean reads") +
  expand_limits(y = 0) +
  theme_bw()


kneaddata.report = kneaddata.data %>%
  # gather(Metric, counts, -c(sample, Study)) %>%
  pivot_longer(
    cols = -c(sample, Study),
    names_to = "Metric",
    values_to = "counts"
  ) %>%
  filter(
    Metric %in%
      c("percent_lowqual", "percent_human", "recovery_rate", "Study")
  )

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



#  kneaddata.data %>% 
#   filter(percent_lowqual > 0.75) %>% 
#   filter(total.reads < 5e+06)

# ERR12405951

# ERR12405935


# kneaddata.data %>% 
#   filter(sample == "ERR12405935")


# Filter samples with low QC. 
kneaddata.data <- kneaddata.data %>% 
  filter(percent_lowqual < 0.75) %>% 
  filter(total.reads > 5e+06)



## Humann ------------------------------------------------------------

# Overview of counts
humann.logs = list.files(
  here("data/studies/gunjur_natmed_2024"),
  pattern = "-report.txt",
  recursive = T,
  full.names = T
)
names(humann.logs) = list.files(
  here("data/studies/gunjur_natmed_2024"),
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

# Filter sample with low reads.

humann.report %>% 
  arrange(total_reads)

humann.report %>%
  arrange(total_species)

humann.report %>%
  arrange(nucleotide_percent)

# ERR12405935 has very low nucleotide_percent assignation and very low number of species detected compared with the rest of the samples (25)
# ERR12405951 has very few reads ~ 8000.


# Load data ------------------------------------------------------------

gunjur.metadata.info = read_xlsx(
  here("data/studies/gunjur_natmed_2024/41591_2024_2823_MOESM3_ESM.xlsx"),
  sheet = 4
)

samplesheet = read_csv(
  here("data/studies/gunjur_natmed_2024/ena_metadata.csv")
) %>% 
  mutate(lane_id = str_replace(run_alias, "SC_RUN_", "")) %>% 
  relocate(lane_id)

gunjur.metadata = read_xlsx(
  here("data/studies/gunjur_natmed_2024/41591_2024_2823_MOESM3_ESM.xlsx"),
  sheet = 4
) %>%
  left_join(samplesheet) %>%
  mutate(
    sample_id = run_accession,
    Sample = run_accession,
    sampleId = run_accession
  ) %>%
  column_to_rownames("sample_id") %>%
  mutate(
    Response = if_else(
      BOR %in% c("CR", "PR"),
      "Responder",
      "Non-responder"
    )
  ) %>%
  mutate(Response = fct_relevel(Response, "Non-responder")) %>%
  mutate(host_disease_status = histology_cohort) %>%
  mutate(Study = "Gunjur") %>%
  mutate(gender = if_else(sex == "F", "female", "male")) %>%
  mutate(age = as.numeric(age)) %>%
  mutate(
    age_bin = cut(age, breaks = seq(10, 100, 5)) %>%
      str_remove("\\(") %>%
      str_remove("\\]") %>%
      str_replace(",", "-")
  ) %>%
  mutate(atb = ifelse(antibiotic_use == TRUE, "yes", "no")) %>%
  mutate(bmi = BMI) %>% 
  mutate(
    pfs_event = ifelse(not_progression_free == TRUE, "yes", "no"),
    pfs_time = PFS,
    os_time = OS,
    os_event = ifelse(deceased == TRUE, "yes", "no")
  ) %>%
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
    bmi,
    atb,
    pfs_event,
    pfs_time,
    os_time,
    os_event
  ) %>%
  mutate(Response = as.character(Response)) %>%
  mutate(timepoint = "PreTreatment") %>%
  filter(!Sample %in% c("ERR12405935", "ERR12405951")) %>% #QC filter
  mutate(therapy = 'anti-PD1 + anti-CTL4') %>%
  inner_join(
    kneaddata.data %>%
      select(sample, clean.reads, human.reads, total.reads) %>%
      rename(Sample = sample)
  )
rownames(gunjur.metadata) <- gunjur.metadata$Sample

# Make metaphlan4 object
gunjur.files = list.files(
  here("data/studies/gunjur_natmed_2024/results/metaphlan/read_stats/"),
  full.names = T
)
names(gunjur.files) = list.files(
  here("data/studies/gunjur_natmed_2024/results/metaphlan/read_stats/")
) %>%
  basename() %>%
  str_remove("-reads.txt") 

gunjur.pseq = make_phyloseq4(
  table_fn = gunjur.files,
  sampledata = gunjur.metadata
)
write_rds(gunjur.pseq, here("output/data/taxaprofiles/gunjur-metaphlan.rda"))


# Import Metacyc abundances
gunjur.metacyc = import_metacyc(
  here("data/studies/gunjur_natmed_2024/results/humann/pathway_abundance.txt"),
  gunjur.metadata,
  taxonomy = F
)
write_rds(gunjur.metacyc, here("output/data/metacyc/gunjur-metacyc.rda"))

# Import Metacyc abundances with taxonomy 
gunjur.metacyc.tax = import_metacyc(
  here("data/studies/gunjur_natmed_2024/results/humann/pathway_abundance.txt"),
  gunjur.metadata,
  taxonomy = T
)
write_rds(gunjur.metacyc.tax, here("output/data/regroup/gunjur-metacyc.tax.rda"))

# Import Metacyc abundances (cpm)
gunjur.metacyc.cpm = import_metacyc(
  here("data/studies/gunjur_natmed_2024/results/humann/pathway_abundance-cpm.txt"),
  gunjur.metadata,
  taxonomy = F
)
write_rds(gunjur.metacyc.cpm, here("output/data/metacyc/gunjur-metacyc-cpm.rda"))

# Import KEGG abundances
gunjur.kegg = import_metacyc(
  here(
    "data/studies/gunjur_natmed_2024/results/humann/regroup/kegg_families.txt"
  ),
  gunjur.metadata,
  type = "kegg",
  taxonomy = F
)
write_rds(gunjur.kegg, here("output/data/regroup/gunjur-kegg.rda"))

# Import KEGG abundances with taxonomy 
gunjur.kegg.tax = import_metacyc(
  here("data/studies/gunjur_natmed_2024/results/humann/regroup/kegg_families.txt"),
  gunjur.metadata,
  type = "kegg",
  taxonomy = T
)
write_rds(gunjur.kegg.tax, here("output/data/regroup/gunjurr-kegg.tax.rda"))

# Import eggNOG abundances
gunjur.eggnog = import_metacyc(
  here(
    "data/studies/gunjur_natmed_2024/results/humann/regroup/eggnog_families.txt"
  ),
  gunjur.metadata,
  type = "kegg",
  taxonomy = F
)
write_rds(gunjur.eggnog, here("output/data/regroup/gunjur-eggnog.rda"))

# Import KEGG-compound abundances
gunjur.compound = import_metacyc(
  here(
    "data/studies/gunjur_natmed_2024/results/humann/regroup/kegg_compounds.txt"
  ),
  gunjur.metadata
)
write_rds(gunjur.compound, here("output/data/regroup/gunjur-compound.rda"))

# Import KEGG-module abundances
gunjur.module = import_metacyc(
  here(
    "data/studies/gunjur_natmed_2024/results/humann/regroup/kegg_modules.txt"
  ),
  gunjur.metadata
)
write_rds(gunjur.module, here("output/data/regroup/gunjur-module.rda"))
