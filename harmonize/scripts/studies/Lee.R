# Process of the Lee dataset
# @ Miguel Parra Martinez

here::here("scripts/studies/Lee.R")
library(here)
source(here("scripts", "helpers.R"))

# QC --------------------------------------------------------------------

## Kneaddata ------------------------------------------------------------

kneaddata.data = list.files(
  here("data/studies/lee_nature_2022/"),
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
  mutate(Study = "Lee") %>%
  mutate(sample = str_replace(sample, "-clean", "")) %>%
  mutate(sample = str_replace(sample, "-kd", "")) %>%
  mutate(clean.reads = final, total.reads = raw) 

# Plot distribution of reads
kd.qc.read_distribution <- kneaddata.data %>% 
  ggplot(aes(y = clean.reads, x = "Lee")) +
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

lee.metadata = read_csv(here("data/studies/lee_nature_2022/SraRunTable.csv"))

lee.metadata = lee.metadata %>%
  janitor::clean_names() %>%
  mutate(sample_id = run) %>%
  column_to_rownames("sample_id") %>%
  mutate(
    PFS12_response = if_else(
      pfs12 %in% c("yes", "no"),
      "Responder",
      "Non-responder"
    )
  ) %>%
  mutate(
    Response = if_else(orr %in% c("CR", "PR"), "Responder", "Non-responder")
  ) %>%
  mutate(host_disease_status = "MELA") %>%
  mutate(Study = "Lee") %>%
  mutate(bmi = host_body_mass_index) %>%
  mutate(BOR = orr) %>% 
  rename(
    gender = host_sex,
    sampleId = sample_name_2,
    age_bin = host_age_5yr_bin,
    Sample = run
  ) %>%
  mutate(age = str_sub(age_bin, 1, 2) %>% as.numeric()) %>%
  select(
    Sample,
    sampleId,
    Study,
    host_disease_status,
    Response,
    BOR,
    cohort,
    gender,
    age,
    age_bin,
    bmi
  ) %>%
  mutate(Response = as.character(Response)) %>%
  mutate(timepoint = "PreTreatment") %>%
  mutate(therapy = 'anti-PD1 / anti-CTL4/ anti-PD1 + anti-CTL4') %>%
  left_join(
    kneaddata.data %>%
      select(sample, clean.reads, human.reads, total.reads) %>%
      rename(Sample = sample)
  )

rownames(lee.metadata) <- lee.metadata$Sample

# Make metaphlan4 object
# lee.files = list.files(
#   here("data/studies/lee_nature_2022/lee_nature_2022_v2/results/metaphlan/read_stats/"),
#   full.names = T
# )
# names(lee.files) = list.files(
#   here("data/studies/lee_nature_2022/lee_nature_2022_v2/results/metaphlan/read_stats/")
# ) %>%
#   basename() %>%
#   str_remove("-reads.txt")
# lee.pseq = make_phyloseq4(table_fn = lee.files, sampledata = lee.metadata)
# write_rds(lee.pseq, here("output/data/taxaprofiles/lee-metaphlan.rda"))

# Make metaphlan4 object
lee.files = list.files(
  here("data/studies/lee_nature_2022/lee_nature_2022_v2/results/metaphlan/read_stats/"),
  full.names = T
)
names(lee.files) = list.files(
  here("data/studies/lee_nature_2022/lee_nature_2022_v2/results/metaphlan/read_stats/")
) %>%
  basename() %>%
  str_remove("-reads.txt")
lee.pseq = make_phyloseq4(table_fn = lee.files, sampledata = lee.metadata)
write_rds(lee.pseq, here("output/data/taxaprofiles/lee-metaphlan.rda"))

# Import Metacyc abundances
lee.metacyc = import_metacyc(
  here("data/studies/lee_nature_2022/lee_nature_2022_v2/results/humann/pathway_abundance.txt"),
  lee.metadata,
  taxonomy = F
)
write_rds(lee.metacyc, here("output/data/metacyc/lee-metacyc.rda"))


# Import Metacyc abundances with taxonomy 
lee.metacyc.tax = import_metacyc(
  here("data/studies/lee_nature_2022/lee_nature_2022_v2/results/humann/pathway_abundance.txt"),
  lee.metadata,
  taxonomy = T
)
write_rds(lee.metacyc.tax, here("output/data/regroup/lee-metacyc.tax.rda"))

# Import Metacyc abundances (cpm)
lee.metacyc.cpm = import_metacyc(
  here("data/studies/lee_nature_2022/lee_nature_2022_v2/results/humann/pathway_abundance-cpm.txt"),
  lee.metadata,
  taxonomy = F
)
write_rds(lee.metacyc.cpm, here("output/data/metacyc/lee-metacyc-cpm.rda"))


# Import KEGG abundances
lee.kegg = import_metacyc(
  here("data/studies/lee_nature_2022/lee_nature_2022_v2/results/humann/regroup/kegg_families.txt"),
  lee.metadata,
  type = "kegg",
  taxonomy = F
)
write_rds(lee.kegg, here("output/data/regroup/lee-kegg.rda"))

# Import KEGG abundances with taxonominy 
lee.kegg.tax = import_metacyc(
  here("data/studies/lee_nature_2022/lee_nature_2022_v2/results/humann/regroup/kegg_families.txt"),
  lee.metadata,
  type = "kegg",
  taxonomy = T
)
write_rds(lee.kegg.tax, here("output/data/regroup/lee-kegg.tax.rda"))

# Import eggNOG abundances
lee.eggnog = import_metacyc(
  here("data/studies/lee_nature_2022/lee_nature_2022_v2/results/humann/regroup/eggnog_families.txt"),
  lee.metadata,
  type = "kegg",
  taxonomy = F
)
write_rds(lee.eggnog, here("output/data/regroup/lee-eggnog.rda"))

# Import KEGG-compound abundances
lee.compound = import_metacyc(
  here("data/studies/lee_nature_2022/lee_nature_2022_v2/results/humann/regroup/kegg_compounds.txt"),
  lee.metadata
)
write_rds(lee.compound, here("output/data/regroup/lee-compound.rda"))

# Import KEGG-module abundances
lee.module = import_metacyc(
  here("data/studies/lee_nature_2022/lee_nature_2022_v2/results/humann/regroup/kegg_modules.txt"),
  lee.metadata
)
write_rds(lee.module, here("output/data/regroup/lee-module.rda"))



