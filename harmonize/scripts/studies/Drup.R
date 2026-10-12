# Process of the DRUP dataset
# @ Miguel Parra Martinez

here::here("scripts/studies/DRUP.R")
library(here)
source(here("scripts", "helpers.R"))

# Load QC for DRUP. ------------------------------------------------------------

kneaddata.data = list.files(
  here("data/studies/drup/biobakery/results/"),
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
  mutate(Study = "DRUP") %>%
  mutate(sample = str_replace(sample, "-clean", "")) %>%
  mutate(sample = str_replace(sample, "-kd", "")) %>%
  mutate(clean.reads = final, total.reads = raw) 

write_csv(kneaddata.data, file = here("output/data/QC/drup.kneaddata.csv"))

# Plot distribution of reads
kd.qc.read_distribution <- kneaddata.data %>% 
  ggplot(aes(y = clean.reads, x = "DRUP")) +
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

# Load TSE object from DRUP-MSI dataset. 
tse_reads <- read_rds("~/projects/gut-microbiome/gut-metabolomics/output/data/metagenomics/drup/tse_reads_metaphlan.rds")


metadata_drup <- coldf(tse_reads) 

tse_drup_coldata_df <- metadata_drup %>%
  mutate(
    sampleId = Sample,
    Study = "DRUP",
    host_disease_status = tumor_class,
    Response = ifelse(Response == "R", "Responder", "Non-responder"),
    gender = sex,
    BOR = BOR_proxy,
    age = Age,
    age_bin = cut(age, breaks = seq(10, 100, 5)) %>%
      str_remove("\\(") %>%
      str_remove("\\]") %>%
      str_replace(",", "-"),
    timepoint = ifelse(timepoint == "T1", "PreTreatment", "PostTreatment"), 
    atb = ifelse(ATB_use == "Yes", "yes", "no"),
    cohort = "MSI",
    batch = "DRUP-MSI", 
    pfs_event = PFS_event,
    pfs_time = PFS_time,
    os_time = OS_time,
    os_event = OS_event,
  ) %>% 
  left_join(drup_bmi) %>% 
    select(
    Sample,
    sampleId,
    Study,
    host_disease_status,
    Response,
    cohort,
    gender,
    BOR,
    age,
    age_bin,
    timepoint,
    bmi,
    atb,
    pfs_time, 
    pfs_event,
    os_time,
    os_event,
  ) %>% 
  mutate(therapy = 'anti-PD1') %>%
  left_join(
    kneaddata.data %>%
      select(sample, clean.reads, human.reads, total.reads) %>%
      rename(Sample = sample)
  ) %>% 
  DataFrame()

rownames(tse_drup_coldata_df) <- tse_drup_coldata_df$Sample

colData(tse_reads) <- tse_drup_coldata_df

drup.pseq <- convertToPhyloseq(tse_reads)

write_rds(drup.pseq, here("output/data/taxaprofiles/drup-metaphlan.rda"))


### Functional data ------------------------------------------------------------

metadata_drup_df <- tse_drup_coldata_df %>% 
  as.data.frame() %>% 
  mutate(Sample2 = Sample)

# Import Metacyc abundances
drup.metacyc = import_metacyc(
  here("data/studies/drup/biobakery/results/humann/pathway_abundance.txt"),
  metadata_drup_df,
  taxonomy = F
)
write_rds(drup.metacyc, here("output/data/metacyc/drup-metacyc.rda"))

# Import Metacyc abundances with taxonomy 
drup.metacyc.tax = import_metacyc(
  here("data/studies/drup/biobakery/results/humann/pathway_abundance.txt"),
  metadata_drup_df,
  taxonomy = T
)
write_rds(drup.metacyc.tax, here("output/data/regroup/drup-metacyc.tax.rda"))


drup.metacyc.cpm = import_metacyc(
  here("data/studies/drup/biobakery/results/humann/pathway_abundance-cpm.txt"),
  metadata_drup_df,
  taxonomy = F
)
write_rds(drup.metacyc.cpm, here("output/data/metacyc/drup-metacyc-cpm.rda"))

# Import KEGG abundances
drup.kegg = import_metacyc(
  here("data/studies/drup/biobakery/results/humann/regroup/kegg_families.txt"),
  metadata_drup_df,
  type = "kegg",
  taxonomy = F
)
write_rds(drup.kegg, here("output/data/regroup/drup-kegg.rda"))

# Import KEGG abundances with taxonomy 
drup.kegg.tax = import_metacyc(
  here("data/studies/drup/biobakery/results/humann/regroup/kegg_families.txt"),
  metadata_drup_df,
  type = "kegg",
  taxonomy = T
)
write_rds(drup.kegg.tax, here("output/data/regroup/drup-kegg.tax.rda"))

# Import KEGG abundances normalized
drup.kegg.relab = import_metacyc(
  here("data/studies/drup/biobakery/results/humann/regroup/kegg_families-relab.txt"),
  metadata_drup_df,
  type = "kegg",
  taxonomy = F
)
write_rds(drup.kegg.relab, here("output/data/regroup/drup-kegg-relab.rda"))

# Import eggNOG abundances
drup.eggnog = import_metacyc(
  here("data/studies/drup/biobakery/results/humann/regroup/eggnog_families.txt"),
  metadata_drup_df,
  type = "kegg",
  taxonomy = F
)
write_rds(drup.eggnog, here("output/data/regroup/drup-eggnog.rda"))

# Import KEGG-compound abundances
drup.compound = import_metacyc(
  here("data/studies/drup/biobakery/results/humann/regroup/kegg_compounds.txt"),
  metadata_drup_df
)
write_rds(drup.compound, here("output/data/regroup/drup-compound.rda"))

# Import KEGG-module abundances
drup.module = import_metacyc(
  here("data/studies/drup/biobakery/results/humann/regroup/kegg_modules.txt"),
  metadata_drup_df
)
write_rds(drup.module, here("output/data/regroup/drup-module.rda"))