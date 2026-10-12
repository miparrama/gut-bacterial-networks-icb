# Process of the Derosa dataset 
# @ Miguel Parra Martinez

here::here("scripts/studies/Derosa.R")
library(here)
source(here("scripts","helpers.R"))


# QC --------------------------------------------------------------------

## Kneaddata ------------------------------------------------------------

kneaddata.data = list.files(
  here("data/studies/derosa_nature_2022"),
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
  mutate(Study = "Derosa") %>%
  mutate(clean.reads = final_single, total.reads = raw_single) %>% 
  filter(!is.na(sample))

# Plot distribution of reads
kd.qc.read_distribution <- kneaddata.data %>% 
  ggplot(aes(y = clean.reads, x = "Derosa")) +
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

# Load data ------------------------------------------------------------

derosa.clinical = read_csv(here("data/studies/derosa_nature_2022/derosa-metadata.csv"))

derosa.metadata = read_csv(here(
  "data/studies/derosa_nature_2022/SraRunTable.csv"
)) %>%
  left_join(derosa.clinical, by = c("Sample Name" = "SampleID")) %>%
  janitor::clean_names() %>%
  mutate(sample_id = run) %>%
  column_to_rownames("sample_id") %>%
  filter(best_response %in% c("CR", "PR", "PD", "SD")) %>%
  mutate(
    Response = if_else(
      best_response %in% c("CR", "PR"),
      "Responder",
      "Non-responder"
    )
  ) %>%
  mutate(Response = fct_relevel(Response, "Non-responder")) %>%
  mutate(atb = if_else(atb == "wnno", "no", atb)) %>%
  mutate(atb = if_else(atb == "wynes", "yes", atb)) %>%
  mutate(host_disease_status = "NSCLC") %>%
  mutate(Study = "Derosa") %>%
  rename(gender = sex, Sample = run, sampleId = sample_name) %>%
  mutate(gender = if_else(gender == "Femme", "female", "male")) %>%
  mutate(age_bin = cut(age, breaks = seq(10, 100, 5))) %>%
  mutate(age_bin = str_remove(age_bin, "\\(")) %>%
  mutate(age_bin = str_remove(age_bin, "\\]")) %>%
  mutate(age_bin = str_replace(age_bin, ",", "-")) %>%
  mutate(os_event = os_evt) %>%
  mutate(BOR = best_response) %>% 
  select(
    Sample,
    sampleId,
    Study,
    host_disease_status,
    Response,
    BOR,
    gender,
    age,
    age_bin,
    atb,
    os_time,
    os_event
  ) %>%
  mutate(Response = as.character(Response)) %>%
  mutate(timepoint = "PreTreatment") %>%
  mutate(therapy = 'anti-PD1 / anti-PDL1') %>%
  left_join(
    kneaddata.data %>%
      select(sample, clean.reads, human.reads, total.reads) %>%
      rename(Sample = sample)
  )
rownames(derosa.metadata) <- derosa.metadata$Sample

# Make metaphlan4 object
derosa.files = list.files(
  here("data/studies/derosa_nature_2022/results/metaphlan/read_stats/"),
  full.names = T
)
names(derosa.files) = list.files(
  here("data/studies/derosa_nature_2022/results/metaphlan/read_stats/")
) %>%
  basename() %>%
  str_remove("-reads.txt")

derosa.pseq = make_phyloseq4(
  table_fn = derosa.files,
  sampledata = derosa.metadata
)
write_rds(derosa.pseq, here("output/data/taxaprofiles/derosa-metaphlan.rda"))

# Import Metacyc abundances
derosa.metacyc = import_metacyc(
  here("data/studies/derosa_nature_2022/results/humann/pathway_abundance.txt"),
  derosa.metadata,
  taxonomy = F
)
write_rds(derosa.metacyc, here("output/data/metacyc/derosa-metacyc.rda"))

# Import Metacyc abundances with taxonomy 
derosa.metacyc.tax = import_metacyc(
  here("data/studies/derosa_nature_2022/results/humann/pathway_abundance.txt"),
  derosa.metadata,
  taxonomy = T
)
write_rds(derosa.metacyc.tax, here("output/data/regroup/derosa-metacyc.tax.rda"))

# Import Metacyc abundances (cpm)
derosa.metacyc.cpm = import_metacyc(
  here("data/studies/derosa_nature_2022/results/humann/pathway_abundance-cpm.txt"),
  derosa.metadata,
  taxonomy = F
)
write_rds(derosa.metacyc.cpm, here("output/data/metacyc/derosa-metacyc-cpm.rda"))

# Import KEGG abundances
derosa.kegg = import_metacyc(
  here("data/studies/derosa_nature_2022/results/humann/regroup/kegg_families.txt"),
  derosa.metadata,
  type = "kegg",
  taxonomy = F
)
write_rds(derosa.kegg, here("output/data/regroup/derosa-kegg.rda"))

# Import KEGG abundances with taxonomy 
derosa.kegg.tax = import_metacyc(
  here("data/studies/derosa_nature_2022/results/humann/regroup/kegg_families.txt"),
  derosa.metadata,
  type = "kegg",
  taxonomy = T
)
write_rds(derosa.kegg.tax, here("output/data/regroup/derosa-kegg.tax.rda"))

# Import eggNOG abundances
derosa.eggnog = import_metacyc(
  here("data/studies/derosa_nature_2022/results/humann/regroup/eggnog_families.txt"),
  derosa.metadata,
  type = "kegg",
  taxonomy = F
)
write_rds(derosa.eggnog, here("output/data/regroup/derosa-eggnog.rda"))

# Import KEGG-compound abundances
derosa.compound = import_metacyc(
  here("data/studies/derosa_nature_2022/results/humann/regroup/kegg_compounds.txt"),
  derosa.metadata
)
write_rds(derosa.compound, here("output/data/regroup/derosa-compound.rda"))

# Import KEGG-module abundances
derosa.module = import_metacyc(
  here("data/studies/derosa_nature_2022/results/humann/regroup/kegg_modules.txt"),
  derosa.metadata
)
write_rds(derosa.module, here("output/data/regroup/derosa-module.rda"))


# RNAseq: ------------------------------------------------------------
# ! The overlap between RNAseq samples and metagenomics samples is not clear. This would have to be checked with the authors of the manuscript.
# ! The matching done here was done by Tom and now it is not clear if it is correct. It is based on the last numbers of the SampleId of the RNAseq data, and the `Sample Name` from the metagenomics data. LUM3_617 would match with S12617. But S39013 can match with LUM2_013 or IMM_013. 

derosa.clinical = read_csv(here("data/studies/derosa_nature_2022/derosa-metadata.csv"))

derosa.counts = read_tsv(here("data/studies/derosa_nature_2022/transcriptome/GSE182328_Gene_counts_matrix.txt")) %>% 
  column_to_rownames("Gene")

derosa.project = read_csv(here("data/studies/derosa_nature_2022/transcriptome/SraRunTable.csv")) %>% 
  janitor::clean_names() %>% 
  select(run, akk_metaominer, disease_state)

# Complete metadata
derosa.temp = derosa.clinical %>% 
  select(SampleID, akksgb, akkmsp) %>% 
  separate(SampleID, into = c("temp", "sample_id"), remove = F) %>% 
  select(-temp) %>% 
  mutate(sample_id = as.numeric(sample_id)) %>% 
  mutate(sample_id = as.factor(sample_id))

derosa.temp %>% count(sample_id, sort = T)

# Aggregate metadata  
derosa.meta = read_tsv(here("data/studies/derosa_nature_2022/transcriptome/sra_explorer_metadata.tsv")) %>% 
  janitor::clean_names() %>% 
  rename(NCBID = accession) %>% 
  separate(sra_nice_filename, into = c("SampleId", "GSM", "Sample", "type", "timepoint")) %>% 
  select(NCBID, Sample) %>% 
  filter(Sample %in% colnames(derosa.counts)) %>% 
  distinct(NCBID, .keep_all = T) %>% 
  left_join(derosa.project, by = c("NCBID" = "run")) %>% 
  mutate(sample_id = substring(Sample, 4)) %>% 
  mutate(sample_id = as.numeric(sample_id)) %>% 
  mutate(sample_id = as.factor(sample_id)) %>% 
  filter(sample_id != 56) %>% 
  left_join(derosa.temp) %>% 
  filter(SampleID != "LUM2_097") %>% 
  filter(SampleID != "LUM3_101")

# Convert gene names to symbols
library(IOBR)
derosa.rna.tpm <- count2tpm(countMat = derosa.counts, idType = "SYMBOL", org = "hsa")

# Get immune signatures
derosa.abis = immunedeconv::deconvolute(derosa.rna.tpm, "abis", tumor = T) %>% 
  pivot_longer(cols = -cell_type, names_to = "sampleId") %>% 
  pivot_wider(names_from = "cell_type", values_from = "value") %>% 
  janitor::clean_names() %>% 
  dplyr::rename(sampleId = sample_id) %>% 
  # left_join(baruch.rna.metadata, by = c("sampleId" = "Sample")) %>% 
  left_join(select(derosa.meta, Sample, SampleID), by = c("sampleId" = "Sample")) %>% 
  write_csv(here("output/data/transcriptome/derosa-abis.csv"))

library(DGEobj.utils)
RawCPM.log2 <- DGEobj.utils::convertCounts(derosa.counts %>% as.matrix(), unit = "CPM", log = F, normalize = "none", prior.count = 0) %>% 
  as.data.frame() 
RawCPM.log2 = log2(RawCPM.log2 + 1)

RawCPM.log2.long = RawCPM.log2 %>%
  rownames_to_column("symbol") %>% 
  pivot_longer(-symbol, names_to = "sampleId", values_to = "cpm") 

mait_yao_et_al = c("SLC4A10", "IL23R", "NCR3", "TMIGD2", "LST1", "COLQ", "ME1", "ZBTB16", "RORC", "TLE1")

# Import marker gene list (from Joris)
marker.gene.signatures = readxl::read_excel(here("data/resources/Danaher_marker_genes_extended.xlsx")) %>% 
  select(Gene, Marker_gene_set) %>% 
  add_row(Gene = mait_yao_et_al, Marker_gene_set = rep("MAIT", length(mait_yao_et_al)) ) 

# Make genes into list
gene.signatures = split(marker.gene.signatures$Gene, marker.gene.signatures[, 'Marker_gene_set'])

# Iterate over signatures and determine significance 
immune.signatures.cpm = gene.signatures %>% 
  map_dfr(.id = "signature", function(x){
    RawCPM.log2.long %>% 
      as_tibble() %>% 
      filter(symbol %in% x) %>% 
      group_by(sampleId) %>% 
      summarise(avg = mean(cpm)) 
  }, .progress = TRUE) %>% 
  pivot_wider(names_from = signature, values_from = avg) %>% 
  janitor::clean_names() %>% 
  rename(sampleId = sample_id) %>% 
  write_csv(here("output/data/transcriptome/derosa-immune-signatures.csv"))


