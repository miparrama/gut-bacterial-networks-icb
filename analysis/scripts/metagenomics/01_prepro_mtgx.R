# Preprocessing of metagenomics data. 
# @ Miguel Parra Martinez
# 2025-01

source("scripts/helpers.R")

# Load metadata -------------------------------------------------------------------

# Import metadata.
drup.gut.metadata <- read_csv(
  "data/metagenomics/DRUP/drup-gut-metadata.csv"
) %>%
  select(
    -any_of(c(
      #Filter old clinical annotations.
      "bor.DRUP",
      "bor.according.to.RECIST",
      "bor.DRUP",
      "drup.clin.benefit",
      "pfsdate",
      "pfsweeks",
      "pfsevent",
      "death.date",
      "osdate",
      "osweeks",
      "osevent",
      "CB",
      "BOR",
      "Response"
    ))
  )

#Import clinical data
clinical_data <- read_csv("data/metagenomics/DRUP/clinical_data_elaborate_20250314.csv") %>% 
  mutate(patientId = str_replace_all(SubjectKey, "-", ""))

#Import dysbiosis data (old dysbiosis data)
functional_dysbiosis <- read_csv ("data/metagenomics/DRUP/drup-functional-dysbiosis.csv") %>% 
  select (Sample, dysbiosis) %>% rename (functional_dysbiosis = dysbiosis)
species_dysbiosis <- read_csv ("data/metagenomics/DRUP/drup-species-dysbiosis.csv") %>% 
  select (Sample, dysbiosis) %>% rename (species_dysbiosis = dysbiosis)

#Import metacyc pathways names
metacyc.names = data.table::fread("data/resources/map_metacyc-pwy_name.txt.gz", header = F) %>% 
  rename(ID = V1, name = V2) %>% 
  as_tibble()

# Generate DRUP metadata -------------------------------------------------------------------

drup.gut.metadata <- drup.gut.metadata %>% 
  left_join(clinical_data %>% select (-Cohort), by = "patientId") %>% 
  left_join(functional_dysbiosis) %>% 
  left_join(species_dysbiosis) %>% 
  mutate (
    Response = case_when(
      BOR_proxy %in% c("CR", "PR") ~ "R", 
      BOR_proxy %in% c("SD", "PD") ~ "NR", 
      .default = "NE"
    ),
    CB = as.factor (CB),
    tumor_class = ifelse (tumor_groups %in% c("CRC", "Other-GI"), 
                          tumor_groups, 
                          "Other-non-GI"),
    ATB_use_unknonw = ifelse (is.na(ATB_use), "Unknown", ATB_use),
    ATB_use = ifelse (ATB_use_unknonw == "Yes", "Yes", "No/Unknown")
  )

write_csv (drup.gut.metadata, "output/data/metagenomics/drup/drup-gut-metadata-metagenomics.csv")

# Filter samples that are not MSI or do not have clinical outcome data and do not have metagenomics data sequenced. 

samples_metagenomics =  list.files("data/metagenomics/DRUP/biobakery/results/kneaddata/") %>% 
  basename() 

drup.gut.metadata.fil <- drup.gut.metadata %>% 
  filter (msStatus != "MSS") %>% 
  filter (!is.na(BOR_proxy)) %>% 
  filter (BOR_proxy != "NE") %>% 
  filter (Sample %in% samples_metagenomics)

write_csv (drup.gut.metadata.fil, "output/data/metagenomics/drup/drup-gut-metadata-metagenomics-fil.csv")

# ------------------------------------------------------------------------ #
# 1. DRUP data --------------------------------------------------------------
# ------------------------------------------------------------------------ #

# Create folder if it does not exist. 
if (!file.exists("output/data/metagenomics/drup/")) {
  dir.create("output/data/metagenomics/drup/", recursive = TRUE)
}

## Metaphlan4 -------------------------------------------------------------
# ------------------------------------------------------------------------ #


files = list.files("/DATA/share/Voesties/data/resources/gut_microbiome/drup_gut/workflows/biobakery/results/metaphlan/read_stats/", full.names = T)
names(files) =  list.files("/DATA/share/Voesties/data/resources/gut_microbiome/drup_gut/workflows/biobakery/results/metaphlan/read_stats/") %>% 
  basename() %>% 
  str_remove("-reads.txt")

files = files[names(files) != "NA"] # Filter negative control sample. 

drup.gut.metadata_df <- drup.gut.metadata.fil %>% 
  filter (Sample %in% names(files)) %>% 
  mutate (Sample2 = Sample) %>% 
  as_data_frame() %>% 
  column_to_rownames(var = "Sample2")

pseq = make_phyloseq4(table_fn = files, sampledata = drup.gut.metadata_df) 

tse_reads = pseq %>% 
  convertFromPhyloseq() %>% 
  mia::transformAssay(assay.type = "counts",  method = "relabundance") %>% 
  mia::transformAssay( assay.type = "relabundance", method = "clr", pseudocount = T)

# Add C.innocuum 
rowData(tse_reads) <- rowtibble(tse_reads) %>% 
  mutate (Species = ifelse (SGB == "SGB4037", "Clostridium_innocuum", Species)) %>% 
  DataFrame()

write_rds(tse_reads, "output/data/metagenomics/drup/tse_reads_metaphlan.rds")

## Write fastspar input
# tse_reads %>%
#   agglomerateByRank(rank= "Species") %>%
#   subsetByPrevalent(prevalence = 0.05) %>%
#   assay("counts") %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "#OTU ID") %>%
#   write_tsv("../fastspar/drup_counts.tsv")

## Pathways --------------------------------------------------------------
# ------------------------------------------------------------------------ #

### Metacyc ###

# Import Metacyc abundances
metacyc.tse = import_metacyc("/DATA/share/Voesties/data/resources/gut_microbiome/drup_gut/workflows/biobakery/results/humann/pathway_abundance.txt", drup.gut.metadata_df, taxonomy = F) %>% 
  convertFromPhyloseq() %>% 
  transformAssay( method = "clr", pseudocount = T)

# Add pathway names to object. 
rowData (metacyc.tse) <- rowtibble(metacyc.tse) %>% 
  rename (ID = Species) %>% 
  left_join (metacyc.names) %>% 
  DataFrame()

# write_rds(metacyc.tse, "output/data/metagenomics/drup/tse_metacyc.rds")

# Import Metacyc abundances with taxa data. 
metacyc.tse.tax = import_metacyc("/DATA/share/Voesties/data/resources/gut_microbiome/drup_gut/workflows/biobakery/results/humann/pathway_abundance.txt",
                                 drup.gut.metadata_df,
                                 taxonomy = T) %>% 
  convertFromPhyloseq() %>% 
  transformAssay( method = "clr", pseudocount = T)

# Add pathway names to object. 
rowData (metacyc.tse.tax) <- rowtibble(metacyc.tse.tax) %>% 
  separate(Species, into = c("ID", "rank"), sep = "\\|") %>% 
  # rename (ID = Species) %>% 
  left_join (metacyc.names) %>% 
  DataFrame()

# write_rds(metacyc.tse.tax, "output/data/metagenomics/drup/tse_metacyc_tax.rds")


### Kegg Orthologs ###

orthologs.tse = import_metacyc("/DATA/share/Voesties/data/resources/gut_microbiome/drup_gut/workflows/biobakery/results/humann/regroup/kegg_families.txt",
                               drup.gut.metadata_df,
                               taxonomy = F,
                               type = "kegg") %>%
  convertFromPhyloseq() %>% 
  transformAssay( method = "clr", pseudocount = T)

# write_rds(orthologs.tse, "output/data/metagenomics/drup/tse_orthologs.rds")


### Eggnog families ###

eggnog.tse = import_metacyc("/DATA/share/Voesties/data/resources/gut_microbiome/drup_gut/workflows/biobakery/results/humann/regroup/eggnog_families.txt",
                               drup.gut.metadata_df,
                               taxonomy = F,
                               type = "eggnog") %>% 
  convertFromPhyloseq() %>% 
  transformAssay( method = "clr", pseudocount = T)

# write_rds(eggnog.tse, "output/data/metagenomics/drup/tse_eggnog_families.rds")


# ------------------------------------------------------------------------ #
# ICB-cohort data --------------------------------------------------------
# ------------------------------------------------------------------------ #

if (!file.exists("output/data/metagenomics/microbeICB/")) {
  dir.create("output/data/metagenomics/microbeICB/", recursive = TRUE)
}

# All data was already merged in the phyloseq format at:
# /DATA/share/Voesties/data/resources/gut_microbiome/microbeICB/integration/

# Create folder if it does not exist. 
if (!file.exists("output/data/metagenomics/microbeICB/")) {
  dir.create("output/data/metagenomics/microbeICB/", recursive = TRUE)
} 

## Metaphlan4 ------------------------------------------------------------------

### SGB level ###

tse_int.SGB <- read_rds("data/metagenomics/int.sgb.rds") %>% 
  convertFromPhyloseq()

# Add bachad adjusted data. 
int.SGB.comp.adj <- read_rds("data/metagenomics/int.sgb.comp.adj.rds") %>% 
  convertFromPhyloseq() %>% 
  assay("counts")

assay(tse_int.SGB, "relabundance_adj") <- int.SGB.comp.adj

# Add batch data:
colData(tse_int.SGB) <- coldf(tse_int.SGB) %>%
  mutate(batch = paste(Study, cohort, sep = ".") %>% str_replace(".NA", "")) %>%
  mutate(n_batch = n(), .by = batch) %>%
  DataFrame()

rowData(tse_int.SGB)$SGB <- rownames(tse_int.SGB)
  
# Save file.
write_rds(tse_int.SGB, "output/data/metagenomics/microbeICB/tse_int.SGB.rds")

### Save counts for fastspar: -----------------------------------------------------

tse_int.SGB <- read_rds("output/data/metagenomics/microbeICB/tse_int.SGB.rds")

# Meta-analaysis
metadata = coldf(tse_int.SGB) %>% 
  mutate(host_disease_status = as.factor(host_disease_status)) %>% 
  mutate(batch = paste(Study, cohort, sep = ".") %>% str_replace(".NA", "")) %>% 
  mutate(Study = ifelse(str_detect(Study, "Kim"), "Kim", Study)) %>% 
  mutate(n = n(), .by = Study)

tmp <- read_rds("output/data/metagenomics/microbeICB/tse_int.SGB.rds") %>%
  agglomerateByRanks("Species") %>%
  filter_se(Sample %in% metadata$Sample) %>% 
  subsetByPrevalent(prevalence = 0.05) # 773 species 1148 (samples)
rownames(tmp) <- rowData(tmp)$Species

tmp %>% assay("counts") %>%
  as.data.frame() %>%
  rownames_to_column(var = "#OTU ID") %>%
  write_tsv("../fastspar/microbe_icb_2_counts.tsv")

## Save SGB and Species names.

rowtibble(tmp) %>% 
  select(Species, SGB) %>% 
  write_csv("output/data/metagenomics/microbeICB/filter_Species_SGB_df.csv")



# -------------------------------------------------------------------------- #
# FG500  -------------------------------------------------------------------
# -------------------------------------------------------------------------- #

if (!file.exists("output/data/metagenomics/healthy/")) {
  dir.create("output/data/metagenomics/healthy/", recursive = TRUE)
}


# Import metadata
healthy.questions = read_csv("../studies/healthy_500fg/FG500_public_data_2022-11-18_13_12_11.csv") %>% 
  janitor::clean_names() %>% 
  rename(age = wat_is_uw_leeftijd_jaar,
         host_subject_id = x500fg_identifier) %>% 
  rename(sex = wat_is_uw_geslacht,
         contraception = gebruikt_u_de_pil_als_anticonceptiemiddel,
         smoking = bent_u_een_roker_of_een_roker_geweest)

# Create metadata table
healthy.metadata = read_csv("../studies/healthy_500fg/SraRunTable.csv") %>% 
  janitor::clean_names() %>% 
  mutate(sample_id = run) %>% 
  dplyr::left_join(healthy.questions) %>% 
  mutate(Study = "Healthy",
         Response = "Healthy",
         timepoint = "Healthy") %>% 
  column_to_rownames("sample_id")


## Metaphlan4 --------------------------------------------------------------

healthy.files = list.files("../studies/healthy_500fg/results/metaphlan/read_stats", full.names = T)
names(healthy.files) =  list.files("../studies/healthy_500fg/results/metaphlan/read_stats") %>% 
  basename() %>% 
  str_remove("-reads.txt")

healthy500fg.pseq = make_phyloseq4(table_fn = healthy.files, sampledata = healthy.metadata, tree = T)

tse_healthy500fg = healthy500fg.pseq %>% 
  convertFromPhyloseq() %>% 
  transformAssay( method = "relabundance") %>% 
  transformAssay( assay.type = "relabundance", method = "clr", pseudocount = T)

# Add C.innocuum 
rowData(tse_healthy500fg) <- rowtibble(tse_healthy500fg) %>% 
  mutate (Species = ifelse (SGB == "SGB4037", "Clostridium_innocuum", Species)) %>% 
  DataFrame()

write_rds (tse_healthy500fg, "output/data/metagenomics/healthy/tse_healthy500fg.rds") 

## Functional ----------------------------------------------------------------

### Eggnog orthologs ### 

healthy500fg.eggnog.tse = import_metacyc("../studies/healthy_500fg/results/humann/eggnog_families.txt",
                                         healthy.metadata,
                                         taxonomy = F,
                                         type = "eggnog") %>% 
  convertFromPhyloseq() %>% 
  transformAssay( method = "clr", pseudocount = T)

write_rds(healthy500fg.eggnog.tse, "output/data/metagenomics/healthy/tse_healthy500fg_eggnog_families.rds")


# ------------------------------------------------------------------------------
# DMP --------------------------------------------------------------------------
# ------------------------------------------------------------------------------


if (!file.exists("output/data/metagenomics/DMP/")) {
  dir.create("output/data/metagenomics/DMP/", recursive = TRUE)
}

# Load data. 

DAG3_merged <- read_csv("data/metagenomics/DMP/DAG3_merged_metaphlan4_Sep2024.csv")
DAG3_merged_phenotype <- read_csv("data/metagenomics/DMP/DAG3_basic_phenos.csv")

# Create phyloseq and TSE objects. 

DAG3_relabundance <- DAG3_merged %>%
  rename(Sample = `...1`) %>%
  select(c(Sample, contains("t__"))) %>%
  mutate(Sample = Sample %>% str_replace("_metaphlan", ""))

tax.table.info = tibble(
  otu_id = colnames(DAG3_relabundance)[-1],
  clade_name = colnames(DAG3_relabundance)[-1],
) %>%
  separate(
    otu_id,
    into = c(
      "Kingdom",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species",
      "SGB"
    ),
    sep = "\\|"
  ) %>%
  mutate(
    Kingdom = str_remove(Kingdom, "k__"),
    Phylum = str_remove(Phylum, "p__"),
    Class = str_remove(Class, "c__"),
    Order = str_remove(Order, "o__"),
    Family = str_remove(Family, "f__"),
    Genus = str_remove(Genus, "g__"),
    Species = str_remove(Species, "s__"),
    SGB = str_remove(SGB, "t__")
  ) %>%
  mutate(SGB = str_remove(SGB, "_group")) %>%
  mutate(
    Phylum = if_else(
      Phylum == "Eukaryota_unclassified",
      "Euryarchaeota",
      Phylum
    )
  ) %>%
  mutate(
    Phylum = if_else(Phylum == "Viruses_unclassified", "Viruses", Phylum)
  ) %>%
  mutate(SGB2 = SGB) %>%
  column_to_rownames(var = "SGB2") %>%
  as.matrix()

otu.table.stats = DAG3_relabundance %>% 
  column_to_rownames(var = "Sample") %>% 
  t() %>% 
  as.data.frame()

rownames(otu.table.stats) <- tax.table.info[,"SGB"] # Change rownames

sample.table <- DAG3_merged_phenotype %>% 
  filter(sampleID %in% colnames(otu.table.stats)) %>% 
  column_to_rownames(var = "sampleID")

otu.table= otu_table(otu.table.stats, taxa_are_rows = T)
taxa.table = tax_table(tax.table.info)
sampledata = sample_data(sample.table)

# Make phyloseq object (counts are relative abundance)
DAG3_pseq = phyloseq::phyloseq(otu.table, taxa.table, sampledata)

DAG3_tse = DAG3_pseq %>% 
  convertFromPhyloseq()

colData(DAG3_tse)$Sample <- colnames(DAG3_tse)

DAG3_tse <- DAG3_tse %>% 
  addPerCellQC() %>% 
  addPrevalence() %>% 
  mia::transformAssay( assay.type = "counts", method = "clr", pseudocount = T) 

DAG3_tse %>% 
  write_rds("output/data/metagenomics/DMP/DAG3_tse.rds")








