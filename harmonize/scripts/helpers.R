# Libraries and functions used in the R project.  
# @ Miguel Parra Martinez

here::i_am("scripts/helpers.R")

# Libraries -------------------------------------------------------------------

library(mia)
library(phyloseq)
library(microbiome)

library(patchwork)
library(ggsci)
library(ggthemes)
library(ggpubr)
library(ggbeeswarm)

library(janitor)
library(readxl)
library(tidyverse)

library(conflicted)
library(here)

# Define conflicts preferences.
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::left_join)
conflicts_prefer(dplyr::inner_join)
conflicts_prefer(dplyr::count)
conflicts_prefer(base::intersect)
conflicts_prefer(microbiome::transform)
conflicts_prefer(janitor::fisher.test)
conflicts_prefer(edgeR::cpm)


# Import resources ------------------------------------------------------------

# Import metacyc Names
metacyc.names = vroom::vroom(
  here("data/resources/map_metacyc-pwy_name.txt.gz"),
  col_names = c("ID", "Name")
)
# Import KEGG modules
module_names = read_csv(
  here("data/resources/module_names.csv")
) %>%
  mutate(short_name = str_remove(Name, " =>.*$"))


# Functions -------------------------------------------------------------------

make_phyloseq4 = function(
  table_fn,
  sampledata,
  tree = T,
  relab = F,
  chocophlan_nwk = here("data/resources/mpa_vOct22_CHOCOPhlAnSGB_202212.nwk")
) {
  # Import metaphlan results (rel w readstats)
  if (relab) {

    metaphlan = table_fn %>%
      map_dfr(.id = "sampleId", .f = function(x) {
        suppressMessages(data.table::fread(x, skip = 4))
      }) %>%
      rename(clade_name = `#clade_name`, reads = relative_abundance) %>%
      filter(grepl("t__", clade_name)) %>%
      select(-additional_species, -NCBI_tax_id) %>%
      separate(
        clade_name,
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
      mutate(SGB = str_remove(SGB, "_group"))

    # Get OTU info (use mp2 read stats)
    otu.table.stats = metaphlan %>%
      pivot_wider(
        names_from = sampleId,
        values_from = reads,
        values_fill = 0
      ) %>%
      select(-Kingdom:-Species) %>%
      column_to_rownames("SGB") %>%
      otu_table(taxa_are_rows = T)
  } else {
    
    metaphlan = table_fn %>%
      map_dfr(.id = "sampleId", .f = function(x) {
        suppressMessages(data.table::fread(x, skip = 5))
      }) %>%
      rename(
        clade_name = `#clade_name`,
        reads = estimated_number_of_reads_from_the_clade
      ) %>%
      filter(grepl("t__", clade_name)) %>%
      select(-coverage, -clade_taxid) %>%
      separate(
        clade_name,
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
      mutate(SGB = str_remove(SGB, "_group"))

    # Get OTU info (use mp2 read stats)
    otu.table.stats = metaphlan %>%
      select(-relative_abundance) %>%
      pivot_wider(
        names_from = sampleId,
        values_from = reads,
        values_fill = 0
      ) %>%
      select(-Kingdom:-Species) %>%
      column_to_rownames("SGB") %>%
      otu_table(taxa_are_rows = T)
  }

  # Get taxa info
  taxa.table = metaphlan %>%
    distinct(SGB, .keep_all = T) %>%
    select(Kingdom:SGB) %>%
    mutate(sp = SGB) %>%
    column_to_rownames("sp") %>%
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
    as.matrix() %>%
    tax_table()

  # Get sample data
  sampledata = sampledata %>%
    rownames_to_column("SAMPLE") %>%
    filter(SAMPLE %in% colnames(otu.table.stats)) %>%
    column_to_rownames("SAMPLE") %>%
    sample_data()

  # Make phyloseq object (w/ estimated read counts)
  pseq.stats = phyloseq(otu.table.stats, taxa.table, sampledata)

  # Get phylogenetic tree
  if (tree) {
    mpa_tree = ape::read.tree(chocophlan_nwk)
    mpa_tree$tip.label = paste0("SGB", mpa_tree$tip.label)
    pseq.tree = merge_phyloseq(pseq.stats, mpa_tree)
    return(pseq.tree)
  } else {
    pseq.stats
  }
}

import_metacyc = function(table, sampledata, taxonomy = T, type = "metacyc", exclude = T) {
  pathways = data.table::fread(table)
  colnames(pathways)[1] = "metacyc"
  
  # Remove taxa
  if(!taxonomy){
    pathways = pathways %>% 
      filter(!(grepl("\\|", metacyc))) 
  }
  
  if(taxonomy){
    pathways = pathways %>% 
      separate(metacyc, into = c("metacyc", "tax"), sep = "\\|")
  }
  
  # Separate names
  pathways = pathways %>% 
    separate(metacyc, into = c("ID", "Name"), sep = ":")
  
  if(exclude) {
    pathways = pathways %>% 
      filter(!(ID %in% c("UNMAPPED", "UNINTEGRATED", "UNGROUPED")))
  }
  
  if(taxonomy){
    pathways = pathways %>% 
      mutate (ID = paste(ID, tax, sep = "|"))
  }
  
  # Fix colnames
  if(type == "metacyc") {
    colnames(pathways) = str_remove(colnames(pathways), "_Abundance")
  } else {
    colnames(pathways) = str_remove(colnames(pathways), "_Abundance-RPKs")
  }
  
  # Fix samples not found
  #sampledata.fil = sampledata %>% 
  #  filter(run %in% colnames(pathways))
  
  # Convert to phyloseq object
  pathways.otu = pathways %>% 
    select(-any_of(c("Name","NA","tax"))) %>% 
    column_to_rownames("ID") %>% 
    otu_table(taxa_are_rows = T) 
  
  pathways.tax = pathways %>% 
    select(-any_of(c("Name","NA"))) %>%  
    column_to_rownames("ID") %>% 
    rownames_to_column("Species") %>% 
    select(Species) %>% 
    mutate(id = Species) %>% 
    column_to_rownames("id") %>% 
    as.matrix() %>% 
    tax_table()
  
  # Return pathways as phyloseq object
  pseq = phyloseq(pathways.otu, pathways.tax, sample_data(sampledata))
  return(pseq)
  
}


# Functions to quickly access the metadata in tse object
rowtibble <- function (tse){
  row_data_tibble <- rowData(tse) %>% as_tibble()
  return(row_data_tibble)
}

coltibble <- function (tse){
  col_data_tibble <- colData(tse) %>% as_tibble()
  return(col_data_tibble)
}

rowdf <- function (tse){
  row_data_df <- rowData(tse) %>% as.data.frame()
  return(row_data_df)
}

coldf <- function (tse){
  col_data_df <- colData(tse) %>% as.data.frame()
  return(col_data_df)
}

