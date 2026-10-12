#-----------------------------------------------------------------------------#
# Project: Bacterial clusters
# Script:  functional_per_bacteria.R
# Purpose: Script to check the functional differences between bacterial clusters of interest
# Author:  Miguel Parra-Martinez
# Date:    2026-03-10
#-----------------------------------------------------------------------------#

# Load variables form 04_bacterial_communities

# Metacyc ---------------------------------------------------------------------------

#Import metacyc pathways names
metacyc.names = data.table::fread("data/resources/map_metacyc-pwy_name.txt.gz", header = F) %>% 
  rename(ID = V1, name = V2) %>% 
  as_tibble()

## Pathways per bacteria ------------------------------------------------------------
metacyc.tse.tax <- read_rds("data/metagenomics/int.metacyc.tax.rda") %>% 
  mia::convertFromPhyloseq() %>%
  transformAssay(method = "relabundance") %>%
  transformAssay(method = "pa") %>% 
  transformAssay(
    method = "clr",
    assay.type = "relabundance",
    pseudocount = T
  )

# Add pathway names to object. 
rowData (metacyc.tse.tax) <- rowdf(metacyc.tse.tax) %>% 
  separate(Species, into = c("ID", "rank"), sep = "\\|") %>% 
  # rename (ID = Species) %>% 
  left_join (metacyc.names) %>% 
  DataFrame()

# Add species column 
rowData(metacyc.tse.tax)$Species = str_extract(rowData(metacyc.tse.tax)$rank, "(?<=s__)[^ ]+")

tax_info_Oct22 <- tax_info_Oct22 %>%
  mutate(
    Species = str_extract(SGBs, "s__[^|]+") %>%
      str_remove("^s__"),
    Species_ncbi = str_split(NCBI_tax, "\\|") %>%
      map_chr(~ .x[7]),
    SGB = str_remove_all(SGBs, "_group") %>% str_extract("(SGB|EUK)\\d+$")
  )

species_cluster_1 <- tax_info_Oct22 %>% 
  filter(SGB %in% sgb_cluster_1) %>% 
  pull(Species) %>% 
  unique()

species_cluster_2 <- tax_info_Oct22 %>% 
  filter(SGB %in% sgb_cluster_2) %>% 
  pull(Species) %>% 
  unique()

metacyc.tax.long <- metacyc.tse.tax %>%
  filter_se_rows(Species %in% c(species_cluster_1, species_cluster_2)) %>%
  meltSE(assay.type = "relabundance", add.row = T) %>%
  mutate(
    cluster = case_when(
      Species %in% species_cluster_1 ~ "1",
      Species %in% species_cluster_2 ~ "2",
      .default = "Other"
    ),
    pa = ifelse(relabundance > 0, 1, 0)
  )

## Difference in abundance between clusters ------------------------------------------------------------ 

ratio_relabundance <- metacyc.tax.long %>%
  group_by(name, cluster, SampleID) %>%
  summarize(relabundance = sum(relabundance)) %>%
  ungroup() %>%
  arrange(SampleID, -relabundance) %>%
  pivot_wider(names_from = cluster, values_from = relabundance, values_fill = 0) %>%
  mutate(ratio = `2` - `1`) %>% 
  arrange(ratio) %>% 
  mutate(
    mean_1 = mean(`1`),
    mean_2 = mean(`2`),
    .by = name
  ) %>% 
  mutate(
    diff = mean_1 - mean_2
  )

relative_diff_pathway <- ratio_relabundance %>%
  distinct(name, mean_1, mean_2, diff) %>%
  mutate(
    percent = abs(diff) / max(mean_1, mean_2),
    .by = name
  ) %>%
  mutate(
    sign_percent = sign(diff) * percent
  )


# Save plot based on difference between pathways of the microbiome. 

pathways_bacteria_contribution_barplot <- ratio_relabundance %>%
  distinct(name, mean_1, mean_2, diff) %>%
  arrange(-diff) %>%
  {
    bind_rows(slice_head(., n = 20), slice_tail(., n = 20))
  } %>%
  left_join(
    fit_lm_metacyc$meta_fits %>%
      mutate(ID = feature) %>%
      left_join(metacyc.names)
  ) %>%
  mutate(sign_pval = -log10(pval) * sign(coef)) %>%
  mutate(class = ifelse(diff > 0, "high in cluster 1", "high in cluster 2")) %>%
  left_join(
    relative_diff_pathway %>% 
      select(name, percent, sign_percent)
  ) %>% 
  filter(!is.na(pval)) %>%
  mutate(name = fct_reorder(str_to_sentence(name), sign_pval)) %>%
  ggplot(aes(x = sign_pval, y = name)) +
  geom_col(aes(fill = sign_percent)) +
  scale_fill_gradient2(
    high = muted("blue"),
    mid = "white",
    low = muted("red"),
    midpoint = 0,
    limits = c(-1, 1),
    guide = guide_colorbar(barwidth = unit(6.5, "cm"))
  ) +
  geom_vline(xintercept = c(-log10(0.05), log10(0.05)), linetype = "dashed") +
  labs(
    y = NULL,
    x = "Response association\n(signed log10 p-value)",
    fill = "Relative abundance difference\n(Cluster 1 - Cluster 2 bacteria)",
    subtitle = "Pathways annotated to Cluster 1 or Cluster 2 bacteria"
  ) +
  theme_bw(base_size = 14) 

pathways_bacteria_contribution_barplot +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom", legend.box = "vertical")

# Correlation plot all pathways
scatter_plot <- fit_lm_community_score %>%
  mutate(sign_pval_community_score = -log10(pval) * sign(coef)) %>%
  select(name, sign_pval_community_score) %>%
  left_join(
    fit_lm_metacyc$meta_fits %>%
      mutate(ID = feature) %>%
      left_join(metacyc.names) %>%
      mutate(sign_pval_response = -log10(pval) * sign(coef)) %>%
      select(name, sign_pval_response)
  ) %>%
  left_join(
    ratio_relabundance %>%
      distinct(name, mean_1, mean_2, diff) %>%
      arrange(-diff)
  ) %>%
  left_join(
    relative_diff_pathway %>%
      select(name, percent, sign_percent)
  ) %>%
  left_join(
    tse_int_metacyc_fil %>%
      transformAssay(method = "relabundance") %>%
      addPerFeatureQC(assay.type = "relabundance") %>%
      rowtibble() %>%
      select(name, mean_relabundance = mean)
  ) %>%
  mutate(name = str_to_sentence(name) %>% str_wrap(width = 30)) %>%
  ggplot(aes(sign_pval_community_score, sign_pval_response)) +
  geom_smooth(method = "lm", color = "black") +
  geom_vline(xintercept = c(-log10(0.05), log10(0.05)), linetype = "dashed") +
  geom_hline(yintercept = c(-log10(0.05), log10(0.05)), linetype = "dashed") +
  geom_label_repel(
    aes(label = name),
    label.size = 0.2,
    size = 3,
    label.padding = 0.15, # shrinks the padding inside each label box (default 0.25 "lines")
    # max.overlaps = Inf, # don't silently drop labels that still overlap
    max.overlaps = 12,
    force = 1, # how strongly labels repel each other (increase if still crowded)
    min.segment.length = 0.2, # shorter threshold before a connector line is drawn
    seed = 43
  ) +
  geom_point(
    data = . %>% filter(is.na(sign_percent)),
    aes(color = sign_percent, size = mean_relabundance)
    # size = 5
  ) +
  # geom_point(
  #   data = . %>% filter(!is.na(sign_percent)),
  #   aes(color = sign_percent, size = mean_relabundance)
  #   # size = 5
  # ) +
  scale_color_gradient2(
    high = muted("blue"),
    mid = "white",
    low = muted("red"),
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  labs(
    x = "Community score association\n(signed log10 p-value, meta-analysis)",
    y = "Response association\n(signed log10 p-value, meta-analysis)",
    color = "Relative abundance difference\n(Cluster 1 - Cluster 2 bacteria)",
    size = "Mean pathway abundance"
  ) +
  # stat_cor(size = 4) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "Pearsons's r = 0.73, p < 0.001, n = 463 ",
    hjust = -0.2,
    vjust = 4,
    size = 4.5
  ) +
  theme_few(base_size = 14) +
  guides(
    color = guide_colorbar(
      order = 2,
      barwidth = 15,
      barheight = 0.8,
      direction = "horizontal",
      title.position = "top"
    ),
    size = guide_legend(
      order = 1,
      direction = "horizontal",
      title.position = "top",
      nrow = 1
    )
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.02),
    legend.justification = c("right", "bottom"),
    legend.box = "vertical",
    legend.box.background = element_rect(
      color = NA,
      fill = alpha("white", 0.8)
    ),
    legend.box.margin = margin(4, 4, 4, 4),
    legend.background = element_blank() # explicitly turn off the per-legend background/box
  )

scatter_plot

ggsave(
  "output/figures/figures_manuscript/metacyc_scatter_all.pdf",
  scatter_plot,
  width = 15,
  height = 10
)

## Save table ------------------------------------------------------------

fit_lm_community_score_table <- fit_lm_community_score %>% 
  mutate(exposure = "Community_score_low") %>% 
  arrange(pval)

fit_lm_response_table <- fit_lm_metacyc$meta_fits %>%
  left_join(metacyc.names, by = c("feature" = "ID")) %>%
  arrange(pval)

relative_diff_pathway_table <- relative_diff_pathway %>% 
  rename(
    mean_abundance_cluster_1 = mean_1,
    mean_abundance_cluster_2 = mean_2,
    difference = diff,
    relative_difference = percent,
    signed_relative_difference = sign_percent
  ) %>%
   left_join(metacyc.names) %>%
  arrange(-relative_difference)

openxlsx::write.xlsx(
  list(
    "Community_score_association" = fit_lm_community_score_table,
    "Response_association" = fit_lm_response_table,
    "Relative_difference_Cluster1_2" = relative_diff_pathway_table
  ),
  file = "output/data/tables_manuscript/metacyc_pathways.xlsx",
  asTable = TRUE
)




