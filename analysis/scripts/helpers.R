# Libraries and functions used in the R project.  
# @ Miguel Parra Martinez

# Libraries ---------------------------------------------------------------

library(janitor)
library(broom.mixed)

library(readxl)

library(ggpubr)
library(ggthemes)
library(ggbeeswarm)
library(ggside)
library(ggnewscale)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggrepel)
library(patchwork)
library(miaViz)
library(scater)
library(ggstats)
library(ggalluvial)
library(ggsci)
library(scales)
library(Cairo)
library(ggvenn)
library(ggnewscale)

library(bluster)

library(survival) 
library(survminer)

library(Maaslin2)
library(MMUPHin)

library(phyloseq)
library(mia)
library(tidyverse)


# Colors -----------------------------------------------------------------

colors_studies <- c(
  "Derosa" = "#FC8D62",
  "DRUP" =  "#E5C494",
  "Gunjur" = "#A6D854", 
  "Kim" = "#FFD92F",
  "Lee" = "#66C2A5",
  "Routy" = "#8DA0CB",
  "Spencer" = "#E78AC3"
)


# Functions ---------------------------------------------------------------

# Map colors to names of brewer paletter. 
map_colors <- function(variable, palette_name = "Set3", remove_na = FALSE) {
  # Check if the input is a factor; if not, convert it to a factor
  if (!is.factor(variable)) {
    variable <- factor(variable)
  }
  
  # Optionally drop NA before computing levels/colors
  if (remove_na) {
    variable <- factor(variable[!is.na(variable)])
  }
  
  # Get the unique levels of the factor
  levels <- levels(variable)
  num_levels <- length(levels)
  
  # Check if the palette has enough colors
  if (num_levels > brewer.pal.info[palette_name, "maxcolors"]) {
    message(sprintf("The selected palette (%s) does not have enough colors for %s levels.",
                     palette_name, 
                     num_levels)
    )
    message("Extending the number of colors in the palette....")
    
    # Generate the color palette with the same number of colors as the number of levels
    color_palette <- colorRampPalette(brewer.pal(8, palette_name))(num_levels)
    
  } else {
    
    # Generate the color palette with the same number of colors as the number of levels
    color_palette <- brewer.pal(num_levels, palette_name)
  }
  
  # Create a named vector that maps each level to a color
  color_mapping <- setNames(color_palette, levels)
  
  return(color_mapping)
}


# Updated survminer::ggforest fixing a bug.
ggforest2 <- function(
  model,
  data = NULL,
  main = "Hazard ratio",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 0.7,
  refLabel = "reference",
  noDigits = 2
) {
  conf.high <- conf.low <- estimate <- NULL
  stopifnot(inherits(model, "coxph"))
  data <- survminer:::.get_data(model, data = data)
  terms <- attr(model$terms, "dataClasses")[-1]
  coef <- as.data.frame(broom::tidy(model, conf.int = TRUE))
  gmodel <- broom::glance(model)
  allTerms <- lapply(seq_along(terms), function(i) {
    var <- names(terms)[i]
    if (terms[i] %in% c("factor", "character")) {
      adf <- as.data.frame(table(data[, var]))
      cbind(var = var, adf, pos = 1:nrow(adf))
    } else if (terms[i] == "numeric") {
      data.frame(var = var, Var1 = "", Freq = nrow(data), pos = 1)
    } else {
      vars <- grep(paste0("^", var, "*."), coef$term, value = TRUE)
      data.frame(
        var = vars,
        Var1 = "",
        Freq = nrow(data),
        pos = seq_along(vars)
      )
    }
  })
  allTermsDF <- do.call(rbind, allTerms)
  colnames(allTermsDF) <- c("var", "level", "N", "pos")
  inds <- apply(allTermsDF[, 1:2], 1, paste0, collapse = "")
  rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
  toShow <- cbind(allTermsDF, coef[inds, ])[, c(
    "var",
    "level",
    "N",
    "p.value",
    "estimate",
    "conf.low",
    "conf.high",
    "pos"
  )]
  toShowExp <- toShow[, 5:7]
  toShowExp[is.na(toShowExp)] <- 0
  toShowExp <- format(exp(toShowExp), digits = noDigits)
  toShowExpClean <- data.frame(
    toShow,
    pvalue = signif(
      toShow[,
        4
      ],
      noDigits + 1
    ),
    toShowExp
  )
  toShowExpClean$stars <- paste0(
    round(toShowExpClean$p.value, noDigits + 1),
    " ",
    ifelse(toShowExpClean$p.value < 0.05, "*", ""),
    ifelse(toShowExpClean$p.value < 0.01, "*", ""),
    ifelse(toShowExpClean$p.value < 0.001, "*", "")
  )
  toShowExpClean$ci <- paste0(
    "(",
    toShowExpClean[, "conf.low.1"],
    " - ",
    toShowExpClean[, "conf.high.1"],
    ")"
  )
  toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] <- refLabel
  toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] <- "<0.001 ***"
  toShowExpClean$stars[is.na(toShowExpClean$estimate)] <- ""
  toShowExpClean$ci[is.na(toShowExpClean$estimate)] <- ""
  toShowExpClean$estimate[is.na(toShowExpClean$estimate)] <- 0
  toShowExpClean$var <- as.character(toShowExpClean$var)
  toShowExpClean$var[duplicated(toShowExpClean$var)] <- ""
  toShowExpClean$N <- paste0("(N=", toShowExpClean$N, ")")
  toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1, ]
  ci_h <- toShowExpClean$conf.high[is.finite(as.numeric(
    toShowExpClean$conf.high.1
  ))]
  ci_l <- toShowExpClean$conf.low[
    !near(as.numeric(toShowExpClean$conf.low.1), 0)
  ]
  toShowExpClean$conf.high[is.infinite(as.numeric(
    toShowExpClean$conf.high.1
  ))] <- max(ci_h, na.rm = TRUE)
  toShowExpClean$conf.low[near(
    as.numeric(toShowExpClean$conf.low.1),
    0
  )] <- min(ci_l, na.rm = TRUE)
  rangeb <- range(
    toShowExpClean$conf.low,
    toShowExpClean$conf.high,
    na.rm = TRUE
  )
  breaks <- axisTicks(rangeb / 2, log = TRUE, nint = 7)
  rangeplot <- rangeb
  rangeplot[1] <- rangeplot[1] - diff(rangeb)
  rangeplot[2] <- rangeplot[2] + 0.15 * diff(rangeb)
  width <- diff(rangeplot)
  y_variable <- rangeplot[1] + cpositions[1] * width
  y_nlevel <- rangeplot[1] + cpositions[2] * width
  y_cistring <- rangeplot[1] + cpositions[3] * width
  y_stars <- rangeb[2]
  x_annotate <- seq_len(nrow(toShowExpClean))
  annot_size_mm <- fontsize *
    as.numeric(grid::convertX(unit(theme_get()$text$size, "pt"), "mm"))
  p <- ggplot(toShowExpClean, aes(seq_along(var), exp(estimate))) +
    geom_rect(aes(
      xmin = seq_along(var) - 0.5,
      xmax = seq_along(var) +
        0.5,
      ymin = exp(rangeplot[1]),
      ymax = exp(rangeplot[2]),
      fill = ordered(seq_along(var) %% 2 + 1)
    )) +
    scale_fill_manual(values = c("#FFFFFF33", "#00000033"), guide = "none") +
    geom_point(pch = 15, size = 4) +
    geom_errorbar(
      aes(ymin = exp(conf.low), ymax = exp(conf.high)),
      width = 0.15
    ) +
    geom_hline(yintercept = 1, linetype = 3) +
    coord_flip(ylim = exp(rangeplot)) +
    ggtitle(main) +
    scale_y_log10(
      name = "",
      labels = sprintf("%g", breaks),
      expand = c(0.02, 0.02),
      breaks = breaks
    ) +
    theme_light() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "none",
      panel.border = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    xlab("") +
    annotate(
      geom = "text",
      x = x_annotate,
      y = exp(y_variable),
      label = toShowExpClean$var,
      fontface = "bold",
      hjust = 0,
      size = annot_size_mm
    ) +
    annotate(
      geom = "text",
      x = x_annotate,
      y = exp(y_nlevel),
      hjust = 0,
      label = toShowExpClean$level,
      vjust = -0.1,
      size = annot_size_mm
    ) +
    annotate(
      geom = "text",
      x = x_annotate,
      y = exp(y_nlevel),
      label = toShowExpClean$N,
      fontface = "italic",
      hjust = 0,
      vjust = ifelse(toShowExpClean$level == "", 0.5, 1.1),
      size = annot_size_mm
    ) +
    annotate(
      geom = "text",
      x = x_annotate,
      y = exp(y_cistring),
      label = toShowExpClean$estimate.1,
      size = annot_size_mm,
      vjust = ifelse(toShowExpClean$estimate.1 == "reference", 0.5, -0.1)
    ) +
    annotate(
      geom = "text",
      x = x_annotate,
      y = exp(y_cistring),
      label = toShowExpClean$ci,
      size = annot_size_mm,
      vjust = 1.1,
      fontface = "italic"
    ) +
    annotate(
      geom = "text",
      x = x_annotate,
      y = exp(y_stars),
      label = toShowExpClean$stars,
      size = annot_size_mm,
      hjust = -0.2,
      fontface = "italic"
    ) +
    annotate(
      geom = "text",
      x = 0.5,
      y = exp(y_variable),
      label = paste0(
        "# Events: ",
        gmodel$nevent,
        "; Global p-value (Log-Rank): ",
        format.pval(gmodel$p.value.log, eps = ".001"),
        " \nAIC: ",
        round(gmodel$AIC, 2),
        "; Concordance Index: ",
        round(gmodel$concordance, 2)
      ),
      size = annot_size_mm,
      hjust = 0,
      vjust = 1.2,
      fontface = "italic"
    )
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  ggpubr::as_ggplot(gt)
}

# Functions to quickly access the rowData and colData in form of a tse
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

min_nonzero <- function(x, na.rm = TRUE) {
  # Remove zeros (and optionally NAs)
  x_nz <- x[x != 0]
  if (na.rm) x_nz <- x_nz[!is.na(x_nz)]
  
  if (length(x_nz) == 0) {
    return(NA)  # or Inf, or stop() if no non-zero values exist
  } else {
    return(min(x_nz))
  }
}

abbrev_species <- function(x) {
  sapply(x, function(name) {
    # Case 1: GGB..._SGB...
    if (str_detect(name, "_SGB")) {
      return(str_extract(name, "SGB\\d+"))
    }
    
    # Case 2: Genus_species
    parts <- str_split_fixed(name, "_", 2)
    genus <- parts[,1]
    rest  <- parts[,2]
    
    genus_abbr <- paste0(substr(genus, 1, 1), ".")
    
    if (rest == "" | is.na(rest)) {
      abbrev = (genus_abbr)
    } else {
      abbrev = (paste0(genus_abbr, rest))
    }

    abbrev = abbrev %>% 
      str_replace_all("_", " ")
    return(abbrev)

  })
}

# Function to filter a SummarizedExperiment by sample metadata
filter_se <- function(se, ...) {
  meta <- as.data.frame(colData(se))
  
  keep <- meta %>% filter(...) %>% rownames()
  
  se[, keep]
}

# Function to filter a SummarizedExperiment by sample metadata
filter_se_rows <- function(se, ...) {
  meta <- as.data.frame(rowData(se))
  
  keep <- meta %>% filter(...) %>% rownames()
  
  se[keep,]
}

# Import humann3 pathway abundances
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
    phyloseq::otu_table(taxa_are_rows = T) 
  
  pathways.tax = pathways %>% 
    select(-any_of(c("Name","NA"))) %>%  
    column_to_rownames("ID") %>% 
    rownames_to_column("Species") %>% 
    select(Species) %>% 
    mutate(id = Species) %>% 
    column_to_rownames("id") %>% 
    as.matrix() %>% 
    phyloseq::tax_table()
  
  # Return pathways as phyloseq object
  pseq = phyloseq(pathways.otu, pathways.tax, sample_data(sampledata))
  return(pseq)
  
}

make_phyloseq4 = function(table_fn, sampledata, tree = T, tree_file = "data/resources/mpa_vOct22_CHOCOPhlAnSGB_202212.nwk", relab = F) {
  
  # Import metaphlan results (rel w readstats)
  if(relab){
    metaphlan = table_fn %>% 
      map_dfr(.id = "sampleId",
              .f = function(x) {suppressMessages(data.table::fread(x, skip = 4))},
              .progress = T) %>% 
      rename(clade_name = `#clade_name`,
             reads = relative_abundance) %>% 
      filter(grepl("t__", clade_name)) %>% # Select only SGB
      select(-additional_species, -NCBI_tax_id)  %>% 
      separate(clade_name, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SGB"), sep = "\\|") %>% 
      mutate(Kingdom = str_remove(Kingdom, "k__"),
             Phylum = str_remove(Phylum, "p__"),
             Class = str_remove(Class, "c__"),
             Order = str_remove(Order, "o__"),
             Family = str_remove(Family, "f__"),
             Genus = str_remove(Genus, "g__"), 
             Species = str_remove(Species, "s__"),
             SGB = str_remove(SGB, "t__")) %>% 
      mutate(SGB = str_remove(SGB, "_group"))
    
    # Get OTU info (use mp2 read stats)
    otu.table.stats = metaphlan %>% 
      pivot_wider(names_from = sampleId, values_from = reads, values_fill = 0) %>% 
      select(-Kingdom:-Species) %>% 
      column_to_rownames("SGB") %>% 
      otu_table(taxa_are_rows = T)
    
  } else{
    metaphlan = table_fn %>% 
      map_dfr(.id = "sampleId",
              .f = function(x) {suppressMessages(data.table::fread(x, skip = 5))},
              .progress = T) %>% 
      rename(clade_name = `#clade_name`,
             reads = estimated_number_of_reads_from_the_clade) %>% 
      filter(grepl("t__", clade_name)) %>% # Select only SGB
      select(-coverage, -clade_taxid)  %>% 
      separate(clade_name, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SGB"), sep = "\\|") %>% 
      mutate(Kingdom = str_remove(Kingdom, "k__"),
             Phylum = str_remove(Phylum, "p__"),
             Class = str_remove(Class, "c__"),
             Order = str_remove(Order, "o__"),
             Family = str_remove(Family, "f__"),
             Genus = str_remove(Genus, "g__"), 
             Species = str_remove(Species, "s__"),
             SGB = str_remove(SGB, "t__")) %>% 
      mutate(SGB = str_remove(SGB, "_group"))
    
    # Get OTU info (use mp2 read stats)
    otu.table.stats = metaphlan %>% 
      select(-relative_abundance) %>% 
      pivot_wider(names_from = sampleId, values_from = reads, values_fill = 0) %>% 
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
    mutate(Phylum = if_else(Phylum == "Eukaryota_unclassified", "Euryarchaeota", Phylum)) %>%
    mutate(Phylum = if_else(Phylum == "Viruses_unclassified", "Viruses", Phylum)) %>% 
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
  if(tree){
    mpa_tree = ape::read.tree(tree_file)
    mpa_tree$tip.label = paste0("SGB", mpa_tree$tip.label)
    pseq.tree = merge_phyloseq(pseq.stats, mpa_tree)
    return(pseq.tree)
  } else{
    pseq.stats
  }
}

correlation_fisher <- function (tse, 
                                prevalence = 0.1,
                                qval_cutoff = 0.1,
                                p_adjust_method = "bonferroni",
                                fraction_sig = 0.05){
  # Function to compute the co-occurrence correlation of bacteria  based on -log10(pval) of a fisher test. 
  
  tse <- tse %>% 
    transformAssay(method = 'pa') %>% 
    subsetByPrevalent(prevalence = prevalence) 
  
  species_fil <- rownames (tse)

  pairwise_combinations <- combn(species_fil, 2, simplify = FALSE)

  presence_mt <- assay(tse, "pa")
  
  interaction_scores_df <- map_dfr(pairwise_combinations, function(pair) {
    mgs1 <- factor(presence_mt[pair[1], ], levels = c("0", "1"))
    mgs2 <- factor(presence_mt[pair[2], ], levels = c("0", "1"))
    table_data <- table(mgs1, mgs2)
    test <- fisher.test(table_data)
    score <- -log10(test$p.value) * sign(test$estimate - 1)
    
    # Return a tibble (data frame) with the desired columns
    tibble(pairs = paste(pair, collapse = " vs "),
           score = score,
           p.value = test$p.value,
           species1 = pair[1],
           species2 = pair[2])
  }, .progress = T)
  
  
  interaction_scores_df$qval <- p.adjust(interaction_scores_df$p.value,
                                         method = p_adjust_method )
  
  df <- tibble(
    pval = interaction_scores_df$p.value, 
    qval = interaction_scores_df$qval
  )
  
  pval_limit <- df %>% 
    arrange(pval) %>% 
    filter (qval < qval_cutoff) %>% 
    tail(1) %>% 
    pull (pval)
  
  
  # Adjacent matrix
  
  sig_species_interactions <- interaction_scores_df %>% 
    pivot_longer(cols = c("species1", "species2"),
                 names_to = "names", values_to = "species") %>% 
    distinct(species) %>% 
    pull(species)
  
  interaction_scores_df_fil <- interaction_scores_df %>% 
    filter ( species2 %in% sig_species_interactions & 
               species1 %in% sig_species_interactions)
  
  
  # Construction of adjacent matrix 
  right_diag <- interaction_scores_df_fil %>% 
    # mutate (score = ifelse (p.value > pval_limit, 0, score)) %>%
    dplyr::select (species1, species2, score ) 
  left_diag <- interaction_scores_df_fil %>%
    # mutate (score = ifelse (p.value > pval_limit, 0, score)) %>%
    mutate(species1 = right_diag$species2,
           species2 = right_diag$species1) %>% 
    dplyr::select (species1, species2, score ) 
  diag <- tibble(
    species1 = sig_species_interactions,
    species2 = sig_species_interactions,
    score = 0
  )
  adj_matrix <- rbind (right_diag, diag, left_diag) %>% 
    pivot_wider(names_from = species2, values_from = score) %>% 
    column_to_rownames(var = "species1") %>%
    as.matrix()
  adj_matrix <- adj_matrix[,sort (colnames (adj_matrix))]
  adj_matrix <- adj_matrix[sort (rownames (adj_matrix)),]
  
  # Construction of pval matrix
  right_diag <- interaction_scores_df_fil %>% 
    # mutate (score = ifelse (p.value > pval_limit, 0, score)) %>%
    dplyr::select (species1, species2, p.value ) 
  left_diag <- interaction_scores_df_fil %>%
    # mutate (score = ifelse (p.value > pval_limit, 0, score)) %>%
    mutate(species1 = right_diag$species2,
           species2 = right_diag$species1) %>% 
    dplyr::select (species1, species2, p.value ) 
  diag <- tibble(
    species1 = sig_species_interactions,
    species2 = sig_species_interactions,
    p.value = 1
  )
  pval_matrix <- rbind (right_diag, diag, left_diag) %>% 
    pivot_wider(names_from = species2, values_from = p.value) %>% 
    column_to_rownames(var = "species1") %>%
    as.matrix()
  
  pval_matrix <- pval_matrix[,sort (colnames (pval_matrix))]
  pval_matrix <- pval_matrix[sort (rownames (pval_matrix)),]
  
  
  # Sparsification of the matrix
  # adj_matrix[pval_matrix > pval_limit] <- 0
  
  # Filter for species with at least X% of significant species. 
  # fraction_significant <- length (species_fil) * (fraction_sig)
  # adj_matrix_fil <- adj_matrix[colSums(adj_matrix != 0) > fraction_significant, colSums(adj_matrix != 0) > fraction_significant]
  
  return (list(adj_matrix, pval_matrix))
}



roc_from_taxa <- function(
  tse,
  group,
  assay_name = "counts",
  response_col = "Response",
  positive_class = "Responder",
  direction_roc = "<",
  batch_col = NULL,
  base_size = 14, 
  make_plot = TRUE
) {
  require(pROC)
  require(ggplot2)
  require(dplyr)

  mat <- assay(tse, assay_name)

  if (!all(group %in% rownames(mat))) {
    message(length(missing), " taxa in 'group' not found in assay rownames")
    # Keep only the variables present
    group <- intersect(group, rownames(mat))
  }

  # Summed abundance for the group
  group_abund <- colSums(mat[group, , drop = FALSE])

  # Metadata
  meta <- as.data.frame(colData(tse))
  meta$abund <- group_abund

  if (!is.null(batch_col)) {
    # --- ROC per batch ---
    # Loop over batches
    roc_dfs <- list()
    auc_table <- data.frame(
      batch = character(),
      auc = numeric(),
      stringsAsFactors = FALSE
    )

    for (b in unique(meta[[batch_col]])) {
      dat <- meta %>% filter(.data[[batch_col]] == b)

      y <- dat[[response_col]]
      if (!is.factor(y)) {
        y <- factor(y)
      }
      y <- relevel(y, ref = setdiff(levels(y), positive_class)[1])
      y_bin <- as.numeric(y == positive_class)

      roc_obj <- roc(
        response = y_bin,
        predictor = dat$abund,
        direction = direction_roc,
        quiet = TRUE
      )

      auc_val <- auc(roc_obj)

      roc_df <- data.frame(
        fpr = 1 - roc_obj$specificities,
        tpr = roc_obj$sensitivities,
        batch = b
      )

      roc_dfs[[b]] <- roc_df
      auc_table <- rbind(
        auc_table,
        data.frame(batch = b, auc = as.numeric(auc_val))
      )
    }

    # Combine all ROC curves into one dataframe
    roc_df_all <- bind_rows(roc_dfs)

    # Plot
    p <- roc_df_all %>%
      arrange(fpr, tpr) %>%
      mutate(batch = fct_relevel(batch, auc_table$batch)) %>%
      ggplot(aes(x = fpr, y = tpr, color = batch)) +
      geom_line(size = 1) +
      geom_abline(linetype = "dashed", color = "gray50") +
      coord_equal() +
      theme_bw(base_size = base_size) +
      labs(
        title = ifelse(
          length(group) < 6,
          paste("ROC curves for", paste(group, collapse = "+")),
          "ROC curve"
        ),
        x = "False Positive Rate",
        y = "True Positive Rate",
        color = "Batch"
      ) +
      scale_color_brewer(
        palette = "Set1",
        labels = paste0(auc_table$batch, " (AUC=", round(auc_table$auc, 3), ")")
      ) +
      scale_x_continuous(expand = c(0.01, 0), limits = c(0, 1)) +
      scale_y_continuous(expand = c(0.01, 0), limits = c(0, 1)) +
      theme(
        legend.position = c(0.95, 0.05), # x = 95% right, y = 5% from bottom
        legend.justification = c("right", "bottom"),
        legend.background = element_rect(
          fill = alpha("white", 0.6),
          color = NA
        ),
        panel.grid = element_blank()
      )
  } else {
    # --- Single ROC across all samples ---
    y <- meta[[response_col]]
    if (!is.factor(y)) {
      y <- factor(y)
    }
    y <- relevel(y, ref = setdiff(levels(y), positive_class)[1])
    y_bin <- as.numeric(y == positive_class)

    roc_obj <- roc(
      response = y_bin,
      predictor = meta$abund,
      direction = direction_roc,
      quiet = TRUE
    )
    auc_val <- auc(roc_obj)

    roc_df <- data.frame(
      fpr = 1 - roc_obj$specificities,
      tpr = roc_obj$sensitivities
    ) %>%
      arrange(fpr, tpr)

    auc_table <- data.frame(batch = "All samples", auc = as.numeric(auc_val))

    # Plot
    p <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
      geom_line(color = "steelblue", size = 1) +
      geom_abline(linetype = "dashed", color = "gray50") +
      coord_equal() +
      theme_bw(base_size = base_size) +
      labs(
        title = ifelse(
          length(group) < 6,
          paste(
            "ROC curve for",
            paste(group, collapse = "+"),
            "(AUC =",
            round(auc_val, 3),
            ")"
          ),
          paste(
            "ROC curve ",
            "(AUC =",
            round(auc_val, 3),
            ")"
          )
        ),
        x = "False Positive Rate",
        y = "True Positive Rate"
      ) +
      scale_x_continuous(expand = c(0.01, 0), limits = c(0, 1)) +
      scale_y_continuous(expand = c(0.01, 0), limits = c(0, 1)) +
      theme(panel.grid = element_blank())
  }
  if(make_plot == TRUE){
    return(list(auc_table = auc_table, plot = p))
  } else {
    return(list(auc_table = auc_table))
  }
  
}

roc_from_col <- function(
  tse,
  metadata_variable,
  assay_name = "counts",
  response_col = "Response",
  positive_class = "Responder",
  direction_roc = "<",
  batch_col = NULL,
  base_size = 14,
  average = FALSE,
  make_plot = TRUE
) {
  require(pROC)
  require(ggplot2)
  require(dplyr)

  # mat <- assay(tse, assay_name)

  # Check if the variable is present.
  if (!metadata_variable %in% colnames(colData(tse))) {
    stop("'", metadata_variable, "' is not a column in colData(tse)")
  }

  group <- colData(tse)[[metadata_variable]]

  # Metadata
  meta <- as.data.frame(colData(tse))

  meta$abund <- group

  if (!is.null(batch_col)) {
    
    batches <- unique(meta[,batch_col])

    # --- ROC per batch ---
    # Loop over batches
    fpr_grid <- seq(0, 1, length.out = 101)

    # --- Updated Loop ---
    roc_list_interpolate <- list()
    roc_dfs <- list()
    auc_table <- data.frame(batch = character(), auc = numeric())

    for (b in unique(meta[[batch_col]])) {
      dat <- meta %>% filter(.data[[batch_col]] == b)

      # Standard ROC calculation
      roc_obj <- roc(
        response = dat[[response_col]],
        predictor = dat$abund,
        direction = direction_roc,
        quiet = TRUE
      )

      # Interpolate TPR onto our fixed FPR grid
      # pROC specificities are usually 1-FPR, so we handle that:
      interp_tpr <- approx(
        x = 1 - roc_obj$specificities,
        y = roc_obj$sensitivities,
        xout = fpr_grid,
        ties = max
      )$y

      roc_list_interpolate[[b]] <- data.frame(
        fpr = fpr_grid,
        tpr = interp_tpr,
        batch = b
      )

      roc_df <- data.frame(
        fpr = 1 - roc_obj$specificities,
        tpr = roc_obj$sensitivities,
        batch = b
      )

      roc_dfs[[b]] <- roc_df

      auc_table <- rbind(
        auc_table,
        data.frame(batch = b, auc = as.numeric(auc(roc_obj)))
      )
    }

    # 2. Combine individual curves
    roc_interpolate_all <- bind_rows(roc_list_interpolate)

    # 3. Calculate the Average (Summary) Curve
    summary_roc <- roc_interpolate_all %>%
      group_by(fpr) %>%
      summarise(
        mean_tpr = mean(tpr, na.rm = TRUE),
        sd_tpr = sd(tpr, na.rm = TRUE),
        se_tpr = sd_tpr / sqrt(n()),
        # lower = max(0, mean_tpr - se_tpr),
        # upper = min(1, mean_tpr + se_tpr)
        lower = max(0, mean_tpr - sd_tpr),
        upper = min(1, mean_tpr + sd_tpr)
      ) %>%
      filter(fpr != 0) %>%
      add_row(
        fpr = 0,
        mean_tpr = 0,
        sd_tpr = 0,
        se_tpr = 0,
        lower = 0,
        upper = 0
      ) %>% 
      filter(fpr != 1) %>%
      add_row(
        fpr = 1,
        mean_tpr = 1,
        sd_tpr = 0,
        se_tpr = 0,
        lower = 0,
        upper = 0
      )
    
    roc_df_all <- bind_rows(roc_dfs)

    # Plot
    if (average == T) {
      p <- roc_df_all %>%
        arrange(fpr, tpr) %>%
        mutate(batch = fct_relevel(batch, auc_table$batch)) %>%
        ggplot(aes(x = fpr, y = tpr, color = batch)) +
        geom_line(size = 1, alpha = 0.75) +
        geom_ribbon(
          data = summary_roc,
          aes(x = fpr, ymin = lower, ymax = upper),
          inherit.aes = F,
          fill = "grey50",
          alpha = 0.3
        ) +
        geom_line(
          data = summary_roc,
          aes(x = fpr, y = mean_tpr),
          color = "black",
          size = 1.5
        ) +
        geom_abline(linetype = "dashed", color = "gray50") +
        coord_equal() +
        theme_bw(base_size = base_size) +
        labs(
          title = metadata_variable,
          x = "False Positive Rate",
          y = "True Positive Rate",
          color = "Batch"
        ) +
        scale_color_manual(
          values = c(map_colors(batches, "Set2")),
          labels = paste0(
            auc_table$batch,
            " (AUC=",
            round(auc_table$auc, 3),
            ")"
          )
        ) +
        scale_x_continuous(expand = c(0.01, 0), limits = c(0, 1)) +
        scale_y_continuous(expand = c(0.01, 0), limits = c(0, 1)) +
        theme(
          legend.position = c(0.98, 0.02),
          legend.justification = c("right", "bottom"),
          legend.background = element_rect(
            fill = alpha("white", 0.6),
            color = NA
          ),
          panel.grid = element_blank()
        )
    } else {
        p <- roc_df_all %>%
        arrange(fpr, tpr) %>%
        mutate(batch = fct_relevel(batch, auc_table$batch)) %>%
        ggplot(aes(x = fpr, y = tpr, color = batch)) +
        geom_line(size = 1) +
        geom_abline(linetype = "dashed", color = "gray50") +
        coord_equal() +
        theme_bw(base_size = base_size) +
        labs(
          title = metadata_variable,
          x = "False Positive Rate",
          y = "True Positive Rate",
          color = "Batch"
        ) +
        scale_color_manual(
          values = c(map_colors(batches, "Set2")),
          labels = paste0(
            auc_table$batch,
            " (AUC=",
            round(auc_table$auc, 3),
            ")"
          )
        ) +
        scale_x_continuous(expand = c(0.01, 0), limits = c(0, 1)) +
        scale_y_continuous(expand = c(0.01, 0), limits = c(0, 1)) +
        theme(
          legend.position = c(0.98, 0.02),
          legend.justification = c("right", "bottom"),
          legend.background = element_rect(
            fill = alpha("white", 0.6),
            color = NA
          ),
          panel.grid = element_blank()
        )
    }

  } else {
    # --- Single ROC across all samples ---
    y <- meta[[response_col]]
    if (!is.factor(y)) {
      y <- factor(y)
    }
    y <- relevel(y, ref = setdiff(levels(y), positive_class)[1])
    y_bin <- as.numeric(y == positive_class)

    roc_obj <- roc(
      response = y_bin,
      predictor = meta$abund,
      direction = direction_roc,
      quiet = TRUE
    )
    auc_val <- auc(roc_obj)

    roc_df <- data.frame(
      fpr = 1 - roc_obj$specificities,
      tpr = roc_obj$sensitivities
    ) %>%
      arrange(fpr, tpr)

    auc_table <- data.frame(batch = "All samples", auc = as.numeric(auc_val))

    # Plot
    p <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
      geom_line(color = "steelblue", size = 1) +
      geom_abline(linetype = "dashed", color = "gray50") +
      coord_equal() +
      theme_bw(base_size = base_size) +
      labs(
        title = metadata_variable,
        x = "False Positive Rate",
        y = "True Positive Rate"
      ) +
      scale_x_continuous(expand = c(0.01, 0), limits = c(0, 1)) +
      scale_y_continuous(expand = c(0.01, 0), limits = c(0, 1)) +
      theme(panel.grid = element_blank())
  }
  if (make_plot == TRUE) {
    return(list(auc_table = auc_table, plot = p))
  } else {
    return(list(auc_table = auc_table))
  }
}


get_abundance_set <- function(tse, bacteria_set, pseudocount = 0.00001) {
  # Function to get the meadian log abundance of a bacterial group. 
  tse <- tse %>% 
    transformAssay(method = "relabundance") %>% 
    transformAssay(method = "log", assay.type = "relabundance", pseudocount = pseudocount)
  
  mat <- assay(tse, "relabundance")
  
  bacteria_present <- bacteria_set[bacteria_set %in% rownames(mat)]
  if(length(bacteria_present) == 0) {
    message(length(bacteria_set) - length(bacteria_set %in% rownames(mat)), " taxa in 'group' not found in assay rownames")
  }
  
  # Subset to selected bacteria
  mat_subset <- mat[bacteria_present, , drop = FALSE]
  
  total_abun <- apply(mat_subset, 2, sum)
  
  mat_subset_counts <- assay(tse, "counts")[bacteria_present, , drop = FALSE]
  
  fraction_present <- apply(mat_subset_counts, 2, function(x) {
    sum(x > 0, na.rm = TRUE) / length(bacteria_present)
  })
  
  return(list(
    "total_abun" = total_abun,
    "prevalence" = fraction_present))
}

plot_survival_quantile <- function(tse, time_col, event_col, var_col, q = 0.5, ...) {
  # tse       : TSE object
  # time_col  : column in colData(tse) with survival time
  # event_col : column in colData(tse) with survival event (0/1)
  # var_col   : continuous variable to stratify samples (colData)
  # q         : quantile cutoff (0 < q < 0.5)
  library(survival)
  library(survminer)
  
  # Extract metadata
  meta <- as.data.frame(colData(tse))
  
  if(!all(c(time_col, event_col, var_col) %in% colnames(meta))) {
    stop("time_col, event_col, or var_col not found in colData(tse)")
  }
  
  # Get values of the variable
  vals <- meta[[var_col]]
  
  # Compute cutoffs
  lower_cut <- quantile(vals, probs = q, na.rm = TRUE)
  upper_cut <- quantile(vals, probs = 1 - q, na.rm = TRUE)
  
  # Assign group: bottom, top, or NA (middle samples)
  meta$group <- NA
  meta$group[vals <= lower_cut] <- paste0("Bottom_", q*100, "%")
  meta$group[vals >= upper_cut] <- paste0("Top_", q*100, "%")
  
  # Keep only defined groups
  df <- meta %>% filter(!is.na(group))
  
  df$event = as.numeric(df[[event_col]])
  df$time = as.numeric(df[[time_col]])
  
  fit <- survfit(Surv(time, event) ~ group, data = df)
  fit
  
  df %>%
    coxph(Surv(time, event) ~ group, data = .) %>%
    broom::tidy()
  
  g = df %>% 
    ggsurvplot(
      fit = survfit(Surv(time, event) ~ group, data = .),
      pval = TRUE,
      risk.table = TRUE,
      ggtheme = theme_bw(), 
      ...
    )
  
  return(g)
}




prevalRatioClusters <- function(
  tse,
  otu_pos,
  otu_neg
) {

  preval_neg= get_median_abundance_set(
    tse,
    otu_neg
  )$prevalence
  preval_pos = get_median_abundance_set(
    tse,
    otu_pos
  )$prevalence

  preval_score = ((preval_pos - preval_neg) + 1) / 2
  return(preval_score)
}


create_forest_plot <- function(
  data,
  name_col = "Name",
  title_annotation_left = "Metabolite protective",
  title_annotation_right = "Metabolite harmful",
  title = "Cox regression results"
) {
  # Prepare the data
  res_plot <- data %>%
    mutate(
      # Format the HR estimates with padding
      est_fmt = str_pad(
        round(estimate, 2),
        width = 4,
        pad = "0",
        side = "right"
      ),
      ci_low_fmt = str_pad(
        round(conf.low, 2),
        width = 4,
        pad = "0",
        side = "right"
      ),
      ci_high_fmt = str_pad(
        round(conf.high, 2),
        width = 4,
        pad = "0",
        side = "right"
      ),
      estimate_lab = paste0(est_fmt, " (", ci_low_fmt, "-", ci_high_fmt, ")"),

      # Format p-values
      p_fmt = case_when(
        p.value < 0.001 ~ "<0.001",
        TRUE ~ str_pad(
          as.character(round(p.value, 3)),
          width = 4,
          pad = "0",
          side = "right"
        )
      ),

      # Rename for plotting
      model = .data[[name_col]]
    ) %>%
    # Add header row
    bind_rows(
      data.frame(
        model = "Metabolite",
        estimate_lab = "Hazard Ratio (95% CI)",
        p_fmt = "p-value",
        estimate = NA,
        conf.low = NA,
        conf.high = NA
      )
    ) %>%
    mutate(model = fct_rev(fct_relevel(model, "Metabolite")))

  # Calculate x-axis limits for the middle plot
  x_min <- min(log(data$conf.low), na.rm = TRUE) - 0.2
  x_max <- max(log(data$conf.high), na.rm = TRUE) + 0.2

  # Middle plot - forest plot on log scale
  p_mid <- res_plot %>%
    filter(model != "Metabolite") %>%
    ggplot(aes(y = model)) +
    theme_classic() +
    geom_point(aes(x = log(estimate)), shape = 15, size = 3) +
    geom_linerange(aes(xmin = log(conf.low), xmax = log(conf.high))) +
    labs(x = "Log Hazard Ratio") +
    coord_cartesian(
      ylim = c(1, nrow(data) + 1),
      xlim = c(x_min, x_max)
    ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    annotate(
      "text",
      x = x_min + 0.1,
      # y = nrow(data) + 1,
       y = 0.5,
      label = title_annotation_left,
      hjust = 0,
      size = 3
    ) +
    annotate(
      "text",
      x = x_max - 0.1,
      # y = nrow(data) + 1,
       y = 0.5,
      label = title_annotation_right,
      hjust = 1,
      size = 3
    ) +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank()
    )

  # Left side - metabolite names only
  p_left1 <- res_plot %>%
    ggplot(aes(y = model)) +
    geom_text(
      aes(x = 0, label = model),
      hjust = 0,
      fontface = ifelse(res_plot$model == "Metabolite", "bold", "plain")
    ) +
    theme_void() +
    coord_cartesian(xlim = c(0, 1))

  # Left side - HR estimates
  p_left2 <- res_plot %>%
    ggplot(aes(y = model)) +
    geom_text(
      aes(x = 0, label = estimate_lab),
      hjust = 0,
      fontface = ifelse(
        res_plot$estimate_lab == "Hazard Ratio (95% CI)",
        "bold",
        "plain"
      )
    ) +
    theme_void() +
    coord_cartesian(xlim = c(0, 1))

  # Right side - p-values
  p_right <- res_plot %>%
    ggplot() +
    geom_text(
      aes(x = 0, y = model, label = p_fmt),
      hjust = 0,
      fontface = ifelse(res_plot$p_fmt == "p-value", "bold", "plain")
    ) +
    theme_void()

  # Layout design with 4 panels
  layout <- c(
    area(t = 0, l = 0, b = 30, r = 4), # metabolite names
    area(t = 0, l = 5, b = 30, r = 8), # HR estimates
    area(t = 1, l = 9, b = 30, r = 14), # forest plot
    area(t = 0, l = 15, b = 30, r = 17) # p-values
  )

  # Combine plots
  final_plot <- p_left1 +
    p_left2 +
    p_mid +
    p_right +
    plot_layout(design = layout) +
    plot_annotation(
      title = title
    )

  # # Combine plots
  # final_plot <- p_left + p_mid + p_right + plot_layout(design = layout)
  
  return(final_plot)
}
