## Co-occurring gut bacterial networks associate with clinical outcome to immune checkpoint blockade

This repository contains the code pertinent to our manuscript. Please refer to the Materials and Methods section of the article for details.
>Miguel Parra-Martinez, Iris L. Mimpen, Catherine Toner-Bartelds, Krijn K. Dijkstra, Lodewyk F.A. Wessels, Emile E. Voest,
Co-occurring gut bacterial networks associate with clinical outcome to immune checkpoint blockade

## Abstract
The gut microbiome is a key regulator of immune function and is linked to immune checkpoint blockade (ICB) outcomes, yet inconsistent signatures across cohorts and limited mechanistic understanding hinder therapeutic translation. Most studies examine individual species, underappreciating that bacteria live in communities that collectively shape metabolic output and host immunity. Using metagenomics data from seven ICB-treated cohorts (n=1,122), unsupervised co-abundance network analysis revealed two anti-correlated clusters (Cluster 1 and Cluster 2) enriched in bacteria associated with favorable and poor response to ICB, respectively. Cluster 1 was composed of common gut commensals, whereas Cluster 2 was enriched for dysbiosis-associated bacteria. A community score based on the two clusters predicted responses to ICB across studies more consistently than single. Microbiomes dominated by either of the two clusters of bacteria showed profound differences in functional pathways and metabolites detected by metagenomics and metabolomics. Collectively, these findings show that no single bacterial species, but two anti-correlated bacterial communities with distinct functions, are associated with ICB response. This provides a foundation for therapies aimed at rebalancing gut microbial communities and their functional output to improve ICB response.

## Analysis Overview and files. 
The analyses performed include:

- `analysis/`: Analysis code of the metagenomic and metabolomic datasets.This includes the R code to generate the figures of the manuscript. 
- `harmonize/`: Code to integrate the different data sources of metagenomic data. 
- `fastspar/`: Co-abundance analysis of co-ocurring bacteria in ICB patients.
- `studies`: Folder were the different studies included in the analysis were processed. 

## Pipeline
This analysis makes use of a Google Cloud (GCP) based Nextflow pipeline to process whole metagenome sequences with the Biobakery tools Humann and Metaphlan.
The pipeline can be found at: https://github.com/miparrama/biobakery-nf

## Contributors

- [Miguel Parra-Martinez](https://github.com/miparrama) - [@miparrama](https://github.com/miparrama)