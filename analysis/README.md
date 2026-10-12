5--- output: html_document ---

# Microbiome data analysis in patients treated with ICB.

This Rproject analyzes the metabolomic and metagenomic data of the MMRd/MSI DRUP cohort treated with immunotherapy.

## Data origin:

-   data/metagenomics/functional_dysbiosis —\> /home/m.p.martinez/projects/gut-microbiome/drup_gut/analysis/output (cancer dyscovery results)

-   data/metabolites/**nkimet_addtional-compounds.xlsx** —\> Uploaded from laptop (obtained from email from Iris?)

-   data/metabolites/**metabolites_hmdb_metaboanalyst.csv**—\> Data obtained from metaboanalist website to annotate metabolites of class 2b with HMDB id. (2024-11-28).

-   data/metagenomics/DRUP/drup-gut-metadata.csv —\> /DATA/share/Voesties/data/resources/gut_microbiome/drup_gut/preprocessing/drup-gut-metadata.csv. This file is generated in generate-metadata.R script (found in Voesties/data/resources/gut-microbiome) and puts together msi-cohort metadata and gut-metagenomics cohort data (questionaire) and timepoints.

-   data/metagenomics/DRUP/clinical_data_elaborate_20250314.csv -\> Clinical metadata of DRUP.

-   data/metagenomics/int.sgb.comp.adj.rds & int.sgb.rds are comming from harmonize/output/data/merge

-   data/metagenomics/int.metacyc.rda -\> harmonize/output/data/merge/int.metacyc.rda

-   data/metagenomics/int.eggnog.rda -\> harmonize/output/data/merge/int.eggnog.rda

-   data/metagenomics/int.metacyc.tax.rda -\> harmonize/output/data/merge/int.metacyc.tax.rda

#