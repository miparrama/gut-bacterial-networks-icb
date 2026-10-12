
# Create conda environment with fastspar
# conda create -n fastspar -c conda-forge -c bioconda fastspar=1.0.0
conda activate fastspar

cohort="microbe_icb_2"
# cohort="$1"

# Create folders
mkdir -p output/covariance
mkdir -p output/correlation
mkdir -p output/pvalues

fastspar --otu_table ${cohort}_counts.tsv \
         --correlation output/correlation/${cohort}_correlation.tsv \
         --covariance output/covariance/${cohort}_covariance.tsv \
         --iterations 50 \
         --threads 10

echo "Fastspar analysis completed."
echo "Starting bootstrap analysis..."

# Generate bootstrap matrices and correlations 
mkdir -p bootstrap/bootstrap_counts
fastspar_bootstrap --otu_table ${cohort}_counts.tsv --number 1000 --prefix bootstrap/bootstrap_counts/${cohort}_fake_data

mkdir -p bootstrap/bootstrap_correlation
mkdir -p bootstrap/bootstrap_covariance
parallel -j 4 fastspar \
    --otu_table {} \
    --correlation bootstrap/bootstrap_correlation/cor_{/} \
    --covariance bootstrap/bootstrap_covariance/cov_{/} -i 5 ::: bootstrap/bootstrap_counts/*

# Compute p-values based on bootstrap results
fastspar_pvalues --otu_table ${cohort}_counts.tsv \
    --correlation output/correlation/${cohort}_correlation.tsv \
    --prefix bootstrap/bootstrap_correlation/cor_${cohort}_fake_data_ \
    --permutations 1000 \
    --outfile output/pvalues/${cohort}_pvalues.tsv

# fastspar --otu_table ${cohort}_counts.txt --correlation median_correlation.tsv --covariance median_covariance.tsv --iterations 50 --threads 10


# Delete bootstrap files to save space
rm -r bootstrap
