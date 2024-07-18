# filter data
## filter snp data using plink and r
### extract low density snps
Rscript 03_filter-data/01_low-density-snps/ld-snps.R
plink --cow --bfile 02_update-genome-assembly/ars --extract 03_filter-data/01_low-density-snps/ld-snps.txt --make-bed --out 03_filter-data/01_low-density-snps/ars-ld

### missing genotype filter to remove individuals with missing ld snps
plink --cow --bfile 03_filter-data/01_low-density-snps/ars-ld --mind 0.95 --make-bed --out 03_filter-data/02_missing-snps/ars-ld-mind
plink --cow --bfile 02_update-genome-assembly/ars --keep 03_filter-data/02_missing-snps/ars-ld-mind.fam --make-bed --out 03_filter-data/02_missing-snps/ars-mind

### identity by state filter to remove duplicate individuals
plink --cow --bfile 03_filter-data/02_missing-snps/ars-mind --distance ibs square --out 03_filter-data/03_identity-by-state/ars-mind
Rscript 03_filter-data/03_identity-by-state/ibs-plot.R
Rscript 03_filter-data/03_identity-by-state/ibs-filter.R
plink --cow --bfile 03_filter-data/02_missing-snps/ars-mind --remove 03_filter-data/03_identity-by-state/ibs-filter.txt --make-bed --out 03_filter-data/03_identity-by-state/ars-mind-ibs
plink --cow --bfile 03_filter-data/02_missing-snps/ars-ld-mind --remove 03_filter-data/03_identity-by-state/ibs-filter.txt --make-bed --out 03_filter-data/03_identity-by-state/ars-ld-mind-ibs

### inbreeding outliers filter to remove outbred individuals from reference breeds
plink --cow --bfile 03_filter-data/03_identity-by-state/ars-mind-ibs --het --out 03_filter-data/04_inbreeding/01_outliers/ars-mind-ibs
Rscript 03_filter-data/04_inbreeding/01_outliers/inbreeding-outliers.R
Rscript 03_filter-data/04_inbreeding/01_outliers/inbreeding-outliers-plot.R
plink --cow --bfile 03_filter-data/03_identity-by-state/ars-mind-ibs --remove 03_filter-data/04_inbreeding/01_outliers/inbreeding-outliers.txt --make-bed --out 03_filter-data/04_inbreeding/01_outliers/ars-mind-ibs-het
plink --cow --bfile 03_filter-data/03_identity-by-state/ars-ld-mind-ibs --remove 03_filter-data/04_inbreeding/01_outliers/inbreeding-outliers.txt --make-bed --out 03_filter-data/04_inbreeding/01_outliers/ars-ld-mind-ibs-het

### inbreeding filter to remove excess individuals from european reference breeds
plink --cow --bfile 03_filter-data/04_inbreeding/01_outliers/ars-mind-ibs-het --het --out 03_filter-data/04_inbreeding/02_filter/ars-mind-ibs-het
Rscript 03_filter-data/04_inbreeding/02_filter/inbreeding-filter.R
plink --cow --bfile 03_filter-data/04_inbreeding/01_outliers/ars-mind-ibs-het --remove 03_filter-data/04_inbreeding/02_filter/inbreeding-filter.txt --make-bed --out 03_filter-data/04_inbreeding/02_filter/ars-mind-ibs-het-filter
plink --cow --bfile 03_filter-data/04_inbreeding/01_outliers/ars-ld-mind-ibs-het --remove 03_filter-data/04_inbreeding/02_filter/inbreeding-filter.txt --make-bed --out 03_filter-data/04_inbreeding/02_filter/ars-ld-mind-ibs-het-filter
plink --cow --bfile 03_filter-data/04_inbreeding/02_filter/ars-ld-mind-ibs-het-filter --het --out 03_filter-data/04_inbreeding/02_filter/ars-ld-mind-ibs-het-filter
Rscript 03_filter-data/04_inbreeding/02_filter/inbreeding-filter-plot.R

### filter snps by autosome, call rate and minor allele frequency
plink --cow --bfile 03_filter-data/04_inbreeding/02_filter/ars-mind-ibs-het-filter --chr 1-29 --geno 0.05 --maf 0.05 --make-bed --out 03_filter-data/05_filter-snps/hd
plink --cow --bfile 03_filter-data/04_inbreeding/02_filter/ars-ld-mind-ibs-het-filter --chr 1-29 --geno 0.05 --maf 0.05 --make-bed --out 03_filter-data/05_filter-snps/ld
