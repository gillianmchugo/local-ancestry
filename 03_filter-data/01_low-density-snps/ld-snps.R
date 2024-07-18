# ld snps
## write list of ld snps
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(readr)) install.packages("readr")

### read schnabel file and write list of ld snps
readr::read_table("01_data-sources/02_downloaded/09_schnabel/UMC_marker_names_180910/9913_ARS1.2_58336_SNP50_marker_name_180910.map",
                  col_names = F) %>%
  dplyr::select(X2) %>%
  readr::write_delim("03_filter-data/01_low-density-snps/ld-snps.txt",
                     col_names = F)