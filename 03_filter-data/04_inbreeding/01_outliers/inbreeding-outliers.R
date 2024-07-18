# inbreeding outliers
## write inbreeding outliers filter to remove outbred individuals from reference breeds
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(readr)) install.packages("readr")
if (!require(tidyr)) install.packages("tidyr")

### read inbreeding data, identify outliers and write filter to remove outbred individuals from reference breeds
readr::read_table("03_filter-data/04_inbreeding/01_outliers/ars-mind-ibs.het") %>%
  dplyr::group_by(FID) %>%
  dplyr::mutate(upper_outlier = dplyr::case_when(F > quantile(F)[4] + 1.5*IQR(F) ~ IID)) %>%
  dplyr::mutate(lower_outlier = dplyr::case_when(F < quantile(F)[2] - 1.5*IQR(F) ~ IID)) %>%
  tidyr::drop_na(lower_outlier) %>%
  dplyr::filter(FID %in% c("HOLS", "ANGU", "JERS",
                           "MUTU", "LAGU", "NDAG",
                           "THAR", "GIR", "NELO")) %>%
  dplyr::select(FID,
                IID) %>%
  readr::write_delim("03_filter-data/04_inbreeding/01_outliers/inbreeding-outliers.txt",
                     col_names = F)