# inbreeding filter
## write inbreeding filter to remove excess individuals from european reference breeds
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(readr)) install.packages("readr")

### read inbreeding data and filter european reference breeds
het <- readr::read_table("03_filter-data/04_inbreeding/02_filter/ars-mind-ibs-het.het") %>%
  dplyr::group_by(FID) %>%
  dplyr::filter(FID %in% c("ANGU", "HOLS", "JERS"))

### identify upper 25 individuals to remove
max <- het %>%
  dplyr::slice_max(F,
                   prop = 0.199) %>%
  dplyr::select(FID,
                IID) 

### identify lower 25 individuals to remove
min <- het %>%
  dplyr::slice_min(F,
                   prop = 0.19) %>%
  dplyr::select(FID,
                IID) 

### combine upper and lower individuals and write filter
dplyr::bind_rows(max,
                 min) %>%
  readr::write_delim("03_filter-data/04_inbreeding/02_filter/inbreeding-filter.txt",
                     col_names = F)