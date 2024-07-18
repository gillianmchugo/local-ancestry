# ibs filter
## write ibs filter to remove duplicate individuals
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(readr)) install.packages("readr")
if (!require(stringr)) install.packages("stringr")

### read ibs results
ibs <- readr::read_table("03_filter-data/03_identity-by-state/ars-mind.mibs",
                         col_names = F) %>%
  as.data.frame()

### read ids
ids <- readr::read_table("03_filter-data/03_identity-by-state/ars-mind.mibs.id",
                         col_names = F) %>%
  dplyr::transmute(population = X1,
                   id = stringr::str_replace_all(X2,
                                                 "\\.",
                                                 "_"))

### set ids as column names
colnames(ibs) <- ids$id

### set ids as row names
rownames(ibs) <- ids$id

### identify triangular matrix
tri <- lower.tri(ibs)

### extract triangular matrix values
ibs_tri <- data.frame(id1 = colnames(ibs)[col(ibs)[tri]],
                      id2 = rownames(ibs)[row(ibs)[tri]],
                      value = ibs[tri])

### rename id columns
ids <- ids %>%
  dplyr::rename(population1 = population,
                id1 = id)

### add population names for id1
ibs_tri <- ibs_tri %>%
  dplyr::inner_join(ids)

### rename id columns
ids <- ids %>%
  dplyr::rename(population2 = population1,
                id2 = id1)

### add population names for id2
ibs_tri <- ibs_tri %>%
  dplyr::inner_join(ids)

### write ibs filter to remove duplicate individuals
ibs_tri %>%
  dplyr::filter(value >= 0.99) %>%
  dplyr::mutate(population = dplyr::case_when(stringr::str_detect(id1,
                                                                  "barbato") ~ population1,
                                              stringr::str_detect(id2,
                                                                  "barbato") ~ population2,
                                              .default = population1)) %>%
  dplyr::mutate(id = dplyr::case_when(stringr::str_detect(id1,
                                                          "barbato") ~ id1,
                                      stringr::str_detect(id2,
                                                          "barbato") ~ id2,
                                      .default = id1)) %>%
  dplyr::select(population,
                id) %>%
  dplyr::distinct() %>%
  readr::write_delim("03_filter-data/03_identity-by-state/ibs-filter.txt",
                     col_names = F)