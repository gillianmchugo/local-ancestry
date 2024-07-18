# local ancestry results
## extract results of local ancestry analyses
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(MOSAIC)) devtools::install_github("mststats/MOSAIC")
if (!require(parallel)) install.packages("parallel")
if (!require(readr)) install.packages("readr")
if (!require(scales)) install.packages("scales")
if (!require(stringr)) install.packages("stringr")
if (!require(tibble)) install.packages("tibble")

### function to read mosaic data and extract mean ancestry results
mosaic_data <- function(population, order = c(1, 2, 3), data = "hd"){
  # list rdata files containing model parameters and local ancestry results
  files <- list.files(paste0("07_local-ancestry-analysis/01_mosaic/02_output/",
                             dplyr::case_when(data == "hd" ~ "01_hd/",
                                              data == "ld" ~ "02_ld/"),
                             dplyr::case_when(population == "ALEN" ~ "01_ALEN",
                                              population == "ANKO" ~ "02_ANKO",
                                              population == "BORA" ~ "03_BORA",
                                              population == "BORG" ~ "04_BORG",
                                              population == "CHIA" ~ "05_CHIA",
                                              population == "EASZ" ~ "06_EASZ",
                                              population == "KARA" ~ "07_KARA",
                                              population == "KETE" ~ "08_KETE",
                                              population == "MARC" ~ "09_MARC",
                                              population == "MARE" ~ "10_MARE",
                                              population == "NDAM" ~ "11_NDAM",
                                              population == "NGAN" ~ "12_NGAN",
                                              population == "ROMA" ~ "13_ROMA",
                                              population == "SHEK" ~ "14_SHEK",
                                              population == "SOMB" ~ "15_SOMB")),
                      full.names = T,
                      pattern = "*.RData")
  # load file 1
  load(files[1])
  # load file 2
  load(files[2])
  # get local ancestry results for snp positions
  localanc_positions <- MOSAIC::grid_to_pos(localanc,
                                            paste0("07_local-ancestry-analysis/01_mosaic/02_output/",
                                                   dplyr::case_when(data == "hd" ~ "01_hd/",
                                                                    data == "ld" ~ "02_ld/"),
                                                   dplyr::case_when(population == "ALEN" ~ "01_ALEN",
                                                                    population == "ANKO" ~ "02_ANKO",
                                                                    population == "BORA" ~ "03_BORA",
                                                                    population == "BORG" ~ "04_BORG",
                                                                    population == "CHIA" ~ "05_CHIA",
                                                                    population == "EASZ" ~ "06_EASZ",
                                                                    population == "KARA" ~ "07_KARA",
                                                                    population == "KETE" ~ "08_KETE",
                                                                    population == "MARC" ~ "09_MARC",
                                                                    population == "MARE" ~ "10_MARE",
                                                                    population == "NDAM" ~ "11_NDAM",
                                                                    population == "NGAN" ~ "12_NGAN",
                                                                    population == "ROMA" ~ "13_ROMA",
                                                                    population == "SHEK" ~ "14_SHEK",
                                                                    population == "SOMB" ~ "15_SOMB")),
                                            g.loc,
                                            chrnos)
  # write files for each chromosome
  for (chromosome in 1:29){
    # read snp information file
    snpinfo <- readr::read_table(paste0("07_local-ancestry-analysis/01_mosaic/02_output/",
                                        dplyr::case_when(data == "hd" ~ "01_hd/",
                                                         data == "ld" ~ "02_ld/"),
                                        dplyr::case_when(population == "ALEN" ~ "01_ALEN",
                                                         population == "ANKO" ~ "02_ANKO",
                                                         population == "BORA" ~ "03_BORA",
                                                         population == "BORG" ~ "04_BORG",
                                                         population == "CHIA" ~ "05_CHIA",
                                                         population == "EASZ" ~ "06_EASZ",
                                                         population == "KARA" ~ "07_KARA",
                                                         population == "KETE" ~ "08_KETE",
                                                         population == "MARC" ~ "09_MARC",
                                                         population == "MARE" ~ "10_MARE",
                                                         population == "NDAM" ~ "11_NDAM",
                                                         population == "NGAN" ~ "12_NGAN",
                                                         population == "ROMA" ~ "13_ROMA",
                                                         population == "SHEK" ~ "14_SHEK",
                                                         population == "SOMB" ~ "15_SOMB"),
                                        "/snpfile.",
                                        chromosome),
                                 col_names = F) %>%
      dplyr::rename(rs = X1,
                    chr = X2,
                    pos = X4) %>%
      dplyr::select(rs,
                    chr,
                    pos)
    # extract ancestries
    localanc_positions_chromosome <- localanc_positions[[chromosome]]
    mean_a <- list()
    mean_b <- list()
    mean_c <- list()
    for (i in 1:dim(localanc_positions_chromosome)[3]){
      mean_a[[i]] <- mean(localanc_positions_chromosome[order[1],
                                                        1:dim(localanc_positions_chromosome)[2],
                                                        i])
      mean_b[[i]] <- mean(localanc_positions_chromosome[order[2],
                                                        1:dim(localanc_positions_chromosome)[2],
                                                        i])
      mean_c[[i]] <- mean(localanc_positions_chromosome[order[3],
                                                        1:dim(localanc_positions_chromosome)[2],
                                                        i])
    }
    mean_a <- unlist(mean_a)
    mean_b <- unlist(mean_b)
    mean_c <- unlist(mean_c)
    # combine ancestries and save file
    tibble::tibble(snpinfo,
                   mean_a,
                   mean_b,
                   mean_c) %>%
      readr::write_csv(paste0("07_local-ancestry-analysis/01_mosaic/02_output/",
                              dplyr::case_when(data == "hd" ~ "01_hd/",
                                               data == "ld" ~ "02_ld/"),
                              dplyr::case_when(population == "ALEN" ~ "01_ALEN",
                                               population == "ANKO" ~ "02_ANKO",
                                               population == "BORA" ~ "03_BORA",
                                               population == "BORG" ~ "04_BORG",
                                               population == "CHIA" ~ "05_CHIA",
                                               population == "EASZ" ~ "06_EASZ",
                                               population == "KARA" ~ "07_KARA",
                                               population == "KETE" ~ "08_KETE",
                                               population == "MARC" ~ "09_MARC",
                                               population == "MARE" ~ "10_MARE",
                                               population == "NDAM" ~ "11_NDAM",
                                               population == "NGAN" ~ "12_NGAN",
                                               population == "ROMA" ~ "13_ROMA",
                                               population == "SHEK" ~ "14_SHEK",
                                               population == "SOMB" ~ "15_SOMB"),
                              "/mosaic-",
                              data,
                              "-",
                              population,
                              "-",
                              stringr::str_pad(chromosome,
                                               2,
                                               pad = "0",
                                               side = "left"),
                              ".csv"))
  }
  # get chromosome 1
  join <- readr::read_csv(paste0("07_local-ancestry-analysis/01_mosaic/02_output/",
                                 dplyr::case_when(data == "hd" ~ "01_hd/",
                                                  data == "ld" ~ "02_ld/"),
                                 dplyr::case_when(population == "ALEN" ~ "01_ALEN",
                                                  population == "ANKO" ~ "02_ANKO",
                                                  population == "BORA" ~ "03_BORA",
                                                  population == "BORG" ~ "04_BORG",
                                                  population == "CHIA" ~ "05_CHIA",
                                                  population == "EASZ" ~ "06_EASZ",
                                                  population == "KARA" ~ "07_KARA",
                                                  population == "KETE" ~ "08_KETE",
                                                  population == "MARC" ~ "09_MARC",
                                                  population == "MARE" ~ "10_MARE",
                                                  population == "NDAM" ~ "11_NDAM",
                                                  population == "NGAN" ~ "12_NGAN",
                                                  population == "ROMA" ~ "13_ROMA",
                                                  population == "SHEK" ~ "14_SHEK",
                                                  population == "SOMB" ~ "15_SOMB"),
                                 "/mosaic-",
                                 data,
                                 "-",
                                 population,
                                 "-01.csv"))
  # add chromosomes 2 - 29
  for(i in 2:29){
    join <- join %>%
      dplyr::full_join(readr::read_csv(paste0("07_local-ancestry-analysis/01_mosaic/02_output/",
                                              dplyr::case_when(data == "hd" ~ "01_hd/",
                                                               data == "ld" ~ "02_ld/"),
                                              dplyr::case_when(population == "ALEN" ~ "01_ALEN",
                                                               population == "ANKO" ~ "02_ANKO",
                                                               population == "BORA" ~ "03_BORA",
                                                               population == "BORG" ~ "04_BORG",
                                                               population == "CHIA" ~ "05_CHIA",
                                                               population == "EASZ" ~ "06_EASZ",
                                                               population == "KARA" ~ "07_KARA",
                                                               population == "KETE" ~ "08_KETE",
                                                               population == "MARC" ~ "09_MARC",
                                                               population == "MARE" ~ "10_MARE",
                                                               population == "NDAM" ~ "11_NDAM",
                                                               population == "NGAN" ~ "12_NGAN",
                                                               population == "ROMA" ~ "13_ROMA",
                                                               population == "SHEK" ~ "14_SHEK",
                                                               population == "SOMB" ~ "15_SOMB"),
                                              "/mosaic-",
                                              data,
                                              "-",
                                              population,
                                              "-",
                                              stringr::str_pad(i,
                                                               2,
                                                               pad = "0",
                                                               side = "left"),
                                              ".csv")))
  }
  # format snp positions, add row numbers and z-scores and save file
  join %>%
    dplyr::mutate(row_num = dplyr::row_number(),
                  z_a = as.numeric(scale(mean_a)),
                  z_b = as.numeric(scale(mean_b)),
                  z_c = as.numeric(scale(mean_c))) %>%
    readr::write_csv(paste0("07_local-ancestry-analysis/01_mosaic/02_output/",
                            dplyr::case_when(data == "hd" ~ "01_hd/",
                                             data == "ld" ~ "02_ld/"),
                            dplyr::case_when(population == "ALEN" ~ "01_ALEN",
                                             population == "ANKO" ~ "02_ANKO",
                                             population == "BORA" ~ "03_BORA",
                                             population == "BORG" ~ "04_BORG",
                                             population == "CHIA" ~ "05_CHIA",
                                             population == "EASZ" ~ "06_EASZ",
                                             population == "KARA" ~ "07_KARA",
                                             population == "KETE" ~ "08_KETE",
                                             population == "MARC" ~ "09_MARC",
                                             population == "MARE" ~ "10_MARE",
                                             population == "NDAM" ~ "11_NDAM",
                                             population == "NGAN" ~ "12_NGAN",
                                             population == "ROMA" ~ "13_ROMA",
                                             population == "SHEK" ~ "14_SHEK",
                                             population == "SOMB" ~ "15_SOMB"),
                            "/mosaic-",
                            data,
                            "-",
                            population,
                            ".csv"))
}

### function to read mosaic data and extract mean ancestry results for easz
mosaic_data_easz <- function(data = "hd"){
  # list rdata files containing model parameters and local ancestry results
  files <- list.files(paste0("07_local-ancestry-analysis/01_mosaic/02_output/",
                             dplyr::case_when(data == "hd" ~ "01_hd/",
                                              data == "ld" ~ "02_ld/"),
                             "06_EASZ"),
                      full.names = T,
                      pattern = "*.RData")
  # load file 1
  load(files[1])
  # load file 3
  load(files[3])
  # get local ancestry results for snp positions for run 1
  run_01_localanc_positions <- MOSAIC::grid_to_pos(localanc,
                                                   paste0("07_local-ancestry-analysis/01_mosaic/02_output/",
                                                          dplyr::case_when(data == "hd" ~ "01_hd/",
                                                                           data == "ld" ~ "02_ld/"),
                                                          "06_EASZ"),
                                                   g.loc,
                                                   chrnos)
  # load file 2
  load(files[2])
  # load file 4
  load(files[4])
  # get local ancestry results for snp positions for run 2
  run_02_localanc_positions <- MOSAIC::grid_to_pos(localanc,
                                                   paste0("07_local-ancestry-analysis/01_mosaic/02_output/",
                                                          dplyr::case_when(data == "hd" ~ "01_hd/",
                                                                           data == "ld" ~ "02_ld/"),
                                                          "06_EASZ"),
                                                   g.loc,
                                                   chrnos)
  # write files for each chromosome
  for (chromosome in 1:29){
    # read snp information file
    snpinfo <- readr::read_table(paste0("07_local-ancestry-analysis/01_mosaic/02_output/",
                                        dplyr::case_when(data == "hd" ~ "01_hd/",
                                                         data == "ld" ~ "02_ld/"),
                                        "06_EASZ/snpfile.",
                                        chromosome),
                                 col_names = F) %>%
      dplyr::rename(rs = X1,
                    chr = X2,
                    pos = X4) %>%
      dplyr::select(rs,
                    chr,
                    pos)
    # extract ancestries for run 1
    run_01_localanc_positions_chromosome <- run_01_localanc_positions[[chromosome]]
    run_01_mean_a <- list()
    run_01_mean_b <- list()
    run_01_mean_c <- list()
    for (i in 1:dim(run_01_localanc_positions_chromosome)[3]){
      run_01_mean_a[[i]] <- mean(run_01_localanc_positions_chromosome[1,
                                                                      1:dim(run_01_localanc_positions_chromosome)[2],
                                                                      i])
      run_01_mean_b[[i]] <- mean(run_01_localanc_positions_chromosome[3,
                                                                      1:dim(run_01_localanc_positions_chromosome)[2],
                                                                      i])
      run_01_mean_c[[i]] <- mean(run_01_localanc_positions_chromosome[2,
                                                                      1:dim(run_01_localanc_positions_chromosome)[2],
                                                                      i])
    }
    run_01_mean_a <- unlist(run_01_mean_a)
    run_01_mean_b <- unlist(run_01_mean_b)
    run_01_mean_c <- unlist(run_01_mean_c)
    # extract ancestries for run 2
    run_02_localanc_positions_chromosome <- run_02_localanc_positions[[chromosome]]
    run_02_mean_a <- list()
    run_02_mean_b <- list()
    run_02_mean_c <- list()
    for (i in 1:dim(run_02_localanc_positions_chromosome)[3]){
      run_02_mean_a[[i]] <- mean(run_02_localanc_positions_chromosome[1,
                                                                      1:dim(run_02_localanc_positions_chromosome)[2],
                                                                      i])
      run_02_mean_b[[i]] <- mean(run_02_localanc_positions_chromosome[3,
                                                                      1:dim(run_02_localanc_positions_chromosome)[2],
                                                                      i])
      run_02_mean_c[[i]] <- mean(run_02_localanc_positions_chromosome[2,
                                                                      1:dim(run_02_localanc_positions_chromosome)[2],
                                                                      i])
    }
    run_02_mean_a <- unlist(run_02_mean_a)
    run_02_mean_b <- unlist(run_02_mean_b)
    run_02_mean_c <- unlist(run_02_mean_c)
    # combine runs
    means <- tibble::tibble(run_01_mean_a,
                            run_01_mean_b,
                            run_01_mean_c,
                            run_02_mean_a,
                            run_02_mean_b,
                            run_02_mean_c) %>%
      dplyr::mutate(run_01_num_samples = 50,
                    run_02_num_samples = 61) %>%
      dplyr::mutate(mean_a = rowSums(across(ends_with("mean_a")) * across(ends_with("num_samples"))) / rowSums(across(ends_with("num_samples")))) %>%
      dplyr::mutate(mean_b = rowSums(across(ends_with("mean_b")) * across(ends_with("num_samples"))) / rowSums(across(ends_with("num_samples")))) %>%
      dplyr::mutate(mean_c = rowSums(across(ends_with("mean_c")) * across(ends_with("num_samples"))) / rowSums(across(ends_with("num_samples")))) %>%
      dplyr::select(mean_a,
                    mean_b,
                    mean_c)
    # combine with snp information and save file
    tibble::tibble(snpinfo,
                   means) %>%
      readr::write_csv(paste0("07_local-ancestry-analysis/01_mosaic/02_output/",
                              dplyr::case_when(data == "hd" ~ "01_hd/",
                                               data == "ld" ~ "02_ld/"),
                              "06_EASZ/mosaic-",
                              data,
                              "-EASZ-",
                              stringr::str_pad(chromosome,
                                               2,
                                               pad = "0",
                                               side = "left"),
                              ".csv"))
  }
  # get chromosome 1
  join <- readr::read_csv(paste0("07_local-ancestry-analysis/01_mosaic/02_output/",
                                 dplyr::case_when(data == "hd" ~ "01_hd/",
                                                  data == "ld" ~ "02_ld/"),
                                 "06_EASZ/mosaic-",
                                 data,
                                 "-EASZ-01.csv"))
  # add chromosomes 2 - 29
  for(i in 2:29){
    join <- join %>%
      dplyr::full_join(readr::read_csv(paste0("07_local-ancestry-analysis/01_mosaic/02_output/",
                                              dplyr::case_when(data == "hd" ~ "01_hd/",
                                                               data == "ld" ~ "02_ld/"),
                                              "06_EASZ/mosaic-",
                                              data,
                                              "-EASZ-",
                                              stringr::str_pad(i,
                                                               2,
                                                               pad = "0",
                                                               side = "left"),
                                              ".csv")))
  }
  # format snp positions, add row numbers and z-scores and save file
  join %>%
    dplyr::mutate(row_num = dplyr::row_number(),
                  z_a = as.numeric(scale(mean_a)),
                  z_b = as.numeric(scale(mean_b)),
                  z_c = as.numeric(scale(mean_c))) %>%
    readr::write_csv(paste0("07_local-ancestry-analysis/01_mosaic/02_output/",
                            dplyr::case_when(data == "hd" ~ "01_hd/",
                                             data == "ld" ~ "02_ld/"),
                            "06_EASZ/mosaic-",
                            data,
                            "-EASZ.csv"))
}

### function to read elai chromosome data and extract mean ancestry results
elai_chromosome_data <- function(chromosome, population, order = c(1, 2, 3), data = "hd"){
  # read mean local ancestry dosage file
  ps21 <- readr::read_table(paste0("07_local-ancestry-analysis/02_elai/02_output/",
                                   dplyr::case_when(data == "hd" ~ "01_hd/",
                                                    data == "ld" ~ "02_ld/"),
                                   dplyr::case_when(population == "ALEN" ~ "01_ALEN",
                                                    population == "ANKO" ~ "02_ANKO",
                                                    population == "BORA" ~ "03_BORA",
                                                    population == "BORG" ~ "04_BORG",
                                                    population == "CHIA" ~ "05_CHIA",
                                                    population == "EASZ" ~ "06_EASZ",
                                                    population == "KARA" ~ "07_KARA",
                                                    population == "KETE" ~ "08_KETE",
                                                    population == "MARC" ~ "09_MARC",
                                                    population == "MARE" ~ "10_MARE",
                                                    population == "NDAM" ~ "11_NDAM",
                                                    population == "NGAN" ~ "12_NGAN",
                                                    population == "ROMA" ~ "13_ROMA",
                                                    population == "SHEK" ~ "14_SHEK",
                                                    population == "SOMB" ~ "15_SOMB"),
                                   "/output/",
                                   data,
                                   "-",
                                   population,
                                   "-",
                                   chromosome,
                                   ".ps21.txt"),
                            col_names = F) %>%
    dplyr::mutate(across(everything(),
                         ~./2))
  # read snp information file
  snpinfo <- readr::read_table(paste0("07_local-ancestry-analysis/02_elai/02_output/",
                                      dplyr::case_when(data == "hd" ~ "01_hd/",
                                                       data == "ld" ~ "02_ld/"),
                                      dplyr::case_when(population == "ALEN" ~ "01_ALEN",
                                                       population == "ANKO" ~ "02_ANKO",
                                                       population == "BORA" ~ "03_BORA",
                                                       population == "BORG" ~ "04_BORG",
                                                       population == "CHIA" ~ "05_CHIA",
                                                       population == "EASZ" ~ "06_EASZ",
                                                       population == "KARA" ~ "07_KARA",
                                                       population == "KETE" ~ "08_KETE",
                                                       population == "MARC" ~ "09_MARC",
                                                       population == "MARE" ~ "10_MARE",
                                                       population == "NDAM" ~ "11_NDAM",
                                                       population == "NGAN" ~ "12_NGAN",
                                                       population == "ROMA" ~ "13_ROMA",
                                                       population == "SHEK" ~ "14_SHEK",
                                                       population == "SOMB" ~ "15_SOMB"),
                                      "/output/",
                                      data,
                                      "-",
                                      population,
                                      "-",
                                      chromosome,
                                      ".snpinfo.txt")) %>%
    dplyr::mutate(chr = chromosome) %>%
    dplyr::select(rs,
                  chr,
                  pos)
  # extract ancestries
  a <- t(ps21[seq(order[1],
                  length(ps21),
                  3)])
  b <- t(ps21[seq(order[2],
                  length(ps21),
                  3)])
  c <- t(ps21[seq(order[3],
                  length(ps21),
                  3)])
  # calculate mean ancestry
  mean_a <- apply(a,
                  1,
                  mean)
  mean_b <- apply(b,
                  1,
                  mean)
  mean_c <- apply(c,
                  1,
                  mean)
  # ensure correct length
  mean_a <- mean_a[1:nrow(snpinfo)]
  mean_b <- mean_b[1:nrow(snpinfo)]
  mean_c <- mean_c[1:nrow(snpinfo)]
  # combine ancestries and save file
  tibble::tibble(snpinfo,
                 mean_a,
                 mean_b,
                 mean_c) %>%
    readr::write_csv(paste0("07_local-ancestry-analysis/02_elai/02_output/",
                            dplyr::case_when(data == "hd" ~ "01_hd/",
                                             data == "ld" ~ "02_ld/"),
                            dplyr::case_when(population == "ALEN" ~ "01_ALEN",
                                             population == "ANKO" ~ "02_ANKO",
                                             population == "BORA" ~ "03_BORA",
                                             population == "BORG" ~ "04_BORG",
                                             population == "CHIA" ~ "05_CHIA",
                                             population == "EASZ" ~ "06_EASZ",
                                             population == "KARA" ~ "07_KARA",
                                             population == "KETE" ~ "08_KETE",
                                             population == "MARC" ~ "09_MARC",
                                             population == "MARE" ~ "10_MARE",
                                             population == "NDAM" ~ "11_NDAM",
                                             population == "NGAN" ~ "12_NGAN",
                                             population == "ROMA" ~ "13_ROMA",
                                             population == "SHEK" ~ "14_SHEK",
                                             population == "SOMB" ~ "15_SOMB"),
                            "/elai-",
                            data,
                            "-",
                            population,
                            "-",
                            stringr::str_pad(chromosome,
                                             2,
                                             pad = "0",
                                             side = "left"),
                            ".csv"))
}

### function to combine elai chromosomes into genome data
elai_genome_data <- function(population, data = "hd"){
  # get chromosome 1
  join <- readr::read_csv(paste0("07_local-ancestry-analysis/02_elai/02_output/",
                                 dplyr::case_when(data == "hd" ~ "01_hd/",
                                                  data == "ld" ~ "02_ld/"),
                                 dplyr::case_when(population == "ALEN" ~ "01_ALEN",
                                                  population == "ANKO" ~ "02_ANKO",
                                                  population == "BORA" ~ "03_BORA",
                                                  population == "BORG" ~ "04_BORG",
                                                  population == "CHIA" ~ "05_CHIA",
                                                  population == "EASZ" ~ "06_EASZ",
                                                  population == "KARA" ~ "07_KARA",
                                                  population == "KETE" ~ "08_KETE",
                                                  population == "MARC" ~ "09_MARC",
                                                  population == "MARE" ~ "10_MARE",
                                                  population == "NDAM" ~ "11_NDAM",
                                                  population == "NGAN" ~ "12_NGAN",
                                                  population == "ROMA" ~ "13_ROMA",
                                                  population == "SHEK" ~ "14_SHEK",
                                                  population == "SOMB" ~ "15_SOMB"),
                                 "/elai-",
                                 data,
                                 "-",
                                 population,
                                 "-01.csv"))
  # add chromosomes 2 - 29
  for(i in 2:29){
    join <- join %>%
      dplyr::full_join(readr::read_csv(paste0("07_local-ancestry-analysis/02_elai/02_output/",
                                              dplyr::case_when(data == "hd" ~ "01_hd/",
                                                               data == "ld" ~ "02_ld/"),
                                              dplyr::case_when(population == "ALEN" ~ "01_ALEN",
                                                               population == "ANKO" ~ "02_ANKO",
                                                               population == "BORA" ~ "03_BORA",
                                                               population == "BORG" ~ "04_BORG",
                                                               population == "CHIA" ~ "05_CHIA",
                                                               population == "EASZ" ~ "06_EASZ",
                                                               population == "KARA" ~ "07_KARA",
                                                               population == "KETE" ~ "08_KETE",
                                                               population == "MARC" ~ "09_MARC",
                                                               population == "MARE" ~ "10_MARE",
                                                               population == "NDAM" ~ "11_NDAM",
                                                               population == "NGAN" ~ "12_NGAN",
                                                               population == "ROMA" ~ "13_ROMA",
                                                               population == "SHEK" ~ "14_SHEK",
                                                               population == "SOMB" ~ "15_SOMB"),
                                              "/elai-",
                                              data,
                                              "-",
                                              population,
                                              "-",
                                              stringr::str_pad(i,
                                                               2,
                                                               pad = "0",
                                                               side = "left"),
                                              ".csv")))
  }
  # format snp positions, add row numbers and z-scores and save file
  join %>%
    dplyr::mutate(row_num = dplyr::row_number(),
                  z_a = as.numeric(scale(mean_a)),
                  z_b = as.numeric(scale(mean_b)),
                  z_c = as.numeric(scale(mean_c))) %>%
    readr::write_csv(paste0("07_local-ancestry-analysis/02_elai/02_output/",
                            dplyr::case_when(data == "hd" ~ "01_hd/",
                                             data == "ld" ~ "02_ld/"),
                            dplyr::case_when(population == "ALEN" ~ "01_ALEN",
                                             population == "ANKO" ~ "02_ANKO",
                                             population == "BORA" ~ "03_BORA",
                                             population == "BORG" ~ "04_BORG",
                                             population == "CHIA" ~ "05_CHIA",
                                             population == "EASZ" ~ "06_EASZ",
                                             population == "KARA" ~ "07_KARA",
                                             population == "KETE" ~ "08_KETE",
                                             population == "MARC" ~ "09_MARC",
                                             population == "MARE" ~ "10_MARE",
                                             population == "NDAM" ~ "11_NDAM",
                                             population == "NGAN" ~ "12_NGAN",
                                             population == "ROMA" ~ "13_ROMA",
                                             population == "SHEK" ~ "14_SHEK",
                                             population == "SOMB" ~ "15_SOMB"),
                            "/elai-",
                            data,
                            "-",
                            population,
                            ".csv"))
}

### function to combine genome data for multiple populations
combine_populations <- function(populations, software = "mosaic", data = "hd"){
  # get population 1 and rename columns
  all_populations <- get(populations)
  join <- readr::read_csv(paste0("07_local-ancestry-analysis/",
                                 dplyr::case_when(software == "mosaic" ~ "01_mosaic/",
                                                  software == "elai" ~ "02_elai/"),
                                 "02_output/",
                                 dplyr::case_when(data == "hd" ~ "01_hd/",
                                                  data == "ld" ~ "02_ld/"),
                                 dplyr::case_when(all_populations[1] == "ALEN" ~ "01_ALEN/",
                                                  all_populations[1] == "ANKO" ~ "02_ANKO/",
                                                  all_populations[1] == "BORA" ~ "03_BORA/",
                                                  all_populations[1] == "BORG" ~ "04_BORG/",
                                                  all_populations[1] == "CHIA" ~ "05_CHIA/",
                                                  all_populations[1] == "EASZ" ~ "06_EASZ/",
                                                  all_populations[1] == "KARA" ~ "07_KARA/",
                                                  all_populations[1] == "KETE" ~ "08_KETE/",
                                                  all_populations[1] == "MARC" ~ "09_MARC/",
                                                  all_populations[1] == "MARE" ~ "10_MARE/",
                                                  all_populations[1] == "NDAM" ~ "11_NDAM/",
                                                  all_populations[1] == "NGAN" ~ "12_NGAN/",
                                                  all_populations[1] == "ROMA" ~ "13_ROMA/",
                                                  all_populations[1] == "SHEK" ~ "14_SHEK/",
                                                  all_populations[1] == "SOMB" ~ "15_SOMB/"),
                                 software,
                                 "-",
                                 data,
                                 "-",
                                 all_populations[1],
                                 ".csv")) %>%
    dplyr::mutate(num_samples = dplyr::case_when(all_populations[1] == "ALEN" ~ 6,
                                                 all_populations[1] == "ANKO" ~ 25,
                                                 all_populations[1] == "BORA" ~ 90,
                                                 all_populations[1] == "BORG" ~ 50,
                                                 all_populations[1] == "CHIA" ~ 19,
                                                 all_populations[1] == "EASZ" ~ 111,
                                                 all_populations[1] == "KARA" ~ 16,
                                                 all_populations[1] == "KETE" ~ 22,
                                                 all_populations[1] == "MARC" ~ 13,
                                                 all_populations[1] == "MARE" ~ 5,
                                                 all_populations[1] == "NDAM" ~ 41,
                                                 all_populations[1] == "NGAN" ~ 27,
                                                 all_populations[1] == "ROMA" ~ 51,
                                                 all_populations[1] == "SHEK" ~ 16,
                                                 all_populations[1] == "SOMB" ~ 23)) %>%
    dplyr::rename_with(~ paste0(all_populations[1],
                                "_",
                                .),
                       - c(rs,
                           pos,
                           chr,
                           row_num))
  # add other populations
  for(p in all_populations[-1]){
    add <- readr::read_csv(paste0("07_local-ancestry-analysis/",
                                  dplyr::case_when(software == "mosaic" ~ "01_mosaic/",
                                                   software == "elai" ~ "02_elai/"),
                                  "02_output/",
                                  dplyr::case_when(data == "hd" ~ "01_hd/",
                                                   data == "ld" ~ "02_ld/"),
                                  dplyr::case_when(p == "ALEN" ~ "01_ALEN/",
                                                   p == "ANKO" ~ "02_ANKO/",
                                                   p == "BORA" ~ "03_BORA/",
                                                   p == "BORG" ~ "04_BORG/",
                                                   p == "CHIA" ~ "05_CHIA/",
                                                   p == "EASZ" ~ "06_EASZ/",
                                                   p == "KARA" ~ "07_KARA/",
                                                   p == "KETE" ~ "08_KETE/",
                                                   p == "MARC" ~ "09_MARC/",
                                                   p == "MARE" ~ "10_MARE/",
                                                   p == "NDAM" ~ "11_NDAM/",
                                                   p == "NGAN" ~ "12_NGAN/",
                                                   p == "ROMA" ~ "13_ROMA/",
                                                   p == "SHEK" ~ "14_SHEK/",
                                                   p == "SOMB" ~ "15_SOMB/"),
                                  software,
                                  "-",
                                  data,
                                  "-",
                                  p,
                                  ".csv")) %>%
      dplyr::mutate(num_samples = dplyr::case_when(p == "ALEN" ~ 6,
                                                   p == "ANKO" ~ 25,
                                                   p == "BORA" ~ 90,
                                                   p == "BORG" ~ 50,
                                                   p == "CHIA" ~ 19,
                                                   p == "EASZ" ~ 111,
                                                   p == "KARA" ~ 16,
                                                   p == "KETE" ~ 22,
                                                   p == "MARC" ~ 13,
                                                   p == "MARE" ~ 5,
                                                   p == "NDAM" ~ 41,
                                                   p == "NGAN" ~ 27,
                                                   p == "ROMA" ~ 51,
                                                   p == "SHEK" ~ 16,
                                                   p == "SOMB" ~ 23)) %>%
      dplyr::rename_with(~ paste0(p,
                                  "_",
                                  .),
                         - c(rs,
                             pos,
                             chr,
                             row_num))
    join <- join %>%
      dplyr::full_join(add)
  }
  # calculate weighted mean ancestries and z-scores and save file
  join <- join %>%
    dplyr::mutate(mean_a = rowSums(across(ends_with("mean_a")) * across(ends_with("num_samples"))) / rowSums(across(ends_with("num_samples")))) %>%
    dplyr::mutate(mean_b = rowSums(across(ends_with("mean_b")) * across(ends_with("num_samples"))) / rowSums(across(ends_with("num_samples")))) %>%
    dplyr::mutate(mean_c = rowSums(across(ends_with("mean_c")) * across(ends_with("num_samples"))) / rowSums(across(ends_with("num_samples")))) %>%
    dplyr::mutate(z_a = as.numeric(scale(mean_a)),
                  z_b = as.numeric(scale(mean_b)),
                  z_c = as.numeric(scale(mean_c))) %>%
    readr::write_csv(paste0("07_local-ancestry-analysis/",
                            dplyr::case_when(software == "mosaic" ~ "01_mosaic/",
                                             software == "elai" ~ "02_elai/"),
                            "02_output/",
                            dplyr::case_when(data == "hd" ~ "01_hd/",
                                             data == "ld" ~ "02_ld/"),
                            "16_combined/",
                            dplyr::case_when(stringr::str_detect(populations,
                                                                 "european") ~ "01_european-hybrids/",
                                             stringr::str_detect(populations,
                                                                 "trypanotolerant_african") ~ "02_trypanotolerant-african-hybrids/",
                                             stringr::str_detect(populations,
                                                                 "trypanosusceptible_african") ~ "03_trypanosusceptible-african-hybrids/",
                                             .default = ""),
                            software,
                            "-",
                            data,
                            "-",
                            stringr::str_replace_all(populations,
                                                     "_",
                                                     "-"),
                            ".csv"))
  # separate combined genomes into chromosomes and save files
  for (i in 1:29){
    join %>%
      dplyr::filter(chr == i) %>%
      readr::write_csv(paste0("07_local-ancestry-analysis/",
                              dplyr::case_when(software == "mosaic" ~ "01_mosaic/",
                                               software == "elai" ~ "02_elai/"),
                              "02_output/",
                              dplyr::case_when(data == "hd" ~ "01_hd/",
                                               data == "ld" ~ "02_ld/"),
                              "16_combined/",
                              dplyr::case_when(stringr::str_detect(populations,
                                                                   "european") ~ "01_european-hybrids/",
                                               stringr::str_detect(populations,
                                                                   "trypanotolerant_african") ~ "02_trypanotolerant-african-hybrids/",
                                               stringr::str_detect(populations,
                                                                   "trypanosusceptible_african") ~ "03_trypanosusceptible-african-hybrids/",
                                               .default = ""),
                              software,
                              "-",
                              data,
                              "-",
                              stringr::str_replace_all(populations,
                                                       "_",
                                                       "-"),
                              "-",
                              stringr::str_pad(i,
                                               2,
                                               pad = "0",
                                               side = "left"),
                              ".csv"))
  }
}

### apply functions to data
#### set lists of analyses
all_hybrids <- c("ROMA", "CHIA", "MARC", "MARE", "ALEN",
                 "NDAM", "BORG", "SOMB", "KETE", "SHEK",
                 "ANKO", "NGAN", "EASZ", "KARA", "BORA")

european_hybrids <- c("ROMA", "CHIA", "MARC", "MARE", "ALEN")

african_hybrids <- c("NDAM", "BORG", "SOMB", "KETE", "SHEK",
                     "ANKO", "NGAN", "EASZ", "KARA", "BORA")

trypanotolerant_african_hybrids <- c("NDAM", "BORG", "SOMB", "KETE", "SHEK")

trypanosusceptible_african_hybrids <- c("ANKO", "NGAN", "EASZ", "KARA", "BORA")

selected_european_hybrids <- c("ROMA", "CHIA")

selected_trypanotolerant_african_hybrids <- c("BORG", "SHEK")

combinations <- c("european_hybrids", "trypanotolerant_african_hybrids", "trypanosusceptible_african_hybrids",
                  "selected_european_hybrids", "selected_trypanotolerant_african_hybrids")

#### read mosaic data and extract mean ancestry results
##### european hybrids
parallel::mclapply(european_hybrids,
                   mosaic_data,
                   mc.cores = 15)
parallel::mclapply(european_hybrids,
                   mosaic_data,
                   data = "ld",
                   mc.cores = 15)

##### african hybrids
parallel::mclapply(african_hybrids,
                   mosaic_data,
                   order = c(1, 3, 2),
                   mc.cores = 15)
parallel::mclapply(african_hybrids,
                   mosaic_data,
                   order = c(1, 3, 2),
                   data = "ld",
                   mc.cores = 15)

mosaic_data_easz()
mosaic_data_easz(data = "ld")

#### read elai chromosome data and extract mean ancestry results
##### european hybrids
for (h in european_hybrids){
  parallel::mclapply(1:29,
                     elai_chromosome_data,
                     population = h,
                     mc.cores = 15)
  parallel::mclapply(1:29,
                     elai_chromosome_data,
                     population = h,
                     data = "ld",
                     mc.cores = 15)
}

##### african hybrids
for (h in african_hybrids){
  parallel::mclapply(1:29,
                     elai_chromosome_data,
                     population = h,
                     order = c(1, 3, 2),
                     mc.cores = 15)
  parallel::mclapply(1:29,
                     elai_chromosome_data,
                     population = h,
                     order = c(1, 3, 2),
                     data = "ld",
                     mc.cores = 15)
}

#### combine elai chromosome data into genome data for all hybrids
parallel::mclapply(all_hybrids,
                   elai_genome_data,
                   mc.cores = 15)
parallel::mclapply(all_hybrids,
                   elai_genome_data,
                   data = "ld",
                   mc.cores = 15)

#### combine genome data for multiple hybrid combinations
parallel::mclapply(combinations,
                   combine_populations,
                   mc.cores = 15)
parallel::mclapply(combinations,
                   combine_populations,
                   data = "ld",
                   mc.cores = 15)

parallel::mclapply(combinations,
                   combine_populations,
                   software = "elai",
                   mc.cores = 15)
parallel::mclapply(combinations,
                   combine_populations,
                   software = "elai",
                   data = "ld",
                   mc.cores = 15)