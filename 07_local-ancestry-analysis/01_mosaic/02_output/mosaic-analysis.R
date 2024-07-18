# mosaic analysis
## local ancestry analysis with mosaic
### install required packages
if (!require(devtools)) install.packages("devtools")
if (!require(dplyr)) install.packages("dplyr")
if (!require(MOSAIC)) devtools::install_github("mststats/MOSAIC")
if (!require(stringr)) install.packages("stringr")

### function to run mosaic analysis
mosaic_analysis <- function(i,
                            data = "hd"){
  # set breeds for analysis
  breed <- c("ALEN", "ANKO", "BORA", "BORG", "CHIA",
             "EASZ", "KARA", "KETE", "MARC", "MARE",
             "NDAM", "NGAN", "ROMA", "SHEK", "SOMB")
  # set sample numbers for analysis
  samples <- c(6, 25, 90, 50, 19,
               50, 16, 22, 13, 5,
               41, 27, 51, 16, 23)
  # run mosaic analysis
  MOSAIC::run_mosaic(A = 3,
                     chrnos = 1:29,
                     datasource = paste0("07_local-ancestry-analysis/01_mosaic/01_input/",
                                         dplyr::case_when(data == "hd" ~ "01_hd",
                                                          data == "ld" ~ "02_ld")),
                     doFst = F,
                     MC = 10,
                     Ne = 400,
                     NUMI = samples[i],
                     pops = c("ANGU", "GIR", "HOLS",
                              "JERS", "LAGU", "MUTU",
                              "NDAG", "NELO", "THAR"),
                     resultsdir = paste0("07_local-ancestry-analysis/01_mosaic/02_output/",
                                         dplyr::case_when(data == "hd" ~ "01_hd/",
                                                          data == "ld" ~ "02_ld/"),
                                         stringr::str_pad(i,
                                                          2,
                                                          side = "left",
                                                          pad = "0"),
                                         "_",
                                         breed[i]),
                     target = breed[i])
}

### apply function to hd and ld data for breeds 1 to 15
lapply(1:15,
       mosaic_analysis)
lapply(1:15,
       mosaic_analysis,
       data = "ld")

### run mosaic for second batch of easz samples for hd and ld data
MOSAIC::run_mosaic(A = 3,
                   chrnos = 1:29,
                   datasource = "07_local-ancestry-analysis/01_mosaic/01_input/01_hd",
                   doFst = F,
                   firstind = 51,
                   MC = 40,
                   Ne = 400,
                   NUMI = 61,
                   pops = c("ANGU", "GIR", "HOLS",
                            "JERS", "LAGU", "MUTU",
                            "NDAG", "NELO", "THAR"),
                   resultsdir = "07_local-ancestry-analysis/01_mosaic/02_output/01_hd/06_EASZ",
                   target = "EASZ")

MOSAIC::run_mosaic(A = 3,
                   chrnos = 1:29,
                   datasource = "07_local-ancestry-analysis/01_mosaic/01_input/02_ld",
                   doFst = F,
                   firstind = 51,
                   MC = 40,
                   Ne = 400,
                   NUMI = 61,
                   pops = c("ANGU", "GIR", "HOLS",
                            "JERS", "LAGU", "MUTU",
                            "NDAG", "NELO", "THAR"),
                   resultsdir = "07_local-ancestry-analysis/01_mosaic/02_output/02_ld/06_EASZ",
                   target = "EASZ")