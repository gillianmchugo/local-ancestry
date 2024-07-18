# recombination rates
## prepare recombination rate files for mosaic analysis
## adapted from https://github.com/jamesandyward/MOSAIC
### read recombination map file
cattle_rmap <- read.table("01_data-sources/02_downloaded/10_ma/doi_10.5061_dryad.q2q84__v1/cattle_rmap.txt",
                          header = T)

### function to write recombination rate files for each chromosome
rates <- function(chr) {
  cattle_rmap$mean_map <- (cattle_rmap$map_f + cattle_rmap$map_m)/2
  chr.rates <- cattle_rmap[cattle_rmap$Chr == chr,]
  chr.rates$cum_rates <- cumsum(chr.rates$mean_map)*100
  sites <- matrix(,
                  nrow = 1,
                  ncol = length(chr.rates$Name))
  pos <- t(chr.rates$Location)
  r_rate <- t(chr.rates$cum_rates)
  d <- list(sites,
            pos,
            r_rate)
  x <- do.call(rbind,
               d)
  x[1,1] <- paste(":sites:",
                  length(chr.rates$Name),
                  sep = "")
  x[is.na(x)] <- " "
  write.table(x,
              quote = F,
              row.names = F,
              col.names = F,
              file = paste0("07_local-ancestry-analysis/01_mosaic/01_input/01_hd/rates.",
                            chr))
  write.table(x,
              quote = F,
              row.names = F,
              col.names = F,
              file = paste0("07_local-ancestry-analysis/01_mosaic/01_input/02_ld/rates.",
                            chr))
}

### apply function
lapply(1:29,
       rates)