# functional enrichment analysis
## functional enrichment analysis of local ancestry results for hd and ld snp data with mosaic and elai
### run functional enrichment analysis
nohup bash -c 'Rscript 08_functional-enrichment-analysis/functional-enrichment.R' > 08_functional-enrichment/log.txt &

### combine enrichment plots
nohup bash -c 'Rscript 08_functional-enrichment-analysis/03_figures/combine-enrichment-plots.R' > /dev/null 2>&1 &
