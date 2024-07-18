# update genome assembly
## update genome assembly of snp data to ars-ucd1.2
### remove "-" from snps without identified reference allele in downloaded schnabel file with sed
sed 's/\t-/\t/g' 01_data-sources/02_downloaded/09_schnabel/UMC_marker_names_180910/9913_ARS1.2_777962_HD_marker_name_180910.map > 01_data-sources/02_downloaded/09_schnabel/9913_ARS1.2_777962_HD_marker_name_180910-edit.map

### update genome assembly using edited file with plink
plink --cow --bfile 01_data-sources/03_merged/merged --update-chr 01_data-sources/02_downloaded/09_schnabel/9913_ARS1.2_777962_HD_marker_name_180910-edit.map 1 2 --update-map 01_data-sources/02_downloaded/09_schnabel/9913_ARS1.2_777962_HD_marker_name_180910-edit.map 4 2 --make-bed --out 02_update-genome-assembly/ars