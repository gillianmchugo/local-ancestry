# phylogenetic analysis
## phylogenetic analysis with treemix
### prepare outgroup data
#### unzip downloaded files with unzip
unzip 01_data-sources/02_downloaded/06_widde/outgroup/cattle_*_PLINK.zip -d 01_data-sources/02_downloaded/06_widde/outgroup

#### convert to binary plink files in forward format with plink and snpchimp
python /path/to/snpchimp/source_codes/iConvert/iConvert.py -p 01_data-sources/02_downloaded/06_widde/outgroup/convert.param
plink --cow --ped 01_data-sources/02_downloaded/06_widde/outgroup/cattle__775585variants__2individuals_updated.ped --map 01_data-sources/02_downloaded/06_widde/outgroup/cattle__775585variants__2individuals.map --make-bed --out 01_data-sources/02_downloaded/06_widde/outgroup/cattle__775585variants__2individuals_updated

#### update ids using edited fam files with plink
##### verdugo
plink --cow --bfile 01_data-sources/02_downloaded/04_verdugo/v_snp_data_3 --update-ids 01_data-sources/02_downloaded/04_verdugo/outgroup/outgroup-verdugo-id.txt --make-bed --out 01_data-sources/02_downloaded/04_verdugo/outgroup/outgroup-verdugo

##### widde
plink --cow --bfile 01_data-sources/02_downloaded/06_widde/outgroup/cattle__775585variants__2individuals_updated --update-ids 01_data-sources/02_downloaded/06_widde/outgroup/outgroup-widde-ids.txt --make-bed --out 01_data-sources/02_downloaded/06_widde/outgroup/outgroup-widde

### merge files with plink
plink --cow --bfile 01_data-sources/02_downloaded/04_verdugo/outgroup/outgroup-verdugo --bmerge 01_data-sources/02_downloaded/06_widde/outgroup/outgroup-widde --make-bed --out 01_data-sources/03_merged/outgroup/outgroup-merged

### update genome assembly using edited file with plink
plink --cow --bfile 01_data-sources/03_merged/outgroup/outgroup-merged --update-chr 01_data-sources/02_downloaded/09_schnabel/9913_ARS1.2_777962_HD_marker_name_180910-edit.map 1 2 --update-map 01_data-sources/02_downloaded/09_schnabel/9913_ARS1.2_777962_HD_marker_name_180910-edit.map 4 2 --make-bed --out 02_update-genome-assembly/outgroup/outgroup-ars

### write list of snps from hd and ld data
plink --cow --bfile 03_filter-data/05_filter-snps/hd --write-snplist --out 03_filter-data/06_outgroup/hd
plink --cow --bfile 03_filter-data/05_filter-snps/ld --write-snplist --out 03_filter-data/06_outgroup/ld

### extract hd and ld snps from outgroup
plink --cow --bfile 02_update-genome-assembly/outgroup/outgroup-ars --extract 03_filter-data/06_outgroup/hd.snplist --make-bed --out 03_filter-data/06_outgroup/outgroup-hd
plink --cow --bfile 02_update-genome-assembly/outgroup/outgroup-ars --extract 03_filter-data/06_outgroup/ld.snplist --make-bed --out 03_filter-data/06_outgroup/outgroup-ld

### merge outgroup with rest of data
plink --cow --bfile 03_filter-data/05_filter-snps/hd --bmerge 03_filter-data/06_outgroup/outgroup-hd --make-bed --out 03_filter-data/06_outgroup/outgroup-hd-merged
plink --cow --bfile 03_filter-data/05_filter-snps/ld --bmerge 03_filter-data/06_outgroup/outgroup-ld --make-bed --out 03_filter-data/06_outgroup/outgroup-ld-merged

### filter snps by call rate and minor allele frequency
plink --cow --bfile 03_filter-data/06_outgroup/outgroup-hd-merged --geno 0.05 --maf 0.05 --make-bed --out 03_filter-data/06_outgroup/outgroup-hd-merged-filter
plink --cow --bfile 03_filter-data/06_outgroup/outgroup-ld-merged --geno 0.05 --maf 0.05 --make-bed --out 03_filter-data/06_outgroup/outgroup-ld-merged-filter

### write gzipped allele frequency cluster files
plink --cow --bfile 03_filter-data/06_outgroup/outgroup-hd-merged-filter --freq gz --family --out 06_phylogenetic-analysis/01_input/hd
plink --cow --bfile 03_filter-data/06_outgroup/outgroup-ld-merged-filter --freq gz --family --out 06_phylogenetic-analysis/01_input/ld

### download plink2treemix python script
wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/treemix/plink2treemix.py -O 06_phylogenetic-analysis/01_input/plink2treemix.py

### use plink2treemix python script to convert files
cd 06_phylogenetic-analysis/01_input
python plink2treemix.py hd.frq.strat.gz hd.gz
python plink2treemix.py ld.frq.strat.gz ld.gz
cd ../../..

### run treemix analysis
nohup bash -c 'for m in {1..16}; do for k in {100..1000..100}; do /path/to/treemix/bin/treemix -i 06_phylogenetic-analysis/01_input/hd.gz -root GAUR -global -m $m -k $k -seed $k -o 06_phylogenetic-analysis/02_output/01_migration-edges/01_hd/hd.$m.$k; done; done' > 06_phylogenetic-analysis/02_output/01_migration-edges/01_hd/log.txt &
nohup bash -c 'for m in {1..16}; do for k in {100..1000..100}; do /path/to/treemix/bin/treemix -i 06_phylogenetic-analysis/01_input/ld.gz -root GAUR -global -m $m -k $k -seed $k -o 06_phylogenetic-analysis/02_output/01_migration-edges/02_ld/ld.$m.$k; done; done' > 06_phylogenetic-analysis/02_output/01_migration-edges/02_ld/log.txt &

### run r script to run optm and plot migration edges
Rscript 06_phylogenetic-analysis/03_figures/migration-edges-plots.R

### download bite package
mkdir /path/to/bite
wget https://github.com/marcomilanesi/BITE/archive/refs/heads/master.zip -O /path/to/bite/master.zip
unzip /path/to/bite/master.zip -d /path/to/bite

### run r script to generate treemix bootstrap script with bite package
Rscript 06_phylogenetic-analysis/01_input/bite-bootstrap.R

### edit treemix bootstrap script from bite package to include location of treemix and global flag
sed 's;\ttreemix;\t\/path\/to\/treemix\/bin\/treemix -global;g' 06_phylogenetic-analysis/01_input/treemix_scripts/Treemix_bootstrap.sh > 06_phylogenetic-analysis/01_input/Treemix_bootstrap-edit.sh
sed 's;parallel;\/path\/to\/parallel\/bin\/parallel;g' 06_phylogenetic-analysis/01_input/Treemix_bootstrap-edit.sh > 06_phylogenetic-analysis/01_input/bite-bootstrap.sh

### change permissions to make treemix bootstrap script from bite package executable
chmod 777 06_phylogenetic-analysis/01_input/bite-bootstrap.sh

### copy reuired files to folders
cp 06_phylogenetic-analysis/01_input/hd.gz 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/01_hd/00
cp 06_phylogenetic-analysis/01_input/hd.gz 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/01_hd/03
cp 06_phylogenetic-analysis/01_input/hd.gz 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/01_hd/12
cp 06_phylogenetic-analysis/01_input/ld.gz 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/02_ld/00
cp 06_phylogenetic-analysis/01_input/ld.gz 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/02_ld/03
cp 06_phylogenetic-analysis/01_input/ld.gz 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/02_ld/11
cp 06_phylogenetic-analysis/01_input/bite-bootstrap.sh 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/01_hd/00
cp 06_phylogenetic-analysis/01_input/bite-bootstrap.sh 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/01_hd/03
cp 06_phylogenetic-analysis/01_input/bite-bootstrap.sh 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/01_hd/12
cp 06_phylogenetic-analysis/01_input/bite-bootstrap.sh 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/02_ld/00
cp 06_phylogenetic-analysis/01_input/bite-bootstrap.sh 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/02_ld/03
cp 06_phylogenetic-analysis/01_input/bite-bootstrap.sh 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/02_ld/11

### run bite treemix bootstrap script
cd 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/01_hd/00
nohup bash -c './bite-bootstrap.sh hd.gz 0 5 500 GAUR 100 /path/to/phylip-3.697/exe/consense hd' > log.txt &
cd ..
cd 03
nohup bash -c './bite-bootstrap.sh hd.gz 3 5 500 GAUR 100 /path/to/phylip-3.697/exe/consense hd' > log.txt &
cd ..
cd 12
nohup bash -c './bite-bootstrap.sh hd.gz 12 15 500 GAUR 100 /path/to/phylip-3.697/exe/consense hd' > log.txt &
cd ../..
cd 02_ld/00
nohup bash -c './bite-bootstrap.sh ld.gz 0 5 500 GAUR 100 /path/to/phylip-3.697/exe/consense ld' > log.txt &
cd ..
cd 03
nohup bash -c './bite-bootstrap.sh ld.gz 3 5 500 GAUR 100 /path/to/phylip-3.697/exe/consense ld' > log.txt &
cd ..
cd 11
nohup bash -c './bite-bootstrap.sh ld.gz 11 15 500 GAUR 100 /path/to/phylip-3.697/exe/consense ld' > log.txt &
cd ../../../../..
rm tmp_NoMatch_*.txt

### unzip files
gunzip -c 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/01_hd/03/hd.treeout.gz > 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/01_hd/03/hd.treeout
gunzip -c 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/01_hd/12/hd.treeout.gz > 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/01_hd/12/hd.treeout
gunzip -c 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/02_ld/03/ld.treeout.gz > 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/02_ld/03/ld.treeout
gunzip -c 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/02_ld/11/ld.treeout.gz > 06_phylogenetic-analysis/02_output/02_bootstrap-replicates/02_ld/11/ld.treeout

### r script to generate figures
Rscript 06_phylogenetic-analysis/03_figures/phylogenetic-plots.R