# local ancestry analysis
## local ancestry analysis of hd and ld snp data with mosaic and elai
### mosaic
#### split data into chromosomes with plink
nohup bash -c 'for i in {1..29}; do plink --cow --bfile 03_filter-data/05_filter-snps/hd --chr $i --make-bed --out 07_local-ancestry-analysis/01_mosaic/01_input/01_hd/hd-$i; done' > 07_local-ancestry-analysis/01_mosaic/01_input/01_hd/log.txt &
nohup bash -c 'for i in {1..29}; do plink --cow --bfile 03_filter-data/05_filter-snps/ld --chr $i --make-bed --out 07_local-ancestry-analysis/01_mosaic/01_input/02_ld/ld-$i; done' > 07_local-ancestry-analysis/01_mosaic/01_input/02_ld/log.txt &

#### phase chromosomes with shapeit
nohup bash -c 'for i in {1..29}; do /path/to/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit --input-bed 07_local-ancestry-analysis/01_mosaic/01_input/01_hd/hd-$i.bed 07_local-ancestry-analysis/01_mosaic/01_input/01_hd/hd-$i.bim 07_local-ancestry-analysis/01_mosaic/01_input/01_hd/hd-$i.fam --effective-size 400 --thread 5 --seed 123456789 --force --output-max 07_local-ancestry-analysis/01_mosaic/01_input/01_hd/hd-$i.phased.haps 07_local-ancestry-analysis/01_mosaic/01_input/01_hd/hd-$i.phased.sample; done' >> 07_local-ancestry-analysis/01_mosaic/01_input/01_hd/log.txt &
mv shapeit* 07_local-ancestry-analysis/01_mosaic/01_input/01_hd/log/
nohup bash -c 'for i in {1..29}; do /path/to/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit --input-bed 07_local-ancestry-analysis/01_mosaic/01_input/02_ld/ld-$i.bed 07_local-ancestry-analysis/01_mosaic/01_input/02_ld/ld-$i.bim 07_local-ancestry-analysis/01_mosaic/01_input/02_ld/ld-$i.fam --effective-size 400 --thread 5 --seed 123456789 --force --output-max 07_local-ancestry-analysis/01_mosaic/01_input/02_ld/ld-$i.phased.haps 07_local-ancestry-analysis/01_mosaic/01_input/02_ld/ld-$i.phased.sample; done' >> 07_local-ancestry-analysis/01_mosaic/01_input/02_ld/log.txt &
mv shapeit* 07_local-ancestry-analysis/01_mosaic/01_input/02_ld/log/

#### download r script and change permissions
wget https://maths.ucd.ie/~mst/MOSAIC/convert_from_haps.R -O 07_local-ancestry-analysis/01_mosaic/01_input/convert_from_haps.R
chmod 777 07_local-ancestry-analysis/01_mosaic/01_input/convert_from_haps.R

#### install required r package
R -e "install.packages('argparser')"

#### copy fam file
cp 03_filter-data/05_filter-snps/ld.fam 07_local-ancestry-analysis/01_mosaic/01_input/02_ld
cp 03_filter-data/05_filter-snps/hd.fam 07_local-ancestry-analysis/01_mosaic/01_input/01_hd

#### convert to mosaic format with r script
nohup bash -c 'for i in {1..29}; do 07_local-ancestry-analysis/01_mosaic/01_input/convert_from_haps.R 07_local-ancestry-analysis/01_mosaic/01_input/01_hd $i hd- .phased.haps hd.fam 07_local-ancestry-analysis/01_mosaic/01_input/01_hd; done' >> 07_local-ancestry-analysis/01_mosaic/01_input/01_hd/log.txt &
nohup bash -c 'for i in {1..29}; do 07_local-ancestry-analysis/01_mosaic/01_input/convert_from_haps.R 07_local-ancestry-analysis/01_mosaic/01_input/02_ld $i ld- .phased.haps ld.fam 07_local-ancestry-analysis/01_mosaic/01_input/02_ld; done' >> 07_local-ancestry-analysis/02_mosaic/01_input/02_ld/log.txt &

#### recombination rates adapted from https://github.com/jamesandyward/MOSAIC
Rscript 07_local-ancestry-analysis/01_mosaic/01_input/recombination-rates.R

#### run mosaic
nohup bash -c 'Rscript 07_local-ancestry-analysis/01_mosaic/02_output/mosaic-analysis.R' > /dev/null 2>&1 &

#### copy snpfiles to plot results
nohup bash -c 'for b in 01_ALEN 02_ANKO 03_BORA 04_BORG 05_CHIA 06_EASZ 07_KARA 08_KETE 09_MARC 10_MARE 11_NDAM 12_NGAN 13_ROMA 14_SHEK 15_SOMB; do cp 07_local-ancestry-analysis/01_mosaic/01_input/01_hd/snpfile.* 07_local-ancestry-analysis/01_mosaic/02_output/01_hd/$b; done' > /dev/null 2>&1 &
nohup bash -c 'for b in 01_ALEN 02_ANKO 03_BORA 04_BORG 05_CHIA 06_EASZ 07_KARA 08_KETE 09_MARC 10_MARE 11_NDAM 12_NGAN 13_ROMA 14_SHEK 15_SOMB; do cp 07_local-ancestry-analysis/01_mosaic/01_input/02_ld/snpfile.* 07_local-ancestry-analysis/01_mosaic/02_output/02_ld/$b; done' > /dev/null 2>&1 &

### elai
#### generate bimbam format files for each breed and chromosome using plink
nohup bash -c 'for b in ALEN ANGU ANKO BORA BORG CHIA EASZ GIR HOLS JERS KARA KETE LAGU MARC MARE MUTU NDAG NDAM NELO NGAN ROMA SHEK SOMB THAR; do for i in {1..29}; do plink --cow --bfile 03_filter-data/05_filter-snps/hd --family --keep-cluster-names $b --chr $i --recode-bimbam --out 07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-$b-$i; done; done' > 07_local-ancestry-analysis/02_elai/01_input/01_hd/log.txt &
nohup bash -c 'for b in ALEN ANGU ANKO BORA BORG CHIA EASZ GIR HOLS JERS KARA KETE LAGU MARC MARE MUTU NDAG NDAM NELO NGAN ROMA SHEK SOMB THAR; do for i in {1..29}; do plink --cow --bfile 03_filter-data/05_filter-snps/ld --family --keep-cluster-names $b --chr $i --recode-bimbam --out 07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-$b-$i; done; done' > 07_local-ancestry-analysis/02_elai/01_input/02_ld/log.txt &

#### run elai analysis
cd 07_local-ancestry-analysis/02_elai/02_output/01_hd/01_ALEN
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-LAGU-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-GIR-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-ALEN-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-ALEN-$i.recode.pos.txt -s 30 -o hd-ALEN-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 02_ANKO
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-LAGU-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-ANKO-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-ANKO-$i.recode.pos.txt -s 30 -o hd-ANKO-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 03_BORA
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-LAGU-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-BORA-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-BORA-$i.recode.pos.txt -s 30 -o hd-BORA-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 04_BORG
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-LAGU-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-BORG-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-BORG-$i.recode.pos.txt -s 30 -o hd-BORG-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 05_CHIA
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-HOLS-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-NDAG-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-GIR-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-CHIA-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-CHIA-$i.recode.pos.txt -s 30 -o hd-CHIA-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 06_EASZ
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-LAGU-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-EASZ-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-EASZ-$i.recode.pos.txt -s 30 -o hd-EASZ-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 07_KARA
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-LAGU-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-KARA-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-KARA-$i.recode.pos.txt -s 30 -o hd-KARA-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 08_KETE
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-LAGU-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-KETE-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-KETE-$i.recode.pos.txt -s 30 -o hd-KETE-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 09_MARC
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-NDAG-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-GIR-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-MARC-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-MARC-$i.recode.pos.txt -s 30 -o hd-MARC-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 10_MARE
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-NDAG-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-GIR-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-MARE-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-MARE-$i.recode.pos.txt -s 30 -o hd-MARE-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 11_NDAM
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-NDAG-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-NDAM-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-NDAM-$i.recode.pos.txt -s 30 -o hd-NDAM-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 12_NGAN
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-HOLS-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-LAGU-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-NGAN-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-NGAN-$i.recode.pos.txt -s 30 -o hd-NGAN-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 13_ROMA
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-LAGU-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-GIR-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-ROMA-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-ROMA-$i.recode.pos.txt -s 30 -o hd-ROMA-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 14_SHEK
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-LAGU-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-SHEK-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-SHEK-$i.recode.pos.txt -s 30 -o hd-SHEK-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 15_SOMB
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-LAGU-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-SOMB-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/01_hd/hd-SOMB-$i.recode.pos.txt -s 30 -o hd-SOMB-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ../..
cd 02_ld/01_ALEN
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-LAGU-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-GIR-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-ALEN-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-ALEN-$i.recode.pos.txt -s 30 -o ld-ALEN-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 02_ANKO
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-LAGU-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-ANKO-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-ANKO-$i.recode.pos.txt -s 30 -o ld-ANKO-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 03_BORA
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-LAGU-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-BORA-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-BORA-$i.recode.pos.txt -s 30 -o ld-BORA-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 04_BORG
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-LAGU-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-BORG-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-BORG-$i.recode.pos.txt -s 30 -o ld-BORG-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 05_CHIA
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-HOLS-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-NDAG-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-GIR-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-CHIA-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-CHIA-$i.recode.pos.txt -s 30 -o ld-CHIA-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 06_EASZ
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-LAGU-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-EASZ-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-EASZ-$i.recode.pos.txt -s 30 -o ld-EASZ-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 07_KARA
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-LAGU-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-KARA-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-KARA-$i.recode.pos.txt -s 30 -o ld-KARA-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 08_KETE
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-LAGU-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-KETE-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-KETE-$i.recode.pos.txt -s 30 -o ld-KETE-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 09_MARC
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-NDAG-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-GIR-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-MARC-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-MARC-$i.recode.pos.txt -s 30 -o ld-MARC-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 10_MARE
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-NDAG-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-GIR-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-MARE-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-MARE-$i.recode.pos.txt -s 30 -o ld-MARE-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 11_NDAM
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-NDAG-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-NDAM-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-NDAM-$i.recode.pos.txt -s 30 -o ld-NDAM-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 12_NGAN
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-HOLS-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-LAGU-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-NGAN-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-NGAN-$i.recode.pos.txt -s 30 -o ld-NGAN-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 13_ROMA
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-LAGU-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-GIR-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-ROMA-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-ROMA-$i.recode.pos.txt -s 30 -o ld-ROMA-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 14_SHEK
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-LAGU-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-SHEK-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-SHEK-$i.recode.pos.txt -s 30 -o ld-SHEK-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ..
cd 15_SOMB
nohup bash -c 'for i in {1..29}; do /path/to/elai/elai-lin -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-ANGU-$i.recode.geno.txt -p 10 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-GIR-$i.recode.geno.txt -p 11 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-LAGU-$i.recode.geno.txt -p 12 -g /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-SOMB-$i.recode.geno.txt -p 1 -pos /home/workspace/gmchugo/local-ancestry/07_local-ancestry-analysis/02_elai/01_input/02_ld/ld-SOMB-$i.recode.pos.txt -s 30 -o ld-SOMB-$i -C 3 -c 15 -mg 200 -R 123456789; done' > log.txt &
cd ../../../../..

### extract results
nohup bash -c 'Rscript 07_local-ancestry-analysis/local-ancestry-results.R' > /dev/null 2>&1 &

### generate local ancestry plots
nohup bash -c 'Rscript 07_local-ancestry-analysis/local-ancestry-plots.R' > 07_local-ancestry-analysis/log.txt &

### combine local ancestry plots
nohup bash -c 'Rscript 07_local-ancestry-analysis/03_comparison/combine-local-ancestry-plots.R' > /dev/null 2>&1 &

### generate correlation plots
nohup bash -c 'Rscript 07_local-ancestry-analysis/03_comparison/correlation-plots.R' > /dev/null 2>&1 &