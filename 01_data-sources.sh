# data sources
## download and prepare data
### download files with wget or web interfaces
#### bahbahani
wget https://datadryad.org/stash/downloads/download_resource/1159 -O 01_data-sources/02_downloaded/01_bahbahani/doi_10.5061_dryad.38jp6__v1.zip

#### barbato
wget https://data.mendeley.com/public-files/datasets/znr5dy4x29/files/1d646dcb-d3ac-4264-a95f-e619ed6f71fa/file_downloaded -O 01_data-sources/02_downloaded/02_barbato/Barbato2020_Dataset.zip

#### upadhyay
wget https://datadryad.org/stash/downloads/download_resource/1793 -O 01_data-sources/02_downloaded/03_upadhyay/doi_10.5061_dryad.f2d1q__v1.zip

#### verdugo
wget https://osf.io/download/r6ntf/ -O 01_data-sources/02_downloaded/04_verdugo/v_snp_data.tar.gz

#### ward
wget https://data.mendeley.com/public-files/datasets/yt3tgpt48d/files/9e6e53d6-8789-44fa-a27d-e23e0b09ecdd/file_downloaded -O 01_data-sources/02_downloaded/05_ward/breed_info.txt
wget https://data.mendeley.com/public-files/datasets/yt3tgpt48d/files/6a99c0de-eb6f-43d8-8911-6dbee105b6b8/file_downloaded -O 01_data-sources/02_downloaded/05_ward/README.txt
wget https://data.mendeley.com/public-files/datasets/yt3tgpt48d/files/86d6bc5e-14fb-44ba-a4e4-f28031af73ca/file_downloaded -O 01_data-sources/02_downloaded/05_ward/source_c.map
wget https://data.mendeley.com/public-files/datasets/yt3tgpt48d/files/fcffbf44-ce5d-49ee-a297-5122d103aa36/file_downloaded -O 01_data-sources/02_downloaded/05_ward/source_c.ped
wget https://data.mendeley.com/public-files/datasets/yt3tgpt48d/files/429ae9bd-0f97-421e-9914-badde24deb59/file_downloaded -O 01_data-sources/02_downloaded/05_ward/source_d.map
wget https://data.mendeley.com/public-files/datasets/yt3tgpt48d/files/1002bbf0-a465-449f-bc1f-59edefcef0ab/file_downloaded -O 01_data-sources/02_downloaded/05_ward/source_d.ped

#### widde
##### select required data from http://widde.toulouse.inra.fr/widde/

#### wragg
wget https://datashare.ed.ac.uk/download/DS_10283_3835.zip -O 01_data-sources/02_downloaded/07_wragg/DS_10283_3835.zip

#### snpchimp
##### select required data from https://webserver.ibba.cnr.it/SNPchimp/index.php/download/download-cow-data

#### schnabel
wget https://www.animalgenome.org/repository/cattle/UMC_bovine_coordinates/UMC_marker_names_180910.zip -O 01_data-sources/02_downloaded/09_schnabel/UMC_marker_names_180910.zip

#### ma
wget https://datadryad.org/stash/downloads/download_resource/1517 -O 01_data-sources/02_downloaded/10_ma/doi_10.5061_dryad.q2q84__v1.zip

### unzip downloaded files with unzip, gunzip or tar
#### bahbahani
unzip 01_data-sources/02_downloaded/01_bahbahani/doi_10.5061_dryad.38jp6__v1.zip -d 01_data-sources/02_downloaded/01_bahbahani/doi_10.5061_dryad.38jp6__v1
unzip 01_data-sources/02_downloaded/01_bahbahani/doi_10.5061_dryad.38jp6__v1/HD_SNP.zip -d 01_data-sources/02_downloaded/01_bahbahani/doi_10.5061_dryad.38jp6__v1

#### barbato
unzip 01_data-sources/02_downloaded/02_barbato/Barbato2020_Dataset.zip -d 01_data-sources/02_downloaded/02_barbato/Barbato2020_Dataset

#### upadhyay
unzip 01_data-sources/02_downloaded/03_upadhyay/doi_10.5061_dryad.f2d1q__v1.zip -d 01_data-sources/02_downloaded/03_upadhyay/doi_10.5061_dryad.f2d1q__v1
unzip 01_data-sources/02_downloaded/03_upadhyay/doi_10.5061_dryad.f2d1q__v1/PrimitiveCattleGenotypes.zip -d 01_data-sources/02_downloaded/03_upadhyay/doi_10.5061_dryad.f2d1q__v1/PrimitiveCattleGenotypes

#### verdugo
tar -xf 01_data-sources/02_downloaded/04_verdugo/v_snp_data.tar.gz -C 01_data-sources/02_downloaded/04_verdugo

#### widde
unzip 01_data-sources/02_downloaded/06_widde/cattle_*_PLINK.zip -d 01_data-sources/02_downloaded/06_widde

#### wragg
unzip 01_data-sources/02_downloaded/07_wragg/DS_10283_3835.zip -d 01_data-sources/02_downloaded/07_wragg/DS_10283_3835

#### snpchimp
gunzip 01_data-sources/02_downloaded/08_snpchimp/SNPchimp_result_*.csv.gz

#### schnabel
unzip 01_data-sources/02_downloaded/09_schnabel/UMC_marker_names_180910.zip -d 01_data-sources/02_downloaded/09_schnabel

#### ma
unzip 01_data-sources/02_downloaded/10_ma/doi_10.5061_dryad.q2q84__v1.zip -d 01_data-sources/02_downloaded/10_ma/doi_10.5061_dryad.q2q84__v1

### convert snp data to binary plink files in forward format with plink and snpchimp
#### bora-ndam
cd 01_data-sources/01_generated/01_bora-ndam
python /path/to/snpchimp/source_codes/PEDDA_ROW/pedda_row.py
plink --cow --file Univ_College_Dublin_BOV770V01_20221212 --make-bed --out Univ_College_Dublin_BOV770V01_20221212
cd ../../..

#### somb
python /path/to/snpchimp/source_codes/iConvert/iConvert.py -p 01_data-sources/01_generated/02_somb/convert.param
plink --cow --ped 01_data-sources/01_generated/02_somb/we_ucd_27042020/PLINK_260420_0702/PLINK_260420_0702/we_ucd_27042020_updated.ped --map 01_data-sources/01_generated/02_somb/we_ucd_27042020/PLINK_260420_0702/PLINK_260420_0702/we_ucd_27042020.map --make-bed --out 01_data-sources/01_generated/02_somb/we_ucd_27042020_updated

#### bahbahani
python /path/to/snpchimp/source_codes/iConvert/iConvert.py -p 01_data-sources/02_downloaded/01_bahbahani/convert.param
plink --cow --ped 01_data-sources/02_downloaded/01_bahbahani/doi_10.5061_dryad.38jp6__v1/HD_SNP_updated.ped --map 01_data-sources/02_downloaded/01_bahbahani/doi_10.5061_dryad.38jp6__v1/HD_SNP.map --make-bed --out 01_data-sources/02_downloaded/01_bahbahani/HD_SNP_updated

#### barbato
plink --cow --bfile 01_data-sources/02_downloaded/02_barbato/Barbato2020_Dataset/Barbato2020_Dataset/Barbato2020_data --recode --out 01_data-sources/02_downloaded/02_barbato/Barbato2020_data
python /path/to/snpchimp/source_codes/iConvert/iConvert.py -p 01_data-sources/02_downloaded/02_barbato/convert.param
plink --cow --ped 01_data-sources/02_downloaded/02_barbato/Barbato2020_data_updated.ped --map 01_data-sources/02_downloaded/02_barbato/Barbato2020_data.map --make-bed --out 01_data-sources/02_downloaded/02_barbato/Barbato2020_data_updated

#### upadhyay
plink --cow --file 01_data-sources/02_downloaded/03_upadhyay/doi_10.5061_dryad.f2d1q__v1/PrimitiveCattleGenotypes/dryad_submission/primitive_cattle --make-bed --out 01_data-sources/02_downloaded/03_upadhyay/primitive_cattle

#### ward
plink --cow --file 01_data-sources/02_downloaded/05_ward/source_c --make-bed --out 01_data-sources/02_downloaded/05_ward/source_c
python /path/to/snpchimp/source_codes/iConvert/iConvert.py -p 01_data-sources/02_downloaded/05_ward/convert.param
plink --cow --ped 01_data-sources/02_downloaded/05_ward/source_d_updated.ped --map 01_data-sources/02_downloaded/05_ward/source_d.map --make-bed --out 01_data-sources/02_downloaded/05_ward/source_d_updated

#### widde
python /path/to/snpchimp/source_codes/iConvert/iConvert.py -p 01_data-sources/02_downloaded/06_widde/convert.param
plink --cow --ped 01_data-sources/02_downloaded/06_widde/cattle__775585variants__263individuals_updated.ped --map 01_data-sources/02_downloaded/06_widde/cattle__775585variants__263individuals.map --make-bed --out 01_data-sources/02_downloaded/06_widde/cattle__775585variants__263individuals_updated

#### wragg
plink --cow --file 01_data-sources/02_downloaded/07_wragg/DS_10283_3835/set1 --make-bed --out 01_data-sources/02_downloaded/07_wragg/set1
python /path/to/snpchimp/source_codes/iConvert/iConvert.py -p 01_data-sources/02_downloaded/07_wragg/convert.param
plink --cow --ped 01_data-sources/02_downloaded/07_wragg/DS_10283_3835/set2_updated.ped --map 01_data-sources/02_downloaded/07_wragg/DS_10283_3835/set2.map --make-bed --out 01_data-sources/02_downloaded/07_wragg/set2_updated

### extract selected samples using edited fam files with plink
#### somb
plink --cow --bfile 01_data-sources/01_generated/02_somb/we_ucd_27042020_updated --keep 01_data-sources/01_generated/02_somb/somb-keep.txt --make-bed --out 01_data-sources/01_generated/02_somb/we_ucd_27042020_updated-keep

#### bahbahani
plink --cow --bfile 01_data-sources/02_downloaded/01_bahbahani/HD_SNP_updated --keep 01_data-sources/02_downloaded/01_bahbahani/bahbahani-keep.txt --make-bed --out 01_data-sources/02_downloaded/01_bahbahani/HD_SNP_updated-keep

#### barbato
plink --cow --bfile 01_data-sources/02_downloaded/02_barbato/Barbato2020_data_updated --keep 01_data-sources/02_downloaded/02_barbato/barbato-keep.txt --make-bed --out 01_data-sources/02_downloaded/02_barbato/Barbato2020_data_updated-keep

#### upadhyay
plink --cow --bfile  01_data-sources/02_downloaded/03_upadhyay/primitive_cattle --keep  01_data-sources/02_downloaded/03_upadhyay/upadhyay-keep.txt --make-bed --out  01_data-sources/02_downloaded/03_upadhyay/primitive_cattle-keep

#### verdugo
plink --cow --bfile  01_data-sources/02_downloaded/04_verdugo/v_snp_data_1 --keep  01_data-sources/02_downloaded/04_verdugo/verdugo-keep.txt --make-bed --out  01_data-sources/02_downloaded/04_verdugo/v_snp_data_1-keep

#### ward
plink --cow --bfile  01_data-sources/02_downloaded/05_ward/source_c --keep  01_data-sources/02_downloaded/05_ward/ward-keep.txt --make-bed --out  01_data-sources/02_downloaded/05_ward/source_c-keep

### update snp data sample ids using edited fam files with plink
#### bora-ndam
plink --cow --bfile 01_data-sources/01_generated/01_bora-ndam/Univ_College_Dublin_BOV770V01_20221212 --update-ids 01_data-sources/01_generated/01_bora-ndam/bora-ndam-ids.txt --make-bed --out 01_data-sources/01_generated/01_bora-ndam/bora-ndam

#### somb
plink --cow --bfile 01_data-sources/01_generated/02_somb/we_ucd_27042020_updated-keep --update-ids 01_data-sources/01_generated/02_somb/somb-ids.txt --make-bed --out 01_data-sources/01_generated/02_somb/somb

#### bahbahani
plink --cow --bfile 01_data-sources/02_downloaded/01_bahbahani/HD_SNP_updated-keep --update-ids 01_data-sources/02_downloaded/01_bahbahani/bahbahani-ids.txt --make-bed --out 01_data-sources/02_downloaded/01_bahbahani/bahbahani

#### barbato
plink --cow --bfile 01_data-sources/02_downloaded/02_barbato/Barbato2020_data_updated-keep --update-ids 01_data-sources/02_downloaded/02_barbato/barbato-ids.txt --make-bed --out 01_data-sources/02_downloaded/02_barbato/Barbato2020_data_updated-keep-ids
plink --cow --bfile 01_data-sources/02_downloaded/02_barbato/Barbato2020_data_updated-keep-ids --update-parents 01_data-sources/02_downloaded/02_barbato/barbato-parent-ids.txt --make-bed --out 01_data-sources/02_downloaded/02_barbato/barbato

#### upadhyay
plink --cow --bfile 01_data-sources/02_downloaded/03_upadhyay/primitive_cattle-keep --update-ids 01_data-sources/02_downloaded/03_upadhyay/upadhyay-ids.txt --make-bed --out 01_data-sources/02_downloaded/03_upadhyay/upadhyay

#### verdugo
plink --cow --bfile 01_data-sources/02_downloaded/04_verdugo/v_snp_data_1-keep --update-ids 01_data-sources/02_downloaded/04_verdugo/verdugo-ids.txt --make-bed --out 01_data-sources/02_downloaded/04_verdugo/verdugo

#### ward
plink --cow --bfile 01_data-sources/02_downloaded/05_ward/source_c-keep --update-ids 01_data-sources/02_downloaded/05_ward/ward-ids-01.txt --make-bed --out 01_data-sources/02_downloaded/05_ward/source_c-keep-ids
plink --cow --bfile 01_data-sources/02_downloaded/05_ward/source_d_updated --update-ids 01_data-sources/02_downloaded/05_ward/ward-ids-02.txt --make-bed --out 01_data-sources/02_downloaded/05_ward/source_d_updated-ids

#### widde
plink --cow --bfile 01_data-sources/02_downloaded/06_widde/cattle__775585variants__263individuals_updated --update-ids 01_data-sources/02_downloaded/06_widde/widde-ids.txt --make-bed --out 01_data-sources/02_downloaded/06_widde/widde

#### wragg
plink --cow --bfile 01_data-sources/02_downloaded/07_wragg/set1 --update-ids 01_data-sources/02_downloaded/07_wragg/wragg-ids-01.txt --make-bed --out 01_data-sources/02_downloaded/07_wragg/set1-ids
plink --cow --bfile 01_data-sources/02_downloaded/07_wragg/set2_updated --update-ids 01_data-sources/02_downloaded/07_wragg/wragg-ids-02.txt --make-bed --out 01_data-sources/02_downloaded/07_wragg/set2_updated-ids

### merge files with plink
#### generated
plink --cow --bfile 01_data-sources/01_generated/01_bora-ndam/bora-ndam --bmerge 01_data-sources/01_generated/02_somb/somb --make-bed --out 01_data-sources/01_generated/generated

#### ward
plink --cow --bfile 01_data-sources/02_downloaded/05_ward/source_c-keep-ids --bmerge  01_data-sources/02_downloaded/05_ward/source_d_updated-ids --make-bed --out 01_data-sources/02_downloaded/05_ward/ward

#### wragg
plink --cow --bfile 01_data-sources/02_downloaded/07_wragg/set1-ids --bmerge  01_data-sources/02_downloaded/07_wragg/set2_updated-ids --make-bed --out 01_data-sources/02_downloaded/07_wragg/wragg

#### merge all
plink --cow --merge-list 01_data-sources/03_merged/merge-list.txt --make-bed --out 01_data-sources/03_merged/merged