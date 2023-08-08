# Ks age distribution and genomic dotplots

## Dataset for Day2
```bash
wget -O Day2.zip https://drive.switch.ch/index.php/s/N97Fpkvpu56bY4Q/download
unzip Day2.zip
ls 2023.07.15_ks

```

## Conda environments

Conda environment for wgd: `conda activate wgd`
Conda environment for ksrates: `conda activate ksrates`
Conda environment for i-adhore: `conda activate wgd`

## WGD and I-adhore
```bash
# Run wgd
conda activate wgd
cd wgd_script/
wgd dmd -I 3 Elaeis.fa -o 01.wgd_dmd
wgd ksd 01.wgd_dmd/Elaeis.fa.mcl Elaeis.fa -o 02.wgd_ksd -n 8
wgd syn -f mRNA -a ID -ks 02.wgd_ksd/Elaeis.fa.ks.tsv Elaeis.gff 01.wgd_dmd/Elaeis.fa.mcl -o 03.wgd_syn
wgd mix -ni 100 --method bgmm -n 1 5 02.wgd_ksd/Elaeis.fa.ks.tsv -o 04.wgd_mix


# Run I-adhore
cd iadhore_script/
conda activate wgd
bash ./running_diamond.shell ../data/Species.info.xls . diamond 4
Rscript install_packages.R
Rscript prepare_i-ADHoRe.v2.R -i ./data/Species.info.xls -o . -c run_iadhore.sh

bash run_iadhore.sh

```
## ksrates

```bash
conda activate ksrates
cd ksrates_script


```