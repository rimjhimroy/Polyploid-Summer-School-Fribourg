source PATH=/home/rstudio/work/polyploid-summer-school-2023/src:$PATH

/home/rstudio/work/polyploid-summer-school-2023/src/bin/wgd dmd -I 3 Elaeis.fa -o 01.wgd_dmd
/home/rstudio/work/polyploid-summer-school-2023/src/bin/wgd ksd 01.wgd_dmd/Elaeis.fa.mcl Elaeis.fa -o 02.wgd_ksd -n 8
/home/rstudio/work/polyploid-summer-school-2023/src/bin/wgd syn -f mRNA -a ID -ks 02.wgd_ksd/Elaeis.fa.ks.tsv Elaeis.gff 01.wgd_dmd/Elaeis.fa.mcl -o 03.wgd_syn
/home/rstudio/work/polyploid-summer-school-2023/src/bin/wgd mix -ni 100 --method bgmm -n 1 5 02.wgd_ksd/Elaeis.fa.ks.tsv -o 04.wgd_mix

