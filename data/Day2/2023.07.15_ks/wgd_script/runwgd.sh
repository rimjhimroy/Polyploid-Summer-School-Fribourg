source PATH=src:$PATH

src/bin/wgd dmd -I 3 data/Day2/2023.07.15_ks/wgd_script/Elaeis.fa -o 01.wgd_dmd
src/bin/wgd ksd 01.wgd_dmd/Elaeis.fa.mcl data/Day2/2023.07.15_ks/wgd_script/Elaeis.fa -o 02.wgd_ksd -n 8
src/bin/wgd syn -f mRNA -a ID -ks 02.wgd_ksd/Elaeis.fa.ks.tsv data/Day2/2023.07.15_ks/wgd_script/Elaeis.gff 01.wgd_dmd/Elaeis.fa.mcl -o 03.wgd_syn
src/bin/wgd mix -ni 100 --method bgmm -n 1 5 02.wgd_ksd/Elaeis.fa.ks.tsv -o 04.wgd_mix

