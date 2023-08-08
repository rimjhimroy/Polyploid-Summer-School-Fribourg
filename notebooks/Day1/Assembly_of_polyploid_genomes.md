# Polyploid Genome Assembly

## Dataset for Day1
```bash
wget -O Day1.zip https://drive.switch.ch/index.php/s/HOb9sSQV7ISIkoR/download
unzip Day1.zip

```
## Conda environments
Conda environment: `conda activate assembly`

## KMC
```bash
cd Day1/ds1
mkdir tmp
kmc -k21 -t16 -m64 -ci1 -cs10000 ds1_hifi.fastq.gz kmcdb tmp
kmc_tools transform kmcdb histogram kmcdb_k21.hist -cx10000
L=$(smudgeplot.py cutoff kmcdb_k21.hist L)
U=$(smudgeplot.py cutoff kmcdb_k21.hist U)
echo $L
echo $U
kmc_tools transform kmcdb -ci"$L" -cx"$U" dump -s kmcdb_L"$L"_U"$U".dump
smudgeplot.py hetkmers -o kmcdb_L"$L"_U"$U" < kmcdb_L"$L"_U"$U".dump

smudgeplot.py plot  kmcdb_L"$L"_U"$U"_coverages.tsv


hifiasm  -o ds1.asm ds1_hifi.fastq.gz # need more memory
```