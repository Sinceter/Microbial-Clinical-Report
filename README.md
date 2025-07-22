# Microbial-Clinical-Report


### Usage
```
Rscript mNGS_Kraken-based.R \
  -f path/to/Kraken-style/reports/folder \
  -c path/to/clinial/table \
  --col_header_NR="Nanopore run" \
  --col_header_Bc="Barcode" \
  --col_header_Cl="Clinical culture result" \
  --read_perct_cutoff="0,10,20,30,40,50" \
  --time="12" <tag labeling > \
  --threads 12 <the threads more, the speed faster>
```

```
# example
Rscript mNGS_Kraken-based.R -f ../KB_BCB_mono160_20250722 -c ../MASTER_BCB_mono160_20250722.xlsx --col_header_NR="Nanopore run" --col_header_Bc="Barcode" --col_header_Cl="Clinical culture result" --read_perct_cutoff="0,10,20,30,40,50" --time="12"  --threads 12
```
