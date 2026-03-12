# Data

Raw data files are **not committed** to this repository (genomic data, large files).

## Dataset: GSE116256 — van Galen et al., Cell 2019

### Option 1 — Manual download from GEO (recommended)

1. Go to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116256
2. Scroll to **Supplementary files**
3. Download the `.h5ad` files (pre-processed by the authors):
   - `GSE116256_RAW.tar` contains per-sample `.h5ad` files
4. Extract and place files in `data/raw/`

### Option 2 — Programmatic download

```python
import GEOparse

gse = GEOparse.get_GEO(geo="GSE116256", destdir="data/raw/")
```

### Option 3 — Direct curl

```bash
# Full supplementary archive (~2GB)
cd data/raw/
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116256/suppl/GSE116256_RAW.tar
tar -xf GSE116256_RAW.tar
```

---

## Expected files after download

```
data/raw/
├── GSM3587923_AML1012-D0.dem.txt.gz      # AML patient samples
├── GSM3587924_AML1012-D0.anno.txt.gz
├── ...                                     # 16 AML patients × 2 files
├── GSM3587999_BM1.dem.txt.gz             # Healthy donor bone marrow
└── GSM3588000_BM1.anno.txt.gz
```

Each patient has two files:
- `.dem.txt.gz` — digital expression matrix (genes × cells)
- `.anno.txt.gz` — cell annotations (cell type labels from van Galen)

---

## Paper reference

van Galen P, Hovestadt V, Wadsworth MH II, et al.  
*Single-Cell RNA-Sequencing Reveals a Hierarchy of Leukemic Cell States and Their Relationship to Normal Hematopoiesis.*  
Cell. 2019;176(6):1265-1281.e24.  
https://doi.org/10.1016/j.cell.2019.01.031
