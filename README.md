# 🧬 AML Single-Cell Landscape

**Single-cell RNA-seq analysis of Acute Myeloid Leukemia — cell type annotation, differentiation trajectories, and malignant cell classification.**

This project is the second part of a portfolio focused on computational oncogenomics applied to AML.
It extends the bulk RNA / mutation analysis from [`aml-tcga-genomics`](https://github.com/YOUR_USERNAME/aml-tcga-genomics) by going to single-cell resolution.

---

## 🔬 Biological question

> *What is the cellular composition of AML tumors, at which differentiation stage are malignant cells blocked, and can we classify them computationally?*

AML is characterized by a differentiation block: hematopoietic progenitors fail to mature and accumulate as blasts. scRNA-seq allows us to map this at single-cell resolution — something bulk RNA-seq cannot do.

---

## 📂 Dataset

**GSE116256 — van Galen et al., Cell 2019**
*"Single-Cell RNA-Sequencing Reveals a Hierarchy of Leukemic Cell States and Their Relationship to Normal Hematopoiesis"*

- 38,410 single cells from **16 AML patients** + 5 healthy donors
- Pre-processed `.h5ad` files available directly on GEO
- Gold standard dataset for AML single-cell research

**Download instructions:** see [`data/README.md`](data/README.md)

---

## 🗂️ Project structure

```
aml-scrna-landscape/
├── README.md
├── environment.yml              # Conda environment
├── data/
│   ├── README.md               # Download instructions (data not committed)
│   └── raw/                    # .h5ad files go here (gitignored)
├── notebooks/
│   ├── 01_data_loading_qc.ipynb
│   ├── 02_preprocessing_normalization.ipynb
│   ├── 03_clustering_annotation.ipynb
│   ├── 04_trajectory_analysis.ipynb
│   └── 05_ml_malignant_classification.ipynb
├── src/
│   └── utils.py                # Reusable functions
└── figures/                    # Exported plots (PNG/SVG)
```

---

## ⚙️ Analysis pipeline

| Notebook | Topic | Key tools |
|---|---|---|
| 01 | Data loading & Quality Control | `scanpy`, `anndata` |
| 02 | Preprocessing & Normalization | `scanpy`, `harmonypy` |
| 03 | Clustering & Cell Type Annotation | `scanpy`, `leiden` |
| 04 | Trajectory / Pseudotime Analysis | `scanpy` PAGA + DPT |
| 05 | ML — Malignant Cell Classification | `scikit-learn`, `xgboost` |

---

## 📊 Key results

*(To be updated as analysis progresses)*

- [ ] UMAP showing AML cell type landscape
- [ ] Marker gene heatmap per cluster
- [ ] PAGA differentiation graph
- [ ] Pseudotime trajectory
- [ ] Malignant vs normal classifier (AUC, feature importance)

---

## 🛠️ Environment setup

This project shares the virtual environment `ven-aml-tcga` with [`aml-tcga-genomics`](https://github.com/YOUR_USERNAME/aml-tcga-genomics).

```bash
# Clone the repo (inside "Projets Github/")
git clone https://github.com/YOUR_USERNAME/aml-scrna-landscape.git

# Activate the shared environment
source ../ven-aml-tcga/bin/activate      # Mac/Linux
# ..\ven-aml-tcga\Scripts\activate       # Windows

# Install the additional scRNA packages (first time only)
pip install scanpy anndata harmonypy leidenalg igraph xgboost shap GEOparse

# Launch Jupyter
cd aml-scrna-landscape
jupyter lab
```

**Directory structure:**
```
Projets Github/
├── ven-aml-tcga/              ← shared virtual environment
├── aml-tcga-genomics/         ← project 1 (bulk RNA / mutations)
└── aml-scrna-landscape/       ← project 2 (this repo)
```

---

## 📚 References

- van Galen P. et al. *Single-Cell RNA-Sequencing Reveals a Hierarchy of Leukemic Cell States and Their Relationship to Normal Hematopoiesis.* Cell, 2019. https://doi.org/10.1016/j.cell.2019.01.031
- Wolf F.A. et al. *SCANPY: large-scale single-cell gene expression data analysis.* Genome Biology, 2018.
- Bergen V. et al. *Generalizing RNA velocity to transient cell states through dynamical modeling.* Nature Biotechnology, 2020.

---

## 👩‍💻 About

Data Scientist specializing in cancer genomics, with a focus on Acute Myeloid Leukemia (AML). I analyze tumor genomic data — mutational profiling, somatic variant analysis, survival modeling — to contribute to a better understanding of cancer biology and treatment response.

Currently working at Airbus Defence & Space, I am in parallel building expertise in computational oncogenomics on TCGA data, with the goal of joining a specialized team in cancer genomics.

I am actively looking for opportunities within oncology research teams — particularly at **l'Oncopole de Toulouse** and affiliated INSERM/CNRS units.

📧 Reach me via [LinkedIn](https://www.linkedin.com/in/jessicalalanne/) | [GitHub](https://github.com/jesslalanne)

---
