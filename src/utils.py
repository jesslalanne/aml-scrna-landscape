"""
utils.py — Reusable functions for AML scRNA-seq analysis pipeline.

Conventions:
    - All variable names and comments in English
    - Functions return AnnData objects or DataFrames, never modify in place silently
    - Each function has a docstring with Parameters / Returns
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
from pathlib import Path
from typing import Optional, List, Dict, Tuple


# ─────────────────────────────────────────────────────────────────────────────
# 1. DATA LOADING
# ─────────────────────────────────────────────────────────────────────────────

def load_van_galen_sample(
    dem_path: str,
    anno_path: str,
    sample_name: str
) -> sc.AnnData:
    """
    Load one AML or healthy donor sample from van Galen 2019 raw files.

    Parameters
    ----------
    dem_path : str
        Path to the .dem.txt.gz file (digital expression matrix, genes x cells).
    anno_path : str
        Path to the .anno.txt.gz file (cell annotations).
    sample_name : str
        Label for this sample (e.g., 'AML1012', 'BM1').

    Returns
    -------
    sc.AnnData
        AnnData object with:
        - .X : count matrix (cells x genes)
        - .obs : cell annotations from van Galen
        - .obs['sample'] : sample label
    """
    # Load expression matrix (genes as rows, cells as columns)
    dem = pd.read_csv(dem_path, sep="\t", index_col=0, compression="gzip")
    # Transpose to cells x genes (standard AnnData convention)
    adata = sc.AnnData(X=dem.T.values, obs=pd.DataFrame(index=dem.columns), var=pd.DataFrame(index=dem.index))

    # Load annotations
    anno = pd.read_csv(anno_path, sep="\t", index_col=0, compression="gzip")
    # Align annotations to cell barcodes
    common_cells = adata.obs.index.intersection(anno.index)
    adata = adata[common_cells].copy()
    adata.obs = anno.loc[common_cells]

    # Tag sample
    adata.obs["sample"] = sample_name
    adata.obs["is_aml"] = not sample_name.startswith("BM")

    print(f"  Loaded {sample_name}: {adata.n_obs} cells, {adata.n_vars} genes")
    return adata


def load_all_samples(data_dir: str) -> sc.AnnData:
    """
    Load and concatenate all samples from van Galen 2019.

    Parameters
    ----------
    data_dir : str
        Path to the folder containing .dem.txt.gz and .anno.txt.gz files.

    Returns
    -------
    sc.AnnData
        Concatenated AnnData with all cells.
    """
    data_path = Path(data_dir)
    dem_files = sorted(data_path.glob("*.dem.txt.gz"))

    adatas = []
    for dem_file in dem_files:
        sample_name = dem_file.stem.replace(".dem.txt", "").split("_")[-1]
        anno_file = dem_file.parent / dem_file.name.replace(".dem.txt.gz", ".anno.txt.gz")

        if not anno_file.exists():
            print(f"  Warning: annotation file not found for {sample_name}, skipping.")
            continue

        adata = load_van_galen_sample(str(dem_file), str(anno_file), sample_name)
        adatas.append(adata)

    combined = sc.concat(adatas, label="sample", keys=[a.obs["sample"].iloc[0] for a in adatas])
    print(f"\nTotal: {combined.n_obs} cells from {len(adatas)} samples")
    return combined


# ─────────────────────────────────────────────────────────────────────────────
# 2. QUALITY CONTROL
# ─────────────────────────────────────────────────────────────────────────────

def compute_qc_metrics(adata: sc.AnnData) -> sc.AnnData:
    """
    Add QC metrics to AnnData: n_genes, total_counts, pct_counts_mt.

    Parameters
    ----------
    adata : sc.AnnData
        Raw AnnData object.

    Returns
    -------
    sc.AnnData
        AnnData with QC columns added to .obs.
    """
    # Flag mitochondrial genes (human: MT- prefix)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True
    )
    return adata


def filter_cells_qc(
    adata: sc.AnnData,
    min_genes: int = 200,
    max_genes: int = 6000,
    max_pct_mt: float = 20.0
) -> sc.AnnData:
    """
    Filter cells based on QC thresholds.

    Parameters
    ----------
    adata : sc.AnnData
    min_genes : int
        Minimum number of genes detected per cell.
    max_genes : int
        Maximum number of genes (removes doublets).
    max_pct_mt : float
        Maximum percentage of mitochondrial counts (removes dying cells).

    Returns
    -------
    sc.AnnData
        Filtered AnnData.
    """
    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs["n_genes_by_counts"] < max_genes].copy()
    adata = adata[adata.obs["pct_counts_mt"] < max_pct_mt].copy()
    n_after = adata.n_obs

    print(f"QC filtering: {n_before} → {n_after} cells retained "
          f"({n_before - n_after} removed, {100*(n_before-n_after)/n_before:.1f}%)")
    return adata


def plot_qc_violin(adata: sc.AnnData, save_path: Optional[str] = None) -> None:
    """
    Plot QC metrics as violin plots.

    Parameters
    ----------
    adata : sc.AnnData
    save_path : str, optional
        If provided, save figure to this path.
    """
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    labels  = ["Genes per cell", "Total counts per cell", "% Mitochondrial counts"]

    for ax, metric, label in zip(axes, metrics, labels):
        ax.violinplot(adata.obs[metric].dropna(), showmedians=True)
        ax.set_title(label, fontsize=11)
        ax.set_xticks([])
        ax.set_ylabel(label)

    plt.suptitle("QC Metrics — van Galen 2019 AML", fontsize=13, fontweight="bold")
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")
    plt.show()


# ─────────────────────────────────────────────────────────────────────────────
# 3. PREPROCESSING
# ─────────────────────────────────────────────────────────────────────────────

def normalize_and_log(adata: sc.AnnData, target_sum: float = 1e4) -> sc.AnnData:
    """
    Normalize counts per cell and apply log1p transformation.

    Parameters
    ----------
    adata : sc.AnnData
    target_sum : float
        Each cell is normalized to this total count (default: 10,000 = CPM-like).

    Returns
    -------
    sc.AnnData
        Normalized AnnData. Raw counts saved in .raw.
    """
    adata.raw = adata  # Save raw counts before normalization
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    print(f"Normalization done: target_sum={int(target_sum)}, log1p applied.")
    return adata


def select_highly_variable_genes(
    adata: sc.AnnData,
    n_top_genes: int = 3000,
    batch_key: Optional[str] = "sample"
) -> sc.AnnData:
    """
    Select highly variable genes (HVGs) for downstream dimensionality reduction.

    Parameters
    ----------
    adata : sc.AnnData
    n_top_genes : int
        Number of HVGs to select.
    batch_key : str, optional
        Column in .obs to use for batch-aware HVG selection.

    Returns
    -------
    sc.AnnData
        AnnData with .var['highly_variable'] column added.
    """
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        batch_key=batch_key,
        flavor="seurat_v3",
        subset=False
    )
    n_hvg = adata.var["highly_variable"].sum()
    print(f"Highly variable genes selected: {n_hvg} / {adata.n_vars}")
    return adata


def run_pca(adata: sc.AnnData, n_comps: int = 50) -> sc.AnnData:
    """
    Run PCA on highly variable genes.

    Parameters
    ----------
    adata : sc.AnnData
    n_comps : int
        Number of principal components.

    Returns
    -------
    sc.AnnData
        AnnData with .obsm['X_pca'] added.
    """
    sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=True)
    print(f"PCA done: {n_comps} components")
    return adata


def run_harmony_batch_correction(
    adata: sc.AnnData,
    batch_key: str = "sample"
) -> sc.AnnData:
    """
    Apply Harmony batch correction on PCA embedding.

    Parameters
    ----------
    adata : sc.AnnData
        Must have .obsm['X_pca'] already computed.
    batch_key : str
        Column in .obs indicating batch (sample / patient).

    Returns
    -------
    sc.AnnData
        AnnData with .obsm['X_pca_harmony'] added.
    """
    import harmonypy as hm

    pca_matrix = adata.obsm["X_pca"].astype(float)
    meta = adata.obs[[batch_key]]

    harmony_out = hm.run_harmony(pca_matrix, meta, batch_key)
    Z = harmony_out.Z_corr
    # Z_corr is (n_components, n_cells) -> transpose to (n_cells, n_components)
    adata.obsm["X_pca_harmony"] = Z.T if Z.shape[1] == adata.n_obs else Z

    print(f"Harmony batch correction done on key='{batch_key}'")
    return adata


# ─────────────────────────────────────────────────────────────────────────────
# 4. CLUSTERING & UMAP
# ─────────────────────────────────────────────────────────────────────────────

def build_neighbors_and_umap(
    adata: sc.AnnData,
    use_rep: str = "X_pca_harmony",
    n_neighbors: int = 15,
    umap_min_dist: float = 0.3
) -> sc.AnnData:
    """
    Build KNN graph and compute UMAP embedding.

    Parameters
    ----------
    adata : sc.AnnData
    use_rep : str
        Key in .obsm to use for neighbor computation.
    n_neighbors : int
        Number of neighbors for KNN graph.
    umap_min_dist : float
        UMAP min_dist parameter (lower = tighter clusters).

    Returns
    -------
    sc.AnnData
        AnnData with .obsm['X_umap'] added.
    """
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors)
    sc.tl.umap(adata, min_dist=umap_min_dist)
    print(f"UMAP computed (n_neighbors={n_neighbors}, min_dist={umap_min_dist})")
    return adata


def run_leiden_clustering(
    adata: sc.AnnData,
    resolution: float = 0.5
) -> sc.AnnData:
    """
    Run Leiden clustering algorithm.

    Parameters
    ----------
    adata : sc.AnnData
    resolution : float
        Higher = more clusters. Typical range: 0.3 – 1.5.

    Returns
    -------
    sc.AnnData
        AnnData with .obs['leiden'] column added.
    """
    sc.tl.leiden(adata, resolution=resolution, key_added="leiden")
    n_clusters = adata.obs["leiden"].nunique()
    print(f"Leiden clustering done: {n_clusters} clusters (resolution={resolution})")
    return adata


# ─────────────────────────────────────────────────────────────────────────────
# 5. CELL TYPE ANNOTATION
# ─────────────────────────────────────────────────────────────────────────────

# AML cell type marker genes — from van Galen 2019 (Table S2)
AML_MARKER_GENES: Dict[str, List[str]] = {
    "HSC":            ["CD34", "SPINK2", "MLLT3", "CRHBP"],
    "Progenitor":     ["CD34", "MPO", "ELANE", "AZU1"],
    "GMP":            ["MPO", "ELANE", "CTSG", "PRTN3"],
    "Promono":        ["LYZ", "MNDA", "CSF1R", "CD14"],
    "Mono":           ["CD14", "LYZ", "S100A8", "S100A9"],
    "cDC":            ["FCER1A", "CD1C", "HLA-DQA1"],
    "pDC":            ["LILRA4", "CLEC4C", "IL3RA"],
    "Erythroid":      ["HBB", "HBA1", "GYPA", "KLF1"],
    "Platelet/MEP":   ["PF4", "PPBP", "GP9", "ITGA2B"],
    "T cell":         ["CD3D", "CD3E", "TRAC"],
    "NK cell":        ["NCAM1", "GNLY", "NKG7"],
    "B cell":         ["CD79A", "MS4A1", "CD19"],
    "Malignant":      ["CD34", "CD33", "CD117"],  # context-dependent
}


def score_cell_types(adata: sc.AnnData) -> sc.AnnData:
    """
    Compute gene set scores for each AML cell type using scanpy's score_genes.

    Adds one column per cell type in .obs (e.g., 'score_HSC', 'score_Mono').

    Parameters
    ----------
    adata : sc.AnnData

    Returns
    -------
    sc.AnnData
    """
    for cell_type, genes in AML_MARKER_GENES.items():
        # Keep only genes present in the dataset
        available_genes = [g for g in genes if g in adata.var_names]
        if len(available_genes) == 0:
            print(f"  Warning: no marker genes found for {cell_type}")
            continue
        sc.tl.score_genes(
            adata,
            gene_list=available_genes,
            score_name=f"score_{cell_type.replace('/', '_').replace(' ', '_')}",
            use_raw=True
        )
    print("Cell type scoring done.")
    return adata


def assign_cell_type_from_scores(adata: sc.AnnData) -> sc.AnnData:
    """
    Assign cell type label based on highest gene score.

    Parameters
    ----------
    adata : sc.AnnData
        Must have score columns (from score_cell_types).

    Returns
    -------
    sc.AnnData
        With .obs['predicted_cell_type'] column.
    """
    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    score_df = adata.obs[score_cols]
    # Get column name with max score per cell
    best_col = score_df.idxmax(axis=1)
    # Clean label: remove 'score_' prefix
    adata.obs["predicted_cell_type"] = best_col.str.replace("score_", "").str.replace("_", " ")
    return adata


# ─────────────────────────────────────────────────────────────────────────────
# 6. TRAJECTORY ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def run_paga(
    adata: sc.AnnData,
    group_key: str = "predicted_cell_type"
) -> sc.AnnData:
    """
    Run PAGA (Partition-based Graph Abstraction) for trajectory inference.

    Parameters
    ----------
    adata : sc.AnnData
    group_key : str
        Column in .obs to use as node labels in the PAGA graph.

    Returns
    -------
    sc.AnnData
        AnnData with .uns['paga'] populated.
    """
    sc.tl.paga(adata, groups=group_key)
    print(f"PAGA computed on '{group_key}'")
    return adata


def run_diffusion_pseudotime(
    adata: sc.AnnData,
    root_cell_type: str = "HSC",
    cell_type_key: str = "predicted_cell_type"
) -> sc.AnnData:
    """
    Compute Diffusion Pseudotime (DPT) starting from a root cell type.

    Parameters
    ----------
    adata : sc.AnnData
    root_cell_type : str
        Cell type label to use as pseudotime origin (typically HSC).
    cell_type_key : str
        Column in .obs containing cell type labels.

    Returns
    -------
    sc.AnnData
        AnnData with .obs['dpt_pseudotime'] column.
    """
    sc.tl.diffmap(adata)

    # Pick a root cell from the specified cell type
    root_candidates = np.where(adata.obs[cell_type_key] == root_cell_type)[0]
    if len(root_candidates) == 0:
        raise ValueError(f"No cells found with cell type '{root_cell_type}'. "
                         f"Available types: {adata.obs[cell_type_key].unique().tolist()}")

    adata.uns["iroot"] = root_candidates[0]
    sc.tl.dpt(adata)
    print(f"DPT computed from root='{root_cell_type}'")
    return adata


# ─────────────────────────────────────────────────────────────────────────────
# 7. ML — MALIGNANT CELL CLASSIFICATION
# ─────────────────────────────────────────────────────────────────────────────

def prepare_ml_dataset(
    adata: sc.AnnData,
    label_key: str = "CellType",
    malignant_label: str = "Malignant",
    use_hvg: bool = True
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Prepare feature matrix and binary labels for malignant cell classification.

    Parameters
    ----------
    adata : sc.AnnData
    label_key : str
        Column in .obs containing cell type labels (van Galen annotation).
    malignant_label : str
        Value in label_key that corresponds to malignant cells.
    use_hvg : bool
        If True, use only highly variable genes as features.

    Returns
    -------
    X : np.ndarray (n_cells, n_features)
    y : np.ndarray (n_cells,) — binary: 1 = malignant, 0 = normal
    feature_names : List[str] — gene names used as features
    """
    if label_key not in adata.obs.columns:
        raise ValueError(f"Column '{label_key}' not found in adata.obs. "
                         f"Available: {adata.obs.columns.tolist()}")

    if use_hvg and "highly_variable" in adata.var.columns:
        gene_mask = adata.var["highly_variable"].values
    else:
        gene_mask = np.ones(adata.n_vars, dtype=bool)

    X = adata.X[:, gene_mask]
    if hasattr(X, "toarray"):
        X = X.toarray()

    y = (adata.obs[label_key] == malignant_label).astype(int).values
    feature_names = adata.var_names[gene_mask].tolist()

    n_malignant = y.sum()
    n_normal = len(y) - n_malignant
    print(f"ML dataset: {len(y)} cells | Malignant: {n_malignant} | Normal: {n_normal}")
    print(f"Features: {len(feature_names)} genes")

    return X, y, feature_names


def plot_feature_importance(
    feature_names: List[str],
    importances: np.ndarray,
    top_n: int = 20,
    title: str = "Top discriminant genes — Malignant vs Normal",
    save_path: Optional[str] = None
) -> None:
    """
    Bar plot of top feature importances from a tree-based classifier.

    Parameters
    ----------
    feature_names : List[str]
    importances : np.ndarray
    top_n : int
        Number of top genes to show.
    title : str
    save_path : str, optional
    """
    idx = np.argsort(importances)[::-1][:top_n]
    top_genes = [feature_names[i] for i in idx]
    top_scores = importances[idx]

    fig, ax = plt.subplots(figsize=(9, 6))
    bars = ax.barh(range(top_n), top_scores[::-1], color="#c0392b", alpha=0.8)
    ax.set_yticks(range(top_n))
    ax.set_yticklabels(top_genes[::-1], fontsize=9)
    ax.set_xlabel("Feature Importance (Gini)", fontsize=11)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.invert_xaxis()

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")
    plt.show()


# ─────────────────────────────────────────────────────────────────────────────
# 8. PLOTTING UTILITIES
# ─────────────────────────────────────────────────────────────────────────────

def set_plot_style() -> None:
    """Apply consistent matplotlib style for the project."""
    matplotlib.rcParams.update({
        "figure.dpi": 120,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.family": "sans-serif",
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
    })


def save_figure(fig: matplotlib.figure.Figure, path: str, dpi: int = 150) -> None:
    """
    Save a matplotlib figure to disk.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    path : str
        Output path (PNG or SVG).
    dpi : int
    """
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    print(f"Saved: {path}")
