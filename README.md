# 10x Genomics Spatial Transcriptomics

## 🧬 Overview

This repository contains four Jupyter notebooks for spatial transcriptomics analysis using **10x Genomics** platforms. The notebooks were developed and run on **Google Colab**, based on the official tutorials from the [Squidpy](https://squidpy.readthedocs.io) and [Scanpy](https://scanpy-tutorials.readthedocs.io) documentation.

Each notebook is self-contained — it installs its own dependencies, downloads or loads its own dataset, runs the full analysis pipeline, and produces all plots inline.

---

## 🔬 Technologies and Platforms

| Technology | Description |
|------------|-------------|
| **10x Xenium** | Sub-cellular resolution in-situ spatial transcriptomics. Assigns individual transcripts to segmented cells. ~480 gene panel. |
| **10x Visium** | Spot-based spatial gene expression. 55 µm spots arranged in a regular grid over tissue sections. Paired with H&E or fluorescence imaging. |
| **Squidpy** | Python library for spatial single-cell analysis — spatial graphs, neighbourhood enrichment, co-occurrence, Moran's I. |
| **Scanpy** | Industry-standard Python library for single-cell RNA-seq — normalisation, PCA, UMAP, Leiden clustering, differential expression. |
| **SpatialData** | Unified multi-modal spatial data framework used for Xenium data loading and organisation. |

---

## 📓 Notebook Descriptions

---

### 1️⃣ Xenium Human Lung (`01_Xenium_Human_Lung/`)

**Source:** https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_xenium.html  
**Dataset:** [Xenium V1 Human Lung Cancer 2-FOV](https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_human_Lung_2fov/Xenium_V1_human_Lung_2fov_outs.zip) — 161,000 cells × 480 genes

#### What this notebook does

The Xenium platform provides true **single-cell spatial resolution** — every transcript is localised to a specific segmented cell, unlike Visium which assigns expression to 55 µm spots that may capture multiple cells.

#### Full Pipeline

**Step 1 — Package Installation**
- Installs `spatialdata`, `spatialdata-io`, `squidpy`, `scanpy`, `seaborn` with no version pins (allows pip to resolve the latest compatible stack and avoids dependency conflicts)

**Step 2 — Data Loading**
- Downloads the Xenium V1 Human Lung zip (~3 GB) from the 10x Genomics public server using `wget`
- Loads with `spatialdata_io.xenium()` into a `SpatialData` object
- Images (`morphology_focus`, `morphology_mip`) are skipped to avoid a known Zarr chunking error (`TypeError: Expected an iterable of integers`) caused by dask loading large TIFFs with irregular chunk shapes
- Works entirely in memory — no `sdata.write()` call

**Step 3 — Extract AnnData Table**
- Gene expression is in `sdata.tables["table"]` as an `AnnData` object
- `adata.obsm["spatial"]` contains (x, y) tissue coordinates for every cell — populated automatically by the Xenium reader

**Step 4 — Quality Control**
- `sc.pp.calculate_qc_metrics()` computes `total_counts` and `n_genes_by_counts` per cell
- Negative control probe percentage and negative decoding codeword percentage are reported (both should be < 1 %)
- Four QC histograms: total transcripts, unique genes, cell area, nucleus/cell area ratio

**Step 5 — Filtering**
- Cells with fewer than 10 total counts removed
- Genes expressed in fewer than 5 cells removed

**Step 6 — Normalisation & Dimensionality Reduction**
- Raw counts preserved in `adata.layers["counts"]`
- `sc.pp.normalize_total()` — library size normalisation
- `sc.pp.log1p()` — log(x + 1) transformation
- `sc.pp.pca()` — top 50 principal components
- `sc.pp.neighbors()` — k-NN graph in PCA space
- `sc.tl.umap()` — 2D non-linear embedding
- `sc.tl.leiden()` — community detection clustering

**Step 7 — Spatial Statistics**
- `sq.gr.spatial_neighbors()` — Delaunay triangulation graph connecting physically proximate cells
- `sq.gr.centrality_scores()` — closeness, clustering coefficient, degree centrality per cluster
- `sq.gr.co_occurrence()` — probability of finding cluster B near cluster A across spatial radii (run on 50 % subsample)
- `sq.gr.nhood_enrichment()` — permutation z-score for each cluster-pair spatial adjacency

**Step 8 — Spatially Variable Genes**
- `sq.gr.spatial_autocorr(mode="moran")` — Moran's I statistic across all genes, identifying genes with spatially structured expression

**Key Outputs**
- UMAP coloured by clusters, total counts, gene counts
- Spatial scatter of Leiden clusters on tissue
- Centrality score bar plots
- Co-occurrence curves
- Neighbourhood enrichment heatmap
- Moran's I ranked gene table and spatial expression plots

---

### 2️⃣ Visium H&E Mouse Brain (`02_Visium_HnE_Mouse_Brain/`)

**Source:** https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html  
**Dataset:** Pre-processed mouse brain coronal section — loaded directly via `sq.datasets.visium_hne_adata()` and `sq.datasets.visium_hne_image()` (no manual download)

#### What this notebook does

This tutorial demonstrates the full Squidpy analysis pipeline on **Visium H&E data**. Visium captures gene expression across a regular grid of 55 µm spots placed on a tissue section stained with Hematoxylin & Eosin (H&E). The key advantage of this tutorial is combining **morphological image features** extracted from the H&E image with gene expression data.

#### Full Pipeline

**Step 1 — Data Loading**
- Built-in dataset from Squidpy — downloads automatically on first run
- `img`: `ImageContainer` holding the high-resolution H&E tissue image
- `adata`: pre-processed mouse brain AnnData with cluster annotations already added

**Step 2 — Spatial Visualisation of Cluster Annotations**
- `sq.pl.spatial_scatter()` shows the tissue section with each spot coloured by its pre-annotated brain region cluster

**Step 3 — Image Feature Extraction**
- `sq.im.calculate_image_features()` extracts **summary statistics** (mean, standard deviation of pixel intensities) from the H&E image within each Visium spot
- Run at two scales: `scale=1.0` (local fine-grained texture) and `scale=2.0` (broader tissue context)
- Features from both scales are concatenated into `adata.obsm["features"]`
- This creates a spot × image-features matrix alongside the existing spot × gene matrix

**Step 4 — Morphology-based Clustering**
- A custom `cluster_features()` function runs PCA → neighbors → Leiden clustering on the image features alone
- This produces a cluster label based purely on tissue morphology visible in H&E, completely independent of gene expression
- `adata.obs["features_cluster"]` stores the morphology-based labels

**Step 5 — Comparison: Gene vs Morphology Clusters**
- Side-by-side `sq.pl.spatial_scatter()` comparing gene-expression Leiden clusters with morphology-based image clusters
- Concordance between the two confirms that H&E tissue structure matches transcriptomic identity

**Step 6 — Neighbourhood Enrichment**
- `sq.gr.spatial_neighbors()` builds the Visium hexagonal grid connectivity graph
- `sq.gr.nhood_enrichment()` tests which brain region clusters are spatially adjacent more than expected
- Results plotted as a colour-coded heatmap (positive z-score = enriched co-localisation)

**Step 7 — Co-occurrence Analysis**
- `sq.gr.co_occurrence()` calculates distance-based spatial probability scores
- Conditioned on the Hippocampus cluster to show which other regions co-occur with it at increasing distances

**Step 8 — Ligand-Receptor Interaction Analysis**
- `sq.gr.ligrec()` runs the CellPhoneDB statistical framework using the OmniPath database
- Tests all annotated ligand-receptor pairs across all cluster combinations using 100 permutations
- Visualised as a dot plot filtered by the Hippocampus → Pyramidal layer axis with strict significance threshold (α = 1e-4)

**Step 9 — Moran's I Spatially Variable Genes**
- `sq.gr.spatial_autocorr(mode="moran")` evaluated on the top 1,000 highly variable genes
- Top spatially variable genes plotted spatially to confirm their patterned expression

---

### 3️⃣ Visium Fluorescence Mouse Brain (`03_Visium_Fluorescence_Mouse_Brain/`)

**Source:** https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_fluo.html  
**Dataset:** Pre-processed mouse brain Visium fluorescence crop — loaded via `sq.datasets.visium_fluo_adata_crop()` and `sq.datasets.visium_fluo_image_crop()`

#### What this notebook does

This tutorial focuses on **image analysis features** specific to fluorescence microscopy data paired with Visium spots. Unlike the H&E tutorial (which uses staining intensity), fluorescence data allows channel-specific analysis of protein markers. The dataset is a smaller crop of the full brain section for faster execution.

#### Full Pipeline

**Step 1 — Data Loading**
- Built-in Squidpy dataset — pre-processed mouse brain with cluster annotations
- `img`: fluorescence `ImageContainer` (multi-channel)
- `adata`: spot-level gene expression with pre-annotated clusters

**Step 2 — Spatial Cluster Visualisation**
- `sq.pl.spatial_scatter()` overlays cluster labels on the fluorescence tissue image

**Step 3 — Image Segmentation**
- `sq.im.segment()` applies a thresholding-based segmentation on the DAPI channel to identify cell nuclei within each Visium spot
- Segmentation results are stored back into the `ImageContainer`

**Step 4 — Segmentation-based Feature Extraction**
- `sq.im.calculate_image_features()` with `features="segmentation"` computes per-spot morphological statistics derived from the segmented nuclei:
  - Number of segmented cells per spot
  - Mean cell area
  - Ratio of tissue covered by cells
- These are stored in `adata.obsm["features_segmentation"]`

**Step 5 — Summary Feature Extraction**
- In addition to segmentation features, `features="summary"` extracts pixel intensity statistics (mean, standard deviation) from the fluorescence channels per spot
- Stored in `adata.obsm["features_summary"]`

**Step 6 — Feature Clustering**
- Features from segmentation and summary statistics are combined
- Leiden clustering on image features produces `adata.obs["features_cluster"]`

**Step 7 — Comparison and Interpretation**
- Side-by-side spatial plots compare the image-based cluster with the gene-expression cluster
- Concordance shows that fluorescence morphology tracks transcriptomic cell identity

---

### 4️⃣ Scanpy Spatial Basic Analysis — Visium (`04_Scanpy_Spatial_Basic_Analysis/`)

**Source:** https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html  
**Dataset:** 10x Visium human lymph node (loaded via `sc.datasets.visium_sge()`) + MERFISH mouse cortex example

#### What this notebook does

This is a **Scanpy-native** tutorial that shows how spatial transcriptomics data fits naturally into the standard Scanpy single-cell analysis workflow. It focuses on Visium data but also includes a MERFISH example. The key message is that spatial context (tissue image + spot coordinates) is stored alongside gene expression in AnnData's `obsm` and `uns` slots, making spatial analysis a seamless extension of standard workflows.

#### Full Pipeline

**Step 1 — Reading the Data**
- `sc.read_visium()` / `sc.datasets.visium_sge()` loads the Visium SpaceRanger output
- Spot coordinates in `adata.obsm["spatial"]`
- H&E image in `adata.uns["spatial"][library_id]["images"]["hires"]`
- Spot size and scale factors stored in `adata.uns["spatial"]`

**Step 2 — QC and Preprocessing**
- `sc.pp.calculate_qc_metrics()` per spot
- Filter spots by minimum gene count and filter genes by minimum spot count
- `sc.pp.normalize_total()` + `sc.pp.log1p()`
- `sc.pp.highly_variable_genes()` selects informative genes for dimensionality reduction

**Step 3 — Manifold Embedding and Clustering**
- `sc.pp.pca()` → `sc.pp.neighbors()` → `sc.tl.umap()` → `sc.tl.leiden()`
- Standard Scanpy single-cell workflow applied to spatial spots

**Step 4 — Spatial Visualisation**
- `sc.pl.spatial()` overlays cluster labels, gene expression, and QC metrics directly onto the H&E image
- Demonstrates how spatial coordinates stored in `adata.obsm["spatial"]` are used for tissue-anchored plotting

**Step 5 — Cluster Marker Genes**
- `sc.tl.rank_genes_groups()` — Wilcoxon rank-sum test for marker genes per cluster
- `sc.pl.rank_genes_groups_dotplot()` — dot plot of top markers per spatial cluster

**Step 6 — MERFISH Example**
- Loads a MERFISH mouse visual cortex dataset
- Shows that the same Scanpy workflow applies to other spatial platforms with minor adjustments

---

## 🛠️ Environment and Dependencies

All notebooks were run on **Google Colab** (Python 3.10). Each notebook installs its own dependencies in Cell 1.

### Core Packages

| Package | Purpose |
|---------|---------|
| `scanpy` | Single-cell analysis: normalisation, PCA, UMAP, Leiden, marker genes |
| `squidpy` | Spatial analysis: graphs, neighbourhood enrichment, Moran's I, image features |
| `spatialdata` | Multi-modal spatial data container (Xenium notebook) |
| `spatialdata-io` | Readers for Xenium, Visium, and other spatial platforms |
| `anndata` | Core data structure (AnnData) for single-cell matrices |
| `dask` | Lazy loading of large spatial arrays |
| `seaborn` | Statistical distribution plots for QC |
| `matplotlib` | All plot rendering |

### Key Technical Notes

**No version pinning:** Packages are installed without pinned version numbers. Manually specifying old versions (e.g. `spatialdata==0.2.5 dask==2023.12.1`) causes `ResolutionImpossible` errors because they are mutually incompatible. Letting pip resolve the latest compatible stack avoids all conflicts.

**Xenium image loading:** The Xenium morphology images (20,000 × 51,000 px, multi-scale TIFFs) are loaded by dask with irregular chunk shapes. Writing these to Zarr causes `TypeError: Expected an iterable of integers. Got ((1,), (3553,), (4096, 1695))`. The fix is to load with `morphology_focus=False, morphology_mip=False` and never call `sdata.write()`.

---

## 📊 Methods Summary

### Spatial Graph Construction
All spatial statistics require a connectivity graph. For **Visium** data, `sq.gr.spatial_neighbors()` uses the regular hexagonal grid. For **Xenium** (arbitrary coordinates), `coord_type="generic", delaunay=True` builds a Delaunay triangulation connecting naturally proximate cells.

### Neighbourhood Enrichment
Compares the observed number of neighbouring spots/cells between each cluster pair against a null distribution from 1,000 random permutations. The resulting z-score indicates whether two clusters are spatially co-localised (positive) or mutually exclusive (negative).

### Co-occurrence Probability
Distance-based score computed directly from (x, y) coordinates:

$$\text{score}(r) = \frac{p(\text{cluster}_B \mid \text{cluster}_A,\, r)}{p(\text{cluster}_B)}$$

Score > 1 means the clusters are more often found together at radius r than expected by chance.

### Moran's I Spatial Autocorrelation
Measures whether a gene's expression is spatially structured. I ≈ +1 means nearby spots share similar expression (spatially clustered); I ≈ 0 means random; I ≈ −1 means dispersed.

### Image Feature Extraction (Visium)
For H&E and fluorescence Visium data, `sq.im.calculate_image_features()` extracts pixel statistics within the tissue image crop surrounding each spot. These features form a second data matrix (spot × image-features) that can be clustered independently to compare morphology-based with transcriptome-based groupings.

## 🔗 Tutorial Links

| Notebook | Official Tutorial |
|----------|------------------|
| Xenium Human Lung | https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_xenium.html |
| Visium H&E Mouse Brain | https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html |
| Visium Fluorescence Mouse Brain | https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_fluo.html |
| Scanpy Spatial Basic Analysis | https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html |
