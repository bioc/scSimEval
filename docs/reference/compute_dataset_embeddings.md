# Compute Low-Dimensional Embeddings for Reference and Simulated Single-Cell Datasets

Projects empirical biological reference and simulated count matrices
into low-dimensional representations (UMAP, t-SNE, or PCA) for
side-by-side and multi-method topological benchmarking. Automatically
performs library-size normalization, variance-based feature selection,
principal component analysis, and optional unsupervised clustering.

## Usage

``` r
compute_dataset_embeddings(
  reference,
  simulated,
  reduction = c("umap", "tsne", "pca"),
  n_pcs = 30,
  perplexity = 30,
  n_neighbors = 15,
  min_dist = 0.3,
  seed = 42,
  cell_types = NULL,
  batch = NULL
)
```

## Arguments

- reference:

  Biological reference `SingleCellExperiment` object or count matrix
  (features x cells).

- simulated:

  A single simulated count matrix/`SingleCellExperiment` or a named list
  of them.

- reduction:

  Character string specifying the dimensionality reduction method:
  `"umap"` (default), `"tsne"`, or `"pca"`.

- n_pcs:

  Integer specifying the number of principal components to calculate
  (default: 30).

- perplexity:

  Numeric perplexity for t-SNE (default: 30; automatically adapted for
  small sample sizes).

- n_neighbors:

  Integer number of nearest neighbors for UMAP (default: 15).

- min_dist:

  Numeric minimum distance parameter for UMAP (default: 0.3).

- seed:

  Random seed for reproducibility (default: 42).

- cell_types:

  Optional factor or character vector of cell type annotations for
  reference cells (or named list if per-dataset).

- batch:

  Optional factor or character vector of batch annotations.

## Value

A tidy `data.frame` containing cell coordinates (`Dim1`, `Dim2`),
`Dataset` name, `Dataset_Type` ("Reference" vs "Simulated"),
`Cell_Type`, `Cluster`, `Library_Size`, `Detected_Features`, and
`Batch`.

## Examples

``` r
data(example_scrna, package = "scSimEval")
emb <- compute_dataset_embeddings(
  reference = example_scrna$ref,
  simulated = list("Splatter" = example_scrna$sim),
  reduction = "umap",
  cell_types = example_scrna$cell_types
)
head(emb)
#>              Cell_ID      Dim1        Dim2   Dataset      Role Dataset_Type
#> 1 GTACTTAAGCTCAATA-1  1.281574  2.05128620 Reference Reference    Reference
#> 2 CTTGACGAGTAGCGGG-1  2.460906  0.01632259 Reference Reference    Reference
#> 3 AAAGCACCAACCCTCC-1  2.786740 -0.33215095 Reference Reference    Reference
#> 4 CGCACAATCTTGATGA-1  2.183438  0.60949848 Reference Reference    Reference
#> 5 AGCTACTAGCGAGTAA-1  1.330224  0.52773355 Reference Reference    Reference
#> 6 CTTACTAGTCTCAGCC-1 -1.490763  1.46527765 Reference Reference    Reference
#>   Cell_Type   Cluster Library_Size Detected_Features   Batch
#> 1  Bud_axis Cluster_1           54                28 Batch 1
#> 2  Bud_axis Cluster_4           53                20 Batch 1
#> 3  Bud_axis Cluster_4           39                21 Batch 1
#> 4  Bud_axis Cluster_4           28                17 Batch 1
#> 5  Bud_axis Cluster_4           79                37 Batch 1
#> 6  Bud_axis Cluster_6           55                24 Batch 1
```
