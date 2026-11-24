test_that("Category III CDI (Fang et al. 2022) handles edge cases and appears as single CDI metric in bubble matrix", {
  data("example_scrna", package = "scSimEval", envir = environment())
  data("example_rice_bud", package = "scSimEval", envir = environment())

  # Continuous small floating-point simulated counts should not collapse to all zeros
  set.seed(11)
  m_float <- matrix(stats::runif(40 * 60, 0, 0.35), nrow = 40, ncol = 60)
  m_float[1:15, 1:30] <- m_float[1:15, 1:30] + 0.25
  cl <- factor(rep(c("TypeA", "TypeB"), each = 30))
  cdi_float <- calc_cdi(m_float, cl)
  expect_true(is.finite(cdi_float$cdi_aic))

  # Multi-batch subsampling edge case where a cluster x batch split drops below 3 cells
  m_mb <- matrix(stats::rpois(30 * 120, lambda = 4), nrow = 30, ncol = 120)
  cl_mb <- factor(rep(c("A", "B", "C"), each = 40))
  bt_mb <- rep("B1", 120)
  bt_mb[c(1:3, 41:60, 81:100)] <- "B2"
  cdi_mb <- calc_cdi(m_mb, cl_mb, batch_label = bt_mb, max_cells = 80L)
  expect_true(is.finite(cdi_mb$cdi_aic))

  # Single CDI metric (AIC-based per Fang et al. 2022, no BIC) appears in bubble matrix
  ev1 <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$sim,
                                      compute_bivariate = FALSE, trajectory_num_trees = 100L, verbose = FALSE)
  ev2 <- evaluate_simulation_accuracy(example_scrna$ref, example_scrna$ref,
                                      compute_bivariate = FALSE, trajectory_num_trees = 100L, verbose = FALSE)
  p <- plot_benchmark_bubble_matrix(list(SimA = ev1, RefSelf = ev2))
  cdi_rows <- p$data[p$data$Metric == "cdi", , drop = FALSE]
  expect_equal(nrow(cdi_rows), 2L)
  expect_true(all(as.character(cdi_rows$Display_Metric) == "CDI"))
  expect_false(any(grepl("BIC", as.character(p$data$Metric), ignore.case = TRUE)))
})

test_that("Category III Generative Precision & Recall are scale-consistent and pass synthetic manifold checks", {
  set.seed(42)
  mk_clusters <- function(n_per = 40L, shift = 0, clusters = 1:4) {
    do.call(cbind, lapply(clusters, function(k) {
      m <- matrix(stats::rnorm(30 * n_per, mean = 0, sd = 1), nrow = 30, ncol = n_per)
      m[(k * 5 - 4):(k * 5), ] <- m[(k * 5 - 4):(k * 5), ] + 5 + shift
      m
    }))
  }

  ref_syn <- mk_clusters(40L, 0, 1:4)
  sim_iid <- mk_clusters(40L, 0, 1:4)
  sim_collapse <- mk_clusters(80L, 0, 1:1)
  sim_disjoint <- mk_clusters(40L, 15, 1:4)

  pr_iid <- calc_generative_precision_recall(ref_syn, sim_iid)
  expect_gt(pr_iid$generative_precision, 0.75)
  expect_gt(pr_iid$generative_recall, 0.75)

  pr_col <- calc_generative_precision_recall(ref_syn, sim_collapse)
  expect_gt(pr_col$generative_precision, 0.70)
  expect_lt(pr_col$generative_recall, 0.45)

  pr_dis <- calc_generative_precision_recall(ref_syn, sim_disjoint)
  expect_lt(pr_dis$generative_precision, 0.05)
  expect_lt(pr_dis$generative_recall, 0.05)

  # Scale invariance: continuous scaled counts vs integer UMI counts should not force (1, 0)
  data("example_rice_bud", package = "scSimEval", envir = environment())
  pr_momtf <- calc_generative_precision_recall(
    example_rice_bud$reference$rna,
    example_rice_bud$simulated$scMoMtF$rna
  )
  expect_gt(pr_momtf$generative_precision, 0.05)
  expect_lt(pr_momtf$generative_precision, 0.95)
  expect_gt(pr_momtf$generative_recall, 0.10)
  expect_lt(pr_momtf$generative_recall, 0.95)
})

test_that("Category V precision, recall, F1, and AUC do not inflate null contrasts to 1.0", {
  set.seed(99)
  n_genes <- 80L
  n_cells <- 90L
  ref_mat <- matrix(stats::rpois(n_genes * n_cells, lambda = 4), nrow = n_genes, ncol = n_cells)
  rownames(ref_mat) <- paste0("Gene_", seq_len(n_genes))
  colnames(ref_mat) <- paste0("Cell_", seq_len(n_cells))
  ct <- factor(rep(c("A", "B", "C"), each = 30L))
  bt <- factor(rep(c("Rep1", "Rep2"), length.out = n_cells))
  # Inject DE signal only between C and {A, B}; A vs B has zero true DEGs
  ref_mat[1:15, 61:90] <- ref_mat[1:15, 61:90] + stats::rpois(15 * 30, lambda = 12)

  # Pure noise simulator: zero true DEGs recovered
  sim_noise <- matrix(stats::rpois(n_genes * n_cells, lambda = 4), nrow = n_genes, ncol = n_cells)
  rownames(sim_noise) <- rownames(ref_mat)
  colnames(sim_noise) <- colnames(ref_mat)

  res_noise <- evaluate_deg_fidelity(
    ref_data = ref_mat, sim_data = sim_noise,
    ref_celltypes = ct, sim_celltypes = ct,
    batch_ref = bt, batch_sim = bt
  )
  tbl <- res_noise$deg_summary_table
  p_val <- tbl$Value[tbl$Metric == "precision"]
  r_val <- tbl$Value[tbl$Metric == "recall"]
  f1_val <- tbl$Value[tbl$Metric == "f1_score"]
  expect_lte(p_val, 0.1)
  expect_lte(r_val, 0.1)
  expect_lte(f1_val, 0.1)
  # Summary F1 must equal exact harmonic mean of pooled precision and recall
  expected_f1 <- if ((p_val + r_val) > 0) 2 * p_val * r_val / (p_val + r_val) else 0
  expect_equal(f1_val, expected_f1, tolerance = 1e-10)
})

test_that("Category VI integration ARI & NMI resolve list cell_types and avoid circular self-clustering", {
  data("example_rice_bud", package = "scSimEval", envir = environment())
  ref_m <- list(RNA = example_rice_bud$reference$rna, ATAC = example_rice_bud$reference$atac)
  sim_m <- list(RNA = example_rice_bud$simulated$MOSim$rna, ATAC = example_rice_bud$simulated$MOSim$atac)
  ct_list <- list(
    SummarizedExperiment::colData(sim_m$RNA)$cell_type,
    SummarizedExperiment::colData(sim_m$ATAC)$cell_type
  )

  # Passing list(ct_m1, ct_m2) with SCEs that contain unimodal PCA must NOT produce fake ARI/NMI = 1.0
  res <- evaluate_multiomics_integration(
    ref_multi = ref_m,
    sim_multi = sim_m,
    cell_types = ct_list,
    verbose = FALSE
  )
  expect_lt(res$metrics$ari, 0.5)
  expect_lt(res$metrics$nmi, 0.6)
})

test_that("Category VII edge_flip, geodesic_correlation, and wcor_features are scale-consistent and accurate", {
  data("example_rice_bud", package = "scSimEval", envir = environment())
  ref_rna <- example_rice_bud$reference$rna
  sim_rna <- example_rice_bud$simulated$scMoMtF$rna
  res <- evaluate_trajectory_metrics(
    ref_data = ref_rna,
    sim_data = sim_rna,
    cell_types_ref = SummarizedExperiment::colData(ref_rna)$cell_type,
    cell_types_sim = SummarizedExperiment::colData(sim_rna)$cell_type,
    num_trees = 200L,
    threads = 2L
  )
  for (m in c("him", "edge_flip", "geodesic_correlation", "F1_branches", "wcor_features")) {
    expect_true(is.finite(res[[m]]))
    expect_gte(res[[m]], 0)
    expect_lte(res[[m]], 1)
  }
  # Normalized edge-flip similarity on two trees sharing 3 of 4 edges equals 0.75
  t1 <- list(milestone_network = data.frame(
    from = c("A", "B", "C", "D"), to = c("B", "C", "D", "E"), stringsAsFactors = FALSE
  ))
  t2 <- list(milestone_network = data.frame(
    from = c("A", "B", "C", "C"), to = c("B", "C", "D", "E"), stringsAsFactors = FALSE
  ))
  expect_equal(calc_dynverse_edge_flip(t1, t2), 0.75, tolerance = 1e-8)
})

test_that("Category I distributional fixes handle scaling, ties, and sparse 2D density cleanly", {
  set.seed(42)
  x <- stats::rnorm(120, mean = 10, sd = 1)
  expect_equal(calc_bhattacharyya(x, x), 0, tolerance = 1e-8)
  expect_gt(calc_bhattacharyya(x, 100 * x), 2.0)
  expect_gt(calc_bhattacharyya(x, x + 50), 2.0)
  expect_equal(calc_bhattacharyya(rep(1, 50), rep(1, 50)), 0)
  expect_gt(calc_bhattacharyya(rep(1, 50), rep(100, 50)), 10)

  # Tie-safe Runs test and NN mismatch on zero-inflated vectors
  x_ties <- c(rep(0, 60), stats::runif(40, 0.1, 1))
  y_ties <- c(rep(0, 60), stats::runif(40, 0.1, 1))
  rt_same <- calc_runs_test(x_ties, y_ties)
  expect_lt(abs(rt_same$runs_statistic), 3.0)
  expect_lt(calc_nn_mismatch(x_ties, x_ties), 0.15)
  expect_gt(calc_nn_mismatch(stats::rnorm(50, 0, 0.1), stats::rnorm(50, 10, 0.1)), 0.90)

  # Sparse 2D EMD handles zero-IQR inputs without returning NA
  m_sp <- cbind(c(rep(0, 80), stats::runif(20, 0.1, 2)), c(rep(1, 80), stats::runif(20, 0.2, 0.9)))
  emd_same <- calc_emd_2d(m_sp, m_sp)
  expect_true(is.finite(emd_same))
  expect_lt(emd_same, 1e-4)
})

test_that("Category VI unpaired integration respects unpaired layout when n1 == n2", {
  data("example_soybean_root", package = "scSimEval", envir = environment())
  res_comb <- evaluate_multiomics_integration(
    ref_multi = example_soybean_root$reference$combined,
    sim_multi = example_soybean_root$simulated$scDesign3$combined,
    pairing = "unpaired",
    verbose = FALSE
  )
  expect_equal(nrow(res_comb$embedding), 200L)
  expect_true(all(vapply(res_comb$metrics, function(v) is.finite(v) && v >= 0 && v <= 1, logical(1))))

  res_list <- evaluate_multiomics_integration(
    ref_multi = list(RNA = example_soybean_root$reference$rna, ATAC = example_soybean_root$reference$atac),
    sim_multi = list(RNA = example_soybean_root$simulated$scDesign3$rna, ATAC = example_soybean_root$simulated$scDesign3$atac),
    pairing = "unpaired",
    verbose = FALSE
  )
  expect_equal(nrow(res_list$embedding), 200L)
  expect_true(all(vapply(res_list$metrics, function(v) is.finite(v) && v >= 0 && v <= 1, logical(1))))
})

test_that("All 4 analysis modes return exact canonical metric counts (44/38/56/56) and DE/DA/DE/DA prefixes", {
  data("example_rice_bud", package = "scSimEval", envir = environment())
  data("example_soybean_root", package = "scSimEval", envir = environment())

  count_canonical <- function(tbl) {
    std <- standardize_benchmark_categories(tbl)
    sub_key <- ifelse(
      grepl("^\\(V\\)", std$Category),
      paste0(std$Category, "||", std$Property, "||", std$Metric),
      paste0(std$Category, "||", std$Metric)
    )
    length(unique(sub_key))
  }

  r_rna <- evaluate_simulation_accuracy(
    example_rice_bud$reference$rna, example_rice_bud$simulated$scDesign3$rna,
    elapsed_time = 87.87, memory_mb = 5713.9, sim_cores = 16L,
    trajectory_num_trees = 100L, verbose = FALSE
  )
  expect_equal(count_canonical(r_rna$metrics_summary_table), 44L)
  p_rna <- plot_benchmark_bubble_matrix(list(scDesign3 = r_rna))
  c5_rna <- unique(as.character(p_rna$data$Display_Metric[grepl("Bio-Signal|\\(V\\)", p_rna$data$Display_Category)]))
  expect_true(all(startsWith(c5_rna, "DE ")))

  r_atac <- evaluate_simulation_accuracy(
    example_rice_bud$reference$atac, example_rice_bud$simulated$scDesign3$atac,
    elapsed_time = 87.87, memory_mb = 5713.9, sim_cores = 16L, verbose = FALSE
  )
  expect_equal(count_canonical(r_atac$metrics_summary_table), 38L)
  p_atac <- plot_benchmark_bubble_matrix(list(scDesign3 = r_atac))
  c5_atac <- unique(as.character(p_atac$data$Display_Metric[grepl("Bio-Signal|\\(V\\)", p_atac$data$Display_Category)]))
  expect_true(all(startsWith(c5_atac, "DA ")))

  r_pair <- evaluate_multiomics_accuracy(
    list(RNA = example_rice_bud$reference$rna, ATAC = example_rice_bud$reference$atac),
    list(RNA = example_rice_bud$simulated$scDesign3$rna, ATAC = example_rice_bud$simulated$scDesign3$atac),
    pairing = "paired", elapsed_time = 87.87, memory_mb = 5713.9, sim_cores = 16L,
    trajectory_num_trees = 100L, verbose = FALSE
  )
  expect_equal(count_canonical(r_pair$benchmark_summary_table), 56L)
  p_pair <- plot_benchmark_bubble_matrix(list(scDesign3 = r_pair))
  expect_equal(length(unique(p_pair$data$Metric_Key)), 56L)
  c5_pair <- unique(as.character(p_pair$data$Display_Metric[grepl("Bio-Signal|\\(V\\)", p_pair$data$Display_Category)]))
  expect_equal(length(c5_pair), 12L)
  expect_equal(sum(startsWith(c5_pair, "RNA DE ")), 6L)
  expect_equal(sum(startsWith(c5_pair, "ATAC DA ")), 6L)

  r_unp <- evaluate_multiomics_accuracy(
    example_soybean_root$reference$combined, example_soybean_root$simulated$scMultiSim$combined,
    pairing = "unpaired", elapsed_time = 107.55, memory_mb = 1354.5, sim_cores = 16L,
    trajectory_num_trees = 100L, verbose = FALSE
  )
  expect_equal(count_canonical(r_unp$benchmark_summary_table), 56L)
})

test_that("Regression tests: normalize_counts, benchmark_resource_usage, multi-batch CMS, clamped ATAC curve, Cat VI disambiguation", {
  # 1. normalize_counts normalizes continuous non-integer non-negative counts
  set.seed(11)
  cont_mat <- Matrix::Matrix(matrix(rgamma(500, shape = 0.5, scale = 0.4), nrow = 25, ncol = 20), sparse = TRUE)
  expect_lt(max(cont_mat), 35)
  norm_cont <- normalize_counts(cont_mat)
  expect_gt(mean(as.numeric(norm_cont)), mean(as.numeric(cont_mat)))

  # 2. benchmark_resource_usage reports realistic peak MiB (not divided by 131072)
  res_mem <- benchmark_resource_usage({ x <- numeric(2e5); sum(x) }, sim_cores = 2L)
  expect_true(is.finite(res_mem$memory_mb) && res_mem$memory_mb > 1)
  expect_equal(res_mem$memory_per_core_mib, res_mem$memory_mb / 2, tolerance = 1e-6)

  # 3. Multi-batch CMS (5 well-mixed batches, 100 cells) yields reasonable mixing > 0.2
  set.seed(42)
  coords_5b <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
  batch_5b  <- factor(rep(paste0("rep", 1:5), each = 20))
  cms_5b <- calc_cms(coords = coords_5b, batch_info = batch_5b, threads = 1L)
  expect_true(is.finite(cms_5b$mean_cms) && cms_5b$mean_cms > 0.2)

  # Single-batch evaluate_batch_metrics returns NA gracefully instead of crashing
  data(example_scrna, package = "scSimEval")
  b1 <- evaluate_batch_metrics(example_scrna$ref, factor(rep("OnlyBatch", ncol(example_scrna$ref))), verbose = FALSE)
  expect_true(is.na(b1$pcr_r2))
  expect_true(is.na(b1$cms))

  # 4. Accessibility sparsity curve RMSE stays bounded in [0, 1]
  data(example_rice_bud, package = "scSimEval")
  ac_res <- evaluate_accessibility_sparsity_curve(
    example_rice_bud$reference$atac,
    example_rice_bud$simulated$scMoMtF$atac
  )
  expect_true(is.finite(ac_res$curve_rmse) && ac_res$curve_rmse >= 0 && ac_res$curve_rmse <= 1)

  # 5. Full 56-metric coverage, Cat V/VI disambiguation, per-metric raw heatmap fill, and raw-only direction arrows
  rds_path <- system.file("extdata", "scSimEval_benchmark_results.rds", package = "scSimEval")
  if (nzchar(rds_path) && file.exists(rds_path)) {
    rds_df <- readRDS(rds_path)
    p_hm <- plot_metric_heatmap(rds_df)
    hm_metrics <- unique(as.character(p_hm$data$Display_Metric))
    expect_equal(length(hm_metrics), 56L)
    expect_false(any(grepl("\u2191|\u2193", hm_metrics)))
    expect_equal(length(unique(p_hm$data$Category)), 8L)

    # Raw mode shows (↑)/(↓) arrows, raw values in cells, and per-metric [0, 1] tile fill
    p_hm_raw <- plot_metric_heatmap(rds_df, scale_fill = "raw")
    hm_raw_metrics <- unique(as.character(p_hm_raw$data$Display_Metric))
    expect_equal(length(hm_raw_metrics), 56L)
    expect_true(all(grepl("\\(\u2191\\)|\\(\u2193\\)$", hm_raw_metrics)))
    expect_true(any(grepl("\\(\u2191\\)$", hm_raw_metrics)))
    expect_true(any(grepl("\\(\u2193\\)$", hm_raw_metrics)))
    expect_true(grepl("Score_Norm", rlang::as_label(p_hm_raw$mapping$fill)))
    expect_equal(p_hm_raw$scales$get_scales("fill")$limits, c(0, 1))
    expect_true(any(as.numeric(p_hm_raw$data$Score_Label[p_hm_raw$data$Metric == "core_seconds"]) > 1, na.rm = TRUE))

    p_mds <- expect_no_warning(plot_metric_mds(rds_df, ordination_by = "summaries"))
    expect_equal(nrow(p_mds$data), 56L)
    expect_false(any(grepl("\u2191|\u2193", p_mds$data$Display_Label)))

    p_pca <- plot_metric_pca(rds_df, panel = "both")
    expect_true(inherits(p_pca, "gg") || inherits(p_pca, "patchwork"))
    p_pca_load <- plot_metric_pca(rds_df, panel = "loadings", top_n_loadings = NULL)
    expect_gte(nrow(p_pca_load$data), 50L)
    expect_false(any(grepl("\u2191|\u2193", p_pca_load$data$Display_Label)))

    p_box <- plot_metric_boxplots(rds_df, score_type = "normalized")
    expect_equal(length(unique(paste(p_box$data$Category, p_box$data$Metric))), 56L)
    expect_false(any(grepl("\u2191|\u2193", as.character(p_box$data$Display_Metric))))
    p_box_raw <- plot_metric_boxplots(rds_df, score_type = "raw")
    expect_true(all(grepl("\\(\u2191\\)|\\(\u2193\\)$", as.character(p_box_raw$data$Display_Metric))))

    p_cat5 <- plot_category_metric_bars(rds_df, category = "(V) Biological Signal & Differential Fidelity")
    expect_equal(length(unique(p_cat5$data$Display_Metric)), 12L)
    expect_false(any(grepl("\u2191|\u2193", as.character(p_cat5$data$Display_Metric))))
    p_cat5_raw <- plot_category_metric_bars(rds_df, category = "(V) Biological Signal & Differential Fidelity", score_type = "raw")
    expect_true(all(grepl("\\(\u2191\\)$", as.character(p_cat5_raw$data$Display_Metric))))

    p_ind <- plot_individual_metric_bar(rds_df, metric = "rna_de_f1")
    expect_false(grepl("\u2191|\u2193", p_ind$labels$title))
    p_ind_raw <- plot_individual_metric_bar(rds_df, metric = "rna_de_f1", score_type = "raw")
    expect_true(grepl("\u2191", p_ind_raw$labels$title))
  }
})
