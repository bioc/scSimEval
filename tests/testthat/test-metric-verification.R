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
