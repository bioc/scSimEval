# ==============================================================================
# (VI) Multiomics Integration & Alignment
# Multiomics Integration & Cross-Modal Alignment Fidelity
# Gold-standard benchmarking metrics extracted from scib, scMultiSim, scCross, Matilda, scMoMtF
# ==============================================================================

#' Calculate Batch Average Silhouette Width (Batch ASW)
#'
#' Evaluates batch mixing within each cell type in integrated or joint multiomics
#' embedding space, following the definitions in \code{scib} and \code{scMultiSim}.
#' For each cell type, the silhouette width of cells with respect to batch labels
#' is computed (\eqn{-1 \le s_i \le 1}), where values near 0 indicate well-mixed batches.
#' The per-cell-type score is \eqn{1 - |s_i|}, averaged across all cells and cell types.
#'
#' @param embed Numeric matrix of cell coordinates (cells x dimensions, e.g. PCA or joint embedding).
#' @param batch Vector or factor of batch or modality labels for each cell.
#' @param cell_type Optional vector or factor of cell-type labels for each cell.
#' @param threads Number of CPU threads (default 1).
#'
#' @return Numeric scalar between 0 and 1, where 1 indicates optimal batch mixing.
#' @references
#' Luecken et al. (2022) Benchmarking atlas-level data integration in single-cell genomics. Nature Methods.
#' Li et al. (2023) scMultiSim: simulation of multi-modality single cell data guided by gene regulatory networks. Cell Systems.
#' @examples
#' set.seed(42)
#' emb <- matrix(rnorm(100 * 10), 100, 10)
#' bt <- rep(c("B1", "B2"), each = 50)
#' ct <- rep(c("T1", "T2"), 50)
#' calc_batch_asw(emb, bt, ct)
#' @keywords internal
#' @export
calc_batch_asw <- function(embed, batch, cell_type = NULL, threads = 1) {
  if (is.null(embed) || nrow(embed) < 3) return(NA_real_)
  n_cells <- nrow(embed)
  batch <- as.character(batch)
  if (length(batch) != n_cells) {
    batch <- if (length(batch) > n_cells) batch[seq_len(n_cells)] else rep(batch, length.out = n_cells)
  }
  if (!is.null(cell_type)) {
    cell_type <- as.character(cell_type)
    if (length(cell_type) != n_cells) {
      cell_type <- if (length(cell_type) > n_cells) cell_type[seq_len(n_cells)] else rep(cell_type, length.out = n_cells)
    }
  }
  if (length(unique(batch)) < 2) return(1.0)
  threads <- if (is.null(threads) || length(threads) == 0L || is.na(threads[1L])) 1L else max(1L, as.integer(threads[1L]))
  has_par_sil <- exists(".parallel_batch_silhouette", mode = "function")

  if (!is.null(cell_type) && length(unique(cell_type)) >= 2) {
    cts <- unique(cell_type)
    valid_groups <- list()
    for (ct in cts) {
      idx <- which(cell_type == ct)
      idx <- idx[idx <= n_cells]
      if (length(idx) >= 4 && length(unique(batch[idx])) >= 2) {
        sub_batch <- as.numeric(as.factor(batch[idx]))
        tb <- table(sub_batch)
        if (all(tb >= 2)) {
          valid_groups[[length(valid_groups) + 1L]] <- list(idx = idx, sub_batch = sub_batch)
        }
      }
    }
    if (length(valid_groups) > 0L) {
      inner_th <- max(1L, as.integer(floor(threads / max(1L, length(valid_groups)))))
      eval_grp <- function(vg, grp_threads) {
        sub_embed <- embed[vg$idx, , drop = FALSE]
        s_vals <- if (has_par_sil && (grp_threads > 1L || !requireNamespace("cluster", quietly = TRUE))) {
          .parallel_batch_silhouette(sub_embed, vg$sub_batch, threads = grp_threads, return_scores = TRUE)
        } else if (requireNamespace("cluster", quietly = TRUE)) {
          d_mat <- stats::dist(sub_embed)
          sil <- tryCatch(cluster::silhouette(vg$sub_batch, d_mat), error = function(e) NULL)
          if (!is.null(sil) && inherits(sil, "silhouette")) sil[, 3] else numeric(0)
        } else {
          numeric(0)
        }
        if (length(s_vals) > 0L) mean(1 - abs(s_vals), na.rm = TRUE) else NA_real_
      }
      asw_per_ct <- if (threads > 1L && length(valid_groups) >= 2L && exists(".sc_parallel_lapply", mode = "function")) {
        unlist(.sc_parallel_lapply(valid_groups, function(vg) eval_grp(vg, inner_th), threads = min(threads, length(valid_groups))), use.names = FALSE)
      } else {
        vapply(valid_groups, function(vg) eval_grp(vg, threads), numeric(1))
      }
      asw_per_ct <- asw_per_ct[!is.na(asw_per_ct)]
      if (length(asw_per_ct) > 0L) {
        return(max(0, min(1, mean(asw_per_ct))))
      }
    }
  }

  # Fallback to global batch silhouette
  b_num <- as.numeric(as.factor(batch))
  s_vals <- if (has_par_sil && ((threads > 1L && n_cells >= 60L) || !requireNamespace("cluster", quietly = TRUE))) {
    .parallel_batch_silhouette(embed, b_num, threads = threads, return_scores = TRUE)
  } else if (requireNamespace("cluster", quietly = TRUE)) {
    d_mat <- stats::dist(embed)
    sil <- tryCatch(cluster::silhouette(b_num, d_mat), error = function(e) NULL)
    if (!is.null(sil) && inherits(sil, "silhouette")) sil[, 3] else numeric(0)
  } else {
    numeric(0)
  }
  if (length(s_vals) > 0L) {
    return(max(0, min(1, mean(1 - abs(s_vals), na.rm = TRUE))))
  }
  1.0
}


#' Calculate Cell-Type Average Silhouette Width (Cell-Type ASW)
#'
#' Evaluates biological conservation and cell-type cluster separation in integrated
#' or joint multiomics space, following \code{scCross} and \code{scib}.
#' The silhouette width of cells with respect to ground-truth cell-type labels is
#' computed and scaled to \eqn{[0, 1]}{[0, 1]} via \eqn{(\mathrm{ASW} + 1) / 2}{(ASW + 1) / 2}.
#'
#' @param embed Numeric matrix of cell coordinates (cells x dimensions).
#' @param cell_type Vector or factor of cell-type labels for each cell.
#' @param threads Number of CPU threads (default 1).
#'
#' @return Numeric scalar between 0 and 1, where 1 indicates distinct, tight cell-type clusters.
#' @references
#' Ding et al. (2023) scCross: deep generative modeling for single-cell cross-modality integration and perturbation analysis.
#' Luecken et al. (2022) Benchmarking atlas-level data integration in single-cell genomics. Nature Methods.
#' @examples
#' set.seed(42)
#' emb <- matrix(rnorm(100 * 10), 100, 10)
#' ct <- rep(c("T1", "T2"), each = 50)
#' calc_cell_type_asw(emb, ct)
#' @keywords internal
#' @export
calc_cell_type_asw <- function(embed, cell_type, threads = 1) {
  if (is.null(embed) || nrow(embed) < 3) return(NA_real_)
  n_cells <- nrow(embed)
  cell_type <- as.character(cell_type)
  if (length(cell_type) != n_cells) {
    cell_type <- if (length(cell_type) > n_cells) cell_type[seq_len(n_cells)] else rep(cell_type, length.out = n_cells)
  }
  if (length(unique(cell_type)) < 2) return(1.0)
  threads <- if (is.null(threads) || length(threads) == 0L || is.na(threads[1L])) 1L else max(1L, as.integer(threads[1L]))

  n <- nrow(embed)
  if (n > 3000) {
    idx <- sample(seq_len(n), 3000)
    embed <- embed[idx, , drop = FALSE]
    cell_type <- cell_type[idx]
  }
  ct_num <- as.numeric(as.factor(cell_type))
  tb <- table(ct_num)
  if (any(tb < 2)) {
    valid_classes <- as.numeric(names(tb)[tb >= 2])
    keep <- ct_num %in% valid_classes
    embed <- embed[keep, , drop = FALSE]
    ct_num <- ct_num[keep]
  }
  if (length(unique(ct_num)) < 2) return(1.0)

  raw_asw <- if (exists(".parallel_batch_silhouette", mode = "function") &&
                 ((threads > 1L && nrow(embed) >= 60L) || !requireNamespace("cluster", quietly = TRUE))) {
    .parallel_batch_silhouette(embed, ct_num, threads = threads)
  } else if (requireNamespace("cluster", quietly = TRUE)) {
    d_mat <- stats::dist(embed)
    sil <- tryCatch(cluster::silhouette(ct_num, d_mat), error = function(e) NULL)
    if (!is.null(sil) && inherits(sil, "silhouette")) mean(sil[, 3], na.rm = TRUE) else NA_real_
  } else {
    NA_real_
  }
  if (!is.na(raw_asw)) {
    return(max(0, min(1, (raw_asw + 1) / 2)))
  }
  1.0
}


#' Calculate KNN Graph Connectivity
#'
#' Evaluates whether cells of the same cell type form a connected component in
#' the k-nearest neighbor (k-NN) graph across modalities/batches, following
#' \code{scCross}, \code{scMultiSim}, and \code{scib}.
#' For each cell type, the relative size of the largest connected component
#' (LCC) in the induced subgraph is calculated: \eqn{|\mathrm{LCC}| / |V_c|}{|LCC| / |V_c|}. The overall score
#' is the mean across all cell types.
#'
#' @param embed Numeric matrix of cell coordinates (cells x dimensions).
#' @param cell_type Vector or factor of cell-type labels for each cell.
#' @param k Number of nearest neighbors to construct the graph (default 15).
#'
#' @return Numeric scalar between 0 and 1, where 1 indicates complete graph connectivity.
#' @references
#' Ding et al. (2023) scCross: deep generative modeling for single-cell cross-modality integration.
#' Li et al. (2023) scMultiSim: simulation of multi-modality single cell data. Cell Systems.
#' Luecken et al. (2022) Benchmarking atlas-level data integration. Nature Methods.
#' @examples
#' set.seed(42)
#' emb <- matrix(rnorm(60 * 5), 60, 5)
#' ct <- rep(c("A", "B"), each = 30)
#' calc_graph_connectivity(emb, ct, k = 10)
#' @keywords internal
#' @export
calc_graph_connectivity <- function(embed, cell_type, k = 15) {
  if (is.null(embed) || nrow(embed) < 2) return(NA_real_)
  n_cells <- nrow(embed)
  cell_type <- as.character(cell_type)
  if (length(cell_type) != n_cells) {
    cell_type <- if (length(cell_type) > n_cells) cell_type[seq_len(n_cells)] else rep(cell_type, length.out = n_cells)
  }
  cts <- unique(cell_type)
  if (length(cts) < 1) return(1.0)
  
  k <- max(1, min(k, n_cells - 1))

  adj <- matrix(FALSE, nrow = n_cells, ncol = n_cells)
  if (requireNamespace("RANN", quietly = TRUE)) {
    nn_mat <- RANN::nn2(as.matrix(embed), k = k + 1L)$nn.idx[, -1L, drop = FALSE]
    for (i in seq_len(n_cells)) {
      nn <- nn_mat[i, ]
      adj[i, nn] <- TRUE
      adj[nn, i] <- TRUE
    }
  } else {
    d_mat <- as.matrix(stats::dist(embed))
    diag(d_mat) <- Inf
    for (i in seq_len(n_cells)) {
      nn <- order(d_mat[i, ])[seq_len(k)]
      adj[i, nn] <- TRUE
      adj[nn, i] <- TRUE
    }
  }
  
  conn_scores <- numeric(length(cts))
  for (ci in seq_along(cts)) {
    c_label <- cts[ci]
    idx <- which(cell_type == c_label)
    idx <- idx[idx <= n_cells]
    if (length(idx) <= 1) {
      conn_scores[ci] <- 1.0
      next
    }
    
    sub_adj <- adj[idx, idx, drop = FALSE]
    n_sub <- length(idx)
    
    visited <- rep(FALSE, n_sub)
    comp_sizes <- integer(0)
    
    for (start_node in seq_len(n_sub)) {
      if (!visited[start_node]) {
        queue <- start_node
        visited[start_node] <- TRUE
        size <- 0
        while (length(queue) > 0) {
          curr <- queue[1]
          queue <- queue[-1]
          size <- size + 1
          nbrs <- which(sub_adj[curr, ] & !visited)
          if (length(nbrs) > 0) {
            visited[nbrs] <- TRUE
            queue <- c(queue, nbrs)
          }
        }
        comp_sizes <- c(comp_sizes, size)
      }
    }
    
    max_comp <- if (length(comp_sizes) > 0) max(comp_sizes) else 1
    conn_scores[ci] <- max_comp / n_sub
  }
  
  mean_conn <- mean(conn_scores, na.rm = TRUE)
  max(0, min(1, mean_conn))
}


#' Safe K-Means Clustering Guarding Unique Row Count
#' @keywords internal
#' @noRd
.safe_kmeans_labels <- function(embed, n_clusters, n_types) {
  k_clust <- if (!is.null(n_clusters)) as.integer(n_clusters) else as.integer(n_types)
  n_uniq <- nrow(unique(embed))
  if (is.na(n_uniq) || n_uniq < 2L) return(NULL)
  k_clust <- max(2L, min(k_clust, n_uniq))
  km <- tryCatch(
    stats::kmeans(embed, centers = k_clust, nstart = 25L, iter.max = 100L),
    error = function(e) {
      tryCatch(
        stats::kmeans(embed, centers = k_clust, nstart = 15L, iter.max = 100L, algorithm = "MacQueen"),
        error = function(e2) NULL
      )
    }
  )
  if (is.null(km)) return(NULL)
  km$cluster
}

#' Calculate Integration Adjusted Rand Index (ARI)
#'
#' Computes the Adjusted Rand Index (ARI) evaluating the concordance between
#' clusters identified in the integrated/joint multiomics embedding and ground-truth
#' cell-type labels, following \code{scib}, \code{scCross}, \code{scMultiSim},
#' \code{Matilda}, and \code{scMoMtF}.
#'
#' @param embed Numeric matrix of cell coordinates (cells x dimensions).
#' @param cell_type Vector or factor of ground-truth cell-type labels.
#' @param n_clusters Number of clusters (defaults to the number of unique cell types).
#'
#' @return Numeric scalar between 0 and 1, where 1 indicates perfect clustering match.
#' @references
#' Hubert & Arabie (1985) Comparing partitions. Journal of Classification.
#' Luecken et al. (2022) Benchmarking atlas-level data integration. Nature Methods.
#' Ding et al. (2023) scCross.
#' Li et al. (2023) scMultiSim. Cell Systems.
#' Yang et al. (2023) Matilda. Nucleic Acids Research.
#' Lan et al. (2024) scMoMtF.
#' @examples
#' set.seed(42)
#' emb <- matrix(rnorm(60 * 4), 60, 4)
#' ct <- rep(c("TypeA", "TypeB"), each = 30)
#' calc_integration_ari(emb, ct)
#' @keywords internal
#' @export
calc_integration_ari <- function(embed, cell_type, n_clusters = NULL) {
  if (is.null(embed) || nrow(embed) < 2) return(NA_real_)
  n_cells <- nrow(embed)
  cell_type <- as.character(cell_type)
  if (length(cell_type) != n_cells) {
    cell_type <- if (length(cell_type) > n_cells) cell_type[seq_len(n_cells)] else rep(cell_type, length.out = n_cells)
  }
  n_types <- length(unique(cell_type))
  if (n_types <= 1) return(1.0)

  pred_cl <- .safe_kmeans_labels(embed, n_clusters, n_types)
  if (is.null(pred_cl)) return(NA_real_)
  
  if (requireNamespace("mclust", quietly = TRUE)) {
    ari <- mclust::adjustedRandIndex(pred_cl, cell_type)
  } else {
    tab <- table(pred_cl, cell_type)
    n <- sum(tab)
    a_i <- rowSums(tab)
    b_j <- colSums(tab)
    sum_n2 <- sum(choose(tab, 2))
    sum_a2 <- sum(choose(a_i, 2))
    sum_b2 <- sum(choose(b_j, 2))
    n2 <- choose(n, 2)
    exp_idx <- (sum_a2 * sum_b2) / n2
    max_idx <- 0.5 * (sum_a2 + sum_b2)
    ari <- if (max_idx == exp_idx) 1.0 else (sum_n2 - exp_idx) / (max_idx - exp_idx)
  }
  max(0, min(1, ari))
}


#' Calculate Integration Normalized Mutual Information (NMI)
#'
#' Computes Normalized Mutual Information (NMI) between clusters in the integrated
#' multiomics space and ground-truth cell-type labels, following \code{scib},
#' \code{scCross}, \code{scMultiSim}, \code{Matilda}, and \code{scMoMtF}.
#'
#' @param embed Numeric matrix of cell coordinates (cells x dimensions).
#' @param cell_type Vector or factor of ground-truth cell-type labels.
#' @param n_clusters Number of clusters (defaults to the number of unique cell types).
#'
#' @return Numeric scalar between 0 and 1, where 1 indicates perfect mutual information.
#' @references
#' Luecken et al. (2022) Benchmarking atlas-level data integration. Nature Methods.
#' Ding et al. (2023) scCross.
#' Li et al. (2023) scMultiSim. Cell Systems.
#' Yang et al. (2023) Matilda. Nucleic Acids Research.
#' Lan et al. (2024) scMoMtF.
#' @examples
#' set.seed(42)
#' emb <- matrix(rnorm(60 * 4), 60, 4)
#' ct <- rep(c("TypeA", "TypeB"), each = 30)
#' calc_integration_nmi(emb, ct)
#' @keywords internal
#' @export
calc_integration_nmi <- function(embed, cell_type, n_clusters = NULL) {
  if (is.null(embed) || nrow(embed) < 2) return(NA_real_)
  n_cells <- nrow(embed)
  cell_type <- as.character(cell_type)
  if (length(cell_type) != n_cells) {
    cell_type <- if (length(cell_type) > n_cells) cell_type[seq_len(n_cells)] else rep(cell_type, length.out = n_cells)
  }
  n_types <- length(unique(cell_type))
  if (n_types <= 1) return(1.0)

  pred_cl <- .safe_kmeans_labels(embed, n_clusters, n_types)
  if (is.null(pred_cl)) return(NA_real_)
  
  tab <- table(pred_cl, cell_type)
  total <- sum(tab)
  p_xy <- tab / total
  p_x <- rowSums(p_xy)
  p_y <- colSums(p_xy)
  
  h_x <- -sum(p_x[p_x > 0] * log(p_x[p_x > 0]))
  h_y <- -sum(p_y[p_y > 0] * log(p_y[p_y > 0]))
  
  p_outer <- outer(p_x, p_y)
  nz <- p_xy > 0
  mi <- sum(p_xy[nz] * log(p_xy[nz] / p_outer[nz]))
  
  denom <- (h_x + h_y) / 2
  if (denom <= 0) return(1.0)
  nmi <- mi / denom
  max(0, min(1, nmi))
}


.fast_mod_pca <- function(m, n_pcs = 15) {
  m_norm <- normalize_counts(m)
  if (nrow(m_norm) > 2000L) {
    gv <- fast_row_vars(m_norm)
    hvg_idx <- order(gv, decreasing = TRUE)[seq_len(2000L)]
    m_hvg <- as.matrix(m_norm[hvg_idx, , drop = FALSE])
  } else if (inherits(m_norm, "Matrix")) {
    m_hvg <- as.matrix(m_norm)
  } else {
    m_hvg <- m_norm
  }
  fast_pca(t(m_hvg), n_pcs = min(n_pcs, ncol(m_hvg) - 1L, nrow(m_hvg) - 1L))
}

#' Evaluate Multiomics Integration & Cross-Modal Alignment Fidelity (Category VI)
#'
#' Evaluates the multiomics integration and cross-modal alignment fidelity of simulated
#' single-cell multiomics data across 5 gold-standard integration metrics extracted from
#' key single-cell benchmarking suites (\code{scib}, \code{scMultiSim}, \code{scCross},
#' \code{Matilda}, and \code{scMoMtF}):
#' \enumerate{
#'   \item \strong{Batch ASW (\code{batch_asw}):} Average Silhouette Width evaluating batch mixing within cell types (scib, scMultiSim).
#'   \item \strong{Cell-Type ASW (\code{cell_type_asw}):} Average Silhouette Width evaluating cell-type bio-conservation in the multiomics space (scCross, scib).
#'   \item \strong{Graph Connectivity (\code{graph_connectivity}):} Relative size of the largest connected component in the k-NN graph for each cell type (scCross, scMultiSim, scib).
#'   \item \strong{ARI (\code{ari}):} Adjusted Rand Index comparing joint multiomics clustering against ground-truth cell types (scib, scCross, scMultiSim, Matilda, scMoMtF).
#'   \item \strong{NMI (\code{nmi}):} Normalized Mutual Information quantifying information concordance between joint clustering and cell types (scib, scCross, scMultiSim, Matilda, scMoMtF).
#' }
#'
#' @param ref_multi Optional reference multiomics list of matrices (e.g. \code{list(rna = ..., atac = ...)}).
#' @param sim_multi Simulated multiomics list of matrices or a matrix representing joint embedding coordinates.
#' @param cell_types Optional cell-type vector or list/data.frame of annotations for each cell.
#' @param batch_info Optional batch/donor/replicate vector or list/data.frame of annotations for each cell.
#' @param embed Optional precomputed joint embedding matrix (cells x dimensions). If NULL, computed automatically from \code{sim_multi}.
#' @param reduction_name Optional name of the precomputed dimensional reduction to extract from \code{SingleCellExperiment} objects.
#' @param n_pcs Number of principal components per modality to retain for joint embedding (default 15).
#' @param k Number of nearest neighbors for graph-based metrics (default 30).
#' @param threads Number of CPU threads for parallel processing (default 1).
#' @param pairing Optional pairing mode (\code{"paired"} or \code{"unpaired"}); auto-detected from cell barcodes if \code{NULL}.
#' @param verbose Logical; whether to print progress messages (default TRUE).
#'
#' @return A list of class \code{scSimEval_integration} containing:
#'   \item{metrics}{Named list of the 5 calculated integration metric values (batch_asw, cell_type_asw, graph_connectivity, ari, nmi).}
#'   \item{integration_summary_table}{A tidy data.frame containing the 5 metrics with columns \code{Category}, \code{Property}, \code{Metric}, \code{Value}, and \code{Modality}.}
#'   \item{embedding}{The joint multiomics embedding matrix used for evaluation.}
#' @examples
#' data(example_multiomics, package = "scSimEval")
#' ct <- rep(c("A", "B"), each = 50)
#' res <- evaluate_multiomics_integration(
#'   sim_multi = example_multiomics$sim_multi,
#'   cell_types = ct
#' )
#' unlist(res[vapply(res, function(x) is.numeric(x) && length(x) == 1L, logical(1))])
#' @export
evaluate_multiomics_integration <- function(
  ref_multi = NULL,
  sim_multi,
  cell_types = NULL,
  batch_info = NULL,
  embed = NULL,
  reduction_name = NULL,
  n_pcs = 15,
  k = 30,
  threads = 1,
  pairing = NULL,
  verbose = TRUE
) {
  threads <- .set_parallel_threads(threads)
  # Enforce SingleCellExperiment input convention
  val_sim <- validate_multiomics_input(sim_multi, "sim_multi", reduction_name = reduction_name)
  sim_multi <- list(val_sim$m1, val_sim$m2)
  if (!is.null(ref_multi)) {
    val_ref <- validate_multiomics_input(ref_multi, "ref_multi", reduction_name = reduction_name)
    ref_multi <- list(val_ref$m1, val_ref$m2)
  }
  
  if (is.null(embed) && !is.null(val_sim$embedding)) {
    embed <- val_sim$embedding
    if (verbose) message("Using precomputed dimensional reduction from sim_multi for Category VI.")
  }
  if (is.null(cell_types) && !is.null(val_sim$cell_types)) {
    cell_types <- val_sim$cell_types
  }
  if (is.null(batch_info) && !is.null(val_sim$batch_info)) {
    batch_info <- val_sim$batch_info
  }

  # 1. Resolve Joint Embedding
  joint_embed <- embed
  ct_vec <- NULL
  bt_vec <- NULL
  p1_cached <- NULL
  p2_cached <- NULL
  is_paired_layout <- FALSE
  n1 <- if (is.list(sim_multi) && length(sim_multi) >= 2L) ncol(sim_multi[[1L]]) else NA_integer_
  n2 <- if (is.list(sim_multi) && length(sim_multi) >= 2L) ncol(sim_multi[[2L]]) else NA_integer_
  c1 <- if (is.list(sim_multi) && length(sim_multi) >= 2L) colnames(sim_multi[[1L]]) else NULL
  c2 <- if (is.list(sim_multi) && length(sim_multi) >= 2L) colnames(sim_multi[[2L]]) else NULL
  detect_paired <- if (!is.null(pairing) && tolower(as.character(pairing)[1L]) == "unpaired") {
    FALSE
  } else if (!is.null(pairing) && tolower(as.character(pairing)[1L]) == "paired") {
    !is.na(n1) && !is.na(n2) && n1 == n2
  } else if (!is.na(n1) && !is.na(n2) && n1 == n2) {
    if (!is.null(c1) && !is.null(c2)) identical(c1, c2) else TRUE
  } else {
    FALSE
  }

  if (is.null(joint_embed)) {
    if (is.matrix(sim_multi) || (is.data.frame(sim_multi) && !is.list(sim_multi))) {
      joint_embed <- as.matrix(sim_multi)
    } else if (is.list(sim_multi) && length(sim_multi) >= 2L) {
      m1 <- sim_multi[[1L]]
      m2 <- sim_multi[[2L]]

      # Paired multiomics (matched barcodes/cells)
      if (detect_paired) {
        is_paired_layout <- TRUE
        p1 <- .fast_mod_pca(m1, n_pcs = n_pcs)
        p2 <- .fast_mod_pca(m2, n_pcs = n_pcs)
        p1_cached <- p1
        p2_cached <- p2
        s1 <- scale(p1)
        s2 <- scale(p2)
        s1[is.na(s1)] <- 0
        s2[is.na(s2)] <- 0
        joint_embed <- cbind(s1, s2)
      } else {
        # Unpaired multiomics: modalities stacked (distinct cell populations)
        p1 <- .fast_mod_pca(m1, n_pcs = n_pcs)
        p2 <- .fast_mod_pca(m2, n_pcs = n_pcs)
        p1_cached <- p1
        p2_cached <- p2
        min_p <- min(ncol(p1), ncol(p2))
        joint_embed <- rbind(scale(p1[, seq_len(min_p), drop = FALSE]),
                             scale(p2[, seq_len(min_p), drop = FALSE]))
        joint_embed[is.na(joint_embed)] <- 0
        bt_vec <- c(rep("Modality_1", n1), rep("Modality_2", n2))
      }
    } else if (is.list(sim_multi) && length(sim_multi) == 1L) {
      joint_embed <- .fast_mod_pca(sim_multi[[1L]], n_pcs = n_pcs)
    }
  } else if (detect_paired && !is.na(n1) && nrow(joint_embed) == n1) {
    is_paired_layout <- TRUE
  } else if (!detect_paired && !is.na(n1) && !is.na(n2) && nrow(joint_embed) == (n1 + n2)) {
    is_paired_layout <- FALSE
    bt_vec <- c(rep("Modality_1", n1), rep("Modality_2", n2))
  }

  if (is.null(joint_embed) || nrow(joint_embed) < 3L) {
    warning("Unable to construct valid joint embedding for multiomics integration evaluation.")
    return(NULL)
  }

  n_cells <- nrow(joint_embed)

  # 2. Always resolve cell_types and batch_info cleanly (including list(ct_m1, ct_m2)),
  #    whether joint_embed was precomputed or freshly constructed!
  .resolve_multi_meta <- function(meta_obj, n_target, paired_flag, n1_v, n2_v) {
    if (is.null(meta_obj)) return(NULL)
    if (is.data.frame(meta_obj)) {
      if (nrow(meta_obj) == n_target) return(as.character(meta_obj[[1L]]))
      return(NULL)
    }
    if (is.list(meta_obj)) {
      if (paired_flag && length(meta_obj) >= 1L && !is.null(meta_obj[[1L]]) && length(meta_obj[[1L]]) == n_target) {
        return(as.character(meta_obj[[1L]]))
      }
      if (paired_flag && length(meta_obj) >= 2L && !is.null(meta_obj[[2L]]) && length(meta_obj[[2L]]) == n_target) {
        return(as.character(meta_obj[[2L]]))
      }
      if (!paired_flag && length(meta_obj) >= 2L && !is.null(meta_obj[[1L]]) && !is.null(meta_obj[[2L]])) {
        comb_v <- c(as.character(meta_obj[[1L]]), as.character(meta_obj[[2L]]))
        if (length(comb_v) == n_target) return(comb_v)
      }
      # Fallback: first element matching target length
      for (el in meta_obj) {
        if (!is.null(el) && length(el) == n_target) return(as.character(el))
      }
      return(NULL)
    }
    if (length(meta_obj) == n_target) return(as.character(meta_obj))
    if (!paired_flag && !is.na(n1_v) && !is.na(n2_v) && length(meta_obj) == n1_v && n_target == (n1_v + n2_v)) {
      return(c(as.character(meta_obj), rep("Unannotated", n2_v)))
    }
    NULL
  }

  ct_vec <- .resolve_multi_meta(cell_types, n_cells, is_paired_layout, n1, n2)
  if (is.null(bt_vec)) {
    bt_vec <- .resolve_multi_meta(batch_info, n_cells, is_paired_layout, n1, n2)
  }

  has_ground_truth_ct <- !is.null(ct_vec) && length(unique(stats::na.omit(ct_vec))) >= 2L

  if (!is.null(bt_vec)) {
    bt_vec <- as.character(bt_vec)
    if (length(bt_vec) != n_cells) {
      bt_vec <- if (length(bt_vec) > n_cells) bt_vec[seq_len(n_cells)] else c(bt_vec, rep("Modality_2", n_cells - length(bt_vec)))
    }
  }

  # If batch is NULL or unimodal, evaluate cross-modal alignment by projecting modalities
  eval_modality_mixing <- is.null(bt_vec) || length(unique(bt_vec)) < 2L

  # For batch_asw / cell_type_asw / graph_connectivity when ct_vec is truly absent,
  # derive modality-1 clusters (not joint_embed clusters) to avoid circularity.
  ct_for_structure <- ct_vec
  if (is.null(ct_for_structure)) {
    k_fallback <- min(5L, max(2L, as.integer(floor(sqrt(n_cells / 5)))))
    if (!is.null(p1_cached) && nrow(p1_cached) == n_cells) {
      cl_fb <- .safe_kmeans_labels(p1_cached, k_fallback, k_fallback)
      ct_for_structure <- if (!is.null(cl_fb)) paste0("Cluster_", cl_fb) else rep("Cluster_1", n_cells)
    } else {
      cl_fb <- .safe_kmeans_labels(joint_embed, k_fallback, k_fallback)
      ct_for_structure <- if (!is.null(cl_fb)) paste0("Cluster_", cl_fb) else rep("Cluster_1", n_cells)
    }
  }

  if (eval_modality_mixing && is.list(sim_multi) && length(sim_multi) >= 2L) {
    m1 <- sim_multi[[1L]]
    m2 <- sim_multi[[2L]]
    p1 <- if (!is.null(p1_cached)) p1_cached else .fast_mod_pca(m1, n_pcs = n_pcs)
    p2 <- if (!is.null(p2_cached)) p2_cached else .fast_mod_pca(m2, n_pcs = n_pcs)
    p1_cached <- p1
    p2_cached <- p2
    min_p <- min(ncol(p1), ncol(p2))
    cross_embed <- rbind(scale(p1[, seq_len(min_p), drop = FALSE]),
                         scale(p2[, seq_len(min_p), drop = FALSE]))
    cross_embed[is.na(cross_embed)] <- 0
    cross_batch <- c(rep("Modality_1", n1), rep("Modality_2", n2))
    cross_ct <- if (length(ct_for_structure) == (n1 + n2)) {
      ct_for_structure
    } else if (length(ct_for_structure) == n1 && n1 == n2) {
      c(ct_for_structure, ct_for_structure)
    } else {
      c(ct_for_structure[seq_len(min(length(ct_for_structure), n1))], rep("Modality_2", n2))
    }
    b_asw <- calc_batch_asw(cross_embed, cross_batch, cross_ct, threads = threads)
  } else {
    bt_target <- if (!is.null(bt_vec)) bt_vec else rep("Batch_1", n_cells)
    b_asw <- calc_batch_asw(joint_embed, bt_target, ct_for_structure, threads = threads)
  }

  ct_asw <- calc_cell_type_asw(joint_embed, ct_for_structure, threads = threads)
  g_conn <- calc_graph_connectivity(joint_embed, ct_for_structure, k = min(15L, n_cells - 1L))

  # Compute ARI and NMI cleanly without circular self-clustering:
  # - When ground-truth cell types exist, cluster joint_embed once and evaluate both ARI and NMI.
  # - When cell types are absent in paired multiomics, measure cross-modal cluster concordance
  #   between modality-1 PCA clusters and modality-2 PCA clusters (never joint vs joint!).
  if (has_ground_truth_ct) {
    n_types <- length(unique(ct_vec))
    pred_cl <- .safe_kmeans_labels(joint_embed, NULL, n_types)
    if (!is.null(pred_cl)) {
      ari_val <- max(0, min(1, calc_ari(pred_cl, ct_vec)))
      nmi_val <- max(0, min(1, calc_nmi(pred_cl, ct_vec)))
    } else {
      ari_val <- NA_real_
      nmi_val <- NA_real_
    }
  } else if (is_paired_layout && is.list(sim_multi) && length(sim_multi) >= 2L) {
    p1 <- if (!is.null(p1_cached)) p1_cached else .fast_mod_pca(sim_multi[[1L]], n_pcs = n_pcs)
    p2 <- if (!is.null(p2_cached)) p2_cached else .fast_mod_pca(sim_multi[[2L]], n_pcs = n_pcs)
    k_xm <- min(5L, max(2L, as.integer(floor(sqrt(n_cells / 5)))))
    cl1 <- .safe_kmeans_labels(p1, k_xm, k_xm)
    cl2 <- .safe_kmeans_labels(p2, k_xm, k_xm)
    if (!is.null(cl1) && !is.null(cl2)) {
      ari_val <- max(0, min(1, calc_ari(cl1, cl2)))
      nmi_val <- max(0, min(1, calc_nmi(cl1, cl2)))
    } else {
      ari_val <- NA_real_
      nmi_val <- NA_real_
    }
  } else {
    ari_val <- NA_real_
    nmi_val <- NA_real_
  }
  
  metrics_list <- list(
    batch_asw = as.numeric(b_asw),
    cell_type_asw = as.numeric(ct_asw),
    graph_connectivity = as.numeric(g_conn),
    ari = as.numeric(ari_val),
    nmi = as.numeric(nmi_val)
  )
  
  summary_df <- data.frame(
    Category = "(VI) Multiomics Integration & Alignment",
    Property = "multiomics_integration",
    Metric = c("batch_asw", "cell_type_asw", "graph_connectivity", "ari", "nmi"),
    Value = as.numeric(c(b_asw, ct_asw, g_conn, ari_val, nmi_val)),
    Modality = "Joint (Cross-Modal)",
    stringsAsFactors = FALSE
  )
  
  res <- list(
    metrics = metrics_list,
    integration_summary_table = summary_df,
    embedding = joint_embed
  )
  class(res) <- c("scSimEval_integration", "list")
  res
}
