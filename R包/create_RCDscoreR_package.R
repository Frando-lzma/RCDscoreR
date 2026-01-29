## ============================================================
## Create a local R package: RCDscoreR
## ============================================================

## ---- 0) packages ----
pkgs <- c("usethis", "devtools", "roxygen2")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(usethis)

## ---- 1) where to create the package ----
# 你可以改成你想放包工程的目录（不要放在 system32 之类路径）
pkg_parent <- "D:/user/Watermelon Sword/(important)biomessage/RCD库开发"
pkg_path <- file.path(pkg_parent, "RCDscoreR")

if (dir.exists(pkg_path)) stop("Package folder already exists: ", pkg_path)

## ---- 2) create package skeleton ----
usethis::create_package(pkg_path, open = FALSE)

## ---- 3) add dependencies ----
usethis::use_package("matrixStats")   # 用于 colSds
usethis::use_package("data.table")    # 可选：如果你后续要放导入 FerrDb 的函数
usethis::use_package("stringr")       # 可选

## ---- 4) add R code file ----
rcd_file <- file.path(pkg_path, "R", "scoring.R")

writeLines(c(
  '`%||%` <- function(x, y) if (is.null(x)) y else x',
  '',
  '#\' Compute housekeeping baseline for each sample/cell',
  '#\'',
  '#\' @param X Numeric matrix, genes x samples/cells. Row names are gene symbols.',
  '#\' @param hk_genes Character vector of housekeeping genes.',
  '#\' @param eps Small constant to avoid zero sd.',
  '#\' @return A list with mu, sd, and housekeeping genes used.',
  '#\' @export',
  'hk_baseline <- function(X, hk_genes, eps = 1e-6) {',
  '  stopifnot(is.matrix(X) || is.data.frame(X))',
  '  X <- as.matrix(X)',
  '  if (is.null(rownames(X))) stop("X must have rownames as gene symbols.")',
  '  hk_genes <- intersect(rownames(X), hk_genes)',
  '  if (length(hk_genes) < 10) stop("Too few housekeeping genes found in X (<10).")',
  '  HK <- X[hk_genes, , drop = FALSE]',
  '  mu <- colMeans(HK, na.rm = TRUE)',
  '  sd <- matrixStats::colSds(as.matrix(HK), na.rm = TRUE)',
  '  sd <- pmax(sd, eps)',
  '  list(mu = mu, sd = sd, hk_genes_used = hk_genes)',
  '}',
  '',
  '#\' Convert expression to positive activity using housekeeping-referenced z-score',
  '#\'',
  '#\' @param X Numeric matrix, genes x samples/cells (log-scale recommended).',
  '#\' @param baseline Output of hk_baseline().',
  '#\' @param z0 Z-score threshold; only (z - z0) > 0 contributes.',
  '#\' @param eps Small constant.',
  '#\' @return A list with Z (z-score matrix) and A (activity matrix).',
  '#\' @export',
  'gene_activity <- function(X, baseline, z0 = 1, eps = 1e-6) {',
  '  X <- as.matrix(X)',
  '  mu <- baseline$mu',
  '  sd <- baseline$sd',
  '  Z <- sweep(X, 2, mu, FUN = "-")',
  '  Z <- sweep(Z, 2, sd + eps, FUN = "/")',
  '  A <- pmax(0, Z - z0)',
  '  list(Z = Z, A = A)',
  '}',
  '',
  '#\' Score one gene set using activity matrix',
  '#\'',
  '#\' @param A Activity matrix (genes x samples/cells).',
  '#\' @param genes Character vector of genes.',
  '#\' @param weights Named numeric vector of gene weights (optional). Names must be gene symbols.',
  '#\' @return Numeric vector (length = ncol(A)).',
  '#\' @export',
  'score_one_set <- function(A, genes, weights = NULL) {',
  '  genes <- intersect(rownames(A), genes)',
  '  if (length(genes) == 0) return(rep(NA_real_, ncol(A)))',
  '  M <- A[genes, , drop = FALSE]',
  '  if (is.null(weights)) {',
  '    return(colMeans(M, na.rm = TRUE))',
  '  } else {',
  '    w <- weights[genes]',
  '    w[is.na(w)] <- 1',
  '    num <- colSums(sweep(M, 1, w, `*`), na.rm = TRUE)',
  '    den <- sum(w, na.rm = TRUE)',
  '    return(num / den)',
  '  }',
  '}',
  '',
  '#\' RCD scoring (Driver / Suppressor / Marker + Net + Conflict)',
  '#\'',
  '#\' This function computes per-sample/per-cell scores for each RCD using',
  '#\' housekeeping-referenced z-score activity.',
  '#\'',
  '#\' @param X Numeric matrix, genes x samples/cells. Row names are gene symbols.',
  '#\'   Recommended: log-scaled (e.g., log1p normalized for scRNA, log2(TPM+1) for bulk).',
  '#\' @param genesets A list with components: driver, suppressor, marker. Each is a named list: RCD -> genes.',
  '#\'   Optionally include conflict (RCD -> overlap genes), but conflict is computed from driver/suppressor anyway.',
  '#\' @param hk_genes Housekeeping genes vector.',
  '#\' @param z0 Z-score threshold for "significant expression". Default 1.',
  '#\' @param weights_driver Optional named numeric vector: gene weights for driver sets.',
  '#\' @param weights_suppressor Optional named numeric vector: gene weights for suppressor sets.',
  '#\' @param weights_marker Optional named numeric vector: gene weights for marker sets.',
  '#\' @param conflict_frac_cutoff Structural conflict cutoff (Jaccard overlap) for flagging.',
  '#\' @return A list with matrices: driver_score, suppressor_score, marker_score, net_tendency,',
  '#\' and vectors: conflict_frac, conflict_flag, hk_used.',
  '#\' @export',
  'rcd_score <- function(X, genesets, hk_genes, z0 = 1,',
  '                      weights_driver = NULL,',
  '                      weights_suppressor = NULL,',
  '                      weights_marker = NULL,',
  '                      conflict_frac_cutoff = 0.05) {',
  '  stopifnot(is.matrix(X) || is.data.frame(X))',
  '  X <- as.matrix(X)',
  '  if (is.null(rownames(X))) stop("X must have rownames as gene symbols.")',
  '  if (is.null(colnames(X))) colnames(X) <- paste0("S", seq_len(ncol(X)))',
  '',
  '  base <- hk_baseline(X, hk_genes = hk_genes)',
  '  act  <- gene_activity(X, baseline = base, z0 = z0)',
  '  A <- act$A',
  '',
  '  all_rcd <- sort(unique(c(names(genesets$driver), names(genesets$suppressor), names(genesets$marker))))',
  '',
  '  drv <- sup <- mkr <- matrix(NA_real_, nrow = length(all_rcd), ncol = ncol(X),',
  '                              dimnames = list(all_rcd, colnames(X)))',
  '',
  '  for (r in all_rcd) {',
  '    drv[r, ] <- score_one_set(A, genesets$driver[[r]] %||% character(0), weights_driver)',
  '    sup[r, ] <- score_one_set(A, genesets$suppressor[[r]] %||% character(0), weights_suppressor)',
  '    mkr[r, ] <- score_one_set(A, genesets$marker[[r]] %||% character(0), weights_marker)',
  '  }',
  '',
  '  net <- drv - sup',
  '',
  '  conflict_frac <- sapply(all_rcd, function(r) {',
  '    g1 <- genesets$driver[[r]] %||% character(0)',
  '    g2 <- genesets$suppressor[[r]] %||% character(0)',
  '    u <- union(g1, g2)',
  '    if (length(u) == 0) return(0)',
  '    length(intersect(g1, g2)) / length(u)',
  '  })',
  '  conflict_flag <- conflict_frac >= conflict_frac_cutoff',
  '',
  '  list(',
  '    driver_score = drv,',
  '    suppressor_score = sup,',
  '    marker_score = mkr,',
  '    net_tendency = net,',
  '    conflict_frac = conflict_frac,',
  '    conflict_flag = conflict_flag,',
  '    hk_used = base$hk_genes_used,',
  '    z0 = z0',
  '  )',
  '}',
  '',
  '#\' Rank top-k RCDs for each sample/cell using a chosen score matrix',
  '#\'',
  '#\' @param score_mat RCD x samples matrix (e.g., res$net_tendency).',
  '#\' @param k Top-k.',
  '#\' @return A data.frame with sample, rank, RCD, score.',
  '#\' @export',
  'rcd_rank_topk <- function(score_mat, k = 3) {',
  '  score_mat <- as.matrix(score_mat)',
  '  if (is.null(colnames(score_mat))) colnames(score_mat) <- paste0("S", seq_len(ncol(score_mat)))',
  '  out <- lapply(seq_len(ncol(score_mat)), function(j) {',
  '    s <- score_mat[, j]',
  '    ord <- order(s, decreasing = TRUE, na.last = NA)',
  '    ord <- head(ord, k)',
  '    data.frame(sample = colnames(score_mat)[j], rank = seq_along(ord),',
  '               RCD = rownames(score_mat)[ord], score = s[ord], row.names = NULL)',
  '  })',
  '  do.call(rbind, out)',
  '}',
  '',
  '#\' Load FerrDb genesets from an .rds file (recommended format)',
  '#\'',
  '#\' @param rds_path Path to FerrDbV3_RCD_genesets.rds you previously saved.',
  '#\' @return genesets list.',
  '#\' @export',
  'load_ferrdb_genesets <- function(rds_path) {',
  '  obj <- readRDS(rds_path)',
  '  if (!is.list(obj) || is.null(obj$genesets)) stop("Invalid rds: expecting a list with $genesets.")',
  '  obj$genesets',
  '}'
), con = rcd_file)

## ---- 5) add README (optional) ----
readme_path <- file.path(pkg_path, "README.md")
writeLines(c(
  "# RCDscoreR",
  "",
  "RCD scoring utilities (Driver/Suppressor/Marker) using housekeeping-referenced z-score activity.",
  "",
  "## Install (local)",
  "```r",
  "devtools::install_local('D:/user/Watermelon Sword/(important)biomessage/RCD库开发/RCDscoreR')",
  "```",
  "",
  "## Main functions",
  "- `load_ferrdb_genesets()`",
  "- `rcd_score()`",
  "- `rcd_rank_topk()`",
  ""
), con = readme_path)

## ---- 6) document + install ----
setwd(pkg_path)
devtools::document()
devtools::install(pkg_path, upgrade = "never")

message("Done. Package created and installed: RCDscoreR")
