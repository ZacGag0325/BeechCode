# 00_master_pipeline.R
############################################################
# (Full file unchanged except the function below)
# Replace the existing make_mll_from_bruvo() helper with this version.
############################################################

### Helper: build deterministic MLL from Bruvo distances + threshold ----
make_mll_from_bruvo <- function(genind_obj, threshold) {
  # Explicit repeat lengths to avoid poppr automatic estimation warning.
  # For this microsatellite panel, all loci are treated as dinucleotide repeats.
  replen_vec <- rep(2, adegenet::nLoc(genind_obj))
  names(replen_vec) <- adegenet::locNames(genind_obj)
  
  bruvo <- poppr::bruvo.dist(genind_obj, replen = replen_vec)
  bruvo_mat <- as.matrix(bruvo)
  
  if (!is.matrix(bruvo_mat) || nrow(bruvo_mat) != ncol(bruvo_mat)) {
    stop("Invalid Bruvo distance matrix: expected square matrix.")
  }
  if (nrow(bruvo_mat) != adegenet::nInd(genind_obj)) {
    stop("Invalid Bruvo distance matrix: nrow does not match nInd(genind).")
  }
  if (any(is.na(bruvo_mat))) {
    stop("Invalid Bruvo distance matrix: contains NA values.")
  }
  if (any(!is.finite(bruvo_mat))) {
    stop("Invalid Bruvo distance matrix: contains non-finite values.")
  }
  
  rn <- rownames(bruvo_mat)
  cn <- colnames(bruvo_mat)
  ids <- adegenet::indNames(genind_obj)
  
  if (is.null(rn) || is.null(cn) || !identical(rn, cn)) {
    stop("Bruvo matrix row/column names are missing or inconsistent.")
  }
  if (!setequal(rn, ids)) {
    stop("Bruvo matrix sample names do not match indNames(genind).")
  }
  
  bruvo_mat <- bruvo_mat[ids, ids, drop = FALSE]
  bruvo_mat <- (bruvo_mat + t(bruvo_mat)) / 2
  diag(bruvo_mat) <- 0
  
  hc <- stats::hclust(stats::as.dist(bruvo_mat), method = "average")
  raw_id <- stats::cutree(hc, h = threshold)
  stable_id <- relabel_mll_ids(raw_id, sample_ids = ids)
  
  list(
    bruvo_mat = bruvo_mat,
    mll_id = stable_id,
    mll_label = paste0("MLL", stable_id)
  )
}