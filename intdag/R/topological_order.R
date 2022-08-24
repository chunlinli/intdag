
topological_order <- function(v, thresh = 0.1) {
    p <- ncol(v)
    q <- nrow(v)
    if (q < p) stop("No sufficient interventions: q < p.")

    v[abs(v) < thresh] <- 0

    an_mat <- matrix(0, p, p)
    in_mat <- matrix(0, q, p)
    iv_mat <- matrix(0, q, p)

    removed_x <- rep(FALSE, q)
    removed_y <- rep(FALSE, p)

    v_abs <- abs(v)
    v_nz <- (v_abs != 0)

    # check if there is a primary variable without intervention.
    while (!all(removed_y)) {

        # leaf-instrument pairs
        iv_targets <- rowSums(as.matrix(v_nz[, !removed_y]))
        one <- min(iv_targets[iv_targets > 0 & !removed_x])
        leaf_iv <- which(!removed_x & iv_targets == one)
        if (length(leaf_iv) == 0) break
        leaf <- rep(NA, length(leaf_iv))
        leaf_iter <- 1
        for (l in leaf_iv) {
            j <- which(v_abs[l, ] == max(v_abs[l, !removed_y]) & !removed_y)[1]
            iv_mat[l, j] <- 1
            in_mat[l, j] <- 1
            leaf[leaf_iter] <- j
            leaf_iter <- leaf_iter + 1
        }
        leaf <- unique(leaf)

        # leaf-noninstrument pairs
        for (j in leaf) {
            leaf_interventions <- which(!removed_x & v_abs[, j] != 0)
            leaf_noniv <- setdiff(leaf_interventions, leaf_iv)
            in_mat[leaf_noniv, j] <- 1
        }

        # ancestral relations
        for (j in leaf) {
            j_instrument <- which(iv_mat[, j] != 0)
            j_descendant <- vector("list", length = length(j_instrument))
            for (l in seq_len(length(j_instrument))) {
                l2 <- j_instrument[l]
                j_descendant[[l]] <- which(removed_y & v_abs[l2, ] != 0)
            }

            j_descendant <- Reduce(union, j_descendant)
            an_mat[j, j_descendant] <- 1
        }

        # peeling-off
        removed_y[leaf] <- TRUE
        removed_x[leaf_iv] <- TRUE
    }

    # reconstruction of topological order
    an_mat <- (solve(diag(p) - an_mat) != 0) - diag(p)
    in_mat <- 1 * (in_mat %*% (diag(p) + an_mat) > 0)

    list(an_mat = an_mat, in_mat = in_mat, iv_mat = iv_mat)
}
