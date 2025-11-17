bc_alpha_beta_table <- function(bgms_fit,
                                var_names = NULL,
                                z_crit = 1.96) {
  # --- basic checks ---
  main_mean <- bgms_fit$posterior_mean_main
  main_summ <- bgms_fit$posterior_summary_main
  
  if (is.null(main_mean) || is.null(main_summ)) {
    stop("Need posterior_mean_main and posterior_summary_main in bgms_fit.")
  }
  
  p <- nrow(main_mean)  # number of variables
  
  if (!all(dim(main_mean) == c(p, 2))) {
    stop("posterior_mean_main must be p x 2 (alpha, beta).")
  }
  if (nrow(main_summ) != 2 * p) {
    stop("posterior_summary_main must have 2*p rows (one alpha row and one beta row per variable).")
  }
  
  # columns in posterior_summary_main
  cn <- colnames(main_summ)
  if (!"sd" %in% cn || !"mean" %in% cn) {
    stop("posterior_summary_main must contain 'mean' and 'sd' columns.")
  }
  
  # row order: rows 1..p = alpha for vars 1..p; rows (p+1)..(2p) = beta
  alpha_rows <- 1:p
  beta_rows  <- (p + 1):(2 * p)
  
  alpha_sd <- main_summ[alpha_rows, "sd"]
  beta_sd  <- main_summ[beta_rows,  "sd"]
  
  alpha_mean <- main_mean[, 1]
  beta_mean  <- main_mean[, 2]
  
  alpha_low  <- alpha_mean - z_crit * alpha_sd
  alpha_high <- alpha_mean + z_crit * alpha_sd
  beta_low   <- beta_mean  - z_crit * beta_sd
  beta_high  <- beta_mean  + z_crit * beta_sd
  
  out <- data.frame(
    var_idx    = seq_len(p),
    alpha_mean = alpha_mean,
    alpha_low  = alpha_low,
    alpha_high = alpha_high,
    beta_mean  = beta_mean,
    beta_low   = beta_low,
    beta_high  = beta_high,
    stringsAsFactors = FALSE
  )
  
  # attach variable names
  if (!is.null(var_names)) {
    if (length(var_names) < p) {
      warning("var_names shorter than number of variables; some names will be NA.")
      tmp <- rep(NA_character_, p)
      tmp[seq_along(var_names)] <- var_names
      var_names <- tmp
    }
    out$var_name <- var_names
  } else if (!is.null(bgms_fit$x) && !is.null(colnames(bgms_fit$x))) {
    out$var_name <- colnames(bgms_fit$x)
  } else {
    out$var_name <- paste0("V", seq_len(p))
  }
  
  out
}


bc_alpha_beta_plot <- function(alpha_beta,
                               label_var = "var_name",
                               core_alpha = 1.5,
                               core_beta  = 1.0) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required.")
  }
  
  df <- alpha_beta
  
  # classify regions (you can tweak core_alpha/core_beta)
  df$core_type <- "peripheral/weak"
  df$core_type[abs(df$alpha_mean) >= core_alpha & df$beta_mean >= core_beta] <- "core"
  df$core_type[abs(df$alpha_mean) <  core_alpha & df$beta_mean >= core_beta] <- "controversial"
  
  ggplot2::ggplot(df,
                  ggplot2::aes(x = alpha_mean, y = beta_mean, color = core_type)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = beta_low, ymax = beta_high),
                           alpha = 0.4, width = 0) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = alpha_low, xmax = alpha_high),
                            alpha = 0.4, height = 0) +
    ggplot2::geom_point(size = 3) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = .data[[label_var]]),
      size = 3,
      max.overlaps = Inf
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
    ggplot2::scale_color_manual(values = c(
      "core"            = "red",
      "controversial"   = "orange",
      "peripheral/weak" = "grey40"
    )) +
    ggplot2::labs(
      x = expression(alpha ~ "(directionality)"),
      y = expression(beta ~ "(neutrality aversion)"),
      color = "Region"
    ) +
    ggplot2::xlim(-5, 5) +
    ggplot2::theme_minimal()
}


# Compute per-participant Blume–Capel-style energy from a bgms fit
#
# E_n = - sum_{i<j} J_ij s_ni s_nj + sum_i D_i s_ni^2 - sum_i h_i s_ni
#
energy_from_bgms <- function(bgms_fit,
                             data,
                             reference_category = 0,
                             na_as_zero = TRUE) {
  J    <- bgms_fit$posterior_mean_pairwise
  main <- bgms_fit$posterior_mean_main
  
  if (is.null(J) || is.null(main)) {
    stop("Need posterior_mean_pairwise and posterior_mean_main in bgms_fit.")
  }
  
  X_raw <- as.matrix(data)
  p <- ncol(X_raw)
  
  if (!all(dim(J) == c(p, p))) {
    stop("posterior_mean_pairwise must be p x p, here p = ", p, ".")
  }
  if (!all(dim(main) == c(p, 2))) {
    stop("posterior_mean_main must be p x 2 (alpha, beta) for p variables.")
  }
  
  # symmetrize J and zero diagonal
  J <- (J + t(J)) / 2
  diag(J) <- 0
  
  # main: col1 = alpha (h), col2 = beta (D)
  h <- main[, 1]
  D <- main[, 2]
  
  # recode data to spins {-1,0,+1}
  S <- X_raw
  S[,] <- NA_real_
  S[X_raw <  reference_category] <- -1
  S[X_raw == reference_category] <-  0
  S[X_raw >  reference_category] <-  1
  if (na_as_zero) S[is.na(S)] <- 0
  
  S_mat <- as.matrix(S)
  
  # pairwise term
  SJ <- S_mat %*% J
  pair_term <- -0.5 * rowSums(SJ * S_mat)
  
  # neutrality term
  s_sq <- S_mat^2
  neut_term  <- rowSums(sweep(s_sq, 2, D, `*`))
  
  # field term
  field_term <- -rowSums(sweep(S_mat, 2, h, `*`))
  
  E <- pair_term + neut_term + field_term
  
  if (!na_as_zero) {
    row_has_na <- apply(is.na(X_raw), 1, any)
    E[row_has_na] <- NA_real_
  }
  if (!is.null(rownames(X_raw))) names(E) <- rownames(X_raw)
  
  E
}

