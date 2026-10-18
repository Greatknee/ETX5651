# Quantile MAR connectedness: heatmap and quantile asymmetry.
# Chen (2021) SVD projection of unrestricted quantile VAR slopes.
# Quantiles: 0.05 and 0.50. Heatmap horizon: 8 quarters.
# Quantile asymmetry: horizons 1 to 20 quarters; lower minus median shares.
#
# Install once: install.packages(c("quantreg", "ggplot2", "readxl"))
# Run: Rscript layer2_qmar.R data.xlsx results_layer2
# Or source in RStudio with the input file in your working directory.

required <- c("quantreg", "ggplot2", "readxl")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(ggplot2))

args <- if (sys.nframe() == 0L) commandArgs(trailingOnly = TRUE) else character()
INPUT <- if (length(args)) args[1] else "data.xlsx"

OUT <- if (length(args) >= 2) args[2] else "results_layer2"
TAUS <- c(0.05, 0.50)
H_MAIN <- 8L
H_GRID <- sort(unique(c(1:20, H_MAIN)))
I <- 2L; K <- 2L; P <- I * K
NODES <- c("US_GDP", "AUS_GDP", "US_NFCI", "AUS_FCI")
COLS <- c("us_gdp_growth_qoq_pct", "aus_gdp_growth_qoq_pct", "us_nfci", "aus_fci")
tau_key <- function(tau) sprintf("%.2f", tau)

# Column stacking: vec(Y_t) = (US_GDP, AUS_GDP, US_NFCI, AUS_FCI)'.
read_data <- function(path) {
  if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl for XLSX input.")
    d <- as.data.frame(readxl::read_excel(path, sheet = "Data", .name_repair = "minimal"))
  } else d <- utils::read.csv(path, check.names = FALSE)
  if (!all(c("quarter", COLS) %in% names(d))) stop("Input columns do not match the supplied data.")
  if (any(!grepl("^[0-9]{4}Q[1-4]$", d$quarter))) stop("Invalid quarter labels.")
  d$quarter_id <- 4L * as.integer(substr(d$quarter, 1, 4)) + as.integer(substr(d$quarter, 6, 6)) - 1L
  d <- d[order(d$quarter_id), ]
  if (anyDuplicated(d$quarter_id) || any(diff(d$quarter_id) != 1)) stop("Duplicate or missing quarters.")
  if (nrow(d) <= P + 1L) stop("Too few observations.")
  if (!all(vapply(d[COLS], is.numeric, logical(1))) || any(!is.finite(as.matrix(d[COLS])))) {
    stop("Missing or nonnumeric observations.")
  }
  d
}

# G(B %x% A) = vec(A) vec(B)'.
rearrange <- function(Phi) {
  G <- matrix(0, I^2, K^2)
  for (s in seq_len(K)) for (r in seq_len(K)) {
    rr <- (r - 1L) * I + seq_len(I)
    ss <- (s - 1L) * I + seq_len(I)
    G[, r + (s - 1L) * K] <- as.vector(Phi[rr, ss, drop = FALSE])
  }
  G
}

project_kronecker <- function(Phi) {
  z <- svd(rearrange(Phi))
  if (z$d[1] <= .Machine$double.eps) stop("Degenerate zero slope matrix.")
  A <- matrix(z$u[, 1], I, I)
  B <- matrix(z$d[1] * z$v[, 1], K, K)
  first <- which(abs(as.vector(A)) > 1e-12)[1]
  if (as.vector(A)[first] < 0) { A <- -A; B <- -B }
  Phi <- kronecker(B, A)
  dimnames(Phi) <- list(NODES, NODES)
  Phi
}

fit_one <- function(L, R, tau) {
  X <- cbind(Intercept = 1, L)
  if (qr(X)$rank != ncol(X)) stop("Rank-deficient design.")
  coefficients <- vapply(seq_len(P), function(q) {
    quantreg::rq.fit.br(X, R[, q], tau = tau)$coefficients
  }, numeric(P + 1L))
  Phi_U <- t(coefficients[-1, , drop = FALSE])
  Phi <- project_kronecker(Phi_U)
  offset <- R - L %*% t(Phi)
  intercept <- apply(offset, 2, stats::quantile, probs = tau, type = 1, names = FALSE)
  residuals <- sweep(offset, 2, intercept, "-")
  centered <- sweep(residuals, 2, colMeans(residuals), "-")
  Omega <- crossprod(centered) / nrow(centered)
  dimnames(Omega) <- list(NODES, NODES)
  if (min(eigen(Omega, symmetric = TRUE, only.values = TRUE)$values) <= 1e-12) {
    stop("Residual covariance is not positive definite.")
  }
  rho <- max(Mod(eigen(Phi, only.values = TRUE)$values))
  if (rho >= 1) warning("Unstable projected matrix at tau = ", tau)
  list(tau = tau, Phi = Phi, Omega = Omega)
}

# Generalised FEVD: recipient rows, source columns; r = 0,...,H-1.
gfevd <- function(Phi, Omega, H) {
  stopifnot(H >= 1, H == as.integer(H), max(abs(Omega - t(Omega))) < 1e-8)
  C <- diag(P)
  numerator <- matrix(0, P, P)
  denominator <- numeric(P)
  for (r in 0:(H - 1L)) {
    Q <- C %*% Omega
    numerator <- numerator + sweep(Q^2, 2, diag(Omega), "/")
    denominator <- denominator + diag(Q %*% t(C))
    C <- C %*% Phi
  }
  raw <- sweep(numerator, 1, denominator, "/")
  theta <- sweep(raw, 1, rowSums(raw), "/")
  dimnames(theta) <- list(NODES, NODES)
  stopifnot(max(abs(rowSums(theta) - 1)) < 1e-10, min(theta) > -1e-12)
  theta
}

matrix_long <- function(M) {
  d <- expand.grid(recipient = rownames(M), source = colnames(M), KEEP.OUT.ATTRS = FALSE)
  d$share_pct <- 100 * as.vector(M)
  d
}

d <- read_data(INPUT)
Y <- as.matrix(d[COLS]); colnames(Y) <- NODES
L <- Y[-nrow(Y), , drop = FALSE]
R <- Y[-1, , drop = FALSE]
fits <- lapply(TAUS, function(tau) fit_one(L, R, tau))

heat <- do.call(rbind, lapply(fits, function(f) {
  data.frame(tau = tau_key(f$tau), matrix_long(gfevd(f$Phi, f$Omega, H_MAIN)))
}))
heat$source <- factor(heat$source, levels = NODES)
heat$recipient <- factor(heat$recipient, levels = rev(NODES))
heat$tau <- factor(heat$tau, levels = c("0.05", "0.50"))

theme_set(theme_minimal(base_size = 12) + theme(
  panel.grid.minor = element_blank(),
  plot.title = element_text(face = "bold"),
  legend.position = "bottom"
))

p_heat <- ggplot(heat, aes(source, recipient, fill = share_pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", share_pct)), size = 3.4) +
  facet_wrap(~tau, nrow = 1, labeller = label_both) +
  scale_fill_gradient(low = "#F1F5F7", high = "#79AEC4", limits = c(0, 100)) +
  labs(x = "Source innovation", y = "Recipient", fill = "Share (%)",
       title = paste("Generalised connectedness at H =", H_MAIN, "quarters")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dir.create(file.path(OUT, "figures"), recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(OUT, "figures", "connectedness_heatmap.pdf"), p_heat,
       width = 12, height = 7, device = grDevices::cairo_pdf, bg = "white")


# Financial-to-growth share: GDP recipient rows, both financial source columns.
# Reuse the fitted models and generalised FEVD used by the heatmap.
C_results <- do.call(rbind, lapply(H_GRID, function(H) {
  do.call(rbind, lapply(fits, function(f) {
    theta <- gfevd(f$Phi, f$Omega, H)
    C <- rowSums(theta[c("US_GDP", "AUS_GDP"), c("US_NFCI", "AUS_FCI"), drop = FALSE])
    data.frame(country = c("US", "Australia"), tau = f$tau, H = H,
               C = unname(C))
  }))
}))

# Match country and horizon before computing D(H) = C_0.05(H) - C_0.50(H).
# C and D are fractions; D_pp is measured in percentage points.
lower <- C_results[C_results$tau == 0.05, c("country", "H", "C")]
median <- C_results[C_results$tau == 0.50, c("country", "H", "C")]
matched <- merge(lower, median, by = c("country", "H"), suffixes = c("_05", "_50"))
D_results <- data.frame(country = matched$country, H = matched$H,
                        D = matched$C_05 - matched$C_50)
D_results$D_pp <- 100 * D_results$D

# Plot 2: D by horizon; zero is equal lower-quantile and median connectedness.
p_D <- ggplot(D_results, aes(H, D_pp, colour = country)) +
  geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c(US = "#245B83", Australia = "#B34B3E")) +
  labs(x = "Forecast horizon H (quarters)", y = "Difference (percentage points)",
       colour = NULL) +
  theme_minimal(base_size = 13) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "figures", "quantile_asymmetry.pdf"), p_D,
       width = 8, height = 4, device = grDevices::cairo_pdf, bg = "white")

if (interactive()) {
  print(p_heat)
  print(p_D)
}

# Positive D indicates a larger financial-to-growth network share under
# lower-quantile dynamics; it does not identify the sign of GDP responses.
